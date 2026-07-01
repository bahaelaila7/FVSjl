# =============================================================================
# keyword_dispatch.jl — process a keyword stream into stand state
#
# Ported from: base/initre.f  (the giant computed-GOTO keyword processor).
#
# INITRE reads keyword records and dispatches each to a handler that mutates the
# stand. We replace the 147-way computed GOTO with named handler functions keyed
# on the keyword. A few keywords pull additional raw lines from the stream
# (TREEFMT reads 2, STDIDENT reads 1); those handlers take the `KeywordReader`.
#
# Scope note: this implements the keywords exercised by the SN test set; others
# are recognized as no-ops (KNOWN_NOOP) so processing doesn't choke. Handlers are
# added here as later chunks need them (thinning, fire, econ, ...).
# =============================================================================

# Variant-AGNOSTIC no-ops: keywords inert in EVERY variant — pure I/O / echo / debug /
# report control, plus report-table requests (no `.sum` effect anywhere; they request a
# C6 DBS/.out table FVSjl doesn't emit yet). Acting on these is a reporting chunk, not a
# simulation change. Variant-SPECIFIC no-op determinations live with each variant (see
# `variant_noop_keywords`) so they get re-verified per variant — the MANAGED bug (a real SN
# growth effect once mislabeled inert) is why that split matters.
const KNOWN_NOOP = Set([
    # pure I/O / echo / debug / report control
    "SCREEN", "NOSCREEN", "STATS", "ECHOSUM", "ECHO", "NOECHO", "NOSUM",
    "NODEBUG", "DEBUG", "CALBSTAT", "REWIND", "ENDFILE", "FVSSTAND",
    # report-table requests (output not yet emitted)
    "TREELIST", "ATRTLIST", "CUTLIST",
    # CCADJ (sstage.f act 444): adjusts CCCOEF/CCCOEF2 used ONLY inside SSTAGE (the structural-
    # stage CLASSIFICATION, Stage et al.) — verified .sum-inert (FVSjl doesn't emit SSTAGE; the
    # coefficient never reaches DGF/CCF growth). Variant-agnostic (sstage is base). Recognize when
    # SSTAGE output is in scope.
    "CCADJ",
    # bare-stand flag (also handled explicitly below; harmless here)
    "NOTREES",
])

"""
    variant_noop_keywords(variant) -> Set{String}

Keywords that are inert *for this variant* but NOT guaranteed inert for others — so each
new variant supplies its own set (default empty). The dispatch treats these as no-ops
alongside the variant-agnostic [`KNOWN_NOOP`]. Re-verify by tracing the keyword's effect
variable when porting another variant.
"""
const EMPTY_STRING_SET = Set{String}()
variant_noop_keywords(::AbstractVariant) = EMPTY_STRING_SET

"Read one raw (un-lexed) line from the keyword stream, advancing the record count."
function read_raw_line!(kr::KeywordReader)
    kr.record_count += 1
    return readline(kr.io)
end

# --- individual keyword handlers (each ported from its initre.f label) --------

# OPTION 52 — CUTEFF (initre.f:5400): the default proportion of selected trees removed/affected,
# used where THINPRSC/THINAUTO cuteff or HTGSTOP/TOPKILL PRB is left blank (replaces the 1.0 default).
function kw_cuteff!(s::StandState, rec::KeywordRecord)
    rec.present[1] && (s.control.cut_eff = Float32(rec.values[1]))
    return
end

# OPTION 12 — TFIXAREA (initre.f:816): total fixed plot area. NOTRE then expands the small-tree
# (DBH < BRK) sample by FP = 1/TFPA instead of the default fixed_plot_inv/π (notre.f:45).
function kw_tfixarea!(s::StandState, rec::KeywordRecord)
    rec.present[1] && (s.plot.total_fixed_plot = Float32(rec.values[1]))
    return
end

# OPTION 10 — DESIGN (initre.f:743): plot design.
function kw_design!(s::StandState, rec::KeywordRecord)
    p, v = s.plot, rec.values
    rec.present[1] && (p.baf = v[1])
    rec.present[2] && (p.fixed_plot_inv = v[2])
    rec.present[3] && (p.min_dbh_var_plot = v[3])
    rec.present[4] && (p.points_inv = nint(v[4]))
    rec.present[5] && (p.nonstockable = nint(v[5]))
    rec.present[6] && (p.sample_weight = v[6])
    (p.sample_weight <= 0f0 && p.points_inv > 0) && (p.sample_weight = Float32(p.points_inv))
    g = v[7]
    (g > 1f0 && g <= 100f0) && (g *= 0.01f0)
    (g > 0f0 && g <= 1f0) && (p.gross_space = g)
    return
end

# OPTION 11 — NUMCYCLE (initre.f:771).
function kw_numcycle!(s::StandState, rec::KeywordRecord)
    n = rec.values[1]
    (n >= 1f0 && nint(n) <= MAXCYC) && (s.control.ncycle = nint(n))
    return
end

# OPTION 2 — TIMEINT (vbase/initre.f:1200): cycle length (IY). Field 1 = cycle index (0/absent =
# all cycles), field 2 = period length in years (default 10). The uniform path sets s.control.year
# (YR/IFINT), the cycle length threaded through the growth models — DDS/HTG scale by FINT/5,
# autcor/year-age by it. A per-cycle index N stores the override in cycle_lengths[N+1] (Fortran
# `I=IABS(ARRAY(1))+1; IY(I)=I2`); build_cycle_schedule! cumulates it into the boundary years.
function kw_timeint!(s::StandState, rec::KeywordRecord)
    v = rec.values; pr = rec.present
    len = pr[2] ? Float32(v[2]) : 10f0
    cyc = pr[1] ? abs(nint(v[1])) : Int32(0)
    if cyc <= 0
        s.control.year = len                                      # all cycles (uniform period)
    else
        idx = Int(cyc) + 1
        idx <= MAXCY1 && (s.control.cycle_lengths[idx] = nint(len))  # cycle N → IY(N+1) = len
    end
    return
end

# OPTION 134 — CYCLEAT (vbase/initre.f:13400): request an extra cycle boundary at a calendar year.
# Collects de-duplicated positive years into control.cycleat_years (Fortran IWORK1); fvs.f /
# build_cycle_schedule! later inserts each as a new boundary strictly inside the run (without
# extending the end or moving the start), increasing the effective cycle count.
function kw_cycleat!(s::StandState, rec::KeywordRecord)
    rec.present[1] || return
    yr = nint(rec.values[1])
    yr <= 0 && return
    yr in s.control.cycleat_years || push!(s.control.cycleat_years, yr)
    return
end

# STRCLASS (ksstag.f): activate the stand structural-stage classification (SSTAGE) + optionally
# override its thresholds. Field 1 = print code (0 = compute but don't print, else print — we always
# compute; the per-cycle `.out` report is Chunk D-report, deferred). Fields 2-7 override
# GAPPCT/SSDBH/SAWDBH/CCMIN/TPAMIN/PCTSMX. The class is then available as the event-monitor
# variables BSCLASS/ASCLASS (structural class), BSTRDBH/ASTRDBH (uppermost-stratum DBH),
# BCANCOV/ACANCOV (canopy cover).
function kw_strclass!(s::StandState, rec::KeywordRecord)
    s.control.strclass_on = true
    v = rec.values; pr = rec.present
    d = s.control.strclass_thresh                       # (gappct, ssdbh, sawdbh, ccmin, tpamin, pctsmx)
    s.control.strclass_thresh = (
        pr[2] ? Float32(v[2]) : d[1], pr[3] ? Float32(v[3]) : d[2], pr[4] ? Float32(v[4]) : d[3],
        pr[5] ? Float32(v[5]) : d[4], pr[6] ? Float32(v[6]) : d[5], pr[7] ? Float32(v[7]) : d[6])
    return
end

# NOHTDREG (sn keyword option 60) — HT-DBH REGRESSION CALIBRATION control (LHTDRG), NOT an establishment flag.
# Per sn keyword option-60 handling: field 2 > 0 INVOKES the per-species height-diameter calibration
# (LHTDRG=.TRUE.); blank/zero SUPPRESSES it (= the SN default, grinit.f:104 LHTDRG=.FALSE.). jl models only the
# default/suppress path: regent.f's HTDBH-inventory-inverse small-tree DBH branch (which is exactly what fires
# when LHTDRG=.FALSE. OR a species is Wykoff-calibrated, IABFLG=1). The INVOKE form would (a) switch IABFLG≠1
# species to regent.f's Wykoff HT-DBH equation branch and (b) run cratet.f's ≥3-obs HT-DBH regression fit —
# both UNPORTED. So the suppress/default form is a faithful no-op; the invoke form is flagged, not silently wrong.
function kw_nohtdreg!(s::StandState, rec::KeywordRecord)
    # initre.f:2605-2674: field 1 = species (SPDECD; <0 group / 0 or blank = all / >0 one), field 2 = invoke flag.
    # field 2 > 0 ⇒ LHTDRG[sp]=TRUE (invoke the HT-DBH calibration); blank/0 ⇒ LHTDRG[sp]=FALSE (= grinit default).
    invoke = length(rec.present) >= 2 && rec.present[2] && Float32(rec.values[2]) > 0f0
    spfield = strip(rec.fields[1]); num = tryparse(Int, spfield)
    setflag(sp) = (1 <= sp <= length(s.control.ht_drag_sp)) && (s.control.ht_drag_sp[sp] = invoke)
    if isempty(spfield) || spfield == "0"
        fill!(s.control.ht_drag_sp, invoke)                 # all species
    elseif num !== nothing && num < 0                       # species group −N
        g = -num
        (1 <= g <= length(s.control.sp_groups)) || return
        for sp in s.control.sp_groups[g]; setflag(sp); end
    else
        idx, _ = resolve_species(spfield, s.variant, s.species, s.coef)
        idx > 0 && setflag(Int(idx))
    end
    return
end

# MORTMSB (base keyword option 137) — ALTERNATE "mature-stand breakup" mortality. Sets the 6 MSB params used
# by the morts.f alternate-mortality block + msbmrt.f kill routine. Field validation + reset-all-on-error mirror
# initre.f:13700-13714 exactly: a bad field dumps an error and restores ALL params to their (MSB-off) defaults.
function kw_mortmsb!(s::StandState, rec::KeywordRecord)
    c = s.control; v = rec.values; pr = rec.present; np = length(pr)
    f(i) = (np >= i && pr[i]) ? Float32(v[i]) : nothing
    reset_msb!() = (c.msb_qmd = 999f0; c.msb_slope = 0f0; c.msb_eff = 0.90f0;
                    c.msb_dlo = 0f0; c.msb_dhi = 999f0; c.msb_flag = Int32(1))
    bad() = (@warn "MORTMSB field out of range — alternate mortality disabled (defaults restored)."; reset_msb!())
    q = f(1); if q !== nothing; (q > 0f0) ? (c.msb_qmd = q) : (return bad()); end          # QMDMSB > 0
    sl = f(2); if sl !== nothing; (-10f0 <= sl <= -1.605f0) ? (c.msb_slope = sl) : (return bad()); end  # SLPMSB ∈ [−10,−1.605]
    ef = f(3); if ef !== nothing; (0f0 < ef <= 1f0) ? (c.msb_eff = ef) : (return bad()); end   # EFFMSB ∈ (0,1]
    dl = f(4); (dl !== nothing && dl >= 0f0) && (c.msb_dlo = dl)                              # DLOMSB ≥ 0
    dh = f(5); (dh !== nothing && dh >= 0f0) && (c.msb_dhi = dh)                              # DHIMSB ≥ 0
    c.msb_dlo >= c.msb_dhi && return bad()                                                    # range must be non-empty
    mf = f(6); if mf !== nothing; (1f0 <= mf <= 3f0) ? (c.msb_flag = Int32(trunc(mf))) : (return bad()); end  # MFLMSB ∈ {1,2,3}
    return
end

# CARBREPT (FFE carbon extension): request the per-cycle Stand Carbon Report. CARBCALC field 1
# selects the carbon-calculation method (0 = FFE fuel-based, 1 = JENKINS national biomass — the
# default and the only one FVSjl's live pools implement via `stand_live_carbon`/`jenkins_biomass`).
kw_carbrept!(s::StandState, ::KeywordRecord) = (s.control.carbon_report_on = true; nothing)
function kw_carbcalc!(s::StandState, rec::KeywordRecord)
    # FLD1 method (0=FFE*, 1=Jenkins), FLD2 units (0=US t/ac*, 1=metric t/ha, 2=metric t/ac). fmin.f opt 46.
    length(rec.present) >= 1 && rec.present[1] && (s.control.carbon_method = Int32(clamp(round(Int, rec.values[1]), 0, 1)))
    length(rec.present) >= 2 && rec.present[2] && (s.control.carbon_units  = Int32(clamp(round(Int, rec.values[2]), 0, 2)))
    return
end

# OPTION 13 — GROWTH (vbase/initre.f:2300): the INPUT growth-data type codes + measurement periods
# used by the LSTART calibration — field 1 = IDG (diameter data type: 0 none/increment, 1/3 = the DG
# field is past DBH → PDBH, 2 = increment), 2 = FINT (DG measurement period), 3 = IHTG (height data
# type), 4 = FINTH (HTG period), 5 = FINTM (mortality period). Defaults IDG/IHTG=0, periods=5 — which
# is FVSjl's current bit-exact behaviour (the DG field is the increment over 5 yr). The captured
# params are stored; the IDG=1/3 past-DBH interpretation + non-default FINT scaling are the WK3
# past-DBH calibration chunk (intree.f:531-537 / dgdriv.f:330) — they CHANGE the calibration, so they
# are deferred until validated against a purpose-built past-DBH / non-5-yr-FINT scenario (the WK3
# residual area, sp33/65), not wired blind.
function kw_growth!(s::StandState, rec::KeywordRecord)
    v = rec.values; pr = rec.present
    pr[1] && (s.control.growth_idg = nint(v[1]))
    pr[2] && v[2] > 0f0 && (s.control.growth_fint = Float32(v[2]))
    pr[3] && (s.control.growth_ihtg = nint(v[3]))
    pr[4] && v[4] > 0f0 && (s.control.growth_finth = Float32(v[4]))
    pr[5] && v[5] > 0f0 && (s.control.growth_fintm = Float32(v[5]))
    return
end

# OPTION 16 — INVYEAR (initre.f:895).
function kw_invyear!(s::StandState, rec::KeywordRecord)
    rec.present[1] && (s.control.cycle_year[1] = nint(rec.values[1]))
    return
end

# OPTION 94 — SITECODE (initre.f:2546): site index for the site species / all.
# Fields: array[1]=species (or 0=all), array[2]=site index.
function kw_sitecode!(s::StandState, rec::KeywordRecord)
    si = rec.values[2]
    isp = nint(rec.values[1])
    if isp <= 0
        fill!(s.plot.sp_site_index, si)
    elseif 1 <= isp <= MAXSP
        s.plot.sp_site_index[isp] = si
        s.plot.site_species = Int32(isp)        # ISISP — the site species
    end
    rec.present[2] && (s.plot.site_index = si)
    return
end

"""
    kw_managed!(s, rec)

MANAGED (initre.f:10000): set the managed-stand flag (`MANAGD`), which turns on the DGF
planted/managed diameter-growth term — `kplant = MANAGD>0 ? 1 : 0`, added as
`dg_planted[sp]·kplant` to ln(DDS) (dgf.f:179/328). With no date (or date ≤ 0): immediate
— an explicit field-2 value of 0 sets unmanaged, anything else (incl. a bare `MANAGED`
card) sets managed. A *dated* MANAGED (field 1 > 0) is scheduled by FVS for that cycle
(OPNEW act 82); that deferred path is not yet ported (rare — MANAGED is normally a
stand-level setup flag).
"""
function kw_managed!(s::StandState, rec::KeywordRecord)
    v = rec.values; p = rec.present
    idt = p[1] ? nint(v[1]) : Int32(0)
    idt > 0 && return                                  # dated MANAGED (OPNEW act 82) — deferred
    s.plot.managed = (p[2] && v[2] == 0f0) ? Int32(0) : Int32(1)
    return
end

"""
    kw_bamax!(s, rec)

BAMAX (initre.f:6800, option 66): pin the stand's maximum basal area. With field 1 > 0,
`ba_max = field 1` (LBAMAX). `site_index_setup!` then derives every species' SDImax
default from it (`sp_sdi_def = BAMAX / (0.5454154·PMSDIU)`, sdical.f:208) instead of the
per-species SDI constants, so the SDImax-driven self-thinning mortality is keyed to the
user's BAMAX. Without the keyword (`ba_max == 0`) the SDImax stays dynamic, as before.
"""
function kw_bamax!(s::StandState, rec::KeywordRecord)
    rec.present[1] && rec.values[1] > 0f0 && (s.control.ba_max = Float32(rec.values[1]))
    return
end

"""
    kw_compress!(s, rec)

COMPRESS (initre.f:8000, option 78): schedule tree-record compression (act 250) — reduce
the record list to `target` classes (field 2, default MAXTRE/2 = 1500), finding `pn1`%
(field 3, default 50) of the classes by Method-1 attribute breaks and the rest by Method-2
principal-component splitting (comprs.f). Field 1 is the date/cycle (default 1).

⚠ The compression *algorithm* (comprs.f, a 1010-line IBM-SSP-eigen PCA clustering) is NOT
yet ported — see `docs/COMPRESS_chunk_plan.md`. This handler RECOGNIZES + schedules the
keyword (icflag 250) so it is no longer silently dropped; `cuts!` currently skips icflag
250 (records pass through uncompressed) until the algorithm lands. A stand using COMPRESS
therefore still diverges from Fortran by the compression — a *tracked*, visible gap.
"""
function kw_compress!(s::StandState, rec::KeywordRecord)
    v = rec.values; p = rec.present
    yr = p[1] ? nint(v[1]) : Int32(1)                  # date/cycle (default 1)
    target = p[2] ? v[2] : Float32(MAXTRE ÷ 2)         # target records (default MAXTRE/2)
    pn1 = p[3] ? v[3] : 50f0                            # % via Method 1 (default 50)
    push!(s.control.schedule, ScheduledActivity(yr, Int32(250), (target, pn1, 0f0, 0f0, 0f0, 0f0)))
    return
end

"""
    kw_nocalib!(s, rec)

NOCALIB (initre.f:5800, option 56): disable diameter-growth self-calibration (LDGCAL) for a
species — field 1 is `0`/blank = all, `−N` = SPGROUP, else a code. Clears `dg_calib_sp[sp]`
so `calibrate_diameter_growth!` skips the COR fit (the species uses its uncalibrated DG).
(SN also clears LHTCAL, but FVSjl does not do large-tree HT self-calibration — `htg_cor` is
0 unless set by HCOR2 — so the HT side is already inert here.)
"""
function kw_nocalib!(s::StandState, rec::KeywordRecord)
    spfield = strip(rec.fields[1]); num = tryparse(Int, spfield)
    if isempty(spfield) || spfield == "0"
        fill!(s.control.dg_calib_sp, false)                 # all species
    elseif num !== nothing && num < 0                       # species group −N
        g = -num
        (1 <= g <= length(s.control.sp_groups)) || return
        for sp in s.control.sp_groups[g]
            (1 <= sp <= length(s.control.dg_calib_sp)) && (s.control.dg_calib_sp[sp] = false)
        end
    else
        idx, _ = resolve_species(spfield, s.variant, s.species, s.coef)
        idx > 0 && (s.control.dg_calib_sp[idx] = false)
    end
    return
end

"""
    kw_sdicalc!(s, rec)

SDICALC (initre.f:14000, option 1400): choose the SDI method + min-DBH thresholds. Field 1 =
`DBHSTAGE` (Reineke min DBH), field 2 = `DBHZEIDE` (Zeide min DBH), field 3 ≥ 1 ⇒ Zeide method
else Reineke. Sets the SHARED `zeide_sdi` (LZEIDE) flag, so BOTH the reported `.sum` SDI column
(`stand_sdi`) AND the SDImax self-thinning mortality (`mortality.jl`, which already reads it)
follow the chosen method — matching Fortran, where SDICALC's LZEIDE drives both.
"""
function kw_sdicalc!(s::StandState, rec::KeywordRecord)
    p = rec.present; v = rec.values
    p[1] && (s.control.dbh_stage = Float32(v[1]))
    p[2] && (s.control.dbh_zeide = Float32(v[2]))
    s.control.zeide_sdi = v[3] >= 1f0                 # field 3 ≥ 1 ⇒ Zeide, else Reineke
    s.control.sdi_method = s.control.zeide_sdi ? "ZEIDE  " : "REINEKE"
    return
end

"""
    kw_resetage!(s, rec)

RESETAGE (initre.f:9500, option 92; applied by resage.f act 443): rebase the stand age so
that at the activity's date the age equals field 2 (`IAGE = age − IDT + IY(1)`). Because
RESAGE runs after DISPLY, the reset year's own `.sum` row keeps the old age and the rebase
takes effect the following row — handled in `summary_row` from `age_reset_year`/`_age`.
Field 1 is the date (a calendar year, or a 1-based cycle number resolved against INVYEAR).
"""
function kw_resetage!(s::StandState, rec::KeywordRecord)
    v = rec.values; p = rec.present
    idt = p[1] ? nint(v[1]) : Int32(1)
    invyr = Int(s.control.cycle_year[1]); period = round(Int, s.control.year)
    period < 1 && (period = 5)
    yr = idt >= Int32(1000) ? Int(idt) : invyr + (Int(idt) - 1) * period   # cycle number → year
    s.control.age_reset_year = Int32(yr)
    s.control.age_reset_age = p[2] ? Int32(nint(v[2])) : Int32(0)
    return
end

"""
    kw_serlcorr!(s, rec)

SERLCORR (initre.f:9300, option 91): set the ARMA(1,1) serial-correlation parameters of the
stochastic diameter growth — field 1 = BJPHI (AR, default 0.74), field 2 = BJTHET (MA,
default 0.42). They define the DGSCOR autocorrelation series `BJRHO` (`autcor`), so a
non-default value changes the per-cycle DG variance/covariance multipliers.
"""
function kw_serlcorr!(s::StandState, rec::KeywordRecord)
    rec.present[1] && (s.control.dg_bjphi = Float32(rec.values[1]))
    rec.present[2] && (s.control.dg_bjthet = Float32(rec.values[2]))
    return
end

# READCOR{D,H,R} read a continuation block of MAXSP per-species correction terms (Fortran 8F10.0,
# initre.f:5650). A blank field reads as 0.0 (Fortran F10.0), and the apply guards `> 0`, so an
# unspecified species gets no correction. Returns the freshly-read array (replacing the default 1s).
function read_species_corr!(kr::KeywordReader, nsp::Integer)
    vals = zeros(Float32, MAXSP)
    @inbounds for ln in 1:cld(nsp, 8)
        line = rpad(read_raw_line!(kr), 80)
        for f in 1:8
            i = (ln - 1) * 8 + f
            i > nsp && break
            field = strip(line[(f-1)*10+1 : f*10])
            isempty(field) || (vals[i] = something(tryparse(Float32, field), 0f0))
        end
    end
    return vals
end

# OPTIONS 54/55 — READCORD/REUSCORD (initre.f:5600/5700): large-tree DIAMETER-growth correction
# COR2 (added as ln(COR2) to DGCON before calibration, dgf.f:1168). READ reloads the terms; REUSE
# re-enables the previously-read terms without reading (multi-stand carry-over).
function kw_readcord!(s::StandState, kr::KeywordReader)
    s.control.dg_cor2 = read_species_corr!(kr, nspecies(s.variant)); s.control.dg_cor2_on = true; return
end
kw_reuscord!(s::StandState) = (s.control.dg_cor2_on = true; return)

# OPTIONS 67/68 — READCORH/REUSCORH (initre.f:6900/7000): large-tree HEIGHT-growth correction
# HCOR2 (added as ln(HCOR2) to HTCON before calibration, htgf.f:332).
function kw_readcorh!(s::StandState, kr::KeywordReader)
    s.control.htg_cor2 = read_species_corr!(kr, nspecies(s.variant)); s.control.htg_cor2_on = true; return
end
kw_reuscorh!(s::StandState) = (s.control.htg_cor2_on = true; return)

# OPTIONS 73/74 — READCORR/REUSCORR (initre.f:7500/7600): small-tree HEIGHT-growth correction
# RCOR2 (the small-tree height constant RHCON = RCOR2, a multiplier on the REGENT con, regent.f:585).
function kw_readcorr!(s::StandState, kr::KeywordReader)
    s.control.regh_cor2 = read_species_corr!(kr, nspecies(s.variant)); s.control.regh_cor2_on = true; return
end
kw_reuscorr!(s::StandState) = (s.control.regh_cor2_on = true; return)

"""
    kw_dgstdev!(s, rec)

DGSTDEV (initre.f:5900, option 57): set DGSD, the number of standard deviations the
stochastic serial-correlation diameter-growth variation is bounded to (default 2.0,
grinit.f:171). `DGSD < 1` turns the random variation OFF (deterministic DG). Consumed by
`dgscor!` + the OLDRN bound in `diameter_growth!`/`small_tree_growth!`.
"""
function kw_dgstdev!(s::StandState, rec::KeywordRecord)
    rec.present[1] && (s.control.dg_stddev_bound = Float32(rec.values[1]))
    return
end

"""
    kw_rannseed!(s, rec)

RANNSEED (initre.f:6300, option 61): reseed the main random-number stream (RANSED). A
present non-zero field 1 installs that seed (forced odd); a present field-1 of 0 means
GETSED — a clock-based seed, intentionally NOT reproduced (non-deterministic); a blank
field restarts the stream from its saved seed `ss`. Set at keyword time, before any draw.
"""
function kw_rannseed!(s::StandState, rec::KeywordRecord)
    if rec.present[1]
        rec.values[1] == 0f0 && return                 # GETSED clock seed — non-deterministic
        ranseed!(s.rng, true, Float32(rec.values[1]))
    else
        ranseed!(s.rng, false, 0f0)                     # restart from the saved seed
    end
    return
end

"""
    kw_sdimax!(s, rec)

SDIMAX (initre.f:3072, option 89): override the maximum stand density index and the self-
thinning bounds. Field 1 is a species (blank/0 = all, −N = SPGROUP, else a code); field 2
(if > 0) sets that species' `sp_sdi_def` (SDIDEF, the per-species SDImax — a flagged user
value that `site_index_setup!` then leaves untouched). Fields 5/6 set the lower/upper
self-thinning *percents* PMSDIL (≥10) / PMSDIU (≤95) — stored here as the fractions the
mortality model uses (`pct_sdimax_mort_lo`/`_hi`, ÷100).
"""
function kw_sdimax!(s::StandState, rec::KeywordRecord)
    p = s.plot; v = rec.values; pr = rec.present
    if pr[2] && v[2] > 0f0                              # per-species max SDI (SDIDEF + MAXSDI)
        spfield = strip(rec.fields[1])
        num = tryparse(Int, spfield)
        targets = Int[]
        if isempty(spfield) || spfield == "0"
            append!(targets, 1:length(p.sp_sdi_def))    # all species
        elseif num !== nothing && num < 0               # species group −N
            g = -num
            (1 <= g <= length(s.control.sp_groups)) && append!(targets, s.control.sp_groups[g])
        else
            idx, _ = resolve_species(spfield, s.variant, s.species, s.coef)
            idx > 0 && push!(targets, idx)
        end
        for sp in targets
            (1 <= sp <= length(p.sp_sdi_def)) && (p.sp_sdi_def[sp] = Float32(v[2]))
        end
    end
    pr[5] && (p.pct_sdimax_mort_lo = max(Float32(v[5]), 10f0) / 100f0)   # PMSDIL (≥10) → fraction
    pr[6] && (p.pct_sdimax_mort_hi = min(Float32(v[6]), 95f0) / 100f0)   # PMSDIU (≤95) → fraction
    return
end

# OPTION 14 — STDINFO (initre.f:808): stand info. Forest/habitat decoding
# (FORKOD/HABTYP) is deferred; the geometry/site fields are set here.
function kw_stdinfo!(s::StandState, rec::KeywordRecord)
    p, v = s.plot, rec.values
    p.user_forest_code = nint(v[1])
    # SN: STDINFO field 2 is the habitat/ecological-unit field, decoded by HABTYP
    # (numeric → index into SNECU; alpha → matched, uppercased) into the PCOM code.
    rec.present[2] && (p.eco_unit = rpad(resolve_eco_unit(rec.fields[2], rec.values[2]), 10))
    rec.present[3] && (p.stand_age = nint(v[3]))
    rec.present[4] && (p.aspect = v[4] * 0.0174533f0)   # degrees → radians (utils.f)
    rec.present[5] && (p.slope  = v[5] / 100f0)         # percent → fraction (utils.f)
    (rec.present[6] && v[6] > 0f0) && (p.elevation = v[6])
    if rec.present[9]
        org = nint(v[9])
        p.stand_origin = (org < 0 || org > 1) ? Int32(0) : org
    end
    # Fort Bragg (forkod.f:137-140): forest 701 (location 701xx) ⇒ IFOR=20 (special
    # longleaf/loblolly DG + bark equations) AND KODFOR is remapped to NC Uwharrie
    # district 81110 (region 8) for downstream FORTYP + VOLEQDEF. Without the KODFOR
    # remap, VOLEQDEF sees region 7 (70106÷10000) and assigns NO R8 Clark equation ⇒
    # every tree gets zero volume.
    if div(p.user_forest_code, 100) == 701
        p.forest_idx = Int32(20)
        p.user_forest_code = Int32(81110)
    end
    # FORKOD phase 3: default lat/long/elev from the forest code (forkod.f:193).
    lat0, long0, elev0 = forest_location(s.coef, div(p.user_forest_code, 100))
    p.latitude  == 0f0 && (p.latitude  = lat0)
    p.longitude == 0f0 && (p.longitude = long0)
    p.elevation == 0f0 && (p.elevation = elev0)
    return
end

# Thinning/harvest keyword → CUTS method code (icflag). Extended per method as the
# cuts! port lands; THINDBH is the first (milestone 1). (cuts.f label dispatch.)
const _THIN_ICFLAG = Dict("THINBTA" => Int32(3), "THINATA" => Int32(4),
                          "THINBBA" => Int32(5), "THINABA" => Int32(6),
                          "THINPRSC" => Int32(7), "THINDBH" => Int32(8),
                          "THINSDI" => Int32(10), "THINHT" => Int32(12),
                          "THINRDEN" => Int32(14), "THINAUTO" => Int32(1),
                          "THINCC" => Int32(11))

# Parse a THIN* activity: field 1 = calendar year, fields 2-7 = the 6 method params.
# Stores a ScheduledActivity for `cuts!` to apply on the matching cycle.
function kw_thin!(s::StandState, rec::KeywordRecord, icflag::Int32)
    v = rec.values
    # Blank date field defaults to cycle 1 (IDT=1, initre.f:1189) — a cycle number, not a
    # year. cuts! interprets dates < 1000 as cycle numbers (FVS cycle = FVSjl cycle + 1).
    yr = rec.present[1] ? nint(v[1]) : Int32(1)
    params = ntuple(i -> Float32(v[i + 1]), 6)
    push!(s.control.schedule, ScheduledActivity(yr, icflag, params))
    return
end

# Decode a species keyword field into the `sp_field_matches` selector: 0 = all, −g = SPGROUP
# group g, or the resolved 1-based species index. (Mirrors SPDECD's IS for a single field.)
function species_selector(s::StandState, spfield::AbstractString)::Int
    f = strip(spfield)
    isempty(f) && return 0
    # SPDECD (spdecd.f:34-114): ISP = IFIX(ARRAY(IPOS)) — a NUMERIC species field is decoded by its
    # numeric value, NOT by an alpha/FIA/PLANTS string match: negative = species GROUP code; 0 = ALL
    # (or alpha code in KARD); a positive integer 1..MAXSP is the species SEQUENCE INDEX directly
    # (ISP=IFIX(ARRAY); KARD=CNSP(ISP)); > MAXSP is an error (ignored here). Only a genuinely alpha
    # field (non-numeric) falls through to the alpha/FIA/PLANTS decode.
    x = tryparse(Float64, f)
    if x !== nothing
        isp = trunc(Int, x)
        isp < 0 && return isp                       # species group code (handled by callers)
        isp == 0 && return 0                        # ALL
        return isp <= nspecies(s.variant) ? isp : 0 # positive = species index (SPDECD); > MAXSP ⇒ ignore
    end
    idx, _ = resolve_species(f, s.variant, s.species, s.coef)
    return idx > 0 ? idx : 0
end

# OPTION 138 — SETSITE (initre.f:13800, act 120): schedule a mid-run site change. Field 1 = date
# (calendar year, or cycle number < 1000), 2 = habitat type, 3 = BAMAX, 4 = species (SPDECD), 5 =
# site index, 6 = site-index flag (0 = direct value, else % change), 7 = SDImax. OPNEW stores
# ARRAY(2..7) as the 6 activity params; `apply_setsite!` enacts them at the matching cycle.
function kw_setsite!(s::StandState, rec::KeywordRecord)
    v = rec.values; pr = rec.present
    date   = pr[1] ? nint(v[1]) : Int32(1)
    ihab   = pr[2] ? Float32(v[2]) : 0f0
    bamax  = pr[3] ? Float32(v[3]) : 0f0
    isp    = Float32(species_selector(s, length(rec.fields) >= 4 ? rec.fields[4] : ""))
    si     = pr[5] ? Float32(v[5]) : 0f0
    siflag = pr[6] ? Float32(v[6]) : 0f0
    sdimax = pr[7] ? Float32(v[7]) : 0f0
    push!(s.control.schedule, ScheduledActivity(date, Int32(120), (ihab, bamax, isp, si, siflag, sdimax)))
    return
end

"""
    apply_setsite!(s) -> Bool

Enact any SETSITE (act 120) scheduled for the current cycle (grincr.f:89-200): set the per-species
site index `sp_site_index` (SITEAR) directly or as a % change (clamped ≥ 1), and optionally BAMAX
(→ `ba_max` + the bamax-derived `sp_sdi_def`) and SDImax (`sp_sdi_def`). Then recompute the
site-dependent DG constants (`dgcons!` = the FVS RCON). Returns whether anything fired. Habitat
(param 1) is not yet wired — SN growth keys off forest type, not the habitat code; a non-zero
habitat is ignored (documented gap). No-op unless a SETSITE is due.
"""
function apply_setsite!(s::StandState)::Bool
    isempty(s.control.schedule) && return false
    p = s.plot; yr = current_cycle_year(s); fvscyc = Int(s.control.cycle) + 1
    applied = false
    for a in s.control.schedule
        a.icflag == Int32(120) || continue
        (Int(a.year) == yr || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)) || continue
        prm = a.params
        bamax = prm[2]; isp = round(Int, prm[3]); si = prm[4]; siflag = prm[5]; sdimax = prm[6]
        bamax > 0f0 && (s.control.ba_max = bamax)
        @inbounds for sp in 1:MAXSP
            sp_field_matches(s.control, isp, sp) || continue
            if siflag == 0f0
                si > 0f0 && (p.sp_site_index[sp] = si)
            else
                si != 0f0 && (p.sp_site_index[sp] += p.sp_site_index[sp] * si / 100f0)
            end
            p.sp_site_index[sp] < 1f0 && (p.sp_site_index[sp] = 1f0)
            bamax > 0f0 && (p.sp_sdi_def[sp] = bamax / (0.5454154f0 * (p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 0.85f0)))
            sdimax > 0f0 && (p.sp_sdi_def[sp] = sdimax)
        end
        applied = true
    end
    applied && dgcons!(s)
    return applied
end

# Growth/mortality multipliers (BAIMULT/HTGMULT/MORTMULT/REGHMULT/REGDMULT — MULTS,
# base/mults.f). Field 1 = date, field 2 = species (0 = all), field 3 = multiplier.
function kw_mult!(s::StandState, rec::KeywordRecord, kind::Symbol)
    v = rec.values; pr = rec.present
    # MORTMULT (kind :mort) additionally carries a DBH window (PRM(3)/PRM(4) → XMDIA1/XMDIA2,
    # morts.f:170-171); the others ignore fields 4-5. Defaults: D1=0, D2=99999 (all trees).
    # FMORTMLT (kind :fmort, fmin.f:2506) SWAPS the species/multiplier fields vs the growth mults:
    # it is (date, MULTIPLIER, species, dbh_lo, dbh_hi) (SPDECD reads species from field 3, ARRAY(2) is the
    # multiplier), whereas MORTMULT etc. are (date, species, multiplier). It is also DBH-windowed.
    windowed = kind === :mort || kind === :fixdg || kind === :fixhtg || kind === :crn || kind === :fmort
    spfield, valfield = kind === :fmort ? (3, 2) : (2, 3)
    d1 = (windowed && pr[4]) ? Float32(v[4]) : 0f0
    d2 = (windowed && pr[5] && Float32(v[5]) > 0f0) ? Float32(v[5]) : 99999f0
    push!(s.control.multipliers,
          GrowthMultiplier(kind, Int32(nint(v[1])), Int32(nint(v[spfield])), Float32(v[valfield]), d1, d2))
    return
end

"""
    active_fmort_mult(control, sp, year, dbh) -> Float32

The FMORTMLT fire-caused-mortality multiplier in effect for species `sp` at cycle `year`, applied only
when the tree's `dbh` is in the keyword's half-open DBH window [d1, d2) (fmin.f:2506; fmeff.f:340
`PMORT = PMORT*FMORTMLT(I)`). 1.0 when none applies. Same most-recent / species-specific precedence.
"""
function active_fmort_mult(c::Control, sp::Integer, year::Integer, dbh::Real)
    isempty(c.multipliers) && return 1f0
    val = 1f0; d1 = 0f0; d2 = 99999f0; bestyr = typemin(Int32); bestspec = false
    @inbounds for m in c.multipliers
        (m.kind === :fmort && sp_field_matches(c, m.species, sp) && m.year <= year) || continue
        spec = m.species != 0
        if m.year > bestyr || (m.year == bestyr && spec && !bestspec)
            val = m.value; d1 = m.d1; d2 = m.d2; bestyr = m.year; bestspec = spec
        end
    end
    return (d1 <= dbh < d2) ? val : 1f0
end

"""
    active_mort_mult(control, sp, year, dbh) -> Float32

The MORTMULT multiplier in effect for species `sp` at cycle `year`, applied only when
the tree's `dbh` falls in the keyword's DBH window [d1, d2) (morts.f:518). Outside the
window — or when no MORTMULT applies — returns 1.0. Same most-recent / species-specific
precedence as `active_multiplier`.
"""
function active_mort_mult(c::Control, sp::Integer, year::Integer, dbh::Real)
    isempty(c.multipliers) && return 1f0
    val = 1f0; d1 = 0f0; d2 = 99999f0; bestyr = typemin(Int32); bestspec = false
    @inbounds for m in c.multipliers
        (m.kind === :mort && sp_field_matches(c, m.species, sp) && m.year <= year) || continue
        spec = m.species != 0
        if m.year > bestyr || (m.year == bestyr && spec && !bestspec)
            val = m.value; d1 = m.d1; d2 = m.d2; bestyr = m.year; bestspec = spec
        end
    end
    return (d1 <= dbh < d2) ? val : 1f0
end

"""
    active_crn_mult(control, sp, year, dbh) -> Float32

The CRNMULT crown-ratio multiplier in effect for species `sp` at cycle `year`, when the
tree's `dbh` is in the keyword's CLOSED DBH window [d1, d2] (sn/crown.f:318). Persists from
the keyword date onward (most recent / species-specific wins). 1.0 when none applies.
"""
function active_crn_mult(c::Control, sp::Integer, year::Integer, dbh::Real)
    isempty(c.multipliers) && return 1f0
    val = 1f0; d1 = 0f0; d2 = 99999f0; bestyr = typemin(Int32); bestspec = false
    @inbounds for m in c.multipliers
        (m.kind === :crn && sp_field_matches(c, m.species, sp) && m.year <= year) || continue
        spec = m.species != 0
        if m.year > bestyr || (m.year == bestyr && spec && !bestspec)
            val = m.value; d1 = m.d1; d2 = m.d2; bestyr = m.year; bestspec = spec
        end
    end
    return (d1 <= dbh <= d2) ? val : 1f0
end

"""
    apply_fix_scalers!(s, stash, kind, fint)

FIXDG/FIXHTG (grincr.f:451-525): a ONE-SHOT per-cycle scaler. In the cycle whose year
range [cyc_start, cyc_start+period) contains the keyword date, multiply DG (kind=:fixdg)
or HTG (:fixhtg) by `value` for every tree of the matching species (0 = all) whose DBH is
in [d1, d2). The tripled upper/lower records (stash dgU/dgL, htgU/htgL) get the same factor
(FVS scales DG(ITFN)/DG(ITFN+1)). Runs after all growth, before mortality (MORTS reads the
scaled DG). Multiple scalers firing the same cycle compound, in keyword order.
"""
function apply_fix_scalers!(s::StandState, stash, kind::Symbol, fint::Float32)
    isempty(s.control.multipliers) && return s
    # cycle window [start,end) from the IY schedule (TIMEINT/CYCLEAT-aware; uniform = +period)
    cyc_start = current_cycle_year(s)
    cyc_end = cycle_year_at(s.control, Int(s.control.cycle) + 1)
    t = s.trees
    nlive = stash === nothing ? t.n : stash.nlive
    isdg = kind === :fixdg
    @inbounds for m in s.control.multipliers
        m.kind === kind || continue
        # one-shot: a fixed-year scaler matches exactly one cycle's [start,end) range. A
        # date before the first cycle fires in cycle 0 (Fortran OPFIND past-date behaviour).
        (cyc_start <= m.year < cyc_end || (m.year < cyc_start && s.control.cycle == 0)) || continue
        for i in 1:nlive
            sp_field_matches(s.control, m.species, t.species[i]) || continue
            d = t.dbh[i]
            (m.d1 <= d < m.d2) || continue
            if isdg
                t.diam_growth[i] *= m.value
                stash !== nothing && (stash.dgU[i] *= m.value; stash.dgL[i] *= m.value)
            else
                t.ht_growth[i] *= m.value
                stash !== nothing && (stash.htgU[i] *= m.value; stash.htgL[i] *= m.value)
            end
        end
    end
    return s
end

# SPGROUP (vbase/initre.f:4726): define a species group. The keyword-line field is the group
# name (optional → auto "GROUPnn"); the NEXT record lists the member species (alpha or numeric
# codes, space-separated, ALL/0/dups skipped, max 90, `&` continues — rare, not handled). Groups
# are numbered in definition order and referenced from a species field by the negative index −N.
# Max 30 groups.
function kw_spgroup!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    line = read_raw_line!(kr)
    length(s.control.sp_groups) >= 30 && return
    members = Int32[]
    for tok in split(line)
        u = uppercase(tok)
        (u == "ALL" || u == "0") && continue
        n = tryparse(Int, tok)
        idx = n !== nothing ? (n > 0 ? Int32(n) : Int32(0)) :
              first(resolve_species(tok, s.variant, s.species, s.coef))
        (idx <= 0 || idx in members) && continue          # skip invalid + duplicates
        push!(members, idx)
        length(members) >= 90 && break
    end
    push!(s.control.sp_groups, members)
    nm = isempty(rec.fields) ? "" : strip(rec.fields[1])
    push!(s.control.sp_group_names,
          isempty(nm) ? "GROUP" * lpad(length(s.control.sp_groups), 2, '0') : uppercase(String(nm)))
    s.control.n_spgroups = Int32(length(s.control.sp_groups))
    return
end

"""
    sp_field_matches(control, isp, sp) -> Bool

Does species `sp` match a keyword species field `isp`? 0 = all species, >0 = that single
species, <0 = every member of SPGROUP group −isp. The one place the ISPCC<0 group branch lives.
"""
@inline function sp_field_matches(c::Control, isp::Integer, sp::Integer)
    isp == 0 && return true
    isp > 0 && return Int(isp) == Int(sp)
    g = -Int(isp)
    return g <= length(c.sp_groups) && (Int32(sp) in c.sp_groups[g])
end

# TREESZCP (base/keywds.f:51 / SIZCAP): a per-species maximum tree size. Field 1 =
# species (0 = all), 2 = cap DBH (SIZCAP[1]), 3 = annual mortality rate applied to
# trees that reach the cap (SIZCAP[2]), 4 = no-mortality flag IDMFLG (SIZCAP[3]),
# 5 = height cap (SIZCAP[4]). Applied immediately to every matching species; the cap
# then governs diameter growth (dgbnd), height growth (htgf) and a size-cap mortality
# floor (morts). No date field — the cap holds for the whole projection.
function kw_treeszcp!(s::StandState, rec::KeywordRecord)
    v = rec.values; pr = rec.present
    isp   = nint(v[1])
    cap   = pr[2] ? Float32(v[2]) : 999f0
    mrate = pr[3] ? Float32(v[3]) : 1f0
    flag  = pr[4] ? Float32(v[4]) : 0f0
    htcap = (pr[5] && Float32(v[5]) > 0f0) ? Float32(v[5]) : 999f0
    sc = s.control.sp_size_cap
    # species selection: 0 = all, >0 = single, <0 = SPGROUP group −isp.
    rng = isp < 0 ? (isp <= -1 && -isp <= length(s.control.sp_groups) ? s.control.sp_groups[-isp] : Int32[]) :
          isp == 0 ? (1:size(sc, 1)) : (isp:isp)
    @inbounds for sp in rng
        sc[sp, 1] = cap; sc[sp, 2] = mrate; sc[sp, 3] = flag; sc[sp, 4] = htcap
    end
    return
end

# HTGSTOP (act 110) / TOPKILL (act 111) top-damage events (htgstp.f). Keyword fields:
# date, species (0/neg = all/group), HT1, HT2, PRB (damage prob), AVEPRB, STDPBR (BACHLO
# mean/sd of the kill proportion). Stored as a ScheduledActivity (icflag = act code).
function kw_htgstp!(s::StandState, rec::KeywordRecord, act::Integer)
    v = rec.values; pr = rec.present
    f(i, dflt) = (pr[i] ? Float32(v[i]) : dflt)
    p = (f(2, 0f0), f(3, 0f0), f(4, 9999f0), f(5, 1f0), f(6, 0f0), f(7, 0f0))
    push!(s.control.htgstp_events, ScheduledActivity(Int32(nint(v[1])), Int32(act), p))
    return
end

# FIXMORT (activity 97, morts.f:781): a one-shot forced-mortality override. Fields: date,
# species (0 = all), rate, DBH window d1/d2, option PRM(5) (0=replace,1=add,2=max,3=mult →
# IP 1/2/3/4), point/size flag PRM(6). Normal (non-concentration) path only — PRM(6)≥10 is
# deferred. The rate is clamped here exactly as Fortran (≤1 when PRM(5) absent; [0,1] when
# PRM(5)<3) so the apply step is a clean per-tree formula.
function kw_fixmort!(s::StandState, rec::KeywordRecord)
    v = rec.values; pr = rec.present
    isp = Float32(nint(v[2]))
    rate = Float32(v[3])
    d1 = (pr[4] ? Float32(v[4]) : 0f0); d1 < 0f0 && (d1 = 0f0)
    d2 = (pr[5] && Float32(v[5]) > 0f0) ? Float32(v[5]) : 999f0
    ip = 1
    if pr[6]                      # PRM(5) given (NP>4)
        p5 = Float32(v[6])
        p5 < 3f0 && (rate = clamp(rate, 0f0, 1f0))
        ip = p5 == 1f0 ? 2 : p5 == 2f0 ? 3 : p5 == 3f0 ? 4 : 1
    else
        rate > 1f0 && (rate = 1f0)
    end
    pflag = pr[7] ? Float32(v[7]) : 0f0
    push!(s.control.fixmort_events,
          ScheduledActivity(Int32(nint(v[1])), Int32(97), (isp, rate, d1, d2, Float32(ip), pflag)))
    return
end

# VOLUME (volkey.f:9915) — override the CUBIC merch standards for a species (0=all, <0=group).
# Keyword fields after the date: 2=species, 3=DBHMIN, 4=TOPD, 5=STMP, 6=FRMCLS, 7=METHC,
# 8=SCFMIND, 9=SCFTOPD, 10=SCFSTMP. FRMCLS/METHC are not used by the SN R8 Clark taper path
# (it has no form-class or method selector), so we carry only the 7 standards the model reads:
# species + DBHMIN/TOPD/STMP/SCFMIND/SCFTOPD in params, SCFSTMP in aux.
function kw_volume!(s::StandState, rec::KeywordRecord)
    v = rec.values
    isp = Float32(nint(v[2]))
    push!(s.control.volume_events,
          ScheduledActivity(Int32(nint(v[1])), Int32(217),
              (isp, Float32(v[3]), Float32(v[4]), Float32(v[5]), Float32(v[8]), Float32(v[9])),
              Float32(v[10])))
    return
end

# BFVOLUME (volkey.f:9905) — override the BOARD-FOOT merch standards for a species.
# Fields after the date: 2=species, 3=BFMIND, 4=BFTOPD, 5=BFSTMP, 6=FRMCLS, 7=METHB
# (FRMCLS/METHB unused by the R8 Clark path).
function kw_bfvolume!(s::StandState, rec::KeywordRecord)
    v = rec.values
    isp = Float32(nint(v[2]))
    push!(s.control.volume_events,
          ScheduledActivity(Int32(nint(v[1])), Int32(218),
              (isp, Float32(v[3]), Float32(v[4]), Float32(v[5]), 0f0, 0f0), 0f0))
    return
end

# Load a per-species defect curve (CFDEFT or BFDEFT rows, sdefet.f:151): the five values at DBH
# 5/10/15/20/25" go to rows 2..6, the 25" value extends flat to rows 7,8,9 (DBH 30/35/40"); row 1
# (DBH 0") stays 0. `isp` selects 0=all / +species / −SPGROUP.
function _set_defect!(c, dmat::Matrix{Float32}, isp::Integer, vals::NTuple{5,Float32})
    setone(sp) = @inbounds begin
        dmat[2, sp] = vals[1]; dmat[3, sp] = vals[2]; dmat[4, sp] = vals[3]
        dmat[5, sp] = vals[4]; dmat[6, sp] = vals[5]
        dmat[7, sp] = vals[5]; dmat[8, sp] = vals[5]; dmat[9, sp] = vals[5]
    end
    if isp == 0
        for sp in 1:size(dmat, 2); setone(sp); end
    elseif isp > 0
        isp <= size(dmat, 2) && setone(isp)
    else
        g = -isp
        (1 <= g <= length(c.sp_groups)) || return
        for sp in c.sp_groups[g]; setone(sp); end
    end
    return
end

_defect_vals(v) = (Float32(v[3]), Float32(v[4]), Float32(v[5]), Float32(v[6]), Float32(v[7]))

# MCDEFECT (sdefet.f, IACTK 215) — per-species CUBIC defect curve (CFDEFT). Fields: 1=date,
# 2=species, 3..7 = defect fractions at DBH 5/10/15/20/25". sdefet.f:84-120: a DATED card is
# scheduled (OPNEW) to take effect at that cycle; an UNDATED card changes the terms now.
# `_sched_defect!` shares that gate with BFDEFECT (IACTK 216).
function _sched_defect!(s::StandState, rec::KeywordRecord, icflag::Int32, immediate_mat)
    vals = _defect_vals(rec.values); isp = nint(rec.values[2])
    if rec.present[1] && rec.values[1] > 0f0          # dated → defer to the cycle (sdefet.f OPNEW)
        push!(s.control.volume_events,
              ScheduledActivity(Int32(nint(rec.values[1])), icflag,
                  (Float32(isp), vals[1], vals[2], vals[3], vals[4], vals[5]), 0f0))
    else                                               # undated → change the terms now
        _set_defect!(s.control, immediate_mat, isp, vals)
    end
end
kw_mcdefect!(s::StandState, rec::KeywordRecord) = _sched_defect!(s, rec, Int32(215), s.control.sp_cf_defect)
# BFDEFECT (sdefet.f, IACTK 216) — per-species BOARD-FOOT defect curve (BFDEFT); reduces board
# feet AND sawtimber cubic by IBDF% in the volume path.
kw_bfdefect!(s::StandState, rec::KeywordRecord) = _sched_defect!(s, rec, Int32(216), s.control.sp_bf_defect)

# Set the log-linear form-model coefficients B0/B1 (CFLA/BFLA) for a species (sdefln.f): field 1 =
# species (alpha/FIA/−group), field 2 = B0 (intercept), field 3 = B1 (slope). `b0`/`b1` are the
# target per-species arrays. Only present fields are written (LNOTBK), preserving the 0/1 defaults.
function _set_form!(s::StandState, b0::Vector{Float32}, b1::Vector{Float32}, rec::KeywordRecord)
    spfield = strip(rec.fields[1])
    has0 = rec.present[2]; has1 = rec.present[3]
    v0 = Float32(rec.values[2]); v1 = Float32(rec.values[3])
    setone(sp) = @inbounds begin
        (1 <= sp <= length(b0)) || return
        has0 && (b0[sp] = v0); has1 && (b1[sp] = v1)
    end
    fnum = tryparse(Float64, spfield)        # numeric species field (may be "0." / "-1")
    if fnum !== nothing
        isp = round(Int, fnum)
        if isp == 0                          # SPDECD IS=0 ⇒ all species
            for sp in 1:length(b0); setone(sp); end
        elseif isp < 0                       # −N ⇒ SPGROUP N
            g = -isp
            (1 <= g <= length(s.control.sp_groups)) || return
            for sp in s.control.sp_groups[g]; setone(sp); end
        else
            setone(isp)
        end
    else                                     # alpha / FIA code
        idx, _ = resolve_species(spfield, s.variant, s.species, s.coef)
        idx > 0 && setone(Int(idx))
    end
    return
end

# FERTILIZE (ffin.f) — schedule a fertilizer application. Fields: 1=date, 2=N (forced to 200, the
# only representable amount), 3=P / 4=K (ignored), 5=efficacy multiplier (default 1). Only the
# efficacy feeds the response (FFERT), so we carry just that (activity 260, params[1]=efficacy).
function kw_fertiliz!(s::StandState, rec::KeywordRecord)
    v = rec.values
    yr = nint(rec.present[1] ? v[1] : 1f0)
    eff = rec.present[5] ? Float32(v[5]) : 1f0
    push!(s.control.fertilize_events,
          ScheduledActivity(Int32(yr), Int32(260), (eff, 0f0, 0f0, 0f0, 0f0, 0f0)))
    return
end

# MCFDLN (sdefln.f, option 39) — set the CUBIC form-model coefficients CFLA0/CFLA1 (activates the
# log-linear merch-cubic form/defect correction in the volume path).
kw_mcfdln!(s::StandState, rec::KeywordRecord) =
    _set_form!(s, s.control.sp_cf_form0, s.control.sp_cf_form1, rec)

# BFFDLN (sdefln.f, option 40) — set the BOARD-FOOT form-model coefficients BFLA0/BFLA1.
kw_bffdln!(s::StandState, rec::KeywordRecord) =
    _set_form!(s, s.control.sp_bf_form0, s.control.sp_bf_form1, rec)

# VOLEQNUM (initre.f:5061) — override the cubic NVEL volume-equation id (VEQNNC) for a species.
# Field 1 = species (alpha code / FIA number / −N SPGROUP), field 2 = the 10-char equation id (e.g.
# "841CLKE318"). Stored and applied AFTER VOLEQDEF assigns the defaults (apply_voleqnum_overrides!).
function kw_voleqnum!(s::StandState, rec::KeywordRecord)
    spfield = strip(rec.fields[1]); eq = strip(rec.fields[2])
    isempty(eq) && return
    eqs = String(eq)
    num = tryparse(Int, spfield)
    if num !== nothing && num < 0                       # species group −N
        g = -num
        (1 <= g <= length(s.control.sp_groups)) || return
        for sp in s.control.sp_groups[g]
            push!(s.control.voleqnum_overrides, (Int32(sp), eqs))
        end
    else                                                # alpha / FIA code
        idx, _ = resolve_species(spfield, s.variant, s.species, s.coef)
        idx > 0 && push!(s.control.voleqnum_overrides, (Int32(idx), eqs))
    end
    return
end

"Apply VOLEQNUM overrides to `species.vol_eq` (after VOLEQDEF has set the defaults)."
function apply_voleqnum_overrides!(s::StandState)
    for (sp, eq) in s.control.voleqnum_overrides
        1 <= sp <= length(s.species.vol_eq) && (s.species.vol_eq[sp] = eq)
    end
    return s
end

"""
    apply_fixmort!(s, killed, n, fint)

FIXMORT (morts.f:781-1042): a one-shot forced-mortality override applied AFTER the BA-check, the
last word on `killed[]`. In the cycle whose range holds the event date, for each of the `n`
original records of the matching species (0 = all) with d1≤DBH<d2, override the kill by the IP
option (1 replace P·rate / 2 add / 3 max / 4 multiply). PRM(6) requests *concentration* of that
mortality (morts.f:838): KBIG by size (10 bottom-up, 20 top-down), KPOINT onto points (1), or both
(11/21). Concentration pools the window's kill into XMORE and re-imposes it whole-record at a time
in the chosen order (size-ranked / point-by-point / size-within-point) until XMORE is spent.
"""
function apply_fixmort!(s::StandState, killed::AbstractVector{Float32}, n::Int, fint::Float32)
    isempty(s.control.fixmort_events) && return
    # cycle window [start,end) from the IY schedule (TIMEINT/CYCLEAT-aware; uniform = +period)
    cyc_start = current_cycle_year(s)
    cyc_end = cycle_year_at(s.control, Int(s.control.cycle) + 1)
    t = s.trees
    iptinv = max(1, Int(s.plot.points_inv))
    for ev in s.control.fixmort_events
        (cyc_start <= ev.year < cyc_end || (ev.year < cyc_start && s.control.cycle == 0)) || continue
        isp = round(Int, ev.params[1]); rate = ev.params[2]
        d1 = ev.params[3]; d2 = ev.params[4]; ip = round(Int, ev.params[5])
        pflag = ev.params[6]
        kbig   = (pflag == 10f0 || pflag == 11f0) ? 1 : (pflag == 20f0 || pflag == 21f0) ? 2 : 0
        kpoint = (pflag == 1f0 || pflag == 11f0 || pflag == 21f0) ? 1 : 0
        inwin(i) = (d1 <= t.dbh[i] < d2) && sp_field_matches(s.control, isp, t.species[i])
        # Concentrate only if KBIG, or KPOINT on a real multi-point inventory (morts.f:845).
        if kbig == 0 && !(kpoint == 1 && iptinv > 1)
            # NORMAL per-tree override (morts.f:1017): set the kill in place per the IP mode.
            @inbounds for i in 1:n
                inwin(i) || continue
                p = t.tpa[i]
                if ip == 1
                    killed[i] = p * rate
                elseif ip == 2
                    killed[i] += max(0f0, p - killed[i]) * rate
                elseif ip == 3
                    killed[i] = max(killed[i], p * rate)
                else
                    killed[i] = min(p, killed[i] * rate)
                end
            end
            continue
        end
        # --- concentration: pool the window's mortality into XMORE per the IP (morts.f:847-883) ---
        xmore = 0f0
        @inbounds for i in 1:n
            inwin(i) || continue
            p = t.tpa[i]
            if ip == 1
                xmore += p * rate; killed[i] = 0f0
            elseif ip == 2
                xmore += max(0f0, p - killed[i]) * rate
            elseif ip == 3
                tmp = max(killed[i], p * rate); tmp > killed[i] && (xmore += tmp - killed[i])
            else
                xmore += killed[i] * rate; killed[i] = 0f0
            end
        end
        # Size rank (only needed when KBIG): WORK3 = ∓(DBH+DG/bark), RDPSRT descending (morts.f:889)
        # — KBIG=1 negates so the smallest grown trees sort first (bottom up), =2 largest first.
        local ord::Vector{Int32}
        if kbig >= 1
            ba = s.calib.bark_a; bb = s.calib.bark_b
            gdbh(i) = t.dbh[i] + t.diam_growth[i] / bark_ratio(ba, bb, t.species[i], t.dbh[i])
            sgn = kbig == 1 ? -1f0 : 1f0
            key = Float32[sgn * gdbh(i) for i in 1:n]
            ord = Vector{Int32}(undef, n); _rdpsrt!(key, ord)
        end
        # Kill whole records in the chosen traversal order until CREDIT reaches XMORE (last partial).
        credit = 0f0
        # returns true once XMORE is fully consumed (stop)
        take!(ix) = begin
            inwin(ix) || return false
            avail = t.tpa[ix] - killed[ix]
            if credit + avail <= xmore + 0.0001f0
                credit += avail; @inbounds killed[ix] = t.tpa[ix]; false
            else
                @inbounds killed[ix] += xmore - credit; credit = xmore; true
            end
        end
        if kpoint == 0                       # SIZE only (morts.f:901)
            for ix in ord; take!(Int(ix)) && break; end
        elseif kbig == 0                     # POINTS only (morts.f:937): point order, record order
            stop = false
            for j in 1:iptinv
                @inbounds for i in 1:n
                    t.plot_id[i] == j || continue
                    take!(i) && (stop = true; break)
                end
                stop && break
            end
        else                                 # SIZE within POINT (morts.f:978): points have priority
            stop = false
            for j in 1:iptinv
                for ix in ord
                    @inbounds (t.plot_id[ix] == j) || continue
                    take!(Int(ix)) && (stop = true; break)
                end
                stop && break
            end
        end
    end
    return
end

"""
    htgstp!(s; fint)

HTGSTOP / TOPKILL top-damage events (htgstp.f). Runs once per cycle, after TRIPLE/MORTS and
before UPDATE applies the increments (gradd.f:158). For each event firing this cycle (one-shot,
matched by date), damage the trees of its species whose height is in (HT1, HT2]:
  * act 110 (HTGSTOP): scale the height increment HTG by PKIL∈[0,1].
  * act 111 (TOPKILL): reduce height to H·(1−PKIL) (PKIL≤0.8); for a tall (H≥25, D≥6) tree
    whose Behre top diameter ≥4 set NORMHT/ITRUNC (permanent broken top); cut the crown ratio.
PKIL = BACHLO(AVEPRB, STDPBR) — deterministic (= AVEPRB, no RNG) when STDPBR≤0; a RANN escape
draw skips a tree when PRB<1. Records are visited in species-sorted (IND1) order so the RNG
stream matches FVS when the event is stochastic.
"""
function htgstp!(s::StandState; fint::Float32 = 5f0)
    isempty(s.control.htgstp_events) && return s
    # cycle window [start,end) from the IY schedule (TIMEINT/CYCLEAT-aware; uniform = +period)
    cyc_start = current_cycle_year(s)
    cyc_end = cycle_year_at(s.control, Int(s.control.cycle) + 1)
    t = s.trees; bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    species_sort!(s)
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    for ev in s.control.htgstp_events
        (cyc_start <= ev.year < cyc_end || (ev.year < cyc_start && s.control.cycle == 0)) || continue
        isp = round(Int, ev.params[1]); ht1 = ev.params[2]; ht2 = ev.params[3]
        prb = ev.params[4]; aveprb = ev.params[5]; stdpbr = ev.params[6]
        act = ev.icflag
        sp1 = isp <= 0 ? 1 : isp; sp2 = isp <= 0 ? MAXSP : isp
        @inbounds for sp in sp1:sp2
            # isp<0 (SPGROUP group) loops all species but keeps only members.
            (isp < 0 && !sp_field_matches(s.control, isp, sp)) && continue
            i1 = isct[sp, 1]; i1 == 0 && continue
            i2 = isct[sp, 2]
            for k in i1:i2
                i = ind1[k]
                t.tpa[i] <= 0f0 && continue
                h = t.height[i]
                (h <= ht1 || h > ht2) && continue
                if prb <= 0.99999f0
                    rann!(s.rng) > prb && continue
                end
                pkil = bachlo(s.rng, aveprb, stdpbr)
                pkil <= 0f0 && continue
                if act == 110
                    pkil > 1f0 && (pkil = 1f0)
                    t.ht_growth[i] *= pkil
                else
                    pkil > 0.8f0 && (pkil = 0.8f0)
                    topk = h * pkil; toph = h - topk
                    itrc2 = floor(Int32, toph * 100f0 + 0.5f0)
                    if t.trunc[i] > 0                       # already topkilled: lower the truncation
                        t.trunc[i] > itrc2 && (t.trunc[i] = itrc2)
                        t.height[i] = toph
                        continue
                    end
                    brk = bark_ratio(bark_a, bark_b, sp, t.dbh[i])
                    d = t.dbh[i] * brk
                    if !(h < 25f0 || d < 6f0)              # tall enough to maybe break permanently
                        af = t.cuft_vol[i] / (0.00545415f0 * d * d * h)
                        af = 0.44244f0 - (0.99167f0 / af) - 1.43237f0 * log(af) +
                             1.68581f0 * sqrt(af) - 0.13611f0 * af * af
                        dtk = topk / h
                        dtk = (dtk / (af * dtk + (1f0 - af))) * d
                        if dtk >= 4f0                       # permanent broken top
                            t.trunc[i] = itrc2
                            t.norm_ht[i] = floor(Int32, h * 100f0 + 0.5f0)
                        end
                    end
                    t.height[i] = toph
                    iod = t.crown_pct[i]                    # crown reduction (skip if bug-adjusted <0)
                    if iod >= 0
                        cn = (Float32(iod) / 100f0) * h - h + toph
                        new = floor(Int32, cn / toph * 100f0 + 0.5f0)
                        new < 5 && (new = Int32(5))
                        t.crown_pct[i] = -new
                    end
                end
            end
        end
    end
    return s
end

"""
    active_multiplier(control, kind, sp, year) -> Float32

The growth/mortality multiplier in effect for species `sp` at the cycle `year`
(MULTS): the most recent matching keyword wins, a species-specific one beating an
all-species (`species==0`) one of the same date. Returns 1.0 when none apply (the
common case — short-circuited).
"""
function active_multiplier(c::Control, kind::Symbol, sp::Integer, year::Integer)
    isempty(c.multipliers) && return 1f0
    val = 1f0; bestyr = typemin(Int32); bestspec = false
    @inbounds for m in c.multipliers
        (m.kind === kind && (m.species == 0 || m.species == sp) && m.year <= year) || continue
        spec = m.species != 0
        if m.year > bestyr || (m.year == bestyr && spec && !bestspec)
            val = m.value; bestyr = m.year; bestspec = spec
        end
    end
    return val
end

# THINQFA (option 141, ICFLAG 17): Q-factor diameter-distribution thin — a TWO-record
# keyword (initre.f:5981). Line 1: [year, loDBH, hiDBH, species, Qfactor, classWidth,
# target]; line 2: a single integer for the target units (≤0 BA, ==1 TPA, >1 SDI → 0/1/2).
# Defaults (initre): loDBH=0, hiDBH=24, Q=1.4, classWidth=2. qfatar → activity aux slot.
function kw_thinqfa!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    v = rec.values; pr = rec.present
    yr     = nint(v[1])
    valmin = pr[2] ? Float32(v[2]) : 0f0
    valmax = pr[3] ? Float32(v[3]) : 24f0
    spec   = pr[4] ? Float32(v[4]) : 0f0
    qfac   = pr[5] ? Float32(v[5]) : 1.4f0
    diacw  = pr[6] ? Float32(v[6]) : 2f0
    tarqfa = pr[7] ? Float32(v[7]) : 0f0
    iset   = something(tryparse(Int, strip(read_raw_line!(kr))), 0)   # 2nd record: units switch
    qfatar = iset <= 0 ? 0f0 : (iset <= 1 ? 1f0 : 2f0)
    push!(s.control.schedule,
          ScheduledActivity(Int32(yr), Int32(17), (valmin, valmax, spec, qfac, diacw, tarqfa), qfatar))
    return
end

# OPTION — ESTAB…END establishment packet (esin.f). ESTAB opens the packet (field 1 =
# date of disturbance); subsequent PLANT(430)/NATURAL(431) cards each schedule a regen
# activity (year, species, TPA, %survival, age, height, shade); END closes it. Other
# establishment keywords (TALLY/SPROUT/…) are recognized and skipped for now. At END a
# TALLY(427) trigger is scheduled at the disturbance date and IDSDAT→-9999 (esin.f:100-117).
function kw_estab!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    s.estab.active = true
    idsdat = rec.present[1] ? nint(rec.values[1]) : Int32(-1)   # ESTAB date field
    sched = s.control.schedule
    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"
            break
        elseif k == "PLANT" || k == "NATURAL"
            ic = k == "PLANT" ? Int32(430) : Int32(431)
            v = r.values
            yr   = r.present[1] ? nint(v[1]) : Int32(1)
            # species (field 2) is SPDECD-decoded: a numeric species index is used directly,
            # but a 2-letter ALPHA code (e.g. LP) must be resolved to its index — the old
            # numeric-only `Float32(v[2])` left alpha codes as 0, so the planting silently
            # created no trees.
            spfld = length(r.fields) >= 2 ? strip(r.fields[2]) : ""
            sp = if isempty(spfld) || tryparse(Int, spfld) !== nothing
                Float32(v[2])                                  # numeric index (or blank → 0)
            else
                Float32(species_selector(s, spfld))            # alpha code → species index
            end
            tpa  = Float32(v[3])
            tpa <= 0f0 && continue                                # esin.f:143 (TPA must be >0)
            surv = (v[4] < 0.001f0 || v[4] > 100f0) ? 100f0 : Float32(v[4])  # esin.f:149
            push!(sched, ScheduledActivity(yr, ic, (sp, tpa, surv, Float32(v[5]), Float32(v[6]), Float32(v[7]))))
        elseif k == "SPROUT"                                  # esin.f opt 26: enable stump sprouting
            # esin.f opt 26 fields: 1=date, 2=species (SPDECD: 0⇒all, <0⇒group), 3=SMULT (sprout-count mult,
            # default 1), 4=HMULT (height mult, default 1), 5=lower DBH (default 0), 6=upper DBH (default 999).
            # Blank species ⇒ no valid sprouting species ⇒ LSPRUT=.FALSE. (esin.f:625). The per-species + stump-
            # DBH-range table is honored in esuckr! via s.control.sprout_overrides (esuckr.f:96-205 activity 450).
            # esin.f sets LSPRUT=.TRUE. first; a BLANK species field reads ARRAY(2)=0 ⇒ SPDECD returns IS=0
            # (ALL species), so a bare `SPROUT <date>` ENABLES all-species sprouting — it does NOT disable.
            # Only the −999 "no species" sentinel (an unrecognized alpha code) flips LSPRUT off (esin.f:625);
            # NOSPROUT is the explicit disable. (Live FVScs: blank `SPROUT 2000.` == `SPROUT 2000. 0.0`.)
            isp = r.present[2] ? Float32(r.values[2]) : 0f0   # SPDECD selector (0/all, <0/group, >0/single)
            if isp == -999f0
                s.control.lsprut = false
            else
                s.control.lsprut = true
                smul  = r.present[3] ? Float32(r.values[3]) : 1f0
                hmul  = r.present[4] ? Float32(r.values[4]) : 1f0
                dmin  = r.present[5] ? Float32(r.values[5]) : 0f0
                dmax  = r.present[6] ? Float32(r.values[6]) : 999f0
                push!(s.control.sprout_overrides, (isp, smul, hmul, dmin, dmax))
                # keep the legacy scalars in sync for any single-species/all default form (back-compat)
                r.present[3] && (s.control.sprout_smult = smul)
                r.present[4] && (s.control.sprout_hmult = hmul)
            end
        elseif k == "NOSPROUT"                                # esin.f opt 27: disable sprouting
            s.control.lsprut = false
        end
        # other establishment keywords (TALLY/…) not yet ported — skipped
    end
    # END processing (esin.f:100-117): schedule the TALLY(427) establishment trigger at
    # the disturbance date, then mark IDSDAT unset so ESNUTR defaults it.
    if idsdat != Int32(-1)
        push!(sched, ScheduledActivity(max(Int32(1), idsdat), Int32(427),
                                       (Float32(idsdat), 0f0, 0f0, 0f0, 0f0, 0f0)))
    end
    s.estab.idsdat = Int32(-9999)
    return
end

"""
    kw_database!(s, rec, kr)

DATABASE (dbsin.f): the DBS output block. Read to END; `DSNOUT` sets the output SQLite file
(its name is on the FOLLOWING line, `READ(IREAD,'(A)')DSNOUT`), `SUMMARY` enables the
FVS_Summary table (ISUMARY=1). Other DBS sub-keywords (SQLOUT/COMPUTDB/TREELIDB/…) are
recognized and skipped — only the Summary table is emitted so far (see `write_dbs_summary!`).
"""
function kw_database!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"
            break
        elseif k == "DSNOUT"
            s.control.dbs_out_file = strip(read_raw_line!(kr))   # filename on the next line
        elseif k == "SUMMARY"
            s.control.dbs_summary = true
        elseif k == "TREELIDB"
            s.control.dbs_treelist = true
        elseif k == "COMPUTDB"
            s.control.dbs_compute = true
        end
    end
    return
end

"""
    kw_fmin!(s, rec, kr)

Parse the FMIN Fire & Fuels keyword block (read to END, like ESTAB). Activates the FFE
and captures the SIMFIRE event (date + fire conditions) and FLAMEADJ (flame multiplier /
crown-fire fraction). The remaining FMIN keywords (snag/fuel setup, reports) are
recognized but not yet ported, so they are skipped.
"""
# Report-only FFE keywords (text reports); recognized but intentionally no-ops here — the equivalent data
# is emitted via the DBS path. Distinguished from genuinely-unported MODEL keywords (which get a warning).
const _FFE_REPORT_KEYWORDS = Set([
    "BURNREPT", "FUELOUT", "SNAGOUT", "SNAGSUM", "MORTREPT", "FUELREPT", "MOREOUT", "LANDOUT",
    "STATFUEL", "MORTCLAS", "SNAGCLAS", "DWDVLOUT", "DWDCVOUT", "FMODLIST", "FUELFOTO", "CANFPROF",
    "SVIMAGES", "CARBCUT"])

# SNAGFALL (fmin.f:464, opt 9): per-species snag fall-rate parameters. Field 1 = species (SPDECD:
# alpha/FIA/0=all/−group), field 2 = FALLX (rate-of-fall correction, clamp ≥ 0.001), field 3 = ALLDWN
# (snag age by which the last 5% fall, clamp ≥ 0). Only present fields are written, preserving the
# fire_species_props.csv defaults; stored sparsely on FFEParams and read by snag_fall_density.
function _snagfall!(s::StandState, rec::KeywordRecord)
    fs = s.fire; fs === nothing && return
    p = fs.params
    spfield = strip(rec.fields[1])
    has_fx = rec.present[2]; has_ad = rec.present[3]
    (has_fx || has_ad) || return
    fx = max(Float32(rec.values[2]), 0.001f0); ad = max(Float32(rec.values[3]), 0f0)
    setone(sp) = @inbounds begin
        (1 <= sp <= length(s.species.class_codes[:, 1])) || return
        has_fx && (p.snag_fallx_ovr[Int32(sp)] = fx)
        has_ad && (p.snag_alldwn_ovr[Int32(sp)] = ad)
    end
    fnum = tryparse(Float64, spfield)
    if fnum !== nothing
        isp = round(Int, fnum)
        if isp == 0                                    # SPDECD IS=0 ⇒ all species
            for sp in 1:length(s.species.class_codes[:, 1]); setone(sp); end
        elseif isp < 0                                 # −N ⇒ SPGROUP N
            g = -isp
            (1 <= g <= length(s.control.sp_groups)) || return
            for sp in s.control.sp_groups[g]; setone(sp); end
        else
            setone(isp)
        end
    else
        idx, _ = resolve_species(spfield, s.variant, s.species, s.coef)
        idx > 0 && setone(Int(idx))
    end
    return
end

# SNAGDCAY (fmin.f:633, opt 11): per-species snag DECAYX (decay-rate multiplier). Field 1 = species
# (SPDECD: alpha/FIA/0=all/−group), field 2 = DECAYX (clamp ≥ 0, fmin.f:643). Only the present field
# is written, preserving the fire_species_props.csv defaults (0.07/0.21/0.35); stored sparsely on
# FFEParams and read by the snag soft-decay transition (snag_summary) + the crown-fall TSOFT (fmscro!).
function _snagdcay!(s::StandState, rec::KeywordRecord)
    fs = s.fire; fs === nothing && return
    p = fs.params
    rec.present[2] || return
    dx = max(Float32(rec.values[2]), 0f0)
    spfield = strip(rec.fields[1])
    setone(sp) = @inbounds begin
        (1 <= sp <= length(s.species.class_codes[:, 1])) || return
        p.snag_decayx_ovr[Int32(sp)] = dx
    end
    fnum = tryparse(Float64, spfield)
    if fnum !== nothing
        isp = round(Int, fnum)
        if isp == 0
            for sp in 1:length(s.species.class_codes[:, 1]); setone(sp); end
        elseif isp < 0
            g = -isp
            (1 <= g <= length(s.control.sp_groups)) || return
            for sp in s.control.sp_groups[g]; setone(sp); end
        else
            setone(isp)
        end
    else
        idx, _ = resolve_species(spfield, s.variant, s.species, s.coef)
        idx > 0 && setone(Int(idx))
    end
    return
end

# SNAGBRK (fmin.f:504, opt 10): per-species snag height-LOSS rates. Field 1 = species (SPDECD), fields
# 2/3 = YRS50 (years to lose 50% height) for HARD/SOFT snags, fields 4/5 = YRS30 (years to reach 30%)
# for HARD/SOFT. Converted to the 4 HTX coefficients FMSNGHT uses (fmin.f:538/546/557/566), with the SN
# constants HTR1=HTR2=0.01 and HTXSFT=2.0. Stored sparsely; empty ⇒ HTX=0 = no height loss (SN default).
function _snagbrk!(s::StandState, rec::KeywordRecord)
    fs = s.fire; fs === nothing && return
    p = fs.params
    HTR = 0.01f0; HTXSFT = 2f0
    yr(i) = rec.present[i] ? max(1, Int(trunc(Float64(rec.values[i])))) : 1   # INT, clamp ≥1 (fmin.f:526…)
    y50h = yr(2); y50s = yr(3); y30h = yr(4); y30s = yr(5)
    # HTX1/3: 50%-loss rate (>0.5·HTD regime); HTX2/4: 30%-loss rate in the <0.5·HTD regime (the 0.3/0.5 step,
    # over YRS30−YRS50 years, bumped to avoid div-0). Only fields actually given are written (else stay 0).
    htx1 = rec.present[2] ? (1f0 - 0.5f0^(1f0 / y50h)) / HTR : 0f0
    htx3 = rec.present[3] ? (1f0 - 0.5f0^(1f0 / y50s)) / (HTR * HTXSFT) : 0f0
    dh = y30h > y50h ? Float32(y30h - y50h) : 0.001f0
    ds = y30s > y50s ? Float32(y30s - y50s) : 0.001f0
    htx2 = rec.present[4] ? (1f0 - 0.6f0^(1f0 / dh)) / HTR : 0f0            # 0.3/0.5 = 0.6
    htx4 = rec.present[5] ? (1f0 - 0.6f0^(1f0 / ds)) / (HTR * HTXSFT) : 0f0
    htx = (htx1, htx2, htx3, htx4)
    spfield = strip(rec.fields[1])
    setone(sp) = @inbounds begin
        (1 <= sp <= length(s.species.class_codes[:, 1])) || return
        p.snag_htx[Int32(sp)] = htx
    end
    fnum = tryparse(Float64, spfield)
    if fnum !== nothing
        isp = round(Int, fnum)
        if isp == 0
            for sp in 1:length(s.species.class_codes[:, 1]); setone(sp); end
        elseif isp < 0
            g = -isp
            (1 <= g <= length(s.control.sp_groups)) || return
            for sp in s.control.sp_groups[g]; setone(sp); end
        else
            setone(isp)
        end
    else
        idx, _ = resolve_species(spfield, s.variant, s.species, s.coef)
        idx > 0 && setone(Int(idx))
    end
    return
end

# Lazily materialise the per-stand DKR override as a copy of the default decay-rate matrix, so a
# FUELMULT/FUELDCAY keyword modifies a private copy and fmcwd! reads it (an empty matrix ⇒ default).
function _ensure_dkr!(p::FFEParams)
    size(p.dkr, 1) == 11 || (p.dkr = copy(_FM_DKR))
    return p.dkr
end

# FUELMULT (fmin.f:1368, opt 29): multiply the total fuel decay rate of every size class by a per-decay-
# class multiplier (fields 1-4 = decay classes 1-4), capped at 1.0. A blank field leaves that class as-is.
function _fuelmult!(s::StandState, rec::KeywordRecord)
    fs = s.fire; fs === nothing && return
    dkr = _ensure_dkr!(fs.params); v = rec.values
    @inbounds for idec in 1:4
        rec.present[idec] || continue
        m = Float32(v[idec])
        for i in 1:11; dkr[i, idec] = min(dkr[i, idec] * m, 1f0); end
    end
    return
end

# FUELDCAY (fmin.f:806, opt 16): set the total decay rate for specific size classes of one decay class.
# Field 1 = decay class ID (clamp 1-4; ID≥5 ⇒ apply class 4's rates to ALL classes), 2 = litter (size 10),
# 3 = duff (11), 4-6 = size classes 1-3, 7 = size classes 4-9. Woody/litter (1-10) capped at 1.0.
function _fueldcay!(s::StandState, rec::KeywordRecord)
    fs = s.fire; fs === nothing && return
    rec.present[1] || return                         # no decay class specified ⇒ keyword ignored
    dkr = _ensure_dkr!(fs.params); v = rec.values
    id = Int(trunc(v[1])); idec = clamp(id, 1, 4)
    rec.present[2] && (dkr[10, idec] = Float32(v[2]))
    rec.present[3] && (dkr[11, idec] = Float32(v[3]))
    rec.present[4] && (dkr[1, idec]  = Float32(v[4]))
    rec.present[5] && (dkr[2, idec]  = Float32(v[5]))
    rec.present[6] && (dkr[3, idec]  = Float32(v[6]))
    rec.present[7] && (@inbounds for j in 4:9; dkr[j, idec] = Float32(v[7]); end)
    if id < 5
        @inbounds for i in 1:10; dkr[i, idec] = min(dkr[i, idec], 1f0); end
    else                                             # ID≥5: copy class-4 rates to all decay classes
        @inbounds for d2 in 1:4, j in 1:11; dkr[j, d2] = dkr[j, 4]; end
        @inbounds for d2 in 1:4, i in 1:10; dkr[i, d2] = min(dkr[i, d2], 1f0); end
    end
    return
end

# FUELINIT (fmin.f:1066, opt 21): initial HARD surface-fuel loadings (tons/ac), 12 params → STFUEL size
# classes (fmcba.f:321-342). −1 / blank = keep default. The PRMS(1) "<1"" field fills sizes 1+2 unless
# they are given explicitly. Stored in FFEParams.stfuel_hard (size 1:11, −1 = no override); fmcba! applies.
function _fuelinit!(s::StandState, rec::KeywordRecord)
    fs = s.fire; fs === nothing && return
    length(fs.params.stfuel_hard) == 11 || (fs.params.stfuel_hard = fill(-1f0, 11))
    sh = fs.params.stfuel_hard; v = rec.values; pr = rec.present
    pr[2]  && (sh[3]  = Float32(v[2]))    # 1-3"
    pr[3]  && (sh[4]  = Float32(v[3]))    # 3-6"
    pr[4]  && (sh[5]  = Float32(v[4]))    # 6-12"
    pr[5]  && (sh[6]  = Float32(v[5]))    # 12-20"
    pr[6]  && (sh[10] = Float32(v[6]))    # litter
    pr[7]  && (sh[11] = Float32(v[7]))    # duff
    pr[8]  && (sh[1]  = Float32(v[8]))    # <.25"
    pr[9]  && (sh[2]  = Float32(v[9]))    # .25-1"
    if pr[1]                              # <1" lumped — split into sizes 1 & 2 (fmcba.f:329-340)
        p1 = Float32(v[1])
        if !pr[8] && !pr[9]
            sh[1] = p1 * 0.5f0; sh[2] = p1 * 0.5f0
        elseif !pr[8] && pr[9]
            sh[1] = max(p1 - Float32(v[9]), 0f0)
        elseif pr[8] && !pr[9]
            sh[2] = max(p1 - Float32(v[8]), 0f0)
        end
    end
    pr[10] && (sh[7] = Float32(v[10]))    # 20-35"
    pr[11] && (sh[8] = Float32(v[11]))    # 35-50"
    pr[12] && (sh[9] = Float32(v[12]))    # >50"
    return
end

# FUELSOFT (fmin.f:2459, opt 53): initial SOFT/rotten surface-fuel loadings, 9 params → size classes 1-9
# directly (fmcba.f:354-362). Stored in FFEParams.stfuel_soft (size 1:11, −1 = no override).
function _fuelsoft!(s::StandState, rec::KeywordRecord)
    fs = s.fire; fs === nothing && return
    length(fs.params.stfuel_soft) == 11 || (fs.params.stfuel_soft = fill(-1f0, 11))
    ss = fs.params.stfuel_soft; v = rec.values; pr = rec.present
    @inbounds for i in 1:9
        pr[i] && (ss[i] = Float32(v[i]))
    end
    return
end

# Lazily materialise the per-stand PRDUFF override as an 11×4 matrix of the uniform 0.02 default, so
# DUFFPROD modifies a private copy and fmcwd! reads it (an empty matrix ⇒ the _FM_PRDUFF default).
function _ensure_prduff!(p::FFEParams)
    size(p.prduff, 1) == 11 || (p.prduff = fill(0.02f0, 11, 4))
    return p.prduff
end

# DUFFPROD (fmin.f:887, opt 17): proportion of decayed material that goes to duff (vs the air), per size
# class of one decay class. Field 1 = decay class ID (clamp 1-4; ID≥5 ⇒ class-4 values to all), 7 = all
# sizes 1-10, 2 = litter(10), 3-5 = sizes 1-3, 6 = sizes 4-9. Clamped to [0,1] (fmin.f:918-925).
function _duffprod!(s::StandState, rec::KeywordRecord)
    fs = s.fire; fs === nothing && return
    rec.present[1] || return                         # no decay class ⇒ ignored
    pd = _ensure_prduff!(fs.params); v = rec.values
    id = Int(trunc(v[1])); idec = clamp(id, 1, 4)
    rec.present[7] && (@inbounds for i in 1:10; pd[i, idec] = Float32(v[7]); end)
    rec.present[2] && (pd[10, idec] = Float32(v[2]))
    rec.present[3] && (pd[1, idec]  = Float32(v[3]))
    rec.present[4] && (pd[2, idec]  = Float32(v[4]))
    rec.present[5] && (pd[3, idec]  = Float32(v[5]))
    rec.present[6] && (@inbounds for j in 4:9; pd[j, idec] = Float32(v[6]); end)
    if id <= 4
        @inbounds for i in 1:10; pd[i, idec] = clamp(pd[i, idec], 0f0, 1f0); end
    else                                             # ID≥5: copy class-4 proportions to all decay classes
        @inbounds for d2 in 1:4, i in 1:10; pd[i, d2] = clamp(pd[i, 4], 0f0, 1f0); end
    end
    return
end

# DEFULMOD (fmin.f:1795, opt 39, act 2539): define/alter a fuel model. Field 2 = model#; fields 3-7 =
# PRMS(2-6) = dead 1/10/100-hr SAV, live SAV, dead 1-hr load; a SUPPLEMENTAL record (7×F10) = PRMS(7-13) =
# dead 10/100-hr load, live-woody load, depth, mext, live-herb SAV, live-herb load. −1/blank = keep the
# standard model's value. Stores the resolved (load, sav, depth, mext) override read by fuel_model_resolved.
function _defulmod!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    fs = s.fire; fs === nothing && return
    rec.present[2] || (read_raw_line!(kr); return)         # no model# ⇒ still consume the supplemental line
    model = Int(trunc(rec.values[2]))
    (1 <= model <= size(s.coef.ffe_fuel_models, 1)) || (read_raw_line!(kr); return)  # invalid model ⇒ skip (still consume the record)
    load, sav, depth, mext = standard_fuel_model(s.coef, model)
    load = copy(load); sav = copy(sav)
    p = fill(-1f0, 13)
    @inbounds for i in 2:6                                  # PRMS(2-6) from fields 3-7
        (rec.present[i+1] && rec.values[i+1] >= 0f0) && (p[i] = Float32(rec.values[i+1]))
    end
    line = read_raw_line!(kr)                               # supplemental record: PRMS(7-13), 7×10-char fields
    @inbounds for j in 0:6
        lo = j * 10 + 1; hi = min((j + 1) * 10, length(line))
        lo <= hi || continue
        val = tryparse(Float32, strip(line[lo:hi]))
        (val !== nothing && val >= 0f0) && (p[7+j] = val)
    end
    p[2] >= 0f0 && (sav[1, 1] = p[2]); p[3] >= 0f0 && (sav[1, 2] = p[3]); p[4] >= 0f0 && (sav[1, 3] = p[4])
    p[5] >= 0f0 && (sav[2, 1] = p[5]); p[12] >= 0f0 && (sav[2, 2] = p[12])
    p[6] >= 0f0 && (load[1, 1] = p[6]); p[7] >= 0f0 && (load[1, 2] = p[7]); p[8] >= 0f0 && (load[1, 3] = p[8])
    p[9] >= 0f0 && (load[2, 1] = p[9]); p[13] >= 0f0 && (load[2, 2] = p[13])
    p[10] >= 0f0 && (depth = p[10]); p[11] >= 0f0 && (mext = p[11])
    fs.defulmod[Int32(model)] = (load, sav, depth, mext)
    return
end

function kw_fmin!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    s.fire === nothing && (s.fire = FireState())
    fs = s.fire
    fs.active = true
    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"
            break
        elseif k == "SIMFIRE"                              # fire event + conditions (fmin.f:292-360)
            v = r.values
            # Field 1 = IDT (DATE or CYCLE), default 1 (fmin.f:309-310 IDT=1; LNOTBK(1)⇒ARRAY(1)). A value
            # ≤ MAXCYC is a 1-based CYCLE number, converted to that cycle's start year (opexpn.f:40-44
            # IY1.LE.MAXCYC ⇒ IY(IY1)); otherwise it is a calendar year. So a no-param SIMFIRE (s10_fire)
            # fires in cycle 1 (= the inventory year), not "never". Each SIMFIRE is its own OPNEW activity,
            # so >1 keyword schedules >1 fire (fire_repeat 2000+2020) — push each onto fire_schedule.
            idt = r.present[1] ? Int(nint(v[1])) : 1
            # FVS IDT cycle numbers are 1-based; jl `cycle_year_at` is 0-based ⇒ cycle IDT → index IDT-1.
            fyear = (1 <= idt <= MAXCYC) ? Int(cycle_year_at(s.control, idt - 1)) : idt
            # Resolve each condition with the FVS default (fmin.f:325-330 PRMS preset), overridden when present.
            swind    = r.present[2] ? Float32(v[2]) : 20f0                          # SWIND  = PRMS(1)
            # FMOIS = INT(PRMS(2)). FVS's FMMOIS only sets the moisture for codes 1..4; for any
            # other code (0, or out-of-range like the 9 in fire_fuel9) it is a NO-OP, leaving the
            # moisture at its last value — which, after the per-cycle PotFire MODERATE pass (FMOIS=3,
            # fmvinit.f:63-66), is dryness model 3. So an invalid code resolves to model 3, NOT a
            # clamp to the very-wet model 4. (Codes 1..4 use the matching FMMOIS table directly.)
            fmois    = r.present[3] ? (local fc = Int32(nint(v[3])); (Int32(1) <= fc <= Int32(4)) ? fc : Int32(3)) : Int32(1)
            atemp    = r.present[4] ? Float32(trunc(v[4])) : 70f0                   # ATEMP  = INT(PRMS(3))
            mortcode = r.present[5] ? clamp(Int32(nint(v[5])), Int32(0), Int32(1)) : Int32(1)   # MKODE = PRMS(4)
            psburn   = r.present[6] ? clamp(Float32(v[6]), 0f0, 100f0) : 100f0      # PSBURN = PRMS(5)
            burnseas = r.present[7] ? clamp(Int32(nint(v[7])), Int32(1), Int32(4)) : Int32(1)   # BURNSEAS = PRMS(6)
            push!(fs.fire_schedule, (Float32(fyear), swind, Float32(fmois), atemp,
                                     Float32(mortcode), psburn, Float32(burnseas)))
            # Keep the scalars in sync with the earliest pending fire for any legacy single-fire read.
            sort!(fs.fire_schedule; by = first)
            ev = fs.fire_schedule[1]
            fs.fire_year = Int32(ev[1]); fs.swind = ev[2]; fs.fmois = Int32(ev[3])
            fs.atemp = ev[4]; fs.mortcode = Int32(ev[5]); fs.psburn = ev[6]; fs.burnseas = Int32(ev[7])
        elseif k == "FLAMEADJ"                             # flame mult + crown fraction (fmburn.f:337-347)
            v = r.values
            r.present[2] && (fs.flmult = Float32(v[2]))                            # FLMULT = FPRMS(1)
            r.present[4] && (fs.crburn = v[4] > -1f0 ? Float32(v[4]) * 0.01f0 : 0f0)  # CRBURN = FPRMS(3)·.01
        elseif k == "CARBREPT"                             # request the FFE Stand Carbon Report (fmcrbout.f)
            s.control.carbon_report_on = true
        elseif k == "POTFIRE" || k == "POTFLAME"           # request the Potential Fire report (fmpofl.f)
            s.control.potfire_report_on = true
        elseif k == "CARBCALC"                             # FLD1 method 0=FFE*/1=Jenkins, FLD2 units 0=US-t/ac*/1=t-ha/2=t-ac
            r.present[1] && (s.control.carbon_method = Int32(clamp(nint(r.values[1]), 0, 1)))
            r.present[2] && (s.control.carbon_units  = Int32(clamp(nint(r.values[2]), 0, 2)))
        elseif k == "FMORTMLT"                             # fire-caused mortality multiplier (fmin.f:2506)
            kw_mult!(s, r, :fmort)
        elseif k == "SNAGFALL"                             # per-species snag fall rates (fmin.f:464, opt 9)
            _snagfall!(s, r)
        elseif k == "SNAGDCAY"                             # per-species snag DECAYX (fmin.f:633, opt 11)
            _snagdcay!(s, r)
        elseif k == "SNAGBRK"                              # per-species snag height-loss HTX (fmin.f:504, opt 10)
            _snagbrk!(s, r)
        elseif k == "FUELMULT"                             # fuel decay-rate multiplier (fmin.f:1368, opt 29)
            _fuelmult!(s, r)
        elseif k == "FUELDCAY"                             # fuel decay-rate set (fmin.f:806, opt 16)
            _fueldcay!(s, r)
        elseif k == "DUFFPROD"                             # proportion of decay → duff (fmin.f:887, opt 17)
            _duffprod!(s, r)
        elseif k == "SNAGPSFT"                             # per-species initial-soft snag fraction (fmin.f:1683, opt 37)
            # field 1 = species (SPDECD), field 2 = PSOFT (proportion soft at creation, clamp [0,1]).
            if r.present[2]
                ps = clamp(Float32(r.values[2]), 0f0, 1f0)
                sel = species_selector(s, r.fields[1])
                if sel > 0
                    fs.params.psoft_ovr[Int32(sel)] = ps
                elseif sel == 0
                    @inbounds for sp in 1:length(s.species.class_codes[:, 1]); fs.params.psoft_ovr[Int32(sp)] = ps; end
                elseif -sel <= length(s.control.sp_groups)
                    for sp in s.control.sp_groups[-sel]; fs.params.psoft_ovr[Int32(sp)] = ps; end
                end
            end
        elseif k == "SALVAGE"                              # remove snags (fmin.f:993, opt 20, act 2520)
            # field 1 = date, 2-7 = min DBH / max DBH / max age / OKSOFT (0=all,1=hard,2=soft) / PROP
            # (fraction removed) / PROPLV (proportion left → down-wood). Defaults match fmin.f:1023-1028.
            v = r.values
            date  = r.present[1] ? nint(v[1]) : Int32(1)
            mindb = r.present[2] ? Float32(v[2]) : 0f0
            maxdb = r.present[3] ? Float32(v[3]) : 999f0
            maxag = r.present[4] ? Float32(v[4]) : 5f0
            oksft = r.present[5] ? Float32(v[5]) : 1f0
            prop  = r.present[6] ? Float32(v[6]) : 0.9f0
            proplv = r.present[7] ? Float32(v[7]) : 0f0
            push!(s.control.schedule,
                  ScheduledActivity(date, Int32(2520), (mindb, maxdb, maxag, oksft, prop, proplv)))
        elseif k == "DEFULMOD"                             # define/alter a fuel model (fmin.f:1795, opt 39, act 2539)
            _defulmod!(s, r, kr)
        elseif k == "FUELMODL"                             # force standard fuel models (fmin.f:1717, opt 38, act 2538)
            # field 1 = date, then up to 3 (model#, weight) pairs (fields 2-7). Valid models 0<m≤MXDFMD;
            # a blank/≤0 weight defaults to 1; weights normalized to sum 1 (fmin.f:1767). No valid model ⇒
            # auto-selection (no override stored).
            v = r.values
            date = r.present[1] ? nint(v[1]) : Int32(1)
            pairs = Tuple{Int32,Float32}[]
            for p in 0:2
                mi = 2 + 2p; wi = 3 + 2p
                r.present[mi] || continue
                m = Int(trunc(v[mi]))
                (m > 0 && m <= 53) || continue            # MXDFMD ≈ 53 standard models
                w = (r.present[wi] && v[wi] > 0f0) ? Float32(v[wi]) : 1f0
                push!(pairs, (Int32(m), w))
            end
            if !isempty(pairs)
                tot = sum(p[2] for p in pairs)
                tot > 0f0 && (pairs = [(p[1], p[2] / tot) for p in pairs])
                push!(fs.fuelmodl, (date, pairs))
            end
        elseif k == "PILEBURN"                             # jackpot/pile burn (fmin.f:1161, opt 23, act 2523)
            # field 1 = date, 2 = type, 3 = AFFECT %, 4 = ATREAT %, 5 = FULCON %, 6 = TRMORT %.
            # Type-1 defaults (fmin.f:1196): 70 / 10 / 80 / 0.
            v = r.values
            date = r.present[1] ? nint(v[1]) : Int32(1)
            typ  = r.present[2] ? Float32(v[2]) : 1f0
            aff  = r.present[3] ? Float32(v[3]) : 70f0
            atr  = r.present[4] ? Float32(v[4]) : 10f0
            ful  = r.present[5] ? Float32(v[5]) : 80f0
            trm  = r.present[6] ? Float32(v[6]) : 0f0
            push!(s.control.schedule, ScheduledActivity(date, Int32(2523), (typ, aff, atr, ful, trm, 0f0)))
        elseif k == "FUELMOVE"                             # transfer fuel between size pools (fmin.f:1515, opt 34, act 2530)
            # field 1 = date, 2 = FROM size (0-11), 3 = TO size, 4 = amount, 5 = proportion, 6 = leave (Z),
            # 7 = target-final (Q). Defaults match fmin.f (6/11/0/0/9999/0).
            v = r.values
            date = r.present[1] ? nint(v[1]) : Int32(1)
            frm  = r.present[2] ? Float32(v[2]) : 6f0
            to   = r.present[3] ? Float32(v[3]) : 11f0
            amt  = r.present[4] ? Float32(v[4]) : 0f0
            prop = r.present[5] ? Float32(v[5]) : 0f0
            leav = r.present[6] ? Float32(v[6]) : 9999f0
            q    = r.present[7] ? Float32(v[7]) : 0f0
            push!(s.control.schedule, ScheduledActivity(date, Int32(2530), (frm, to, amt, prop, leav, q)))
        elseif k == "SALVSP"                               # salvage species cut/leave list (fmin.f:149, opt 1, act 2501)
            # field 1 = date, 2 = species (SPDECD: 0=all/idx/−group), 3 = flag (<1 cut-list, ≥1 leave-list).
            date = r.present[1] ? nint(r.values[1]) : Int32(1)
            sel  = species_selector(s, length(r.fields) >= 2 ? r.fields[2] : "")
            mode = (r.present[3] && r.values[3] >= 1f0) ? 1f0 : 0f0
            push!(s.control.schedule,
                  ScheduledActivity(date, Int32(2501), (Float32(sel), mode, 0f0, 0f0, 0f0, 0f0)))
        elseif k == "FUELTRET"                             # fuel-treatment depth adjustment (fmin.f:1264, opt 25, act 2525)
            # field 1 = date, 2 = treatment type (0-2), 3 = harvest type (1-3), 4 = depth mult (−1 = use the
            # DPMULT table). DPMULT(HARTYP, FTREAT+1): FTREAT 0 → 1.0/1.3/1.6 by harvest, FTREAT 1 → 0.83, 2 → 0.75.
            v = r.values
            date = r.present[1] ? nint(v[1]) : Int32(1)
            ftreat = r.present[2] ? clamp(Int(trunc(v[2])), 0, 2) : 0
            hartyp = r.present[3] ? clamp(Int(trunc(v[3])), 1, 3) : 1
            dpmod = (r.present[4] && v[4] >= 0f0) ? Float32(v[4]) :
                    (ftreat == 0 ? (1f0, 1.3f0, 1.6f0)[hartyp] : ftreat == 1 ? 0.83f0 : 0.75f0)
            push!(fs.fueltret, (date, dpmod))
        elseif k == "FIRECALC"                             # fire-calc method (fmin.f:2293, opt 49, act 2549)
            # SN default IFLOGIC=0 (OLD FM logic, fminit.f:824) — exactly jl's FMCFMD path. The USAV/UBD/heat
            # overrides (fields 4-9) apply ONLY to method 1 (new FM logic) / 2 (modelled loads → FM89), which
            # are alternative fire-behavior models not ported. So method 0 is a faithful no-op; warn otherwise.
            (r.present[2] && Int(trunc(r.values[2])) != 0) &&
                @warn "FIRECALC method 1/2 (new FM logic / modelled loads) not ported — defaulting to the faithful old-FM-logic (method 0) path"
        elseif k == "DROUGHT"
            # SN no-op: DROUGHT sets IDRYB/IDRYE (drought years), but those affect the fuel model only in the
            # UT/CR/LS variants — "not used in OZ-FFE" (the Southern FFE; fmvinit.f:1113). Recognized.
        elseif k == "CANCALC"
            # SN no-op: CANCALC sets canopy base-height / bulk-density options for the CROWN-fire model
            # (FMCFIR), which the SN variant does not run (potential_fire/fmburn skip crown fire). Recognized.
        elseif k == "SOILHEAT"
            # report-only: SOILHEAT requests the soil-heating report when a fire occurs; jl emits no
            # soil-heating report, so this is a recognized no-op (like the other FFE report keywords).
        elseif k == "FUELPOOL"                             # per-species fuel decay class (fmin.f:967, opt 19)
            # field 1 = species (SPDECD: alpha/FIA/0=all/−group), field 2 = decay class 1-4.
            if r.present[2]
                idec = Int(trunc(r.values[2]))
                if 1 <= idec <= 4
                    sel = species_selector(s, r.fields[1])
                    if sel > 0
                        fs.params.dkrcls_ovr[Int32(sel)] = Int32(idec)
                    elseif sel == 0                        # 0 ⇒ all species
                        @inbounds for sp in 1:length(s.species.class_codes[:, 1])
                            fs.params.dkrcls_ovr[Int32(sp)] = Int32(idec)
                        end
                    elseif -sel <= length(s.control.sp_groups)  # −g ⇒ SPGROUP g
                        for sp in s.control.sp_groups[-sel]
                            fs.params.dkrcls_ovr[Int32(sp)] = Int32(idec)
                        end
                    end
                end
            end
        elseif k == "FUELINIT"                             # initial hard fuel loadings (fmin.f:1066, opt 21)
            _fuelinit!(s, r)
        elseif k == "FUELSOFT"                             # initial soft fuel loadings (fmin.f:2459, opt 53)
            _fuelsoft!(s, r)
        elseif k == "SNAGINIT"                             # add user snags (fmin.f:1119, opt 22)
            # field 1 = species (SPDECD), 2 = DBH at death, 3 = ht at death (HTDEAD, taper+death-vol),
            # 4 = CURRENT ht (HTIH — the snag's present top that truncates the fall cone in CWD1), 5 = age,
            # 6 = density stems/ac. 0/−999 species ⇒ ignored (fmin.f:1142). A pre-broken SNAGINIT snag
            # (`… 50 40 …` ⇒ HTDEAD 50 / HTIH 40) falls as the fat LOWER 40ft of a 50ft cone ⇒ more sz5.
            sp = species_selector(s, r.fields[1])
            if sp > 0
                v = r.values
                d   = r.present[2] ? Float32(v[2]) : -1f0
                htd = r.present[3] ? Float32(v[3]) : -1f0
                htc = r.present[4] ? Float32(v[4]) : -1f0
                age = r.present[5] ? Float32(v[5]) : -1f0
                den = r.present[6] ? Float32(v[6]) : -1f0
                push!(fs.snaginit, (Float32(sp), d, htd, htc, age, den))
            end
        elseif k == "POTFMOIS"                             # PotFire moisture (fmin.f:1391, opt 30)
            # field 1 = IFIRE (1=severe/2=moderate), fields 2-8 = 7 moisture % for that scenario; a blank
            # field uses the FMMOIS default for the scenario, blank herb (8) ⇒ the (resolved) woody value.
            v = r.values
            ifire = r.present[1] ? clamp(Int(trunc(v[1])), 1, 2) : 1
            def = fuel_moisture(ifire == 1 ? 1 : 3)        # severe→model 1, moderate→model 3
            dpct = (def[1,1], def[1,2], def[1,3], def[1,4], def[1,5], def[2,1], def[2,2]) .* 100f0
            m = Float32[r.present[i+1] ? Float32(v[i+1]) : dpct[i] for i in 1:7]
            r.present[8] || (m[7] = m[6])                  # herb defaults to woody (fmin.f:1421)
            fs.params.potf[ifire].mois = (m[1],m[2],m[3],m[4],m[5],m[6],m[7])
        elseif k == "POTFWIND"                             # PotFire wind (fmin.f:1658, opt 35) — sev/mod
            r.present[1] && (fs.params.potf[1].wind = Float32(r.values[1]))
            r.present[2] && (fs.params.potf[2].wind = Float32(r.values[2]))
        elseif k == "POTFTEMP"                             # PotFire temperature (fmin.f:1671, opt 36)
            r.present[1] && (fs.params.potf[1].temp = Float32(r.values[1]))
            r.present[2] && (fs.params.potf[2].temp = Float32(r.values[2]))
        elseif k == "POTFSEAS"                             # PotFire season (fmin.f:2000, opt 41)
            r.present[1] && (fs.params.potf[1].season = Int32(trunc(r.values[1])))
            r.present[2] && (fs.params.potf[2].season = Int32(trunc(r.values[2])))
        elseif k == "POTFPAB"                              # PotFire % area burned (fmin.f:2012, opt 42)
            r.present[1] && (fs.params.potf[1].pab = Float32(r.values[1]))
            r.present[2] && (fs.params.potf[2].pab = Float32(r.values[2]))
        elseif k == "MOISTURE"                             # fuel-moisture override (fmin.f:237, opt 5)
            v = r.values
            idt = r.present[1] ? Int32(trunc(v[1])) : Int32(1)            # IDT = date/cycle (default 1)
            # 7 fuel-moisture % (1hr/10hr/100hr/3+/duff/live-woody/live-herb); blank → 0, except a blank
            # live-herb (field 8) defaults to the live-woody value (fmin.f:277-279).
            pr = ntuple(i -> r.present[i+1] ? Float32(v[i+1]) : 0f0, 7)
            (!r.present[8]) && (pr = (pr[1], pr[2], pr[3], pr[4], pr[5], pr[6], pr[6]))
            push!(fs.moisture_ovr, (idt, pr))
        elseif k == "SNAGPBN"                              # post-burn snag-fall params (fmin.f:1233, opt 24)
            v = r.values; p = fs.params
            r.present[1] && (p.pb_soft = clamp(Float32(v[1]), 0f0, 1f0))   # PBSOFT, clamp [0,1]
            r.present[2] && (p.pb_smal = clamp(Float32(v[2]), 0f0, 1f0))   # PBSMAL, clamp [0,1]
            r.present[3] && (p.pb_time = max(Float32(v[3]), 1f0))          # PBTIME, min 1
            r.present[4] && (p.pb_size = max(Float32(v[4]), 0f0))          # PBSIZE, min 0
            r.present[5] && (p.pb_scor = max(Float32(v[5]), 0f0))          # PBSCOR, min 0
        elseif k in _FFE_REPORT_KEYWORDS
            # report-only FFE keywords (BURNREPT/FUELOUT/SNAGSUM/…): the text reports aren't emitted; the
            # equivalent data is available via the DBS path. Recognized, intentionally a no-op here.
        else
            # Transparency: jl handles only a subset of the ~53 FMIN/FFE keywords (fmin.f TABLE). The rest
            # are model/override keywords (MOISTURE/SNAGPBN/FUELMODL/FUELINIT/…) that change the simulation
            # when used — silently skipping them would produce wrong results with no signal. Warn so the run
            # is not silently unfaithful. (See docs/audit/INDEX.md "SYSTEMATIC GAP — FMIN handler".)
            @warn "FMIN/FFE keyword not yet ported — IGNORED (using defaults; result may diverge from FVS)" keyword=k
        end
    end
    return
end

"""
    kw_econ!(s, rec, kr)

Parse the ECON economic-analysis keyword block (read to END). Captures the management
costs (ANNUCST), variable harvest costs by DBH class (HRVVRCST) and harvest revenues by
species + DBH (HRVRVN) into `EconState` for the discounting core; other ECON keywords
are recognized but not yet ported.
"""
function kw_econ!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    s.econ === nothing && (s.econ = EconState())
    ec = s.econ
    ec.active = true
    while true
        r = read_keyword!(kr)
        (r.status == KW_EOF || r.status == KW_STOP) && break
        k = strip(r.name)
        isempty(k) && continue
        if k == "END"
            break
        elseif k == "STRTECON"                             # ecin.f: field1=start year/delay, field2=DISCOUNT RATE (%),
            # field3=known SEV, field4=compute-SEV flag. jl honors the DISCOUNT RATE (the audit GAP); the start-year
            # delay and SEV are not modeled. Rate is a PERCENT (eccalc.f:91 rate=discountRate/100). Default 0 (ecinit.f:15).
            r.present[2] && (ec.discount_rate = Float32(r.values[2]) / 100f0)
        elseif k == "ANNUCST"                              # annual management cost ($/ac/yr)
            r.present[1] && (ec.ann_cost += Float32(r.values[1]))
        elseif k == "HRVVRCST"                             # variable harvest cost by DBH class
            r.present[1] || continue
            hi = (r.present[4] && r.values[4] > 0f0) ? Float32(r.values[4]) : 999f0
            push!(ec.hrv_cost, EconCostRev(Float32(r.values[1]), Int32(nint(r.values[2])),
                                           Float32(r.values[3]), hi, Int32(0)))
        elseif k == "HRVRVN"                               # harvest revenue by species + DBH (field 4 = species)
            r.present[1] || continue
            # units (field 2): 1=TPA, 2=BF_1000(whole-tree bd ft), 3=FT3_100(whole-tree cu ft), 4=BF_1000_LOG,
            # 5=FT3_100_LOG (ECNCOM.F77:19). The LOG units 4/5 are REPORT-ONLY (echarv.f populates the
            # FVS_EconHarvestValue detail table; they do NOT flow to FVS_EconSummary.Revenue/PNV — confirmed 0
            # live even with periods after the harvest). Unit 4 (BF_1000_LOG) is now PORTED: the R9LOGS full-stem
            # log-bucking + DIB-grading (r8clark_vol.jl `_r8_scribner_bf_by_dib` → econ.jl `accrue_log_grade!`)
            # reproduces FVS_EconHarvestValue bit-exact (econ_strtecon.key HRVRVN 300 4 10.0 ALL + THINSDI 2000:
            # SM=16bf/$5, HI=9/$3, AB=41/$12, SK=5/$1). Unit 5 (FT3_100_LOG) is now ALSO ported: the per-log
            # CUBIC bucking R9LGCFT (r8clark_vol.jl `_r8_cuft_by_dib` → econ.jl `accrue_log_grade_cuft!`) fills
            # the Ft3_Removed/Ft3_Value columns. Both units flow through `econ_harvest_value_rows`.
            code = length(r.fields) >= 4 ? strip(uppercase(r.fields[4])) : "ALL"
            sp = code == "ALL" || isempty(code) ? Int32(0) :
                 Int32(first(resolve_species(code, s.variant, s.species, s.coef)))
            sp < 0 && (sp = Int32(0))
            push!(ec.hrv_rev, EconCostRev(Float32(r.values[1]), Int32(nint(r.values[2])),
                                          Float32(r.values[3]), 999f0, sp))
        end
    end
    return
end

# OPTION 15 — STDIDENT (initre.f:862): stand id from the next raw line.
function kw_stdident!(s::StandState, kr::KeywordReader)
    record = rpad(read_raw_line!(kr), 250)[1:250]
    i1 = findfirst(c -> c != ' ', record[1:26])
    if i1 === nothing
        s.plot.stand_id = rpad(" ", 26)
        s.control.title = rpad("", 72)
        return
    end
    i2 = something(findnext(==(' '), record, i1), 27)
    i2 = min(i2, 27)
    nplt = record[i1:min(i2, 26)]
    s.plot.stand_id = length(nplt) <= 26 ? rpad(nplt, 26) : nplt[1:26]
    rest = i2 + 1 <= length(record) ? record[i2+1:end] : ""
    s.control.title = length(rest) >= 72 ? rest[1:72] : rpad(rest, 72)
    return
end

# OPTION 5 — TREEFMT (initre.f:646): custom tree-record format, 2×80 raw lines.
function kw_treefmt!(s::StandState, kr::KeywordReader)
    line1 = read_raw_line!(kr)
    line2 = read_raw_line!(kr)
    s.control.tree_format = rpad(line1, 80)[1:80] * rpad(line2, 80)[1:80]
    return
end

# OPTION — IF / THEN / ENDIF (event monitor, evmon.f). Reads the condition expression
# (raw lines up to THEN), then the activity keywords up to ENDIF, into a
# ConditionalActivity whose acts fire only in cycles where the condition is true.
# The IF-block activities carry NO date (field 1 blank) — their year is filled in at
# fire time — so they parse with the same (year=v[1]=0, params=v[2:7]) layout.
# OPTION 33 — COMPUTE (vbase/initre.f:1266 → EVUSRV): an event-monitor user-variable block.
# Field 1 = start year (IDT, default 1). The following raw lines are `NAME = expression`
# assignments, re-evaluated each cycle from the start year (a later def may use an earlier one);
# the block ends at END. The values become variables readable by IF/THEN conditions.
function kw_compute!(s::StandState, rec::KeywordRecord, kr::KeywordReader)
    date = rec.present[1] ? nint(rec.values[1]) : Int32(1)
    while !eof(kr.io)
        t = strip(read_raw_line!(kr))
        isempty(t) && continue
        uppercase(t) == "END" && break
        eq = findfirst('=', t)
        eq === nothing && continue
        name = uppercase(strip(t[1:eq-1])); expr = strip(t[eq+1:end])
        (isempty(name) || isempty(expr)) && continue
        push!(s.control.compute_defs, (date, String(name), parse_event_condition(String(expr))))
    end
    return
end

function kw_if!(s::StandState, kr::KeywordReader)
    cond = ""
    while true
        line = strip(read_raw_line!(kr))
        (isempty(line) || startswith(uppercase(line), "THEN")) && (isempty(line) ? continue : break)
        cond *= " " * line
        uppercase(line) == "THEN" && break
    end
    acts = ScheduledActivity[]
    while true
        rec = read_keyword!(kr)
        (rec.status == KW_EOF || rec.status == KW_STOP) && break
        kw = strip(rec.name); isempty(kw) && continue
        uppercase(kw) == "ENDIF" && break
        v = rec.values
        if haskey(_THIN_ICFLAG, kw)
            push!(acts, ScheduledActivity(nint(v[1]), _THIN_ICFLAG[kw],
                                          ntuple(i -> Float32(v[i + 1]), 6)))
        elseif kw == "SPECPREF"
            push!(acts, ScheduledActivity(nint(v[1]), Int32(201),
                                          ntuple(i -> Float32(v[i + 1]), 6)))
        end
        # unknown keywords inside IF (e.g. COMPUTE) are skipped for now
    end
    push!(s.control.conditionals,
          ConditionalActivity(parse_event_condition(cond), acts, strip(cond)))
    return
end

"""
    process_keywords!(state, kr, base_path) -> Symbol

Process keyword records from `kr` until PROCESS / STOP / EOF. `base_path` is the
keyword file path with the extension stripped (used to locate the `.tre` file for
TREEDATA). Returns the terminating reason (:process, :stop, :eof).
"""
function process_keywords!(s::StandState, kr::KeywordReader, base_path::AbstractString)
    trees_loaded = false   # an explicit TREEDATA was processed
    notrees      = false   # NOTREES suppresses the default tree read
    nkw          = 0       # real keywords seen (0 ⇒ bare STOP/EOF, not a stand)
    # INITRE end (initre.f:334-336): if no TREEDATA keyword ran, read the tree file once
    # anyway — so a stand without TREEDATA (e.g. snt01's FFE stand, which only REWINDs the
    # shared data unit) still gets its trees. NOTREES (bare stands) suppresses this, and a
    # bare terminator (STOP/EOF with no keywords) is not a stand at all (nkw==0).
    finish(reason) = (nkw > 0 && !trees_loaded && !notrees &&
                      load_trees!(s, base_path * ".tre"); reason)
    while true
        rec = read_keyword!(kr)
        rec.status == KW_EOF && return finish(:eof)
        rec.status == KW_STOP && return finish(:stop)
        kw = strip(rec.name)
        isempty(kw) && continue                      # blank-line record
        nkw += 1
        if     kw == "DESIGN";   kw_design!(s, rec)
        elseif kw == "TFIXAREA"; kw_tfixarea!(s, rec)      # total fixed plot area (notre.f:45)
        elseif kw == "CUTEFF";   kw_cuteff!(s, rec)        # default cut/affect proportion EFF (initre.f:5400)
        elseif kw == "NUMCYCLE"; kw_numcycle!(s, rec)
        elseif kw == "TIMEINT";  kw_timeint!(s, rec)      # cycle length (period); default 5
        elseif kw == "GROWTH";   kw_growth!(s, rec)        # input growth-data type codes + measurement periods (initre.f:2300)
        elseif kw == "INVYEAR";  kw_invyear!(s, rec)
        elseif kw == "SITECODE"; kw_sitecode!(s, rec)
        elseif kw == "STDINFO";  kw_stdinfo!(s, rec)
        elseif kw == "MANAGED";  kw_managed!(s, rec)       # managed-stand flag → DGF kplant term (dgf.f:179)
        elseif kw == "BAMAX";    kw_bamax!(s, rec)         # max basal area → SDImax override (initre.f:6800)
        elseif kw == "SDIMAX";   kw_sdimax!(s, rec)        # per-species SDImax + PMSDIL/PMSDIU (initre.f:3072)
        elseif kw == "RANNSEED"; kw_rannseed!(s, rec)      # reseed the main RNG stream (initre.f:6300)
        elseif kw == "DGSTDEV";  kw_dgstdev!(s, rec)       # DGSD bound on stochastic DG variation (initre.f:5900)
        elseif kw == "NOCALIB";  kw_nocalib!(s, rec)       # disable DG self-calibration per species (initre.f:5800)
        elseif kw == "SERLCORR"; kw_serlcorr!(s, rec)      # ARMA(1,1) DGSCOR phi/theta (initre.f:9300)
        elseif kw == "READCORD"; kw_readcord!(s, kr)       # read large-tree DG correction COR2 → DGCON (initre.f:5600)
        elseif kw == "REUSCORD"; kw_reuscord!(s)           # reuse prior COR2 (initre.f:5700)
        elseif kw == "READCORH"; kw_readcorh!(s, kr)       # read large-tree HTG correction HCOR2 → HTCON (initre.f:6900)
        elseif kw == "REUSCORH"; kw_reuscorh!(s)           # reuse prior HCOR2 (initre.f:7000)
        elseif kw == "READCORR"; kw_readcorr!(s, kr)       # read small-tree HTG correction RCOR2 → RHCON (initre.f:7500)
        elseif kw == "REUSCORR"; kw_reuscorr!(s)           # reuse prior RCOR2 (initre.f:7600)
        elseif kw == "RESETAGE"; kw_resetage!(s, rec)      # rebase stand age at a date (resage.f act 443)
        elseif kw == "SDICALC";  kw_sdicalc!(s, rec)       # SDI method + thresholds (LZEIDE drives report+mortality, initre.f:14000)
        elseif kw == "COMPRESS"; kw_compress!(s, rec)      # schedule record compression (initre.f:8000; algorithm TODO)
        elseif kw == "STDIDENT"; kw_stdident!(s, kr)
        elseif kw == "TREEFMT";  kw_treefmt!(s, kr)
        elseif kw == "IF";       kw_if!(s, kr)
        elseif kw == "ESTAB";    kw_estab!(s, rec, kr)
        elseif kw == "DATABASE"; kw_database!(s, rec, kr)  # DBS output block (DSNOUT/SUMMARY → SQLite)
        elseif kw == "FMIN";     kw_fmin!(s, rec, kr)      # Fire & Fuels Extension block (SIMFIRE/FLAMEADJ)
        elseif kw == "ECON";     kw_econ!(s, rec, kr)      # ECON economic-analysis block (ANNUCST/HRVVRCST/HRVRVN)
        # SPROUT/NOSPROUT are establishment-extension sub-keywords (read by ESIN inside an
        # ESTAB…END block, esin.f opt 26/27) — NOT top-level base keywords (the Fortran base
        # processor rejects a standalone SPROUT with "INVALID KEYWORD"). Handled in kw_estab!.
        elseif kw == "TREEDATA"; load_trees!(s, base_path * ".tre"); trees_loaded = true
        elseif kw == "NOTREES";  notrees = true       # bare stand — no tree-data read
        elseif kw == "NOAUTOES"; s.control.lsprut = false  # disable automatic establishment (incl. auto stump-sprouting); estab.f LFLAG
        elseif kw == "NOSPROUT"; s.control.lsprut = false  # disable stump sprouting (esin.f opt 27), standalone form
        elseif kw == "THINQFA"; kw_thinqfa!(s, rec, kr)   # 2-record keyword
        elseif kw == "SPGROUP"; kw_spgroup!(s, rec, kr)   # species group: name + next-record species list
        elseif kw == "COMPUTE"; kw_compute!(s, rec, kr)   # event-monitor variable block (NAME = expr … END)
        elseif kw == "NOTRIPLE"; s.control.icl4 = Int32(0)            # disable record tripling (initre.f:5500)
        elseif kw == "NUMTRIP";  rec.present[1] && (s.control.icl4 = Int32(nint(rec.values[1])))  # set ICL4 (initre.f:2709)
        elseif haskey(_THIN_ICFLAG, kw); kw_thin!(s, rec, _THIN_ICFLAG[kw])
        elseif kw == "SPECPREF"; kw_thin!(s, rec, Int32(201))   # cut modifier: species preference
        elseif kw == "MINHARV";  kw_thin!(s, rec, Int32(200))   # cut modifier: minimum-harvest thresholds
        elseif kw == "SPLEAVE";  kw_thin!(s, rec, Int32(206))   # cut modifier: leave named species
        elseif kw == "TCONDMLT"; kw_thin!(s, rec, Int32(202))   # cut modifier: tree-condition weight (TCWT)
        elseif kw == "LEAVESP";  kw_thin!(s, rec, Int32(206))   # alias for SPLEAVE
        elseif kw == "YARDLOSS"                                 # yarding loss. FVS fields (initre.f:3637-45):
            # field1 = DATE, field2 = PRLOST, field3 = PRDSNG, field4 = PRCRWN. jl previously read PRLOST from
            # field1 (the DATE!) ⇒ YARDLOSS was silently INACTIVE (prlost=0). PRLOST = proportion of the
            # harvest lost in yarding (left on site); of that LOSS, PRDSNG is downed + (1−PRDSNG) standing snags.
            rec.present[2] && (s.control.yardloss_prlost = clamp(Float32(rec.values[2]), 0f0, 1f0))
            rec.present[3] && (s.control.yardloss_prdsng = clamp(Float32(rec.values[3]), 0f0, 1f0))
        elseif kw == "SALVAGE"                                  # ABANDONED in Fortran (cuts.f:103) — recognized
                                                                # no-op so the keyword doesn't fall through silently
        elseif kw == "SETPTHIN"; kw_thin!(s, rec, Int32(248))   # point-thin prescription (point, metric)
        elseif kw == "THINPT";   kw_thin!(s, rec, Int32(15))    # point thin (residual + class + dir)
        elseif kw == "BAIMULT";  kw_mult!(s, rec, :bai)   # diameter-growth (BA-increment) multiplier
        elseif kw == "HTGMULT";  kw_mult!(s, rec, :htg)   # height-growth multiplier
        elseif kw == "MORTMULT"; kw_mult!(s, rec, :mort)  # mortality-rate multiplier
        elseif kw == "REGHMULT"; kw_mult!(s, rec, :regh)  # regen height-growth multiplier
        elseif kw == "REGDMULT"; kw_mult!(s, rec, :regd)  # regen diameter-growth multiplier
        elseif kw == "TREESZCP"; kw_treeszcp!(s, rec)     # per-species size cap (SIZCAP)
        elseif kw == "FIXDG";    kw_mult!(s, rec, :fixdg)  # one-shot DG scaler (grincr.f:451)
        elseif kw == "FIXHTG";   kw_mult!(s, rec, :fixhtg) # one-shot HTG scaler (grincr.f:492)
        elseif kw == "HTGSTOP";  kw_htgstp!(s, rec, 110)   # top-damage: scale height growth
        elseif kw == "TOPKILL";  kw_htgstp!(s, rec, 111)   # top-damage: top-kill (htgstp.f)
        elseif kw == "FIXMORT";  kw_fixmort!(s, rec)       # forced-mortality override (morts.f:781)
        elseif kw == "VOLUME";   kw_volume!(s, rec)        # cubic merch-standard override (volkey.f)
        elseif kw == "BFVOLUME"; kw_bfvolume!(s, rec)      # board-foot merch-standard override (volkey.f)
        elseif kw == "MCDEFECT"; kw_mcdefect!(s, rec)      # cubic-volume defect curve (sdefet.f/vols.f)
        elseif kw == "BFDEFECT"; kw_bfdefect!(s, rec)      # board-foot defect curve (sdefet.f/vols.f)
        elseif kw == "VOLEQNUM"; kw_voleqnum!(s, rec)      # cubic volume-equation override (initre.f:5061)
        elseif kw == "FERTILIZ"; kw_fertiliz!(s, rec)      # fertilizer growth response (ffin.f/ffert.f)
        elseif kw == "MCFDLN";   kw_mcfdln!(s, rec)        # cubic form-model coefs CFLA0/CFLA1 (sdefln.f)
        elseif kw == "BFFDLN";   kw_bffdln!(s, rec)        # board form-model coefs BFLA0/BFLA1 (sdefln.f)
        elseif kw == "CRNMULT";  kw_mult!(s, rec, :crn)    # crown-ratio-change multiplier (crown.f:319)
        elseif kw == "FMORTMLT"; kw_mult!(s, rec, :fmort)  # FFE fire-caused mortality multiplier (fmeff.f:340)
        elseif kw == "CYCLEAT";  kw_cycleat!(s, rec)       # extra cycle-boundary year (initre.f opt 134)
        elseif kw == "SETSITE";  kw_setsite!(s, rec)       # scheduled mid-run site-index/BAMAX/SDImax change (act 120)
        elseif kw == "CUTLIST";  s.control.dbs_cutlist = true  # emit the FVS_CutList DBS table (dbscuts.f, ICUTLIST)
        elseif kw == "STRCLASS"; kw_strclass!(s, rec)      # activate SSTAGE structural-stage classification (ksstag.f)
        elseif kw == "CARBREPT"; kw_carbrept!(s, rec)      # request the FFE Stand Carbon Report (fmcrbout.f)
        elseif kw == "CARBCALC"; kw_carbcalc!(s, rec)      # carbon method 0=FFE / 1=JENKINS
        elseif kw == "NOHTDREG"; kw_nohtdreg!(s, rec)      # HT-DBH (LHTDRG) calibration control: suppress=no-op, invoke=warn
        elseif kw == "MORTMSB";  kw_mortmsb!(s, rec)       # alternate "mature-stand breakup" mortality (msbmrt.f)
        elseif kw == "PROCESS";  return finish(:process)
        elseif kw in KNOWN_NOOP || kw in variant_noop_keywords(s.variant)
            # recognized no-op — variant-agnostic, or inert for this variant
        else
            # unrecognized keyword — record it so it can't hide as a SILENT gap (YARDLOSS was one such
            # gap in snt01.key/sn.key). Surfaced via s.control.unrecognized_keywords for tests/diagnostics.
            # Only plausible keyword NAMES (alphabetic first char) — not stray numeric data-record tokens.
            (!isempty(kw) && isletter(first(kw))) && push!(s.control.unrecognized_keywords, kw)
        end
    end
end
