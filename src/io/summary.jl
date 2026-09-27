# =============================================================================
# summary.jl — the `.sum` stand-summary table writer (SUMOUT)
#
# Ported from: base/sumout.jl (SUMOUT). Writes the machine-readable summary file:
# a `-999` header line per stand, then one fixed-format row per simulation period.
#
# Column order of the 29-field row (sumout.f / IOSUM, with the before/after-
# treatment split):
#   year age TPA  BA SDI CCF TopHt QMD                         (start of period)
#   Tcuft Mcuft Scuft Bdft                                     (start volumes)
#   remTPA remTcuft remMcuft remScuft remBdft                  (removals)
#   atBA atSDI atCCF atTopHt atQMD                             (after treatment)
#   period accretion mortality MAI  fortype sizecls stockcls   (growth + class)
# =============================================================================

const _SUM_ROW_FMT = Printf.Format(
    "%s%s%s%s%s%s%s%5.1f" *
    "%s%s%s%s%s%s%s%s%s" *
    "%s%s%s%s%5.1f  %s%s%s  %6.1f %s %s%s\n")

# Fortran `Iw` integer edit: right-justified in width `w`, but ALL `w` chars become `*` when the
# value (incl. any sign) doesn't fit — exactly what gfortran prints on integer overflow (e.g. a CCF
# ≥10000 in the sumout.f I4 field ⇒ `****`). Byte-identical to `%wd` for values that fit, so the
# integer columns pre-format through this before the `%s`-based row FORMAT below.
@inline _fi(x::Integer, w::Int) = (s = string(x); length(s) <= w ? lpad(s, w) : "*"^w)

# METRIC variants (BC, Canada/ON) use metric/vbase/sumout.f, whose row FORMAT 20 writes 7I6 volume
# integers (IOSUM 4-10 = total, merch, board, remTrees, remTotal, remMerch, remBoard) — it DROPS the
# sawlog cubic (scuft) + removed-sawlog columns the imperial 9I6 carries. Same leading/trailing fields.
const _SUM_ROW_FMT_METRIC = Printf.Format(
    "%s%s%s%s%s%s%s%5.1f" *
    "%s%s%s%s%s%s%s" *
    "%s%s%s%s%5.1f  %s%s%s  %6.1f %s %s%s\n")

"""
    SummaryRow

One `.sum` period row. Mirrors the SUMOUT/IOSUM columns; the `at_*`
(after-treatment) fields equal the start values when there is no cut in the
period, and the removal fields are 0.
"""
Base.@kwdef mutable struct SummaryRow
    year::Int; age::Int; tpa::Int
    ba::Int; sdi::Int; ccf::Int; topht::Int; qmd::Float64
    cuft::Int; mcuft::Int; scuft::Int; bdft::Int
    rem_tpa::Int = 0; rem_cuft::Int = 0; rem_mcuft::Int = 0; rem_scuft::Int = 0; rem_bdft::Int = 0
    at_ba::Int; at_sdi::Int; at_ccf::Int; at_topht::Int; at_qmd::Float64
    period::Int = 0; accretion::Int = 0; mortality::Int = 0; mai::Float64 = 0.0
    fortype::Int; sizecls::Int; stockcls::Int
end

"Write one SUMOUT-format period row."
function write_sum_row(io::IO, r::SummaryRow; metric::Bool = false)
    if metric
        # metric/vbase/sumout.f FORMAT 20: 7I6 volume block (total, merch, board, remTrees, remTotal,
        # remMerch, remBoard) — DROP the imperial-only sawlog columns scuft + rem_scuft.
        Printf.format(io, _SUM_ROW_FMT_METRIC,
            _fi(r.year,4), _fi(r.age,4), _fi(r.tpa,6), _fi(r.ba,4), _fi(r.sdi,5), _fi(r.ccf,4), _fi(r.topht,4), r.qmd,
            _fi(r.cuft,6), _fi(r.mcuft,6), _fi(r.bdft,6),
            _fi(r.rem_tpa,6), _fi(r.rem_cuft,6), _fi(r.rem_mcuft,6), _fi(r.rem_bdft,6),
            _fi(r.at_ba,4), _fi(r.at_sdi,5), _fi(r.at_ccf,4), _fi(r.at_topht,4), r.at_qmd,
            _fi(r.period,6), _fi(r.accretion,5), _fi(r.mortality,6), r.mai,
            _fi(r.fortype,3), _fi(r.sizecls,1), _fi(r.stockcls,1))
    else
        Printf.format(io, _SUM_ROW_FMT,
            _fi(r.year,4), _fi(r.age,4), _fi(r.tpa,6), _fi(r.ba,4), _fi(r.sdi,5), _fi(r.ccf,4), _fi(r.topht,4), r.qmd,
            _fi(r.cuft,6), _fi(r.mcuft,6), _fi(r.scuft,6), _fi(r.bdft,6),
            _fi(r.rem_tpa,6), _fi(r.rem_cuft,6), _fi(r.rem_mcuft,6), _fi(r.rem_scuft,6), _fi(r.rem_bdft,6),
            _fi(r.at_ba,4), _fi(r.at_sdi,5), _fi(r.at_ccf,4), _fi(r.at_topht,4), r.at_qmd,
            _fi(r.period,6), _fi(r.accretion,5), _fi(r.mortality,6), r.mai,
            _fi(r.fortype,3), _fi(r.sizecls,1), _fi(r.stockcls,1))
    end
    return io
end

# =============================================================================
# CSV form of the .sum — the modern, named-column, machine-readable summary (the
# output analog of the input .tre→.csv / .key→.yaml modernization). Same data as the
# fixed-column .sum (columns documented in docs/FORMATS.md §4), flattened across stands
# with a leading StandID/MgmtID so it loads straight into pandas/R/a spreadsheet.
# =============================================================================
"Named-column header for the summary CSV (StandID/MgmtID/Title + the SUMOUT columns)."
const SUM_CSV_HEADER = [
    "StandID", "MgmtID", "Title", "Year", "Age", "Tpa", "BA", "SDI", "CCF", "TopHt", "QMD",
    "TCuFt", "MCuFt", "SCuFt", "BdFt", "RTpa", "RTCuFt", "RMCuFt", "RSCuFt", "RBdFt",
    "ATBA", "ATSDI", "ATCCF", "ATTopHt", "ATQMD", "PrdLen", "Accret", "Mort", "MAI",
    "ForType", "SizeCls", "StkCls"]

_sumf1(x) = @sprintf("%.1f", Float64(x))   # the float columns (QMD/ATQMD/MAI) print like the .sum
# Minimal CSV quoting: wrap in double-quotes (and double any embedded quote) iff the field
# carries a comma / quote / newline. The Title (the STDIDENT description) can contain commas.
_csvq(s) = (t = String(s); occursin(r"[\",\n]", t) ? '"' * replace(t, '"' => "\"\"") * '"' : t)

"Render one `SummaryRow` as a CSV row (SUM_CSV_HEADER order) for stand `sid`/`mid`/`title`."
function sum_csv_row(sid::AbstractString, mid::AbstractString, title::AbstractString, r::SummaryRow)
    join((_csvq(sid), _csvq(mid), _csvq(title),
          r.year, r.age, r.tpa, r.ba, r.sdi, r.ccf, r.topht, _sumf1(r.qmd),
          r.cuft, r.mcuft, r.scuft, r.bdft, r.rem_tpa, r.rem_cuft, r.rem_mcuft, r.rem_scuft, r.rem_bdft,
          r.at_ba, r.at_sdi, r.at_ccf, r.at_topht, _sumf1(r.at_qmd),
          r.period, r.accretion, r.mortality, _sumf1(r.mai), r.fortype, r.sizecls, r.stockcls), ',')
end

"""
    write_sum_csv(io, stands)

Write the summary CSV: the header then, for each `(stand_id, mgmt_id, title, rows)` in `stands`,
one line per `SummaryRow`. `stands` is a vector of such 4-tuples (one per projected stand). The
`Title` (the STDIDENT description that the fixed-column .sum drops) is carried as its own column.
"""
function write_sum_csv(io::IO, stands)
    println(io, join(SUM_CSV_HEADER, ','))
    for (sid, mid, title, rows) in stands
        for r in rows
            println(io, sum_csv_row(String(sid), String(mid), String(title), r))
        end
    end
    return io
end

"""
    write_sum_header(io, nperiods, stand_id, mgmt_id, sample_wt, variant, date, time, nplots)

Write the per-stand `-999` header line that precedes the period rows.
"""
# Fortran E15.7 edit descriptor: value = 0.DDDDDDD·10^p (0.1≤|m|<1), 7 mantissa digits, 2-digit exponent,
# right-justified in width 15. FVS's .sum header sample-weight uses this — NOT C %15.7E (which prints
# D.DDDDDDD·10^(p-1), 8 sig digits): 11.0 → Fortran "0.1100000E+02" vs C "1.1000000E+01". (io-serialization #2)
function _fortran_e15_7(x::Real)
    x = Float64(x)
    x == 0.0 && return "  0.0000000E+00"
    neg = x < 0; a = abs(x)
    p = floor(Int, log10(a)) + 1
    mi = round(Int, a / 10.0^p * 1.0e7)
    mi >= 10_000_000 && (mi = 1_000_000; p += 1)
    lpad(@sprintf("%s0.%07dE%s%02d", neg ? "-" : "", mi, p < 0 ? "-" : "+", abs(p)), 15)
end

function write_sum_header(io::IO, nperiods::Integer, stand_id::AbstractString,
                          mgmt_id::AbstractString, sample_wt::Real, variant::AbstractString,
                          date::AbstractString, time::AbstractString, nplots::Integer)
    @printf(io, "-999%5d %-26s %-4s%s %-2s %-10s %-8s %-10s %-11s%3d\n",
        nperiods, stand_id, mgmt_id, _fortran_e15_7(Float32(sample_wt)), variant, date, time,
        "          ", "           ", nplots)
    return io
end

"""
    write_sum_file(io, state; period=5, stand_id="", mgmt_id="NONE",
                   sample_wt=nothing, variant="SN", date="", time="", header=true)

Run the full projection and write the complete `.sum` table: the `-999` header
(when `header`) plus one period row per cycle. Each row's start-of-period stats
come from `summary_row` at that cycle; the growth columns (accretion/mortality)
come from `grow_cycle!` advancing to the next period. The final cycle has no
growth period. Cumulative removed merch volume feeds MAI. Requires `state` set up
through `setup_growth!` + `compute_forest_type!` + `compute_volumes!`.
"""
# FVS forms the .sum volume totals by turning each per-tree volume into a per-acre value IN PLACE
# (CFV(I)=CFV(I)*PROB(I), fvs.f:221 at the inventory / gradd.f:303 each cycle), PCTILE-ing those, then dividing
# back (CFV(I)=CFV(I)/PROB(I), fvs.f:269 / gradd.f:350). In REAL*4 (v·p)/p is not always v, so every downstream
# reader of the per-tree arrays before the next VOLS — the FVS_TreeList/CutList rows, the next cycle's CUTS
# removals, ECHARV — sees the round-tripped value (e.g. MCFV 25.0 → 24.999998). gradd.f divides only when
# PROB>0 and round-trips the biomass/carbon arrays only under LFIANVB; the inventory pass (fvs.f) does all of them.
function _vol_prob_roundtrip!(s::StandState, cycle0::Bool)
    t = s.trees
    rt(v::Float32, p::Float32) = (v * p) / p
    bio = cycle0 || s.control.fia_nvb
    @inbounds for i in 1:t.n
        p = t.tpa[i]
        p > 0f0 || continue
        t.cuft_vol[i] = rt(t.cuft_vol[i], p);         t.bdft_vol[i] = rt(t.bdft_vol[i], p)
        t.merch_cuft_vol[i] = rt(t.merch_cuft_vol[i], p); t.saw_cuft_vol[i] = rt(t.saw_cuft_vol[i], p)
        if bio
            t.abvgrd_bio[i] = rt(t.abvgrd_bio[i], p);   t.merch_bio[i] = rt(t.merch_bio[i], p)
            t.cubsaw_bio[i] = rt(t.cubsaw_bio[i], p);   t.foliage_bio[i] = rt(t.foliage_bio[i], p)
            t.abvgrd_carb[i] = rt(t.abvgrd_carb[i], p); t.merch_carb[i] = rt(t.merch_carb[i], p)
            t.cubsaw_carb[i] = rt(t.cubsaw_carb[i], p); t.foliage_carb[i] = rt(t.foliage_carb[i], p)
        end
    end
    return s
end

function write_sum_file(io::IO, s::StandState; period::Int = 5,
                        stand_id::AbstractString = "", mgmt_id::AbstractString = "NONE",
                        sample_wt = nothing, variant::AbstractString = "SN",
                        date::AbstractString = "", time::AbstractString = "", header::Bool = true,
                        collect_rows::Union{Nothing,Vector} = nothing, cycle_hook = nothing,
                        compute_collect::Union{Nothing,Vector} = nothing,
                        cutlist_collect::Union{Nothing,Vector} = nothing,
                        atrtlist_collect::Union{Nothing,Vector} = nothing,
                        carbon_collect::Union{Nothing,Vector} = nothing,
                        potfire_collect::Union{Nothing,Vector} = nothing,
                        hrvcarbon_collect::Union{Nothing,Vector} = nothing,
                        climate_collect::Union{Nothing,Vector} = nothing,
                        canprof_collect::Union{Nothing,Vector} = nothing,
                        strclass_collect::Union{Nothing,Vector} = nothing,
                        dm_collect::Union{Nothing,Vector} = nothing,
                        dm_top4::Vector{Int} = Int[],
                        wwpb_barrier::Union{Nothing,Function} = nothing)
    build_cycle_schedule!(s)                 # ensure the IY boundary-year array is current (idempotent)
    # Accretion/mortality (IOSUM 15/16): disply.f stores INT(OACC(7)/GROSPC+0.5) — grow_cycle! returns the
    # imperial per-area values; the metric variants convert in `metric_sumout` (sumout.f INT(IOSUM·FT3pACRtoM3pHA)).
    _acc_mort(x) = trunc(Int, x + 0.5f0)
    met = _metric_variant(s.variant)
    ncyc = Int(s.control.ncycle_eff)         # rows = ncyc + 1 (inventory + each cycle, post-CYCLEAT)
    ncyc < 1 && (ncyc = Int(s.control.ncycle))
    # FFE Stand Carbon Report (CARBREPT): carbon_on gates only the REPORT-row collection. The per-cycle
    # fuel DYNAMICS (below) run for any FFE-active stand.
    carbon_on = carbon_collect !== nothing && s.fire !== nothing && s.fire.active
    # FFE surface-fuel dynamics (decay + snag falldown + crown-lift → down wood) run EVERY cycle for any
    # FFE-active stand, exactly as FVS fmmain.f does — NOT only when the Carbon report is on. The down-
    # wood pools (FireState.cwd) feed FMCFMD's (SMALL,LARGE) fuel-model selection, so a SIMFIRE stand's
    # fire behavior depends on them having evolved since inventory. The prior carbon-only gate froze cwd
    # at the inventory value (fire_fuel9 2005 sm=7.02 == 1990, vs FVS's accumulated 9.19), which selected
    # fuel model 5 over 12 and gave byram 2905 vs FVS's 4194 (~6× low on the FM5 component).
    # FFE dynamics require the variant's fuel tables to be ported; a variant still building its growth port
    # (e.g. NC/Klamath) has no ffe_fuel_live yet ⇒ FFE stays inert (no snag seeding / fuel dynamics) rather
    # than erroring. Ported variants all have fuel tables ⇒ unaffected.
    # `ffe_fuel_live` is the EASTERN (SN/NE/CS/LS) live-fuel table; the WESTERN variants (NC and the CR
    # family) carry their live+dead fuel in their own cover-type/top-2 loaders, so their ffe_fuel_live is
    # empty — but their FFE fuel loop MUST still run. Without NC here, ffe_on was false ⇒ the per-cycle
    # ffe_fuel_update! + fire_smlg were SKIPPED ⇒ the fire sampled an unaccumulated ~empty down-wood pool
    # (fire_smlg=(0.37,0) vs the real cwd ~10) ⇒ wrong fuel-model weights / under-fire. (The CR-family
    # variants likely share this latent run_keyfile gap — validate + fold them in separately.)
    # EastCascades (like NC/the CR family) carries its live+dead fuel in its own single-cover-type loaders
    # (ec_live/dead_fuel_loading), so its ffe_fuel_live is empty — but its FFE fuel loop MUST still run, else
    # the per-cycle ffe_fuel_update! + fire_smlg stash are skipped ⇒ the SIMFIRE samples a (0,0) down-wood
    # point ⇒ FMDYN drops the natural-fuel model (ect01_ffe FMD-set {6,9} vs live {9,10,6}) ⇒ under-fire.
    # (Unlike WC/PN/CA/WS whose crown fire dominates and masks the fuel-model weights, EC's is a SURFACE fire.)
    # IE/KT (InlandEmpire/Kootenai) carry their live+dead fuel in their own COVTYP loaders (ie_live/dead_
    # fuel_loading), so ffe_fuel_live is empty — but the FFE fuel loop MUST still run, else the per-cycle
    # ffe_fuel_update! + fire_smlg stash are skipped ⇒ the SIMFIRE samples an empty (0,0) down-wood point ⇒
    # FMCFMD falls back to fuel model 8 (byram ~370 vs the accumulated-fuel ~637 on crown-prone dense stands)
    # ⇒ under-fire. Folded in with the decay-table fix (ie/fmvinit.f DKR is decay-class-INDEPENDENT, _FM_DKR_NR)
    # so the accumulated fine fuel matches the oracle rather than over-retaining ~2.2×. (The rest of the CR
    # family — CR/EM/CI/TT/UT/WC/PN/CA/WS — shares this latent gap; fold them in + validate separately.)
    # BM (bm/fmcba.f bm_live/dead_fuel_loading + bm/fmvinit.f DKR, both ported) folded in: without it the BM
    # SIMFIRE sampled a (0,0) down-wood point ⇒ FMDYN picked model 8 alone (live 10+12 at SMALL 1.58/LARGE
    # 29.87) ⇒ flame 0.8 ft SURFACE vs live 11.5 ft PASSIVE ⇒ catastrophic under-kill (41134029010497:
    # 1078→222 TPA vs live 1078→1).
    ffe_on = s.fire !== nothing && s.fire.active &&
             (!isempty(s.coef.ffe_fuel_live) || s.variant isa Klamath || s.variant isa EastCascades ||
              s.variant isa SouthCentralOregon || s.variant isa OregonCoast || s.variant isa Olympic ||
              s.variant isa InlandEmpire || s.variant isa Kootenai ||
              s.variant isa BlueMountains || _ffe_west_vol(s.variant) ||
              s.variant isa WestCascades || s.variant isa PacificNorthwest)   # OC/OP live fuel in fire_fuel_covtype_live.csv (non-reserved) ⇒ ffe_fuel_live empty
    # WC/PN (wc/pn fmmain.f run FMSDIT + the annual FMSNAG/FMCWD/FMCADD loop like every FFE variant) were off this
    # list ⇒ no inventory snags, no fuel dynamics: pnt01 stand 4's 2003 SIMFIRE sampled SMALL/LARGE 0.46/0.00 vs live
    # 3.79/11.65 ⇒ fuel models 2/5 instead of live 5/10 (FMDYN 0.66/0.34), flame 4.33 vs 6.09 ft, 2013 TPA 127 vs 71.
    # CI/TT/UT/CR (ci/tt/ut/cr fmsdit.f + FMSNAG/FMCWD/FMCADD, the same annual loop as IE/EM) were off this list ⇒ no
    # inventory snags and no fuel dynamics: S248112 no-fire DDW 2.40 flat vs live 4.65→4.34 (CI), Standing_Dead 0.
    if ffe_on
        ffe_seed_input_snags!(s)             # inventory snags from the input dead records (FMSADD ITYP=3)
        fill!(s.fire.crown_lift_annual, 0f0)
        snapshot_ffe_oldcrown!(s)            # FMOLDC at inventory: gives the 1st cycle's crown-lift a valid
                                             # OLDCRW (else the 1st cycle's fine down-wood is lost; DDW gap)
    end
    # SNAGINIT (act 2522, FMSADD) is a SCHEDULED activity that FVS runs in FMMAIN DURING cycle 1, AFTER the
    # inventory carbon report — NOT at inventory like the input dead-tree snags. So the user snags first appear
    # in cycle-1's report (live net01 SnagDet is empty at the inventory year; the SNAGINIT cohort shows up the
    # next cycle). Defer the add to the start of the first growing cycle.
    snaginit_pending = ffe_on
    g = s.plot.gross_space
    # .sum header sample weight = SAMWT (s.plot.sample_weight), NOT gross_space. (gross_space
    # is the non-stockable expansion used for per-acre normalization below; SAMWT is the
    # stand's sampling weight, which FVS prints in the -999 header — e.g. 11, not 1.1.)
    sw = sample_wt === nothing ? s.plot.sample_weight : sample_wt
    # FVS SAMWT default is 1.0 (a blank DESIGN sample-weight field ⇒ 1.0). When the points_inv
    # fallback can't fire (points_inv 0), sample_weight stays 0; emit the 1.0 default so the
    # header matches FVS.
    sw <= 0f0 && (sw = 1f0)
    if header
        write_sum_header(io, ncyc + 1, stand_id, mgmt_id, sw, variant, date, time, Int(s.plot.pi))
    end
    # MAI (evtstv.f:348-440 per cycle, disply.f:391 final row) — FVS's state machine, not a plain (merch+removed)/age:
    # MAIFLG=1 (MAI off for good) when the inventory age is 0 on a stocked stand or a later row finds the age reset to
    # 0 (ZERO=AGE−(IY(ICYC)−IY(ICYC−1))=0) without a new stand; NEWSTD=1 after a clearcut or on bare ground. TOTREM
    # adds the previous cycle's INTEGER merch removal IOSUM(9) only when the age advanced and the row's MAI is on,
    # restarts at an age decrease or a clearcut; the final row reads TOTREM without the last cycle's removal.
    maiflg = false; newstd = false; agelst = 0; totrem = 0f0
    cc_prev = false        # the last growing cycle's after-treatment values show a clearcut (evtstv.f:426-427)
    prev_increment = 0f0   # IOSUM(9) of the most recent growing cycle
    cover_year0 = 0        # COVER: inventory year (IY(1)) for ICVAGE offset
    di(x) = trunc(Int, x + 0.5)
    prev_rem_scuft = 0                          # last growing cycle's sawlog-cubic removal (IOSUM(22) carry)
    for c in 0:ncyc
        # The inventory row / cycle-0 TreeList read the density CRATET's final DENSE left (cratet.f:578, over CRATET's
        # IND — IND1-seeded RDPSRT(.FALSE.), or the :257 identity RDPSRT(.TRUE.) with dead records), the same state the
        # first grow_cycle! rebuilds before DGDRIV. The setup pass used a fresh sort, so equal-DBH trees swapped their
        # PCT/PTBAL (MEASURED FVSem_g16 196420598020004 2012: 33 BAPctile/PtBAL cells permuted among tied records).
        c == 0 && compute_density!(s; cratet_ind = _fvs_ind_lifecycle(s.variant))
        compute_forest_type!(s)
        last = c == ncyc
        per = last ? 0 : cycle_period_at(s.control, c)   # THIS cycle's length (varies w/ TIMEINT/CYCLEAT)
        # main columns reflect the start-of-cycle (pre-thin) stand.
        # MAI terminal-row quirk (evtstv.f:414 + disply.f:392): intermediate rows accumulate
        # removed merch with a one-cycle lag, but the FINAL row's MAI is loaded from the
        # un-incremented TOTREM — i.e. it EXCLUDES the last growing cycle's removal.
        age_c = _summary_age(s); mai_zero = false; zero_c = -1
        if c == 0                                                   # evtstv.f:356-377 (ICYC=1)
            if age_c <= 0
                mai_zero = true
                stand_tpa(s) == 0f0 ? (newstd = true) : (maiflg = true)
            end
        elseif !last                                                # evtstv.f:384-400
            zero_c = age_c - cycle_period_at(s.control, c - 1)
            age_c < agelst && (totrem = 0f0)
            if zero_c == 0 || (maiflg && !newstd)
                mai_zero = true
                newstd || (maiflg = true)
            elseif age_c > agelst
                totrem += prev_increment
            end
        else                                                        # disply.f:391 (IOSUM(2)=AGE)
            mai_zero = maiflg || age_c <= 0
        end
        r = summary_row(s; period = per, total_removed_merch = totrem, cycle0 = c == 0,
                        final_row = last, mai_zero = mai_zero)
        if c > 0 && !last                                           # evtstv.f:426-434
            cc_prev && (newstd = true; totrem = 0f0)
            zero_c == 0 && stand_tpa(s) == 0f0 && (newstd = true)
        end
        last || (agelst = age_c)
        _vol_prob_roundtrip!(s, c == 0)   # fvs.f:221/269 (cycle 0) / gradd.f:303/350: per-tree V·PROB ... /PROB
        # per-cycle hook (DBS TreeList): the start-of-cycle (pre-thin) tree list at year r.year.
        # `c` is the cycle index (0 = inventory) — dbstrls.f emits input dead records only at cycle 0.
        # dbstrls.f binds PrdLen = IFINT: at the end of FVS cycle ICYC (year IY(ICYC+1) = jl cycle c's start) that is
        # the FINT of the cycle just grown (jl c−1) — so the terminal list reports the last cycle's length, not the
        # .sum final row's 0; the inventory list (c=0, before any growth) carries the first cycle's FINT.
        cycle_hook === nothing || cycle_hook(s, r.year, c == 0 ? per : cycle_period_at(s.control, c - 1), c)
        # Test-only observer (tiered suite bit-identity snapshots, test/harness/tiered/snapshot.jl): a callback in the
        # CURRENT TASK's local storage sees the same start-of-cycle state. Task-local ⇒ safe when stands run on
        # parallel tasks; absent ⇒ one Dict lookup per summary row, no effect on the simulation.
        let snap = get(task_local_storage(), :fvsjl_snapshot_hook, nothing)
            snap === nothing || snap(s, r.year, c)
        end
        # COVER report-only accumulator (CVCNOP): the canopy statistics of the start-of-cycle
        # (pre-thin) stand at year r.year → slot IP1=c+1. Gated on the COVER activity 900.
        if s.cover !== nothing && s.cover.active
            compute_density!(s)
            c == 0 && (cover_year0 = Int(r.year))
            # FINT as CVCBMS sees it (DDS=(2·D·DG+DG²)/FINT): at the inventory row (fvs.f:301, before any
            # GRINCR) FINT is still the GROWTH/DG_MEASURE measurement period (initre.f:831, dbsstandin.f:701;
            # default grinit.f) — not the first cycle length. Later rows: the grown cycle's IY(ICYC+1)-IY(ICYC)
            # (grincr.f:65-66).
            cover_fint = c == 0 ? s.control.growth_fint : cycle_period_at(s.control, c - 1)
            # Inventory row: CVCBMS reads FVS's post-calibration DG array (dgdriv.f DO 220 dubs every record
            # without a measured increment), not the raw input increments jl keeps in diam_growth.
            if c == 0 && s.variant isa BlueMountains
                dg_in = copy(view(s.trees.diam_growth, 1:s.trees.n))
                s.trees.diam_growth[1:s.trees.n] .= bm_cycle0_dg(s)
                cover_accumulate!(s.cover, s, r.year, cover_year0, cover_fint)
                s.trees.diam_growth[1:s.trees.n] .= dg_in
            else
                cover_accumulate!(s.cover, s, r.year, cover_year0, cover_fint)
            end
        end
        # FFE Stand Carbon Report row (FMCRBOUT, fmmain.f:206) — sampled at the FVS phase: AFTER FMBURN
        # (fire kill + snag booking + consumption) but BEFORE UPDATE grows the stand. For a non-fire cycle
        # that phase equals the cycle-top, pre-growth stand (sampled here). For a SIMFIRE cycle the row
        # must be POST-fire, so defer it to grow_cycle!'s `carbon_hook` (fired right after the fire, on the
        # post-fire cycle-start-size stand) — else the fire's snag/AGL/Released effects surface one row late.
        # Carbon-Released-from-Fire: 0 unless a SIMFIRE burned in r.year (fmburn! records it in burn_reports);
        # convert tons-C/ac → the report units (same factor as stand_carbon_report's pools, carbon.jl).
        _carb_push(st; vtrip::Bool = false) = begin
            rel = 0f0; tcon = 0f0
            if st.fire !== nothing
                @inbounds for br in st.fire.burn_reports
                    br.year == Int(r.year) && (rel = br.released; tcon = get(br, :totcon, 0f0)::Float32)
                end
            end
            uf = st.control.carbon_units == 1 ? 0.90718f0 / 0.4046945f0 :     # METRIC.F77 TItoTM / ACRtoHA
                 st.control.carbon_units == 2 ? 0.90718f0 : 1f0
            # FVS_Fuels Consumed = NINT(TOTCON) of the fire burned in this FMDOUT year (fmdout.f:269/403)
            fl = merge(ffe_fuel_loadings(st), (consumed = tcon,))
            push!(carbon_collect, (r.year, stand_carbon_report(st; vtrip = vtrip), fl,
                                   snag_summary(st), ffe_down_wood(st), rel * uf, snag_detail(st)))
        end
        # A SIMFIRE cycle: the fire (inside grow_cycle!'s mortality_and_fire!) must consume + snag the
        # START-of-cycle fuels, so this cycle's pre-grow ffe_fuel_update! is WITHHELD and its period handed
        # to grow_cycle! (run post-fire, FVS FMBURN→annual-loop order). The carbon row is likewise deferred
        # to the post-fire `carbon_hook`. (`fire_cycle` adds the carbon-report gate.)
        fire_this_cycle = !last && _fire_due(s) && per > 1   # OPCYCL: cycle range contains fire_year
        fire_cycle = carbon_on && fire_this_cycle
        # FFE reports (FMMAIN, gradd.f:118): FMCBA → FMSSUM/FMPOCR/FMPOFL/FMDOUT/FMCRBOUT run inside GRADD, AFTER this
        # cycle's CUTS (grincr.f:292 → cuts.f:1823 FMSCUT slash + YARDLOSS snags), so a thinned cycle reports the POST-
        # cut stand and the one-time dead-fuel load (fmcba.f:457) reads the post-cut FMTBA/PERCOV. Growing cycles run
        # this right after cuts! below; the final (post-projection) row keeps the cycle-top call (no cut happens).
        carb_v3_pending = nothing   # (carbon_collect index, vtrip) of this cycle's pre-growth FMCRBOUT row
        _ffe_reports! = function ()
        # FMMAIN (and its FMCRBOUT/FMDOUT/FMSSUM reports) runs once per projection CYCLE (fvs.f cycle loop), never
        # for the post-projection final row — live FVS_Carbon/Fuels/SnagSum carry NUMCYCLE rows, no final year.
        if carbon_on && !fire_cycle && !last
            compute_density!(s)
            fmcba!(s)
            _vt = _fm_will_triple(s)
            _carb_push(s; vtrip = _vt)   # FMMAIN runs on the tripled list in a tripling cycle
            carb_v3_pending = (length(carbon_collect), _vt)   # V(3) re-derived at the FMMAIN seam (grow_cycle!)
        end
        # FVS_PotFire: the potential-fire behavior under fixed severe/moderate weather (FMPOFL), per cycle
        if potfire_collect !== nothing && s.fire !== nothing && s.fire.active && !isempty(s.coef.ffe_fuel_live)
            compute_density!(s)
            pfr = potential_fire_report(s)
            pfr !== nothing && push!(potfire_collect, (r.year, pfr))
        end
        # FVS_CanProfile (fmpocr.f mode 2, fmmain.f:188): the PRE-growth (cycle-start inventory) canopy crown-fuel
        # profile — reported alongside FMPOFL, BEFORE this cycle's growth (distinct from the post-growth carbon path).
        if canprof_collect !== nothing && s.fire !== nothing && s.fire.active
            compute_density!(s)
            push!(canprof_collect, (Int(r.year), canopy_crfill(s)))
        end
        end
        last && _ffe_reports!()
        # FVS_StrClass (sstage.f → dbsstrclass.f): the SSTAGE structure classification, BEFORE-thin (Removal_Code
        # 0) at the cycle-top stand. The AFTER-thin (cd=1) row is captured post-cuts! below (non-last cycles).
        if strclass_collect !== nothing
            compute_density!(s)
            push!(strclass_collect, (Int(r.year), 0, structure_report(s)))
        end
        # FVS_DM_* dwarf-mistletoe summary (misprt.f MISPRT): the START-OF-CYCLE stand state; the DM
        # mortality columns project the cycle's mortality (MISMRT) over this cycle's length. NAGE advances
        # with the report year (IAGE + year − inventory year). Collected only on DM-infected stands.
        if dm_collect !== nothing && _dm_report_active(s)
            compute_density!(s)
            nage = Int(s.plot.stand_age) + (Int(r.year) - Int(s.control.cycle_year[1]))
            # MISMRT projects the cycle's mortality; the final report row (per==0) reuses the
            # previous cycle's length (matching FVS's non-zero final-row DM mortality).
            perdm = last ? cycle_period_at(s.control, max(c - 1, 0)) : per
            push!(dm_collect, (Int(r.year),
                  mistletoe_report(s; fint = Float32(perdm), top4 = dm_top4, nage = nage)))
        end
        if !last
            # FMSDIT (grincr.f:227, before CUTS): FMCROW's height percentiles for this cycle's CROWNW.
            ffe_on && ffe_snapshot_hpct!(s)
            # DBS FVS_Compute: snapshot the active COMPUTE variables at this (growing) cycle's
            # start — only the growing cycles get a row (the event monitor runs during growth).
            compute_collect === nothing ||
                push!(compute_collect, (r.year, snapshot_compute!(s, r.year, c)))
            # apply this cycle's scheduled thin (CUTS) BEFORE growth; report the removed
            # + after-treatment columns on THIS row (matching the Fortran .sum). cuts! is
            # idempotent, so grow_cycle!'s own cuts! call below is then a no-op.
            # FVS_CutList: arm the per-record cut sink for this (real) thin, then stash + disarm.
            cutlist_collect === nothing || (s.control.cutlist_capture = Any[])
            atrtlist_collect === nothing || (s.control.atrtlist_capture = Any[])
            econ_cycle_start!(s)   # ECON ECSETP/ECSTATUS(…,0) precede CUTS (grincr.f:273) — ECHARV needs the start year
            rem = cuts!(s; fint = Float32(per))
            if cutlist_collect !== nothing
                push!(cutlist_collect, (r.year, per, s.control.cutlist_capture))
                s.control.cutlist_capture = nothing
            end
            if atrtlist_collect !== nothing
                push!(atrtlist_collect, (r.year, per, s.control.atrtlist_capture))
                s.control.atrtlist_capture = nothing
            end
            # FVS_StrClass AFTER-thin row (Removal_Code 1), post-cuts! (identical to the cd=0 row on a no-thin cycle).
            if strclass_collect !== nothing
                compute_density!(s)
                push!(strclass_collect, (Int(r.year), 1, structure_report(s)))
            end
            if rem.tpa > 0f0
                compute_density!(s)
                r.rem_tpa  = di(rem.tpa / g);  r.rem_cuft  = di(rem.cuft / g)
                r.rem_mcuft = di(rem.mcuft / g); r.rem_scuft = di(rem.scuft / g)
                # metric vols.f never loads MCFV (see summary_row), so cuts.f's CMCUT = Σ PREM·MCFV stays 0 and
                # OMCREM(7)/OSCREM(7) are never accumulated (cuts.f:1779-1783): removed merch/sawlog cubic are 0.
                met && (r.rem_mcuft = 0; r.rem_scuft = 0)
                r.rem_bdft = di(rem.bdft / g)
                r.at_ba = di(stand_ba(s) / g);  r.at_sdi = di(stand_sdi(s) / g)
                r.at_ccf = di(stand_ccf(s) / g); r.at_topht = di(stand_top_height(s))
                r.at_qmd = Float64(stand_qmd(s))   # QDBHAT=ATAVD (REAL): FVS_Summary binds it whole; the .sum prints F5.1
                # TOTREM: evtstv.f:414 (CASE DEFAULT — every variant; the eastern CASE is commented out) accumulates
                # the INTEGER IOSUM(9,ICYC-1)=INT(OMCREM(7)/GROSPC+0.5) (disply.f:341), and BCYMAI=(TOTREM+CURVOL)/AGE.
                # The per-acre float removal agreed at the .sum's F5.1 but not in FVS_Summary's REAL MAI (MEASURED
                # FVSem_g16 196378260020004 thinbba 2032: live 10.668750 = 1707/160, jl 10.670339 = 1707.25/160).
                prev_increment = Float32(r.rem_mcuft)
                # evtstv.f:506-529 TSTV2(5)=BA, (6)=RELDEN, (14)=SDI after treatment all 0 with trees removed
                cc_prev = stand_ba(s) == 0f0 && stand_ccf(s) == 0f0 && stand_sdi(s) == 0f0
            else
                prev_increment = 0f0   # this growing cycle had no removal
                cc_prev = false
            end
            _ffe_reports!()          # FMMAIN reports on the post-cut stand (see _ffe_reports! above)
            # SNAGINIT is processed at the top of FMSNAG (fmsnag.f:88-100), i.e. inside the annual loop that follows the
            # reports, so its snags first surface in the NEXT cycle's report.
            if snaginit_pending
                ffe_add_snaginit!(s); snaginit_pending = false
            end
            # FVS_Hrv_Carbon: collect AFTER this cycle's cut so year r.year's harvest is booked (KYR=1),
            # matching FMCHRVOUT. Habitat default 1 (SC); SE (2) applies to specific national forest codes.
            hrvcarbon_collect === nothing || s.fire === nothing || !s.fire.active ||
                push!(hrvcarbon_collect, (r.year, harvested_carbon_report(s, r.year, 1)))
            # FFE annual fuel loop BEFORE growth (report→fuel→grow). FVS interleaves the annual fuel steps
            # with FMBURN, so a fire fires on the cycle-start + fire-year's single annual step — NOT the
            # period-end fuel. When a SIMFIRE burns this cycle, split the loop: advance 1 year, stash the
            # (SMALL,LARGE) the fire burns on, then advance the rest. Non-fire cycles run the full loop once.
            ffe_defer_init = false
            r6_defer_fuel = false
            if ffe_on
                # FMMAIN runs FMBURN (the fire, fmmain.f:170) BEFORE the annual fuel loop (FMSNAG/FMCWD/
                # FMCADD, fmmain.f:228), so the fire samples the START-OF-CYCLE down wood. Stash (SMALL,LARGE)
                # for the fire basis here. For a fire cycle the annual loop is DEFERRED into grow_cycle!
                # (run post-fire); non-fire cycles advance the pools the full period now (report→fuel→grow).
                # A fire in the FIRST FFE cycle (s10_fire SIMFIRE@cycle1) burns before any prior ffe_fuel_update!
                # loaded the dead-fuel pools, so init them now (FVS FMCBA precedes the first FMBURN). Otherwise
                # the stash would read zero cwd ⇒ low-fuel model ⇒ under-kill (jl 119 vs live 57 TPA). For
                # cycle≥2 fires the pools are already loaded (fuels_init), so fire_carbon stays bit-exact.
                fire_this_cycle && !s.fire.fuels_init && (compute_density!(s); fmcba!(s))
                fire_this_cycle && (s.fire.fire_smlg = _small_large_fuel(s.fire))
                # Defer the VERY FIRST (init-year) dead-fuel load into grow_cycle! post-cuts! (FVS loads FMCBA
                # after the cut phase). Later cycles (already init) run the annual loop pre-grow as before.
                # CR-ONLY: the eastern (SN/NE/CS/LS) FFE + its PotFIRE/carbon reports are separately validated and
                # their init-year has no thin (pre==post), so keep their pre-grow load path untouched (doctrine #5).
                if fire_this_cycle
                elseif r6_ffe_code(s.variant) !== :none && s.fire.fuels_init
                    # R6 variants (BM/EC/SO/PN/OP/WC): FMSNAG's FMR6HTLS height-loss draws start from the main RNG
                    # state at FMMAIN (gradd.f:118, after this cycle's DGDRIV/REGENT/MORTS/TRIPLE/MISTOE draws) and are
                    # rolled back each year (fmsnag.f:113-116/290-293). Pre-grow, jl held the PREVIOUS cycle's FMMAIN
                    # state (live BM 22960873010497 cycle-3 draws showed up in jl's cycle 4). Run the annual loop at the
                    # FMMAIN point inside grow_cycle! (mortality_and_fire!'s post_fire seam), like a fire cycle.
                    r6_defer_fuel = true
                elseif s.fire.fuels_init || !(s.variant isa CentralRockies)
                    ffe_fuel_update!(s, per)
                else
                    ffe_defer_init = true
                end
            end
            # FMMAIN (fmmain.f:139-206): FMCBA runs once at the year start; between FMBURN and FMCRBOUT only SN re-runs it
            # (sn/fmburn.f:589, the post-burn FULIV2 shrub age), so EM's fire-year FLIVE is the PRE-fire load (live
            # 196378260020004 2022 Shrub_Herb 0.263 = ½·(0.285+0.242); a post-fire FMCBA gave 0.392).
            # Only SN re-runs it (fmburn.f:588 `IF (BURNYR.EQ.IYR .AND. VARACD.EQ.'SN')`); no western variant does (EC ect01 2003
            # Surface_Shrub post-fire FMCBA 1.649 vs live pre-fire 0.420). CS/LS/NE keep the re-run pending a live check.
            _refmcba = st -> (st.variant isa Southern || st.variant isa CentralStates || st.variant isa LakeStates ||
                              st.variant isa Northeast)
            chook = fire_cycle ? (st -> (compute_density!(st); _refmcba(st) && fmcba!(st); _carb_push(st))) : nothing
            _v3p = carb_v3_pending; carb_v3_pending = nothing
            fhook = _v3p === nothing ? nothing :
                    ((st, stash) -> (e = carbon_collect[_v3p[1]];
                                     carbon_collect[_v3p[1]] = Base.setindex(e, carbon_report_fmmain_v3(e[2], st, stash, _v3p[2]), 2)))
            gr = grow_cycle!(s; fint = Float32(per), carbon_hook = chook, fmmain_hook = fhook,
                             fuel_period = (fire_this_cycle || r6_defer_fuel) ? per : nothing,
                             ffe_init_period = ffe_defer_init ? per : nothing,
                             wwpb_barrier = wwpb_barrier)   # advances cycle (PPE mode-2 LIVE seam)
            r.accretion = _acc_mort(gr.accretion)
            r.mortality = _acc_mort(gr.mortality)
            # Climate-FVS report (clauestb.f/DBSCLSUM): collected POST-growth, labeled with the cycle-START year,
            # viability sampled at report_year+fint/2 — the offset FVS uses (see climate_report). c.spmort1/2 were
            # just populated by grow_cycle!'s apply_climate_mort!.
            climate_collect === nothing || (s.climate !== nothing && s.climate.active) &&
                push!(climate_collect, (Int(r.year), something(s.climate.pending_report,
                                                              climate_report(s; report_year = Int(r.year), fint = per))))
            if ffe_on                                   # crown-lift from THIS growth (FMSDIT) + FMOLDC snapshot
                compute_crown_lift!(s, per); snapshot_ffe_oldcrown!(s)   # OLDCRW from last FMSDIT's CROWNW
                ffe_snapshot_hpct!(s)                   # then FMCROW (fmsdit.f:128) re-ranks the grown stand
            end
        end
        # (FMCHRVOUT runs only inside FMMAIN's cycle loop, so no FVS_Hrv_Carbon row for the terminal .sum year —
        #  live Hrv_Carbon has exactly the FVS_Carbon years in all 96 tiered fixtures.)
        # disply.f:382-387 zeroes the FINAL row's removal columns IOSUM(7..10) (and 14..16) but NOT IOSUM(22), the
        # later-added sawlog-cubic removal (disply.f:342 INT(OSCREM(7)/GROSPC+.5)); CUTS zeroes OSCREM only at its
        # own entry (cuts.f:323-329) and fvs.f:432 resets only ONTREM(7) ⇒ the final row carries the LAST growing
        # cycle's sawlog removal (live econ_strtecon 2005: SCuFt removed 23 = the 2000 thin's).
        last ? (r.rem_scuft = prev_rem_scuft) : (prev_rem_scuft = r.rem_scuft)
        # metric/vbase/sumout.f: BC/ON print (and hand DBSSUMRY) the per-ha metric conversion of the IOSUM row.
        rout = met ? metric_sumout(r) : r
        write_sum_row(io, rout; metric = met)
        collect_rows === nothing || push!(collect_rows, rout)
    end
    return io
end

# The two metric variants (canada BC / ON): FVS compiles metric/vbase/{disply,sumout}.f and the metric dbsqlite writers.
_metric_variant(v) = v isa BritishColumbia || v isa Ontario

"""
    metric_sumout(r) -> SummaryRow

metric/vbase/sumout.f:328-352 (identical in the BC and ON builds): the per-ha metric row FVS prints to the .sum
and passes to DBSSUMRY, formed from disply.f's IMPERIAL per-acre integers in `r` with REAL arithmetic truncated by
INT(): trees/SDI `INT(IOSUM/ACRtoHA)`, BA `INT(IOLDBA·FT2pACRtoM2pHA)`, heights `INT(IBTAVH·FTtoM)`, volumes /
accretion / mortality `INT(IOSUM·FT3pACRtoM3pHA)`; CCF, year, age, period and the classes pass through. QMD/ATQMD
= QSDBT·INtoCM / QDBHAT·INtoCM and MAI = BCYMAI·FT3pACRtoM3pHA stay unrounded REALs (F5.1/F6.1 on print).
"""
function metric_sumout(r::SummaryRow)
    ha(x)  = trunc(Int, Float32(x) / ACRE_TO_HA)
    ba(x)  = trunc(Int, Float32(x) * FT2PACRE_TO_M2PHA)
    ht(x)  = trunc(Int, Float32(x) * FT_TO_M)
    vol(x) = trunc(Int, Float32(x) * FT3PACRE_TO_M3PHA)
    cm(x)  = Float64(Float32(x) * IN_TO_CM)
    SummaryRow(year = r.year, age = r.age, tpa = ha(r.tpa), ba = ba(r.ba), sdi = ha(r.sdi), ccf = r.ccf,
        topht = ht(r.topht), qmd = cm(r.qmd),
        cuft = vol(r.cuft), mcuft = vol(r.mcuft), scuft = vol(r.scuft), bdft = vol(r.bdft),
        rem_tpa = ha(r.rem_tpa), rem_cuft = vol(r.rem_cuft), rem_mcuft = vol(r.rem_mcuft),
        rem_scuft = vol(r.rem_scuft), rem_bdft = vol(r.rem_bdft),
        at_ba = ba(r.at_ba), at_sdi = ha(r.at_sdi), at_ccf = r.at_ccf, at_topht = ht(r.at_topht),
        at_qmd = cm(r.at_qmd), period = r.period, accretion = vol(r.accretion), mortality = vol(r.mortality),
        mai = Float64(Float32(r.mai) * FT3PACRE_TO_M3PHA),
        fortype = r.fortype, sizecls = r.sizecls, stockcls = r.stockcls)
end

"""
    summary_row(state) -> SummaryRow

Build the start-of-period `.sum` row from the current (already-grown-to-cycle)
stand state: per-acre TPA/BA/SDI/CCF/top-height/QMD, the four stand volumes, and
the forest-type / size / stocking classes. Removal, after-treatment and growth
(accretion/mortality/MAI) fields are filled by the cycle driver. The integer
columns use FVS's truncate-after-+0.5 rounding (`_dtrunc`)."""
# TSTV1(2)/IOSUM(2): the stand age at the current cycle's start (evtstv.f:260 IAGE+IY(ICYC)−IY(1)), rebased after a
# RESETAGE (resage.f runs after DISPLY, so the reset row itself keeps the old age).
function _summary_age(s::StandState)::Int
    yr = cycle_year_at(s.control, Int(s.control.cycle))
    age = Int(s.plot.stand_age) + (yr - Int(s.control.cycle_year[1]))
    ry = Int(s.control.age_reset_year)
    (ry >= 0 && yr > ry) && (age = Int(s.control.age_reset_age) + (yr - ry))
    return age
end

function summary_row(s::StandState; period::Int = 0, total_removed_merch::Real = 0,
                     accretion::Real = 0, mortality::Real = 0, cycle0::Bool = false,
                     final_row::Bool = false, mai_zero::Bool = false)
    g = s.plot.gross_space
    dt(x) = trunc(Int, x + 0.5f0)
    # Metric variants (BC, Ontario): this row holds disply.f's IMPERIAL per-acre IOSUM stage (the same integers
    # every variant stores); metric/vbase/sumout.f converts them to per-ha metric only when printing / calling
    # DBSSUMRY — see `metric_sumout`. QMD and MAI are kept UNROUNDED here (QSDBT / BCYMAI are REALs).
    met  = _metric_variant(s.variant)
    # BM: FVS's .sum TPA and volume totals are PCTILE totals (gradd.f:289-322 / cratet.f:682) — a Float32
    # cumulative sum walking IND BACKWARDS (smallest DBH first, pctile.f), over PROB and over CFV·PROB etc.
    # formed in Float32 — not a record-order sum. IND = CRATET's order on the cycle-0 row (bm_cratet_ind!),
    # gradd.f:186's fresh RDPSRT(.TRUE.) after. Record order flipped knife-edge rows by ±1 (41134819010497:
    # per-record TPA bit-identical, record-order Σ 1331.49988 → 1331 vs live 1332). BM-gated (the base
    # gradd.f is shared; other variants not yet re-validated on this order).
    bm_ind = nothing
    if s.variant isa BlueMountains && s.trees.n > 0
        bm_ind = Vector{Int32}(undef, s.trees.n)
        cycle0 ? bm_cratet_ind!(s, bm_ind) : _rdpsrt!(view(s.trees.dbh, 1:s.trees.n), bm_ind)
    end
    pctile_tot(w) = (acc = 0f0; @inbounds(for k in length(bm_ind):-1:1; acc += w(Int(bm_ind[k])); end); acc)
    tpa  = bm_ind === nothing ? dt(stand_tpa(s) / g) :
           dt(pctile_tot(i -> s.trees.tpa[i]) / g)
    ba   = dt(stand_ba(s) / g)
    sdi  = dt(stand_sdi(s) / g)
    ccf  = dt(stand_ccf(s) / g)
    # BM cycle-0 row: FVS's AVH (DENSE at cratet.f:692 / AVHT40 :624) walks the IND CRATET left — the IND1-seeded
    # RDPSRT(.FALSE.) of cratet.f:166 when no dead were deleted (:197 skips :270), else :270's fresh sort
    # (bm_cratet_ind!). A fresh sort here broke 40-TPA-cutoff DBH ties (23900114010900 PP/GF 8.3": 45 vs live 46).
    toph = dt(stand_top_height(s; cratet_ind = cycle0 && _fvs_ind_lifecycle(s.variant)))
    qmd  = Float64(stand_qmd(s))   # QSDBT=ORMSQD (REAL): sumout.f:209 binds it whole to FVS_Summary; the .sum prints F5.1
    t = s.trees
    # STRICTLY SEQUENTIAL Float32 accumulation (ACC += VOL[i]·PROB[i], i=1..n) to match FVS's DISPLY DO-loop
    # order — Julia's `sum(generator)` may use PAIRWISE reduction, which reorders the Float32 adds and flips
    # the rendered integer by 1 on knife-edge rows (the non-associative tree-SUM residual).
    function vtot(f)
        # `f` is a runtime Symbol ⇒ getfield(t, f) infers as Any, boxing `fld[i]` on every add (measured
        # ~50 KB/cycle + a type-instability). All vtot fields are Vector{Float32}, so assert it: concrete
        # `fld` ⇒ allocation-free, type-stable, and the sequential Float32 accumulation order is unchanged.
        fld = getfield(t, f)::Vector{Float32}; acc = 0f0
        if bm_ind !== nothing
            acc = pctile_tot(i -> fld[i] * t.tpa[i])        # CFV(I)=CFV(I)*PROB(I) then PCTILE (gradd.f:288-322)
        else
            @inbounds for i in 1:t.n
                acc += fld[i] * t.tpa[i]
            end
        end
        # disply.f: IOSUM(k)=INT(O..CUR(7)/GROSPC+0.5) — the imperial per-acre integer (metric variants too;
        # sumout.f's m³/ha conversion is the second stage, in `metric_sumout`).
        return dt(acc / g)
    end
    # Year/age come from the cycle-boundary schedule (IY, build_cycle_schedule!): the calendar
    # year at this cycle's start, and the age advanced by the elapsed years from the inventory.
    # For uniform cycles this is exactly cycle_year[1] + cyc·per (bit-exact); non-uniform TIMEINT
    # and CYCLEAT-inserted boundaries make the steps vary.
    cyc = Int(s.control.cycle)
    yr = cycle_year_at(s.control, cyc)
    age = Int(s.plot.stand_age) + (yr - Int(s.control.cycle_year[1]))
    # RESETAGE (resage.f): after the reset year (run after DISPLY, so the reset row itself keeps
    # the old age) the stand age is rebased — age(Y) = age_reset_age + (Y − reset_year).
    ry = Int(s.control.age_reset_year)
    (ry >= 0 && yr > ry) && (age = Int(s.control.age_reset_age) + (yr - ry))
    mcuft = vtot(:merch_cuft_vol)
    # BC (metric): vols.f computes merch into WK1 but NEVER loads the MCFV summary array (imperial vols.f
    # does: MCFV=MCF), so live's OMCCUR(7)=PCTILE(MCFV)=0 ⇒ the .sum merch column (IOSUM 5) AND the
    # merch-based MAI (BCYMAI) are structurally 0 every cycle. Mirror that: zero the summary merch (the
    # per-tree merch_cuft_vol stays intact for other consumers). Board is 0 too ("not computed in this variant").
    met && (mcuft = 0)
    # MAI (BCYMAI, disply.f:383): (merch cuft + cumulative removed merch) / age.
    # `total_removed_merch` carries cross-cycle removals (0 at the inventory).
    # Computed in Float32 to match FVS REAL*4 (the %.1f rounding differs from Float64).
    # Metric final row (disply.f:516-523): the EASTERN metric variant ON loads BCYMAI from IOSUM(4) (total
    # cubic), BC from IOSUM(5) (merch, structurally 0); the per-cycle rows (evtstv.f:360-418 CASE DEFAULT)
    # use IOSUM(5) for both. BCYMAI = (IOSUM+TOTREM)/AGE in REAL, left imperial here.
    mvol = (met && final_row && s.variant isa Ontario) ? vtot(:cuft_vol) : mcuft
    # After a RESETAGE that rebased the age to ZERO, FVS shuts off MAI (disply.f:391-394 BCYMAI=0 when
    # MAIFLG≠0; evtstv.f:396 sets it when ZERO=age−period==0, i.e. the age was reset to 0, and persists it).
    # A RESETAGE to a NON-zero age keeps MAI on (e.g. s17_managed resets to 40 → MAI stays 62.5). Non-RESETAGE
    # runs have ry<0 ⇒ untouched (bit-exact); bare-ground (NEWSTD) has no RESETAGE ⇒ its own MAI path unchanged.
    mai = (mai_zero || (ry >= 0 && yr > ry && Int(s.control.age_reset_age) == 0)) ? 0f0 :
          (age > 0 ? (met ? (Float32(mvol) + Float32(total_removed_merch)) / Float32(age) :
                            Float32(mcuft + total_removed_merch) / Float32(age)) : 0f0)
    SummaryRow(
        year = yr, age = age, tpa = tpa,
        ba = ba, sdi = sdi, ccf = ccf, topht = toph, qmd = qmd,
        cuft = vtot(:cuft_vol), mcuft = mcuft,
        # ON metric .sum reports NMV (Mowraski net merch, IOSUM 6) in the board column via
        # bdft_vol; BC leaves board 0 ("not computed in this variant"). Both keep merch (IOSUM 5)
        # 0 at inventory (metric vols.f never loads the MCFV summary array — the shared `met` quirk).
        scuft = met ? 0 : vtot(:saw_cuft_vol),
        bdft = (s.variant isa Ontario) ? vtot(:bdft_vol) : (met ? 0 : vtot(:bdft_vol)),
        at_ba = ba, at_sdi = sdi, at_ccf = ccf, at_topht = toph, at_qmd = qmd,
        period = period, mai = mai,
        accretion = trunc(Int, accretion + 0.5), mortality = trunc(Int, mortality + 0.5),
        fortype = Int(s.plot.forest_type), sizecls = Int(s.plot.size_class),
        stockcls = Int(s.plot.stocking_class))
end
