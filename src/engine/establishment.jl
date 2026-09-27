# =============================================================================
# establishment.jl — regeneration / establishment (ESNUTR → ESTAB)
#
# Ported from: base/esnutr.f (cycle hook) + base/estab.f (tree creation) +
# base/estab_helpers.f (ESSUBH/ESTIME) + base/esinit.f.
#
# SN's PARTIAL (keyword-driven) establishment model: no auto-ingrowth. When an
# ESTAB packet scheduled PLANT(430)/NATURAL(431) activities are due, ESTAB creates
# the regen trees. A single bare plot (NPTIDS=1) is replicated MINREP=50 times: each
# replicate independently draws an established HEIGHT per species (ESSUBH height-at-
# age + a BACHLO draw on the establishment RNG), and contributes one record per
# species carrying plantedTPA·survival/100 / dupnpt TPA (400/50 = 8) ⇒ 50×2 = 100
# records, 800 TPA. The trees enter AFTER growth+mortality (GRADD order) so they're
# fresh (full TPA) this period; their DBH is derived from the established height.
# =============================================================================

# =============================================================================
# es_pasmax_xcsmax — PASSALL / PASMAX excess-tree cap (estab.f:1288,1318-1321).
#
# When a plot's per-species EXCESS regeneration count (integer-valued, incremented by
# 1.0 per non-"best" tree) is passed to the tree list it is broken into IBRKUP records
# (≈1 record per 5 excess trees), each carrying the SAME per-record excess XCSMAX. The
# PASSALL keyword's PASMAX caps that per-record count so the TOTAL excess passed for the
# species cannot exceed PASMAX (esin.f echo: "MAXIMUM NUMBER OF EXCESS TREES PASSED PER
# PLOT PER SPECIES"). The record PROB scales LINEARLY with XCSMAX, so this is a purely
# DETERMINISTIC post-draw scaling — it consumes NO RNG draw.
#
#   IBRKUP = INT(EXCESS/5.0 + 1.0)          (estab.f:1288)
#   XCSMAX = EXCESS/BRKUP                    (estab.f:1318)
#   FTEMP  = PASMAX/BRKUP                    (estab.f:1319)
#   IF (XCSMAX > FTEMP) XCSMAX = FTEMP       (estab.f:1320-1321)  ← the cap
#
# Returns (xcsmax, ibrkup). Faithful transcription; VALIDATED BIT-EXACT against the live
# instrumented FVSie oracle (single-.o estab.f dump swap) across EXCESS 1..6 × PASMAX
# {1,5}: e.g. (EXCESS=6,PASMAX=5)→BRK=2,XCS=2.5 (uncapped 3.0 → capped), (EXCESS=5,
# PASMAX=1)→BRK=2,XCS=0.5, (EXCESS=4,PASMAX=5)→BRK=1,XCS=4.0 (uncapped, 4<5). See
# test/unit/test_estab_passall.jl. Default PASMAX=5 (CONFID, esinit.f:50); the AUTOES
# ingrowth path overrides PASMAX=15 (estab.f:249, NTALLY==99). This kernel is the ported,
# oracle-validated primitive; the live IE tally seam is not yet wired to it (the excess
# TPA effect is cornered by the placeholder-height / probabilistic-excess establishment
# approximation — see the PASSALL branch in keyword_dispatch.jl:kw_estab!).
# =============================================================================
function es_pasmax_xcsmax(excess::Real, pasmax::Real)
    exc = Float32(excess)
    ibrkup = trunc(Int, exc / 5f0 + 1f0)          # INT(EXCESS/5.0 + 1.0)
    brk = Float32(ibrkup)
    xcsmax = exc / brk                            # EXCESS/BRKUP
    ftemp = Float32(pasmax) / brk                 # PASMAX/BRKUP
    xcsmax > ftemp && (xcsmax = ftemp)            # cap
    return xcsmax, ibrkup
end

const _ES_MINREP = 50          # MINREP: target plot replication (esinit.f) — the DEFAULT for
                               # Establishment.minrep (state.jl); the live value is per-stand via the
                               # MINPLOTS keyword (esin.f opt 20), read at both idup call sites as est.minrep.

# XMIN: per-species establishment min height (blkdat.f) lives in
# data/southern/species_coefficients.csv as the `estab_min_ht` column.
# HHTMAX: per-species max establishment height (blkdat.f).
const _ES_HHTMAX = Float32[23.0,27.0,21.0,21.0,22.0,20.0,24.0,18.0,18.0,17.0,22.0,
    (20.0 for _ in 12:90)...]

# NE HHTMAX: per-species MAX established height (ne/blkdat.f DATA HHTMAX/, 108 values). A HARD cap on the
# grown establishment height (not the soft site-curve HTMAX) — e.g. YB(30)=22, WO(55)=16 are reached exactly
# by live (all trees clamped). The SN _ES_HHTMAX above is wrong for NE (it uses a 20-ft fill for sp≥12).
const _NE_ES_HHTMAX = Float32[
    20,24,18,16,18,16,16,18,20,14, 14,16,16,16,16,16,16,16,14,14,
    16,18,12,20,16,20,16,16,18,22, 20,18,18,18,14,14,14,14,18,14,
    24,24,18,24,28,24,18,20,20,24, 24,20,20,26,16,14,12,12,16,16,
    14,16,14,16,16,12,20,16,16,14, 16,12,12,12,12,12,12,18,20,12,
    20,20,20,20,16,16,16,24,14,24, 32,18,16,16,16,16,12,10,16,18,
    30,20,20,18,16,20,20,30]

# NE ESSUBH per-species reference age CARAGE (essubh.f DATA MAPNE/, 108 values — DISTINCT from the htcalc
# curve-index MAPNE). The planted base height is (NC-128 height at this age / this age) · min(5, TIME−DELAY).
const _NE_ESSUBH_REFAGE = Int[
    20,10,15,20,15,20,20,20, 5,20, 15,20,20,10,20,20,20,20,20,10,
    10,15,20,15,10,20,20,20,20,20, 20,20,20,20,20,20,20,20,20,20,
    20,20,20,35,35,20,10,20,20,20, 10,20,15,20,10,10,10,10,10,10,
    10,30,10,10,30,30,20,10,10,20, 10,10,10,20,10,10,10,20,20,10,
    25,25,10,25,25,10,10,20,10,10, 10,10,20,20,20,20,20,10,10,10,
    10,10,10,10,10,10,10,10]

# CS ESSUBH per-species reference age CARAGE (cs/essubh.f DATA MAPCS/, 96 values — DISTINCT from the
# htcalc curve-index MAPCS). Same NE-style base height = (NC-128 height at CARAGE / CARAGE)·min(5, TIME−DELAY).
const _CS_ESSUBH_REFAGE = Int[
    10,10,10,15,20,10,5,20,20,25, 25,25,25,20,20,20,20,20,20,20,
    20,10,20,20,20,35,20,15,20,10, 15,20,20,20,10,20,20,20,20,20,
    20,20,20,20,20,35,10,20,10,10, 10,10,10,30,10,30,10,10,10,10,
    20,10,35,30,10,10,20,10,10,20, 20,10,10,20,20,20,10,20,20,20,
    20,10,10,10,10,10,10,10,10,20, 10,20,25,20,10,10]

# CS planted/regen height cap HHTMAX (cs/blkdat.f DATA HHTMAX/, 96 species). Clamps the REPORTED
# seedling height (estab.f:496 / esgent.f:64); the DBH is taken from the UNCAPPED grown height.
const _CS_ES_HHTMAX = Float32[
    16,27,14,14,14,16,20,20,18,16, 20,20,16,14,14,14,18,14,14,14,
    14,14,14,14,18,28,20,24,20,16, 18,26,16,14,12,20,16,20,12,20,
    24,16,16,24,24,24,16,20,16,16, 16,20,12,16,14,12,12,20,16,20,
    14,14,20,16,20,14,20,20,18,20, 20,12,20,24,20,20,24,20,24,20,
    18,18,20,32,10,20,20,18,16,20, 12,20,20,20,20,16]

# LS planted/regen height cap HHTMAX (ls/blkdat.f DATA HHTMAX/, 68 species) + ESSUBH reference age CARAGE
# (ls/essubh.f DATA MAPLS/, 68 species — DISTINCT from the htcalc MAPLS curve map).
const _LS_ES_HHTMAX = Float32[14,20,18,18,20,18,18,20,16,24,16,16,16,16,18,24,24,18,20,26,16,12,20,22,16,16,16,14,24,16,16,14,12,20,16,20,20,14,14,20,20,24,18,20,18,20,20,24,10,16,18,20,20,20,12,18,16,20,16,24,30,20,20,20,32,20,18,20]
const _LS_ESSUBH_REFAGE = Int[20,15,20,20,5,15,15,20,20,10,20,20,10,10,20,35,15,15,20,20,20,20,20,20,20,20,20,20,20,10,30,10,10,20,10,10,20,20,20,20,20,20,20,20,20,20,10,10,10,10,10,10,10,20,10,10,10,10,25,20,10,10,10,10,10,10,10,10]

# Establishment min-height (XMIN) + max seedling height (HHTMAX) per species — from each variant's
# blkdat.f (VERIFIED: IE blkdat.f:62 XMIN == _IE_ES_XMIN). EM/BM/UT/CI had no establishment.jl ⇒ the dispatch
# fell to the missing `:estab_min_ht` coef ⇒ KeyError crash on ESTAB/PLANT-keyword stands (full utt01/emt01/
# bmt01). NOTE: blkdat.f XMIN (establishment), NOT regent.f XMIN (small-tree/regen — a DIFFERENT array).
const _EM_ES_XMIN   = Float32[1,1,1,1,0.5,0.5,1,0.5,0.5,1,3,6,3,3,3,3,6,0.5,3]
const _EM_ES_HHTMAX = Float32[23,27,21,27,18,6,24,18,18,17,16,16,16,16,16,16,16,22,16]
const _BM_ES_XMIN   = Float32[0.9,1.7,1,1,0.5,0.5,1.3,0.5,0.5,1,1,1,1,1,6,1,1,1]
const _BM_ES_HHTMAX = Float32[23,27,21,21,22,6,24,18,18,17,23,9,20,20,16,20,17,20]
const _UT_ES_XMIN   = Float32[1,1,1,0.5,0.5,6,1,0.5,0.5,1,0.5,0.5,0.5,0.5,0.5,0.5,0.5,3,3,0.5,0.5,3,0.5,0.5]
const _UT_ES_HHTMAX = Float32[9,9,10,7,7,16,10,7,7,10,6,6,10,6,6,6,9,16,16,6,6,16,9,10]
const _CI_ES_XMIN   = Float32[1,1,1,0.5,0.5,0.5,1,0.5,0.5,1,1,1,6,0.5,0.5,1,3,0.5,3]
const _CI_ES_HHTMAX = Float32[23,27,21,21,22,20,24,18,18,17,27,27,16,6,6,27,16,22,16]
# EC (East Cascades, MAXSP=32) establishment min height XMIN / max seedling height HHTMAX — ec/blkdat.f:140,153.
const _EC_ES_XMIN   = Float32[1,1,1,0.5,0.5,0.5,1,0.5,0.5,1,1,0.5,1,1,1,0.5,1.5,1,1,1,1,1,1,1,1,1,1,1,1,1,0.5,1]
const _EC_ES_HHTMAX = Float32[23,27,21,21,22,20,24,18,18,17,20,22,20,20,20,20,20,20,20,20,20,50,20,20,20,20,20,20,20,20,22,20]
# NC (Klamath) establishment per-species min height (XMIN) / max sprout height (HHTMAX) — nc/blkdat.f:73,80.
const _NC_ES_XMIN   = Float32[1,1,1,0.5,1,0.5,0.5,1,0.5,1,1,1]
const _NC_ES_HHTMAX = Float32[27,31,25,25,26,24,28,20,20,18,26,25]
# NC subsequent/planted base height (nc/essubh.f): a FIXED per-species table (no age/site/EMSQR),
# clamped [XMIN,HHTMAX] by the shared engine (like UT/TT).
const _NC_ESSUBH_HHT = Float32[1,1,1,1,7,1,7,7,1,0.8,7,2]
# OP (Olympic, MAXSP=39) establishment min height XMIN / max sprout height HHTMAX — op/blkdat.f:69-72
# (DATA XMIN/ 1.0, 2*1.5, 7*1.0, 1.4, 3*1.0, 1.3, 1.5, 13*1.0, 1.5, 9*1.0 /  and  HHTMAX/ 21*20.0,50.0,17*20.0 /).
# These are the ESCOMN establishment arrays (NOT regent.f XMIN). OP essubh.f builds its planted-seedling
# height from op/smhgdg.f (Gould–Harrington small-tree height growth), so there is no fixed ESSUBH_HHT table.
const _OP_ES_XMIN   = Float32[1.0, 1.5, 1.5, 1,1,1,1,1,1,1, 1.4, 1,1,1, 1.3, 1.5,
                              1,1,1,1,1,1,1,1,1,1,1,1,1, 1.5, 1,1,1,1,1,1,1,1,1]
const _OP_ES_HHTMAX = Float32[fill(20f0, 21)..., 50f0, fill(20f0, 17)...]

# BM (Blue Mountains) estab-model species → SUMSP slot map (esaddt.f:150-167): the .es1 payload
# reports per-species TPA<1" for these eight species. Non-BM variants just report 0 in every slot
# (the payload only matters to a real external BM model; the bridge's A/B rides on the .es2).
const _ES_ADDT_BMMAP = Dict{Int,Int}(4=>1, 9=>2, 2=>3, 7=>4, 1=>5, 8=>6, 10=>7, 5=>8)

# esaddt.f:126-176 — write the `.es1` stand summary the external regeneration model reads. The
# header (stand id, planting year IPYR, ADDTREES field-2 offset) is exact; the density block is a
# best-effort port of the Blue-Mountains payload (habitat/slope/aspect/elevation, pre/post SDI & BA
# over GROSPC, stand TPA>1", per-species TPA<1"). Fortran right-justifies A30/I30/F30.1 fields.
function write_es1(s::StandState, path::AbstractString, sid::AbstractString, ipyr::Int, iyr1::Int32)
    t = s.trees; p = s.plot
    grospc = p.gross_space > 0f0 ? p.gross_space : 1f0
    sumsp = zeros(Float32, 8); sum1 = 0f0
    @inbounds for i in 1:t.n
        if t.dbh[i] > 1f0
            sum1 += t.tpa[i]
        else
            slot = get(_ES_ADDT_BMMAP, Int(t.species[i]), 0)
            slot > 0 && (sumsp[slot] += t.tpa[i])
        end
    end
    open(path, "w") do io
        println(io, lpad(strip(String(sid)), 30))              # (1) StandID (A30)
        println(io, lpad(string(ipyr), 30))                    # (2) planting year IPYR (I30)
        println(io, lpad(string(Int(iyr1)), 30))               # (3) ADDTREES fld2 (I30)
        println(io, lpad(string(Int(p.habitat_code)), 30))     # (4) hab code KODTYP
        println(io, lpad(string(Int(p.slope_raw)), 30))        # (5) slope ISLOP
        println(io, lpad(string(Int(p.aspect_deg)), 30))       # (6) aspect IASPEC
        _es1_f(io, p.elevation * 100f0)                        # (7) elev*100
        _es1_f(io, p.sdi_before_cut / grospc)                  # (8) bsdi
        _es1_f(io, p.sdi_after_cut / grospc)                   # (9) asdi
        _es1_f(io, p.old_ba / grospc)                          # (10) ba
        _es1_f(io, p.at_ba / grospc)                           # (11) aba
        _es1_f(io, sum1)                                       # (12) stand TPA >1"
        for k in 1:8; _es1_f(io, sumsp[k]); end                # (13-20) spp TPA <1"
    end
    return nothing
end
_es1_f(io, x::Real) = println(io, lpad(string(round(Float32(x); digits=1)), 30))

# esaddt.f:181-195 + base/oprdat.f — read the external model's `.es2` activity block and OPADD-
# schedule its activities. Layout: line 1 = IKEEP (esaddt reads it, I10); then OPRDAT scans for a
# line equal to the stand id (NPLT), and reads `IACTK IDT NPRMS PRMS(1..NPRMS)` records (free-form)
# until an `End` line. Records with IDT ≥ the current cycle-start year are scheduled (oprdat.f:55).
# 431/430 (NATURAL/PLANT) route into the already-validated regen path (the `due` filter in
# establish!). Returns IKEEP so the caller can honor the delete/keep of the file (esaddt.f:191-195).
function read_es2!(s::StandState, path::AbstractString, sid::AbstractString)::Int
    lines = readlines(path)
    isempty(lines) && return 0
    ikeep = something(tryparse(Int, strip(lines[1])), 0)
    target = strip(String(sid))
    icyc_year = Int(current_cycle_year(s))                     # IY(ICYC)
    i = 2; n = length(lines); nadd = 0
    while i <= n
        if strip(lines[i]) == target
            i += 1
            while i <= n && !startswith(lstrip(lines[i]), "End")
                toks = split(strip(lines[i]))
                if length(toks) >= 3
                    iactk = tryparse(Int, toks[1]); idt = tryparse(Int, toks[2])
                    nprms = tryparse(Int, toks[3])
                    if iactk !== nothing && idt !== nothing && nprms !== nothing
                        prms = Float32[something(tryparse(Float32, toks[3+j]), 0f0)
                                       for j in 1:min(nprms, length(toks) - 3)]
                        # oprdat.f:55 IF (IDT.GE.IY(ICYC)) — schedule only current-or-future dates
                        # (cycle-number dates <1000 are relative, always kept for the due filter).
                        if idt >= icyc_year || (0 < idt < 1000)
                            pr = ntuple(k -> k <= length(prms) ? prms[k] : 0f0, 6)
                            push!(s.control.schedule, ScheduledActivity(Int32(idt), Int32(iactk), pr))
                            nadd += 1
                        end
                    end
                end
                i += 1
            end
            break
        end
        i += 1
    end
    return ikeep
end

# esnutr.f:59 CALL ESADDT(1) — the ADDTREES external-regeneration bridge (estb/esaddt.f). For each
# scheduled activity-432 due this cycle: write the `.es1` summary, run the external command, read the
# `.es2` activity block, and OPADD-schedule its activities. Runs at the top of the ESNUTR seam so the
# `due` filter in establish! creates the returned PLANT/NATURAL regen the SAME cycle.
function addtrees_bridge!(s::StandState, yr::Int32, per::Int)
    fvscyc = Int(s.control.cycle) + 1
    kdt = cycle_year_at(s.control, Int(s.control.cycle) + 1) - 1   # IY(ICYC+1)-1 (end of cycle)
    stem = replace(s.control.keyword_file, " " => "")
    sid  = strip(s.plot.stand_id)
    var  = replace(String(s.control.variant_code), " " => "")
    for at in s.estab.addtrees
        at.fired && continue
        idt = Int(at.idt)
        # OPFIND mode-1: the 432 fires in the cycle its date falls in (calendar year OR cycle number).
        due = (yr <= idt < yr + per) || (0 < idt < 1000 && idt == fvscyc)
        due || continue
        at.fired = true
        isempty(strip(at.cmdln)) && continue
        ipyr = Int(at.iyr1) + kdt
        # esaddt.f:73-107 — filename stem <KWDFIL>_<NPLT>_<KDT>_<VARACD>, spaces removed from each part.
        base = string(stem, "_", replace(String(sid), " " => ""), "_", kdt, "_", var)
        es1 = base * ".es1"; es2 = base * ".es2"
        write_es1(s, es1, sid, ipyr, at.iyr1)
        # esaddt.f:116-124,177 — the command is the CMDLN with the .es1 filename appended, then SYSTEM.
        cmd = string(at.cmdln, " ", es1)
        try
            run(pipeline(Cmd(String.(split(cmd))); stdout=devnull, stderr=devnull))
        catch
            # A failed external command mirrors FVS's CALL SYSTEM returning nonzero: the .es2 may be
            # absent (nothing scheduled). Do not abort the run.
        end
        isfile(es2) || continue
        ikeep = read_es2!(s, es2, sid)
        ikeep == 1 || (try; rm(es2); catch; end)   # esaddt.f:191-195 IKEEP≠1 ⇒ delete
    end
    return nothing
end

"""
    snapshot_esb_inputs!(state)

Freeze the AUTOES ESB1 BAAOLD at stand SETUP, mirroring FVS `base/fvs.f:201 CALL ESFLTR`
(which runs once, before any growth cycle): the inventory per-point OVERSTORY (D≥REGNBK=2.999)
basal area BAAINV(1), used as ESB1 = ESTOCK(BAAOLD=BAAINV). ESFLTR is never re-called after a
harvest, so this stays the INVENTORY value even when the first regen tally fires after a cycle-2
thin — the key to matching the post-thin re-stocking magnitude. (ESB's small-tree TPACRE is NOT
frozen — estab.f:301-322 reads the live tally-time small-tree TPA.) IE/EM only; no-op otherwise.
Idempotent. The per-point BA uses the same PTBAA scale as `point_basal_area!` (validated vs live
BAAA), filtered to point-1 overstory records.
"""
function snapshot_esb_inputs!(s::StandState)
    (s.variant isa InlandEmpire || s.variant isa EasternMontana) || return s
    isnan(s.estab.inv_baaold) || return s              # snapshot once (setup)
    p, t = s.plot, s.trees
    scale = p.gross_space > 0f0 ? p.pi / p.gross_space : 1f0   # PTBAA scale (= point_basal_area!)
    # PER-INVENTORY-POINT OVERSTORY BAAINV(NNID) (esfltr.f:67, D≥REGNBK): each stockable point accumulates its own
    # frozen inventory overstory BA. ESB1(NCOUNT)=ESTOCK(BAAINV(NNID)) is per-point in estab.f, so the ESB−ESB1
    # actual-vs-predicted stocking correction is per-point — a single stand/point-1 value over-corrects the open
    # points of a heterogeneous multi-point stand (dense point 1 ⇒ +4 logit ⇒ the open points' ingrowth PROB1
    # saturates ⇒ ~4× AUTOES over-production on M333 subalpine stands). MEASURED FVSie_g16 1856003217290487.
    nptids = max(1, Int(p.points_inv) - Int(p.nonstockable))
    ptbaold = zeros(Float32, nptids)
    # esfltr.f:61-68 (record order): PIX=PI−FLOAT(NONSTK); BAAINV(N)=BAAINV(N)+0.005454154*D*D*ZPROB*PIX, evaluated
    # left to right (((c·D)·D)·PROB)·PIX — MEASURED FVSie_g16 3356357010690 point 4 BAAINV 4397AE6E; the former
    # P·c·D·D·(PI/GROSPC) form rounded ESB1(4) 4 ULP off (C0084734 vs C0084738) ⇒ PROB1 ⇒ every cohort record's TPA.
    pix = p.pi - Float32(p.nonstockable)
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d >= 2.999f0 || continue                              # OVERSTORY (D≥REGNBK)
        pid = Int(t.plot_id[i])
        (1 <= pid <= nptids) && (ptbaold[pid] += 0.005454154f0 * d * d * t.tpa[i] * pix)
    end
    s.estab.inv_point_baaold = ptbaold
    s.estab.inv_baaold = ptbaold[1]                           # point-1 value (scalar path / EM)
    return s
end

"""
    establish!(state; fint=5f0) -> Bool

Create scheduled PLANT/NATURAL regen for the current cycle (ESNUTR/ESTAB). Runs at
the end of `grow_cycle!` (GRADD order). Idempotent per year. Returns whether any
tree was created. No-op unless an ESTAB packet is active.
"""
# estb_planted_height — IE/EM (estb/estab.f) DO 322 height of ONE PLANT/NATURAL keyword tree (estab.f:1010-1041):
# ESSUBH base height from the plot's EMSQR and DILATE=FIRST(2,sp), then a user height TREEHT≥0.1 replaces it with a
# lognormal BACHLO(ln TREEHT, 0.5) redrawn until within [0.5,2]·TREEHT (:1026-1034, :estab draws from the CURRENT
# s.rng.es0 — the caller positions it at the plot's post-ESAVE state), +HTADJ, floor 0.05 (else +HTADJ, floor XMIN),
# cap HHTMAX. Shared by the AUTOES tally (NBEST ranks the planted trees at these heights, :1090-1144) and by
# establish! (which books them), so both see identical values.
function estb_planted_height(s::StandState, a, per::Int, yr::Int, emsqr::Float32, dil::Float32)::Float32
    sp = round(Int, a.params[1])
    pyr = (0 < Int(a.year) < 1000) ? Int(cycle_year_at(s.control, Int(a.year) - 1)) : Int(a.year)
    delay = pyr - yr
    gentim = max(per - 5, 0)
    trage = a.params[4] < 0.5f0 ? 2f0 : a.params[4]; trage > 10f0 && (trage = 10f0)
    age = Float32(per) - Float32(delay) - Float32(gentim) + trage; age < 1f0 && (age = 1f0)
    slo = s.plot.slope
    hht = if s.variant isa InlandEmpire
        iage = trunc(Int, (Float32(per) - Float32(delay) - Float32(gentim)) + 0.5f0)
        iage < 1 && (iage = 1); iage > 20 && (iage = 20)
        ie_essubh(sp, age, clamp(s.plot.basal_area, 1f0, 400f0), em_ihtser(Int(s.plot.habitat_code)), 1, 3,
                  slo*cos(s.plot.aspect), slo*sin(s.plot.aspect), slo, s.plot.elevation, emsqr * dil * _IE_ES_BNORML[iage])
    else
        em_essubh_hht(sp, flog(age), clamp(s.plot.basal_area, 1f0, 400f0), slo*cos(s.plot.aspect), slo*sin(s.plot.aspect),
                      slo, s.plot.elevation, em_ihtser(Int(s.plot.habitat_code)), 3, 1)
    end
    hadj = isempty(s.estab.ht_adj) ? 0f0 : get(s.estab.ht_adj, Int32(sp), 0f0)
    treeht = a.params[5]
    if treeht >= 0.1f0
        hht = treeht; xh = flog(hht)
        while true
            xxh = fexp(bachlo(s.rng, xh, 0.5f0; stream = :estab))
            (0.5f0 * hht <= xxh <= 2f0 * hht) && (hht = xxh; break)
        end
        hht += hadj; hht < 0.05f0 && (hht = 0.05f0)
    else
        hht += hadj
        xmn = s.variant isa InlandEmpire ? _IE_ES_XMIN[sp] : _EM_ES_XMIN[sp]
        hht < xmn && (hht = xmn)
    end
    hmx = s.variant isa InlandEmpire ? _IE_ES_HHTMAX[sp] : _EM_ES_HHTMAX[sp]
    hht > hmx && (hht = hmx)
    return hht
end

# strp estab.f:174-270 (TT/UT) head of an ESTAB call with tally number `ntally`: NTALLY==1 draws a fresh ESDRAW off the
# running ESRANN stream (:174-178), every call reseeds ESRNSD(ESDRAW) (:179) and fills WK6 (:202-205). On NTALLY==1 the
# per-plot site prep IPPREP (:207-270) is sampled without replacement off WK6 — MECHPREP/BURNPREP %-of-plots (esetpr.f)
# or, with neither keyword, strp/esprep.f's flat 0.75/0.20/0.05 — and later tallies reuse it. Each point's replicates
# are processed grouped by prep type (DO 202 ITYPEP), so replicate slot k of a point takes the k-th smallest prep of its
# block (ie_esetpr_sample). Only TT's PP essubh reads it (UPRE(IPREP), tt/essubh.f:57,123).
function _strp_estab_rng!(s::StandState, ntally::Integer, nptids::Integer, idup::Integer, yr::Integer, per::Integer)
    ntally == 1 && (s.estab.es_seed = floor(esrann!(s.rng) * 100000f0 + 0.5f0))
    esd = s.estab.es_seed
    (esd % 2f0 == 0f0) && (esd += 1f0)                                   # ESRNSD odd-force (esrann.f:56)
    s.rng.es0 = Float64(esd)
    wk6 = Float32[esrann!(s.rng) for _ in 1:(nptids * idup)]
    if ntally == 1 && s.variant isa Teton
        pmech_pct = nothing; pburn_pct = nothing
        for a in s.control.schedule
            (a.icflag == Int32(493) || a.icflag == Int32(491)) || continue
            ay = Int(a.year)
            idt = (0 < ay < 1000) ? Int(cycle_year_at(s.control, ay - 1)) : ay
            (Int(yr) <= idt < Int(yr) + per) || continue
            a.icflag == Int32(493) ? (pmech_pct = a.params[2]) : (pburn_pct = a.params[2])
        end
        sumup = if pmech_pct !== nothing || pburn_pct !== nothing
            es = ie_esetpr(pmech_pct, pburn_pct)
            ie_esetpr_normalize(0f0, es.pmech, es.pburn, es.ialn2, es.ialn3)
        else
            ie_esetpr_normalize(0.75f0, 0.20f0, 0.05f0, 0, 0)
        end
        s.estab.es_ipprep = Int32.(ie_esetpr_sample(sumup, wk6, nptids, idup))
    end
    return s
end

# estab.f:1496-1505 (IE/EM, after ESGENT): ABIRTH(I)=ABIRTH(I)+GENTIM for every record established this cycle,
# GENTIM being estab.f's final value this cycle (the last PLANT's :1053-1058 reset, else :448's FINT−5).
function esgent_add_gentim!(s::StandState, nstart::Int, fint::Float32)
    gentim = s.estab.gentim_cyc == Int32(s.control.cycle) ? s.estab.gentim_post : max(fint - 5f0, 0f0)
    @inbounds for i in (nstart+1):s.trees.n; s.trees.birth_age[i] += gentim; end
    return s
end

function establish!(s::StandState; fint::Float32 = 5f0)::Bool
    s.estab.active || return false
    # AK (SoutheastAlaska) has its own full establishment model (ak/estab.f, 2477 lines, + ak/esblkd.f habitat
    # tables, ak/essubh.f, ak/esgent.f) that jl has not ported; the shared path would read eastern `estab_min_ht`
    # AK doesn't carry (akt01's PLANT stand). Leave ESTAB/PLANT inert for AK, loudly, rather than erroring.
    if s.variant isa SoutheastAlaska
        s.estab.active = false
        @warn "AK establishment (ak/estab.f) is not ported — ESTAB/PLANT/NATURAL are inert; the stand diverges from FVSak"
        return false
    end
    t = s.trees; sd = s.coef.species
    es_xmin = s.variant isa CentralRockies ? _CR_ES_XMIN :
              s.variant isa InlandEmpire ? _IE_ES_XMIN :
              s.variant isa Teton ? _TT_ES_XMIN :
              s.variant isa EasternMontana ? _EM_ES_XMIN :
              s.variant isa BlueMountains ? _BM_ES_XMIN :
              s.variant isa Utah ? _UT_ES_XMIN :
              s.variant isa CentralIdaho ? _CI_ES_XMIN :
              s.variant isa EastCascades ? _EC_ES_XMIN :
              s.variant isa Klamath ? _NC_ES_XMIN :
              s.variant isa Olympic ? _OP_ES_XMIN :
              (s.variant isa WestCascades || s.variant isa PacificNorthwest) ? _OP_ES_XMIN :   # wc/pn blkdat.f:70-71 XMIN = op/blkdat.f DATA
              sd[:estab_min_ht]   # per-species establishment min height (eastern SN/NE/CS/LS have this column)
    es_hhtmax = s.variant isa Northeast ? _NE_ES_HHTMAX :
                s.variant isa CentralStates ? _CS_ES_HHTMAX :
                s.variant isa LakeStates ? _LS_ES_HHTMAX :
                s.variant isa CentralRockies ? _CR_ES_HHTMAX :
                s.variant isa InlandEmpire ? _IE_ES_HHTMAX :
                s.variant isa Teton ? _TT_ES_HHTMAX :
                s.variant isa EasternMontana ? _EM_ES_HHTMAX :
                s.variant isa BlueMountains ? _BM_ES_HHTMAX :
                s.variant isa Utah ? _UT_ES_HHTMAX :
                s.variant isa CentralIdaho ? _CI_ES_HHTMAX :
                s.variant isa EastCascades ? _EC_ES_HHTMAX :
                s.variant isa Klamath ? _NC_ES_HHTMAX :
                s.variant isa Olympic ? _OP_ES_HHTMAX :
                (s.variant isa WestCascades || s.variant isa PacificNorthwest) ? _OP_ES_HHTMAX : _ES_HHTMAX   # wc/pn blkdat.f:73 = op
    per = round(Int, fint)
    yr = Int32(current_cycle_year(s))   # IY schedule; yr+per below = next boundary (fint is per-cycle)
    # esnutr.f:59 CALL ESADDT(1): the ADDTREES external-regen bridge runs at the TOP of the ESNUTR
    # seam, before the tally. It may OPADD-schedule PLANT/NATURAL activities (from the external
    # model's .es2) that the `due` filter below then creates THIS cycle. Runs before the years_done
    # guard so a 432 fires exactly once per cycle it is due, independent of the regen idempotency.
    isempty(s.estab.addtrees) || addtrees_bridge!(s, yr, per)
    yr in s.estab.years_done && return false
    # PLANT/NATURAL dates < 1000 are CYCLE NUMBERS (FVS 1-based), not calendar years — the same OPNEW/OPFIND
    # convention cuts! applies (cuts.jl:203-208). Without the cycle-number clause a `PLANT 2 ...` (cycle 2) never
    # matched `yr <= a.year` (2016 <= 2 is false) ⇒ planting silently never fired (bit vs live only on real FIA
    # stands scheduled by cycle; thin/salvage already resolved this way, ESTAB was the omission).
    fvscyc = Int(s.control.cycle) + 1
    due = [a for a in s.control.schedule
           if (a.icflag == Int32(430) || a.icflag == Int32(431)) &&
              ((yr <= a.year < yr + per) || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc))]

    # NPTIDS = IPTINV − NONSTK (esplt2.f:74): the STOCKABLE inventory points, not the raw
    # plot count. Driving DUPNPT/IDUP and so the regen record count + its per-record RNG draws.
    nptids = max(1, Int(s.plot.points_inv) - Int(s.plot.nonstockable))
    # estab.f:199-207: IDUP = smallest I with NPTIDS·I ≥ MINREP = CEIL(MINREP/NPTIDS) (not floor); the
    # MAXPLT cap doesn't bind for the divergent 1<NPTIDS<MINREP cases. NPTIDS=1 ⇒ ceil=floor=50 (BARE stand).
    idup   = max(1, cld(Int(s.estab.minrep), nptids))   # MINPLOTS keyword (esin.f MINREP; default 50)
    dupnpt = Float32(nptids * idup)
    # strp ESNUTR/ESTAB (TT/UT): ESTAB runs whenever esnutr.f schedules a tally (estab_prep_esnutr! →
    # est.cyc_ntally), not only when a PLANT/NATURAL is due. A new ESTAB date restarts NTALLY=1 (fresh ESDRAW off the
    # running ESRANN stream + fresh site preps); the ≤19-yr continuation (NTALLY+1) with nothing to plant still
    # reseeds with ESDRAW and burns the WK6 fill and every plot's EMSQR/ESAVE draws (estab.f:174-205,394-417,528),
    # which moves the stream the NEXT NTALLY=1 draws its ESDRAW from. MEASURED FVStt_g16 (ESTAB 1992 + ESTAB 2012
    # PLANTs): live ESTAB 2012 runs NTALLY=1 with ESDRAW 73599; jl kept counting (NTALLY=3, the 1992 seed) ⇒ every
    # 2012 seedling height off.
    strp = s.variant isa Teton || s.variant isa Utah
    strp_nt = (strp && s.estab.cyc_ntally_year == yr) ? Int(s.estab.cyc_ntally) : -1
    if isempty(due)
        if strp_nt > 0
            _strp_estab_rng!(s, strp_nt, nptids, idup, yr, per)
            for _ in 1:(nptids * idup)                               # estab.f:394-417 EMSQR ×2, ESAVE, then :528 reseed
                esrann!(s.rng); esrann!(s.rng)
                _esave = floor(esrann!(s.rng) * 100000f0 + 0.5f0)
                (_esave % 2f0 == 0f0) && (_esave += 1f0); s.rng.es0 = Float64(_esave)
            end
            s.estab.ntally = Int32(strp_nt)
            push!(s.estab.years_done, yr)
        end
        return false
    end
    # ESSUBH base height from age uses the variant's site-curve: SN Chapman-Richards (ht_curve_b*),
    # NE NC-128 (ne_htcalc_height). bc is SN-only (NE has no ht_curve_b* coefs).
    bc = (s.variant isa Northeast || s.variant isa CentralStates || s.variant isa LakeStates ||
          s.variant isa CentralRockies || s.variant isa InlandEmpire || s.variant isa Teton ||
          s.variant isa EasternMontana || s.variant isa BlueMountains || s.variant isa Utah ||
          s.variant isa CentralIdaho || s.variant isa EastCascades ||
          s.variant isa Klamath || s.variant isa Olympic || s.variant isa WestCascades ||
          s.variant isa PacificNorthwest) ? nothing :   # western variants use a fixed/XMIN base, not the SN ht-curve
         (sd[:ht_curve_b1], sd[:ht_curve_b2], sd[:ht_curve_b3], sd[:ht_curve_b4], sd[:ht_curve_b5])
    montane = !isempty(s.plot.eco_unit) && s.plot.eco_unit[1] == 'M'
    ifor = Int(s.plot.forest_idx)
    # Natural-height random-draw acceptance window (estab.f:482-483 SN/CS vs :489-490 NE). FVS draws
    # RAN~N(0.5,0.25) and REDRAWS until RAN falls in the window; the window is VARIANT-SPECIFIC:
    # NE accepts [-2.5, 2.5]; SN and CS accept [0.0, 1.5]. The narrower SN/CS window truncates the
    # tails (no negative RAN ⇒ fewer trees pinned to the XMIN floor, a longer upper tail) AND rejects
    # more draws, so it also changes how many :estab draws each replicate consumes — a different window
    # desyncs the whole establishment RNG stream (and, downstream, the shared small-tree growth RANN
    # stream), which is the bare_natural sawtimber-tail divergence (D10). jl previously hardcoded the
    # NE window on the shared path.
    # Establishment default-height RAN acceptance window (estab.f:483/490): SN = [0,1.5] (sn/estab.f:486);
    # NE, CS, AND LS all = [-2.5,2.5] (ne/cs/ls estab.f:490). The old `Northeast ? … : (0,1.5)` wrongly gave
    # CS AND LS the SN window [0,1.5], which REJECTS the low tail (RAN<0) ⇒ biased the planted-seedling
    # heights HIGH (esp. the smallest, whose small-RAN draws live accepts) — the BARE-PLANT over-sizing.
    ran_lo, ran_hi = (s.variant isa Southern || s.variant isa CentralRockies || s.variant isa InlandEmpire || s.variant isa Teton || s.variant isa Utah || s.variant isa EastCascades || s.variant isa Olympic || s.variant isa WestCascades || s.variant isa PacificNorthwest || s.variant isa BlueMountains) ? (0f0, 1.5f0) : (-2.5f0, 2.5f0)   # CR/IE/TT/UT/EC/OP/WC/PN/BM = SN window (cr/estab.f:486; ec/estab.f:486; op estab.f:486; BM strp/estab.f:486 RAN∈[0,1.5])
    # gentim/delay/trage timing (esnutr/estab/essubh): age = FINT − delay − gentim + trage.
    # estab.f:448-449 — GENTIM = FINT−5 (clamped ≥0), depends ONLY on FINT, never IDSDAT/calendar
    # year. (Was `yr − idsdat`, a confirmed bandaid B5; masked today by the es_xmin height floor.)
    gentim = max(per - 5, 0)
    # Each new regen tree's crown ratio uses the per-point CCF computed by DENSE from the EXISTING (pre-regen)
    # overstory: regent.f:178 `CR=0.89722−0.0000461·PCCF(IPCCF)` with `IPCCF=ITRE(I)` (the tree's point). We now
    # carry that exact per-point value (`density.point_ccf`, filled by `point_density!` at start-of-cycle) and
    # index it by each record's point below — replacing the prior whole-stand `stand_ccf` approximation. The
    # coefficient is tiny (4.6e-5), so a bare/sparse stand (CCF≈0) is unchanged to print resolution.
    created = false
    nstart = t.n        # tree count before establishment (phase-2 crown pass starts here)
    # IPTIDS (esplt2.f:77-131): the STOCKABLE inventory point indices, indexed by the estab outer loop
    # (estab.f:313 `ITRE=IPTIDS[nn]`) — NOT the raw loop counter nn. A nonstockable plot (its `mort_code==8`
    # ".tre" record is skipped in treeinput.jl:91) has NO stored tree, so it is absent from the overstory
    # plot_ids; the stockable points are exactly the distinct plot_ids that DO carry a record. Using nn
    # directly put regen on the nonstockable point and skipped a stockable one, reading the wrong
    # `density.point_ccf[plot_id]` ⇒ the estab_pccf 7-tree crown residual (plant_stocked point 7 nonstockable).
    # FALLBACK to nn when the count doesn't match NPTIDS (bare stands: no overstory ⇒ empty ⇒ identity nn),
    # so every no-nonstockable-point scenario stays bit-exact.
    iptids = sort(unique(Int(t.plot_id[i]) for i in 1:nstart))
    use_iptids = length(iptids) == nptids
    # REGENT-LESTB's BALMOD competition uses the PRE-establishment density — the new seedlings do NOT compete
    # in their own creation cycle (live FVSne debug: GMOD=1.0 / AVH=0 for a BARE stand; the DENSE/BAL the cycle
    # uses predates the regen). Snapshot the BAL over the existing overstory (1:nstart) NOW, before any seedling
    # is added; computing it AFTER (over the cohort, the old code) over-counted the seedlings' own BA and
    # under-grew the established cohort ~4% (dbh 1.12 vs live 1.17 ⇒ the cyc-1 SDI/CCF deficit).
    ebau_pre = zeros(Float32, 50)
    s.variant isa Northeast && ne_badist!(ebau_pre, s)
    # LS ls_balmod reads RMSQD from the DENSE common, which for the establishment cohort is the
    # PRE-establishment stand QMD (live FMEFF stamp: BARE stand → RMSQD=0, so ls_balmod takes the
    # rmsqd≤0 → omega=b4 branch, GM 0.745 for jack pine). Snapshot it BEFORE the seedlings are added;
    # `stand_qmd(s)` recomputed after would include the cohort (0.626) and flip ls_balmod to the else
    # branch (GM 1.0) ⇒ the cohort over-grows (the BARE-PLANT seedling over-sizing).
    rmsqd_pre = stand_qmd(s)
    # CS REGENT-LESTB BALMOD needs the POST-growth OVERSTORY BA/AVH (pre-seedling) — same snapshot pattern as
    # ebau_pre/rmsqd_pre above. plot.basal_area/avg_height are STALE here (set pre-growth at cycle start;
    # compute_density! refreshes them only AFTER establish!, simulate.jl:463), and stand_ba/stand_top_height
    # computed later in the phase-2 loop would over-count the new seedlings' own BA (breaks the BARE-GROUND case).
    # Recompute NOW over the overstory (1:nstart, no seedlings yet): cs_estab overstory 134.1/70.2; bare stand 0/0.
    ov_ba_pre = stand_ba(s)
    ov_avh_pre = stand_top_height(s)
    # estab.f:175-205 — pre-replicate :estab RNG setup, consumed BEFORE any per-replicate height draw.
    # On the FIRST tally (NTALLY==1) FVS draws once to derive ESDRAW = INT(DRAW·1e5+0.5) and SAVEs it;
    # every tally then reseeds the establishment stream with ESRNSD(.TRUE.,ESDRAW) (odd-forced) and
    # consumes IDUP·NPTIDS draws filling the WK6 site-prep vector. Both advance the :estab stream ahead
    # of the height draws, so jl MUST consume them or every replicate's BACHLO height is off — which
    # (D10) shifts the sp3 seedling sizes, hence the cycle each crosses 3" DBH into the large-tree DGF,
    # which desyncs the sp13 DGSCOR serial-correlation stream and spreads the sawtimber tail.
    # IE/EM (estb/estab.f): the PLANT/NATURAL trees are booked INSIDE the tally's per-plot loop (DO 322, estab.f:975),
    # each plot NCOUNT drawing its planted heights from the post-ESAVE state of its OWN body (:967) before the
    # ESRNSD(ESAVE) reseed (:1075). ie_autoes_establish! (same ESNUTR seam, just before) computed those states from
    # the tally's seed chain; use them verbatim and leave ESS0 at the post-tally ESAVE it set. No ESTAB call this
    # cycle (no states) ⇒ the replicate chain below.
    es_ps = s.estab.es_plot_state
    use_ps = (s.variant isa InlandEmpire || s.variant isa EasternMontana) &&
             s.estab.es_plot_year == yr && length(es_ps) == nptids * idup
    es0_post_tally = s.rng.es0
    pl_plot = Int32[]                  # plot NCOUNT of each record this pass creates (use_ps mode)
    # The tally already ran ESTAB's preamble (NTALLY bookkeeping, ESDRAW, WK6 fill) for this same ESTAB call; the
    # keyword trees are part of that call, not a second one ⇒ skip it (an extra NTALLY increment would e.g. turn a
    # LONE tally's NTALLY=0 into 1 and fire a spurious 20-yr continuation next cycle).
    if !use_ps
        if strp
            _strp_estab_rng!(s, strp_nt > 0 ? strp_nt : 1, nptids, idup, yr, per)
            s.estab.ntally = Int32(strp_nt > 0 ? strp_nt : 1)
        else
            s.estab.ntally += Int32(1)
            if s.estab.ntally == Int32(1)
                s.estab.es_seed = floor(esrann!(s.rng) * 100000f0 + 0.5f0)   # fresh ESDRAW (NTALLY==1)
            end
            esd = s.estab.es_seed
            (esd % 2f0 == 0f0) && (esd += 1f0)                               # ESRNSD odd-force (esrann.f:56)
            s.rng.es0 = Float64(esd)
            for _ in 1:(nptids * idup)                                       # WK6 site-prep fill (estab.f:202-205)
                esrann!(s.rng)
            end
        end
    end
    _strp_ipprep = (s.variant isa Teton && length(s.estab.es_ipprep) == nptids * idup) ? s.estab.es_ipprep : Int32[]
    # estab.f outer loop: `for nn in 1:NPTIDS` (each inventory point) × `idup` replicates
    # → NPTIDS·idup records total. For a BARE stand every point is identical (BAAA=0,
    # uniform slope/aspect/habitat), so the per-point variables don't vary; only the
    # record count and the ESRANN draw count scale with NPTIDS. (ptree already divides by
    # dupnpt = NPTIDS·idup, so the planted TPA is conserved across all the records.)
    _ie_first2 = Dict{Int,Float32}()   # IE PLANT-height DILATE = FIRST(2,sp) accumulator (estab.f:181 init 0.1, :1039 sqrt); IE branch only
    # strp estab.f (TT/UT): GENTIM=FINT−5 once before the plot loop (:317); each PLANT resets it to FINT−DELAY−5
    # (:508-512) AFTER its ESSUBH call, so the next PLANT's ESSUBH AGE (essubh.f:70) reads the previous one's value.
    _gchain = Float32(gentim)
    _tt_first2 = Dict{Int,Float32}()   # TT FIRST(2,sp) DILATE accumulator (estab.f:104 init 0.1, :467/:491 sqrt)
    @inbounds for nn in 1:nptids, rep in 1:idup
        # per-replicate establishment RNG draws (estab.f:216-221): two for emsqr
        # (previously discarded on the no-treeht path), one for esdraw (the re-seed value).
        # estab.f:646-650 EMSQR = ±DRAW2 (sign from DRAW1<0.5). These two draws align with the live oracle only
        # for the FIRST plot; from plot 2 on, live's per-plot ESRANN count is inflated by the AUTOES natural-regen
        # tally (STOADJ block + species tally, estab.f:651+) that jl does not model (that is #143). The IE PLANT
        # branch below now USES these two draws to form EMSQR (RNG-NEUTRAL — the draws were already taken); the CI
        # essubh disp stays 0 (#143-entangled, .sum-inert on cit01). The essubh MEAN (PN) is bit-exact vs live.
        _emd1 = esrann!(s.rng); _emd2 = esrann!(s.rng)
        _emsqr = (_emd1 < 0.5f0 ? -1f0 : 1f0) * _emd2          # estab.f:646-650 (consumed here whether or not used)
        esdraw = floor(esrann!(s.rng) * 100000f0 + 0.5f0)
        if use_ps                                              # plot NCOUNT = (nn-1)·IDUP + rep (estab.f NCOUNT order)
            ncount = (nn - 1) * idup + rep
            s.rng.es0 = es_ps[ncount]; _emsqr = s.estab.es_plot_emsqr[ncount]
        end
        _kph = 0                                               # PLANT index within this plot (DO 322 order)
        for a in due
            sp = round(Int, a.params[1]); (1 <= sp <= MAXSP) || continue
            ptree = a.params[2] * (a.params[3] / 100f0) / dupnpt
            ptree <= 0f0 && continue
            _kph += 1
            # a cycle-number date (<1000) resolves to the calendar year at that cycle (cycle_year_at) before the
            # DELAY offset — else `delay = 2 - 2016 = -2014` ⇒ age≈2019 ⇒ grossly over-sized "seedlings". A
            # calendar-year date carries its own sub-cycle offset unchanged. The date is FVS's 1-BASED cycle number
            # (the `due` filter matches a.year == control.cycle+1) while cycle_year_at takes the 0-based cycle ⇒
            # a.year−1. Passing a.year gave the NEXT cycle's start year (DELAY=FINT ⇒ AGE clamped to 1 ⇒ ~1-ft
            # seedlings); live FVSbm_g16 gives `PLANT 2.0` output identical to the calendar `PLANT <IY(2)>`.
            pyr    = (0 < Int(a.year) < 1000) ? Int(cycle_year_at(s.control, Int(a.year) - 1)) : Int(a.year)
            delay  = pyr - Int(yr)
            trage  = a.params[4] < 0.5f0 ? 2f0 : a.params[4]; trage > 10f0 && (trage = 10f0)
            age = Float32(per) - Float32(delay) - Float32(gentim) + trage; age < 1f0 && (age = 1f0)
            si  = s.plot.sp_site_index[sp]
            # ESSUBH base height (essubh.f:73-82). NE uses its OWN formula — NOT the site-curve height at the
            # tree age: a per-species reference age CARAGE (essubh.f MAPNE, distinct from the htcalc curve map),
            # H = NC-128 site-curve height at CARAGE, then HHT = (H/CARAGE)·min(5, TIME−DELAY) (avg juvenile rate
            # × available time). The `age` above is FVS's REGENT-start AGE (essubh.f:93), used by growth, not the
            # planted height. SN keeps the Curtis-Arney htcalc_height(age).
            _dil_last = get(_ie_first2, sp, 0.1f0)                 # FIRST(2,sp) this tree reads (IE branch advances it)
            hht = if s.variant isa Northeast
                carage = Float32(_NE_ESSUBH_REFAGE[sp])
                (ne_htcalc_height(sp, si, carage) / carage) * min(5f0, Float32(per) - Float32(delay))
            elseif s.variant isa CentralStates
                # CS ESSUBH (cs/essubh.f:72-81): identical NE-style base height with the CS refage map +
                # the CS NC-128 forward curve — H at CARAGE, then (H/CARAGE)·min(5, TIME−DELAY).
                carage = Float32(_CS_ESSUBH_REFAGE[sp])
                (cs_htcalc_height(sp, si, carage) / carage) * min(5f0, Float32(per) - Float32(delay))
            elseif s.variant isa LakeStates
                # LS ESSUBH (ls/essubh.f:68-77): identical NE/CS-style base height — carage from ls/essubh.f's
                # own MAPLS map, H via the LS NC-128 curve (htcalc IVAR=1), then (H/CARAGE)·min(5, TIME−DELAY).
                carage = Float32(_LS_ESSUBH_REFAGE[sp])
                (ls_htcalc_height(sp, si, carage) / carage) * min(5f0, Float32(per) - Float32(delay))
            elseif s.variant isa CentralRockies
                _CR_ESSUBH_HHT[sp]        # cr/essubh.f: a FIXED per-species base height (not a height-at-age curve)
            elseif s.variant isa InlandEmpire
                # IE NATURAL/PLANT base height (ie/essubh.f) — the subsequent/planted-tree height model
                # HHT = EXP(PN + EMSQR·DILATE·BNORM·SIG). IHTSER from the shared estab MYGRUP→MYHTS bracket
                # (em_ihtser); IPREP=1 (NONE) / IPHY=3 defaults (esplt2.f:191-192); BAA = overstory competition
                # clamp[1,400]; XCOS/XSIN = cos/sin(aspect)·slope. disp = EMSQR·DILATE·BNORM (estab.f:646-650 +
                # essubh.f:118), the log-normal dispersion realized per PLANT record — same form the IE tally path
                # uses (inlandempire/establishment.jl:1002-1005). Restores the planted-height SPREAD: the median
                # EXP(PN) (disp=0) was biased LOW (Jensen: E[exp(disp·σ)]>1, plus the XMIN floor lifts the low
                # tail) ⇒ the cohort out-grows/crosses breast height too slowly ⇒ the BA/CCF/SDI under-bias.
                # LIMITATION (distribution-level, NOT per-record bit-exact): FVS draws EMSQR once per plot and
                # advances FIRST(2,sp) in the tally block BEFORE the PLANT block WITHIN the same plot; jl runs the
                # whole tally then the whole PLANT block, so neither the per-plot EMSQR nor the FIRST(2,sp)
                # accumulator is threaded between the two — this reproduces the aggregate (seed-invariant) spread,
                # not the exact per-record heights. Used ONLY when the keyword gives no AVE.HEIGHT (the treeht≥0.1
                # branch below overrides it with the user height + lognormal BACHLO draw). XMIN floor applied below.
                let _slo = s.plot.slope
                    _iage = trunc(Int, (Float32(per) - Float32(delay) - Float32(gentim)) + 0.5f0)   # essubh.f AGE=TIME-DELAY-GENTIM; IAGE=INT(AGE+.5)
                    _iage < 1 && (_iage = 1); _iage > 20 && (_iage = 20)
                    _dil = get(_ie_first2, sp, 0.1f0); _ie_first2[sp] = sqrt(_dil)                   # FIRST(2,sp): 0.1 → sqrt per use (estab.f:1039)
                    ie_essubh(sp, age, clamp(s.plot.basal_area, 1f0, 400f0),
                              em_ihtser(Int(s.plot.habitat_code)), 1, 3,
                              _slo*cos(s.plot.aspect), _slo*sin(s.plot.aspect), _slo, s.plot.elevation,
                              _emsqr * _dil * _IE_ES_BNORML[_iage])
                end
            elseif s.variant isa Teton
                # tt/essubh.f: a fixed per-species base height, except PP (10) which takes the CI PP subsequent-height
                # model (essubh.f:123-125) HHT=EXP(PN+EMSQR·DILATE·BNORM·0.49076), PN=−1.99480+1.53946·ln(AGE)
                # −0.00402·BAA −0.14710 +UPRE(IPREP) −0.01155·ELEV; AGE=TIME−DELAY−GENTIM+TRAGE (TIME=FINT, the
                # rounded/clamped DELAY, the GENTIM chain), IAGE=INT(AGE+.5). DILATE=FIRST(2,sp) (estab.f:467), then
                # FIRST(2,sp)=SQRT(DILATE) (:491) for every PLANT species. IPREP=1 (UPRE=0). Clamped [XMIN,HHTMAX].
                _dil_tt = get(_tt_first2, sp, 0.1f0); _tt_first2[sp] = sqrt(_dil_tt)
                if sp == 10
                    _pdt = Float32(clamp(delay, -3, per))
                    _aget = Float32(per) - _pdt - _gchain + trage; _aget < 1f0 && (_aget = 1f0)
                    _iaget = clamp(trunc(Int, _aget + 0.5f0), 1, length(_IE_ES_BNORML))
                    _ipr = isempty(_strp_ipprep) ? 1 : Int(_strp_ipprep[(nn - 1) * idup + rep])
                    _upre = _ipr == 2 ? 0.20729f0 : _ipr == 3 ? 0.18491f0 : _ipr == 4 ? 0.11864f0 : 0f0
                    _pnt = -1.99480f0 + 1.53946f0 * flog(_aget) - 0.00402f0 * clamp(s.plot.basal_area, 1f0, 400f0) -
                           0.14710f0 + _upre - 0.01155f0 * s.plot.elevation
                    fexp(_pnt + _emsqr * _dil_tt * _IE_ES_BNORML[_iaget] * 0.49076f0)
                else
                    _TT_ESSUBH_HHT[sp]
                end
            elseif s.variant isa CentralIdaho
                # CI subsequent/planted base height (ci/essubh.f, HHT=EXP(PN + disp·SIG)). IHTSER from the
                # shared estb habitat-bracket chain (em_ihtser == the shared estab MYGRUP→MYHTS map, estab.f:493);
                # IPREP=1/IPHY=3 defaults (esplt2.f:191-192). BAA=overstory competition clamp[1,400]; XCOS/XSIN=
                # cos/sin(aspect)·slope (estab.f:480). disp = EMSQR·DILATE·BNORM (per-stand 2-draw EMSQR × per-species
                # sqrt-shrink DILATE × deterministic BNORML[IAGE]); MEASUREMENT PASS uses disp=0 (deterministic mean),
                # validated .sum-inert on cit01 (planted seedlings stay sub-threshold, never enter the summary TPA).
                let _slo = s.plot.slope
                    ci_essubh(sp, age, clamp(s.plot.basal_area, 1f0, 400f0),
                              em_ihtser(Int(s.plot.habitat_code)), 1, 3,
                              _slo*cos(s.plot.aspect), _slo*sin(s.plot.aspect), _slo, s.plot.elevation, 0f0)
                end
            elseif s.variant isa BlueMountains
                # BM base height (bm/essubh.f): HHT = SMHTGF(sp, MODE=0, DTIME=AGE) = the small-tree height-at-total-
                # age curve. Deterministic — NO EMSQR/DILATE/ELEV (bm/essubh.f discards them). SI = per-species SITEAR.
                bm_essubh_hht(sp, si, age)
            elseif s.variant isa Utah
                # UT base height (ut/essubh.f): a FIXED per-species table (no age/site/EMSQR), clamped [XMIN,HHTMAX]
                # by the shared engine. (The UT stump-sprout subsystem — the strp/esuckr.f + vstrp/essprt.f
                # CASE('UT') Crouch aspen model — is ported and dispatched in sprout.jl; utt01 runs end-to-end
                # with no crash and is bit-exact-or-cornered vs FVSut_g16, the residual being the pre-existing
                # UT growth straddle #206, which is present even in the unthinned control stand.)
                _UT_ESSUBH_HHT[sp]
            elseif s.variant isa EasternMontana
                # EM subsequent/planted base height (em/essubh.f, deterministic EXP(PN)). IHTSER from the habitat
                # code bracket search; IPHY=3 / IPREP=1 defaults (esplt2.f). BAA=overstory competition BA clamp[1,400].
                _slo = s.plot.slope
                em_essubh_hht(sp, flog(age), clamp(s.plot.basal_area, 1f0, 400f0),
                              _slo*cos(s.plot.aspect), _slo*sin(s.plot.aspect), _slo, s.plot.elevation,
                              em_ihtser(Int(s.plot.habitat_code)), 3, 1)
            elseif s.variant isa Klamath
                _NC_ESSUBH_HHT[sp]        # nc/essubh.f fixed per-species base height; clamped [XMIN,HHTMAX]
            elseif s.variant isa EastCascades
                # EC base height (ec/essubh.f → ec/smhtgf.f MODE=0): the small-tree height-at-total-age
                # curve HHT = SMHTGF(sp, AGE) with SI = the species' SITEAR. Deterministic (no EMSQR/DILATE/
                # ELEV). The PLANT-no-height branch below then adds the [0,1.5] RAN draw (ec/estab.f:485).
                ec_essubh_hht(sp, si, age)
            elseif s.variant isa WestCascades || s.variant isa PacificNorthwest
                # WC/PN base height (wc/essubh.f == pn/essubh.f == op/essubh.f): a 1.0-ft, 0.1" seedling grown 5 yr by SMHGDG MODE=0
                # (CR=0.5, PTBAL=PTBA=0, AVHT=(5/FINT)·AVH+((FINT−5)/FINT)·ATAVH, smhgdg.f:182-190); HHT=1+HG5,
                # ×FINT/5 when FINT<5 (essubh.f:75-91). Redwood (sp 17) is HHT=2.0 flat. Deterministic.
                if sp == 17
                    2f0
                else
                    avht = (5f0 / Float32(per)) * s.plot.avg_height +
                           ((Float32(per) - 5f0) / Float32(per)) * s.plot.at_avg_ht
                    hg5, _ = _rg_smhgdg(s.variant, sp, 1f0, 0.1f0, 0.5f0, 0f0, 0f0, si, avht)   # PN: no DF King SI
                    hht_wc = 1f0 + hg5
                    per < 5 && (hht_wc *= Float32(per) / 5f0)
                    hht_wc
                end
            elseif s.variant isa Olympic
                # OP base height (op/essubh.f): a 1.0-ft seedling grown 5 yr by op/smhgdg.f (Gould–Harrington
                # small-tree HG5). ESSUBH-mode (MODE=0) forces CR=0.5, PTBA=PTBAL=0, RELHT=H/AVHT with
                # AVHT=(5/FINT)·AVH+((FINT−5)/FINT)·ATAVH (smhgdg.f:234-242). HHT=1+HG5, scaled ·FINT/5 when FINT<5
                # (essubh.f:75-91). Redwood (sp 17) is assigned HHT=2.0 flat (no SMHGDG). Deterministic — no draw.
                if sp == 17
                    2f0
                else
                    avht = (5f0 / Float32(per)) * s.plot.avg_height +
                           ((Float32(per) - 5f0) / Float32(per)) * s.plot.at_avg_ht
                    hg5, _ = op_smhgdg(sp, 1f0, 0.1f0, 0.5f0, 0f0, 0f0, avht, si)
                    hht_op = 1f0 + hg5
                    per < 5 && (hht_op *= Float32(per) / 5f0)
                    hht_op
                end
            else
                htcalc_height(bc, sp, si, age, montane)
            end
            treeht = a.params[5]
            # HTADJ (esin.f opt 15 → esnutr.f 442): per-species height adjustment added to HHT BEFORE the
            # XMIN/0.05 floor and HHTMAX clamp (estab.f:932/1033/1036). Default 0 (empty dict) ⇒ inert.
            hadj = isempty(s.estab.ht_adj) ? 0f0 : get(s.estab.ht_adj, Int32(sp), 0f0)
            _nph = Int(s.estab.es_plot_nph)
            if use_ps                                               # IE/EM: the tally's exact DO 322 height
                _ncnt = (nn - 1) * idup + rep
                _dl = (_nph > 0 && length(s.estab.es_plot_dil) == nptids * idup * _nph && _kph <= _nph) ?
                      s.estab.es_plot_dil[(_ncnt - 1) * _nph + _kph] :
                      _dil_last          # no tally dilations (no-stocking branch): planted-only FIRST(2) chain
                hht = estb_planted_height(s, a, per, Int(yr), _emsqr, _dl)
            elseif treeht >= 0.1f0                                  # PLANT specified a height
                hht = treeht; xh = flog(hht)
                while true
                    xxh = fexp(bachlo(s.rng, xh, 0.5f0; stream = :estab))
                    (0.5f0 * hht <= xxh <= 2f0 * hht) && (hht = xxh; break)
                end
                hht += hadj                                        # estab.f:1033 HHT=HHT+HTADJ (before the 0.05 floor)
                hht < 0.05f0 && (hht = 0.05f0)                      # PLANT floor 0.05 (estab.f:1034)
            elseif s.variant isa EasternMontana || s.variant isa CentralIdaho ||
                   s.variant isa Klamath ||
                   s.variant isa InlandEmpire
                # (BM is NOT in this group: FVSbm is built from strp/estab.f, whose no-user-height PLANT path
                # (estab.f:485-489) DOES draw RAN=BACHLO(0.5,0.25) in [0,1.5] and adds it — live FVSbm_g16
                # debug: ESSUBH 7.805 → HHT 8.41 for WL. BM takes the default RAN branch below.)
                # Shared estb/estab.f:1035-1037 PLANT (no user height): HHT = essubh + HTADJ(default 0), floor XMIN —
                # NO RAN draw. Only the user-specified-height branch (treeht≥0.1, estab.f:1026-1034) draws the lognormal
                # BACHLO perturbation. jl already consumes the per-replicate EMSQR/ESDRAW draws (line ~218) for stream
                # sync, so skipping this extra RAN keeps the :estab stream aligned vs the live oracle. EM validated;
                # CI added #154 (was wrongly taking the else RAN-branch below → +~0.5 ft spurious height + a stream
                # desync). NOTE (CR/IE/TT): same shared-source no-draw applies, latent behind their essubh branches
                # (their validation used no-PLANT DB stands); fold them in when a PLANT .key is validated per variant.
                hht += hadj                                        # estab.f:1036 HHT=HHT+HTADJ (before the XMIN floor)
                hht < es_xmin[sp] && (hht = es_xmin[sp])
            else                                                   # default: RAN~N(0.5,0.25), accept RAN∈[ran_lo,ran_hi]
                while true
                    ran = bachlo(s.rng, 0.5f0, 0.25f0; stream = :estab)
                    (ran_lo <= ran <= ran_hi) && (hht += ran; break)  # estab.f:483/490 (variant-specific window)
                end
                hht += hadj                                        # estab.f:932 HEIGHT(N)=HHT+HTADJ (before the XMIN floor)
                hht < es_xmin[sp] && (hht = es_xmin[sp])           # default/natural floor XMIN (estab.f:1037)
            end
            hht > es_hhtmax[sp] && (hht = es_hhtmax[sp])
            if s.variant isa Teton || s.variant isa Utah               # strp estab.f:508-512 GENTIM reset (chain)
                _pdg = Float32(clamp(delay, -3, per))
                _gchain = (Float32(per) - _pdg) < 5f0 ? 0f0 : Float32(per) - _pdg - 5f0
            end
            ibrkup = floor(Int, ptree / 10f0 + 1f0); brk = Float32(ibrkup)
            # Establishment DBH from the grown seedling height (esgent.f:55-62). A seedling still BELOW
            # breast height (HT < 4.5 ft) has no real DBH — FVS assigns the nominal `DBH = 0.1 + 0.001·HT`
            # (esgent.f:56), NOT the HTDBH⁻¹ inverse. jl previously ran HTDBH⁻¹ for every seedling, which
            # over-sized sub-breast-height regen (bare_natural: DBH 0.225 vs live 0.10 at HT~3.4 ft),
            # inflating stand BA ~0.26% and biasing large-tree DGF growth (D10). Only HT ≥ 4.5 uses the
            # inverse, floored to the species min DIAM + the height-proportional add.
            if s.variant isa BlueMountains || s.variant isa Teton || s.variant isa Utah
                # strp/estab.f:626 DBH(ITRN)=0.1 for every new record regardless of height; REGENT(LESTB) (bm_esgent!)
                # then assigns the dubbed DK / D+0.001·HK. The HTDBH inverse here gave 1.3"/2.7" planted WL/PP at
                # birth, which fed the wrong D into the birth-cycle REGENT. TT/UT build the same strp estab.f
                # (tt|ut/estab.f:626); a ≥4.5-ft TT/UT seedling fell to the shared HTDBH inverse, which has no TT/UT
                # :htdbh_p2 coefficients (KeyError crash on any tall PLANT, e.g. TT NC HHT 10+RAN).
                dbh = 0.1f0
            elseif hht < 4.5f0
                dbh = 0.1f0 + 0.001f0 * hht
            elseif s.variant isa EastCascades || s.variant isa Olympic || s.variant isa WestCascades ||
                   s.variant isa PacificNorthwest   # pn/estab.f:626 == wc's
                # ec/estab.f:626 and op/estab.f:626 (wc/estab.f:626 is the same source) both assign the establishment DBH = 0.1 flat; their
                # esgent.f only recomputes DBH when WK4<1 (a partial birth cycle). A full-birth-cycle
                # PLANT/NATURAL tree (WK4=1) keeps DBH=0.1 even after its height exceeds breast height —
                # height grows, DBH does not. (EC has its own ec_htdbh_dbh; OP is ORGANON-volume — neither
                # has the shared :htdbh_* coef arrays, so the shared inverse both mis-modeled them and crashed.)
                dbh = 0.1f0
            else
                dbh = _htdbh_dbh(sd, sp, hht, ifor; isne = s.variant isa Northeast); dbh < 0.1f0 && (dbh = 0.1f0)
                dbh += 0.001f0 * hht
            end
            for _ in 1:ibrkup
                n = t.n + 1; n + Int(t.ndead) > length(t.dbh) && break   # leave room for the dead block (t.n+1…t.n+ndead); else the volume loop `1:(t.n+ndead)` overruns the MAXTRE arrays (intermittent SIGSEGV on dense ESTAB stands with inventory dead records)
                t.n = n
                use_ps && push!(pl_plot, Int32((nn - 1) * idup + rep))
                t.iestat[n]      = Int32(0)  # estab.f:1438 PLANT/NATURAL records: IESTAT=0 (slot may be reused)
                t.zrand[n]       = -999f0    # estab.f:1424 ZRAND(ITRN)=-999.
                if s.variant isa InlandEmpire || s.variant isa EasternMontana
                    # estb/estab.f:1427-1439: DG=HTG=0, OLDPCT=OLDRN=0, WK1=WK2=0, MISPUTZ(ITRN,0) — clear a reused slot.
                    t.diam_growth[n] = 0f0; t.ht_growth[n] = 0f0
                    t.old_crown_pct[n] = 0f0; t.old_random[n] = 0f0
                    t.dg_prev[n] = 0f0; t.mort_pa[n] = 0f0; t.dmr[n] = Int32(0)
                end
                t.tree_id[n]     = Int32(10000000 + (Int(s.control.cycle) + 1) * 10000 + n)   # IDTREE=IDCMP1+ICYC*10000+ITRN (estab.f:164-165,1440) ⇒ TreeList "ES" id
                # IMC (TreeVal): estb/estab.f:1385-1386 — 1, but 2 for a planted tree NOT ranked best (NOTE≠1) while
                # STOADJ>0; strp/estab.f:600 always 1. NOTE comes from the tally's NBEST pass (es_plot_note, plot-major).
                _note = (use_ps && _nph > 0 && length(s.estab.es_plot_note) == nptids * idup * _nph && _kph <= _nph) ?
                        s.estab.es_plot_note[((nn - 1) * idup + rep - 1) * _nph + _kph] : 1
                t.mort_code[n]   = (_note != 1 && s.estab.stoadj > 0f0) ? Int32(2) : Int32(1)
                t.species[n]     = Int32(sp)
                t.dbh[n]         = dbh
                t.height[n]      = hht
                t.tpa[n]         = ptree / brk
                # WK4 = HTIMLT (estb/estab.f:1436). IE/EM (estb): DO 322 calls ESSUBH (TIME=FINT) BEFORE the HTIMLT
                # block, and ie|em/essubh.f MODIFIES its arguments: DELAY → NINT, floored at −3, capped at TIME
                # (:49-53), and TRAGE → TIME−DELAY (:61). Then estab.f:1055-1063 GENTIM = FINT−DELAY<5 ? 0 :
                # FINT−DELAY−5, HTIMLT = min(TRAGE,GENTIM)/(GENTIM+1e-4) — so a standard PLANT gets 5/5.0001 =
                # 0.99998 (MEASURED live FVSie esgent WK4), NOT 1.0: <1 sends a sub-breast-height seedling down
                # esgent.f:60-62 (DBH=0.1+0.001·HT, DG=0) instead of REGENT's 0.1+DIAM·0.01+0.001·HK.
                # Other variants keep the full birth-cycle HTG (1.0; guards slot reuse). #193
                # TT/UT (strp estab.f:506-516, same shape after tt|ut/essubh.f:64-69 rounds/clamps DELAY and sets
                # TRAGE=TIME−DELAY): HTIMLT = min(TRAGE,GENTIM)/(GENTIM+1e-4) = 0.99998 for a start-of-cycle PLANT.
                t.htimlt[n]      = if s.variant isa InlandEmpire || s.variant isa EasternMontana ||
                                      s.variant isa Teton || s.variant isa Utah ||
                                      s.variant isa WestCascades || s.variant isa PacificNorthwest ||   # wc/pn estab.f:508-516
                                      s.variant isa Olympic   # op/estab.f == wc's (0.99998 for a PLANT); inert until OP ESGENT exists
                    _pd = Float32(clamp(delay, -3, per))
                    _pgen = (Float32(per) - _pd) < 5f0 ? 0f0 : Float32(per) - _pd - 5f0
                    min(Float32(per) - _pd, _pgen) / (_pgen + 0.0001f0)
                else
                    1.0f0
                end
                # ABIRTH = AGEPL + GENTIM (estab.f:628/707) — the REGENT-start `age` already computed above IS
                # FVS's tree age (essubh.f:93). jl left birth_age=0 ⇒ established trees ran ~AGEPL+GENTIM (=7 for a
                # default PLANT) years too YOUNG ⇒ htgf's even-aged site curve (steeper when young) over-predicted
                # height growth as planted stands approached the site asymptote (late-cycle TopHt jl-high). CR-gated:
                # the eastern variants share this latent gap but are separately validated (avoid unvalidated churn).
                # CR/TT: western even-aged htgf reads birth_age. IE joins for CLIMATE ONLY: IE growth does NOT
                # read birth_age (only the Climate-FVS BIRTHYR does — climate_growth_wk4! + inlandempire/regent.jl
                # clim_treemult), so setting it is byte-identical for climate-off IE but gives established regen the
                # correct BIRTHYR=THISYR-ABIRTH so the Leites XDF/XPP/XWL transfer distance is nonzero (was: birth_age
                # =0 ⇒ BIRTHYR=now ⇒ XRELGR≡1 ⇒ under-grown diameter/volume under CLIMATE — matches oracle ABIRTH 5-8).
                # EM builds the identical estab.f (FVSem_buildDir/estab.f == FVSie's): ABIRTH=AGADSB/AGEXC/AGEPL (:1235/
                # :1324/:1414) + GENTIM (:1504). EM reads it in Climate-FVS BIRTHYR and the aspen REGENT (HITE1=f(ABIRTH));
                # birth_age=0 gave EM AutoEstb regen BIRTHYR=THISYR ⇒ DF GrowthMult 1.133 vs live 1.045 (stand 5352355010661).
                # #252: every western variant whose estab.f carries ABIRTH(ITRN)=AGEPL… + ABIRTH(I)=ABIRTH(I)+GENTIM
                # (the strp estab.f of UT/NC/PN/WC/EC/SO/CA/WS/OC/OP, the IE-family one of CI/KT, and BC's) — needed by
                # Climate-FVS BIRTHYR (clgmult/clmorts) now that those variants are climate-wired. BM keeps its own
                # measured AGEPL form below; IE/EM take the per-record AGEPL form below.
                (s.variant isa CentralRockies ||
                 s.variant isa CentralIdaho || s.variant isa Kootenai || s.variant isa Klamath ||
                 s.variant isa PacificNorthwest || s.variant isa WestCascades || s.variant isa EastCascades ||
                 s.variant isa SouthCentralOregon || s.variant isa CentralCalifornia || s.variant isa WestSierra ||
                 s.variant isa OregonCoast || s.variant isa Olympic || s.variant isa BritishColumbia) &&
                    (t.birth_age[n] = age)   # ABIRTH=AGEPL+GENTIM (estab.f:628/707)
                # IE/EM: ABIRTH is AGEPL = FINT−DELAY+TRAGE at creation (estab.f:1064/1414, the ESSUBH-clamped DELAY
                # and the ORIGINAL TRAGE of :990) — what REGENT(LESTB) reads for aspen SITAGE (ie/regent.f:572; em
                # HITE1=f(ABIRTH)) — and esgent_add_gentim! adds the cycle's final GENTIM afterwards (estab.f:1504).
                # jl had stored AGE−GENTIM(FINT−5): IE planted aspen then grew from SITAGE 2 instead of live's 7
                # (DBH 0.7728 vs 0.7812).
                if s.variant isa InlandEmpire || s.variant isa EasternMontana
                    _pdi = Float32(clamp(delay, -3, per))
                    t.birth_age[n] = Float32(per) - _pdi + trage
                    s.estab.gentim_post = (Float32(per) - _pdi) < 5f0 ? 0f0 : Float32(per) - _pdi - 5f0
                    s.estab.gentim_cyc = Int32(s.control.cycle)
                end
                # BM (strp/estab.f:517,628): ABIRTH = AGEPL = FINT−DELAY+TRAGE. Read by the birth-cycle aspen REGENT
                # (bm/regent.f:319 LESTB ⇒ SITAGE=ABIRTH) and Climate-FVS BIRTHYR.
                s.variant isa BlueMountains && (t.birth_age[n] = Float32(per) - Float32(delay) + trage)
                # TT/UT (strp estab.f:517,628): ABIRTH = AGEPL = FINT−DELAY+TRAGE with the essubh-rounded DELAY and the
                # ORIGINAL TRAGE (:438) — read by REGENT(LESTB) as SITAGE for aspen/MM (tt/smhtgf.f:45, tt|ut/regent.f
                # LESTB). esgent_add_gentim! adds the call's final GENTIM after ESGENT (estab.f:707).
                if s.variant isa Teton || s.variant isa Utah
                    _pdi = Float32(clamp(delay, -3, per))
                    t.birth_age[n] = Float32(per) - _pdi + trage
                    s.estab.gentim_post = (Float32(per) - _pdi) < 5f0 ? 0f0 : Float32(per) - _pdi - 5f0
                    s.estab.gentim_cyc = Int32(s.control.cycle)
                end
                # Records go on inventory point `nn` (estab.f:313 ITRE=IPTIDS[nn]).
                # point_ba scales each point's raw BA by PI/GROSPC with PI=NPTIDS, so with
                # the planted TPA spread evenly over NPTIDS points each point_ba comes back
                # to the full stand BA — matching the oracle's pba=ba_v fallback (PTBAA≤0)
                # for fresh establishment, for any NPTIDS (NPTIDS=1 ⇒ this is point 1).
                t.plot_id[n]     = use_iptids ? Int32(iptids[nn]) : Int32(nn)   # IPTIDS[nn] = nn-th STOCKABLE point
                t.crown_pct[n]   = Int32(0)            # crown set in phase 2 (REGENT lestb)
                t.crown_ratio[n] = 0f0
                t.norm_ht[n]     = Int32(0)
                # estab.f:604-623 clears these on every new record; the slot may hold a deleted/dead record
                # (a stale ITRUNC>0 with NORMHT=0 sent r4_topkill's cftopk/bftopk to NaN board feet). Same
                # reset the sprout path does (sprout.jl). All are 0 on a never-used slot ⇒ inert there.
                t.trunc[n]       = Int32(0)
                t.defect[n]      = Int32(0); t.special[n] = Int32(0)
                t.cull[n]        = 0f0; t.decay_code[n] = Int32(0); t.woodland_stems[n] = Int32(0)
                t.sort_key[n]    = Float64(n)
                # Newly-established trees carry NO volume in their birth cycle (FVS: VOLS runs BEFORE the new
                # records are inserted — see grow_cycle! note "regen first gets volume from the next cycle's
                # VOLS"). Zero the volume fields so the post-ESTAB .sum does not sum a STALE value inherited from
                # this array slot's prior occupant — a DEAD inventory tree (compute_volumes! writes cuft_vol over
                # 1:(n+ndead), so slots n+1…n+ndead hold real dead-tree volume) or a removed record. Without this,
                # an all-dead-inventory (HISTORY=6) + heavy-ESTAB stand (EM/IE) over-reports TCuFt/MCuFt/BdFt by
                # the reused slot's volume (a 0.1" seedling summing 22 cuft ⇒ 2600 cuft/acre). Correct volume is
                # filled by next cycle's compute_volumes!.
                t.cuft_vol[n]       = 0f0; t.merch_cuft_vol[n] = 0f0
                t.saw_cuft_vol[n]   = 0f0; t.bdft_vol[n]       = 0f0
                t.merch_top_cf[n] = 0f0; t.merch_top_bf[n] = 0f0   # estab.f HT2TD(ITRN,1..2)=0 (a reused slot kept its old top heights)
                created = true
            end
        end
        es = esdraw; (es % 2f0 == 0f0) && (es += 1f0); s.rng.es0 = Float64(es)  # ESRNSD(true,esdraw)
    end
    use_ps && (s.rng.es0 = es0_post_tally)             # ESS0 = the tally's post-ESAVE reseed (estab.f:1075, last plot)
    # FVS books all of a plot's records inside the per-plot loop — best (estab.f:1255) → excess (:1344) → PLANT/
    # NATURAL (:1436) — before the next plot. jl's AUTOES tally booked its naturals (all plots) just before this
    # pass; interleave so plot k's planted records follow plot k's naturals. Storage order drives the REGENT(LESTB)
    # crown-dub draw order (regent.f:301, I=ITRNIN..ITRN) and the SPESRT species-major IND1 used for every later
    # per-record draw, so the block must match the oracle's layout, not just its contents.
    aplot = s.estab.es_aut_plot
    created_pos = collect((nstart + 1):t.n)            # storage positions of the records THIS pass created (tracked
                                                       # through the interleave below — PHASE 2 must crown-dub these)
    if use_ps && !isempty(pl_plot) && !isempty(aplot) && Int(s.estab.es_aut_first) + length(aplot) - 1 == nstart &&
       all(>(0), aplot)
        lo = Int(s.estab.es_aut_first)
        srcs = vcat(collect(lo:nstart), collect((nstart + 1):t.n))
        tags = vcat(aplot, pl_plot)
        cat_ = vcat(zeros(Int, length(aplot)), ones(Int, length(pl_plot)))   # 0 = natural, 1 = planted
        order = sortperm(collect(1:length(srcs)); by = k -> (tags[k], cat_[k], k))
        perm = srcs[order]
        if perm != srcs
            permute_records!(t, lo, perm)
            @inbounds for i in lo:t.n; t.sort_key[i] = Float64(i); end
            # FVS stamps IDTREE=IDCMP1+ICYC*10000+ITRN at the slot where it CREATES each record (estab.f:1260/1350/
            # 1440); after re-laying jl's block into that slot order, re-stamp so the TreeList "ES" ids match.
            icyc_id = (Int(s.control.cycle) + 1) * 10000
            @inbounds for i in lo:t.n; t.tree_id[i] = Int32(10000000 + icyc_id + i); end
            # the created (PLANT/NATURAL keyword) records now sit wherever the interleave put them — not at
            # nstart+1:t.n (those slots may now hold the tally's naturals). Follow them so EM's PHASE-2 crown dub
            # hits the new keyword records, as before the interleave (IE's own esgent crown-dubs every new record).
            created_pos = [lo + k - 1 for k in eachindex(perm) if perm[k] > nstart]
        end
    end
    # PHASE 2 — ESGENT → REGENT(lestb): assign each new tree its open-grown crown in
    # SPESRT (species-then-record) order (regent.f:107-116). cr = 0.89722 −
    # 0.0000461·PCCF + 0.07985·N(0,1)[±1], clamp [0.20,0.90]; the crown draw uses the
    # MAIN RANN stream (separate from the ESRANN heights). The per-cycle CROWN
    # (crown_ratio_update!, run after) then applies its ±1%/yr change limit (~85).
    if created
        newidx = sort(created_pos; by = i -> (Int(t.species[i]), i))
        # NE only: REGENT(LESTB) also GROWS each new seedling its creation cycle (esgent.f:48). SN's
        # essubh assigns the full height-at-age directly, so SN needs no growth here; NE's essubh gives
        # a BASE height that this grows to the cycle-end height (the BARE-stand TopHt fix). XWT=0 for LESTB.
        # CS shares NE's REGENT(LESTB) shape: ESSUBH gives a BASE height that this grows to the cycle-end
        # height, with the CS NC-128 increment + CS balmod (cs/regent.f:118-340, FNT=FINT−5).
        ne_estab = s.variant isa Northeast
        cs_estab = s.variant isa CentralStates
        ls_estab = s.variant isa LakeStates
        local ebau_e, b3_e, avh_e, scale_e, rdiam_e, rnd_e
        local cb1_e, cb2_e, cb3_e, ba_e
        local lcheck_e, lb1_e, lb2_e, lb3_e, lb4_e, lc1_e, lc2_e, lbamax_e, rmsqd_e
        if ne_estab
            ebau_e = ebau_pre                              # PRE-establishment BAL (snapshot above), not the cohort's
            b3_e = sd[:dg_b3]; avh_e = s.plot.avg_height
            # REGENT LESTB period: FNT = FINT−5 (regent.f:118-124; LSKIPH ⇒ no ht growth when FINT≤5).
            scale_e = per > 5 ? Float32(per - 5) / NE_REGENT_REGYR : 0f0   # CON=HGADJ=XRHGRO=1
            rdiam_e = sd[:regent_min_diam]
            rnd_e = s.control.dg_stddev_bound >= 1f0        # DGSD random ±10%
        elseif cs_estab
            cb1_e = sd[:balmod_b1]; cb2_e = sd[:balmod_b2]; cb3_e = sd[:balmod_b3]
            avh_e = ov_avh_pre                             # POST-growth overstory (pre-seedling snapshot above) —
            ba_e = ov_ba_pre                               # NOT the stale plot.* (regent.f end-of-cycle regen timing)
            scale_e = per > 5 ? Float32(per - 5) / 10f0 : 0f0   # FNT/REGYR, REGYR=10 (CON=HGADJ=XRHGRO=1)
            rdiam_e = sd[:htdbh_db]                         # DIAM floor (= cs/regent.f DIAM == htdbh_db)
            rnd_e = s.control.dg_stddev_bound >= 1f0        # DGSD random ±10%
        elseif ls_estab
            # LS shares NE/CS's REGENT(LESTB) shape: ESSUBH gives a BASE height (5-yr) that this grows to the
            # cycle-end height via the LS NC-128 increment (MAPLS) + ls_balmod (ls/regent.f). BARE stand: BA≈0,
            # RMSQD≈0 ⇒ ls_balmod omega=b4/gm≈1 and AVH=0 ⇒ no competition suppression (live FVSls GMOD=1).
            lcheck_e = sd[:balmod_check]; lb1_e = sd[:balmod_b1]; lb2_e = sd[:balmod_b2]; lb3_e = sd[:balmod_b3]
            lb4_e = sd[:balmod_b4]; lc1_e = sd[:balmod_c1]; lc2_e = sd[:balmod_c2]; lbamax_e = sd[:balmod_bamax1]
            avh_e = s.plot.avg_height; ba_e = s.plot.basal_area; rmsqd_e = rmsqd_pre  # pre-establishment RMSQD (DENSE)
            scale_e = per > 5 ? Float32(per - 5) / 10f0 : 0f0   # FNT/REGYR, REGYR=10 (CON=HGADJ=XRHGRO=1)
            rdiam_e = sd[:htdbh_db]
            rnd_e = s.control.dg_stddev_bound >= 1f0
        end
        # IE has its OWN esgent (ie_esgent!, run after this) which does the REGENT(LESTB) crown-dub in the
        # faithful ie/regent.f DO-13 STORAGE order for ALL new records (PLANT + AUTOES). Drawing it here too
        # (in SPESRT species order) would DOUBLE-consume the main stream for IE PLANT/NATURAL records and use
        # the wrong (eastern) order. So skip the crown-dub draw for IE — ie_esgent! is the sole IE crown-dub.
        # BM likewise: bm/regent.f LESTB draws the crown RAN (regent.f:257-264) and the height ZZRAN (:358-360)
        # INTERLEAVED per record on the main stream, so bm_esgent! owns the crown draw too.
        # EM likewise: em/regent.f LESTB draws each new record's crown in STORAGE order (DO 13) before its ZRANDs.
        # TT (tt/regent.f:302-316 DO 13, storage order, before the SMHTGF ZRANDs) and UT (ut/regent.f:230-241, inside
        # the species-major record loop, interleaved with ZZRAN) likewise draw the crown in their own esgent.
        _ie_own_esgent = s.variant isa InlandEmpire || s.variant isa BlueMountains || s.variant isa EasternMontana ||
                         s.variant isa EastCascades || s.variant isa Teton || s.variant isa Utah ||
                         s.variant isa WestCascades ||
                         s.variant isa PacificNorthwest   # WC/PN: regent.f LESTB draws the crown (wc_esgent!)
        @inbounds for i in newidx
            _ie_own_esgent && continue
            ran_cr = 0f0
            while true
                ran_cr = bachlo(s.rng, 0f0, 1f0)
                -1f0 <= ran_cr <= 1f0 && break
            end
            pccf = s.density.point_ccf[Int(t.plot_id[i])]      # PCCF(IPCCF), IPCCF=ITRE(I) (regent.f:160,178)
            cr = clamp(0.89722f0 - 0.0000461f0 * pccf + 0.07985f0 * ran_cr, 0.20f0, 0.90f0)
            icr0 = floor(Int32, cr * 100f0 + 0.5f0)
            t.crown_pct[i]   = icr0
            t.crown_ratio[i] = Float32(icr0)
            if ne_estab                                        # REGENT(LESTB) height growth + new DBH
                sp = Int(t.species[i]); h = t.height[i]; si = s.plot.sp_site_index[sp]
                # XRHGRO = REGHMULT (regent.f HTGR = HTCALC·CON·SCALE·HGADJ·XRHGRO); the LESTB path must apply
                # it too (was hardcoded 1 ⇒ REGHMULT ignored for the establishment cohort, mult_reghmult diverged).
                xrhgro = active_multiplier(s.control, :regh, sp, Int(yr))
                if ne_htcalc_htmax(sp, si) - h <= 1f0
                    htgr = 0.1f0
                else
                    # regent.f:224 HTGR = HTCALC·CON·SCALE·HGADJ·XRHGRO — the LESTB path must apply CON =
                    # exp(htg_cor_small) (= RHCON·exp(HCOR)) too, as NE small_tree_growth.jl:48 does. It was OMITTED
                    # here ⇒ planted seedlings over-grew (WP: CON=0.914, live rawHTGR 7.98 vs jl 8.73; live-stamped).
                    htgr = ne_htcalc_incr(sp, si, ne_htcalc_age(sp, si, h)) *
                           fexp(s.calib.htg_cor_small[sp]) * scale_e * xrhgro
                end
                gmod = ne_balmod(b3_e[sp], ebau_e, t.dbh[i])
                relht = avh_e > 0f0 ? min(h / avh_e, 1f0) : 0f0
                htgr = max(htgr * (1f0 - (1f0 - gmod) * (1f0 - relht)), 0.1f0)
                if rnd_e
                    rh = 0f0
                    while true; rh = bachlo(s.rng, 0f0, 1f0); -1f0 <= rh <= 1f0 && break; end
                    htgr = max(htgr + rh * 0.1f0 * htgr, 0.1f0)
                end
                hk = h + htgr
                # DBH is derived from the UNCAPPED grown height (the HHTMAX clamp below only bounds the
                # REPORTED height, not the diameter — live YB: dbh from the grown ~23.5 ⇒ 1.8, height clamped
                # to HHTMAX 22). Computing dbh from the clamped height under-sized it (SDI/CCF dropped).
                if hk <= 4.5f0                       # regent.f:290-293: DG=0, DBH=D+0.001·HK (no Wykoff inverse)
                    t.dbh[i] = t.dbh[i] + 0.001f0 * hk
                else
                    dnew = _htdbh_dbh(sd, sp, hk, ifor; isne = s.variant isa Northeast); dnew < 0.1f0 && (dnew = 0.1f0)
                    dnew < rdiam_e[sp] && (dnew = rdiam_e[sp])
                    t.dbh[i] = dnew + 0.001f0 * hk
                end
                hk > _NE_ES_HHTMAX[sp] && (hk = _NE_ES_HHTMAX[sp])   # HARD HHTMAX clamp on the REPORTED height
                t.height[i] = hk
            elseif cs_estab                                    # CS REGENT(LESTB) height growth + new DBH
                sp = Int(t.species[i]); h = t.height[i]; si = s.plot.sp_site_index[sp]
                xrhgro = active_multiplier(s.control, :regh, sp, Int(yr))   # REGHMULT (was hardcoded 1)
                if cs_htcalc_htmax(sp, si) - h <= 1f0          # HTMAX−H ≤ 1 ⇒ HTG=0.1 (cs/regent.f:206)
                    htgr = 0.1f0
                else
                    htgr = cs_htcalc_incr(sp, si, cs_htcalc_age(sp, si, h)) * scale_e * xrhgro
                end
                # regent.f:156 BAL=(1-PCT/100)·BA, PCT=BA-percentile. A new seedling (D=0.1) is the smallest ⇒
                # PCT≈0 ⇒ BAL≈full overstory BA (live BAL=134.1=BA). crown_ratio was a wrong proxy (gave ~0.5·BA).
                bal = ba_e
                gmod = cs_balmod(cb1_e[sp], cb2_e[sp], cb3_e[sp], bal, ba_e, t.dbh[i])
                relht = avh_e > 0f0 ? min(h / avh_e, 1f0) : 0f0
                htgr = max(htgr * (1f0 - (1f0 - gmod) * (1f0 - relht)), 0.1f0)
                if rnd_e
                    rh = 0f0
                    while true; rh = bachlo(s.rng, 0f0, 1f0); -1f0 <= rh <= 1f0 && break; end
                    htgr = max(htgr + rh * 0.1f0 * htgr, 0.1f0)
                end
                hk = h + htgr
                # CS LESTB dbh (cs/regent.f:338-341): DBH = htdbh⁻¹(hk), floored to DIAM (or DIAM if hk<4.5),
                # THEN + 0.001·hk. (The htdbh inverse needs hk>4.5; for hk<4.5 the DIAM floor applies first.)
                if hk < 4.5f0
                    dbhk = rdiam_e[sp]
                else
                    dbhk = _htdbh_dbh(sd, sp, hk, ifor; isne = s.variant isa Northeast)
                    dbhk < rdiam_e[sp] && (dbhk = rdiam_e[sp])
                end
                t.dbh[i] = dbhk + 0.001f0 * hk
                hk > _CS_ES_HHTMAX[sp] && (hk = _CS_ES_HHTMAX[sp])   # HARD HHTMAX clamp on the REPORTED height
                t.height[i] = hk
            elseif ls_estab                                    # LS REGENT(LESTB) height growth + new DBH
                sp = Int(t.species[i]); h = t.height[i]; si = s.plot.sp_site_index[sp]
                xrhgro = active_multiplier(s.control, :regh, sp, Int(yr))
                if ls_htcalc_htmax(sp, si) - h <= 1f0
                    htgr = 0.1f0
                else
                    # regent.f:224 CON = exp(htg_cor_small) (= RHCON·exp(HCOR)), as LS small_tree_growth.jl:40 applies.
                    # Was omitted here (inert for JP where CON≈1, but a latent bug for CON≠1 species — cf. the NE fix).
                    htgr = ls_htcalc_incr(sp, si, ls_htcalc_age(sp, si, h)) *
                           fexp(s.calib.htg_cor_small[sp]) * scale_e * xrhgro
                end
                gmod = ls_balmod(sp, t.dbh[i], ba_e, rmsqd_e, lcheck_e, lb1_e, lb2_e, lb3_e, lb4_e, lc1_e, lc2_e, lbamax_e)
                relht = avh_e > 0f0 ? min(h / avh_e, 1f0) : 0f0
                htgr = max(htgr * (1f0 - (1f0 - gmod) * (1f0 - relht)), 0.1f0)
                if rnd_e
                    rh = 0f0
                    while true; rh = bachlo(s.rng, 0f0, 1f0); -1f0 <= rh <= 1f0 && break; end
                    htgr = max(htgr + rh * 0.1f0 * htgr, 0.1f0)
                end
                hk = h + htgr
                if hk < 4.5f0                       # ls/regent.f LESTB dbh: DIAM floor (or htdbh⁻¹), + 0.001·hk
                    dbhk = rdiam_e[sp]
                else
                    dbhk = _htdbh_dbh(sd, sp, hk, ifor; isne = s.variant isa Northeast)
                    dbhk < rdiam_e[sp] && (dbhk = rdiam_e[sp])
                end
                t.dbh[i] = dbhk + 0.001f0 * hk
                hk > _LS_ES_HHTMAX[sp] && (hk = _LS_ES_HHTMAX[sp])   # HARD HHTMAX clamp on the REPORTED height
                t.height[i] = hk
            end
        end
        # ESGENT calls SPESRT to RE-ESTABLISH the species-order sort after adding
        # regen (esgent.f:41-44). SPESRT/LNKCHN visit records in ascending-record
        # order, so reset the lineage key to the physical record position: otherwise
        # stale TRIPLE lineage keys (3·K+offset) from earlier cycles, which are never
        # reconciled without a thinning compaction, scramble the post-establishment
        # species_sort! order and desync the per-tree DGSCOR RNG stream from FVS.
        @inbounds for i in 1:t.n
            t.sort_key[i] = Float64(i)
        end
        compute_density!(s)
    end
    push!(s.estab.years_done, yr)
    return created
end

# =============================================================================
# ESTAB site-prep activity bookkeeping (BURNPREP 491 / MECHPREP 493) — the IACT(,4) status that ECON's
# MECHCST/BURNCST reads (eccalc.f OPSTUS/OPGET3). Pure status bookkeeping: it never changes establishment.
# =============================================================================

"Resolved date of a scheduled activity (a date in 1..MAXCYC is a cycle number ⇒ IY(date), OPEXPN)."
_prep_date(s::StandState, a) = (0 < Int(a.year) < 1000) ? Int(cycle_year_at(s.control, Int(a.year) - 1)) : Int(a.year)

"""
    esetpr_mark!(s, idsdat, kdt)

ESETPR (esetpr.f) status effects: for BURNPREP then MECHPREP, walk the PENDING activities (OPGET2: IACT(,4)=0)
dated in [IDSDAT, KDT] in date order (IOPSRT; ties = input order); the LAST one found is marked done in its own
year (OPDON2(…,IDT,…,I)), every earlier one is deleted (OPDON2(…,−1,…)). The "no-parameter BURNPREP ⇒ METH=3,
RETURN before MECHPREP" branch is inert: ESPRIN always stores the %-of-plots parameter.
"""
function esetpr_mark!(s::StandState, idsdat::Integer, kdt::Integer)
    st = s.estab.prep_status
    for code in (Int32(491), Int32(493))
        pend = Tuple{Int,Int}[]
        for (i, a) in enumerate(s.control.schedule)
            a.icflag == code || continue
            haskey(st, i) && continue
            d = _prep_date(s, a)
            (idsdat <= d <= kdt) && push!(pend, (d, i))
        end
        isempty(pend) && continue
        sort!(pend)
        for (k, (d, i)) in enumerate(pend)
            st[i] = k == length(pend) ? Int32(d) : Int32(-1)
        end
    end
    return
end

"""
    estab_prep_cancel_cycle!(s)

estab.f `IF(NTALLY.GT.1) CALL OPFIND(2,MYACTS(3),NTODO) … OPDEL1`: a continuation / ingrowth tally deletes the
pending BURNPREP/MECHPREP activities due THIS cycle.
"""
function estab_prep_cancel_cycle!(s::StandState)
    y1 = Int(current_cycle_year(s)); y2 = Int(cycle_year_at(s.control, Int(s.control.cycle) + 1))
    st = s.estab.prep_status
    for (i, a) in enumerate(s.control.schedule)
        (a.icflag == Int32(491) || a.icflag == Int32(493)) || continue
        haskey(st, i) && continue
        d = _prep_date(s, a)
        (y1 <= d < y2) && (st[i] = Int32(-1))
    end
    return
end

"""
    estab_prep_tally!(s, ntally, idsdat, kdt)

The site-prep effect of one ESTAB call (estab.f): NTALLY>1 (continuation, or ESTB/AK ingrowth 99) cancels this
cycle's preps; NTALLY=1 runs ESETPR over [IDSDAT, KDT]. Idempotent per cycle (ESNUTR runs once per cycle in FVS).
"""
function estab_prep_tally!(s::StandState, ntally::Integer, idsdat::Integer, kdt::Integer)
    yr = Int32(current_cycle_year(s))
    yr in s.estab.prep_years_done && return
    push!(s.estab.prep_years_done, yr)
    if ntally > 1
        estab_prep_cancel_cycle!(s)
    elseif ntally == 1
        esetpr_mark!(s, idsdat, kdt)
    end
    return
end

"""
    estab_prep_esnutr!(s)

ESNUTR's ESTAB-call decision (strp/ls esnutr.f:145-290) for the variants whose establishment is not jl's IE/EM
AUTOES tally — replicated ONLY to drive the site-prep status (ESETPR / NTALLY>1 cancel). Per cycle, KDT=IY(ICYC+1)−1:
  1. a TALLYONE(428) — else TALLYTWO(429) — due this cycle (the last one, OPGET(NTODO)): NTALLY=IACTK−427,
     IDSDAT=PRMS(1), LONE=.TRUE.; >20 yr stale ⇒ canceled (no ESTAB); a TALLYTWO whose TALLYONE was not done after
     IDSDAT becomes NTALLY=1 (OPSTUS(428,IDSDAT,KDT,0,…));
  2. else a TALLY(427) due this cycle: the latest disturbance date PRMS(1) ⇒ IDSDAT, NTALLY=1 (stale ⇒ canceled);
  3. else the 20-yr continuation: NTALLY>0 ∧ KDT−IDSDAT≤19 ⇒ NTALLY+1;
  4. else a PLANT/NATURAL due this cycle: NTALLY=1 (STRP/LS) or 99 (ESTB CI/KT, AK), IDSDAT=IY(ICYC+1)−20.
Then ESTAB (site-prep effect) and `IF (LONE) NTALLY=0`; the chosen tally is OPDONE'd at KDT, the other tallies due
this cycle are deleted (esnutr.f:300-310). IDSDAT −9999 ⇒ IY(1)−20 (esnutr.f:100).
"""
function estab_prep_esnutr!(s::StandState)
    sched = s.control.schedule
    # TT/UT also take the tally number from here (est.cyc_ntally) to drive the strp ESTAB RNG (establish!).
    (s.variant isa Teton || s.variant isa Utah ||
     any(a -> a.icflag == Int32(491) || a.icflag == Int32(493), sched)) || return
    est = s.estab
    yr = Int32(current_cycle_year(s))
    yr in est.prep_years_done && return
    icyc = Int(s.control.cycle) + 1
    y2 = Int(cycle_year_at(s.control, icyc)); kdt = y2 - 1
    est.prep_idsdat == Int32(-99999) && (est.prep_idsdat = est.idsdat)
    est.prep_idsdat == Int32(-9999) && (est.prep_idsdat = Int32(Int(cycle_year_at(s.control, 0)) - 20))
    ydate(d) = (1 <= d < 1000) ? Int(cycle_year_at(s.control, d - 1)) : d
    due(i, a) = !haskey(est.tally_status, i) && _ec_cycle_of(s, _prep_date(s, a)) == icyc
    tallies(code) = [i for (i, a) in enumerate(sched) if a.icflag == code && due(i, a)]
    ntally = 0; lone = false; chosen = 0
    t1 = tallies(Int32(428)); isempty(t1) && (t1 = tallies(Int32(429)))
    if !isempty(t1)
        chosen = t1[end]; a = sched[chosen]
        ntally = Int(a.icflag) - 427
        est.prep_idsdat = Int32(ydate(round(Int, a.params[1]))); lone = true
        if kdt + 1 - Int(est.prep_idsdat) > 20
            ntally = 0; chosen = 0
        elseif ntally == 2
            done1 = [(Int(sched[i].year), st) for (i, st) in est.tally_status
                     if sched[i].icflag == Int32(428) && Int(est.prep_idsdat) <= _prep_date(s, sched[i]) <= kdt]
            ist = isempty(done1) ? 0 : last(sort(done1))[2]
            (isempty(done1) || ist <= est.prep_idsdat) && (ntally = 1)
        end
    else
        t0 = tallies(Int32(427))
        if !isempty(t0)
            best = -1
            for i in t0
                d = round(Int, sched[i].params[1])
                d > best && (best = d; chosen = i)
            end
            est.prep_idsdat = Int32(ydate(best)); ntally = 1
            kdt + 1 - Int(est.prep_idsdat) > 20 && (ntally = 0; chosen = 0)
        elseif est.prep_ntally > 0 && kdt - Int(est.prep_idsdat) <= 19
            ntally = Int(est.prep_ntally) + 1
        else
            y1 = Int(current_cycle_year(s))
            npnats = count(a -> (a.icflag == Int32(430) || a.icflag == Int32(431)) &&
                                ((y1 <= Int(a.year) < y2) || (0 < Int(a.year) < 1000 && Int(a.year) == icyc)), sched)
            if npnats > 0
                estb_like = s.variant isa CentralIdaho || s.variant isa Kootenai || s.variant isa SoutheastAlaska
                ntally = estb_like ? 99 : 1
                est.prep_idsdat = Int32(y2 - 20)
            end
        end
    end
    est.cyc_ntally = Int32(ntally); est.cyc_ntally_year = yr
    if ntally > 0
        est.prep_ntally = Int32(ntally)
        estab_prep_tally!(s, ntally, Int(est.prep_idsdat), kdt)
        lone && (est.prep_ntally = Int32(0))
    else
        push!(est.prep_years_done, yr)
    end
    # esnutr.f:300-310: the executed tally is done at KDT; every other tally due this cycle is deleted
    for code in (Int32(427), Int32(428), Int32(429)), i in tallies(code)
        est.tally_status[i] = i == chosen ? Int32(kdt) : Int32(-1)
    end
    return
end

"IE/EM with the AUTOES tally switched off: the PLANT/NATURAL catch-all ESTAB call is NTALLY=99 (estb esnutr.f:355)."
function estab_prep_npnats_estb!(s::StandState)
    any(a -> a.icflag == Int32(491) || a.icflag == Int32(493), s.control.schedule) || return
    icyc = Int(s.control.cycle) + 1
    y1 = Int(current_cycle_year(s)); y2 = Int(cycle_year_at(s.control, Int(s.control.cycle) + 1))
    any(a -> (a.icflag == Int32(430) || a.icflag == Int32(431)) &&
             ((y1 <= Int(a.year) < y2) || (0 < Int(a.year) < 1000 && Int(a.year) == icyc)), s.control.schedule) || return
    estab_prep_tally!(s, 99, y2 - 20, y2 - 1)
    return
end
