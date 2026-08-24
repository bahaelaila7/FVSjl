# =============================================================================
# PPE (Parallel Processing Extension) — LANDSCAPE / multi-stand orchestration.
# =============================================================================
# FAITHFUL PORT of the PPMAIN master-cycle control structure (recovered source:
# scratchpad/ppe/recovered/ppbase/ppmain.f + alstd1/alstd2/ppmclk/splaex, from FVS
# rev bc6e2377^; the "no source" premise was wrong — see [[fvsjl-wwpb-ppe-oracle-absent]]).
#
# ---------------------------------------------------------------------------
# What PPMAIN actually does (ppmain.f, distilled):
#   1. INITIALIZE every stand up to the master starting year (INSTND → GRSTND with
#      ICYC=1 → PUTSTD), then ALINIT.
#   2. MASTER-CYCLE LOOP (label 40..210), one iteration per master cycle:
#        a. PPMCLK      — increment the master clock; LEORUN when MICYC > MNCYC.
#        b. ALSTD1      — all-stand start-of-cycle review: GPLOAD (mode-2 point-group
#                         interaction list) + BMSETP (beetle dispersal run lists) +
#                         HVALOC (multistand harvest allocation, only if LHVALL).
#        c. C11SRT      — order the stands for processing (CISNUM → ISNSRT).
#        d. STAND LOOP (label 50): GETSTD each stand, set its processing mode via
#           PPINTR:
#             mode 1 → GRSTND projects it independently to the master-cycle end.
#             mode 2 → GRINCR grows it only up to mortality; BMSDIT seeds the beetle.
#        e. ALSTD2 (label 150, only if any mode-2 stand): BMDRV — the WWPB multi-stand
#           beetle DISPERSAL across the landscape. Then a phase-2 stand loop applies
#           the result per stand: BMKILL → GRADD → GRCEND (finish the cycle's growth).
#        f. End-of-run: emit the landscape reports (ESOUT/MPBOUT/…/GENPRT).
#   3. PPAUSE / stop-restart between master cycles (LPAUSE).
#
# GETSTD/PUTSTD/SPARLS are Fortran-FVS's mechanism for swapping ONE stand's global
# COMMON in/out of a direct-access file (FVS holds a single stand in COMMON at a
# time). In Julia this is STRUCTURALLY UNNECESSARY: each stand is its own state
# object, and the landscape is just a Vector of them — the DA file IS the Vector.
#
# ---------------------------------------------------------------------------
# MODE-1 EQUIVALENCE (the proof this port rests on, since mode-1 needs no oracle):
# In mode 1 the stands DO NOT interact within a master cycle — each is projected
# independently by GRSTND to the master-cycle boundary, and the landscape statistics
# are SPLAEX/ALSTD2's area-weighted aggregate Σ(value·area)/Σ(area) over the stands.
# Therefore the mode-1 landscape result is IDENTICAL to: project each stand fully
# (all its cycles) with the already-oracle-validated per-stand engine `run_keyfile`,
# then area-weight-aggregate the per-report-year rows. The master-cycle stepping and
# the C11SRT processing ORDER are behaviorally INERT for mode-1 output (no coupling ⇒
# order-independent), so we run each stand once and aggregate — bit-identical to the
# stepped form, and every per-stand value is oracle-validated. The C11SRT order is
# still computed (via the bit-exact ppe_index_qsort!) so the port is faithful and the
# kernel participates; it becomes LOAD-BEARING only for mode-2 / harvest scheduling.
#
# MODE-2 (interstand beetle dispersal) + MXHRVP (multistand harvest) require TRUE
# master-cycle stepping — pausing every stand mid-projection to exchange landscape
# state — which is exactly getstd/putstd's job and is NOT exposed by run_keyfile
# (FVSjl projects a whole keyfile per call). Those are documented composition SEAMS
# below (ppe_neighbors + the wwpb outbreak), left as honest stubs rather than faked:
# the beetle biology itself is already bit-exact (src/engine/wwpb*.jl) and runs
# per-stand via the DISPERSE keyword; only the cross-stand spatial redistribution of
# beetle pressure is the missing landscape coupling.
# =============================================================================

"""
One member stand of a PPE landscape: a keyfile to project and its area weight
(acres, PPE's per-stand contribution to TOTALWT). `col`/`row` are optional hex-grid
coordinates for the mode-2 spatial-neighbor machinery (SPNB*/HXINDX); when unset the
stand has no spatial identity and only mode-1 (independent) aggregation applies.
"""
struct PPEStand
    keyfile::String
    area::Float64
    col::Int          # hex-grid column (0 = unplaced)
    row::Int          # hex-grid row    (0 = unplaced)
end
PPEStand(keyfile::AbstractString; area::Real = 1.0, col::Integer = 0, row::Integer = 0) =
    PPEStand(String(keyfile), Float64(area), Int(col), Int(row))

"""
The area-weighted PPE landscape aggregate for one reporting/master-cycle year
(PPEXCM PTSTV1). Fields map to the documented PTSTV1(1..9) slots; `nstands` is how
many member stands reported that year (aligned on the `.sum` year).
"""
struct PPEAggregate
    year::Int
    avbtpa::Float64      # PTSTV1(1)  average before-thin trees/acre
    avbtcuft::Float64    # PTSTV1(2)  average before-thin total cubic volume
    avbmcuft::Float64    # PTSTV1(3)  average before-thin merch cubic volume
    avbbdft::Float64     # PTSTV1(4)  average before-thin board-foot volume
    avbacc::Float64      # PTSTV1(5)  last-cycle accretion (cuft/acre/yr)
    avbmort::Float64     # PTSTV1(6)  last-cycle mortality  (cuft/acre/yr)
    msperiod::Int        # PTSTV1(7)  master-cycle period length (yrs to next report)
    totalwt::Float64     # PTSTV1(8)  Σ area over reporting stands
    oldtarg::Float64     # PTSTV1(9)  last master-cycle harvest target (0 = none set)
    nstands::Int
end

# Parse the per-cycle before-thin fields out of one `.sum` data row (US/imperial
# 28-token layout). tpa=col3, cuft=col9, mcuft=col10, bdft=col12 from the FRONT;
# period-length/accretion/mortality read from the END (end-5 / end-4 / end-3) so they
# are robust to the volume-block width. Returns nothing for header/blank/non-data
# lines. `prd` is IOSUM(14) (PERIOD LENGTH, YEARS) — the CMADDS accretion/mortality
# weight (see the aggregation note below). Returns nothing for header/blank lines.
function _ppe_parse_sum_row(line::AbstractString)
    f = split(strip(line))
    length(f) < 20 && return nothing
    all(c -> isdigit(c) || c == '-', f[1]) || return nothing        # year is an integer
    try
        year = parse(Int, f[1])
        tpa  = parse(Float64, f[3])
        cuft = parse(Float64, f[9])
        mcuft = parse(Float64, f[10])
        bdft = parse(Float64, f[12])
        prd  = parse(Float64, f[end-5])
        acc  = parse(Float64, f[end-4])
        mort = parse(Float64, f[end-3])
        return (year = year, tpa = tpa, cuft = cuft, mcuft = mcuft, bdft = bdft,
                prd = prd, acc = acc, mort = mort)
    catch
        return nothing
    end
end

# The CMADDS/CMPRT2 composite-yield aggregation (ppbase/cmadds.f + cmprt2.f, the
# routine that builds the "COMPOSITE YIELD STATISTICS" table — validated bit-exact
# against the historical FVSppe oracle, scratchpad/ppe/recovered). `stand_rows` is one
# (area, rows) pair per landscape stand, `rows` a vector of `_ppe_parse_sum_row`
# NamedTuples. For each report year, CMADDS accumulates the weighted sums and CMPRT
# averages them (IFIX(x+.5) round-to-nearest, kept as Float here — callers round on
# compare):
#   TREES/CU FT/MERCH (IOSUM 3..6): area-weighted mean  Σ(v·w) / Σ(w)         (÷ PRBSUM)
#   ACCRETION/MORTALITY (IOSUM 15,16): weighted by area·PERIOD, per cmadds.f
#       ALL(15)+=v·(w·prd); IALL(15)=ALL(15)/ALL(14) where ALL(14)=Σ(w·prd)   (÷ Σ w·prd)
#   PERIOD (IOSUM 14): IALL(14)=Σ(prd·w)/Σ(w)                                 (÷ PRBSUM)
#   TOTAL SAMPLE WEIGHT (IOSUM 17): Σ w
# (BA/CCF/DOM HT (IOSUM 11..13) also divide by PRBSUM but use each stand's AFTER-thin
# residual value, and REMOVALS (IOSUM 7..10) divide by the treated weight HRVSUM — both
# outside PTSTV1(1..9) and not produced here; see the CHECK-1 note in the PPE tests.)
function _ppe_aggregate(stand_rows)
    wtpa = Dict{Int,Float64}(); wcuft = Dict{Int,Float64}(); wmcuft = Dict{Int,Float64}()
    wbdft = Dict{Int,Float64}(); wacc = Dict{Int,Float64}(); wmort = Dict{Int,Float64}()
    wsum = Dict{Int,Float64}(); wprd = Dict{Int,Float64}(); ncnt = Dict{Int,Int}()
    for (area, rows) in stand_rows
        for r in rows
            a = area; wp = a * r.prd
            wtpa[r.year]  = get(wtpa, r.year, 0.0)  + r.tpa  * a
            wcuft[r.year] = get(wcuft, r.year, 0.0) + r.cuft * a
            wmcuft[r.year]= get(wmcuft, r.year, 0.0)+ r.mcuft* a
            wbdft[r.year] = get(wbdft, r.year, 0.0) + r.bdft * a
            wacc[r.year]  = get(wacc, r.year, 0.0)  + r.acc  * wp   # ALL(15) += v·(w·prd)
            wmort[r.year] = get(wmort, r.year, 0.0) + r.mort * wp   # ALL(16) += v·(w·prd)
            wprd[r.year]  = get(wprd, r.year, 0.0)  + wp            # ALL(14) = Σ(w·prd)
            wsum[r.year]  = get(wsum, r.year, 0.0)  + a             # PRBSUM
            ncnt[r.year]  = get(ncnt, r.year, 0)    + 1
        end
    end
    out = PPEAggregate[]
    for yr in sort(collect(keys(wsum)))
        w = wsum[yr]
        w <= 0.0 && continue
        wp = wprd[yr]
        # ACC/MOR divide by Σ(w·prd); zero when the growth period is zero (final row).
        acc  = wp > 0.0 ? wacc[yr]  / wp : 0.0
        mort = wp > 0.0 ? wmort[yr] / wp : 0.0
        # PTSTV1(7) MSPERIOD == IALL(14) = Σ(prd·w)/Σw (the master-cycle period length).
        msperiod = floor(Int, wp / w + 0.5)
        push!(out, PPEAggregate(yr, wtpa[yr]/w, wcuft[yr]/w, wmcuft[yr]/w, wbdft[yr]/w,
                                acc, mort, msperiod, w,
                                0.0,               # PTSTV1(9) OLDTARG: no landscape harvest target set
                                ncnt[yr]))
    end
    return out
end

"""
    ppe_processing_order(stands) -> Vector{Int}

PPMAIN's C11SRT master stand ordering (ppmain.f:170 `CALL C11SRT(NOSTND,CISNUM,
ISNSRT,.FALSE.)`): stands are processed each master cycle in ascending order of their
internal stand identifier. Uses the bit-exact `ppe_index_qsort!` (the ported C11SRT)
over each stand's identity key (keyfile basename, the Julia analog of CISNUM). Returns
a permutation of `1:length(stands)`. Behaviorally INERT for mode-1 output (stands do
not interact) but load-bearing for mode-2 / harvest; wired so the port is faithful.
"""
function ppe_processing_order(stands::AbstractVector{PPEStand})
    n = length(stands)
    n == 0 && return Int[]
    keys = [basename(s.keyfile) for s in stands]          # CISNUM analog
    index = collect(1:n)                                   # ISNSRT pre-loaded 1:n (PPMAIN passes LSEQ=.FALSE.)
    ppe_index_qsort!(index, keys; lseq = false)            # ascending character sort of the INDEX
    return index
end

"""
    ppe_neighbors(stand, stands; nstnd = length(stands)) -> Vector{Int}

The SPNB* / HXINDX spatial-neighbor lookup used by the mode-2 interstand disturbance
spread (ALSTD2/BMDRV). Given a placed `stand` (col/row set), returns the landscape
indices of its up-to-6 hex neighbors (N/NE/SE/S/SW/NW) via the bit-exact `ppe_hxindx`.
Stands with no coordinates (col/row 0) have no neighbors. This is the wired seam for
mode-2; the cross-stand beetle-pressure redistribution that consumes it is not yet
composed (see the mode-2 note in the file header).
"""
function ppe_neighbors(stand::PPEStand, stands::AbstractVector{PPEStand};
                       nstnd::Integer = length(stands))
    (stand.col == 0 && stand.row == 0) && return Int[]
    # placed stands carry a 1-based hex index = (row-1)*ncols + col; recover it and
    # map each hex-neighbor pointer back to the stand occupying that hex.
    hexof(s) = (s.col == 0 && s.row == 0) ? 0 :
               begin
                   nr = trunc(Int, sqrt(Float32(nstnd))); nc = nr
                   nc * nr < nstnd && (nr += 1)
                   (s.row - 1) * nc + s.col
               end
    self = hexof(stand)
    self == 0 && return Int[]
    byhex = Dict(hexof(s) => i for (i, s) in enumerate(stands) if hexof(s) != 0)
    out = Int[]
    for nei in 1:6
        p = ppe_hxindx(nei, self, Int(nstnd))
        p != 0 && haskey(byhex, p) && push!(out, byhex[p])
    end
    return out
end

"""
    ppe_landscape_mode2!(members, w, coeffs; iyr1, iyr2, dispersal-params...) -> (kills, ls)

The MODE-2 interstand-beetle-dispersal composition (ppmain.f label-150 ALSTD2 phase +
the phase-2 BMKILL stand loop). `members` is a vector of NamedTuples
`(trees=<FVS treelist>, area=<acres>, xloc=<m>, yloc=<m>, stock=<Bool>,
sp_alpha=<sp->code>)` — one per landscape stand, each already GROWN to mortality
(GRINCR, mode-2). This runs the LANDSCAPE step faithfully:

  1. bmsdit! bins each stand's treelist into its WwpbStand (the FVS→BM bridge),
     seeds rvdsc=1 (drought off) and any inventory beetle damage.
  2. bmdrv_multi! runs the master per-year loop, redistributing BKP ACROSS the
     placed stands (bmatct_multi!, golden-validated) — the piece mode-1 cannot do.
  3. bmkill! hands each stand's beetle mortality back to its treelist (returned as
     `kills[i]`, a per-record WK2 mortality vector; the caller applies it, mirroring
     ppmain.f's `CALL BMKILL` before GRADD/GRCEND).

Returns `(kills, ls)` where `ls::WwpbLandscape` carries the per-stand BKPOUT/BKPIN/
SELFBKP dispersal ledger. Every kernel here is bit-exact (bmsdit!/bmcgrf!/bmcbkp!/
bmcnum!/bmatct_multi!/bmistd!/bmmort!/bmkill!, driver-golden-validated); the
composition is the faithful bmdrv.f/ppmain.f control structure. At one stand /
OUTOFF it collapses to the single-stand DISPERSE path (wwpb_apply!).

RESIDUAL SEAM (measured architectural block): feeding this LIVE, mid-projection —
i.e. pausing every member stand's FVSjl projection AT the outbreak cycle boundary to
exchange landscape BKP, then resuming — is exactly getstd/putstd's master-cycle
stepping, which `run_keyfile` (whole-keyfile-per-call) does not expose. So this
function takes the already-grown per-stand treelists as input rather than driving the
engine; wiring it into a live per-cycle FVSjl landscape stepper is the one remaining
coupling. It is ADDITIVE + INERT (no simulate.jl seam; mode-1 stays byte-identical).
"""
function ppe_landscape_mode2!(members, w, coeffs; iyr1::Integer, iyr2::Integer,
                              ipson::Bool=false, usera::NTuple{3,Float32}=(1.0f0,1.0f0,1.0f0),
                              selfa::NTuple{3,Float32}=(1.0f0,1.0f0,1.0f0),
                              userc::NTuple{3,Float32}=(100.0f0,100.0f0,100.0f0),
                              urmax::NTuple{3,Float32}=(15.0f0,15.0f0,15.0f0),
                              outoff::Bool=true, ufloat::Float32=-1.0f0, sdd::Float32=0.0f0,
                              rvod::Float32=1.0f0, stocko::Float32=1.0f0,
                              seed_pbkill=nothing)
    n = length(members)
    stands = WwpbStand[]
    area  = Float32[]; xloc = Float32[]; yloc = Float32[]; stock = Bool[]
    for (i, m) in enumerate(members)
        st = WwpbStand()
        bmsdit!(st, m.trees, w, m.sp_alpha)
        fill!(st.rvdsc, 1.0f0)                    # drought model not run ⇒ neutral
        if seed_pbkill !== nothing
            sp = seed_pbkill[i]
            sp !== nothing && (st.pbkill .= sp)
        end
        push!(stands, st)
        push!(area, Float32(m.area)); push!(xloc, Float32(m.xloc)); push!(yloc, Float32(m.yloc))
        push!(stock, get(m, :stock, true))
    end
    ls = WwpbLandscape(xloc, yloc, area, stock)
    bmdrv_multi!(ls, stands, w, coeffs; area=area, iyr1=Int(iyr1), iyr2=Int(iyr2),
                 ipson=ipson, usera=usera, selfa=selfa, userc=userc, urmax=urmax,
                 outoff=outoff, ufloat=ufloat, sdd=sdd, rvod=rvod, stocko=stocko)
    # phase-2 handback: BMKILL per stand (ppmain.f label-150 second loop).
    kills = Vector{Vector{Float32}}(undef, n)
    for (i, m) in enumerate(members)
        t = m.trees
        wk2 = zeros(Float32, t.n)                 # no pre-existing FVS mortality in this composition
        bmkill!(stands[i], w, t, wk2, m.sp_alpha)
        kills[i] = wk2
    end
    return kills, ls, stands
end

"""
    ppe_run_landscape(stands; variant, keyargs...) -> Vector{PPEAggregate}

Faithful port of the PPMAIN mode-1 master-cycle landscape run. Orders the member
stands by the C11SRT processing order (`ppe_processing_order`), projects each with the
oracle-validated per-stand engine `run_keyfile` (GRSTND-equivalent), then runs the
bit-exact CMADDS/CMPRT2 composite-yield aggregation (`_ppe_aggregate`) per report year:
the before-thin TREES/volume columns are area-weighted means `Σ(v·w)/Σw`, while
ACC/MOR are area·period-weighted `Σ(v·w·prd)/Σ(w·prd)` (they coincide when every stand
shares one period). `variant` and any extra keyword args forward to `run_keyfile`.

Mode-1 (no interstand interaction) is the reachable regime and is bit-identical to
the master-cycle-stepped form (see the equivalence proof in the file header); every
per-stand value is oracle-validated. Mode-2 (interstand beetle) + MXHRVP are the
documented seams (`ppe_neighbors`, the wwpb outbreak) — not applied here.
"""
function ppe_run_landscape(stands::AbstractVector{PPEStand}; variant, keyargs...)
    isempty(stands) && return PPEAggregate[]
    order = ppe_processing_order(stands)              # C11SRT (inert for mode-1 output)
    # project each stand and collect its (area, parsed .sum rows), visiting stands in the
    # PPE processing order (order-independent for mode-1, but faithful), then run the
    # bit-exact CMADDS/CMPRT2 composite-yield aggregation.
    stand_rows = Tuple{Float64,Vector{Any}}[]
    for si in order
        st = stands[si]
        sumtext = run_keyfile(st.keyfile; variant = variant, keyargs...)
        rows = Any[]
        for line in split(sumtext, '\n')
            r = _ppe_parse_sum_row(line)
            r === nothing || push!(rows, r)
        end
        push!(stand_rows, (st.area, rows))
    end
    return _ppe_aggregate(stand_rows)
end

# =============================================================================
# MODE-2 LIVE interstand-beetle coupling (the master-cycle lockstep barrier).
# =============================================================================
# `ppe_run_landscape_live!` steps N stands IN LOCKSTEP: every member is projected by
# the ordinary per-stand engine (each_stand → setup_growth! → write_sum_file →
# grow_cycle!) but PAUSES at the WWPB seam (simulate.jl, post-MORTS/pre-GRADD, records
# un-tripled) via the `wwpb_barrier` hook. When all stands have reached that seam for a
# master cycle, ONE landscape `bmdrv_multi!` dispersal runs across the placed neighbors
# (redistributing beetle pressure over the SPLALO geometry / HXINDX), and each stand's
# beetle kill is handed back (`bmkill!` → WK2 → t.tpa) before the stand resumes GRADD.
# This is exactly what Fortran-FVS needed getstd/putstd for (swap one stand's COMMON in/
# out of a DA file to pause it); FVSjl's per-stand StandState makes it natural.
#
# The lockstep barrier is implemented with cooperative Julia Tasks + a Channel rendezvous
# (NOT OS threads): each stand runs on an `@async` Task; at the seam it deposits its
# freshly-binned WwpbStand and `put!`s an arrival token, then blocks on `take!` of its
# release channel. The coordinator (this function's own task) gathers all N arrivals, runs
# the ONE `bmdrv_multi!` cascade, then releases every task. Because @async is cooperatively
# scheduled (one task runs at a time, yielding only at the Channel ops) and the BMRANN RNG
# is advanced ONLY inside `bmdrv_multi!`/`bmistd!` (the coordinator, in fixed stand order —
# `bmsdit!`/`bmkill!` draw no random numbers), the run is fully DETERMINISTIC regardless of
# task scheduling. `wwpb_barrier === nothing` on every ordinary run keeps the single-stand
# path byte-identical (the guard test_multicycle stays 339/11).
#
# CONSTRAINTS (deadlock-safety of the fixed-N barrier): every member must share the same
# NUMCYCLE (so all reach exactly `ncyc` barriers in lockstep) and must NOT schedule a
# SIMFIRE (a fire cycle triples inside mortality_and_fire! and skips the WWPB seam ⇒ that
# stand would never arrive). Both hold for a WWPB dispersal landscape.
#
# At MXSTND=1 the barrier reproduces the single-stand DISPERSE path (`wwpb_apply!`) BYTE-
# IDENTICAL — `bmdrv_multi!` collapses to `bmatct_single!` (verified by test), and the
# bmsdit!/seed/bmkill! wrapping mirrors `wwpb_outbreak_cycle!`. See test_ppe_landscape_live.jl.

"""
One member of a LIVE (mode-2) PPE landscape: the keyfile to project, its area/spatial
location (SPLALO xloc/yloc in METERS, area in acres, STOCK flag), and its per-cycle beetle
seed (`seed_class`/`seed_tpa` — the synthetic inventory-damage kick-off, 0 = a pure neighbor
that only RECEIVES dispersed pressure). The outbreak year-window is landscape-global
(`iyr1`/`iyr2` on `ppe_run_landscape_live!`).
"""
struct PPELiveMember
    keyfile::String
    area::Float32
    xloc::Float32
    yloc::Float32
    stock::Bool
    seed_class::Int
    seed_tpa::Float32
end
PPELiveMember(keyfile::AbstractString; area::Real=1.0, xloc::Real=0.0, yloc::Real=0.0,
              stock::Bool=true, seed_class::Integer=0, seed_tpa::Real=0.0) =
    PPELiveMember(String(keyfile), Float32(area), Float32(xloc), Float32(yloc), stock,
                  Int(seed_class), Float32(seed_tpa))

# Build one member's projectable StandState (the run_keyfile prelude: first stand of the
# keyfile, notre! + setup_growth! + compute_volumes!). Returns (state, stand_id, mgmt_id).
function _ppe_member_state(keyfile::AbstractString; variant, faithful::Bool=true)
    local st = nothing
    for s in each_stand(keyfile; variant = variant, faithful = faithful)
        st = s; break                          # one stand per member keyfile
    end
    st === nothing && error("PPE live: keyfile $keyfile produced no stand")
    notre!(st); setup_growth!(st); compute_volumes!(st)
    sid = strip(st.plot.stand_id)
    mid = strip(st.plot.mgmt_id); mid = isempty(mid) ? "NONE" : String(mid)
    return (st, String(sid), mid)
end

# The species→alpha-code closure a stand's bmsdit!/bmkill! need (mirrors wwpb_apply!).
function _ppe_spalpha(s)
    code = s.coef.code_alpha
    return sp::Int -> (1 <= sp <= length(code)) ? code[sp] : ""
end

# The bmsdit! treelist input, snapshotted at the barrier (a mktrees-shaped NamedTuple over
# the CURRENT live records — after old_tpa restore, so tpa == cycle-start TPA). Reused
# verbatim by the standalone-decision recompute and read only by bmsdit!/bmkill!.
function _ppe_tree_snapshot(t)
    n = t.n
    return (n = n,
            species   = Int32[t.species[i]        for i in 1:n],
            dbh       = Float32[t.dbh[i]           for i in 1:n],
            tpa       = Float32[t.tpa[i]           for i in 1:n],
            height    = Float32[t.height[i]        for i in 1:n],
            crown_pct = Int32[t.crown_pct[i]       for i in 1:n],
            ht_growth = Float32[t.ht_growth[i]     for i in 1:n],
            cuft_vol  = Float32[t.cuft_vol[i]      for i in 1:n])
end

_ppe_seedvec(class::Integer, tpa::Real) = begin
    v = zeros(Float32, WWPB_NSCL)
    (1 <= class <= WWPB_NSCL) && (v[class] = Float32(tpa))
    v
end

"""
    ppe_run_landscape_live!(members; variant, iyr1, iyr2, period=5, dispersal-params...)
        -> NamedTuple

Run the LIVE mode-2 landscape (the lockstep master-cycle barrier described above). Returns
a NamedTuple:

  * `sums`         — per-stand `.sum` text (Vector{String}), landscape index order.
  * `kills`        — `kills[i][cyc]` the per-record beetle-reconciled WK2 the barrier applied
                     to stand `i` at master cycle `cyc` (the in-flight dispersal decision).
  * `snaps`        — `snaps[i][cyc]` the bmsdit! treelist snapshot at that barrier (the input
                     the standalone cascade is re-fed for the equivalence proof).
  * `fvsmort`      — `fvsmort[i][cyc]` the FVS density/background mortality WK2 baseline the
                     beetle kill is MAX-combined with (bmkill! seed).
  * `landscapes`   — per-cycle `WwpbLandscape` (BKPOUT/BKPIN/SELFBKP dispersal ledger).
  * `stand_ids`    — per-stand id; `ncyc` — barriers per stand; `w0`/`coeffs` — the shared
                     dispersal state's ORIGINAL parameters (for the standalone recompute).

Every dispersal kernel (`bmsdit!`/`bmdrv_multi!`/`bmatct_multi!`/`bmkill!`) is bit-exact vs
pristine Fortran; this function is only the orchestration plumbing, proven zero-divergence by
`ppe_landscape_replay` + `ppe_standalone_decisions` (see the equivalence test).
"""
function ppe_run_landscape_live!(members::AbstractVector{PPELiveMember}; variant,
                                 iyr1::Integer, iyr2::Integer, period::Integer = 5,
                                 faithful::Bool = true,
                                 date::AbstractString = "01-01-2026",
                                 time::AbstractString = "00:00:00",
                                 usera::NTuple{3,Float32}=(1.0f0,1.0f0,1.0f0),
                                 selfa::NTuple{3,Float32}=(1.0f0,1.0f0,1.0f0),
                                 userc::NTuple{3,Float32}=(100.0f0,100.0f0,100.0f0),
                                 urmax::NTuple{3,Float32}=(15.0f0,15.0f0,15.0f0),
                                 outoff::Bool=true)
    n = length(members)
    n == 0 && error("PPE live: empty landscape")

    # shared landscape dispersal state (one BMRANN stream, matching bmdrv_multi!/PPMAIN).
    w = wwpb_defaults!(variant)
    w.pbspec = Int32(1); w.iyr1 = Int32(iyr1); w.iyr2 = Int32(iyr2)
    coeffs = wwpb_init_coeffs(w.upsiz)
    w0_seed = w.rng_s0                                # ORIGINAL RNG (for the standalone recompute)

    # build every member's StandState + landscape geometry.
    states  = Vector{Any}(undef, n); stand_ids = Vector{String}(undef, n)
    mgmt_ids = Vector{String}(undef, n); spαs = Vector{Any}(undef, n)
    xloc = Float32[]; yloc = Float32[]; area = Float32[]; stock = Bool[]; seed = Vector{Float32}[]
    for (i, m) in enumerate(members)
        st, sid, mid = _ppe_member_state(m.keyfile; variant = variant, faithful = faithful)
        states[i] = st; stand_ids[i] = sid; mgmt_ids[i] = mid; spαs[i] = _ppe_spalpha(st)
        # sarea == the single-stand harness's sarea (gross_space, else 1) so MXSTND=1 collapses.
        sar = st.plot.gross_space > 0.0f0 ? st.plot.gross_space : 1.0f0
        push!(area, sar); push!(xloc, m.xloc); push!(yloc, m.yloc); push!(stock, m.stock)
        push!(seed, _ppe_seedvec(m.seed_class, m.seed_tpa))
    end
    ncyc = Int(states[1].control.ncycle_eff); ncyc < 1 && (ncyc = Int(states[1].control.ncycle))

    # per-cycle records + rendezvous channels.
    kills   = [Vector{Vector{Float32}}() for _ in 1:n]
    snaps   = [Vector{Any}()             for _ in 1:n]
    fvsmort = [Vector{Vector{Float32}}() for _ in 1:n]
    landscapes = WwpbLandscape[]
    dep = Vector{WwpbStand}(undef, n)
    arrivals = Channel{Int}(n)
    releases = [Channel{Nothing}(1) for _ in 1:n]

    make_barrier(i) = (s, old_tpa, fint) -> begin
        t = s.trees; nn = t.n
        fmort = Float32[old_tpa[k] - t.tpa[k] for k in 1:nn]
        @inbounds for k in 1:nn; t.tpa[k] = old_tpa[k]; end   # beetle works on cycle-start stand
        push!(snaps[i], _ppe_tree_snapshot(t)); push!(fvsmort[i], copy(fmort))
        stw = WwpbStand()
        if nn > 0
            bmsdit!(stw, t, w, spαs[i]); fill!(stw.rvdsc, 1.0f0); stw.pbkill .= seed[i]
        else
            fill!(stw.rvdsc, 1.0f0)
        end
        dep[i] = stw
        put!(arrivals, i); take!(releases[i])                 # rendezvous: coordinator runs bmdrv_multi!
        wk2 = copy(fmort)
        nn > 0 && bmkill!(dep[i], w, t, wk2, spαs[i])
        @inbounds for k in 1:nn; t.tpa[k] = max(old_tpa[k] - wk2[k], 0.0f0); end
        push!(kills[i], wk2)
    end

    tasks = Vector{Task}(undef, n)
    for i in 1:n
        tasks[i] = @async begin
            io = IOBuffer()
            write_sum_file(io, states[i]; period = Int(period), stand_id = stand_ids[i],
                           mgmt_id = mgmt_ids[i], variant = variant_code(states[i].variant),
                           date = date, time = time, wwpb_barrier = make_barrier(i))
            String(take!(io))
        end
    end

    for _cyc in 1:ncyc                                        # coordinator (lockstep master cycles)
        for _ in 1:n; take!(arrivals); end                   # wait for every stand at the seam
        ls = WwpbLandscape(xloc, yloc, area, stock)
        bmdrv_multi!(ls, dep, w, coeffs; area = area, iyr1 = Int(iyr1), iyr2 = Int(iyr2),
                     usera = usera, selfa = selfa, userc = userc, urmax = urmax, outoff = outoff)
        push!(landscapes, ls)
        for i in 1:n; put!(releases[i], nothing); end         # release all → resume GRADD
    end

    sums = String[fetch(tasks[i]) for i in 1:n]
    return (sums = sums, kills = kills, snaps = snaps, fvsmort = fvsmort,
            landscapes = landscapes, stand_ids = stand_ids, ncyc = ncyc,
            seed = seed, area = area, xloc = xloc, yloc = yloc, stock = stock,
            spαs = spαs, w0_seed = w0_seed, iyr1 = Int(iyr1), iyr2 = Int(iyr2),
            usera = usera, selfa = selfa, userc = userc, urmax = urmax, outoff = outoff)
end

"""
    ppe_standalone_decisions(res; variant) -> kills

Recompute the per-cycle dispersal decisions OUTSIDE the engine, using the ALREADY-VALIDATED
standalone `bmdrv_multi!` cascade fed the exact per-cycle barrier treelists `res.snaps` that
the LIVE run recorded. A fresh landscape BMRANN stream is stepped from the same initial seed;
each cycle rebins the recorded snapshots (`bmsdit!`), runs the one landscape dispersal, and
hands back `bmkill!` on the same FVS-mortality baseline (`res.fvsmort`). This is the "premade
decisions" side of the equivalence proof: identical kernels, but the decisions are computed
standalone from the barrier treelists rather than in-flight. Returns `kills[i][cyc]`.
"""
function ppe_standalone_decisions(res; variant)
    n = length(res.snaps); ncyc = res.ncyc
    w = wwpb_defaults!(variant)
    w.pbspec = Int32(1); w.iyr1 = Int32(res.iyr1); w.iyr2 = Int32(res.iyr2)
    w.rng_s0 = res.w0_seed
    coeffs = wwpb_init_coeffs(w.upsiz)
    kills = [Vector{Vector{Float32}}() for _ in 1:n]
    for cyc in 1:ncyc
        stands = WwpbStand[]
        for i in 1:n
            snap = res.snaps[i][cyc]; stw = WwpbStand()
            if snap.n > 0
                bmsdit!(stw, snap, w, _ppe_snap_alpha(res, i)); fill!(stw.rvdsc, 1.0f0)
                stw.pbkill .= res.seed[i]
            else
                fill!(stw.rvdsc, 1.0f0)
            end
            push!(stands, stw)
        end
        ls = WwpbLandscape(res.xloc, res.yloc, res.area, res.stock)
        bmdrv_multi!(ls, stands, w, coeffs; area = res.area, iyr1 = res.iyr1, iyr2 = res.iyr2,
                     usera = res.usera, selfa = res.selfa, userc = res.userc, urmax = res.urmax,
                     outoff = res.outoff)
        for i in 1:n
            snap = res.snaps[i][cyc]; wk2 = copy(res.fvsmort[i][cyc])
            snap.n > 0 && bmkill!(stands[i], w, snap, wk2, _ppe_snap_alpha(res, i))
            push!(kills[i], wk2)
        end
    end
    return kills
end

# the alpha-code closure for a recorded member (rebuilt from the member's own StandState is
# unavailable post-run; the landscape shares one species→code map per variant, so recover it
# from the live states via a stored closure). Stored on the result to keep decisions faithful.
_ppe_snap_alpha(res, i) = res.spαs[i]

"""
    ppe_landscape_replay(members, kills; variant, period=5, ...) -> Vector{String}

The PREMADE-decision projection: re-project each member stand INDEPENDENTLY (a plain
sequential loop — NO Tasks, NO channel, NO cross-stand coupling), injecting the fixed
per-cycle beetle kill `kills[i][cyc]` at the WWPB seam instead of running any dispersal.
The per-stand growth path is the identical engine (write_sum_file → grow_cycle!); the ONLY
thing supplied from outside is the pre-computed WK2. If the LIVE concurrent barrier had
perturbed any stand's projection beyond delivering its kill vector (a race, a mis-routed
kill, an RNG/state leak, an ordering effect), this independent replay would DIVERGE. Byte-
identical `sums` therefore prove the lockstep plumbing (barrier + state extraction + kill
handback) adds zero divergence — mirroring the mode-1 stepped==independent equivalence.
"""
function ppe_landscape_replay(members::AbstractVector{PPELiveMember}, kills; variant,
                              period::Integer = 5, faithful::Bool = true,
                              date::AbstractString = "01-01-2026",
                              time::AbstractString = "00:00:00")
    n = length(members)
    sums = String[]
    for i in 1:n
        st, sid, mid = _ppe_member_state(members[i].keyfile; variant = variant, faithful = faithful)
        cyc = Ref(0); ki = kills[i]
        barrier = (s, old_tpa, fint) -> begin
            cyc[] += 1
            t = s.trees; wk2 = ki[cyc[]]
            @inbounds for k in 1:t.n; t.tpa[k] = max(old_tpa[k] - wk2[k], 0.0f0); end
        end
        io = IOBuffer()
        write_sum_file(io, st; period = Int(period), stand_id = sid, mgmt_id = mid,
                       variant = variant_code(st.variant), date = date, time = time,
                       wwpb_barrier = barrier)
        push!(sums, String(take!(io)))
    end
    return sums
end
