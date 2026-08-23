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
# accretion/mortality read from the END (end-4 / end-3) so they are robust to the
# volume-block width. Returns nothing for header/blank/non-data lines.
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
        acc  = parse(Float64, f[end-4])
        mort = parse(Float64, f[end-3])
        return (year = year, tpa = tpa, cuft = cuft, mcuft = mcuft, bdft = bdft, acc = acc, mort = mort)
    catch
        return nothing
    end
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
    ppe_run_landscape(stands; variant, keyargs...) -> Vector{PPEAggregate}

Faithful port of the PPMAIN mode-1 master-cycle landscape run. Orders the member
stands by the C11SRT processing order (`ppe_processing_order`), projects each with the
oracle-validated per-stand engine `run_keyfile` (GRSTND-equivalent), then forms the
area-weighted PPE landscape aggregates (PPEXCM PTSTV1(1..9)) per report year —
SPLAEX/ALSTD2's aggregation, each `AVB* = Σ(value·area)/Σ(area)`. `variant` and any
extra keyword args forward to `run_keyfile`.

Mode-1 (no interstand interaction) is the reachable regime and is bit-identical to
the master-cycle-stepped form (see the equivalence proof in the file header); every
per-stand value is oracle-validated. Mode-2 (interstand beetle) + MXHRVP are the
documented seams (`ppe_neighbors`, the wwpb outbreak) — not applied here.
"""
function ppe_run_landscape(stands::AbstractVector{PPEStand}; variant, keyargs...)
    isempty(stands) && return PPEAggregate[]
    order = ppe_processing_order(stands)              # C11SRT (inert for mode-1 output)
    # accumulate area-weighted sums per year (SPLAEX aggregation), visiting stands in
    # the PPE processing order (order-independent for mode-1, but faithful).
    wtpa = Dict{Int,Float64}(); wcuft = Dict{Int,Float64}(); wmcuft = Dict{Int,Float64}()
    wbdft = Dict{Int,Float64}(); wacc = Dict{Int,Float64}(); wmort = Dict{Int,Float64}()
    wsum = Dict{Int,Float64}(); ncnt = Dict{Int,Int}()
    for si in order
        st = stands[si]
        sumtext = run_keyfile(st.keyfile; variant = variant, keyargs...)
        for line in split(sumtext, '\n')
            r = _ppe_parse_sum_row(line)
            r === nothing && continue
            a = st.area
            wtpa[r.year]  = get(wtpa, r.year, 0.0)  + r.tpa  * a
            wcuft[r.year] = get(wcuft, r.year, 0.0) + r.cuft * a
            wmcuft[r.year]= get(wmcuft, r.year, 0.0)+ r.mcuft* a
            wbdft[r.year] = get(wbdft, r.year, 0.0) + r.bdft * a
            wacc[r.year]  = get(wacc, r.year, 0.0)  + r.acc  * a
            wmort[r.year] = get(wmort, r.year, 0.0) + r.mort * a
            wsum[r.year]  = get(wsum, r.year, 0.0)  + a
            ncnt[r.year]  = get(ncnt, r.year, 0)    + 1
        end
    end
    yrs = sort(collect(keys(wsum)))
    out = PPEAggregate[]
    for (i, yr) in enumerate(yrs)
        w = wsum[yr]
        w <= 0.0 && continue
        # PTSTV1(7) MSPERIOD = the master-cycle period = years to the next report row
        # (0 on the final row). Landscape-level constant in PPE; here derived from the
        # reporting cadence the member stands share.
        msperiod = i < length(yrs) ? (yrs[i+1] - yr) : 0
        push!(out, PPEAggregate(yr, wtpa[yr]/w, wcuft[yr]/w, wmcuft[yr]/w, wbdft[yr]/w,
                                wacc[yr]/w, wmort[yr]/w, msperiod, w,
                                0.0,               # PTSTV1(9) OLDTARG: no landscape harvest target set
                                ncnt[yr]))
    end
    return out
end
