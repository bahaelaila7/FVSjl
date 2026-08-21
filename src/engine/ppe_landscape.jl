# =============================================================================
# PPE (Parallel Processing Extension) — LANDSCAPE / multi-stand harness.
# =============================================================================
# USER-directed 2026-08-21 ("reconstruct a PPE harness like WWPB"). The PPE outer
# routines (PPMAIN / ALSTD2 / SPLAEX) have NO source in the FVS tree — only the
# COMMON-block DECLARATIONS survive (archive/PPEcommons/PPEPRM.F77 + PPEXCM.F77).
# So, exactly as the WWPB landscape harness (src/engine/wwpb_landscape.jl) did for
# the absent PPMAIN/ALSTD2/SPLAEX beetle orchestration, this is a FAITHFUL
# RECONSTRUCTION of PPE's documented behavior, NOT a bit-exact port (there is no
# oracle to run — the code is absent). It composes the ALREADY-VALIDATED per-stand
# projection (run_keyfile) into a landscape.
#
# PPE data model (PPEPRM.F77): a landscape holds up to MXSTND=10000 stands, each
# area-weighted (TOTALWT), stepped through master cycles (MSPERIOD). PPE's
# "PPE-defined variables" (PPEXCM.F77 PTSTV1) are the AREA-WEIGHTED landscape
# aggregates evaluated each master cycle — the reconstructed core here:
#   (1) AVBTPA   average before-thin trees/acre
#   (2) AVBTCUFT average before-thin total cubic volume
#   (3) AVBMCUFT average before-thin merch cubic volume
#   (4) AVBBDFT  average before-thin board-foot volume
#   (5) AVBACC   last-cycle accretion (cuft/acre/yr)
#   (6) AVBMORT  last-cycle mortality (cuft/acre/yr)
#   (7) MSPERIOD master cycle period length (years)
#   (8) TOTALWT  total area (weight)
# Each AVB* = Σ_stands(value · area) / Σ_stands(area), i.e. the area-weighted mean
# of the per-stand `.sum` before-thin values (PPE's SPLAEX/ALSTD2 aggregation).
# ADDITIVE + INERT: a new engine file with no simulate.jl seam — a single stand
# projects byte-identically (gate: multicycle 339/11). Multi-stand harvest policies
# (MXHRVP) + spatial neighbors (MXNBS/SPNBCM) are documented PPE extension points,
# left as reconstruction stubs (they need their own absent-source reconstruction).
# =============================================================================

"""
One member stand of a PPE landscape: a keyfile to project and its area weight
(acres). `area` is PPE's per-stand contribution to TOTALWT (area-weighting).
"""
struct PPEStand
    keyfile::String
    area::Float64
end
PPEStand(keyfile::AbstractString; area::Real = 1.0) = PPEStand(String(keyfile), Float64(area))

"""
The area-weighted PPE landscape aggregate for one reporting year (PPEXCM PTSTV1).
`nstands` is how many member stands reported that year (aligned on `.sum` year).
"""
struct PPEAggregate
    year::Int
    avbtpa::Float64      # PTSTV1(1)
    avbtcuft::Float64    # PTSTV1(2)
    avbmcuft::Float64    # PTSTV1(3)
    avbbdft::Float64     # PTSTV1(4)
    avbacc::Float64      # PTSTV1(5)
    avbmort::Float64     # PTSTV1(6)
    totalwt::Float64     # PTSTV1(8)  Σ area over reporting stands
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
    ppe_run_landscape(stands; variant, keyargs...) -> Vector{PPEAggregate}

Reconstructed PPE landscape run: project each member `PPEStand` with `run_keyfile`
(the validated per-stand engine), then form the area-weighted PPE landscape
aggregates (PPEXCM PTSTV1) per reporting year — SPLAEX/ALSTD2's job. `variant` and
any extra keyword args are forwarded to `run_keyfile`. Behavior-faithful (there is
no PPE oracle, its source being absent); the per-stand values ARE oracle-validated,
and the aggregate is the documented Σ(value·area)/Σ(area).
"""
function ppe_run_landscape(stands::AbstractVector{PPEStand}; variant, keyargs...)
    isempty(stands) && return PPEAggregate[]
    # accumulate area-weighted sums per year (SPLAEX aggregation)
    wtpa = Dict{Int,Float64}(); wcuft = Dict{Int,Float64}(); wmcuft = Dict{Int,Float64}()
    wbdft = Dict{Int,Float64}(); wacc = Dict{Int,Float64}(); wmort = Dict{Int,Float64}()
    wsum = Dict{Int,Float64}(); ncnt = Dict{Int,Int}()
    for st in stands
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
    out = PPEAggregate[]
    for yr in sort(collect(keys(wsum)))
        w = wsum[yr]
        w <= 0.0 && continue
        push!(out, PPEAggregate(yr, wtpa[yr]/w, wcuft[yr]/w, wmcuft[yr]/w, wbdft[yr]/w,
                                wacc[yr]/w, wmort[yr]/w, w, ncnt[yr]))
    end
    return out
end
