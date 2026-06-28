# =============================================================================
# econ.jl — economic analysis core (C8, ECON extension)
#
# Ported from: bin/FVSsn_buildDir/eccalc.f (the ECON discounting / present-value core).
#
# The ECON extension values a stand's management over time: costs (fixed annual + per-
# harvest variable) and revenues (per harvest, by product/diameter) are discounted to a
# present net value (PNV), with a benefit/cost ratio and a realizable rate of return.
# These are the self-contained financial kernels; the keyword layer (ANNUCST/HRVVRCST/
# HRVRVN that build the per-year cost/revenue streams from the harvested volumes) wraps
# them. Values use Float32 to match the Fortran REAL arithmetic.
# =============================================================================

"""
    econ_present_value(amt, time, rate) -> Float32

Present value of `amt` discounted `time` years at annual `rate` (computePV, eccalc.f:999):
`amt / (1+rate)^time`. Non-positive amounts return 0.
"""
@inline econ_present_value(amt::Float32, time::Integer, rate::Float32)::Float32 =
    amt <= 0f0 ? 0f0 : amt / (1f0 + rate)^time

"""
    econ_pnv(undisc_cost, undisc_rev, rate) -> (; pnv, disc_cost, disc_rev)

Present net value over the analysis horizon (computePNV, eccalc.f:1006). `undisc_cost`
and `undisc_rev` are the undiscounted cost/revenue in each year 1…endTime. Costs accrue
at the *start* of each year (discounted `i−1` years), revenues at the *end* (discounted
`i` years); `pnv = disc_rev − disc_cost`.
"""
function econ_pnv(undisc_cost::AbstractVector{Float32}, undisc_rev::AbstractVector{Float32},
                  rate::Float32)
    disc_cost = 0f0; disc_rev = 0f0
    n = min(length(undisc_cost), length(undisc_rev))
    @inbounds for i in 1:n
        disc_cost += econ_present_value(undisc_cost[i], i - 1, rate)
        disc_rev  += econ_present_value(undisc_rev[i],  i,     rate)
    end
    return (; pnv = disc_rev - disc_cost, disc_cost = disc_cost, disc_rev = disc_rev)
end

# FVS computes the B/C ratio + RRR only when the discounted cost (and, for RRR, revenue) exceeds NEAR_ZERO
# (eccalc.f:58 `NEAR_ZERO=0.01`, :681/:685), else leaves them blank — guarding against tiny-cost blow-ups.
const ECON_NEAR_ZERO = 0.01f0

"Benefit/cost ratio = discounted revenue / discounted cost (eccalc.f:683); 0 unless cost > NEAR_ZERO (0.01)."
@inline econ_bc_ratio(disc_rev::Float32, disc_cost::Float32)::Float32 =
    disc_cost > ECON_NEAR_ZERO ? disc_rev / disc_cost : 0f0

"""
    econ_rate_of_return(disc_rev, disc_cost, endtime, rate) -> Float32

Realizable rate of return (%) (eccalc.f:686): `100·((revDisc/costDisc)^(1/endTime)·
(1+rate) − 1)`. Zero unless both discounted streams exceed NEAR_ZERO (0.01, eccalc.f:681/685).
"""
@inline function econ_rate_of_return(disc_rev::Float32, disc_cost::Float32,
                                     endtime::Integer, rate::Float32)::Float32
    (disc_cost > ECON_NEAR_ZERO && disc_rev > ECON_NEAR_ZERO) || return 0f0
    return 100f0 * ((disc_rev / disc_cost)^(1f0 / endtime) * (1f0 + rate) - 1f0)
end

"""
    econ_sev(net_rotation, rate, rotation) -> Float32

Soil expectation value (Faustmann land value): the present value of an infinite series
of identical rotations, each returning a present net value `net_rotation` over `rotation`
years at annual `rate` — `net_rotation·(1+rate)^rotation / ((1+rate)^rotation − 1)`.
"""
@inline function econ_sev(net_rotation::Float32, rate::Float32, rotation::Integer)::Float32
    f = (1f0 + rate)^rotation
    return f > 1f0 ? net_rotation * f / (f - 1f0) : 0f0
end

"""
    econ_forest_value(pnv, sev_input, rate, endtime) -> (; forest_value, reprod_value)

Forest and reproduction (bare-land) value at the end of the analysis given a known
end-of-horizon land value `sev_input` (eccalc.f:649-655): the SEV is discounted back
`endtime` years and added to the management `pnv`; the reproduction value subtracts the
starting land value.
"""
@inline function econ_forest_value(pnv::Float32, sev_input::Float32, rate::Float32, endtime::Integer)
    disc_sev = econ_present_value(sev_input, endtime, rate)
    fv = pnv + disc_sev
    return (; forest_value = fv, reprod_value = fv - sev_input)
end

# ECON quantity-units of measure (ECNCOM.F77:19).
const ECON_TPA      = Int32(1)   # per tree
const ECON_BF_1000  = Int32(2)   # per thousand board feet (whole tree)
const ECON_FT3_100  = Int32(3)   # per hundred cubic feet (whole tree)
const ECON_BF_1000_LOG = Int32(4)   # per MBF, log-graded by inside-bark DIB class (echarv.f BF_1000_LOG)
const ECON_FT3_100_LOG = Int32(5)   # per CCF, log-graded by DIB class (FT3_100_LOG; needs per-log cubic)

"""
    harvest_value(records, sp, dbh, tpa, cuft, bdft) -> Float32

Cost or revenue from harvesting one tree record (species `sp`, DBH `dbh`, `tpa` trees/ac,
`cuft` cubic and `bdft` board feet per tree) against a list of `EconCostRev` (echarv.f /
eccalc.f): each record whose species (0 = all) and DBH range `[lo, hi)` match the tree
contributes `amount × volume`, where the volume is per-tree (TPA), per-MBF
(BF_1000 = bdft·tpa/1000), or per-CCF (FT3_100 = cuft·tpa/100). The LOG-graded units 4/5
(BF_1000_LOG / FT3_100_LOG) need the ecvol.f per-log bucking + echarv.f DIB-grade bucketing
subsystem — UNPORTED, so they fall through to `vol=0` here (kw_econ! warns at parse time).
Inert for sn.key: live FVS also yields 0 revenue there (the harvests don't qualify).
"""
function harvest_value(records::AbstractVector{EconCostRev}, sp::Integer, dbh::Float32,
                       tpa::Float32, cuft::Float32, bdft::Float32)::Float32
    total = 0f0
    @inbounds for r in records
        ((r.sp == 0 || r.sp == sp) && r.dbh_lo <= dbh < r.dbh_hi) || continue
        vol = r.unit == ECON_TPA     ? tpa :
              r.unit == ECON_BF_1000 ? bdft * tpa / 1000f0 :
              r.unit == ECON_FT3_100 ? cuft * tpa / 100f0  : 0f0
        total += r.amount * vol
    end
    return total
end

"""
    _log_dia_grp(records, sp, value) -> Union{Int32,Nothing}

ECHARV `getDiaGrp` (echarv.f:146): among the BF_1000_LOG (unit 4) HRVRVN records that apply to
species `sp` (an explicit species record, or `sp==0` = ALL), return the DIB class of the largest
class-lower-bound `≤ value`, encoded as `Int32(round(lo·10))`; `nothing` when `value` is below every
class. Mirrors getDiaGrp's "largest class ≤ value" pick over the descending-sorted class list.
"""
function _log_dia_grp(records::AbstractVector{EconCostRev}, sp::Integer, value::Float32,
                      unit::Int32 = ECON_BF_1000_LOG)
    best = -Inf32
    @inbounds for r in records
        (r.unit == unit && (r.sp == 0 || r.sp == sp)) || continue
        (r.dbh_lo <= value && r.dbh_lo > best) && (best = r.dbh_lo)
    end
    best == -Inf32 && return nothing
    return Int32(round(best * 10f0))
end

"""
    accrue_log_grade!(ec, sp, year, harvTpa, netTreeBf, by_dib)

Accumulate one removed tree's log-graded board feet into `ec.log_grade_rev` (echarv.f:70-100, the
BF_1000_LOG path). `by_dib` is the tree's per-log-DIB gross Scribner BF (from `_r8_scribner_bf_by_dib`,
stashed in `compute_volumes!`). Each log's gross BF is scaled by the whole-tree defect proportion
`defProp = ΣgrossLogBf / netTreeBf` and the removed `harvTpa`, then bucketed by the log-end DIB's class
(`_log_dia_grp`). Keyed `(year, species, classLo·10)`. REPORT-ONLY — feeds FVS_EconHarvestValue, not PNV.
"""
function accrue_log_grade!(ec::EconState, sp::Integer, year::Int32, harvTpa::Float32,
                           netTreeBf::Float32, by_dib::Dict{Int,Float32})
    (netTreeBf > 0f0 && !isempty(by_dib)) || return
    treeVol = 0f0
    for v in values(by_dib); treeVol += v; end
    treeVol > 0f0 || return
    defProp = treeVol / netTreeBf
    @inbounds for (idib, gbf) in by_dib
        cls = _log_dia_grp(ec.hrv_rev, sp, Float32(idib))
        cls === nothing && continue
        key = (year, Int32(sp), cls)
        ec.log_grade_rev[key] = get(ec.log_grade_rev, key, 0f0) + gbf * harvTpa / defProp
    end
    return
end

"""
    accrue_log_grade_cuft!(ec, sp, year, harvTpa, netTreeFt3, by_dib)

Cubic (FT3_100_LOG, unit 5) analog of `accrue_log_grade!` (echarv.f:103-132). `by_dib` is the tree's
per-log-DIB gross cubic feet (from `_r8_cuft_by_dib`); `netTreeFt3` is the tree's net merchantable cubic
(echarv.f `ft3PerTree`). Each log's gross cuft is scaled by `defProp = ΣgrossLogFt3 / netTreeFt3` and the
removed `harvTpa`, bucketed by the log-end DIB class, keyed `(year, species, classLo·10)` in `log_grade_ft3`.
"""
function accrue_log_grade_cuft!(ec::EconState, sp::Integer, year::Int32, harvTpa::Float32,
                                netTreeFt3::Float32, by_dib::Dict{Int,Float32})
    (netTreeFt3 > 0f0 && !isempty(by_dib)) || return
    treeVol = 0f0
    for v in values(by_dib); treeVol += v; end
    treeVol > 0f0 || return
    defProp = treeVol / netTreeFt3
    @inbounds for (idib, gft3) in by_dib
        cls = _log_dia_grp(ec.hrv_rev, sp, Float32(idib), ECON_FT3_100_LOG)
        cls === nothing && continue
        key = (year, Int32(sp), cls)
        ec.log_grade_ft3[key] = get(ec.log_grade_ft3, key, 0f0) + gft3 * harvTpa / defProp
    end
    return
end

"""
    econ_harvest_value_rows(ec, coef) -> Vector{<:NamedTuple}

Build the FVS_EconHarvestValue rows (eccalc.f:760-852) from the accumulated log-graded revenue —
both board (`ec.log_grade_rev`, unit 4 BF_1000_LOG) and cubic (`ec.log_grade_ft3`, unit 5 FT3_100_LOG).
One row per (year, species, unit, DIB class) in FVS's species → unit (board before cubic) → ascending-DIB
order: `Min_DIB` = the class lower bound, `Max_DIB` = the next-larger class's lower bound (or 999.9 for the
top class, eccalc.f:776). Board rows fill `Board_Ft_Removed = nint(Σbf)`, `Board_Ft_Value = Total_Value =
nint(Σbf/1000·price)`; cubic rows fill `Ft3_Removed = nint(Σft3)`, `Ft3_Value = Total_Value =
nint(Σft3/100·price)` (the /100 mirrors board's /1000, eccalc.f:801/813/819). The unused unit's volume
columns + the TPA/Tons/DBH columns are left null. `price` is the class's HRVRVN amount.
"""
function econ_harvest_value_rows(ec::EconState, coef)
    rows = NamedTuple[]
    (isempty(ec.log_grade_rev) && isempty(ec.log_grade_ft3)) && return rows
    # one emitter shared by both units; `scale` is the per-unit price divisor (1000 BF / 100 CuFt).
    function _emit!(acc, unit::Int32, scale::Float32, cubic::Bool)
        for key in sort!(collect(keys(acc)))               # (year, speciesIdx, classLo·10) → species-index order
            yr, sp, clsx10 = key
            vol = acc[key]
            lo = Float32(clsx10) / 10f0
            los = sort!(unique(Float32[r.dbh_lo for r in ec.hrv_rev
                        if r.unit == unit && (r.sp == 0 || r.sp == sp)]))
            nexti = findfirst(>(lo), los)
            maxdib = nexti === nothing ? 999.9f0 : los[nexti]
            price = 0f0
            for r in ec.hrv_rev
                (r.unit == unit && r.dbh_lo == lo && (r.sp == sp || r.sp == 0)) || continue
                price = r.amount
                r.sp == sp && break                         # species-specific price overrides ALL
            end
            spi = Int(sp)
            val = round(Int, vol / scale * price, RoundNearestTiesAway)
            push!(rows, (year = Int(yr),
                         sp_fvs = String(strip(coef.code_alpha[spi])),
                         sp_plants = String(strip(coef.code_plants[spi])),
                         sp_fia = String(strip(coef.code_fia[spi])),
                         min_dib = Float64(lo), max_dib = Float64(maxdib),
                         ft3_removed = cubic ? round(Int, vol, RoundNearestTiesAway) : missing,
                         ft3_value   = cubic ? val : missing,
                         bf_removed  = cubic ? missing : round(Int, vol, RoundNearestTiesAway),
                         bf_value    = cubic ? missing : val,
                         total_value = val))
        end
    end
    _emit!(ec.log_grade_rev, ECON_BF_1000_LOG, 1000f0, false)   # FVS unit order: board (4) before cubic (5)
    _emit!(ec.log_grade_ft3, ECON_FT3_100_LOG, 100f0, true)
    return rows
end

"""
    econ_value_harvest(ec, sp, dbh, tpa, cuft, bdft) -> (; cost, revenue)

Total variable harvest cost and revenue for a set of cut trees (parallel arrays of
species / DBH / removed TPA / per-tree cubic & board feet) against the stand's ECON
cost/revenue tables (`ec.hrv_cost`, `ec.hrv_rev`). These are the per-harvest cash flows
that accumulate into the undiscounted streams the discounting core (`econ_pnv`) values.
"""
function econ_value_harvest(ec::EconState, sp::AbstractVector, dbh::AbstractVector{Float32},
                            tpa::AbstractVector{Float32}, cuft::AbstractVector{Float32},
                            bdft::AbstractVector{Float32})
    cost = 0f0; revenue = 0f0
    @inbounds for i in eachindex(sp)
        tpa[i] > 0f0 || continue
        s = Int(sp[i])
        cost    += harvest_value(ec.hrv_cost, s, dbh[i], tpa[i], cuft[i], bdft[i])
        revenue += harvest_value(ec.hrv_rev,  s, dbh[i], tpa[i], cuft[i], bdft[i])
    end
    return (; cost, revenue)
end

"""
    record_harvest!(ec, year, sp, dbh, removed_tpa, cuft, bdft)

Value a cycle's harvest (the per-tree removed TPA from `cuts!`) against the ECON tables
and append the `(year, cost, revenue)` cash flow to `ec.harvests` (no-op if nothing was
cut or no cost/revenue resulted). The discounting core (`econ_stand_pnv`) values the
accumulated flows.
"""
function record_harvest!(ec::EconState, year::Integer, sp::AbstractVector,
                         dbh::AbstractVector{Float32}, removed_tpa::AbstractVector{Float32},
                         cuft::AbstractVector{Float32}, bdft::AbstractVector{Float32})
    v = econ_value_harvest(ec, sp, dbh, removed_tpa, cuft, bdft)
    (v.cost > 0f0 || v.revenue > 0f0) && push!(ec.harvests, (Float32(year), v.cost, v.revenue))
    return
end

"""
    econ_stand_pnv(ec, end_year) -> (; pnv, disc_cost, disc_rev)

Present net value of the stand's management from `ec.base_year` to `end_year` (eccalc.f):
the annual cost `ann_cost` accrues every year (discounted at the year's start) and each
recorded harvest's cost/revenue is discounted from its year, at `ec.discount_rate`.
"""
function econ_stand_pnv(ec::EconState, end_year::Integer)
    rate = ec.discount_rate
    base = Int(ec.base_year)
    (base < 0 || end_year < base) && return (; pnv = 0f0, disc_cost = 0f0, disc_rev = 0f0)
    disc_cost = 0f0; disc_rev = 0f0
    # annual management cost: accrues at the start of each year (time t = 0 … end-base-1)
    if ec.ann_cost > 0f0
        for t in 0:(end_year - base - 1)
            disc_cost += econ_present_value(ec.ann_cost, t, rate)
        end
    end
    # harvest cash flows (revenue at end of its year, cost at start)
    for (yr, cost, rev) in ec.harvests
        t = Int(yr) - base
        # eccalc.f:114-117: cost accrues at START of year → discounted beginTime-1 = (yr-startYear) = t yrs;
        # revenue accrues at END of year → discounted beginTime = (yr-startYear+1) = t+1 yrs. (base = startYear.)
        disc_cost += econ_present_value(cost, max(0, t),     rate)
        disc_rev  += econ_present_value(rev,  max(0, t + 1), rate)
    end
    return (; pnv = disc_rev - disc_cost, disc_cost, disc_rev)
end
