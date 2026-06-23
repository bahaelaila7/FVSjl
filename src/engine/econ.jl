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

"Benefit/cost ratio = discounted revenue / discounted cost (eccalc.f:683)."
@inline econ_bc_ratio(disc_rev::Float32, disc_cost::Float32)::Float32 =
    disc_cost > 0f0 ? disc_rev / disc_cost : 0f0

"""
    econ_rate_of_return(disc_rev, disc_cost, endtime, rate) -> Float32

Realizable rate of return (%) (eccalc.f:686): `100·((revDisc/costDisc)^(1/endTime)·
(1+rate) − 1)`. Zero unless both discounted streams are positive.
"""
@inline function econ_rate_of_return(disc_rev::Float32, disc_cost::Float32,
                                     endtime::Integer, rate::Float32)::Float32
    (disc_cost > 0f0 && disc_rev > 0f0) || return 0f0
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

"""
One ECON harvest cost or revenue record (HRVVRCST / HRVRVN): `amount` per `unit`, applied
to harvested trees with DBH in `[dbh_lo, dbh_hi)`.
"""
struct EconCostRev
    amount::Float32
    unit::Int32
    dbh_lo::Float32
    dbh_hi::Float32
end

"""
    harvest_value(records, dbh, tpa, cuft, bdft) -> Float32

Cost or revenue from harvesting one tree record (DBH `dbh`, `tpa` trees/ac, `cuft` cubic
and `bdft` board feet per tree) against a list of `EconCostRev` (echarv.f / eccalc.f):
each record whose DBH range contains the tree contributes `amount × volume`, where the
volume is per-tree (TPA), per-MBF (BF_1000 = bdft·tpa/1000), or per-CCF
(FT3_100 = cuft·tpa/100). Log-graded units are handled by the log-bucking layer (TODO).
"""
function harvest_value(records::AbstractVector{EconCostRev}, dbh::Float32, tpa::Float32,
                       cuft::Float32, bdft::Float32)::Float32
    total = 0f0
    @inbounds for r in records
        (r.dbh_lo <= dbh < r.dbh_hi) || continue
        vol = r.unit == ECON_TPA     ? tpa :
              r.unit == ECON_BF_1000 ? bdft * tpa / 1000f0 :
              r.unit == ECON_FT3_100 ? cuft * tpa / 100f0  : 0f0
        total += r.amount * vol
    end
    return total
end
