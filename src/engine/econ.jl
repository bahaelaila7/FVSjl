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
