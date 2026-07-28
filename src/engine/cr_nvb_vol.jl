# =============================================================================
# cr_nvb_vol.jl — CR volume, NVB method (NVEL nsvb.f, National-Scale Volume & Biomass).
#
# Total cubic (TCF) = Vtotib (total stem INSIDE-bark cubic), from Table S1 coefficients
# (data/centralrockies/nvb/s1_ib.csv) via CalcVOLWT eqn forms 1-5, looked up by (SPCD from
# VOLEQ(8:10), DIVISION from the eco-province). VALIDATED per-tree vs live (crt01 WF sp5:
# jl vib == live tcf 3.706/1.854/18.389/3.027, division 0).
#
# Merch cubic (MCF) = Vtotib · Rmrch, Rmrch = (1-(1-h1/H)^a)^b (Table S5 ratio, eqn 6), h1 =
# height to the merch top diameter — TODO (needs the NVBC taper HTmrch). SCF/BF also TODO.
# For now TCF is exact; MCF/SCF/BF return 0 (partial — improves the .sum TCuFt column).
# =============================================================================

const _NVB_S1 = Ref{Union{Nothing,Dict{Tuple{Int,Int},NTuple{10,Float32}}}}(nothing)

function _nvb_load_s1()
    _NVB_S1[] !== nothing && return _NVB_S1[]
    d = Dict{Tuple{Int,Int},NTuple{10,Float32}}()
    path = joinpath(@__DIR__, "..", "..", "data", "centralrockies", "nvb", "s1_ib.csv")
    for (li, line) in enumerate(eachline(path))
        li == 1 && continue                                   # header
        f = split(line, ',')
        length(f) < 13 && continue
        spcd = parse(Int, strip(f[1])); div = parse(Int, strip(f[2])); eqn = round(Int, parse(Float32, f[4]))
        coefs = (Float32(eqn), (parse(Float32, f[i]) for i in 5:13)...)  # eqn, a,a0,a1,b,b0,b1,b2,c,c1
        d[(spcd, div)] = coefs
    end
    _NVB_S1[] = d
    return d
end

"CalcVOLWT (nsvb.f:862) — volume/weight from the equation form + coefficients."
@inline function _nvb_calcvolwt(spcd::Int, d::Float32, h::Float32, eqn::Int,
                                a, a0, a1, b, b0, b1, b2, c, c1, wdsg::Float32)::Float32
    eqn <= 0 && return 0f0
    if eqn == 1
        return a * fpow(d, b) * fpow(h, c)
    elseif eqn == 2
        k = spcd < 300 ? 9f0 : 11f0
        return d < k ? a0 * fpow(d, b0) * fpow(h, c) : a0 * fpow(k, b0 - b1) * fpow(d, b1) * fpow(h, c)
    elseif eqn == 3
        return a * fpow(d, a1 * fpow(1f0 - fexp(-b1 * d), c1)) * fpow(h, c)
    elseif eqn == 4
        return a * fpow(d, b) * fpow(h, c) * fexp(-(b2 * d))
    elseif eqn == 5
        return a * fpow(d, b) * fpow(h, c) * wdsg / 62.4f0
    end
    return 0f0
end

"CR NVB total inside-bark cubic (Vtotib). `voleq`=NVB eq id, `division`=eco-province division (0 default).
Returns a 15-vec with VOL[1]=Vtotib=TCF. MCF/SCF/BF TODO."
function cr_nvb_vol(voleq::AbstractString, d::Float32, h::Float32; division::Int = 0, wdsg::Float32 = 0f0)
    vol = zeros(Float32, 15)
    (d < 1f0 || h < 5f0) && return vol
    spcd = tryparse(Int, voleq[8:10]); spcd === nothing && return vol
    tbl = _nvb_load_s1()
    row = get(tbl, (spcd, division), get(tbl, (spcd, 0), nothing))
    row === nothing && return vol
    eqn = round(Int, row[1])
    vib = _nvb_calcvolwt(spcd, d, h, eqn, row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], wdsg)
    vol[1] = vib
    return vol
end
