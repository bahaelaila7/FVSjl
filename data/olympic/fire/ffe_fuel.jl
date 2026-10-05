# =============================================================================
# ffe_fuel.jl (olympic) — OP FFE initial fuel loading (op/fmcba.f) + op_cwcalc crown width for PERCOV.
# Mirror of oregoncoast/fire/ffe_fuel.jl for Olympic's 39 NWO species. Fuel tables from the committed OP
# CSVs (fire_fuel_covtype_*); op/cwcalc.f OPMAP crown width, forest 708→606 Mt Hood BF.
# =============================================================================

const _OP_FFE_FUELDIR = normpath(joinpath(@__DIR__, ".."))

function _op_read_fuel_csv(fname, ncol)
    rows = Vector{Vector{Float32}}()
    for (li, line) in enumerate(eachline(joinpath(_OP_FFE_FUELDIR, fname)))
        li == 1 && continue
        isempty(strip(line)) && continue
        toks = split(line, ",")
        push!(rows, Float32[parse(Float32, strip(toks[1 + k])) for k in 1:ncol])
    end
    M = zeros(Float32, length(rows), ncol)
    @inbounds for i in 1:length(rows), j in 1:ncol; M[i, j] = rows[i][j]; end
    M
end

const _OP_FUINIE = _op_read_fuel_csv("fire_fuel_covtype_dead_est.csv", 11)
const _OP_FUINII = _op_read_fuel_csv("fire_fuel_covtype_dead_init.csv", 11)
const _OP_FULIV  = _op_read_fuel_csv("fire_fuel_covtype_live.csv", 4)
const _OP_FULIVI = _OP_FULIV[:, 1:2]
const _OP_FULIVE = _OP_FULIV[:, 3:4]

"op/fmcba.f live herb/shrub fuel — top-2 cover est↔init blend by PERCOV (ALGSLP, XCOV=[10,60])."
@inline function op_live_fuel_loading(covca::NTuple{2,Int}, covcawt::NTuple{2,Float32}, percov::Float32)::NTuple{2,Float32}
    herb = 0f0; shrub = 0f0
    @inbounds for j in 1:2
        c = covca[j]; (1 <= c <= 39) || continue; wt = covcawt[j]
        herb  += _cr_algslp2(percov, 10f0, 60f0, _OP_FULIVI[c, 1] * wt, _OP_FULIVE[c, 1] * wt)
        shrub += _cr_algslp2(percov, 10f0, 60f0, _OP_FULIVI[c, 2] * wt, _OP_FULIVE[c, 2] * wt)
    end
    return (herb, shrub)
end

"op/fmcba.f initial dead surface fuel (11 size classes) — top-2 cover est↔init blend by PERCOV."
function op_dead_fuel_loading(covca::NTuple{2,Int}, covcawt::NTuple{2,Float32}, percov::Float32)::Vector{Float32}
    out = zeros(Float32, 11)
    @inbounds for isz in 1:11, j in 1:2
        c = covca[j]; (1 <= c <= 39) || continue; wt = covcawt[j]
        out[isz] += _cr_algslp2(percov, 10f0, 60f0, _OP_FUINII[c, isz] * wt, _OP_FUINIE[c, isz] * wt)
    end
    return out
end

# op/cwcalc.f OPMAP (39 NWO species → CWEQN = FIA + model#). Forest 708 (BLM Salem) → 606 Mt Hood BF.
const _OP_CWMAP = ("01105","01505","01703","01905","02006","09805","02206","04205","08105","09305",
                   "10805","11605","11705","11905","12205","20205","21104","24205","26305","26403",
                   "31206","35106","36102","63102","63102","74605","74705","81505","06405","07204",
                   "10105","10305","23104","35106","35106","35106","31206","12205","12205")

# op/cwcalc.f = the shared western cwcalc.f library: `_cwcalc_national` with the KODFOR Region-6 BF (cwcalc.f CASE(609,800)
# Olympic, CASE(606,708) Mt Hood/BLM Salem, …). The former per-equation copy hard-wired the 606 table and errored on the
# other OPMAP codes (26305 WH, 35106 RA, 24205 RC, 01105 SF — jl CRASH on the tiered OP fixture).
function op_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32;
                   kodfor::Int = 708)::Float32
    (1 <= sp <= 39) || return 0f0
    eqn = _OP_CWMAP[sp]
    bf = (601 <= kodfor < 1000) ? get(_R6_CWBF, (kodfor, eqn[1:3]), 1f0) : 1f0
    return _cwcalc_national(eqn, d, h, cr, barea, el, hi; bf = bf)
end
