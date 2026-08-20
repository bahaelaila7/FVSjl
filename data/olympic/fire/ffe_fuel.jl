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

"op/cwcalc.f forest-grown crown width for FFE PERCOV. Crookston-R6 model 2 (_cr_r6m2, BF folded into the
leading coef). 606 Mt Hood BF: WF(015)=1.130, LP(108)=0.944, others not in the 606 table ⇒ BF=1.0.
Incremental — errors on un-ported CWEQN (the opt01 ref-stand species: WF/ES/LP/SP/PP/DF, all model 05)."
function op_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32)::Float32
    (1 <= sp <= 39) || return 0f0
    eqn = _OP_CWMAP[sp]
    cl = cr * h * 0.01f0
    ba1 = (barea < 1f0 ? 1f0 : barea) + 1f0
    if     eqn == "20205"; return _cr_r6m2(6.0227f0*1.000f0, 0.54361f0,-0.20669f0,0.20395f0,-0.00644f0,-0.00378f0, d,h,cl,ba1,el,  1f0,75f0,80f0)  # DF (BF 202=1.0)
    elseif eqn == "01505"; return _cr_r6m2(5.0312f0*1.130f0, 0.53680f0,-0.18957f0,0.16199f0, 0.04385f0,-0.00651f0, d,h,cl,ba1,el,  2f0,75f0,35f0)  # WF (BF 015=1.130)
    elseif eqn == "10805"; return _cr_r6m2(6.6941f0*0.944f0, 0.81980f0,-0.36992f0,0.17722f0,-0.01202f0,-0.00882f0, d,h,cl,ba1,el,  1f0,79f0,40f0)  # LP (BF 108=0.944)
    elseif eqn == "11705"; return _cr_r6m2(3.5930f0*1.000f0, 0.63503f0,-0.22766f0,0.17827f0, 0.04267f0,-0.00290f0, d,h,cl,ba1,el,  5f0,75f0,56f0)  # SP (BF 117=1.0)
    elseif eqn == "12205"; return _cr_r6m2(4.7762f0*1.000f0, 0.74126f0,-0.28734f0,0.17137f0,-0.00602f0,-0.00209f0, d,h,cl,ba1,el, 13f0,75f0,50f0)  # PP (BF 122=1.0)
    elseif eqn == "09305"; return _cr_r6m2(6.7575f0*1.000f0, 0.55048f0,-0.25204f0,0.19002f0, 0.0f0,    -0.00313f0, d,h,cl,ba1,el,  1f0,85f0,40f0)  # ES (no BAREA term; BF 093=1.0)
    else
        error("op_cwcalc: crown-width equation $eqn (species $sp) not yet ported — add its op/cwcalc.f CASE.")
    end
end
