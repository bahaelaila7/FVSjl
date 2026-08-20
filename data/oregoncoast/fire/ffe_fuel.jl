# =============================================================================
# ffe_fuel.jl (oregoncoast) — OC FFE initial fuel loading (fmcba.f), per-species COVTYP with the
# top-2-cover-type est↔init ALGSLP blend. Structurally identical to CA/NC (ca/fmcba.f), reusing the
# shared _cr_algslp2 piecewise-linear interpolator (XCOV = [10%, 60%] cover).
#
# The four DATA tables are loaded from the committed CSVs (single source of truth, parser-verified vs
# oc/fmcba.f):  fire_fuel_covtype_dead_est.csv = FUINIE (established, ~60% cover),
# fire_fuel_covtype_dead_init.csv = FUINII (initiating, ~10%), fire_fuel_covtype_live.csv = FULIVI+FULIVE.
# 50 ORGANON species × 11 dead size classes (<.25,.25-1,1-3,3-6,6-12,12-20,20-35,35-50,>50,LITT,DUFF).
# =============================================================================

const _OC_FFE_FUELDIR = normpath(joinpath(@__DIR__, ".."))

function _oc_read_fuel_csv(fname, ncol)
    rows = Vector{Vector{Float32}}()
    for (li, line) in enumerate(eachline(joinpath(_OC_FFE_FUELDIR, fname)))
        li == 1 && continue                       # header
        isempty(strip(line)) && continue
        toks = split(line, ",")
        # column 1 is species_index; take the next ncol numeric columns
        push!(rows, Float32[parse(Float32, strip(toks[1 + k])) for k in 1:ncol])
    end
    M = zeros(Float32, length(rows), ncol)
    @inbounds for i in 1:length(rows), j in 1:ncol; M[i, j] = rows[i][j]; end
    M
end

const _OC_FUINIE = _oc_read_fuel_csv("fire_fuel_covtype_dead_est.csv", 11)   # [50 × 11] established dead
const _OC_FUINII = _oc_read_fuel_csv("fire_fuel_covtype_dead_init.csv", 11)  # [50 × 11] initiating dead
# live CSV columns: herb_init, shrub_init, herb_est, shrub_est
const _OC_FULIV  = _oc_read_fuel_csv("fire_fuel_covtype_live.csv", 4)        # [50 × 4]
const _OC_FULIVI = _OC_FULIV[:, 1:2]    # initiating (herb, shrub)
const _OC_FULIVE = _OC_FULIV[:, 3:4]    # established (herb, shrub)

"oc/fmcba.f live herb/shrub fuel — top-2 cover types blended est↔init by PERCOV (ALGSLP, XCOV=[10,60])."
@inline function oc_live_fuel_loading(covca::NTuple{2,Int}, covcawt::NTuple{2,Float32}, percov::Float32)::NTuple{2,Float32}
    herb = 0f0; shrub = 0f0
    @inbounds for j in 1:2
        c = covca[j]; (1 <= c <= 50) || continue; wt = covcawt[j]
        herb  += _cr_algslp2(percov, 10f0, 60f0, _OC_FULIVI[c, 1] * wt, _OC_FULIVE[c, 1] * wt)
        shrub += _cr_algslp2(percov, 10f0, 60f0, _OC_FULIVI[c, 2] * wt, _OC_FULIVE[c, 2] * wt)
    end
    return (herb, shrub)
end

"oc/fmcba.f initial dead surface fuel (11 size classes) — top-2 cover est↔init blend by PERCOV (STFUEL)."
function oc_dead_fuel_loading(covca::NTuple{2,Int}, covcawt::NTuple{2,Float32}, percov::Float32)::Vector{Float32}
    out = zeros(Float32, 11)
    @inbounds for isz in 1:11, j in 1:2
        c = covca[j]; (1 <= c <= 50) || continue; wt = covcawt[j]
        out[isz] += _cr_algslp2(percov, 10f0, 60f0, _OC_FUINII[c, isz] * wt, _OC_FUINIE[c, isz] * wt)
    end
    return out
end

# =============================================================================
# oc_cwcalc — OC forest-grown crown width (oc/cwcalc.f CRWDTH) for FFE PERCOV. Distinct from the CCF's
# open-grown R5CRWD (oc_tree_ccf): this is the Bechtold/Crookston forest-grown model, keyed by OCMAP →
# CWEQN (FIA code + model#: 03=Crookston R1, 04=Crookston R6 m1, 05=Crookston R6 m2). Forest 711 (BLM
# Medford) uses the 610 Rogue River BF table (cwcalc.f CASE(610,710,711)). BF multiplies the leading coef
# before the cap, so it folds into `a` of _cr_r6m2. Incremental like ca_cwcalc — errors on un-ported CWEQN.
const _OC_CWMAP = ("04105","08105","24205","01703","02006","02105","20205","26305","26403","10105",
                   "10305","10805","10805","11301","11605","11705","11905","12205","12702","12702",
                   "06405","09204","21104","23104","11605","80102","80502","80702","80702","81505",
                   "81802","82102","83902","31206","31206","35106","36102","63102","35106","31206",
                   "31206","63102","63102","74605","74705","31206","98102","98102","31206","21104")

"oc/cwcalc.f Crookston-R1 crown width: a·exp(k + bcl·ln(CL) + bd·ln(D)), OMIND=1 small-tree ×D, capped."
@inline function _oc_crookr1(a::Float32, k::Float32, bcl::Float32, bd::Float32, d::Float32, cl::Float32, cap::Float32)::Float32
    dm = d < 1f0 ? 1f0 : d
    cw = a * fexp(k + bcl * flog(cl) + bd * flog(dm))
    d < 1f0 && (cw *= d)
    cw > cap ? cap : cw
end

function oc_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32)::Float32
    (1 <= sp <= 50) || return 0f0
    eqn = _OC_CWMAP[sp]
    cl = cr * h * 0.01f0
    ba1 = (barea < 1f0 ? 1f0 : barea) + 1f0     # cwcalc.f:859 BAREA=max(BA,1); formula uses (BAREA+1)
    # Crookston-R6 model 2 (BF folded into the leading coef; BF from the 610 Rogue River table).
    if     eqn == "20205"; return _cr_r6m2(6.0227f0*1.000f0, 0.54361f0,-0.20669f0,0.20395f0,-0.00644f0,-0.00378f0, d,h,cl,ba1,el,  1f0,75f0,80f0)  # DF (BF 202=1.0)
    elseif eqn == "10805"; return _cr_r6m2(6.6941f0*0.944f0, 0.81980f0,-0.36992f0,0.17722f0,-0.01202f0,-0.00882f0, d,h,cl,ba1,el,  1f0,79f0,40f0)  # LP (BF 108=0.944)
    elseif eqn == "11705"; return _cr_r6m2(3.5930f0*1.048f0, 0.63503f0,-0.22766f0,0.17827f0, 0.04267f0,-0.00290f0, d,h,cl,ba1,el,  5f0,75f0,56f0)  # SP (BF 117=1.048)
    elseif eqn == "12205"; return _cr_r6m2(4.7762f0*0.918f0, 0.74126f0,-0.28734f0,0.17137f0,-0.00602f0,-0.00209f0, d,h,cl,ba1,el, 13f0,75f0,50f0)  # PP (BF 122=0.918)
    elseif eqn == "01703"; return _oc_crookr1(1.0303f0, 1.14079f0, 0.20904f0, 0.38787f0, d, cl, 40f0)   # grand fir (Crookston R1, no BF)
    elseif eqn == "09204"; return _nc_donnelly(2.8232f0, 0.66326f0, d, 38f0)   # Brewer spruce (Donnelly D-power, no BF)
    else
        error("oc_cwcalc: crown-width equation $eqn (species $sp) not yet ported — add its oc/cwcalc.f CASE.")
    end
end
