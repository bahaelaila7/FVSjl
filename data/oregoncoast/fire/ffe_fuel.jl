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
# CWEQN (FIA code + model#: 03=Crookston R1, 04=Crookston R6 m1, 05=Crookston R6 m2), evaluated by the shared
# `_cwcalc_national` library with the KODFOR Region-6 BF (see oc_cwcalc below).
const _OC_CWMAP = ("04105","08105","24205","01703","02006","02105","20205","26305","26403","10105",
                   "10305","10805","10805","11301","11605","11705","11905","12205","12702","12702",
                   "06405","09204","21104","23104","11605","80102","80502","80702","80702","81505",
                   "81802","82102","83902","31206","31206","35106","36102","63102","35106","31206",
                   "31206","63102","63102","74605","74705","31206","98102","98102","31206","21104")


# oc/cwcalc.f IS the shared western cwcalc.f library (VIE $Id): one SELECT CASE(CWEQN) — `_cwcalc_national` — with the
# Region-6 forest bias factor BF chosen by the post-FORKOD KODFOR (cwcalc.f:468-876; CASE(610,710,711) Rogue River,
# CASE(611,712) Siskiyou, …) on the leading coefficient of the Crookston-R6 model-2 ('…05') forms. The former
# per-equation copy here hard-wired the 610 table and errored on every other OCMAP code (36102 CO, 63102 TO, 11605 KP,
# 80502 CY, 81802 BL, 04105 PC — jl CRASH on 34 of the 100 tiered OC cases).
function oc_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32;
                   kodfor::Int = 610)::Float32
    (1 <= sp <= 50) || return 0f0
    eqn = _OC_CWMAP[sp]
    bf = (601 <= kodfor < 1000) ? get(_R6_CWBF, (kodfor, eqn[1:3]), 1f0) : 1f0
    return _cwcalc_national(eqn, d, h, cr, barea, el, hi; bf = bf)
end
