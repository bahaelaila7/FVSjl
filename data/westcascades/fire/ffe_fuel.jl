# =============================================================================
# data/westcascades/fire/ffe_fuel.jl — WC (WestCascades) FFE initial surface-fuel loading.
#
# Ported from wc/fmcba.f DATA FULIVE/FULIVI/FUINIE/FUINII (39 species). UNLIKE the California-CWHR
# variants (NC/WS/CA) which interpolate the top-two cover species, WC uses the SINGLE dominant cover
# type COVTYP (wc/fmcba.f:476-480, 528-533) — the same structure as the interior western variants
# (CR/IE/EM/…). Live (herb, shrub) + dead (11 size classes) pools interpolated by PERCOV over
# XCOV=[10,60] between the initiating (10% cover) and established (60% cover) tables. Values from
# J. Brown / Ottmar 2000, transcribed verbatim from wc/fmcba.f. Reuses `_cr_algslp2`.
#
# wc_cwcalc (crown width for FMCBA's PERCOV): wct01 is forest 618 (Willamette) = R6 IFOR=6. WC is NOT
# in cwcalc.f's R5CRWD gate (only SO/WS, NC IFOR≤3|5, CA/OC IFOR≤5 branch there) — so WC always uses
# the Crookston R6 CAMAP path (wc/cwcalc.f WCMAP). The forest-618 BF (bark/crown factor) multiplies the
# leading coefficient only, so we reuse `_cr_r6m2(a·BF, …)` (bit-exact: Fortran computes a·BF left-assoc
# in Float32). BF[618]: 017=0.972, 093=0.857, 108=0.903, 117=1.097, 122=1.070 (015/202 absent ⇒ BF=1).
# =============================================================================

# wc/fmcba.f DATA FULIVE — established (60% cover) (herb, shrub) per species 1..39.
const _WC_FULIVE = Float32[
    0.15 0.10;  0.15 0.10;  0.15 0.10;  0.15 0.10;  0.15 0.10;  0.00 0.00;  0.15 0.10;
    0.20 0.20;  0.20 0.20;  0.30 0.20;  0.20 0.10;  0.20 0.25;  0.20 0.25;  0.15 0.10;
    0.20 0.25;  0.20 0.20;  0.20 0.20;  0.20 0.20;  0.20 0.20;  0.15 0.10;  0.20 0.20;
    0.20 0.20;  0.20 0.20;  0.20 0.20;  0.25 0.25;  0.25 0.25;  0.25 0.25;  0.23 0.22;
    0.14 0.35;  0.20 0.20;  0.20 0.10;  0.20 0.10;  0.20 0.20;  0.20 0.20;  0.25 0.25;
    0.25 0.25;  0.25 0.25;  0.00 0.00;  0.25 0.25]
# wc/fmcba.f DATA FULIVI — initiating (10% cover) (herb, shrub) per species 1..39.
const _WC_FULIVI = Float32[
    0.30 2.00;  0.30 2.00;  0.30 2.00;  0.30 2.00;  0.30 2.00;  0.00 0.00;  0.30 2.00;
    0.40 2.00;  0.40 2.00;  0.30 2.00;  0.40 1.00;  0.25 0.10;  0.25 0.10;  0.30 2.00;
    0.25 0.10;  0.40 2.00;  0.40 2.00;  0.40 2.00;  0.40 2.00;  0.30 2.00;  0.40 2.00;
    0.40 2.00;  0.40 2.00;  0.40 2.00;  0.18 2.00;  0.18 1.32;  0.18 1.32;  0.55 0.35;
    0.10 2.06;  0.40 2.00;  0.40 1.00;  0.40 1.00;  0.40 2.00;  0.40 2.00;  0.18 1.32;
    0.18 1.32;  0.18 1.32;  0.00 0.00;  0.18 1.32]

# wc/fmcba.f DATA FUINIE — established (60% cover), 11 size classes
# (<.25, .25-1, 1-3, 3-6, 6-12, 12-20, 20-35, 35-50, >50, Litter, Duff) per species 1..39.
const _WC_FUINIE = Float32[
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 1 Pacific silver fir
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 2 white fir
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 3 grand fir
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 4 subalpine fir
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 5 California/Shasta red fir
    0.0 0.0 0.0  0.0  0.0  0.0 0.0 0.0 0.0 0.0  0.0;  # 6 ---
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 7 noble fir
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 8 Alaska cedar/western larch
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 9 incense-cedar
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 10 Engelmann/Sitka spruce
    0.9 0.9 1.2  7.0  8.0  0.2 0.0 0.0 0.0 0.6 30.0;  # 11 lodgepole pine
    0.7 0.7 1.6  2.5  2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 12 Jeffrey pine
    0.7 0.7 1.6  2.5  2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 13 sugar pine
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 14 western white pine
    0.7 0.7 1.6  2.5  2.5  0.0 0.0 0.0 0.0 1.4  5.0;  # 15 ponderosa pine
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 16 Douglas-fir
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 17 coast redwood
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 18 western redcedar
    0.7 0.7 3.0  7.0  7.0 10.0 0.0 0.0 0.0 1.0 35.0;  # 19 western hemlock
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 20 mountain hemlock
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 21 bigleaf maple
    0.7 0.7 1.6  2.5  2.5  5.0 0.0 0.0 0.0 0.8 30.0;  # 22 red alder
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 23 white alder/madrone
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 24 paper birch
    0.7 0.7 0.8  1.2  1.2  0.5 0.0 0.0 0.0 1.4  0.0;  # 25 giant chinkapin/tanoak
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 26 quaking aspen
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 27 black cottonwood
    0.7 0.7 0.8  1.2  1.2  0.5 0.0 0.0 0.0 1.4  0.0;  # 28 Or white/Cal black oak
    0.1 0.2 0.4  0.5  0.8  1.0 0.0 0.0 0.0 0.1  0.0;  # 29 juniper
    0.9 0.9 1.6  3.5  3.5  0.0 0.0 0.0 0.0 0.6 10.0;  # 30 subalpine larch
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 31 whitebark pine
    0.9 0.9 1.2  7.0  8.0  0.2 0.0 0.0 0.0 0.6 30.0;  # 32 knobcone pine
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 33 Pacific yew
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 34 Pacific dogwood
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 35 hawthorn
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 36 bitter cherry
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;  # 37 willow
    0.0 0.0 0.0  0.0  0.0  0.0 0.0 0.0 0.0 0.0  0.0;  # 38 ---
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8]  # 39 other
const _WC_FUINII = Float32[    # initiating (10% cover)
    0.7 0.7 1.6 4.0 4.0 6.0 0.0 0.0 0.0 0.3 12.0;  # 1 Pacific silver fir
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 2 white fir
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 3 grand fir
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 4 subalpine fir
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 5 California/Shasta red fir
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  0.0;  # 6 ---
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 7 noble fir
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 8 Alaska cedar/western larch
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 9 incense-cedar
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 10 Engelmann/Sitka spruce
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3 12.0;  # 11 lodgepole pine
    0.1 0.1 0.2 0.5 0.5 0.0 0.0 0.0 0.0 0.5  0.8;  # 12 Jeffrey pine
    0.1 0.1 0.2 0.5 0.5 0.0 0.0 0.0 0.0 0.5  0.8;  # 13 sugar pine
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 14 western white pine
    0.1 0.1 0.2 0.5 0.5 0.0 0.0 0.0 0.0 0.5  0.8;  # 15 ponderosa pine
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 16 Douglas-fir
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 17 coast redwood
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 18 western redcedar
    0.5 0.5 2.0 2.8 2.8 6.0 0.0 0.0 0.0 0.5 12.0;  # 19 western hemlock
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 20 mountain hemlock
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 21 bigleaf maple
    0.1 0.1 0.2 0.5 0.5 3.0 0.0 0.0 0.0 0.4 12.0;  # 22 red alder
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 23 white alder/madrone
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 24 paper birch
    0.1 0.1 0.1 0.2 0.2 0.0 0.0 0.0 0.0 0.5  0.0;  # 25 giant chinkapin/tanoak
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 26 quaking aspen
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 27 black cottonwood
    0.1 0.1 0.1 0.2 0.2 0.0 0.0 0.0 0.0 0.5  0.0;  # 28 Or white/Cal black oak
    0.2 0.4 0.2 0.0 0.0 0.0 0.0 0.0 0.0 0.2  0.0;  # 29 juniper
    0.5 0.5 1.0 1.4 1.4 0.0 0.0 0.0 0.0 0.3  5.0;  # 30 subalpine larch
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 31 whitebark pine
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3 12.0;  # 32 knobcone pine
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 33 Pacific yew
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 34 Pacific dogwood
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 35 hawthorn
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 36 bitter cherry
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6;  # 37 willow
    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  0.0;  # 38 ---
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8  5.6]  # 39 other

# WC live herb/shrub fuel (wc/fmcba.f:476-480): interpolate INITIATING↔ESTABLISHED by PERCOV, keyed by
# the SINGLE dominant-species cover type COVTYP (NOT the NC/WS/CA top-two).
@inline function wc_live_fuel_loading(covtyp::Int, percov::Float32)::NTuple{2,Float32}
    (covtyp < 1 || covtyp > 39) && (covtyp = 16)     # DF default (wc/fmcba.f:462)
    herb  = _cr_algslp2(percov, 10f0, 60f0, _WC_FULIVI[covtyp, 1], _WC_FULIVE[covtyp, 1])
    shrub = _cr_algslp2(percov, 10f0, 60f0, _WC_FULIVI[covtyp, 2], _WC_FULIVE[covtyp, 2])
    return (herb, shrub)
end

# WC dead surface-fuel loading (wc/fmcba.f:528-533): interpolate FUINII↔FUINIE by PERCOV per size class.
function wc_dead_fuel_loading(covtyp::Int, percov::Float32)::Vector{Float32}
    (covtyp < 1 || covtyp > 39) && (covtyp = 16)
    out = Vector{Float32}(undef, 11)
    @inbounds for isz in 1:11
        out[isz] = _cr_algslp2(percov, 10f0, 60f0, _WC_FUINII[covtyp, isz], _WC_FUINIE[covtyp, isz])
    end
    return out
end

# =============================================================================
# wc_cwcalc — WC crown width (ft) for FMCBA's PERCOV. Crookston R6 model-2 (wc/cwcalc.f WCMAP), forest 618.
# The 5-char WCMAP code = FIA(3) + eqn#(2); eqn 05 = Crookston R6 model-2 (a·BF·D^b·H^c·CL^dd·(BA+1)^e·
# exp(EL)^f), eqn 03 = the incense-cedar-family log form. Only the 6 wct01 species are ported; others error
# (a follow-up crown-width chunk, mirroring the CA F4a scope). BF folded into the leading coefficient (a·BF).
const _WC_CWMAP = ("01105","01505","01703","01905","02006","09805","02206","04205","08105","09305",
                   "10805","11605","11705","11905","12205","20205","21104","24205","26305","26403",
                   "31206","35106","31206","37506","63102","74605","74705","81505","06405","07204",
                   "10105","10305","23104","35106","35106","35106","31206","12205","12205")

function wc_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32)::Float32
    (1 <= sp <= 39) || return 0f0
    eqn = _WC_CWMAP[sp]; cl = cr * h * 0.01f0; ba1 = barea + 1f0
    if     eqn == "01505"; return _cr_r6m2(5.0312f0*1.0f0,   0.53680f0,-0.18957f0,0.16199f0, 0.04385f0,-0.00651f0, d,h,cl,ba1,el, 2f0,75f0,35f0)   # WF 015 BF=1
    elseif eqn == "10805"; return _cr_r6m2(6.6941f0*0.903f0, 0.81980f0,-0.36992f0,0.17722f0,-0.01202f0,-0.00882f0, d,h,cl,ba1,el, 1f0,79f0,40f0)   # LP 108 BF=0.903
    elseif eqn == "11705"; return _cr_r6m2(3.5930f0*1.097f0, 0.63503f0,-0.22766f0,0.17827f0, 0.04267f0,-0.00290f0, d,h,cl,ba1,el, 5f0,75f0,56f0)   # SP 117 BF=1.097
    elseif eqn == "12205"; return _cr_r6m2(4.7762f0*1.070f0, 0.74126f0,-0.28734f0,0.17137f0,-0.00602f0,-0.00209f0, d,h,cl,ba1,el, 13f0,75f0,50f0)  # PP 122 BF=1.070
    elseif eqn == "20205"; return _cr_r6m2(6.0227f0*1.0f0,   0.54361f0,-0.20669f0,0.20395f0,-0.00644f0,-0.00378f0, d,h,cl,ba1,el, 1f0,75f0,80f0)   # DF 202 BF=1
    elseif eqn == "09305"; return _cr_r6m2(6.7575f0*0.857f0, 0.55048f0,-0.25204f0,0.19002f0, 0f0,     -0.00313f0, d,h,cl,ba1,el, 1f0,85f0,40f0)   # ES 093 BF=0.857
    elseif eqn == "01703"                                                                          # GF 017 log form (no BF, no H/BA/EL)
        dm = d >= 1f0 ? d : 1f0
        v = 1.0303f0 * fexp(1.14079f0 + 0.20904f0*flog(cl) + 0.38787f0*flog(dm))
        d < 1f0 && (v *= d); v > 40f0 && (v = 40f0); return v
    else
        error("wc_cwcalc: crown-width equation $(eqn) (WC species $(sp)) not yet ported — wct01 exercises " *
              "only WF/GF/LP/SP/PP/DF/ES; the remaining WCMAP equations are a follow-up crown-width chunk.")
    end
end
