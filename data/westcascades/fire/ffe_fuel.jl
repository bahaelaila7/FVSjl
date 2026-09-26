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
# wc_cwcalc — WC crown width (ft), the one CRWDTH (FMCBA PERCOV, FVS_TreeList/CutList CrWidth, …). wc/cwcalc.f
# CASE('WC') CWEQN=WCMAP(ISPC) evaluated by the shared national library `_cwcalc_national`; WC is outside the
# R5CRWD gate. BF = the Region-6 forest bias factor (cwcalc.f:477-876, `_R6_CWBF` by KODFOR×FIASP) for KODFOR in
# 601..999 — KODFOR is wc/forkod.f's JFOR(IFOR) (613 → 605). jl had hard-coded forest 618 and 7 species; the other
# 32 WCMAP species errored ("02206 not yet ported"), so any TreeList on a stand with e.g. noble fir crashed.
const _WC_CWMAP = ("01105","01505","01703","01905","02006","09805","02206","04205","08105","09305",
                   "10805","11605","11705","11905","12205","20205","21104","24205","26305","26403",
                   "31206","35106","31206","37506","63102","74605","74705","81505","06405","07204",
                   "10105","10305","23104","35106","35106","35106","31206","12205","12205")

function wc_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32;
                   kodfor::Int = 618)::Float32
    (1 <= sp <= 39) || return 0f0
    eq = _WC_CWMAP[sp]
    bf = (601 <= kodfor < 1000) ? get(_R6_CWBF, (kodfor, eq[1:3]), 1f0) : 1f0
    return _cwcalc_national(eq, d, h, cr, barea, el, hi; bf = bf)
end
