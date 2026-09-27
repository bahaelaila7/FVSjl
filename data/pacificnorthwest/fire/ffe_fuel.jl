# =============================================================================
# data/pacificnorthwest/fire/ffe_fuel.jl — PN (PacificNorthwest) FFE initial surface-fuel loading.
#
# Ported from pn/fmcba.f DATA FULIVE/FULIVI/FUINIE/FUINII (39 species). Like WC (and the interior
# western variants), PN uses the SINGLE dominant cover type COVTYP (pn/fmcba.f) — NOT the top-2 COVCA
# interpolation of the California-CWHR variants (NC/WS/CA). Live (herb, shrub) + dead (11 size classes)
# pools interpolated by PERCOV over XCOV=[10,60] between the initiating (10%) and established (60%)
# tables. Transcribed VERBATIM from pn/fmcba.f. Reuses `_cr_algslp2`.
#
# PN differs from WC (near-clone) only in coefficient DATA: slot 6 = SS (Sitka spruce; WC blank), and
# the FUINIE(sp10 ES → DF table)/FUINII(sp8/9/10/16/17/18/21 first-two classes 1.1→1.6, sp1 12-20"
# 6.0→0.0) tables carry genuine per-species differences (pn/fmcba.f vs wc/fmcba.f).
#
# pn_cwcalc (crown width for FMCBA's PERCOV): pnt01 is forest 612 (SIUSLAW) = R6 IFOR. PN is NOT in
# cwcalc.f's R5CRWD gate (only SO/NC/CA/OC variants branch there) — so PN always uses the Crookston R6
# CAMAP path. pn/cwcalc.f is BYTE-IDENTICAL to wc/cwcalc.f (md5 2bace281): SAME eqn codes + SAME base
# coefficients; only the forest-specific BF differs (forest 612 SIUSLAW: DF 202=0.977, RC 242=0.905,
# WH 263=0.924; all other species BF=1.0 — vs WC's forest-618 Willamette factors).
# =============================================================================

# pn/fmcba.f DATA FULIVE — established (60% cover) (herb, shrub) per species 1..39.
const _PN_FULIVE = Float32[
    0.15 0.10;  0.15 0.10;  0.15 0.10;  0.15 0.10;  0.15 0.10;  0.30 0.20;  0.15 0.10;
    0.20 0.20;  0.20 0.20;  0.30 0.20;  0.20 0.10;  0.20 0.25;  0.20 0.25;  0.15 0.10;
    0.20 0.25;  0.20 0.20;  0.20 0.20;  0.20 0.20;  0.20 0.20;  0.15 0.10;  0.20 0.20;
    0.20 0.20;  0.20 0.20;  0.20 0.20;  0.25 0.25;  0.25 0.25;  0.25 0.25;  0.23 0.22;
    0.14 0.35;  0.20 0.20;  0.20 0.10;  0.20 0.10;  0.20 0.20;  0.20 0.20;  0.25 0.25;
    0.25 0.25;  0.25 0.25;  0.00 0.00;  0.25 0.25]
# pn/fmcba.f DATA FULIVI — initiating (10% cover) (herb, shrub) per species 1..39.
const _PN_FULIVI = Float32[
    0.30 2.00;  0.30 2.00;  0.30 2.00;  0.30 2.00;  0.30 2.00;  0.30 2.00;  0.30 2.00;
    0.40 2.00;  0.40 2.00;  0.30 2.00;  0.40 1.00;  0.25 0.10;  0.25 0.10;  0.30 2.00;
    0.25 0.10;  0.40 2.00;  0.40 2.00;  0.40 2.00;  0.40 2.00;  0.30 2.00;  0.40 2.00;
    0.40 2.00;  0.40 2.00;  0.40 2.00;  0.18 2.00;  0.18 1.32;  0.18 1.32;  0.55 0.35;
    0.10 2.06;  0.40 2.00;  0.40 1.00;  0.40 1.00;  0.40 2.00;  0.40 2.00;  0.18 1.32;
    0.18 1.32;  0.18 1.32;  0.00 0.00;  0.18 1.32]

# pn/fmcba.f DATA FUINIE — established (60% cover), 11 size classes
# (<.25, .25-1, 1-3, 3-6, 6-12, 12-20, 20-35, 35-50, >50, Litter, Duff) per species 1..39.
const _PN_FUINIE = Float32[
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 1 Pacific silver fir
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 2 white fir
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 3 grand fir
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;  # 4 subalpine fir
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 5 California/Shasta red fir
    0.7 0.7 3.0  7.0  7.0 10.0 0.0 0.0 0.0 1.0 35.0;  # 6 Sitka spruce
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;  # 7 noble fir
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 8 Alaska cedar/western larch
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 9 incense-cedar
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 10 Engelmann spruce (PN: uses DF table)
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
# pn/fmcba.f DATA FUINII — initiating (10% cover).
const _PN_FUINII = Float32[
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 1 Pacific silver fir (PN: 12-20"=0.0)
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 2 white fir
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 3 grand fir
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 4 subalpine fir
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 5 California/Shasta red fir
    0.5 0.5 2.0 2.8 2.8 6.0 0.0 0.0 0.0 0.5 12.0;  # 6 Sitka spruce
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 7 noble fir
    1.6 1.6 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 8 Alaska cedar/western larch
    1.6 1.6 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 9 incense-cedar
    1.6 1.6 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 10 Engelmann spruce
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3 12.0;  # 11 lodgepole pine
    0.1 0.1 0.2 0.5 0.5 0.0 0.0 0.0 0.0 0.5  0.8;  # 12 Jeffrey pine
    0.1 0.1 0.2 0.5 0.5 0.0 0.0 0.0 0.0 0.5  0.8;  # 13 sugar pine
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 14 western white pine
    0.1 0.1 0.2 0.5 0.5 0.0 0.0 0.0 0.0 0.5  0.8;  # 15 ponderosa pine
    1.6 1.6 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 16 Douglas-fir
    1.6 1.6 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 17 coast redwood
    1.6 1.6 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 18 western redcedar
    0.5 0.5 2.0 2.8 2.8 6.0 0.0 0.0 0.0 0.5 12.0;  # 19 western hemlock
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 20 mountain hemlock
    1.6 1.6 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 21 bigleaf maple
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

# PN live herb/shrub fuel (pn/fmcba.f): interpolate INITIATING↔ESTABLISHED by PERCOV, keyed by the
# SINGLE dominant-species cover type COVTYP (NOT the NC/WS/CA top-two).
@inline function pn_live_fuel_loading(covtyp::Int, percov::Float32)::NTuple{2,Float32}
    (covtyp < 1 || covtyp > 39) && (covtyp = 16)     # DF default (pn/fmcba.f)
    herb  = _cr_algslp2(percov, 10f0, 60f0, _PN_FULIVI[covtyp, 1], _PN_FULIVE[covtyp, 1])
    shrub = _cr_algslp2(percov, 10f0, 60f0, _PN_FULIVI[covtyp, 2], _PN_FULIVE[covtyp, 2])
    return (herb, shrub)
end

# PN dead surface-fuel loading (pn/fmcba.f): interpolate FUINII↔FUINIE by PERCOV per size class.
function pn_dead_fuel_loading(covtyp::Int, percov::Float32)::Vector{Float32}
    (covtyp < 1 || covtyp > 39) && (covtyp = 16)
    out = Vector{Float32}(undef, 11)
    @inbounds for isz in 1:11
        out[isz] = _cr_algslp2(percov, 10f0, 60f0, _PN_FUINII[covtyp, isz], _PN_FUINIE[covtyp, isz])
    end
    return out
end

# =============================================================================
# pn_cwcalc — PN crown width (ft). pn/cwcalc.f is BYTE-IDENTICAL to wc/cwcalc.f (md5 2bace281): the same WCMAP
# codes, the national equations and the Region-6 bias factor BF by KODFOR×FIASP (cwcalc.f:188-876). So PN
# evaluates exactly wc_cwcalc with its own post-FORKOD KODFOR (609/612/800/708/709/712). The old PN port covered
# only the 7 pnt01 species at forest 612 and errored on the rest (a PN TreeList / FFE run on species 8+ crashed).
# =============================================================================
pn_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32;
          kodfor::Int = 612)::Float32 = wc_cwcalc(sp, d, h, cr, barea, el, hi; kodfor = kodfor)
