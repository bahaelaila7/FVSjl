# =============================================================================
# data/eastcascades/fire/ffe_fuel.jl — EC (EastCascades) FFE initial surface-fuel loading + crown width.
#
# Ported from ec/fmcba.f DATA FULIVE/FULIVI/FUINIE/FUINII (32 species) + ec/cwcalc.f (ECMAP). Like WC/PN and
# the interior western variants, EC uses the SINGLE dominant cover type COVTYP (ec/fmcba.f) — NOT the top-2
# COVCA interpolation of the California-CWHR variants (NC/WS/CA). Live (herb, shrub) + dead (11 size classes)
# pools interpolated by PERCOV over XCOV=[10,60] between the initiating (10%) and established (60%) tables.
# Transcribed VERBATIM from ec/fmcba.f. Reuses `_cr_algslp2` + `_cr_r6m2`.
#
# ec_cwcalc (crown width for FMCBA's PERCOV): ect01 is forest 608 (OKANOGAN, R6 IFOR≤4). EC's base cwidth.f
# fills CRWDTH via CWCALC (cwcalc.f), which for VARACD='EC' selects the Crookston R6 model-2 equations via
# ECMAP + the forest-608 BF (bark/crown factor). The 5-char ECMAP code = FIA(3)+eqn#(2); eqn 05 = Crookston
# R6 model-2, eqn 03 = the western-larch log form. Only the 6 ect01 species (DF/ES/LP/PP/SF/WL) are ported;
# the remaining ECMAP equations are a follow-up crown-width chunk (mirrors the WC/CA F4a scope).
# =============================================================================

# ec/fmcba.f DATA FULIVE — established (60% cover) (herb, shrub) per species 1..32.
const _EC_FULIVE = Float32[
    0.15 0.10;  0.20 0.20;  0.20 0.20;  0.15 0.10;  0.20 0.20;  0.15 0.10;  0.20 0.10;
    0.15 0.20;  0.15 0.20;  0.20 0.25;  0.20 0.20;  0.15 0.20;  0.20 0.20;  0.20 0.10;
    0.15 0.10;  0.15 0.10;  0.20 0.20;  0.20 0.20;  0.14 0.35;  0.20 0.20;  0.20 0.20;
    0.20 0.20;  0.20 0.20;  0.25 0.25;  0.20 0.20;  0.25 0.25;  0.25 0.25;  0.23 0.22;
    0.25 0.25;  0.25 0.25;  0.15 0.20;  0.25 0.25]
# ec/fmcba.f DATA FULIVI — initiating (10% cover) (herb, shrub) per species 1..32.
const _EC_FULIVI = Float32[
    0.30 2.00;  0.40 2.00;  0.40 2.00;  0.30 2.00;  0.40 2.00;  0.30 2.00;  0.40 1.00;
    0.30 2.00;  0.30 2.00;  0.25 0.10;  0.40 2.00;  0.30 2.00;  0.40 2.00;  0.40 1.00;
    0.30 2.00;  0.30 2.00;  0.40 2.00;  0.40 2.00;  0.10 2.06;  0.40 2.00;  0.40 2.00;
    0.40 2.00;  0.40 2.00;  0.18 2.00;  0.40 2.00;  0.18 1.32;  0.18 1.32;  0.55 0.35;
    0.18 1.32;  0.18 1.32;  0.30 2.00;  0.18 1.32]

# ec/fmcba.f DATA FUINIE — established (60% cover), 11 size classes
# (<.25, .25-1, 1-3, 3-6, 6-12, 12-20, 20-35, 35-50, >50, Litter, Duff) per species 1..32.
const _EC_FUINIE = Float32[
    1.0 1.0 1.6 10.0 10.0 10.0 0.0 0.0 0.0 0.8 30.0;   # 1 wp
    0.9 0.9 1.6  3.5  3.5  0.0 0.0 0.0 0.0 0.6 10.0;   # 2 wl
    0.9 0.9 1.6  3.5  3.5  0.0 0.0 0.0 0.0 0.6 10.0;   # 3 df
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;   # 4 sf
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;   # 5 rc
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;   # 6 gf
    0.9 0.9 1.2  7.0  8.0  0.0 0.0 0.0 0.0 0.6 15.0;   # 7 lp
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;   # 8 es
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;   # 9 af
    0.7 0.7 1.6  2.5  2.5  0.0 0.0 0.0 0.0 1.4  5.0;   # 10 pp
    0.7 0.7 3.0  7.0  7.0 10.0 0.0 0.0 0.0 1.0 35.0;   # 11 wh
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;   # 12 mh
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;   # 13 py
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;   # 14 wb
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;   # 15 nf
    0.7 0.7 3.0  7.0  7.0  0.0 0.0 0.0 0.0 0.6 25.0;   # 16 wf
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;   # 17 ll
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;   # 18 yc
    0.1 0.2 0.4  0.5  0.8  1.0 0.0 0.0 0.0 0.1  0.0;   # 19 wj
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;   # 20 bm
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;   # 21 vn
    0.7 0.7 1.6  2.5  2.5  5.0 0.0 0.0 0.0 0.8 30.0;   # 22 ra
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;   # 23 pb
    0.7 0.7 0.8  1.2  1.2  0.5 0.0 0.0 0.0 1.4  0.0;   # 24 gc
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;   # 25 dg
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;   # 26 as
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;   # 27 cw
    0.7 0.7 0.8  1.2  1.2  0.5 0.0 0.0 0.0 1.4  0.0;   # 28 wo
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;   # 29 pl
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8;   # 30 wi
    1.1 1.1 2.2 10.0 10.0  0.0 0.0 0.0 0.0 0.6 30.0;   # 31 os
    0.2 0.6 2.4  3.6  5.6  0.0 0.0 0.0 0.0 1.4 16.8]   # 32 oh

# ec/fmcba.f DATA FUINII — initiating (10% cover), 11 size classes per species 1..32.
const _EC_FUINII = Float32[
    0.6 0.6 0.8  6.0  6.0  6.0 0.0 0.0 0.0 0.4 12.0;   # 1 wp
    0.5 0.5 1.0  1.4  1.4  0.0 0.0 0.0 0.0 0.3  5.0;   # 2 wl
    0.5 0.5 1.0  1.4  1.4  0.0 0.0 0.0 0.0 0.3  5.0;   # 3 df
    0.5 0.5 2.0  2.8  2.8  0.0 0.0 0.0 0.0 0.3 12.0;   # 4 sf
    1.6 1.6 3.6  6.0  8.0  6.0 0.0 0.0 0.0 0.5 12.0;   # 5 rc
    0.5 0.5 2.0  2.8  2.8  0.0 0.0 0.0 0.0 0.3 12.0;   # 6 gf
    0.6 0.7 0.8  2.8  3.2  0.0 0.0 0.0 0.0 0.3  7.0;   # 7 lp
    0.7 0.7 1.6  4.0  4.0  0.0 0.0 0.0 0.0 0.3 12.0;   # 8 es
    0.7 0.7 1.6  4.0  4.0  0.0 0.0 0.0 0.0 0.3 12.0;   # 9 af
    0.1 0.1 0.2  0.5  0.5  0.0 0.0 0.0 0.0 0.5  0.8;   # 10 pp
    0.5 0.5 2.0  2.8  2.8  6.0 0.0 0.0 0.0 0.5 12.0;   # 11 wh
    0.7 0.7 1.6  4.0  4.0  0.0 0.0 0.0 0.0 0.3 12.0;   # 12 mh
    1.1 1.1 3.6  6.0  8.0  6.0 0.0 0.0 0.0 0.5 12.0;   # 13 py
    0.7 0.7 1.6  4.0  4.0  0.0 0.0 0.0 0.0 0.3 12.0;   # 14 wb
    0.5 0.5 2.0  2.8  2.8  0.0 0.0 0.0 0.0 0.3 12.0;   # 15 nf
    0.5 0.5 2.0  2.8  2.8  0.0 0.0 0.0 0.0 0.3 12.0;   # 16 wf
    0.5 0.5 1.0  1.4  1.4  0.0 0.0 0.0 0.0 0.3  5.0;   # 17 ll
    1.1 1.1 3.6  6.0  8.0  6.0 0.0 0.0 0.0 0.5 12.0;   # 18 yc
    0.2 0.4 0.2  0.0  0.0  0.0 0.0 0.0 0.0 0.2  0.0;   # 19 wj
    1.1 1.1 3.6  6.0  8.0  6.0 0.0 0.0 0.0 0.5 12.0;   # 20 bm
    1.1 1.1 3.6  6.0  8.0  6.0 0.0 0.0 0.0 0.5 12.0;   # 21 vn
    0.1 0.1 0.2  0.5  0.5  3.0 0.0 0.0 0.0 0.4 12.0;   # 22 ra
    1.1 1.1 3.6  6.0  8.0  6.0 0.0 0.0 0.0 0.5 12.0;   # 23 pb
    0.1 0.1 0.1  0.2  0.2  0.0 0.0 0.0 0.0 0.5  0.0;   # 24 gc
    1.1 1.1 3.6  6.0  8.0  6.0 0.0 0.0 0.0 0.5 12.0;   # 25 dg
    0.1 0.4 5.0  2.2  2.3  0.0 0.0 0.0 0.0 0.8  5.6;   # 26 as
    0.1 0.4 5.0  2.2  2.3  0.0 0.0 0.0 0.0 0.8  5.6;   # 27 cw
    0.1 0.1 0.1  0.2  0.2  0.0 0.0 0.0 0.0 0.5  0.0;   # 28 wo
    0.1 0.4 5.0  2.2  2.3  0.0 0.0 0.0 0.0 0.8  5.6;   # 29 pl
    0.1 0.4 5.0  2.2  2.3  0.0 0.0 0.0 0.0 0.8  5.6;   # 30 wi
    0.7 0.7 1.6  4.0  4.0  0.0 0.0 0.0 0.0 0.3 12.0;   # 31 os
    0.1 0.4 5.0  2.2  2.3  0.0 0.0 0.0 0.0 0.8  5.6]   # 32 oh

# EC live herb/shrub surface-fuel loading (ec/fmcba.f:443-449): interpolate FULIVI↔FULIVE by PERCOV.
@inline function ec_live_fuel_loading(covtyp::Int, percov::Float32)::NTuple{2,Float32}
    (covtyp < 1 || covtyp > 32) && (covtyp = 3)      # DF default (ec/fmcba.f:431)
    herb  = _cr_algslp2(percov, 10f0, 60f0, _EC_FULIVI[covtyp, 1], _EC_FULIVE[covtyp, 1])
    shrub = _cr_algslp2(percov, 10f0, 60f0, _EC_FULIVI[covtyp, 2], _EC_FULIVE[covtyp, 2])
    return (herb, shrub)
end

# EC dead surface-fuel loading (ec/fmcba.f:497-504): interpolate FUINII↔FUINIE by PERCOV per size class.
function ec_dead_fuel_loading(covtyp::Int, percov::Float32)::Vector{Float32}
    (covtyp < 1 || covtyp > 32) && (covtyp = 3)
    out = Vector{Float32}(undef, 11)
    @inbounds for isz in 1:11
        out[isz] = _cr_algslp2(percov, 10f0, 60f0, _EC_FUINII[covtyp, isz], _EC_FUINIE[covtyp, isz])
    end
    return out
end

# ec/fmcba.f DATA MAPDRY(155) — each EC habitat code (ITYPE) → dry(0)/mesic(1)/moist(2). Used by FMCFMD via
# the ECMOIST entry: the moist-mixed metagroup (MMIXCT) is selected iff MAPDRY==1, else dry-mixed (DMIXCT).
const _EC_MAPDRY = Int[
    0,1,2,2,1,2,2,1,0,0, 0,0,0,0,0,1,0,0,0,1, 0,0,0,0,1,1,0,0,0,0,
    0,0,0,1,0,0,0,0,0,1, 1,1,1,2,1,1,1,2,2,2,
    2,2,2,2,1,1,2,2,2,1, 1,1,1,1,1,1,1,2,2,1, 1,1,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,2,2,1,1, 2,2,2,2,2,2,2,1,2,2,
    2,2,2,2,2,2,2,2,2,0, 0,0,0,0,2,2,1,2,2,2, 2,0,0,0,0,1,2,2,2,2,
    1,1,1,0,2,0,0,0,0,0, 2,0,2,0,0,1,2,2,2,1, 2,0,2,0,0]

# ec/fmcba.f ENTRY ECMOIST — the moist/dry habitat flag for FMCFMD (1 ⇒ moist-mixed MMIXCT, else dry-mixed).
@inline ec_moist(itype::Integer)::Int = (1 <= itype <= length(_EC_MAPDRY)) ? _EC_MAPDRY[itype] : 0

# ec/cwcalc.f ECMAP — EC crown width (ft) for FMCBA's PERCOV. Crookston R6 model-2 (a·BF·D^b·H^c·CL^dd·
# (BA+1)^e·exp(EL)^f), forest-608 (OKANOGAN) BF applied to the leading coefficient. WL (073, eqn 03) uses the
# western-larch log form (no BF). Only the 6 ect01 species are ported; others error (follow-up chunk).
function ec_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32)::Float32
    (1 <= sp <= 32) || return 0f0
    cl = cr * h * 0.01f0; ba1 = barea + 1f0
    if     sp == 3   # DF 20205 BF=1
        return _cr_r6m2(6.0227f0,       0.54361f0,-0.20669f0,0.20395f0,-0.00644f0,-0.00378f0, d,h,cl,ba1,el, 1f0,75f0,80f0)
    elseif sp == 8   # ES 09305 BF=1 (093 not in forest-608 BF list); no (BA+1) term
        return _cr_r6m2(6.7575f0,       0.55048f0,-0.25204f0,0.19002f0, 0f0,      -0.00313f0, d,h,cl,ba1,el, 1f0,85f0,40f0)
    elseif sp == 7   # LP 10805 BF=1.114
        return _cr_r6m2(6.6941f0*1.114f0, 0.81980f0,-0.36992f0,0.17722f0,-0.01202f0,-0.00882f0, d,h,cl,ba1,el, 1f0,79f0,40f0)
    elseif sp == 10  # PP 12205 BF=1
        return _cr_r6m2(4.7762f0,       0.74126f0,-0.28734f0,0.17137f0,-0.00602f0,-0.00209f0, d,h,cl,ba1,el, 13f0,75f0,50f0)
    elseif sp == 4   # SF 01105 BF=1; EL floor 4
        return _cr_r6m2(4.4799f0,       0.45976f0,-0.10425f0,0.11866f0, 0.06762f0,-0.00715f0, d,h,cl,ba1,el, 4f0,72f0,33f0)
    elseif sp == 2   # WL 07303 log form (no BF), uses ALOG(BAREA) directly
        dm = d >= 1f0 ? d : 1f0
        v = 1.02478f0 * fexp(0.99889f0 + 0.19422f0*flog(cl) + 0.59423f0*flog(dm) -
                             0.09078f0*flog(h) - 0.02341f0*flog(barea))
        d < 1f0 && (v *= d); v > 40f0 && (v = 40f0); return v
    else
        error("ec_cwcalc: crown-width equation for EC species $(sp) not yet ported — ect01 exercises only " *
              "DF/ES/LP/PP/SF/WL; the remaining ECMAP equations are a follow-up crown-width chunk.")
    end
end
