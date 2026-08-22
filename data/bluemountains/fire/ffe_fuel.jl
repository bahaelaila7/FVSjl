# BM FFE surface-fuel loading (bm/fmcba.f FULIVE/FULIVI/FUINIE/FUINII), crown-biomass species map
# (bm/fmcrow.f ISPMAP), and a crown-width remap (bm/cwcalc.f BMMAP → the CR species carrying the same
# crown-width equation, so the western cr_cwcalc library serves BM). Reuses _cr_algslp2.
#
# Live/dead loading mirrors CI: FLIVE(I)=ALGSLP(PERCOV,[10,60],[FULIVI,FULIVE]); STFUEL(ISZ,2)=
# ALGSLP(PERCOV,[10,60],[FUINII,FUINIE]) (bm/fmcba.f:319-382). MAXSP=18; covtyp>18 clamps to DF(3).

const _BM_FULIVE = Float32[
    0.15 0.1;  # 1 WP
    0.2  0.2;  # 2 WL
    0.2  0.2;  # 3 DF
    0.15 0.1;  # 4 GF
    0.2  0.2;  # 5 MH
    0.14 0.35; # 6 WJ
    0.2  0.1;  # 7 LP
    0.15 0.2;  # 8 ES
    0.15 0.2;  # 9 AF
    0.2  0.25; # 10 PP
    0.2  0.1;  # 11 WB
    0.2  0.1;  # 12 LM
    0.2  0.2;  # 13 PY
    0.2  0.2;  # 14 YC
    0.25 0.25; # 15 AS
    0.25 0.25; # 16 CW
    0.2  0.25; # 17 OS
    0.25 0.25; # 18 OH
]

const _BM_FULIVI = Float32[
    0.30 2.0;   # 1 WP
    0.4  2.0;   # 2 WL
    0.4  2.0;   # 3 DF
    0.30 2.0;   # 4 GF
    0.4  2.0;   # 5 MH
    0.10 2.06;  # 6 WJ
    0.4  1.0;   # 7 LP
    0.30 2.0;   # 8 ES
    0.30 2.0;   # 9 AF
    0.25 0.10;  # 10 PP
    0.4  1.0;   # 11 WB
    0.4  1.0;   # 12 LM
    0.4  2.0;   # 13 PY
    0.4  2.0;   # 14 YC
    0.18 1.32;  # 15 AS
    0.18 1.32;  # 16 CW
    0.25 0.10;  # 17 OS
    0.18 1.32;  # 18 OH
]

const _BM_FUINIE = Float32[
    1.0 1.0 1.6 10.0 10.0 10.0 0.0 0.0 0.0 0.8 30.0;  # 1 WP
    0.9 0.9 1.6 3.5  3.5  0.0  0.0 0.0 0.0 0.6 10.0;  # 2 WL
    0.9 0.9 1.6 3.5  3.5  0.0  0.0 0.0 0.0 0.6 10.0;  # 3 DF
    0.7 0.7 3.0 7.0  7.0  0.0  0.0 0.0 0.0 0.6 25.0;  # 4 GF
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 5 MH
    0.1 0.2 0.4 0.5  0.8  1.0  0.0 0.0 0.0 0.1 0.0;   # 6 WJ
    0.9 0.9 1.2 7.0  8.0  0.0  0.0 0.0 0.0 0.6 15.0;  # 7 LP
    1.1 1.1 2.2 10.0 10.0 0.0  0.0 0.0 0.0 0.6 30.0;  # 8 ES
    1.1 1.1 2.2 10.0 10.0 0.0  0.0 0.0 0.0 0.6 30.0;  # 9 AF
    0.7 0.7 1.6 2.5  2.5  0.0  0.0 0.0 0.0 1.4 5.0;   # 10 PP
    0.9 0.9 1.2 7.0  8.0  0.0  0.0 0.0 0.0 0.6 15.0;  # 11 WB
    0.9 0.9 1.2 7.0  8.0  0.0  0.0 0.0 0.0 0.6 15.0;  # 12 LM
    0.9 0.9 1.6 3.5  3.5  0.0  0.0 0.0 0.0 0.6 10.0;  # 13 PY
    2.2 2.2 5.2 15.0 20.0 15.0 0.0 0.0 0.0 1.0 35.0;  # 14 YC
    0.2 0.6 2.4 3.6  5.6  0.0  0.0 0.0 0.0 1.4 16.8;  # 15 AS
    0.2 0.6 2.4 3.6  5.6  0.0  0.0 0.0 0.0 1.4 16.8;  # 16 CW
    0.7 0.7 1.6 2.5  2.5  0.0  0.0 0.0 0.0 1.4 5.0;   # 17 OS
    0.2 0.6 2.4 3.6  5.6  0.0  0.0 0.0 0.0 1.4 16.8;  # 18 OH
]

const _BM_FUINII = Float32[
    0.6 0.6 0.8 6.0 6.0 6.0 0.0 0.0 0.0 0.4 12.0;  # 1 WP
    0.5 0.5 1.0 1.4 1.4 0.0 0.0 0.0 0.0 0.3 5.0;   # 2 WL
    0.5 0.5 1.0 1.4 1.4 0.0 0.0 0.0 0.0 0.3 5.0;   # 3 DF
    0.5 0.5 2.0 2.8 2.8 0.0 0.0 0.0 0.0 0.3 12.0;  # 4 GF
    1.6 1.6 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 5 MH
    0.2 0.4 0.2 0.0 0.0 0.0 0.0 0.0 0.0 0.2 0.0;   # 6 WJ
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3 7.0;   # 7 LP
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 8 ES
    0.7 0.7 1.6 4.0 4.0 0.0 0.0 0.0 0.0 0.3 12.0;  # 9 AF
    0.1 0.1 0.2 0.5 0.5 0.0 0.0 0.0 0.0 0.5 0.8;   # 10 PP
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3 7.0;   # 11 WB
    0.6 0.7 0.8 2.8 3.2 0.0 0.0 0.0 0.0 0.3 7.0;   # 12 LM
    0.5 0.5 1.0 1.4 1.4 0.0 0.0 0.0 0.0 0.3 5.0;   # 13 PY
    1.1 1.1 3.6 6.0 8.0 6.0 0.0 0.0 0.0 0.5 12.0;  # 14 YC
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8 5.6;   # 15 AS
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8 5.6;   # 16 CW
    0.1 0.1 0.2 0.5 0.5 0.0 0.0 0.0 0.0 0.5 0.8;   # 17 OS
    0.1 0.4 5.0 2.2 2.3 0.0 0.0 0.0 0.0 0.8 5.6;   # 18 OH
]

@inline function bm_live_fuel_loading(covtyp::Int, percov::Float32)::NTuple{2,Float32}
    (covtyp < 1 || covtyp > 18) && (covtyp = 3)
    (_cr_algslp2(percov, 10f0, 60f0, _BM_FULIVI[covtyp, 1], _BM_FULIVE[covtyp, 1]),
     _cr_algslp2(percov, 10f0, 60f0, _BM_FULIVI[covtyp, 2], _BM_FULIVE[covtyp, 2]))
end

function bm_dead_fuel_loading(covtyp::Int, percov::Float32)::Vector{Float32}
    (covtyp < 1 || covtyp > 18) && (covtyp = 3)
    o = Vector{Float32}(undef, 11)
    @inbounds for isz in 1:11
        o[isz] = _cr_algslp2(percov, 10f0, 60f0, _BM_FUINII[covtyp, isz], _BM_FUINIE[covtyp, isz])
    end
    o
end

# bm/fmcrow.f ISPMAP — FFE crown-biomass group (FMCROWE eastern index / FMCROWW western dispatch).
# Species 15/16/18 (aspen/cottonwood/other-hwd) use FMCROWE (Jenkins); all others use FMCROWW (= CR path).
const _BM_ISPMAP = Int[15, 8, 3, 4, 24, 16, 11, 18, 1, 13, 14, 11, 7, 8, 41, 17, 13, 41]
@inline bm_uses_fmcrowe(spiw::Integer) = spiw == 15 || spiw == 16 || spiw == 18

# bm/cwcalc.f BMMAP crown-width equation per species → the CR species (cr_cwcalc _CR_CWMAP) that carries
# the SAME 5-char crown-width equation, so the western Crookston/Bechtold library (cr_cwcalc) reproduces
# BM's CRWDTH. bmt01's species (DF/GF/WL/LP/ES/MH) all map exactly. The BM-unique equations (06405 WJ,
# 23104 PY, 04205 YC, 74705 CW, 31206 OH) have no CR carrier => nearest-genus fallback (cornered for
# non-bmt01 stands).
const _BM_TO_CR_CWSP = Int[15, 8, 3, 4, 6, 16, 11, 17, 1, 13, 14, 10, 7, 7, 19, 19, 13, 19]
# Region-6 forest bias factor BF (cwcalc.f CASE(614) UMATILLA, bmt01's forest) — per-FIASP, applied ONLY on the
# R6-Model-2 (·BF·) eqns, NOT the log-form ones (07303 WL / 01703 GF have no BF in cwcalc.f). WP(119)=1.128,
# DF(202)=1.055, LP(108)=1.244, ES(093)=1.137, AF(019)=1.110, PP(122)=1.035; WL/GF and the rest = 1.0 (log-form /
# not in the 614 table). Folded into cr_cwcalc's leading coef via its `bf` kwarg ⇒ FVS_TreeList CrWidth bit-exact.
const _BM_CWBF = Float32[1.128, 1.0, 1.055, 1.0, 1.0, 1.0, 1.244, 1.137, 1.110, 1.035, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# `forest_bf` = apply the R6 forest BF (TreeList/CutList forest-grown CRWDTH, IWHO=0). The FFE PERCOV path (fmcba,
# fuel_model) and StrClass call with the default forest_bf=false ⇒ BF=1, matching the oracle FFE crown-biomass
# (bmt01_fire is BF-FREE — measured: with-BF regressed 2030 BA 131→165). Only _forest_crwdth opts in.
@inline bm_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32, el::Float32, hi::Float32;
                  forest_bf::Bool = false)::Float32 =
    (1 <= sp <= 18) ? cr_cwcalc(_BM_TO_CR_CWSP[sp], d, h, cr, barea, el, hi;
                                bf = forest_bf ? _BM_CWBF[sp] : 1f0) : 0f0
