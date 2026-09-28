# =============================================================================
# data/southeastalaska/fire/ffe_fuel.jl — AK (SoutheastAlaska) FFE data from fire/ak/*.f (the six AK-specific
# FFE sources; every other FFE routine FVSak links is byte-identical to FVSpn's, see the AK notes in the engine).
#
#  * ak/fmcba.f   FULIVE/FULIVI (established/initiating live herb+shrub by the dominant-BA cover species,
#                 ALGSLP over PERCOV 10..60) and FUINI (the initial dead fuel by the INVENTORY forest type,
#                 FTDEADFU 1-15; hard pool only). Bare stand ⇒ COVTYP 11 (western hemlock) the first year.
#  * ak/fmbrkt.f  FOFEM bark thickness B1 (FMBRKT = DBH·B1).
#  * ak/fmvinit.f DKR (4 decay classes × 11 sizes), HTX (1.0 except YC/RC = 0), HTR1=HTR2=0.02, TFALL.
#  * the vbase fmr6sdcy/fmr6fall AK branches live in the shared R6 layer (engine/fire/r6_snag_fall.jl, code :AK).
# =============================================================================

# ak/fmcba.f DATA FULIVE / FULIVI — (herb, shrub) per species 1..23.
const _AK_FULIVE = Float32[
    0.15 0.10; 0.15 0.10; 0.20 0.20; 0.30 0.20; 0.30 0.20; 0.30 0.20; 0.30 0.20; 0.30 0.20; 0.20 0.10; 0.20 0.20;
    0.20 0.20; 0.15 0.20; 0.30 0.20; 0.20 0.20; 0.20 0.20; 0.20 0.20; 0.20 0.20; 0.25 0.25; 0.25 0.25; 0.25 0.25;
    0.25 0.25; 0.25 0.25; 0.25 0.25]
const _AK_FULIVI = Float32[
    0.30 2.00; 0.30 2.00; 0.40 2.00; 0.30 2.00; 0.30 2.00; 0.30 2.00; 0.30 2.00; 0.30 2.00; 0.40 1.00; 0.40 2.00;
    0.40 2.00; 0.30 2.00; 0.30 2.00; 0.40 2.00; 0.40 2.00; 0.40 2.00; 0.40 2.00; 0.18 1.32; 0.18 1.32; 0.18 1.32;
    0.18 1.32; 0.18 1.32; 0.18 1.32]
# ak/fmcba.f DATA FUINI(11 size classes, 15 dead-fuel forest types).
const _AK_FUINI = Float32[
    0.06 0.47 0.82 0.60 2.66 2.72 0.0 0.0 0.0  3.75 49.75;   # 1  122 white spruce
    0.10 0.19 0.88 0.02 0.17 0.00 0.0 0.0 0.0  4.75 49.75;   # 2  125 black spruce
    0.05 0.16 0.50 0.18 0.40 1.01 0.0 0.0 0.0  1.62 55.33;   # 3  270 mountain hemlock
    0.03 0.08 0.40 0.26 1.25 3.89 0.0 0.0 0.0  1.45 76.93;   # 4  271 Alaska yellow cedar
    0.02 0.07 0.14 0.32 0.82 0.28 0.0 0.0 0.0  1.76 55.33;   # 5  281 lodgepole pine
    0.09 0.36 1.21 0.63 1.98 9.26 0.0 0.0 0.0  2.00 59.80;   # 6  301 western hemlock
    0.05 0.25 0.41 0.31 1.32 1.26 0.0 0.0 0.0  2.43 79.86;   # 7  304 western redcedar
    0.07 0.28 1.07 0.83 1.67 8.00 0.0 0.0 0.0  2.42 49.61;   # 8  305 Sitka spruce (also 264, 268)
    0.02 0.17 0.42 0.09 0.83 0.00 0.0 0.0 0.0  8.08 47.03;   # 9  703 cottonwood
    0.04 0.49 2.11 0.08 0.54 0.00 0.0 0.0 0.0  3.54 52.49;   # 10 709/704 cottonwood-willow
    0.04 0.04 0.44 0.05 0.39 0.00 0.0 0.0 0.0 16.58 67.23;   # 11 901 aspen
    0.09 0.43 1.64 0.17 1.76 0.19 0.0 0.0 0.0 13.82 61.51;   # 12 902 paper birch
    0.15 1.50 2.68 0.15 5.62 5.39 0.0 0.0 0.0  2.78 59.80;   # 13 911 red alder
    0.02 0.12 0.82 0.56 1.74 0.81 0.0 0.0 0.0  7.28 46.41;   # 14 999 nonstocked
    0.22 0.47 1.62 0.62 2.88 7.99 0.0 0.0 0.0  5.35 46.41]   # 15 others

"ak/fmcba.f live herb/shrub: ALGSLP(PERCOV, [10,60], [FULIVI, FULIVE]) of the cover species COVTYP."
@inline function ak_live_fuel_loading(covtyp::Int, percov::Float32)::NTuple{2,Float32}
    (covtyp < 1 || covtyp > 23) && (covtyp = 11)
    herb  = _cr_algslp2(percov, 10f0, 60f0, _AK_FULIVI[covtyp, 1], _AK_FULIVE[covtyp, 1])
    shrub = _cr_algslp2(percov, 10f0, 60f0, _AK_FULIVI[covtyp, 2], _AK_FULIVE[covtyp, 2])
    return (herb, shrub)
end

"ak/fmcba.f:270-319 FTDEADFU: the initial dead-fuel row by the forest type IFORTP at the inventory year."
function ak_ftdeadfu(ifortp::Integer)::Int
    ifortp == 122 && return 1
    ifortp == 125 && return 2
    ifortp == 270 && return 3
    ifortp == 271 && return 4
    ifortp == 281 && return 5
    ifortp == 301 && return 6
    ifortp == 304 && return 7
    (ifortp == 305 || ifortp == 264 || ifortp == 268) && return 8
    ifortp == 703 && return 9
    (ifortp == 709 || ifortp == 704) && return 10
    ifortp == 901 && return 11
    ifortp == 902 && return 12
    ifortp == 911 && return 13
    ifortp == 999 && return 14
    return 15
end
ak_dead_fuel_loading(ifortp::Integer)::Vector{Float32} = Float32[_AK_FUINI[ak_ftdeadfu(ifortp), k] for k in 1:11]

# ak/fmbrkt.f DATA B1 — FOFEM bark thickness (FMBRKT = DBH·B1).
const _AK_FM_BARK_B1 = Float32[0.047, 0.041, 0.022, 0.031, 0.025, 0.025, 0.032, 0.027, 0.028, 0.035, 0.040, 0.040,
                               0.025, 0.026, 0.026, 0.027, 0.027, 0.040, 0.044, 0.044, 0.041, 0.041, 0.044]

# ak/fmvinit.f DKR(size 1..11, decay class 1..4).
const _FM_DKR_AK = Float32[
    0.052 0.061 0.073 0.098;   # <0.25"
    0.052 0.061 0.073 0.098;   # 0.25-1"
    0.052 0.061 0.073 0.098;   # 1-3"
    0.012 0.025 0.041 0.077;   # 3-6"
    0.012 0.025 0.041 0.077;   # 6-12"
    0.009 0.018 0.031 0.058;   # 12-20"
    0.009 0.018 0.031 0.058;   # 20-35"
    0.009 0.018 0.031 0.058;   # 35-50"
    0.009 0.018 0.031 0.058;   # >50"
    0.35  0.4   0.45  0.5  ;   # litter
    0.002 0.002 0.003 0.003]   # duff

# ak/fmvinit.f HTX(I,1:4): 1.0 for every species except Alaska-cedar (3) and western redcedar (10) = 0.0.
const _AK_FM_HTX = NTuple{4,Float32}[(sp == 3 || sp == 10) ? (0f0, 0f0, 0f0, 0f0) : (1f0, 1f0, 1f0, 1f0) for sp in 1:23]

# ak/fmvinit.f TFALL(I,0:5) (foliage, <0.25", 0.25-1", 1-3", 3-6", 6-12"); TFALL(I,0)=MIN(TFALL(I,0),LEAFLF),
# TFALL(I,2)=MIN(TFALL(I,2),TFALL(I,3)), TFALL(I,5)=TFALL(I,4).
function _ak_tfall_row(sp::Int)
    t = (sp in (5, 6, 7, 13)) ? (2f0, 5f0, 5f0, 10f0, 50f0) :
        (sp in (3, 10))       ? (5f0, 15f0, 15f0, 30f0, 55f0) :
        (sp in (4, 11, 12))   ? (1f0, 5f0, 5f0, 15f0, 50f0) :
        (sp in (1, 2, 8, 9))  ? (2f0, 5f0, 5f0, 15f0, 50f0) :
                                (1f0, 10f0, 15f0, 15f0, 50f0)          # 14:23
    leaflf = sp == 1 || sp == 2 ? 7f0 : sp in (3, 8, 10, 11) ? 5f0 : sp == 4 ? 1f0 : sp in (5, 6, 7, 13) ? 6f0 :
             sp == 9 ? 3f0 : sp == 12 ? 4f0 : 1f0
    f0 = min(t[1], leaflf); f2 = min(t[3], t[4])
    return (f0, t[2], f2, t[4], t[5], t[5])
end
const _FM_TFALL_AK = let m = zeros(Float32, 23, 6)
    for sp in 1:23, k in 1:6
        m[sp, k] = _ak_tfall_row(sp)[k]
    end
    m
end

