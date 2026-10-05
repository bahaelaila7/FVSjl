# fire/fuel_decay.jl — FFE down-wood / litter / duff decay (FMCWD, fire/base/fmcwd.f).
#
# Each cycle the surface-fuel pools `fire.cwd[size 1:11, soft(1)/hard(2), decay-class 1:4]` decay by
# size- and decay-class-specific annual rates (DKR), with a fraction (PRDUFF) of the decayed woody
# material routed to duff, and a hard→soft transfer for the woody classes. Soft material decays 10%
# faster than hard (the ·1.1). This is the *decay* half of the FFE fuel dynamics; the *additions*
# half (litterfall + crown breakage + snag falldown, FMCADD) is a companion routine — and the two are
# COUPLED: litter (size 10) decays at 0.65/yr (≈ gone in one 5-yr cycle), so the forest-floor pool only
# holds up because annual litterfall replenishes it. So this routine is faithful on its own but the
# grown-cycle Stand Carbon Report (DDW/Floor) only validates once FMCADD lands too — see
# docs/FFE_FUEL_DYNAMICS_chunk_plan.md.

# DKR — annual decay rate by [size class 1:11, decay class 1:4] (sn/fmvinit.f:70-104). Woody classes
# 1-9 use 0.11 (decay class 1) or 0.11/0.11/0.09/0.07… (classes 2-4, which copy class 2); litter
# (10) = 0.65, duff (11) = 0.002.
const _FM_DKR = Float32[
    0.11 0.11 0.11 0.11     # 1  (<0.25")
    0.11 0.11 0.11 0.11     # 2  (0.25-1")
    0.11 0.09 0.09 0.09     # 3  (1-3")
    0.11 0.07 0.07 0.07     # 4  (3-6")
    0.11 0.07 0.07 0.07     # 5  (6-12")
    0.11 0.07 0.07 0.07     # 6  (12-20")
    0.11 0.07 0.07 0.07     # 7  (20-35")
    0.11 0.07 0.07 0.07     # 8  (35-50")
    0.11 0.07 0.07 0.07     # 9  (>50")
    0.65 0.65 0.65 0.65     # 10 litter
    0.002 0.002 0.002 0.002 # 11 duff
]
# Lake States annual decay rates (ls/fmvinit.f:72-95). LS is decay-class-INDEPENDENT (DKR(I,J)=DKR(I,1))
# and its litter loss is 0.31/yr — NOT the SN 0.65 — so LS litter equilibrates ~2× higher, which is the
# down-wood loading FMDYN weights the fire fuel model on. Applying the SN table to LS decayed the litter ~2×
# too fast ⇒ SMALL down-wood ~1.47× low ⇒ FMDYN under-weighted the hot model ⇒ under-scorch ⇒ fire under-kill
# (S82-S86; ls_simfire 2020 TPA). Woody classes are also slower than SN (0.06/0.02 vs 0.07). Duff 0.002 == SN.
const _FM_DKR_LS = Float32[
    0.11 0.11 0.11 0.11     # 1
    0.11 0.11 0.11 0.11     # 2
    0.09 0.09 0.09 0.09     # 3
    0.06 0.06 0.06 0.06     # 4
    0.06 0.06 0.06 0.06     # 5
    0.02 0.02 0.02 0.02     # 6
    0.02 0.02 0.02 0.02     # 7
    0.02 0.02 0.02 0.02     # 8
    0.02 0.02 0.02 0.02     # 9
    0.31 0.31 0.31 0.31     # 10 litter (ls/fmvinit.f:94)
    0.002 0.002 0.002 0.002 # 11 duff
]
# Northeast annual decay rates (ne/fmvinit.f:73-96). Decay-class-INDEPENDENT; litter 0.40/yr (NOT SN 0.65);
# woody 0.19/0.19/0.11/0.07/0.03… (Fahey/Arthur/Foster-Lang). Same SN-mismatch class as LS.
const _FM_DKR_NE = Float32[
    0.19 0.19 0.19 0.19     # 1
    0.19 0.19 0.19 0.19     # 2
    0.11 0.11 0.11 0.11     # 3
    0.07 0.07 0.07 0.07     # 4
    0.03 0.03 0.03 0.03     # 5
    0.03 0.03 0.03 0.03     # 6
    0.03 0.03 0.03 0.03     # 7
    0.03 0.03 0.03 0.03     # 8
    0.03 0.03 0.03 0.03     # 9
    0.40 0.40 0.40 0.40     # 10 litter (ne/fmvinit.f:96)
    0.002 0.002 0.002 0.002 # 11 duff
]
# Central States annual decay rates (cs/fmvinit.f:70-93). Litter 0.65 (== SN) but woody is decay-class-
# INDEPENDENT: classes 3-9 use 0.09/0.07 for ALL decay classes (vs SN's 0.11 at decay class 1). Duff 0.002.
const _FM_DKR_CS = Float32[
    0.11 0.11 0.11 0.11     # 1
    0.11 0.11 0.11 0.11     # 2
    0.09 0.09 0.09 0.09     # 3
    0.07 0.07 0.07 0.07     # 4
    0.07 0.07 0.07 0.07     # 5
    0.07 0.07 0.07 0.07     # 6
    0.07 0.07 0.07 0.07     # 7
    0.07 0.07 0.07 0.07     # 8
    0.07 0.07 0.07 0.07     # 9
    0.65 0.65 0.65 0.65     # 10 litter (cs/fmvinit.f:90, == SN)
    0.002 0.002 0.002 0.002 # 11 duff
]
# Central Rockies annual decay rates (cr/fmvinit.f:80-112). Decay class 1: woody 0.12/0.12/0.09 (classes 1-3),
# 0.015 (classes 4-9 = LARGE wood, "agrees with published values for large CWD, Brown et al. 1998"); litter 0.5,
# duff 0.002. Decay classes 2-4 are the class-1 rate × 0.45 (the CR/UT modifier, fmvinit.f:108-112). The SN
# default decayed CR's LARGE down-wood ~7× too fast (0.11 vs 0.015 at decay class 1) ⇒ the down-wood pool
# collapsed over cycles (crt01 large 16.7→1.6) ⇒ FMDYN under-weighted the heavy model on non-MCCT fire stands
# (SFCT/WSCT/ASCT model 8-vs-10) + the crt01 byram residual. This is the CR analogue of the LS DKR fix above.
const _FM_DKR_CR = Float32[
    0.12   0.054   0.054   0.054      # 1  (<0.25")
    0.12   0.054   0.054   0.054      # 2  (0.25-1")
    0.09   0.0405  0.0405  0.0405     # 3  (1-3")
    0.015  0.00675 0.00675 0.00675    # 4  (3-6")   — LARGE CWD 0.015/yr (Brown et al. 1998)
    0.015  0.00675 0.00675 0.00675    # 5  (6-12")
    0.015  0.00675 0.00675 0.00675    # 6  (12-20")
    0.015  0.00675 0.00675 0.00675    # 7  (20-35")
    0.015  0.00675 0.00675 0.00675    # 8  (35-50")
    0.015  0.00675 0.00675 0.00675    # 9  (>50")
    0.5    0.225   0.225   0.225      # 10 litter (cr/fmvinit.f:99)
    0.002  0.0009  0.0009  0.0009     # 11 duff (cr/fmvinit.f:100)
]
# Oregon Coast R6-Oregon decay rates (oc/fmcba.f:610-671, Kim Mellen-McLean CWD model, forest 711 = R6 Oregon
# KODFOR≥600 branch). Base DKRT by size×decay-class (small woody 1-3" = 0.076/0.081/0.090/0.113; large 4-9" =
# 0.019/0.025/0.033/0.058; litter 0.5; duff 0.002), then ×DKRADJ(TEMP,MOIST,K) habitat adjustment. Applied at
# the FFE stand's default ecoclass CWC221 = a MESIC Douglas-fir site ⇒ TEMP=2/MOIST=2 ⇒ DKRADJ(2,2,K) = (1,
# 1.7, 1) for K=1(size≤3)/2(size4-5)/3(size≥6) — so the 3-12" classes decay ×1.7. (Exact per-stand ITYPE→
# CAHMC/CAWMD habitat resolution is a refinement; the mesic default fits the CWC221 reference stand.) The SN
# default `_FM_DKR` (woody 0.07-0.11) decayed OC's LARGE down-wood ~4-6× too fast ⇒ under-loaded 2003 surface
# fuel ⇒ the fire under-carried. Same class as the CR/NC DKR fixes above.
const _FM_DKR_OC = Float32[
    0.076  0.081  0.090  0.113     # 1  (<0.25")   K=1 ×1.0
    0.076  0.081  0.090  0.113     # 2  (0.25-1")  K=1 ×1.0
    0.076  0.081  0.090  0.113     # 3  (1-3")     K=1 ×1.0
    0.0323 0.0425 0.0561 0.0986    # 4  (3-6")     K=2 ×1.7
    0.0323 0.0425 0.0561 0.0986    # 5  (6-12")    K=2 ×1.7
    0.019  0.025  0.033  0.058     # 6  (12-20")   K=3 ×1.0
    0.019  0.025  0.033  0.058     # 7  (20-35")   K=3 ×1.0
    0.019  0.025  0.033  0.058     # 8  (35-50")   K=3 ×1.0
    0.019  0.025  0.033  0.058     # 9  (>50")     K=3 ×1.0
    0.5    0.5    0.5    0.5        # 10 litter (not habitat-adjusted; fmcba.f DO I=1,9)
    0.002  0.002  0.002  0.002     # 11 duff
]

# Variant-default DKR: SN uses `_FM_DKR`; LS/NE/CS/CR use their own faithful tables (…/fmvinit.f).
_fm_dkr_default(::AbstractVariant) = _FM_DKR
_fm_dkr_default(::OregonCoast) = _FM_DKR_OC      # oc/fmcba.f R6-Oregon DKRT × DKRADJ(mesic CWC221)
_fm_dkr_default(::LakeStates) = _FM_DKR_LS
_fm_dkr_default(::Northeast) = _FM_DKR_NE
_fm_dkr_default(::CentralStates) = _FM_DKR_CS
_fm_dkr_default(::CentralRockies) = _FM_DKR_CR
# North-Rockies decay table (ie/kt/em/ci fmvinit.f): the SAME base woody/litter/duff rates as CR
# (0.12/0.12/0.09 fine, 0.015 coarse, 0.5 litter, 0.002 duff), but DECAY-CLASS-INDEPENDENT — these four
# variants set `DKR(I,J)=DKR(I,1)` for all decay classes 2-4 and DO NOT apply the CR/UT/TT `×0.45`
# decay-class modifier (verified: ie/fmvinit.f:77-88, kt/:73-83, em/:87-96, ci/:118-128 all lack the
# `MODIFER FOR CR/UT` block that CR/UT/TT carry, cr/fmvinit.f:97-113). The earlier assignment to
# `_FM_DKR_CR` was a MISVERIFICATION (it checked fmcwd.f's decay ALGORITHM, identical across the family,
# not fmvinit.f's DKR DATA + the ×0.45): it over-retained the fine dead fuel ~2.2× in the dominant decay
# class 4 (0.054 vs 0.12 woody, 0.225 vs 0.5 litter) ⇒ the IE surface fire over-carried on stale fine
# fuel and needed the fuel loop disabled to compensate. This table restores the faithful IE/KT/EM/CI rates.
const _FM_DKR_NR = Float32[
    0.12  0.12  0.12  0.12     # 1  (<0.25")
    0.12  0.12  0.12  0.12     # 2  (0.25-1")
    0.09  0.09  0.09  0.09     # 3  (1-3")
    0.015 0.015 0.015 0.015    # 4  (3-6")
    0.015 0.015 0.015 0.015    # 5  (6-12")
    0.015 0.015 0.015 0.015    # 6  (12-20")
    0.015 0.015 0.015 0.015    # 7  (20-35")
    0.015 0.015 0.015 0.015    # 8  (35-50")
    0.015 0.015 0.015 0.015    # 9  (>50")
    0.5   0.5   0.5   0.5       # 10 litter (ie/fmvinit.f:86)
    0.002 0.002 0.002 0.002    # 11 duff   (ie/fmvinit.f:87)
]
_fm_dkr_default(::InlandEmpire) = _FM_DKR_NR     # ie/fmvinit.f:77-88 — DKR(I,J)=DKR(I,1), NO ×0.45
_fm_dkr_default(::Kootenai) = _FM_DKR_NR         # kt/fmvinit.f:73-83 — NO ×0.45
_fm_dkr_default(::EasternMontana) = _FM_DKR_NR   # em/fmvinit.f:87-96 — NO ×0.45
_fm_dkr_default(::CentralIdaho) = _FM_DKR_NR     # ci/fmvinit.f:118-128 — NO ×0.45
# canada/fire/bc/fmvinit.f:66-93: DKR(1:3,1) 0.12/0.12/0.09, the 3"+ woody classes 0.03 (Klenner 2021, vs the N-Rockies
# 0.015), DKR(I,J)=DKR(I,1); litter DKR(10,J)=0.355 and duff DKR(11,J)=0.0033 (Kurz et al. 2009) — no ×0.45.
const _FM_DKR_BC = let m = copy(_FM_DKR_NR)
    m[4:9, :] .= 0.03f0; m[10, :] .= 0.355f0; m[11, :] .= 0.0033f0
    m
end
_fm_dkr_default(::BritishColumbia) = _FM_DKR_BC
# TT (tt/fmvinit.f:72-92): DKR(I,J) = DKR(I,1)·0.45 for the woody classes I = 1..9 only, THEN litter/duff DKR(10,J) = 0.5,
# DKR(11,J) = 0.002 for every decay class — unlike CR/UT, whose ×0.45 loop runs to MXFLCL after the litter/duff set
# (cr|ut/fmvinit.f), giving 0.225/0.0009. jl gave TT the CR table ⇒ litter decayed at 0.225 not 0.5 (MEASURED FVStt_g16
# 2780339010690 FVS_Fuels litter 2010: live 1.363 t/ac, jl 2.832).
const _FM_DKR_TT = let m = copy(_FM_DKR_CR)
    m[10, :] .= 0.5f0; m[11, :] .= 0.002f0
    m
end
_fm_dkr_default(::Teton) = _FM_DKR_TT
_fm_dkr_default(::Utah) = _FM_DKR_CR             # ut/fmvinit.f — DKR ×0.45 (== CR)
_fm_dkr_default(::SoutheastAlaska) = _FM_DKR_AK   # ak/fmvinit.f DKR(11,4) (data/southeastalaska/fire/ffe_fuel.jl)
# NC (Klamath) base decay table (nc/fmvinit.f:70-92) — decay-class-INDEPENDENT and MUCH slower than the SN
# default (woody 0.0125-0.025 vs SN 0.07-0.11); subsequently ×DCYMLT (nc/fmcba.f:405, Dunning-code/site index).
# Without this NC fell through to the SN `_FM_DKR` ⇒ LARGE down-wood decayed ~3.7× too fast ⇒ low fuel-model
# weight on the hot model ⇒ weak surface fire. Used as the base that the first-FFE-year DCYMLT then scales.
const _FM_DKR_NC = Float32[
    0.025  0.025  0.025  0.025     # 1  (<0.25")
    0.025  0.025  0.025  0.025     # 2  (0.25-1")
    0.025  0.025  0.025  0.025     # 3  (1-3")
    0.0125 0.0125 0.0125 0.0125    # 4  (3-6")
    0.0125 0.0125 0.0125 0.0125    # 5  (6-12")
    0.0125 0.0125 0.0125 0.0125    # 6  (12-20")
    0.0125 0.0125 0.0125 0.0125    # 7  (20-35")
    0.0125 0.0125 0.0125 0.0125    # 8  (35-50")
    0.0125 0.0125 0.0125 0.0125    # 9  (>50")
    0.5    0.5    0.5    0.5        # 10 litter
    0.002  0.002  0.002  0.002     # 11 duff
]
_fm_dkr_default(::Klamath) = _FM_DKR_NC          # fallback; the first FFE year replaces it with the ×DCYMLT table

# DCYMLT (nc/fmcba.f:395-405 + dunn.f GETDUNN = DUNN50): the site index of the site species → a Dunning
# code (0-7, site↑ ⇒ code↓) → the decay-rate multiplier (Dunning 0-1 → 1.5, 2-4 → 1.0, 5-7 → 0.5).
const _NC_DUNN_XSR = (23f0, 31f0, 39f0, 49f0, 56f0, 75f0, 90f0, 106f0)   # DUNN50, ascending
const _NC_DUNN_YSR = Float32[7, 6, 5, 4, 3, 2, 1, 0]                     # Dunning codes
const _NC_DCY_XD   = (0f0, 1f0, 2f0, 3f0, 4f0, 5f0, 6f0, 7f0)
const _NC_DCY_YD   = Float32[1.5, 1.5, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5]
@inline function nc_dcymlt(si::Float32)::Float32
    dun = _ffe_algslp(si, _NC_DUNN_XSR, _NC_DUNN_YSR)
    return _ffe_algslp(dun, _NC_DCY_XD, _NC_DCY_YD)
end
# NC decay table scaled by DCYMLT (all classes; nct01 has no FUELDCAY ⇒ every SETDECAY<0 ⇒ all scaled).
@inline nc_adjusted_dkr(si::Float32)::Matrix{Float32} = _FM_DKR_NC .* nc_dcymlt(si)

# WS (ws/fmvinit.f:120-150): the NC table (0.025 / 0.0125 woody, litter 0.5, duff 0.002, PRDUFF 0.02), scaled at the
# first FFE year by the same Dunning DCYMLT (ws/fmcba.f:541-581; ws/dunn.f GETDUNN = DUNN50) on SITEAR(ISISP) — RF's
# SITEAR(7) when the site species is GB(21)/MC(41), whose site indices are 100-year Dunning values.
_fm_dkr_default(::WestSierra) = _FM_DKR_NC
function ws_adjusted_dkr(sitear::AbstractVector{Float32}, isisp::Int)::Matrix{Float32}
    si = (isisp == 21 || isisp == 41) ? sitear[7] : ((1 <= isisp <= length(sitear)) ? sitear[isisp] : 0f0)
    return _FM_DKR_NC .* nc_dcymlt(si)
end

# CA (ca/fmcba.f:579-745, the decay rates CA sets at the first FFE year instead of fmvinit.f): KODFOR 500-599 (R5) the
# California table (= the NC values), else the Oregon base table (= the EC/SO base) × DKRADJ(CAHMC(ITYPE),CAWMD(ITYPE),K)
# capped at 1 with the small-wood bump (r6_adjusted_dkr); litter 0.5, duff 0.002; then, for KODFOR < 600, ×DCYMLT on
# SITEAR(ISISP) (GETDUNN = DUNN50, as NC/WS). ITYPE is HABTYP's PCOML index (46 = CWC221 when undecoded, ca/sitset.f).
const _FM_CAHMC = Int8[2,2,2,2,2,1,2,1,2,2, 2,2,2,1,2,1,1,2,2,2, 2,2,2,2,2,2,2,3,2,2, 2,2,2,2,3,3,2,2,2,3,
                       2,2,2,2,2,2,2,1,1,1, 3,2,2,2,2,3,3,3,2,3, 2,2,3,3,2,2,2,2,2,2, 2,2,2,2,1,1,2,2,2,2,
                       2,2,2,2,1,2,2,2,1,2]
const _FM_CAWMD = Int8[2,3,3,2,2,3,3,3,2,2, 2,2,2,3,2,3,3,2,3,1, 1,2,2,1,2,2,2,2,3,2, 3,3,3,3,2,1,2,2,2,2,
                       1,2,2,2,2,2,3,3,3,3, 1,1,2,2,2,2,2,2,2,1, 2,2,2,1,2,2,3,2,2,2, 2,2,2,2,3,3,1,2,2,2,
                       2,2,2,2,3,2,2,2,3,2]
_fm_dkr_default(::CentralCalifornia) = _FM_DKR_NC
function ca_adjusted_dkr(kodfor::Int, itype::Int, sitear::AbstractVector{Float32}, isisp::Int)::Matrix{Float32}
    if 500 <= kodfor < 600
        dkr = copy(_FM_DKR_NC)
    else
        it = (1 <= itype <= 90) ? itype : 46
        dkr = r6_adjusted_dkr(_FM_DKR_EC, Int(_FM_CAHMC[it]), Int(_FM_CAWMD[it]))
        dkr[10, :] .= 0.5f0; dkr[11, :] .= 0.002f0
    end
    if kodfor < 600
        si = (1 <= isisp <= length(sitear)) ? sitear[isisp] : 0f0
        dkr .*= nc_dcymlt(si)
    end
    return dkr
end

# BM has its OWN base decay table (bm/fmvinit.f:68-113) — NOT the CR table. Faster litter (0.65 vs CR 0.5)
# and different woody rates; used as the base the habitat DKRADJ then scales (see bm_adjusted_dkr).
const _FM_DKR_BM = Float32[
    0.076  0.081  0.090  0.113      # 1  (<0.25")
    0.076  0.081  0.090  0.113      # 2  (0.25-1")
    0.076  0.081  0.090  0.113      # 3  (1-3")
    0.019  0.025  0.033  0.058      # 4  (3-6")
    0.019  0.025  0.033  0.058      # 5  (6-12")
    0.019  0.025  0.033  0.058      # 6  (12-20")
    0.019  0.025  0.033  0.058      # 7  (20-35")
    0.019  0.025  0.033  0.058      # 8  (35-50")
    0.019  0.025  0.033  0.058      # 9  (>50")
    0.65   0.65   0.65   0.65       # 10 litter (bm/fmvinit.f:111)
    0.002  0.002  0.002  0.002      # 11 duff  (bm/fmvinit.f:112)
]
_fm_dkr_default(::BlueMountains) = _FM_DKR_BM     # bm/fmvinit.f — distinct from CR; DKRADJ-scaled at 1st yr

# EC has its OWN base decay table (ec/fmvinit.f:79-124) — the SAME woody rates as BM (0.076-0.113 fine,
# 0.019-0.058 coarse) but litter = 0.50/yr (ec/fmvinit.f:123, NOT BM 0.65). Used as the base the habitat
# DKRADJ then scales (ec_adjusted_dkr). Without it EC fell through to the SN _FM_DKR ⇒ coarse down-wood
# decayed ~3× too fast ⇒ LARGE down-wood ~1.8× low ⇒ FMDYN dropped the natural-fuel model 10 ⇒ under-fire.
const _FM_DKR_EC = Float32[
    0.076  0.081  0.090  0.113      # 1  (<0.25")
    0.076  0.081  0.090  0.113      # 2  (0.25-1")
    0.076  0.081  0.090  0.113      # 3  (1-3")
    0.019  0.025  0.033  0.058      # 4  (3-6")
    0.019  0.025  0.033  0.058      # 5  (6-12")
    0.019  0.025  0.033  0.058      # 6  (12-20")
    0.019  0.025  0.033  0.058      # 7  (20-35")
    0.019  0.025  0.033  0.058      # 8  (35-50")
    0.019  0.025  0.033  0.058      # 9  (>50")
    0.50   0.50   0.50   0.50       # 10 litter (ec/fmvinit.f:123)
    0.002  0.002  0.002  0.002      # 11 duff  (ec/fmvinit.f:124)
]
_fm_dkr_default(::EastCascades) = _FM_DKR_EC      # ec/fmvinit.f — DKRADJ-scaled at 1st yr

# WC/PN/OP base decay table (wc/pn/op fmvinit.f:68-113, identical): Mellen-McLean CWD workshop rates by decay class,
# litter 0.5, duff 0.002; scaled by DKRADJ(TEMP,MOIST,K) from WCWMC/WCWMD (wc/fmcba.f:490) or PNWMC/PNWMD
# (pn,op/fmcba.f:464) at the first FFE year. jl ran WC/PN on the SN table (coarse wood 0.07/yr vs 0.012-0.077,
# litter 0.65) — pnt01/wct01 stand 4 LARGE down wood at the 2003 fire ~35% low.
const _FM_DKR_WC = Float32[
    0.069  0.081  0.097  0.131      # 1  (<0.25")
    0.069  0.081  0.097  0.131      # 2  (0.25-1")
    0.069  0.081  0.097  0.131      # 3  (1-3")
    0.012  0.025  0.041  0.077      # 4  (3-6")
    0.012  0.025  0.041  0.077      # 5  (6-12")
    0.012  0.025  0.041  0.077      # 6  (12-20")
    0.012  0.025  0.041  0.077      # 7  (20-35")
    0.012  0.025  0.041  0.077      # 8  (35-50")
    0.012  0.025  0.041  0.077      # 9  (>50")
    0.50   0.50   0.50   0.50       # 10 litter
    0.002  0.002  0.002  0.002      # 11 duff
]
_fm_dkr_default(::Union{WestCascades,PacificNorthwest,Olympic}) = _FM_DKR_WC

# wc/fmcba.f and pn,op/fmcba.f carry their OWN DKRADJ(TEMP,MOIST,K) tables (DATA (((DKRADJ(I,J,K),K=1,3),J=1,3),I=1,3)),
# not the BM/EC/SO one (_FM_DKRADJ). WC stand 4 (ITYPE 52, TEMP 3 / MOIST 1): BM/EC table ×1.21/1.79/1.21 vs WC's own
# ×0.75/1.0/0.75 — the down wood decayed ~2× too fast (no-fire 2033 Surface_ge3 6.8 vs live 11.0).
function _dkradj_from(vals::Vector{Float32})
    a = Array{Float32}(undef, 3, 3, 3); n = 0
    for i in 1:3, j in 1:3, k in 1:3
        n += 1; a[i, j, k] = vals[n]
    end
    a
end
const _FM_DKRADJ_WC = _dkradj_from(Float32[
    1.35, 2, 1.35,  1.49, 2, 1.49,  1.7, 2, 1.7,
    0.875, 1.5, 0.875,  1, 2, 1,  1.21, 2, 1.21,
    0.75, 1, 0.75,  0.825, 1.3, 0.825,  0.875, 1.5, 0.875])
const _FM_DKRADJ_PN = _dkradj_from(Float32[      # pn/fmcba.f == op/fmcba.f
    1.21, 2, 1.21,  1.35, 2, 1.35,  1.7, 2, 1.7,
    0.825, 1.3, 0.825,  1, 2, 1,  1.35, 2, 1.35,
    0.75, 1, 0.75,  0.75, 1, 0.75,  0.925, 1.7, 0.925])

"""
    r6_adjusted_dkr(base, temp, moist) -> Matrix{Float32}

The R6 first-year habitat decay adjustment shared by bm/ec/wc/pn/op fmcba.f: woody classes 1-9 scaled by
DKRADJ(TEMP,MOIST,K) (K = 1 for sizes 1-3, 2 for 4-5, 3 for 6-9) capped at 1, then (sizes 9→2) a smaller class
decaying slower than the next larger one is bumped up to it. Litter/duff keep the base.
"""
function r6_adjusted_dkr(base::Matrix{Float32}, temp::Integer, moist::Integer;
                         adj::Array{Float32,3} = _FM_DKRADJ)::Matrix{Float32}
    dkrtab = adj
    dkr = copy(base)
    @inbounds for i in 1:9
        k = i <= 3 ? 1 : (i <= 5 ? 2 : 3)
        adj = dkrtab[temp, moist, k]
        for j in 1:4
            v = dkr[i, j] * adj
            dkr[i, j] = v > 1f0 ? 1f0 : v
        end
    end
    @inbounds for i in 9:-1:2, j in 1:4
        (dkr[i, j] - dkr[i-1, j]) > 0f0 && (dkr[i-1, j] = dkr[i, j])
    end
    return dkr
end

# SO (SouthCentralOregon) Oregon base decay table (so/fmcba.f:770-808) — BYTE-IDENTICAL to the EC/BM woody
# rates (0.076-0.113 fine, 0.019-0.058 coarse) with litter 0.50/yr (so/fmcba.f:844) and duff 0.002 — i.e. the
# same matrix as _FM_DKR_EC. Reused via the alias. (SO's California-forest branch uses a flat 0.025/0.0125
# table — a KODFOR 500-599/701 path not exercised by the R6 Deschutes reference stand.) Scaled by the habitat
# DKRADJ(TEMP,MOIST,K) at the first FFE year (so_adjusted_dkr).
const _FM_DKR_SO = _FM_DKR_EC
_fm_dkr_default(::SouthCentralOregon) = _FM_DKR_SO

# SO habitat → temperature (SOHMC) / moisture (SOWMD) class (so/fmcba.f:93-119, from FMR6SDCY). Same DKRADJ
# table as BM/EC (_FM_DKRADJ). 92 plant-association codes; SOHMC 1=hot/2=mod/3=cold, SOWMD 1=wet/2=mesic/3=dry.
const _FM_SOHMC = Int8[
    2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
    3, 1, 2, 2, 2, 2, 2, 2, 2, 3,
    3, 3, 1, 3, 1, 2, 2, 2, 1, 1,
    1, 1, 1, 1, 2, 2, 1, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 1, 3, 2, 2, 2, 1,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2]
const _FM_SOWMD = Int8[
    2, 2, 2, 1, 1, 2, 1, 1, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 2,
    1, 2, 1, 1, 2, 1, 1, 2, 1, 3,
    1, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 2, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 2, 3, 3, 3, 2, 3, 2, 2,
    2, 2, 3, 2, 3, 1, 2, 2, 2, 3,
    1, 3, 3, 3, 3, 2, 3, 2, 2, 2,
    1, 2]

"""
    so_adjusted_dkr(itype) -> Matrix{Float32}

SO habitat-conditioned decay rates (so/fmcba.f:764-836): the SO base DKR scaled by `DKRADJ(TEMP,MOIST,K)`
(the SAME table BM/EC use) for the stand's habitat `itype` (SOHMC/SOWMD), capped at 1.0, then a second pass
(size 9→2) bumps any size class decaying slower than the next-larger class up to the larger's rate. Only
woody classes 1-9; litter (10)/duff (11) keep the base. Applied once at the first FFE year (no FuelDcay).
"""
# so/fmcba.f:746-761,833-837: on the California forests (KODFOR 500-599, 701) the SO decay rates are the flat CA set —
# 0.025/yr for the <3" classes, 0.0125 for 3"+, every decay class, litter 0.5, duff 0.002 — with NO habitat DKRADJ
# (that adjustment is the Oregon branch only). MEASURED FVSso_g16 15364795010497 (forest 505) FVS_Fuels 2020: surface
# <3" 7.97 live vs 5.83 jl on the Oregon matrix.
function so_california_dkr()::Matrix{Float32}
    dkr = copy(_FM_DKR_SO)
    @inbounds for i in 1:9, j in 1:4
        dkr[i, j] = i <= 3 ? 0.025f0 : 0.0125f0
    end
    @inbounds for j in 1:4
        dkr[10, j] = 0.5f0; dkr[11, j] = 0.002f0
    end
    return dkr
end

function so_adjusted_dkr(itype::Integer)::Matrix{Float32}
    dkr = copy(_FM_DKR_SO)
    (itype < 1 || itype > length(_FM_SOHMC)) && return dkr
    temp = Int(_FM_SOHMC[itype]); moist = Int(_FM_SOWMD[itype])
    @inbounds for i in 1:9
        k = i <= 3 ? 1 : (i <= 5 ? 2 : 3)
        adj = _FM_DKRADJ[temp, moist, k]
        for j in 1:4
            v = dkr[i, j] * adj
            dkr[i, j] = v > 1f0 ? 1f0 : v
        end
    end
    @inbounds for i in 9:-1:2, j in 1:4
        (dkr[i, j] - dkr[i-1, j]) > 0f0 && (dkr[i-1, j] = dkr[i, j])
    end
    return dkr
end

# EC habitat → temperature (ECHMC) / moisture (ECWMD) class (ec/fmcba.f:82-124, from FMR6SDCY). Same DKRADJ
# table as BM (_FM_DKRADJ). 155 habitat codes; 1=hot/2=mod/3=cold (ECHMC), 1=wet/2=mesic/3=dry (ECWMD).
const _FM_ECHMC = Int8[
    3,3,2,2,2,2,2,2,2,2, 2,2,2,2,1,2,2,2,2,2, 2,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,2,2,2,2, 2,2,2,2,2,2,3,3,3,3,
    3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3, 3,3,2,2,2,3,3,2,2,3,
    2,2,2,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,3,3,
    3,3,3,3,3,3,3,3,3,1, 1,1,2,1,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2]
const _FM_ECWMD = Int8[
    3,3,2,2,2,2,1,2,2,2, 3,3,3,3,3,2,3,3,3,3, 2,2,3,3,2,2,3,2,2,2,
    2,3,2,2,2,3,3,3,3,3, 2,2,2,2,2,2,3,2,2,2,
    1,2,2,2,3,3,1,2,2,3, 1,2,2,2,2,3,3,2,2,3, 3,3,2,2,2,2,2,1,1,1,
    2,2,2,2,2,2,2,2,3,2, 2,2,2,2,2,2,1,2,3,2,
    3,2,1,1,2,2,2,1,2,3, 3,3,3,3,2,2,2,2,2,2, 2,3,3,3,3,3,2,2,2,2,
    2,2,2,2,2,2,2,3,2,3, 2,3,2,3,3,2,2,2,2,2,
    2,3,2,3,2]

"""
    ec_adjusted_dkr(itype) -> Matrix{Float32}

EC habitat-conditioned decay rates (ec/fmcba.f:457-491): the EC base DKR scaled by `DKRADJ(TEMP,MOIST,K)`
(the SAME table BM uses) for the stand's habitat `itype` (ECHMC/ECWMD), capped at 1.0, then a second pass
(size 9→2) bumps any size class decaying slower than the next-larger class up to the larger's rate. Only
woody classes 1-9; litter (10)/duff (11) keep the base. Applied once at the first FFE year (no FuelDcay).
"""
function ec_adjusted_dkr(itype::Integer)::Matrix{Float32}
    dkr = copy(_FM_DKR_EC)
    (itype < 1 || itype > length(_FM_ECHMC)) && return dkr
    temp = Int(_FM_ECHMC[itype]); moist = Int(_FM_ECWMD[itype])
    @inbounds for i in 1:9
        k = i <= 3 ? 1 : (i <= 5 ? 2 : 3)
        adj = _FM_DKRADJ[temp, moist, k]
        for j in 1:4
            v = dkr[i, j] * adj
            dkr[i, j] = v > 1f0 ? 1f0 : v
        end
    end
    @inbounds for i in 9:-1:2, j in 1:4
        (dkr[i, j] - dkr[i-1, j]) > 0f0 && (dkr[i-1, j] = dkr[i, j])
    end
    return dkr
end
const _FM_PRDUFF = 0.02f0   # proportion of decayed woody material that becomes duff (fmvinit.f:112)

# ── BM decay-rate habitat adjustment (bm/fmcba.f:67-113, 333-368) ──────────────────────────────────
# BM is the lone western variant that DOES NOT use the CR decay table verbatim: at the first FFE year it
# multiplies the base DKR by a habitat-conditioned factor DKRADJ(TEMP,MOIST,K) (from FMR6SDCY). Each BM
# habitat type (ITYPE, 1-92 = the decoded KODTYP / PCOML index) maps to a temperature class BMHMC (1=hot,
# 2=moderate, 3=cold) and a moisture class BMWMD (1=wet, 2=mesic, 3=dry); K is a size-class group
# (1: size 1-3 / <3", 2: size 4-5 / 3-12", 3: size 6-9 / >12"). Without this, BM down-wood decays at the
# generic CR rate and OVER-ACCUMULATES over the projection (measured: bmt01 SMALL/LARGE ~2.4× live), which
# pushes FMDYN onto the hot fuel model → spurious crown fire → mortality over-kill.
const _FM_BMHMC = Int8[   # temperature class by habitat type (bm/fmcba.f:73-84)
    3,3,2,2,2,2,2,2,2,2, 2,2,2,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3,
    2,3,3,2,3,2,2,2,3,3, 3,2,2,2,2,2,3,3,1,1,
    1,1,2,2,2,1,2,2,1,2, 1,1,1,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,2,2,
    2,2,2,2,2,2,3,2,2,2, 2,1]
const _FM_BMWMD = Int8[   # moisture class by habitat type (bm/fmcba.f:88-99)
    3,3,3,3,3,2,3,3,3,3, 3,2,3,2,1,2,3,1,1,2, 1,1,2,2,2,2,2,3,2,3,
    2,3,3,1,1,1,2,1,2,3, 3,3,2,2,2,2,3,3,3,3,
    3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,1,1,2, 2,2,2,1,1,1,3,3,3,2,
    2,2,3,3,2,1,3,2,1,2, 2,2]
# DKRADJ[TEMP, MOIST, K] (bm/fmcba.f:101 — DATA order fills K innermost, then MOIST(J), then TEMP(I)).
const _FM_DKRADJ = let a = Array{Float32}(undef, 3, 3, 3)
    vals = Float32[
        1.7,  2.0,  1.7,   1.49, 1.91, 1.49,  0.75, 0.85, 0.75,    # TEMP=1 (hot):  MOIST 1,2,3
        1.35, 1.85, 1.35,  1.0,  1.7,  1.0,   0.875, 1.2, 0.875,   # TEMP=2 (mod):  MOIST 1,2,3
        1.21, 1.79, 1.21,  1.14, 1.76, 1.14,  0.75, 0.85, 0.75]    # TEMP=3 (cold): MOIST 1,2,3
    n = 0
    for i in 1:3, j in 1:3, k in 1:3
        n += 1; a[i, j, k] = vals[n]
    end
    a
end

"""
    bm_adjusted_dkr(itype) -> Matrix{Float32}

BM habitat-conditioned decay rates (bm/fmcba.f:333-368): the BM base DKR scaled by `DKRADJ(TEMP,MOIST,K)`
for the stand's habitat type `itype`, capped at 1.0, then a second pass (size 9→2) bumps any size class that
would decay slower than the next-larger class up to the larger class's rate. Only the woody classes 1-9 are
adjusted; litter (10) and duff (11) keep the CR default. Applied once, at the first FFE year, when the user
has not overridden decay with FuelDcay.
"""
function bm_adjusted_dkr(itype::Integer)::Matrix{Float32}
    dkr = copy(_FM_DKR_BM)
    (itype < 1 || itype > length(_FM_BMHMC)) && return dkr
    temp = Int(_FM_BMHMC[itype]); moist = Int(_FM_BMWMD[itype])
    @inbounds for i in 1:9
        k = i <= 3 ? 1 : (i <= 5 ? 2 : 3)
        adj = _FM_DKRADJ[temp, moist, k]
        for j in 1:4
            v = dkr[i, j] * adj
            dkr[i, j] = v > 1f0 ? 1f0 : v
        end
    end
    # bump smaller wood up to larger wood's rate where it would otherwise decay more slowly (fmcba.f:359-368)
    @inbounds for i in 9:-1:2, j in 1:4
        (dkr[i, j] - dkr[i-1, j]) > 0f0 && (dkr[i-1, j] = dkr[i, j])
    end
    return dkr
end

"""
    apply_fuelmove!(s) -> Bool

FUELMOVE (act 2530, fmtret.f:203-368): transfer surface fuel between size categories at a scheduled cycle.
Each due activity moves XGET = max(amount, proportion·source, source−leave, target−current) tons/ac (capped
at the available source) from size class FROM to TO; size class 0 is the import (FROM=0) / export (TO=0)
sink. The per-size-class totals are then written back by scaling each class's cwd sub-pools (soft/hard ×
decay) by new/old, or dumping into the hard/fast pool if the class was empty. No-op without a due FUELMOVE.
"""
function apply_fuelmove!(s::StandState)::Bool
    fs = s.fire
    (fs === nothing || !fs.active || isempty(s.control.schedule)) && return false
    yr = Int(current_cycle_year(s)); fvscyc = Int(s.control.cycle) + 1
    cwd = fs.cwd
    # FORG/FSRC/FTRG by size class 0:11 (0 = outside sink); +1 offset so idx 1 = class 0, idx j+1 = class j.
    forg = zeros(Float32, 12)
    @inbounds for j in 1:11; forg[j+1] = sum(@view cwd[j, :, :]); end
    fsrc = copy(forg); ftrg = zeros(Float32, 12); altered = false
    for a in s.control.schedule
        a.icflag == Int32(2530) || continue
        (Int(a.year) == yr || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc)) || continue
        ifrm = Int(round(a.params[1])); ito = Int(round(a.params[2]))
        x = a.params[3]; y = a.params[4]; z = a.params[5]; q = a.params[6]
        (0 <= ifrm <= 11 && 0 <= ito <= 11 && ifrm != ito && x >= 0f0 && 0f0 <= y <= 1f0 && z >= 0f0) || continue
        fi = ifrm + 1; ti = ito + 1
        xget = 0f0
        if ifrm > 0
            fsrc[fi] <= 0f0 && continue
            xget = q >= 0f0 ? max(x, y * fsrc[fi], fsrc[fi] - z, q - fsrc[ti]) : max(x, y * fsrc[fi], fsrc[fi] - z)
            xget > fsrc[fi] && (xget = fsrc[fi])
            fsrc[fi] -= xget
        else
            xget = q >= 0f0 ? max(x, q - fsrc[ti]) : x
        end
        ftrg[ti] += xget
        xget > 0f0 && (altered = true)
    end
    altered || return false
    @inbounds for j in 1:11
        ft = fsrc[j+1] + ftrg[j+1]
        abs(forg[j+1] - ft) >= 1f-6 || continue
        if forg[j+1] <= 1f-6
            cwd[j, 2, 3] = ft                              # empty class → all to hard/fast (CWD(1,J1,2,3))
        else
            sc = ft / forg[j+1]
            for k in 1:2, l in 1:4; cwd[j, k, l] *= sc; end
        end
    end
    return true
end

"""
    fmcwd!(s, nyrs) -> StandState

Apply `nyrs` years of FFE surface-fuel decay to `fire.cwd` (FMCWD, fmcwd.f:78-134). Duff is decayed
first (so woody-decay duff additions land after), then each woody/litter class 1-10: a PRDUFF fraction
of the decayed amount is moved to duff, the pool is reduced by `(1−DKR·{1.1 soft})^nyrs`, and (woody
classes < 10) a `nyrs·ln(1−DKR)/ln(0.64)` fraction of the hard pool transfers to soft. No-op unless
FFE is active.
"""
function fmcwd!(s::StandState, nyrs::Integer)
    fs = s.fire
    (fs === nothing || !fs.active) && return s
    cwd = fs.cwd; n = Float32(nyrs); ni = Int(nyrs)   # (1-DKR)**NYRS is REAL**INTEGER (__powisf2, fpowi)
    # FUELMULT/FUELDCAY override the DKR matrix; DUFFPROD overrides PRDUFF — both fall back to defaults
    # when unset (size 0×0).
    dkr = size(fs.params.dkr, 1) == 11 ? fs.params.dkr : _fm_dkr_default(s.variant)
    has_pd = size(fs.params.prduff, 1) == 11; pdm = fs.params.prduff
    @inbounds for L in 1:4
        # duff (size 11) first, so woody decay can add to it below
        cwd[11, 1, L] *= fpowi(1f0 - dkr[11, L] * 1.1f0, ni)
        cwd[11, 2, L] *= fpowi(1f0 - dkr[11, L], ni)
        cwd[11, 1, L] < 0f0 && (cwd[11, 1, L] = 0f0)
        cwd[11, 2, L] < 0f0 && (cwd[11, 2, L] = 0f0)
        for J in 1:10
            dk = dkr[J, L]
            pd = has_pd ? pdm[J, L] : _FM_PRDUFF       # DUFFPROD-overridable proportion-to-duff
            # amount decayed this cycle → a PRDUFF fraction becomes duff (added to the hard duff pool)
            amt = cwd[J, 1, L] - cwd[J, 1, L] * fpowi(1f0 - dk * 1.1f0, ni)
            amt < 1f-9 && (amt = 0f0); cwd[11, 2, L] += amt * pd
            amt = cwd[J, 2, L] - cwd[J, 2, L] * fpowi(1f0 - dk, ni)
            amt < 1f-9 && (amt = 0f0); cwd[11, 2, L] += amt * pd
            # decrease the pools
            cwd[J, 1, L] *= fpowi(1f0 - dk * 1.1f0, ni); cwd[J, 1, L] < 1f-9 && (cwd[J, 1, L] = 0f0)
            cwd[J, 2, L] *= fpowi(1f0 - dk, ni);        cwd[J, 2, L] < 1f-9 && (cwd[J, 2, L] = 0f0)
            # hard → soft transfer (woody classes only)
            if J < 10
                # fmcwd.f:117 TOSOFT = NYRS*(LOG(1-DKR(J,L)))/(LOG(0.64)) — REAL LOG = glibc logf (flog); Julia's own
                # Float32 log differs by an ULP on some DKR (MEASURED FVSbm_g16 22960873010497 salvage 2007: CWD(1,4,1,2)
                # 2.55130029 live vs 2.5513005 jl after the first FMCWD).
                tosoft = clamp(n * flog(1f0 - dk) / flog(0.64f0), 0f0, 1f0) * cwd[J, 2, L]
                cwd[J, 1, L] += tosoft; cwd[J, 2, L] -= tosoft
                cwd[J, 1, L] < 1f-9 && (cwd[J, 1, L] = 0f0)
                cwd[J, 2, L] < 1f-9 && (cwd[J, 2, L] = 0f0)
            end
        end
    end
    return s
end
