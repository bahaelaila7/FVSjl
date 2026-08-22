# =============================================================================
# IE (Inland Empire) FVS_TreeList crown width — national cwcalc.f, IEMAP dispatch.
#
# cwcalc.f is the shared NATIONAL crown-width routine: one flat SELECT CASE(CWEQN)
# keyed by a 5-char code (3-char FIA sp + 2-char model), and each variant supplies a
# per-species map (IEMAP for IE) that picks the code.  IE is Region 1/4 so KODFOR<601
# ⇒ the forest bias-factor section is skipped and BF=1.0 (vie/cwcalc.f:477 GO TO 10).
#
# The equation FORMS (_em_r1 Crookston-R1 log / _em_powf power / _em_r6m2 Crookston-R6
# model-2 / _em_bech1 Bechtold-2004 model-1 / _em_bech2 model-2) are the shared national
# forms already ported + validated bit-exact for EM (easternmontana/crown.jl, 353/353 vs
# FVSem_g16); IE reuses them verbatim.  Only IE's map + the 8 codes EM does not carry
# (11903 01703 26303 24203 10602 23104 32102 12205) are added here.
#
# _IE_CWMAP is IEMAP (vie/cwcalc.f:184) verbatim, indexed by the jl IE species index —
# which equals the vie ISPC index (proven by the bit-exact per-species DGF port).  Two
# FIA-122 species (PP idx10 '12203' log, OS idx23 '12205' R6M2) need the INDEX, not the
# FIA code, to disambiguate — hence the index-keyed map.
#
# REPORTING-ONLY (FVS_TreeList CrWidth); INERT for growth/mortality/CCF/volume.
# Validated: FVSie_clean ietl FVS_TreeList — the 7 emitted species (all Crookston-R1
# log-form: WP WL DF GF WH? RC? — present set WP/WL/DF/GF/ES/AF/LP) bit-exact.  The
# other 15 codes are transcription-faithful national cwcalc.f (un-exercised on this
# stand; the shared forms are the EM-validated ones).
# =============================================================================
const _IE_CWMAP = ("11903","07303","20203","01703","26303","24203","10803","09303",
                   "01903","12203","26405","10105","11301","07204","10602","06602",
                   "23104","74605","74902","32102","37506","74902","12205")

function ie_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                   el::Float32, hi::Float32)::Float32
    (1 <= sp <= 23) || return 0.5f0
    barea <= 1f0 && (barea = 1f0)          # cwcalc.f: IF(BAREA.LE.1.) BAREA=1.
    cl  = cr * h * 0.01f0                   # CL = CR*H*0.01 (CR is crown-ratio percent)
    ba1 = barea + 1f0
    eqn = _IE_CWMAP[sp]
    cw = if eqn == "11903";     _em_r1(1.0405f0,1.2799f0,0.11941f0,0.42745f0,0f0,-0.07182f0, 1f0,35f0, d,h,cl,barea)
    elseif eqn == "07303";      _em_r1(1.02478f0,0.99889f0,0.19422f0,0.59423f0,-0.09078f0,-0.02341f0, 1f0,40f0, d,h,cl,barea)
    elseif eqn == "20203";      _em_r1(1.01685f0,1.48372f0,0.27378f0,0.49646f0,-0.18669f0,-0.01509f0, 1f0,80f0, d,h,cl,barea)
    elseif eqn == "01703";      _em_r1(1.0303f0,1.14079f0,0.20904f0,0.38787f0,0f0,0f0, 1f0,40f0, d,h,cl,barea)
    elseif eqn == "26303";      _em_r1(1.02460f0,1.3522f0,0.24844f0,0.412117f0,-0.104357f0,0.03538f0, 0.1f0,54f0, d,h,cl,barea)
    elseif eqn == "24203";      _em_r1(1.03597f0,1.46111f0,0.26289f0,0.18779f0,0f0,0f0, 1f0,45f0, d,h,cl,barea)
    elseif eqn == "10803";      _em_r1(1.03992f0,1.58777f0,0.30812f0,0.64934f0,-0.38964f0,0f0, 0.7f0,40f0, d,h,cl,barea)
    elseif eqn == "09303";      _em_r1(1.02687f0,1.28027f0,0.2249f0,0.47075f0,-0.15911f0,0f0, 0.1f0,40f0, d,h,cl,barea)
    elseif eqn == "01903";      _em_r1(1.02886f0,1.01255f0,0.30374f0,0.37093f0,-0.13731f0,0f0, 0.1f0,30f0, d,h,cl,barea)
    elseif eqn == "12203";      _em_r1(1.02687f0,1.49085f0,0.1862f0,0.68272f0,-0.28242f0,0f0, 2f0,46f0, d,h,cl,barea)
    elseif eqn == "26405";      _em_r6m2(3.7854f0,0.54684f0,-0.12954f0,0.16151f0,0.03047f0,-0.00561f0, d,h,cl,ba1,el,10f0,79f0,45f0)
    elseif eqn == "10105";      _em_r6m2(2.2354f0,0.66680f0,-0.11658f0,0.16927f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,40f0)
    elseif eqn == "11301";      _em_bech1(4.0181f0,0.8528f0, d,25f0)
    elseif eqn == "07204";      _em_powf(2.2586f0,0.68532f0,33f0, d)
    elseif eqn == "10602";      _em_bech2(-5.4647f0,1.9660f0,-0.0395f0,0.0427f0,-0.0259f0, d,cr,hi,-40f0,11f0,25f0,true)
    elseif eqn == "06602";      _em_bech2(-4.1599f0,1.3528f0,-0.0233f0,0.0633f0,-0.0423f0, d,cr,hi,-37f0,19f0,29f0,true)
    elseif eqn == "23104";      _em_powf(6.1297f0,0.45424f0,30f0, d)
    elseif eqn == "74605";      _em_r6m2(4.7961f0,0.64167f0,-0.18695f0,0.18581f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,45f0)
    elseif eqn == "74902";      _em_bech2(4.1687f0,1.5355f0,0f0,0f0,0.1275f0, d,cr,hi,-26f0,-2f0,35f0,false)
    elseif eqn == "32102";      _em_bech2(5.9765f0,0.8648f0,0f0,0.0675f0,0f0, d,cr,hi,-9.9f9,9.9f9,39f0,false)
    elseif eqn == "37506";      _em_powf(5.8980f0,0.4841f0,25f0, d)
    elseif eqn == "12205";      _em_r6m2(4.7762f0,0.74126f0,-0.28734f0,0.17137f0,-0.00602f0,-0.00209f0, d,h,cl,ba1,el,13f0,75f0,50f0)
    else 0f0 end
    # cwcalc.f final CRWDTH clamp (after label 9000).
    cw < 0.5f0 && (cw = 0.5f0)
    cw > 99.9f0 && (cw = 99.9f0)
    return cw
end
