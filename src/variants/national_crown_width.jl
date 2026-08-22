# =============================================================================
# National FVS_TreeList crown-width dispatch — shared cwcalc.f SELECT CASE(CWEQN).
#
# cwcalc.f is ONE national routine: a flat SELECT CASE over a 5-char code (3-char FIA sp
# + 2-char model), and each variant supplies a per-species map ({V}MAP) that picks the
# code.  _cwcalc_national is that SELECT CASE, keyed by the code string; each variant's
# {v}_cwcalc = _cwcalc_national(_V_CWMAP[sp], …).  The equation FORMS are the shared
# national forms ported + validated bit-exact for EM (easternmontana/crown.jl, 353/353 vs
# FVSem_g16): _em_r1 (Crookston-R1 log), _em_powf (power a·D^b), _em_r6m2 (Crookston-R6
# model-2), _em_bech1 (Bechtold-2004 model-1), _em_bech2 (model-2).  Interior variants are
# Region 1/4 (KODFOR<601 ⇒ BF=1.0 skipped) so BF=1 throughout here.
#
# Codes present cover IE/KT/CI (and the shared subset of EM); TT/UT/BC add their own codes
# as they are validated.  REPORTING-ONLY (FVS_TreeList CrWidth); INERT for growth/vol/CCF.
# =============================================================================
# AK equation '08' (Region 10): general product a·D^b·H^ch·CL^ccl·(BAREA+1)^cba·EXP(EL)^cel.
# D<0.1 ⇒ flat CW=0.5 (NO small-tree ×(D/floor) scaling, unlike the R1/R6 forms).  Absent
# terms carry a 0 exponent ⇒ ^0 = ×1.0 (bit-exact no-op).  EL clamp [ello,elhi] (wide when
# the code has no EXP(EL) factor).  Transcendentals via glibc powf/expf (as the EM forms).
@inline function _cw08(a::Float32, b::Float32, ch::Float32, ccl::Float32, cba::Float32,
                       cel::Float32, ello::Float32, elhi::Float32, cap::Float32,
                       d::Float32, h::Float32, cl::Float32, ba1::Float32, el::Float32)::Float32
    d < 0.1f0 && return 0.5f0
    elc = el < ello ? ello : (el > elhi ? elhi : el)
    cw = a * _emcw_pow(d, b) * _emcw_pow(h, ch) * _emcw_pow(cl, ccl) *
         _emcw_pow(ba1, cba) * _emcw_pow(_emcw_exp(elc), cel)
    cw > cap && (cw = cap)
    return cw
end

function _cwcalc_national(eqn::AbstractString, d::Float32, h::Float32, cr::Float32,
                          barea::Float32, el::Float32, hi::Float32)::Float32
    barea <= 1f0 && (barea = 1f0)          # cwcalc.f: IF(BAREA.LE.1.) BAREA=1.
    cl  = cr * h * 0.01f0                   # CL = CR*H*0.01 (CR is crown-ratio percent)
    ba1 = barea + 1f0
    cw = if     eqn == "11903"; _em_r1(1.0405f0,1.2799f0,0.11941f0,0.42745f0,0f0,-0.07182f0, 1f0,35f0, d,h,cl,barea)
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
    # --- CI-added national codes (WH/AF R6-model-2, WJ R6-model-2 no-elev, MC Bechtold-m2) ---
    elseif eqn == "26305";      _em_r6m2(6.0384f0,0.51581f0,-0.21349f0,0.17468f0,0.06143f0,-0.00571f0, d,h,cl,ba1,el,1f0,72f0,54f0)
    elseif eqn == "01905";      _em_r6m2(5.8827f0,0.51479f0,-0.21501f0,0.17916f0,0.03277f0,-0.00828f0, d,h,cl,ba1,el,10f0,85f0,30f0)
    elseif eqn == "06405";      _em_r6m2(5.1486f0,0.73636f0,-0.46927f0,0.39114f0,-0.05429f0,0f0, d,h,cl,ba1,el,1f0,999f0,36f0)
    elseif eqn == "47502";      _em_bech2(4.0105f0,0.8611f0,0f0,0f0,-0.0431f0, d,cr,hi,-37f0,27f0,29f0,false)
    # --- TT/UT-added national codes (R6-model-2 DF/ES/LP/WF, Donnelly BI, Bechtold GO/GB) ---
    elseif eqn == "20205";      _em_r6m2(6.0227f0,0.54361f0,-0.20669f0,0.20395f0,-0.00644f0,-0.00378f0, d,h,cl,ba1,el,1f0,75f0,80f0)
    elseif eqn == "09305";      _em_r6m2(6.7575f0,0.55048f0,-0.25204f0,0.19002f0,0f0,-0.00313f0, d,h,cl,ba1,el,1f0,85f0,40f0)
    elseif eqn == "10805";      _em_r6m2(6.6941f0,0.81980f0,-0.36992f0,0.17722f0,-0.01202f0,-0.00882f0, d,h,cl,ba1,el,1f0,79f0,40f0)
    elseif eqn == "01505";      _em_r6m2(5.0312f0,0.53680f0,-0.18957f0,0.16199f0,0.04385f0,-0.00651f0, d,h,cl,ba1,el,2f0,75f0,35f0)
    elseif eqn == "31206";      _em_powf(7.5183f0,0.4461f0,30f0, d)
    elseif eqn == "81402";      _em_bech2(0.3309f0,0.8918f0,0f0,0.0510f0,0f0, d,cr,hi,-9.9f9,9.9f9,19f0,false)
    elseif eqn == "10201";      _em_bech1(7.4251f0,0.8991f0, d,25f0)
    # --- AK-added national codes (R10): R6-model-2 SF/YC/SS/RC, Donnelly RA, and the '08' form ---
    elseif eqn == "01105";      _em_r6m2(4.4799f0,0.45976f0,-0.10425f0,0.11866f0,0.06762f0,-0.00715f0, d,h,cl,ba1,el,4f0,72f0,33f0)
    elseif eqn == "04205";      _em_r6m2(3.3756f0,0.45445f0,-0.11523f0,0.22547f0,0.08756f0,-0.00894f0, d,h,cl,ba1,el,16f0,62f0,59f0)
    elseif eqn == "09805";      _em_r6m2(8.48f0,0.70692f0,-0.38812f0,0.17127f0,0f0,0f0, d,h,cl,ba1,el,1f0,999f0,50f0)
    elseif eqn == "24205";      _em_r6m2(6.2382f0,0.29517f0,-0.10673f0,0.23219f0,0.05341f0,-0.00787f0, d,h,cl,ba1,el,1f0,72f0,45f0)
    elseif eqn == "35106";      _em_powf(7.0806f0,0.4771f0,35f0, d)
    elseif eqn == "09508";      _cw08(3.391358f0,0.638945f0,-0.395285f0,0.264254f0,0f0,0f0, -9.9f9,9.9f9,16f0, d,h,cl,ba1,el)
    elseif eqn == "09408";      _cw08(8.515316f0,0.630576f0,-0.867757f0,0.477791f0,0.10021f0,-0.015034f0, 1f0,85f0,40f0, d,h,cl,ba1,el)
    elseif eqn == "74708";      _cw08(0.790658f0,0.551987f0,0.446434f0,0f0,0f0,-0.048415f0, -9.9f9,9.9f9,56f0, d,h,cl,ba1,el)
    elseif eqn == "37508";      _cw08(2.725006f0,0.53601f0,0f0,0.196372f0,-0.015305f0,0f0, -9.9f9,9.9f9,53f0, d,h,cl,ba1,el)
    elseif eqn == "74608";      _cw08(2.386015f0,0.63014f0,-0.147121f0,0.274356f0,0f0,0f0, -9.9f9,9.9f9,48f0, d,h,cl,ba1,el)
    else 0f0 end
    cw < 0.5f0 && (cw = 0.5f0)              # cwcalc.f final CRWDTH clamp [0.5, 99.9]
    cw > 99.9f0 && (cw = 99.9f0)
    return cw
end
