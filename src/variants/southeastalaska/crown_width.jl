# =============================================================================
# AK (SoutheastAlaska) FVS_TreeList crown width — national cwcalc.f, AKMAP dispatch.
#
# AKMAP (vie/cwcalc.f:91) for AK's 23 species (SF AF YC TA WS LS BE SS LP RC WH MH OS AD RA
# PB AB BA AS CW WI SU OH); AK is Region 10 (KODFOR≥1000 ⇒ BF=1.0).  ak_cwcalc =
# _cwcalc_national(AKMAP[sp], …).  AK contributed the '08' form (_cw08: a·D^b·H^ch·CL^ccl·
# (BA+1)^cba·EXP(EL)^cel, D<0.1⇒flat 0.5) via codes 09508/09408/74708/37508/74608, plus
# R6-m2 SF/YC/SS/RC (01105/04205/09805/24205) and Donnelly RA (35106); the rest shared.
#
# REPORTING-ONLY (FVS_TreeList CrWidth); INERT for growth/mortality/CCF/volume.
# Validated vs FVSak_clean aktl FVS_TreeList: function-level A/B bit-exact on the emitted
# species.  AK is the 18th variant past the eastern 0.5 default.
# =============================================================================
const _AK_CWMAP = ("01105","01905","04205","09508","09408","09408","09508","09805",
                   "10805","24205","26305","26405","09408","74708","35106","37508",
                   "37508","74708","74608","74708","74708","74708","74708")

@inline function ak_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                           el::Float32, hi::Float32)::Float32
    (1 <= sp <= 23) || return 0.5f0
    _cwcalc_national(_AK_CWMAP[sp], d, h, cr, barea, el, hi)
end
