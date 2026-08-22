# =============================================================================
# UT (Utah) FVS_TreeList crown width — national cwcalc.f, UTMAP dispatch.
#
# UTMAP (vie/cwcalc.f:262) for UT's 24 species (WB LM DF WF BS AS LP ES AF PP PI WJ GO PM
# RM UJ GB NC FC MC BI BE OS OH); UT is Region 4 ⇒ BF=1.0.  ut_cwcalc =
# _cwcalc_national(UTMAP[sp], …); all codes shared with the national dispatch (UT
# contributed 01505/81402/10201).  Note UT's PP (idx10) and OS (idx23) both map to the
# R6-model-2 '12205' (unlike CI/IE PP '12203'); UTMAP is index-keyed so this is exact.
#
# REPORTING-ONLY (FVS_TreeList CrWidth); INERT for growth/mortality/CCF/volume.
# Validated vs FVSut_clean uttl FVS_TreeList: function-level A/B bit-exact on the emitted
# species.  UT is the 17th variant past the eastern 0.5 default.
# =============================================================================
const _UT_CWMAP = ("10105","11301","20205","01505","09305","74605","10805","09305",
                   "01905","12205","10602","06405","81402","10602","06405","06405",
                   "10201","74902","74902","47502","31206","74902","12205","81402")

@inline function ut_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                           el::Float32, hi::Float32)::Float32
    (1 <= sp <= 24) || return 0.5f0
    _cwcalc_national(_UT_CWMAP[sp], d, h, cr, barea, el, hi)
end
