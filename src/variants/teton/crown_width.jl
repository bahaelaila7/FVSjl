# =============================================================================
# TT (Teton) FVS_TreeList crown width — national cwcalc.f, TTMAP dispatch.
#
# TTMAP (vie/cwcalc.f:253) for TT's 18 species (WB LM DF PM BS AS LP ES AF PP UJ RM BI MM
# NC MC OS OH); TT is Region 4 ⇒ BF=1.0.  tt_cwcalc = _cwcalc_national(TTMAP[sp], …), all
# codes shared with the national dispatch (TT contributed 20205/09305/10805/31206).
#
# REPORTING-ONLY (FVS_TreeList CrWidth); INERT for growth/mortality/CCF/volume.
# Validated vs FVStt_clean tttl FVS_TreeList: function-level A/B bit-exact on the emitted
# species.  TT is the 16th variant past the eastern 0.5 default.
# =============================================================================
const _TT_CWMAP = ("10105","11301","20205","10602","09305","74605","10805","09305",
                   "01905","12203","06405","06405","31206","32102","74902","47502",
                   "12205","74902")

@inline function tt_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                           el::Float32, hi::Float32)::Float32
    (1 <= sp <= 18) || return 0.5f0
    _cwcalc_national(_TT_CWMAP[sp], d, h, cr, barea, el, hi)
end
