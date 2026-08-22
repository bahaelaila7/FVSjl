# =============================================================================
# CI (CentralIdaho) FVS_TreeList crown width — national cwcalc.f, CIMAP dispatch.
#
# CIMAP (vie/cwcalc.f:129) for CI's 19 species (WP WL DF GF WH RC LP ES AF PP WB PY AS WJ
# MC LM CW OS OH); CI is Region 4 ⇒ BF=1.0.  ci_cwcalc = _cwcalc_national(CIMAP[sp], …).
# CI contributed 4 codes to the shared national dispatch: 26305 (WH R6-m2), 01905 (AF
# R6-m2), 06405 (WJ R6-m2 no-elev), 47502 (MC Bechtold-m2); the rest are shared with
# EM/IE.  Note CIMAP disambiguates two FIA-122 species by INDEX (PP idx10 '12203' log vs
# OS idx18 '12205' R6-m2).
#
# REPORTING-ONLY (FVS_TreeList CrWidth); INERT for growth/mortality/CCF/volume.
# Validated vs FVSci_clean citl FVS_TreeList: function-level A/B bit-exact on the emitted
# species.  CI is the 15th variant past the eastern 0.5 default.
# =============================================================================
const _CI_CWMAP = ("11903","07303","20203","01703","26305","24203","10803","09303",
                   "01905","12203","10105","23104","74605","06405","47502","11301",
                   "74902","12205","74902")

@inline function ci_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                           el::Float32, hi::Float32)::Float32
    (1 <= sp <= 19) || return 0.5f0
    _cwcalc_national(_CI_CWMAP[sp], d, h, cr, barea, el, hi)
end
