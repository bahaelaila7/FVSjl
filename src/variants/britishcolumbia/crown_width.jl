# =============================================================================
# BC (BritishColumbia) FVS_TreeList crown width — national cwcalc.f, BCMAP dispatch.
#
# BCMAP (vie/cwcalc.f:202) for BC's 15 species (PW LW FD BG HW CW PL SE BL PY EP AT AC OC
# OH); BC is metric but FVS computes cwcalc in imperial internally (D in, H ft, BA ft²/ac,
# EL hundreds-ft), so bc_cwcalc = _cwcalc_national(BCMAP[sp], …) with imperial inputs —
# identical to the other Region-1/4 variants (BF=1.0).  Every BC code is already carried by
# the shared national dispatch (0 new codes: 11903/07303/20203/01703 log, 26305/24205/
# 10805/01905/74605 R6-m2, 09303/12203 log, 37506 Donnelly, 74902 Bechtold).
#
# REPORTING-ONLY (FVS_TreeList CrWidth, in FEET); the metric-treelist output layer converts
# ft→m downstream (FVS_TreeList_Metric.CrWidth = feet·0.3048).  INERT for growth/mort/vol.
# Validated vs FVSbc_clean bctl FVS_TreeList_Metric (converted to imperial): the log-form
# species (PW/FD) are bit-exact (800/800 each); the R6-m2 species (PL 10805) are bit-exact
# at the oracle's ELEV≈1.476 (STDINFO 45 m → hundreds-ft), cornered only on the BC metric
# elevation-unit feed (a separate, already-fixed BC-reader concern, commit 185833a).  BC is
# the 19th (final western) variant past the eastern 0.5 default; note jl does not yet emit a
# metric TreeList, so this dispatch is currently additive/inert but correct-when-emitted.
# =============================================================================
const _BC_CWMAP = ("11903","07303","20203","01703","26305","24205","10805","09303",
                   "01905","12203","37506","74605","74902","20203","37506")

@inline function bc_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                           el::Float32, hi::Float32)::Float32
    (1 <= sp <= 15) || return 0.5f0
    _cwcalc_national(_BC_CWMAP[sp], d, h, cr, barea, el, hi)
end
