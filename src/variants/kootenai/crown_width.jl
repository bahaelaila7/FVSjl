# =============================================================================
# KT (Kootenai) FVS_TreeList crown width — national cwcalc.f, KTMAP dispatch.
#
# KTMAP (vie/cwcalc.f:211) is IEMAP's first 11 codes verbatim (KT MAXSP=11: WP WL DF GF
# WH RC LP ES AF PP MH); KT is Region 1 ⇒ BF=1.0.  Every KT code is already carried by
# ie_cwcalc's dispatch (which reuses the EM-353/353-validated national form helpers), so
# kt_cwcalc delegates to ie_cwcalc (identical codes, same forms, BF=1.0).
#
# REPORTING-ONLY (FVS_TreeList CrWidth); INERT for growth/mortality/CCF/volume.
# Validated vs FVSkt_clean kttl FVS_TreeList: function-level A/B bit-exact on all emitted
# species (WP/LP/ES log-forms + MH/OT R6M2 26405).  14th variant past the eastern 0.5.
# =============================================================================
const _KT_CWMAP = ("11903","07303","20203","01703","26303","24203","10803","09303",
                   "01903","12203","26405")

@inline function kt_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                           el::Float32, hi::Float32)::Float32
    (1 <= sp <= 11) || return 0.5f0
    ie_cwcalc(sp, d, h, cr, barea, el, hi)   # KTMAP[sp] == IEMAP[sp] for sp 1..11
end
