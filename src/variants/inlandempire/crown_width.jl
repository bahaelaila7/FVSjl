# =============================================================================
# IE (Inland Empire) FVS_TreeList crown width — national cwcalc.f, IEMAP dispatch.
#
# cwcalc.f is the shared NATIONAL crown-width routine: one flat SELECT CASE(CWEQN)
# keyed by a 5-char code (3-char FIA sp + 2-char model), and each variant supplies a
# per-species map (IEMAP for IE) that picks the code.  IE is Region 1/4 so KODFOR<601
# ⇒ the forest bias-factor section is skipped and BF=1.0 (vie/cwcalc.f:477 GO TO 10).
#
# ie_cwcalc = _cwcalc_national(IEMAP[sp], …); the code→equation SELECT CASE lives in the
# shared national_crown_width.jl (reuses the EM-353/353-validated national form helpers).
# IE contributed the 8 codes EM does not carry (11903 01703 26303 24203 10602 23104 32102
# 12205) to that shared dispatch.
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

@inline function ie_cwcalc(sp::Int, d::Float32, h::Float32, cr::Float32, barea::Float32,
                           el::Float32, hi::Float32)::Float32
    (1 <= sp <= 23) || return 0.5f0
    _cwcalc_national(_IE_CWMAP[sp], d, h, cr, barea, el, hi)
end
