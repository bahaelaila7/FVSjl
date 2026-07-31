# =============================================================================
# inlandempire.jl — the IE (Inland Empire) variant singleton + registration.
#
# Ported from: ie/*.f (FVS "INLAND EMPIRE PROGNOSIS"). IE is the 3rd western variant
# (CR hub + KT done) and REUSES the KT western-Wykoff-DDS engine at a discount:
#   • MORTALITY = IDENTICAL Hamilton (ie/morts.f:278 RIP=2.76253+0.222310·√D−…) — reuse
#     the KT mortality! form with IE coefficients (POT/IPDG/PMSC).
#   • CROWN / REGENT / VOLUME (Flewelling FW2, geocode 'I'/INGY already in _fw2_jsp) reuse.
#   • DG (ie/dgf.f) + HEIGHT (ie/htgf.f) share the Wykoff family but have LARGER equations —
#     the IE DDS term set + height form MUST be VERIFIED vs live FVSie (instrument-replay).
#
# IE infra (ie/blkdat.f, ie/PRGPRM.F77): VARACD='IE', MAXSP=23, IFINT=10 (YR=10), seed 55329,
# habitat OCURHT(16,MAXSP). Species JSP (23): WP WL DF GF WH RC LP ES AF PP MH  WB LM LL PM RM
# PY AS CO MM PB OH OS (first 11 ≈ KT; +12 more). Oracle: /workspace/.iework/FVSie_clean
# (relink_ie.sh mirrors relink_kt.sh). PORT IN PROGRESS — hooks error loudly until implemented (doctrine #5).
# =============================================================================

"""
    InlandEmpire <: AbstractVariant

The FVS Inland Empire variant (VARACD "IE", MAXSP = 23) — 3rd western variant, reusing the
KT engine at a discount. PORT IN PROGRESS (chunk 0 scaffold).
"""
struct InlandEmpire <: AbstractVariant end

variant_code(::InlandEmpire) = "IE"
nspecies(::InlandEmpire) = 23
htg_period(::InlandEmpire) = 10f0

const IE_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "inlandempire"))

coefficients(::InlandEmpire) = cached_coefficients(() -> load_species_coefficients(IE_DATADIR), "IE")
