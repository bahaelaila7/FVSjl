# =============================================================================
# centralidaho.jl — the CI (Central Idaho) variant singleton + registration.
#
# Ported from: ci/*.f (FVS "CENTRAL IDAHO PROGNOSIS"). CI is a Tier-1 discount of the
# western Wykoff-DDS cluster — closest template is IE (shared Northern-Rockies conifer set,
# ci/dgf.f is the standard western Wykoff DDS reading DGLD/DGCR/DGDBAL/DGHAB, NO GEMDG).
# Distinguishing infra (MEASURED, ci/blkdat.f + ci/grinit.f):
#   • VARACD='CI', MAXSP=19, seed 55329, IFINT=10 (10-yr cycle).
#   • LZEIDE=.TRUE.  ⇒ ZEIDE SDI (like UT/TT — NOT Stage like IE/EM/BM). ← key difference from IE.
#   • DGSD=1.7 ; LHTDRG=.TRUE. for all except sp15 (MC).
#   • Habitat: ci/habtyp.f ITYPE indexes the NI (Northern Idaho) 30 habitat types → OCURHT(16,MAXSP)
#     habitat-group dimension (same shape as IE/KT).
# Species JSP (19): WP WL DF GF WH RC LP ES AF PP WB PY AS WJ MC LM CW OS OH
# FIAJSP:          119 073 202 017 263 242 108 093 019 122 101 231 746 064 475 113 747 299 998
# Oracle: /workspace/.ciwork/FVSci_clean (relink from bin/FVSci_buildDir/*.o + isoc23 shim).
# Canonical test stand: cit01 (tests/FVSci/cit01.key; forest 412, habitat 520).
#
# PORT IN PROGRESS — chunk 0 scaffold. Growth/volume hooks have NO CI method yet, so dispatching
# on CentralIdaho() errors loudly (doctrine #5) until each chunk lands + is validated vs FVSci_clean.
# =============================================================================

"""
    CentralIdaho <: AbstractVariant

The FVS Central Idaho variant (VARACD "CI", MAXSP = 19) — a Tier-1 discount of the western
Wykoff cluster (IE template + UT-style Zeide SDI). PORT IN PROGRESS (chunk 0 scaffold).
"""
struct CentralIdaho <: AbstractVariant end

variant_code(::CentralIdaho) = "CI"
nspecies(::CentralIdaho) = 19
htg_period(::CentralIdaho) = 10f0
