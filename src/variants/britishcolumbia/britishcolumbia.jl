# =============================================================================
# britishcolumbia.jl — the BC (British Columbia) variant singleton + registration.
#
# Ported from: bc/*.f (FVS "BRITISH COLUMBIA"). BC is a Tier-1 discount of the western Wykoff-DDS
# cluster — closest template **IE** (Stage SDI, standard western Wykoff dgf reading DGLD; even
# cleaner than CI which is Zeide). Distinguishing infra (MEASURED, bc/blkdat.f + bc/grinit.f):
#   • VARACD='BC', MAXSP=15, seed 55329, IFINT=10, DGSD=2.0.
#   • LZEIDE=.FALSE. ⇒ STAGE SDI (like IE/EM/BM).  LHTDRG=.TRUE. all species.
#   • Bark: bc/bratio.f BARK2=0, IMAP=2 for all ⇒ constant BRATIO=BARK1 (IE-style, no branches).
# Species JSP (15): PW LW FD BG HW CW PL SE BL PY EP AT AC OC OH (BC codes; first 10 FIA = IE conifers).
# FIAJSP: 119 073 202 017 263 242 108 093 019 122 375 746 747 202 998
# Oracle: /workspace/.bcwork/FVSbc_clean (relink_bc.sh + bc_stubs.o for 4 missing DBS output symbols).
#
# PORT IN PROGRESS — chunk 0 scaffold. Un-ported hooks error loudly (doctrine #5).
# =============================================================================

"""
    BritishColumbia <: AbstractVariant

The FVS British Columbia variant (VARACD "BC", MAXSP = 15) — a Tier-1 discount of the western
Wykoff cluster (IE template + Stage SDI). PORT IN PROGRESS (chunk 0 scaffold).
"""
struct BritishColumbia <: AbstractVariant end

variant_code(::BritishColumbia) = "BC"
nspecies(::BritishColumbia) = 15
htg_period(::BritishColumbia) = 10f0

const BC_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "britishcolumbia"))

coefficients(::BritishColumbia) = cached_coefficients(() -> load_species_coefficients(BC_DATADIR), "BC")
