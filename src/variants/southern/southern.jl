# =============================================================================
# southern.jl — the Southern (SN) variant
#
# Ported from: sn/  (sn/blkdat.f, sn/dgf.f, sn/htgf.f, sn/morts.f, ...)
#
# Declares the variant singleton and its identity. The growth/mortality/regen
# methods that specialize the engine hooks (from src/variants/variant.jl) are
# added in later chunks under this directory (diameter_growth.jl, etc.).
# =============================================================================

"""
    Southern <: AbstractVariant

The FVS Southern variant (variant code "SN"), 90 species. Singleton — carried as
the type parameter of `StandState` so all variant dispatch is devirtualized.
"""
struct Southern <: AbstractVariant end

variant_code(::Southern) = "SN"
nspecies(::Southern) = 90

# Keywords inert *for SN specifically* — verified against SN code, NOT guaranteed for other
# variants (re-check the effect variable when porting one). Treated as no-ops by the keyword
# dispatch alongside the variant-agnostic `KNOWN_NOOP`.
#   FIXCW  (sn/cwidth.f): crown-width override — verified .sum-inert in SN (CRWDTH is output-
#          only; a live-Fortran SN FIXCW run is byte-identical). Other variants may feed crown
#          width into cover/growth.
#   FIAVBC: switches volume/biomass to the FIA National Volume Library — FVSjl has only the R8
#          Clark equations (SN default LFIANVB=.FALSE.). A variant on the NVB path would diverge.
#   NOAUTOES / AUTOES / NOHTDREG: establishment-control flags — establishment is variant-specific;
#          their inertness is a current-state assumption (AUTOES would need the auto-establishment
#          trigger ported — a latent gap of the MANAGED kind).
const KNOWN_NOOP_SN = Set([
    "FIXCW", "FIAVBC", "NOAUTOES", "AUTOES", "NOHTDREG",
])

# Insect/disease EXTENSIONS that FVS itself STUBS OUT of the SN variant — the SN build links only the
# `ex*.f` no-op stubs (exmpb/exdfb/exdftm/exmist/exbm/exrd …), and the keyword handler emits
# `FVS11 ERROR: REQUESTED EXTENSION IS NOT PART OF THIS PROGRAM` then leaves the simulation unchanged
# (verified: DFB on snt01 ⇒ byte-identical .sum). These are NOT unported omissions — they are not part
# of SN FVS at all (mostly western insects/diseases with no Southern hosts), so the faithful SN
# behaviour is exactly this: recognize, do nothing. (Block keywords like MISTOE/RDIN read sub-keywords
# in full FVS; the stub doesn't consume them, and FVSjl likewise ignores both the keyword and its
# orphaned sub-keyword lines — same nil .sum.)
const SN_UNSUPPORTED_EXTENSIONS = Set([
    "MPB", "DFB", "DFTM", "WSBW", "BWE", "BRUST", "MISTOE", "RDIN", "ANIN", "RRIN",
    "RDBBMORT", "RDSUM", "RDDETAIL",   # root-disease report sub-keywords (also ex* in SN)
])

const _SN_NOOP_ALL = union(KNOWN_NOOP_SN, SN_UNSUPPORTED_EXTENSIONS)
variant_noop_keywords(::Southern) = _SN_NOOP_ALL
