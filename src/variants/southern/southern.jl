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
