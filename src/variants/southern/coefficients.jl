# =============================================================================
# coefficients.jl (southern) — locate and cache the Southern variant's CSV data.
#
# The generic container + loader live in core/coefficients.jl; here we just point
# at the Southern data directory and cache the load.
# =============================================================================

"Directory holding the Southern variant's coefficient CSVs."
const SN_DATADIR = normpath(joinpath(@__DIR__, "..", "..", "..", "data", "southern"))

"""
    coefficients(::Southern) -> SpeciesCoefficients

Cached Southern coefficient load (first call reads `data/southern/*.csv`).
"""
coefficients(::Southern) = cached_coefficients(() -> load_species_coefficients(SN_DATADIR), "SN")
