# =============================================================================
# forest_location.jl — default stand lat/long/elevation by national forest
#
# Ported from: sn/forkod.f (phase 3, the forest → location table).
#
# When a stand doesn't supply latitude/longitude (no LOCATE keyword), FVS fills
# defaults from the national forest code. These feed the Hopkins bioclimatic index
# used by several crown-width and growth equations, so getting them right matters.
# The forest number is KODFOR ÷ 100 (e.g. 80106 → 801).
# =============================================================================

# The forest → (latitude, longitude, elevation·100ft) table (forkod.f:195+) is loaded
# from data/southern/forest_locations.csv into `coef.forest_location`.

"Default (lat, long, elev) for a forest number; (0,0,0) if unknown."
forest_location(c::SpeciesCoefficients, forest::Integer) =
    get(c.forest_location, Int(forest), (0f0, 0f0, 0f0))
