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

"forest number → (latitude, longitude, elevation·100ft) defaults (forkod.f:195+)."
const SN_FOREST_LOCATION = Dict{Int,NTuple{3,Float32}}(
    801 => (32.37f0, 86.30f0,  7f0), 802 => (37.99f0, 84.18f0, 12f0),
    803 => (34.30f0, 83.82f0, 17f0), 804 => (35.16f0, 84.88f0, 22f0),
    805 => (30.44f0, 84.28f0,  1f0), 806 => (31.32f0, 92.43f0,  2f0),
    807 => (33.31f0, 89.17f0,  3f0), 808 => (37.27f0, 79.94f0, 21f0),
    809 => (34.50f0, 93.06f0,  9f0), 810 => (35.28f0, 93.13f0, 13f0),
    811 => (35.60f0, 82.55f0, 25f0), 812 => (34.00f0, 81.04f0,  4f0),
    813 => (31.34f0, 94.73f0,  3f0), 701 => (35.60f0, 82.55f0, 25f0),
    905 => (37.95f0, 91.77f0, 10f0), 908 => (37.74f0, 88.54f0,  4f0),
)

"Default (lat, long, elev) for a forest number; (0,0,0) if unknown."
forest_location(forest::Integer) = get(SN_FOREST_LOCATION, Int(forest), (0f0, 0f0, 0f0))
