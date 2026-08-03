# =============================================================================
# crown.jl (centralidaho) — CI per-tree CCF (ci/ccfcal.f MODE=1). Chunk 5 (CCF portion).
# The stand CCF = Σ CCFT·TPA is the RELDEN the DG (dgf!) CONSPP term reads. Same western
# polynomial as KT/IE (Bush/Crookston): D≥10 → RD1+D·RD2+D²·RD3; D<10 → RDA·D^RDB. CI's RD*
# coefficients (ci/ccfcal.f DATA) match KT for the shared N-Rockies conifers.
# (Crown-WIDTH MODE=2 B1..B6 + crown-ratio model land separately with crown_ratio_update!.)
# =============================================================================

const CI_RD1 = Float32[0.03,0.02,0.11,0.04,0.03,0.03,0.01925,0.03,0.03,0.03,0.01925,0.01925,0.03,0.01925,0.0204,0.01925,0.03,0.03,0.03]
const CI_RD2 = Float32[0.0167,0.0148,0.0333,0.027,0.0215,0.0238,0.01676,0.0173,0.0216,0.018,0.01676,0.01676,0.0238,0.01676,0.0246,0.01676,0.0215,0.0215,0.0215]
const CI_RD3 = Float32[0.0023,0.00338,0.00259,0.00405,0.00363,0.0049,0.00365,0.00259,0.00405,0.00281,0.00365,0.00365,0.0049,0.00365,0.0074,0.00365,0.00363,0.00363,0.00363]
const CI_RDA = Float32[0.009884,0.007244,0.017299,0.015248,0.011109,0.008915,0.009187,0.007875,0.011402,0.007813,0.009187,0.009187,0.008915,0.009187,0.0,0.009187,0.011109,0.011109,0.011109]
const CI_RDB = Float32[1.6667,1.8182,1.5571,1.7333,1.725,1.78,1.76,1.736,1.756,1.768,1.76,1.76,1.78,1.76,0.0,1.76,1.725,1.725,1.725]

"""Per-tree CCF contribution (ci/ccfcal.f MODE=1), excluding the ×TPA factor."""
@inline function ci_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    if d >= 10f0
        return CI_RD1[sp] + d * CI_RD2[sp] + d * d * CI_RD3[sp]
    else
        return CI_RDA[sp] * Float32(d) ^ CI_RDB[sp]
    end
end
