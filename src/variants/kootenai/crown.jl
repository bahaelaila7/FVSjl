# =============================================================================
# crown.jl (kootenai) — KT per-tree CCF (kt/ccfcal.f MODE=1). The stand CCF = Σ CCFT is the RELDEN
# that the DG (dgf!) and height (htgf) chunks read. UNLIKE the eastern/CR stand_ccf (crown-width → area),
# KT's CCF is a DIRECT per-species polynomial (Bush/Crookston, kt/ccfcal.f):
#   D >= 10:  CCFT = RD1(sp) + D*RD2(sp) + D^2*RD3(sp)
#   D <  10:  CCFT = RDA(sp) * D^RDB(sp)
#   CCFT *= P (tree TPA);  stand CCF = Σ CCFT.
# Crown-WIDTH (ccfcal MODE=2, B1..B6) is the crown-ratio/crown model — separate; lands with crown_ratio!.
# =============================================================================

const KT_RD1 = Float32[0.03,0.02,0.11,0.04,0.03,0.03,0.01925,0.03,0.03,0.03,0.03]
const KT_RD2 = Float32[0.0167,0.0148,0.0333,0.0270,0.0215,0.0238,0.01676,0.0173,0.0216,0.0180,0.0215]
const KT_RD3 = Float32[0.00230,0.00338,0.00259,0.00405,0.00363,0.00490,0.00365,0.00259,0.00405,0.00281,0.00363]
const KT_RDA = Float32[0.009884,0.007244,0.017299,0.015248,0.011109,0.008915,0.009187,0.007875,0.011402,0.007813,0.011109]
const KT_RDB = Float32[1.6667,1.8182,1.5571,1.7333,1.7250,1.7800,1.7600,1.7360,1.7560,1.7680,1.7250]

"""Per-tree CCF contribution (kt/ccfcal.f MODE=1), excluding the ×TPA factor."""
@inline function kt_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    if d >= 10f0
        return KT_RD1[sp] + d * KT_RD2[sp] + d * d * KT_RD3[sp]
    else
        return KT_RDA[sp] * d ^ KT_RDB[sp]
    end
end
