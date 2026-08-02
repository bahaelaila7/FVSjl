# =============================================================================
# crown.jl (teton) — TT per-tree CCF (tt/ccfcal.f MODE=1). Chunk 3 prerequisite
# (RELDEN = stand CCF + PCCF feed the DG), also reused by chunk 5 (crown).
#
# tt/ccfcal.f: BREAK=10 for PP(10)/NC(15)/OH(18), else 1. Per-tree CCF (before ×P):
#   D ≥ BREAK           : RD1 + D·RD2 + D²·RD3                          (large-tree poly)
#   0.1 < D < BREAK     : sp13/16 → D·(RD1+RD2+RD3) ; else RDA·D^RDB   (small-tree)
#   D ≤ 0.1             : sp13/16 → D·(RD1+RD2+RD3) ; sp10 → RDA·D^RDB ; else 0.001
# Returns the per-tree CCF (the caller multiplies by tpa for the stand/point CCF sum),
# matching em_tree_ccf. Verified vs raw tt/ccfcal.f DATA (WB RD1=.01925, DF RD1=.11).
# =============================================================================

const TT_RD1 = Float32[0.01925, 0.01925, 0.11, 0.01925, 0.03, 0.03, 0.01925, 0.03, 0.03, 0.03, 0.01925, 0.01925, 0.0204, 0.03, 0.03, 0.0204, 0.01925, 0.03]
const TT_RD2 = Float32[0.0168, 0.0168, 0.0333, 0.01676, 0.0173, 0.0238, 0.0168, 0.0173, 0.0216, 0.018, 0.01676, 0.01676, 0.0246, 0.0238, 0.0215, 0.0246, 0.0168, 0.0215]
const TT_RD3 = Float32[0.00365, 0.00365, 0.00259, 0.00365, 0.00259, 0.0049, 0.00365, 0.00259, 0.00405, 0.00281, 0.00365, 0.00365, 0.0074, 0.0049, 0.00363, 0.0074, 0.00365, 0.00363]
const TT_RDA = Float32[0.009187, 0.009187, 0.017299, 0.009187, 0.007875, 0.008915, 0.009187, 0.007875, 0.011402, 0.007813, 0.009187, 0.009187, 0.0, 0.008915, 0.011109, 0.0, 0.009187, 0.011109]
const TT_RDB = Float32[1.76, 1.76, 1.5571, 1.76, 1.736, 1.78, 1.76, 1.736, 1.756, 1.768, 1.76, 1.76, 0.0, 1.78, 1.725, 0.0, 1.76, 1.725]

@inline function tt_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    poly()  = TT_RD1[sp] + d * TT_RD2[sp] + d * d * TT_RD3[sp]
    small() = TT_RDA[sp] * d ^ TT_RDB[sp]
    brk = (sp == 10 || sp == 15 || sp == 18) ? 10f0 : 1f0
    if d >= brk
        return poly()
    elseif d > 0.1f0
        return (sp == 13 || sp == 16) ? d * (TT_RD1[sp] + TT_RD2[sp] + TT_RD3[sp]) : small()
    else
        (sp == 13 || sp == 16) && return d * (TT_RD1[sp] + TT_RD2[sp] + TT_RD3[sp])
        sp == 10 && return small()
        return 0.001f0
    end
end
