# =============================================================================
# crown.jl (utah) — per-tree CCF (ut/ccfcal.f MODE=1). Chunk 5 (partial: the CCF path
# that feeds dgf!/htgf RELDEN + PCCF; crown-ratio dub + crown-width come later).
#
# UT CCF is the direct per-species polynomial (NI INT-133 / Paine-Hann form), summed over
# the tree list = stand CCF (RELDEN) + per-point PCCF (like EM/KT/IE/TT, not the CR crown-
# width→area path):
#   large (D≥BREAK):  CCFT = RD1 + D·RD2 + D²·RD3
#   small (D>0.1):    sp 20/21 (MC/BI): D·(RD1+RD2+RD3);  else RDA·D^RDB
#   tiny  (D≤0.1):    sp 20/21: D·(RD1+RD2+RD3);          else 0.001
#   BREAK = 10 for sp 17:19,22 (GB/NC/FC/BE); else 1.
# Coefficients verbatim from ut/ccfcal.f DATA (species 1=WB..24=OH).
# =============================================================================

const UT_RD1 = Float32[0.01925, 0.01925, 0.11, 0.04, 0.03, 0.03, 0.01925, 0.03, 0.03, 0.03,
                       0.01925, 0.01925, 0.03, 0.01925, 0.01925, 0.01925, 0.01925, 0.03, 0.03, 0.0204,
                       0.0204, 0.03, 0.01925, 0.03]
const UT_RD2 = Float32[0.01676, 0.01676, 0.0333, 0.0270, 0.0173, 0.0238, 0.01676, 0.0173, 0.0216, 0.0180,
                       0.01676, 0.01676, 0.0215, 0.01676, 0.01676, 0.01676, 0.01676, 0.0215, 0.0215, 0.0246,
                       0.0246, 0.0215, 0.01676, 0.0215]
const UT_RD3 = Float32[0.00365, 0.00365, 0.00259, 0.00405, 0.00259, 0.00490, 0.00365, 0.00259, 0.00405, 0.00281,
                       0.00365, 0.00365, 0.00363, 0.00365, 0.00365, 0.00365, 0.00365, 0.00363, 0.00363, 0.0074,
                       0.0074, 0.00363, 0.00365, 0.00363]
const UT_RDA = Float32[0.009187, 0.009187, 0.017299, 0.015248, 0.007875, 0.008915, 0.009187, 0.007875, 0.011402, 0.007813,
                       0.009187, 0.009187, 0.011109, 0.009187, 0.009187, 0.009187, 0.009187, 0.011109, 0.011109, 0.0,
                       0.0, 0.011109, 0.009187, 0.011109]
const UT_RDB = Float32[1.7600, 1.7600, 1.5571, 1.7333, 1.7360, 1.7800, 1.7600, 1.7360, 1.7560, 1.7680,
                       1.7600, 1.7600, 1.7250, 1.7600, 1.7600, 1.7600, 1.7600, 1.7250, 1.7250, 0.0,
                       0.0, 1.7250, 1.7600, 1.7250]

# ut/ccfcal.f MODE=1: per-tree CCF (before the ×P expansion the caller applies).
@inline function ut_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    dd = Float32(d)
    brk = (17 <= sp <= 19 || sp == 22) ? 10f0 : 1f0
    if dd >= brk
        return UT_RD1[sp] + dd * UT_RD2[sp] + dd * dd * UT_RD3[sp]
    elseif dd > 0.1f0
        return (sp == 20 || sp == 21) ? dd * (UT_RD1[sp] + UT_RD2[sp] + UT_RD3[sp]) : UT_RDA[sp] * dd ^ UT_RDB[sp]
    else
        return (sp == 20 || sp == 21) ? dd * (UT_RD1[sp] + UT_RD2[sp] + UT_RD3[sp]) : 0.001f0
    end
end
