# =============================================================================
# crown.jl (easternmontana) — per-tree CCF (em/ccfcal.f MODE=1).
#
# EM CCF is the direct per-species polynomial (Paine-Hann / NI INT-133 form), summed
# over the tree list = stand CCF (like KT/IE ccfcal, not the CR crown-width→area path).
#   large:  CCFT = RD1 + D·RD2 + D²·RD3
#   small:  CCFT = RDA · D^RDB
# The large/small DBH threshold is species-group dependent (em/ccfcal.f SELECT CASE):
#   sp 5 (LL):                 D≥10 → poly, else RDA·D^RDB           (no 0.001 floor)
#   sp 1-4,6-10,12,17,18:      D≥1  → poly; D>0.1 → RDA·D^RDB; else 0.001
#   sp 11,13-16,19:            D≥10 → poly; D>0.1 → RDA·D^RDB; else 0.001
# Coefficients transcribed verbatim from em/ccfcal.f DATA (species 1=WB..19=OH).
# =============================================================================

const EM_RD1 = Float32[0.0186, 0.0392, 0.0388, 0.01925, 0.03, 0.01925, 0.01925, 0.03, 0.0172,
                       0.0219, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.0204, 0.03]
const EM_RD2 = Float32[0.0146, 0.0180, 0.0269, 0.01676, 0.0216, 0.01676, 0.01676, 0.0173, 0.00876,
                       0.0169, 0.0215, 0.0238, 0.0215, 0.0215, 0.0215, 0.0215, 0.0238, 0.0246, 0.0215]
const EM_RD3 = Float32[0.00288, 0.00207, 0.00466, 0.00365, 0.00405, 0.00365, 0.00365, 0.00259, 0.00112,
                       0.00325, 0.00363, 0.00490, 0.00363, 0.00363, 0.00363, 0.00363, 0.00490, 0.0074, 0.00363]
const EM_RDA = Float32[0.009884, 0.007244, 0.017299, 0.009187, 0.011402, 0.009187, 0.009187, 0.007875,
                       0.011402, 0.007813, 0.011109, 0.008915, 0.011109, 0.011109, 0.011109, 0.011109,
                       0.008915, 0.011109, 0.011109]
const EM_RDB = Float32[1.6667, 1.8182, 1.5571, 1.7600, 1.7560, 1.7600, 1.7600, 1.7360, 1.7560,
                       1.7780, 1.7250, 1.7800, 1.7250, 1.7250, 1.7250, 1.7250, 1.7800, 1.7250, 1.7250]

# em/ccfcal.f MODE=1: per-tree CCF (before the ×P expansion the caller applies).
@inline function em_tree_ccf(sp::Integer, d::Real)::Float32
    d <= 0f0 && return 0f0
    poly()  = EM_RD1[sp] + d * EM_RD2[sp] + d * d * EM_RD3[sp]
    small() = EM_RDA[sp] * d ^ EM_RDB[sp]
    if sp == 5                                   # LL: D≥10 poly, else small (no floor)
        return d >= 10f0 ? poly() : small()
    elseif sp == 11 || (13 <= sp <= 16) || sp == 19   # GA/CW/BA/PW/NC/OH: D≥10 threshold
        return d >= 10f0 ? poly() : (d > 0.1f0 ? small() : 0.001f0)
    else                                          # 1-4,6-10,12,17,18: D≥1 threshold
        return d >= 1f0 ? poly() : (d > 0.1f0 ? small() : 0.001f0)
    end
end
