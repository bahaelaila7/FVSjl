# =============================================================================
# units.jl — unit-conversion constants
#
# Ported from: common/METRIC.F77
#
# FVS is an imperial-unit model internally; these factors are used only at the
# I/O edges when metric output is requested (the per-stand `metric` flag lives on
# Control, not here — these are pure constants).
# =============================================================================

# Metric → US
const CM_TO_IN  = 0.3937f0
const CM_TO_FT  = 0.0328084f0
const M_TO_IN   = 39.37f0
const M_TO_FT   = 3.28084f0
const KM_TO_MI  = 0.6214f0
const M2_TO_FT2 = 10.763867f0
const HA_TO_ACRE = 2.471f0
const M3_TO_FT3 = 35.314455f0
const KG_TO_LB  = 2.2046226f0
const TONNE_TO_TON = 1.102311f0
const KJ_TO_BTU = 0.9478171f0

# US → Metric
const IN_TO_CM  = 2.54f0
const FT_TO_CM  = 30.48f0
const IN_TO_M   = 0.0254001f0
const FT_TO_M   = 0.3048f0
const MI_TO_KM  = 1.609f0
const FT2_TO_M2 = 0.0929034f0
const ACRE_TO_HA = 0.4046945f0
const FT3_TO_M3 = 0.028317f0
const LB_TO_KG  = 0.4535924f0
const TON_TO_TONNE = 0.90718f0
const BTU_TO_KJ = 1.0550559f0

# Composite
const M2PHA_TO_FT2PACRE = 4.3560773f0

@inline celsius_to_f(c::Real) = 1.8f0 * Float32(c) + 32.0f0
@inline f_to_celsius(f::Real) = 0.554f0 * Float32(f) - 17.7f0
