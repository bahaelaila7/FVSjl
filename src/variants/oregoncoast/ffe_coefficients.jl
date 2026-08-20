# =============================================================================
# ffe_coefficients.jl (oregoncoast) — OC FFE (Fire & Fuels) per-species coefficients.
#
# The ORGANON variants (OC/OP) ported growth+volume only; their coef.species lacks the FFE columns jl's
# shared FFE reads (fmcba/crown_biomass/fuel). This is the OC crown-biomass chunk (extracted 2026-08-20
# from bin/FVSoc_buildDir), the first of the NC-FFE-style OC/OP FFE port (see scratchpad/ocffe/FFE_PORT_HANDOFF.md).
#
# OC_FFE_V2T   — fmvinit.f V2T (crown-biomass, LB/CUFT by species, from NI/SO variants). jl's crown_biomass
#                multiplies by _FM_P2T=0.0005 (lb→tons, =1/2000), so these are the RAW V2T values.
# OC_FFE_ISPMAP — fmcrow.f DATA ISPMAP (SPIW→SPIE): OC species → FFE crown-biomass group (fmcrow.f:151).
# OC_FFE_DBHMIN — the crown-biomass merch DBH gate = OC's BLM merch min (sp_dbh_min = 7.0).
# =============================================================================

# fmvinit.f SELECT CASE (I) → V2T(I), OC species 1..50 (LB/CUFT).
const OC_FFE_V2T = Float32[
    24.3, 21.8, 19.3, 21.8, 22.5, 22.5, 28.7, 26.2, 26.2, 22.5,   # 1-10  (7=DF 28.7)
    23.7, 23.7, 23.7, 22.5, 21.2, 21.2, 22.5, 23.7, 23.7, 23.7,   # 11-20
    34.9, 20.6, 21.2, 26.2, 23.7, 49.9, 49.9, 37.4, 37.4, 37.4,   # 21-30
    34.9, 37.4, 49.9, 27.4, 37.4, 23.1, 36.2, 36.2, 27.4, 31.2,   # 31-40
    31.8, 36.2, 28.7, 21.8, 19.3, 22.5, 34.9, 36.2, 34.9, 21.2]   # 41-50 (50=redwood 21.2)

# fmcrow.f DATA ISPMAP / OC species (SPIW) → FFE crown-biomass group (SPIE).
const OC_FFE_ISPMAP = Int32[
    7, 20, 7, 4, 4, 4, 3, 6, 24, 14,
    11, 11, 11, 11, 15, 15, 15, 13, 13, 11,
    16, 18, 19, 7, 11, 17, 17, 21, 17, 21,
    21, 21, 17, 5, 44, 23, 10, 17, 56, 29,
    46, 17, 60, 41, 17, 64, 17, 17, 21, 19]

# OC BLM merchantable DBH min (sitset.f BLM; also OC's sp_dbh_min) = the crown-biomass merch gate.
const OC_FFE_DBHMIN = 7.0f0
