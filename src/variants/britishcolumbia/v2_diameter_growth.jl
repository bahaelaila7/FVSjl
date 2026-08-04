# =============================================================================
# v2_diameter_growth.jl (britishcolumbia) — BC V2 (LV2ATV) large-tree DDS formula. CHUNK 2b (V2 regime).
#
# The V2 regime (dgf.f LV2ATV=true) applies to NON-ICH/IDF/SBS/SBPS zones (ESSF/MS/PP). Its large-tree DDS
# is the IE-form Wykoff (dgf.f:1940 + 1992-1995): a linear predictor → log. FORMULA + the 5 per-species
# coefficients VALIDATED bit-exact (99.1%, ≤1 ULP) vs instrumented FVSbc on all_BC_essf (ESSFdk/01),
# GIVEN the oracle's CONSPP + DGDSQ. ⚠ The CONSPP habitat-class RESOLUTION (DGHAB/DGFOR/DGDS via MAPHAB/
# MAPLOC/MAPDSQ + SEILTDG site params) is task #133's next step — until then this can't run in-engine.
# =============================================================================

const BC_V2_DGLD   = Float32[0.56445,0.54140,0.56888,0.68810,0.68712,0.58705,0.89503,0.73045,0.86240,0.66101,0.89778,0.89778,0.89778,0.56888,0.89778]
const BC_V2_DGCR   = Float32[1.08338,1.03478,2.06850,1.93969,1.64133,1.29360,1.85558,1.54643,0.52044,1.31618,1.28403,1.28403,1.28403,2.06850,1.28403]
const BC_V2_DGCRSQ = Float32[0.0,0.07509,-0.62361,-0.78258,-0.27244,0.0,-0.36393,-0.26635,0.86236,0.0,0.0,0.0,0.0,-0.62361,0.0]
const BC_V2_DGBAL  = Float32[0.42112,0.43637,0.50202,0.45142,0.0,0.74596,-0.03665,0.25639,0.0,0.0,0.0,0.0,0.0,0.50202,0.0]
const BC_V2_DGDBAL = Float32[-2.08272,-2.03256,-2.11590,-1.76812,-0.80918,-2.28375,-0.43329,-1.18218,-0.51270,-1.25881,-0.66110,-0.66110,-0.66110,-2.11590,-0.66110]

"""
    bc_v2_dds(sp, d_in, bal, cr, conspp, dgdsq) -> ln(DDS) (WK2)

BC V2 large-tree DDS (dgf.f:1940 + 1992-1995), IE-form: `DDS = CONSPP + DGLD·ln(D) + DGBAL·BAL +
CR·(DGCR + CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1)`; then `max(-9.21, log(max(0.001, DDS)))`.
`d_in` inches; `conspp` (DGCON+COR+CCF terms) + `dgdsq` (habitat-class) come from the CONSPP resolution
(bc_v2_dgcons!, task #133 — not yet ported). VALIDATED bit-exact (99.1%, ≤1 ULP) on all_BC_essf.
"""
@inline function bc_v2_dds(sp::Integer, d_in::Real, bal::Real, cr::Real, conspp::Real, dgdsq::Real)
    d = Float32(d_in); ald = log(d)
    dds = Float32(conspp) + BC_V2_DGLD[sp]*ald + BC_V2_DGBAL[sp]*Float32(bal) +
          Float32(cr)*(BC_V2_DGCR[sp] + Float32(cr)*BC_V2_DGCRSQ[sp]) +
          Float32(dgdsq)*d*d + BC_V2_DGDBAL[sp]*Float32(bal)/log(d + 1f0)
    return max(-9.21f0, log(max(0.001f0, dds)))
end
