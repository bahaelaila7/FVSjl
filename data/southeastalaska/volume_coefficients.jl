# =============================================================================
# volume_coefficients.jl (southeastalaska) — NVEL Region-10 F32 Flewelling profile
# coefficients (SHP_AK / FDBT_AK), transcribed from volume/NVEL/f_alaska.f.
#
# AK volume (ak/sitset.f VOLEQDEF, IREGN=10) assigns three NVEL families:
#   A00F32W###  → PROFILE (Flewelling 2-pt, inside-bark) — the R10 conifers (SF/AF/YC/SS/LP/RC/WH/MH)
#   A00DVEW###  → DVEST woodland (TA/WS/LS/BE/OS/PB/AB/BA/AS/CW/WI/SU/OH)  [NOT yet ported — see volume.jl]
#   A32CURW###  → PROFILE R10 CUR (AD/RA hardwoods)                        [NOT yet ported]
# akt01 is 100% F32 species, so the F32 profile is the cyc0-validation path.
#
# FWINIT (fwinit.f, GEOCODE 'A'): FIA→internal Flewelling JSP:
#   042(YC)→31, 242(RC)→32, 098(SS)→33, 260/263/264(WH/MH/LP/SF/AF/hemlock)→34.
# SHP_AK (f_alaska.f) uses JRSP=JSP-30 with hemlock(4) folding onto spruce(3)'s F column:
#   jsp31→F1(YC), jsp32→F2(RC), jsp33→F3(spruce), jsp34→F3(spruce/hemlock share the shape).
# The SHP kernel is IDENTICAL to the shared _fw2_shp_core (SHP_OT) — only the F table differs.
# FDBT (bark→DBHIB) is NOT shared: it is per-JSP (spruce jsp33 ≠ hemlock jsp34), see _ak_fdbt below.
# =============================================================================

# F(10:47) per Flewelling species column (f[k] == Fortran F(k+9)); 38 values each.
# f_alaska.f DATA (F(I,1)=YC/jsp31, F(I,2)=RC/jsp32, F(I,3)=spruce&hemlock/jsp33,34).
const _AK_VOL_F1 = Float64[   # ALASKA YELLOW CEDAR (jsp 31)
    0.36073443e+01, 0.28808170e+00, 0.36880717e-03, -0.59261770e+00, -0.26176422e+00,
    0.21473459e+00, 0.10424453e-04, 0.00000000e+00, 0.12753659e+05, -0.29418473e+04,
   -0.59464947e-01, -0.52356508e+00, 0.31885182e+01, 0.44300416e+00, -0.11222337e+00,
    0.62259640e+01, -0.24733297e+01, -0.80698148e+01, 0.23730057e+01, 0.69999943e+01,
   -0.97392805e+01, -0.17557578e+01, 0.30946244e+01, -0.88493982e-01, 0.23735651e+01,
    0.12393680e+01, -0.74687251e+00, -0.11540908e+00, 0.41429549e+01, 0.18799039e+01,
   -0.70710595e+00, -0.59270766e-01, -0.17168611e+04, 0.35569956e+03, 0.00000000e+00,
    0.10969542e+01, -0.28489405e-01, -0.87696686e-01]

const _AK_VOL_F2 = Float64[   # WESTERN RED CEDAR (jsp 32)
    0.98477853e+00, 0.58511441e+00, 0.30909028e-03, -0.89204162e+00, -0.23390452e+00,
   -0.49998897e-02, 0.54604625e-01, 0.00000000e+00, 0.12769089e+05, -0.28750172e+04,
    0.23208522e+01, 0.97290535e+00, 0.21836353e+01, 0.96352866e-01, -0.60131935e-01,
    0.46427702e+00, -0.78654704e+00, 0.30346464e+01, -0.50430550e+00, -0.25389446e+01,
   -0.21025226e+02, 0.28035137e+01, 0.87909814e+01, -0.86687109e+00, 0.73773176e+00,
   -0.13553797e+01, -0.24498232e+00, 0.50321072e+00, 0.17521623e+01, 0.35018471e+01,
    0.14592921e+00, -0.12726548e+00, -0.17269922e+04, 0.35994895e+03, 0.00000000e+00,
    0.88681845e+01, 0.58686630e+02, 0.51018047e+00]

const _AK_VOL_F3 = Float64[   # SPRUCE & HEMLOCK (jsp 33 or 34)
    0.49336868e+01, 0.40232019e-01, 0.83224395e-03, 0.90107598e+00, -0.90527505e+00,
    0.89268065e+00, 0.14535135e+00, 0.00000000e+00, 0.78000000e+01, 0.00000000e+00,
    0.00000000e+00, 0.99000000e+00, 0.50129548e+01, -0.32927757e+00, -0.82895182e+00,
    0.46568444e+01, -0.16673502e+01, -0.54307696e+01, 0.12439151e+01, 0.13910769e+02,
   -0.56519913e+01, 0.11187845e+02, 0.69281311e+01, 0.23890379e+02, 0.96523635e+01,
   -0.95228393e+01, -0.20216598e+01, 0.19022197e+01, 0.63832200e+01, 0.46545050e+01,
    0.81554530e+01, 0.14744247e+02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
    0.14356636e+01, -0.50828256e+01, -0.21749222e+02]
