# BM (BlueMountains) FVS_TreeList CrWidth via bm_cwcalc + folded forest-614 BF (regression, 2026-08-22).
# BM's crown width reuses the CR Crookston library (bm_cwcalc → cr_cwcalc via _BM_TO_CR_CWSP), and bmt01's forest 614
# (Umatilla) applies a per-FIASP BF on the R6-Model-2 eqns: WP(119)=1.128, DF(202)=1.055, LP(108)=1.244, ES(093)=1.137,
# AF(019)=1.110, PP(122)=1.035 (WL/GF use log-form eqns which have NO BF; the rest = 1.0). The BF is folded into
# cr_cwcalc's leading coef via its `bf` kwarg. IMPORTANT: the BF is a TreeList/forest-grown (IWHO=0) thing — the FFE
# PERCOV path (fmcba) and StrClass are BF-FREE, so bm_cwcalc defaults to forest_bf=false and only _forest_crwdth opts
# in (bm_cwcalc(...; forest_bf=true)). VALIDATED vs FVSbm_clean bmt01 FVS_TreeList: all 29 inventory rows match jl's
# bm_cwcalc(forest_bf=true) at full precision (0 mismatches). BM is the 10th variant wired past the eastern 0.5 default.

using Test
using FVSjl

@testset "BM FVS_TreeList CrWidth = bm_cwcalc + forest-614 BF (forest_bf), bit-exact vs FVSbm_clean" begin
    ba = 85.13126f0; el = 45.0f0; hi = -323.6175f0   # bmt01 stand context (forest 614, el 45; inline ⇒ Hopkins -323.6175)
    cwT(sp, d, h, cr) = clamp(FVSjl.bm_cwcalc(sp, Float32(d), Float32(h), Float32(cr), ba, el, hi; forest_bf = true),  0.5f0, 99.9f0)
    cwF(sp, d, h, cr) = clamp(FVSjl.bm_cwcalc(sp, Float32(d), Float32(h), Float32(cr), ba, el, hi; forest_bf = false), 0.5f0, 99.9f0)
    # BM species indices: DF=3, ES=8. Oracle FVS_TreeList CrWidth (FVSbm_clean bmt01, inventory 1990, with 614 BF):
    @test isapprox(cwT(3, 1.2, 11.0, 55), 5.058; atol = 0.005)   # DF 1.2"/11'/cr55 (BF 1.055)
    @test isapprox(cwT(3, 4.0, 20.0, 25), 8.273; atol = 0.005)   # DF 4.0"/20'/cr25 (BF 1.055)
    @test isapprox(cwT(8, 3.2, 17.0, 45), 9.125; atol = 0.005)   # ES 3.2"/17'/cr45 (BF 1.137)
    # the FFE path (forest_bf=false) is BF-free ⇒ ~5% smaller for a BF>1 species (must NOT equal the TreeList value):
    @test cwF(3, 1.2, 11.0, 55) < cwT(3, 1.2, 11.0, 55)
    @test isapprox(cwF(3, 1.2, 11.0, 55) * 1.055f0, cwT(3, 1.2, 11.0, 55); atol = 0.005)   # BF is exactly ×1.055
end
