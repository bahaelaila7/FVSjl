# BM (BlueMountains) CRWDTH via bm_cwcalc + the Region-6 forest BF (cwcalc.f) — regression.
# BM's crown width reuses the CR Crookston library (bm_cwcalc → cr_cwcalc via _BM_TO_CR_CWSP) × a per-forest BF keyed
# by the post-FORKOD KODFOR × FIASP=CWEQN(1:3), applied only on the R6 ·BF· equation forms (WL 07303 / GF 01703 /
# MH 26403 are log-form: no BF), then the cwcalc.f [0.5, 99.9] clamp. FVS has ONE CRWDTH array (CWIDTH→CWCALC IWHO=0)
# that FVS_TreeList, FMCBA (PERCOV) and BMSTAGE all read — the earlier "FFE path is BF-free" belief was REFUTED by
# FVSbm_g16 DEBUG CWIDTH on Malheur(604) stand 41137075010497: BF-free DF was ×0.9452 (=1/1.058, the 604 DF BF), a
# 614-baked PP ×1.035 (604 has no PP BF); with the forest BF all 29 inventory trees match and PERCOV 53.74 = live.
# The 614 rows below: FVSbm_clean bmt01 FVS_TreeList CrWidth (inventory 1990).

using Test
using FVSjl

@testset "BM CRWDTH = bm_cwcalc × R6 forest BF (per KODFOR) + cwcalc clamp" begin
    ba = 85.13126f0; el = 45.0f0; hi = -323.6175f0   # bmt01 stand context (forest 614, el 45; inline ⇒ Hopkins -323.6175)
    cw(sp, d, h, cr; kf = 614) = FVSjl.bm_cwcalc(sp, Float32(d), Float32(h), Float32(cr), ba, el, hi; kodfor = kf)
    # BM species indices: DF=3, ES=8. Oracle FVS_TreeList CrWidth (FVSbm_clean bmt01, forest 614):
    @test isapprox(cw(3, 1.2, 11.0, 55), 5.058; atol = 0.005)   # DF 1.2"/11'/cr55 (614 BF 1.055)
    @test isapprox(cw(3, 4.0, 20.0, 25), 8.273; atol = 0.005)   # DF 4.0"/20'/cr25 (614 BF 1.055)
    @test isapprox(cw(8, 3.2, 17.0, 45), 9.125; atol = 0.005)   # ES 3.2"/17'/cr45 (614 BF 1.137)
    # per-forest BF lookup (cwcalc.f SELECT CASE (KODFOR)):
    @test FVSjl.bm_cw_bf(3, 614) == 1.055f0      # DF Umatilla
    @test FVSjl.bm_cw_bf(3, 604) == 1.058f0      # DF Malheur (the measured 1/0.9452)
    @test FVSjl.bm_cw_bf(10, 614) == 1.035f0     # PP Umatilla
    @test FVSjl.bm_cw_bf(10, 604) == 1f0         # PP has no Malheur entry
    @test FVSjl.bm_cw_bf(4, 614) == 1f0          # GF 01703 is log-form ⇒ no BF even though 614 lists 017
    @test FVSjl.bm_cw_bf(3, 1000) == 1f0         # KODFOR ≥ 1000 skips the BF section
    # the BF scales the whole width: Malheur/Umatilla DF ratio == 1.058/1.055 (same tree, same stand context)
    @test isapprox(cw(3, 4.0, 20.0, 25; kf = 604) / cw(3, 4.0, 20.0, 25), 1.058f0 / 1.055f0; rtol = 1f-5)
    # cwcalc.f:2391-2392 final clamp: a 0.1" seedling floors at 0.5 (live DEBUG CWIDTH 0.500)
    @test cw(4, 0.1, 1.0, 95) == 0.5f0
end
