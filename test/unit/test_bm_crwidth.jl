# BM (Blue Mountains) CRWDTH — the ONE crown-width implementation `bm_cwcalc` (cwidth.f → cwcalc.f IWHO=0): BMMAP
# (cwcalc.f:102-106) through the national cwcalc library × the Region-6 forest BF for the post-FORKOD KODFOR
# (cwcalc.f:477-860), then the [0.5, 99.9] clamp. FVS keeps one CRWDTH array, read by FMCBA (PERCOV), BMSTAGE,
# FVS_TreeList/CutList, THINCC (cuts.f:1228 via CCCLS), COVER CVCW (cvcw.f:80), SSTAGE/StrClass (sstage.f:238) and new
# ESUCKR sprouts — all of which call bm_cwcalc.
# VALIDATED vs live FVSbm_g16: FVS DEBUG CWIDTH on Malheur(604) 41137075010497 (all inventory trees; PERCOV 53.74 = live);
# bm_sample120 FIA TreeList CrWidth 1261/1261 cycle-0 rows; bmt01 (FVSbm_clean, forest 614) TreeList rows.

using Test
using FVSjl

@testset "BM CRWDTH = bm_cwcalc (BMMAP × R6 forest BF + clamp), bit-exact vs FVSbm" begin
    # --- bmt01 inventory anchors (FVSbm_clean FVS_TreeList, forest 614, BA 85.13126, EL 45, Hopkins -323.6175)
    ba = 85.13126f0; el = 45.0f0; hi = -323.6175f0
    cw(sp, d, h, cr; kf = 614) = FVSjl.bm_cwcalc(sp, Float32(d), Float32(h), Float32(cr), ba, el, hi; kodfor = kf)
    @test isapprox(cw(3, 1.2, 11.0, 55), 5.058; atol = 0.005)   # DF 1.2"/11'/cr55 (614 BF 1.055)
    @test isapprox(cw(3, 4.0, 20.0, 25), 8.273; atol = 0.005)   # DF 4.0"/20'/cr25
    @test isapprox(cw(8, 3.2, 17.0, 45), 9.125; atol = 0.005)   # ES 3.2"/17'/cr45 (614 BF 1.137)

    # --- live FIA TreeList anchors whose equations have no stand-BA / elevation term (exact without stand context)
    cwx(sp, d, h, cr) = FVSjl.bm_cwcalc(sp, Float32(d), Float32(h), Float32(cr), 50f0, 45f0, 0f0; kodfor = 614)
    @test cwx(4, 1.0, 9.0, 80) == 4.871013f0                      # GF 01703 (R1 log form)
    @test cwx(4, 1.4, 11.0, 75) == 5.710271f0
    @test cwx(18, 5.4, 14.0, 25) == 15.952932f0                   # OH 31206 (power form, BM-unique)
    @test cwx(18, 6.9, 17.0, 30) == 17.79632f0
    @test cwx(18, 4.0, 16.0, 30) == 13.953997f0
    @test cwx(15, 2.1848745, 15.396993, 81) == 7.591541f0         # AS 74605 (FIA 746 has no BM BF)
    @test cwx(15, 0.1, 1.01, 90) == 0.5f0                         # cwcalc.f final clamp CW ≥ 0.5
    @test cw(4, 0.1, 1.0, 95) == 0.5f0                            # (live DEBUG CWIDTH 0.500)

    # --- BMMAP (cwcalc.f:102-106): WP WL DF GF MH WJ LP ES AF PP WB LM PY YC AS CW OS OH
    @test FVSjl._BM_CWEQN == ("11905", "07303", "20205", "01703", "26403", "06405", "10805", "09305", "01905",
                              "12205", "10105", "11301", "23104", "04205", "74605", "74705", "12205", "31206")

    # --- per-forest BF lookup (cwcalc.f SELECT CASE (KODFOR)); only the R6 model-2 ('…05') forms carry BF
    @test FVSjl.bm_cw_bf(3, 614) == 1.055f0      # DF Umatilla
    @test FVSjl.bm_cw_bf(3, 604) == 1.058f0      # DF Malheur (the measured 1/0.9452)
    @test FVSjl.bm_cw_bf(3, 616) == 1f0          # Wallowa-Whitman has no DF entry
    @test FVSjl.bm_cw_bf(2, 607) == 1f0          # WL 07303 is log-form ⇒ no BF even though 607 lists 073
    @test FVSjl.bm_cw_bf(10, 614) == 1.035f0     # PP Umatilla
    @test FVSjl.bm_cw_bf(10, 604) == 1f0         # PP has no Malheur entry
    @test FVSjl.bm_cw_bf(4, 614) == 1f0          # GF 01703 is log-form ⇒ no BF even though 614 lists 017
    @test FVSjl.bm_cw_bf(3, 1000) == 1f0         # KODFOR ≥ 1000 skips the BF section
    @test FVSjl.bm_cw_bf(3, 500) == 1f0          # KODFOR < 601 skips it too

    # the BF scales the whole width of a large tree: Malheur/Umatilla DF ratio == 1.058/1.055
    @test isapprox(cw(3, 4.0, 20.0, 25; kf = 604) / cw(3, 4.0, 20.0, 25), 1.058f0 / 1.055f0; rtol = 1f-5)
    @test cw(4, 12.0, 60.0, 40; kf = 614) == cw(4, 12.0, 60.0, 40; kf = 616)          # GF: BF-free everywhere

    # cwcalc.f CASE('20205') small-tree branch is `6.0227*1.0*…` — the only '…05' code whose D<1 branch drops BF
    # (live: BM FIA DF seedlings D=0.1 carry the no-BF value, e.g. 0.5238 where BF would give 0.5526).
    @test cw(3, 0.5, 3.0, 95; kf = 614) == cw(3, 0.5, 3.0, 95; kf = 616)
    @test cw(1, 0.5, 3.0, 95; kf = 614) > cw(1, 0.5, 3.0, 95; kf = 500)                # WP 11905 keeps BF at D<1
end
