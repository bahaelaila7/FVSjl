# BM (Blue Mountains) CRWDTH — cwcalc.f BMMAP + Region-6 forest BF via the national cwcalc library (2026-09-19).
# CRWDTH feeds the FVS_TreeList/CutList CrWidth column, THINCC canopy-cover thinning (cuts.f:1228 via CCCLS), the
# COVER CVCW canopy-cover report (cvcw.f:80), SSTAGE/FVS_StrClass (sstage.f:238) and new ESUCKR sprouts.
# VALIDATED vs live FVSbm_g16 on the 120-stand bm_sample120 FIA set: TreeList CrWidth 1261/1261 cycle-0 rows exact
# (10 species) + 87/87 on every identical-state (stand,year); COVER canopy-cover row 120/120 cycle-0 exact; THINCC
# cut selection + per-record PREM bit-exact vs the cuts.f DEBUG stream. Anchors below are live TreeList rows whose
# equations carry no stand-BA / elevation term, so they pin the kernel without stand context.

using Test
using FVSjl

@testset "BM CRWDTH — cwcalc.f BMMAP + R6 BF, bit-exact vs FVSbm_g16" begin
    cw(sp, d, h, cr; kodfor = 614) = FVSjl.bm_crwdth(sp, Float32(d), Float32(h), Float32(cr), 50f0, 45f0, 0f0, kodfor)

    # cwcalc.f:102-106 BMMAP, WP WL DF GF MH WJ LP ES AF PP WB LM PY YC AS CW OS OH
    @test FVSjl._BM_CWMAP == ("11905", "07303", "20205", "01703", "26403", "06405", "10805", "09305", "01905",
                              "12205", "10105", "11301", "23104", "04205", "74605", "74705", "12205", "31206")

    # live FVS_TreeList anchors (BA/elevation-free equations ⇒ exact without stand context)
    @test cw(4, 1.0, 9.0, 80) == 4.871013f0                       # GF 01703 (R1 log form)
    @test cw(4, 1.4, 11.0, 75) == 5.710271f0
    @test cw(18, 5.4, 14.0, 25) == 15.952932f0                    # OH 31206 (power form, D only)
    @test cw(18, 6.9, 17.0, 30) == 17.79632f0
    @test cw(18, 4.0, 16.0, 30) == 13.953997f0
    @test cw(15, 2.1848745, 15.396993, 81) == 7.591541f0          # AS 74605 (R6 model 2; FIA 746 has no BM BF)
    @test cw(15, 0.1, 1.01, 90) == 0.5f0                          # cwcalc.f final clamp CW ≥ 0.5

    # R6 bias factor table (cwcalc.f:554/613/699/751) — per forest × FIASP; absent pairs keep 1.0
    @test FVSjl.bm_cw_bf(614, "202") == 1.055f0
    @test FVSjl.bm_cw_bf(604, "202") == 1.058f0
    @test FVSjl.bm_cw_bf(616, "202") == 1f0                       # Wallowa-Whitman has no DF entry
    @test FVSjl.bm_cw_bf(616, "264") == 1.077f0
    @test FVSjl.bm_cw_bf(607, "073") == 0.879f0
    @test FVSjl.bm_cw_bf(500, "202") == 1f0                       # cwcalc.f:477 — KODFOR < 601 ⇒ no BF

    # BF scales only the '…05' R6 model-2 forms: DF 20205 large tree scales by BF; GF 01703 (R1 log) is BF-free
    @test cw(3, 12.0, 60.0, 40; kodfor = 614) ≈ 1.055f0 * cw(3, 12.0, 60.0, 40; kodfor = 616) rtol = 1f-6
    @test cw(4, 12.0, 60.0, 40; kodfor = 614) == cw(4, 12.0, 60.0, 40; kodfor = 616)

    # cwcalc.f CASE('20205') small-tree branch is `6.0227*1.0*…` — the only '…05' code whose D<1 branch drops BF.
    # (Live: BM FIA DF seedlings D=0.1 carry the no-BF value, e.g. 0.5238 where BF would give 0.5526.)
    @test cw(3, 0.5, 3.0, 95; kodfor = 614) == cw(3, 0.5, 3.0, 95; kodfor = 616)
    @test cw(1, 0.5, 3.0, 95; kodfor = 614) > cw(1, 0.5, 3.0, 95; kodfor = 500)   # WP 11905 keeps BF when D<1
end
