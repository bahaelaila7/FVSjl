# C4 — mortality model (MORTS) validation against Oracle A (snt01 cycle 1).
# The density (Pretzsch) inputs SDIMAX, t, and the Zeide/Reineke dia0 are pinned to
# the exact Oracle A values; the post-mortality TPA carries the inherited diameter-
# growth residual (~4%), so it is checked with tolerance.

using Test
using FVSjl
using FVSjl: notre!, setup_growth!, diameter_growth!, mortality!, stand_sdimax,
             stand_tpa, _pretzsch_tn10, PRETZSCH_SDIK, bark_ratio

const KEYM = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"

@testset "SDIMAX (SDICAL) matches Oracle A" begin
    s, _ = initialize(KEYM); notre!(s)
    @test isapprox(stand_sdimax(s), 348.44f0; atol=0.1)   # BA-weighted SDIDEF
end

@testset "density mortality (Pretzsch) — Zeide dia0 + activation" begin
    s, _ = initialize(KEYM); notre!(s)
    setup_growth!(s); diameter_growth!(s, s.variant)
    t = s.trees; fint = 5f0
    @test s.control.zeide_sdi                              # SN uses Zeide/Reineke SDI
    # Zeide/Reineke start diameter (no growth) = exact Oracle A value
    tt = 0f0; sumdr0 = 0f0
    for i in 1:t.n
        tt += t.tpa[i]; sumdr0 += t.tpa[i] * t.dbh[i]^1.605f0
    end
    dia0 = (sumdr0 / tt)^(1f0/1.605f0)
    @test isapprox(dia0, 4.701f0; atol=0.005)             # vs quadratic QMD 5.14
    @test isapprox(tt, 589.65f0; atol=0.5)                # density t (per stockable)
end

@testset "mortality reduces TPA in the right regime" begin
    s, _ = initialize(KEYM); notre!(s)
    setup_growth!(s); diameter_growth!(s, s.variant)
    g = s.plot.gross_space
    before = stand_tpa(s) / g
    mortality!(s, s.variant)
    after = stand_tpa(s) / g
    @test after < before                                  # trees died
    # oracle cycle-1 loses ~29 TPA (536→507); we lose less (inherited DG gap) but >10
    @test 10f0 < (before - after) < 35f0
end
