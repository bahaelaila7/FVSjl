# UT (Utah) density-mortality distribution = UTMRT (ut/utmrt.f), regression 2026-08-25.
#
# The western full-population FIA sweep + a mature over-max-SDI sign-tally flagged UT: on MATURE over-max-SDI
# stands the port UNDER-killed cycle-1 TPA (both-sides-traced CN 559750145126144: live 426→245, old jl 426→346;
# CN 4809005010690: live 471→322, old jl 471→391 — SAME direction). Root cause: the port applied a UNIFORM
# per-tree self-thinning rate, but UT's morts.f distributes the density (or background) mortality TOTAL across
# records by UTMRT (ut/utmrt.f, CALL at ut/morts.f:553) — a BA-percentile/shade-tolerance geometric progression
# that CONCENTRATES kills on suppressed (low-percentile = small) trees — and then runs the QMD-convergence loop
# (ut/morts.f:257-603, GO TO 10). Uniform matches the TOTAL first-pass kill but not the distribution; feeding the
# BAMAX residual-BA cap a uniform (BA-heavy) kill removes BA too fast per tree, leaving far too much TPA.
#
# Self-contained unit test (no DB / no oracle) of the two invariants the fix restores:
#   (1) `_varmrt_efftr!(::Utah)` = PEFF·VARADJ·0.01 with PEFF DECREASING in the BA percentile (ut/utmrt.f:97-131);
#       species 17-19,22 (GB/NC/FC/BE) additionally scale by (100−CRI)/100 (the CR-variant form).
#   (2) `_varmrt!` distributes a fixed TOKILL CONCENTRATED on the low-percentile (small) trees, NOT uniformly.
using Test
using FVSjl

@testset "UT (Utah) UTMRT mortality distribution" begin
    # UT VARADJ (ut/utmrt.f): DF(sp3)=0.55, AS(sp6)=1.00, ES(sp8)=0.50, PP(sp10)=0.85.
    peff(pct) = clamp(0.84525f0 - 0.01074f0 * pct + 0.0000002f0 * Float32(pct)^3, 0.01f0, 1.00f0)

    @testset "efftr = PEFF·VARADJ·0.01, decreasing in BA percentile (utmrt.f:97-131)" begin
        t = FVSjl.TreeList(3)
        for (i, pct) in enumerate((10f0, 50f0, 100f0))
            t.species[i] = Int32(3); t.dbh[i] = Float32(2 + i); t.tpa[i] = 100f0
            t.crown_ratio[i] = pct                       # crown_ratio holds FVS's PCT (stand_pct!)
            t.crown_pct[i] = Int32(40)
        end
        t.n = 3
        efftr = zeros(Float32, 3)
        pass1 = FVSjl._varmrt_efftr!(efftr, nothing, FVSjl.Utah(), t, 3)
        # exact UTMRT efficiency for DF (VARADJ 0.55, ·0.01 — NOT the NC ·0.1)
        @test efftr[1] ≈ peff(10f0)  * 0.55f0 * 0.01f0
        @test efftr[2] ≈ peff(50f0)  * 0.55f0 * 0.01f0
        @test efftr[3] ≈ peff(100f0) * 0.55f0 * 0.01f0
        # concentrated on the suppressed (small) tree; dominant hits the 0.01 PEFF floor
        @test efftr[1] > efftr[2] > efftr[3]
        @test efftr[3] ≈ 0.01f0 * 0.55f0 * 0.01f0        # PEFF floored at PCT=100
        @test pass1 ≈ sum(t.tpa[i] * efftr[i] for i in 1:3)
        # species tolerance scales it: an intolerant AS (VARADJ 1.0) is killed harder than a tolerant ES (0.50)
        t.species[1] = Int32(6);  efa = zeros(Float32, 3); FVSjl._varmrt_efftr!(efa, nothing, FVSjl.Utah(), t, 3)
        t.species[1] = Int32(8);  efe = zeros(Float32, 3); FVSjl._varmrt_efftr!(efe, nothing, FVSjl.Utah(), t, 3)
        @test efa[1] > efe[1]
        @test efa[1] ≈ peff(10f0) * 1.00f0 * 0.01f0
        @test efe[1] ≈ peff(10f0) * 0.50f0 * 0.01f0
    end

    @testset "species 17-19,22 use the CR crown-ratio form PEFF·((100−CRI)/100)·VARADJ·0.01" begin
        t = FVSjl.TreeList(2)
        for i in 1:2; t.dbh[i] = 5f0; t.tpa[i] = 100f0; t.crown_ratio[i] = 30f0; t.crown_pct[i] = Int32(40); end
        t.species[1] = Int32(18); t.species[2] = Int32(3)   # sp18 NC (crown form) vs sp3 DF (plain)
        t.n = 2
        efftr = zeros(Float32, 2); FVSjl._varmrt_efftr!(efftr, nothing, FVSjl.Utah(), t, 2)
        @test efftr[1] ≈ peff(30f0) * ((100f0 - 40f0)/100f0) * 0.90f0 * 0.01f0   # sp18 VARADJ=0.90, ×0.6 crown
        @test efftr[2] ≈ peff(30f0) * 0.55f0 * 0.01f0                            # sp3 plain
    end

    @testset "_varmrt! distributes TOKILL concentrated on small trees, total preserved" begin
        t = FVSjl.TreeList(4)
        for (i, pct) in enumerate((8f0, 30f0, 60f0, 100f0))
            t.species[i] = Int32(3); t.dbh[i] = Float32(i); t.tpa[i] = 100f0
            t.crown_ratio[i] = pct; t.crown_pct[i] = Int32(40)
        end
        t.n = 4
        killed = zeros(Float32, 4); efftr = zeros(Float32, 4); temwk2 = zeros(Float32, 4)
        tokill = 120f0
        sumkil = FVSjl._varmrt!(killed, efftr, temwk2, nothing, FVSjl.Utah(), t, 4, tokill)
        @test sumkil ≈ tokill rtol = 1e-3
        @test sum(killed) ≈ tokill rtol = 1e-3
        @test all(killed[i] <= t.tpa[i] + 1f-4 for i in 1:4)
        # strictly monotone-decreasing kill down the size gradient — the whole point of the fix
        @test killed[1] > killed[2] > killed[3] > killed[4]
        @test killed[1] > 3 * killed[4]
    end
end
