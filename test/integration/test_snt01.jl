# test_snt01.jl — C5 .sum regression: the snt01 cycle-0 row is bit-exact on every
# stand-state/volume/class field, and cycle-1 tracks the Fortran baseline closely.
# Locks in the validated volume + DG-calibration + tripling + summary pipeline.

using Test
using FVSjl

const _SNT01_KEY = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"

@testset "snt01 .sum cycle-0 (bit-exact start-state) + cycle-1 tracking" begin
    if !isfile(_SNT01_KEY)
        @test_skip "snt01.key not available"
    else
        s, _ = initialize(_SNT01_KEY)
        notre!(s)
        FVSjl.setup_growth!(s)
        FVSjl.compute_forest_type!(s)
        FVSjl.compute_volumes!(s)

        r0 = FVSjl.summary_row(s; period = 5)
        # Fortran snt01.sum baseline, 1990 row (start-of-period columns + classes).
        @test r0.year == 1990
        @test r0.age == 60
        @test r0.tpa == 536
        @test r0.ba == 77
        @test r0.sdi == 160
        @test r0.ccf == 218
        @test r0.topht == 63
        @test isapprox(r0.qmd, 5.1; atol = 0.05)
        @test r0.cuft == 1368        # total cubic
        @test r0.mcuft == 1149       # merch cubic
        @test r0.scuft == 68         # sawtimber cubic
        @test r0.bdft == 285         # board feet
        @test r0.fortype == 520
        @test r0.sizecls == 2
        @test r0.stockcls == 2
        @test isapprox(r0.mai, 19.1; atol = 0.05)

        # Advance one cycle; 1995 row should track the baseline (507/103/6.1) closely.
        FVSjl.grow_cycle!(s)
        FVSjl.compute_forest_type!(s)
        g = s.plot.gross_space
        @test isapprox(stand_tpa(s) / g, 507; atol = 3)        # baseline 507
        @test isapprox(stand_ba(s)  / g, 103; atol = 4)        # baseline 103
        @test isapprox(stand_qmd(s), 6.1; atol = 0.15)         # baseline 6.1
        @test isapprox(stand_top_height(s), 63; atol = 1)      # baseline 63
    end
end
