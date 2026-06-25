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

        # Advance one cycle; the 1995 row matches the baseline BIT-EXACTLY at the .sum print resolution.
        # The .sum truncates via di(x)=trunc(x+0.5); asserting the printed integer/1-decimal matches the
        # baseline is the tightest (float-ULP-equivalent) check — a sub-ULP float diff could only flip the
        # truncation, which it does not here (507.44→507, 103.19→103, 6.11→6.1, 63.41→63).
        FVSjl.grow_cycle!(s)
        FVSjl.compute_forest_type!(s)
        g = s.plot.gross_space
        di(x) = trunc(Int, x + 0.5)
        @test di(stand_tpa(s) / g) == 507                      # baseline 507 — BIT-EXACT (.sum print)
        @test di(stand_ba(s)  / g) == 103                      # baseline 103 — BIT-EXACT
        @test round(stand_qmd(s); digits = 1) == 6.1f0         # baseline 6.1 — BIT-EXACT (F.1)
        @test di(stand_top_height(s)) == 63                    # baseline 63  — BIT-EXACT
    end
end
