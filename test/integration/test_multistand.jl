# test_multistand.jl — multi-stand driver (each_stand). A keyword file holds several
# stands separated by PROCESS and ended by STOP. FVS re-runs INITRE per stand (ITRN
# resets) but the tree-record format (TREFMT) persists in COMMON across stands, and a
# stand with no TREEDATA still reads the shared tree file (initre.f:334 default INTREE)
# unless NOTREES. This guards three regressions that all surfaced together:
#   (1) stands 2+ fell back to the DEFAULT tree format → misparsed → garbage TPA,
#   (2) the FFE stand (no TREEDATA, only REWIND) loaded zero trees,
#   (3) the trailing STOP after the last PROCESS produced a phantom 6th stand.

using Test, FVSjl

const _MS_KEY = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"

@testset "multi-stand driver (each_stand) — snt01 5 stands" begin
    if !isfile(_MS_KEY)
        @test_skip "snt01.key not available"
    else
        stands = each_stand(_MS_KEY)
        @test length(stands) == 5                  # exactly 5 (no phantom terminator stand)

        # Stands 1-4 all read the SAME snt01.tre (27 trees) — via TREEDATA (1-3) or the
        # default INTREE after REWIND (4) — so every one starts at TPA 536 bit-exact,
        # which only holds if TREFMT persisted across stands.
        for i in 1:4
            s = stands[i]
            notre!(s)
            g = s.plot.gross_space
            @test s.trees.n == 27
            @test isapprox(stand_tpa(s) / g, 536.0; atol = 0.5)
            @test isapprox(stand_ba(s)  / g,  77.0; atol = 1.0)
        end
        # Stand 4 (FFE) has its own inventory year 1993 and NO TREEDATA keyword.
        @test Int(stands[4].control.cycle_year[1]) == 1993

        # Stand 5 is the bare-ground PLANT stand: NOTREES ⇒ no tree-file read, ESTAB active.
        notre!(stands[5])
        @test stands[5].trees.n == 0
        @test stands[5].estab.active
    end
end
