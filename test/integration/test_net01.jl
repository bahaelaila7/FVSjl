# test_net01.jl — NE variant net01 validation vs the live FVSne oracle (tests/FVSne/net01.sum.save).
# FVSjulia/Oracle-A has NO NE, so the committed .sum.save is the sole ground truth. Cycle-0 stand state
# (TPA/BA/QMD/TopHt) needs only tree-parse + the shared density — the first bit-exact NE gate.
using Test
using FVSjl

const _NET01_KEY = "/workspace/ForestVegetationSimulator/tests/FVSne/net01.key"

@testset "net01 (NE) cycle-0 stand state — bit-exact vs FVSne oracle" begin
    if !isfile(_NET01_KEY)
        @test_skip "net01.key not available"
    else
        n = first(FVSjl.each_stand(_NET01_KEY; variant = Northeast()))
        FVSjl.notre!(n)
        g = n.plot.gross_space
        di(x) = trunc(Int, x + 0.5)
        # net01.sum.save stand-1 (UNTHINNED) 1990 row: 536 77 160 176 63 5.1 ...
        @test di(stand_tpa(n) / g) == 536               # TPA   — BIT-EXACT
        @test di(stand_ba(n) / g) == 77                 # BA    — BIT-EXACT
        @test round(stand_qmd(n); digits = 1) == 5.1f0  # QMD   — BIT-EXACT
        @test di(stand_top_height(n)) == 63             # TopHt — BIT-EXACT
        # NE-specific parse correctness: first tree is Jack Pine (JP=19), dbh 11.5
        @test n.trees.species[1] == 19
        @test n.trees.dbh[1] ≈ 11.5f0
    end
end
