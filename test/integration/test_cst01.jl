# =============================================================================
# test_cst01.jl — Central States (CS) variant, cst01 cycle-0 stand columns
#
# Validates the CS variant's inventory-cycle (1990) stand statistics against the
# live FVScs .sum (tests/FVScs/cst01.sum). All six geometric columns are bit-exact
# (TPA/BA/SDI/CCF reported per-gross-acre = stockable ÷ GROSPC; QMD/TopHt unscaled).
# This stand exercises the GROSPC<1 path (11 plots, 1 non-stockable ⇒ GROSPC 0.909)
# that the SN/NE test stands never hit — the .sum writer's ÷GROSPC scale-back.
# =============================================================================

using Test
using FVSjl

@testset "CS cst01 cycle-0 stand columns (vs live FVScs)" begin
    key = "/workspace/ForestVegetationSimulator/tests/FVScs/cst01.key"
    if !isfile(key)
        @info "cst01.key not present; skipping CS cycle-0 test"
    else
        s, reason = FVSjl.initialize(key; variant = CentralStates())
        @test reason == :process
        @test s.trees.n == 27
        FVSjl.notre!(s)                    # NOTRE — BAF/fixed-plot expansion × GROSPC
        FVSjl.setup_growth!(s)             # CRATET dub + cs_dgcons! bark copy (CFTOPK broken-top)
        FVSjl.compute_volumes!(s)          # eastern R9 Clark cubic + R9LOGS board feet, CS merch
        r = FVSjl.summary_row(s; period = 0)

        # Live FVScs cst01.sum, 1990 inventory row (per-gross-acre):
        @test r.tpa   == 536
        @test r.ba    == 77
        @test r.sdi   == 160
        @test r.ccf   == 169
        @test r.topht == 63
        @test round(r.qmd, digits = 1) == 5.1

        # Volume columns (per-gross-acre) — all bit-exact. Sawtimber-cubic + board feet ride the
        # CS merch standards (_cs_merch); total/merch cubic depend on the CFTOPK broken-top path
        # (cs_dgcons! bark copy — without it the two broken-top trees read ~4 cuft/acre low).
        @test r.cuft  == 1517
        @test r.mcuft == 1300
        @test r.scuft == 497
        @test r.bdft  == 2903
    end
end
