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
        FVSjl.compute_forest_type!(s)      # FORTYP/STKVAL (shared, FIA-keyed; CS reuses NE stocking CSVs)
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

        # Trailing classification fields (FORTYP/STKVAL): 503 = W.OAK-R.OAK-HICKORY, size/stock
        # class 22. The .sum growth columns (period/accretion/mortality) need the cycle 0→1 DG
        # projection (cs/dgf.f, chunk 3) and are intentionally not asserted here.
        @test r.fortype  == 503
        @test r.sizecls  == 2
        @test r.stockcls == 2
    end
end

@testset "CS cst01 cycle-1 growth (vs live FVScs)" begin
    key = "/workspace/ForestVegetationSimulator/tests/FVScs/cst01.key"
    if !isfile(key)
        @info "cst01.key not present; skipping CS cycle-1 test"
    else
        # Full growth spine end-to-end: cs_dgf! DG + calibration, cs htgf height growth, cs regent
        # small-tree growth, TWIGS crown, varmrt + background mortality, R9 volume. The DG calibration
        # uses the CS GST DBH≥5 floor (cs/dgdriv.f:380) — with the SN/NE ≥3 floor, WO would spuriously
        # calibrate (debug-stamped: live FN<5 ⇒ COR=0 for all species) and over-grow, inflating BA/SDI.
        s = first(FVSjl.each_stand(key; variant = CentralStates()))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s); FVSjl.compute_forest_type!(s)
        g = s.plot.gross_space; di(x) = trunc(Int, x + 0.5)
        FVSjl.grow_cycle!(s; fint = 10f0)
        # Live FVScs cst01.sum 2000 row — all six stand columns BIT-EXACT:
        @test di(FVSjl.stand_tpa(s) / g)        == 518
        @test di(FVSjl.stand_ba(s) / g)         == 99
        @test di(FVSjl.stand_sdi(s) / g)        == 196
        @test di(FVSjl.stand_ccf(s) / g)        == 202
        @test di(FVSjl.stand_top_height(s))     == 68
        @test round(FVSjl.stand_qmd(s); digits = 1) == 5.9f0

        # cycle-1 grown-stand volume (live 2000: Tcuft 2110 / Mcuft 1887 / Scuft 886 / Bdft 5084).
        # Tcuft + Scuft are bit-exact; Mcuft/Bdft land within 1 unit (~0.05% grown-stand rounding floor).
        FVSjl.compute_volumes!(s)
        v = FVSjl.summary_row(s; period = 0)
        @test v.cuft  == 2110
        @test v.scuft == 886
        @test abs(v.mcuft - 1887) <= 1
        @test abs(v.bdft  - 5084) <= 1
    end
end
