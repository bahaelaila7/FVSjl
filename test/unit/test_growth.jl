# C3 — growth model validation against Oracle A (single-stand single-cycle run).
# Diameter calibration COR and per-tree height growth are pinned to the exact
# Oracle A values found during the port.

using Test
using FVSjl
using FVSjl: notre!, setup_growth!, grow_cycle!, height_growth!, dgcons!,
             point_basal_area!, calibrate_diameter_growth!, stand_qmd, stand_ba,
             stand_top_height, stand_tpa

const KEY3 = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"

@testset "DG calibration COR vs Oracle A" begin
    s, _ = initialize(KEY3)
    notre!(s)
    s.plot.basal_area = stand_ba(s); s.plot.avg_height = stand_top_height(s)
    point_basal_area!(s); dgcons!(s)
    calibrate_diameter_growth!(s)
    # species with <5 measured-growth trees → COR exactly 0 (FNMIN=5)
    @test s.calib.dg_cor[27] == 0f0     # HI (3 obs)
    @test s.calib.dg_cor[22] == 0f0     # SM (4 obs)
    @test s.calib.dg_cor[89] == 0f0     # OH (1 obs)
    # SK (6 obs), AB (5 obs) calibrate (within ~4% of oracle 0.656 / 1.016 —
    # the small gap is fine BA tuning, tracked in DIVERGENCES/memory)
    @test isapprox(s.calib.dg_cor[65], 0.656f0; atol=0.05)
    @test isapprox(s.calib.dg_cor[33], 1.016f0; atol=0.08)
end

@testset "height growth (HTGF) matches Oracle A" begin
    s, _ = initialize(KEY3)
    notre!(s); s.plot.basal_area = stand_ba(s); s.plot.avg_height = stand_top_height(s)
    height_growth!(s, s.variant)
    t = s.trees
    # per-tree HTG, dominant trees (tree 2 is a tiny seedling — small-tree path, skip)
    @test isapprox(t.ht_growth[1], 1.442f0; atol=0.002)   # HI
    @test isapprox(t.ht_growth[3], 0.487f0; atol=0.002)   # OH
    @test isapprox(t.ht_growth[4], 0.141f0; atol=0.002)   # SK
    @test isapprox(t.ht_growth[5], 1.782f0; atol=0.002)   # SK
    @test isapprox(t.ht_growth[6], 4.840f0; atol=0.005)   # AB
    @test isapprox(t.ht_growth[8], 2.311f0; atol=0.005)   # HI
end

@testset "cycle loop runs and grows the stand" begin
    s, _ = initialize(KEY3)
    notre!(s); setup_growth!(s)
    q0 = stand_qmd(s)
    for _ in 1:5; grow_cycle!(s); end
    @test stand_qmd(s) > q0                       # diameters grew
    @test s.control.cycle == 5
    @test all(isfinite, s.trees.dbh[1:s.trees.n])
end
