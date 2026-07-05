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
    # SK (6 obs), AB (5 obs) calibrate. RE-GROUNDED vs LIVE FVSsn (debug-stamped dgdriv.f COR(ISPC) at
    # the WC·CORI calibration point, snt01): live COR SK=0.700993, AB=1.085818 — jl matches BIT-EXACT
    # (Float32). The old targets 0.656/1.016 were Oracle A's (FVSjulia) values, ~7% wrong; the loose
    # atol 0.05/0.08 existed only to absorb that OA error. (jl's snt01.sum is bit-exact vs live all cycles.)
    # Float32-ULP floor: jl matches the 6-decimal live stamp to |Δ|=1.19e-7 (one Float32 ULP at this
    # magnitude — the stamp literal rounds to a neighbouring Float32). atol 2f-7 corners that single op
    # (was 1f-4, ~1000× padded; the old 0.05/0.08 absorbed Oracle-A's ~7% error, since removed).
    @test isapprox(s.calib.dg_cor[65], 0.700993f0; atol=2f-7)
    @test isapprox(s.calib.dg_cor[33], 1.085818f0; atol=2f-7)
end

@testset "height growth (HTGF) — per-tree HTG vs the 3-decimal reference" begin
    s, _ = initialize(KEY3)
    notre!(s); s.plot.basal_area = stand_ba(s); s.plot.avg_height = stand_top_height(s)
    height_growth!(s, s.variant)
    t = s.trees
    # per-tree HTG, dominant trees (tree 2 is a tiny seedling — small-tree path, skip). jl's HTGF is the
    # live-validated computation (snt01.sum is bit-exact vs live FVSsn every cycle); it RENDERS to the
    # recorded 3-decimal reference EXACTLY (measured |Δ| ≤ 4.3e-4, all round to the field). Compare the
    # rounded value `==` — the recorded values happen to match Oracle A but are here anchored to jl's own
    # live-validated result at print resolution (was atol 0.002/0.005, ~5-10× the 3-decimal half-width).
    @test round(Float64(t.ht_growth[1]); digits=3) == 1.442   # HI
    @test round(Float64(t.ht_growth[3]); digits=3) == 0.487   # OH
    @test round(Float64(t.ht_growth[4]); digits=3) == 0.141   # SK
    @test round(Float64(t.ht_growth[5]); digits=3) == 1.782   # SK
    @test round(Float64(t.ht_growth[6]); digits=3) == 4.840   # AB
    @test round(Float64(t.ht_growth[8]); digits=3) == 2.311   # HI
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
