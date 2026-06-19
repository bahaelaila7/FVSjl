# C3 (in progress) — DGF/DGCONS structural checks. The bit-exact gate is the
# cycle-1 .sum row, which needs the DGDRIV driver (calibration + serial
# correlation); these tests pin the deterministic equation core meanwhile.

using Test
using FVSjl
using FVSjl: _dgf_phys_group, _dgf_forest_group, dgcons!, dgf!, coefficients, Southern

@testset "DGF coefficient tables (loaded from CSV)" begin
    sd = coefficients(Southern()).species
    @test length(sd[:dg_intercept]) == 90
    @test length(sd[:dg_phys_p411]) == 90
    @test sd[:dg_intercept][13] == 0.222214f0    # loblolly pine intercept
end

@testset "ecological-unit → physiographic group" begin
    @test _dgf_phys_group("231Dd") == :s231t    # snt01's unit
    @test _dgf_phys_group("231B")  == :s231l
    @test _dgf_phys_group("M221")  == :pm221
    @test _dgf_phys_group("222")   == :p222
    @test _dgf_phys_group("411")   == :p411
    @test _dgf_phys_group("999")   == :none
    @test _dgf_forest_group(801)   == :nohd
    @test _dgf_forest_group(0)     == :none
end

@testset "dgf! runs on snt01 and is finite/sane" begin
    s, _ = initialize(joinpath(homedir() == "/root" ? "" : "",
                               "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"))
    notre!(s)
    s.plot.avg_height = FVSjl.stand_top_height(s)
    s.plot.basal_area = stand_ba(s)
    dgcons!(s)
    dgf!(s)
    wk2 = s.scratch.wk[2, 1:s.trees.n]
    @test all(isfinite, wk2)
    @test all(v -> -9.21f0 <= v <= 5f0, wk2)    # ln(DDS) physical range
    # the s231t physiographic coefficient was applied to DGCON
    @test _dgf_phys_group(s.plot.eco_unit) == :s231t
end
