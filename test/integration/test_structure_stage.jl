# test_structure_stage.jl — SSTAGE stand structural-stage class (structure_stage.jl).
#
# `structure_class` classifies the stand each cycle (1=SI, 2=SE, 3=UR, 4/5/6) via canopy
# height-gap stratification + per-stratum cover + the dominant-cohort 70th-percentile DBH
# (SSTGHP). The CLASS is SSTAGE's primary discrete output and is validated bit-exact against the
# Fortran "Structural statistics" report's Struct-Class column (produced by the STRCLASS keyword,
# which works even in the stripped ground-truth binary). Two fire-free + one FFE stand confirm it.
# (The full per-stratum REPORT columns need the exact CRWDTH source — see SSTAGE_chunk_plan.md.)

using Test, FVSjl

# Step the projection and collect the structural class at each cycle's start.
function _ss_classes(stand, ncyc)
    classes = Int[]
    FVSjl.notre!(stand); FVSjl.setup_growth!(stand); FVSjl.compute_volumes!(stand)
    for c in 0:ncyc
        FVSjl.compute_density!(stand)
        push!(classes, FVSjl.structure_class(stand).class)
        c < ncyc && FVSjl.grow_cycle!(stand; fint = 5f0)
    end
    return classes
end

@testset "SSTAGE — structural-stage class vs Fortran" begin
    # COVOLP cover helper: cap at 100, monotone in crown area.
    @test FVSjl._ss_cover([1e9], 1:1, 1.0) == 100.0
    @test 0.0 < FVSjl._ss_cover([5000.0, 5000.0], 1:2, 1.0) < 100.0

    sc = joinpath(@__DIR__, "..", "harness", "scenarios", "fire_early.key")
    if !isfile(sc)
        @test_skip "fire_early scenario not available"
    else
        s, _ = FVSjl.initialize(sc)
        got = _ss_classes(s, 5)
        # Fortran Structural statistics Struct-Class column (fire_early, FFE stand):
        # 1990=3(UR), 1995..2015=2(SE)
        @test got == [3, 2, 2, 2, 2, 2]
    end

    snt = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"
    if !isfile(snt)
        @test_skip "snt01.key not available"
    else
        s1 = first(FVSjl.each_stand(snt))
        got = _ss_classes(s1, 10)
        # Fortran snt01 stand-1 (S248112) Struct-Class: UR@1990, SE@1995..2040
        @test got == vcat([3], fill(2, 10))
    end
end
