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

    # STRCLASS keyword: activates SSTAGE + overrides thresholds (sawdbh = field 4).
    sk = FVSjl.StandState(FVSjl.Southern())
    rec = FVSjl.KeywordRecord("STRCLASS", "", ["1", "", "", "30", fill("", 8)...],
        Float32[1, 0, 0, 30, zeros(Float32, 8)...], [true, false, false, true, falses(8)...],
        12, FVSjl.KW_OK, 0)
    FVSjl.kw_strclass!(sk, rec)
    @test sk.control.strclass_on
    @test sk.control.strclass_thresh[3] == 30f0           # SAWDBH override applied

    # event-monitor integration: BSCLASS drives an IF condition (fires at UR, not SE).
    if isfile(joinpath(@__DIR__, "..", "harness", "scenarios", "fire_early.key"))
        st, _ = FVSjl.initialize(joinpath(@__DIR__, "..", "harness", "scenarios", "fire_early.key"))
        FVSjl.notre!(st); FVSjl.setup_growth!(st); FVSjl.compute_volumes!(st); FVSjl.compute_density!(st)
        ctx = FVSjl.EventCtx(1, 1990, st)
        @test FVSjl._event_var("BSCLASS", ctx) == 3f0     # UR (matches structure_class)
        @test FVSjl._event_var("BSTRDBH", ctx) > 0f0      # uppermost-stratum DBH
        cond = FVSjl.parse_event_condition("BSCLASS EQ 3")
        @test FVSjl.eval_event(cond, ctx) != 0f0          # fires at 1990 (UR)
        FVSjl.grow_cycle!(st; fint = 5f0); FVSjl.compute_density!(st)
        @test FVSjl.eval_event(cond, FVSjl.EventCtx(2, 1995, st)) == 0f0   # not at 1995 (SE)
    end
end
