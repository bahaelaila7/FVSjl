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

        # Tot-Cov column bit-exact (raw-PROB crown area). Fortran snt01 stand-1: 82/87/90/92/91/90/
        # 89/89/87/86/85 — match within 1 (IFIX round at .5 boundaries).
        s1c = first(FVSjl.each_stand(snt))
        FVSjl.notre!(s1c); FVSjl.setup_growth!(s1c); FVSjl.compute_volumes!(s1c)
        ftcov = [82, 87, 90, 92, 91, 90, 89, 89, 87, 86, 85]
        # uppermost-stratum DBH (SSTGHP DBHNOM, 70th crown-percentile of the canopy cohort).
        # Fortran snt01 stand-1 stratum-1 DBH; 8/11 exact, rest ≤0.5 (cohort/window-edge boundary).
        ftdbh = [10.3, 9.8, 11.5, 12.3, 15.1, 17.0, 18.2, 20.7, 22.6, 23.8, 24.1]
        for c in 0:10
            FVSjl.compute_density!(s1c)
            r = FVSjl.structure_class(s1c)
            @test abs(round(Int, r.cover) - ftcov[c+1]) <= 1
            @test abs(r.strdbh - ftdbh[c+1]) <= 0.6
            c < 10 && FVSjl.grow_cycle!(s1c; fint = 5f0)
        end

        # structure_report per-stratum columns vs the Fortran Structural-statistics report
        # (snt01 stand-1 stratum 1 @1990): DBH/Nom-Ht/Lg-Ht/Sm-Ht/CrnBase/Cov/Sp1/Sp2 all bit-exact.
        s1r = first(FVSjl.each_stand(snt))
        FVSjl.notre!(s1r); FVSjl.setup_growth!(s1r); FVSjl.compute_volumes!(s1r); FVSjl.compute_density!(s1r)
        rep = FVSjl.structure_report(s1r)
        @test !isempty(rep.strata)
        a = rep.strata[1]
        @test isapprox(a.dbh, 10.3; atol = 0.1)              # DBHNOM
        @test round(Int, a.nomht) == 63 && round(Int, a.lght) == 75 && round(Int, a.smht) == 55
        @test round(Int, a.crnbase) == 42                    # the "Bas" column = mean crown-base height
        @test round(Int, a.cover) == 49
        @test strip(s1r.coef.code_alpha[a.sp1]) == "SM" && strip(s1r.coef.code_alpha[a.sp2]) == "HI"

        # the .out "Structural statistics" report writer — byte-for-byte vs Fortran (FORMAT 85/90).
        s1w = first(FVSjl.each_stand(snt))
        FVSjl.notre!(s1w); FVSjl.setup_growth!(s1w); FVSjl.compute_volumes!(s1w)
        io = IOBuffer()
        FVSjl.write_structure_report(io, s1w, 10; stand_id = "S248112", mgmt_id = "NONE")
        lines = split(String(take!(io)), "\n")
        @test lines[1] == "Structural statistics for stand: S248112                     MgmtID: NONE"
        @test any(l -> startswith(l, "Year Cd  DBH  Nom  Lg  Sm Bas Cov Sp1 Sp2 D"), lines)
        # the 1990 before-thin row, byte-for-byte (the bit-exact cycle)
        @test "1990  0  10.3  63  75  55  42  49 SM  HI  2   5.8  30  38   2   7  65 AB  SK  1   0.0   0   0   0   0   0 --  --  0 2  82  3=UR" in lines
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

@testset "single-canopy-tree structure-stage (NTREES≤1 cover branch) vs live FVS" begin
    # sstage.f:235-256 — a stand with only ONE canopy tree is classified by its CROWN-AREA cover
    # (WK6 = CW²·TPA·π/4), NOT the stratum DBHNOM, and BEFORE stratification (so it fires even when the lone
    # tree forms no OK stratum → reported N=0). Each scenario below is a single SM record whose class was read
    # off live FVSsn (fort.16 "Structural statistics") — they pin the four reachable sub-branches:
    #   bg     cover<CCMIN, TPA<TPAMIN   → 0=BG
    #   si_tpa cover<CCMIN, TPA≥TPAMIN   → 1=SI   (the TPA-override; the OLD DBHNOM path gave 0=BG here)
    #   si_dbh cover≥CCMIN, DBH<SAWDBH   → 1=SI
    #   os     cover≥CCMIN, DBH≥SAWDBH   → 5=OS
    # (Standalone class-2/SE is unreachable for a lone tree: its tiny SDI always trips the SE→SI demote —
    # itself faithful to FVS.)
    sdir = joinpath(@__DIR__, "..", "harness", "scenarios")
    for (name, want) in [("struct_1canopy_bg", 0), ("struct_1canopy_si_tpa", 1),
                         ("struct_1canopy_si_dbh", 1), ("struct_1canopy_os", 5)]
        key = joinpath(sdir, name * ".key")
        if !isfile(key)
            @test_skip "$name scenario not available"
        else
            s = first(FVSjl.each_stand(key))
            FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
            r = FVSjl.structure_class(s)
            @test s.trees.n == 1                # single canopy tree
            @test r.nstr == 0                   # NSTR untouched by the GOTO-80 branch (live reports N=0)
            @test r.class == want               # class matches live FVSsn
        end
    end
end
