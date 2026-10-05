# test_tiered_east.jl — live-fixture regression tests for the fixes that brought CS / LS / NE / OC / OP / ON into the
# tiered suite (test/fixtures/tiered/<v>/, goldens from the live oracles named in each PROVENANCE.toml). Each testset runs
# one fixture case through the committed engine path (tiered_runner.run_case → run_keyfile) and asserts that the cells
# the fix touches now equal the live golden exactly (compare_case: printed .sum fields, DBS cells numeric-exact).
#   julia --project=. -e 'using Test, FVSjl; include("test/integration/test_tiered_east.jl")'
using Test, FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

_case(v, cn, r) = (d = mktempdir(); (txt, db, crashed, err) = run_case(v, cn, r; dir = d);
                   (crashed ? error("jl crashed on $v/$cn/$r: $err") : compare_case(v, cn, r, txt, db)))
_cells(ms, file; col = nothing, year = nothing) =
    [m for m in ms if m.file == file && (col === nothing || m.col == col) && (year === nothing || m.year == string(year))]

@testset "tiered EAST/ORGANON fixes vs live goldens" begin
    @testset "CS 66519757010661 none — R9 InvReference, cycle-0 dead rows, REGENT direct DBH, BALMOD/BAGE5/VDBH" begin
        ms = _case("CS", "66519757010661", "none")
        # voleqdef.f R9_EQN: VEQNNC = '900CLKE'+FIA (was blank)
        @test isempty(_cells(ms, "FVS_InvReference"))
        # crown.f DO 79 dead-crown dub, ptbal.f CASE('CS',…) PTBALT=0, cratet DENSE dead PCT (dead row TreeIndex 3000, 2002)
        @test isempty(_cells(ms, "FVS_TreeList"; year = 2002))
        # cs/regent.f:291-294 DBH(K)=D+0.001*HK before MORTS; bark (0+b·d)/d → b; cs/dgf.f BAGE5 IND1 order + VDBHC*1./D;
        # cs/balmod.f full-precision B2/B3 — the whole projection is now within one TCuFt ULP cell
        @test isempty(_cells(ms, "sum"))
        @test length(ms) <= 1
    end
    @testset "CS 1813567613290487 none — CRATET dead-inclusive DENSE BA for the crown dub (cratet.f:168-170)" begin
        ms = _case("CS", "1813567613290487", "none")
        @test isempty(_cells(ms, "FVS_TreeList"; year = 2024))
    end
    @testset "LS 68756497010661 none — htcalc/balmod libm + association, regent powf, habtyp FVS14" begin
        ms = _case("LS", "68756497010661", "none")
        @test isempty(_cells(ms, "FVS_Error"))                 # ls/habtyp.f:126 FVS14 (initre.f:384-387 default call)
        @test isempty(_cells(ms, "sum"))
        @test isempty(_cells(ms, "FVS_TreeList"; col = "HtG"))  # htcalc.f:394 ((H-BH)/B1)/SI**B2, powf/expf/logf
        @test length(ms) <= 3
    end
    @testset "LS 21073815010661 none — ls/sitset.f:262-267 FVS54 (blank site species/index ⇒ RN/60 default)" begin
        ms = _case("LS", "21073815010661", "none")
        @test isempty(_cells(ms, "FVS_Error"))                 # FVS08, FVS14, FVS54 rows in live order
    end
    @testset "NE 259381087489998 none — crown dub before calibration, dead-inclusive CRATET BA" begin
        ms = _case("NE", "259381087489998", "none")
        @test isempty(_cells(ms, "FVS_TreeList"; year = 2013))
    end
    @testset "OC 645142535126144 none — R6 VEQNNC volumes (R6_EQN westside), ca/sitset.f SDIDEF, national cwcalc" begin
        ms = _case("OC", "645142535126144", "none")
        @test isempty(_cells(ms, "sum"))                       # cycle-0 TCuFt 18506 → live 12329
        @test isempty(_cells(ms, "FVS_InvReference"))          # CFVolEq/BFVolEq 616BEHW…/F06FW2W…, SDIMax
    end
    @testset "OP 1543978576290487 none — R6 VEQNNC volumes" begin
        ms = _case("OP", "1543978576290487", "none")
        @test isempty(_cells(ms, "sum"))
    end
    @testset "OP — Flewelling F_HG for WH (organon/whphg.f SITECV_F) and national cwcalc (were jl crashes)" begin
        for (cn, r) in (("12973618010497", "plant_cal"), ("15072502010497", "none"))
            d = mktempdir(); txt, db, crashed, err = run_case("OP", cn, r; dir = d)
            @test !crashed
        end
        # SITECV_F/F_HG kernel spot values (Float32 recomputation of the Fortran expressions)
        @test FVSjl.op_sitecv_f(30f0, 50f0) > 0f0
        g, phg = FVSjl.op_f_hg(120f0, 80f0, 5f0)
        @test 0f0 < g < 500f0 && phg > 0f0
    end
    @testset "ON 21197892010661 none — merch standards (grinit/sitset) + InvReference LocationCode" begin
        ms = _case("ON", "21197892010661", "none")
        @test isempty(_cells(ms, "FVS_InvReference"))
    end
end
