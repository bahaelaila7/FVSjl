# test_dbs_invref.jl — C6 DBS FVS_InvReference table (dbsinvref.f) via the DATABASE block.
#
# FVS_InvReference is a once-per-stand dump of the variant's species master list: FVS/PLANTS/FIA
# codes, the SDI method + per-species SDImax and site index, and the cubic/board volume-equation
# ids + merch specs. All data the engine already holds (after compute_volumes!). Validated
# **bit-exact across all 90 species × 19 columns** vs live Fortran's FVSOut.db; this test checks the
# row count, schema, and a few reference species + the merch-spec constants.

using Test, FVSjl, SQLite, DBInterface

@testset "C6 DBS — FVS_InvReference table" begin
    tre = joinpath(@__DIR__, "..", "harness", "scenarios", "dbs_compute.tre")
    if !isfile(tre)
        @test_skip "dbs scenario not available"
    else
        dir = mktempdir()
        db = joinpath(dir, "out.db")
        cp(tre, joinpath(dir, "ref.tre"); force = true)
        key = joinpath(dir, "ref.key")
        open(key, "w") do io
            print(io, """
STDIDENT
REFDB
STDINFO        80106   231Dd        60.0     315.0      30.0       7.0
INVYEAR       1990.0
NUMCYCLE         1.0
SITECODE          63      60.
DESIGN                                        11.0       1.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
DATABASE
DSNOUT
$db
SUMMARY
END
TREEDATA
PROCESS
STOP
""")
        end
        FVSjl.run_keyfile(key; faithful = true)
        @test isfile(db)
        d = SQLite.DB(db)
        try
            rows = [NamedTuple(r) for r in DBInterface.execute(d,
                "SELECT * FROM FVS_InvReference ORDER BY SpeciesNum")]
            @test length(rows) == 90                              # one row per SN master-list species
            r1 = rows[1]
            @test strip(r1.SpeciesFVS) == "FR" && strip(r1.SpeciesPlants) == "ABIES"
            @test strip(r1.SpeciesFIA) == "010"
            @test strip(r1.SDIType) == "ZEIDE"                    # SN default LZEIDE
            @test r1.SDIMax == 655 && r1.SiteIndex == 58          # bit-exact vs Fortran
            @test strip(r1.CFVolEq) == "841CLKE261" && strip(r1.CFCruiseType) == "FVS"
            # species 1 (FR) merch specs: CF 4/4/0.5, saw 10/7/1, board 10/7/1
            @test r1.CFMinDBH == 4.0 && r1.CFTopDia == 4.0 && r1.CFStump == 0.5
            @test r1.CFSawMinDBH == 10.0 && r1.CFSawTopDia == 7.0 && r1.CFSawStump == 1.0
            @test r1.BFMinDBH == 10.0 && r1.BFTopDia == 7.0 && r1.BFStump == 1.0
            # every species has a populated (non-zero) top diameter + stump (specs filled, not 0)
            @test all(r.CFTopDia > 0 && r.CFStump > 0 && r.BFMinDBH > 0 for r in rows)
            # the woodland/special species carry CFMinDBH = 6 (matched bit-exact vs Fortran)
            @test count(r -> r.CFMinDBH == 6.0, rows) == 9
            @test rows[2].SDIMax == 354 && strip(rows[2].SpeciesFVS) == "JU"   # 2nd species spot-check
        finally
            SQLite.close(d)
        end
        rm(dir; recursive = true, force = true)
    end
end
