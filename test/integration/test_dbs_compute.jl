# test_dbs_compute.jl — C6 DBS FVS_Compute table (dbscmpu.f) via the DATABASE/COMPUTDB block.
#
# COMPUTDB sends the per-cycle COMPUTE event-monitor variables to a dynamic-schema SQLite table
# (one REAL column per variable). Rows are written for the GROWING cycles only (the event monitor
# runs during growth), matching Fortran. The values are the start-of-cycle event-monitor variables
# — here MYBA=BBA and MYSDI=BSDI, both validated bit-exact (to Float32) vs live Fortran's FVSOut.db
# (this also exercises the BSDI = raw Reineke SDIBC fix end-to-end through the DBS path).

using Test, FVSjl, SQLite, DBInterface

@testset "C6 DBS — FVS_Compute table (COMPUTDB)" begin
    tre = joinpath(@__DIR__, "..", "harness", "scenarios", "dbs_compute.tre")
    if !isfile(tre)
        @test_skip "dbs_compute scenario not available"
    else
        dir = mktempdir()
        db = joinpath(dir, "out.db")
        cp(tre, joinpath(dir, "dbs_compute.tre"); force = true)
        key = joinpath(dir, "dbs_compute.key")
        open(key, "w") do io
            print(io, """
STDIDENT
CMPDB
STDINFO        80106   231Dd        60.0     315.0      30.0       7.0
INVYEAR       1990.0
NUMCYCLE         3.0
SITECODE          63      60.
DESIGN                                        11.0       1.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
DATABASE
DSNOUT
$db
COMPUTDB
END
COMPUTE            0
MYBA = BBA
MYSDI = BSDI
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
            cols = [r.name for r in DBInterface.execute(d, "PRAGMA table_info(FVS_Compute)")]
            @test "MYBA" in cols && "MYSDI" in cols          # dynamic schema: one col per COMPUTE var
            got = [(Int(r.Year), Float64(r.MYBA), Float64(r.MYSDI))
                   for r in DBInterface.execute(d, "SELECT Year,MYBA,MYSDI FROM FVS_Compute ORDER BY Year")]
            # rows only for the GROWING cycles (1990/1995/2000, NOT the final 2005)
            @test length(got) == 3
            @test [g[1] for g in got] == [1990, 1995, 2000]
            # values bit-exact (Float32) vs live Fortran FVS_Compute (BSDI fix → MYSDI)
            want = [(77.39207, 202.93901), (103.19379, 252.89455), (126.26046, 292.84149)]
            for (g, w) in zip(got, want)
                @test isapprox(g[2], w[1]; atol = 0.01)      # MYBA = BBA
                @test isapprox(g[3], w[2]; atol = 0.05)      # MYSDI = BSDI (raw Reineke)
            end
        finally
            SQLite.close(d)
        end
        rm(dir; recursive = true, force = true)
    end
end
