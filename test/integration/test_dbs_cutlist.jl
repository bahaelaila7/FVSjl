# test_dbs_cutlist.jl — C6 DBS FVS_CutList table (dbscuts.f) via the CUTLIST keyword.
#
# CUTLIST sends the per-cycle REMOVED records to a SQLite table with the FVS_TreeList per-tree
# columns, but TPA = removed trees/acre. The records are captured non-invasively by `_log_cut!`
# (a gated observer in the thinning path — zero effect when off). This SN Fortran build writes the
# cut list only to a TEXT dataset (not FVS_CutList), so there is no Fortran table to diff; instead
# we validate that the CutList RECONSTRUCTS the `.sum` removed columns (RTpa/RTCuFt/RMCuFt) — which
# are themselves bit-exact vs Fortran — i.e. Σ(TPA)=RTpa and Σ(TPA·vol)=removed volume.

using Test, FVSjl, SQLite, DBInterface

@testset "C6 DBS — FVS_CutList table (CUTLIST)" begin
    tre = joinpath(@__DIR__, "..", "harness", "scenarios", "dbs_compute.tre")
    if !isfile(tre)
        @test_skip "dbs scenario not available"
    else
        dir = mktempdir()
        db = joinpath(dir, "out.db")
        cp(tre, joinpath(dir, "cut.tre"); force = true)
        key = joinpath(dir, "cut.key")
        thin = rpad("THINBTA", 10) * lpad("1995", 10) * lpad("80", 10)   # residual BA 80 @ 1995
        open(key, "w") do io
            print(io, """
STDIDENT
CUTDB
STDINFO        80106   231Dd        60.0     315.0      30.0       7.0
INVYEAR       1990.0
NUMCYCLE         3.0
SITECODE          63      60.
DESIGN                                        11.0       1.0
$thin
CUTLIST
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
        sumtxt = FVSjl.run_keyfile(key; faithful = true)
        # parse the .sum removed columns at 1995: RTpa=f[13], RTCuFt=f[14], RMCuFt=f[15]
        r1995 = nothing
        for ln in split(sumtxt, "\n")
            f = split(ln); length(f) >= 15 && f[1] == "1995" && (r1995 = f)
        end
        @test r1995 !== nothing
        rtpa = parse(Int, r1995[13]); rtcuft = parse(Int, r1995[14]); rmcuft = parse(Int, r1995[15])
        @test rtpa > 0                                      # a partial thin actually happened

        @test isfile(db)
        d = SQLite.DB(db)
        try
            cols = [r.name for r in DBInterface.execute(d, "PRAGMA table_info(FVS_CutList)")]
            @test "TPA" in cols && "Species" in cols && "BAPctile" in cols
            recs = [NamedTuple(r) for r in DBInterface.execute(d,
                "SELECT TPA,TCuFt,MCuFt FROM FVS_CutList WHERE Year=1995")]
            @test !isempty(recs)
            stpa  = sum(Float64(r.TPA) for r in recs)
            stcuft = sum(Float64(r.TPA) * Float64(coalesce(r.TCuFt, 0.0)) for r in recs)
            smcuft = sum(Float64(r.TPA) * Float64(coalesce(r.MCuFt, 0.0)) for r in recs)
            # the cut records reconstruct the .sum removed aggregates: DBS full-precision Σ vs the RENDERED-INTEGER
            # .sum removed cols (parse(Int,·)) → irreducible width = the PRINT HALF-WIDTH 0.5 (category-2). Was ≤1 (2× pad).
            @test round(Int, stpa)  == round(Int, rtpa)     # Σ removed TPA renders to the .sum integer (was ≤0.5)
            @test round(Int, stcuft) == round(Int, rtcuft)  # Σ removed total cubic renders to the .sum integer
            @test round(Int, smcuft) == round(Int, rmcuft)  # Σ removed merch cubic renders to the .sum integer
        finally
            SQLite.close(d)
        end
        rm(dir; recursive = true, force = true)
    end
end
