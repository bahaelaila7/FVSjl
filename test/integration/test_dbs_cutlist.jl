# test_dbs_cutlist.jl — C6 DBS FVS_CutList table (dbscuts.f) via the CUTLIST keyword.
#
# CUTLIST sends the per-cycle REMOVED records to a SQLite table with the FVS_TreeList per-tree
# columns, but TPA = removed trees/acre. The records are captured non-invasively by `_log_cut!`
# (a gated observer in the thinning path — zero effect when off). As in FVS the table needs BOTH a
# scheduled CUTLIST activity (initre.f opt 92, act 199 → PRTRLS(2) from cuts.f:1738, only on a cycle
# that actually removed trees) AND the DATABASE CUTLIDB flag (dbsin.f opt 17); the earlier form of this
# test (bare CUTLIST, no CUTLIDB) produced no FVS_CutList in live FVS. We validate that the CutList
# RECONSTRUCTS the `.sum` removed columns (RTpa/RTCuFt/RMCuFt) — themselves bit-exact vs Fortran — and
# that the table is absent when either the activity or the DBS flag is missing.

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
$(rpad("CUTLIST", 10) * lpad("0", 10))
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
DATABASE
DSNOUT
$db
SUMMARY
CUTLIDB
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
            @test "TPA" in cols && "SpeciesFVS" in cols && "SpeciesFIA" in cols && "BAPctile" in cols
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
            # TruncHt (dbscuts.f (ITRUNC+5)/100, feet) — was the raw hundredths ITRUNC (100× too large); assert feet-scale.
            th = [Int(r.TruncHt) for r in DBInterface.execute(d, "SELECT TruncHt FROM FVS_CutList")]
            @test all(<(1000), th)                          # feet, not hundredths (the 100× bug)
            # CrWidth (dbscuts.f CW=CRWDTH(I), the shared _forest_crwdth dispatch) — was t.crown_width[i]=0 for most
            # variants; SN (eastern) now gets the open-grown crown_width. Assert real widths present (not all 0).
            cw = [Float64(r.CrWidth) for r in DBInterface.execute(d, "SELECT CrWidth FROM FVS_CutList") if r.CrWidth !== missing]
            @test !isempty(cw) && maximum(cw) > 1.0
        finally
            SQLite.close(d)
        end
        rm(dir; recursive = true, force = true)

        # PRTRLS gating (prtrls.f:84-178 + dbsin.f opt 17): FVS writes FVS_CutList only for a SCHEDULED CUTLIST
        # activity AND the DBS CUTLIDB flag — live FVSsn_g16 creates no table if either is missing.
        for (nm, cutkw, dbkw) in (("no-CUTLIDB", rpad("CUTLIST", 10) * lpad("0", 10), ""),
                                  ("no-CUTLIST", "", "CUTLIDB\n"),
                                  ("CUTLIST-cycle1-no-cut", "CUTLIST", "CUTLIDB\n"))   # blank date ⇒ cycle 1; thin is in 1995
            gdir = mktempdir(); gdb = joinpath(gdir, "out.db")
            cp(tre, joinpath(gdir, "cut.tre"); force = true)
            gkey = joinpath(gdir, "cut.key")
            write(gkey, "STDIDENT\nCUTDB\nSTDINFO        80106   231Dd        60.0     315.0      30.0       7.0\n" *
                  "INVYEAR       1990.0\nNUMCYCLE         3.0\nSITECODE          63      60.\n" *
                  "DESIGN                                        11.0       1.0\n$thin\n" *
                  (isempty(cutkw) ? "" : cutkw * "\n") *
                  "TREEFMT\n(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,\n" *
                  "T52,I2,T66,5I1,T54,7I1,T75,F3.0)\nDATABASE\nDSNOUT\n$gdb\nSUMMARY\n$(dbkw)END\nTREEDATA\nPROCESS\nSTOP\n")
            FVSjl.run_keyfile(gkey; faithful = true)
            d = SQLite.DB(gdb)
            try
                @test !("FVS_CutList" in [r.name for r in SQLite.tables(d)])
            finally
                SQLite.close(d)
            end
            rm(gdir; recursive = true, force = true)
        end
    end
end
