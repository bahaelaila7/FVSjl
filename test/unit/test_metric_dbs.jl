# test_metric_dbs.jl — BC/ON metric DBS tables (metric/dbsqlite) vs LIVE FVSbc_clean / FVSon_g16.
#
# The BC and ON builds link metric/dbsqlite + metric/vbase: DBSSUMRY/DBSTRLS/DBSCUTS/DBSATRTLS always write *_Metric
# tables (East naming for VARACD CS/LS/NE/SN/ON, dbssumry.f:91-99, dbstrls.f:121-136), and the summary row is
# sumout.f:328-352's per-ha conversion of disply.f's imperial IOSUM integers. Fixtures (test/fixtures/metric_dbs):
# 5-tree BC and 8-tree ON (ont01) stands, TREELIST every cycle + a THINDBH 30% cut at the inventory year with
# CUTLIST/ATRTLIST, goldens = the live DBs (bcmet.live.db / onmet.live.db) and the live .out summary rows.
#
# Asserted exact against live: table names + column names/types, the inventory .sum row and FVS_Summary_Metric row,
# the ON Summary_East_Metric being EMPTY (its INSERT names MCuM,SCuM,NCuM,… that its CREATE lacks ⇒ prepare fails,
# dbssumry.f:185-204), and every inventory-year TreeList/CutList/ATRTList column except the per-tree volumes (the BC/ON
# volume kernels carry a ≤3e-5 relative residual — checked at rtol 5e-5) and the Cut/ATRT DG (open: DBSCUTS binds the
# post-calibration DG(I), jl keeps the input increment). Later cycles carry BC/ON growth residuals and are not asserted.
using Test, FVSjl, SQLite, DBInterface

_mcols(db, t) = [(r.name, lowercase(r.type)) for r in DBInterface.execute(db, "PRAGMA table_info($t)")]
_mrows(db, t, yr) = Dict(r[:TreeIndex] => r for r in
    (Dict(pairs(NamedTuple(x))) for x in DBInterface.execute(db, "SELECT * FROM $t WHERE Year = $yr")))

@testset "BC/ON metric DBS tables vs live" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "metric_dbs")
    for (stem, variant, east, yr0) in (("bcmet", FVSjl.BritishColumbia(), false, 1992),
                                       ("onmet", FVSjl.Ontario(), true, 2004))
        @testset "$stem" begin
            dir = mktempdir()
            for f in ("$stem.key", "$stem.tre"); cp(joinpath(fx, f), joinpath(dir, f)); end
            sumtxt = cd(() -> FVSjl.run_keyfile("$stem.key"; variant = variant, output = :sum), dir)
            live = SQLite.DB(joinpath(fx, "$stem.live.db")); jl = SQLite.DB(joinpath(dir, "$stem.db"))
            try
                sfx = east ? "_East_Metric" : "_Metric"
                tabs = ["FVS_Summary", "FVS_TreeList", "FVS_CutList", "FVS_ATRTList"] .* sfx
                jtabs = [t.name for t in SQLite.tables(jl)]
                for t in tabs
                    @test t in jtabs
                    @test _mcols(jl, t) == _mcols(live, t)                     # names, order, declared types
                end
                @test !any(t -> t in jtabs, ("FVS_Summary", "FVS_TreeList", "FVS_CutList", "FVS_ATRTList"))

                # .sum inventory row (sumout.f metric FORMAT 20) — byte-for-byte fields vs the live .out table
                lrow = split(first(readlines(joinpath(fx, "$(stem)_live_sum.txt"))))
                jrow = split(only(l for l in split(sumtxt, '\n') if startswith(l, string(yr0))))
                @test jrow == lrow

                srows(db) = [collect(values(NamedTuple(x)))[2:end] for x in DBInterface.execute(db,
                             "SELECT * FROM FVS_Summary$sfx WHERE Year = $yr0")]
                if east
                    @test isempty(srows(live)) && isempty(srows(jl))            # FVS INSERT/CREATE mismatch ⇒ 0 rows
                else
                    @test srows(jl) == srows(live)
                end

                vol = ("TCuM", "MCuM", "CCum", "CCuM")
                for t in tabs[2:4]
                    L = _mrows(live, t, yr0); J = _mrows(jl, t, yr0)
                    @test sort(collect(keys(J))) == sort(collect(keys(L)))
                    for (i, lr) in L, (c, lv) in lr
                        (c === :CaseID || !haskey(J, i)) && continue
                        jv = J[i][c]
                        if String(c) in vol
                            @test isapprox(jv, lv; rtol = 5e-5, atol = 1e-6)
                        elseif c === :DG && t != tabs[2]
                            continue                                         # open item (post-calibration DG)
                        else
                            @test jv == lv
                        end
                    end
                end
            finally
                SQLite.close(live); SQLite.close(jl)
            end
        end
    end
end
