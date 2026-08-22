# FVS_TreeList DG/HtG for OC (regression for the ORGANON diam_growth/ht_growth clobber, 2026-08-22).
# OC's diameter_growth!(::OregonCoast) applies DBH/HT INLINE (ORGANON does growth+mort+crown in one EXECUTE) and
# used to ZERO t.diam_growth/t.ht_growth so the shared grow_cycle! apply-loop wouldn't double-apply — but that
# left the FVS_TreeList DG/HtG columns 0 for every PROJECTED cycle (the inventory cycle was fine). CR (Wykoff) was
# always bit-exact, so it was OC-specific. Fix: the shared apply-loop SKIPS OC (it already added 0 there, so this is
# .sum-inert) and OC KEEPS the applied increment for the report. This test asserts the projected-cycle DG/HtG are
# now POPULATED (the bug was exactly 0); the per-tree values are bit-exact-or-cornered on the OC broken-top knife-edge.

using Test
using FVSjl
using SQLite, DBInterface

@testset "FVS_TreeList DG/HtG populated for OC projected cycles (was 0)" begin
    fix = joinpath(@__DIR__, "..", "fixtures", "octreelist")
    dir = mktempdir()
    cp(joinpath(fix, "octreelist.tre"), joinpath(dir, "octreelist.tre"))
    key = joinpath(dir, "octreelist.key")
    outdb = joinpath(dir, "out.db")
    open(key, "w") do io
        for ln in readlines(joinpath(fix, "octreelist.key"))
            println(io, strip(ln) == "octreelist_oracle.db" ? outdb : ln)
        end
    end
    FVSjl.run_keyfile(key; variant = FVSjl.variant_from_code("OC"))
    db = SQLite.DB(outdb)
    @test "FVS_TreeList" in [t.name for t in SQLite.tables(db)]

    dghtg(yr) = begin
        s = 0.0
        for r in DBInterface.execute(db, "SELECT DG, HtG FROM FVS_TreeList WHERE Year=$yr")
            s += abs(Float64(r.DG)) + abs(Float64(r.HtG))
        end
        s
    end
    # inventory cycle was always populated; the REGRESSION is the projected cycles being 0
    @test dghtg(1998) > 100.0     # oracle ≈ 171.8 (was 0.0 before the fix)
    @test dghtg(2003) > 100.0     # oracle ≈ 172.5 (was 0.0 before the fix)
end
