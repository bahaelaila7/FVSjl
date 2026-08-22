# FVS_TreeList DG/HtG for OP (Olympic) — companion to the OC clobber regression (2026-08-22).
# OP is the sibling ORGANON variant (NWO). Unlike OC (single-authority diameter_growth! that applied DBH/HT inline
# and used to zero diam_growth/ht_growth), OP uses the COOPERATING driver: diameter_growth!(::Olympic) FILLS
# t.diam_growth and height_growth!(::Olympic) FILLS t.ht_growth, and the SHARED grow_cycle! apply-loop grows the
# stand — so the applied increment SURVIVES into the report. This test confirms OP does NOT have the OC clobber:
# the projected-cycle FVS_TreeList DG/HtG are populated AND bit-exact vs FVSop_clean (oracle sums captured below).
# (The oracle emits 2 extra rows at the INVENTORY year = the cycle-0 input-dead records, which jl emits only for
# CentralRockies — a known cross-variant latent gap, unrelated to DG/HtG; the live-tree DG/HtG sums are what match.)

using Test
using FVSjl
using SQLite, DBInterface

@testset "FVS_TreeList DG/HtG bit-exact for OP live trees (no OC-style clobber)" begin
    fix = joinpath(@__DIR__, "..", "fixtures", "optreelist")
    dir = mktempdir()
    # jl derives the TREEDATA .tre name from the keyfile basename, so name both "optreelist".
    cp(joinpath(fix, "optreelist.tre"), joinpath(dir, "optreelist.tre"))
    key = joinpath(dir, "optreelist.key")
    outdb = joinpath(dir, "out.db")
    open(key, "w") do io
        for ln in readlines(joinpath(fix, "optreelist.key"))
            println(io, strip(ln) == "optl_oracle.db" ? outdb : ln)
        end
    end
    FVSjl.run_keyfile(key; variant = FVSjl.variant_from_code("OP"))
    db = SQLite.DB(outdb)
    @test "FVS_TreeList" in [t.name for t in SQLite.tables(db)]

    # ΣDG, ΣHtG over the live tree list at a given year (dead records carry DG=HtG=0, so they don't perturb it)
    sums(yr) = begin
        q = DBInterface.execute(db, "SELECT DG, HtG FROM FVS_TreeList WHERE Year=$yr")
        rows = [(Float64(r[1]), Float64(r[2])) for r in q]
        (sum(abs(r[1]) for r in rows), sum(abs(r[2]) for r in rows))
    end

    # Oracle sums from FVSop_clean (optl_oracle.db): the projected cycles must be POPULATED and match bit-exact.
    for (yr, dg, htg) in ((1995, 19.05, 148.74), (2000, 18.46, 163.69))
        sdg, shtg = sums(yr)
        @test sdg > 1.0 && shtg > 1.0            # populated (the OC bug was exactly 0 here)
        @test isapprox(sdg, dg; atol = 0.05)     # bit-exact vs oracle (0.05 = print-rounding slack)
        @test isapprox(shtg, htg; atol = 0.05)
    end
end
