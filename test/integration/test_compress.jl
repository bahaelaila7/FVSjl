# test_compress.jl — COMPRESS keyword recognition + scheduling (initre.f:8000, option 78).
#
# COMPRESS schedules act 250 (reduce the tree list to NCLAS classes). The clustering
# ALGORITHM (comprs.f, 1010-line IBM-SSP-eigen PCA) is a tracked chunk (docs/COMPRESS_
# chunk_plan.md) and NOT yet ported — so this tests only that the keyword is RECOGNIZED
# and scheduled (no longer silently dropped) and that a COMPRESS stand still runs cleanly
# (records pass through uncompressed until the algorithm lands). Do NOT mark COMPRESS done
# on the strength of this test — it does not exercise the compression.

using Test, FVSjl

@testset "COMPRESS keyword recognition + scheduling" begin
    mkrec(fields, vals, present) =
        FVSjl.KeywordRecord("COMPRESS", "", fields, vals, present, 12, FVSjl.KW_OK, 0)

    # explicit params: date 1990, target 200, pn1 60
    s = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_compress!(s, mkrec(["1990", "200", "60", fill("", 9)...],
                                Float32[1990, 200, 60, zeros(Float32, 9)...],
                                [true, true, true, falses(9)...]))
    a = last(s.control.schedule)
    @test a.icflag == Int32(250)
    @test a.year == Int32(1990)
    @test a.params[1] == 200f0 && a.params[2] == 60f0

    # defaults: date 1, target MAXTRE/2 = 1500, pn1 50
    s2 = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_compress!(s2, mkrec(fill("", 12), zeros(Float32, 12), falses(12)))
    d = last(s2.control.schedule)
    @test d.year == Int32(1) && d.params[1] == 1500f0 && d.params[2] == 50f0

    # COMPRESS is no longer in the no-op set (it is a real, now-recognized effect)
    @test !("COMPRESS" in FVSjl.KNOWN_NOOP)
    @test !("COMPRESS" in FVSjl.variant_noop_keywords(FVSjl.Southern()))

    # a stand carrying COMPRESS still runs cleanly (the scheduled act 250 is skipped by
    # cuts! PASS 2 — records pass through uncompressed until the algorithm is ported)
    dir = joinpath(@__DIR__, "..", "harness", "scenarios")
    key = joinpath(dir, "fire_early.key")
    if isfile(key)
        off = FVSjl.run_keyfile(key; faithful = true)
        cmpkey = joinpath(dir, "_compress_on.key")
        cp(joinpath(dir, "fire_early.tre"), joinpath(dir, "_compress_on.tre"); force = true)
        open(cmpkey, "w") do io
            for l in eachline(key)
                println(io, l)
                strip(l) == "INVYEAR       1990.0" && println(io, "COMPRESS        1990      100.")
            end
        end
        try
            on = FVSjl.run_keyfile(cmpkey; faithful = true)
            @test length(split(on, "\n")) == length(split(off, "\n"))   # runs, same shape
        finally
            rm(cmpkey; force = true); rm(joinpath(dir, "_compress_on.tre"); force = true)
        end
    end
end
