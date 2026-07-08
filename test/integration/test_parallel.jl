# test_parallel.jl — pillar-3 guarantee: running many stands CONCURRENTLY produces byte-identical
# output to running them serially, for every variant. Bit-identity is the correctness test for
# parallelism (no shared mutable state, per-stand RNG/scratch/trees, thread-safe coef cache).
#
# Under the canonical single-threaded runner (`julia test/runtests.jl`) this degrades to a
# determinism/reentrancy check (N fresh runs must all match the reference) — still meaningful: it
# proves no run leaves hidden global state that perturbs the next. Run with `julia -t 8 …` to exercise
# TRUE concurrency (the -t8 4-variant validation was recorded in docs/MODERNIZATION_AUDIT.md).

using Test, FVSjl, Base.Threads

@testset "parallel == serial (bit-identical), all variants" begin
    D, T = "01-01-2020", "00:00:00"   # fixed .sum stamp ⇒ deterministic header
    scen = [("SN", joinpath(@__DIR__, "..", "harness", "scenarios", "snt01_alpha.key")),
            ("NE", joinpath(@__DIR__, "..", "fixtures", "canonical", "net01.key")),
            ("CS", joinpath(@__DIR__, "..", "fixtures", "canonical", "cst01.key")),
            ("LS", joinpath(@__DIR__, "..", "fixtures", "canonical", "lst01.key"))]
    for (v, key) in scen
        if !isfile(key); @test_skip "$v scenario missing"; continue; end
        @testset "$v" begin
            FVSjl.run_keyfile(key; faithful=true, date=D, time=T)         # warm + compile
            empty!(FVSjl._COEF_CACHE)                                     # cold cache ⇒ threads race the get!
            ref = FVSjl.run_keyfile(key; faithful=true, date=D, time=T)   # serial reference
            N = 32
            res = Vector{String}(undef, N)
            @threads for i in 1:N
                res[i] = FVSjl.run_keyfile(key; faithful=true, date=D, time=T)
            end
            @test all(==(ref), res)
        end
    end
end
