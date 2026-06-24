# test_rannseed.jl — RANNSEED keyword (option 61) → reseed the main RNG stream vs Fortran.
#
# RANNSEED installs a new seed for the main random stream (RANSED), which drives the
# stochastic DGSCOR serial-correlation diameter growth + mortality. Before the fix the
# keyword was unrecognized (ignored) so the stand kept the default 55329 seed. Checks:
#   1. SETS — the handler reseeds rng.ss (forced odd, RANSED);
#   2. FIRES — the run differs from the same stand at the default seed;
#   3. FORTRAN — TPA/BA/SDI/CCF/TopHt/QMD match live Fortran (± DGSCOR/Scribner noise).

using Test, FVSjl

const _RS_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_rs_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_rs_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]

@testset "RANNSEED → reseed main RNG stream" begin
    # 1. SETS — a present non-zero field installs the (odd-forced) seed; blank restarts
    s = FVSjl.StandState(FVSjl.Southern())
    rec(v, present) = FVSjl.KeywordRecord("RANNSEED", "", [string(v), fill("", 11)...],
                                          Float32[v, zeros(Float32, 11)...],
                                          [present, falses(11)...], 12, FVSjl.KW_OK, 0)
    FVSjl.kw_rannseed!(s, rec(12345f0, true))
    @test s.rng.ss == 12345f0          # odd already → unchanged
    FVSjl.kw_rannseed!(s, rec(1000f0, true))
    @test s.rng.ss == 1001f0           # even → forced odd (RANSED)

    key = joinpath(_RS_DIR, "rannseed.key"); sav = joinpath(_RS_DIR, "rannseed.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "rannseed scenario not available"
    else
        jl = _rs_rows(FVSjl.run_keyfile(key; faithful = true))
        # 2. FIRES — dropping RANNSEED (default 55329 seed) changes the stochastic stream
        offkey = joinpath(_RS_DIR, "_rann_off.key")
        cp(joinpath(_RS_DIR, "rannseed.tre"), joinpath(_RS_DIR, "_rann_off.tre"); force = true)
        open(offkey, "w") do io
            for l in eachline(key); startswith(strip(l), "RANNSEED") || println(io, l); end
        end
        try
            off = _rs_rows(FVSjl.run_keyfile(offkey; faithful = true))
            @test any(jl[i][10] != off[i][10] for i in eachindex(jl))   # MerchCuFt differs somewhere
        finally
            rm(offkey; force = true); rm(joinpath(_RS_DIR, "_rann_off.tre"); force = true)
        end
        # 3. FORTRAN — TPA/BA/SDI/CCF/TopHt/QMD match live Fortran
        ft = _rs_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            for c in (3, 4, 5, 6, 7); @test parse(Int, j[c]) == parse(Int, f[c]); end
            @test abs(parse(Float64, j[8]) - parse(Float64, f[8])) <= 0.05
        end
    end
end
