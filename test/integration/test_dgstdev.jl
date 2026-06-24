# test_dgstdev.jl — DGSTDEV keyword (option 57) → DGSD bound on stochastic DG variation vs Fortran.
#
# DGSTDEV sets DGSD, the std-dev bound on the serial-correlation diameter-growth variation
# (default 2.0; DGSD<1 turns the random variation OFF → deterministic DG). Before the fix
# the keyword was unrecognized and DGSD was hardcoded 2.0. The scenario sets DGSTDEV 0
# (deterministic). Checks: SETS control.dg_stddev_bound; FIRES (≠ default); FORTRAN bit-exact.

using Test, FVSjl

const _DG_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_dg_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_dg_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]

@testset "DGSTDEV → DGSD stochastic-DG bound" begin
    # 1. SETS — and default is 2.0 (grinit.f)
    s = FVSjl.StandState(FVSjl.Southern())
    @test s.control.dg_stddev_bound == 2f0
    FVSjl.kw_dgstdev!(s, FVSjl.KeywordRecord("DGSTDEV ", "", ["0", fill("", 11)...],
                      zeros(Float32, 12), [true, falses(11)...], 12, FVSjl.KW_OK, 0))
    @test s.control.dg_stddev_bound == 0f0

    key = joinpath(_DG_DIR, "dgstdev.key"); sav = joinpath(_DG_DIR, "dgstdev.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "dgstdev scenario not available"
    else
        jl = _dg_rows(FVSjl.run_keyfile(key; faithful = true))
        # 2. FIRES — DGSTDEV 0 (deterministic) differs from the same stand with default DGSD=2
        offkey = joinpath(_DG_DIR, "_dgstdev_off.key")
        cp(joinpath(_DG_DIR, "dgstdev.tre"), joinpath(_DG_DIR, "_dgstdev_off.tre"); force = true)
        open(offkey, "w") do io
            for l in eachline(key); startswith(strip(l), "DGSTDEV") || println(io, l); end
        end
        try
            off = _dg_rows(FVSjl.run_keyfile(offkey; faithful = true))
            @test any(jl[i][10] != off[i][10] for i in eachindex(jl))
        finally
            rm(offkey; force = true); rm(joinpath(_DG_DIR, "_dgstdev_off.tre"); force = true)
        end
        # 3. FORTRAN — deterministic DG matches live Fortran exactly (no RNG noise)
        ft = _dg_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            for c in (3, 4, 5, 6, 7); @test parse(Int, j[c]) == parse(Int, f[c]); end
            @test abs(parse(Float64, j[8]) - parse(Float64, f[8])) <= 0.05
        end
    end
end
