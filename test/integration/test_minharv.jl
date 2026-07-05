# test_minharv.jl — MINHARV minimum-harvest thresholds (cuts.f:400/1556) vs live Fortran.
#
# MINHARV sets per-cycle harvest minimums (BA / total / merch / sawlog cubic / board feet). After a
# cycle's thinning is computed, if the TOTAL removal falls below ANY threshold the whole cut is
# CANCELED — the pre-thin stand is restored. The scenario schedules a THINBBA (leave 80 ft² BA, which
# removes < 100) together with MINHARV BAMIN=100, so the thin is canceled. Checks:
#   1. CANCELED — the run is identical to the same stand with the THINBBA removed (no thin at all);
#   2. NON-VACUOUS — without the MINHARV line the THINBBA does change the stand (so it really fired);
#   3. FORTRAN — matches live Fortran on TPA/BA/cubic columns (board feet within Scribner noise).

using Test, FVSjl

const _MH_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_mh_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_mh_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_mhcol(r, c) = parse(Float64, r[c])

@testset "MINHARV minimum-harvest gate vs Fortran" begin
    key = joinpath(_MH_DIR, "minharv.key")
    sav = joinpath(_MH_DIR, "minharv.sum.save")
    tre = joinpath(_MH_DIR, "minharv.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "minharv scenario not available"
    else
        jl = _mh_rows(FVSjl.run_keyfile(key; faithful = true))
        lines = readlines(key)

        # 1. CANCELED: dropping the THINBBA too gives the same result (the gate canceled the cut).
        noth = tempname() * ".key"
        write(noth, join(filter(l -> !startswith(l, "MINHARV") && !startswith(l, "THINBBA"), lines), "\n") * "\n")
        cp(tre, replace(noth, ".key" => ".tre"); force = true)
        nothrows = _mh_rows(FVSjl.run_keyfile(noth; faithful = true))
        @test jl == nothrows

        # 2. NON-VACUOUS: without MINHARV, the THINBBA actually changes the stand.
        nomh = tempname() * ".key"
        write(nomh, join(filter(l -> !startswith(l, "MINHARV"), lines), "\n") * "\n")
        cp(tre, replace(nomh, ".key" => ".tre"); force = true)
        nomhrows = _mh_rows(FVSjl.run_keyfile(nomh; faithful = true))
        @test jl != nomhrows

        # 3. matches live Fortran (cubic columns bit-exact; board feet within Scribner noise).
        ft = _mh_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                @test _mhcol(jl[i], 3) == _mhcol(ft[i], 3)   # TPA — BIT-EXACT
                @test _mhcol(jl[i], 4) == _mhcol(ft[i], 4)   # BA  — BIT-EXACT
                @test _mhcol(jl[i], 11) == _mhcol(ft[i], 11)   # SCuFt — BIT-EXACT (measured Δ0 all cycles)
            end
            # TCuFt (9) / MCuFt (10): bit-exact bar a print-boundary ULP — residual is the non-associative
            # Float32 tree-SUM accumulation order (doctrine #9: exposed, not a passing ≤1).
            @test_broken all(_mhcol(jl[i], 9)  == _mhcol(ft[i], 9)  for i in 1:length(jl))  # TCuFt — non-associative tree-SUM order
            @test_broken all(_mhcol(jl[i], 10) == _mhcol(ft[i], 10) for i in 1:length(jl))  # MCuFt — non-associative tree-SUM order
            # Board feet (12): BIT-EXACT every cycle bar a single print-boundary ULP (the per-acre Scribner
            # sum lands within one ULP of the +0.5 integer-render knife-edge); residual is the non-associative
            # Float32 tree-SUM accumulation order (doctrine #9: exposed, not a passing ≤1).
            @test_broken all(_mhcol(jl[i], 12) == _mhcol(ft[i], 12) for i in 1:length(jl))  # BdFt — non-associative tree-SUM order
        end
    end
end
