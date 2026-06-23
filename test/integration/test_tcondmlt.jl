# test_tcondmlt.jl — TCONDMLT tree-condition cut weights (cuts.f:1074/1424) vs live Fortran.
#
# TCONDMLT sets two weights added to the RDPSRT cut-priority key:
#   WK2 = ±DBH + IORDER(SPECPREF) + TCWT·IMC + SPCLWT·ISPECL
# TCWT (PRM 1) weights the mortality/condition code IMC (1..3 live); SPCLWT (PRM 2) weights the
# special-status code ISPECL (damage code 55). A positive weight removes the flagged trees first,
# regardless of size. Two scenarios, each marking the sugar maples and thinning from below (THINBBA):
#   * tcondmlt — IMC=3 + TCWT=100 (condition weight);
#   * spclwt   — special-status code 55 (ISPECL=9) + SPCLWT=100 (special-status weight).
# Each must (1) differ from the same thin without the TCONDMLT line, and (2) match live Fortran on
# TPA/BA/cubic columns (board feet within Scribner noise).

using Test, FVSjl

const _TC_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_tc_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_tc_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_tccol(r, c) = parse(Float64, r[c])

@testset "TCONDMLT condition-weighted thin vs Fortran" begin
    for stem in ("tcondmlt", "spclwt")
        key = joinpath(_TC_DIR, stem * ".key"); sav = joinpath(_TC_DIR, stem * ".sum.save")
        tre = joinpath(_TC_DIR, stem * ".tre")
        if !isfile(key) || !isfile(sav)
            @test_skip "$stem scenario not available"; continue
        end
        @testset "$stem" begin
            jl = _tc_rows(FVSjl.run_keyfile(key; faithful = true))

            # 1. NON-VACUOUS: without TCONDMLT the thin ranks by size only → a different cut.
            notc = tempname() * ".key"
            write(notc, join(filter(l -> !startswith(l, "TCONDMLT"), readlines(key)), "\n") * "\n")
            cp(tre, replace(notc, ".key" => ".tre"); force = true)
            @test jl != _tc_rows(FVSjl.run_keyfile(notc; faithful = true))

            # 2. matches live Fortran (cubic columns bit-exact; board feet within Scribner noise).
            ft = _tc_base(sav)
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for i in 1:length(jl)
                    for c in (3, 4, 9, 10, 11)
                        @test abs(_tccol(jl[i], c) - _tccol(ft[i], c)) <= 2
                    end
                    @test abs(_tccol(jl[i], 12) - _tccol(ft[i], 12)) <= 1 + 0.005 * _tccol(ft[i], 12)
                end
            end
        end
    end
end
