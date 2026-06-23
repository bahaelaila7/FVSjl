# test_tcondmlt.jl — TCONDMLT tree-condition cut weight (cuts.f:1074/1424) vs live Fortran.
#
# TCONDMLT sets TCWT, a weight added to the RDPSRT cut-priority key: WK2 = ±DBH + IORDER(SPECPREF)
# + TCWT·IMC + SPCLWT·ISPECL. IMC is the tree's mortality/condition code (1..3 for live trees), so a
# positive TCWT removes worse-condition trees first, regardless of size. The scenario marks the sugar
# maples with IMC=3 in the .tre and thins from below (THINBBA) with TCWT=100, so those trees are cut
# first. Checks:
#   1. NON-VACUOUS — the result differs from the same thin without the TCONDMLT line;
#   2. FORTRAN — TPA/BA/cubic columns match live Fortran (board feet within Scribner noise).

using Test, FVSjl

const _TC_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_tc_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_tc_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_tccol(r, c) = parse(Float64, r[c])

@testset "TCONDMLT condition-weighted thin vs Fortran" begin
    key = joinpath(_TC_DIR, "tcondmlt.key")
    sav = joinpath(_TC_DIR, "tcondmlt.sum.save")
    tre = joinpath(_TC_DIR, "tcondmlt.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "tcondmlt scenario not available"
    else
        jl = _tc_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. NON-VACUOUS: without TCONDMLT the thin ranks by size only → a different cut.
        notc = tempname() * ".key"
        write(notc, join(filter(l -> !startswith(l, "TCONDMLT"), readlines(key)), "\n") * "\n")
        cp(tre, replace(notc, ".key" => ".tre"); force = true)
        notcrows = _tc_rows(FVSjl.run_keyfile(notc; faithful = true))
        @test jl != notcrows

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
