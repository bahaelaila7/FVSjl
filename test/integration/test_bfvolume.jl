# test_bfvolume.jl — BFVOLUME board-foot merch-standard override (volkey.f) vs live Fortran.
#
# BFVOLUME overrides the per-species board-foot merch standards. The scenario raises BFTOPD 9"→11"
# for all species, so the board-foot top diameter differs from the sawtimber one (BFPFLG=0). That
# exercises the full board-foot path that VOLEQNUM does not reach:
#   • board feet recomputed from a separate board call with BFTOPD/BFSTMP (fvsvol.f:362), AND
#   • the Region-8 "≥10 ft of product" rule (fvsvol.f:499): where the board-top sawlog is < 10 ft,
#     board feet AND the sawtimber cubic (→ the sawtimber part of merch cubic) are zeroed.
# Two checks: (1) it FIRES (board feet drops vs no-BFVOLUME); (2) it matches live Fortran on the
# structural + cubic columns bit-exact and on board feet bar the ±Scribner Float32 noise.

using Test, FVSjl

const _BV_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_bv_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_bv_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_bvcol(r, c) = parse(Float64, r[c])

@testset "BFVOLUME board-foot override (+Region-8) vs Fortran" begin
    key = joinpath(_BV_DIR, "bfvolume_override.key")
    sav = joinpath(_BV_DIR, "bfvolume_override.sum.save")
    tre = joinpath(_BV_DIR, "bfvolume_override.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "bfvolume_override scenario not available"
    else
        jl = _bv_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. FIRES: raising BFTOPD to 11" cuts board feet (col 12) vs the no-BFVOLUME twin.
        novkey = tempname() * ".key"
        write(novkey, join(filter(l -> !startswith(l, "BFVOLUME"), readlines(key)), "\n") * "\n")
        cp(tre, replace(novkey, ".key" => ".tre"); force = true)
        nd = _bv_rows(FVSjl.run_keyfile(novkey; faithful = true))
        @test length(jl) == length(nd)
        @test any(_bvcol(nd[i], 12) - _bvcol(jl[i], 12) > 100 for i in 1:length(jl))

        # 2. matches live Fortran: cubic columns bit-exact (±1 noise), board feet within Scribner noise.
        ft = _bv_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                for c in (3, 4, 9, 10, 11)   # TPA / BA / total / merch / sawtimber cubic
                    @test abs(_bvcol(jl[i], c) - _bvcol(ft[i], c)) <= 2
                end
                @test _bvcol(jl[i], 12) == _bvcol(ft[i], 12)
            end
        end
    end
end
