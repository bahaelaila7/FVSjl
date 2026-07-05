# test_sprout_regen.jl — ESUCKR stump/root-sprout regeneration vs live Fortran.
#
# The scenario heavily thins (THINBBA, leave 40 ft² BA) a stand containing sugar maple
# (a sprouting species) with sprouting enabled for SM via an ESTAB…END / SPROUT block.
# After the 2000 thin the cut stumps regenerate as sprouts (esuckr.f), so the next
# cycle's TPA jumps back up. Checks:
#   1. NON-VACUOUS — removing the SPROUT line gives a different (much lower) post-thin
#      stand, i.e. the sprouts are actually being created;
#   2. FORTRAN — matches live Fortran on TPA/BA/cubic columns (board feet within
#      Scribner noise), including the sprout-regeneration cycle.

using Test, FVSjl

const _SP_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_sp_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_sp_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_spcol(r, c) = parse(Float64, r[c])

@testset "ESUCKR stump-sprout regeneration vs Fortran" begin
    key = joinpath(_SP_DIR, "sprout.key")
    sav = joinpath(_SP_DIR, "sprout.sum.save")
    tre = joinpath(_SP_DIR, "sprout.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "sprout scenario not available"
    else
        jl = _sp_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. NON-VACUOUS: without the SPROUT line the cut stumps don't regenerate.
        nosp = tempname() * ".key"
        write(nosp, join(filter(l -> !startswith(strip(l), "SPROUT"), readlines(key)), "\n") * "\n")
        cp(tre, replace(nosp, ".key" => ".tre"); force = true)
        nosprows = _sp_rows(FVSjl.run_keyfile(nosp; faithful = true))
        @test jl != nosprows
        # the post-thin (2005) TPA must be much higher with sprouting than without
        post(rows) = (i = findfirst(r -> r[1] == "2005", rows); i === nothing ? 0.0 : _spcol(rows[i], 3))
        @test post(jl) > post(nosprows) + 100      # sprouts add hundreds of TPA

        # 2. matches live Fortran — BIT-EXACT on every structural + volume column (measured Δ=0
        # across all 11 rows for TPA/BA/cubic AND board feet; the old `<=1.5` was pure padding —
        # the sprout-regen cycle reproduces live's records exactly, so the .sum renders identically).
        ft = _sp_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                for c in (3, 4, 5, 6, 9, 10, 12)    # TPA, QMD-cols, BA, cubic volumes, board feet
                    @test _spcol(jl[i], c) == _spcol(ft[i], c)
                end
            end
        end
    end
end
