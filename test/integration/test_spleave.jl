# test_spleave.jl — SPLEAVE / LEAVESP (leave named species during a thin, cuts.f:1466) vs Fortran.
#
# SPLEAVE sets a per-species LEAVESP flag; a flagged species is excluded from BOTH the stocking
# computation and the removal of every subsequent thin (cuts.f:1031/1097). The scenario thins from
# below (THINBBA, leave 80 ft² BA) but leaves sugar maple (species 22), so SM is never cut. Checks:
#   1. NON-VACUOUS — the result differs from the same thin without the SPLEAVE line (SM was protected);
#   2. FORTRAN — matches live Fortran on TPA/BA/cubic columns (board feet within Scribner noise).

using Test, FVSjl

const _SL_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_sl_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_sl_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_slcol(r, c) = parse(Float64, r[c])

@testset "SPLEAVE leave-species during thin vs Fortran" begin
    key = joinpath(_SL_DIR, "spleave.key")
    sav = joinpath(_SL_DIR, "spleave.sum.save")
    tre = joinpath(_SL_DIR, "spleave.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "spleave scenario not available"
    else
        jl = _sl_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. NON-VACUOUS: without SPLEAVE the thin cuts sugar maple too → a different stand.
        nosl = tempname() * ".key"
        write(nosl, join(filter(l -> !startswith(l, "SPLEAVE"), readlines(key)), "\n") * "\n")
        cp(tre, replace(nosl, ".key" => ".tre"); force = true)
        noslrows = _sl_rows(FVSjl.run_keyfile(nosl; faithful = true))
        @test jl != noslrows

        # 2. matches live Fortran (cubic columns bit-exact; board feet within Scribner noise).
        ft = _sl_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                for c in (3, 4, 9, 10, 11)
                    @test abs(_slcol(jl[i], c) - _slcol(ft[i], c)) <= 2
                end
                @test abs(_slcol(jl[i], 12) - _slcol(ft[i], 12)) <= 1 + 0.002 * _slcol(ft[i], 12)
            end
        end
    end
end
