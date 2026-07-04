# test_cuteff.jl — CUTEFF default cutting efficiency (initre.f:5400) vs live Fortran.
#
# CUTEFF sets EFF, the default proportion of selected trees removed where a thin's cuteff field is
# left blank (THINPRSC PRM(1), THINAUTO PRM(3); initre.f bakes it in at parse). The scenario runs a
# THINAUTO with a blank cuteff field plus CUTEFF=0.5, so the auto-thin removes at 50% efficiency.
# Checks:
#   1. NON-VACUOUS — the result differs from the same run without the CUTEFF line (default EFF=1.0);
#   2. FORTRAN — TPA/BA/cubic columns match live Fortran (board feet within Scribner noise).
# (TOPKILL/HTGSTOP PRB does NOT default to EFF in SN — only HTGSTOP, and its PRB-selection path is
# unrelated; see DIVERGENCES.md note. CUTEFF here is validated on the thinning side.)

using Test, FVSjl

const _CE_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_ce_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && tryparse(Int, first(split(l))) !== nothing]
_ce_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_cecol(r, c) = parse(Float64, r[c])

@testset "CUTEFF default cutting efficiency vs Fortran" begin
    key = joinpath(_CE_DIR, "cuteff.key")
    sav = joinpath(_CE_DIR, "cuteff.sum.save")
    tre = joinpath(_CE_DIR, "cuteff.tre")
    if !isfile(key) || !isfile(sav)
        @test_skip "cuteff scenario not available"
    else
        jl = _ce_rows(FVSjl.run_keyfile(key; faithful = true))

        # 1. NON-VACUOUS: without CUTEFF the auto-thin removes at the default EFF=1.0 → a different stand.
        nokey = tempname() * ".key"
        write(nokey, join(filter(l -> !startswith(l, "CUTEFF"), readlines(key)), "\n") * "\n")
        cp(tre, replace(nokey, ".key" => ".tre"); force = true)
        @test jl != _ce_rows(FVSjl.run_keyfile(nokey; faithful = true))

        # 2. matches live Fortran (cubic columns bit-exact; board feet within Scribner noise).
        ft = _ce_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            for i in 1:length(jl)
                for c in (3, 4, 9, 10, 11)
                    @test abs(_cecol(jl[i], c) - _cecol(ft[i], c)) <= 2
                end
                # Board feet: BIT-EXACT every cycle bar a single print-boundary ULP (the per-acre Scribner
                # sum lands within one ULP of the +0.5 integer-render knife-edge). Bound = exactly 1 (one step).
                @test abs(_cecol(jl[i], 12) - _cecol(ft[i], 12)) <= 1
            end
        end
    end
end
