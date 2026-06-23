# test_multipliers.jl — growth/mortality keyword multipliers (MULTS) vs live Fortran.
#
# BAIMULT (diameter-growth), HTGMULT (height-growth), and MORTMULT (background mortality)
# are per-species dated multipliers (base/mults.f). Each scenario applies a 1.5–2.0×
# multiplier from the inventory year; the .sum.save baselines are live-Fortran output for
# the exact same keys. FVSjl must reproduce them, AND the multiplier must visibly act
# (guards against a silent no-op).
#   * mult_htgmult / mult_baimult — base stand (s29 records, no thin), bit-exact every cycle.
#   * mult_mortmult — a bare PLANT stand whose early cycles are background-mortality (where
#     MORTMULT applies); validated on those cycles (later cycles carry the regen density tail).
#   * mult_mortmult_win — same stand with a DBH-windowed MORTMULT (3× for DBH<4", morts.f:518):
#     the multiplier kills the small planted trees, then stops once they grow past 4" (≈2012),
#     so its TPA diverges from the windowless 3× — guards the per-tree window logic.

using Test, FVSjl

const _MULT_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_mult_rows(txt) = [split(l) for l in split(txt, "\n")
                   if length(split(l)) >= 11 && tryparse(Int, first(split(l))) !== nothing]
_mult_base(path) = [split(l) for l in eachline(path)
                    if length(split(l)) >= 11 &&
                       (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]

@testset "growth/mortality multipliers (MULTS) vs Fortran" begin
    # (scenario, cycles validated bit-exact) — cols 3 TPA / 4 BA / 7 TopHt / 8 QMD.
    # mult_reghmult/mult_regdmult are bare-stand regen multipliers: validated on the early
    # cycles where the small-tree REGENT model (and thus XRHGRO/XRDGRO) is in effect.
    for (nm, ncyc) in (("mult_htgmult", 11), ("mult_baimult", 11), ("mult_mortmult", 6),
                       ("mult_mortmult_win", 7), ("mult_reghmult", 5), ("mult_regdmult", 5))
        key  = joinpath(_MULT_DIR, nm * ".key")
        base = joinpath(_MULT_DIR, nm * ".sum.save")
        if !isfile(key) || !isfile(base)
            @test_skip "$nm scenario not available"; continue
        end
        @testset "$nm" begin
            jl = _mult_rows(FVSjl.run_keyfile(key; faithful = true))
            ft = _mult_base(base)
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for i in 1:min(ncyc, length(jl)), c in (3, 4, 7, 8)
                    @test abs(parse(Float64, jl[i][c]) - parse(Float64, ft[i][c])) <= 1
                end
            end
        end
    end

    # the multipliers must visibly act (base stand reaches ~TopHt 79 / QMD 15 by 2040)
    if isfile(joinpath(_MULT_DIR, "mult_htgmult.key"))
        h = _mult_rows(FVSjl.run_keyfile(joinpath(_MULT_DIR, "mult_htgmult.key"); faithful = true))
        @test parse(Float64, h[end][7]) >= 88     # HTGMULT 1.5 ⇒ TopHt ~90
        b = _mult_rows(FVSjl.run_keyfile(joinpath(_MULT_DIR, "mult_baimult.key"); faithful = true))
        @test parse(Float64, b[end][8]) >= 22     # BAIMULT 1.5 ⇒ QMD ~23
    end
end
