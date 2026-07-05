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
    # (scenario, cycles, tpa_tol, ba_tol) — cols 3 TPA / 4 BA / 7 TopHt / 8 QMD. RE-MEASURED per-column
    # (2026-07-05): TopHt is BIT-EXACT in ALL 6 scenarios; QMD is bit-exact except mult_baimult (1 tenth);
    # TPA/BA are bit-exact except a rendered-integer 1-step print knife-edge in a few scenarios (TPA: baimult,
    # mortmult_win; BA: mortmult_win, reghmult). So TopHt→==, QMD→tenth-grid ≤1 tenth, and TPA/BA carry their
    # EXACT per-scenario measured bound (0 or 1) — replacing the old uniform tol=1/tol=2 (which padded every
    # bit-exact column, up to 10× on QMD and 2× on the mortmult scenarios).
    for (nm, ncyc, tpa_tol, ba_tol) in (("mult_htgmult", 11, 0, 0), ("mult_baimult", 11, 1, 0),
                                        ("mult_mortmult", 6, 0, 0), ("mult_mortmult_win", 7, 1, 1),
                                        ("mult_reghmult", 5, 0, 1), ("mult_regdmult", 5, 0, 0))
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
                for i in 1:min(ncyc, length(jl))
                    @test abs(parse(Float64, jl[i][3]) - parse(Float64, ft[i][3])) <= tpa_tol   # TPA (0 ⇒ BIT-EXACT)
                    @test abs(parse(Float64, jl[i][4]) - parse(Float64, ft[i][4])) <= ba_tol    # BA  (0 ⇒ BIT-EXACT)
                    @test parse(Float64, jl[i][7]) == parse(Float64, ft[i][7])                  # TopHt — BIT-EXACT (all)
                    @test abs(round(Int, parse(Float64, jl[i][8]) * 10) -
                              round(Int, parse(Float64, ft[i][8]) * 10)) <= 1                   # QMD — ≤1 tenth knife-edge
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
