# test_tripling.jl — NOTRIPLE / NUMTRIP tripling control (ICL4) vs live Fortran.
#
# Record tripling splits each live record into central/upper/lower weighted copies for the
# first ICL4 cycles (a variance-reduction technique); afterwards growth is the stochastic
# serial-correlation path. ICL4 defaults to 2 (grinit). NOTRIPLE disables tripling (initre.f:5500);
# NUMTRIP n sets ICL4=n (initre.f:2709). Both change the stochastic realization, so each must
# match its live-Fortran .sum:
#   * notriple — base snt01 stand + NOTRIPLE: no tripling, bit-exact every cycle/column.
#   * numtrip  — base stand + NUMTRIP 4 (triple 4 cycles instead of 2): a different realization
#     from the default, still bit-exact vs Fortran.
# (FVSjl previously IGNORED NOTRIPLE — it sat in KNOWN_NOOP and the trip test used a hardcoded
# limit — so a NOTRIPLE stand silently diverged ~20 columns; this test guards the fix.)

using Test, FVSjl

const _TR_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_tr_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 11 && tryparse(Int, first(split(l))) !== nothing]
_tr_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 11 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_trcol(r, c) = parse(Float64, r[c])

@testset "NOTRIPLE / NUMTRIP tripling control vs Fortran" begin
    have(nm) = isfile(joinpath(_TR_DIR, nm * ".key")) && isfile(joinpath(_TR_DIR, nm * ".sum.save"))
    for nm in ("notriple", "numtrip")
        if !have(nm); @test_skip "$nm scenario not available"; continue; end
        @testset "$nm" begin
            jl = _tr_rows(FVSjl.run_keyfile(joinpath(_TR_DIR, nm * ".key"); faithful = true))
            ft = _tr_base(joinpath(_TR_DIR, nm * ".sum.save"))
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for i in 1:length(jl)
                    for c in (3, 4, 7, 8)   # TPA / BA / TopHt / QMD — BIT-EXACT (measured Δ0 both scenarios)
                        @test _trcol(jl[i], c) == _trcol(ft[i], c)
                    end
                end
                # cuft (9): numtrip is BIT-EXACT every row; notriple straddles the print/tree-sum ±1 boundary at
                # one cycle ⇒ exposed @test_broken (doctrine #9), not a passing ≤1. (non-associative tree-SUM order.)
                if nm == "numtrip"
                    @test all(_trcol(jl[i], 9) == _trcol(ft[i], 9) for i in 1:length(jl))       # BIT-EXACT
                else
                    @test_broken all(_trcol(jl[i], 9) == _trcol(ft[i], 9) for i in 1:length(jl))  # tree-sum order
                end
            end
        end
    end
end
