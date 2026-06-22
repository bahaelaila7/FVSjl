# test_cuts_coverage.jl — coverage for the CUTS methods/modifiers the decision flow
# flags as unported (see DECISION_FLOW_DETAILED.md §CUTS). The headline snt01.key/
# sn.key exercise SPECPREF/THINPRSC/YARDLOSS/THINSDI only in their thinned stands,
# which the suite never validated. These single-stand scenarios isolate each semantic.
#
# Scenario keys are committed (test/harness/scenarios/cut_*.key); the .tre is a copy
# of snt01.tre (gitignored, so we materialise it here). Golden numbers are inline —
# the values Oracle A (FVSjulia) produces for the affected .sum column at the cut year
# (verified 2026-06-22). Each unported semantic is a `@test_broken`: it documents the
# gap, keeps the suite green, and FLIPS to a failure when the feature starts matching
# the oracle (i.e. when it gets ported). The THINBTA blank-dbhhi fix is a hard `@test`.

using Test, FVSjl

const _HARNESS  = joinpath(@__DIR__, "..", "harness", "scenarios")
const _SNT01_TRE = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.tre"

# Run a single-stand scenario and return its .sum data rows (Float64 vectors).
function _sum_rows(name)
    key = joinpath(_HARNESS, name * ".key")
    tre = joinpath(_HARNESS, name * ".tre")
    isfile(tre) || cp(_SNT01_TRE, tre)            # materialise the (gitignored) .tre
    s, _ = initialize(key)
    FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
    io = IOBuffer(); FVSjl.write_sum_file(io, s)
    rows = Vector{Vector{Float64}}()
    for ln in eachline(IOBuffer(String(take!(io))))
        occursin("-999", ln) && continue
        toks = split(strip(ln)); isempty(toks) && continue
        v = tryparse.(Float64, toks)
        all(!isnothing, v) && push!(rows, Float64.(v))
    end
    return rows
end

# value at the given calendar year + column index (1-based .sum field).
_at(rows, yr, col) = (r = findfirst(x -> x[1] == yr, rows); r === nothing ? NaN : rows[r][col])

@testset "CUTS coverage (decision-flow gap tracker)" begin
    if !isfile(joinpath(_HARNESS, "cut_specpref.key")) || !isfile(_SNT01_TRE)
        @test_skip "cut_* scenarios / snt01.tre not available"
    else
        # cols: 1=yr … 13=rem_tpa … 7=TopHt … 24=MORT-ish (yarding-affected)
        sp = _sum_rows("cut_specpref")
        pr = _sum_rows("cut_thinprsc")
        yl = _sum_rows("cut_yardloss")
        sd = _sum_rows("cut_thinsdi")

        # Regression (FIXED): THINBTA with blank dbhhi must actually thin.
        @testset "THINBTA fires with blank dbhhi (cut_specpref)" begin
            @test _at(sp, 2000.0, 13) > 0          # rem_tpa non-zero at the cut year
        end

        # Unported semantics — SHOULD match the oracle, currently do not. These flip
        # to passing (alerting) when the keyword is ported.
        @testset "SPECPREF reorders the cut (TopHt 2005 == oracle 59)" begin
            @test_broken isapprox(_at(sp, 2005.0, 7), 59.0; atol = 1)
        end
        @testset "THINPRSC removes the prescribed TPA (rem_tpa 2000 == 259)" begin
            @test_broken isapprox(_at(pr, 2000.0, 13), 259.0; atol = 2)
        end
        @testset "YARDLOSS affects removed-volume accounting (col24 2000 == 124)" begin
            @test_broken isapprox(_at(yl, 2000.0, 24), 124.0; atol = 1)
        end
        @testset "THINSDI thins to target SDI (rem_tpa 2000 == 20)" begin
            @test_broken isapprox(_at(sd, 2000.0, 13), 20.0; atol = 2)
        end
    end
end
