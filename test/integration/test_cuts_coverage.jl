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
        ht = _sum_rows("cut_thinht")
        cc = _sum_rows("cut_thincc")
        rd = _sum_rows("cut_thinrden")
        au = _sum_rows("cut_thinauto")
        qf = _sum_rows("cut_thinqfa")

        # PORTED + validated: THINBTA blank-dbhhi fires, and SPECPREF reorders the cut
        # so the year-2000 removal is bit-exact to Oracle A (rem_tpa 327, at_BA 96).
        @testset "SPECPREF cut is bit-exact at the cut year (cut_specpref)" begin
            @test _at(sp, 2000.0, 13) > 0                          # THINBTA fires
            @test isapprox(_at(sp, 2000.0, 13), 327.0; atol = 1)   # rem_tpa
            @test isapprox(_at(sp, 2000.0, 18),  96.0; atol = 1)   # at-treatment BA
        end
        # STILL OPEN (downstream of the cut): post-thin DGSCOR growth — the 2005 top
        # height / sawtimber volume of a thinned stand. Shared by s29 + snt01 stands 3-4.
        @testset "post-thin growth matches (cut_specpref TopHt 2005 == 59)" begin
            @test_broken isapprox(_at(sp, 2005.0, 7), 59.0; atol = 1)
        end
        @testset "THINPRSC removes the prescribed TPA (rem_tpa 2000 == 259)" begin
            @test isapprox(_at(pr, 2000.0, 13), 259.0; atol = 2)   # ported (cut-code marked)
        end
        @testset "YARDLOSS affects removed-volume accounting (col24 2000 == 124)" begin
            @test_broken isapprox(_at(yl, 2000.0, 24), 124.0; atol = 1)
        end
        @testset "THINSDI thins to target SDI (rem_tpa 2000 == 20)" begin
            @test isapprox(_at(sd, 2000.0, 13), 20.0; atol = 2)   # ported (Zeide SDI + proportional)
        end
        # ── remaining label_325/400 thins: live-Fortran ground truth (TPA @2005). Each
        #    is ⛔ until ported; the expected value is the Fortran result so the test
        #    validates the port (flip @test_broken→@test when it lands). NOT tuned to FVSjl.
        @testset "THINHT thins a height class (TPA 2005 == 133)" begin
            @test isapprox(_at(ht, 2005.0, 3), 133.0; atol = 2)   # ported (label_325 on height)
        end
        @testset "THINCC thins to residual crown cover (TPA 2005 == 71)" begin
            @test isapprox(_at(cc, 2005.0, 3), 71.0; atol = 2)         # ported (CCCLS crown cover)
        end
        @testset "THINRDEN thins to Curtis RD (TPA 2005 == 418)" begin
            @test isapprox(_at(rd, 2005.0, 3), 418.0; atol = 2)        # ported (Curtis RD + proportional)
        end
        @testset "THINAUTO auto-thins on stocking (TPA 2005 == 231)" begin
            @test isapprox(_at(au, 2005.0, 3), 231.0; atol = 2)        # ported (AUTSTK + recurring from-below)
        end
        @testset "THINQFA thins to a Q-factor distribution (TPA 2005 == 89)" begin
            @test isapprox(_at(qf, 2005.0, 3), 89.0; atol = 2)         # ported (CUTQFA Q-factor + 2-record)
        end
    end
end
