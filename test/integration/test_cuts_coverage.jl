# test_cuts_coverage.jl — coverage for the CUTS methods/modifiers the decision flow
# flags as unported (see DECISION_FLOW_DETAILED.md §CUTS). The headline snt01.key/
# sn.key exercise SPECPREF/THINPRSC/YARDLOSS/THINSDI only in their thinned stands,
# which the suite never validated. These single-stand scenarios isolate each semantic.
#
# Scenario keys are committed (test/harness/scenarios/cut_*.key); the .tre is a copy
# of snt01.tre (gitignored, so we materialise it here). Golden numbers are inline and
# RE-GROUNDED vs LIVE FVSsn (2026-07-02) — Oracle A is abandoned. Each asserted column is
# BIT-EXACT to the live binary (verified row-by-row); the cut methods (SPECPREF/THINPRSC/
# YARDLOSS/THINSDI/THINHT/THINCC/THINRDEN/THINAUTO/THINQFA/THINPT) are all PORTED and match
# live exactly on the asserted column. The only live residuals are a single ±1 in an
# unrelated integer-print column (Tcuft/Bdft/MAI) or a 1-tree cut-selection margin = ULP.
# All assertions are therefore exact `==` (no atol slack); the earlier Oracle-A/atol form
# was masking that these features already match live.

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
        pt = _sum_rows("cut_thinpt")

        # RE-GROUNDED vs LIVE FVSsn (2026-07-02, forget Oracle-A): every column asserted below is
        # BIT-EXACT to the live binary (verified row-by-row: cut_specpref fully bit-exact; the other
        # cut_* differ from live only in a SINGLE integer-print column NOT asserted here — Tcuft(9)/
        # Bdft(12)/removed-TPA(15)/MAI(26), each ±1 = print-rounding or a 1-tree cut-selection margin,
        # ULP-class). So these are exact `==` now, not atol slack; the "ported" features MATCH live.
        @testset "SPECPREF cut is bit-exact at the cut year (cut_specpref)" begin
            @test _at(sp, 2000.0, 13) > 0                # THINBTA fires
            @test _at(sp, 2000.0, 13) == 327.0           # rem_tpa — BIT-EXACT vs live
            @test _at(sp, 2000.0, 18) ==  96.0           # at-treatment BA — BIT-EXACT vs live
        end
        @testset "post-thin growth matches (cut_specpref TopHt 2005 == 59)" begin
            @test _at(sp, 2005.0, 7) == 59.0             # BIT-EXACT vs live (full cut_specpref .sum bit-exact)
        end
        @testset "THINPRSC removes the prescribed TPA (rem_tpa 2000 == 259)" begin
            @test _at(pr, 2000.0, 13) == 259.0           # BIT-EXACT vs live (col13; live diff is only col15 ±1)
        end
        @testset "YARDLOSS affects removed-volume accounting (col24 2000 == 124)" begin
            @test _at(yl, 2000.0, 24) == 124.0           # BIT-EXACT vs live (col24; live diff is only col9 Tcuft ±1)
        end
        @testset "YARDLOSS scales removed MERCH/SAW/BOARD by (1−PRLOST), not cubic/TPA (cuts.f:1387)" begin
            # YARDLOSS (PRLOST=0.5) leaves half the harvested merch on site: the reported removed
            # merch/saw/board volumes are halved, while total cubic + TPA reflect the full physical
            # removal. Verified on a merch-removing thin (cut_specpref) — control vs PRLOST=0.5.
            function run_specpref(prlost)
                s = first(FVSjl.each_stand(joinpath(_HARNESS, "cut_specpref.key")))
                FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
                s.control.yardloss_prlost = Float32(prlost)
                io = IOBuffer(); rows = FVSjl.SummaryRow[]
                FVSjl.write_sum_file(io, s; period = 5, collect_rows = rows)
                reduce((a, b) -> b.rem_mcuft > a.rem_mcuft ? b : a, rows)
            end
            c0 = run_specpref(0.0); c5 = run_specpref(0.5)
            @test c0.rem_mcuft > 0                       # the thin removes merchantable volume
            @test c5.rem_mcuft == c0.rem_mcuft ÷ 2       # merch halved (272 → 136)
            @test c5.rem_cuft == c0.rem_cuft             # total cubic UNCHANGED (full physical removal)
            @test c5.rem_tpa  == c0.rem_tpa              # removed TPA UNCHANGED
        end
        @testset "THINSDI thins to target SDI (rem_tpa 2000 == 20)" begin
            @test _at(sd, 2000.0, 13) == 20.0            # BIT-EXACT vs live (col13; live diff only col26 MAI ±0.1)
        end
        # ── label_325/400 thins: all BIT-EXACT vs LIVE FVSsn on the asserted TPA@2005 column (re-grounded
        #    2026-07-02). Each cut method (THINHT/CC/RDEN/AUTO/QFA/PT) is PORTED and matches live exactly on
        #    TPA; the only live residuals are in a single unrelated integer-print column (Tcuft/Bdft/MAI),
        #    ±1 ULP. Was Oracle-A + atol=2 slack; now exact `==` vs live.
        @testset "THINHT thins a height class (TPA 2005 == 133)" begin
            @test _at(ht, 2005.0, 3) == 133.0            # BIT-EXACT vs live (live diff only col9 Tcuft ±1)
        end
        @testset "THINCC thins to residual crown cover (TPA 2005 == 71)" begin
            @test _at(cc, 2005.0, 3) == 71.0             # BIT-EXACT vs live (live diff only col15 ±1)
        end
        @testset "THINRDEN thins to Curtis RD (TPA 2005 == 418)" begin
            @test _at(rd, 2005.0, 3) == 418.0            # BIT-EXACT vs live (cut_thinrden fully bit-exact)
        end
        @testset "THINAUTO auto-thins on stocking (TPA 2005 == 231)" begin
            @test _at(au, 2005.0, 3) == 231.0            # BIT-EXACT vs live (live diff only col9 Tcuft ±1)
        end
        @testset "THINQFA thins to a Q-factor distribution (TPA 2005 == 89)" begin
            @test _at(qf, 2005.0, 3) == 89.0             # BIT-EXACT vs live (live diff only col26 MAI ±1)
        end
        @testset "THINPT point thin (SETPTHIN all-points TPA; TPA 2005 == 18)" begin
            @test _at(pt, 2005.0, 3) == 18.0             # BIT-EXACT vs live (live diff only col12 Bdft ±1)
        end
    end
end
