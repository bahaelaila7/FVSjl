# PPE MXHRVP — DIFFERENTIAL end-to-end: the SELECTED-changes-the-yield path + the partial-cut
# (status-4/HVPART) branch, vs the FVSppe oracle golden (test/fixtures/ppe/mxhrvp/msp_diff.golden.txt).
#
# msp_diff.key: CREDIT = SELECTED * BBA (so HVYLDS[selected]=BBA ≠ HVTHIN[not]=0 — a genuine
# differential, unlike the degenerate msp.key), TARGET=150. HVSEL then selects a SUBSET: stand 1
# FULL (status 3), stand 2 PARTIAL (status 4, HVPART), stand 3 NO (status 2) — SELECTED RESOURCE
# fills to exactly the target every master cycle. Oracle partials: 1990 93.8% · 2000 36.2% · 2010 2.3%.
#
# This exercises "the test must exercise the semantic": HVYLDS≠HVTHIN driving the greedy subset
# selection + the partial-cut fraction — the paths the degenerate all-selected msp.key never hit.
# CORNERING: cyc0 BIT-EXACT (BBA=77.39206 ⇒ HVPART=0.938 == oracle 93.8%); cyc1+ HVPART is a
# function of BBA, so it inherits the KNOWN EC-variant before-thin-BA growth straddle (37.8% vs
# 36.2%, 4.0% vs 2.3%) — a cornered residual, NOT a MXHRVP defect. The SELECT PATTERN + the
# fill-to-target (SELECTED RESOURCE=150) are bit-exact every cycle.
using Test
using FVSjl: ppe_run_landscape_harvest!, PPEStand, EastCascades

@testset "PPE MXHRVP — DIFFERENTIAL selection + partial cut vs FVSppe oracle" begin
    kf  = joinpath(@__DIR__, "..", "fixtures", "ppe", "mxhrvp", "stand.key")
    st  = PPEStand(kf; area = 11.0)
    res = ppe_run_landscape_harvest!([st, st, st]; variant = EastCascades(),
              labels = ["ALL", "ALL", "ALL"], mslabel = "ALL",
              target_expr = "150", priority_expr = "BBA", credit_expr = "SELECTED * BBA",
              master_years = [1990, 2000, 2010], lprtct = true)   # EXACT (partial cut) is the oracle default
    @test length(res) == 3
    r1, r2, r3 = res

    # --- the DIFFERENTIAL: HVYLDS(selected) ≠ HVTHIN(not selected) every cycle ---
    for r in res
        @test all(r.credit_notsel .== 0f0)          # HVTHIN = SELECTED(0)*BBA = 0
        @test all(r.credit_sel .> 0f0)              # HVYLDS = SELECTED(1)*BBA = BBA > 0
        @test any(r.credit_sel .!= r.credit_notsel) # genuinely differential (not the degenerate case)
    end

    # --- SELECT PATTERN (tie-break-invariant): exactly one FULL(3) + one PARTIAL(4) + one NO(2) ---
    #     The 3 stands are physically IDENTICAL (same keyfile ⇒ exactly equal BBA), so WHICH stand
    #     gets full/partial/not is an RDPSRT unstable-sort tie-break on equal priorities (#206-class,
    #     a named primitive) — the aggregate (2 selected filling to target) is identical either way.
    #     FVSjl's tie-break puts the partial on a different physical stand than the oracle; harmless.
    for r in res
        @test sort(r.status) == [2, 3, 4]                         # one of each, regardless of assignment
        @test count(r.selected) == 2                              # two selected (one full + one partial)
        @test isapprox(r.selected_resource, 150f0; atol = 1f-3)   # fills to target (oracle 0.1500000E+03)
        @test isapprox(r.pct_of_target, 100f0; atol = 0.05f0)     # oracle PERCENT OF TARGET = 100.0
    end

    # --- cyc0 (1990): the partial fraction is BIT-EXACT (BBA=77.39206 ⇒ HVPART=0.93818 == 93.8%) ---
    @test isapprox(r1.hvpart, 0.93818f0; atol = 0.0005f0)         # oracle **YES= 93.8 PERCENT
    @test isapprox(r1.credit_sel[1], 77.39f0; atol = 0.005f0)     # HVYLDS = BBA = oracle 0.7739E+02

    # --- cyc1/cyc2: HVPART cornered to the EC before-thin-BA growth straddle (documented) ---
    #     FVSjl BBA 108.87 ⇒ HVPART 37.8% (oracle 36.2%); BBA 144.29 ⇒ 4.0% (oracle 2.3%).
    @test isapprox(r2.hvpart, 0.378f0; atol = 0.002f0)            # FVSjl EC growth (oracle 0.362)
    @test isapprox(r3.hvpart, 0.040f0; atol = 0.003f0)            # FVSjl EC growth (oracle 0.023)
    @test abs(r2.hvpart - 0.362f0) < 0.02                         # within the EC-growth corner vs oracle
end
