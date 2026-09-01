# PPE MXHRVP — END-TO-END: the ppe_run_landscape_harvest! coordinator (hvaloc.f) vs the
# proven FVSppe oracle golden (test/fixtures/ppe/mxhrvp/msp.golden.txt).
#
# The oracle golden (msp.key: 3 EC S248112 stands, one policy MSPLABEL=ALL / TARGET=1000 /
# PRIORITY=BBA / CREDIT=BBA) has HVALOC/HVSEL select ALL 3 stands every master cycle (target
# 1000 ≫ available), with per-stand PRIORITY=CREDIT=BBA and SELECTED RESOURCE = Σ BBA:
#   1990: BBA 77.39, SELECTED 232.18 (23.2%)   2000: 110.2, 330.48 (33.0)   2010: 146.6, 439.90 (44.0)
#
# This reconstructs that landscape from FVSjl standalone stands and runs the coordinator.
# CORNERING: FVSjl's before-thin BA (=BBA) is BIT-EXACT at cyc0 (77.39206 == oracle 0.7739E+02)
# but under-grows beyond cyc0 (108.87 vs 110.2, 144.29 vs 146.6) — the KNOWN EC-variant growth
# residual (#207 EC growth is cyc0-bit-exact; beyond that a cornered straddle), INHERITED by
# CREDIT=BBA, NOT a MXHRVP defect. So the SELECTION LOGIC is validated bit-exact (every SELECT
# outcome + the cyc0 numerics), and the cyc1+ numeric residual is cornered to EC growth.
using Test
using FVSjl: ppe_run_landscape_harvest!, PPEStand, EastCascades

@testset "PPE MXHRVP — end-to-end coordinator vs FVSppe oracle (msp landscape)" begin
    kf  = joinpath(@__DIR__, "..", "fixtures", "ppe", "mxhrvp", "stand.key")
    st  = PPEStand(kf; area = 11.0)                       # SAMPLE WEIGHT 33 = 3 × 11
    res = ppe_run_landscape_harvest!([st, st, st]; variant = EastCascades(),
              labels = ["ALL", "ALL", "ALL"], mslabel = "ALL",
              target_expr = "1000", priority_expr = "BBA", credit_expr = "BBA",
              master_years = [1990, 2000, 2010])

    @test length(res) == 3
    r1, r2, r3 = res

    # --- SELECTION OUTCOMES: bit-exact vs oracle (all 3 stands YES every master cycle) ---
    @test all(r.selected == [true, true, true] for r in res)
    @test [r.year for r in res] == [1990, 2000, 2010]
    @test all(r.target == 1000f0 for r in res)

    # --- cyc0 (1990): BIT-EXACT numerics vs oracle (FVSjl BBA=77.39206 rounds to 0.7739E+02) ---
    @test all(isapprox(c, 77.39f0; atol = 0.005f0) for c in r1.credit)      # BBA = oracle 0.7739E+02
    @test all(isapprox(p, 77.39f0; atol = 0.005f0) for p in r1.priority)    # PRIORITY = CREDIT = BBA
    @test isapprox(r1.selected_resource, 232.176f0; atol = 0.01f0)          # oracle 0.2321762E+03
    @test isapprox(r1.pct_of_target, 23.2f0; atol = 0.05f0)                 # oracle PERCENT OF TARGET
    @test r1.nonselected_supply < 1f-3                                      # ~1.5e-5 running-subtraction residual

    # --- cyc1/cyc2: SELECT bit-exact; CREDIT cornered to the EC growth residual (documented) ---
    #     FVSjl 2000: BBA 108.87 (oracle 110.2); 2010: 144.29 (oracle 146.6). Under-grows ≤1.6%.
    @test round(r2.selected_resource, digits = 1) ≈ 326.6 atol = 0.2  # = 3 × 108.87 (FVSjl EC growth)
    @test round(r3.selected_resource, digits = 1) ≈ 432.9 atol = 0.2  # = 3 × 144.29
    @test abs(r2.selected_resource - 330.48) / 330.48 < 0.02          # within the EC-growth corner vs oracle
    @test abs(r3.selected_resource - 439.90) / 439.90 < 0.02
end
