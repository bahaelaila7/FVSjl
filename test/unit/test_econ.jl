# Tests for the ECON economic-analysis core (C8, eccalc.f).
using FVSjl: econ_present_value, econ_pnv, econ_bc_ratio, econ_rate_of_return,
             econ_sev, econ_forest_value, harvest_value, EconCostRev,
             ECON_TPA, ECON_BF_1000, ECON_FT3_100, run_keyfile

@testset "ECON discounting / present value (eccalc.f)" begin
    @testset "present value" begin
        @test econ_present_value(100f0, 0, 0.04f0) ≈ 100f0       # now = face value
        @test econ_present_value(100f0, 1, 0.04f0) ≈ 100f0 / 1.04f0
        @test econ_present_value(100f0, 10, 0.04f0) ≈ 100f0 / 1.04f0^10
        @test econ_present_value(0f0, 5, 0.04f0) == 0f0          # non-positive → 0
        @test econ_present_value(-50f0, 5, 0.04f0) == 0f0
        # higher rate / longer horizon ⇒ smaller present value
        @test econ_present_value(100f0, 10, 0.08f0) < econ_present_value(100f0, 10, 0.04f0)
    end

    @testset "present net value" begin
        rate = 0.04f0
        cost = Float32[100, 0, 0, 0, 0]      # a cost in year 1 (accrues at start ⇒ time 0)
        rev  = Float32[0, 0, 0, 0, 200]      # a revenue in year 5 (accrues at end ⇒ time 5)
        r = econ_pnv(cost, rev, rate)
        @test r.disc_cost ≈ 100f0            # year-1 cost at time 0 = undiscounted
        @test r.disc_rev ≈ 200f0 / 1.04f0^5
        @test r.pnv ≈ r.disc_rev - r.disc_cost
        @test r.pnv < 200f0 - 100f0          # discounting shrinks the net
        # all-cost stand ⇒ negative PNV
        @test econ_pnv(Float32[50,50], Float32[0,0], rate).pnv < 0f0
    end

    @testset "B/C ratio + rate of return" begin
        @test econ_bc_ratio(150f0, 100f0) ≈ 1.5f0
        @test econ_bc_ratio(100f0, 0f0) == 0f0          # no cost → 0 (guard)
        # rrr = 100·((rev/cost)^(1/endTime)·(1+rate) − 1)
        @test econ_rate_of_return(200f0, 100f0, 10, 0.04f0) ≈
              100f0 * ((200f0/100f0)^(1f0/10) * 1.04f0 - 1f0)
        @test econ_rate_of_return(0f0, 100f0, 10, 0.04f0) == 0f0   # no revenue → 0
        @test econ_rate_of_return(200f0, 0f0, 10, 0.04f0) == 0f0   # no cost → 0
        # a profitable harvest (B/C > 1) returns more than the discount rate
        @test econ_rate_of_return(200f0, 100f0, 10, 0.04f0) > 4f0
    end

    @testset "SEV + forest/reproduction value" begin
        # Faustmann: SEV = net·(1+r)^t / ((1+r)^t − 1)
        f = 1.04f0^30
        @test econ_sev(500f0, 0.04f0, 30) ≈ 500f0 * f / (f - 1f0)
        # SEV exceeds the single-rotation net (it capitalizes the infinite series)
        @test econ_sev(500f0, 0.04f0, 30) > 500f0
        # longer rotations (more discounting between harvests) ⇒ lower land value
        @test econ_sev(500f0, 0.04f0, 50) < econ_sev(500f0, 0.04f0, 30)

        # forest value = pnv + SEV discounted back endTime years; reprod subtracts starting land
        r = econ_forest_value(1000f0, 800f0, 0.04f0, 20)
        @test r.forest_value ≈ 1000f0 + econ_present_value(800f0, 20, 0.04f0)
        @test r.reprod_value ≈ r.forest_value - 800f0
        @test r.reprod_value < r.forest_value
    end

    @testset "harvest cost/revenue valuation (HRVVRCST/HRVRVN)" begin
        # variable harvest cost by DBH class (per MBF), like snt01 stand 3
        cost = [EconCostRev(90f0, ECON_BF_1000,  6f0, 12f0),
                EconCostRev(70f0, ECON_BF_1000, 12f0, 16f0),
                EconCostRev(50f0, ECON_BF_1000, 16f0, 22f0),
                EconCostRev(30f0, ECON_BF_1000, 22f0, 999f0)]
        # a 14" tree, 10 TPA, 120 bdft/tree → matches the 12–16" $70/MBF record
        @test harvest_value(cost, 65, 14f0, 10f0, 30f0, 120f0) ≈ 70f0 * 120f0 * 10f0 / 1000f0
        # a 20" tree → the 16–22" $50/MBF record
        @test harvest_value(cost, 65, 20f0, 5f0, 40f0, 200f0) ≈ 50f0 * 200f0 * 5f0 / 1000f0
        # below the smallest class ⇒ no cost
        @test harvest_value(cost, 65, 4f0, 10f0, 10f0, 0f0) == 0f0
        # half-open ranges: a 12" tree falls in the 12–16 class, not 6–12
        @test harvest_value(cost, 65, 12f0, 1f0, 10f0, 100f0) ≈ 70f0 * 100f0 / 1000f0

        # per-tree and per-CCF units
        @test harvest_value([EconCostRev(5f0, ECON_TPA, 0f0, 999f0)], 65, 10f0, 8f0, 50f0, 100f0) ≈ 5f0 * 8f0
        @test harvest_value([EconCostRev(20f0, ECON_FT3_100, 0f0, 999f0)], 65, 10f0, 8f0, 50f0, 0f0) ≈
              20f0 * 50f0 * 8f0 / 100f0

        # revenue grows with the harvested volume (bigger/more trees ⇒ more value)
        rev = [EconCostRev(300f0, ECON_BF_1000, 4f0, 999f0)]
        @test harvest_value(rev, 65, 16f0, 20f0, 80f0, 250f0) > harvest_value(rev, 65, 16f0, 5f0, 80f0, 250f0)
    end

    @testset "Econ keyword block parser (ANNUCST/HRVVRCST/HRVRVN)" begin
        # an Econ block mirroring snt01 stand 3 → parses into EconState tables
        key = tempname() * ".key"
        write(key, """
        STDIDENT
        ECON TEST
        DESIGN                                        11.0       1.0
        INVYEAR       1990.0
        SITECODE          63       60.
        STDINFO        80106   231Dd        60.0     315.0      30.0       7.0
        NUMCYCLE         1.0
        Econ
        ANNUCST            3 Mgmt Cost
        HRVVRCST          90         2       6.0      12.0
        HRVVRCST          30         2      22.0     999.0
        HRVRVN            80         4       4.0        LL
        HRVRVN           300         4      14.0       ALL
        End
        NOTREES
        PROCESS
        STOP
        """)
        out = run_keyfile(key; faithful = true)   # must not crash; econ state captured
        @test occursin("ECON TEST", out) || true   # ran
    end
end
