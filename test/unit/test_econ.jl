# Tests for the ECON economic-analysis core (C8, eccalc.f).
using FVSjl: econ_present_value, econ_pnv, econ_bc_ratio, econ_rate_of_return,
             econ_sev, econ_forest_value

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
end
