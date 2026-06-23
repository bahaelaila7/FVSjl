# Tests for the standing live-tree carbon pools (FFE F8, fmcrbout.f).
# Carbon = Jenkins biomass × TPA × 0.5, summed over the tree list.
using FVSjl: stand_live_carbon, jenkins_biomass, StandState, init_blockdata!,
             init_merch_standards!, coefficients, Southern

@testset "standing live-tree carbon (FMCRBOUT)" begin
    s = StandState(Southern()); init_blockdata!(s, s.variant); init_merch_standards!(s)
    coef = coefficients(Southern())
    t = s.trees; t.n = 3
    t.species[1]=Int32(65); t.dbh[1]=14f0; t.tpa[1]=40f0
    t.species[2]=Int32(33); t.dbh[2]= 8f0; t.tpa[2]=30f0
    t.species[3]=Int32(22); t.dbh[3]= 4f0; t.tpa[3]=20f0

    c = stand_live_carbon(s)

    # equals the hand-summed Jenkins biomass × TPA × 0.5
    above = 0f0; merch = 0f0; root = 0f0
    for i in 1:3
        a, m, r = jenkins_biomass(coef, t.species[i], t.dbh[i])
        above += a*t.tpa[i]; merch += m*t.tpa[i]; root += r*t.tpa[i]
    end
    @test c.aboveground ≈ above * 0.5f0
    @test c.merch ≈ merch * 0.5f0
    @test c.belowground ≈ root * 0.5f0

    # physical sanity: positive pools, aboveground > belowground > 0, merch ≤ aboveground
    @test c.aboveground > 0f0 && c.belowground > 0f0
    @test c.aboveground > c.belowground
    @test c.merch <= c.aboveground

    # scales linearly with TPA
    for i in 1:3; t.tpa[i] *= 2f0; end
    c2 = stand_live_carbon(s)
    @test c2.aboveground ≈ 2f0 * c.aboveground

    # zero-TPA trees contribute nothing
    for i in 1:3; t.tpa[i] = 0f0; end
    c0 = stand_live_carbon(s)
    @test c0.aboveground == 0f0 && c0.belowground == 0f0
end
