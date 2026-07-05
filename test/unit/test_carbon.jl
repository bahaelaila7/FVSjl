# Tests for the standing live-tree carbon pools (FFE F8, fmcrbout.f).
# Carbon = Jenkins biomass × TPA × 0.5, summed over the tree list.
using FVSjl: stand_live_carbon, standing_dead_carbon, down_wood_carbon, forest_floor_carbon,
             shrub_herb_carbon, stand_carbon,
             jenkins_biomass, StandState, init_blockdata!, init_merch_standards!,
             coefficients, Southern, FireState, add_snag!, fmcba!

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
    @test c.aboveground == above * 0.5f0
    @test c.merch == merch * 0.5f0
    @test c.belowground == root * 0.5f0

    # physical sanity: positive pools, aboveground > belowground > 0, merch ≤ aboveground
    @test c.aboveground > 0f0 && c.belowground > 0f0
    @test c.aboveground > c.belowground
    @test c.merch <= c.aboveground

    # scales linearly with TPA
    for i in 1:3; t.tpa[i] *= 2f0; end
    c2 = stand_live_carbon(s)
    @test c2.aboveground == 2f0 * c.aboveground

    # zero-TPA trees contribute nothing
    for i in 1:3; t.tpa[i] = 0f0; end
    c0 = stand_live_carbon(s)
    @test c0.aboveground == 0f0 && c0.belowground == 0f0
end

@testset "dead + total stand carbon pools (FMCRBOUT)" begin
    s = StandState(Southern()); init_blockdata!(s, s.variant); init_merch_standards!(s)
    s.plot.forest_type = Int32(520)
    s.plot.latitude = 35f0; s.plot.longitude = -80f0; s.plot.elevation = 10f0
    coef = coefficients(Southern())
    t = s.trees; t.n = 2
    t.species[1]=Int32(65); t.dbh[1]=14f0; t.height[1]=72f0; t.tpa[1]=40f0; t.crown_pct[1]=Int32(40)
    t.species[2]=Int32(33); t.dbh[2]= 8f0; t.height[2]=55f0; t.tpa[2]=30f0; t.crown_pct[2]=Int32(45)

    # no FFE → all dead pools zero
    @test standing_dead_carbon(s) == 0f0 && down_wood_carbon(s) == 0f0

    s.fire = FireState(); s.fire.active = true
    fmcba!(s)                                         # populate the cwd / flive surface-fuel pools
    @test down_wood_carbon(s) > 0f0                   # down-wood (woody size classes 1–9) carbon
    # fmcrbout.f carbon fractions: down wood (classes 1–9) at 0.5, forest floor (litter+duff,
    # classes 10–11) at 0.37 (Smith & Heath), shrub/herb (FLIVE) at 0.5 — NOT a flat 0.5.
    @test down_wood_carbon(s)    == sum(@view s.fire.cwd[1:9, :, :]) * 0.5f0
    @test forest_floor_carbon(s) == sum(@view s.fire.cwd[10:11, :, :]) * 0.37f0
    @test shrub_herb_carbon(s)   == (s.fire.flive[1] + s.fire.flive[2]) * 0.5f0

    # snag carbon = Jenkins aboveground × standing density × 0.5
    add_snag!(s.fire, 65, 14f0, 40f0, 2003)
    a, _, _ = jenkins_biomass(coef, 65, 14f0)
    @test standing_dead_carbon(s) == a * 40f0 * 0.5f0

    # the combined report: total = live(above+below) + snag + down wood + forest floor + shrub/herb
    c = stand_carbon(s)
    @test c.standing_dead == standing_dead_carbon(s)
    @test c.down_wood == down_wood_carbon(s)
    @test c.forest_floor == forest_floor_carbon(s) && c.shrub_herb == shrub_herb_carbon(s)
    @test c.total == c.live_above + c.live_below + c.standing_dead + c.down_wood + c.forest_floor + c.shrub_herb
    @test c.total > 0f0
end
