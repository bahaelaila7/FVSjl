# Tests for fire fuel consumption + carbon release (FFE F7/F8, FMCONS).
using FVSjl: fire_consumption_fractions, apply_fire_consumption!, fuel_moisture,
             FireState, StandState, init_blockdata!, init_merch_standards!, fmcba!, fmburn!, Southern

@testset "fire fuel consumption (FMCONS)" begin
    @testset "consumption fractions" begin
        dry = fire_consumption_fractions(fuel_moisture(1))   # very dry
        @test length(dry) == 11
        @test dry[1] == 0.9f0 && dry[2] == 0.9f0             # <1" classes
        @test dry[3] == 0.65f0                                # 1–3"
        @test dry[10] == 1.0f0                                # litter fully consumed
        @test all(0f0 .<= dry .<= 1f0)
        # >3" classes consume a positive fraction when dry
        @test all(dry[i] > 0f0 for i in 4:9)
        # duff: drier ⇒ more consumed than wet
        wet = fire_consumption_fractions(fuel_moisture(4))   # very wet
        @test dry[11] > wet[11]
        # the large-fuel diameter-reduction also burns less when wetter
        @test dry[4] >= wet[4]
    end

    @testset "apply consumption + carbon release" begin
        s = StandState(Southern()); init_blockdata!(s, s.variant); init_merch_standards!(s)
        s.plot.forest_type = Int32(520)
        s.plot.latitude = 35f0; s.plot.longitude = -80f0; s.plot.elevation = 10f0
        s.fire = FireState(); s.fire.active = true
        s.trees.n = 0
        fmcba!(s)                                            # populate the down-wood pools
        before = sum(s.fire.cwd)
        released = apply_fire_consumption!(s.fire, fuel_moisture(1))
        after = sum(s.fire.cwd)
        @test before > 0f0
        @test after < before                                # fuels were consumed
        @test released ≈ (before - after) * 0.5f0           # carbon = consumed biomass × 0.5
        @test released > 0f0
    end

    @testset "fmburn! reports carbon release" begin
        s = StandState(Southern()); init_blockdata!(s, s.variant); init_merch_standards!(s)
        s.plot.forest_type = Int32(520)
        s.plot.latitude = 35f0; s.plot.longitude = -80f0; s.plot.elevation = 10f0
        t = s.trees; t.n = 2
        t.species[1]=Int32(65); t.dbh[1]=14f0; t.height[1]=72f0; t.tpa[1]=30f0; t.crown_pct[1]=Int32(40)
        t.species[2]=Int32(33); t.dbh[2]= 8f0; t.height[2]=55f0; t.tpa[2]=30f0; t.crown_pct[2]=Int32(45)
        s.fire = FireState(); s.fire.active = true
        res = fmburn!(s; wind = 20f0, fmois = 1, year = 2003)
        @test res.carbon_released > 0f0                      # the fire released carbon from fuel
    end
end
