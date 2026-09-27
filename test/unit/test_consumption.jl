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
        before_cwd = copy(s.fire.cwd)
        released = apply_fire_consumption!(s.fire, fuel_moisture(1))
        after = sum(s.fire.cwd)
        @test before > 0f0
        @test after < before                                # fuels were consumed
        # Carbon release splits by pool (fmcrbout.f:151 `V(11)=BIOCON(1)·0.37 + BIOCON(2)·0.50`, with
        # fmdout.f:286-287 BIOCON(1)=consumed litter(10)+duff(11), BIOCON(2)=consumed woody): forest floor
        # converts at 0.37 (Smith & Heath NE-722), all woody classes (1–9) at 0.50 — NOT a uniform 0.5.
        woody_c = sum(before_cwd[1:9,  :, :] .- s.fire.cwd[1:9,  :, :])
        ff_c    = sum(before_cwd[10:11, :, :] .- s.fire.cwd[10:11, :, :])
        # Released = woody·0.5 + forest-floor·0.37 + live herb/shrub burned·0.5 (FMCONS BURNLV: herb 1.0,
        # shrub 0.6; in BIOCON(2), fmdout.f). The live-fuel term closed the fire_carbon released gap (5.13→5.5).
        live_c = s.fire.flive[1] + 0.6f0 * s.fire.flive[2]
        # the pools lose CWD·PRBURN (fmcons.f:309 CWD·(1−PRBURN)), so the pool difference matches to rounding
        @test released ≈ woody_c * 0.5f0 + ff_c * 0.37f0 + live_c * 0.5f0 rtol = 1f-6
        @test ff_c > 0f0                                     # litter/duff WAS consumed (so the split matters)
        @test released < (before - after) * 0.5f0           # the 0.37 forest-floor pulls it below uniform-0.5
        @test released > 0f0
    end

    @testset "FMCONS record: PSBURN scaling, empty 1-3\" class, EXPOSR, smoke" begin
        mois = fuel_moisture(1)
        fs = FireState(); fs.active = true; fs.flive = (0.2f0, 0.5f0)
        fs.cwd[1, 1, 1] = 1f0; fs.cwd[2, 1, 1] = 1f0; fs.cwd[10, 1, 1] = 2f0; fs.cwd[11, 1, 1] = 4f0
        c = FVSjl.fire_consumption!(fs, mois; psburn = 50f0, burncr = 0.3f0)
        # 1-3" class empty ⇒ the <1" classes burn 100% (fmcons.f:135), ×PSBURN/100
        @test c.burned[1] == 1f0 * (1f0 * 50f0 / 100f0) && c.burned[2] == c.burned[1]
        @test fs.cwd[1, 1, 1] == 1f0 * (1f0 - 1f0 * 50f0 / 100f0)
        @test c.burnlv == ((1f0 * 50f0 / 100f0) * 0.2f0, (0.6f0 * 50f0 / 100f0) * 0.5f0)
        prduf = max(0f0, 83.7f0 - 0.426f0 * mois[1, 5] * 100f0)
        @test c.exposr == (-8.98f0 + 0.899f0 * prduf) * 50f0 / 100f0
        # smoke: dead classes × EMMFAC(IM,IL,1,IPM) + live/crown × EMFACL, lb → tons (FMFOUT SMOKE·P2T)
        im = mois[1, 4] <= 0.20f0 ? 3 : mois[1, 4] <= 0.375f0 ? 2 : 1
        ts = sum(c.burned[il] * FVSjl._FM_EMMFAC[im, il, 1] for il in 1:11)
        @test c.smoke[1] ≈ (ts + (c.burnlv[1] + c.burnlv[2] + 0.3f0) * 21.3f0) * 0.0005f0 rtol = 1f-6
        fs2 = FireState(); fs2.cwd[3, 1, 1] = 1f0; fs2.cwd[1, 1, 1] = 1f0
        @test FVSjl.fire_consumption!(fs2, mois).burned[1] == 0.9f0   # 1-3" fuel present, PRBURN(3)=0.65 ⇒ 0.9
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
