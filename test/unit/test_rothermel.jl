# Tests for the Rothermel surface-fire model (FFE F5, FMFINT).
#
# Like the other inert FFE pieces, exact bit-for-bit agreement can only be confirmed
# against live Fortran once the fuel-model construction (F5b) + burn integration are in
# (FMFINT has no standalone .sum output). These tests pin the physical invariants, the
# Byram→flame relationship, and the too-moist cutoff — the behaviors a correct Rothermel
# implementation must satisfy — plus determinism.
using FVSjl: rothermel_surface_fire, fuel_moisture, fire_wind_reduction,
             scorch_height, crown_volume_scorched, fire_tree_mortality, coefficients, Southern,
             build_dynamic_fuel_model, StandState, init_blockdata!, init_merch_standards!, FireState, fmcba!,
             standard_fuel_model

@testset "Rothermel surface fire (FMFINT)" begin
    # a single dead 1-hr fuel (grass-like): load, SAV, depth, Mx, moisture
    function grass(; wind = 0f0, slope = 0f0, m = 0.06f0)
        load = zeros(Float32, 2, 4); sav = zeros(Float32, 2, 4); mois = zeros(Float32, 2, 4)
        load[1, 1] = 0.034f0; sav[1, 1] = 3500f0; mois[1, 1] = m
        rothermel_surface_fire(load, sav, 1.0f0, 0.12f0, mois; wind = wind, slope_tan = slope)
    end

    r0 = grass(); r5 = grass(wind = 5f0); r10 = grass(wind = 10f0)

    @testset "physical invariants" begin
        # a dry grass fire spreads and has positive intensity
        @test r0.spread > 0f0 && r0.byram > 0f0 && r0.flame > 0f0
        @test r0.sigma ≈ 3500f0                       # characteristic SAV = the single class's SAV
        # wind strictly increases spread / intensity
        @test r10.spread > r5.spread > r0.spread
        @test r10.byram  > r5.byram  > r0.byram
        # slope increases spread
        @test grass(slope = 0.5f0).spread > r0.spread
        # wetter fuel spreads slower
        @test grass(m = 0.10f0).spread < grass(m = 0.04f0).spread
        # fuel wetter than its moisture of extinction (0.12) cannot carry fire
        rwet = grass(m = 0.30f0)
        @test rwet.byram == 0f0 && rwet.flame == 0f0 && rwet.spread == 0f0
    end

    @testset "Byram → flame relationship" begin
        # flame length is Byram's 0.45·(I/60)^0.46
        for r in (r0, r5, r10)
            @test r.flame ≈ 0.45f0 * (r.byram / 60f0)^0.46f0
        end
    end

    @testset "fuel structure effects" begin
        # a coarser timber-litter fuel (lower SAV, deeper bed) spreads slower than fine grass
        load = zeros(Float32, 2, 4); sav = zeros(Float32, 2, 4); mois = zeros(Float32, 2, 4)
        load[1, 1] = 0.10f0; sav[1, 1] = 2000f0
        load[1, 2] = 0.10f0; sav[1, 2] = 109f0       # 10-hr
        load[1, 3] = 0.15f0; sav[1, 3] = 30f0        # 100-hr (dropped: SAV<16? no, 30≥16 kept)
        mois[1, 1] = 0.06f0; mois[1, 2] = 0.07f0; mois[1, 3] = 0.08f0
        timber = rothermel_surface_fire(load, sav, 0.3f0, 0.25f0, mois; wind = 5f0)
        @test timber.spread > 0f0
        @test timber.spread < r5.spread               # grass spreads faster than timber litter
        # determinism
        @test timber.byram == rothermel_surface_fire(load, sav, 0.3f0, 0.25f0, mois; wind = 5f0).byram
    end

    @testset "empty fuel" begin
        z = zeros(Float32, 2, 4)
        r = rothermel_surface_fire(z, z, 1f0, 0.12f0, z; wind = 5f0)
        @test r.byram == 0f0 && r.flame == 0f0
    end

    @testset "standard fuel models (Anderson 13)" begin
        coef = coefficients(Southern())
        # database values (fminit.f) for the models snt01 stand 4 uses
        l10, s10, d10, m10 = standard_fuel_model(coef, 10)   # timber litter + understory
        @test l10[1, 1] ≈ 0.13820f0 && l10[1, 2] ≈ 0.09183f0 && l10[1, 3] ≈ 0.23003f0
        @test l10[2, 1] ≈ 0.09183f0                          # live woody
        @test s10[1, 1] == 2000f0 && s10[1, 2] == 109f0 && s10[1, 3] == 30f0
        @test d10 == 1.0f0 && m10 ≈ 0.25f0
        l5, _, d5, m5 = standard_fuel_model(coef, 5)         # brush
        @test l5[1, 1] ≈ 0.04591f0 && d5 == 2.0f0 && m5 ≈ 0.20f0

        # Under the very-dry FMMOIS model 1, the brush model 5 (live woody, drying out) carries a HOTTER
        # fire than timber-litter model 10 — verified vs live FVSsn FMFINT at fmois=1: FM5 byram 8988 >
        # FM10 6519 (FWIND=1); at wind 2, FVSjl FM5 10021 > FM10 6608. Both carry fire.
        mois = fuel_moisture(1)
        r10 = rothermel_surface_fire(l10, s10, d10, m10, mois; wind = 2f0)
        r5  = rothermel_surface_fire(l5, standard_fuel_model(coef, 5)[2], d5, m5, mois; wind = 2f0)
        @test r10.byram > 0f0 && r10.flame > 0f0
        @test r5.byram > r10.byram > 0f0

        # the weighted 10(96%)+5(4%) blend under a canopy-reduced wind reproduces the
        # Fortran PotFire flame magnitude (~3–5 ft), not the dynamic model's low ~2 ft.
        fwind = 20f0 * fire_wind_reduction(87.6f0)
        flame = 0f0
        for (fm, wt) in ((10, 0.96f0), (5, 0.04f0))
            ld, sv, dp, mx = standard_fuel_model(coef, fm)
            flame += rothermel_surface_fire(ld, sv, dp, mx, mois; wind = fwind).flame * wt
        end
        @test 3.5f0 < flame < 7f0
    end

    @testset "fuel moisture scenario (FMMOIS)" begin
        vd = fuel_moisture(1); vw = fuel_moisture(4)
        @test vd[1, 1] ≈ 0.05f0 && vd[1, 5] ≈ 0.40f0 && vd[2, 1] ≈ 0.55f0   # very dry
        @test vw[1, 1] ≈ 0.16f0 && vw[2, 1] ≈ 1.5f0                          # very wet
        # 1-hr moisture rises monotonically from very dry → very wet
        @test fuel_moisture(1)[1,1] < fuel_moisture(2)[1,1] < fuel_moisture(3)[1,1] < fuel_moisture(4)[1,1]
    end

    @testset "canopy wind reduction (CANCLS/CORFAC)" begin
        @test fire_wind_reduction(2f0)   ≈ 0.5f0     # below first class
        @test fire_wind_reduction(17.5f0) ≈ 0.3f0    # exact class breakpoint
        @test fire_wind_reduction(37.5f0) ≈ 0.2f0
        @test fire_wind_reduction(90f0)  ≈ 0.1f0     # above last class
        @test fire_wind_reduction(11.25f0) ≈ 0.4f0   # midway 5↔17.5 ⇒ midway 0.5↔0.3
        @test fire_wind_reduction(10f0) > fire_wind_reduction(50f0)  # denser canopy shelters more
    end

    @testset "input layer → Rothermel → mortality chain" begin
        coef = coefficients(Southern())
        # very-dry fuel moisture, a grassy/timber surface fuel, 20-mph open wind under 50% canopy
        mois = fuel_moisture(1)
        load = zeros(Float32, 2, 4); sav = zeros(Float32, 2, 4)
        load[1, 1] = 0.05f0; sav[1, 1] = 2000f0
        load[1, 2] = 0.10f0; sav[1, 2] = 109f0
        fwind = 20f0 * fire_wind_reduction(50f0)
        r = rothermel_surface_fire(load, sav, 0.4f0, 0.20f0, mois; wind = fwind)
        @test r.byram > 0f0 && r.flame > 0f0
        # carry the Byram intensity through scorch → CSV → mortality (a thin-barked oak)
        sch = scorch_height(r.byram, 77f0, fwind)
        csv = crown_volume_scorched(sch, 50f0, 40)
        pmort = fire_tree_mortality(coef, 64, 6f0, r.flame, csv)   # small scarlet oak
        @test 0f0 <= pmort <= 1f0
        @test pmort > fire_tree_mortality(coef, 64, 20f0, r.flame, csv)  # big tree survives better
    end

    @testset "dynamic fuel model (FMCFMD3) → fire, full integration" begin
        s = StandState(Southern()); init_blockdata!(s, s.variant); init_merch_standards!(s)
        s.plot.forest_type = Int32(520)
        s.plot.latitude = 35f0; s.plot.longitude = -80f0; s.plot.elevation = 10f0
        t = s.trees; t.n = 4
        for (i, (sp, d, h, ic)) in enumerate(((65,14f0,72f0,40),(33,8f0,55f0,45),
                                              (22,4f0,18f0,50),(65,2f0,8f0,55)))
            t.species[i]=Int32(sp); t.dbh[i]=d; t.height[i]=h; t.tpa[i]=30f0; t.crown_pct[i]=Int32(ic)
        end
        s.fire = FireState(); s.fire.active = true
        fmcba!(s)                                       # populate cwd / flive / cover
        mois = fuel_moisture(1)                          # very dry
        load, sav, depth, mext = build_dynamic_fuel_model(s, mois)

        # dead 1-hr load = 0–.25" + litter pools (lb/ft²); depth and Mx are positive
        currcwd1 = sum(@view s.fire.cwd[1, :, :]) * 0.04591f0
        currcwd10 = sum(@view s.fire.cwd[10, :, :]) * 0.04591f0
        @test load[1, 1] ≈ currcwd1 + currcwd10
        @test depth > 0f0 && 0f0 < mext < 1f0
        # the understory trees (≤ 6 ft: the 2-ft sapling) put live-woody load in via crown biomass
        @test load[2, 1] > 0f0
        @test sav[1, 1] == 2000f0 && sav[1, 2] == 109f0  # USAV / 10-hr SAV
        # dry herb ⇒ part of the live herb moves to the dead-herb class
        @test load[1, 4] >= 0f0

        # the constructed model carries a fire under very-dry/windy conditions
        fwind = 20f0 * fire_wind_reduction(s.fire.percov)
        r = rothermel_surface_fire(load, sav, depth, mext, mois; wind = fwind)
        @test r.byram > 0f0 && r.flame > 0f0
        # determinism
        l2, s2, d2, m2 = build_dynamic_fuel_model(s, mois)
        @test l2 == load && d2 == depth && m2 == mext
    end
end
