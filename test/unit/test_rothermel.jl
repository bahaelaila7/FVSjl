# Tests for the Rothermel surface-fire model (FFE F5, FMFINT).
#
# Like the other inert FFE pieces, exact bit-for-bit agreement can only be confirmed
# against live Fortran once the fuel-model construction (F5b) + burn integration are in
# (FMFINT has no standalone .sum output). These tests pin the physical invariants, the
# Byram→flame relationship, and the too-moist cutoff — the behaviors a correct Rothermel
# implementation must satisfy — plus determinism.
using FVSjl: rothermel_surface_fire, fuel_moisture, fire_wind_reduction,
             scorch_height, crown_volume_scorched, fire_tree_mortality, coefficients, Southern

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
end
