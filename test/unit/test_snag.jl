# Tests for the snag falldown + decay dynamics (FFE F7, FMSFALL + DECAYX).
using FVSjl: snag_fall_density, snag_decay_fraction, coefficients, Southern, coef_col,
             add_snag!, update_snags!, snag_standing_density, FireState, StandState,
             init_blockdata!, init_merch_standards!, fmburn!

@testset "snag falldown + decay (FMSFALL)" begin
    coef = coefficients(Southern())

    # reference mirroring fmsfall.f for a record at full density
    function ref_fall(ksp, d, origden, denttl)
        base = max(0.01, -0.001679*d + 0.064311)
        modrate = min(1.0, base * coef_col(coef, :snag_fallx)[ksp])
        if d < 12 && ksp != 2
            return modrate * origden
        end
        alldwn = coef_col(coef, :snag_alldwn)[ksp]
        x = (0.05 - 1) / (-modrate)
        fallm2 = alldwn <= x ? 2.0 : 0.05 / (alldwn - x)
        if denttl <= 0.05*origden
            return fallm2 * origden
        end
        dfalln = modrate * origden
        if denttl < dfalln + 0.05*origden
            dfalln = denttl - origden*(0.05 - fallm2)
        end
        return dfalln
    end

    @testset "small-snag linear fall" begin
        for (sp, d) in ((5, 8f0), (65, 8f0), (33, 6f0))
            @test isapprox(snag_fall_density(coef, sp, d, 100f0, 100f0), Float32(ref_fall(sp, Float64(d), 100, 100)); atol = 2f-6)   # jl F32 vs F64 ref_fall: few-ULP piecewise-formula rounding (measured max 1.9e-6)
        end
        # fast-falling pine (snag class 1) falls more than average-class oak
        @test snag_fall_density(coef, 5, 8f0, 100f0, 100f0) > snag_fall_density(coef, 65, 8f0, 100f0, 100f0)
        # bigger snags stand longer ⇒ fewer fall per year
        @test snag_fall_density(coef, 65, 4f0, 100f0, 100f0) > snag_fall_density(coef, 65, 10f0, 100f0, 100f0)
    end

    @testset "large-snag last-5% logic" begin
        # a large oak snag (≥12") uses the ALLDWN ramp
        @test isapprox(snag_fall_density(coef, 65, 16f0, 100f0, 100f0), Float32(ref_fall(65, 16.0, 100, 100)); atol = 2f-6)   # jl F32 vs F64 ref_fall: few-ULP piecewise-formula rounding (measured max 1.9e-6)
        # redcedar (sp 2) keeps the last-5% logic even when small
        @test isapprox(snag_fall_density(coef, 2, 8f0, 100f0, 100f0), Float32(ref_fall(2, 8.0, 100, 100)); atol = 2f-6)   # jl F32 vs F64 ref_fall: few-ULP piecewise-formula rounding (measured max 1.9e-6)
        # at/below 5% remaining, the final fall rate clears the remainder
        f = snag_fall_density(coef, 65, 16f0, 100f0, 4f0)
        @test isapprox(f, Float32(ref_fall(65, 16.0, 100, 4)); atol = 2f-6)   # jl F32 vs F64 ref_fall: few-ULP piecewise-formula rounding (measured max 1.9e-6)
        @test f >= 0f0
    end

    @testset "decay fraction (DECAYX)" begin
        @test snag_decay_fraction(coef, 5)  == 0.07f0   # fast snag class 1
        @test snag_decay_fraction(coef, 65) == 0.21f0   # average snag class 2
        @test snag_decay_fraction(coef, 2)  == 0.35f0   # slow snag class 3 (redcedar)
    end

    @testset "snag list: creation + per-cycle aging" begin
        s = StandState(Southern()); init_blockdata!(s, s.variant)
        s.fire = FireState()
        # add two cohorts; a zero-density add is a no-op
        add_snag!(s.fire, 65, 14f0, 40f0, 2003)
        add_snag!(s.fire, 5,  10f0, 25f0, 2003)
        add_snag!(s.fire, 33,  8f0,  0f0, 2003)        # no-op
        @test length(s.fire.snags.sp) == 2
        @test snag_standing_density(s.fire) ≈ 65f0
        @test all(s.fire.snags.den_soft .== 0f0)       # new snags start hard

        # age 5 years: some fall (transfer to CWD), some hard→soft. update_snags! advances each snag by
        # its OWN age (current cycle year − death year, capped at nyears) — FMSNAG ages by death year, so
        # set the current year 5 yrs past the 2003 deaths.
        s.control.cycle_year[1] = Int32(2008)
        fell = update_snags!(s, 5)
        @test fell > 0f0
        @test snag_standing_density(s.fire) < 65f0
        @test snag_standing_density(s.fire) ≈ 65f0 - fell
        # Snags created HARD stay in the HARD pool for the FALL — FVS's DENIH/DENIS are the snag's INITIAL
        # state at creation (these were created hard via add_snag!), and the per-snag HARD decay flag that
        # flips at DKTIME does NOT move the fall density (verified vs instrumented FMSNAG: DFIS=0). So the
        # SCNV 0.80 soft factor applies only to snags SEEDED soft, never to softened mortality snags.
        @test all(s.fire.snags.den_soft .== 0f0)
        # the fast-falling pine cohort loses a larger fraction than the oak
        sn = s.fire.snags
        oak_left = (sn.den_hard[1] + sn.den_soft[1]) / sn.origden[1]
        pine_left = (sn.den_hard[2] + sn.den_soft[2]) / sn.origden[2]
        @test pine_left < oak_left

        # age well past DKTIME: still HARD (no density transition) — the bole keeps falling to the hard pool
        s.control.cycle_year[1] = Int32(2050)
        update_snags!(s, 5)
        @test all(s.fire.snags.den_soft .== 0f0)

        # the fallen snags transferred biomass into the down-wood (CWD) pools
        @test sum(s.fire.cwd) > 0f0
        # a 14" oak snag → the 12–20" down-wood size class (6)
        @test sum(@view s.fire.cwd[6, :, :]) > 0f0

        # eventually nearly all fall (advance to 25 yrs past death)
        s.control.cycle_year[1] = Int32(2028)
        update_snags!(s, 20)
        @test snag_standing_density(s.fire) < 1f0
    end

    @testset "fire creates snags (fmburn! integration)" begin
        s = StandState(Southern()); init_blockdata!(s, s.variant); init_merch_standards!(s)
        s.plot.forest_type = Int32(520)
        s.plot.latitude = 35f0; s.plot.longitude = -80f0; s.plot.elevation = 10f0
        t = s.trees; t.n = 2
        t.species[1]=Int32(65); t.dbh[1]=14f0; t.height[1]=72f0; t.tpa[1]=30f0; t.crown_pct[1]=Int32(40)
        t.species[2]=Int32(22); t.dbh[2]= 4f0; t.height[2]=18f0; t.tpa[2]=30f0; t.crown_pct[2]=Int32(50)
        s.fire = FireState(); s.fire.active = true
        res = fmburn!(s; wind = 20f0, fmois = 1, year = 2003)
        @test res.killed > 0f0
        # the killed TPA became standing snags
        @test snag_standing_density(s.fire) ≈ res.killed
        @test all(s.fire.snags.year .== 2003)
    end
end

# LS snag dynamics differ from SN/CS: FMSFALL uses the "new equation" BASE = −0.006·d+0.18 (fmsfall.f:128,
# a MUCH faster fall than SN's −0.001679·d+0.064311), the small-snag linear-fall breakpoint is 18" (12" for
# cedar/tamarack ksp 10,11,14), and snags LOSE HEIGHT over time via non-zero default HTX per snag class
# (fmvinit.f:823-875) — SN/NE keep HTX=0. These drive the LS Stand-Dead snag pool.
@testset "LS snag dynamics (FMSFALL new-equation + FMSNGHT height loss)" begin
    using FVSjl: LakeStates, ffe_snag_height_loss!, snag_bole_carbon, KeywordReader,
                 KeywordRecord, kw_fmin!, each_stand, notre!, setup_growth!, compute_volumes!
    lc = coefficients(LakeStates())

    @testset "LS FMSFALL base rate (−0.006·d+0.18) vs SN" begin
        # sp 10 (snag class 6, FALLX 0.53), DBH 11 (small → linear): MODRATE = (−0.006·11+0.18)·0.53
        modrate = (-0.006f0*11f0 + 0.18f0) * coef_col(lc, :snag_fallx)[10]
        @test snag_fall_density(lc, 10, 11f0, 50f0, 50f0; variant = LakeStates()) ≈ modrate * 50f0
        # the LS rate is far faster than the SN formula would give for the same snag
        sn_modrate = min(1f0, max(0.01f0, -0.001679f0*11f0 + 0.064311f0) * coef_col(lc, :snag_fallx)[10])
        @test snag_fall_density(lc, 10, 11f0, 50f0, 50f0; variant = LakeStates()) > 2f0 * sn_modrate * 50f0
        # LS linear-fall breakpoint is 18" (a 15" non-cedar snag still falls linearly, unlike SN's 12")
        @test snag_fall_density(lc, 15, 15f0, 50f0, 50f0; variant = LakeStates()) ≈
              clamp((-0.006f0*15f0 + 0.18f0) * coef_col(lc, :snag_fallx)[15], 0.01f0, 1f0) * 50f0
    end

    @testset "LS default HTX seeded by FMIN (fmvinit height-loss classes)" begin
        # LS default HTX by snag class (fmvinit.f:831-873): class1=3.0, 2=1.0, 3/4=0, 5=0.65, 6=0.45, hemlock=0
        htx = coef_col(lc, :snag_htx)
        @test htx[6] == 3.0f0    # jack-pine group (class 1)
        @test htx[1] == 1.0f0    # class 2
        @test htx[15] ≈ 0.65f0  # maple/ash group (class 5)
        @test htx[10] ≈ 0.45f0  # cedar/oak group (class 6)
        @test htx[12] ≈ 0.0f0   # hemlock exception (fmvinit.f:873)
    end
end
