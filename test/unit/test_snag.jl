# Tests for the snag falldown + decay dynamics (FFE F7, FMSFALL + DECAYX).
using FVSjl: snag_fall_density, snag_decay_fraction, coefficients, Southern, coef_col

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
            @test snag_fall_density(coef, sp, d, 100f0, 100f0) ≈ Float32(ref_fall(sp, Float64(d), 100, 100))
        end
        # fast-falling pine (snag class 1) falls more than average-class oak
        @test snag_fall_density(coef, 5, 8f0, 100f0, 100f0) > snag_fall_density(coef, 65, 8f0, 100f0, 100f0)
        # bigger snags stand longer ⇒ fewer fall per year
        @test snag_fall_density(coef, 65, 4f0, 100f0, 100f0) > snag_fall_density(coef, 65, 10f0, 100f0, 100f0)
    end

    @testset "large-snag last-5% logic" begin
        # a large oak snag (≥12") uses the ALLDWN ramp
        @test snag_fall_density(coef, 65, 16f0, 100f0, 100f0) ≈ Float32(ref_fall(65, 16.0, 100, 100))
        # redcedar (sp 2) keeps the last-5% logic even when small
        @test snag_fall_density(coef, 2, 8f0, 100f0, 100f0) ≈ Float32(ref_fall(2, 8.0, 100, 100))
        # at/below 5% remaining, the final fall rate clears the remainder
        f = snag_fall_density(coef, 65, 16f0, 100f0, 4f0)
        @test f ≈ Float32(ref_fall(65, 16.0, 100, 4))
        @test f >= 0f0
    end

    @testset "decay fraction (DECAYX)" begin
        @test snag_decay_fraction(coef, 5)  ≈ 0.07f0   # fast snag class 1
        @test snag_decay_fraction(coef, 65) ≈ 0.21f0   # average snag class 2
        @test snag_decay_fraction(coef, 2)  ≈ 0.35f0   # slow snag class 3 (redcedar)
    end
end
