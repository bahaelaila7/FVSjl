# Unit tests for the stump-sprout sub-routines (NSPREC / SPRTHT / ESSPRT, SN).
# Expected values are hand-computed directly from essprt.f's SN SELECT CASE blocks.
using FVSjl: nsprec_sn, sprtht_sn, essprt_sn, coefficients, Southern

@testset "sprout sub-routines (NSPREC/SPRTHT/ESSPRT)" begin
    coef = coefficients(Southern())

    @testset "NSPREC sprout count" begin
        @test nsprec_sn(5, 6.0f0) == 1          # sp5: DSTMP<7 → 1
        @test nsprec_sn(5, 7.0f0) == 0          # sp5: DSTMP≥7 → 0
        @test nsprec_sn(33, 4.0f0) == 1         # oak grp: <5 → 1
        @test nsprec_sn(33, 5.0f0) == 1         # NINT(-1+0.4·5)=NINT(1.0)=1
        @test nsprec_sn(33, 7.5f0) == 2         # NINT(2.0)=2
        @test nsprec_sn(33, 8.75f0) == 3        # NINT(2.5)=3  (ties away from zero)
        @test nsprec_sn(33, 10.0f0) == 3        # NINT(3.0)=3
        @test nsprec_sn(33, 11.0f0) == 3        # >10 → 3
        @test nsprec_sn(82, 7.5f0) == 2         # same oak/sweetgum group
        @test nsprec_sn(65, 8.0f0) == 1         # non-grouped species → default 1
        @test nsprec_sn(99, 3.0f0) == 1         # out-of-range default
    end

    @testset "SPRTHT sprout height" begin
        @test sprtht_sn(33, 70.0f0, 5) ≈ (0.1f0 + 70f0/50f0) * 5f0   # curve sp
        @test sprtht_sn(59, 70.0f0, 5) ≈ (0.1f0 + 70f0/50f0) * 5f0   # upper range
        @test sprtht_sn(58, 70.0f0, 5) ≈ 0.5f0 + 0.5f0 * 5f0         # gap → default
        @test sprtht_sn(89, 70.0f0, 5) ≈ 0.5f0 + 0.5f0 * 5f0         # >87 → default
        @test sprtht_sn(5,  60.0f0, 3) ≈ (0.1f0 + 60f0/50f0) * 3f0   # explicit sp5
    end

    @testset "ESSPRT survival multiplier" begin
        # constant-multiplier species
        @test essprt_sn(coef, 5,  1.0f0, 6.0f0, 520) ≈ 1.0f0 * 0.42f0
        @test essprt_sn(coef, 33, 1.0f0, 6.0f0, 520) ≈ 1.0f0 * 0.93f0
        @test essprt_sn(coef, 22, 2.0f0, 6.0f0, 520) ≈ 2.0f0 * 0.73f0
        # logistic species (sp20: a=4.1975, b=-0.1821)
        @test essprt_sn(coef, 20, 1.0f0, 8.0f0, 520) ≈
              1f0 / (1f0 + exp(-(4.1975f0 - 0.1821f0 * 8f0)))
        # default logistic (sp65 falls through to a=2.7386, b=-0.1076)
        @test essprt_sn(coef, 65, 1.0f0, 8.0f0, 520) ≈
              1f0 / (1f0 + exp(-(2.7386f0 - 0.1076f0 * 8f0)))
        # forest-special sp64: special forest 809 uses the cubic-poly form
        @test essprt_sn(coef, 64, 1.0f0, 10.0f0, 809) ≈ (57.3f0 - 0.0032f0*1000f0)/100f0
        # forest-special sp64 in a common forest uses the ELSE logistic (3.8897, -0.2260)
        @test essprt_sn(coef, 64, 1.0f0, 10.0f0, 520) ≈
              1f0 / (1f0 + exp(-(3.8897f0 - 0.2260f0 * 10f0)))
        # sp77 special-forest inverse-logistic form
        @test essprt_sn(coef, 77, 1.0f0, 10.0f0, 905) ≈
              1f0 / (1f0 + exp(-(-2.8058f0 + 22.6839f0 *
                                  (1f0 / ((10f0/0.7788f0) - 0.4403f0)))))
    end
end
