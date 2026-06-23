# Unit tests for fire-caused tree mortality (FFE F6, FMEFF + FMBRKT).
# Expected values come from a reference transcription of fmeff.f / fmbrkt.f.
using FVSjl: fire_bark_thickness, fire_mortality_group, fire_tree_mortality,
             fire_mortality_adjust, scorch_height, crown_volume_scorched,
             coefficients, Southern, coef_col

@testset "fire mortality (FMEFF) + bark thickness (FMBRKT)" begin
    coef = coefficients(Southern())
    b1 = [0.019,0.022,0.024,0.025,0.026,0.027,0.028,0.029,0.030,0.031,0.032,0.033,
          0.034,0.035,0.036,0.037,0.038,0.039,0.040,0.041,0.042,0.043,0.044,0.045,
          0.046,0.047,0.048,0.049,0.050,0.052,0.055,0.057,0.059,0.060,0.062,0.063,
          0.068,0.072,0.081]
    mb0 = [1.0229,0.1683,1.2165,0.8221,2.7750]
    mb1 = [-0.2646,-0.1332,-0.4758,-0.4098,-1.1224]
    mb2 = [2.6232,3.4152,6.0415,8.4682,2.8312]

    @testset "bark thickness" begin
        # general species: DBH · B1[EQNUM[sp]]
        for (sp, d) in ((1, 10f0), (65, 12f0), (33, 8f0), (90, 20f0))
            eq = Int(coef_col(coef, :bark_eqnum)[sp])
            @test fire_bark_thickness(coef, sp, Float32(d)) ≈ Float32(d * b1[eq])
        end
        # shortleaf pine (5): Harmon quadratic
        d = 10f0
        ref = max(0f0, (0.07f0 + 0.09f0*d*2.54f0 - 0.0001f0*(d*2.54f0)^2)/2.54f0)
        @test fire_bark_thickness(coef, 5, d) ≈ ref
    end

    @testset "mortality groups" begin
        @test fire_mortality_group(63) == 1 && fire_mortality_group(74) == 1
        @test fire_mortality_group(64) == 2 && fire_mortality_group(78) == 2
        @test fire_mortality_group(27) == 3
        @test fire_mortality_group(20) == 4
        @test fire_mortality_group(54) == 5
        @test fire_mortality_group(65) == 6 && fire_mortality_group(1) == 6
    end

    # reference mortality (mirrors fmeff.f)
    function ref_mort(sp, d, flame, csv)
        g = fire_mortality_group(sp)
        if 1 <= g <= 5
            charht = flame*0.7
            xm = -(mb0[g] + mb1[g]*d*2.54 + mb2[g]*charht/3.28)
            mn = log(1/0.000001 - 1)
            return xm >= mn ? 0.0 : 1/(1+exp(xm))
        else
            eq = Int(coef_col(coef, :bark_eqnum)[sp])
            bt = sp == 5 ? max(0.0,(0.07+0.09*d*2.54-0.0001*(d*2.54)^2)/2.54) : d*b1[eq]
            xm = exp(-1.941 + 6.316*(1-exp(-bt)) - 0.000535*csv*csv)
            return 1/(1+xm)
        end
    end

    @testset "tree mortality" begin
        for (sp, d, fl, csv) in ((64, 10f0, 8f0, 50f0),   # scarlet oak, group 2
                                 (27, 14f0, 6f0, 30f0),   # hickory, group 3
                                 (20,  6f0, 10f0, 80f0),  # red maple, group 4
                                 (54, 12f0, 5f0, 40f0),   # black gum, group 5
                                 (65, 10f0, 8f0, 50f0),   # northern oak, group 6 (Reinhardt)
                                 (1,  18f0, 4f0, 20f0),   # fir, group 6
                                 (5,  10f0, 8f0, 60f0))   # shortleaf pine, group 6 + special bark
            @test fire_tree_mortality(coef, sp, Float32(d), Float32(fl), Float32(csv)) ≈
                  Float32(ref_mort(sp, d, fl, csv))
        end
        # monotonicity: bigger flame ⇒ higher mortality (group-2 oak)
        @test fire_tree_mortality(coef, 64, 10f0, 12f0, 50f0) >
              fire_tree_mortality(coef, 64, 10f0, 4f0, 50f0)
        # bigger DBH (thicker bark) ⇒ lower mortality (same group-2 oak, same flame)
        @test fire_tree_mortality(coef, 64, 20f0, 8f0, 50f0) <
              fire_tree_mortality(coef, 64, 5f0, 8f0, 50f0)
        # probabilities stay in [0,1]
        for sp in (64, 27, 20, 54, 65, 1, 5), fl in (2f0, 8f0, 20f0)
            p = fire_tree_mortality(coef, sp, 10f0, fl, 50f0)
            @test 0f0 <= p <= 1f0
        end
    end

    @testset "scorch height (Van Wagner) + crown volume scorched" begin
        # SCH = (63/(140-ATEMP)) · BYRAM'^(7/6) / sqrt(BYRAM' + FWIND^3), BYRAM'=BYRAM/60
        sch_ref(byram, atemp, fwind) = (b = byram/60; (63/(140-atemp)) * (b^(7/6)/sqrt(b + fwind^3)))
        for (by, at, wd) in ((3000f0, 77f0, 5f0), (12000f0, 90f0, 10f0), (500f0, 60f0, 2f0))
            @test scorch_height(by, at, wd) ≈ Float32(sch_ref(by, at, wd))
        end
        @test scorch_height(12000f0, 90f0, 5f0) > scorch_height(3000f0, 90f0, 5f0)  # hotter ⇒ taller scorch

        # CSV: crown length CRL = HT·CR; scorch length SL = SCH-(HT-CRL), clamped to [0,CRL]
        function csv_ref(sch, ht, cr)
            crl = ht*(cr/100); crl <= 0 && return 100.0
            sl = sch - (ht - crl); sl = clamp(sl, 0.0, crl)
            100*(sl*(2*crl - sl)/(crl*crl))
        end
        for (sch, ht, cr) in ((40f0, 60f0, 40), (80f0, 60f0, 40), (10f0, 60f0, 40), (100f0, 50f0, 50))
            @test crown_volume_scorched(sch, ht, cr) ≈ Float32(csv_ref(sch, ht, cr))
        end
        @test crown_volume_scorched(5f0, 60f0, 0) == 100f0    # no crown → fully scorched
        @test crown_volume_scorched(0f0, 60f0, 40) == 0f0     # scorch below crown → 0
        @test crown_volume_scorched(100f0, 60f0, 40) == 100f0 # scorch past crown top → 100%
        # full chain: byram → scorch → CSV → mortality (group-6 species)
        by = 4000f0
        sch = scorch_height(by, 77f0, 5f0)
        csv = crown_volume_scorched(sch, 60f0, 40)
        flame = 0.45f0 * (by/60f0)^0.46f0
        p = fire_tree_mortality(coef, 1, 12f0, flame, csv)
        @test 0f0 <= p <= 1f0
    end

    @testset "SN has no FMEFF season/size adjustments" begin
        # The maple-<4″, hardwood-≤1″ and early-season reductions in fmeff.f are gated by
        # IF (VARACD .EQ. 'LS'/'ON'/'NE') and do NOT apply to SN — fire_mortality_adjust is
        # a no-op for every species/size/season (fmeff.f:278-326). The only universal SN
        # rule (dbh≤1 & csv>50 ⇒ 1.0, fmeff.f:330) lives in fmburn!, not here.
        for (sp, d, seas) in ((26, 3.0f0, 3), (65, 1.0f0, 3), (30, 6.0f0, 1),
                              (30, 2.0f0, 1), (65, 8.0f0, 1), (5, 8.0f0, 1), (26, 3.0f0, 1))
            @test fire_mortality_adjust(0.3f0, sp, d, seas) == 0.3f0
        end
    end
end
