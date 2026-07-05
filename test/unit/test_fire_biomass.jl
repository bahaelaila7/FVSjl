# Unit tests for the Jenkins tree-biomass model (FFE F1, FMCBIO).
# Expected values hand-computed directly from fmcbio.f's equations + coefficients.
using FVSjl: jenkins_biomass, coefficients, Southern, coef_col,
             crown_biomass, StandState, init_blockdata!, init_merch_standards!,
             ffe_dead_fuel_type, ffe_live_fuel_type, ffe_dead_fuel_loading, ffe_live_fuel_loading,
             ffe_forest_type, fmcba!, FireState

@testset "Jenkins tree biomass (FMCBIO)" begin
    coef = coefficients(Southern())
    intocm = 2.54f0; ktoti = 1.102311f0 / 1000f0

    # reference implementation mirroring fmcbio.f for a given (group, jgrp, dbh)
    function ref(igrp, jgrp, dbh, dbhmin)
        b0a = (-2.0336f0,-2.2304f0,-2.5384f0,-2.5356f0,-2.0773f0,
               -2.2094f0,-1.9123f0,-2.4800f0,-2.0127f0,-0.7152f0)[igrp]
        b1a = (2.2592f0,2.4435f0,2.4814f0,2.4349f0,2.3323f0,
               2.3867f0,2.3651f0,2.4835f0,2.4342f0,1.7029f0)[igrp]
        b0m = (-0.3737f0,-0.3065f0)[jgrp]; b1m = (-1.8055f0,-5.4240f0)[jgrp]
        b0b = (-1.5619f0,-1.6911f0)[jgrp]; b1b = (0.6614f0,0.8160f0)[jgrp]
        dcm = dbh * intocm
        if dcm >= 2.5f0
            a = exp(b0a + b1a*log(dcm)); r = a*exp(b0b + b1b/dcm)
        else
            a = exp(b0a + b1a*log(2.5f0))*(dcm/2.5f0); r = a*exp(b0b + b1b/2.5f0)
        end
        m = dbh >= dbhmin ? a*exp(b0m + b1m/dcm) : 0f0
        (a*ktoti, m*ktoti, r*ktoti)
    end

    # sp1 = group 3 (softwood, jgrp 1), dbh_min 4.0
    @test all(jenkins_biomass(coef, 1, 10.0f0) .== ref(3, 1, 10.0f0, 4.0f0))
    # sp65 = group 9 (hardwood, jgrp 2)
    dm65 = coef_col(coef, :dbh_min)[65]
    @test all(jenkins_biomass(coef, 65, 8.0f0) .== ref(9, 2, 8.0f0, dm65))
    # sp88 = group 1 (softwood)
    dm88 = coef_col(coef, :dbh_min)[88]
    @test all(jenkins_biomass(coef, 88, 15.0f0) .== ref(1, 1, 15.0f0, dm88))

    # small-tree (< 2.5 cm) scaling branch
    @test all(jenkins_biomass(coef, 1, 0.5f0) .== ref(3, 1, 0.5f0, 4.0f0))

    # merch is zero below the species merch-DBH limit, non-zero above it
    @test jenkins_biomass(coef, 1, 3.0f0)[2] == 0f0          # 3 < 4.0 dbh_min
    @test jenkins_biomass(coef, 1, 12.0f0)[2] > 0f0
    # aboveground monotonic increasing in DBH; root < aboveground
    a5  = jenkins_biomass(coef, 65, 5.0f0)
    a15 = jenkins_biomass(coef, 65, 15.0f0)
    @test a15[1] > a5[1]
    @test a15[3] < a15[1]
    # zero / non-positive DBH → all zero
    @test jenkins_biomass(coef, 1, 0.0f0) == (0f0, 0f0, 0f0)
end

@testset "FFE species properties (fmvinit.f / fmcrow.f ISPMAP)" begin
    coef = coefficients(Southern())
    # SN→LS species map (ISPMAP, fmcrow.f:148) — used to pick the crown-biomass equation set
    @test coef_col(coef, :ls_spi)[1]  == 8f0     # fir → LS 8 (true fir/hemlock)
    @test coef_col(coef, :ls_spi)[2]  == 14f0    # redcedar → LS 14
    @test coef_col(coef, :ls_spi)[22] == 26f0    # sugar maple → LS 26
    @test coef_col(coef, :ls_spi)[65] == 34f0    # black willow group → LS 34 (red oak set)
    # V2T wood specific gravity (lb/cuft, fmvinit.f) — the SG passed to FMCROWE
    @test coef_col(coef, :v2t)[1]  == 20.6f0
    @test coef_col(coef, :v2t)[2]  == 27.4f0
    @test coef_col(coef, :v2t)[22] == 34.9f0
    # snag decay/fall classes (fmvinit.f) — feed the snag dynamics (F3/F7)
    @test coef_col(coef, :dkr_cls)[1]  == 4f0
    @test coef_col(coef, :snag_cls)[2] == 3f0
    @test coef_col(coef, :tfall_cls)[2] == 1f0
    # snag fall/decay dynamics by snag class (fmvinit.f:1060-1086), with species overrides
    @test coef_col(coef, :snag_decayx)[5] == 0.07f0  # shortleaf pine = fast snag class 1
    @test coef_col(coef, :snag_fallx)[5]  == 7.17f0
    @test coef_col(coef, :snag_alldwn)[5] == 50.0f0  # pine ALLDWN override
    @test coef_col(coef, :snag_alldwn)[2] == 100.0f0 # redcedar ALLDWN override
    @test coef_col(coef, :snag_decayx)[65] == 0.21f0 # oak = average snag class 2
    @test coef_col(coef, :snag_alldwn)[65] == 15.0f0
    # every species populated (1..90, no gaps)
    @test all(coef_col(coef, :v2t)[s] > 0f0 for s in 1:90)
    @test all(1 <= coef_col(coef, :ls_spi)[s] <= 68 for s in 1:90)
end

@testset "crown biomass by size class (FMCROWE) — structural" begin
    # FMCROWE has no .sum-visible output and feeds canopy bulk density (validated only
    # at F5/F6), so these tests pin structure/invariants, not absolute bit-exactness.
    s = StandState(Southern()); init_blockdata!(s, s.variant); init_merch_standards!(s)
    cb(sp, d, h, ic) = crown_biomass(s, sp, Float32(d), Float32(h), ic)

    # degenerate trees → all zero
    @test cb(65, 0.0, 70.0, 40) == (0f0, 0f0, 0f0, 0f0, 0f0, 0f0)
    @test cb(65, 12.0, 0.0, 40) == (0f0, 0f0, 0f0, 0f0, 0f0, 0f0)

    for (sp, d, h, ic) in ((65, 12.0, 70.0, 40), (1, 10.0, 60.0, 50),
                           (22, 8.0, 55.0, 45), (33, 15.0, 72.0, 55))
        xv = cb(sp, d, h, ic)
        @test length(xv) == 6
        @test xv[1] > 0f0                       # XV[0]: foliage is non-trivial
        @test xv[6] == 0f0                      # XV[5]: unused
        @test all(x -> x >= -1f-3, xv)          # size classes non-negative (bar Float32 noise)
        @test sum(xv) > 0f0                      # total crown weight positive
        @test xv == cb(sp, d, h, ic)            # deterministic
    end

    # species map to different proportion forms (oak vs pine vs maple) ⇒ distinct splits
    oak  = cb(65, 12.0, 70.0, 40)               # ls_spi 34 → red-oak form
    pine = cb(1,  12.0, 70.0, 40)               # ls_spi 8  → shortleaf-pine form
    @test oak != pine
end

@testset "FFE surface fuel loading (FMCBA / FUINI / FULIV)" begin
    coef = coefficients(Southern())
    # dead-fuel forest-type classification (fmcba.f:260)
    @test ffe_dead_fuel_type(101) == 1          # eastern white pine
    @test ffe_dead_fuel_type(161) == 3          # loblolly–shortleaf
    @test ffe_dead_fuel_type(520) == 6          # oak–hickory (snt01's forest type)
    @test ffe_dead_fuel_type(805) == 9          # maple–beech–birch
    @test ffe_dead_fuel_type(999) == 6          # default → oak–hickory
    # live-fuel classification (fmcba.f:161)
    @test ffe_live_fuel_type(4) == 1            # pines
    @test ffe_live_fuel_type(7) == 3            # redcedar
    @test ffe_live_fuel_type(6) == 4            # oak savannah
    @test ffe_live_fuel_type(2) == 2            # default → hardwoods
    # dead-fuel loadings (FUINI) — oak–hickory (type 6): litter 4.28, duff 5.91
    fd = ffe_dead_fuel_loading(coef, 520)
    @test length(fd) == 11
    @test fd[1]  == 0.13f0                        # <0.25" class
    @test fd[3]  == 1.93f0                        # 1–3" class
    @test fd[10] == 4.28f0                        # litter
    @test fd[11] == 5.91f0                        # duff
    # white pine (type 1): duff is the heaviest at 12.52
    @test ffe_dead_fuel_loading(coef, 103)[11] == 12.52f0
    # live herb/shrub loadings (FULIV)
    @test ffe_live_fuel_loading(coef, 4) == (0.1f0, 0.25f0)     # pines
    @test ffe_live_fuel_loading(coef, 7) == (1.0f0, 5.0f0)      # redcedar
end

@testset "FFE forest-type classifier (FMSNFT)" begin
    s = StandState(Southern()); init_blockdata!(s, s.variant)
    ft(code) = (s.plot.forest_type = Int32(code); ffe_forest_type(s))
    s.trees.n = 0
    @test ft(520) == 1          # oak-hickory hardwood (snt01 stand 4)
    @test ft(181) == 7          # eastern redcedar
    @test ft(605) == 8          # St. Francis type
    @test ft(999) == 9          # nonstocked
    @test ft(888) == 1          # default → hardwood
    @test ft(161) == 0          # pine type but no trees → IFFEFT unset (0)
    # pine-fraction split for a loblolly-shortleaf (161) stand
    t = s.trees; s.plot.forest_type = Int32(161)
    t.n = 2
    t.species[1] = Int32(13); t.dbh[1] = 12.0f0; t.tpa[1] = 100.0f0   # loblolly pine
    t.species[2] = Int32(65); t.dbh[2] = 12.0f0; t.tpa[2] = 10.0f0    # an oak
    @test ffe_forest_type(s) == 4                                     # >70% pine BA → pine
    t.tpa[2] = 100.0f0                                                # 50/50 pine/hardwood
    @test ffe_forest_type(s) == 2                                     # ≤50% pine → hardwood/pine
end

@testset "FFE per-cycle fuel/cover update (FMCBA)" begin
    function build()
        s = StandState(Southern()); init_blockdata!(s, s.variant)
        s.plot.forest_type = Int32(520)
        s.plot.latitude = 35f0; s.plot.longitude = -80f0; s.plot.elevation = 10f0
        t = s.trees; t.n = 3
        t.species[1]=Int32(65); t.dbh[1]=14f0; t.height[1]=72f0; t.tpa[1]=40f0; t.crown_pct[1]=Int32(40)
        t.species[2]=Int32(33); t.dbh[2]= 8f0; t.height[2]=55f0; t.tpa[2]=30f0; t.crown_pct[2]=Int32(45)
        t.species[3]=Int32(22); t.dbh[3]= 6f0; t.height[3]=45f0; t.tpa[3]=20f0; t.crown_pct[3]=Int32(50)
        s.fire = FireState(); s.fire.active = true
        s
    end
    s = build(); fmcba!(s); fs = s.fire

    @test fs.covtyp == 65                       # oak carries the most basal area
    @test fs.bigdbh == 14f0
    @test 0f0 < fs.percov <= 100f0
    # forest 520 → hardwood (IFFEFT 1) → hardwood live fuels FULIV[2]
    @test fs.flive == (0.01f0, 0.03f0)
    # the BA-weighted decay-class split conserves each size class's FUINI total. The split multiplies fd[isz]
    # by per-class Float32 fractions that sum back to 1 only to Float32 precision, so Σcwd == fd·(Σfrac) differs
    # from fd by a NON-ASSOCIATIVE SUM-ORDER ULP (measured max|Δ|=1.19e-7 = 1 Float32 eps). atol 2f-7 = that
    # sum-order width (was the loose ≈ default rtol≈3.4e-4).
    fd = ffe_dead_fuel_loading(coefficients(Southern()), 520)
    for isz in 1:11
        @test isapprox(sum(@view fs.cwd[isz, 2, :]), fd[isz]; atol = 2f-7)
    end
    # determinism
    s2 = build(); fmcba!(s2)
    @test s2.fire.covtyp == fs.covtyp && s2.fire.percov == fs.percov

    # no-op when FFE is inactive
    s3 = build(); s3.fire.active = false; fmcba!(s3)
    @test s3.fire.covtyp == 0 && s3.fire.percov == 0f0 && all(s3.fire.cwd .== 0f0)
end
