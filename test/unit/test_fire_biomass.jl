# Unit tests for the Jenkins tree-biomass model (FFE F1, FMCBIO).
# Expected values hand-computed directly from fmcbio.f's equations + coefficients.
using FVSjl: jenkins_biomass, coefficients, Southern, coef_col

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
    @test all(jenkins_biomass(coef, 1, 10.0f0) .≈ ref(3, 1, 10.0f0, 4.0f0))
    # sp65 = group 9 (hardwood, jgrp 2)
    dm65 = coef_col(coef, :dbh_min)[65]
    @test all(jenkins_biomass(coef, 65, 8.0f0) .≈ ref(9, 2, 8.0f0, dm65))
    # sp88 = group 1 (softwood)
    dm88 = coef_col(coef, :dbh_min)[88]
    @test all(jenkins_biomass(coef, 88, 15.0f0) .≈ ref(1, 1, 15.0f0, dm88))

    # small-tree (< 2.5 cm) scaling branch
    @test all(jenkins_biomass(coef, 1, 0.5f0) .≈ ref(3, 1, 0.5f0, 4.0f0))

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
    @test coef_col(coef, :v2t)[1]  ≈ 20.6f0
    @test coef_col(coef, :v2t)[2]  ≈ 27.4f0
    @test coef_col(coef, :v2t)[22] ≈ 34.9f0
    # snag decay/fall classes (fmvinit.f) — feed the snag dynamics (F3/F7)
    @test coef_col(coef, :dkr_cls)[1]  == 4f0
    @test coef_col(coef, :snag_cls)[2] == 3f0
    @test coef_col(coef, :tfall_cls)[2] == 1f0
    # every species populated (1..90, no gaps)
    @test all(coef_col(coef, :v2t)[s] > 0f0 for s in 1:90)
    @test all(1 <= coef_col(coef, :ls_spi)[s] <= 68 for s in 1:90)
end
