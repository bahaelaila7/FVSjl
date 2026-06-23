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
