# IE AUTOES ESTOCK P(stocking) — bit-exact vs live FVSie (task #143 chunk A1).
# Oracle: FVSie_clean on iet01 stand-4 (SHELTERWOOD WITH AUTO REGENERATION) with a DEBUG keyword
# dumps `PN FOR STOCKING= 0.2116` (estab.f:538) alongside its inputs (IHAB=10, IPREP=1, SLO=0.30,
# ASPECT=5.498, ELEV=34, BAA=1, TIME=1, SQREGT=1). IHAB=10 → IEQ=3 (cedar/hemlock series).
# See docs/AUTOES_CHUNK_PLAN.md for the full measured target set.
using Test
using FVSjl

@testset "IE ESTOCK P(stocking) — task #143 chunk A1" begin
    # measured iet01 stand-4 inputs → oracle PN = 0.2116 (cedar/hemlock series, IEQ=3)
    pn = FVSjl.ie_estock(10, 1, 0.30f0, cos(5.498f0), sin(5.498f0), 34.0f0,
                         1.00f0, log(1.00f0), 1.0f0, 1.0f0, 0.0f0, 0.0f0, 4)
    @test isapprox(pn, 0.2116f0; atol = 1f-3)
    # the caller forms the inventory stocking prob FTEMP = 1/(1+exp(-PN)) = 0.5527 (estab.f:540)
    ftemp = 1f0 / (1f0 + exp(-pn))
    @test isapprox(ftemp, 0.5527f0; atol = 1f-3)

    # IEQ dispatch coverage: IHAB=3→DF(1), 6→GF(2), 10→cedar/hemlock(3), 12→subalpine(4).
    # (regression guard: each branch must return a finite logit for nominal inputs)
    for ihab in (3, 6, 10, 12)
        p = FVSjl.ie_estock(ihab, 1, 0.30f0, cos(5.498f0), sin(5.498f0), 34.0f0,
                            50.0f0, log(50.0f0), 1.0f0, 1.0f0, 0.0f0, 0.0f0, 4)
        @test isfinite(p)
    end
    # IPREP>3 "all roads" branch must also compute
    @test isfinite(FVSjl.ie_estock(10, 4, 0.30f0, 0.7f0, -0.7f0, 34.0f0, 1.0f0, 0.0f0, 1.0f0, 1.0f0, 0.0f0, 0.0f0, 4))
end
