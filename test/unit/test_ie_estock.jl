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

@testset "IE ESNSPE P(#species) — task #143 chunk A2a" begin
    # iet01 stand-4 plot-1: ISER=4 (WH), ITPP=2, TPP=2, TPPLN=ln2, BAA=1, ELEV=34, REGT=1, BWAF=0.
    # XCOS=cos(asp)·SLO, XSIN=sin(asp)·SLO (SLO-weighted aspect). Oracle PSPE=(0.543,0.393,0,0,0,0).
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    psp = FVSjl.ie_esnspe(4, 2, 2.0f0, log(2.0f0), 1.0f0, 34.0f0, 1.0f0, 0.0f0, xcos, xsin, slo)
    @test isapprox(psp[1], 0.543f0; atol = 1f-3)
    @test isapprox(psp[2], 0.393f0; atol = 1f-3)
    @test psp[3] == 0f0 && psp[4] == 0f0 && psp[5] == 0f0 && psp[6] == 0f0  # ITPP=2 gates ≥3 off
    # ITPP=6 exercises all six count-logits (regression guard: finite & in (0,1))
    psp6 = FVSjl.ie_esnspe(4, 6, 10.0f0, log(10.0f0), 50.0f0, 34.0f0, 1.0f0, 0.0f0, xcos, xsin, slo)
    @test all(0f0 .< collect(psp6) .< 1f0)
end

@testset "IE ESPADV P(advance species) — task #143 chunk A2b" begin
    # iet01 stand-4 plot-1: IHAB=10, IPREP=1, IFO=4, IPHY≠1, TIME=1 (measured), BAA=1, ELEV=34, REGT=1,
    # BWAF=BWB4=0, occupancy=1 (OCURHT/OCURNF/XESMLT all 1 for sp1-9), OVER<9.95 (no bumps). XCOS/XSIN SLO-weighted.
    slo = 0.30f0; xcos = cos(5.498f0) * slo; xsin = sin(5.498f0) * slo
    # measured OCURHT(grp10,·)=1.0 for sp1-9, 0 for sp10+ (PP) — the occupancy that zeroes PP.
    occ = Float32[1, 1, 1, 1, 1, 1, 1, 1, 1, 0]; over = zeros(Float32, 10)
    padv = FVSjl.ie_espadv(10, 1, 4, 3, xcos, xsin, slo, 1.0f0, 1.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ, over)
    # oracle PADV (ie/espadv.f dump, 3 dp): WP WL DF GF WH RC LP ES AF PP
    oracle = (0.062f0, 0.005f0, 0.048f0, 0.485f0, 0.283f0, 0.122f0, 0.001f0, 0.014f0, 0.039f0, 0.0f0)
    for i in 1:10
        @test isapprox(padv[i], oracle[i]; atol = 6f-4)
    end
    # occupancy zeroes the species out
    occ0 = zeros(Float32, 10)
    padv0 = FVSjl.ie_espadv(10, 1, 4, 3, xcos, xsin, slo, 1.0f0, 1.0f0, 34.0f0, 1.0f0, 0.0f0, 0.0f0, occ0, over)
    @test all(collect(padv0) .== 0f0)
end
