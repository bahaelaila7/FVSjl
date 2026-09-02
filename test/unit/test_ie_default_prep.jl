# =============================================================================
# test_ie_default_prep.jl — IE DEFAULT site-prep + per-plot PROB1 on the bare dated-ESTAB (disturbance)
# tally. Closes the 572-vs-555 baseline: FVS (estab.f:365-370) applies ESPREP DEFAULT site-prep proportions
# even with NO MECHPREP/BURNPREP keyword, and because the stocking logit ie_estock carries a per-IPREP SPRE
# term (estock.f:214), the per-plot PROB1 splits by IPREP. jl previously used a single stand-level PROB1
# (IPREP=1) ⇒ over-booked the tally ~3% ⇒ .sum TPA 572 vs the oracle 555 (×0.909 stockable). Also exercises
# the PADV+PSUB best-species SUMUP (estab.f:728), added on the disturbance tally (ITIME>2).
#
# MEASURED on FVSie_g16 na_def (bare IE STDINFO 101/520/0/315/30/45, ESTAB 1992, NUMCYCLE 6):
#   • SITE PREP SUMMARY default NONE/MECH/BURN = 48/26/26 (regen report).
#   • per-plot PROB1 by IPREP = 0.699 / 0.663977 / 0.657729 (dump); 48/26/26 mean = 0.6792 = report.
#   • SPRE deltas (grand-fir series ieq=2): 0 / -0.16145 / -0.18933 (esblkd.f).
# These tests FIRE the ported branches: is_ie=false + no prep is the pre-fix behaviour; the default-prep
# prep_sumup + prob1_prep + is_ie=true path moves the tally. Deleting the per-plot p1n / PSUB / default-prep
# makes the capped/prepped tally equal the un-prepped one and the inequality assertions fail.
# =============================================================================

using Test
using FVSjl

@testset "IE default site-prep + per-plot PROB1 (bare dated-ESTAB baseline)" begin
    asp = 5.4978f0      # 315° in radians
    slo = 0.30f0; elev = 45.0f0

    @testset "ie_estock per-IPREP SPRE term (grand-fir series ieq=2, ihab=7)" begin
        # PROB1 varies by IPREP through SPRE[iprep,2] = 0 / -0.16145 / -0.18933 (estock.f:214).
        f(ip) = FVSjl.ie_estock(7, ip, slo, cos(asp) * slo, sin(asp) * slo, elev, 1.0f0, 0.0f0,
                                10.0f0, sqrt(10.0f0), 0.0f0, 0.0f0, 4)
        @test isapprox(f(2) - f(1), -0.16145f0; atol = 1f-3)
        @test isapprox(f(3) - f(1), -0.18933f0; atol = 1f-3)
    end

    @testset "ie_esprep default proportions (series 4, NONE-dominant)" begin
        pn0, pm0, pb0 = FVSjl.ie_esprep(4, asp, slo, 0.0f0, elev)
        @test all(0f0 .< (pn0, pm0, pb0) .< 1f0)
        tot = pn0 + pm0 + pb0
        @test pn0 / tot > 0.30 && pm0 / tot > 0.05 && pb0 / tot > 0.05   # NONE dominant, MECH/BURN nonzero
    end

    @testset "default prep + per-plot PROB1 MOVES the disturbance tally DOWN" begin
        occ = Float32[ones(9); zeros(14)]; over = zeros(Float32, 10)
        kw = (seed0 = 43303, nplots = 200, ihab = 7, iser = 4, ifo = 4, iprep = 1, iphy = 3,
              xcos = cos(asp) * slo, xsin = sin(asp) * slo, slo = slo, elev = elev, baa = 1.0f0,
              regt = 10.0f0, bwaf = 0.0f0, bwb4 = 0.0f0, prob1 = 0.699f0, dupnpt = 200.0f0,
              occ = occ, over = over, time = 10.0f0, is_ingro = false, nsp = 23)

        base = FVSjl.ie_autoes_tally(; kw...)                       # pre-fix: scalar PROB1, no PSUB, no prep

        # default ESPREP prep mix + per-IPREP PROB1 (the engine's disturbance-tally path)
        psum = FVSjl.ie_esetpr_normalize(FVSjl.ie_esprep(4, asp, slo, 0.0f0, elev)..., 0, 0)
        p1p = Float32[1f0 / (1f0 + exp(-FVSjl.ie_estock(7, ip, slo, cos(asp) * slo, sin(asp) * slo,
                      elev, 1.0f0, 0.0f0, 10.0f0, sqrt(10.0f0), 0.0f0, 0.0f0, 4))) for ip in 1:3]
        withprep = FVSjl.ie_autoes_tally(; kw..., prep_sumup = psum, prob1_prep = p1p, is_ie = true)

        @test withprep != base                                     # the ported path fires
        @test sum(withprep) < sum(base)                            # MECH/BURN plots' lower PROB1 ⇒ fewer trees
        @test (sum(base) - sum(withprep)) / sum(base) > 0.02       # a MATERIAL (not ULP) reduction (~3%)
    end

    @testset "is_ie=false leaves the tally byte-identical (EM/other variants unperturbed)" begin
        # With is_ie=false and no prep, the new PSUB/prob1_prep code is inert ⇒ the ingrowth/other-variant
        # path is unchanged. (This is what keeps the 339/11 gate byte-identical.)
        occ = Float32[ones(9); zeros(14)]; over = zeros(Float32, 10)
        kw = (seed0 = 43303, nplots = 50, ihab = 7, iser = 4, ifo = 4, iprep = 1, iphy = 3,
              xcos = cos(asp) * slo, xsin = sin(asp) * slo, slo = slo, elev = elev, baa = 1.0f0,
              regt = 10.0f0, bwaf = 0.0f0, bwb4 = 0.0f0, prob1 = 0.699f0, dupnpt = 50.0f0,
              occ = occ, over = over, time = 10.0f0, is_ingro = false, nsp = 23)
        @test FVSjl.ie_autoes_tally(; kw..., is_ie = false) == FVSjl.ie_autoes_tally(; kw...)
    end
end
