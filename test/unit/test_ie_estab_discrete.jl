# =============================================================================
# test_ie_estab_discrete.jl — IE discrete-establishment PASSALL/PASMAX excess-tree cap wired into the live
# AUTOES tally (ie_autoes_tally), gated to InlandEmpire. estab.f:1079-1145 (NOTE best/excess split) +
# 1288-1370 (excess pass, es_pasmax_xcsmax). Ported 2026-09-02 on branch estab-ie-discrete.
#
# SEMANTIC (mapped from estab.f, MEASURED on the instrumented FVSie_g16 dump — na_def/na_p1, single-.o swap):
#   • Per plot, the excess regen count per species is EXCESS(I) = (#trees of I) − (best trees of I). The best
#     count NBEST = min(ITP, max(4, NUMSPE)) (measured bit-exact, 54/54 plots) and the NBEST trees are the
#     NUMSPE best-species trees + the (NBEST−NUMSPE) TALLEST excess trees.
#   • The excess pass books FTEMP2·300·min(EXCESS(I),PASMAX)/DUPNPT per species (FTEMP2=SUMESP/(EXCESS+1e-5));
#     PASMAX replaces EXCESS with min(EXCESS,PASMAX) — the es_pasmax_xcsmax kernel (BRKUP·XCSMAX=min(EXCESS,PASMAX)).
#   • So PASMAX only REDUCES a species' excess TPA when EXCESS(I) > PASMAX; otherwise the tally is byte-identical.
#     Default PASMAX=5 (CONFID); the AUTOES ingrowth path forces 15; PASSALL n sets it to n.
#
# These tests FIRE the ported branch: `pasmax=Inf32` (the default for every non-IE caller) is the pre-wiring
# behaviour, and a finite `pasmax` that is smaller than some plot's EXCESS moves the tally. Deleting the cap
# block in ie_autoes_tally makes the capped tally equal the uncapped one and the inequality assertions fail.
# =============================================================================

using Test
using FVSjl

@testset "IE discrete-establishment PASMAX excess cap (ie_autoes_tally, IE-gated)" begin
    # Dense iet01-style AUTOES tally (ihab=10, 50 plots) — the same fixture the #143 species split validates,
    # which has plots whose per-species EXCESS exceeds 1, so a low PASMAX bites.
    occ = Float32[ones(9); zeros(14)]; over = zeros(Float32, 10); slo = 0.30f0
    kw = (seed0 = 43303, nplots = 50, ihab = 10, iser = 4, ifo = 4, iprep = 1, iphy = 3,
          xcos = cos(5.498f0) * slo, xsin = sin(5.498f0) * slo, slo = slo, elev = 34.0f0, baa = 1.0f0,
          regt = 1.0f0, bwaf = 0.0f0, bwb4 = 0.0f0, prob1 = 0.5527f0, dupnpt = 50.0f0, occ = occ, over = over)

    base = FVSjl.ie_autoes_tally(; kw...)                        # pasmax=Inf32 default ⇒ uncapped (pre-wiring)
    cap1 = FVSjl.ie_autoes_tally(; kw..., pasmax = 1f0)          # PASSALL 1 ⇒ cap EXCESS→1 per species
    cap5 = FVSjl.ie_autoes_tally(; kw..., pasmax = 5f0)          # default CONFID=5
    big  = FVSjl.ie_autoes_tally(; kw..., pasmax = 1000f0)       # ≥ any plot's EXCESS (ITPP≤MAXTPP=25) ⇒ uncapped

    @testset "uncapped path is byte-identical (Inf == a PASMAX larger than any EXCESS)" begin
        # The cap block only fires when EXCESS(I) > PASMAX; with PASMAX≥1000 no species is ever capped, so the
        # tally must equal the Inf-default (this is what keeps the 339/11 gate + other variants byte-identical).
        @test big == base
    end

    @testset "a low PASMAX moves the tally DOWN (the cap fires)" begin
        @test sum(cap1) < sum(base)                              # PASSALL 1 caps ⇒ fewer excess trees
        @test sum(cap1) <= sum(cap5) <= sum(base)                # monotone in PASMAX
        @test all(cap1 .<= base .+ 1f-6)                         # every species non-increasing under the cap
        @test (sum(base) - sum(cap1)) / sum(base) > 0.05         # a MATERIAL (not ULP) reduction
    end

    # =========================================================================
    # Establishment-cohort HEIGHT-CLASS / WK4 emission (estab.f DO 99 + DO 33/228). The `emit` path books the
    # per-tree records FVS creates — advance (ESADVH, WK4≈0.60), subsequent (ESSUBH, WK4≈0.20/0.00) first trees
    # + tripled excess (ESXCSH mean height, WK4=STOMLT) — instead of one collapsed WK4=0.60 record per species.
    # The collapse over-projected the cohort (all TPA grown at the max multiplier). Two invariants:
    #   (1) TPA is CONSERVED: Σ emit-record TPA == Σ tally (the distribution repartitions, it does not add/remove).
    #   (2) The distribution is HETEROGENEOUS: WK4 spans {≈0.60 advance, ≈0.20/0.00 subsequent}, not all 0.60.
    # =========================================================================
    @testset "emit path: WK4 height-class distribution (advance/subsequent/excess)" begin
        emit = NTuple{5,Float64}[]
        # iet01 ingrowth tally (ihab=10 ⇒ IHTSER=MYHTS(10)=4), GENTIM=5 (FINT=10), ESPADV called (ingrowth).
        tally = FVSjl.ie_autoes_tally(; kw..., pasmax = 15f0, emit = emit, ihtser = 4, gentim = 5f0,
                                      call_espadv = true)
        @test !isempty(emit)
        # (1) TPA conservation: the emit records repartition exactly the same total TPA as the tally.
        emit_tot = sum(r[5] for r in emit)
        @test isapprox(emit_tot, sum(tally); rtol = 0.02)
        # (2) Heterogeneous WK4 — the whole point of the fix. Bucket the per-record WK4.
        wk4s = round.([r[4] for r in emit]; digits = 2)
        adv_tpa = sum(r[5] for r in emit if round(r[4]; digits = 2) >= 0.55)   # advance (WK4≈0.60)
        sub_tpa = sum(r[5] for r in emit if round(r[4]; digits = 2) <  0.55)   # subsequent/excess (WK4≤0.20)
        @test any(w -> w >= 0.55, wk4s)                          # some advance records (WK4≈0.60)
        @test any(w -> w <  0.55, wk4s)                          # some subsequent/excess (WK4≈0.20/0.00)
        @test sub_tpa > 0.10 * emit_tot                          # the low-WK4 tail is MATERIAL (not a rounding sliver)
        # Sanity: every record carries a valid species/height/positive TPA.
        @test all(1 <= Int(r[1]) <= 23 && r[3] > 0.0 && r[5] > 0.0 for r in emit)
    end
end
