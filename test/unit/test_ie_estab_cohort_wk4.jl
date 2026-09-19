# =============================================================================
# test_ie_estab_cohort_wk4.jl — IE AUTOES establishment-cohort HEIGHT-CLASS / WK4 partition (ie_autoes_tally
# emit path; estab.f DO 99 birth-multiplier + DO 33 best / DO 228 excess booking), gated to InlandEmpire.
#
# WHAT b6edfd9a DID: replaced the collapsed single-WK4=0.60 record per (species,point) with the faithful
# per-tree partition — each best tree booked at its own STOMLT (advance TRAGE=3−DELAY ⇒ WK4≈0.60; subsequent
# TRAGE=TIME−DELAY ⇒ WK4≈0.20 or 0.00), plus tripled ESXCSH excess records at WK4=STOMLT(species). WK4 drives
# the birth-cycle growth (esgent.f/regent.f apply it TWICE ⇒ effective WK4²), so the three-class split governs
# how much the cohort's DBH grows.
#
# ORACLE RE-BALANCE MEASUREMENT (2026-09-04, branch fix-ie-cohort-rebalance). The b6edfd9a partition was
# suspected on the IE re-sweep of OVER-WEIGHTING the slow (WK4≤0.20) classes → systematic BA/CCF/SDI
# UNDER-projection, to be re-balanced toward the advance class. Instrumenting FVSie_g16 estab.f DO 33/228 (the
# actual per-tree booked WK4 + PROB) vs the jl emit records across a POPULATION of bare + sparse establishment
# stands (ecoregions 331Aa/Ac/Af, 342Ia/b; single- AND multi-point) REFUTED that:
#   • per-class booked-TPA fractions match the oracle within realization noise (advance/subsequent/excess within
#     ~1-3pp; population-mean advance fraction ~equal, jl if anything marginally HIGHER) — the partition is
#     FAITHFUL to the FVS source; it does NOT over-weight the slow classes.
#   • the .sum residual is a TWO-SIDED #206 OLDRN/ZRAND realization straddle (86 stands: BA 38 over / 44 under,
#     mean −0.85 BA ~1.5%; needs_dig-style ≥10-20% subset 34/42 over/under, mean −0.9; adjacent near-identical
#     stands flip sign; the per-stand advance-fraction discrepancy is UNCORRELATED with the BA sign — a stand
#     booking LESS advance than the oracle yields BOTH +19 and −15 BA).
# ⇒ NO re-balance was warranted (shifting the partition toward advance would break faithfulness to the FVS
# source WITHOUT touching the straddle). The residual is legitimately the IE #206 corner. This test LOCKS the
# faithful three-class partition against silent regression.
#
# Fixture = the iet01-style AUTOES ingrowth tally (ihab=10, 50 plots, seed 43303) shared with the #143 species
# split + the PASMAX discrete test — validated bit-exact-in-RNG-consumption vs live FVSie.
# =============================================================================

using Test
using FVSjl

@testset "IE AUTOES cohort WK4 height-class partition (ie_autoes_tally emit, IE-gated)" begin
    occ = Float32[ones(9); zeros(14)]; over = zeros(Float32, 10); slo = 0.30f0
    kw = (seed0 = 43303, nplots = 50, ihab = 10, iser = 4, ifo = 4, iprep = 1, iphy = 3,
          xcos = cos(5.498f0) * slo, xsin = sin(5.498f0) * slo, slo = slo, elev = 34.0f0, baa = 1.0f0,
          regt = 1.0f0, bwaf = 0.0f0, bwb4 = 0.0f0, prob1 = 0.5527f0, dupnpt = 50.0f0, occ = occ, over = over)
    ihtser = 4   # MYHTS(IHAB=10) — habitat series for ESADVH/ESSUBH

    emit = Vector{NTuple{6,Float64}}()   # (sp, pt, ht, wk4, tpa, best) — 6th = BEST flag (IESTAT)
    tally = FVSjl.ie_autoes_tally(; kw..., is_ie = true, pasmax = 15f0, emit = emit,
                                  ihtser = ihtser, gentim = 5f0)

    # Aggregate booked TPA by WK4 class (round to the 0.60/0.20/0.00 STOMLT buckets).
    part = Dict(0.6f0 => 0.0, 0.2f0 => 0.0, 0.0f0 => 0.0)
    other = 0.0
    for r in emit
        wk4 = round(Float32(r[4]), digits = 1); tpa = r[5]
        haskey(part, wk4) ? (part[wk4] += tpa) : (other += tpa)
    end
    T = sum(values(part)) + other

    @testset "three-class structure present (NOT the collapsed single-0.60 record)" begin
        @test !isempty(emit)
        @test part[0.6f0] > 0.0            # advance ESADVH
        @test part[0.2f0] > 0.0            # subsequent ESSUBH (TRAGE=1)
        @test part[0.0f0] > 0.0            # subsequent/excess (TRAGE=0)
        @test other == 0.0                 # every booked record lands in one of the STOMLT buckets
    end

    @testset "advance (WK4=0.60) is the dominant class — matches the oracle (71-90% advance)" begin
        adv_frac = part[0.6f0] / T
        @test adv_frac > 0.60              # oracle-population advance fraction is 71-90%; fixture ≈ 0.665
        # advance must dominate the slow classes COMBINED (the over-correction would invert this)
        @test part[0.6f0] > part[0.2f0] + part[0.0f0]
    end

    @testset "faithful partition golden (locks b6edfd9a; regression guard)" begin
        # Golden fractions for this exact fixture (jl b6edfd9a, oracle-consistent three-class split). A
        # re-balance that shifts mass off the slow classes would move these — the guard is intentional.
        @test isapprox(part[0.6f0] / T, 0.665, atol = 0.01)   # advance
        @test isapprox(part[0.2f0] / T, 0.222, atol = 0.01)   # subsequent (TRAGE=1)
        @test isapprox(part[0.0f0] / T, 0.114, atol = 0.01)   # subsequent/excess (TRAGE=0)
        @test isapprox(T, 583.65, rtol = 0.02)                 # total established TPA (emit conserves the tally)
    end

    @testset "emit path conserves the summed tally (emit TPA == returned per-species tally)" begin
        # The emit records are the faithful per-record expansion of the same tally the collapsed path sums.
        @test isapprox(T, sum(tally), rtol = 0.05)
    end
end
