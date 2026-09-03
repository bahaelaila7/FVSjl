# test_tt_mortality.jl — TT (Teton) density-mortality D10-recalibration loop (tt/morts.f label 10), 2026-09-03.
#
# tt/morts.f wraps the whole TN10/RN/kill/TTMRT block in a QMD-convergence loop (label 10, IPASS≤10): the
# selective (percentile) self-thinning raises the post-kill Reineke self-thin diameter D10N above the
# uniform-growth D10 estimate, which tightens the SDI density target (TMD10=CONST·D10^−1.605) and kills MORE.
# FVS re-runs with D10←D10N until |D10−D10N|≤0.1 (or D10N≤DIA0, or IPASS=10). The port UNDER-self-thinned dense
# conifer stands because it stopped at the first pass (both-sides g16-traced CN 388908802489998 cyc1: oracle
# iterates 2 passes 2967→2682; old jl held at the pass-1 value 2967). The sibling Zeide-SDI variant UT already
# ports the identical loop (test_ut_mortality.jl); this guards the TT copy.
#
# Self-contained (no DB / no oracle): a dense DF cohort well above the SDImax 85% line, with a wide diameter
# spread so the percentile kill lifts QMD past the 0.1 convergence threshold ⇒ the loop takes >1 pass.
using Test, FVSjl

@testset "TT (Teton) morts D10-recalibration self-thin loop" begin
    # Build a dense single-species (DF, sp3) stand: 30 records 1.0–6.0", 250 tpa each (7500 tpa), 0.5" DG,
    # crown_ratio (=FVS PCT percentile) spread 0–100, SDImax seeded to 300 so tt ≫ the 85% self-thin line.
    build = function ()
        s = FVSjl.StandState(FVSjl.Teton()); t = s.trees; p = s.plot; n = 30
        for i in 1:n
            t.species[i]     = Int32(3)
            t.dbh[i]         = Float32(1.0 + 5.0 * (i - 1) / (n - 1))
            t.tpa[i]         = 250f0
            t.diam_growth[i] = 0.5f0
            t.crown_ratio[i] = Float32((i - 1) * 100 / (n - 1))   # PCT (BA percentile) — TTMRT allocator input
            t.crown_pct[i]   = Int32(40)
            t.dg_prev[i]     = 0f0
        end
        t.n = n
        p.sp_sdi_def[3] = 300f0                                    # DF SDImax (stand_sdimax uses sp_sdi_def)
        ba = 0f0; for i in 1:n; ba += 0.0054542f0 * t.dbh[i]^2 * t.tpa[i]; end
        p.basal_area = ba
        s
    end

    s = build(); t = s.trees
    tpa0 = sum(t.tpa[i] for i in 1:t.n)
    @test tpa0 ≈ 7500f0
    FVSjl.mortality!(s, FVSjl.Teton(); book_snags = false)
    residual = sum(t.tpa[i] for i in 1:t.n)

    # FIRES — the dense over-max-SDI stand self-thins hard (>90% killed).
    @test (tpa0 - residual) / tpa0 > 0.9

    # ITERATES — the D10-recalibration loop drives the residual BELOW the first-pass (single-pass) target.
    # Measured: with the loop residual ≈ 528.8 tpa; a single pass (loop removed) leaves 635.8 tpa. The loop's
    # extra ~107 tpa of kill is the whole point of the fix, so residual must sit well under the single-pass value.
    @test residual < 600f0                       # single-pass would leave 635.8 — this FAILS without the loop
    @test isapprox(residual, 528.8f0; atol = 3f0)   # pinned multi-pass residual (regression guard)
end
