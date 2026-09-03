# =============================================================================
# test_op_native_mortality_rip.jl — OP (Olympic) FVS-native no-big-6 RIP mortality
# (op/morts.f BM0..BM5 logistic + Gould–Harrington small-tree model) vs live FVSop_clean.
#
# When a stand has NO ORGANON big-6 record (no GF sp3 / DF sp16), ORGANON never runs
# (op_org_ran=false) and op/morts.f falls back to its FVS-native per-species RIP mortality.
# The validated cooperating stand S248112 (op2c) HAS a DF ⇒ it exercises the MORTEXP-for-all
# path; this test swaps DF/GF → WF in that stand (opnb6.tre) so ORGANON is disabled and the
# RIP path fires. Oracle FVSop_clean shows RIP declining TPA 536→344 over 10 cycles.
#
# Reference: test/harness/scenarios/opnb6.sum.save (FVSop_clean, opnb6.key/opnb6.tre, DGSD=0
# ⇒ deterministic ⇒ a TRUE bit-exact bar on the mortality trajectory). VALIDATED BIT-EXACT on
# TPA / BA / SDI / accretion-mortality growth every cycle. The only residuals are cubic/board-
# foot VOLUME columns from cycle 1+ (the pre-existing OC/OP ORGANON GMV merch-cubic volume
# corner — identical surviving population, so NOT the mortality port) and one QMD rounding
# straddle (2020: 9.9 vs 10.0); those columns are not asserted here.
#
# Guards the RIP path from silently breaking: without the port, mortality!(::Olympic) on a
# no-big-6 stand used to `error(...)`; a regression to under/over-kill moves the TPA trajectory.
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "OP no-big-6 RIP mortality vs FVSop_clean (opnb6, NUMCYCLE=10)" begin
    key = joinpath(@__DIR__, "..", "harness", "scenarios", "opnb6.key")
    if !isfile(key) || !isfile(joinpath(@__DIR__, "..", "harness", "scenarios", "opnb6.tre"))
        @test_skip "opnb6 scenario not available"
    else
        rows = F.SummaryRow[]
        for s in F.each_stand(key; variant = F.Olympic())
            F.notre!(s); F.setup_growth!(s); F.compute_volumes!(s)
            # This stand has NO big-6 (all WF/LP/PP/SP): ORGANON must be OFF ⇒ RIP path taken.
            @test s.calib.op_org_ran == false
            push!(rows, F.summary_row(s))
            for _ in 1:Int(s.control.ncycle)
                F.grow_cycle!(s; fint = 5f0)
                F.compute_volumes!(s)
                push!(rows, F.summary_row(s))
            end
            break
        end
        @test length(rows) == 11

        # Oracle TPA/BA/SDI trajectory (RIP firing) — bit-exact every cycle.
        oracle_tpa = [536, 511, 488, 466, 445, 425, 406, 389, 373, 358, 344]
        oracle_ba  = [77,  97,  120, 144, 169, 194, 219, 244, 268, 291, 313]
        oracle_sdi = [184, 220, 257, 295, 332, 369, 403, 435, 465, 493, 519]
        oracle_yr  = [1990,1995,2000,2005,2010,2015,2020,2025,2030,2035,2040]
        for (k, r) in enumerate(rows)
            @test r.year == oracle_yr[k]
            @test r.tpa  == oracle_tpa[k]
            @test r.ba   == oracle_ba[k]
            @test r.sdi  == oracle_sdi[k]
        end

        # The RIP path must actually REMOVE trees (guard against a no-op regression):
        # cumulative kill of ~192 TPA (536→344) over 10 cycles.
        @test rows[1].tpa - rows[end].tpa == 192
    end
end
