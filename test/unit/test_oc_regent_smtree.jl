# =============================================================================
# test_oc_regent_smtree.jl — OC REGENT small-tree height model (smhtgf.f) vs live FVSoc_clean.
#
# OC/OP routed ALL IORG=0 trees to the large-tree oc_htgf_native, which drops the RELHT competition
# suppression when PCCF<100 and so over-grew sub-breast-height seedlings ~7ft/cycle (tn2 DF: oracle
# 2.0→3.8 vs jl 2.0→10.7). FVS routes trees with DBH<XMAX(sp) to REGENT/SMHTGF instead. OC/OP are
# DGSD=0 ⇒ the small-tree path is deterministic ⇒ bit-exact.
#
# Ground truth = the live FVSoc_clean scoped-DEBUG REGENT dump (oct01 stand-1 cyc0, scratchpad/oc/):
# DF (firs eq, CON=1.0) HTGRR=HTGR=1.8110 for (D=0.1, H=2.0, ICR=55, SI=92, RELHT=H/AVH). Validated
# 2026-08-20. HCOR (small-tree height calibration) is a follow-on (only GF calibrated here, ~+1.9%).
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "OC REGENT smhtgf — DF firs equation bit-exact vs live dump" begin
    # DF = sp 7 (firs, MAPSP=2). Inputs from the live dump / jl setup for oct01 stand-1 tree 2.
    si = 92.0f0
    h = 2.0f0
    avh = 63.43875f0
    relht = min(h / avh, 1.05f0)
    cr_regent = 55.0f0 / 10.0f0                 # ICR/10, the regent CR scale
    bal = 0.0f0                                 # firs eq does not use BAL
    htgr = F.oc_smhtgf(7, 0.1f0, h, cr_regent, 184.0f0, bal, si, relht)
    @test isapprox(htgr, 1.8110f0; atol = 1.0f-3)   # live REGENT dump HTGRR

    # firs branch is BAL-independent (only SI/CR/RELHT); GF (sp 4) is also firs.
    @test F.oc_smhtgf(4, 0.1f0, h, cr_regent, 184.0f0, 0f0, si, relht) ≈ htgr atol=1f-4

    # XMAX routing thresholds (regent.f DATA): 4.0" default, 10.0 for GS(23)/RW(50).
    @test F.OC_REG_XMAX[7] == 4.0f0
    @test F.OC_REG_XMAX[23] == 10.0f0 && F.OC_REG_XMAX[50] == 10.0f0
end

@testset "OC REGENT — small-tree height end-to-end (heights match live, was +7ft off)" begin
    fx = joinpath(@__DIR__, "..", "fixtures", "oregoncoast")
    key = joinpath(fx, "ocgro.key")
    if !isfile(key)
        @test_skip "oregoncoast small-tree fixture not present"
    else
        state = cd(fx) do
            snaps = Dict{Int,Tuple{Float32,Float32}}()
            for s in F.each_stand("ocgro.key"; variant = F.variant_from_code("OC"))
                F.notre!(s); F.setup_growth!(s); F.compute_volumes!(s)
                hook = (st, yr, pl, cy) -> (cy == 1 && foreach(
                    i -> (snaps[Int(st.trees.tree_id[i])] = (st.trees.dbh[i], st.trees.height[i])), 1:st.trees.n))
                io = IOBuffer(); F.write_sum_file(io, s; variant = "OC", cycle_hook = hook)
                break
            end
            snaps
        end
        # HEIGHT (smhtgf): tn2 (DF seedling) live cyc1 h=3.80 (bug gave 10.7); tn15 (GF) h=5.10 (was 12.1).
        @test haskey(state, 2)  && abs(state[2][2]  - 3.80f0) < 0.15f0
        @test haskey(state, 15) && abs(state[15][2] - 5.10f0) < 0.15f0
        # DBH-from-height (HTDBH Curtis/Arney): tn15 GF cyc1 live d=0.5 (the bug gave d=2.3 via /bark;
        # the raw-HT1 formula gave 1.3). Must track the oracle within a NINT.
        @test haskey(state, 15) && abs(state[15][1] - 0.5f0) < 0.15f0
    end
end
