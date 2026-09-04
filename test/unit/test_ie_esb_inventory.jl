# IE/EM AUTOES ESB inventory-calibration (D1 post-thin re-stocking fix).
#
# D1 — the dominant post-disturbance tree-COUNT partition bug (EM/IE, Wykoff-general): after a cycle-2 thin
# opens the canopy, the AUTOES re-stocking cohort is mis-PRODUCED (BA converges, TPA/QMD diverge, seed-invariant).
# ROOT CAUSE (MEASURED vs FVS{ie,em}_g16 estab.f):
#   (1) ESB1 = ESTOCK(BAAOLD) must read the INVENTORY per-point OVERSTORY BA (ESFLTR BAAINV, base/fvs.f:201,
#       frozen at setup — never recomputed after a harvest), NOT the CURRENT post-thin per-point BA. jl formerly
#       used the thinned BA (e.g. 40 vs inventory-overstory 133) ⇒ ESB1 +2.4 instead of −0.97 ⇒ esb_shift
#       wrong-SIGNED ⇒ PROB1 0.18 vs oracle 0.87 ⇒ ~10× under-production of the disturbance re-stocking cohort.
#   (2) The disturbance CONTINUATION tally (NTALLY=2) must REUSE the stored ESB−ESB1 (estab.f:264/319/511 skip
#       the reset+recompute for NTALLY≠1, so estab.f:579 keeps applying the persisted ESB1(NCOUNT)). jl formerly
#       zeroed it (the `_ntally==1||99` gate excluded 2) ⇒ the continuation lost the stocking shift.
#
# MEASURED oracle values (FVSie_g16, IE stand 2978686010690, cyc-2 disturbance tally NTALLY=1):
#   TPACRE=1703 → ESB=1.160 ; BAAOLD=BAAINV=132.74 (inventory overstory) → ESB1=−0.965 ; esb_shift=+2.125
#   ⇒ PROB1 = logistic(PN=−0.260 + 2.125) = 0.866.  The post-thin BA (≈40) would give ESB1≈+2.4 ⇒ esb_shift≈−1.2
#   ⇒ PROB1≈0.18 (the pre-fix behaviour). Aggregate A/B (32 IE stands, THINBBA regime): sum|last-cycle TPA gap|
#   3126→2541 (10 improved / 4 regressed); EM unchanged (no regression).
using Test
using FVSjl

@testset "IE AUTOES ESB1 BAAOLD sensitivity (D1 esb_shift correction)" begin
    # ESB1 = ESTOCK(BAAOLD, TIME=0). Only the BA-carrying habitat series (IEQ=3 cedar/hemlock, e.g. IHAB=9)
    # make ESB1 depend on BAAOLD — exactly the stands where the inventory-vs-post-thin BA choice moves the answer.
    # Measured ie_estock (IHAB=9, IFO=4, ELEV=58.5): BAAOLD=132.74 → −2.460 ; BAAOLD=40 → −2.600.
    esb1_inv  = FVSjl.ie_estock(9, 1, 0f0, 1f0, 0f0, 58.5f0, 132.74f0, log(132.74f0), 0f0, 0f0, 0f0, 0f0, 4)
    esb1_thin = FVSjl.ie_estock(9, 1, 0f0, 1f0, 0f0, 58.5f0, 40f0,     log(40f0),     0f0, 0f0, 0f0, 0f0, 4)
    # The inventory (higher) overstory BA yields a HIGHER ESB1 than the thinned BA. Since esb_shift = ESB − ESB1,
    # the inventory BAAOLD produces a SMALLER shift ⇒ LOWER re-stocking PROB1 — the correction that pulls the
    # post-thin over-production down toward the oracle (A/B: 11802509 974→678 [oracle 700]; 3115334 250→154 [162]).
    @test esb1_inv > esb1_thin                    # higher inventory BA ⇒ higher ESB1
    # For a grand-fir series (IEQ=2, IHAB=4) ESB1 carries NO BA term — the BAAOLD choice is inert there, so those
    # stands are untouched by this fix (their D1 residual is the separate advance-regen/TPACRE cohort issue).
    a = FVSjl.ie_estock(4, 1, 0f0, 1f0, 0f0, 0f0, 132.74f0, log(132.74f0), 0f0, 0f0, 0f0, 0f0, 10)
    b = FVSjl.ie_estock(4, 1, 0f0, 1f0, 0f0, 0f0, 40f0,     log(40f0),     0f0, 0f0, 0f0, 0f0, 10)
    @test a == b
end

@testset "snapshot_esb_inputs! freezes inventory overstory BAAINV (ESFLTR)" begin
    # A synthetic IE stand: one small-tree (D<REGNBK) record + two overstory records on point 1. snapshot_esb_inputs!
    # must accumulate ONLY the D≥REGNBK point-1 basal area with the PTBAA scale (= point_basal_area!), and leave a
    # non-IE/EM variant untouched.
    s = FVSjl.StandState(FVSjl.InlandEmpire())
    t = s.trees
    t.n = 3
    t.dbh[1] = 1.0f0;  t.tpa[1] = 500f0; t.plot_id[1] = Int32(1)   # small tree — excluded from BAAINV
    t.dbh[2] = 10.0f0; t.tpa[2] = 30f0;  t.plot_id[2] = Int32(1)   # overstory, point 1
    t.dbh[3] = 12.0f0; t.tpa[3] = 20f0;  t.plot_id[3] = Int32(1)   # overstory, point 1
    s.plot.pi = 1f0; s.plot.gross_space = 1f0                      # PTBAA scale = 1
    @test isnan(s.estab.inv_baaold)
    FVSjl.snapshot_esb_inputs!(s)
    expected = 30f0 * 0.005454154f0 * 100f0 + 20f0 * 0.005454154f0 * 144f0   # D≥2.999 point-1 BA only
    @test isapprox(s.estab.inv_baaold, expected; rtol = 1f-5)
    # Idempotent: a second call must not change the frozen value.
    prev = s.estab.inv_baaold
    t.dbh[2] = 20f0                                                # would change BA if recomputed
    FVSjl.snapshot_esb_inputs!(s)
    @test s.estab.inv_baaold == prev

    # Non-IE/EM variant: no-op (field stays NaN).
    s2 = FVSjl.StandState(FVSjl.CentralRockies())
    s2.trees.n = 1; s2.trees.dbh[1] = 10f0; s2.trees.tpa[1] = 10f0; s2.trees.plot_id[1] = Int32(1)
    FVSjl.snapshot_esb_inputs!(s2)
    @test isnan(s2.estab.inv_baaold)
end
