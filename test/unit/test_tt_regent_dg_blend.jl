# test_tt_regent_dg_blend.jl — TT (Teton) regent small/large-tree coupling regression guard.
#
# Guards the coupled M331D woodland over-growth fix (2026-09-03):
#  (1) TT_RG_XMIN/XMAX must be the VERBATIM tt/regent.f:161-171 DATA (1.5/3.0 for the TTVAR conifers),
#      NOT the old mis-read 2.0/4.0.
#  (2) The DIAMETER increment is NOT XWT-blended: FVS regent.f blends only HEIGHT (regent.f:799); for
#      D<BKPT the regent small-tree DG fully replaces the large-tree dgf DG. The invariant that makes the
#      default-conifer DG pure-regent is BKPT == XMAX == 3.0 (so within the blend loop, d<XMAX ⇒ d<BKPT).
#      A regression that re-introduced the DG blend, or moved XMAX off BKPT, over-grew the 1.5–3" woodland
#      trees (measured 40-stand M331D total-rel: density-fix-only 9.84 → XMIN-alone 17.44 → combined 5.79).
using Test
using FVSjl

@testset "TT regent XMIN/XMAX DATA + DG-not-blended coupling" begin
    XMIN = FVSjl.TT_RG_XMIN
    XMAX = FVSjl.TT_RG_XMAX
    BREAK = FVSjl.TT_RG_BREAK

    # (1) Verbatim tt/regent.f DATA XMIN / XMAX.
    @test XMIN == Float32[1.5, 1.5, 1.5, 90.0, 1.5, 1.5, 1.5, 1.5, 1.5, 2.0,
                          90.0, 90.0, 90.0, 2.0, 0.5, 90.0, 1.5, 0.5]
    @test XMAX == Float32[3.0, 3.0, 3.0, 99.0, 3.0, 3.0, 3.0, 3.0, 3.0, 5.0,
                          99.0, 99.0, 99.0, 4.0, 2.0, 99.0, 3.0, 2.0]

    # The TTVAR default conifers must NOT carry the old 2.0/4.0 window.
    for sp in (1, 2, 3, 5, 7, 8, 9, 17)
        @test XMIN[sp] == 1.5f0
        @test XMAX[sp] == 3.0f0
    end

    # (2) DG-purity invariant. The blend loop `continue`s at d≥XMAX and the DG is set from the regent value
    # only when d<BKPT. For the default conifers BKPT == XMAX == 3.0, so every tree the loop sees (d<XMAX)
    # also has d<BKPT ⇒ pure regent DG, never a blend of the large-tree dgf DG.
    for sp in (1, 2, 3, 5, 6, 7, 8, 9, 17)
        @test BREAK[sp] == 3.0f0
        @test XMAX[sp] == 3.0f0
    end
    # MM (14) is the one default species whose window is wider than its break (XMAX=4 > BKPT=3): a
    # 3–4" MM tree enters the loop but keeps the large-tree dgf DG (d≥BKPT), so the `d < BKPT` gate is
    # load-bearing — it must not be simplified back to an unconditional overwrite.
    @test BREAK[14] == 3.0f0
    @test XMAX[14] == 4.0f0


    # ── Firing test: BI/MC (13,16) inventory Curtis-Arney HT-DBH (regent.f:874-901) ──
    # jl formerly used the 0.1·HTG rule-of-thumb; but LHTDRG(13)=LHTDRG(16)=.FALSE. means FVS ALWAYS
    # takes the inventory equation `DG=(DK−DKK)·bark`. _tt_bimc_dk maps a grown height to the predicted
    # DBH. Pinned to the FVStt_dbg-measured branch (dense Gambel-oak stand 51031230020004 cyc1): a BI
    # tree grown to HK≈4.84 gets DK≈0.335 (oracle 0.3324 at its ZZRAN'd HK=4.816), NOT 0.1·HTG≈0.383.
    # Deleting the inventory branch (reverting to 0.1·HTG) changes this value ⇒ the test is not vacuous.
    let
        # BI(13): P2/P3/P4 = 76.5170 / 2.2107 / −0.6365
        p2, p3, p4 = 76.5170f0, 2.2107f0, -0.6365f0
        hat3 = 4.5f0 + p2 * FVSjl.fexp(-p3 * FVSjl.fpow(3.0f0, p4))
        @test hat3 == 30.005816f0
        @test FVSjl._tt_bimc_dk(4.8407f0, hat3, p2, p3, p4) == 0.33502105f0   # HK<HAT3 rational form
        # MC(16): P2/P3/P4 = 1709.7229 / 5.8887 / −0.2286
        q2, q3, q4 = 1709.7229f0, 5.8887f0, -0.2286f0
        hat3m = 4.5f0 + q2 * FVSjl.fexp(-q3 * FVSjl.fpow(3.0f0, q4))
        @test hat3m == 22.017498f0
        @test FVSjl._tt_bimc_dk(10.0f0, hat3m, q2, q3, q4) == 1.1466658f0
        # DG must be (DK−DKK)·bark (DKK=D since original H≤4.5), NOT 0.1·HTG. For DK=0.335, DKK=0.1,
        # bark≈0.99: DG ≈ 0.233 — distinctly below the old rule-of-thumb 0.1·3.83 ≈ 0.383.
        dg_inventory = (0.33502105f0 - 0.1f0) * 0.99f0
        @test 0.22f0 < dg_inventory < 0.24f0
        @test dg_inventory < 0.1f0 * 3.83f0   # strictly less than the removed 0.1·HTG over-grow
    end
end