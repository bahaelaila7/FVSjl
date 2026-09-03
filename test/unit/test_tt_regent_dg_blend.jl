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
end
