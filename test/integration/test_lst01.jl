# =============================================================================
# test_lst01.jl — Lake States (LS) variant, lst01 cycle-0 stand columns
#
# Validates the LS variant's inventory-cycle (1990) stand statistics against a
# freshly-relinked live FVSls .sum (tests/FVSls/lst01.key, stand 248112). All six
# geometric columns are bit-exact (TPA/BA/SDI/CCF reported per-gross-acre =
# stockable ÷ GROSPC; QMD/TopHt unscaled). Same inventory stand as snt01/cst01,
# so the BAF/fixed-plot expansion × GROSPC path is exercised, driven by the LS
# 68-species data (roster/translation/coeffs), LS SITSET (68×68 SICOEF fan-out),
# LS htdbh (Wykoff+Curtis-Arney height dub), and the shared eastern crown-width
# library (LS crown_width_species map + NE equation coefficients) for CCF.
#
# Live FVSls lst01.sum 1990 row: 536 77 160 171 63 5.1  (vol 1551 1338 480 1887).
# Volume columns wait on CHUNK 2 (LS = R9 DVEE Gevorkiantz + Clark mix, not yet ported).
# =============================================================================

using Test
using FVSjl

@testset "LS lst01 cycle-0 stand columns (vs live FVSls)" begin
    key = "/workspace/ForestVegetationSimulator/tests/FVSls/lst01.key"
    if !isfile(key)
        @info "lst01.key not present; skipping LS cycle-0 test"
    else
        s, reason = FVSjl.initialize(key; variant = LakeStates())
        @test reason == :process
        @test s.trees.n == 27
        FVSjl.notre!(s)                    # NOTRE — BAF/fixed-plot expansion × GROSPC
        FVSjl.setup_growth!(s)             # CRATET dub (LS htdbh) + crown init
        FVSjl.compute_forest_type!(s)      # FORTYP/STKVAL (shared, FIA-keyed)
        r = FVSjl.summary_row(s; period = 0)

        # Live FVSls lst01.sum, 1990 inventory row (per-gross-acre) — all six BIT-EXACT:
        @test r.tpa   == 536
        @test r.ba    == 77
        @test r.sdi   == 160
        @test r.ccf   == 171
        @test r.topht == 63
        @test round(r.qmd, digits = 1) == 5.1

        # Volume (per-gross-acre). LS rides the shared eastern R9 Clark cubic (default METHC=6) but its
        # `.sum` BdFt is SCRIBNER (vol2, ls/vols.f:348), not the International ¼" (vol10) NE/CS report —
        # ported via `_r9_scribner_bf` (scribner_factor.csv). SCuFt + BdFt are BIT-EXACT vs live FVSls.
        FVSjl.compute_volumes!(s)
        v = FVSjl.summary_row(s; period = 0)
        @test v.scuft == 480
        @test v.bdft  == 1887
        # TCuFt 1551 / MCuFt 1338: jl 1546 / 1333 (Δ5, ~0.3% — open, tracked in the ledger).
    end
end

@testset "LS lst01 cycle-1 growth spine (vs live FVSls)" begin
    key = "/workspace/ForestVegetationSimulator/tests/FVSls/lst01.key"
    if !isfile(key)
        @info "lst01.key not present; skipping LS cycle-1 test"
    else
        # Full LS growth spine: ls_dgf! DG + shared DGDRIV calibration, LS htgf (NC-128 MAPLS + ls_balmod),
        # LS regent small-tree, TWIGS crown (LS BCR), varmrt + LS background mortality (PMSC/PMD via IMAPLS).
        s = first(FVSjl.each_stand(key; variant = LakeStates()))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s); FVSjl.compute_forest_type!(s)
        g = s.plot.gross_space; di(x) = trunc(Int, x + 0.5)
        FVSjl.grow_cycle!(s; fint = 10f0)
        # Live FVSls lst01.sum 2000 row: 524/104/203/210/64/6.0 — ALL columns now BIT-EXACT. The former
        # BA/SDI/CCF Δ1 was TWO small-tree bugs, both fixed: (1) the HTDBH mode-1 DB floor (htdbh.f:343,
        # missing for LS/NE/CS — a sugar-maple seedling got a negative Wykoff-inverse DBH → grew to 0.68″
        # vs the floored 0.20″); (2) the crown-ratio BA basis (crown.f uses RAW per-acre BA, jl used
        # basal_area/gross_space) which shifted 2 un-crown-capped trees' CR by 2% → their dgf DDS. With
        # both, lst01 base stand is bit-exact 1990-2040 (50 yrs); the 2050+ tail (TPA Δ1-5, BA/SDI exact)
        # is the tripling-spread DGSCOR (NOTRIPLE is bit-exact all 100 yrs).
        @test di(FVSjl.stand_tpa(s) / g)    == 524
        @test di(FVSjl.stand_top_height(s)) == 64
        @test round(FVSjl.stand_qmd(s); digits = 1) == 6.0f0
        @test di(FVSjl.stand_ba(s) / g)  == 104
        @test di(FVSjl.stand_sdi(s) / g) == 203
        @test di(FVSjl.stand_ccf(s) / g) == 210
    end
end
