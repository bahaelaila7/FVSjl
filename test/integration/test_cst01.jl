# =============================================================================
# test_cst01.jl — Central States (CS) variant, cst01 cycle-0 stand columns
#
# Validates the CS variant's inventory-cycle (1990) stand statistics against the
# live FVScs .sum (tests/FVScs/cst01.sum). All six geometric columns are bit-exact
# (TPA/BA/SDI/CCF reported per-gross-acre = stockable ÷ GROSPC; QMD/TopHt unscaled).
# This stand exercises the GROSPC<1 path (11 plots, 1 non-stockable ⇒ GROSPC 0.909)
# that the SN/NE test stands never hit — the .sum writer's ÷GROSPC scale-back.
# =============================================================================

using Test
using FVSjl

@testset "CS cst01 cycle-0 stand columns (vs live FVScs)" begin
    key = "/workspace/ForestVegetationSimulator/tests/FVScs/cst01.key"
    if !isfile(key)
        @info "cst01.key not present; skipping CS cycle-0 test"
    else
        s, reason = FVSjl.initialize(key; variant = CentralStates())
        @test reason == :process
        @test s.trees.n == 27
        FVSjl.notre!(s)                    # NOTRE — BAF/fixed-plot expansion × GROSPC
        FVSjl.setup_growth!(s)             # CRATET dub + cs_dgcons! bark copy (CFTOPK broken-top)
        FVSjl.compute_volumes!(s)          # eastern R9 Clark cubic + R9LOGS board feet, CS merch
        FVSjl.compute_forest_type!(s)      # FORTYP/STKVAL (shared, FIA-keyed; CS reuses NE stocking CSVs)
        r = FVSjl.summary_row(s; period = 0)

        # Live FVScs cst01.sum, 1990 inventory row (per-gross-acre):
        @test r.tpa   == 536
        @test r.ba    == 77
        @test r.sdi   == 160
        @test r.ccf   == 169
        @test r.topht == 63
        @test round(r.qmd, digits = 1) == 5.1

        # Volume columns (per-gross-acre) — all bit-exact. Sawtimber-cubic + board feet ride the
        # CS merch standards (_cs_merch); total/merch cubic depend on the CFTOPK broken-top path
        # (cs_dgcons! bark copy — without it the two broken-top trees read ~4 cuft/acre low).
        @test r.cuft  == 1517
        @test r.mcuft == 1300
        @test r.scuft == 497
        @test r.bdft  == 2903

        # Trailing classification fields (FORTYP/STKVAL): 503 = W.OAK-R.OAK-HICKORY, size/stock
        # class 22. The .sum growth columns (period/accretion/mortality) need the cycle 0→1 DG
        # projection (cs/dgf.f, chunk 3) and are intentionally not asserted here.
        @test r.fortype  == 503
        @test r.sizecls  == 2
        @test r.stockcls == 2
    end
end

@testset "CS cst01 cycle-1 growth (vs live FVScs)" begin
    key = "/workspace/ForestVegetationSimulator/tests/FVScs/cst01.key"
    if !isfile(key)
        @info "cst01.key not present; skipping CS cycle-1 test"
    else
        # Full growth spine end-to-end: cs_dgf! DG + calibration, cs htgf height growth, cs regent
        # small-tree growth, TWIGS crown, varmrt + background mortality, R9 volume. The DG calibration
        # uses the CS GST DBH≥5 floor (cs/dgdriv.f:380) — with the SN/NE ≥3 floor, WO would spuriously
        # calibrate (debug-stamped: live FN<5 ⇒ COR=0 for all species) and over-grow, inflating BA/SDI.
        s = first(FVSjl.each_stand(key; variant = CentralStates()))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s); FVSjl.compute_forest_type!(s)
        g = s.plot.gross_space; di(x) = trunc(Int, x + 0.5)
        FVSjl.grow_cycle!(s; fint = 10f0)
        # Live FVScs cst01.sum 2000 row — all six stand columns BIT-EXACT:
        @test di(FVSjl.stand_tpa(s) / g)        == 518
        @test di(FVSjl.stand_ba(s) / g)         == 99
        @test di(FVSjl.stand_sdi(s) / g)        == 196
        @test di(FVSjl.stand_ccf(s) / g)        == 202
        @test di(FVSjl.stand_top_height(s))     == 68
        @test round(FVSjl.stand_qmd(s); digits = 1) == 5.9f0

        # cycle-1 grown-stand volume (live 2000: Tcuft 2110 / Mcuft 1887 / Scuft 886 / Bdft 5084).
        # Tcuft + Scuft are bit-exact; Mcuft/Bdft land within 1 unit (~0.05% grown-stand rounding floor).
        FVSjl.compute_volumes!(s)
        v = FVSjl.summary_row(s; period = 0)
        @test v.cuft  == 2110
        @test v.scuft == 886
        @test abs(v.mcuft - 1887) <= 1
        @test abs(v.bdft  - 5084) <= 1
    end
end

@testset "CS cst01 multi-cycle projection (vs live FVScs, 10 cycles)" begin
    key = "/workspace/ForestVegetationSimulator/tests/FVScs/cst01.key"
    if !isfile(key)
        @info "cst01.key not present; skipping CS multi-cycle test"
    else
        # First stand (S248112 UNTHINNED CONTROL, NOAUTOES) projected 10×10yr cycles. Validates the
        # full growth spine PAST the bit-exact cycles 0–2 into the single-precision (Float32 ULP) tail.
        # Oracle = live FVScs relinked from bin/FVScs_buildDir (captured 2026-06-30). Stand columns are
        # bit-exact through cycle 2 (1990/2000/2010); cycles 3–10 accumulate a single-precision drift
        # (TPA ±3, SDI ±1, QMD ±0.1, ~1–1.5% on board feet which amplifies DBH steps via Scribner) —
        # the accepted SN-COMPRESS-class ULP floor (RNG S0 bit-aligned every cycle; logic verified).
        # Live (per-gross-acre): {Year => (TPA,BA,SDI,CCF,TopHt,QMD)}
        live = Dict(
            1990 => (536,  77, 160, 169, 63, 5.1), 2000 => (518,  99, 196, 202, 68, 5.9),
            2010 => (476, 122, 231, 234, 70, 6.9), 2020 => (440, 145, 264, 264, 75, 7.8),
            2030 => (310, 147, 255, 253, 80, 9.3), 2040 => (233, 145, 240, 237, 83, 10.7),
            2050 => (179, 143, 227, 223, 86, 12.1), 2060 => (136, 142, 216, 208, 89, 13.8),
            2070 => (110, 142, 209, 200, 91, 15.4), 2080 => ( 93, 143, 205, 195, 95, 16.8),
            2090 => ( 81, 144, 202, 194, 99, 18.0))
        s = first(FVSjl.each_stand(key; variant = CentralStates()))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        g = s.plot.gross_space; di(x) = trunc(Int, x + 0.5)
        for yr in 1990:10:2090
            tpa, ba, sdi, ccf, topht, qmd =
                di(FVSjl.stand_tpa(s)/g), di(FVSjl.stand_ba(s)/g), di(FVSjl.stand_sdi(s)/g),
                di(FVSjl.stand_ccf(s)/g), di(FVSjl.stand_top_height(s)), round(FVSjl.stand_qmd(s); digits=1)
            L = live[yr]
            if yr <= 2010                                   # cycles 0–2: BIT-EXACT
                @test (tpa, ba, sdi, ccf, topht) == (L[1], L[2], L[3], L[4], L[5])
                @test round(Float32(qmd); digits=1) == Float32(L[6])
            else                                            # cycles 3–10: Float32 ULP floor
                @test abs(tpa   - L[1]) <= 3
                @test abs(ba    - L[2]) <= 1
                @test abs(sdi   - L[3]) <= 1
                @test abs(ccf   - L[4]) <= 2
                @test abs(topht - L[5]) <= 2
                @test abs(qmd   - L[6]) <= 0.15
            end
            yr < 2090 && FVSjl.grow_cycle!(s; fint = 10f0)
        end
    end
end
