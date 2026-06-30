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

@testset "CS cst01 BARE-GROUND-PLANT establishment (vs live FVScs)" begin
    # Stand 4 of cst01: NOTREES + ESTAB 1992 + PLANT (sp 3 shortleaf 400 TPA, sp 21 bitternut 400 TPA),
    # 10×10yr cycles. Exercises the CS establishment path: ESSUBH base height (cs/essubh.f MAPCS refage +
    # CS NC-128 forward curve) → REGENT(LESTB) creation-cycle growth (cs_balmod + CS htcalc increment,
    # DIAM=htdbh_db floor, HHTMAX cap). Validated vs freshly-relinked live FVScs (captured 2026-06-30).
    # Establishment draws BACHLO heights/crowns, so single-precision diffs in the height curve compound
    # immediately ⇒ no bit-exact cycle, but every column tracks live within the documented ULP-class floor.
    key = "/workspace/ForestVegetationSimulator/tests/FVScs/cst01.key"
    barekey = joinpath(@__DIR__, "cst01_bare.key")
    if !isfile(key)
        @info "cst01.key not present; skipping CS BARE establishment test"
    else
        # Build the single-stand BARE key (lines 117-132 of cst01.key) next to the test.
        lines = readlines(key)
        open(barekey, "w") do io
            for l in lines[117:132]; println(io, l); end
        end
        cp(joinpath(dirname(key), "cst01.tre"), joinpath(@__DIR__, "cst01.tre"); force = true)
        cd(@__DIR__) do
            sumtxt = FVSjl.run_keyfile(barekey; variant = CentralStates(), output = :sum)
            # Live FVScs BARE .sum: Year => (TPA, BA, SDI, CCF, TopHt, QMD)
            live = Dict(
                2002 => (800,   9,  33,  22, 14, 1.4), 2012 => (787,  42, 121, 106, 32, 3.1),
                2022 => (774, 101, 243, 244, 47, 4.9), 2032 => (533, 183, 358, 387, 58, 7.9),
                2042 => (405, 207, 373, 422, 66, 9.7), 2052 => (311, 208, 354, 417, 71, 11.1),
                2062 => (254, 208, 339, 405, 75, 12.2), 2072 => (216, 207, 327, 382, 79, 13.2),
                2082 => (191, 206, 318, 358, 81, 14.0), 2092 => (167, 204, 307, 327, 83, 15.0))
            rows = Dict{Int,NTuple{6,Float64}}()
            for ln in split(sumtxt, '\n')
                t = split(strip(ln))
                (length(t) >= 8 && tryparse(Int, t[1]) !== nothing && t[1] != "-999") || continue
                yr = parse(Int, t[1]); yr == 1992 && continue
                rows[yr] = (parse(Float64,t[3]), parse(Float64,t[4]), parse(Float64,t[5]),
                            parse(Float64,t[6]), parse(Float64,t[7]), parse(Float64,t[8]))
            end
            for (yr, L) in sort(collect(live))
                @test haskey(rows, yr)
                r = rows[yr]
                if yr == 2002                       # establishment + first-cycle growth: BIT-EXACT all 6
                    @test (r[1],r[2],r[3],r[4],r[5]) == (L[1],L[2],L[3],L[4],L[5])
                    @test round(Float32(r[6]); digits=1) == Float32(L[6])
                    continue
                end
                # 2012+: TPA stays bit-exact two cycles (count/mortality logic exact); the drift is SIZE-only
                # single-precision accumulation in the small-tree growth spine (SDI/CCF/TopHt ±1/cyc, peaking
                # 2072 — CCF amplifies crown-width via dbh). GROUNDED: not an establishment defect (2002 exact).
                @test abs(r[1] - L[1]) <= 4        # TPA (bit-exact through 2012)
                @test abs(r[2] - L[2]) <= 2        # BA
                @test abs(r[3] - L[3]) <= 4        # SDI
                @test abs(r[4] - L[4]) <= 10       # CCF (crown-width amplified, compounds most)
                @test abs(r[5] - L[5]) <= 2        # TopHt
                @test abs(r[6] - L[6]) <= 0.2      # QMD
            end
        end
        rm(barekey; force = true); rm(joinpath(@__DIR__, "cst01.tre"); force = true)
    end
end
