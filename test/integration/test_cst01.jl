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
        # ALL FOUR are BIT-EXACT (re-measured: Mcuft=1887, Bdft=5084 exactly — the old ≤1 was over-cautious).
        FVSjl.compute_volumes!(s)
        v = FVSjl.summary_row(s; period = 0)
        @test v.cuft  == 2110
        @test v.scuft == 886
        @test v.mcuft == 1887   # BIT-EXACT (was over-cautious ≤1)
        @test v.bdft  == 5084   # BIT-EXACT (was over-cautious ≤1)
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
        # cycles 3–10 drift rows: (tpa,ba,sdi,ccf,topht,qmd, L1..L6) collected for the after-loop @test_broken
        drift = NTuple{12,Float64}[]
        for yr in 1990:10:2090
            tpa, ba, sdi, ccf, topht, qmd =
                di(FVSjl.stand_tpa(s)/g), di(FVSjl.stand_ba(s)/g), di(FVSjl.stand_sdi(s)/g),
                di(FVSjl.stand_ccf(s)/g), di(FVSjl.stand_top_height(s)), round(FVSjl.stand_qmd(s); digits=1)
            L = live[yr]
            if yr <= 2010                                   # cycles 0–2: BIT-EXACT
                @test (tpa, ba, sdi, ccf, topht) == (L[1], L[2], L[3], L[4], L[5])
                @test round(Float32(qmd); digits=1) == Float32(L[6])
            else                                            # cycles 3–10: CS height-path Float32 residual (much reduced)
                # ROOT (docs/TOLERANCE_AUDIT.md, CORRECTED 2026-07-05): NOT a "transcendental floor" — it was
                # two REAL op-order bugs in the CS HEIGHT path (balmod is height-only ⇒ DBH/BA stay bit-exact):
                #   (1) cs/balmod.f:67 `(1.-TEMBA/210.)**.5` — FVS `**0.5` lowers to powf, jl used sqrt (differ
                #       ~0.05%); fixed to fpow(x,0.5).  (2) htcalc.f:394 age inversion `(H-BH)/B1/SI**B2` +
                #       `1./B4/SI**B5` are SEQUENTIAL divisions ((a/b)/c), jl divided by the product a/(b·c);
                #       fixed to match. Together these took TopHt from drifting every grown cycle to Δ0 on
                #       2020–2060+2080 (residual Δ1@2070, Δ-2@2090). The remaining TPA Δ1 is a regen-tree VARMRT
                #       KNIFE-EDGE (BA BIT-EXACT every cycle ⇒ only a <1-TPA regen tree flips). CORNERED BY
                #       ELIMINATION (2026-07-05uuu corner-campaign): the regen height = htg1·(1+OLDRN)·gmod feeds
                #       relht=min(hti/AVH,1) → varmrt EFFTR → the discrete deletion that flips ONE <1-TPA regen tree.
                #       Every DECONFOUNDABLE factor was FFI-routed to gfortran AND tested INERT on this flip:
                #       (1) PEFF `RELHTA**3` (sn mortality.jl:138/ne:23 → fpow); (2) htg1 NC-128 SI**B2/EXP/ALOG/**EX
                #       (ne_htcalc_age/incr → fpow/fexp/flog); (3) gmod balmod EXP(-B3·BAL) (ne_balmod → fexp);
                #       (4) AVH — a +1 Float32 ULP bump to stand_top_height left TPA UNCHANGED. ★ LIVE-CONFIRMED
                #       (2026-07-06, FVS_TreeList DBS per-tree differential across cycles): the per-tree grown DBH
                #       SET is BIT-EXACT at cycle 1 (2000: 0/81 mismatch, max|Δ|=0), then diverges AT TRIPLING —
                #       2010 (cycle 2) 3/243 (max 0.009"), 2020 56/243 (max 0.027"). BA still rounds bit-exact
                #       because the sub-0.03" DBH diffs are small vs BA~145. NOTRIPLE isolates it FURTHER — it is
                #       NOT tripling-caused: under NOTRIPLE exactly ONE record diverges at cycle 2 (2010: 1/27,
                #       matched BY TreeId = id 11, FIA 400 hickory, DG Δ0.009"); tripling merely MULTIPLIES that
                #       one record ×3 (→3/243). Cycle 1 bit-exact ⇒ its inputs entering cycle 2 are bit-exact, so
                #       with bit-exact predictors the ONLY cycle-2-specific stochastic factor is the DGSCOR serial-
                #       correlation (dgscor.f, first applies at cycle 2 via the carried OLDRN). Its EXP is already
                #       fexp-routed. ★★ DGSCOR REFUTED by a live base/dgscor.f STAMP (2026-07-06): dumped FVS's
                #       per-tree pre-exp FRM for all 27 NOTRIPLE cycle-2 records and compared to jl's — 0/27
                #       mismatch, the DGSCOR frm is BIT-EXACT (incl. id 11: FVS 0.0673524439335 == jl). So the
                #       "irreducible DGSCOR" label was WRONG (the re-trace discipline paying off). The REAL seed is
                #       UPSTREAM: at cycle 1 (2000) id 11's DBH (9.71887) and Ht (76.0995) are BIT-EXACT but its
                #       CROWN differs — jl crown_pct 50 vs live PctCr 49 (crown CHANGE 1990→2000: jl +5 vs live +4).
                #       The 2000 crown feeds the 2000→2010 dgf! `crwn·cr` term ⇒ wk2 ⇒ dds5 ⇒ the 2010 DBH diverges
                #       (11.38251 vs 11.37337). VERDICT: NOT a primitive — a real, potentially FIXABLE CS crown-ratio
                #       model divergence (task: trace CS crown_ratio_update! vs FVS cratet/CS crown change). The
                #       libm routings (peff/htcalc/balmod/AVH/DGF-exp/log/bachlo) were all inert because the seed is
                #       the CROWN, not any transcendental. This is a lead to CLOSE the class, not accept it.
                @test ba == L[2]        # BA — BIT-EXACT every grown cycle
                push!(drift, (Float64(tpa), Float64(ba), Float64(sdi), Float64(ccf), Float64(topht), Float64(qmd),
                              Float64(L[1]), Float64(L[2]), Float64(L[3]), Float64(L[4]), Float64(L[5]), Float64(L[6])))
            end
            yr < 2090 && FVSjl.grow_cycle!(s; fint = 10f0)
        end
        # ★ ALL BIT-EXACT (2026-07-06): the crown-ratio raw-BA fix (crown_ratio.jl — the /gross_space band-aid was
        # the root, proven by the live dgscor.f + FVS_TreeList stamps above) closed every drift column. Was
        # @test_broken (mislabeled "height-transcendental/DGSCOR floor"); now green ==.
        @test all(d[1] == d[7]  for d in drift)                                # TPA — BIT-EXACT
        @test all(d[3] == d[9]  for d in drift)                                # SDI — BIT-EXACT
        @test all(d[4] == d[10] for d in drift)                                # CCF — BIT-EXACT
        @test all(d[5] == d[11] for d in drift)                                # TopHt — BIT-EXACT
        @test all(round(Int, d[6]*10) == round(Int, d[12]*10) for d in drift)  # QMD — BIT-EXACT (rendered tenth)
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
                # BIT-EXACT all cycles: the QMDGE5 cumulative-cap fix (this planted stand seeds species 3, a
                # CS cap-13 species, so the earlier local-per-tree cap gave the wrong QMDGE5 → biased the small-
                # tree RDBH/RDBHSQ growth terms). With the cumulative species-order cap, jl == live FVScs bit-
                # exact across every cycle (integer cols by `==`; QMD by rendered 1-decimal `==`). The old
                # "SIZE-only single-precision accumulation" comment was the QMDGE5 bug in disguise.
                @test r[1] == L[1]        # TPA — BIT-EXACT
                @test r[2] == L[2]        # BA
                @test r[3] == L[3]        # SDI
                @test r[4] == L[4]        # CCF
                @test r[5] == L[5]        # TopHt
                @test round(Float32(r[6]); digits=1) == Float32(L[6])   # QMD rendered to the .sum's 1 decimal
            end
        end
        rm(barekey; force = true); rm(joinpath(@__DIR__, "cst01.tre"); force = true)
    end
end

@testset "CS cst01 thinning prescriptions BIT-EXACT (vs live FVScs)" begin
    # Full multi-stand cst01.key: the THINDBH (stand 2, every 3rd cycle) and shelterwood THINPRSC/
    # SPECPREF/THINBTA (stands 3-4) prescriptions. The thinning SELECTION + the post-cut state are
    # bit-exact vs live (proving the cut logic is faithful for CS); the deep-thinned tails accumulate
    # the documented single-precision floor amplified at discrete thin/classification thresholds.
    key = "/workspace/ForestVegetationSimulator/tests/FVScs/cst01.key"
    if !isfile(key)
        @info "cst01.key not present; skipping CS thinning test"
    else
        cp(joinpath(dirname(key), "cst01.tre"), joinpath(@__DIR__, "cst01.tre"); force = true)
        cd(@__DIR__) do
            sumtxt = FVSjl.run_keyfile(key; variant = CentralStates(), output = :sum)
            stands = Vector{Vector{NTuple{7,Float64}}}()   # per stand: (yr,TPA,BA,SDI,CCF,TopHt,QMD)
            cur = nothing
            for ln in split(sumtxt, '\n')
                if startswith(ln, "-999"); cur = NTuple{7,Float64}[]; push!(stands, cur); continue; end
                t = split(strip(ln))
                (length(t) >= 8 && tryparse(Int, t[1]) !== nothing && cur !== nothing) || continue
                push!(cur, (parse(Float64,t[1]), parse(Float64,t[3]), parse(Float64,t[4]),
                            parse(Float64,t[5]), parse(Float64,t[6]), parse(Float64,t[7]), parse(Float64,t[8])))
            end
            @test length(stands) == 5
            # Stand 2 (THINDBH) and stand 3 (shelterwood THINPRSC): the cut + post-cut growth is BIT-EXACT
            # through cycle 2. Live values (per-gross-acre): (Year, TPA, BA, SDI, CCF, TopHt, QMD).
            s_thindbh = Dict(1990=>(536,77,160,169,63,5.1), 2000=>(518,99,196,202,68,5.9), 2010=>(476,122,231,234,70,6.9))
            s_shelter = Dict(1990=>(536,77,160,169,63,5.1), 2000=>(235,83,161,163,68,8.0), 2010=>(226,103,190,188,70,9.1))
            for (idx, expect) in ((2, s_thindbh), (3, s_shelter))
                rows = Dict(Int(r[1]) => r for r in stands[idx])
                for (yr, L) in expect
                    @test haskey(rows, yr)
                    r = rows[yr]
                    @test (Int(r[2]),Int(r[3]),Int(r[4]),Int(r[5]),Int(r[6])) == (L[1],L[2],L[3],L[4],L[5])
                    @test round(Float32(r[7]); digits=1) == Float32(L[6])
                end
            end
        end
        rm(joinpath(@__DIR__, "cst01.tre"); force = true)
    end
end

@testset "CS stump sprouting (SPROUT + clearcut, vs live FVScs)" begin
    # CS natural sprouting (SPROUT/ESUCKR): clearcut (THINBTA→20 BA) at 2000 with all-species SPROUT,
    # then ESUCKR creates stump sprouts the next cycle. The 2010 sprout-regeneration cycle is BIT-EXACT
    # vs live FVScs (cs/essprt.f CASE('CS') PREM/NSPREC/SPRTHT + the Wykoff sprout DBH). Validates the
    # last CS natural-process model. (cst01.tre = the CS inventory; site/forest as cst01.)
    tre = "/workspace/ForestVegetationSimulator/tests/FVScs/cst01.tre"
    if !isfile(tre)
        @info "cst01.tre not present; skipping CS sprout test"
    else
        key = joinpath(@__DIR__, "cs_sprout.key")
        cp(tre, joinpath(@__DIR__, "cs_sprout.tre"); force = true)   # TREEDATA reads <keystem>.tre
        open(key, "w") do io
            print(io, """
SCREEN
STATS
STDIDENT
CSSPROUT SPROUT TEST.
DESIGN                                        11.0       1.0
STDINFO          905                60.0     315.0      30.0      10.0       40.
SITECODE          19        60
INVYEAR       1990.0
NUMCYCLE        10.0
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,
T52,I2,T66,5I1,T54,7I1,T75,F3.0)
TREEDATA
ECHOSUM
ESTAB          2000.
SPROUT         2000.       0.0
END
THINBTA        2000.      20.0
PROCESS
STOP
""")
        end
        cd(@__DIR__) do
            sumtxt = FVSjl.run_keyfile(key; variant = CentralStates(), output = :sum)
            row2010 = nothing
            for ln in split(sumtxt, '\n')
                t = split(strip(ln))
                if length(t) >= 8 && t[1] == "2010"; row2010 = t; break; end
            end
            @test row2010 !== nothing
            if row2010 !== nothing
                # Live FVScs 2010 sprout-regen row: TPA 434, BA 23, SDI 44, CCF 40, TopHt 44, QMD 3.1 — bit-exact.
                @test parse(Int, row2010[3]) == 434       # TPA (the sprout cohort)
                @test parse(Int, row2010[4]) == 23        # BA
                @test parse(Int, row2010[5]) == 44        # SDI
                @test parse(Int, row2010[6]) == 40        # CCF
                @test parse(Int, row2010[7]) == 44        # TopHt
                @test round(parse(Float32, row2010[8]); digits = 1) == 3.1f0  # QMD
            end
        end
        rm(key; force = true); rm(joinpath(@__DIR__, "cs_sprout.tre"); force = true)
    end
end
