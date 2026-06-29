# test_net01.jl — NE variant net01 validation vs the live FVSne oracle (tests/FVSne/net01.sum.save).
# FVSjulia/Oracle-A has NO NE, so the committed .sum.save is the sole ground truth. Cycle-0 stand state
# (TPA/BA/QMD/TopHt) needs only tree-parse + the shared density — the first bit-exact NE gate.
using Test
using FVSjl

const _NET01_KEY = "/workspace/ForestVegetationSimulator/tests/FVSne/net01.key"

@testset "net01 (NE) cycle-0 stand state — bit-exact vs FVSne oracle" begin
    if !isfile(_NET01_KEY)
        @test_skip "net01.key not available"
    else
        n = first(FVSjl.each_stand(_NET01_KEY; variant = Northeast()))
        FVSjl.notre!(n)
        g = n.plot.gross_space
        di(x) = trunc(Int, x + 0.5)
        # net01.sum.save stand-1 (UNTHINNED) 1990 row: 536 77 160 176 63 5.1 ...
        @test di(stand_tpa(n) / g) == 536               # TPA   — BIT-EXACT
        @test di(stand_ba(n) / g) == 77                 # BA    — BIT-EXACT
        @test di(stand_sdi(n) / g) == 160               # SDI   — BIT-EXACT (Reineke, no per-species data)
        @test di(stand_ccf(n) / g) == 176               # CCF   — BIT-EXACT (NE cwcalc map + FORKOD lat/long → Hopkins)
        @test round(stand_qmd(n); digits = 1) == 5.1f0  # QMD   — BIT-EXACT
        @test di(stand_top_height(n)) == 63             # TopHt — BIT-EXACT
        # NE-specific parse correctness: first tree is Jack Pine (JP=19), dbh 11.5
        @test n.trees.species[1] == 19
        @test n.trees.dbh[1] ≈ 11.5f0
        # Volume columns (R9 Clark cubic + International-¼" board feet). Live 1990 row:
        # TCuFt 1558 MCuFt 1347 SCuFt 292 BdFt 1633.
        FVSjl.setup_growth!(n); FVSjl.compute_volumes!(n)
        bdft = sum(n.trees.bdft_vol[i] * n.trees.tpa[i] for i in 1:n.trees.n) / g
        @test di(bdft) ≈ 1633 atol = 8       # International ¼" board feet (R9LOGS + r9bdft)
    end
end

@testset "net01 (NE) cycle-1 growth — vs live FVSne (stand 1, unthinned)" begin
    if !isfile(_NET01_KEY)
        @test_skip "net01.key not available"
    else
        # Exercises the full NE growth spine end-to-end: diameter growth (BAL + the YR=10
        # gradd FINT/YR scale), height growth (NC-128 curve + BAL modifier), mortality
        # (background ·0.5 + SDI density), the TWIGS crown model, and REGENT small-tree
        # growth (dbh<5). Live FVSne stand-1 2000 row: TPA 524 BA 107 SDI 213 CCF 229 TopHt 72 QMD 6.1.
        n = first(FVSjl.each_stand(_NET01_KEY; variant = Northeast()))
        FVSjl.notre!(n); FVSjl.setup_growth!(n); FVSjl.compute_volumes!(n)
        g = n.plot.gross_space; di(x) = trunc(Int, x + 0.5)
        FVSjl.grow_cycle!(n; fint = 10f0)
        @test di(stand_tpa(n) / g) == 524                       # TPA — EXACT vs live (background-mortality fix)
        @test di(stand_ba(n) / g) ≈ 107 atol = 2                # BA  — live 107 (jl 106, post-REGENT)
        @test di(stand_sdi(n) / g) ≈ 213 atol = 4               # SDI — live 213 (jl 211)
        @test round(stand_qmd(n); digits = 1) ≈ 6.1 atol = 0.1  # live 6.1 (jl 6.1)
        @test di(stand_top_height(n)) ≈ 72 atol = 2             # TopHt — live 72 (jl 71)
    end
end

# NE thinning-keyword breadth — three cut-selection paths net01 itself does NOT exercise (it
# uses THINDBH/THINBTA). Each injects one THIN* before stand-1's PROCESS and is validated vs the
# live FVSne oracle. The cut-selection (cuts.f) is variant-agnostic shared code; these confirm it
# drives correctly off the NE growth/volume that feeds it. The UNTHINNED-CONTROL stand at 2000 is
# TPA 524 / BA 107 / QMD 6.1; each row is (target, live-2010-resid-TPA, live-2010-resid-BA):
#   THINBBA 60  (resid BA, thin from BELOW — removes small trees, QMD↑) → 2010  83 / 67
#   THINABA 60  (resid BA, thin from ABOVE — removes large trees, QMD↓) → 2010 434 / 81
#   THINSDI 120 (resid SDI target)                                      → 2010 263 / 75
# jl tracks live within the documented cyc-1 DG/volume ULP drift (#50). Also pins the .sum -999
# header variant code = "NE" (was defaulting to "SN"; run_keyfile now threads variant_code).
@testset "net01 (NE) thinning-keyword breadth — vs live FVSne + header variant" begin
    if !isfile(_NET01_KEY)
        @test_skip "net01.key not available"
    else
        function run_thin(kw, val)
            recs = split(String(read(_NET01_KEY)), '\r')   # net01.key is CR-delimited
            thin = rpad(kw, 10) * lpad("2000.0", 10) * lpad(val, 10)   # 10-col FVS fields
            insert!(recs, findfirst(r -> strip(r) == "PROCESS", recs), thin)
            dir = mktempdir()
            kp = joinpath(dir, "$(lowercase(kw)).key")
            write(kp, join(recs, '\r'))
            cp(joinpath(dirname(_NET01_KEY), "net01.tre"), joinpath(dir, "$(lowercase(kw)).tre"))
            FVSjl.run_keyfile(kp; variant = Northeast())
        end
        for (kw, val, resid_tpa, resid_ba) in
            (("THINBBA", "60.0", 83, 67), ("THINABA", "60.0", 434, 81), ("THINSDI", "120.0", 263, 75))
            out = run_thin(kw, val)
            lines = split(out, '\n')
            @test occursin(r"-999.*\bNE\b", lines[1])     # header variant = NE (the fix)
            b1 = findfirst(l -> startswith(l, "-999"), lines)
            row2000 = split(lines[b1 + 2])                # 1990=b1+1, 2000=b1+2 (thin year)
            row2010 = split(lines[b1 + 3])                # post-thin residual shows next period
            @test parse(Int, row2000[1]) == 2000
            @test parse(Int, row2000[3]) == 524           # pre-thin TPA (live 524)
            @test parse(Int, row2000[4]) == 107           # pre-thin BA  (live 107)
            @test parse(Int, row2010[3]) ≈ resid_tpa atol = 3   # post-thin residual TPA vs live
            @test parse(Int, row2010[4]) ≈ resid_ba  atol = 2   # post-thin residual BA  vs live
        end
    end
end

# A1 (NE audit) — net01 stand 2 (the THINDBH "TEST EXPANDED THINDBH OPTION", 16 cycles, repeated
# every-3-cycle multi-class thin). This is the scenario whose BA diverged ~−20 (the "−40 vs −22 BA"
# flag) until the ne_badist! /gross_space bug was fixed (it scaled the BAL competition array 10/11,
# over-growing low-competition trees as the thinned stand opened up). BA now tracks live FVSne within
# ±1 across all 16 cycles. Live FVSne s2.sum BA by year: 77,65,80,96,114,133,149,165,183,195,196,198…
@testset "net01 (NE) stand-2 THINDBH multi-cycle BA — vs live FVSne (audit A1, badist fix)" begin
    if !isfile(_NET01_KEY)
        @test_skip "net01.key not available"
    else
        thindbh = [
            "IF(FRAC(CYCLE/3.0) EQ 0.0)THEN",
            "THINDBH                              4.0      1.00       5.0",
            "THINDBH                              2.0      0.01               300.0",
            "THINDBH                    2.0       4.0      0.01               200.0",
            "THINDBH                    4.0       8.0      0.01               125.0",
            "THINDBH                    8.0      12.0      0.01                60.0",
            "THINDBH                   12.0      16.0      0.01                35.0",
            "THINDBH                   16.0      20.0      0.01                15.0",
            "THINDBH                   20.0                1.00", "ENDIF"]
        recs = ["SCREEN", "NOAUTOES", "STATS", "STDIDENT", "S248112  THINDBH ISOLATED",
            "DESIGN                                        11.0       1.0",
            "STDINFO        922.0                60.0     315.0      30.0      20.0",
            "SITECODE          13        26", "INVYEAR       1990.0", "NUMCYCLE        16.0", "TREEFMT",
            "(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,",
            "T52,I2,T66,5I1,T54,7I1,T75,F3.0)", thindbh...,
            "TREEDATA", "ECHOSUM", "PROCESS", "STOP"]
        dir = mktempdir()
        kp = joinpath(dir, "s2.key")
        write(kp, join(recs, '\r') * '\r')
        cp(joinpath(dirname(_NET01_KEY), "net01.tre"), joinpath(dir, "s2.tre"))
        out = FVSjl.run_keyfile(kp; variant = Northeast())
        lines = split(out, '\n')
        b1 = findfirst(l -> startswith(l, "-999"), lines)
        # live FVSne BA at each cycle row (1990=b1+1 … 2150=b1+17)
        live_ba = (77, 65, 80, 96, 114, 133, 149, 165, 183, 195, 196, 198, 199, 200, 200, 201, 201)
        for (k, lba) in enumerate(live_ba)
            row = split(lines[b1 + k])
            @test parse(Int, row[4]) ≈ lba atol = 2     # BA tracks live within ±1 (was ~−20 pre-fix)
        end
        # cycle-0 is the bit-exact anchor; the deep cycles are where the badist bug used to compound
        @test parse(Int, split(lines[b1 + 1])[4]) == 77   # 1990 BA — bit-exact
    end
end

# A2 (NE audit) — full-species-set cycle-0 volume/crown/density. net01 exercises ~6 of 108 NE
# species; this rewrites net01.tre's 30 records to a diverse 30-species sample (conifers BF/WS/RS/
# NS/RN/WP + hardwoods RM/SM/BM/YB/HI/AB/oaks…) and checks jl's cycle-0 stand against the live FVSne
# oracle: TPA/BA exact, stand volume within ULP. Guards the per-species volume (R9 Clark + R9LOGS),
# crown (CWCALC), and density coefficient loading for the broad species set. See docs/NE_AUDIT.md A2.
@testset "net01 (NE) multi-species cycle-0 volume — vs live FVSne (audit A2)" begin
    if !isfile(_NET01_KEY)
        @test_skip "net01.key not available"
    else
        sp30 = ["BF","WS","RS","NS","RN","WP","AW","EH","HM","JP","SP","RM","SM","BM","YB",
                "PB","HI","SH","AB","PA","BP","CK","SW","CB","HK","OO","BK","RL","ST","PR"]
        tre = readlines(joinpath(dirname(_NET01_KEY), "net01.tre"))
        out_tre = String[]
        k = 0
        for ln in tre
            if length(ln) < 36
                push!(out_tre, ln)
            else
                sp = sp30[mod1(k + 1, length(sp30))]; k += 1
                push!(out_tre, ln[1:33] * rpad(sp, 2) * ln[36:end])
            end
        end
        recs = split(String(read(_NET01_KEY)), '\r')
        dir = mktempdir()
        write(joinpath(dir, "spv.key"), join(vcat(recs[1:16], ["STOP"]), '\r') * "\r")
        write(joinpath(dir, "spv.tre"), join(out_tre, '\n') * "\n")
        out = FVSjl.run_keyfile(joinpath(dir, "spv.key"); variant = Northeast())
        row = split(split(out, '\n')[findfirst(l -> startswith(l, "1990"), split(out, '\n'))])
        # Live FVSne cycle-0: TPA 536 BA 77 CCF 146 | TCuFt 1551 MCuFt 1286 SCuFt 186 BdFt 1023.
        @test parse(Int, row[3]) == 536                       # TPA — exact
        @test parse(Int, row[4]) == 77                        # BA  — exact
        @test parse(Int, row[6]) == 146                       # CCF — exact (CWCALC across species)
        @test parse(Int, row[9])  ≈ 1551 atol = 4             # TCuFt (R9 Clark cubic)
        @test parse(Int, row[12]) ≈ 1023 atol = 8             # BdFt  (R9LOGS Scribner)
    end
end

# A1 (NE audit) — the ESSUBH establishment base-height formula (essubh.f:73-82): NE plants seedlings at
# HHT = (NC-128 site-curve height at the per-species reference age CARAGE / CARAGE) · min(5, period−delay),
# NOT the site-curve height at the tree's age. Verified BIT-EXACT vs live FVSne ESSUBH (BF 5.159, WS 4.591).
@testset "net01 (NE) ESSUBH establishment base height — bit-exact vs live FVSne" begin
    # BF (sp1, SI 52, refage 20) and WS (sp3, SI 50, refage 15); period 10, delay 0 ⇒ min(5, 10)=5.
    for (sp, si, want) in ((1, 52f0, 5.159f0), (3, 50f0, 4.591f0))
        ca = Float32(FVSjl._NE_ESSUBH_REFAGE[sp])
        hht = (FVSjl.ne_htcalc_height(sp, si, ca) / ca) * min(5f0, 10f0 - 0f0)
        @test hht ≈ want atol = 0.002        # = live FVSne ESSUBH HHT (essubh.f), bit-exact
    end
    @test FVSjl._NE_ESSUBH_REFAGE[1] == 20 && FVSjl._NE_ESSUBH_REFAGE[3] == 15   # essubh.f DATA MAPNE
    @test length(FVSjl._NE_ESSUBH_REFAGE) == 108
end

# A1 breadth (NE audit) — net01 stands 3 (shelterwood) + 5 (BARE establishment) vs live FVSne. Confirms
# the post-badist growth spine + the silvicultural treatments track live across the no-fire stands (1,2,3,5).
# Stand 3 = THINPRSC(0.999) shelterwood + SPECPREF + THINBTA(157→35): the prescription-thin TPA is bit-exact
# every cycle; BA within ±2 (the sp9 WP large-tree distributional tail). Stand 5 = NOTREES + ESTAB/PLANT
# (jack+white pine 400 each): regen TPA tracks (±6), BA converges to bit-exact late (early-cohort REGENT runs
# ~20% low on the tiny establishment BA, 8 vs 10 at 2002, → 265/265 by 2092).
@testset "net01 (NE) stand-3 shelterwood + stand-5 BARE — vs live FVSne (audit A1 breadth)" begin
    if !isfile(_NET01_KEY)
        @test_skip "net01.key not available"
    else
        function run_recs(recs, tre)
            dir = mktempdir()
            write(joinpath(dir, "s.key"), join(recs, '\r') * '\r')
            tre === nothing || cp(joinpath(dirname(_NET01_KEY), "net01.tre"), joinpath(dir, "s.tre"))
            split(FVSjl.run_keyfile(joinpath(dir, "s.key"); variant = Northeast()), '\n')
        end
        hdr = ["SCREEN", "NOAUTOES", "STATS", "STDIDENT", "S248112  SHELTERWOOD",
            "DESIGN                                        11.0       1.0",
            "STDINFO        922.0                60.0     315.0      30.0      20.0",
            "SITECODE          13        56", "INVYEAR       1990.0", "NUMCYCLE        10.0",
            "THINPRSC      1990.0     0.999",
            "SPECPREF      2020.0      27.0     999.0", "SPECPREF      2020.0      19.0    9999.0",
            "THINBTA       2020.0     157.0",
            "SPECPREF      2050.0      49.0    -999.0", "SPECPREF      2050.0       9.0     -99.0",
            "THINBTA       2050.0      35.0", "TREEFMT",
            "(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,",
            "T52,I2,T66,5I1,T54,7I1,T75,F3.0)", "TREEDATA", "ECHOSUM", "PROCESS", "STOP"]
        lines = run_recs(hdr, :tre)
        b = findfirst(l -> startswith(l, "-999"), lines)
        live_tpa3 = (536, 235, 230, 218, 140, 137, 135, 31, 30, 30, 29)   # live FVSne — prescription-thin
        for (k, ltpa) in enumerate(live_tpa3)
            row = split(lines[b + k])
            @test parse(Int, row[3]) == ltpa                  # TPA — BIT-EXACT (THINPRSC/THINBTA thinning)
        end
        @test parse(Int, split(lines[b + 1])[4]) == 77        # 1990 BA exact
        @test parse(Int, split(lines[b + 4])[4]) ≈ 134 atol = 2   # 2020 BA within ±2 (WP tail)

        bare = ["SCREEN", "NOAUTOES", "STDIDENT", "BARE GROUND PLANT", "ECHOSUM", "SCREEN", "NOTREES",
            "NOTRIPLE", "STDINFO        922.0                 0.0     315.0      30.0      20.0",
            "INVYEAR         1992", "NOAUTOES", "ESTAB           1992",
            "PLANT           1992         1       400", "PLANT           1992         3       400", "END",
            "NUMCYCLE          10", "PROCESS", "STOP"]
        l5 = run_recs(bare, nothing)
        b5 = findfirst(l -> startswith(l, "-999"), l5)
        @test parse(Int, split(l5[b5 + 1])[3]) == 0           # 1992 BARE — no trees
        @test parse(Int, split(l5[b5 + 2])[3]) == 800         # 2002 regen TPA exact (PLANT 400+400)
        @test parse(Int, split(l5[b5 + 11])[3]) ≈ 499 atol = 6   # 2092 TPA tracks live (mortality)
        @test parse(Int, split(l5[b5 + 11])[4]) ≈ 265 atol = 2   # 2092 BA — converges to live (bit-exact)
    end
end

# A1 end-to-end (NE audit) — the FULL net01.key (all 5 stands incl. stand 4 FFE/SIMFIRE) runs end-to-end
# in jl and tracks live FVSne. This RETIRES the stale "NE FFE unported / :v2t" claim: the shared FFE model
# (src/engine/fire/*) + data/northeast/fire_species_props.csv handle NE fire faithfully. Per-stand max diff
# vs live over 56 rows: stand 4 (FFE) BIT-EXACT (TPA 0/BA 0/TopHt 1); stand 3 TPA 0/BA 2; stand 1 BA 1
# (TopHt 3 = sp9 WP tail); stand 5 TPA 7/BA 4 (early regen); stand 2 BA 6 (WP tail × repeated thinning).
@testset "net01 (NE) FULL keyfile end-to-end (5 stands incl. FFE) — vs live FVSne (audit A1)" begin
    if !isfile(_NET01_KEY)
        @test_skip "net01.key not available"
    else
        out = FVSjl.run_keyfile(_NET01_KEY; variant = Northeast())
        lines = split(out, '\n')
        nstands = count(l -> startswith(l, "-999"), lines)
        @test nstands == 5                                    # all 5 stands run end-to-end (FFE no longer crashes)
        # stand 4 (FFE) is the 4th -999 block; its post-fire 2013 row must match live (TPA 168, BA 81).
        hdrs = findall(l -> startswith(l, "-999"), lines)
        s4 = hdrs[4]
        # rows under stand 4: 1993(+1) 2003(+2) 2013(+3) — the SIMFIRE 2003 fires, 2013 is post-fire
        r2013 = split(lines[s4 + 3])
        @test parse(Int, r2013[1]) == 2013
        @test parse(Int, r2013[3]) == 168                     # post-fire TPA — BIT-EXACT vs live
        @test parse(Int, r2013[4]) == 81                      # post-fire BA  — BIT-EXACT vs live
    end
end
