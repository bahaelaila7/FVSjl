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
