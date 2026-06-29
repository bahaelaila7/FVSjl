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

# THINBBA (thin-from-below to a residual basal area) — a thinning path net01 itself does NOT
# exercise (net01 uses THINDBH/THINBTA). Built by injecting a THINBBA into net01 stand-1 and
# validated vs the live FVSne oracle: at 2000 the thin takes the UNTHINNED-CONTROL stand
# (TPA 524 / BA 107 / QMD 6.1) down to residual BA ~55 (target 60), removing the small trees
# so TPA→83 and QMD→10.8. Also pins the .sum -999 header variant code = "NE" (was defaulting
# to "SN" — run_keyfile now threads variant_code(s.variant)). See docs/NE_PORT_STATUS.md.
@testset "net01 (NE) THINBBA thin-from-below — vs live FVSne + header variant" begin
    if !isfile(_NET01_KEY)
        @test_skip "net01.key not available"
    else
        # Inject a THINBBA (residual BA 60) before stand-1's PROCESS; run through run_keyfile.
        raw = read(_NET01_KEY)                       # net01.key is CR-delimited
        recs = split(String(raw), '\r')
        thin = rpad("THINBBA", 10) * lpad("2000.0", 10) * lpad("60.0", 10)   # 10-col FVS fields
        pidx = findfirst(r -> strip(r) == "PROCESS", recs)
        insert!(recs, pidx, thin)
        dir = mktempdir()
        kp = joinpath(dir, "thinbba.key")
        write(kp, join(recs, '\r'))
        cp(joinpath(dirname(_NET01_KEY), "net01.tre"), joinpath(dir, "thinbba.tre"))
        out = FVSjl.run_keyfile(kp; variant = Northeast())
        # Header variant must read NE (the fix), not the SN default.
        @test occursin(r"-999.*\bNE\b", first(split(out, '\n')))
        # Parse stand-1's 2000 row (2nd data row of the 1st -999 block).
        lines = split(out, '\n')
        b1 = findfirst(l -> startswith(l, "-999"), lines)
        row2000 = split(lines[b1 + 2])                # 1990 = b1+1, 2000 = b1+2 (thin year)
        row2010 = split(lines[b1 + 3])                # post-thin stand shows at the next period
        @test parse(Int, row2000[1]) == 2000
        # 2000 col3/col4 = start-of-period (pre-thin) TPA/BA. Live: 524 / 107.
        @test parse(Int, row2000[3]) == 524           # pre-thin TPA (live 524)
        @test parse(Int, row2000[4]) == 107           # pre-thin BA  (live 107)
        # 2010 col3/col4 = the THINNED residual grown one cycle. Live: TPA 83 / BA 67.
        # Confirms the thin-from-below removed the small trees (524→83) to residual BA ~55.
        @test parse(Int, row2010[3]) ≈ 83 atol = 2    # post-thin residual TPA (live 83)
        @test parse(Int, row2010[4]) ≈ 67 atol = 2    # post-thin residual BA  (live 67)
    end
end
