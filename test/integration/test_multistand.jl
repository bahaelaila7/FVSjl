# test_multistand.jl — multi-stand driver (each_stand). A keyword file holds several
# stands separated by PROCESS and ended by STOP. FVS re-runs INITRE per stand (ITRN
# resets) but the tree-record format (TREFMT) persists in COMMON across stands, and a
# stand with no TREEDATA still reads the shared tree file (initre.f:334 default INTREE)
# unless NOTREES. This guards three regressions that all surfaced together:
#   (1) stands 2+ fell back to the DEFAULT tree format → misparsed → garbage TPA,
#   (2) the FFE stand (no TREEDATA, only REWIND) loaded zero trees,
#   (3) the trailing STOP after the last PROCESS produced a phantom 6th stand.

using Test, FVSjl

const _MS_KEY = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"

@testset "multi-stand driver (each_stand) — snt01 5 stands" begin
    if !isfile(_MS_KEY)
        @test_skip "snt01.key not available"
    else
        stands = each_stand(_MS_KEY)
        @test length(stands) == 5                  # exactly 5 (no phantom terminator stand)

        # Stands 1-4 all read the SAME snt01.tre (27 trees) — via TREEDATA (1-3) or the
        # default INTREE after REWIND (4) — so every one starts at TPA 536 bit-exact,
        # which only holds if TREFMT persisted across stands.
        for i in 1:4
            s = stands[i]
            notre!(s)
            g = s.plot.gross_space
            @test s.trees.n == 27
            # BIT-EXACT vs live: the .sum prints trunc(perAcre+0.5); snt01 cyc0 = TPA 536 / BA 77 (live).
            # (Internal per-acre stand_tpa/g=536.048, stand_ba/g=77.39 → the .sum integers 536/77.) Was
            # atol=0.5/1.0 slack covering the pre-rounding value; now assert the live .sum integer exactly.
            @test trunc(Int, stand_tpa(s) / g + 0.5) == 536   # .sum TPA == live
            @test trunc(Int, stand_ba(s)  / g + 0.5) ==  77   # .sum BA == live
        end
        # Stand 4 (FFE) has its own inventory year 1993 and NO TREEDATA keyword.
        @test Int(stands[4].control.cycle_year[1]) == 1993

        # Stand 5 is the bare-ground PLANT stand: NOTREES ⇒ no tree-file read, ESTAB active.
        notre!(stands[5])
        @test stands[5].trees.n == 0
        @test stands[5].estab.active
    end
end

@testset "full multi-cycle driver (run_keyfile) — snt01" begin
    if !isfile(_MS_KEY)
        @test_skip "snt01.key not available"
    else
        txt = run_keyfile(_MS_KEY)
        rows = [strip(l) for l in split(txt, '\n') if !isempty(strip(l))]
        hdrs = filter(l -> startswith(l, "-999"), rows)
        @test length(hdrs) == 5                      # one -999 header per stand
        @test [split(h)[3] for h in hdrs] == ["S248112","S248112","S248112","FFE","BARE"]

        # Stand 1 (unthinned control) is bit-exact modulo Float32 ulp — assert the first
        # block's data rows match the Fortran baseline to ±2 on TPA/BA/cuft.
        savef = _MS_KEY[1:end-4] * ".sum.save"
        data = filter(l -> !startswith(l, "-999"), rows)[1:11]   # stand-1 block (11 cycles)
        base = isfile(savef) ?
            filter(l -> !startswith(l, "-999"),
                   [strip(l) for l in eachline(savef)])[1:11] : String[]
        isempty(base) && @test_skip "snt01.sum.save baseline not available"
        for (m, o) in zip(data, base)
            mf = split(m); of = split(o)
            @test mf[1] == of[1]                                       # year
            @test parse(Int, mf[3]) == parse(Int, of[3])               # TPA — BIT-EXACT
            @test parse(Int, mf[4]) == parse(Int, of[4])               # BA  — BIT-EXACT
        end
        # total cuft — BIT-EXACT bar a print-boundary ULP; the residual is the non-associative Float32
        # TREE-SUM accumulation order (doctrine #9: exposed, not a passing ±1).
        @test_broken all(parse(Int, split(m)[9]) == parse(Int, split(o)[9]) for (m, o) in zip(data, base))  # total cuft (col 9) — tree-sum order
        # NOTE: stands 2 (THINDBH in IF/THEN event monitor), 3 (THINPRSC), and 4 (FFE
        # fire) require management subsystems still being ported — the driver runs them
        # but their thinning/fire is not yet applied. Tracked, not asserted here.
    end
end
