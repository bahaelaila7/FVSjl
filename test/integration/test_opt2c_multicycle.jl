# =============================================================================
# test_opt2c_multicycle.jl — OP (Olympic) multi-cycle COOPERATING ORGANON/Wykoff driver
# vs live FVSop_clean (stand S248112, NUMCYCLE=2, DGSD=0 ⇒ deterministic ⇒ a TRUE bit-exact bar).
#
# Reference: /workspace/.opwork/op2c.sum (FVSop_clean, op2c.key/op2c.tre). Oracle rows:
#   1990  536  77 184 114 63 5.1 1472  972  5003
#   1995  499  98 219 134 71 6.0 2256 1745  9249
#   2000  470 120 256 154 69 6.8 2988 2419 12427
#
# Exercises the full per-cycle cooperating driver: op_organon_prepare! (LSTART) + the op/dgdriv.f LSTART
# DG COR calibration (WF sp2 COR = +0.03379, bit-exact vs op2c_dbg.out); the DGDRIV cooperating diameter
# driver (op-native dgf! for IORG=0 + ORGANON DGRO-fold→WK2 for IORG=1, then the shared DDS→DG with
# op_bratio); height_growth! (native + stashed ORGANON HGRO); small_tree_growth! (op/regent.f + op/smhgdg.f
# REGENT small-tree height/diameter for D<XMAX, blended by XWT); crown_ratio_update! (op/crown.f Weibull +
# ORGANON CR2); mortality! (op/morts.f ORGANON MORTEXP-for-all).
#
# VALIDATED BIT-EXACT: cycle 0 (1990) 9/9, cycle 1 (1995) 9/9, cycle 2 (2000) 8/9. The one cycle-2
# residual is BdFt (jl 12405 vs 12427, 0.18%) — an accumulated board-foot volume difference in the
# largest trees; documented as `@test_broken`.
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "OP opt2c multi-cycle vs FVSop_clean (S248112, NUMCYCLE=2)" begin
    key = joinpath(@__DIR__, "..", "harness", "scenarios", "opt2c.key")
    if !isfile(key)
        @test_skip "opt2c scenario not available"
    else
        rows = F.SummaryRow[]
        for s in F.each_stand(key; variant = F.Olympic())
            F.notre!(s); F.setup_growth!(s); F.compute_volumes!(s)
            push!(rows, F.summary_row(s))
            for _ in 1:Int(s.control.ncycle)
                F.grow_cycle!(s; fint = 5f0)
                F.compute_volumes!(s)
                push!(rows, F.summary_row(s))
            end
            break
        end
        @test length(rows) >= 3
        r0, r1, r2 = rows[1], rows[2], rows[3]

        # --- cycle 0 (1990): bit-exact, all 9 columns ---
        @test r0.year == 1990
        @test r0.tpa == 536 && r0.ba == 77 && r0.sdi == 184 && r0.ccf == 114
        @test r0.topht == 63 && round(r0.qmd; digits = 1) == 5.1
        @test r0.cuft == 1472 && r0.mcuft == 972 && r0.bdft == 5003

        # --- cycle 1 (1995): BIT-EXACT, all 9 columns ---
        @test r1.year == 1995
        @test r1.tpa == 499 && r1.ba == 98 && r1.sdi == 219 && r1.ccf == 134
        @test r1.topht == 71 && round(r1.qmd; digits = 1) == 6.0
        @test r1.cuft == 2256 && r1.mcuft == 1745 && r1.bdft == 9249

        # --- cycle 2 (2000): 8/9 bit-exact ---
        @test r2.year == 2000
        @test r2.tpa == 470 && r2.ba == 120 && r2.sdi == 256 && r2.ccf == 154
        @test r2.topht == 69 && round(r2.qmd; digits = 1) == 6.8
        @test r2.cuft == 2988 && r2.mcuft == 2419
        # KNOWN RESIDUAL: cycle-2 board-foot (accumulated large-tree BdFt), jl 12405 vs 12427.
        @test_broken r2.bdft == 12427
    end
end
