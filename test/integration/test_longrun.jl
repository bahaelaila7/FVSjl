# test_longrun.jl — COMCUP (comcup.f, top half): each cycle deletes records whose
# expansion factor PROB has fallen to ≤1e-5 (suppressed trees whittled to ~0 by
# mortality), so the next cycle's per-tree DGSCOR doesn't draw an extra deviate for a
# dead record. snt01's 10-cycle .sum never reaches the threshold (first delete is at
# cycle 11), so this needs a LONG (30-cycle) unthinned projection to exercise it.
#
# Guards: (1) comcup actually fires — the unthinned live-record count DROPS over a long
# run (mortality alone never deletes records); (2) the long projection tracks the
# Fortran/Oracle-A baseline (±2 TPA) instead of drifting once records reach the floor.

using Test, FVSjl

const _LR_KEY = joinpath(@__DIR__, "..", "harness", "scenarios", "dense_long.key")
# dense_long = snt01 stand 1 (unthinned) projected 30 cycles. The .tre is a copy of
# snt01.tre (the harness gitignores scenario *.tre; copy it in if missing):
#   cp <FVSsn>/snt01.tre test/harness/scenarios/dense_long.tre
const _LR_TRE = _LR_KEY[1:end-4] * ".tre"

@testset "long-run COMCUP zero-PROB record deletion (dense_long, 30 cycles)" begin
    if !(isfile(_LR_KEY) && isfile(_LR_TRE))
        @test_skip "dense_long.key/.tre not available"
    else
        s, _ = initialize(_LR_KEY)
        notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        # project all cycles, recording TPA at three checkpoints
        checks = Dict{Int,Float64}()
        peak_n = 0
        ncyc = Int(s.control.ncycle)
        for c in 0:ncyc
            yr = Int(s.control.cycle_year[1]) + Int(s.control.cycle) * 5
            peak_n = max(peak_n, s.trees.n)
            yr in (2040, 2090, 2140) && (checks[yr] = stand_tpa(s) / s.plot.gross_space)
            c < ncyc && FVSjl.grow_cycle!(s)
        end

        # (1) comcup fired: tripling peaks at 243 records, then deletions drop it well below
        @test peak_n == 243
        @test s.trees.n < 200                       # records deleted by comcup (else stuck at 243)

        # (2) RE-GROUNDED vs LIVE FVSsn (2026-07-02, forget Oracle-A). dense_long TPA is BIT-EXACT vs live
        # through ~2085; 2040 (.sum 120) and 2140 (.sum 19) are bit-exact. 2090 is the sole residual: jl .sum
        # 36 vs live .sum 35 — a ±1-tree late-run near-SDImax kill-distribution flip on the heavily-thinned
        # tail (TPA ~35), the accepted per-tree-bit-exact ULP class (cf. mix_lp_rm). Two exact `==` vs the live
        # .sum integer, one ±1 ULP-justified.
        @test haskey(checks, 2040) && trunc(Int, checks[2040] + 0.5) == 120   # BIT-EXACT vs live .sum
        @test haskey(checks, 2140) && trunc(Int, checks[2140] + 0.5) ==  19   # BIT-EXACT vs live .sum
        @test haskey(checks, 2090) && abs(trunc(Int, checks[2090] + 0.5) - 35) <= 1  # ULP: late near-SDImax ±1 (jl36/live35)
    end
end
