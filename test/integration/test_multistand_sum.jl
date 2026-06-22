# test_multistand_sum.jl — multi-stand .sum parity vs the live-Fortran baseline.
#
# The earlier multi-stand test only checked structure + cycle-0 TPA, and the sweep
# harness emitted only stand 1 — so stands 2-5 of a multi-stand key were never diffed
# against Fortran. This test closes that gap: run_keyfile must emit ALL 5 stand blocks
# (55 rows) and the two "clean" stands must track Fortran every cycle. Stand 3 in
# particular inherits cross-stand state (TREFMT persists, RNG carries) from stands 1-2,
# so its agreement is the real guard against a multi-stand state-carry bug.
#
# snt01's five stands and their expected behaviour vs Fortran:
#   1 unthinned          → ulp tail (validated here, ±1 TPA/BA, ±8 cuft)
#   2 IF/THEN THINDBH     → known 3rd-thin (2030) class-boundary residual (~187 vs 158)
#   3 THINPRSC shelterwd  → ulp tail (validated here) — the state-carry guard
#   4 FFE fire            → diverges (C7 fire not built) — only structure checked
#   5 BARE PLANT regen    → regen tail — only structure checked

using Test, FVSjl

const _MSS_KEY  = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"
const _MSS_BASE = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.sum.save"

_mss_rows(txt) = [split(l) for l in split(txt, "\n")
                  if !occursin("-999", l) && length(split(l)) >= 11]

@testset "multi-stand .sum parity vs Fortran (run_keyfile, all 5 stands)" begin
    if !isfile(_MSS_KEY) || !isfile(_MSS_BASE)
        @test_skip "snt01 key/baseline not available"
    else
        jl = _mss_rows(FVSjl.run_keyfile(_MSS_KEY; faithful = true))
        ft = _mss_rows(read(_MSS_BASE, String))

        # every stand emitted (regression: harness once produced only stand 1)
        @test length(jl) == 55
        @test length(ft) == 55

        if length(jl) == 55 && length(ft) == 55
            tpa(r) = parse(Float64, r[3]); ba(r) = parse(Float64, r[4]); cuft(r) = parse(Float64, r[9])
            # stands 1 (unthinned) and 3 (THINPRSC, inherits cross-stand state) must
            # track Fortran every cycle to within the single-precision tail.
            for st in (1, 3), c in 1:11
                i = (st - 1) * 11 + c
                @test abs(tpa(jl[i])  - tpa(ft[i]))  <= 1
                @test abs(ba(jl[i])   - ba(ft[i]))   <= 1
                @test abs(cuft(jl[i]) - cuft(ft[i])) <= 8
            end
            # stand 2 is thinned (IF/THEN THINDBH at 2000/2015/2030) — it must show a
            # sharp thinning drop (>20% in one cycle, far above ~5-15% self-mortality),
            # and the managed regime ends with MORE trees than stand 1's unthinned stand,
            # which self-thins harder via SDI mortality (182 vs 120).
            s2 = [parse(Float64, jl[11 + c][3]) for c in 1:11]
            @test any(c -> s2[c] < 0.80 * s2[c-1], 2:11)
            @test tpa(jl[2 * 11]) > tpa(jl[1 * 11])
            # stand 5 (BARE) starts from zero trees and regenerates
            @test tpa(jl[4 * 11 + 1]) == 0.0
            @test tpa(jl[4 * 11 + 5]) > 0.0
        end
    end
end
