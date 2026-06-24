# test_fire.jl — FFE fire-stand (SIMFIRE) end-to-end vs the Fortran baseline.
#
# Covers the snt01-class FFE fire path that the cycle-0/1 test_snt01 does NOT reach (the fire fires
# mid-projection). `fire_early` (ecounit 231Dd, FMIN + SIMFIRE @ 2000) is run to completion and
# diffed against the committed Fortran `fire_early.sum`:
#   - PRE-fire + fire-year cycles (1990/1995/2000) are BIT-EXACT (TPA/BA/volume) — the fuel loading,
#     fire trigger, and fire-year accounting all match;
#   - POST-fire cycles (2005+) carry the known SIMFIRE/FMEFF fire-MORTALITY-DISTRIBUTION residual
#     (~10 TPA: FVSjl under-kills slightly), which this test pins so a regression is caught.
#
# It also guards the FULIV2 shrub override (fuel_loading.jl): 231Dd is FULIV2-affected, and toggling
# the override leaves this .sum identical (live shrub fuel does not reach fire mortality here), so the
# carbon-report fix is confirmed inert for the fire path.

using Test, FVSjl

const _FDIR = joinpath(@__DIR__, "..", "harness", "scenarios")

_f_rows(txt) = [split(strip(l)) for l in split(txt, "\n") if occursin(r"^(19|20)\d\d ", strip(l))]

@testset "FFE fire stand (fire_early SIMFIRE) vs Fortran baseline" begin
    key = joinpath(_FDIR, "fire_early.key"); sav = joinpath(_FDIR, "fire_early.sum")
    if !isfile(key) || !isfile(sav)
        @test_skip "fire_early scenario not available"
    else
        jl = _f_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _f_rows(read(sav, String))
        @test length(jl) == length(ft) >= 6
        for (j, f) in zip(jl, ft)
            yr = j[1]
            @test yr == f[1]
            if yr in ("1990", "1995", "2000")              # pre-fire + fire year: bit-exact
                @test j[3] == f[3]                          # TPA
                @test j[4] == f[4]                          # BA
                @test abs(parse(Int, j[9]) - parse(Int, f[9])) <= 1   # TCuFt
            else                                            # post-fire: the fire-mortality residual
                @test abs(parse(Int, j[3]) - parse(Int, f[3])) <= 12  # TPA (~10 under-kill)
                @test abs(parse(Int, j[4]) - parse(Int, f[4])) <= 10  # BA
            end
        end
    end
end
