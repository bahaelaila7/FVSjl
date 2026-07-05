# =============================================================================
# test_simfire_schedule.jl — SIMFIRE date/cycle default + multiple-fire scheduling (D9).
#
# Two faithfulness bugs the divergence sweep surfaced (s10_fire 789% / fire_repeat 288% TPA):
#
#  (1) A no-parameter SIMFIRE must fire in CYCLE 1, not "never". FVS fmin.f:309 defaults the date
#      field IDT to 1; opexpn.f:40-44 treats a value ≤ MAXCYC as a 1-based CYCLE number and converts
#      it to that cycle's start year. jl had left fire_year=0 (no fire) when the field was blank.
#  (2) Each SIMFIRE keyword is its own OPNEW activity, so >1 keyword schedules >1 fire. jl stored a
#      single scalar fire_year that the second SIMFIRE overwrote — only the last fire fired. Now a
#      `fire_schedule` list holds every event; each fires in the cycle whose range contains its year.
#
# Plus the cycle-1 fuel-init fix: a fire in the FIRST FFE cycle burns before any ffe_fuel_update!
# loaded the dead-fuel pools, so summary.jl now runs FMCBA before stashing the fire's fuel basis
# (else zero cwd ⇒ low-fuel model ⇒ under-kill: 119 vs live 57 TPA).
#
# Golden = live FVSsn .sum. We assert the FIRE-YEAR rows bit-exact (the rows immediately at/after each
# fire — the fire timing + severity). Later post-fire cycles carry the separately-documented post-fire
# DG residual (fire_burn/early ~4% Bdft), so they are not asserted here.
# =============================================================================

using Test
using FVSjl

# pull a single stand's .sum rows as Dict(year => Vector{Int} of the first 12 columns)
function _simfire_rows(key)
    txt = FVSjl.run_keyfile(key; variant = FVSjl.Southern(), output = :sum)
    rows = Dict{Int,Vector{Int}}()
    for ln in split(txt, '\n')
        t = split(strip(ln))
        length(t) >= 12 || continue
        y = tryparse(Int, t[1])
        (y === nothing || y < 1900) && continue
        rows[y] = [parse(Int, t[i]) for i in 1:12 if tryparse(Int, t[i]) !== nothing][1:min(12, end)]
    end
    return rows
end

@testset "SIMFIRE date-default + multi-fire scheduling (D9) vs live FVSsn" begin
    sdir = joinpath(@__DIR__, "..", "harness", "scenarios")

    # (1) s10_fire — bare `SIMFIRE` (no params) ⇒ fires in cycle 1 (1990→1995). The 1995 row is the
    #     fire-year row; live = TPA 57, BA 33, SDI 59, CCF 64, TopHt 63→66, QMD 10.3, Tcuft 777.
    k1 = joinpath(sdir, "s10_fire.key")
    if isfile(k1)
        r = _simfire_rows(k1)
        @test haskey(r, 1995)
        if haskey(r, 1995)
            @test r[1995][3] == 57     # TPA — the fire killed the stand to 57 (was 507 unburned)
            @test r[1995][4] == 33     # BA
            @test r[1995][5] == 59     # SDI
            @test r[1995][6] == 64     # CCF
        end
    else
        @test_skip "s10_fire.key not available"
    end

    # (2) fire_repeat — SIMFIRE 2000 + SIMFIRE 2020 ⇒ BOTH fire. The 2005 row (after the 2000 fire) is
    #     bit-exact vs live: TPA 113, BA 73, SDI 126, CCF 139. (The 2020 fire also fires — TPA drops at
    #     2025 — but 15 yr of post-fire DG drift makes that row carry the documented residual.)
    k2 = joinpath(sdir, "fire_repeat.key")
    if isfile(k2)
        r = _simfire_rows(k2)
        @test haskey(r, 2005)
        if haskey(r, 2005)
            @test r[2005][3] == 113    # TPA after the FIRST (2000) fire — bit-exact
            @test r[2005][4] == 73     # BA
            @test r[2005][5] == 126    # SDI
            @test r[2005][6] == 139    # CCF
        end
        # the SECOND fire (2020) must also fire: 2025 TPA drops well below the unburned ~108 trajectory.
        @test haskey(r, 2025)
        # 2020 fire fired: jl 2025 TPA 64 vs CONFIRMED live 66 (ran /tmp/FVSsn_new on fire_repeat.key, 2026-07-05;
        # 2005 bit-exact 113). Δ2 = the 15-yr post-fire DG-drift + fire-kill-distribution residual on a rendered
        # integer. Was a loose `<= 70` regime threshold (would pass anything below the ~108 unburned trajectory).
        haskey(r, 2025) && @test abs(r[2025][3] - 66) <= 2
    else
        @test_skip "fire_repeat.key not available"
    end
end
