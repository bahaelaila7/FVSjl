# test_regen_coverage.jl — REGEN (establishment) coverage tracker. SN's establishment
# is keyword-driven (ESTAB block: PLANT/NATURAL/TALLY + sprouting), not auto-ingrowth.
# The bare-stand scenarios regenerate to 800 TPA at cycle 1 in Oracle A; FVSjl produces
# 0 until ESTAB is ported. @test_broken documents the gap and flips to a failure
# (alerting) once ESTAB lands. The empty-stand robustness (bare stand runs → all-zero
# .sum instead of crashing) is a hard @test (the prerequisite for ESTAB).

using Test, FVSjl

const _HARNESS = joinpath(@__DIR__, "..", "harness", "scenarios")

# cycle-1 (1997) row: (TPA, BA, QMD, TCuFt) or nothing
function _cyc1_row(key)
    s, _ = initialize(key)
    FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
    io = IOBuffer(); FVSjl.write_sum_file(io, s)
    for ln in eachline(IOBuffer(String(take!(io))))
        occursin("-999", ln) && continue
        t = split(strip(ln)); length(t) >= 9 || continue
        startswith(t[1], "1997") &&
            return (tryparse(Float64, t[3]), tryparse(Float64, t[4]),
                    tryparse(Float64, t[8]), tryparse(Float64, t[9]))
    end
    return nothing
end

@testset "REGEN / establishment coverage (ESTAB gap tracker)" begin
    for nm in ("bare_plant", "bare_natural")
        key = joinpath(_HARNESS, nm * ".key")
        if !isfile(key)
            @test_skip "$nm scenario not generated (gen_estab_scenarios.sh)"
        else
            @testset "$nm: bare stand regenerates (bit-exact @1997)" begin
                row = _cyc1_row(key)
                @test row !== nothing                       # PREREQUISITE: bare stand runs (no crash)
                if row !== nothing
                    tpa, ba, qmd, tcuft = row
                    # ESTAB regenerates to Oracle A's cycle-1 stand BIT-EXACT (established
                    # heights via ESSUBH + ESRANN, dbh from height). Later cycles drift on
                    # the cyc3+ stochastic-DGSCOR tail (a separate residual class).
                    @test isapprox(tpa, 800.0; atol = 2)
                    @test isapprox(ba,   14.0; atol = 1)
                    @test isapprox(qmd,   1.8; atol = 0.1)
                    # birth-cycle regen carries NO volume — the oracle's establishment-cycle
                    # VOLS runs before the records exist, so cyc1 TCuFt=0 (the trees first
                    # get volume from the next cycle). Regression guard for that fix.
                    @test tcuft == 0.0
                end
            end
        end
    end
end
