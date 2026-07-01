# =============================================================================
# test_growth_fint.jl — GROWTH FINT≠5 diameter-measurement period (D2) vs live FVSsn.
#
# The GROWTH keyword's field-1 FINT is the period (years) the input diameter increment was measured
# over (default 5). It feeds the DG calibration SCALE (= 5/FINT) AND the first projection cycle's
# serial-correlation `old` period (dgdriv PVMLT basis): AUTCOR(new = cycle length, old = FINT).
#
# jl previously hardcoded the first-cycle `old` to the variant native YR (htg_period), so a GROWTH 10
# run used AUTCOR(5,5) CORR=0.3196 instead of the correct AUTCOR(5,10) CORR=0.3906 — a ~0.4% cuft /
# ~1.2% bdft under-growth (the D2 residual). A live debug-stamp of dgdriv proved CORR=0.3906 and the
# calibration (COR/OLDRN/VARDG) was already bit-exact, localizing the miss to the projection `old` period.
#
# Fix (diameter_growth.jl): the first-cycle `old` = growth_fint when it is overridden from the universal
# 5-yr default, else htg_period. Default runs (growth_fint=5) are unchanged in BOTH variants. Golden =
# live FVSsn: growth_fint10 is now BIT-EXACT (1995 Tcuft 2848 / Bdft 11115, 2000 3308 / 13836).
# =============================================================================

using Test
using FVSjl

@testset "GROWTH FINT≠5 measurement period (D2) vs live FVSsn" begin
    key = joinpath(@__DIR__, "..", "harness", "scenarios", "growth_fint10.key")
    if !isfile(key)
        @test_skip "growth_fint10.key not available"
    else
        txt = FVSjl.run_keyfile(key; variant = FVSjl.Southern(), output = :sum)
        rows = Dict{Int,Vector{String}}()
        for ln in split(txt, '\n')
            t = split(strip(ln))
            length(t) >= 12 && (y = tryparse(Int, t[1]); y !== nothing && 1900 < y < 2100) && (rows[y] = t)
        end
        # cyc0 was already bit-exact; the fix makes the projection cycles bit-exact too.
        golden = Dict(1990 => (2295, 8566), 1995 => (2848, 11115), 2000 => (3308, 13836))
        for (yr, (tc, bd)) in sort(collect(golden))
            @test haskey(rows, yr)
            if haskey(rows, yr)
                @test parse(Int, rows[yr][9])  == tc    # Tcuft — was 2835 (0.46% low) at 1995 before the fix
                @test parse(Int, rows[yr][12]) == bd    # Bdft  — was 10977 (1.24% low) at 1995 before the fix
            end
        end
    end
end
