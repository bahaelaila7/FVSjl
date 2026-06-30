# =============================================================================
# test_r8clark_special.jl — R8 Clark COEFFSO%DIB17 special species (D7).
#
# Baldcypress (221), pondcypress (222), and green ash (544) are the three species
# for which live FVS r8prep.f SKIPS the form-class `(FCLSS−AFI)/BFI` step (the
# `IF(SPEC.NE.221.AND..NE.222.AND..NE.544)` gate), so their secondary-coefficient
# DIB17 falls through to the `:507` floor = COEFFS%DIB17 (the raw dib17). jl had
# computed `(dib17−AFI)/BFI` for all species → over-extracted merch/saw/board cubic
# (all_GA cyc0 Bdft 223 vs live 174). The fix (r8clark_vol.jl) makes Mcuft/Scuft/Bdft
# bit-exact for these species while leaving every other species unchanged (the :507
# floor is a no-op when BFI<1). Golden = live FVSsn (homogeneous all-one-species stands).
# =============================================================================

using Test
using FVSjl

@testset "R8 Clark COEFFSO%DIB17 — cypress/green-ash merch volume (D7) vs live FVSsn" begin
    sdir = joinpath(@__DIR__, "..", "harness", "scenarios")
    # cyc0 (1990) volume columns: Tcuft, Mcuft, Scuft, Bdft — bit-exact vs live.
    golden = Dict(
        "all_GA" => (1253, 900,  47,  174),   # green ash (544)
        "all_PC" => (1600, 1026, 287, 861),   # pondcypress (222)
        "all_BY" => (1466, 1129, 377, 1362),  # baldcypress (221)
    )
    for (scn, (tc, mc, sc, bd)) in sort(collect(golden))
        key = joinpath(sdir, scn * ".key")
        tre = joinpath(sdir, scn * ".tre")
        if !isfile(key) || !isfile(tre)
            @info "$scn scenario not generated (gen_species_scenarios.sh); skipping"
            continue
        end
        sumtxt = FVSjl.run_keyfile(key; variant = Southern(), output = :sum)
        row = nothing
        for ln in split(sumtxt, '\n')
            t = split(strip(ln))
            if length(t) >= 12 && t[1] == "1990"; row = t; break; end
        end
        @test row !== nothing
        if row !== nothing
            @test parse(Int, row[9])  == tc    # Tcuft (was already bit-exact)
            @test parse(Int, row[10]) == mc    # Mcuft
            @test parse(Int, row[11]) == sc    # Scuft
            @test parse(Int, row[12]) == bd    # Bdft
        end
    end
end
