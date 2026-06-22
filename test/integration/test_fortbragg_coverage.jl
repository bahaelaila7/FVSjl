# test_fortbragg_coverage.jl — Fort Bragg (forest 701) coverage. STDINFO 701xx maps
# to IFOR=20 (special longleaf/loblolly DG + bark) AND remaps KODFOR to NC Uwharrie
# 81110 (region 8) so VOLEQDEF assigns the R8 Clark equations. Without the KODFOR
# remap, VOLEQDEF sees region 7 (70106÷10000) and gives every tree ZERO volume — this
# guards that regression (the .sum cubic columns were all 0).

using Test, FVSjl

const _FB_HARNESS = joinpath(@__DIR__, "..", "harness", "scenarios")

# cycle-0 (inventory year) row: total cubic foot volume (TCuFt), or nothing
function _fb_cyc0_tcuft(key, invyear)
    s, _ = initialize(key)
    FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
    io = IOBuffer(); FVSjl.write_sum_file(io, s)
    for ln in eachline(IOBuffer(String(take!(io))))
        occursin("-999", ln) && continue
        t = split(strip(ln)); length(t) >= 9 || continue
        startswith(t[1], invyear) && return tryparse(Float64, t[9])
    end
    return nothing
end

@testset "Fort Bragg (forest 701) volume coverage" begin
    # forest 701 remaps to Uwharrie 81110 (region 8) for VOLEQDEF
    s, _ = initialize(joinpath(_FB_HARNESS, "s30_fortbragg_ll.key"))
    @test s.plot.forest_idx == 20                  # IFOR=20: special LL/LP DG + bark
    @test s.plot.user_forest_code == 81110         # KODFOR remapped to region 8

    for (nm, invyear, want) in (("s30_fortbragg_ll", "1990", 1421.0),
                                ("s31_fortbragg_lp", "1990", 1401.0))
        key = joinpath(_FB_HARNESS, nm * ".key")
        if !isfile(key)
            @test_skip "$nm scenario not generated"
        else
            tcuft = _fb_cyc0_tcuft(key, invyear)
            @test tcuft !== nothing
            # the headline guard: NONZERO volume (was 0 before the KODFOR remap)
            tcuft !== nothing && @test tcuft > 0.0
            # and close to Oracle A's cycle-0 TCuFt (R8 Clark cubic)
            tcuft !== nothing && @test isapprox(tcuft, want; atol = 5)
        end
    end
end
