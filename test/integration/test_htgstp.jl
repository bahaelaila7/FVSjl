# test_htgstp.jl — HTGSTOP / TOPKILL top-damage events (htgstp.f) vs live Fortran.
#
# Both keywords damage trees in a height window (HT1, HT2] with kill proportion
# PKIL = BACHLO(AVEPRB, STDPBR). The scenarios use STDPBR=0 and PRB=1, which makes the
# subsystem DETERMINISTIC (bachlo returns AVEPRB without drawing, and PRB≥1 skips the RANN
# escape) — so they validate the damage mechanics independent of RNG ordering:
#   * htgstop_det — HTGSTOP (act 110): scale the height increment by 0.5 in the 1995 cycle.
#   * topkill_det — TOPKILL (act 111): top-kill trees taller than 30' by 0.5 (height → TOPH,
#     NORMHT/ITRUNC for permanent broken tops, crown ratio cut). Recovery over later cycles
#     exercises the negative-ICR ("crown adjusted elsewhere") bypass in crown_ratio_update!.
# The stochastic path (STDPBR>0 / PRB<1) walks records in species-sorted IND1 order so the
# RANN/BACHLO stream lines up with FVS; that path is in place but validated separately.

using Test, FVSjl

const _HG_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_hg_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 11 && tryparse(Int, first(split(l))) !== nothing]
_hg_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 11 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_hcol(r, c) = parse(Float64, r[c])

@testset "HTGSTOP / TOPKILL top damage vs Fortran" begin
    have(nm) = isfile(joinpath(_HG_DIR, nm * ".key")) && isfile(joinpath(_HG_DIR, nm * ".sum.save"))
    runjl(nm) = (_hg_rows(FVSjl.run_keyfile(joinpath(_HG_DIR, nm * ".key"); faithful = true)),
                 _hg_base(joinpath(_HG_DIR, nm * ".sum.save")))

    for nm in ("htgstop_det", "topkill_det")
        if !have(nm); @test_skip "$nm scenario not available"; continue; end
        @testset "$nm" begin
            jl, ft = runjl(nm)
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for i in 1:length(jl), c in (3, 4, 7, 8)   # TPA / BA / TopHt / QMD
                    @test abs(_hcol(jl[i], c) - _hcol(ft[i], c)) <= 1
                end
            end
        end
    end

    # TOPKILL must visibly act: the 1995-cycle top-kill drops TopHt far below the base ≈ 64.
    if have("topkill_det")
        jl, _ = runjl("topkill_det")
        r2000 = findfirst(r -> r[1] == "2000", jl)
        @test r2000 !== nothing && _hcol(jl[r2000], 7) < 40   # post-topkill TopHt (base ≈ 64)
    end
end
