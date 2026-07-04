# test_fix_scalers.jl — FIXDG / FIXHTG one-shot growth scalers vs live Fortran.
#
# FIXDG/FIXHTG (grincr.f:451-525) multiply the diameter / height increment by a factor for
# trees of a given species (0 = all) in a DBH window [d1, d2), in the SINGLE cycle whose
# year range contains the keyword date (one-shot — OPDONE marks it accomplished). The scaler
# also hits the tripled upper/lower records (DG(ITFN)/DG(ITFN+1)), so the port must scale the
# stashed dgU/dgL (htgU/htgL) too. It runs after all growth and before MORTS, so the reduced
# DG feeds the (D+G) the mortality model sees. Each .sum.save is live-Fortran output.
#   * fixdg_all  — FIXDG 0.3× all species, all DBH, in the 1995 cycle.
#   * fixdg_win  — FIXDG 0.5× only for 5"≤DBH<10" (exercises the window).
#   * fixhtg_all — FIXHTG 0.4× all species, all DBH.

using Test, FVSjl

const _FIX_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_fix_rows(txt) = [split(l) for l in split(txt, "\n")
                  if length(split(l)) >= 11 && tryparse(Int, first(split(l))) !== nothing]
_fix_base(path) = [split(l) for l in eachline(path)
                   if length(split(l)) >= 11 &&
                      (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_fcol(r, c) = parse(Float64, r[c])

@testset "FIXDG / FIXHTG one-shot scalers vs Fortran" begin
    have(nm) = isfile(joinpath(_FIX_DIR, nm * ".key")) && isfile(joinpath(_FIX_DIR, nm * ".sum.save"))
    runjl(nm) = (_fix_rows(FVSjl.run_keyfile(joinpath(_FIX_DIR, nm * ".key"); faithful = true)),
                 _fix_base(joinpath(_FIX_DIR, nm * ".sum.save")))

    for nm in ("fixdg_all", "fixdg_win", "fixhtg_all")
        if !have(nm); @test_skip "$nm scenario not available"; continue; end
        @testset "$nm" begin
            jl, ft = runjl(nm)
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for i in 1:length(jl), c in (3, 4, 7, 8)   # TPA / BA / TopHt / QMD
                    @test _fcol(jl[i], c) == _fcol(ft[i], c)
                end
            end
        end
    end

    # FIXDG 0.3× must visibly suppress growth (the 1995-cycle QMD lands below the base ≈ 7.0).
    if have("fixdg_all")
        jl, _ = runjl("fixdg_all")
        @test _fcol(jl[2], 8) < 6.8     # 2000 QMD, base ≈ 7.0
    end
end
