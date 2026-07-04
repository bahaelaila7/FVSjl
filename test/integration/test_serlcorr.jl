# test_serlcorr.jl — SERLCORR keyword (option 91) → ARMA(1,1) DGSCOR phi/theta vs Fortran.
#
# SERLCORR sets BJPHI (AR, field 1) / BJTHET (MA, field 2) of the stochastic diameter-growth
# serial correlation (default 0.74 / 0.42). They define the BJRHO autocorrelation series that
# autcor uses for the per-cycle DG variance/covariance multipliers. Before the fix the keyword
# was unrecognized and BJPHI/BJTHET were compile-time consts. Checks: SETS the control fields;
# FIRES (≠ default); FORTRAN bit-exact; and the default path keeps the precomputed const.

using Test, FVSjl

const _SL_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_sl_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_sl_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]

@testset "SERLCORR → ARMA(1,1) DGSCOR phi/theta" begin
    # 1. SETS — defaults match grinit.f (0.74 / 0.42); the keyword overrides both
    s = FVSjl.StandState(FVSjl.Southern())
    @test s.control.dg_bjphi == 0.74f0 && s.control.dg_bjthet == 0.42f0
    FVSjl.kw_serlcorr!(s, FVSjl.KeywordRecord("SERLCORR", "", ["0.50", "0.30", fill("", 10)...],
                       Float32[0.5, 0.3, zeros(Float32, 10)...], [true, true, falses(10)...], 12, FVSjl.KW_OK, 0))
    @test s.control.dg_bjphi == 0.5f0 && s.control.dg_bjthet == 0.3f0
    # default series equals the precomputed const; a custom one differs
    @test FVSjl.dg_bjrho_series() == FVSjl.DG_BJRHO
    @test FVSjl.dg_bjrho_series(0.5f0, 0.3f0) != FVSjl.DG_BJRHO

    key = joinpath(_SL_DIR, "serlcorr.key"); sav = joinpath(_SL_DIR, "serlcorr.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "serlcorr scenario not available"
    else
        jl = _sl_rows(FVSjl.run_keyfile(key; faithful = true))
        # 2. FIRES — SERLCORR changes the stochastic projection
        offkey = joinpath(_SL_DIR, "_serl_off.key")
        cp(joinpath(_SL_DIR, "serlcorr.tre"), joinpath(_SL_DIR, "_serl_off.tre"); force = true)
        open(offkey, "w") do io
            for l in eachline(key); startswith(strip(l), "SERLCORR") || println(io, l); end
        end
        try
            off = _sl_rows(FVSjl.run_keyfile(offkey; faithful = true))
            @test any(jl[i][10] != off[i][10] for i in eachindex(jl))
        finally
            rm(offkey; force = true); rm(joinpath(_SL_DIR, "_serl_off.tre"); force = true)
        end
        # 3. FORTRAN — TPA/BA/SDI/CCF/TopHt/QMD match live Fortran
        ft = _sl_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            for c in (3, 4, 5, 6, 7); @test parse(Int, j[c]) == parse(Int, f[c]); end
            @test parse(Float64, j[8]) == parse(Float64, f[8])
        end
    end
end
