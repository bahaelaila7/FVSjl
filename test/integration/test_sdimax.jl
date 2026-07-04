# test_sdimax.jl — SDIMAX keyword (option 89) → per-species SDImax + self-thinning bounds vs Fortran.
#
# SDIMAX sets a species' max SDI (SDIDEF, field 2) and/or the lower/upper self-thinning
# percents PMSDIL/PMSDIU (fields 5/6, stored as fractions ÷100). The SDImax-driven mortality
# (mortality.jl) then caps the stand to that density. Before the fix SDIMAX was unrecognized
# (silently ignored). Checks: SETS the plot fields; FIRES (≠ no-SDIMAX); FORTRAN bit-exact.

using Test, FVSjl

const _SX_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_sx_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_sx_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]

@testset "SDIMAX → per-species SDImax + self-thinning bounds" begin
    # 1. SETS — field 2 sets sp_sdi_def for all species (field 1 = 0); fields 5/6 → fractions
    s = FVSjl.StandState(FVSjl.Southern())
    flds = ["0", "300.", "", "", "20.", "90.", fill("", 6)...]
    vals = Float32[0, 300, 0, 0, 20, 90, zeros(Float32, 6)...]
    prs  = [true, true, false, false, true, true, falses(6)...]
    FVSjl.kw_sdimax!(s, FVSjl.KeywordRecord("SDIMAX  ", "", flds, vals, prs, 12, FVSjl.KW_OK, 0))
    @test all(==(300f0), s.plot.sp_sdi_def)
    @test s.plot.pct_sdimax_mort_lo ≈ 0.20f0   # PMSDIL 20% → fraction
    @test s.plot.pct_sdimax_mort_hi ≈ 0.90f0   # PMSDIU 90% → fraction

    key = joinpath(_SX_DIR, "sdimax.key"); sav = joinpath(_SX_DIR, "sdimax.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "sdimax scenario not available"
    else
        jl = _sx_rows(FVSjl.run_keyfile(key; faithful = true))
        # 2. FIRES — removing SDIMAX lets the stand grow denser
        offkey = joinpath(_SX_DIR, "_sdimax_off.key")
        cp(joinpath(_SX_DIR, "sdimax.tre"), joinpath(_SX_DIR, "_sdimax_off.tre"); force = true)
        open(offkey, "w") do io
            for l in eachline(key); startswith(strip(l), "SDIMAX  ") || println(io, l); end
        end
        try
            off = _sx_rows(FVSjl.run_keyfile(offkey; faithful = true))
            @test maximum(parse(Int, r[5]) for r in off) > maximum(parse(Int, r[5]) for r in jl)
        finally
            rm(offkey; force = true); rm(joinpath(_SX_DIR, "_sdimax_off.tre"); force = true)
        end
        # 3. FORTRAN — TPA/BA/SDI/CCF/TopHt/QMD match live Fortran (board feet ± Scribner noise)
        ft = _sx_base(sav)
        @test length(jl) == length(ft)
        for (j, f) in zip(jl, ft)
            for c in (3, 4, 5, 6, 7); @test parse(Int, j[c]) == parse(Int, f[c]); end
            @test parse(Float64, j[8]) == parse(Float64, f[8])
        end
    end
end
