# test_multicycle.jl — multi-cycle regression vs a Fortran-oracle golden fixture.
#
# Closes the "automated suite only checks cyc0/cyc1" blindspot: the deep cycle-by-
# cycle behaviour (where BAMAX, self-thinning and the DGSCOR tail live) is now in
# Pkg.test(). Tolerances are loose enough to absorb the known untripled-cycle
# serial-correlation tail (~1-2% on volume) but tight enough to catch a gross
# regression like the pre-BAMAX ~20% TPA/BA overshoot. Golden values are the
# FVSjulia oracle's per-cycle .sum rows (golden_multicycle.csv), which match the
# Fortran baseline on these standard/numeric scenarios.
using Test
using FVSjl

const _GOLD = joinpath(@__DIR__, "golden_multicycle.csv")
const _SCEN = joinpath(@__DIR__, "..", "harness", "scenarios")

# scenario => Vector of (cycle,tpa,ba,sdi,qmd,tcuft)
function _load_golden(path)
    g = Dict{String,Vector{NTuple{6,Float64}}}()
    for (i, ln) in enumerate(eachline(path))
        i == 1 && continue
        f = split(strip(ln), ',')
        isempty(f[1]) && continue
        push!(get!(g, f[1], Vector{NTuple{6,Float64}}()),
              (parse.(Float64, f[2:7])...,))
    end
    return g
end

@testset "multi-cycle regression vs Fortran-oracle golden" begin
    if !isfile(_GOLD) || !isdir(_SCEN)
        @test_skip "golden fixture or scenarios not available"
    else
        gold = _load_golden(_GOLD)
        for (scn, rows) in sort(collect(gold); by = first)
            key = joinpath(_SCEN, scn * ".key")
            isfile(key) || (@test_skip "$scn.key missing"; continue)
            s, _ = initialize(key)
            notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_forest_type!(s)
            FVSjl.compute_volumes!(s)
            g = s.plot.gross_space
            # Tolerances set to the MEASURED jl-vs-golden floor across ALL scenarios (incl. the
            # THINBTA/ATA/BBA/ABA thins, whose RDPSRT residual is now ≤0.6 TPA — the prior 6%/atol=10
            # carve-out was stale, and the "≤0.6 TPA / ≤0.06% cuft" claim under-stated the real cuft 0.4%).
            # Measured maxima: TPA 1.02 (0.74%), BA 0.69, SDI 0.56, QMD 0.09, cuft 18 (0.4%). Bounds ≈1.5×
            # that — single-precision floor; a real regression now fails.
            tT, rT = 1.6, 0.011
            tQ     = 0.13
            tC, rC = 27.0, 0.007
            @testset "$scn" begin
                for (cyc, tpa, ba, sdi, qmd, tcuft) in rows
                    FVSjl.compute_forest_type!(s)
                    mtpa = stand_tpa(s) / g; mba = stand_ba(s) / g
                    msdi = stand_sdi(s) / g; mqmd = stand_qmd(s)
                    mtcuft = FVSjl.summary_row(s; period = 0).cuft
                    @test isapprox(mba,  ba;  atol = 1.1, rtol = 0.011)
                    @test isapprox(mtpa, tpa; atol = tT, rtol = rT)
                    @test isapprox(msdi, sdi; atol = 1.0, rtol = 0.009)
                    @test isapprox(mqmd, qmd; atol = tQ)
                    @test isapprox(mtcuft, tcuft; atol = tC, rtol = rC)
                    Int(cyc) < 10 && FVSjl.grow_cycle!(s)
                end
            end
        end
    end
end
