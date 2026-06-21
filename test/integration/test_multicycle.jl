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
            @testset "$scn" begin
                for (cyc, tpa, ba, sdi, qmd, tcuft) in rows
                    FVSjl.compute_forest_type!(s)
                    mtpa = stand_tpa(s) / g; mba = stand_ba(s) / g
                    msdi = stand_sdi(s) / g; mqmd = stand_qmd(s)
                    mtcuft = FVSjl.summary_row(s; period = 0).cuft
                    # BA is BAMAX-capped ⇒ kept TIGHT: it is the regression sentinel
                    # (a self-thinning/BAMAX regression blows BA). TPA/QMD/cuft carry
                    # the known untripled-cycle serial-correlation tail (up to ~9% on
                    # mixed stands at late cycles), so they are looser — still well
                    # below a gross regression (the pre-BAMAX overshoot was ~20%+).
                    @test isapprox(mba,  ba;  atol = 3, rtol = 0.05)   # sentinel
                    @test isapprox(mtpa, tpa; atol = 4, rtol = 0.10)
                    @test isapprox(msdi, sdi; atol = 6, rtol = 0.08)
                    @test isapprox(mqmd, qmd; atol = 0.7, rtol = 0.05)  # rtol for large late-cycle QMD (sparse post-thin stands)
                    @test isapprox(mtcuft, tcuft; atol = 40, rtol = 0.08)
                    Int(cyc) < 10 && FVSjl.grow_cycle!(s)
                end
            end
        end
    end
end
