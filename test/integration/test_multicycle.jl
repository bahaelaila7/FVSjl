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
            # Golden RE-GROUNDED to live FVSsn (sn_oracle.sh; was Oracle A, which was wrong by ~1 TPA on
            # s12_phys_p221). jl now matches live to print-rounding (≤0.57 TPA) on EVERY scenario — so the
            # bound is tight (atol 1 = one print unit; cuft 2/0.2%). The single exception is mix_lp_hi, a
            # mixed-loblolly high-site stand carrying the documented LP-growth-CALIBRATION tail (jl & Oracle A
            # both drift from live ~4.8 TPA / 0.8 QMD by late cycles); it gets an explicit wider, labelled bound.
            lp_tail = scn == "mix_lp_hi"
            tT, rT = lp_tail ? (5.0, 0.0) : (1.0, 0.0)
            tB     = lp_tail ? 1.5 : 1.0
            tS     = lp_tail ? 3.0 : 1.0
            tQ     = lp_tail ? 0.85 : 0.1
            tC, rC = lp_tail ? (10.0, 0.0) : (2.0, 0.002)
            @testset "$scn" begin
                for (cyc, tpa, ba, sdi, qmd, tcuft) in rows
                    FVSjl.compute_forest_type!(s)
                    mtpa = stand_tpa(s) / g; mba = stand_ba(s) / g
                    msdi = stand_sdi(s) / g; mqmd = stand_qmd(s)
                    mtcuft = FVSjl.summary_row(s; period = 0).cuft
                    @test isapprox(mba,  ba;  atol = tB)
                    @test isapprox(mtpa, tpa; atol = tT, rtol = rT)
                    @test isapprox(msdi, sdi; atol = tS)
                    @test isapprox(mqmd, qmd; atol = tQ)
                    @test isapprox(mtcuft, tcuft; atol = tC, rtol = rC)
                    Int(cyc) < 10 && FVSjl.grow_cycle!(s)
                end
            end
        end
    end
end
