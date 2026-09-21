# test_multicycle.jl — multi-cycle regression vs a live-FVSsn golden fixture (10 SN scenarios, 11 cycles).
#
# Closes the "automated suite only checks cyc0/cyc1" blindspot for SN: BA/SDI/QMD every cycle plus the TPA and cuft
# series, all compared at PRINT PRECISION (the rendered .sum value; no tolerances). Golden values are live FVSsn
# per-cycle .sum rows (golden_multicycle.csv). This is a narrow SN smoke gate; the multi-variant / multi-regime /
# multi-output gate is test_tiered.jl (test/harness/tiered/README.md).
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
            # Golden RE-GROUNDED to live FVSsn (sn_oracle.sh). Every column is compared at PRINT PRECISION (the
            # rendered .sum value): BA/SDI/TPA/cuft as rendered integers, QMD at 1 decimal. No tolerances.
            @testset "$scn" begin
                tpa_pairs  = Tuple{Float64,Float64}[]   # (measured per-acre TPA, golden) per cycle
                cuft_pairs = Tuple{Float64,Float64}[]   # (measured cuft, golden) per cycle
                for (cyc, tpa, ba, sdi, qmd, tcuft) in rows
                    FVSjl.compute_forest_type!(s)
                    mtpa = stand_tpa(s) / g; mba = stand_ba(s) / g
                    msdi = stand_sdi(s) / g; mqmd = stand_qmd(s)
                    mtcuft = FVSjl.summary_row(s; period = 0).cuft
                    # BA + SDI: jl RENDERS to the golden's print-rounded integer exactly (measured di(jl)==golden
                    # every scenario/cycle) — compare the rendered integer `==` (doctrine's preferred form).
                    @test trunc(Int, mba + 0.5) == trunc(Int, ba + 0.5)     # BA — rendered-integer BIT-EXACT
                    @test trunc(Int, msdi + 0.5) == trunc(Int, sdi + 0.5)   # SDI — rendered-integer BIT-EXACT
                    @test round(Float64(mqmd); digits = 1) == qmd   # QMD — rendered 1-dec BIT-EXACT
                    # TPA + cuft at PRINT PRECISION (the .sum renders both as integers, like BA/SDI above): the
                    # golden holds the printed integer, so compare the rendered integers — an unrounded-vs-printed
                    # float compare was never satisfiable and parked 10 scenarios as @test_broken spuriously.
                    push!(tpa_pairs,  (Float64(trunc(Int, mtpa + 0.5)), Float64(trunc(Int, tpa + 0.5))))
                    push!(cuft_pairs, (Float64(trunc(Int, Float64(mtcuft) + 0.5)), Float64(trunc(Int, tcuft + 0.5))))
                    Int(cyc) < 10 && FVSjl.grow_cycle!(s)
                end
                # Rendered-integer TPA/cuft must match every cycle (anything that PRINTS differently fails).
                @test all(a == b for (a, b) in tpa_pairs)    # TPA  — rendered-integer
                @test all(a == b for (a, b) in cuft_pairs)   # cuft — rendered-integer
            end
        end
    end
end
