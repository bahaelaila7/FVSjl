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
            # Golden RE-GROUNDED to live FVSsn (sn_oracle.sh; was Oracle A, wrong by ~1 TPA on s12_phys_p221).
            # jl matches live to print-rounding on EVERY scenario incl. mix_lp_hi (measured maxima with the
            # per-cycle compute_forest_type!: TPA 0.57, BA 0.48, SDI 0.49, QMD 0.05, cuft 1.0). Uniform tight
            # bound = one print unit. (The earlier "LP-calibration tail" was a measurement artifact — a
            # tolerance-probe loop that omitted the per-cycle FORTYP recompute, which feeds diameter growth.)
            # atol = 1 print unit (golden fields are print-rounded integers; jl's internal float is exact,
            # so the max faithful gap is the print half-width + sub-ULP — measured ≤0.57/0.48/0.49/0.05/1.0
            # for TPA/BA/SDI/QMD/cuft). cuft tightened to 1 print unit (rtol dropped): the old 2.0+0.2% was
            # over-caution — the measured cuft max is exactly 1.0 across every scenario.
            tT, rT = 1.0, 0.0
            tB     = 1.0
            tS     = 1.0
            tQ     = 0.1
            tC, rC = 1.0, 0.0
            @testset "$scn" begin
                for (cyc, tpa, ba, sdi, qmd, tcuft) in rows
                    FVSjl.compute_forest_type!(s)
                    mtpa = stand_tpa(s) / g; mba = stand_ba(s) / g
                    msdi = stand_sdi(s) / g; mqmd = stand_qmd(s)
                    mtcuft = FVSjl.summary_row(s; period = 0).cuft
                    # BA + SDI: jl RENDERS to the golden's print-rounded integer exactly (measured di(jl)==golden
                    # every scenario/cycle) — compare the rendered integer `==` (doctrine's preferred form,
                    # stronger than the old atol=1.0 float bound). TPA stays a float knife-edge (di can differ by 1
                    # where the per-acre value straddles the +0.5 boundary — the growth-transcendental).
                    @test trunc(Int, mba + 0.5) == trunc(Int, ba + 0.5)     # BA — rendered-integer BIT-EXACT
                    @test isapprox(mtpa, tpa; atol = tT, rtol = rT)
                    @test trunc(Int, msdi + 0.5) == trunc(Int, sdi + 0.5)   # SDI — rendered-integer BIT-EXACT
                    @test isapprox(mqmd, qmd; atol = tQ)
                    @test isapprox(mtcuft, tcuft; atol = tC, rtol = rC)
                    Int(cyc) < 10 && FVSjl.grow_cycle!(s)
                end
            end
        end
    end
end
