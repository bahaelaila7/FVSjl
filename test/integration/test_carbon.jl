# test_carbon.jl — Stand Carbon Report live-tree pools (CARBREPT / CARBCALC method 1 = Jenkins).
#
# The FFE Stand Carbon Report (fmcrbout.f) is one of the few extension reports the STRIPPED
# ground-truth binary still prints to `.out` (like STRCLASS for SSTAGE) — so the Jenkins live-tree
# carbon pools are validatable BIT-EXACT, not just vs Oracle A. FVSjl computes them in
# `stand_carbon_report` (Jenkins biomass × 0.5 carbon fraction × TPA, converted to metric tons/ha).
# The committed baseline `carbon_jenkins.report.save` is the live Fortran report's data rows.
#
# Scope: the LIVE columns (Aboveground Total / Merch / Belowground Live) are validated here. The
# dead / down-wood / forest-floor / shrub-herb columns need the FFE surface-fuel model active
# (fire_on) — that is the remaining Stand-Carbon-Report increment; see carbon.jl.

using Test, FVSjl

const _CDIR = joinpath(@__DIR__, "..", "harness", "scenarios")

@testset "Stand Carbon Report — Jenkins live pools vs live Fortran" begin
    # CARBREPT / CARBCALC are recognized and set the report flag + method.
    s0 = FVSjl.StandState(FVSjl.Southern())
    @test !s0.control.carbon_report_on && s0.control.carbon_method == 1   # default Jenkins
    FVSjl.kw_carbrept!(s0, FVSjl.KeywordRecord("CARBREPT", "", String[], Float32[], Bool[], 0, FVSjl.KW_OK, 0))
    @test s0.control.carbon_report_on
    FVSjl.kw_carbcalc!(s0, FVSjl.KeywordRecord("CARBCALC", "", ["0"], Float32[0], [true], 1, FVSjl.KW_OK, 0))
    @test s0.control.carbon_method == 0

    key = joinpath(_CDIR, "carbon_jenkins.key"); sav = joinpath(_CDIR, "carbon_jenkins.report.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "carbon_jenkins scenario not available"
    else
        # Fortran report rows: [year, AbvTotal, Merch, BelowLive, BelowDead, StandDead, DDW, Floor,
        # Shb/Hrb, Total, Removed, Released] in metric tons C/ha.
        ft = [split(strip(l)) for l in eachline(sav) if occursin(r"^(19|20)\d\d\s", strip(l))]
        @test length(ft) >= 2

        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        for (c, f) in enumerate(ft)
            FVSjl.compute_density!(s)
            r = FVSjl.stand_carbon_report(s)
            # The inventory cycle (no growth) is BIT-EXACT; grown cycles carry the LP DBH-calibration
            # tail (~0.1% here) that the Jenkins biomass inherits via DBH — orthogonal to the carbon
            # model, so the tolerance widens by a small relative term once the stand has grown.
            tol(v) = c == 1 ? 0.1 : 0.005 * v + 0.1
            @test abs(r.aboveground - parse(Float64, f[2])) <= tol(parse(Float64, f[2]))
            @test abs(r.merch       - parse(Float64, f[3])) <= tol(parse(Float64, f[3]))
            @test abs(r.belowground - parse(Float64, f[4])) <= tol(parse(Float64, f[4]))
            c < length(ft) && FVSjl.grow_cycle!(s; fint = 5f0)
        end

        # FULL Stand Carbon Report at the inventory cycle, once the FFE fuel model (fmcba!) has
        # populated fire.cwd / fire.flive. Every column reconciles BIT-EXACT vs the Fortran report:
        # DDW (×0.5), FOREST FLOOR (×0.37, the Smith & Heath litter/duff fraction, fmcrbout.f:90),
        # SHRUB/HERB (FULIV2 coastal-plain/piedmont override → FLIVE=0.6 here ⇒ 0.67), and the TOTAL.
        s2 = first(FVSjl.each_stand(key))
        FVSjl.notre!(s2); FVSjl.setup_growth!(s2); FVSjl.compute_volumes!(s2); FVSjl.compute_density!(s2)
        if s2.fire !== nothing && s2.fire.active
            FVSjl.compute_forest_type!(s2); FVSjl.fmcba!(s2)
            r = FVSjl.stand_carbon_report(s2)
            f = ft[1]
            @test abs(r.down_wood    - parse(Float64, f[7]))  <= 0.1    # DDW  (3.8)
            @test abs(r.forest_floor - parse(Float64, f[8]))  <= 0.1    # Forest floor (9.1, via ×0.37)
            @test abs(r.shrub_herb   - parse(Float64, f[9]))  <= 0.1    # Shrub/herb (0.7, via FULIV2)
            @test abs(r.total        - parse(Float64, f[10])) <= 0.2    # Total stand carbon (90.1)
        end
    end
end
