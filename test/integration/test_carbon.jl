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

            # GROWN-cycle FOREST FLOOR + DDW via the FFE annual fuel loop (FMCWD decay + FMCADD
            # litterfall + woody breakage, NYRS=1 per year, crown held at the cycle's start). BOTH
            # reconcile vs the Fortran 1995 report row — validating the decay/litterfall/breakage
            # coupling AND crown_biomass (foliage + the V2T/2000 woody fix, checked vs a Fortran XV dump).
            for _ in 1:5
                FVSjl.fmcwd!(s2, 1); FVSjl.fmcadd_litterfall!(s2); FVSjl.fmcadd_woody!(s2)
            end
            r95 = FVSjl.stand_carbon_report(s2)
            @test abs(r95.forest_floor - parse(Float64, ft[2][8])) <= 0.1   # 1995 Floor = 6.6
            @test abs(r95.down_wood    - parse(Float64, ft[2][7])) <= 0.1   # 1995 DDW   = 2.5
        end
    end
end

@testset "Stand Carbon Report — FFE fuel driver across grown cycles" begin
    # `ffe_fuel_update!` (the per-cycle FFE fuel driver: fmcba! + the annual decay/litterfall/breakage
    # loop) evolves the surface-fuel pools across grown cycles. Validated vs the multi-cycle Fortran
    # report (carbon_jenkins, 4 cycles): the live Jenkins pools track within the LP growth tail, and
    # the FOREST FLOOR reconciles every cycle. DOWN DEAD WOOD reconciles until tree mortality begins
    # (2000+, where the report's Stand-Dead/Below-Dead columns turn on) — the snag-debris falldown
    # (CWD2B) that feeds DDW from dying trees is the remaining FMCADD term.
    key = joinpath(_CDIR, "carbon_jenkins.key"); sav = joinpath(_CDIR, "carbon_jenkins.report.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "carbon_jenkins scenario not available"
    else
        ft = [split(strip(l)) for l in eachline(sav) if occursin(r"^(19|20)\d\d\s", strip(l))]
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        s.fire !== nothing && s.fire.active && FVSjl.compute_forest_type!(s)
        s.fire !== nothing && s.fire.active && FVSjl.fmcba!(s)
        for (c, f) in enumerate(ft)
            FVSjl.compute_density!(s)
            r = FVSjl.stand_carbon_report(s)
            # live Jenkins pools within the LP growth-calibration tail
            @test abs(r.aboveground - parse(Float64, f[2])) <= 0.02 * parse(Float64, f[2]) + 0.2
            @test abs(r.belowground - parse(Float64, f[4])) <= 0.02 * parse(Float64, f[4]) + 0.2
            # forest floor reconciles every cycle (decay + litterfall)
            @test abs(r.forest_floor - parse(Float64, f[8])) <= 0.2
            # Below-Dead (dead coarse roots, BIOROOT) reconciles bit-exact every cycle
            @test abs(r.belowground_dead - parse(Float64, f[5])) <= 0.1
            # DDW reconciles bit-exact before mortality (1990/1995). After, the snag-bole + CWD2B-crown
            # falldown feed it, but the net DDW carries a ~1.9 within-cycle residual (a missing addition
            # source, e.g. live-crown breakage from the growing canopy — previously masked by the
            # jenkins whole-tree snag double-count), so post-mortality DDW is bounded, not bit-exact.
            f[1] in ("1990", "1995") && @test abs(r.down_wood - parse(Float64, f[7])) <= 0.1
            f[1] in ("2000", "2005") && @test r.down_wood <= parse(Float64, f[7]) + 0.2  # never over the report
            # STAND-DEAD reconciles BIT-EXACT via the faithful snag STEM-VOLUME bole (cuft·V2T) + the
            # CWD2B crown-still-in-waiting (not whole-tree Jenkins) — validated vs the instrumented
            # Fortran SNGBOLE/SNGTOT dump (bole 3.72/3.26, crown 1.46/1.19, sum = the report's 5.2/4.5).
            f[1] in ("1990",) && @test r.standing_dead == 0f0
            f[1] in ("2000", "2005") && @test abs(r.standing_dead - parse(Float64, f[6])) <= 0.1
            TO = 0.90718474 / 0.40468564
            f[1] == "2000" && @test abs(FVSjl.snag_bole_carbon(s) * TO - 3.72) <= 0.05
            f[1] == "2005" && @test abs(FVSjl.snag_bole_carbon(s) * TO - 3.28) <= 0.05
            f[1] == "2000" && @test abs(FVSjl.snag_crown_carbon(s) * TO - 1.46) <= 0.05
            if c < length(ft)
                # evolve the fuels with the START-of-cycle crown (FVS records the crown at the END of
                # each cycle for the NEXT cycle's litterfall, fmmain.f:264), THEN grow the trees.
                FVSjl.ffe_fuel_update!(s, 5)
                FVSjl.grow_cycle!(s; fint = 5f0)
            end
        end
    end
end
