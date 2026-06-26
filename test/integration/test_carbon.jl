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
            @test abs(r.down_wood    - parse(Float64, f[7]))  <= 0.05   # DDW  (3.8)        — BIT-EXACT
            @test abs(r.forest_floor - parse(Float64, f[8]))  <= 0.05   # Forest floor (9.1) — BIT-EXACT
            @test abs(r.shrub_herb   - parse(Float64, f[9]))  <= 0.05   # Shrub/herb (0.7)  — BIT-EXACT
            @test abs(r.total        - parse(Float64, f[10])) <= 0.05   # Total stand carbon — BIT-EXACT

            # GROWN-cycle FOREST FLOOR + DDW via the FFE annual fuel loop (FMCWD decay + FMCADD
            # litterfall + woody breakage, NYRS=1 per year, crown held at the cycle's start). BOTH
            # reconcile vs the Fortran 1995 report row — validating the decay/litterfall/breakage
            # coupling AND crown_biomass (foliage + the V2T/2000 woody fix, checked vs a Fortran XV dump).
            for _ in 1:5
                FVSjl.fmcwd!(s2, 1); FVSjl.fmcadd_litterfall!(s2); FVSjl.fmcadd_woody!(s2)
            end
            r95 = FVSjl.stand_carbon_report(s2)
            @test abs(r95.forest_floor - parse(Float64, ft[2][8])) <= 0.05  # 1995 Floor = 6.6 — BIT-EXACT
            @test abs(r95.down_wood    - parse(Float64, ft[2][7])) <= 0.05  # 1995 DDW   = 2.5 — BIT-EXACT
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
        # NB carbon_jenkins is a NON-bit-exact-GROWTH fixture (a synthetic LP stand with a diameter-growth
        # calibration tail). The bit-exact carbon validation is owned by carbon_snt (bit-exact growth, all
        # live pools Δ=0); here the LIVE pools track only within the growth tail, so the live-pool agreement
        # is NOT asserted numerically (it would not be ULP — the divergence is the growth, not the carbon).
        # What IS bit-exact here and asserted at ULP: the dead/floor pools that the FFE fuel model drives,
        # plus the snag-bole/crown split vs an instrumented Fortran dump.
        for (c, f) in enumerate(ft)
            FVSjl.compute_density!(s)
            r = FVSjl.stand_carbon_report(s)
            # forest floor + below-dead (dead coarse roots, BIOROOT) reconcile at print resolution
            @test abs(r.forest_floor     - parse(Float64, f[8])) <= 0.06   # tiny litterfall growth-tail effect
            @test abs(r.belowground_dead - parse(Float64, f[5])) <= 0.05   # BIT-EXACT
            # DDW: BIT-EXACT before mortality (1990/1995); the post-mortality dead-pool flow has the same
            # known crown-lift-timing gap as carbon_snt (~1.9) — tracked honestly as @test_broken, not hidden.
            f[1] in ("1990", "1995") && @test abs(r.down_wood - parse(Float64, f[7])) <= 0.05
            f[1] in ("2000", "2005") && @test_broken abs(r.down_wood - parse(Float64, f[7])) <= 0.05
            # STAND-DEAD: 0 before mortality; BIT-EXACT after, via the faithful snag STEM-VOLUME bole
            # (cuft·V2T) + the CWD2B crown-in-waiting — validated vs the instrumented Fortran SNGBOLE/SNGTOT.
            f[1] in ("1990",) && @test r.standing_dead == 0f0
            f[1] in ("2000", "2005") && @test abs(r.standing_dead - parse(Float64, f[6])) <= 0.05
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

@testset "Stand Carbon Report — .out writer (CARBREPT) byte-exact vs Fortran" begin
    # The FFE Stand Carbon Report .out writer (write_carbon_report, fmcrbout.f FORMATs 700-800): the
    # fixed header block + the per-row FORMAT must be byte-for-byte vs the Fortran .out, and the pool
    # values are the validated metric-tons/ha pools. The INVENTORY row (1990) is bit-exact in every
    # column; grown rows track within the LP growth tail (aboveground ~0.5%) with the post-mortality
    # DDW/Shrub bounded (the documented crown-lift/FLIVE residual).
    key = joinpath(_CDIR, "carbon_jenkins.key"); sav = joinpath(_CDIR, "carbon_jenkins.report.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "carbon_jenkins scenario not available"
    else
        ft = [strip(l) for l in eachline(sav) if occursin(r"^(19|20)\d\d\s", strip(l))]
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        io = IOBuffer()
        FVSjl.write_carbon_report(io, s, length(ft) - 1; stand_id = "CARBJENK", mgmt_id = "NONE")
        out = split(String(take!(io)), '\n')

        # header block byte-exact (title, units, stand-id, the three column-header lines, separators)
        @test out[1] == "-"^110
        @test any(l == "                                         STAND CARBON REPORT (BASED ON STOCKABLE AREA)" for l in out)
        @test any(l == "STAND ID: CARBJENK                      MGMT ID: NONE" for l in out)
        @test any(l == "YEAR    Total    Merch     Live     Dead     Dead      DDW    Floor  Shb/Hrb   Carbon   Carbon  from Fire" for l in out)

        rows = [l for l in out if occursin(r"^(19|20)\d\d ", l)]
        @test length(rows) == length(ft)
        # INVENTORY row byte-for-byte identical to the Fortran report row (every column)
        @test strip(rows[1]) == ft[1]                       # INVENTORY row byte-exact — BIT-EXACT
        # Grown rows: the live aboveground/merch agreement is the GROWTH tail of this synthetic fixture
        # (NOT a carbon property — carbon_snt validates those bit-exact), so it is not asserted at ULP here.
        # DDW carries the same known post-mortality dead-pool gap as carbon_snt → @test_broken, not hidden.
        for (i, r) in enumerate(rows)
            i == 1 && continue
            mv = parse.(Float64, split(strip(r))); fv = parse.(Float64, split(ft[i]))
            if mv[1] < 2000                                 # pre-mortality: DDW is BIT-EXACT
                @test abs(mv[7] - fv[7]) <= 0.05
            else                                            # post-mortality: known dead-pool flow gap
                @test_broken abs(mv[7] - fv[7]) <= 0.05
            end
        end
    end
end

@testset "Stand Carbon Report — LIVE pools bit-exact on a bit-exact-growth FFE stand (carbon_snt)" begin
    # carbon_snt = snt01_alpha's bit-exact species + FMIN/CARBREPT (no LP growth tail). On bit-exact
    # growth the LIVE carbon pools (Aboveground Total/Merch, Belowground-Live) and the Forest Floor
    # reconcile BIT-EXACT vs live Fortran every cycle — proving the live-carbon model is correct (the
    # carbon_jenkins ~0.5% residuals were purely its growth tail, not a carbon bug). The DEAD columns
    # (Below-Dead, Stand-Dead, DDW) and Shrub/Herb differ by the three remaining FFE pieces: FMSSEE
    # (input-snag seeding → the inventory Stand-Dead/Below-Dead), the crown-lift (DDW), and FLIVE
    # live-fuel growth (Shb/Hrb) — see FFE_FUEL_DYNAMICS_chunk_plan.md.
    key = joinpath(_CDIR, "carbon_snt.key"); sav = joinpath(_CDIR, "carbon_snt.report.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "carbon_snt scenario not available"
    else
        ft = [parse.(Float64, split(strip(l))) for l in eachline(sav) if occursin(r"^(19|20)\d\d\s", strip(l))]
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        io = IOBuffer()
        FVSjl.write_carbon_report(io, s, length(ft) - 1; stand_id = "S248112", mgmt_id = "NONE")
        rows = [parse.(Float64, split(strip(l)))
                for l in split(String(take!(io)), '\n') if occursin(r"^(19|20)\d\d ", l)]
        @test length(rows) == length(ft)
        # The report prints F7.1 (one decimal), so a bit-exact column matches to the print resolution 0.05.
        for (mv, fv) in zip(rows, ft)
            @test mv[2] ≈ fv[2] atol = 0.05    # Aboveground Total — BIT-EXACT
            @test mv[3] ≈ fv[3] atol = 0.05    # Merch             — BIT-EXACT
            @test mv[4] ≈ fv[4] atol = 0.05    # Belowground Live  — BIT-EXACT
            @test mv[8] ≈ fv[8] atol = 0.05    # Forest Floor      — BIT-EXACT
            @test mv[9] ≈ fv[9] atol = 0.05    # Shrub/Herb (FLIVE) — BIT-EXACT (post-grow FLIVE refresh)
        end
        # DEAD POOLS (BelowD/StandD/DDW) are NOT yet bit-exact — tracked honestly as @test_broken on the
        # max residual across cycles (not hidden behind a loose passing tolerance). The dead-pool flow has a
        # known intermediate-cycle gap (the crown-lift is applied in the next cycle's fuel loop while FVS
        # applies it same-cycle; the inventory + final cycles ARE bit-exact). These @test_broken will flip
        # to a (passing) failure if the dead pools ever reconcile — see docs/FFE_FUEL_DYNAMICS_chunk_plan.md.
        maxd(c) = maximum(abs(rows[i][c] - ft[i][c]) for i in 1:length(ft))
        @test maxd(5) <= 0.05              # Belowground Dead — BIT-EXACT (input-snag root XDCAY = (1−CRDCAY)^10)
        @test_broken maxd(6) <= 0.05       # Stand Dead       (max Δ ≈ 0.6, CWD2B crown-flow phasing)
        @test_broken maxd(7) <= 0.05       # DDW              (max Δ ≈ 1.3, CWD2B/decay phasing)
    end
end

@testset "Stand Carbon Report — emitted by the LIVE run_keyfile path (CARBREPT integration)" begin
    # The CARBREPT keyword (inside the FMIN block) must drive the carbon report from the SAME main
    # simulation (write_sum_file) — not a separate re-simulation — and append it to the .out. This proves
    # the drop-in path: a CARBREPT .key produces the report in a normal run. Values must equal the
    # standalone write_carbon_report (same single simulation) and track the Fortran .out.
    key = joinpath(_CDIR, "carbon_snt.key"); sav = joinpath(_CDIR, "carbon_snt.report.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "carbon_snt scenario not available"
    else
        out = FVSjl.run_keyfile(key)
        @test occursin("STAND CARBON REPORT", out)                 # the report block is in the .out
        @test occursin("YEAR    Total    Merch", out)              # column header present
        ft = Dict(parse(Int, split(strip(l))[1]) => parse.(Float64, split(strip(l)))
                  for l in eachline(sav) if occursin(r"^(19|20)\d\d\s", strip(l)))
        # slice from the carbon-report header so we don't pick up the .sum table rows (also year-led)
        olines = split(out, '\n')
        ci = findfirst(l -> occursin("STAND CARBON REPORT", l), olines)
        rows = [parse.(Float64, split(strip(l)))
                for l in olines[ci:end] if occursin(r"^(19|20)\d\d ", strip(l)) &&
                    haskey(ft, parse(Int, split(strip(l))[1])) && length(split(strip(l))) >= 11]
        @test !isempty(rows)
        for mv in rows
            fv = ft[Int(mv[1])]
            @test mv[2] ≈ fv[2] atol = 0.05         # Aboveground Total — BIT-EXACT
            @test mv[4] ≈ fv[4] atol = 0.05         # Belowground Live  — BIT-EXACT
        end
        # DDW post-mortality dead-pool flow gap (same as the writer test) — broken, not hidden behind slack.
        @test_broken maximum(abs(mv[7] - ft[Int(mv[1])][7]) for mv in rows) <= 0.05
    end
end

@testset "CARBCALC=0 FFE-fuel live carbon method (semantic; live oracle binary-blocked)" begin
    # The FFE carbon method (CARBCALC=0) computes Above = crown(foliage+woody) + stem, Merch = stem
    # (fmcrbout.f:120-141 + fmdout.f:225-258), vs JENKINS national biomass. The live FFE oracle is
    # degenerate in the stripped validation binary (zero live pools), so this is a SOURCE-FAITHFUL port
    # validated by its semantics: Above = crown+stem > Merch = stem > 0, in the Jenkins ballpark, and the
    # method only changes Above/Merch (roots + all dead pools are identical to JENKINS, fmcrbout.f:144).
    key = joinpath(_CDIR, "carbon_snt.key")
    if !isfile(key)
        @test_skip "carbon_snt scenario not available"
    else
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        FVSjl.compute_forest_type!(s); FVSjl.fmcba!(s)
        ffe = FVSjl.ffe_live_carbon(s); jen = FVSjl.stand_live_carbon(s)
        @test ffe.merch > 0f0                                  # stem biomass present
        @test ffe.aboveground > ffe.merch                     # Above = crown + stem  >  Merch = stem
        @test 0.5 < ffe.aboveground / jen.aboveground < 2.0   # same order of magnitude as Jenkins
        # the report branch: method 0 changes Above/Merch but NOT belowground / dead pools
        s.control.carbon_method = Int32(0); r0 = FVSjl.stand_carbon_report(s)
        s.control.carbon_method = Int32(1); r1 = FVSjl.stand_carbon_report(s)
        @test r0.aboveground != r1.aboveground                # method actually switches the basis
        @test r0.belowground ≈ r1.belowground                 # roots unchanged (Jenkins in both)
        @test r0.standing_dead ≈ r1.standing_dead             # dead pools unchanged
        @test r0.down_wood ≈ r1.down_wood
    end
end

using SQLite, DBInterface

@testset "FVS_Carbon DBS table writer (dbsfmcrpt.f schema + round-trip)" begin
    # write_dbs_carbon! writes the FFE Stand Carbon Report pools to the FVS_Carbon SQLite table (the same
    # metric-tons/ha values as the .out report). Validate the schema + the report→column mapping round-trip.
    key = joinpath(_CDIR, "carbon_snt.key")
    if !isfile(key)
        @test_skip "carbon_snt scenario not available"
    else
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        FVSjl.compute_forest_type!(s); FVSjl.fmcba!(s)
        rep = FVSjl.stand_carbon_report(s)
        rows = [(1990, rep), (1995, rep)]
        dbpath = joinpath(mktempdir(), "carb.db")
        FVSjl.write_dbs_carbon!(dbpath, "CASE1", "STAND1", rows)
        db = SQLite.DB(dbpath)
        res = [(; Year = r.Year, CaseID = r.CaseID, StandID = r.StandID,
                Above = r.Aboveground_Total_Live, SD = r.Standing_Dead,
                DDW = r.Forest_Down_Dead_Wood, Total = r.Total_Stand_Carbon)
               for r in DBInterface.execute(db, "SELECT * FROM FVS_Carbon ORDER BY Year")]
        SQLite.close(db)
        @test length(res) == 2                                            # two cycles inserted
        @test res[1].Year == 1990 && res[2].Year == 1995
        @test res[1].CaseID == "CASE1" && res[1].StandID == "STAND1"
        @test res[1].Above ≈ Float64(rep.aboveground)                     # report→column mapping
        @test res[1].SD ≈ Float64(rep.standing_dead)
        @test res[1].DDW ≈ Float64(rep.down_wood)
        @test res[1].Total ≈ Float64(rep.total)
    end
end

@testset "FVS_Fuels DBS table — loadings grounded in validated DDW/Floor + round-trip (dbsfuels.f)" begin
    # ffe_fuel_loadings (DBSFUELS inputs, tons/ac biomass) is composed of the SAME FFE pools that give the
    # validated carbon report: surface woody lt3+ge3 = the down-wood biomass (= down_wood_carbon/0.5),
    # litter+duff = the forest-floor biomass. So the loadings are value-grounded; here we check that
    # relationship + the FVS_Fuels schema/round-trip.
    key = joinpath(_CDIR, "carbon_snt.key")
    if !isfile(key)
        @test_skip "carbon_snt scenario not available"
    else
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        FVSjl.compute_forest_type!(s); FVSjl.fmcba!(s)
        f = FVSjl.ffe_fuel_loadings(s)
        # down-wood biomass (lt3+ge3) reconciles with the validated DDW carbon pool (DDW = biomass × 0.5)
        @test (f.lt3 + f.ge3) ≈ FVSjl.down_wood_carbon(s) / 0.5f0 rtol = 1f-4
        # litter+duff reconcile with the forest-floor pool (floor carbon = biomass × 0.37)
        @test (f.litter + f.duff) ≈ FVSjl.forest_floor_carbon(s) / 0.37f0 rtol = 1f-4
        @test f.s3to6 + f.s6to12 + f.ge12 ≈ f.ge3 rtol = 1f-4      # ge3 size split is consistent
        @test f.surf_total > 0f0
        # DBS round-trip (collection is (year, carbon_report, fuel) 3-tuples)
        rows = [(1990, FVSjl.stand_carbon_report(s), f)]
        dbpath = joinpath(mktempdir(), "fuels.db")
        FVSjl.write_dbs_fuels!(dbpath, "C1", "S1", rows)
        db = SQLite.DB(dbpath)
        res = [(; Year = r.Year, Litter = r.Surface_Litter, lt3 = r.Surface_lt3,
                ST = r.Surface_Total, Tot = r.Total_Biomass)
               for r in DBInterface.execute(db, "SELECT * FROM FVS_Fuels")]
        SQLite.close(db)
        @test length(res) == 1 && res[1].Year == 1990
        @test res[1].lt3 ≈ Float64(f.lt3)
        @test res[1].ST ≈ Float64(f.surf_total)
    end
end

@testset "FVS_SnagSum DBS table — snag density by hard/soft × DBH class (dbsfmssnag.f)" begin
    # snag_summary maps FVSjl's per-record den_hard/den_soft into the FMSSUM cumulative DBH classes
    # (SNPRCL 0/12/18/24/30/36). class1 (≥0) equals the total; the DBS round-trip preserves it.
    key = joinpath(_CDIR, "carbon_snt.key")
    if !isfile(key)
        @test_skip "carbon_snt scenario not available"
    else
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        s.fire !== nothing && s.fire.active && FVSjl.ffe_seed_input_snags!(s)
        sg = FVSjl.snag_summary(s)
        @test sg.hard[1] ≈ sg.hard[7]                         # class 1 (≥0") == hard total
        @test sg.soft[1] ≈ sg.soft[7]                         # class 1 == soft total
        @test sg.hard[7] + sg.soft[7] ≈ FVSjl.snag_standing_density(s.fire) rtol = 1f-4
        @test sg.hard[2] <= sg.hard[1]                        # cumulative: ≥12" ⊆ ≥0"
        rows = [(1990, FVSjl.stand_carbon_report(s), FVSjl.ffe_fuel_loadings(s), sg)]
        dbpath = joinpath(mktempdir(), "snag.db")
        FVSjl.write_dbs_snagsum!(dbpath, "C1", "S1", rows)
        db = SQLite.DB(dbpath)
        res = [(; Year = r.Year, H1 = r.Hard_snags_class1, HT = r.Hard_snags_total,
                Tot = r.Hard_soft_snags_total)
               for r in DBInterface.execute(db, "SELECT * FROM FVS_SnagSum")]
        SQLite.close(db)
        @test length(res) == 1 && res[1].Year == 1990
        @test res[1].H1 ≈ Float64(sg.hard[1])
        @test res[1].Tot ≈ Float64(sg.hard[7] + sg.soft[7])
    end
end

@testset "FVS_Down_Wood_Vol/Cov DBS tables — volume & cover from validated cwd (dbsfmdwvol/cov.f)" begin
    # ffe_down_wood derives volume (cwd biomass·2000/density) + cover (a·vol^b) from the SAME validated cwd
    # pools. Check: total volume reconciles with the down-wood biomass; cover follows the power law; round-trip.
    key = joinpath(_CDIR, "carbon_snt.key")
    if !isfile(key)
        @test_skip "carbon_snt scenario not available"
    else
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        FVSjl.compute_forest_type!(s); FVSjl.fmcba!(s)
        dw = FVSjl.ffe_down_wood(s)
        @test dw.vol_hard[8] ≈ sum(dw.vol_hard[1:7]) rtol = 1f-4      # total = sum of bins
        @test dw.cov_hard[7] ≈ sum(dw.cov_hard[1:6]) rtol = 1f-4
        # hard total volume = hard down-wood biomass × 2000 / 24.96 (SG 0.4)
        hard_bio = sum(@view s.fire.cwd[1:9, 2, :])
        @test dw.vol_hard[8] ≈ hard_bio * 2000f0 / 24.96f0 rtol = 1f-3
        # cover bin 1 (3-6", = size class 4 = volume bin 2) follows the power law
        @test dw.cov_hard[1] ≈ 0.0166f0 * dw.vol_hard[2]^0.8715f0 rtol = 1f-4
        rows = [(1990, FVSjl.stand_carbon_report(s), FVSjl.ffe_fuel_loadings(s), FVSjl.snag_summary(s), dw)]
        dbpath = joinpath(mktempdir(), "dwd.db")
        FVSjl.write_dbs_dwd_vol!(dbpath, "C1", "S1", rows)
        FVSjl.write_dbs_dwd_cov!(dbpath, "C1", "S1", rows)
        db = SQLite.DB(dbpath)
        v = [(; T = r.DWD_Volume_Total_Hard) for r in DBInterface.execute(db, "SELECT * FROM FVS_Down_Wood_Vol")]
        c = [(; T = r.DWD_Cover_Total_Hard) for r in DBInterface.execute(db, "SELECT * FROM FVS_Down_Wood_Cov")]
        SQLite.close(db)
        @test length(v) == 1 && v[1].T ≈ Float64(dw.vol_hard[8])
        @test length(c) == 1 && c[1].T ≈ Float64(dw.cov_hard[7])
    end
end

@testset "potential_fire — dual-scenario surface fire behavior (FVS_PotFire core, fmpofl.f)" begin
    # potential_fire computes the SEVERE (fmois1/20mph/70F) + MODERATE (fmois3/8mph/60F) surface fire
    # behavior WITHOUT applying mortality — the value-grounded core of FVS_PotFire (SN skips crown fire).
    s = FVSjl.StandState(FVSjl.Southern())
    FVSjl.init_blockdata!(s, s.variant); FVSjl.init_merch_standards!(s)
    s.plot.forest_type = Int32(520); s.plot.latitude = 35f0; s.plot.longitude = -80f0; s.plot.elevation = 10f0
    t = s.trees; t.n = 2
    t.species[1] = Int32(65); t.dbh[1] = 14f0; t.height[1] = 72f0; t.tpa[1] = 30f0; t.crown_pct[1] = Int32(40)
    t.species[2] = Int32(22); t.dbh[2] = 4f0;  t.height[2] = 18f0; t.tpa[2] = 30f0; t.crown_pct[2] = Int32(50)
    s.fire = FVSjl.FireState(); s.fire.active = true
    pf = FVSjl.potential_fire(s)
    @test pf !== nothing
    @test pf.severe.flame >= pf.moderate.flame                # severe weather → larger fire
    @test pf.severe.scorch >= pf.moderate.scorch
    @test pf.severe.ba_kill >= pf.moderate.ba_kill >= 0f0     # more mortality under severe
    @test pf.severe.smoke > 0f0 && pf.moderate.smoke > 0f0
    @test !isempty(pf.severe.models)
    # non-mutating: the stand TPA is unchanged after a potential-fire computation
    @test t.tpa[1] == 30f0 && t.tpa[2] == 30f0
end

@testset "FVS_Hrv_Carbon — harvested-wood-products carbon fate (fmscut.f/fmchrvout.f, FAPROP)" begin
    # accrue_harvest_carbon! buckets cut merch biomass by product (DBH vs CDBRK) × group (biogrp>5);
    # harvested_carbon_report distributes it by the FAPROP year-since-harvest decay curves into
    # Products/Landfill/Energy/Emissions. The live oracle is binary-blocked, so validate the model semantics.
    s = FVSjl.StandState(FVSjl.Southern())
    FVSjl.init_blockdata!(s, s.variant); FVSjl.init_merch_standards!(s)
    s.fire = FVSjl.FireState(); s.fire.active = true
    FVSjl.accrue_harvest_carbon!(s, 11, 14f0, 2.0f0, 2000)   # pine 14" > 9 → sawtimber softwood
    FVSjl.accrue_harvest_carbon!(s, 65, 6f0, 1.0f0, 2000)    # oak 6" < 11 → pulpwood hardwood
    @test haskey(s.fire.hwp_fate, (2000, 2, 1))              # saw-softwood bucket
    @test haskey(s.fire.hwp_fate, (2000, 1, 2))              # pulp-hardwood bucket
    r0 = FVSjl.harvested_carbon_report(s, 2000, 1)
    r30 = FVSjl.harvested_carbon_report(s, 2030, 1)
    @test r0.products > 0f0
    @test r0.stored ≈ r0.products + r0.landfill              # stored = in-use + landfill
    @test r0.removed ≈ r0.energy + r0.emissions + r0.stored  # removed = energy+emissions+stored
    @test r30.removed ≈ r0.removed rtol = 1f-4               # total removed is fixed at harvest
    @test r30.products < r0.products                         # in-use wood decays over time
    @test r30.landfill >= r0.landfill                        # landfill accumulates
    # DBS round-trip
    dbpath = joinpath(mktempdir(), "hrv.db")
    FVSjl.write_dbs_hrvcarbon!(dbpath, "C1", "S1", [(2000, r0), (2030, r30)])
    db = SQLite.DB(dbpath)
    res = [(; Year = x.Year, P = x.Products, St = x.Merch_Carbon_Stored, Rm = x.Merch_Carbon_Removed)
           for x in DBInterface.execute(db, "SELECT * FROM FVS_Hrv_Carbon ORDER BY Year")]
    SQLite.close(db)
    @test length(res) == 2 && res[1].Year == 2000
    @test res[1].P ≈ Float64(r0.products)
    @test res[1].Rm ≈ Float64(r0.removed)
end

@testset "FVS_PotFire DBS table — potential fire report + writer (fmpofl.f / dbsfmpf.f)" begin
    # potential_fire_report bundles the dual-scenario surface fire + canopy bulk density + torching
    # probability; write_dbs_potfire! writes the 27-col FVS_PotFire. SN: Tot_Flame = Surf_Flame, indices −1.
    s = FVSjl.StandState(FVSjl.Southern())
    FVSjl.init_blockdata!(s, s.variant); FVSjl.init_merch_standards!(s)
    s.plot.forest_type = Int32(520); s.plot.latitude = 35f0; s.plot.longitude = -80f0; s.plot.elevation = 10f0
    t = s.trees; t.n = 3
    t.species[1] = Int32(65); t.dbh[1] = 14f0; t.height[1] = 72f0; t.tpa[1] = 80f0;  t.crown_pct[1] = Int32(30)
    t.species[2] = Int32(22); t.dbh[2] = 4f0;  t.height[2] = 18f0; t.tpa[2] = 200f0; t.crown_pct[2] = Int32(60)
    t.species[3] = Int32(65); t.dbh[3] = 8f0;  t.height[3] = 40f0; t.tpa[3] = 120f0; t.crown_pct[3] = Int32(45)
    s.fire = FVSjl.FireState(); s.fire.active = true
    FVSjl.fmcba!(s)
    r = FVSjl.potential_fire_report(s)
    @test r !== nothing
    @test r.tot_flame_sev == r.surf_flame_sev                 # SN: no crown fire → total = surface
    @test r.torch_index == -1f0 && r.crown_index == -1f0      # FMCFIR skipped in SN
    @test 0f0 <= r.canopy_density <= 0.35f0                   # CBD capped at 0.35
    @test r.surf_flame_sev >= r.surf_flame_mod                # severe weather → larger fire
    @test 0f0 <= r.ptorch_sev <= 1f0 && 0f0 <= r.ptorch_mod <= 1f0
    dbpath = joinpath(mktempdir(), "pf.db")
    FVSjl.write_dbs_potfire!(dbpath, "C1", "S1", [(2003, r)])
    db = SQLite.DB(dbpath)
    res = [(; Year = x.Year, SF = x.Surf_Flame_Sev, CH = x.Canopy_Ht, CD = x.Canopy_Density,
            TI = x.Torch_Index, FM1 = x.Fuel_Mod1)
           for x in DBInterface.execute(db, "SELECT * FROM FVS_PotFire")]
    SQLite.close(db)
    @test length(res) == 1 && res[1].Year == 2003
    @test res[1].SF ≈ Float64(r.surf_flame_sev)
    @test res[1].CD ≈ Float64(r.canopy_density)
    @test res[1].TI == -1.0
end

@testset "FVS_BurnReport DBS table — captured from a SIMFIRE event (dbsfmburn.f)" begin
    # fmburn! captures a burn-event record (moistures, wind, flame, scorch, weighted fuel models) into
    # fire.burn_reports; write_dbs_burnreport! writes the FVS_BurnReport row. Validate capture + round-trip.
    s = FVSjl.StandState(FVSjl.Southern())
    FVSjl.init_blockdata!(s, s.variant); FVSjl.init_merch_standards!(s)
    s.plot.forest_type = Int32(520); s.plot.latitude = 35f0; s.plot.longitude = -80f0; s.plot.elevation = 10f0
    t = s.trees; t.n = 2
    t.species[1] = Int32(65); t.dbh[1] = 14f0; t.height[1] = 72f0; t.tpa[1] = 30f0; t.crown_pct[1] = Int32(40)
    t.species[2] = Int32(22); t.dbh[2] = 4f0;  t.height[2] = 18f0; t.tpa[2] = 30f0; t.crown_pct[2] = Int32(50)
    s.fire = FVSjl.FireState(); s.fire.active = true
    res = FVSjl.fmburn!(s; wind = 20f0, fmois = 1, year = 2003)
    @test length(s.fire.burn_reports) == 1                    # one fire event captured
    b = s.fire.burn_reports[1]
    @test b.year == 2003
    @test b.flame ≈ res.flame && b.scorch ≈ res.scorch         # matches the FireResult
    @test b.killed ≈ res.killed
    @test !isempty(b.models)                                   # at least one weighted fuel model
    dbpath = joinpath(mktempdir(), "burn.db")
    FVSjl.write_dbs_burnreport!(dbpath, "C1", "S1", s.fire.burn_reports)
    db = SQLite.DB(dbpath)
    rows = [(; Year = r.Year, Flame = r.Flame_length, Scorch = r.Scorch_height,
             M1 = r.One_Hr_Moisture, FM1 = r.FuelModl1, Slope = r.Slope)
            for r in DBInterface.execute(db, "SELECT * FROM FVS_BurnReport")]
    SQLite.close(db)
    @test length(rows) == 1 && rows[1].Year == 2003
    @test rows[1].Flame ≈ Float64(b.flame)
    @test rows[1].M1 ≈ Float64(b.mois[1,1]) * 100              # moisture reported as percent
    @test rows[1].Slope == 0                                   # SN surface-fire path: no slope term

    # FVS_Mortality: killed vs total TPA by DBH class (Total = killed + remaining, pre-fire)
    @test b.killed_ba > 0f0
    @test sum(b.clskil) ≈ b.killed rtol = 1f-4
    @test all(b.totcls[c] >= b.clskil[c] for c in 1:7)        # total ≥ killed in each class
    FVSjl.write_dbs_mortality!(dbpath, "C1", "S1", s.fire.burn_reports)
    FVSjl.write_dbs_consumption!(dbpath, "C1", "S1", s.fire.burn_reports)
    db2 = SQLite.DB(dbpath)
    mort = [(; K3 = r.Killed_class3, T3 = r.Total_class3, BA = r.Bakill)
            for r in DBInterface.execute(db2, "SELECT * FROM FVS_Mortality")]
    cons = [(; ST = r.Surface_Total) for r in DBInterface.execute(db2, "SELECT * FROM FVS_Consumption")]
    SQLite.close(db2)
    @test length(mort) == 1 && mort[1].BA ≈ Float64(b.killed_ba)
    @test mort[1].K3 ≈ Float64(b.clskil[3])                   # 14" tree → class 3 (10-20")
    @test length(cons) == 1 && cons[1].ST ≈ Float64(b.consumed.surf_total)
    @test b.consumed.surf_total >= 0f0                        # fire consumes (≥0) surface fuel
end

@testset "Input-snag seeding — inventory Stand-Dead from input dead records (FMSADD ITYP=3)" begin
    # ffe_seed_input_snags! seeds FFE snags from the input dead-tree records (carbon_snt has sp65 d34.6
    # hist=8 and sp27 d7.2 hist=6). The snag STEM-volume bole (local height-dub + R8 Clark volume × V2T,
    # since the dead partition has no input height/volume) reproduces the Fortran inventory Stand-Dead
    # (3.8 mt/ha). Validated as a standalone mechanism here; wiring it into the grown-cycle report needs
    # the age-dependent snag falldown + the crown-lift (else grown-cycle DDW/Stand-Dead overshoot) —
    # see FFE_FUEL_DYNAMICS_chunk_plan.md.
    key = joinpath(_CDIR, "carbon_snt.key")
    if !isfile(key)
        @test_skip "carbon_snt scenario not available"
    else
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        @test s.trees.ndead == 2                        # the two input dead records were read
        s.fire !== nothing && s.fire.active && FVSjl.compute_forest_type!(s)
        s.fire !== nothing && s.fire.active && FVSjl.fmcba!(s)
        FVSjl.ffe_seed_input_snags!(s)
        TO = 0.90718474 / 0.40468564
        # input snags carry no crown (history ≥7 ⇒ crown fallen), so Stand-Dead = the stem bole
        @test FVSjl.snag_crown_carbon(s) == 0f0
        # input-snag bole = 3.77 vs the Fortran inventory Stand-Dead 3.8 (within print resolution): FVS's
        # snag bole is the MERCHANTABLE cubic (FMDOUT→FMSVOL→CFVOL, v[4]), not the gross total stem (v[1]);
        # the <top-dia tip is a large fraction for small snags (sp27 d7.2: v[1]=5.2 vs v[4]=4.8 = FVS).
        @test abs(FVSjl.snag_bole_carbon(s) * TO - 3.8) <= 0.05
        @test FVSjl.snag_standing_density(s.fire) > 0f0            # snags actually created
    end
end

@testset "Snag falldown — age-aware Stand-Dead tracks the Fortran (carbon_snt, deaths-spreading fix)" begin
    # update_snags! ages each snag by its OWN (current year − death year), not a blanket nyears (FMSNAG):
    # a cycle's fresh mortality is dated ~at the boundary, so it must not fall a full cycle. With the
    # blanket-nyears falldown the Stand-Dead COLLAPSED (3.9→1.9); age-aware it INCREASES and tracks the
    # Fortran (3.8→4.4→5.4→9.5). This drives the snag-bole half of Stand-Dead; the small residual is the
    # un-flowed mortality crown (cwd2b) + the pre-inventory input-snag age. Demonstrates the dynamics fix.
    key = joinpath(_CDIR, "carbon_snt.key")
    if !isfile(key)
        @test_skip "carbon_snt scenario not available"
    else
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        s.fire !== nothing && s.fire.active && FVSjl.compute_forest_type!(s)
        s.fire !== nothing && s.fire.active && FVSjl.fmcba!(s)
        FVSjl.ffe_seed_input_snags!(s)
        TO = 0.90718474 / 0.40468564
        fF = [3.8, 4.4, 5.4, 9.5]
        prev = 0.0; maxresid = 0.0
        for c in 1:4
            FVSjl.compute_density!(s)
            sd = FVSjl.standing_dead_carbon(s) * TO
            c >= 2 && @test sd > prev        # INCREASES every cycle (was collapsing before the fix) — the
                                             # semantic this test exists for (age-aware snag falldown)
            maxresid = max(maxresid, abs(sd - fF[c]))
            prev = sd
            if c < 4
                FVSjl.ffe_fuel_update!(s, 5)      # cwd2b crown flow + decay (as the report does)
                # grow_cycle! ALREADY ages the snags (update_snags!, simulate.jl:211, run before
                # mortality creates the cycle's new snags) — do NOT call it again or snags fall twice.
                FVSjl.grow_cycle!(s; fint = 5f0)
            end
        end
        # Stand-Dead tracks the Fortran but is NOT bit-exact (the ~0.7 dead-pool / height-dub residual) —
        # @test_broken on the max, honest rather than hidden behind a 0.8 tolerance.
        @test_broken maxresid <= 0.05
    end
end
