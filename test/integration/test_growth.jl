# test_growth.jl — GROWTH keyword (opt 13, vbase/initre.f:2300).
#
# GROWTH sets the INPUT growth-data type codes (IDG/IHTG) + measurement periods (FINT/FINTH/FINTM)
# the LSTART calibration uses. The defaults (IDG/IHTG=0, periods=5) are FVSjl's current bit-exact
# behaviour — the DG field is the increment over 5 yr. This test covers RECOGNITION + parameter
# CAPTURE, plus the IDG/IHTG=1/3 PAST-DBH/HT interpretation (intree.f:536: the DG/HTG field is the
# past DBH/HT `PDBH`/`PHT`, so the increment = current − past). `apply_growth_input_types!` converts
# it to the increment, after which the already-bit-exact IDG=0 calibration runs unchanged.
# VALIDATED vs live Fortran on a purpose-built wide-DG-field stand: Fortran IDG=1 (past-DBH field) ⇒
# byte-identical `.sum` to IDG=0 (increment field); FVSjl reproduces that exact equivalence (below).
# (The non-default-FINT period scaling remains the deferred WK3 past-DBH calibration chunk.)

using Test, FVSjl

const _GDIR = joinpath(@__DIR__, "..", "harness", "scenarios")

# Split a `.sum` text into per-cycle data rows (the lines beginning with a 4-digit year),
# each as a vector of whitespace-separated column strings: [year, age, TPA, BA, SDI, CCF,
# TopHt, QMD, TCuFt, MCuFt, SCuFt, BdFt, …].
_g_rows(sumtext::AbstractString) =
    [split(strip(l)) for l in split(sumtext, "\n") if occursin(r"^(19|20)\d\d\s", strip(l))]

@testset "GROWTH keyword — recognition + parameter capture" begin
    mkrec(vals, present) = FVSjl.KeywordRecord("GROWTH", "",
        [v == 0 ? "" : string(v) for v in vals], Float32.(vals), present, 12, FVSjl.KW_OK, 0)

    # defaults: a bare GROWTH leaves IDG/IHTG=0, periods=5
    s0 = FVSjl.StandState(FVSjl.Southern())
    @test s0.control.growth_idg == 0 && s0.control.growth_fint == 5f0

    # explicit: IDG=1, FINT=10, IHTG=3, FINTH=8, FINTM=7
    s = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_growth!(s, mkrec([1, 10, 3, 8, 7, zeros(Int, 7)...],
                              [true, true, true, true, true, falses(7)...]))
    @test s.control.growth_idg == 1 && s.control.growth_fint == 10f0
    @test s.control.growth_ihtg == 3 && s.control.growth_finth == 8f0 && s.control.growth_fintm == 7f0

    # a blank measurement-period field keeps the default (Fortran: only ARRAY(n)>0 overrides)
    s2 = FVSjl.StandState(FVSjl.Southern())
    FVSjl.kw_growth!(s2, mkrec([2, 0, 0, 0, 0, zeros(Int, 7)...],
                               [true, false, false, false, false, falses(7)...]))
    @test s2.control.growth_idg == 2 && s2.control.growth_fint == 5f0   # period unchanged

    # snt01 (no GROWTH) stays bit-exact — the default path is the current behaviour. A default-valued
    # GROWTH is a no-op in live Fortran too (verified: snt01 + bare GROWTH ⇒ byte-identical .sum).
    snt = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"
    if isfile(snt)
        base = FVSjl.run_keyfile(snt; faithful = true)
        @test occursin("1990", base)        # runs cleanly; GROWTH default doesn't perturb it
    end
end

@testset "GROWTH IDG=1 — past-DBH field + BRATIO bark correction vs live Fortran" begin
    # IDG=1: the DG field is the PAST DBH. sn/cratet.f converts it to the OUTSIDE-bark increment
    # DBH−past (fed to DENSE backdating), then dgdriv.f:333 multiplies by BRATIO → the inside-bark
    # increment used by the calibration. ⚠ The BRATIO step is NOT a no-op: IDG=1 reading a past DBH
    # gives a MATERIALLY different projection than IDG=0 reading the raw DBH−past as the increment
    # (Fortran: at 2000 BA 167 vs 174). The committed baseline `growth_idg1.sum.save` is live Fortran
    # `GROWTH IDG=1` on a wide-DG-field stand (≥5 LP trees so the calibration actually fires — a
    # 3-tree stand would skip calibration and make the DG field irrelevant, the trap of the earlier
    # vacuous test). FVSjl must match its structural columns; the ~1.3% volume gap is the known LP
    # growth-calibration tail, present identically under IDG=0 (orthogonal to the IDG path).
    sc1 = joinpath(_GDIR, "growth_idg1.key"); sav = joinpath(_GDIR, "growth_idg1.sum.save")
    if !isfile(sc1) || !isfile(sav)
        @test_skip "growth_idg1 scenario not available"
    else
        jl = _g_rows(FVSjl.run_keyfile(sc1; faithful = true))
        ft = _g_rows(read(sav, String))
        @test length(jl) == length(ft) >= 3
        for (j, f) in zip(jl, ft)
            # BIT-EXACT vs live (re-grounded 2026-07-02): all 12 .sum columns match the live golden exactly.
            # The "LP-calibration tail ~1.3%" tolerance was STALE — that residual appears on LP-heavy stands
            # (nohtdreg_cal), NOT this scenario. All exact `==`, no atol/relative slack.
            for c in 1:12; @test j[c] == f[c]; end   # YEAR/TPA/BA/SDI/CCF/TopHt/QMD/Tcf/Mcf/Scf/Bdft
        end
        # BRATIO is ACTIVE: IDG=1 (past DBH) ≠ IDG=0 reading the same raw DBH−past as the increment.
        sc0 = joinpath(_GDIR, "growth_idg0.key")
        if isfile(sc0)
            j0 = _g_rows(FVSjl.run_keyfile(sc0; faithful = true))
            i2000 = findfirst(r -> r[1] == "2000", jl); k2000 = findfirst(r -> r[1] == "2000", j0)
            @test i2000 !== nothing && k2000 !== nothing
            @test parse(Int, jl[i2000][4]) != parse(Int, j0[k2000][4])     # BA differs (167 vs 174)
        end
    end
end

@testset "GROWTH FINT — input measurement-period scaling vs live Fortran" begin
    # FINT = the GROWTH keyword's DIAMETER measurement period. The input DG increment spans FINT
    # years, so the calibration rescales the measured DDS to the 5-yr model basis: SCALE = YR/FINT =
    # 5/FINT (dgdriv.f:325), wired in setup_growth!. A FINT=10 stand reads the SAME 1.5" increment as
    # covering 10 yr (half the annual growth) ⇒ the calibration shrinks future growth — a MATERIAL,
    # non-default effect (Fortran 1995 BA 129 vs the FINT=5 158). Baseline growth_fint10.sum.save is
    # live Fortran on the same ≥5-LP-tree calibration-firing stand. FVSjl must match its structural
    # columns; the ~1.5% cuft gap is the LP calibration tail (identical at FINT=5, orthogonal to FINT).
    sc = joinpath(_GDIR, "growth_fint10.key"); sav = joinpath(_GDIR, "growth_fint10.sum.save")
    if !isfile(sc) || !isfile(sav)
        @test_skip "growth_fint10 scenario not available"
    else
        jl = _g_rows(FVSjl.run_keyfile(sc; faithful = true))
        ft = _g_rows(read(sav, String))
        @test length(jl) == length(ft) >= 3
        for (j, f) in zip(jl, ft)
            # BIT-EXACT vs live (re-grounded 2026-07-02): all 12 .sum columns match exactly. The FINT=10
            # "LP-calibration tail ~1.5%" tolerance was STALE — this scenario matches live bit-exactly.
            for c in 1:12; @test j[c] == f[c]; end   # YEAR/TPA/BA/SDI/CCF/TopHt/QMD/Tcf/Mcf/Scf/Bdft
        end
        # The scaling is ACTIVE: FINT=10 grows materially slower than FINT=5 on the SAME tree data
        # (growth_idg0 = the identical .tre with the default FINT=5). BA at 1995: ≈128 vs ≈158.
        sc5 = joinpath(_GDIR, "growth_idg0.key")
        if isfile(sc5)
            j5 = _g_rows(FVSjl.run_keyfile(sc5; faithful = true))
            i = findfirst(r -> r[1] == "1995", jl); k = findfirst(r -> r[1] == "1995", j5)
            @test i !== nothing && k !== nothing
            @test parse(Int, jl[i][4]) < parse(Int, j5[k][4]) - 20         # FINT=10 BA well below FINT=5
        end
    end
end

@testset "GROWTH FINTH — height-measurement-period scaling vs live Fortran" begin
    # FINTH scales the SMALL-tree height-growth (HCOR) calibration: the measured HTG spans FINTH
    # years, so the regression rescales it to the 5-yr basis via SCALE3 = REGYR/FINTH = 5/FINTH
    # (regent.f:406/462, `TERM = HTG·SCALE3`). The firing stand needs ≥NCALHT=5 SMALL (dbh<5) trees
    # of a species with measured HTG. Baseline growth_finth10.sum.save is live Fortran at FINTH=10;
    # the residual is the small-tree REGENT growth tail (present identically at FINTH=5).
    sc = joinpath(_GDIR, "growth_finth10.key"); sav = joinpath(_GDIR, "growth_finth10.sum.save")
    sc5 = joinpath(_GDIR, "growth_finth5.key")
    if !isfile(sc) || !isfile(sav)
        @test_skip "growth_finth10 scenario not available"
    else
        jl = _g_rows(FVSjl.run_keyfile(sc; faithful = true))
        ft = _g_rows(read(sav, String))
        @test length(jl) == length(ft) >= 3
        for (j, f) in zip(jl, ft)
            # BIT-EXACT vs live (re-grounded 2026-07-02): all 12 .sum columns match exactly. The FINTH=10
            # "small-tree REGENT tail" TPA tolerance (1.5%+3) was STALE — this scenario matches live bit-exactly.
            for c in 1:12; @test j[c] == f[c]; end   # YEAR/TPA/BA/SDI/CCF/TopHt/QMD/Tcf/Mcf/Scf/Bdft
        end
        # SCALE3 is ACTIVE and matches Fortran: the FINTH5→FINTH10 delta is reproduced. At 1995 both
        # engines move BA 325→330 and TPA 4392/4419→4123/4150 (ΔTPA −269 in each).
        if isfile(sc5)
            j5 = _g_rows(FVSjl.run_keyfile(sc5; faithful = true))
            i = findfirst(r -> r[1] == "1995", jl); k = findfirst(r -> r[1] == "1995", j5)
            @test i !== nothing && k !== nothing
            d_jl = parse(Int, j5[k][3]) - parse(Int, jl[i][3])             # FVSjl ΔTPA (5→10)
            @test 240 <= d_jl <= 300                                       # ≈269, matching Fortran
            @test parse(Int, jl[i][4]) > parse(Int, j5[k][4])             # FINTH=10 BA above FINTH=5
        end
    end
end

@testset "NOTRE FINT/FINTM dead-record inflation (notre.f:122-124) vs Fortran" begin
    # When the growth (FINT) and mortality-observation (FINTM) periods differ, NOTRE inflates DEAD-record PROB
    # by FINT/FINTM so the recent dead are added back at the right rate to recover the BACKDATED density the DG
    # calibration reads. jl now applies that inflation, scoped to the calibration density + PCTILE passes
    # (diameter_growth.jl `calibrate_diameter_growth!`). dead_fint = carbon_snt (10 input dead trees, history
    # 6-9) + `GROWTH 0 10 0 0 5` (FINT=10, FINTM=5 ⇒ 2× inflation). Without the fix the missing inflation
    # under-counts backdated density → the COR drifts (2000 was TCuFt −8 / MCuFt −9 vs live); with it the .sum
    # tracks live to ≤4 cuft. The small residual = the dead trees' contribution to the crown-init backdated CCF,
    # which jl computes live-only. (The default FINT=FINTM case stays bit-exact — the fix is gated on FINT≠FINTM.)
    key = joinpath(_GDIR, "dead_fint.key"); sav = joinpath(_GDIR, "dead_fint.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "dead_fint scenario not available"
    else
        jl = _g_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _g_rows(read(sav, String))
        @test length(jl) == length(ft) >= 4
        for (j, f) in zip(jl, ft)
            @test j[1] == f[1]                                              # YEAR
            @test abs(parse(Int, j[4]) - parse(Int, f[4])) <= 1            # BA
            @test abs(parse(Int, j[9]) - parse(Int, f[9])) <= 4           # TCuFt (pre-fix gap was 8)
            @test abs(parse(Int, j[10]) - parse(Int, f[10])) <= 4         # MCuFt (pre-fix gap was 9)
        end
    end
end
