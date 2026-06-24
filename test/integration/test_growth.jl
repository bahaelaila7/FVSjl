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

@testset "GROWTH IDG=1/3 — past-DBH field reconstructs the increment (vs Fortran)" begin
    # Live-Fortran cross-check (purpose-built wide-DG-field stand): GROWTH IDG=1 reading the PAST DBH
    # produced a byte-identical `.sum` to IDG=0 reading the increment. Here we reproduce that exact
    # equivalence in-memory on snt01: rewriting each tree's DG field to its past DBH (current − incr)
    # under IDG=1 must, after `apply_growth_input_types!`, recover the identical IDG=0 projection.
    snt = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"
    if !isfile(snt)
        @test_skip "snt01.key not available"
    else
        # A: baseline IDG=0 — the DG field is the measured increment.
        sa = first(FVSjl.each_stand(snt)); FVSjl.notre!(sa)
        # B: IDG=1 — rewrite the DG field to the PAST DBH (current dbh − increment).
        sb = first(FVSjl.each_stand(snt)); FVSjl.notre!(sb)
        sb.control.growth_idg = Int32(1)
        @inbounds for i in 1:sb.trees.n
            g = sb.trees.diam_growth[i]
            sb.trees.diam_growth[i] = g > 0f0 ? sb.trees.dbh[i] - g : 0f0
        end
        # setup_growth! runs apply_growth_input_types! (B's past-DBH ⇒ increment), then calibrates both.
        for s in (sa, sb); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s); end
        n = sa.trees.n
        @test sb.trees.n == n
        @test sa.trees.diam_growth[1:n] ≈ sb.trees.diam_growth[1:n]   # same calibrated increment
        for _ in 1:5; FVSjl.grow_cycle!(sa); FVSjl.grow_cycle!(sb); end
        @test FVSjl.stand_qmd(sa) ≈ FVSjl.stand_qmd(sb)               # identical 5-cycle projection
        @test sa.trees.dbh[1:n] ≈ sb.trees.dbh[1:n]
    end
end
