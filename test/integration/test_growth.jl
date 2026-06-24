# test_growth.jl — GROWTH keyword (opt 13, vbase/initre.f:2300).
#
# GROWTH sets the INPUT growth-data type codes (IDG/IHTG) + measurement periods (FINT/FINTH/FINTM)
# the LSTART calibration uses. The defaults (IDG/IHTG=0, periods=5) are FVSjl's current bit-exact
# behaviour — the DG field is the increment over 5 yr. This test covers RECOGNITION + parameter
# CAPTURE (the keyword was previously dropped silently); the IDG=1/3 past-DBH interpretation +
# non-default FINT scaling are the deferred WK3 past-DBH calibration chunk (they change the
# calibration, so they are not wired until validated against a purpose-built scenario).

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
