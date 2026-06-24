# test_crown_lift.jl — FFE crown-LIFT rate X (FMSDIT, fmsdit.f:103-117).
#
# As a tree's crown base rises over a cycle, the lower crown dies into down wood at the annual fraction
# X = (NEWBOT − OLDBOT)/OLDCRL/CYCLEN. This is the dominant post-mortality down-wood (DDW) addition —
# its absence was the last carbon-report residual (FVSjl 2.1 vs Fortran 3.8 @2000). This unit test pins
# the X EQUATION + its edge cases against the Fortran, so the validated formula is locked ahead of the
# tree-record plumbing that feeds it the previous-cycle crown (see FFE_FUEL_DYNAMICS_chunk_plan.md).

using Test, FVSjl
using FVSjl: crown_lift_rate

@testset "FFE crown-lift rate X (FMSDIT)" begin
    # Worked example: crown base rises 32 → 35 ft over a 5-yr cycle, old crown length 8 ft.
    #   OLDBOT = oldht − oldcrl = 40 − 8 = 32 ; NEWBOT = ht − ht·cr% = 50 − 50·0.30 = 35
    #   rise = 3 ; X = 3 / 8 / 5 = 0.075 /yr
    @test crown_lift_rate(40f0, 8f0, 50f0, 30f0, 5) ≈ 0.075f0

    # No base rise (crown grew downward / static) ⇒ X = 0 (fmsdit.f:106 guard).
    # oldht=40 oldcrl=8 ⇒ OLDBOT=32 ; ht=44 cr%=30 ⇒ NEWBOT = 44−13.2 = 30.8 < 32 ⇒ no lift.
    @test crown_lift_rate(40f0, 8f0, 44f0, 30f0, 5) == 0f0

    # Degenerate old crown (OLDCRL ≤ 0.001) ⇒ X = 0, no divide-by-zero (fmsdit.f:106).
    @test crown_lift_rate(40f0, 0f0, 50f0, 30f0, 5) == 0f0
    @test crown_lift_rate(40f0, 0.0005f0, 50f0, 30f0, 5) == 0f0

    # Inversely proportional to cycle length (same geometry, longer cycle ⇒ smaller annual fraction).
    @test crown_lift_rate(40f0, 8f0, 50f0, 30f0, 10) ≈ 0.0375f0           # double cyclen ⇒ half X
    # NB OLDCRL enters BOTH OLDBOT and the denominator, so it is not a clean inverse knob: with
    # oldcrl=4 the old base 36 > the new base 35 ⇒ no rise ⇒ X=0 (the fmsdit.f:106 guard).
    @test crown_lift_rate(40f0, 4f0, 50f0, 30f0, 5) == 0f0
    # A larger base rise (taller tree, higher crown-ratio base) scales X up proportionally.
    @test crown_lift_rate(40f0, 8f0, 60f0, 40f0, 5)  ≈ ((60f0-60f0*0.4f0)-32f0)/8f0/5f0
end
