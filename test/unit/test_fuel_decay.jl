# test_fuel_decay.jl — FFE surface-fuel decay (FMCWD, fuel_decay.jl).
#
# The decay half of the FFE fuel dynamics: cwd[size, soft/hard, decay] *= (1−DKR·{1.1 soft})^nyrs
# per size/decay class, a PRDUFF fraction of the decayed woody → duff, and a hard→soft transfer.
# End-to-end grown-cycle carbon validation needs the FMCADD additions too (litter decays away here);
# this unit test pins the decay EQUATION + the DKR/PRDUFF constants against the Fortran (fmcwd.f /
# fmvinit.f), so the faithful transliteration is locked before the additions chunk lands.

using Test, FVSjl
using FVSjl: FireState, fmcwd!, _FM_DKR, _FM_PRDUFF, StandState, Southern

@testset "FMCWD surface-fuel decay" begin
    # DKR constants (sn/fmvinit.f:70-104): woody class-1 rate 0.11; class-2+ 0.07 for big wood,
    # 0.09 for 1-3"; litter 0.65; duff 0.002. PRDUFF 0.02.
    @test _FM_DKR[1, 1] == 0.11f0 && _FM_DKR[3, 2] == 0.09f0 && _FM_DKR[4, 2] == 0.07f0
    @test _FM_DKR[10, 1] == 0.65f0 && _FM_DKR[11, 1] == 0.002f0
    @test _FM_PRDUFF == 0.02f0

    s = StandState(Southern())
    s.fire = FireState(); s.fire.active = true
    cwd = s.fire.cwd; fill!(cwd, 0f0)
    # seed a single decay class L=1: a hard woody 1-3" pool, hard litter, hard duff
    cwd[3, 2, 1] = 10f0      # woody 1-3" (DKR 0.11 in class 1)
    cwd[10, 2, 1] = 5f0      # litter (DKR 0.65)
    cwd[11, 2, 1] = 6f0      # duff   (DKR 0.002)

    fmcwd!(s, 5)             # five years

    # duff persists (0.002/yr) and gains the PRDUFF transfer from EVERY decayed class 1-10 (here the
    # woody 1-3" pool AND the litter pool both shed 0.02 of what they decayed into duff)
    woody_decayed  = 10f0 - 10f0 * (1f0 - 0.11f0)^5
    litter_decayed =  5f0 -  5f0 * (1f0 - 0.65f0)^5
    duff_expected = 6f0 * (1f0 - 0.002f0)^5 + (woody_decayed + litter_decayed) * _FM_PRDUFF
    @test cwd[11, 2, 1] == duff_expected              # BIT-EXACT: jl computes (1−DKR)^nyrs as a single power, == the closed form (was rtol 1e-4 padding)
    @test cwd[11, 2, 1] > 5.9f0                       # essentially unchanged + small gain

    # litter crashes at 0.65/yr ⇒ ~5·0.35^5 ≈ 0.026 (needs litterfall to hold up — the coupling)
    @test cwd[10, 2, 1] < 0.05f0

    # the woody hard pool decays AND sheds to soft (J<10 hard→soft transfer), conserving woody mass
    # minus what decayed/duffed: hard+soft ≈ 10·(1−0.11)^5 (the decay), split across soft/hard
    woody_after = cwd[3, 1, 1] + cwd[3, 2, 1]
    @test woody_after == 10f0 * (1f0 - 0.11f0)^5      # BIT-EXACT: hard+soft conserves to the closed-form decay (was rtol 1e-3 padding)
    @test cwd[3, 1, 1] > 0f0                           # some moved hard → soft

    # inactive FFE ⇒ no-op
    s2 = StandState(Southern())
    @test fmcwd!(s2, 5) === s2                         # s2.fire === nothing
end
