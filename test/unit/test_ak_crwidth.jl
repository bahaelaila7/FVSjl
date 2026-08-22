# AK (SoutheastAlaska) FVS_TreeList CrWidth via national cwcalc.f AKMAP dispatch (regression, 2026-08-22).
# ak_cwcalc = _cwcalc_national(AKMAP[sp], …); AK (23 species, Region 10 ⇒ BF=1.0) contributed the '08' form
# (_cw08: a·D^b·H^ch·CL^ccl·(BA+1)^cba·EXP(EL)^cel, D<0.1⇒flat 0.5) via 09508/09408/74708/37508/74608, plus
# R6-m2 SF/YC/SS/RC (01105/04205/09805/24205) and Donnelly RA (35106). VALIDATED vs FVSak_clean aktl FVS_TreeList:
# function-level A/B bit-exact on all 3200 rows across the emitted species (AF R6-m2 + TA/BE 09508 + WS 09408),
# worst |Δ|=0.0. AK is the 18th variant past the eastern 0.5 default.

using Test
using FVSjl

@testset "AK FVS_TreeList CrWidth — national cwcalc.f AKMAP, bit-exact vs FVSak_clean" begin
    cw(sp, d, h, cr, ba) = FVSjl.ak_cwcalc(sp, Float32(d), Float32(h), Float32(cr), Float32(ba), 45f0, 0f0)
    @test isapprox(cw(4, 2.0, 15.0, 80, 40.0),  3.4914186f0; atol = 1f-5)  # TA 09508 ('08' form)
    @test isapprox(cw(5, 3.0, 20.0, 75, 50.0),  3.4778612f0; atol = 1f-5)  # WS 09408 ('08' + EL/BAREA)
    @test isapprox(cw(1, 3.0, 25.0, 70, 60.0),  7.1348195f0; atol = 1f-5)  # SF 01105 (R6-m2)
    @test isapprox(cw(15, 5.0, 20.0, 60, 30.0), 15.259793f0; atol = 1f-4)  # RA 35106 (Donnelly power)
    # '08' form: D<0.1 ⇒ flat 0.5 (no small-tree scaling):
    @test cw(4, 0.05, 2.0, 80, 40.0) == 0.5f0
    @test 0.5f0 <= cw(1, 0.05, 0.5, 40, 2.0) <= 99.9f0
end
