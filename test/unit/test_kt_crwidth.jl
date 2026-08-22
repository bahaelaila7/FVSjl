# KT (Kootenai) FVS_TreeList CrWidth via national cwcalc.f KTMAP dispatch (regression, 2026-08-22).
# KTMAP (vie/cwcalc.f) = IEMAP's first 11 codes verbatim (KT MAXSP=11, Region 1, BF=1.0); every code is
# already carried by ie_cwcalc's national-form dispatch, so kt_cwcalc delegates to ie_cwcalc. VALIDATED vs
# FVSkt_clean kttl FVS_TreeList: function-level A/B bit-exact on all 6132 rows across the 6 emitted species
# (WP/LP/ES/AF/DF log-forms + MH/OT R6M2 26405), worst |Δ|=0.0. KT is the 14th variant past the eastern 0.5.

using Test
using FVSjl

@testset "KT FVS_TreeList CrWidth — national cwcalc.f KTMAP → ie_cwcalc, bit-exact vs FVSkt_clean" begin
    cw(sp, d, h, cr, ba) = FVSjl.kt_cwcalc(sp, Float32(d), Float32(h), Float32(cr), Float32(ba), 45f0, 0f0)
    @test isapprox(cw(1, 0.7, 5.0, 90, 3.0),   2.896856f0;  atol = 1f-5)   # WP 11903 (log, BAREA)
    @test isapprox(cw(7, 1.0, 8.0, 85, 30.0),  4.0850024f0; atol = 1f-5)   # LP 10803 (log)
    @test isapprox(cw(11, 2.0, 12.0, 80, 50.0), 5.0578794f0; atol = 1f-5)  # MH/OT 26405 (R6M2, elev)
    # KTMAP[sp] == IEMAP[sp] for sp 1..11 — kt_cwcalc delegates:
    @test FVSjl._KT_CWMAP == FVSjl._IE_CWMAP[1:11]
    @test cw(1, 0.7, 5.0, 90, 3.0) == FVSjl.ie_cwcalc(1, 0.7f0, 5.0f0, 90f0, 3.0f0, 45f0, 0f0)
    @test 0.5f0 <= cw(1, 0.05, 0.5, 40, 2.0) <= 99.9f0
end
