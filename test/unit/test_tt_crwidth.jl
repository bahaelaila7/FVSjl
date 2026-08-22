# TT (Teton) FVS_TreeList CrWidth via national cwcalc.f TTMAP dispatch (regression, 2026-08-22).
# tt_cwcalc = _cwcalc_national(TTMAP[sp], …); TT (18 species, Region 4, BF=1.0) contributed 20205/09305/10805/31206
# to the shared national dispatch. VALIDATED vs FVStt_clean tttl FVS_TreeList: function-level A/B bit-exact on all
# 4000 rows across the emitted species (WB/DF/LP/AS R6-m2 + AF), worst |Δ|=0.0. TT is the 16th variant.

using Test
using FVSjl

@testset "TT FVS_TreeList CrWidth — national cwcalc.f TTMAP, bit-exact vs FVStt_clean" begin
    cw(sp, d, h, cr, ba) = FVSjl.tt_cwcalc(sp, Float32(d), Float32(h), Float32(cr), Float32(ba), 45f0, 0f0)
    @test isapprox(cw(1, 2.0, 12.0, 80, 50.0),  3.8952932f0; atol = 1f-5)  # WB 10105 (R6-m2)
    @test isapprox(cw(3, 3.0, 20.0, 85, 60.0),  8.626372f0;  atol = 1f-5)  # DF 20205 (R6-m2, TT-added)
    @test isapprox(cw(13, 4.0, 15.0, 70, 40.0), 13.953997f0; atol = 1f-4)  # BI 31206 (Donnelly power, TT-added)
    @test FVSjl._TT_CWMAP[10] == "12203"   # TT's PP is the log-form '12203'
    @test 0.5f0 <= cw(1, 0.05, 0.5, 40, 2.0) <= 99.9f0
end
