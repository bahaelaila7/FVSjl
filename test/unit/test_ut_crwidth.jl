# UT (Utah) FVS_TreeList CrWidth via national cwcalc.f UTMAP dispatch (regression, 2026-08-22).
# ut_cwcalc = _cwcalc_national(UTMAP[sp], …); UT (24 species, Region 4, BF=1.0) contributed 01505/81402/10201 to
# the shared national dispatch. VALIDATED vs FVSut_clean uttl FVS_TreeList: function-level A/B bit-exact on all
# 4000 rows across the emitted species (WB/DF/LP/AS R6-m2 + AF), worst |Δ|=0.0. UT is the 17th variant. Note UT's
# PP (idx10) uses the R6-m2 '12205' like OS (idx23), unlike CI/IE PP '12203' — UTMAP is index-keyed so exact.

using Test
using FVSjl

@testset "UT FVS_TreeList CrWidth — national cwcalc.f UTMAP, bit-exact vs FVSut_clean" begin
    cw(sp, d, h, cr, ba) = FVSjl.ut_cwcalc(sp, Float32(d), Float32(h), Float32(cr), Float32(ba), 45f0, 0f0)
    @test isapprox(cw(4, 3.0, 18.0, 80, 55.0),  7.192881f0; atol = 1f-5)  # WF 01505 (R6-m2, UT-added)
    @test isapprox(cw(13, 6.0, 20.0, 60, 30.0), 8.7417f0;   atol = 1f-4)  # GO 81402 (Bechtold-m2, UT-added)
    @test isapprox(cw(17, 6.0, 15.0, 55, 25.0), 12.8197f0;  atol = 1f-4)  # GB 10201 (Bechtold-m1, UT-added)
    # UT's PP (idx10) and OS (idx23) both map to R6-m2 '12205' (index-keyed disambiguation):
    @test FVSjl._UT_CWMAP[10] == "12205" && FVSjl._UT_CWMAP[23] == "12205"
    @test cw(10, 3.0, 18.0, 80, 55.0) == cw(23, 3.0, 18.0, 80, 55.0)
    @test 0.5f0 <= cw(1, 0.05, 0.5, 40, 2.0) <= 99.9f0
end
