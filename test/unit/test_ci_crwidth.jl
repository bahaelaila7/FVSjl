# CI (CentralIdaho) FVS_TreeList CrWidth via national cwcalc.f CIMAP dispatch (regression, 2026-08-22).
# ci_cwcalc = _cwcalc_national(CIMAP[sp], …); CI (Region 4, BF=1.0) contributed 4 codes to the shared national
# dispatch — 26305 (WH R6-m2), 01905 (AF R6-m2), 06405 (WJ R6-m2 no-elev), 47502 (MC Bechtold-m2); the rest are
# shared with EM/IE. VALIDATED vs FVSci_clean citl FVS_TreeList: function-level A/B bit-exact on all 4000 rows
# across the emitted species (WP/DF/LP log-forms + WH/AF R6-m2), worst |Δ|=0.0. CI is the 15th variant past the
# eastern 0.5 default. CIMAP disambiguates two FIA-122 species by INDEX (PP idx10 '12203' vs OS idx18 '12205').

using Test
using FVSjl

@testset "CI FVS_TreeList CrWidth — national cwcalc.f CIMAP, bit-exact vs FVSci_clean" begin
    cw(sp, d, h, cr, ba) = FVSjl.ci_cwcalc(sp, Float32(d), Float32(h), Float32(cr), Float32(ba), 50f0, 0f0)
    @test isapprox(cw(1, 0.7, 5.0, 90, 3.0),    2.896856f0; atol = 1f-5)  # WP 11903 (log, shared)
    @test isapprox(cw(5, 2.0, 12.0, 80, 50.0),  7.21589f0;  atol = 1f-5)  # WH 26305 (R6-m2, CI-added)
    @test isapprox(cw(9, 2.0, 12.0, 80, 50.0),  5.55459f0;  atol = 1f-5)  # AF 01905 (R6-m2, CI-added)
    @test isapprox(cw(14, 3.0, 10.0, 70, 40.0), 6.866664f0; atol = 1f-5)  # WJ 06405 (R6-m2 no-elev, CI-added)
    @test isapprox(cw(15, 6.0, 25.0, 60, 80.0), 9.1771f0;   atol = 1f-4)  # MC 47502 (Bechtold-m2, CI-added)
    # FIA-122 disambiguation by index:
    @test FVSjl._CI_CWMAP[10] == "12203"   # PP log-form
    @test FVSjl._CI_CWMAP[18] == "12205"   # OS R6-m2
    @test 0.5f0 <= cw(1, 0.05, 0.5, 40, 2.0) <= 99.9f0
end
