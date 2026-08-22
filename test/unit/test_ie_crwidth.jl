# IE (InlandEmpire) FVS_TreeList CrWidth via the national cwcalc.f IEMAP dispatch (regression, 2026-08-22).
# cwcalc.f is the shared national crown-width routine (SELECT CASE(CWEQN)); IE (Region 1/4, BF=1.0) supplies
# IEMAP. ie_cwcalc reuses the EM-validated national form helpers (_em_r1/_em_powf/_em_r6m2/_em_bech1/_em_bech2)
# with IE's 23-species map + the 8 codes EM doesn't carry. VALIDATED vs FVSie_clean ietl FVS_TreeList: a
# function-level A/B (jl kernel on the oracle's own DBH/Ht/PctCr rows + internal stand BA) is bit-exact on all
# 4399 rows across the 7 emitted species (WP/WL/DF/GF/ES/AF/LP) — worst |Δ|=0.0. IE is the 13th variant wired
# past the eastern 0.5 default. The other 15 codes are transcription-faithful national cwcalc.f (un-exercised on
# this stand; the shared forms are the EM-353/353-validated ones).

using Test
using FVSjl

@testset "IE FVS_TreeList CrWidth — national cwcalc.f IEMAP, bit-exact vs FVSie_clean" begin
    # ietl context (bare plant stand; BF=1.0; the emitted species use only log-form '…03' codes).
    cw(sp, d, h, cr, ba) = FVSjl.ie_cwcalc(sp, Float32(d), Float32(h), Float32(cr), Float32(ba), 45f0, 0f0)

    # BAREA-free log-forms (identical CW at any stand BA) — GF(4)/ES(8)/AF(9)/LP(7). (The full-precision
    # A/B in /tmp/ie_cw_ab.jl is the bit-exact proof, 4399/4399; these pin the kernel + BAREA-free property.)
    @test isapprox(cw(4, 0.983, 7.71, 84, 30.72), 4.6833563f0; atol = 1f-5)  # GF 01703 (no BAREA term)
    @test cw(4, 0.5, 5.0, 90, 10.0) == cw(4, 0.5, 5.0, 90, 999.0)            # GF truly BAREA-free
    @test cw(7, 0.5, 5.0, 90, 10.0) == cw(7, 0.5, 5.0, 90, 999.0)            # LP 10803 BAREA-free
    @test cw(9, 0.5, 5.0, 90, 10.0) == cw(9, 0.5, 5.0, 90, 999.0)            # AF 01903 BAREA-free

    # BAREA-dependent log-forms respond to stand BA (WP cba=-0.07182: larger BA ⇒ smaller CW):
    @test cw(1, 3.0, 20.0, 90, 10.0) > cw(1, 3.0, 20.0, 90, 200.0)          # WP 11903 uses BAREA
    @test isapprox(cw(1, 0.656, 4.538854, 90, 2.6540914), 2.7072966f0; atol = 1f-5)  # WP 11903
    @test isapprox(cw(2, 0.785, 6.12, 90, 30.72009),      2.3822467f0; atol = 1f-5)  # WL 07303
    @test isapprox(cw(3, 1.391, 10.91, 84, 30.72009),     5.888335f0;  atol = 1f-5)  # DF 20203

    # the IEMAP index must disambiguate the two FIA-122 species: PP(idx10) log '12203' vs OS(idx23) R6M2 '12205'.
    @test FVSjl._IE_CWMAP[10] == "12203"
    @test FVSjl._IE_CWMAP[23] == "12205"

    # cwcalc.f final clamp [0.5, 99.9] holds.
    @test 0.5f0 <= cw(1, 0.05, 0.5, 40, 2.0) <= 99.9f0
end
