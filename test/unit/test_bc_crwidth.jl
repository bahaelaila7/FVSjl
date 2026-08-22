# BC (BritishColumbia) FVS_TreeList CrWidth via national cwcalc.f BCMAP dispatch (regression, 2026-08-22).
# BC is metric but cwcalc runs in imperial internally; bc_cwcalc = _cwcalc_national(BCMAP[sp], …) with imperial
# inputs (BF=1.0). BCMAP needs 0 new codes (all shared with EM/IE/CI/AK). VALIDATED vs FVSbc_clean bctl
# FVS_TreeList_Metric (converted metric→imperial): log-form species PW/FD bit-exact (800/800 each); R6-m2 species
# (PL 10805) bit-exact at the oracle's ELEV≈1.476 (STDINFO 45 m → hundreds-ft). The metric output layer converts
# feet→m downstream (CrWidth·0.3048). BC is the 19th (final western) variant; jl does not yet emit a metric
# TreeList so this dispatch is additive/inert but correct-when-emitted.

using Test
using FVSjl

@testset "BC FVS_TreeList CrWidth — national cwcalc.f BCMAP, bit-exact-or-cornered vs FVSbc_clean" begin
    cw(sp, d, h, cr, ba, el) = FVSjl.bc_cwcalc(sp, Float32(d), Float32(h), Float32(cr), Float32(ba), Float32(el), 0f0)
    @test isapprox(cw(1, 2.0, 15.0, 80, 10.0, 1.476), 5.738732f0; atol = 1f-5)  # PW 11903 (log, BAREA) — bit-exact class
    @test isapprox(cw(3, 3.0, 20.0, 85, 30.0, 1.476), 9.124165f0; atol = 1f-5)  # FD 20203 (log, BAREA) — bit-exact class
    @test isapprox(cw(7, 3.0, 18.0, 75, 25.0, 1.476), 8.514159f0; atol = 1f-5)  # PL 10805 (R6-m2, EL) — cornered on EL feed
    # BCMAP: 0 new codes; index-keyed (OC idx14 reuses DF's log '20203'):
    @test FVSjl._BC_CWMAP[14] == "20203"
    @test 0.5f0 <= cw(1, 0.05, 0.5, 40, 2.0, 1.476) <= 99.9f0
end
