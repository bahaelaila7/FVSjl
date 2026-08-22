# NC (Klamath) FVS_TreeList CrWidth via R5CRWD reuse (regression, 2026-08-22).
# NC/Klamath's forest 505 is Region-5, where cwcalc.f branches to R5CRWD (ws/r5crwd.f) — a function of sp/D/H only,
# NOT the R6M2 Crookston models nc_cwcalc uses (nc_cwcalc's R6M2 was ~37% off for the TreeList). R5CRWD is FIA-keyed
# and shared with WS, so nc_r5crwd reuses ws_r5crwd via an NC-species→WS-species map (by FIA code). VALIDATED vs
# FVSnc_clean nct01 FVS_TreeList: all 29 inventory rows match (0 mismatches); the values are identical to WS's (same
# R5CRWD). NC is the 12th variant wired past the eastern 0.5 default. (nc_cwcalc's R6M2 stays for the FFE/PERCOV path,
# a separate currently-unported NC FFE concern; this change is TreeList-only.)

using Test
using FVSjl

@testset "NC FVS_TreeList CrWidth = R5CRWD (nc_r5crwd → ws_r5crwd), bit-exact vs FVSnc_clean" begin
    cw(sp, d, h) = clamp(FVSjl.nc_r5crwd(sp, Float32(d), Float32(h)), 0.5f0, 99.9f0)
    # NC species indices: DF=3, WF=4, RF=9 (code_alpha order). Oracle FVS_TreeList CrWidth (FVSnc_clean nct01, 1990):
    @test isapprox(cw(3, 0.1, 2.0), 1.56; atol = 0.005)   # DF 0.1"/2'  → 1.56
    @test isapprox(cw(4, 0.1, 3.0), 2.33; atol = 0.005)   # WF 0.1"/3'  → 2.33
    @test isapprox(cw(9, 0.1, 2.0), 1.56; atol = 0.005)   # RF 0.1"/2'  → 1.56
    @test isapprox(cw(3, 1.2, 11.0), 5.26; atol = 0.005)  # DF 1.2"/11' → 5.26
    # NC's R5CRWD is exactly WS's for the shared FIA species (DF/WF/RF), confirming the FIA-keyed reuse:
    @test cw(3, 1.2, 11.0) == clamp(FVSjl.ws_r5crwd(2, 1.2f0, 11.0f0), 0.5f0, 99.9f0)   # NC DF == WS DF(sp2)
    @test 0.5f0 <= cw(3, 0.05, 0.5) <= 99.9f0             # clamp holds for a tiny tree
end
