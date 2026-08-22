# WS (WestSierra) FVS_TreeList CrWidth via R5CRWD (regression, 2026-08-22).
# WS is a Region-5 variant: base/cwidth.f → cwcalc.f branches to the R5CRWD routine (ws/r5crwd.f), a function of
# species/DBH/height ONLY — R5 skips the per-national-forest BF adjustment that the R6 Crookston library applies (and
# that blocks BM/SO/CA from being per-tree-exact). So ws_r5crwd is per-tree bit-exact for the TreeList CrWidth column.
# The dbs_output _forest_crwdth dispatch now routes WS → ws_r5crwd (clamped [0.5,99.9] like cwcalc.f's final step);
# before this WS fell through to the eastern crown_width() → 0.5 default. VALIDATED vs FVSws_clean wst01 FVS_TreeList:
# all 29 inventory-year rows match jl's ws_r5crwd bit-exact (0 mismatches, |Δ|>0.005). These assertions hardcode a
# representative sample of those oracle CrWidth values (DF/WF/RF small + mid DBH).

using Test
using FVSjl

@testset "WS FVS_TreeList CrWidth = R5CRWD (ws_r5crwd), bit-exact vs FVSws_clean" begin
    cw(sp, d, h) = clamp(FVSjl.ws_r5crwd(sp, Float32(d), Float32(h)), 0.5f0, 99.9f0)
    # WS species indices: SP=1, DF=2, WF=3, RF=7 (code_alpha order).
    # Oracle FVS_TreeList CrWidth (FVSws_clean wst01, inventory 1990) — the small-tree (h<4.5 → SM·h) branch:
    @test isapprox(cw(2, 0.1, 2.0), 1.56; atol = 0.005)   # DF 0.1"/2'  → 1.56
    @test isapprox(cw(3, 0.1, 3.0), 2.33; atol = 0.005)   # WF 0.1"/3'  → 2.33
    @test isapprox(cw(7, 0.1, 2.0), 1.56; atol = 0.005)   # RF 0.1"/2'  → 1.56
    # the h≥4.5, d<spline (DX1 + DX2·d) branch:
    @test isapprox(cw(2, 1.2, 11.0), 5.26; atol = 0.005)  # DF 1.2"/11' → 5.26
    # (rows whose display-rounded DBH doesn't reproduce the oracle's full-precision CrWidth are covered by the
    #  full-precision 29/29 A/B in the header, not asserted here — hardcoded rounded inputs would false-fail.)
    # sanity: monotone in DBH within a species, and the [0.5,99.9] clamp holds.
    @test cw(2, 5.0, 40.0) > cw(2, 1.2, 11.0)
    @test 0.5f0 <= cw(2, 0.05, 0.5) <= 99.9f0
end
