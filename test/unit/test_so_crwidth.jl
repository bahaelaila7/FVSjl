# SO (SouthCentralOregon) FVS_TreeList CrWidth — forest-601 BF (already in so_cwcalc) + the Hopkins-index fix
# (regression, 2026-08-22). so_cwcalc already folds the forest-601 (Deschutes) BF into its R6-Model-2 coefs, so the
# R6M2 species were bit-exact; but the 3 Bechtold-2004 species (GC/MC/MB) use the Hopkins bioclimatic index HI, which
# jl computed from lat/lon=0 (HI=-323.6, clamped) instead of so/grinit.f's default TLAT=42/TLONG=121 (HI≈-4.37) —
# so their FVS_TreeList CrWidth was off 1-2.6. FIX: so_grinit! now defaults the stand lat/lon to 42/121 when unset,
# feeding the correct HI. VALIDATED vs FVSso_clean sot01 FVS_TreeList: all 33 inventory species bit-exact (0
# mismatches). Gate 339/11 (lat/lon inert on the growth .sum — DGF uses elevation/site, not HI). SO is the 11th
# variant wired past the eastern 0.5 default.

using Test
using FVSjl

@testset "SO FVS_TreeList CrWidth — forest-601 BF + Hopkins-index fix, bit-exact vs FVSso_clean" begin
    # sot01 context (forest 601 Deschutes; lat/lon default 42/121 ⇒ HI=-4.367; el 45; BA 101.13126).
    ba = 101.13126f0; el = 45.0f0; hi = -4.367f0
    cw(sp, d, h, cr) = clamp(FVSjl.so_cwcalc(sp, Float32(d), Float32(h), Float32(cr), ba, el, hi), 0.5f0, 99.9f0)
    # the 3 Bechtold-2004 species the Hopkins fix corrected (SO species indices GC=29, MC=30, MB=31):
    @test isapprox(cw(29, 6.6, 30.0, 65), 12.451; atol = 0.005)   # GC (63102, HI clamp [-55,15])
    @test isapprox(cw(30, 5.8, 28.0, 65),  9.193; atol = 0.005)   # MC (47502, HI clamp [-37,27])
    @test isapprox(cw(31, 5.0, 25.0, 25),  8.504; atol = 0.005)   # MB (→MC 47502)

    # the Hopkins default is wired in so_grinit!: a fresh SO stand carries lat/lon 42/121 (not 0).
    s = FVSjl.StandState(FVSjl.SouthCentralOregon())
    FVSjl.init_blockdata!(s, s.variant)
    FVSjl.so_grinit!(s)
    @test s.plot.latitude == 42f0 && s.plot.longitude == 121f0
end
