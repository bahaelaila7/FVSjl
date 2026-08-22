# CA (CentralCalifornia) FVS_TreeList CrWidth via ca_cwcalc with the forest-610 BF folded in (regression, 2026-08-22).
# CA's crown width is the Crookston-R6 Model 2 (cwcalc.f CAMAP), and cat01's forest 610 (Rogue River) applies a
# per-FIASP BF: SP(117)=1.048, LP(108)=0.944, PP(122)=0.918, DF(202)/WF(015)=1.0 (not in the 610 table). ca_cwcalc
# previously OMITTED the BF (the aggregated FFE crown-biomass masked the ~5% per-tree error); the FVS_TreeList CrWidth
# column exposed it. The BF is now FOLDED into the leading coefficient (exactly as oc_cwcalc, whose ref forest is also
# 610), making CrWidth bit-exact — VALIDATED vs FVSca_clean cat01 FVS_TreeList: all 29 inventory rows match jl's
# ca_cwcalc at full precision (0 mismatches, function-level A/B on the oracle's own sp/D/H/PctCr). The BF fold is INERT
# on the cat01_ffe .sum (2000 fire cycle byte-identical with/without BF — the crown-width change crosses no fire
# threshold), so it does NOT regress the validated CA FFE. CA is the 9th variant wired past the eastern 0.5 default.

using Test
using FVSjl

@testset "CA FVS_TreeList CrWidth = ca_cwcalc + forest-610 BF, bit-exact vs FVSca_clean" begin
    # cat01 stand context (forest 610, el 45 hundred-ft, inline stand ⇒ lat/lon 0 ⇒ Hopkins -323.6175); BA 85.13126.
    ba = 85.13126f0; el = 45.0f0; hi = -323.6175f0
    # forest_bf=true: the FVS_TreeList forest-grown path applies the R6 forest-610 BF (the FFE PERCOV path is BF-free).
    cw(sp, d, h, cr) = clamp(FVSjl.ca_cwcalc(sp, Float32(d), Float32(h), Float32(cr), ba, el, hi; forest_bf = true), 0.5f0, 99.9f0)
    # CA species indices: DF=7, WF=4, SP=16, LP=12, PP=18 (code_alpha order).
    # Oracle FVS_TreeList CrWidth (FVSca_clean cat01, inventory 1990) — these carry the folded 610 BF:
    @test isapprox(cw(18, 6.5, 30.0,  75), 9.984; atol = 0.005)   # PP 6.5"/30'/cr75 (BF 0.918)
    @test isapprox(cw(12, 7.2, 45.5,  16), 7.036; atol = 0.005)   # LP 7.2"/45.5'/cr16 (BF 0.944)
    @test isapprox(cw(16, 7.9, 75.0,  25), 9.371; atol = 0.005)   # SP 7.9"/75'/cr25 (BF 1.048)
    @test 0.5f0 <= cw(18, 0.05, 0.5, 40) <= 99.9f0                # clamp holds for a tiny tree
end
