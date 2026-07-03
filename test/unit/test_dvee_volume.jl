# D35: the R9 Gevorkiantz '900DVEE' volume model (VOLUME keyword METHC=5, CS).
# Ground truth = live FVScs R9VOL debug-stamp on the cst01_method5 BARE-GROUND-PLANT
# stand (shortleaf pine, fia 110), captured 2026-07-03. Total cubic + pulp merch cubic
# must reproduce live to Float32 precision; the merch cubic must gate off HT2PRD (0 for
# stems below the merch top).
@testset "DVEE (R9 Gevorkiantz) volume vs live stamp" begin
    # (fia, dbh, httot, live V1, live V4); si=60, ba=90, iforst=5 (CS forest 905)
    cases = ((110, 4.582f0, 38.40f0, 1.8470f0, 0.7682f0),
             (110, 3.588f0, 29.41f0, 0.8674f0, 0.0f0),     # below merch top ⇒ V4=0
             (110, 4.670f0, 38.71f0, 1.9344f0, 0.7984f0))
    for (fia, d, h, v1, v4) in cases
        tcf, mcf, scf, bf = FVSjl.r9vol_gevorkiantz(fia, d, h, 5)  # iforst=5 (CS); SI/BA default like fvsvol
        @test isapprox(tcf, v1; atol = 5f-3)     # form-factor total cubic 0.42π·D²·H/576
        @test isapprox(mcf, v4; atol = 5f-3)     # pulp merch cubic (0 when HT2PRD=0)
        @test scf == 0f0 && bf == 0f0            # sawtimber/board not ported (sub-9" range)
    end
    # HT2PRD (R9_MHTS) merch-height log count: 0 for the small stem, 1 for the merch one
    @test FVSjl._r9_mhts_ht2prd(110, 3.588f0, 29.41f0, 60, 90, "CS") == 0
    @test FVSjl._r9_mhts_ht2prd(110, 4.582f0, 38.40f0, 60, 90, "CS") == 1
end
