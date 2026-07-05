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
    # atol 6f-4 = the STAMP-PRECISION floor (was padded 5f-3, ~8× the real residual). The live R9VOL
    # ground truth was captured as a 3-4 decimal debug print; measured jl residual ≤5.94e-4 (tcf@case3
    # 1.933898/1.9344, mcf@case2b 10.395594/10.395). That residual = the stamp's own 3-decimal round
    # (±5e-4) + a sub-ULP Float32 op-order tail in the form-factor (0.42π·D²·H/576) / Gevorkiantz
    # polynomial. Irreducible vs a printed stamp (no raw Fortran bits to compare `==`).
    for (fia, d, h, v1, v4) in cases
        tcf, mcf, scf, bf = FVSjl.r9vol_gevorkiantz(fia, d, h, 5)  # iforst=5 (CS); SI/BA default like fvsvol
        @test isapprox(tcf, v1; atol = 6f-4)     # form-factor total cubic 0.42π·D²·H/576
        @test isapprox(mcf, v4; atol = 6f-4)     # pulp merch cubic (0 when HT2PRD=0)
        @test scf == 0f0 && bf == 0f0            # sawtimber/board not ported (sub-9" range)
    end
    # HT2PRD (R9_MHTS) merch-height log count: 0 for the small stem, 1 for the merch one
    @test FVSjl._r9_mhts_ht2prd(110, 3.588f0, 29.41f0, 60, 90, "CS") == 0
    @test FVSjl._r9_mhts_ht2prd(110, 4.582f0, 38.40f0, 60, 90, "CS") == 1

    # SAWTIMBER merch = VOL(4)_saw + VOL(7)=PT·GCB (fvsvol.f MCF=TVOL(4)+TVOL(7)); ground truth from the
    # fvsvol.f MCF stamp. (fia 110, iforst=5). Both the '912' sawlog polynomial and the PT·GCB topwood.
    for (d, h, mcf_live) in ((9.004f0, 50.98f0, 8.689f0), (9.006f0, 70.36f0, 10.395f0))
        _, mcf, _, _ = FVSjl.r9vol_gevorkiantz(110, d, h, 5)
        @test isapprox(mcf, mcf_live; atol = 6f-4)   # stamp-precision floor (see above)
    end
end
