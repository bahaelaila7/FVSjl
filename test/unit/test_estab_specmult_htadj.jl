# SPECMULT (esin.f opt 7 → esnutr.f XESMLT) + HTADJ (opt 15 → esnutr.f HTADJ) ESTAB-packet keywords
# (2026-08-22). SPECMULT sets a per-species establishment-occupancy multiplier (occ = OCURHT·XESMLT·OCURNF);
# HTADJ sets a per-species height adjustment added to the established-tree height before the XMIN/HHTMAX clamps.
# Both are parsed INSIDE the ESTAB packet (kw_estab!), default inert (XESMLT 1.0 / HTADJ 0.0).
#
# HTADJ validated bit-exact-or-cornered vs FVSci_clean (cihtadj vs cibase FVS_TreeList): the targeted species
# (DF, HTADJ +5.0) shifts +5.0 ft (jl Δmean 4.993, max 5.0 = oracle 4.9976/5.0); the untargeted species (WP,
# hadj=0) is unchanged to a 0.2% growth-realization straddle (jl WP Δ 0.015, oracle 0.0 — cornered class).
# SPECMULT threading is proven deterministic here: ie_espadv[sp] (advance-regen prob) scales linearly with
# occ[sp] (XESMLT), ratio 2.0. Its end-to-end established-TPA effect is the STOCKADJ-class regen-regime straddle
# (advance-regen-dominated stands; masked on ingrowth-dominated fixtures — the documented establishment trap).

using Test
using FVSjl
using FVSjl: KeywordReader, read_keyword!, kw_estab!, StandState, InlandEmpire, init_blockdata!,
             ie_estab_indices, ie_espadv, ie_ocurht, autoes_ocurnf

# Build a fixed 10-col FVS keyword line: keyword rpad-10, then each field rpad-10.
kwline(kw, fields...) = rpad(kw, 10) * join(rpad.(string.(fields), 10)) * "\n"

@testset "SPECMULT + HTADJ ESTAB-packet keywords" begin
    @testset "parsing: single / all / group species selectors" begin
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        body = "ESTAB     1992\n" *
               kwline("SPECMULT", 1992, 3, 2.0) *      # DF (sp 3) ×2.0
               kwline("HTADJ",    1992, 3, 5.0) *      # DF (sp 3) +5.0 ft
               kwline("SPECMULT", 1992, 0, 1.5) *      # ALL species ×1.5 (sel 0)
               "END\n"
        kr = KeywordReader(IOBuffer(body))
        rec = read_keyword!(kr)                        # ESTAB
        kw_estab!(s, rec, kr)
        # SPECMULT 0 (all) ran last ⇒ every species = 1.5, incl DF (overwrites the earlier 2.0 for sp 3):
        @test s.estab.spec_mult[Int32(3)] == 1.5f0
        @test s.estab.spec_mult[Int32(10)] == 1.5f0
        @test length(s.estab.spec_mult) == length(s.coef.code_alpha)
        # HTADJ only DF (sp 3):
        @test s.estab.ht_adj[Int32(3)] == 5.0f0
        @test !haskey(s.estab.ht_adj, Int32(1))
    end

    @testset "SPECMULT threads XESMLT into the advance-regen occupancy (occ = OCURHT·XESMLT·OCURNF)" begin
        idx = ie_estab_indices(260, 116)              # ihab 3, ifo 16
        occ1 = Float32[i == 3 ? 1f0 : (i == 10 ? 1f0 : 0f0) for i in 1:23]
        occ2 = copy(occ1); occ2[3] = 2f0              # XESMLT(DF)=2.0
        over = zeros(Float32, 10)
        p1 = ie_espadv(idx.ihab, idx.iprep, idx.ifo, idx.iphy, 1f0, 0f0, 0.3f0, 1f0, 5f0, 45f0, 1f0, 0f0, 0f0, occ1, over)
        p2 = ie_espadv(idx.ihab, idx.iprep, idx.ifo, idx.iphy, 1f0, 0f0, 0.3f0, 1f0, 5f0, 45f0, 1f0, 0f0, 0f0, occ2, over)
        @test p1[3] > 0f0
        @test isapprox(p2[3] / p1[3], 2.0f0; atol = 1f-4)   # padv[DF] scales with XESMLT
    end

    @testset "defaults inert (empty dicts)" begin
        s = StandState(InlandEmpire()); init_blockdata!(s, s.variant)
        @test isempty(s.estab.spec_mult) && isempty(s.estab.ht_adj)
    end
end
