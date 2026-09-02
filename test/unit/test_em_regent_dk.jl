# test_em_regent_dk.jl — EM CRVAR/UTVAR small-tree height→diameter DK model (em/regent.f DO-25, lines
# 914-1012). Regression guard for the dense-aspen/cottonwood under-growth fix: seedlings that grow past
# breast height (hk≥4.5') must derive DIAMETER from the inverse-Wykoff HD curve, NOT keep the ~0 large-tree
# DG that froze them (measured stand 488929762126144: BA 62→244-class after the fix).
using Test
using FVSjl

@testset "EM regent CRVAR/UTVAR height→diameter DK model" begin
    dk = FVSjl._em_crut_dg
    # Aspen (sp 12): d=1.8", h=12', blended new height hk=15.09' (htg=3.09), bark=1.0, dgmx=2.5, scale=1.
    # Hand-computed from DK=(HT2/(ln(hk-4.5)-HT1))-1 with HT1(12)=4.4421, HT2(12)=-6.5405:
    #   DK=2.1412, DKK(h=12)=1.6947, DGK=0.44655 → DDS transform → DG≈0.4468.
    dg = dk(12, 1.8f0, 12.0f0, 15.09f0, 3.09f0, 0.0f0, 1.0f0, 2.5f0, 1.0f0)
    @test isapprox(dg, 0.4468f0; atol = 2.0f-3)      # matches the FVSem_g16 per-tree dump magnitude (~0.43)

    # Frozen-diameter regression guard: a tall seedling (0.1" dbh grown to ~14') must get SUBSTANTIAL DG,
    # never the ~0 it had before the fix. (bark=0.9 realistic.)
    dg_seed = dk(12, 0.1f0, 10.0f0, 14.0f0, 4.0f0, 0.0f0, 0.9f0, 2.5f0, 1.0f0)
    @test dg_seed > 0.3f0

    # Monotone in new height: taller → more catch-up diameter.
    @test dk(12, 1.5f0, 11.0f0, 18.0f0, 7.0f0, 0.0f0, 0.9f0, 2.5f0, 1.0f0) >
          dk(12, 1.5f0, 11.0f0, 13.0f0, 2.0f0, 0.0f0, 0.9f0, 2.5f0, 1.0f0)

    # DGMAX cap honoured (tiny cap clamps DGK before the DDS transform ⇒ small DG).
    @test dk(12, 0.1f0, 10.0f0, 30.0f0, 20.0f0, 0.0f0, 0.9f0, 0.01f0, 1.0f0) < 0.2f0

    # Juniper (sp 6) SITEAR-linear branch: DK=(hk-4.5)*10/(sitear-4.5); positive growth for a tall seedling.
    @test dk(6, 0.2f0, 8.0f0, 12.0f0, 4.0f0, 30.0f0, 0.9f0, 2.0f0, 1.0f0) > 0.0f0
end
