# test_west_smhtgf_live.jl — the western small-tree height increment (SMHTGF) in the build's own REAL*4 arithmetic.
#
# gfortran compiles EXP/ALOG/** of REAL*4 to glibc expf/logf/powf (FMath fexp/flog/fpow), and the statements are
# evaluated as written (ca/smhtgf.f:152-155 SMHMOD=1.016605*CRMOD*RHMOD; HTGR=DOMHTGR*SMHMOD). Goldens: HTGRR as
# REGENT(LESTB) receives it from SMHTGF in the live oracles (private FVSca_g16 / FVSws_g16 builds with a WRITE after
# the SMHTGF call), for the planted DF/WF seedlings of the tiered plant_cyc stands CA 15320267010497 and
# WS 15353585010497 at their planting cycle; inputs as the oracle printed them.

using Test, FVSjl
const _SM = FVSjl

@testset "CA/WS SMHTGF == live (expf/logf/powf, statement order)" begin
    # CA sp 3 (MAPSP 2, the fir/DF Ritchie–Hann form): (H, CR=ICR/10, RELHT) at BA=BAL=24.483538, SI=72
    for (h, cr, relht, live) in ((2.6090896f0, 8.4f0, 0.082533956f0, 1.9879019f0),
                                 (2.3077304f0, 8.4f0, 0.07300099f0,  1.9076823f0),   # native exp/^: 1.9076821
                                 (2.4261236f0, 9.0f0, 0.076746151f0, 1.9522065f0))   # native exp/^: 1.9522064
        @test _SM.ca_smhtgf(3, 0.1f0, h, cr, 24.483538f0, 24.483538f0, 72f0, relht) === live
    end
    # WS sp 3 (the fir group): CR=ICR/10 at BA=BAL=227.47774, SI=54
    for (cr, live) in ((8.3f0, 3.1913865f0), (9.0f0, 3.4446633f0), (8.6f0, 3.2937050f0))
        @test _SM.ws_smhtgf(3, 0.1f0, cr, 227.47774f0, 227.47774f0, 54f0, 2.5f0) === live
    end
end
