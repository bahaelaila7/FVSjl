# NC (Klamath) Region-6 (SISKIYOU, IFOR 4, forest 611) volume — the VEQNNC forest/region fix (2026-09-05).
#
# NC's VEQNNC is forest-dependent (sitset.f → VOLEQDEF(VAR='NC',IREGN=KODFOR/100,FORST)). Region-5 forests use
# the 500WO2W/500DVEW table; SISKIYOU (a Region-6 forest) instead gets westside Flewelling F06FW2W202 (DF), INGY
# I00FW2W093/073 (WF/PP), and Region-6 Behre 616BEHW<fia> (all else). jl previously hardcoded the R5 table for
# ALL forests, so on Siskiyou stands large sound DF was computed with the fatter R5 R5TAP taper — TOTAL CUBIC
# over-predicted +7% at D=20 → +30% at D=54 (hidden: nct01 is a Region-5 stand with D≤35). This fix routes
# IFOR-4 stands through the R6 kernels (already ported for WC: wc_fw2_westside_vol; BM/WC: bm_r6vol3 Behre).
# Oracle values captured by instrumenting FVSnc_g16's TCUBIC on forest 611 (per-tree TVOL1).

using Test
using Printf
using FVSjl
const _M = FVSjl

@testset "NC Region-6 (Siskiyou) volume — total cubic bit-exact vs FVSnc_g16 forest 611" begin
    # --- Westside Flewelling DF (F06FW2W202): the large-tree over-prediction fix ---
    # Per-tree TVOL1 (total cubic), sound DF, oracle FVSnc_g16 forest 611. bark = printed BARK (4-dp).
    df = [ (19.8f0, 121f0, 0.8286f0, 81.2f0),
           (33.3f0, 195f0, 0.8308f0, 329.9f0),
           (40.0f0, 212f0, 0.8313f0, 489.9f0),
           (54.4f0, 221f0, 0.8320f0, 827.4f0) ]
    for (d, h, bark, orc) in df
        v = _M.wc_fw2_westside_vol("F06FW2W202", d, h, bark; topd = 6.0f0, bftopd = 6.0f0)[1]
        @test isapprox(v, orc; atol = 0.15f0)              # ≤ one 0.1-rounding tick (giant-tree bark knife-edge)
    end
    # The R6 Flewelling REPLACES the R5 R5TAP over-prediction: the old path was systematically OVER, the new one
    # converges to the oracle across the whole D range (ratio 1.074→1.300 OLD ⇒ 1.000 NEW).
    for (d, h, bark, orc) in df
        old = _M.nc_wo2w_vol("500WO2W202", d, h)[1]        # the buggy R5TAP path
        new = _M.wc_fw2_westside_vol("F06FW2W202", d, h, bark; topd = 6.0f0, bftopd = 6.0f0)[1]
        @test old > orc * 1.05f0                            # OLD over-predicts (≥ +5%, up to +30% at D=54)
        @test abs(new - orc) <= 0.15f0                      # NEW is bit-exact-or-one-tick
        @test new < old                                     # NEW strictly reduces the fat R5TAP stem
    end

    # --- Region-6 Behre (616BEHW), hardwoods — total cubic bit-exact (bm_r6vol3 + NC SISKFC form class) ---
    beh = [ (5,  6.2f0, 42.0f0,  0.9794f0, 5.995f0),        # madrone
            (5,  14.5f0, 57.56f0, 0.9806f0, 35.098f0),      # madrone
            (11, 22.4f0, 72.0f0,  0.9457f0, 81.444f0),      # laurel (other hardwoods)
            (11, 12.9f0, 45.0f0,  0.9369f0, 21.881f0),      # laurel
            (11, 4.4f0, 27.0f0,   0.8967f0, 1.146f0) ]      # laurel small
    for (sp, d, h, bark, orc) in beh
        v = _M.nc_behre_vol(sp, d, h, bark)[1]
        @test isapprox(v, orc; atol = 0.05f0)              # bit-exact
    end

    # --- Form class table (formcl.f SISKFC): DBH-class lookup ---
    @test _M.nc_siskfc(3, 5.0f0) == 90                      # DF class 1 (D<11)
    @test _M.nc_siskfc(3, 25.0f0) == 81                     # DF class 3 (21≤D<31)
    @test _M.nc_siskfc(3, 54.4f0) == 80                     # DF class 5 (D>40.9)
    @test _M.nc_siskfc(5, 6.2f0) == 98                      # MA class 1

    # --- Region-5 path UNCHANGED (no regression): the WO2W R5TAP DF stays as-is on R5 forests ---
    @test _M.nc_wo2w_vol("500WO2W202", 20.9f0, 130f0)[1] > 0f0   # R5 path still produces volume
end
