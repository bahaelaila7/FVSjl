# SO (SouthCentralOregon) habitat-driven SDIMAX — alpha plant-association PV_CODE decode, regression 2026-08-25.
#
# Same alpha-PV_CODE decode-bug class as EC/PN/WC/CA: the FIA DB stores the stand plant association as an
# ALPHA ecoclass code (FVS_STANDINIT_COND.PV_CODE, e.g. SO "CWS313"). SO's R6 forests (IFOR≤3 or 10) take
# the ECOCLS PA-specific per-species SDImax; the FIA reader never decoded PV_CODE, so so_sitset! seeded every
# stand's site-species SDImax from the CPS111 default RSDI (285) regardless of the real PA. On a non-default
# PA the SDIMAX ran wrong (CWS313 → 810, not 285 = ~2.8× low) and the density self-thin fired at the wrong
# level. Most SO alpha stands crosswalk to CPS111 (= the default) so were already right — only the minority
# on other PAs were affected.
#
# MEASURED vs live FVSso_clean, CN 24397775010900 (PV_CODE CWS313, no ref, LOCATION 602 = IFOR 2 = R6):
#   live "PLANT COMMUNITY CODE USED IS CWS313"; live .out "SDI MAX" = 549 549 850 810 850 549 665 810 850
#   527 615 810 ... — jl AFTER the decode reproduces that array bit-exact (site species WF(4) SDImax
#   285 → 810). The default CPS111 stand (CN 1127523558290487, CPS211/653 → no crosswalk → default) is
#   UNCHANGED at 285 (site species PP(10)) — decode is inert for default-PA stands.
#
# Self-contained unit test of the decode invariants (no DB / no oracle):
using Test
using FVSjl

@testset "SO habitat PV_CODE decode → SO_PCOML KODTYP (fixes non-default-PA SDIMAX)" begin
    # No reference code ⇒ HBDECD matches the PV_CODE directly against PCOML.
    @test FVSjl.so_habitat_kodtyp("CWS313", "") == 89   # PCOML[89]; ecocls site sp WF(4) RSDI 810
    @test FVSjl.so_habitat_kodtyp("CPS111", "") == 49   # the R6 default PA
    # PVREF6 crosswalk: full (PV_CODE, PV_REF_CODE) match → HABPVR → HBDECD.
    @test FVSjl.so_habitat_kodtyp("CPS211", "605") == 52  # (CPS211,605) → CPS211 identity → PCOML[52]
    # A reference present but no full-match pair ⇒ KARD2 blank ⇒ 0 ⇒ CPS111 default (matches live CPS111).
    @test FVSjl.so_habitat_kodtyp("CPS211", "653") == 0
    @test FVSjl.so_habitat_kodtyp("ZZZZZZ", "999") == 0
    @test FVSjl.so_habitat_kodtyp("", "") == 0

    # habitat_code → PA feeds so_sitset!'s ECOCLS SDImax seed.
    @test FVSjl.so_habtyp(89) == "CWS313"
    @test FVSjl.so_habtyp(0)  == "CPS111"               # unresolved ⇒ default
    @test FVSjl.so_ecocls("CWS313")[1].sdimx == 810f0   # the non-default PA SDImax the decode restores
    @test FVSjl.so_ecocls("CWS313")[1].fvsseq == 4      # site species WF
    @test FVSjl.so_ecocls("CPS111")[1].sdimx == 285f0   # the wrong default it replaces (site species PP)
    @test FVSjl.so_ecocls("CPS111")[1].fvsseq == 10
end
