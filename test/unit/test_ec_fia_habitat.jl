# EC (EastCascades) FIA-DB alpha plant-association PV_CODE decode + SDIDEF(ISEQ) seed, regression 2026-08-25.
#
# EC's STDINFO-KEYWORD alpha decode was fixed in 70535053 (ec_hbdecd), but the FIA-DB PV_CODE column had no
# reader branch ⇒ habitat_code stayed 0 ⇒ ec_sitset! grew every EC FIA stand on the poor default site. This
# wires ec_habitat_kodtyp (ec/habtyp.f PVREF6 crosswalk + ec_hbdecd PCOML match) into the FIA reader, and
# restores the omitted ec/sitset.f:155-156 SDIDEF(ISEQ)=RSDI seed (surfaced once a real FIA PA whose ecoclass
# site species (DF) differs from the stand's DB site species (PP) was decoded).
#
# MEASURED vs live FVSec_clean, CN 1143082551290487 (PV_CODE CDG131 / ref 626, LOCATION 613 = R6):
#   live "HABITAT/PV_CODE: CDG131"; live .out "SDI MAX" = 530 532 530 629 629 629 554 629 575 530 (per
#   species); jl BEFORE = uniform 331 (habitat_code 0 default), AFTER = 530 532 530 629 629 629 554 629 575
#   530 — BIT-EXACT (habitat_code 11, site species PP(10), DF(3) SDImax 530 seeded, others fanned/capped).
#
# Self-contained unit test of the decode invariants (no DB / no oracle):
using Test
using FVSjl

@testset "EC FIA-DB habitat PV_CODE decode → EC_PCOML KODTYP + SDIDEF(ISEQ) seed" begin
    # PVREF6 crosswalk (ref present) then ec_hbdecd PCOML match.
    @test FVSjl.ec_habitat_kodtyp("CDG131", "626") == 11    # → CDG131 (PCOML[11]), ecocls DF SDImax 530
    @test FVSjl.ec_habitat_kodtyp("CDS715", "627") == 39    # → CDS715 (PCOML[39])
    # No reference code ⇒ HBDECD matches the PV_CODE directly (the STDINFO-keyword code path).
    @test FVSjl.ec_habitat_kodtyp("CWS222", "") == 129      # the 70535053 keyword-fix PA
    @test FVSjl.ec_habitat_kodtyp("CDG131", "") == 11
    # Unrecognized / blank ⇒ 0 ⇒ ec_sitset! CPS241 default.
    @test FVSjl.ec_habitat_kodtyp("ZZZZZZ", "999") == 0
    @test FVSjl.ec_habitat_kodtyp("", "") == 0

    # The decoded index feeds ec_sitset! via ec_habtyp/ec_ecocls.
    @test FVSjl.ec_habtyp(11) == "CDG131"
    @test FVSjl.ec_ecocls("CDG131")[1].sdimx == 530f0
    @test FVSjl.ec_ecocls("CDG131")[1].fvsseq == 3          # ecoclass site species DF (≠ a PP-site FIA stand)
end
