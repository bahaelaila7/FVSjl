# CA (CentralCalifornia) habitat-driven SDIMAX — alpha plant-association PV_CODE decode, regression 2026-08-25.
#
# The alpha-PV_CODE decode-bug class (already fixed for EC/PN/WC): several Region-5/6 variants store the stand
# plant association as an ALPHA ecoclass code in the FIA DB (FVS_STANDINIT_COND.PV_CODE, e.g. CA "CDH524"). The
# FIA reader never decoded CA's PV_CODE into the KODTYP index into CA_PCOML, so ca_sitset! fell back to the
# CWC221 default (ecocls SDImx 815) instead of the stand's own ecoclass (CDH524/641 → CDS511 → SDImx 635),
# inflating the BA-weighted stand SDIMAX ~1.3× so the Wykoff density self-thin under-fired on dense stands.
#
# MEASURED vs live FVSca_clean, CN 1127513363290487 (PV_CODE CDH524, ref 641, LOCATION 611 = R6):
#   live .out "SDI MAX" = 635 (all species); jl BEFORE = 815 (CWC221 default), AFTER = 635 (bit-exact).
#   cyc0 .sum was already bit-exact (SDI 444); by 2071 the un-decoded jl kept SDI 709 vs live's self-thinned
#   544 (TPA 884 vs 563) — the classic under-kill; decoding CDH524→CDS511 fires the self-thin toward live.
# Only R6 forests (KODFOR≥600, ~99.7% of CA's FIA population) go through ca/habtyp.f's PVREF6+PCOML path.
#
# Self-contained unit test of the decode invariants (no DB / no oracle):
using Test
using FVSjl

@testset "CA habitat PV_CODE decode → CA_PCOML KODTYP (fixes SDIMAX under-thin)" begin
    # PVREF6 crosswalk: (PV_CODE, PV_REF_CODE) full match → HABPVR → HBDECD PCOML index.
    @test FVSjl.ca_habitat_kodtyp("CDH524", "641") == 18   # → CDS511 (pcoml[18]), ecocls SDImx 635
    @test FVSjl.ca_habitat_kodtyp("CDS215", "641") == 18   # → CDS511 too
    @test FVSjl.ca_habitat_kodtyp("CDS121", "641") == 18
    @test FVSjl.ca_habitat_kodtyp("CPS211", "610") == 31   # → CPG141 (pcoml[31])
    @test FVSjl.ca_habitat_kodtyp("CDC521", "641") == 7    # identity crosswalk → CDC521 (pcoml[7])
    # No reference code ⇒ HBDECD matches the PV_CODE directly against PCOML.
    @test FVSjl.ca_habitat_kodtyp("CDS511", "") == 18
    # Unrecognized / blank ⇒ 0 ⇒ ca_habtyp default (CWC221). A reference present but no full-match pair
    # leaves KARD2 blank in ca/pvref6.f ⇒ 0 (default), NOT the raw PV_CODE.
    @test FVSjl.ca_habitat_kodtyp("HTS121", "641") == 0    # (HTS121,641) → blank HABPVR ⇒ default
    @test FVSjl.ca_habitat_kodtyp("ZZZZZZ", "999") == 0
    @test FVSjl.ca_habitat_kodtyp("", "") == 0

    # The decoded index feeds ca_sitset! via ca_habtyp: CDS511 (SDImx 635) replaces the CWC221 (815) default.
    @test FVSjl.ca_habtyp(18) == "CDS511"
    @test FVSjl.ca_habtyp(0)  == "CWC221"                  # undecodable ⇒ default
    @test FVSjl.ca_ecocls("CDS511")[1].sdimx == 635f0
    @test FVSjl.ca_ecocls("CWC221")[1].sdimx == 815f0      # the wrong default it replaces
end
