# NC (Klamath) habitat-driven SDIMAX — alpha plant-association PV_CODE decode, regression 2026-08-25.
#
# Same alpha-PV_CODE decode-bug class as EC/PN/WC/CA/SO. NC's R6 forests (611 Siskiyou = IFOR 4, 712 BLM
# Coos Bay = IFOR 7) take the ECOCLS PA-specific per-species SDImax; the FIA reader never decoded NC's ALPHA
# PV_CODE and nc_sitset! was a first-pass that seeded a provisional uniform default (~720) instead of the
# stand's ecoclass. So the SDIMAX was wrong on EVERY R6 NC stand — even the default. This wires both: the
# reader decode (PV_CODE → KODTYP into NC_PCOML) and the faithful nc/sitset.f ECOCLS PA seed + C6-ratio fan
# (cap FORMAX 850).
#
# MEASURED vs live FVSnc_clean (LOCATION 611 = IFOR 4 = R6):
#   default CN 1127525637290487 (PV_CODE HTS121 not in PCOML ⇒ CWC221 default): live/jl "SDI MAX" =
#     850 850 815 850 850 850 569 850 850 850 (site species DF(3) 815, C6-fan capped 850) — BIT-EXACT
#     (was jl uniform 720 before the fix).
#   non-default CN 41135977010497 (PV_CODE HTS221, no ref, PCOML[81]): live "PLANT COMMUNITY CODE USED IS
#     HTS221"; live/jl "SDI MAX" = 850 850 830 850 850 850 580 850 850 850 (DF(3) 830) — BIT-EXACT.
#
# Self-contained unit test of the decode invariants (no DB / no oracle):
using Test
using FVSjl

@testset "NC habitat PV_CODE decode → NC_PCOML KODTYP (fixes R6 SDIMAX)" begin
    # No reference code ⇒ HBDECD matches the PV_CODE directly against PCOML.
    @test FVSjl.nc_habitat_kodtyp("HTS221", "") == 81   # PCOML[81]; ecocls site sp DF(3) RSDI 830
    @test FVSjl.nc_habitat_kodtyp("CWC221", "") == 46   # the R6 default PA
    # A PV_CODE not in PCOML (e.g. HTS121) ⇒ 0 ⇒ CWC221 default (matches live's CWC221).
    @test FVSjl.nc_habitat_kodtyp("HTS121", "")    == 0
    @test FVSjl.nc_habitat_kodtyp("HTS121", "641") == 0
    @test FVSjl.nc_habitat_kodtyp("ZZZZZZ", "999") == 0
    @test FVSjl.nc_habitat_kodtyp("", "") == 0

    # habitat_code → PA feeds nc_sitset!'s ECOCLS SDImax seed + C6 fan.
    @test FVSjl.nc_habtyp(81) == "HTS221"
    @test FVSjl.nc_habtyp(46) == "CWC221"
    @test FVSjl.nc_habtyp(0)  == "CWC221"               # unresolved ⇒ default
    @test FVSjl.nc_ecocls("HTS221")[1].sdimx == 830f0   # non-default PA SDImax the decode restores
    @test FVSjl.nc_ecocls("CWC221")[1].sdimx == 815f0   # the default the fan starts from (site species DF)
    @test FVSjl.nc_ecocls("CWC221")[1].fvsseq == 3
    # C6-ratio fan invariant: species 7 (BO) = DF SDImax · C6(7)/C6(3), the rest cap at FORMAX 850.
    @test round(815f0 * FVSjl.NC_C6[7] / FVSjl.NC_C6[3]) == 569f0
    @test round(830f0 * FVSjl.NC_C6[7] / FVSjl.NC_C6[3]) == 580f0
end
