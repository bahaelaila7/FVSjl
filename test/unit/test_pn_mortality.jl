# PN (PacificNorthwest) + WC (WestCascades) habitat-driven SDIMAX / density self-thin, regression 2026-08-25.
#
# The western full-population FIA sweep flagged PN with the density-mortality under-kill phenotype: on dense
# over-max-SDI FIA stands the port UNDER-killed cycle-1 TPA ~3-4× (both-sides-traced vs live FVSpn_clean,
# CN 451024282489998: inventory 4997 → live 4997→1050 vs old jl →4178, live SDI 406 vs jl 777). PN/WC do NOT
# use the shared VARMRT/NWCMRT distributor — vwc/morts.f is the ORGANON per-tree logistic-RIP model with an
# integer-PASS SDI self-thin (scale every tree's kill by the smallest PASS that brings post-mortality SDIA <
# SDIMAX). The bug was NOT in mortality.jl: SDICAL's weighted SDIMAX ran ~2× high because the FIA reader
# never decoded the ALPHA plant-association PV_CODE ("CHS512") into the KODTYP index into PN_PCOML, so
# pn_sitset! fell back to the CHS133 default (SDIDEF 1606→FORMAX 950) instead of the stand's ecoclass SDIDEF
# (CHS512 → 485). With SDIMAX 2× high the PASS loop never fires and the stand keeps ~4× too much TPA.
#
# A SECOND latent bug surfaced once the real PA was decoded: pn/wc sitset.f:91-96 (SDIDEF(ISISP)=RSDI) was
# omitted, so on a PA whose ecoclass site species (e.g. WH) differs from the stand's site species (DF) the
# site species' SDIDEF stayed 0 and the DO-80 fill propagated 0 to every species ⇒ stand_sdimax→0 ⇒ morts.f
# "SDIMAX<5 ⇒ kill ALL trees" wiped the stand. Both are fixed (habitat decode + SDIDEF(isisp) seed).
#
# Signed sign-tally (cycle-1 TPA, sign of jl−live) vs live FVSpn_clean, 12 dense over-max-SDI FIA stands:
#   BEFORE (both bugs): one-directional (jl under-kills the dense stands / kills-all the WH-site ones)
#   AFTER : bit-exact-or-cornered (the lone residual is a ~0.1% grown-QMD straddle at the density cap).
#
# Self-contained unit test of the two decode invariants the fix restores (no DB / no oracle):
using Test
using FVSjl

@testset "PN/WC habitat PV_CODE decode → PCOML KODTYP (fixes SDIMAX under-thin)" begin
    @testset "PN pn_habitat_kodtyp (pn/habtyp.f + PVREF6 + HBDECD)" begin
        # Direct alpha plant-association string → PN_PCOML index (HBDECD path, no reference code).
        @test FVSjl.pn_habitat_kodtyp("CHS512", "") == 60     # pcoml[60] = CHS512, ecocls SDIDEF 485
        @test FVSjl.pn_habitat_kodtyp("CHS324", "") == 51     # WH-site PA (ecocls SDIDEF 865)
        @test FVSjl.pn_habitat_kodtyp("CHS136", "") == 42
        @test FVSjl.pn_habitat_kodtyp("CHS422", "") == 58
        # PVREF6 crosswalk: (PV_CODE, PV_REF_CODE) present ⇒ full match → HABPVR → HBDECD. (CHS512,618)→CHS512.
        @test FVSjl.pn_habitat_kodtyp("CHS512", "618") == 60
        # Unrecognized / blank ⇒ pn/habtyp.f ITYPE=40 default (PCOML[40] = CHS133).
        @test FVSjl.pn_habitat_kodtyp("CHS6", "")  == FVSjl.PN_HAB_DEFAULT == 40
        @test FVSjl.pn_habitat_kodtyp("", "")      == 40
        # A reference code present but no full-match pair ⇒ default (any partial match → ITYPE=40).
        @test FVSjl.pn_habitat_kodtyp("ZZZZZZ", "999") == 40
        # Numeric sequence-number fallback (HBDECD: IFIX(ARRAY2) in 1..NPA).
        @test FVSjl.pn_habitat_kodtyp("60", "") == 60
        # The decoded index feeds pn_sitset!: CHS512 → ecocls SDIDEF 485 (< the CHS133 950 default it replaces).
        @test FVSjl.pn_ecocls("CHS512")[1].sdimx == 485f0
        @test FVSjl.pn_ecocls("CHS133")[1].sdimx == 1606f0    # the wrong default (capped to FORMAX 950)
    end

    @testset "WC wc_habitat_kodtyp (wc/habtyp.f)" begin
        @test FVSjl.wc_habitat_kodtyp("CFS551", "") == 52     # WC default PA (pcoml[52])
        @test FVSjl.wc_habitat_kodtyp("", "")       == FVSjl.WC_HAB_DEFAULT == 52
        # A real WC alpha PA resolves to its own PCOML index (not the default) when present in the table.
        let idx = findfirst(!=("CFS551"), FVSjl.WC_PCOML)
            pa = FVSjl.WC_PCOML[idx]
            @test FVSjl.wc_habitat_kodtyp(pa, "") == idx
        end
    end
end
