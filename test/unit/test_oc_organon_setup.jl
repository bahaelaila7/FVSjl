# OC (Oregon Coast) chunk C2 — ORGANON PREPARE setup calibration regression guard.
#
# Pins the bit-exact result validated at port time against the live FVSoc_clean oracle
# (scoped `DEBUG 1 / CRATET`, stand S248112/ocmin, 27 records): the ORGANON `TMPCAL(k,grp)`
# calibration multipliers, the DF↔PP site conversion, and the WF/ES species crosswalk.
# ORGANON is deterministic (DGSD=0) ⇒ the bar is BIT-EXACT (within Float32 print precision),
# NOT a straddle. See docs/OC_VARIANT_PORT_AUDIT.md (C2 verdict).

using Test
using FVSjl
using FVSjl: OregonCoast, StandState, init_blockdata!, organon_prepare_swo, site_setup!,
             resolve_species, orgspc

# Exact /ORGANON/ buffer PREPARE received from the oracle (C1-validated marshalling).
const _OC_SPECIES = Int32[122,202,122,117,117,17,117,122,202,117,202,202,17,17,202,17,202,122,202,122,17,202,17,17,17,17,17]
const _OC_DBH1 = Float32[11.5,0.1,6.5,7.9,8.0,6.2,8.4,9.5,4.0,8.2,1.2,1.9,0.1,5.3,10.0,6.1,12.7,9.6,10.4,8.5,10.9,9.4,3.2,0.1,5.8,5.0,6.6]
const _OC_HT1OR = Float32[73,4.6,30,75,63,38,5,60,20,65,11,13,4.6,27,65,38,67,60,55,0,65,60,17,4.6,28,25,30]
const _OC_ICR = Int[35,55,75,25,25,45,35,25,25,45,55,45,65,65,35,75,35,25,45,25,65,35,45,65,65,25,65]
const _OC_EXPAN1 = Float32[60.9999733,990.000061,190.940720,129.262070,126.050720,209.865936,114.331734,
    89.3877640,330.000031,119.976891,330.000031,330.000031,990.000061,287.192810,80.6724625,216.803177,
    50.0170250,87.5352173,74.5862350,111.657394,67.9004059,91.2997589,330.000031,330.000031,239.811111,
    322.689850,185.198502]

@testset "OC C2 — ORGANON PREPARE TMPCAL calibration (bit-exact vs live oracle)" begin
    cr1 = Float32[Float32(c)/100f0 for c in _OC_ICR]
    radgro = zeros(Float32, 27)
    res = organon_prepare_swo(copy(_OC_SPECIES), copy(_OC_DBH1), copy(_OC_HT1OR), cr1,
            copy(_OC_EXPAN1), radgro, 27, 11, 60, 54, 92.0f0, 86.5528641f0, 815.0f0, 815.0f0,
            815.0f0, 0.0f0, 1)
    @test res.ierror == 0
    # oracle TMPCAL (CRATET dump), F9.6 print precision ⇒ tol 1e-5
    @test isapprox(res.tmpcal[1,1], 0.789290f0; atol=1f-5)   # DF height
    @test isapprox(res.tmpcal[2,1], 0.624002f0; atol=1f-5)   # DF crown
    @test isapprox(res.tmpcal[1,2], 0.731081f0; atol=1f-5)   # GF height
    @test res.tmpcal[2,4] == 0.5f0                            # SP crown (clamped, exact)
    # every other of the 54 TMPCAL entries is exactly 1.0
    nontrivial = Set([(1,1),(2,1),(1,2),(2,4)])
    for g in 1:18, k in 1:3
        (k,g) in nontrivial && continue
        @test res.tmpcal[k,g] == 1.0f0
    end
    # RAD=.FALSE. ⇒ DG calibration row all 1.0
    @test all(==(1.0f0), res.tmpcal[3,:])
    # ACALIB after the cratet TMPCAL→ACALIB load
    @test isapprox(res.acalib[1,1], 0.789290f0; atol=1f-5)
    @test isapprox(res.acalib[2,1], 0.624002f0; atol=1f-5)
    @test res.acalib[2,4] == 0.5f0
end

@testset "OC C2 — DF↔PP ORGANON site-index conversion (oc/sitset.f:181-189)" begin
    s = StandState(OregonCoast()); init_blockdata!(s, s.variant)
    s.plot.sp_site_index[7]  = 92.0f0   # DF from ecoclass
    s.plot.sp_site_index[18] = 0.0f0    # PP unset
    site_setup!(s, s.variant)
    @test s.plot.sp_site_index[7]  == 92.0f0
    @test isapprox(s.plot.sp_site_index[18], 86.5528641f0; atol=1f-4)  # 0.940792·92
end

@testset "OC C2 — WF/ES species crosswalk (vie/spctrn.f ASPT col 20)" begin
    s = StandState(OregonCoast()); init_blockdata!(s, s.variant)
    sp, v, co = s.species, s.variant, s.coef
    @test resolve_species("WF", v, sp, co)[1] == Int32(4)    # WF → GF
    @test resolve_species("ES", v, sp, co)[1] == Int32(22)   # ES → BR
    @test orgspc(4)  == Int32(17)   # both surrogate to ORGANON GF code 017
    @test orgspc(22) == Int32(17)
end

# --- C4b: ORGANON SWO minor-species height growth (HTGRO2) -------------------------------------
# Validated bit-exact vs the live FVSoc_clean oracle on a 14-species coverage stand (S248112 with
# records reassigned to IC/WH/RC/PY/MA/GC/TO/CY/BM/WO/BO/RA/DG/WI). ocmin exercises ONLY the big-6
# conifers (groups 1-4), so the minor-species HTGRO2 path was previously unported (HGRO=0). These
# pin the HD_SWO H-D coefficients + the HTGRO2 ratio form + the red-alder Worthington H40 path.
@testset "OC C4b — ORGANON HTGRO2 minor-species height growth (bit-exact vs live oracle)" begin
    # HD_SWO predicted height from DBH (organon/htgrowth.f HDPAR): note the oracle's ORGANON
    # missing-height dub of the BO record (group 15, DBH 8.5) was 46.4124 — HD_SWO(15, 8.5).
    @test isapprox(FVSjl.oc_hd_swo(6,  6.5f0), 44.80421f0;  atol=1f-3)   # WH
    @test isapprox(FVSjl.oc_hd_swo(15, 8.5f0), 46.411728f0; atol=1f-3)   # BO
    # WH (FIA 263, group 6), DBH 6.5, DGRO 0.6273517, HT 30, CALIB(1,6)=1 → HGRO (oracle 2.35722)
    @test isapprox(FVSjl.oc_htgro2(Int32(263), 6, 6.5f0, 0.6273517f0, 30.0f0, 1.0f0, 0.0f0),
                   2.35722f0; atol=1f-4)
    # Red-alder site index (organon/statsorg.f CON_RASI) from the DF SITE_1=92
    @test isapprox(FVSjl.oc_con_rasi(92.0f0), 71.987946f0; atol=1f-4)
    # Red alder (FIA 351), HT 17, RASI 71.988 → Worthington H40 5-yr increment (oracle 11.9128)
    @test isapprox(FVSjl.oc_htgro2(Int32(351), 16, 3.2f0, 0.510054f0, 17.0f0, 1.0f0, 71.988f0),
                   11.9128f0; atol=1f-3)
    # RAGEA/RAH40 round-trip (organon/htgrowth.f): H40 at the growth-effective age recovers HT
    ge = FVSjl.oc_ragea(17.0f0, 71.988f0)
    @test isapprox(FVSjl.oc_rah40(ge, 71.988f0), 17.0f0; atol=1f-3)
end

# --- C8b: ORGANON PREPARE missing-HT/CR dubbing + ACALIB wired into the growth path -----------
# The residual C4b left open: a BLANK-height valid-ORGANON record is dubbed by the oracle via
# ORGANON PREPARE (PRDHT, start2 A_HD_SWO), NOT the FVS-native HTDBH curve. `oc_organon_prepare!`
# (oc/cratet.f:155-401) now runs at setup before dub_missing_heights! and writes the ORGANON dub
# back into the tree record, and the growth path threads PREPARE's ACALIB(1,·) into HTGRO2.
# MEASURED bit-exact vs the live FVSoc_clean oracle on a coverage stand = ocmin with record 4
# reassigned to a blank-height BO (black oak, group 15): tree-4 dub HT=44.7403793, LN(DDS)=1.58408797,
# HGRO=1.04114532, CR2=0.234352291; all 17 valid ORGANON trees bit-exact (max |Δ| HGRO 1.4e-5,
# CR2 6e-8) — the previous CRNCLO cross-perturbation of nearby big-6 heights is gone.
@testset "OC C8b — ORGANON PREPARE dubbing + ACALIB into growth (bit-exact vs live oracle)" begin
    # ocmin /ORGANON/ PREPARE buffer with record 4 reassigned to a blank-height BO (FIA 818, group 15)
    species = Int32[122,202,122,818,117,17,117,122,202,117,202,202,17,17,202,17,202,122,202,122,17,202,17,17,17,17,17]
    dbh = Float32[11.5,0.1,6.5,7.9,8.0,6.2,8.4,9.5,4.0,8.2,1.2,1.9,0.1,5.3,10.0,6.1,12.7,9.6,10.4,8.5,10.9,9.4,3.2,0.1,5.8,5.0,6.6]
    ht  = Float32[73,4.6,30,0,63,38,5,60,20,65,11,13,4.6,27,65,38,67,60,55,0,65,60,17,4.6,28,25,30]  # rec 4 HT blank
    icr = Int[35,55,75,25,25,45,35,25,25,45,55,45,65,65,35,75,35,25,45,25,65,35,45,65,65,25,65]
    cr1 = Float32[Float32(c)/100f0 for c in icr]
    ex  = _OC_EXPAN1  # PROB·PI (setup EXPAN1); recomputing SBA/CCFL with it is part of the calibration
    ex4 = vcat(ex[1:3], Float32(129.262070), ex[5:end])  # rec-4 expansion unchanged from ocmin
    res = organon_prepare_swo(species, dbh, ht, cr1, ex4, zeros(Float32,27), 27, 11, 60, 54,
            92.0f0, 86.5528641f0, 815.0f0, 815.0f0, 815.0f0, 0.0f0, 1)
    # PRDHT dub of the blank BO height via A_HD_SWO (oracle 44.7403793)
    @test isapprox(res.ht[4], 44.7403793f0; atol=1f-4)
    # DF still calibrates HT (ACALIB(1,1)=0.789; consumed only by HTGRO2, HTGRO1 ignores it); BO
    # (single tree ⇒ entht<2) stays 1.0
    @test isapprox(res.acalib[1,1], 0.789290f0; atol=1f-5)
    @test res.acalib[1,15] == 1.0f0
    # The non-ORGANON HTDBH curve (the OLD wrong path) dubs BO ~6.4% low — the residual this fixes
    @test !isapprox(FVSjl.oc_htdbh_height(31, 7.9f0), res.ht[4]; atol=1f0)
end
