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
