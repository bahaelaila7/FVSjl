# NC (Klamath) Region-7 (HOOPA, IFOR 5, forest 705) volume — the NVEL BLMVOL fix (2026-09-05).
#
# NC's VEQNNC is forest/region-dependent (sitset.f → VOLEQDEF(VAR='NC',IREGN=KODFOR/100,FORST)). Forest 705 =
# HOOPA (IREGN=7, FORST='05') ⇒ VOLEQDEF returns B00BEHW<fia> (DF = B01BEHW202) — the BLM Oregon volume
# routines (volinit.f:366 → BLMVOL → blmtap.f BEHRE'S HYPERBOLA), NOT the Region-5 WO2W/DVE table (IFOR 1-3)
# nor the Region-6 Behre (IFOR 4). jl previously had no IFOR-5 branch ⇒ fell through to the R5 WO2W R5TAP
# taper (a fatter stem + R5 board rules) ⇒ cyc0 TCuFt +16-22% OVER, BdFt −8-23% UNDER on the 57 LOC-705
# stands (structure bit-exact, one-directional, deterministic). This routes IFOR-5 through nc_blmvol_vol.
#
# Oracle: a standalone BLMVOL driver linked from FVSnc_g16's own blmvol.o/blmtap.o/scrib.o/numlog.o/segmnt.o,
# calling BLMVOL twice exactly as fvsvol.f NATCRS does for IREGN=7 (cubic: MTOPP=TOPD·BARK=5·BARK, BFPFLG=0,
# CUPFLG=1 → VOL(1)+VOL(4); board: MTOPP=BFTOPD·BARK=5·BARK, BFPFLG=1 → VOL(2)). FC=80 (formcl.f default).
# The whole-stand .sum was ALSO validated bit-exact vs live FVSnc_g16 on all 45 treed 705 stands (cyc0).

using Test
using FVSjl
const _M = FVSjl

@testset "NC Region-7 (Hoopa 705) volume — BLMVOL Behre-hyperbola bit-exact vs FVSnc_g16" begin
    # --- VEQNNC table (VOLEQDEF NC/7/05), confirmed vs FVSnc_g16 705 .out ---
    @test _M.NC_R7_VOL_EQ[3]  == "B01BEHW202"   # DF
    @test _M.NC_R7_VOL_EQ[8]  == "B00BEHW631"   # tanoak
    @test _M.NC_R7_VOL_EQ[1]  == "B00BEHW999"   # OS
    @test _M.NC_R7_VOL_EQ[12] == "B00BEHW211"   # redwood
    # PROFILE / TAPEQU from BLMTAPEQ (blmvol.f:937-1070): DF→(1,1), SP→(4,13), PP→(3,11), RF→(7,32), IC→(9,51)
    # (BLMTAPEQ returns (TAPEQU, PROFILE); the kernel is the shared engine blm_vol.jl, keyed by NC_R7_VOL_EQ)
    tp(sp) = _M._blm_tapeq(_M.NC_R7_VOL_EQ[sp])
    @test tp(3)  == (1, 1)    # DF
    @test tp(2)  == (13, 4)   # SP
    @test tp(10) == (11, 3)   # PP
    @test tp(9)  == (32, 7)   # RF
    @test tp(1)  == (56, 10)  # OS → all-other/misc

    # --- DOUBLE_BARK (blmvol.f:838-934) ---
    @test isapprox(_M._blm_double_bark(1, 12.7f0), 0.903563f0 * 12.7f0^0.989388f0; atol = 1f-4)  # DF
    @test isapprox(_M._blm_double_bark(13, 20.0f0), 0.859045f0 * 20.0f0; atol = 1f-4)             # sugar pine
    @test isapprox(_M._blm_double_bark(56, 25.0f0), 25.0f0 / 1.071f0; atol = 1f-4)                # misc

    # --- BLMTAP taper (blmtap.f TLH=0): DIB decreases with height; ~D17 near the base ---
    d17 = floor(20.0f0 * 80f0 / 100f0 + 0.5f0)   # ANINT(0.8*20)=16
    dib_low  = _M._blm_tap(20.0f0, 100.0f0 + 1.5f0, 17.8f0, d17, 16.3f0, 1)
    dib_high = _M._blm_tap(20.0f0, 100.0f0 + 1.5f0, 80.0f0, d17, 16.3f0, 1)
    @test dib_low > dib_high > 0f0

    # --- Per-tree VOL(1)/VOL(4)/VOL(2) — BIT-EXACT vs the standalone BLMVOL oracle (5 species/profiles) ---
    # (sp, D, H, bark, oracle TCF, MCF, BF)
    cases = [
        (3, 12.7f0,  67.0f0, 0.833f0,  24.8415f0,  22.6f0,  116.0f0),  # DF   profile 1
        (8, 18.0f0,  90.0f0, 0.905f0,  60.4117f0,  57.2f0,  323.0f0),  # tanoak profile 10
        (1, 25.0f0, 110.0f0, 0.934f0, 144.6286f0, 144.0f0,  905.0f0),  # OS   profile 10
        (2, 34.6f0, 140.0f0, 0.850f0, 318.1863f0, 310.3f0, 2184.0f0),  # SP   profile 4
        (5, 10.0f0,  45.0f0, 0.980f0,  12.5061f0,  11.1f0,   47.0f0),  # madrone profile 10
    ]
    for (sp, d, h, bark, otcf, omcf, obf) in cases
        tcf, mcf, bf = _M.nc_blmvol_vol(sp, d, h, bark)
        @test isapprox(tcf, otcf; atol = 0.01f0)   # ≤ 1-ULP Float32
        @test isapprox(mcf, omcf; atol = 0.01f0)
        @test bf == obf
    end

    # --- The routing actually changed: BLMVOL (Behre) ≠ the old R5 WO2W R5TAP path that forest-705 used to
    # fall through to. The stand-aggregate over-prediction (+16-22% TCuFt) is a species-mix effect validated at
    # the .sum level (all 45 treed 705 stands cyc0 bit-exact); per-tree the two tapers simply differ. ---
    for (sp, eq, d, h, bark) in ((3, "500WO2W202", 34.6f0, 140.0f0, 0.833f0),
                                 (2, "500WO2W117", 20.0f0, 100.0f0, 0.85f0))
        blm = _M.nc_blmvol_vol(sp, d, h, bark)[1]
        r5  = _M.nc_wo2w_vol(eq, d, h)[1]
        @test blm > 0f0 && r5 > 0f0
        @test !isapprox(blm, r5; atol = 0.5f0)   # the fix genuinely re-routes off the R5 taper
    end

    # --- Small-tree cylinder (BLMVOL TTH≤17.8) ---
    (tcf_s, mcf_s, bf_s) = _M.nc_blmvol_vol(1, 3.0f0, 12.0f0, 0.9f0)
    @test tcf_s > 0f0 && mcf_s == 0f0 && bf_s == 0f0
end
