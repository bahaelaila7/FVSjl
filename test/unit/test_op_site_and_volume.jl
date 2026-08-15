# =============================================================================
# test_op_site_and_volume.jl — OP (Olympic) chunk-2 A/B vs the live FVSop_clean oracle:
#   (1) the ecoclass SITE-INDEX FAN (op/forkod+habtyp+ecocls+sichg+htcalc+sitset), and
#   (2) the BLM Behre-taper cubic/merch/board VOLUME (op/voleqdef+formcl BLM708 + blmvol/blmtap).
#
# Reference stand S248112 (opdbg.key): forest 708 (BLM Salem, IFOR=4), habitat 40 → plant
# association CHS133 → DF(16) site species, SITEAR(16)=98.
#
# (1) SITEAR fan validated vs the FVSop_clean "SITECODE" dump (the 18 printed species):
#     WF=97.98189, ES=139.15092, LP=98.22826, SP=139.15092, PP=139.15092, DF=98.00, WH=87.67.
# (2) Volume validated PER-TREE vs the live oracle, INSTRUMENTED: a WRITE was added to a temp copy
#     of fvsvol.f (buildDir kept pristine, marker count 0 afterward), FVSop_dbg relinked, and the
#     cyc0 per-tree D/HT/TCF/MCF/BF captured (fort.16). All 27 tree records reproduce the oracle
#     TCF (total cuft), MCF (merch cuft, D≥7), and Scribner BF (D≥7) to the F10.3 print. The two
#     dubbed heights (66.31 DF, 62.39 LP) come from the CRATET/ORGANON height-dub (chunk 3).
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "OP site-index fan (S248112, forest 708 / habitat 40)" begin
    s = F.StandState(F.Olympic())
    p = s.plot
    p.user_forest_code = Int32(708)
    p.habitat_code = Int32(40)
    p.site_species = Int32(0)
    fill!(p.sp_site_index, 0f0); fill!(p.sp_sdi_def, 0f0)
    F.op_site_index_setup!(s)
    si = p.sp_site_index

    @test Int(p.forest_idx) == 4          # forkod: 708 → IFOR 4 (BLM Salem)
    @test Int(p.site_species) == 16       # ecocls CHS133 site species = DF

    # SITEAR fan bit-exact vs the oracle SITECODE dump (to the King/Cochran curve full precision).
    @test si[2]  ≈ 97.98189f0  atol=1f-3   # WF
    @test si[10] ≈ 139.15092f0 atol=1f-3   # ES
    @test si[11] ≈ 98.22826f0  atol=1f-3   # LP
    @test si[13] ≈ 139.15092f0 atol=1f-3   # SP
    @test si[15] ≈ 139.15092f0 atol=1f-3   # PP
    @test si[16] ≈ 98.0f0      atol=1f-3   # DF (site species)
    @test si[19] ≈ 87.67f0     atol=1f-2   # WH (DF→WH Nigh conversion)
    # whole 18-species rounded fan matches the dump (SF..RC)
    @test [round(Int, si[i]) for i in 1:18] ==
          [139,98,98,139,98,98,139,139,139,139,98,139,139,139,139,98,139,98]
    @test Int(round(p.sp_sdi_def[16])) == 950   # SDIDEF(DF) = min(1606, FORMAX 950)
end

@testset "OP BLM volume per-tree (S248112 cyc0, vs instrumented FVSop_dbg)" begin
    # (ISPC, D, HT, oracle TCF, oracle MCF, oracle BF) — 27 cyc0 records from fort.16 (fvsvol dump).
    oracle = [
        (11,11.500f0,73.00f0,16.879f0,16.200f0,73f0),(16,0.100f0,2.00f0,0f0,0f0,0f0),
        (15,6.500f0,30.00f0,2.533f0,0f0,0f0),(13,7.900f0,75.00f0,9.179f0,7.500f0,44f0),
        (13,8.000f0,63.00f0,8.254f0,6.400f0,37f0),(2,6.200f0,38.00f0,3.391f0,0f0,0f0),
        (13,8.400f0,5.00f0,0.923f0,0f0,0f0),(11,9.500f0,60.00f0,8.309f0,7.800f0,37f0),
        (16,4.000f0,20.00f0,0.744f0,0f0,0f0),(13,8.200f0,65.00f0,8.403f0,6.400f0,37f0),
        (16,1.200f0,11.00f0,0.040f0,0f0,0f0),(16,1.900f0,13.00f0,0.115f0,0f0,0f0),
        (2,0.100f0,3.00f0,0f0,0f0,0f0),(2,5.300f0,27.00f0,1.788f0,0f0,0f0),
        (16,10.000f0,65.00f0,15.601f0,14.400f0,73f0),(2,6.100f0,38.00f0,3.283f0,0f0,0f0),
        (16,12.700f0,67.00f0,24.841f0,22.600f0,116f0),(11,9.600f0,60.00f0,11.310f0,11.100f0,57f0),
        (16,10.400f0,66.31f0,15.810f0,14.400f0,73f0),(11,8.500f0,62.39f0,8.522f0,7.100f0,37f0),
        (2,10.900f0,65.00f0,19.707f0,18.100f0,86f0),(16,9.400f0,60.00f0,14.802f0,13.000f0,65f0),
        (10,3.200f0,17.00f0,0.409f0,0f0,0f0),(10,0.100f0,2.00f0,0f0,0f0,0f0),
        (10,5.800f0,28.00f0,2.286f0,0f0,0f0),(10,5.000f0,25.00f0,1.507f0,0f0,0f0),
        (2,6.600f0,30.00f0,3.065f0,0f0,0f0)]
    for (sp,d,h,otcf,omcf,obf) in oracle
        tcf = F.op_tree_cuft(sp, d, h)
        v4, v2 = F.op_tree_mvol(sp, d, h)
        mcf = d >= 7f0 ? v4 : 0f0
        bf  = d >= 7f0 ? v2 : 0f0
        @test tcf ≈ otcf atol=0.05     # F10.3 total cubic print
        @test mcf ≈ omcf atol=0.05     # merch cubic (D≥DBHMIN 7)
        @test bf  ≈ obf  atol=0.5      # Scribner board (D≥BFMIND 7)
    end
end
