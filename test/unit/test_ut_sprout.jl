# =============================================================================
# test_ut_sprout.jl — Utah stump/root-sprout subsystem (strp/esuckr.f +
# vstrp/essprt.f CASE('UT'), Crouch aspen ASSPTN).
#
# Guards the #155 UT sprout dispatch (the esuckr! coordinator used to lack a UT
# path and crashed / mis-dispatched into the SN essprt model). Expected values
# are hand-computed straight from the Fortran CASE('UT') SELECT blocks:
#   NSPREC  vstrp/essprt.f:1168  — CASE(6)=2; CASE(13,22)=NINT(0.2·D); CASE(18,19)=NINT(-1+0.4·D)
#   ESSPRT  vstrp/essprt.f:627   — CASE(22)=0.50; CASE(21)=0.70; CASE(13,18,19)=0.80; ELSE 1
#   SPRTHT  vstrp/essprt.f:1429  — CASE(13,18,19,21,22)=(0.1+SI/100)·AGE; CASE(6)=(0.1+SI/80)·AGE
#   ISPSPE  blkdat.f:115         — the six UT sprouters {6,13,18,19,21,22}
# These are NON-VACUOUS: deleting the UT branch in sprout.jl / mis-setting the
# is_sprouting flags flips these assertions.
# =============================================================================
using Test
using FVSjl: nsprec_ut, essprt_ut, sprtht_ut, esuckr!, coefficients, coef_col,
             Utah, StandState, init_blockdata!, _log_cut!

@testset "UT stump-sprout subsystem (CASE('UT'))" begin

    @testset "NSPREC UT sprout count (vstrp/essprt.f:1168)" begin
        @test nsprec_ut(6, 8.0f0)  == 2          # aspen → always 2
        @test nsprec_ut(6, 20.0f0) == 2
        # sp13/22: <5→1; 5..10→NINT(0.2·D); >10→2
        @test nsprec_ut(13, 4.0f0)  == 1
        @test nsprec_ut(13, 5.0f0)  == 1         # NINT(1.0)=1
        @test nsprec_ut(13, 8.0f0)  == 2         # NINT(1.6)=2
        @test nsprec_ut(22, 8.0f0)  == 2
        @test nsprec_ut(13, 11.0f0) == 2         # >10 → 2
        # sp18/19: <5→1; 5..10→NINT(-1+0.4·D); >10→3
        @test nsprec_ut(18, 4.0f0)  == 1
        @test nsprec_ut(18, 7.5f0)  == 2         # NINT(2.0)=2
        @test nsprec_ut(19, 8.75f0) == 3         # NINT(2.5)=3 (ties away from zero)
        @test nsprec_ut(18, 10.0f0) == 3         # NINT(3.0)=3
        @test nsprec_ut(18, 11.0f0) == 3         # >10 → 3
        @test nsprec_ut(21, 8.0f0)  == 1         # sp21 not in NSPREC CASE → default 1
        @test nsprec_ut(3,  8.0f0)  == 1         # non-sprouter default 1
    end

    @testset "ESSPRT UT survival multiplier (vstrp/essprt.f:627)" begin
        @test essprt_ut(22, 1.0f0, 8.0f0) == 0.50f0
        @test essprt_ut(21, 1.0f0, 8.0f0) == 0.70f0
        @test essprt_ut(13, 2.0f0, 8.0f0) == 2.0f0 * 0.80f0
        @test essprt_ut(18, 1.0f0, 8.0f0) == 0.80f0
        @test essprt_ut(19, 1.0f0, 8.0f0) == 0.80f0
        @test essprt_ut(6,  3.0f0, 8.0f0) == 3.0f0   # aspen → ELSE ×1 (TPA set by ASSPTN, not ESSPRT)
    end

    @testset "SPRTHT UT sprout height (vstrp/essprt.f:1429)" begin
        @test sprtht_ut(13, 70.0f0, 5) == (0.1f0 + 70f0/100f0) * 5f0
        @test sprtht_ut(21, 70.0f0, 5) == (0.1f0 + 70f0/100f0) * 5f0
        @test sprtht_ut(22, 70.0f0, 5) == (0.1f0 + 70f0/100f0) * 5f0
        @test sprtht_ut(6,  80.0f0, 3) == (0.1f0 + 80f0/80f0)  * 3f0   # aspen uses /80
        @test sprtht_ut(3,  70.0f0, 5) == 0.5f0 + 0.5f0 * 5f0          # non-sprouter default
    end

    @testset "UT sprouter gate = ISPSPE {6,13,18,19,21,22} (blkdat.f:115)" begin
        coef = coefficients(Utah())
        spr = coef_col(coef, :is_sprouting)
        for sp in (6, 13, 18, 19, 21, 22)
            @test spr[sp] == 1f0
        end
        # sp20 (MC) is NOT a UT sprouter — was wrongly flagged before the fix.
        @test spr[20] == 0f0
        for sp in (1, 3, 7, 10, 14, 20, 23, 24)
            @test spr[sp] == 0f0
        end
    end

    @testset "ESUCKR UT dispatch — aspen cut sprouts, no crash (#155)" begin
        # Build a UT stand, cut an aspen (sp6), and drive the sprout coordinator.
        # Before the #155 fix this path crashed (esuckr!→essprt_sn dispatch gap).
        s = StandState(Utah()); init_blockdata!(s, s.variant)
        t = s.trees; t.n = 2
        t.species[1] = Int32(6);  t.dbh[1] = 8.0f0; t.plot_id[1] = Int32(1)  # aspen (sprouter)
        t.species[2] = Int32(3);  t.dbh[2] = 9.0f0; t.plot_id[2] = Int32(1)  # DF (not a sprouter)
        s.plot.sp_site_index[6] = 60f0                                        # SITEAR(aspen)
        s.control.lsprut = true; s.plot.cycle_length = 10f0
        empty!(s.control.cut_log)
        _log_cut!(s, t, 1, 30.0f0)   # cut 30 TPA of aspen
        _log_cut!(s, t, 2, 10.0f0)   # DF cut → NOT logged (non-sprouter)
        @test length(s.control.cut_log) == 1          # only the aspen cut logged
        n_before = t.n
        created = esuckr!(s; fint = 10f0)             # must NOT throw
        @test created
        @test t.n > n_before                          # sprout records added
        # aspen → NSPREC=2 records, each IMC=2 (sprout regeneration)
        newrecs = (n_before + 1):t.n
        @test length(newrecs) == 2
        for n in newrecs
            @test t.mort_code[n] == Int32(2)          # IMC=2
            @test t.species[n]   == Int32(6)          # aspen
            @test t.tpa[n]        > 0f0               # Crouch ASSPTN sucker TPA
            @test t.birth_age[n] == 10f0              # ABIRTH = ISHAG
        end
    end
end
