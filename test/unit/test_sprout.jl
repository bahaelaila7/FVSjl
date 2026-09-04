# Unit tests for the stump-sprout sub-routines (NSPREC / SPRTHT / ESSPRT, SN).
# Expected values are hand-computed directly from essprt.f's SN SELECT CASE blocks.
using FVSjl: nsprec_sn, sprtht_sn, essprt_sn, sprout_dbh, coefficients, Southern,
             StandState, init_blockdata!, _log_cut!,
             nsprec_ie, essprt_ie, sprtht_ie, nsprec_em, essprt_em, sprtht_em,
             InlandEmpire, EasternMontana

@testset "sprout sub-routines (NSPREC/SPRTHT/ESSPRT)" begin
    coef = coefficients(Southern())

    @testset "NSPREC sprout count" begin
        @test nsprec_sn(5, 6.0f0) == 1          # sp5: DSTMP<7 → 1
        @test nsprec_sn(5, 7.0f0) == 0          # sp5: DSTMP≥7 → 0
        @test nsprec_sn(33, 4.0f0) == 1         # oak grp: <5 → 1
        @test nsprec_sn(33, 5.0f0) == 1         # NINT(-1+0.4·5)=NINT(1.0)=1
        @test nsprec_sn(33, 7.5f0) == 2         # NINT(2.0)=2
        @test nsprec_sn(33, 8.75f0) == 3        # NINT(2.5)=3  (ties away from zero)
        @test nsprec_sn(33, 10.0f0) == 3        # NINT(3.0)=3
        @test nsprec_sn(33, 11.0f0) == 3        # >10 → 3
        @test nsprec_sn(82, 7.5f0) == 2         # same oak/sweetgum group
        @test nsprec_sn(65, 8.0f0) == 1         # non-grouped species → default 1
        @test nsprec_sn(99, 3.0f0) == 1         # out-of-range default
    end

    @testset "SPRTHT sprout height" begin
        @test sprtht_sn(33, 70.0f0, 5) == (0.1f0 + 70f0/50f0) * 5f0   # curve sp
        @test sprtht_sn(59, 70.0f0, 5) == (0.1f0 + 70f0/50f0) * 5f0   # upper range
        @test sprtht_sn(58, 70.0f0, 5) == 0.5f0 + 0.5f0 * 5f0         # gap → default
        @test sprtht_sn(89, 70.0f0, 5) == 0.5f0 + 0.5f0 * 5f0         # >87 → default
        @test sprtht_sn(5,  60.0f0, 3) == (0.1f0 + 60f0/50f0) * 3f0   # explicit sp5
    end

    @testset "ESSPRT survival multiplier" begin
        # constant-multiplier species
        @test essprt_sn(coef, 5,  1.0f0, 6.0f0, 520) == 1.0f0 * 0.42f0
        @test essprt_sn(coef, 33, 1.0f0, 6.0f0, 520) == 1.0f0 * 0.93f0
        @test essprt_sn(coef, 22, 2.0f0, 6.0f0, 520) == 2.0f0 * 0.73f0
        # logistic species (sp20: a=4.1975, b=-0.1821)
        @test essprt_sn(coef, 20, 1.0f0, 8.0f0, 520) ==
              1f0 / (1f0 + exp(-(4.1975f0 - 0.1821f0 * 8f0)))
        # default logistic (sp65 falls through to a=2.7386, b=-0.1076)
        @test essprt_sn(coef, 65, 1.0f0, 8.0f0, 520) ==
              1f0 / (1f0 + exp(-(2.7386f0 - 0.1076f0 * 8f0)))
        # forest-special sp64: special forest 809 uses the cubic-poly form
        @test essprt_sn(coef, 64, 1.0f0, 10.0f0, 809) == (57.3f0 - 0.0032f0*1000f0)/100f0
        # forest-special sp64 in a common forest uses the ELSE logistic (3.8897, -0.2260)
        @test essprt_sn(coef, 64, 1.0f0, 10.0f0, 520) ==
              1f0 / (1f0 + exp(-(3.8897f0 - 0.2260f0 * 10f0)))
        # sp77 special-forest inverse-logistic form
        @test essprt_sn(coef, 77, 1.0f0, 10.0f0, 905) ==
              1f0 / (1f0 + exp(-(-2.8058f0 + 22.6839f0 *
                                  (1f0 / ((10f0/0.7788f0) - 0.4403f0)))))
    end

    @testset "ESUCKR sprout DBH (Wykoff H-D inverse, IABFLG=1)" begin
        # DBH = HT2/(ln(HT-4.5) - HT1) - 1, floored at 0.1; HT<=4.5 -> 0.1.
        # sp65: HT1=4.5142, HT2=-5.2205
        @test sprout_dbh(coef, 65, 15.0f0) == (-5.2205f0 / (log(15f0 - 4.5f0) - 4.5142f0) - 1f0)
        @test sprout_dbh(coef, 33, 12.0f0) == (-4.7206f0 / (log(12f0 - 4.5f0) - 4.4772f0) - 1f0)
        @test sprout_dbh(coef, 22,  4.5f0) == 0.1f0     # HT == 4.5 -> floor
        @test sprout_dbh(coef, 22,  3.0f0) == 0.1f0     # HT < 4.5  -> floor
        @test sprout_dbh(coef, 65, 15.0f0) >= 0.1f0     # never below floor
    end

    @testset "ESTUMP cut-log fidelity (_log_cut!)" begin
        s = StandState(Southern()); init_blockdata!(s, s.variant)
        t = s.trees; t.n = 3
        t.species[1] = Int32(22); t.dbh[1] = 8.0f0;  t.plot_id[1] = Int32(2)  # sprouting
        t.species[2] = Int32(17); t.dbh[2] = 9.0f0;  t.plot_id[2] = Int32(1)  # NOT sprouting
        t.species[3] = Int32(65); t.dbh[3] = 12.0f0; t.plot_id[3] = Int32(3)  # sprouting
        s.control.lsprut = true; s.plot.cycle_length = 10f0
        empty!(s.control.cut_log)
        _log_cut!(s, t, 1, 5.0f0)
        _log_cut!(s, t, 2, 4.0f0)   # non-sprouting species 17 → skipped
        _log_cut!(s, t, 3, 6.0f0)
        _log_cut!(s, t, 1, 0.0f0)   # zero PREM → skipped
        log = s.control.cut_log
        @test length(log) == 2                       # only the two sprouting cuts logged
        @test log[1] == (species = Int32(22), dstmp = 8.0f0,  prem = 5.0f0, plot = Int32(2), ishag = Int32(10))
        @test log[2] == (species = Int32(65), dstmp = 12.0f0, prem = 6.0f0, plot = Int32(3), ishag = Int32(10))
        # lsprut off ⇒ nothing logged at all
        s.control.lsprut = false; empty!(s.control.cut_log)
        _log_cut!(s, t, 1, 5.0f0)
        @test isempty(s.control.cut_log)
    end

    # IE/EM stump & root sprouting (ie/em essprt.f CASE('IE')/CASE('EM')). Coeffs verbatim from essprt.f.
    @testset "IE sprout (NSPREC/ESSPRT/SPRTHT CASE('IE'))" begin
        # NSPREC (essprt.f:221): aspen(18)→2; CO(19) DSTMP-branched; else 1.
        @test nsprec_ie(18, 8.0f0) == 2
        @test nsprec_ie(19, 4.0f0) == 1                       # DSTMP<5 → 1
        @test nsprec_ie(19, 7.5f0) == 2                       # NINT(-1+0.4·7.5)=NINT(2.0)=2
        @test nsprec_ie(19, 10.0f0) == 3                      # ≥10 → 3
        @test nsprec_ie(17, 8.0f0) == 1                       # PY → default 1
        # ESSPRT (essprt.f:70): PY(17)×0.40, CO(19)×0.90, MM/PB(20,21)×0.70, aspen(18)/default ×1.
        @test essprt_ie(17, 1.0f0, 6.0f0) == 0.40f0
        @test essprt_ie(19, 1.0f0, 6.0f0) == 0.90f0
        @test essprt_ie(20, 1.0f0, 6.0f0) == 0.70f0
        @test essprt_ie(21, 1.0f0, 6.0f0) == 0.70f0
        @test essprt_ie(18, 1.0f0, 6.0f0) == 1.0f0            # aspen handled by ASSPTN, not ESSPRT
        # SPRTHT (essprt.f:329): aspen(18) SI/80; PY/CO/MM/PB(17,19,20,21) SI/100; else 0.5+0.5·IAG.
        @test sprtht_ie(18, 48.0f0, 10) == (0.1f0 + 48f0/80f0) * 10f0
        @test sprtht_ie(21, 50.0f0, 10) == (0.1f0 + 50f0/100f0) * 10f0
        @test sprtht_ie(9,  50.0f0, 10) == 0.5f0 + 0.5f0 * 10f0
    end

    @testset "EM sprout (NSPREC/ESSPRT/SPRTHT CASE('EM'))" begin
        # NSPREC (essprt.f:205): aspen(12)→2; CW/BA/PW/NC(13:16) DSTMP-branched; else 1.
        @test nsprec_em(12, 8.0f0) == 2
        @test nsprec_em(14, 4.0f0) == 1
        @test nsprec_em(14, 7.5f0) == 2
        @test nsprec_em(16, 10.0f0) == 3
        @test nsprec_em(11, 8.0f0) == 1                       # GA → default 1
        # ESSPRT (essprt.f:92): GA(11) DSTMP-split 0.80/0.50; CW/PW(13,15)×0.90; BA(14) 0.80/0.50; NC(16)×0.80; PB(17)×0.70.
        @test essprt_em(11, 1.0f0, 12.0f0) == 0.80f0
        @test essprt_em(11, 1.0f0, 13.0f0) == 0.50f0
        @test essprt_em(13, 1.0f0, 6.0f0) == 0.90f0
        @test essprt_em(14, 1.0f0, 25.0f0) == 0.80f0
        @test essprt_em(14, 1.0f0, 26.0f0) == 0.50f0
        @test essprt_em(16, 1.0f0, 6.0f0) == 0.80f0
        @test essprt_em(17, 1.0f0, 6.0f0) == 0.70f0
        @test essprt_em(12, 1.0f0, 6.0f0) == 1.0f0            # aspen via ASSPTN
        # SPRTHT (essprt.f:316): AS/PB(12,17) SI/80; GA/CW/BA/PW/NC(11,13:16) SI/100; else default.
        @test sprtht_em(12, 48.0f0, 10) == (0.1f0 + 48f0/80f0) * 10f0
        @test sprtht_em(11, 50.0f0, 10) == (0.1f0 + 50f0/100f0) * 10f0
        @test sprtht_em(5,  50.0f0, 10) == 0.5f0 + 0.5f0 * 10f0
    end

    # ESTUMP cut-log marks IE {17,18,19,20,21} and EM {11..17} as sprouters (is_sprouting CSV column).
    @testset "IE/EM ESTUMP cut-log sprouter set" begin
        for (variant, sprouters, nonsp) in ((InlandEmpire(), (17,18,19,20,21), (1,10,15)),
                                            (EasternMontana(), (11,12,13,14,15,16,17), (1,10,18)))
            s = StandState(variant); init_blockdata!(s, s.variant)
            t = s.trees; s.control.lsprut = true; s.plot.cycle_length = 10f0
            allsp = (sprouters..., nonsp...)
            t.n = length(allsp)
            for (i, sp) in enumerate(allsp); t.species[i] = Int32(sp); t.dbh[i] = 8.0f0; t.plot_id[i] = Int32(1); end
            empty!(s.control.cut_log)
            for i in 1:t.n; _log_cut!(s, t, i, 3.0f0); end
            logged = Set(Int(r.species) for r in s.control.cut_log)
            @test logged == Set(sprouters)          # exactly the ISPSPE set is logged
        end
    end
end
