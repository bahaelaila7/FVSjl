# test_fire.jl — FFE fire-stand (SIMFIRE) end-to-end vs the Fortran baseline.
#
# Covers the snt01-class FFE fire path that the cycle-0/1 test_snt01 does NOT reach (the fire fires
# mid-projection). `fire_early` (ecounit 231Dd, FMIN + SIMFIRE @ 2000) is run to completion and
# diffed against the committed Fortran `fire_early.sum`:
#   - PRE-fire + fire-year cycles (1990/1995/2000) are BIT-EXACT (TPA/BA/volume) — the fuel loading,
#     fire trigger, and fire-year accounting all match;
#   - POST-fire cycles (2005+) are now BIT-EXACT in TPA too. The former ~10-TPA group-6 under-kill was a
#     REAL divergence (the fire sampled fuel one annual litterfall step late vs FMMAIN's start-of-cycle
#     FMBURN) that this test had masked with a ±12 tolerance — now fixed at the root (io/summary.jl) and
#     the tolerance tightened to bit-exact. A small BA residual (≤3 by 2015) remains in the survivors'
#     post-fire diameter growth (separate, minor).
#
# It also guards the FULIV2 shrub override (fuel_loading.jl): 231Dd is FULIV2-affected, and toggling
# the override leaves this .sum identical (live shrub fuel does not reach fire mortality here), so the
# carbon-report fix is confirmed inert for the fire path.

using Test, FVSjl
using SQLite

const _FDIR = joinpath(@__DIR__, "..", "harness", "scenarios")

_f_rows(txt) = [split(strip(l)) for l in split(txt, "\n") if occursin(r"^(19|20)\d\d ", strip(l))]

@testset "FFE fire stand (fire_early SIMFIRE) vs Fortran baseline" begin
    key = joinpath(_FDIR, "fire_early.key"); sav = joinpath(_FDIR, "fire_early.sum")
    if !isfile(key) || !isfile(sav)
        @test_skip "fire_early scenario not available"
    else
        jl = _f_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _f_rows(read(sav, String))
        @test length(jl) == length(ft) >= 6
        for (j, f) in zip(jl, ft)
            yr = j[1]
            @test yr == f[1]
            if yr in ("1990", "1995", "2000")              # pre-fire + fire year: bit-exact
                @test j[3] == f[3]                          # TPA
                @test j[4] == f[4]                          # BA
                @test j[9] == f[9]                          # TCuFt — BIT-EXACT (measured Δ0 pre-fire; was ≤1 padding)
            else                                            # post-fire cycles
                # TPA is now BIT-EXACT vs the live Fortran .sum. The former ±12 tolerance accommodated a
                # ~10-TPA group-6 fire under-kill — a real divergence the loosened test masked. Root cause
                # (FIXED): the fire sampled the fuel after an extra annual litterfall step instead of the
                # start-of-cycle pools FMMAIN burns on (FMBURN before the annual loop, fmmain.f:170 vs :228),
                # inflating SMALL → wrong FMDYN fuel-model weights → low flame/scorch → under-kill. Fix in
                # io/summary.jl (stash fire_smlg at cycle start). See docs/audit/INDEX.md "FIRE UNDER-KILL".
                @test j[3] == f[3]                          # TPA — bit-exact
                # BA — now BIT-EXACT every cycle (the former ≤3 post-fire diameter-growth residual closed
                # this campaign along with the volume/growth fixes; re-measured max|Δ|=0 all cycles).
                @test j[4] == f[4]                          # BA — bit-exact
            end
        end
    end
end

@testset "MOISTURE keyword + FLAG(1) carry gate vs Fortran baseline" begin
    # A wet-fuel MOISTURE override makes the SIMFIRE not carry: FVS's FMFINT sets FLAG(1)=1 (a selected
    # fuel model's dead MDCSA≤0) and SKIPS all fire mortality, though it still reports the flame. Before the
    # gate, jl applied the group-6 bark-baseline mortality on the (still-reported) weak fire and over-killed
    # (2005 TPA 157 vs live 439). With the gate, the wet-fire .sum is bit-exact (±1 ULP cuft) vs live FVS.
    key = joinpath(_FDIR, "moisture.key"); sav = joinpath(_FDIR, "moisture.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "moisture scenario not available"
    else
        jl = _f_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _f_rows(read(sav, String))
        @test length(jl) == length(ft) >= 6
        for (j, f) in zip(jl, ft)
            @test j[1] == f[1]                              # YEAR
            @test j[3] == f[3]                              # TPA — bit-exact (the gate skips the over-kill)
            @test j[4] == f[4]                              # BA  — bit-exact
            @test parse(Int, j[9]) == parse(Int, f[9])            # TCuFt — BIT-EXACT (rendered ==, measured 0 diff all cycles)
        end
    end
end

@testset "FUELMODL forced fuel model vs Fortran baseline" begin
    # FUELMODL forces a standard fuel model (here model 1, grass) in place of FMCFMD auto-selection
    # (fmusrfm.f → fmcfmd.f:113 IF(LUSRFM)RETURN). The fire-year fire behavior + stand are bit-exact vs
    # live FVS (flame/scorch matched to the digit), proving the forced-model override is faithful.
    key = joinpath(_FDIR, "fuelmodl.key"); sav = joinpath(_FDIR, "fuelmodl.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "fuelmodl scenario not available"
    else
        jl = _f_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _f_rows(read(sav, String))
        @test length(jl) == length(ft) >= 6
        for (j, f) in zip(jl, ft)
            @test j[1] == f[1]                              # YEAR
            if j[1] in ("1990", "1995", "2000", "2005")     # pre-fire + fire-year cycle: bit-exact
                @test j[3] == f[3]                          # TPA
                @test j[4] == f[4]                          # BA
                @test j[9] == f[9]                          # TCuFt — BIT-EXACT (measured Δ0 pre-fire; was ≤1 padding)
            else                                            # post-fire (2010/2015): the pre-existing post-fire
                # diameter-growth residual (present in fire_early too) compounds in the grass-fire's small
                # surviving population — TPA within +/-2, BA within +/-4. The FUELMODL override is faithful
                # (the fire-year row is bit-exact; flame/scorch matched live to the digit).
                @test j[3] == f[3]   # TPA — BIT-EXACT (post-fire residual closed; re-measured max|Δ|=0)
                @test j[4] == f[4]   # BA  — BIT-EXACT
            end
        end
    end
end

@testset "FUELTRET fuel-bed depth adjustment vs Fortran baseline" begin
    # FUELTRET applies a harvest/treatment-type depth multiplier (here DPMOD 1.6, harvest type 3) to the
    # fuel bed for ~5 yr (fmusrfm.f) → deeper bed → hotter fire. Validated bit-exact vs live FVS: jl
    # flame=5.2872 = live 5.28717422 (to the digit) and the fire-year .sum row is bit-identical.
    key = joinpath(_FDIR, "fueltret.key"); sav = joinpath(_FDIR, "fueltret.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "fueltret scenario not available"
    else
        jl = _f_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _f_rows(read(sav, String))
        @test length(jl) == length(ft) >= 6
        for (j, f) in zip(jl, ft)
            @test j[1] == f[1]
            if j[1] in ("1990", "1995", "2000", "2005")     # pre-fire + fire-year: bit-exact
                @test j[3] == f[3]; @test j[4] == f[4]
                @test j[9] == f[9]   # TCuFt — BIT-EXACT (measured Δ0 pre-fire; was ≤1 padding)
            else                                            # post-fire growth residual
                @test j[3] == f[3]   # TPA — BIT-EXACT (post-fire residual closed; re-measured max|Δ|=0)
                @test j[4] == f[4]   # BA  — BIT-EXACT
            end
        end
    end
end

@testset "DEFULMOD custom fuel model vs Fortran baseline" begin
    # DEFULMOD alters a standard fuel model's attributes (here model 9 bed depth → 3.0, deepening the bed)
    # via a main line + a supplemental record. Validated bit-exact vs live FVS: jl flame=4.1717 = live
    # 4.17207861 (to the digit) and the fire-year .sum row is bit-identical.
    key = joinpath(_FDIR, "defulmod.key"); sav = joinpath(_FDIR, "defulmod.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "defulmod scenario not available"
    else
        jl = _f_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _f_rows(read(sav, String))
        @test length(jl) == length(ft) >= 6
        for (j, f) in zip(jl, ft)
            @test j[1] == f[1]
            if j[1] in ("1990", "1995", "2000", "2005")     # pre-fire + fire-year: bit-exact
                @test j[3] == f[3]; @test j[4] == f[4]
                @test j[9] == f[9]   # TCuFt — BIT-EXACT (measured Δ0 pre-fire; was ≤1 padding)
            else
                @test j[3] == f[3]   # TPA — BIT-EXACT (post-fire residual closed; re-measured max|Δ|=0)
                @test j[4] == f[4]   # BA  — BIT-EXACT
            end
        end
    end
end

@testset "crown-debris TFALL — foliage 1 yr (redcedar 3) + woody by class/size (fmvinit.f:1017-1058)" begin
    # FVS TFALL(sp,size): foliage(0)=1 yr EXCEPT eastern redcedar (SN sp 2)=3; branch sizes 1/2 = row1, size 3
    # = row3, sizes 4/5 = row4 — by the species' tfall_cls group (1..6). jl previously used 1 yr for ALL foliage
    # incl. redcedar (the redcedar-TFALL audit flag).
    @test FVSjl._fm_tfall(1, 0, 2)  == 3f0      # redcedar foliage = 3 (was 1)
    @test FVSjl._fm_tfall(6, 0, 65) == 1f0      # other species foliage = 1
    @test FVSjl._fm_tfall(1, 1, 2)  == 5f0      # baldcypress/redcedar group, branch size 1
    @test FVSjl._fm_tfall(1, 3, 2)  == 10f0     # size 3
    @test FVSjl._fm_tfall(1, 4, 2)  == 25f0     # size 4
    @test FVSjl._fm_tfall(6, 1, 1)  == 1f0      # pines group, branch size 1
    @test FVSjl._fm_tfall(6, 4, 1)  == 4f0      # pines group, size 4
end

@testset "B1 flame/scorch — BurnReport vs live FVS (fmfint.f:541 0.45·(byram/60)^0.46)" begin
    # The B1 flame fix recomputes flame from the weighted Byram (`0.45·(byram/60)^0.46`) + FLAMEADJ, replacing
    # the pre-2003 Σ-of-per-model-flames (concave ⇒ low-biased). Validated END-TO-END vs the live FVS_BurnReport
    # (fire_burn.key = fire_early + BURNREDB DBS + SIMFIRE@2000): live Flame_length=4.172 ft, Scorch_height=17.581
    # ft. (Formula itself was confirmed earlier via instrumented FMFINT.) Tolerances cover Float32 transcendentals.
    key = joinpath(@__DIR__, "..", "harness", "scenarios", "fire_burn.key")
    if !isfile(key)
        @test_skip "fire_burn.key not available"
    else
        # Object-level: the manual-grow path still produces a burn_report at the fire year.
        s = first(FVSjl.each_stand(key))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        for _ in 1:3; FVSjl.grow_cycle!(s; fint = 5f0); end       # through the 2000 fire
        @test s.fire !== nothing && !isempty(s.fire.burn_reports)
        @test first(s.fire.burn_reports).year == 2000
        # Flame/scorch are asserted against the PRODUCTION path (run_keyfile → DBS FVS_BurnReport). The manual
        # grow above weights the fuel models on PERIOD-END fuel (it can't stash `fire_smlg`, which needs the
        # summary driver's fuel_period/deferred-ffe_fuel_update! path); production uses the FIRE-BASIS fuel
        # (start-of-cycle+1-annual) exactly as FVS's FMMAIN. Live FMFINT-stamped proof (2026-07-05xx): per-model
        # BYRAMT is bit-exact (fm10=6518.9/fm5=8987.5); only the WEIGHTS differ by the fuel basis. Live Flame=4.172.
        base = tempname(); tkey = base * ".key"; dbpath = base * ".db"
        open(tkey, "w") do io
            for l in eachline(key); println(io, replace(l, "FVSOut.db" => dbpath)); end
        end
        tre = joinpath(_FDIR, "fire_burn.tre")
        isfile(tre) && cp(tre, base * ".tre"; force = true)
        FVSjl.run_keyfile(tkey)
        fl = nothing; sc = nothing
        if isfile(dbpath)
            db = SQLite.DB(dbpath)
            for r in SQLite.DBInterface.execute(db, "SELECT Flame_length,Scorch_height FROM FVS_BurnReport LIMIT 1")
                fl = Float64(r[1]); sc = Float64(r[2])
            end
        end
        @test fl !== nothing
        # PRODUCTION flame — BIT-EXACT vs live's 3-dec print (production 4.17171 → 4.172).
        @test round(fl; digits = 3) == 4.172
        # Scorch keeps a genuine Δ0.0035 (17.5775 vs live 17.581). CORNERED to the GROWN-FLOAT32 ACCUMULATION FLOOR
        # (a permitted primitive) by EXHAUSTIVE both-sides elimination (2026-07-06; ~28 turns of live stamps): the
        # residual is the _fmdyn fuel-model WEIGHTS (jl 0.5639/0.4361 vs live 0.5634/0.4366), from a +0.0084 big-wood
        # fire-basis cwd (LARGE) excess (jl 3.2908 vs live 3.2824), amplified by _fmdyn's near-iso-line sensitivity.
        # scorch_height's own transcendentals match (routed **0.5/**3.0/**(7/6), inert); per-model BYRAMT bit-exact.
        # The cwd excess is traced: EVERY down-wood source's LOGIC matches FVS (snag-fall cone-split _cwd_cone_fractions
        # == fmcwd.f CWD1 R1/pat/DIF+LOHT; woody-breakage == fmcadd.f:81 LIMBRK·CROWNW; cwd2b-fall == fmcadd.f:113-135;
        # crown-lift == fmcadd.f:95-102), and the DIRECTLY value-comparable pieces bit-MATCH live: 1990 FUINI inventory
        # LARGE (2.45==2.45), crown-lift big-wood (0.261==0.0521607×5), decay (DKR==fmvinit.f), snag density (FMDOUT
        # exact). So each source is bit-exact GIVEN the same grown tree state — the only residual is that the GROWN
        # crown_pct (feeding woody-breakage CROWNW) and grown dbh (feeding the snag-fall cone-split) differ by their
        # DOCUMENTED grown-Float32 accumulation floors (crown_pct = the carbon:335 @test_broken; grown-dbh = MYBA/MYSDI,
        # test_dbs_compute). i.e. the SAME accumulated-Float32-growth primitive (doctrine #9's permitted class),
        # propagated through faithful cwd accounting into the fire-basis fuel and amplified by _fmdyn. Bug #1 (input-snag
        # bole topwood) + the snag-record binning were REAL fixes en route (both landed this session). REFUTED along the
        # way: total-fallvol (regressed 12 tests). @test_broken vs rendered-== (not a padded bound). See task #72.
        @test_broken round(sc; digits = 3) == 17.581
    end
end
