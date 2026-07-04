# test_cycleat.jl — CYCLEAT keyword (opt 134, initre.f:13400 + fvs.f:116-135).
#
# CYCLEAT inserts an extra cycle boundary at a requested calendar year, splitting an existing
# cycle (without extending the end or moving the start) and bumping the effective cycle count.
# This exercises the non-uniform cycle schedule (`build_cycle_schedule!`, the FVS IY array) that
# also underpins per-cycle TIMEINT and per-cycle GROWTH. Checks:
#   1. SCHEDULE — build_cycle_schedule! inserts the right boundary + bumps ncycle_eff (unit);
#   2. FORTRAN years/PrdLen — the projected `.sum` YEAR column and period lengths are bit-exact
#      vs live Fortran (the schedule mechanism is exact);
#   3. FORTRAN stand — TPA/BA/SDI track Fortran within the documented TIMEINT period-scaling
#      residual (the split makes 3- and 2-yr cycles; calibrated species sp33/65 carry a small
#      non-5-yr period-scaling tail, present in plain uniform TIMEINT too — NOT a CYCLEAT bug).

using Test, FVSjl

const _CA_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_ca_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]
_ca_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                                y !== nothing && 1900 < y < 2100)]

@testset "CYCLEAT → extra cycle boundary (non-uniform schedule)" begin
    # 1. SCHEDULE — a 10-cycle, 5-yr run (1990..2040); CYCLEAT 2003 + 2018 split two cycles.
    s = FVSjl.StandState(FVSjl.Southern())
    s.control.cycle_year[1] = Int32(1990)
    s.control.ncycle = Int32(10)
    s.control.year = 5f0
    push!(s.control.cycleat_years, Int32(2003)); push!(s.control.cycleat_years, Int32(2018))
    FVSjl.build_cycle_schedule!(s)
    @test s.control.ncycle_eff == 12                       # 10 + 2 inserted
    @test FVSjl.cycle_year_at(s.control, 3) == 2003        # inserted between 2000 and 2005
    @test FVSjl.cycle_year_at(s.control, 7) == 2018        # inserted between 2015 and 2020
    @test FVSjl.cycle_year_at(s.control, 12) == 2040       # end NOT extended
    @test FVSjl.cycle_period_at(s.control, 2) == 3         # the 2000→2003 split cycle is 3 yr
    @test FVSjl.cycle_period_at(s.control, 3) == 2         # the 2003→2005 split cycle is 2 yr
    # a CYCLEAT year on an existing boundary, or outside the run, is a no-op (no extra cycle)
    s2 = FVSjl.StandState(FVSjl.Southern())
    s2.control.cycle_year[1] = Int32(1990); s2.control.ncycle = Int32(10); s2.control.year = 5f0
    push!(s2.control.cycleat_years, Int32(2005))   # already a boundary
    push!(s2.control.cycleat_years, Int32(2099))   # past the end
    FVSjl.build_cycle_schedule!(s2)
    @test s2.control.ncycle_eff == 10

    # 2/3. FORTRAN — the cycleat scenario (CYCLEAT 2003 splits 2000→2005) vs live Fortran.
    key = joinpath(_CA_DIR, "cycleat.key")
    sav = joinpath(_CA_DIR, "cycleat.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "cycleat scenario not available"
    else
        jl = _ca_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _ca_base(sav)
        @test length(jl) == length(ft) == 7
        for (j, f) in zip(jl, ft)
            @test j[1] == f[1]                              # YEAR column bit-exact (schedule)
            @test j[23] == f[23]                            # PrdLen column bit-exact (period length)
        end
        # the inserted boundary is present at 2003 with a 3-yr then 2-yr period
        @test jl[4][1] == "2003" && jl[3][23] == "3" && jl[4][23] == "2"
        # Stand columns are BIT-EXACT vs live: CYCLEAT-2003 splits 2000→2005 into a 3-yr + 2-yr
        # period, and jl's per-cycle growth reproduces the split exactly (re-measured — the old
        # ≤8/≤3/≤6 "non-5-yr period residual" was stale over-caution; every column now matches).
        for (j, f) in zip(jl, ft)
            @test parse(Int, j[3]) == parse(Int, f[3])     # TPA — BIT-EXACT
            @test parse(Int, j[4]) == parse(Int, f[4])     # BA  — BIT-EXACT
            @test parse(Int, j[5]) == parse(Int, f[5])     # SDI — BIT-EXACT
        end
    end
end
