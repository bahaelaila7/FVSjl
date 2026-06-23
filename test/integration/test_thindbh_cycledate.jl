# test_thindbh_cycledate.jl — a thin scheduled with a cycle-number date (a small date, FVS
# 1-based) must fire in the matching cycle. Scoped from initre.f:1189 (blank date field →
# IDT=1, a cycle number) + OPNEW/OPFIND (dates < 1000 are cycle numbers).
#
# Regression for the snt01 stand-4 divergence: `THINDBH 3.` parsed its blank date as year 0,
# which never matched a calendar-year cycle, so the thin silently never fired. With the fix
# (blank date → cycle 1; cuts! treats dates < 1000 as cycle numbers, FVS cycle = FVSjl
# cycle + 1) it fires in cycle 1, and snt01 stand 4 is bit-exact vs Fortran through 2003.

using Test, FVSjl
using FVSjl: StandState, Southern, init_blockdata!, init_merch_standards!, cuts!, ScheduledActivity

@testset "cycle-number-dated thin fires in the matching cycle" begin
    function stand()
        s = StandState(Southern()); init_blockdata!(s, s.variant); init_merch_standards!(s)
        s.control.cycle_year[1] = Int32(1990)
        t = s.trees; t.n = 4
        for (i, d) in enumerate((16f0, 12f0, 8f0, 4f0))
            t.species[i] = Int32(65); t.dbh[i] = d; t.height[i] = 70f0; t.tpa[i] = 40f0
            t.crown_pct[i] = Int32(45)
        end
        s
    end

    # THINBTA (icflag 3) to 40 ft² BA, scheduled with the cycle-NUMBER date 1 (= FVS cycle 1 =
    # FVSjl cycle 0). It must fire in cycle 0 even though the calendar year (1990) ≠ 1.
    s = stand(); s.control.cycle = Int32(0)
    push!(s.control.schedule, ScheduledActivity(Int32(1), Int32(3), (40f0,1f0,0f0,0f0,0f0,0f0)))
    @test cuts!(s; fint = 5f0).tpa > 0f0                  # fired in cycle 0

    # the same date does NOT fire in a later cycle (cycle-number 1 ≠ FVS cycle 2)
    s2 = stand(); s2.control.cycle = Int32(1)
    push!(s2.control.schedule, ScheduledActivity(Int32(1), Int32(3), (40f0,1f0,0f0,0f0,0f0,0f0)))
    @test cuts!(s2; fint = 5f0).tpa == 0f0

    # cycle-number 2 fires in FVSjl cycle 1 (= FVS cycle 2), not cycle 0
    s3 = stand(); s3.control.cycle = Int32(0)
    push!(s3.control.schedule, ScheduledActivity(Int32(2), Int32(3), (40f0,1f0,0f0,0f0,0f0,0f0)))
    @test cuts!(s3; fint = 5f0).tpa == 0f0
    s3.control.cycle = Int32(1); empty!(s3.control.years_cut)
    @test cuts!(s3; fint = 5f0).tpa > 0f0

    # a real calendar-year date (≥ 1000) still matches by year, not as a cycle number
    s4 = stand(); s4.control.cycle = Int32(2)            # year = 1990 + 2·5 = 2000
    push!(s4.control.schedule, ScheduledActivity(Int32(2000), Int32(3), (40f0,1f0,0f0,0f0,0f0,0f0)))
    @test cuts!(s4; fint = 5f0).tpa > 0f0
end
