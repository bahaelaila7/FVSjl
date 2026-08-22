# OC ongoing-mortality snag recruitment (regression for the ORGANON mortality→snag gap, 2026-08-22).
# OC's mortality is applied INLINE in diameter_growth!(::OregonCoast) (ORGANON does growth+mort+crown in one
# EXECUTE), which left mortality!(::OregonCoast) a no-op — so the shared book_mortality_snags! path never ran and
# the ongoing (non-fire) ORGANON mortality never reached the FFE standing-snag pool. OC then booked snags only from
# input-dead + fire, so snag_summary decayed to ~0.1x the FVSoc_clean oracle in the tail and missed the first
# cycle's mortality pulse (ocsnag 1993 14.76 vs oracle 55.26; oracle 1998 92.99, pre-fix jl 75.43 -> post-fix 92.99
# BIT-EXACT). Fix: mortality!(::OregonCoast) books the ORGANON-killed density (t.mort_pa) as snags, exactly as OP's
# mortality! and every Wykoff variant do. This test runs a NO-FIRE OC projection and asserts the snag pool RECRUITS
# ongoing mortality (grows well past the input-dead-only baseline); pre-fix it stayed ~14.76 and decayed.

using Test
using FVSjl
using SQLite, DBInterface

@testset "OC books ongoing ORGANON mortality into the FFE snag pool (was input+fire only)" begin
    fix = joinpath(@__DIR__, "..", "fixtures", "oregoncoast")
    dir = mktempdir()
    cp(joinpath(fix, "ocsnag.tre"), joinpath(dir, "ocsnag.tre"))
    key = joinpath(dir, "ocsnag.key")
    outdb = joinpath(dir, "out.db")
    open(key, "w") do io
        for ln in readlines(joinpath(fix, "ocsnag.key"))
            println(io, strip(ln) == "OCSNAGTEST.db" ? outdb : ln)
        end
    end
    FVSjl.run_keyfile(key; variant = FVSjl.variant_from_code("OC"))
    db = SQLite.DB(outdb)
    @test "FVS_SnagSum" in [t.name for t in SQLite.tables(db)]
    tot = Dict{Int,Float64}()
    for r in DBInterface.execute(db, "SELECT Year, Hard_soft_snags_total t FROM FVS_SnagSum ORDER BY Year")
        tot[r.Year] = r.t
    end
    SQLite.close(db)

    # cycle 0 (inventory 1993) = input-dead only (2 records, LP+SP died 1988) = 14.76, matching the oracle.
    @test isapprox(tot[1993], 14.76; atol = 0.05)
    # projected cycles RECRUIT ongoing ORGANON mortality: the pool must climb well past the input baseline.
    # (pre-fix it stayed <= 14.76 and decayed via fall-down; here it grows to ~70 by 2008.)
    @test tot[1998] > 30.0
    @test tot[2008] > tot[1998] > tot[1993]      # monotone recruitment while ongoing mortality > fall-down
    @test tot[2008] > 4 * tot[1993]              # strong recruitment signal (~70 vs 14.76)
end
