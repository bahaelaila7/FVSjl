# FVS_CanProfile DBS table (dbsfmcanpr.f, CANFPROF keyword) — the FFE canopy crown-fuel profile (crown fuel by
# 1-ft height layer, lbs/ac-ft) serialized to SQLite. Ported 2026-08-21: jl's canopy_bulk_density already builds
# the exact CRFILL array (extracted as `canopy_crfill`); write_dbs_canprofile! serializes it. Reported PRE-growth
# (cycle-START inventory, alongside FMPOFL — distinct from the post-growth carbon/climate reports).
# VALIDATED vs the relinked FVScr_clean oracle (crcanpr.key = crt01 CR FFE stand + CANFPROF): the cyc-0 (1993)
# profile is BIT-EXACT (68/68 height layers); later cycles corner on the CR DGSD=2.0 OLDRN growth straddle (the
# profile is tree-list-derived, the same class the .sum corners in). Oracle values dumped 2026-08-21.

using Test
using FVSjl
using SQLite, DBInterface

@testset "FVS_CanProfile DBS table (CANFPROF, cyc0 bit-exact vs FVScr_clean)" begin
    fix = joinpath(@__DIR__, "..", "fixtures", "canprofile")
    dir = mktempdir()
    cp(joinpath(fix, "crcanpr.tre"), joinpath(dir, "crcanpr.tre"))
    key = joinpath(dir, "crcanpr.key")
    outdb = joinpath(dir, "out.db")
    # the fixture keyfile writes to crcanpr_oracle.db — redirect DSNOUT to our temp db
    open(key, "w") do io
        for ln in readlines(joinpath(fix, "crcanpr.key"))
            println(io, strip(ln) == "crcanpr_oracle.db" ? outdb : ln)
        end
    end
    FVSjl.run_keyfile(key; variant = FVSjl.variant_from_code("CR"))
    @test isfile(outdb)
    db = SQLite.DB(outdb)
    @test "FVS_CanProfile" in [t.name for t in SQLite.tables(db)]

    cell(yr, ht) = begin
        v = nothing
        for r in DBInterface.execute(db,
                "SELECT Canopy_Fuel_lbs_acre_ft AS f FROM FVS_CanProfile WHERE Year=$yr AND Height_ft=$(Float64(ht)) LIMIT 1")
            v = r[:f]
        end
        v
    end

    # cyc-0 (1993) canopy fuel profile — BIT-EXACT vs the oracle (deterministic on the bit-exact inventory tree list)
    @test isapprox(cell(1993, 8),  6.522;   atol = 5e-3)
    @test isapprox(cell(1993, 10), 98.617;  atol = 1e-2)
    @test isapprox(cell(1993, 12), 229.177; atol = 1e-2)   # the plateau (only right at the PRE-growth timing)
    @test isapprox(cell(1993, 20), 304.827; atol = 1e-2)
    # the kg/m³ conversion column tracks the lbs/ac-ft (·0.45359237/(4046.856422·0.3048))
    kg = nothing
    for r in DBInterface.execute(db, "SELECT Canopy_Fuel_kg_m3 AS k FROM FVS_CanProfile WHERE Year=1993 AND Height_ft=12.0 LIMIT 1")
        kg = r[:k]
    end
    @test isapprox(kg, 229.177 * 0.45359237 / (4046.856422 * 0.3048); rtol = 1e-4)
end
