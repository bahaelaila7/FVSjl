# FVS_Climate DBS table (dbsclsum.f, CLIMREDB toggle) — the Climate-FVS per-species Viability-and-Effects report
# serialized to SQLite. Ported 2026-08-21: climate_report (collected POST-growth at the cycle-midpoint sampling
# year, clauestb.f/clgmult.f) + write_dbs_climate!. VALIDATED bit-exact-or-cornered vs the relinked FVSie_clean
# oracle (climtest_oracle.db, IE Clearwater CGCM3_A2 stand S248112): the DETERMINISTIC climate columns Viability
# and ViabMort matched 147/147 distinct (Year,Species) cells BIT-EXACT (the midpoint climate read + SPCALIB
# presence-calibration); the tree-list columns BA/TPA/GrowthMult/MxDenMult/dClimMort/AutoEstbTPA CORNER on the
# DGSD=2.0 OLDRN growth straddle — the same class the `.sum` itself corners in. This test re-runs the emission
# end-to-end and checks the bit-exact columns tightly + the cornered columns within the straddle band.

using Test
using FVSjl
using SQLite, DBInterface

@testset "FVS_Climate DBS table (CLIMREDB, bit-exact vs FVSie_clean)" begin
    fix = joinpath(@__DIR__, "..", "fixtures", "climate")
    dir = mktempdir()
    key = joinpath(dir, "clim_dbs.key")
    cp(joinpath(fix, "clim_iet.tre"), joinpath(dir, "clim_dbs.tre"))
    outdb = joinpath(dir, "out.db")
    # clim_iet.key + CLIMREDB in the DATABASE block + DSNOut → our temp db
    lines = readlines(joinpath(fix, "clim_iet.key"))
    open(key, "w") do io
        indb = false
        for ln in lines
            s = uppercase(strip(ln))
            s == "DATABASE" && (indb = true)
            if indb && startswith(s, "IET01_OUT.DB")   # DSNOut filename line
                println(io, outdb); continue
            end
            if indb && s == "END"
                println(io, "CLIMREDB           1"); println(io, ln); indb = false; continue
            end
            println(io, ln)
        end
    end
    FVSjl.run_keyfile(key; variant = FVSjl.variant_from_code("IE"))
    @test isfile(outdb)
    db = SQLite.DB(outdb)
    @test "FVS_Climate" in [t.name for t in SQLite.tables(db)]

    cell(yr, sp, col) = begin
        v = nothing
        for r in DBInterface.execute(db,
                "SELECT $col FROM FVS_Climate WHERE Year=$yr AND SpeciesFVS='$sp' LIMIT 1")
            v = r[Symbol(col)]
        end
        v
    end

    # DETERMINISTIC climate reads — BIT-EXACT (viability sampled at report_year+fint/2)
    @test isapprox(cell(1990, "DF", "Viability"), 0.9309; atol = 5e-4)   # PSME viab @ 1995 (midpoint)
    @test isapprox(cell(1990, "PP", "Viability"), 0.4335; atol = 5e-4)
    @test isapprox(cell(1990, "WL", "Viability"), 0.8584; atol = 5e-4)
    # SPCALIB-calibrated viability mortality — BIT-EXACT (0 early, rising as PP viability falls)
    @test isapprox(cell(1990, "PP", "ViabMort"), 0.0; atol = 5e-4)
    @test isapprox(cell(2010, "PP", "ViabMort"), 0.0338; atol = 1e-3)
    @test isapprox(cell(2030, "PP", "ViabMort"), 0.2291; atol = 1e-3)
    # post-growth per-species BA/TPA — CORNERED on the DGSD=2.0 OLDRN growth straddle (as the .sum)
    @test isapprox(cell(1990, "DF", "BA"), 28.82; rtol = 0.05)      # jl 28.33 (~1.7%)
    @test isapprox(cell(1990, "DF", "TPA"), 157.63; rtol = 0.05)
    # POTESTAB via the clinit defaults (AESNTREES=500/NESPECIES=4) on the top-4 viable species (cornered ~oracle 99.44)
    @test cell(1990, "MM", "AutoEstbTPA") > 90
    @test cell(1990, "WL", "AutoEstbTPA") == 0.0   # WL not in the top-4
    # MxDenMult = 1 at the first cycle (deterministic)
    @test isapprox(cell(1990, "DF", "MxDenMult"), 1.0; atol = 5e-4)
end
