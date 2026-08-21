# FVS_StrClass DBS table (dbsstrclass.f, STRCLSDB toggle) — the SSTAGE stand-structure classification (up to 3
# height strata × DBHNOM/heights/crown-base/cover/dominant-species/status + strata count, cover, class) serialized
# to SQLite. Ported 2026-08-21: jl's `structure_report` already computes the per-stratum data (validated vs the
# sstage.f `.out` report); write_dbs_strclass! serializes it, two rows/cycle (Removal_Code 0=before / 1=after thin).
# Along the way fixed a REAL gap: `_ss_strata` didn't dispatch OC to `oc_cwcalc` (fell to the generic crown_width →
# 0 cover → wrong strata); added the OregonCoast branch. VALIDATED vs the relinked FVSoc_clean (OC 30-tree stand):
# cyc-0 (1993) + 1998 + 2008 ALL 43 columns BIT-EXACT; later cycles (2003 DBHNOM ~1%, 2013 ht/cover ±1-NINT) corner
# on the OC multi-cycle growth straddle (the structure faithfully reflects the cornered tree list). Oracle 2026-08-21.

using Test
using FVSjl
using SQLite, DBInterface

@testset "FVS_StrClass DBS table (STRCLSDB, cyc0 bit-exact vs FVSoc_clean)" begin
    fix = joinpath(@__DIR__, "..", "fixtures", "strclass")
    dir = mktempdir()
    cp(joinpath(fix, "strclass.tre"), joinpath(dir, "strclass.tre"))
    key = joinpath(dir, "strclass.key")
    outdb = joinpath(dir, "out.db")
    open(key, "w") do io
        for ln in readlines(joinpath(fix, "strclass.key"))
            println(io, strip(ln) == "strclass_oracle.db" ? outdb : ln)
        end
    end
    FVSjl.run_keyfile(key; variant = FVSjl.variant_from_code("OC"))
    @test isfile(outdb)
    db = SQLite.DB(outdb)
    @test "FVS_StrClass" in [t.name for t in SQLite.tables(db)]

    # materialize the cyc-0 (1993, Removal_Code 0) row
    row = nothing
    for r in DBInterface.execute(db, "SELECT * FROM FVS_StrClass WHERE Year=1993 AND Removal_Code=0 LIMIT 1")
        row = Dict{String,Any}(String(k) => getproperty(r, k) for k in propertynames(r))
    end
    @test row !== nothing

    # cyc-0 stratum-1 (dominant): BIT-EXACT vs the oracle (deterministic on the bit-exact inventory tree list)
    @test isapprox(Float64(row["Stratum_1_DBH"]), 7.169552803039551; atol = 1e-4)
    @test row["Stratum_1_Nom_Ht"] == 37
    @test row["Stratum_1_Lg_Ht"] == 75
    @test row["Stratum_1_Sm_Ht"] == 2
    @test row["Stratum_1_Crown_Base"] == 13
    @test row["Stratum_1_Crown_Cover"] == 50
    @test strip(String(row["Stratum_1_SpeciesFVS_1"])) == "GF"    # grand fir dominant (crown-area)
    @test strip(String(row["Stratum_1_SpeciesFVS_2"])) == "DF"
    @test strip(String(row["Stratum_1_SpeciesPLANTS_1"])) == "ABGR"
    @test strip(String(row["Stratum_1_SpeciesFIA_1"])) == "017"
    @test row["Stratum_1_Status_Code"] == 2                       # dominant stratum
    # single-stratum stand: stratum 2/3 absent → "--" / 0
    @test strip(String(row["Stratum_2_SpeciesFVS_1"])) == "--"
    @test Float64(row["Stratum_2_DBH"]) == 0.0
    # whole-stand summary
    @test row["Number_of_Strata"] == 1
    @test row["Total_Cover"] == 50
    @test strip(String(row["Structure_Class"])) == "1=SI"          # stand initiation
end
