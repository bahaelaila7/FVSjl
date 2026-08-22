# FVS_CalibStats DBS table (dbscalib.f, CALBSTDB toggle) — the DGSCOR large-tree DG-calibration sample statistics
# per calibrated species (NUMCAL / STDRAT / WC / CORTEM), serialized to SQLite. Ported 2026-08-22: jl already
# computes every quantity in the shared calibration (calibrate_diameter_growth!, dgdriv.f framework) — it just
# discarded them; captured cal_ntree/cal_stdrat/cal_wci/cal_cortem there + write_dbs_calibstats!.
# VALIDATED vs the relinked FVScr_clean (CR crt01 stand): ALL 9 data columns BIT-EXACT for both calibrated species
# (WF n=5, ES n=6). ScaleFactor = CORTEM = exp(COR) at CALIBRATION time (before the CORMLT re-scale), NOT exp(the
# final dg_cor). Oracle values dumped 2026-08-22.

using Test
using FVSjl
using SQLite, DBInterface

@testset "FVS_CalibStats DBS table (CALBSTDB, bit-exact vs FVScr_clean)" begin
    fix = joinpath(@__DIR__, "..", "fixtures", "calibstats")
    dir = mktempdir()
    cp(joinpath(fix, "calibstats.tre"), joinpath(dir, "calibstats.tre"))
    key = joinpath(dir, "calibstats.key")
    outdb = joinpath(dir, "out.db")
    open(key, "w") do io
        for ln in readlines(joinpath(fix, "calibstats.key"))
            println(io, strip(ln) == "calibstats_oracle.db" ? outdb : ln)
        end
    end
    FVSjl.run_keyfile(key; variant = FVSjl.variant_from_code("CR"))
    @test isfile(outdb)
    db = SQLite.DB(outdb)
    @test "FVS_CalibStats" in [t.name for t in SQLite.tables(db)]

    rows = Dict{String,Dict{String,Any}}()
    for r in DBInterface.execute(db, "SELECT * FROM FVS_CalibStats")
        rows[strip(String(r.SpeciesFVS))] = Dict{String,Any}(String(k) => getproperty(r, k) for k in propertynames(r))
    end
    @test length(rows) == 2                                # WF + ES calibrated

    wf = rows["WF"]
    @test strip(String(wf["TreeSize"])) == "LG"
    @test strip(String(wf["SpeciesPLANTS"])) == "ABCO"
    @test strip(String(wf["SpeciesFIA"])) == "015"
    @test wf["NumTrees"] == 5
    @test isapprox(Float64(wf["ScaleFactor"]),   1.4653407; atol = 1e-4)   # CORTEM = exp(COR)
    @test isapprox(Float64(wf["StdErrRatio"]),   2.0792110; atol = 1e-4)
    @test isapprox(Float64(wf["WeightToInput"]), 0.8554077; atol = 1e-4)   # WC
    @test isapprox(Float64(wf["ReadCorMult"]),   1.5631036; atol = 1e-4)   # exp(log(CORTEM)/WC)

    es = rows["ES"]
    @test es["NumTrees"] == 6
    @test isapprox(Float64(es["ScaleFactor"]),   0.8065776; atol = 1e-4)
    @test isapprox(Float64(es["StdErrRatio"]),   1.5855684; atol = 1e-4)
    @test isapprox(Float64(es["WeightToInput"]), 0.7788284; atol = 1e-4)
end
