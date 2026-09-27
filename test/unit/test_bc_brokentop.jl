# test_bc_brokentop.jl — BC .tre record layout + broken-top (ITRUNC/NORMHT) volume, vs live FVSbc_clean.
#
# bc/blkdat.f TREFMT is the metric layout (F5.1 DBH, F4.1 DG, 2F5.1 HT/THT). jl used the SN default, which
# read HT as blank and the real height as the top-kill height, so every BC tree became top-killed at its own
# height. And bc/vols.f:137-138 volumes a top-killed tree at NORMHT truncated at ITRUNC.
# Fixture: 5 trees (3 broken tops, THT 60/70/50 m), cycle-0 total/merch cubic from FVS_TreeList_Metric.
using FVSjl, Test, SQLite, DBInterface

@testset "BC tree format + broken-top volume vs FVSbc_clean" begin
    src = joinpath(@__DIR__, "..", "fixtures", "bc_brokentop")
    dir = mktempdir()
    for f in ("bctop.key", "bctop.tre"); cp(joinpath(src, f), joinpath(dir, f)); end
    cd(dir) do
        FVSjl.run_keyfile("bctop.key"; variant = FVSjl.BritishColumbia(), output = :sum)
    end
    db = SQLite.DB(joinpath(dir, "bctop.db"))
    rows = [(r.TreeId, r.TCuFt, r.MCuFt, r.TruncHt) for r in
            DBInterface.execute(db, "SELECT TreeId, TCuFt, MCuFt, TruncHt FROM FVS_TreeList WHERE Year = 1992 ORDER BY TreeId")]
    SQLite.close(db)
    live_tcum = Dict("1" => 0.4899158477783203, "2" => 0.6418090462684631, "3" => 1.223409652709961,
                     "4" => 0.21288253366947174, "5" => 0.31424522399902344)      # FVSbc_clean TCuM (m³)
    live_mcum = Dict("3" => 1.1557263135910034)
    live_trunc_m = Dict("1" => 0, "2" => 60, "3" => 70, "4" => 0, "5" => 50)
    @test length(rows) == 5
    for (id, tcf, mcf, trh) in rows
        id = strip(id)
        @test tcf * 0.0283168 ≈ live_tcum[id] rtol = 5e-5
        @test mcf * 0.0283168 ≈ get(live_mcum, id, 0.0) atol = 1e-4
        @test (trh > 0) == (live_trunc_m[id] > 0)          # only the THT trees are top-killed
    end
end
