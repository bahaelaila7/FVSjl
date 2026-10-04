# test_bc_db_live.jl — BC DATABASE-input stand vs LIVE FVSbc.
#
# bc_skyranch.db = the FVS tests/FVSbc YSM-SkyRanch "SkyRanch-Control" stand (978 records on 3 plots; BASAL_AREA_FACTOR
# −20, INV_PLOT_SIZE 20, no NUM_PLOTS). The stock FVSbc_clean SIGSEGVs reading FVS_TreeInit (intree.f passes DBSTREESIN
# 30 of its 34 arguments — the same stale call ON's intree.f had); the golden comes from the same objects relinked with
# that call completed (private scratch build). Goldens: the inventory-year FVS_TreeList_Metric (TPH, DBH, Ht per record)
# and the cycle-0 .sum row.

using Test, FVSjl, SQLite, DBInterface
const _BDB = FVSjl
const _BDB_FX = joinpath(@__DIR__, "..", "fixtures", "bc_db")

const _BDB_ROWS, _BDB_TL = let dir = mktempdir()
    for f in ("bc_skyranch.key", "bc_skyranch.db"); cp(joinpath(_BDB_FX, f), joinpath(dir, f)); end
    txt = cd(() -> _BDB.run_keyfile("bc_skyranch.key"; variant = _BDB.BritishColumbia(), output = :sum), dir)
    db = SQLite.DB(joinpath(dir, "bc_skyranch.out.db"))
    ([split(l) for l in split(txt, '\n') if occursin(r"^\d{4} ", l)],
     Dict(Int(r.TreeIndex) => (Float64(r.TPH), Float64(r.DBH), Float64(r.Ht))
          for r in DBInterface.execute(db, "SELECT TreeIndex,TPH,DBH,Ht FROM FVS_TreeList_Metric WHERE Year=2019")))
end

# (1) canada/bc dbsstandin.f:357-376 (the METRIC reader, as ON's): BASAL_AREA_FACTOR <0 ⇒ /HAtoACR, INV_PLOT_SIZE /HAtoACR,
# BRK_DBH ·CMtoIN; TREE_COUNT raw. (2) initre.f:329-333: with no NUM_PLOTS/NONSTK_PLOTS (grinit IPTINV=NONSTK=−9999)
# IPTINV/NONSTK are the COUNTED plots IPTKNT / non-stockable NSTKNT (intree.f:324-362) — here 3 plots, so PI=3. jl
# pre-scaled TREE_COUNT ×0.40468564 and kept PI=1: every record 19.9996 TPH vs live 6.6667.
@testset "BC SkyRanch DB 2019: TPH/DBH/Ht per record + cycle-0 .sum row == live (metric design factors, IPTINV=IPTKNT)" begin
    n = 0
    for ln in eachline(joinpath(_BDB_FX, "bc_skyranch_2019_tl_live.csv"))
        (startswith(ln, "#") || startswith(ln, "Year")) && continue
        f = split(ln, ','); k = parse(Int, f[2]); n += 1
        @test (k, get(_BDB_TL, k, (NaN, NaN, NaN))) == (k, (parse(Float64, f[3]), parse(Float64, f[4]), parse(Float64, f[5])))
    end
    @test n == length(_BDB_TL) == 978
    @test first(_BDB_ROWS) == split(first(readlines(joinpath(_BDB_FX, "bc_skyranch_live.rows"))))
end
