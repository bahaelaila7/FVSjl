# test_list_ht2td.jl — FVS_TreeList / FVS_CutList Ht2TDCF & Ht2TDBF (merch-top heights, HT2TD(:,2)/(:,1)) vs live.
#
# FVS fills HT2TD in NATCRS/FVSVOL (fvsvol.f:337-339 cubic call, :484-487 board call) from NVEL's HT1PRD: for
# Flewelling FW2 equations HT1PRD = MERLEN's LMERCH+STUMP (profile.f:335-342, SF_HS at the product top), for
# region-6 Behre 1+ΣXLEN (r6vol.f:121-125); VOLS zeroes the array first (vols.f:86-90) and dead cycle-0 records
# are filled too (IPASS=2). The goldens are live FVSbm_g16 output for the two BM establishing stands of the
# bm_seedling fixture (bit-exact every cycle), with TREELIST 0 + a cycle-5 THINBBA + CUTLIST 0.
using FVSjl, Test, SQLite, DBInterface

const _HT_DIR = joinpath(@__DIR__, "..", "fixtures", "list_ht2td")
const _HT_DB  = joinpath(@__DIR__, "..", "fixtures", "bm_seedling", "bm_seedling.db")

function _ht_key(cn, outdb)
    join(["STDIDENT", cn, "DATABASE", "DSNOUT", outdb, "TREELIDB", "CUTLIDB", "DSNin", abspath(_HT_DB),
          "StandSQL", "SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL",
          "TreeSQL", "SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'", "EndSQL", "END",
          "NUMCYCLE         10", "TREELIST         0.0", "THINBBA          5.0      80.0", "CUTLIST          0.0",
          "ECHOSUM", "PROCESS", "STOP"], "\n") * "\n"
end

function _ht_rows(db)
    out = Dict{Tuple{String,Int,String,Int},NTuple{4,Float64}}()
    d = SQLite.DB(db)
    try
        names = [r.name for r in SQLite.tables(d)]
        for tbl in ("FVS_TreeList", "FVS_CutList")
            tbl in names || continue
            for r in DBInterface.execute(d, "SELECT Year, TreeId, TreeIndex, DBH, Ht, Ht2TDCF, Ht2TDBF FROM $tbl")
                out[(tbl, Int(r.Year), string(r.TreeId), Int(r.TreeIndex))] =
                    (Float64(r.DBH), Float64(r.Ht), Float64(r.Ht2TDCF), Float64(r.Ht2TDBF))
            end
        end
    finally
        SQLite.close(d)
    end
    out
end

@testset "FVS_TreeList/CutList Ht2TDCF & Ht2TDBF vs live (BM)" begin
    if !isfile(_HT_DB)
        @test_skip "bm_seedling fixture db missing"
    else
        for cn in ("449747082489998", "15144796010497")
            gold = joinpath(_HT_DIR, "$cn.live.csv")
            isfile(gold) || (@test_skip "$cn golden missing"; continue)
            live = Dict{Tuple{String,Int,String,Int},NTuple{4,Float64}}()
            for (k, ln) in enumerate(eachline(gold))
                k == 1 && continue
                f = split(ln, ',')
                live[(f[1], parse(Int, f[2]), f[3], parse(Int, f[4]))] =
                    (parse(Float64, f[5]), parse(Float64, f[6]), parse(Float64, f[7]), parse(Float64, f[8]))
            end
            dir = mktempdir(); outdb = joinpath(dir, "out.db"); kf = joinpath(dir, "s.key")
            write(kf, _ht_key(cn, outdb))
            FVSjl.run_keyfile(kf; variant = FVSjl.BlueMountains())
            jl = _ht_rows(outdb)
            @test Set(keys(jl)) == Set(keys(live))             # same rows written (gating + ids)
            @test count(v -> v[3] != 0, values(live)) > 0      # the stands really exercise merch-sized trees
            # HT2TD is a pure function of the tree's (DBH, HT, species, bark): assert it bit-exact on every row whose
            # per-tree DBH and HT equal live. Rows whose DBH/HT themselves differ by ~1 ULP carry the known upstream
            # per-tree growth drift (hidden at .sum precision; BM none-residual campaign) — tracked, not asserted here.
            ngeomdiff = 0
            for (k, l) in live
                j = get(jl, k, nothing); j === nothing && continue
                if j[1] == l[1] && j[2] == l[2]
                    @test (j[3], j[4]) == (l[3], l[4])
                else
                    ngeomdiff += 1
                end
            end
            @test_broken ngeomdiff == 0                          # upstream per-tree DBH/HT ULP drift (open)
        end
    end
end
