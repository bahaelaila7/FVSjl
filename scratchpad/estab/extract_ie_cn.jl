using SQLite, DBInterface
src = SQLite.DB("/workspace/SQLite_FIADB_ENTIRE.db")
fix = ARGS[1]; sid = ARGS[2]
isfile(fix) && rm(fix)
dst = SQLite.DB(fix)
for tbl in ("FVS_STANDINIT_COND", "FVS_TREEINIT_COND")
    schsql = ""
    for r in DBInterface.execute(src, "SELECT sql FROM sqlite_master WHERE type='table' AND name='$tbl'")
        schsql = r.sql
    end
    DBInterface.execute(dst, schsql)
    stmt = nothing; n = 0
    for r in DBInterface.execute(src, "SELECT * FROM $tbl WHERE STAND_CN='$sid'")
        nt = NamedTuple(r)
        if stmt === nothing
            cols = join(string.(keys(nt)), ",")
            qs = join(fill("?", length(keys(nt))), ",")
            stmt = DBInterface.prepare(dst, "INSERT INTO $tbl ($cols) VALUES ($qs)")
        end
        DBInterface.execute(stmt, collect(values(nt)))
        n += 1
    end
    println("$tbl: $n rows")
end
