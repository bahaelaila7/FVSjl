# build_subdb.jl — Pillar-1 scale enabler. Copy the FVS-ready rows (FVS_STANDINIT_COND + FVS_TREEINIT_COND)
# for a stratified sample of N stands of one VARIANT from the read-only master DB into a SMALL, INDEXED
# working DB. Running live-FVS + FVSjl against this sub-DB is ~100× faster than the master (whose STAND_CN
# columns are unindexed ⇒ every per-stand query full-scans the 2.2M/8M-row tables). Master never modified.
#
# The copy is done in SQLite's C engine (ATTACH + CREATE TABLE AS SELECT ... WHERE STAND_CN IN (...)) — one
# C-speed table pass, NOT millions of per-row Julia calls.
#
# Usage: julia --project=. test/harness/fia/build_subdb.jl <VARIANT> <N> <out.db>

import SQLite, DBInterface
const MASTER = "/workspace/SQLite_FIADB_ENTIRE.db"

function build(variant, n, outdb)
    src = SQLite.DB(MASTER)
    cns = String[]
    for r in DBInterface.execute(src, "SELECT STAND_CN FROM FVS_STANDINIT_COND WHERE VARIANT='$variant' AND STAND_CN IS NOT NULL ORDER BY ECOREGION, LOCATION, STAND_CN")
        push!(cns, string(r.STAND_CN))
    end
    SQLite.close(src)
    total = length(cns); n = min(n, total)
    stride = total / n
    idx = unique(clamp.(round.(Int, (0.5:1.0:n) .* stride), 1, total))
    sample = cns[idx]
    inlist = join(["'" * replace(s, "'" => "''") * "'" for s in sample], ",")

    isfile(outdb) && rm(outdb)
    dst = SQLite.DB(outdb)
    DBInterface.execute(dst, "ATTACH DATABASE '$MASTER' AS m")
    for tbl in ("FVS_STANDINIT_COND", "FVS_TREEINIT_COND")
        DBInterface.execute(dst, "CREATE TABLE $tbl AS SELECT * FROM m.$tbl WHERE STAND_CN IN ($inlist)")
        DBInterface.execute(dst, "CREATE INDEX ix_$(tbl)_cn ON $tbl(STAND_CN)")
        cnt = 0; for r in DBInterface.execute(dst, "SELECT COUNT(*) c FROM $tbl"); cnt = r.c; end
        println(stderr, "  $tbl: $cnt rows")
    end
    DBInterface.execute(dst, "DETACH DATABASE m")
    open(outdb * ".stands", "w") do io
        for s in sample; println(io, s, '\t', variant); end
    end
    println(stderr, "VARIANT=$variant population=$total sampled=$(length(sample)) → $outdb (+ .stands)")
    return length(sample)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 3 || error("usage: build_subdb.jl <VARIANT> <N> <out.db>")
    build(ARGS[1], parse(Int, ARGS[2]), ARGS[3])
end
