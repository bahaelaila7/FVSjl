# Extract the sampled stands into a small INDEXED subset DB (fast repeated validation).
# ATTACHes the 66GB source READ-ONLY and copies only the matching STANDINIT_COND +
# TREEINIT_COND rows into a fresh DB, then indexes STAND_CN. Usage:
#   julia --project=. extract_subset.jl <sampledir> <out.db>
using SQLite
const SRC = "/workspace/SQLite_FIADB_ENTIRE.db"
function main(sampledir, outdb)
    cns = String[]
    for v in ("sn","ne","cs","ls")
        f = joinpath(sampledir, "$(v)_sample.txt"); isfile(f) || continue
        for l in eachline(f); s=strip(l); isempty(s) || push!(cns, split(s,'\t')[1]); end
    end
    println("sampled STAND_CN: ", length(cns))
    isfile(outdb) && rm(outdb)
    db = SQLite.DB(outdb)
    DBInterface.execute(db, "ATTACH DATABASE 'file:$SRC?mode=ro&immutable=1' AS src")
    # values-list of quoted CNs
    inlist = join(["'" * replace(c, "'"=>"''") * "'" for c in cns], ",")
    for t in ("FVS_STANDINIT_COND","FVS_TREEINIT_COND")
        DBInterface.execute(db, "CREATE TABLE $t AS SELECT * FROM src.$t WHERE STAND_CN IN ($inlist)")
        n = first(DBInterface.execute(db, "SELECT COUNT(*) c FROM $t"))[1]
        DBInterface.execute(db, "CREATE INDEX idx_$(t)_cn ON $t(STAND_CN)")
        println("  $t: $n rows copied + indexed")
    end
    println("subset DB written: $outdb")
end
main(ARGS[1], ARGS[2])
