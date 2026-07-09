# expand_batch.jl — emit the next batch of FVS-ready STAND_CNs for a variant, in the SAME deterministic order
# as extract_sample.jl (ECOREGION, LOCATION, STAND_CN), sliced [offset, offset+limit). Read-only on the master.
# This drives the full-population coverage sweep (goal: eventually cover ALL FVS-ready FIA for the 4 ported
# variants SN/NE/CS/LS). A cursor per variant (test/harness/fia/expand/<v>.cursor) records the next offset.
#
# Usage: julia --project=. test/harness/fia/expand_batch.jl <VARIANT> <offset> <limit> <out.stands>
#   writes up to <limit> lines "STAND_CN<TAB>VARIANT"; prints the population total + emitted count to stderr.
import SQLite, DBInterface
const DB = "/workspace/SQLite_FIADB_ENTIRE.db"
function main(v, offset, limit, out)
    db = SQLite.DB(DB)
    tot = 0
    for r in DBInterface.execute(db, "SELECT COUNT(*) AS n FROM FVS_STANDINIT_COND WHERE VARIANT='$v' AND STAND_CN IS NOT NULL")
        tot = Int(r.n)
    end
    q = """SELECT STAND_CN FROM FVS_STANDINIT_COND
           WHERE VARIANT='$v' AND STAND_CN IS NOT NULL
           ORDER BY ECOREGION, LOCATION, STAND_CN
           LIMIT $limit OFFSET $offset"""
    n = 0
    open(out, "w") do io
        for r in DBInterface.execute(db, q)
            println(io, string(r.STAND_CN), '\t', v); n += 1
        end
    end
    println(stderr, "VARIANT=$v population=$tot offset=$offset emitted=$n remaining=$(max(0, tot - offset - n))")
    println(n)  # stdout: emitted count (for the orchestrator)
end
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 4 || error("usage: expand_batch.jl <VARIANT> <offset> <limit> <out.stands>")
    main(ARGS[1], parse(Int, ARGS[2]), parse(Int, ARGS[3]), ARGS[4])
end
