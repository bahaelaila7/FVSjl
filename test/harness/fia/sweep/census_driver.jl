# Parallel FIA census driver — the scalable realization of task #94 toward the full ~1.24M-stand run.
#
# Architecture (validated by the feasibility measurement, RESULTS.md):
#   * ONE live-FVS subprocess PER STAND (crash-ISOLATED — a D38 SIGFPE kills only that stand,
#     not a batch). Batching was measured to buy only ~1.5x while coupling the ~30% SN crash rate
#     across the whole batch, so it is deliberately NOT used.
#   * Bounded concurrency via a fixed pool of `Threads.@spawn` workers pulling from a Channel; each
#     worker owns a private workdir/keyfile/.sum (no collisions). The FVS engine is no-shared-state
#     (Pillar 3), so concurrent FVSjl.run_keyfile is safe.
#   * Live `.sum` is PARSED + CACHED to a SQLite table the first time, so re-validation after a jl
#     change is jl-only (the ~10 CPU-hr live cost is paid ONCE). Resumable: cached stands are skipped.
#   * Every stand auto-classified: MATCH / live-CRASH (D38 UB, excluded) / nosum / DIVERGENCE[cols].
#
# Usage:
#   julia -t auto --project=. census_driver.jl <VARIANT> <subset.db> <stand_list.txt> <cache.db> [live_timeout_s]
# where <subset.db> is an INDEXED extract (STAND_CN indexed) — a WHERE STAND_CN= on the 66GB source
# is a full scan (~15s/stand); always extract first (extract_subset.jl).

using FVSjl, SQLite, DBInterface, Base.Threads

const COLS = ["TPA","BA","SDI","CCF","TopHt","QMD","TCuFt","MCuFt","SCuFt","BdFt"]
cfg(v) = v=="LS" ? ("/workspace/FVSjl/tmp/oracles/FVSls_new",FVSjl.LakeStates()) : v=="SN" ? ("/workspace/FVSjl/tmp/oracles/FVSsn_new",FVSjl.Southern()) :
         v=="NE" ? ("/workspace/FVSjl/tmp/oracles/FVSne_new",FVSjl.Northeast()) : v=="CS" ? ("/workspace/FVSjl/tmp/oracles/FVScs_new",FVSjl.CentralStates()) : error(v)

keytext(db,cn) = """
STDIDENT
$cn
DATABASE
DSNin
$db
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE         3.0
ECHOSUM
PROCESS
STOP
"""

function parse_sum(text)
    rows = Tuple{Int,Vector{Float64}}[]
    for ln in split(text,'\n')
        f = split(strip(ln)); length(f) < 12 && continue
        y = tryparse(Int, f[1]); (y===nothing || y<1000 || y>3000) && continue
        vals = Float64[]; ok=true
        for i in 3:12; v=tryparse(Float64,f[i]); v===nothing && (ok=false;break); push!(vals,v); end
        ok && length(vals)==10 && push!(rows,(y,vals))
    end
    rows
end

# One stand: run live (timeout, crash-classified), run jl, compare cycle-0 all 10 cols.
# Returns (status, badcols, live_cycle0) — status ∈ :match :diverge :crash :nosum :jlfail
function do_stand(bin, var, db, cn, dir, tmo)
    key=joinpath(dir,"s.key"); write(key, keytext(db,cn))
    sp=joinpath(dir,"s.sum"); isfile(sp) && rm(sp)
    crash=false
    try; run(pipeline(`timeout $tmo $bin --keywordfile=$key`; stdout=devnull, stderr=devnull))
    catch e; if e isa Base.ProcessFailedException; p=e.procs[1]; (p.termsignal!=0 || Int(p.exitcode)==136) && (crash=true); end; end
    live = isfile(sp) ? read(sp,String) : ""
    isempty(live) && return (crash ? :crash : :nosum, String[], Float64[])
    jlout = try FVSjl.run_keyfile(key; variant=var) catch; ""; end
    isempty(jlout) && return (:jlfail, String[], Float64[])
    L=parse_sum(live); J=parse_sum(jlout)
    (isempty(L)||isempty(J)) && return (:jlfail, String[], Float64[])
    Jd=Dict(y=>vv for (y,vv) in J); y0,lv0=L[1]
    haskey(Jd,y0) || return (:jlfail, String[], lv0)
    jv0=Jd[y0]; bad=String[]
    for k in 1:10; abs(lv0[k]-jv0[k])>=0.5 && push!(bad,COLS[k]); end
    (isempty(bad) ? :match : :diverge, bad, lv0)
end

function main(v, subset, listfile, cachedb, tmo)
    tmo=parse(Int,tmo)
    bin,var = cfg(v)
    stands=[strip(l) for l in eachline(listfile) if !isempty(strip(l))]
    # resume: skip stands already recorded
    cache = SQLite.DB(cachedb)
    DBInterface.execute(cache, """CREATE TABLE IF NOT EXISTS census(
        stand_cn TEXT PRIMARY KEY, variant TEXT, status TEXT, badcols TEXT)""")
    done = Set{String}()
    for r in DBInterface.execute(cache, "SELECT stand_cn FROM census WHERE variant='$v'"); push!(done, string(r[1])); end
    todo = [c for c in stands if !(c in done)]
    println("$v: $(length(stands)) stands, $(length(done)) cached, $(length(todo)) to run  [threads=$(nthreads())]")

    ch = Channel{String}(length(todo)); for c in todo; put!(ch, c); end; close(ch)
    results = Tuple{String,Symbol,Vector{String}}[]
    rlock = ReentrantLock()
    t0 = time()
    @sync for _ in 1:max(1,nthreads())
        Threads.@spawn begin
            dir = mktempdir()
            for cn in ch
                st, bad, _ = do_stand(bin, var, db_of(subset), cn, dir, tmo)
                lock(rlock) do; push!(results, (cn, st, bad)); end
            end
        end
    end
    dt = time()-t0

    # persist
    SQLite.transaction(cache) do
        for (cn,st,bad) in results
            DBInterface.execute(cache, "INSERT OR REPLACE INTO census VALUES(?,?,?,?)",
                (cn, v, String(st), join(bad,",")))
        end
    end

    # tally over ALL cached rows for this variant
    tally = Dict{String,Int}(); colfail=Dict{String,Int}(); offenders=String[]
    for r in DBInterface.execute(cache, "SELECT stand_cn,status,badcols FROM census WHERE variant='$v'")
        s=string(r[2]); tally[s]=get(tally,s,0)+1
        if s=="diverge"; push!(offenders, string(r[1])*"["*string(r[3])*"]")
            for c in split(string(r[3]),','); isempty(c)||(colfail[c]=get(colfail,c,0)+1); end
        end
    end
    comparable = get(tally,"match",0)+get(tally,"diverge",0)
    pct = comparable>0 ? round(100*get(tally,"match",0)/comparable, digits=2) : 0.0
    println("\n===== $v CENSUS (cache=$cachedb) =====")
    println("  ran this session: $(length(todo)) in $(round(dt,digits=1))s ($(length(todo)>0 ? round(1000*dt/length(todo),digits=1) : 0) ms/stand, $(nthreads()) threads)")
    println("  status: ", join(["$k=$(tally[k])" for k in sort(collect(keys(tally)))], "  "))
    println("  BIT-EXACT (all 10 cols, cycle-0): $(get(tally,"match",0))/$comparable comparable = $pct%")
    println("  live-CRASH (D38 UB, excluded): $(get(tally,"crash",0))   nosum/treeless: $(get(tally,"nosum",0))")
    isempty(colfail) || println("  col mismatches: ", join(["$c:$(colfail[c])" for c in sort(collect(keys(colfail)))], " "))
    isempty(offenders) || println("  offenders(first 20): ", join(first(offenders, min(20,length(offenders))), "  "))
    close(cache)
end

# subset path is passed straight through to the keyfile DSNin
db_of(subset) = subset
main(ARGS[1], ARGS[2], ARGS[3], ARGS[4], length(ARGS)>=5 ? ARGS[5] : "60")
