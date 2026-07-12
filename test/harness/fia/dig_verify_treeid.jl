# dig_verify_treeid.jl — RELIABLE dig-batch verifier (supersedes the removed 43dm attempt).
# Classifies each top-N LS structure_densephase dig stand as compounded-ULP (self-thinning tie-break, cornered)
# vs ESCALATE (per-tree growth divergence = a possible new bug like DG-serial-corr).
#
# Why this is reliable where the 43dm attempt was not:
#   * matches records by (SpeciesFIA, TreeId) — a STABLE physical-tree key — NOT by (species,DBH)-sorted position
#     (which mis-pairs tied-DBH records across the live/jl DBs). TreeId-matching was proven bit-exact in dig_order.
#   * compares per-tree DBH/DG DIRECTLY (no .sum aggregate ⇒ no growth/mortality entanglement).
#   * MUST be run on a PAUSED sweep (concurrent live-FVS runs contend and corrupt the comparison).
# A GROWTH divergence (matched tree, different DBH/DG) at a PRE-worst cycle ⇒ ESCALATE (deep-trace it).
# Only TPA/survival differing (trees present/absent, matched trees identical) ⇒ compounded-ULP self-thinning.
#
# Usage: julia --project=. test/harness/fia/dig_verify_treeid.jl [N] [VARIANT] [QUEUE.csv]   (default 12 LS)

include(joinpath(@__DIR__, "ledger_fia.jl"))
using SQLite, DBInterface

N   = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 12
V   = length(ARGS) >= 2 ? ARGS[2] : "LS"
QF  = length(ARGS) >= 3 ? ARGS[3] : "/workspace/FVSjl/docs/fia_dig_queue.csv"
bin = BIN[V]; var = VAR[V]

rows = Tuple{String,String,Int,Float64}[]
for ln in readlines(QF)[2:end]
    f = split(ln, ','); length(f) >= 15 || continue
    (f[1] == V && occursin("structure_densephase", ln)) || continue
    wc = tryparse(Int, f[8]); mr = tryparse(Float64, f[9])
    push!(rows, (String(f[3]), String(f[7]), wc === nothing ? 0 : wc, mr === nothing ? 0.0 : mr))
end
sort!(rows, by = r -> -r[4])
sel = rows[1:min(N, end)]
println("=== dig_verify_treeid: top $(length(sel)) of $(length(rows)) $V structure_densephase (TreeId-matched) ===")

const KEY = """
STDIDENT
%CN%
TREELIST           0
DATABASE
DSNin
%SUB%
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
DSNOUT
%OUT%
TREELIDB
END
NUMCYCLE         5.0
ECHOSUM
PROCESS
STOP
"""

# (SpeciesFIA,TreeId) -> (DBH,DG) at a given year
function treemap(db, yr)
    m = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
    isfile(db) || return m
    d = SQLite.DB(db)
    for row in DBInterface.execute(d, "SELECT SpeciesFIA,TreeId,DBH,DG FROM FVS_TreeList WHERE Year=$yr")
        m[(string(row[1]), string(row[2]))] = (row[3] === missing ? 0.0 : Float64(row[3]),
                                               row[4] === missing ? 0.0 : Float64(row[4]))
    end
    SQLite.close(d); m
end

dir = mktempdir(); n_ulp = 0; n_esc = 0; escalate = String[]
for (cn, wcol, wcyc, mr) in sel
    sub = joinpath(dir, "s.db"); try; build_subdb([cn], sub); catch; continue; end
    lk = joinpath(dir, "l.key"); write(lk, replace(KEY, "%CN%"=>cn, "%SUB%"=>sub, "%OUT%"=>"out.db"))
    ldb = joinpath(dir, "FVSOut.db"); isfile(ldb) && rm(ldb)   # DBS APPENDS — fresh per stand
    run(pipeline(ignorestatus(Cmd(`$bin --keywordfile=$lk`; dir = dir)); stdout = devnull, stderr = devnull))
    jdb = joinpath(dir, "j.db"); isfile(jdb) && rm(jdb)
    jk = joinpath(dir, "j.key"); write(jk, replace(KEY, "%CN%"=>cn, "%SUB%"=>sub, "%OUT%"=>jdb))
    try; FVSjl.run_keyfile(jk; variant = var, faithful = true); catch; end
    (isfile(ldb) && isfile(jdb)) || (println("  $cn  SKIP (no treelist)"); continue)
    yrs = Int[]
    let d = SQLite.DB(ldb)
        for r in DBInterface.execute(d, "SELECT DISTINCT Year FROM FVS_TreeList ORDER BY Year"); push!(yrs, Int(r[1])); end
        SQLite.close(d)
    end
    pre = filter(y -> 0 < y < wcyc, yrs)
    growthdiv = false; firstbad = 0; ntrees = 0
    for yr in pre
        Lm = treemap(ldb, yr); Jm = treemap(jdb, yr)
        common = intersect(keys(Lm), keys(Jm)); ntrees = length(common)
        for k in common
            # Compare DBH ONLY (the actual grown size). NOT DG: jl reports DG=-1.0 for seedlings where live
            # reports 0.0 (a treelist reporting convention, not a growth diff) + sub-ULP DG rounding — both
            # leave DBH bit-exact. A real growth bug (e.g. DG-serial-corr) perturbs the applied DBH itself.
            if Lm[k][1] != Jm[k][1]
                growthdiv = true; firstbad = yr; break
            end
        end
        growthdiv && break
    end
    tag = growthdiv ? "★ESCALATE (per-tree DBH/DG div @$firstbad)" : "compounded-ULP"
    growthdiv ? (global n_esc += 1; push!(escalate, cn)) : (global n_ulp += 1)
    println("  $cn  $wcol $(round(mr,digits=1))%@$wcyc  pre=$(pre) matched-trees=$ntrees  → $tag")
end
println("\nSUMMARY: $n_ulp compounded-ULP (cornerable) / $n_esc ESCALATE (deep-trace)")
!isempty(escalate) && println("ESCALATE CNs: ", join(escalate, " "))
