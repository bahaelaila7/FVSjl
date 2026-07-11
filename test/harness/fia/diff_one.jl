# diff_one.jl — per-cycle side-by-side live-vs-jl .sum dump for ONE stand under a regime (Pillar-3 dig helper).
# Usage: julia --project=. test/harness/fia/diff_one.jl <CN> <SN|NE|CS|LS> [regime]
using FVSjl
import SQLite, DBInterface
include(joinpath(@__DIR__, "manage_fia.jl"))   # reuse keytext/regime_block/build_subdb/parse_sum/BIN/VAR

function dump(cn, v, regime)
    bin = BIN[v]; var = VAR[v]; dir = mktempdir()
    db = joinpath(dir,"sub.db"); build_subdb([cn], db)
    live = run_live(bin, cn, db, regime, dir)
    keyf = joinpath(dir,"jl.key"); write(keyf, keytext(cn, db, regime))
    jlout = FVSjl.run_keyfile(keyf; variant=var)
    L = parse_sum(live); J = Dict(y=>vv for (y,vv) in parse_sum(jlout))
    cols = ["TPA","BA","SDI","CCF","TopHt","QMD"]
    println("CN=$cn  regime=$regime  variant=$v")
    println("year  ", join([rpad(c,20) for c in cols]))
    for (y,lv) in L
        haskey(J,y) || (println("$y  (jl missing)"); continue)
        jv = J[y]
        cells = String[]
        for k in 1:6
            mark = lv[k]==jv[k] ? " " : "*"
            push!(cells, rpad(string(mark, lv[k], "/", jv[k]), 20))
        end
        println("$y  ", join(cells))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    dump(ARGS[1], ARGS[2], length(ARGS)>=3 ? ARGS[3] : "thinbba")
end
