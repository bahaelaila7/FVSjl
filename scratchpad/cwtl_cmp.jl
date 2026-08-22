using SQLite, DBInterface
function cw(p, yr)
    db = SQLite.DB(p)
    ("FVS_TreeList" in [t.name for t in SQLite.tables(db)]) || return Dict{String,Tuple{Float64,String,Float64}}()
    d = Dict{String,Tuple{Float64,String,Float64}}()
    for r in DBInterface.execute(db, "SELECT TreeId,CrWidth,SpeciesFVS,DBH FROM FVS_TreeList WHERE Year=$yr")
        d[string(r[1])] = (Float64(r[2]), string(r[3]), Float64(r[4]))
    end
    d
end
function main()
    op, jp, yr = ARGS[1], ARGS[2], parse(Int, ARGS[3])
    o = cw(op, yr); j = cw(jp, yr)
    if isempty(o) || isempty(j)
        println("  EMPTY: oracle=$(length(o)) jl=$(length(j))"); return
    end
    sh = intersect(keys(o), keys(j)); w = 0.0; n = 0; worst = ""
    for id in sh
        dd = abs(o[id][1] - j[id][1]); dd > 0.05 && (n += 1); dd > w && (w = dd; worst = id)
    end
    print("  trees o=$(length(o)) j=$(length(j)) shared=$(length(sh)) | CrWidth worst|Δ|=$(round(w,digits=4)) diffs=$n")
    n > 0 && print(" [worst id=$worst sp=$(o[worst][2]) o=$(round(o[worst][1],digits=3)) j=$(round(j[worst][1],digits=3)) DBH=$(round(o[worst][3],digits=2))]")
    println()
    c = 0
    for id in sort(collect(sh), by=x->parse(Int,x))
        dd = abs(o[id][1] - j[id][1])
        if dd > 0.05 && c < 5
            println("    diff id=$id sp=$(o[id][2]) o=$(round(o[id][1],digits=3)) j=$(round(j[id][1],digits=3)) DBH=$(round(o[id][3],digits=2))"); c += 1
        end
    end
end
main()
