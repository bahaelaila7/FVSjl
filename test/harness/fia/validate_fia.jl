# FIA validation harness: FVSjl (native reader) vs live FVS, per stand, diff .sum.
using FVSjl

const DB = "/workspace/SQLite_FIADB_ENTIRE.db"

# variant → (live binary, FVSjl variant)
variant_cfg(v) = v == "LS" ? ("/tmp/FVSls_new", FVSjl.LakeStates()) :
                 v == "SN" ? ("/tmp/FVSsn_new", FVSjl.Southern()) :
                 v == "NE" ? ("/tmp/FVSne_new", FVSjl.Northeast()) :
                 v == "CS" ? ("/tmp/FVScs_new", FVSjl.CentralStates()) : error("variant $v")

keyfile_text(cn) = """
STDIDENT
$cn
DATABASE
DSNin
$DB
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE         5.0
ECHOSUM
PROCESS
STOP
"""

# parse .sum data rows → Vector of (year, [tpa,ba,sdi,ccf,topht,qmd])
function parse_sum(text)
    rows = Tuple{Int,Vector{Float64}}[]
    for ln in split(text, '\n')
        f = split(strip(ln))
        (length(f) < 8) && continue
        y = tryparse(Int, f[1]); (y === nothing || y < 1000 || y > 3000) && continue
        vals = Float64[]
        ok = true
        for i in 3:8
            v = tryparse(Float64, f[i]); v === nothing && (ok = false; break); push!(vals, v)
        end
        ok && push!(rows, (y, vals))
    end
    return rows
end

function run_live(bin, cn, dir)
    key = joinpath(dir, "s.key")
    write(key, keyfile_text(cn))
    for f in ("s.sum","s.out"); fp=joinpath(dir,f); isfile(fp) && rm(fp); end
    try
        run(pipeline(`$bin --keywordfile=$key`; stdout=devnull, stderr=devnull))
    catch; end
    sp = joinpath(dir, "s.sum")
    isfile(sp) ? read(sp, String) : ""
end

function main(listfile, variant)
    bin, var = variant_cfg(variant)
    stands = [split(strip(l), '\t') for l in eachline(listfile) if !isempty(strip(l))]
    dir = mktempdir()
    # aggregate: per-cycle abs-rel-diff of the 6 cols; cycle-0 exact count
    n_ok = 0; n_run = 0; c0_exact = 0
    col_names = ["TPA","BA","SDI","CCF","TopHt","QMD"]
    # accumulate mean |rel diff| by cycle index
    sumrel = Dict{Int,Vector{Float64}}(); cnt = Dict{Int,Int}()
    worst = Tuple{Float64,String}[]
    for (cn, tag) in stands
        n_run += 1
        print(stderr, "[$n_run/$(length(stands))] $cn[$tag] live..."); flush(stderr)
        live = run_live(bin, cn, dir)
        print(stderr, isempty(live) ? "NOSUM " : "ok jl..."); flush(stderr)
        isempty(live) && (println(stderr); continue)
        keyf = joinpath(dir, "jl.key"); write(keyf, keyfile_text(cn))
        jlout = try FVSjl.run_keyfile(keyf; variant=var) catch e; println(stderr, "JLERR:$e"); ""; end
        println(stderr, isempty(jlout) ? "nojl" : "ok"); flush(stderr)
        (isempty(jlout)) && continue
        L = parse_sum(live); J = parse_sum(jlout)
        (isempty(L) || isempty(J)) && continue
        n_ok += 1
        # align by year
        Jd = Dict(y=>v for (y,v) in J)
        # cycle-0 exact?
        y0, lv0 = L[1]
        if haskey(Jd, y0)
            jv0 = Jd[y0]
            if all(abs.(lv0 .- jv0) .< 0.5)   # .sum is integer-rounded (QMD 1 dp) ⇒ <0.5 = same printed
                c0_exact += 1
            end
        end
        maxrel = 0.0
        for (ci, (y, lv)) in enumerate(L)
            haskey(Jd, y) || continue
            jv = Jd[y]
            rel = [lv[k]==0 ? (jv[k]==0 ? 0.0 : 1.0) : abs(lv[k]-jv[k])/abs(lv[k]) for k in 1:6]
            get!(sumrel, ci, zeros(6)); sumrel[ci] .+= rel; cnt[ci] = get(cnt,ci,0)+1
            maxrel = max(maxrel, maximum(rel))
        end
        push!(worst, (maxrel, "$cn[$tag]"))
    end
    println("\n===== FIA validation: $variant =====")
    println("stands run=$n_run  both-produced-sum=$n_ok  cycle0-printed-identical=$c0_exact / $n_ok")
    println("\nMean |rel diff| by cycle (over $n_ok stands):")
    println("  cyc  ", join(lpad.(col_names,7), " "))
    for ci in sort(collect(keys(sumrel)))
        m = sumrel[ci] ./ cnt[ci]
        println("  ", lpad(ci-1,3), "  ", join([lpad(string(round(m[k]*100,digits=2))*"%",7) for k in 1:6], " "))
    end
    sort!(worst, rev=true)
    println("\nWorst 8 stands by max rel diff (any cycle/col):")
    for (r,s) in worst[1:min(8,end)]; println("  ", round(r*100,digits=1), "%  ", s); end
end

main("/tmp/fia_val/ls_stands.txt", "LS")
