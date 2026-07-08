# manage_fia.jl — Pillar-3: management-scenario compatibility on real FIA plots.
# Same driver as validate_fia.jl but injects a silvicultural KEYWORD BLOCK before PROCESS, so the
# projection runs UNDER management. Diffs the 6-col .sum (TPA/BA/SDI/CCF/TopHt/QMD) live-FVS vs FVSjl
# across all cycles — validating that the keyword removes the same trees and the residual stand projects
# identically. Regime is selectable so thin/plant/fire can each be swept.
#
# Usage: FIA_DB=<subdb> julia --project=. test/harness/fia/manage_fia.jl <standlist> <SN|NE|CS|LS> [regime]
#   regime ∈ {thinbba, thinbta, thindbh}  (default thinbba)  — cycle-2 thin, universal (scales to stand).

using FVSjl
const DB = get(ENV, "FIA_DB", "/workspace/SQLite_FIADB_ENTIRE.db")
const BIN = Dict("SN"=>"/tmp/FVSsn_new","NE"=>"/tmp/FVSne_new","CS"=>"/tmp/FVScs_new","LS"=>"/tmp/FVSls_new")
const VAR = Dict("SN"=>FVSjl.Southern(),"NE"=>FVSjl.Northeast(),"CS"=>FVSjl.CentralStates(),"LS"=>FVSjl.LakeStates())

# Regime keyword blocks. Scheduled by CYCLE (field 1 = 2 ⇒ cycle 2, aligns across stands regardless of
# inventory year — calendar years don't align and live/jl interpret them differently). Residual chosen low
# so the thin FIRES on most real stands (no-op stands still validate).
# FVS reads keyword records FIXED-FORMAT (A10 keyword, then F10.0 fields) — columns MUST align to 10-char
# boundaries or the values scramble. kwrec() enforces that: keyword left-justified in 10, each field
# right-justified in 10.
kwrec(kw, fields...) = rpad(kw, 10) * join(lpad(string(f), 10) for f in fields)
regime_block(r) =
    r == "thinbta" ? kwrec("THINBTA", "2.0", "40.0") :               # thin from above to residual BA 40
    r == "thindbh" ? kwrec("THINDBH", "2.0", "0.0", "99.0", "0.5") : # cut 50% across all DBH
    r == "simfire" ? "FMIn\n" * kwrec("SIMFIRE", "2.0", "10.00", "1", "50.0") * "\nEnd" : # prescribed fire cyc2 (FFE)
    r == "plant"   ? "ESTAB\n" * kwrec("PLANT", "2.0", "3", "400") * "\nEnd" :          # plant sp3 400tpa cyc2 (ESTAB/regen)
    r == "salvage" ? kwrec("SALVAGE", "2.0", "0.0", "999.0", "0.9") :                   # salvage 90% of dead cyc2
                     kwrec("THINBBA", "2.0", "40.0")                  # thin from below to residual BA 40 (default)

keytext(cn, regime) = """
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
$(regime_block(regime))
ECHOSUM
PROCESS
STOP
"""

function parse_sum(text)
    rows = Tuple{Int,Vector{Float64}}[]
    for ln in split(text, '\n')
        f = split(strip(ln)); length(f) < 8 && continue
        y = tryparse(Int, f[1]); (y===nothing || y<1000 || y>3000) && continue
        vals = Float64[]; ok = true
        for i in 3:8; v = tryparse(Float64, f[i]); v===nothing && (ok=false; break); push!(vals, v); end
        ok && push!(rows, (y, vals))
    end
    rows
end

function run_live(bin, cn, regime, dir)
    key = joinpath(dir, "s.key"); write(key, keytext(cn, regime))
    for f in ("s.sum","s.out"); fp=joinpath(dir,f); isfile(fp) && rm(fp); end
    try; run(pipeline(`$bin --keywordfile=$key`; stdout=devnull, stderr=devnull)); catch; end
    sp = joinpath(dir,"s.sum"); isfile(sp) ? read(sp,String) : ""
end

function main(listfile, v, regime)
    bin = BIN[v]; var = VAR[v]; dir = mktempdir()
    stands = [split(strip(l), '\t')[1] for l in eachline(listfile) if !isempty(strip(l))]
    n_run=0; n_ok=0; n_pass=0; n_thinned=0; failures=String[]
    worst = Tuple{Float64,String}[]
    fired_pass=0; fired_fail=0; noop_pass=0; noop_fail=0   # split pass rate by whether the thin FIRED
    for cn in stands
        n_run += 1
        print(stderr, "[$n_run/$(length(stands))] $cn live..."); flush(stderr)
        live = run_live(bin, cn, regime, dir)
        print(stderr, isempty(live) ? "NOSUM " : "ok jl..."); flush(stderr)
        isempty(live) && (println(stderr); continue)
        keyf = joinpath(dir,"jl.key"); write(keyf, keytext(cn, regime))
        jlout = try FVSjl.run_keyfile(keyf; variant=var) catch e; println(stderr, "JLERR:$e"); ""; end
        println(stderr, isempty(jlout) ? "nojl" : "ok"); flush(stderr)
        isempty(jlout) && continue
        L = parse_sum(live); J = parse_sum(jlout)
        (isempty(L) || isempty(J)) && continue
        n_ok += 1
        Jd = Dict(y=>vv for (y,vv) in J)
        # detect whether the thin actually fired: BA drops from cycle 1 → cycle 2 on the LIVE run
        fired = length(L) >= 3 && L[3][2][2] < L[2][2][2] - 0.5
        fired && (n_thinned += 1)
        maxrel = 0.0
        for (y, lv) in L
            haskey(Jd, y) || continue
            jv = Jd[y]
            for k in 1:6
                r = lv[k]==0 ? (jv[k]==0 ? 0.0 : 1.0) : abs(lv[k]-jv[k])/abs(lv[k])
                maxrel = max(maxrel, r)
            end
        end
        push!(worst, (maxrel, cn))
        stand_pass = length(J)==length(L) &&
                     all(haskey(Jd,y) && all(lv[k]==Jd[y][k] for k in 1:6) for (y,lv) in L)
        stand_pass ? (n_pass += 1) : push!(failures, cn)
        fired ? (stand_pass ? (fired_pass += 1) : (fired_fail += 1)) :
                (stand_pass ? (noop_pass  += 1) : (noop_fail  += 1))
    end
    if haskey(ENV,"FIA_FAILOUT") && !isempty(failures)
        mag = Dict(s=>r for (r,s) in worst)
        open(ENV["FIA_FAILOUT"],"w") do io
            for cn in failures; println(io, cn, '\t', v, '\t', round(get(mag,cn,0.0)*100,digits=1)); end
        end
    end
    println("\n===== MANAGEMENT ($regime): $v =====")
    println("stands run=$n_run  both-sum=$n_ok  thin-fired(live)=$n_thinned")
    println("BIT-EXACT (all cycles, 6 cols): $n_pass / $n_ok    FAIL: $(length(failures))")
    println("  thin FIRED: bit-exact $fired_pass / $(fired_pass+fired_fail)   |   thin NO-OP: bit-exact $noop_pass / $(noop_pass+noop_fail) (= growth-only, cf. Pillar-2)")
    buckets = ["<1%"=>0,"1-2%"=>0,"2-5%"=>0,"5-10%"=>0,">10%"=>0]
    for (r,_) in worst; p=r*100; k = p<1 ? 1 : p<2 ? 2 : p<5 ? 3 : p<10 ? 4 : 5; buckets[k]=buckets[k].first=>buckets[k].second+1; end
    println("Worst-rel-diff histogram: ", join(["$(b.first):$(b.second)" for b in buckets], "  "))
    sort!(worst, rev=true)
    println("Worst 8: ", join(["$(round(r*100,digits=1))%:$s" for (r,s) in worst[1:min(8,end)]], "  "))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: manage_fia.jl <standlist> <variant> [regime]")
    main(ARGS[1], ARGS[2], length(ARGS) >= 3 ? ARGS[3] : "thinbba")
end
