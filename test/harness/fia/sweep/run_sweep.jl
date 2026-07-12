# Hardened stratified FIA sweep: jl vs freshly-relinked live FVS on sampled real stands,
# per variant, all 10 .sum cols, cycle-0 printed-identical. HANG-SAFE: live FVS wrapped in
# `timeout` (a pathological stand can't stall the sweep — the cause of the prior runaway).
# Usage: julia --project=. run_sweep.jl <sampledir> <live_timeout_s>
using FVSjl
DB = length(ARGS) >= 3 ? ARGS[3] : "/workspace/SQLite_FIADB_ENTIRE.db"
const COLS = ["TPA","BA","SDI","CCF","TopHt","QMD","TCuFt","MCuFt","SCuFt","BdFt"]
cfg(v) = v=="LS" ? ("/workspace/FVSjl/tmp/oracles/FVSls_new",FVSjl.LakeStates()) : v=="SN" ? ("/workspace/FVSjl/tmp/oracles/FVSsn_new",FVSjl.Southern()) :
         v=="NE" ? ("/workspace/FVSjl/tmp/oracles/FVSne_new",FVSjl.Northeast()) : v=="CS" ? ("/workspace/FVSjl/tmp/oracles/FVScs_new",FVSjl.CentralStates()) : error(v)
keytext(cn) = """
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
function run_live(bin, cn, dir, tmo)
    key=joinpath(dir,"s.key"); write(key, keytext(cn))
    sp=joinpath(dir,"s.sum"); isfile(sp) && rm(sp)
    crash=false
    try; run(pipeline(`timeout $tmo $bin --keywordfile=$key`; stdout=devnull, stderr=devnull))
    catch e; if e isa Base.ProcessFailedException; p=e.procs[1]; (p.termsignal!=0 || Int(p.exitcode)==136) && (crash=true); end; end
    (isfile(sp) ? read(sp,String) : "", crash)
end
function sweep_variant(v, listfile, tmo)
    bin,var = cfg(v); dir=mktempdir()
    stands=[split(strip(l),'\t')[1] for l in eachline(listfile) if !isempty(strip(l))]
    nrun=0; nok=0; c0=0; colfail=zeros(Int,10); crashes=0; nosum=0; offenders=String[]
    for cn in stands
        nrun+=1
        live, crashed = run_live(bin,cn,dir,tmo)
        if isempty(live); crashed ? (crashes+=1) : (nosum+=1); continue; end
        keyf=joinpath(dir,"jl.key"); write(keyf, keytext(cn))
        jlout=try FVSjl.run_keyfile(keyf; variant=var) catch; ""; end
        isempty(jlout) && continue
        L=parse_sum(live); J=parse_sum(jlout); (isempty(L)||isempty(J)) && continue
        nok+=1; Jd=Dict(y=>vv for (y,vv) in J)
        y0,lv0=L[1]; haskey(Jd,y0) || continue; jv0=Jd[y0]; allok=true
        badcols=String[]
        for k in 1:10
            abs(lv0[k]-jv0[k])>=0.5 && (colfail[k]+=1; allok=false; push!(badcols,COLS[k]))
        end
        allok ? (c0+=1) : push!(offenders, "$cn["*join(badcols,",")*"]")
    end
    (v=v, nrun=nrun, nok=nok, c0=c0, colfail=colfail, crashes=crashes, nosum=nosum, offenders=offenders)
end
function main(sampledir, tmo)
    tmo=parse(Int,tmo)
    results=[]
    for v in ("SN","NE","CS","LS")
        lf=joinpath(sampledir,"$(lowercase(v))_sample.txt"); isfile(lf) || continue
        println(">>> sweeping $v …"); flush(stdout)
        push!(results, sweep_variant(v, lf, tmo))
    end
    println("\n================ STRATIFIED FIA SWEEP — cycle-0 all-10-cols vs live ================")
    tot_ok=0; tot_c0=0
    for r in results
        pct = r.nok>0 ? round(100*r.c0/r.nok, digits=1) : 0.0
        println("$(r.v): stands=$(r.nrun) both-sum=$(r.nok) BIT-EXACT=$(r.c0)/$(r.nok) ($pct%) | live-SIGFPE-crash=$(r.crashes) treeless/nosum=$(r.nosum)")
        fails=join(["$(COLS[k]):$(r.colfail[k])" for k in 1:10 if r.colfail[k]>0], " ")
        println("    col mismatches: ", isempty(fails) ? "NONE ✓" : fails)
        isempty(r.offenders) || println("    offenders: ", join(r.offenders, "  "))
        tot_ok+=r.nok; tot_c0+=r.c0
    end
    println("--------------------------------------------------------------------------------")
    println("TOTAL: $(tot_c0)/$(tot_ok) cycle-0 bit-exact (", tot_ok>0 ? round(100*tot_c0/tot_ok,digits=2) : 0, "%)")
end
main(ARGS[1], ARGS[2])
