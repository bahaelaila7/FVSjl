# Enhanced western FIA sweep — TREED-only cycle-0 all-10-col vs freshly-relinked live FVS.
# Beyond run_sweep.jl: (1) EXCLUDES treeless stands from the denominator — western live FVS emits an
# all-zero .sum for treeless FIA conditions (unlike eastern, which emits none), so counting them dilutes
# the rate with trivially-exact zeros; (2) splits GROWTH cols (TPA/BA/SDI/CCF/TopHt/QMD) vs VOLUME cols
# (TCuFt/MCuFt/SCuFt/BdFt) exactness; (3) buckets volume-delta MAGNITUDE (<2% / 2-10% / >10%) to separate
# cornered thresholds from real divergence; (4) captures jl errors (crash/KeyError) with a repro STAND_CN.
# Merges multiple sample dirs (batch-1 + augmentation), deduped by STAND_CN.
# Usage: julia --project=. run_sweep_western.jl <tmo_s> <dir1>=<db1> [<dir2>=<db2> ...]
using FVSjl
const COLS = ["TPA","BA","SDI","CCF","TopHt","QMD","TCuFt","MCuFt","SCuFt","BdFt"]
const GROWTH = 1:6; const VOL = 7:10
cfg(v) = v=="CR" ? ("/workspace/.crwork/FVScr_clean",FVSjl.CentralRockies()) : v=="IE" ? ("/workspace/.iework/FVSie_clean",FVSjl.InlandEmpire()) :
         v=="EM" ? ("/workspace/.emwork/FVSem_clean",FVSjl.EasternMontana()) : v=="UT" ? ("/workspace/.utwork/FVSut_clean",FVSjl.Utah()) :
         v=="CI" ? ("/workspace/.ciwork/FVSci_clean",FVSjl.CentralIdaho()) : v=="TT" ? ("/workspace/.ttwork/FVStt_clean",FVSjl.Teton()) : error(v)
keytext(cn,db) = "STDIDENT\n$cn\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nNUMCYCLE         3.0\nECHOSUM\nPROCESS\nSTOP\n"
function parse_sum(text)
    rows = Tuple{Int,Vector{Float64}}[]
    cols = ((9,14),(15,18),(19,23),(24,27),(28,31),(32,36),(37,42),(43,48),(49,54),(55,60))
    for ln in split(text,'\n')
        s = rstrip(ln); length(s) < 60 && continue
        y = tryparse(Int, strip(s[1:4])); (y===nothing || y<1000 || y>3000) && continue
        vals = Float64[]; ok=true
        for (a,b) in cols; v=tryparse(Float64,strip(s[a:b])); v===nothing && (ok=false;break); push!(vals,v); end
        ok && length(vals)==10 && push!(rows,(y,vals))
    end
    rows
end
function run_live(bin, cn, db, dir, tmo)
    key=joinpath(dir,"s.key"); write(key, keytext(cn,db))
    sp=joinpath(dir,"s.sum"); isfile(sp) && rm(sp)
    crash=false
    try; run(pipeline(`timeout $tmo $bin --keywordfile=$key`; stdout=devnull, stderr=devnull))
    catch e; if e isa Base.ProcessFailedException; p=e.procs[1]; (p.termsignal!=0 || Int(p.exitcode)==136) && (crash=true); end; end
    (isfile(sp) ? read(sp,String) : "", crash)
end
treeless(lv0) = lv0[1] < 0.5 && lv0[2] < 0.5    # year-0 TPA & BA both ~0
reldelta(l,j) = abs(l-j)/max(abs(l),1.0)
function sweep_variant(v, cndb, tmo)
    bin,var = cfg(v); dir=mktempdir()
    treed=0; c0=0; gexact=0; vexact=0; colfail=zeros(Int,10)
    crashes=0; treelessn=0; jlerr=String[]; offenders=String[]
    volbucket=[0,0,0]   # <2%, 2-10%, >10% (max vol reldelta among the divergent vol cols)
    for (cn,db) in cndb
        live, crashed = run_live(bin,cn,db,dir,tmo)
        if isempty(live); crashed && (crashes+=1); continue; end
        L=parse_sum(live); isempty(L) && continue
        y0,lv0=L[1]
        if treeless(lv0); treelessn+=1; continue; end       # exclude treeless from denominator
        keyf=joinpath(dir,"jl.key"); write(keyf, keytext(cn,db))
        jlout=""; err=nothing
        try; jlout=FVSjl.run_keyfile(keyf; variant=var); catch e; err=e; end
        if err!==nothing; length(jlerr)<3 && push!(jlerr, "$cn: "*first(split(sprint(showerror,err),'\n'))); continue; end
        isempty(jlout) && continue
        J=parse_sum(jlout); isempty(J) && continue
        Jd=Dict(y=>vv for (y,vv) in J); haskey(Jd,y0) || continue; jv0=Jd[y0]
        treed+=1
        badcols=String[]; gbad=false; vbad=false; maxvrel=0.0
        for k in 1:10
            if abs(lv0[k]-jv0[k])>=0.5
                colfail[k]+=1; push!(badcols,COLS[k])
                k in GROWTH ? (gbad=true) : (vbad=true; maxvrel=max(maxvrel, reldelta(lv0[k],jv0[k])))
            end
        end
        gbad || (gexact+=1); vbad || (vexact+=1)
        if isempty(badcols); c0+=1
        else
            push!(offenders, "$cn["*join(badcols,",")*(vbad ? "|Δv="*string(round(100*maxvrel,digits=1))*"%" : "")*"]")
            vbad && (maxvrel<0.02 ? (volbucket[1]+=1) : maxvrel<0.10 ? (volbucket[2]+=1) : (volbucket[3]+=1))
        end
    end
    (; v, treed, c0, gexact, vexact, colfail, crashes, treelessn, jlerr, offenders, volbucket)
end
function main(args)
    tmo=parse(Int,args[1])
    pairs=[(split(a,'=')[1],split(a,'=')[2]) for a in args[2:end]]  # (dir,db)
    # collect per-variant unique CNs with their db (first dir wins on dup)
    order=String[]; byvar=Dict{String,Vector{Tuple{String,String}}}(); seen=Dict{String,Set{String}}()
    for (dir,db) in pairs
        for f in sort(readdir(dir))
            endswith(f,"_sample.txt") || continue
            vv=uppercase(replace(f,"_sample.txt"=>""))
            vv in ("CR","IE","EM","UT","CI","TT") || continue
            haskey(byvar,vv) || (byvar[vv]=Tuple{String,String}[]; seen[vv]=Set{String}(); push!(order,vv))
            for l in eachline(joinpath(dir,f)); s=strip(l); isempty(s)&&continue; cn=split(s,'\t')[1]
                cn in seen[vv] && continue; push!(seen[vv],cn); push!(byvar[vv],(cn,db)); end
        end
    end
    results=[]
    for v in order
        println(">>> sweeping $v ($(length(byvar[v])) unique stands) …"); flush(stdout)
        push!(results, sweep_variant(v, byvar[v], tmo))
    end
    println("\n===== WESTERN FIA SWEEP — TREED-only cycle-0 all-10-col vs live (treeless excluded) =====")
    tt=0; tc=0
    for r in results
        pct = r.treed>0 ? round(100*r.c0/r.treed,digits=1) : 0.0
        gp = r.treed>0 ? round(100*r.gexact/r.treed,digits=1) : 0.0
        vp = r.treed>0 ? round(100*r.vexact/r.treed,digits=1) : 0.0
        println("$(r.v): treed=$(r.treed) ALL10=$(r.c0)/$(r.treed) ($pct%) | GROWTH-exact=$(r.gexact) ($gp%) VOL-exact=$(r.vexact) ($vp%) | live-crash=$(r.crashes) treeless-excl=$(r.treelessn)")
        fails=join(["$(COLS[k]):$(r.colfail[k])" for k in 1:10 if r.colfail[k]>0], " ")
        println("    col mismatches: ", isempty(fails) ? "NONE ✓" : fails, "   | vol-Δ buckets [<2%,2-10%,>10%]=", r.volbucket)
        isempty(r.jlerr) || println("    jl ERRORS: ", join(r.jlerr, " || "))
        isempty(r.offenders) || println("    offenders: ", join(r.offenders, "  "))
        tt+=r.treed; tc+=r.c0
    end
    println("-"^80)
    println("TOTAL TREED: $tc/$tt all-10-col bit-exact (", tt>0 ? round(100*tc/tt,digits=2) : 0, "%)")
end
main(ARGS)
