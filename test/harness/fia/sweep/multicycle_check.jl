# Multi-cycle FIA differential — compares jl vs live per-cycle (NOT just cycle-0), across all matching years.
# Closes the blind spot that hid the BM multi-cycle collapse (the western sweeps are all year-0 = pre-growth).
# Reports, per variant: multi-cycle GROWTH bit-exact rate + explicit COLLAPSE detection (jl→0 where live>0).
# Usage: julia --project=. multicycle_check.jl <ncycles> <tmo_s> <dir>=<db> [<dir2>=<db2> ...]
using FVSjl
const COLS = ["TPA","BA","SDI","CCF","TopHt","QMD"]     # growth cols (1:6 of the .sum)
cfg(v) = v=="CR" ? ("/workspace/.crwork/FVScr_clean",FVSjl.CentralRockies()) : v=="IE" ? ("/workspace/.iework/FVSie_clean",FVSjl.InlandEmpire()) :
         v=="EM" ? ("/workspace/.emwork/FVSem_clean",FVSjl.EasternMontana()) : v=="UT" ? ("/workspace/.utwork/FVSut_clean",FVSjl.Utah()) :
         v=="CI" ? ("/workspace/.ciwork/FVSci_clean",FVSjl.CentralIdaho()) : v=="TT" ? ("/workspace/.ttwork/FVStt_clean",FVSjl.Teton()) :
         v=="BM" ? ("/workspace/.bmwork/FVSbm_clean",FVSjl.BlueMountains()) : error(v)
function keytext(cn, db, nc)
    "STDIDENT\n$cn\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\n" *
    "EndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nNUMCYCLE $nc.0\nECHOSUM\nPROCESS\nSTOP\n"
end
function parse_sum(text)
    rows = Tuple{Int,Vector{Float64}}[]
    cols = ((9,14),(15,18),(19,23),(24,27),(28,31),(32,36))
    for ln in split(text,'\n')
        s = rstrip(ln); length(s) < 40 && continue
        y = tryparse(Int, strip(s[1:4])); (y===nothing || y<1000 || y>3000) && continue
        vals = Float64[]; ok=true
        for (a,b) in cols; v=tryparse(Float64,strip(s[a:b])); v===nothing && (ok=false;break); push!(vals,v); end
        ok && length(vals)==6 && push!(rows,(y,vals))
    end
    rows
end
function run_live(bin, cn, db, dir, nc, tmo)
    key=joinpath(dir,"mc.key"); write(key, keytext(cn,db,nc))
    sp=joinpath(dir,"mc.sum"); isfile(sp) && rm(sp)
    try; run(pipeline(`timeout $tmo $bin --keywordfile=$key`; stdout=devnull, stderr=devnull)); catch; end
    isfile(sp) ? read(sp,String) : ""
end
treeless(v) = v[1] < 0.5 && v[2] < 0.5
function check_variant(v, cndb, nc, tmo)
    bin,var = cfg(v); dir=mktempdir()
    treed=0; mcexact=0; collapse=0; drift=0; collapse_cns=String[]; drift_cns=String[]
    for (cn,db) in cndb
        live = run_live(bin,cn,db,dir,nc,tmo); isempty(live) && continue
        L = parse_sum(live); isempty(L) && continue
        treeless(L[1][2]) && continue
        keyf=joinpath(dir,"jl.key"); write(keyf, keytext(cn,db,nc))
        local jl=""; try; r=FVSjl.run_keyfile(keyf; variant=var); jl = r isa AbstractString ? r : string(r); catch; continue; end
        J = parse_sum(jl); isempty(J) && continue
        Jd = Dict(y=>vv for (y,vv) in J)
        treed += 1
        allexact = true; iscollapse = false
        for (y,lv) in L
            haskey(Jd,y) || (allexact=false; continue)
            jv = Jd[y]
            # collapse: live has trees, jl zeroed them
            if lv[1] > 0.5 && jv[1] < 0.5; iscollapse = true; end
            for k in 1:6
                abs(lv[k]-jv[k]) >= 0.5 && (allexact = false)
            end
        end
        if iscollapse; collapse += 1; length(collapse_cns)<6 && push!(collapse_cns, cn)
        elseif allexact; mcexact += 1
        else; drift += 1; length(drift_cns)<6 && push!(drift_cns, cn); end
    end
    (; v, treed, mcexact, collapse, drift, collapse_cns, drift_cns)
end
function main(args)
    nc = parse(Int, args[1]); tmo = parse(Int, args[2])
    pairs=[(split(a,'=')[1],split(a,'=')[2]) for a in args[3:end]]
    order=String[]; byvar=Dict{String,Vector{Tuple{String,String}}}(); seen=Dict{String,Set{String}}()
    for (dir,db) in pairs, f in sort(readdir(dir))
        endswith(f,"_sample.txt") || continue
        vv=uppercase(replace(f,"_sample.txt"=>"")); vv in ("CR","IE","EM","UT","CI","TT","BM") || continue
        haskey(byvar,vv) || (byvar[vv]=Tuple{String,String}[]; seen[vv]=Set{String}(); push!(order,vv))
        for l in eachline(joinpath(dir,f)); s=strip(l); isempty(s)&&continue; cn=split(s,'\t')[1]
            cn in seen[vv] && continue; push!(seen[vv],cn); push!(byvar[vv],(cn,db)); end
    end
    println("===== MULTI-CYCLE FIA DIFFERENTIAL (NUMCYCLE $nc, all cycles vs live) =====")
    for v in order
        r = check_variant(v, byvar[v], nc, tmo)
        pct = r.treed>0 ? round(100*r.mcexact/r.treed,digits=1) : 0.0
        println("$(r.v): treed=$(r.treed) MULTI-CYCLE-exact=$(r.mcexact) ($pct%) | drift=$(r.drift) COLLAPSE=$(r.collapse)")
        r.collapse>0 && println("    ⚠ COLLAPSE stands: ", join(r.collapse_cns, " "))
        r.drift>0 && println("    drift stands: ", join(r.drift_cns, " "))
    end
end
main(ARGS)
