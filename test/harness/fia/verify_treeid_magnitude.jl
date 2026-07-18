include(joinpath("/workspace/FVSjl/test/harness/fia","ledger_fia.jl"))
using SQLite, DBInterface
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
function treemap(db, yr)
    m=Dict{Tuple{String,String,String},Float64}(); isfile(db)||return m
    d=SQLite.DB(db)
    for r in DBInterface.execute(d,"SELECT SpeciesFIA,PtIndex,TreeId,DBH FROM FVS_TreeList WHERE Year=$yr")
        m[(string(r[1]),string(r[2]),string(r[3]))] = r[4]===missing ? 0.0 : Float64(r[4])
    end
    SQLite.close(d); m
end
years(db)=(ys=Int[]; d=SQLite.DB(db); for r in DBInterface.execute(d,"SELECT DISTINCT Year FROM FVS_TreeList ORDER BY Year"); push!(ys,Int(r[1])); end; SQLite.close(d); ys)
THRESH = length(ARGS)>=3 ? parse(Float64,ARGS[3]) : 1.0   # % relative DBH
lines=filter(!isempty, strip.(readlines(ARGS[1]))); out=open(ARGS[2],"w")
for (i,ln) in enumerate(lines)
    p=split(ln); v=p[1]; cn=p[2]; bin=BIN[v]; var=VAR[v]
    dir=mktempdir(); sub=joinpath(dir,"s.db")
    try build_subdb([cn],sub) catch; println(out,"$v $cn BUILDERR"); continue end
    lk=joinpath(dir,"l.key"); write(lk, replace(KEY,"%CN%"=>cn,"%SUB%"=>sub,"%OUT%"=>"out.db"))
    ldb=joinpath(dir,"FVSOut.db"); isfile(ldb)&&rm(ldb)
    run(pipeline(ignorestatus(Cmd(`$bin --keywordfile=$lk`; dir=dir)); stdout=devnull,stderr=devnull))
    jdb=joinpath(dir,"j.db"); isfile(jdb)&&rm(jdb)
    jk=joinpath(dir,"j.key"); write(jk, replace(KEY,"%CN%"=>cn,"%SUB%"=>sub,"%OUT%"=>jdb))
    try FVSjl.run_keyfile(jk; variant=var, faithful=true) catch end
    (isfile(ldb)&&isfile(jdb))||(println(out,"$v $cn NOTREE"); continue)
    ys=years(ldb); worst_rel=0.0; worst_abs=0.0; worst_cyc=0; worst_when_tpamatch=0.0; esc_cyc=0
    for (ci,yr) in enumerate(ys)
        Lm=treemap(ldb,yr); Jm=treemap(jdb,yr)
        tpa_match = keys(Lm)==keys(Jm)
        for k in intersect(keys(Lm),keys(Jm))
            a=Lm[k]; b=Jm[k]; a==b && continue
            rel = 100*abs(a-b)/max(abs(a),1e-6); ab=abs(a-b)
            if rel>worst_rel; worst_rel=rel; worst_abs=ab; worst_cyc=ci; end
            if tpa_match && rel>worst_when_tpamatch; worst_when_tpamatch=rel; end
            if tpa_match && rel>THRESH && esc_cyc==0; esc_cyc=ci; end
        end
    end
    verdict = esc_cyc>0 ? "ESCALATE" : "ULP"
    println(out, "$v $cn $verdict worstRel%=$(round(worst_rel,digits=3)) worstAbs=$(round(worst_abs,digits=4)) worstRel_tpamatch%=$(round(worst_when_tpamatch,digits=3)) esc_cyc=$esc_cyc"); flush(out)
    i%20==0 && println(stderr,"  $i/$(length(lines))")
end
close(out)
