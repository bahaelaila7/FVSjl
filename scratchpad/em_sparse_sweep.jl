using FVSjl, SQLite
db="/workspace/SQLite_FIADB_ENTIRE.db"; DB=SQLite.DB(db)
var=FVSjl.EasternMontana(); bin="/workspace/.emwork/FVSem_clean"
fb(t)=(l="";for ln in split(t,'\n');occursin(r"^\s*\d{4}",ln)&&(l=strip(ln));end;isempty(l)&&return(0,0.0);f=split(l);length(f)<4&&return(0,0.0);(parse(Int,f[3]),parse(Float64,f[4])))
key(cn)="STDIDENT\n$cn\nNUMCYCLE 5.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
live(cn)=(d=mktempdir();kf=joinpath(d,"k.key");write(kf,key(cn));sp=joinpath(d,"k.sum");cd(d);try;run(pipeline(`timeout 60 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull));catch;end;isfile(sp) ? read(sp,String) : "")
jlrun(cn)=(d=mktempdir();kf=joinpath(d,"k.key");write(kf,key(cn));FVSjl.run_keyfile(kf;variant=var))
# SPARSE mature EM: large avg DBH (>13"), LOW BA (sparse), few records → the low-competition regime
q="""SELECT t.STAND_CN cn FROM FVS_TREEINIT_COND t JOIN FVS_STANDINIT_COND s ON s.STAND_CN=t.STAND_CN
     WHERE s.VARIANT='EM' AND t.DIAMETER>0 GROUP BY t.STAND_CN
     HAVING SUM(0.005454154*t.DIAMETER*t.DIAMETER*t.TREE_COUNT) BETWEEN 30 AND 90 AND AVG(t.DIAMETER)>13 AND COUNT(*) BETWEEN 2 AND 12
     ORDER BY t.STAND_CN LIMIT 6"""
stands=[r.cn for r in DBInterface.execute(DB,q)]
println("EM SPARSE-mature stands: ", length(stands))
for cn in stands
  lt=live(cn); isempty(lt) && (println("  NOLIVE/CRASH cn=$cn"); continue)
  jT,jB=fb(jlrun(cn)); lT,lB=fb(lt); (jB==0&&lB==0)&&continue
  d=lB>0 ? round((jB-lB)/lB*100,digits=1) : 0.0
  println("  cn=$cn jl BA=$jB TPA=$jT | live BA=$lB TPA=$lT (Δ$d%)")
end
println("done")
