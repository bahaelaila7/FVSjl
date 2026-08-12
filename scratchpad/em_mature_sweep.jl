using FVSjl, SQLite
db="/workspace/SQLite_FIADB_ENTIRE.db"; DB=SQLite.DB(db)
var=FVSjl.EasternMontana(); bin="/workspace/.emwork/FVSem_clean"
fb(t)=(l="";for ln in split(t,'\n');occursin(r"^\s*\d{4}",ln)&&(l=strip(ln));end;isempty(l)&&return(0,0.0);f=split(l);length(f)<4&&return(0,0.0);(parse(Int,f[3]),parse(Float64,f[4])))
key(cn)="STDIDENT\n$cn\nNUMCYCLE 5.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
live(cn)=(d=mktempdir();kf=joinpath(d,"k.key");write(kf,key(cn));sp=joinpath(d,"k.sum");cd(d);try;run(pipeline(`timeout 90 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull));catch;end;isfile(sp) ? read(sp,String) : "")
jlrun(cn)=(d=mktempdir();kf=joinpath(d,"k.key");write(kf,key(cn));FVSjl.run_keyfile(kf;variant=var))
# mature EM stands, moderate density (avoid the SDI>1500 instability), few records for speed
q="""SELECT t.STAND_CN cn FROM FVS_TREEINIT_COND t JOIN FVS_STANDINIT_COND s ON s.STAND_CN=t.STAND_CN
     WHERE s.VARIANT='EM' AND t.DIAMETER>0 GROUP BY t.STAND_CN
     HAVING SUM(0.005454154*t.DIAMETER*t.DIAMETER*t.TREE_COUNT) BETWEEN 120 AND 350 AND AVG(t.DIAMETER)>7 AND COUNT(*) BETWEEN 3 AND 25
     ORDER BY t.STAND_CN LIMIT 10"""
stands=[r.cn for r in DBInterface.execute(DB,q)]
println("EM mature stands: ", length(stands))
ne=0;corn=0;dv=0;n=0
for cn in stands
  lt=live(cn); isempty(lt) && (println("  NOLIVE/CRASH cn=$cn"); continue)
  jT,jB=fb(jlrun(cn)); lT,lB=fb(lt); (jB==0&&lB==0)&&continue
  global n+=1; d=lB>0 ? abs(jB-lB)/lB : 0.0
  sign = jB>lB ? "+" : "-"
  if jB==lB&&jT==lT; global ne+=1
  elseif d<=0.03; global corn+=1
  else; global dv+=1; println("  DIV cn=$cn jl BA=$jB TPA=$jT | live BA=$lB TPA=$lT ($sign$(round(d*100,digits=1))%)"); end
end
println("EM mature: n=$n bit-exact=$ne cornered=$corn diverged=$dv")
