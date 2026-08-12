using FVSjl, SQLite
db="/workspace/SQLite_FIADB_ENTIRE.db"; var=FVSjl.EasternMontana(); bin="/workspace/.emwork/FVSem_clean"
# full per-cycle TPA/BA from .sum (all year rows)
rows(t)=[(parse(Int,split(strip(l))[1]),parse(Int,split(strip(l))[3]),parse(Float64,split(strip(l))[4])) for l in split(t,'\n') if occursin(r"^\s*\d{4}\s",l) && length(split(strip(l)))>=4]
key(cn)="STDIDENT\n$cn\nNUMCYCLE 4.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
live(cn)=(d=mktempdir();kf=joinpath(d,"k.key");write(kf,key(cn));sp=joinpath(d,"k.sum");cd(d);try;run(pipeline(`timeout 90 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull));catch;end;isfile(sp) ? read(sp,String) : "")
jlrun(cn)=(d=mktempdir();kf=joinpath(d,"k.key");write(kf,key(cn));FVSjl.run_keyfile(kf;variant=var))
for cn in ["5332701010661"]
  nrec=0
  for r in DBInterface.execute(SQLite.DB(db),"SELECT COUNT(*) c FROM FVS_TREEINIT_COND WHERE STAND_CN=? AND DIAMETER>0",[cn]); nrec=r.c; end
  println("=== stand $cn (start tree-records=$nrec) ===")
  lt=live(cn); jr=rows(jlrun(cn)); lr=rows(lt)
  println("  year |   jl TPA/BA   |  live TPA/BA")
  for (y,jT,jB) in jr
    m=filter(x->x[1]==y,lr); isempty(m)&&continue; (_,lT,lB)=m[1]
    println("  $y |  $jT / $jB  |  $lT / $lB")
  end
end
