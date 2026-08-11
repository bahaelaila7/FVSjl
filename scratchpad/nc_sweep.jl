using FVSjl, SQLite
db="/workspace/SQLite_FIADB_ENTIRE.db"
DB=SQLite.DB(db)
var=FVSjl.Klamath(); bin="/workspace/.ncwork/FVSnc_clean"
function finalba(txt::String)
    last=""
    for ln in split(txt,'\n'); occursin(r"^\s*\d{4}",ln) && (last=strip(ln)); end
    isempty(last) && return (0,0.0)
    f=split(last); length(f)<4 && return (0,0.0); (parse(Int,f[3]), parse(Float64,f[4]))
end
function runlive(bin,cn)
    kt="STDIDENT\n$cn\nNUMCYCLE 5.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
    d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,kt); sp=joinpath(d,"k.sum"); cd(d)
    try; run(pipeline(`timeout 120 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull)); catch; end
    isfile(sp) ? read(sp,String) : ""
end
q="SELECT t.STAND_CN FROM FVS_STANDINIT_COND s JOIN (SELECT DISTINCT STAND_CN FROM FVS_TREEINIT_COND) t ON s.STAND_CN=t.STAND_CN WHERE s.VARIANT='NC' LIMIT 12"
stands=[r.STAND_CN for r in DBInterface.execute(DB,q)]
println("NC stands found: ", length(stands))
ne=0; corn=0; div=0; crash=0; n=0; nolive=0
for cn in stands
    kt="STDIDENT\n$cn\nNUMCYCLE 5.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
    d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,kt)
    jt=""; try; jt=FVSjl.run_keyfile(kf; variant=var); catch e; global crash+=1; continue; end
    lt=runlive(bin,cn); isempty(lt) && (global nolive+=1; continue)
    jT,jB=finalba(jt); lT,lB=finalba(lt)
    (jB==0 && lB==0) && continue
    global n+=1
    dba = lB>0 ? abs(jB-lB)/lB : 0.0
    if jB==lB && jT==lT; global ne+=1
    elseif dba<=0.03; global corn+=1
    else; global div+=1; println("  DIV cn=$cn jl BA=$jB TPA=$jT | live BA=$lB TPA=$lT (Δ$(round(dba*100,digits=1))%)"); end
end
println("NC: n=$n bit-exact=$ne cornered(≤3%)=$corn diverged(>3%)=$div jl-crash=$crash no-live=$nolive")
