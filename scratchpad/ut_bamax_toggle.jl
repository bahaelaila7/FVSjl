# Toggle test: run the densest UT stand(s) with the BAMAX cap ACTIVE, and also with it
# NEUTRALIZED (env FVSJL_NO_BAMAX=1 → the cap code checks this and skips), vs live.
# Proves the cap fires: no-cap should over-grow (diverge high), cap should match live.
# NOTE requires a temporary env hook in utah/mortality.jl (added for this test only).
using FVSjl, SQLite
db="/workspace/SQLite_FIADB_ENTIRE.db"; DB=SQLite.DB(db)
var=FVSjl.Utah(); bin="/workspace/.utwork/FVSut_clean"
function finalba(txt)
    last=""; for ln in split(txt,'\n'); occursin(r"^\s*\d{4}",ln) && (last=strip(ln)); end
    isempty(last) && return (0,0.0); f=split(last); length(f)<4 && return (0,0.0); (parse(Int,f[3]), parse(Float64,f[4]))
end
function runlive(cn)
    kt="STDIDENT\n$cn\nNUMCYCLE 5.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
    d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,kt); sp=joinpath(d,"k.sum"); cd(d)
    try; run(pipeline(`timeout 120 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull)); catch; end
    isfile(sp) ? read(sp,String) : ""
end
function runjl(cn)
    kt="STDIDENT\n$cn\nNUMCYCLE 5.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
    d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,kt); FVSjl.run_keyfile(kf; variant=var)
end
# stand CNs passed as ARGS (the divergent/densest ones from the sweep)
for cn in ARGS
    lt=runlive(cn); _,lB=finalba(lt)
    ENV["FVSJL_NO_BAMAX"]="1"; _,jB_off=finalba(runjl(cn))
    delete!(ENV,"FVSJL_NO_BAMAX"); _,jB_on=finalba(runjl(cn))
    println("cn=$cn  live BA=$lB  jl(no-cap)=$jB_off  jl(cap)=$jB_on")
end
