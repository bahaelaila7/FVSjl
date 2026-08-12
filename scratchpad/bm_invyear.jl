using FVSjl
db="/workspace/SQLite_FIADB_ENTIRE.db"; var=FVSjl.BlueMountains(); bin="/workspace/.bmwork/FVSbm_clean"
# INVYEAR after STDIDENT, NUMCYCLE after DATABASE END
key(cn)="STDIDENT\n$cn\nINVYEAR 2018\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nNUMCYCLE 5.0\nECHOSUM\nPROCESS\nSTOP\n"
rows(txt)=[strip(ln) for ln in split(txt,'\n') if occursin(r"^\s*\d{4}\s",ln)]
function runlive(cn); d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,key(cn)); sp=joinpath(d,"k.sum"); cd(d)
    try; run(pipeline(`timeout 120 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull)); catch; end
    isfile(sp) ? read(sp,String) : ""; end
function runjl(cn); d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,key(cn)); cd(d)
    try; r=FVSjl.run_keyfile(kf; variant=var); return r isa AbstractString ? r : string(r); catch e; println("jlerr ",e); return ""; end; end
for cn in ARGS
    lr=rows(runlive(cn)); jr=rows(runjl(cn))
    println("cn=$cn  LIVE cycles=$(length(lr))  JL cycles=$(length(jr))")
    println("  live years: ", join([split(r)[1] for r in lr]," "))
    println("  jl   years: ", join([split(r)[1] for r in jr]," "))
end
