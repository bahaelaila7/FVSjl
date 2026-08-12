using FVSjl, SQLite
db="/workspace/SQLite_FIADB_ENTIRE.db"; var=FVSjl.EasternMontana(); bin="/workspace/.emwork/FVSem_clean"
key(cn)="STDIDENT\n$cn\nNUMCYCLE 6.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
cn="5332701010661"
d=mktempdir();kf=joinpath(d,"k.key");write(kf,key(cn));sp=joinpath(d,"k.sum");cd(d)
try;run(pipeline(`timeout 90 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull));catch;end
println("=== LIVE .sum (year rows) ==="); for l in split(read(sp,String),'\n'); occursin(r"^\s*\d{4}\s",l) && println("  ",strip(l)[1:min(60,end)]); end
d2=mktempdir();kf2=joinpath(d2,"k.key");write(kf2,key(cn))
jt=FVSjl.run_keyfile(kf2;variant=var)
println("=== JL .sum (year rows) ==="); for l in split(jt,'\n'); occursin(r"^\s*\d{4}\s",l) && println("  ",strip(l)[1:min(60,end)]); end
