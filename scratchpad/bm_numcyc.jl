using FVSjl
db="/workspace/SQLite_FIADB_ENTIRE.db"; var=FVSjl.BlueMountains(); bin="/workspace/.bmwork/FVSbm_clean"
# test 3 NUMCYCLE formats for column-sensitivity
fmts=Dict("freeform_1space"=>"NUMCYCLE 5.0","aligned_col11"=>"NUMCYCLE          5","int_2space"=>"NUMCYCLE  5")
key(cn,nc)="STDIDENT\n$cn\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\n$nc\nECHOSUM\nPROCESS\nSTOP\n"
rows(txt)=[strip(ln) for ln in split(txt,'\n') if occursin(r"^\s*\d{4}\s",ln)]
function runlive(cn,nc); d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,key(cn,nc)); sp=joinpath(d,"k.sum"); cd(d)
    try; run(pipeline(`timeout 120 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull)); catch; end
    isfile(sp) ? read(sp,String) : ""; end
function runjl(cn,nc); d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,key(cn,nc)); cd(d)
    try; r=FVSjl.run_keyfile(kf; variant=var); return r isa AbstractString ? r : string(r); catch; return ""; end; end
cn="645068643126144"
for (name,nc) in fmts
    println("--- $name : '$nc' ---")
    println("  LIVE cycles=$(length(rows(runlive(cn,nc))))   JL cycles=$(length(rows(runjl(cn,nc))))")
end
