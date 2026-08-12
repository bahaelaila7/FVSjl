using FVSjl
db="/workspace/SQLite_FIADB_ENTIRE.db"; var=FVSjl.EasternMontana(); bin="/workspace/.emwork/FVSem_clean"
cn="11864108010690"
function key(noae)
    ae = noae ? "NOAUTOES\n" : ""
    "STDIDENT\n$cn\n$ae" * "DATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nNUMCYCLE" * lpad("5",12) * "\nECHOSUM\nPROCESS\nSTOP\n"
end
tpa(txt)=[(strip(ln)[1:4], strip(split(strip(ln))[3])) for ln in split(txt,'\n') if occursin(r"^\s*\d{4}\s",ln)]
function runlive(noae); d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,key(noae)); sp=joinpath(d,"k.sum"); cd(d)
    try; run(pipeline(`timeout 120 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull)); catch; end
    isfile(sp) ? read(sp,String) : ""; end
function runjl(noae); d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,key(noae)); cd(d)
    try; r=FVSjl.run_keyfile(kf; variant=var); return r isa AbstractString ? r : string(r); catch e; return "ERR $e"; end; end
println("TPA per cycle (year:tpa):")
println("live default : ", tpa(runlive(false)))
println("live NOAUTOES: ", tpa(runlive(true)))
println("jl   default : ", tpa(runjl(false)))
println("jl   NOAUTOES: ", tpa(runjl(true)))
