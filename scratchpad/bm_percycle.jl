# Per-cycle live-vs-jl side-by-side for BM multi-cycle DRIFT stands (#140/#169).
# Reuses the DATABASE+NUMCYCLE keyfile pattern; prints EVERY summary cycle row from
# both live FVSbm_clean and jl so the drift magnitude / first-divergent-cycle / sign
# (jl-HIGH BA => under-thin per #140) is directly readable.
using FVSjl, SQLite
db="/workspace/SQLite_FIADB_ENTIRE.db"
var=FVSjl.BlueMountains(); bin="/workspace/.bmwork/FVSbm_clean"
# NUMCYCLE value COLUMN-ALIGNED (cols 11-20) — free-form "NUMCYCLE 5.0" is mis-parsed by
# FVS's fixed-column reader → only 1 projection cycle. Aligned integer → true 5 projections.
key(cn)="STDIDENT\n$cn\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nNUMCYCLE          5\nECHOSUM\nPROCESS\nSTOP\n"
# grab summary rows: lines beginning with a 4-digit year
rows(txt)=[strip(ln) for ln in split(txt,'\n') if occursin(r"^\s*\d{4}\s",ln)]
function runlive(cn)
    d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,key(cn)); sp=joinpath(d,"k.sum"); cd(d)
    try; run(pipeline(`timeout 120 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull)); catch; end
    isfile(sp) ? read(sp,String) : ""
end
function runjl(cn)
    d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,key(cn)); cd(d)
    try; r=FVSjl.run_keyfile(kf; variant=var); return r isa AbstractString ? r : string(r)
    catch e; println("  jl err: ",e); return ""; end
end
for cn in ARGS
    println("\n########## BM stand $cn ##########")
    lr=rows(runlive(cn)); jr=rows(runjl(cn))
    println("--- cols: YEAR AGE TPA BA SDI CCF TOPHT QMD ... (fixed-width; read YEAR/TPA/BA) ---")
    n=max(length(lr),length(jr))
    for i in 1:n
        println("L$(i>length(lr) ? "  (none)" : ": "*lr[i])")
        println("J$(i>length(jr) ? "  (none)" : ": "*jr[i])")
    end
end
