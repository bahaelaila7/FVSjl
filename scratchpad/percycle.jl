# General per-cycle live-vs-jl dumper for any western variant (true 5-projection, col-aligned NUMCYCLE).
# Usage: julia --project=. scratchpad/percycle.jl <VAR> <cn> [cn...]
# Prints per-cycle YEAR/TPA/BA/SDI/CCF/TopHt/QMD for live oracle and jl → classify DG(BA) vs AUTOES(TPA).
using FVSjl
db="/workspace/SQLite_FIADB_ENTIRE.db"
cfg=Dict("BM"=>("/workspace/.bmwork/FVSbm_clean",FVSjl.BlueMountains()),
         "EM"=>("/workspace/.emwork/FVSem_clean",FVSjl.EasternMontana()),
         "IE"=>("/workspace/.iework/FVSie_clean",FVSjl.InlandEmpire()),
         "UT"=>("/workspace/.utwork/FVSut_clean",FVSjl.Utah()),
         "TT"=>("/workspace/.ttwork/FVStt_clean",FVSjl.Teton()),
         "CI"=>("/workspace/.ciwork/FVSci_clean",FVSjl.CentralIdaho()),
         "CR"=>("/workspace/.crwork/FVScr_clean",FVSjl.CentralRockies()))
V=ARGS[1]; bin,var=cfg[V]
key(cn)="STDIDENT\n$cn\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nNUMCYCLE" * lpad("5",12) * "\nECHOSUM\nPROCESS\nSTOP\n"
# fixed-col parse: YEAR(1-4) TPA(9-14) BA(15-18) SDI(19-23) CCF(24-27) TopHt(28-31) QMD(32-36)
function prows(txt)
    r=[]
    for ln in split(txt,'\n'); s=rstrip(ln); length(s)<40 && continue
        y=tryparse(Int,strip(s[1:4])); (y===nothing||y<1000||y>3000) && continue
        push!(r,(y,strip(s[9:14]),strip(s[15:18]),strip(s[19:23]),strip(s[24:27]),strip(s[28:31]),strip(s[32:36])))
    end; r
end
function runlive(cn); d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,key(cn)); sp=joinpath(d,"k.sum"); cd(d)
    try; run(pipeline(`timeout 120 $bin --keywordfile=$kf`;stdout=devnull,stderr=devnull)); catch; end
    isfile(sp) ? read(sp,String) : ""; end
function runjl(cn); d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,key(cn)); cd(d)
    try; r=FVSjl.run_keyfile(kf; variant=var); return r isa AbstractString ? r : string(r); catch e; println("jlerr ",e); return ""; end; end
for cn in ARGS[2:end]
    println("\n######### $V stand $cn #########  YEAR  TPA   BA  SDI  CCF  TopHt QMD")
    L=prows(runlive(cn)); J=prows(runjl(cn)); Jd=Dict(r[1]=>r for r in J)
    for lr in L
        y=lr[1]; jr=get(Jd,y,nothing)
        println("  L $y  ", join(lr[2:end],"  "))
        println("  J $y  ", jr===nothing ? "(none)" : join(jr[2:end],"  "))
    end
end
