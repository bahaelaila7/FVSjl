using FVSjl, SQLite
db="/workspace/SQLite_FIADB_ENTIRE.db"
# (CN, live_BA, live_TPA) from the sweep
divs=[("850447806290487",227.0,1526),("850447807290487",59.0,2365),("248613816489998",317.0,819),
      ("248615564489998",252.0,3269),("1288130126290487",167.0,1474),("1123874220290487",468.0,625),
      ("1123874415290487",251.0,592)]
function finalba(txt)
    last=""; for ln in split(txt,'\n'); occursin(r"^\s*\d{4}",ln) && (last=strip(ln)); end
    isempty(last) && return (0,0.0); f=split(last); (parse(Int,f[3]), parse(Float64,f[4]))
end
for (cn,lB,lT) in divs
    kt="STDIDENT\n$cn\nNUMCYCLE 5.0\nDATABASE\nDSNin\n$db\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nECHOSUM\nPROCESS\nSTOP\n"
    d=mktempdir(); kf=joinpath(d,"k.key"); write(kf,kt)
    jt=""; try; jt=FVSjl.run_keyfile(kf; variant=FVSjl.Klamath()); catch e; println("$cn CRASH: ",sprint(showerror,e)[1:60]); continue; end
    jT,jB=finalba(jt)
    d2 = lB>0 ? abs(jB-lB)/lB*100 : 0
    println("$cn: jl BA=$jB TPA=$jT | live BA=$lB TPA=$lT (Δ$(round(d2,digits=1))%)")
end
