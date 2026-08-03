using FVSjl
const M = FVSjl
const DB = "/workspace/SQLite_FIADB_ENTIRE.db"
keytext(cn) = "STDIDENT\n$cn\nDATABASE\nDSNin\n$DB\nStandSQL\nSELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nTreeSQL\nSELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\nNUMCYCLE         3.0\nECHOSUM\nPROCESS\nSTOP\n"
function parse_sum(text)
    rows = Tuple{Int,Vector{Float64}}[]
    cols = ((9,14),(15,18),(19,23),(24,27),(28,31),(32,36),(37,42),(43,48),(49,54),(55,60))
    for ln in split(text,'\n')
        s = rstrip(ln); length(s) < 60 && continue
        y = tryparse(Int, strip(s[1:4])); (y===nothing || y<1000 || y>3000) && continue
        vals = Float64[]; ok=true
        for (a,b) in cols; v=tryparse(Float64,strip(s[a:b])); v===nothing && (ok=false;break); push!(vals,v); end
        ok && length(vals)==10 && push!(rows,(y,vals))
    end
    rows
end
const COLS = ["TPA","BA","SDI","CCF","TopHt","QMD","TCuFt","MCuFt","SCuFt","BdFt"]
cases = [("CI","/workspace/.ciwork/FVSci_clean",M.CentralIdaho(),["51052888020004","3200106010690"]),
         ("UT","/workspace/.utwork/FVSut_clean",M.Utah(),["198895166020004"]),
         ("EM","/workspace/.emwork/FVSem_clean",M.EasternMontana(),["2342318010690"])]
dir = mktempdir()
for (v,bin,var,cns) in cases, cn in cns
    key = joinpath(dir,"k.key"); write(key, keytext(cn))
    sp = joinpath(dir,"k.sum"); isfile(sp) && rm(sp)
    try; run(pipeline(`$bin --keywordfile=$key`; stdout=devnull, stderr=devnull)); catch; end
    L = parse_sum(read(sp,String))
    jl = M.run_keyfile(key; variant=var)
    J = parse_sum(jl)
    println("\n### $v $cn")
    Jd = Dict(y=>vv for (y,vv) in J)
    for (y,lv) in L
        haskey(Jd,y) || continue
        jv = Jd[y]
        diffs = [ "$(COLS[k]) L=$(lv[k]) J=$(jv[k])" for k in 7:10 if abs(lv[k]-jv[k])>=0.5 ]
        println("  yr $y  TCuFt L=$(lv[7]) J=$(jv[7]) | MCuFt L=$(lv[8]) J=$(jv[8]) | SCuFt L=$(lv[9]) J=$(jv[9]) | BdFt L=$(lv[10]) J=$(jv[10])",
                isempty(diffs) ? "  ✓VOL-EXACT" : "  ✗ "*join(diffs,"; "))
    end
end
