# signature.jl — for each stand in a list, run live FVS + FVSjl on the (indexed) sub-DB and report the
# FIRST diverging .sum column + cycle + direction (jl HIGH/LOW vs live). Grouping these signatures reveals
# whether a batch of failures shares a cause. Reuses the FIA_DB sub-DB + live binary.
#
# Usage: FIA_DB=<subdb> julia --project=. test/harness/fia/signature.jl <standlist> <SN|NE|CS|LS>

using FVSjl
const DB = get(ENV, "FIA_DB", "/workspace/SQLite_FIADB_ENTIRE.db")
const BIN = Dict("SN"=>"/workspace/FVSjl/tmp/oracles/FVSsn_new","NE"=>"/workspace/FVSjl/tmp/oracles/FVSne_new","CS"=>"/workspace/FVSjl/tmp/oracles/FVScs_new","LS"=>"/workspace/FVSjl/tmp/oracles/FVSls_new")
const VAR = Dict("SN"=>FVSjl.Southern(),"NE"=>FVSjl.Northeast(),"CS"=>FVSjl.CentralStates(),"LS"=>FVSjl.LakeStates())
const COLS = ["TPA","BA","SDI","CCF","TopHt","QMD"]

keytext(cn) = """
STDIDENT
$cn
DATABASE
DSNin
$DB
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE         5.0
ECHOSUM
PROCESS
STOP
"""

# parse .sum → year => [TPA,BA,SDI,CCF,TopHt,QMD]
# FVS .sum is FIXED-WIDTH (sumout.f FMT 9014, eastern SN/NE/CS/LS); split() misparses on field overflow
# (CCF≥1000 abuts SDI). Parse fixed cols: TPA[9:14] BA[15:18] SDI[19:23] CCF[24:27] TopHt[28:31] QMD[32:36].
function psum(text)
    d = Dict{Int,Vector{Float64}}()
    fcols = ((9,14),(15,18),(19,23),(24,27),(28,31),(32,36))
    for ln in split(text, '\n')
        s = rstrip(ln); length(s) < 36 && continue
        y = tryparse(Int, strip(s[1:4])); (y===nothing || y<1000 || y>3000) && continue
        v = Float64[]; ok = true
        for (a,b) in fcols; x = tryparse(Float64, strip(s[a:b])); x===nothing && (ok=false; break); push!(v, x); end
        ok && (d[y] = v)
    end
    d
end

function sig(cn, var, bin, dir)
    key = joinpath(dir, "s.key"); write(key, keytext(cn))
    for f in ("s.sum","s.out"); fp=joinpath(dir,f); isfile(fp) && rm(fp); end
    try; run(pipeline(`$bin --keywordfile=$key`; stdout=devnull, stderr=devnull)); catch; end
    sp = joinpath(dir,"s.sum"); L = isfile(sp) ? psum(read(sp,String)) : Dict{Int,Vector{Float64}}()
    isempty(L) && return "LIVE_NOSUM"
    J = try psum(FVSjl.run_keyfile(key; variant=var)) catch e; return "JLERR"; end
    isempty(J) && return "JL_NOSUM"
    for y in sort(collect(keys(L)))
        haskey(J, y) || return "yr$y:MISSING"
        for c in 1:6
            if L[y][c] != J[y][c]
                dir_s = J[y][c] > L[y][c] ? "HIGH" : "LOW"
                return "$(COLS[c])@$y:$dir_s(jl$(J[y][c])/lv$(L[y][c]))"
            end
        end
    end
    return "PASS"
end

function main(listfile, v)
    bin = BIN[v]; var = VAR[v]; dir = mktempdir()
    for l in eachline(listfile)
        isempty(strip(l)) && continue
        cn = split(strip(l), '\t')[1]
        println(cn, "  ", sig(cn, var, bin, dir))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: signature.jl <standlist> <variant>")
    main(ARGS[1], ARGS[2])
end
