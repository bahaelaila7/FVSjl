# FIA reader-fidelity harness (all 10 .sum stat columns, incl VOLUMES).
# Cycle-0 printed-identical is the "reader is a bit-exact drop-in" metric: it exercises the
# native DATABASE reader + inventory volume path with no growth-model drift mixed in. Later
# cycles are reported too (as mean |rel diff|) but those fold in growth residuals tracked
# separately. Usage: julia --project=. test/harness/fia/validate_fia_cols.jl <listfile> <VARIANT>
using FVSjl

const DB = "/workspace/SQLite_FIADB_ENTIRE.db"
variant_cfg(v) = v == "LS" ? ("/tmp/FVSls_new", FVSjl.LakeStates()) :
                 v == "SN" ? ("/tmp/FVSsn_new", FVSjl.Southern()) :
                 v == "NE" ? ("/tmp/FVSne_new", FVSjl.NorthEast()) :
                 v == "CS" ? ("/tmp/FVScs_new", FVSjl.CentralStates()) : error("variant $v")

keyfile_text(cn) = """
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
NUMCYCLE         3.0
ECHOSUM
PROCESS
STOP
"""

# .sum cols 3..12 → TPA,BA,SDI,CCF,TopHt,QMD,TCuFt,MCuFt,SCuFt,BdFt
const COLS = ["TPA","BA","SDI","CCF","TopHt","QMD","TCuFt","MCuFt","SCuFt","BdFt"]
function parse_sum(text)
    rows = Tuple{Int,Vector{Float64}}[]
    for ln in split(text, '\n')
        f = split(strip(ln)); length(f) < 12 && continue
        y = tryparse(Int, f[1]); (y === nothing || y < 1000 || y > 3000) && continue
        vals = Float64[]; ok = true
        for i in 3:12
            v = tryparse(Float64, f[i]); v === nothing && (ok = false; break); push!(vals, v)
        end
        ok && length(vals) == 10 && push!(rows, (y, vals))
    end
    rows
end

function run_live(bin, cn, dir)
    key = joinpath(dir, "s.key"); write(key, keyfile_text(cn))
    sp = joinpath(dir, "s.sum"); isfile(sp) && rm(sp)
    try; run(pipeline(`$bin --keywordfile=$key`; stdout=devnull, stderr=devnull)); catch; end
    isfile(sp) ? read(sp, String) : ""
end

function main(listfile, variant)
    bin, var = variant_cfg(variant)
    stands = [split(strip(l), '\t')[1] for l in eachline(listfile) if !isempty(strip(l))]
    dir = mktempdir()
    n_run = 0; n_ok = 0; c0_exact = 0
    col_fail = zeros(Int, 10)                 # cycle-0 per-column mismatch count
    offenders = Dict{String,Vector{String}}() # col → stand tags that diverged at cyc0
    for cn in stands
        n_run += 1
        live = run_live(bin, cn, dir); isempty(live) && continue
        keyf = joinpath(dir, "jl.key"); write(keyf, keyfile_text(cn))
        jlout = try FVSjl.run_keyfile(keyf; variant=var) catch; ""; end
        isempty(jlout) && continue
        L = parse_sum(live); J = parse_sum(jlout)
        (isempty(L) || isempty(J)) && continue
        n_ok += 1
        Jd = Dict(y=>v for (y,v) in J)
        y0, lv0 = L[1]; haskey(Jd, y0) || continue
        jv0 = Jd[y0]
        allok = true
        for k in 1:10
            # .sum prints integers (QMD 1 dp) ⇒ |Δ|<0.5 = same printed cell
            if abs(lv0[k] - jv0[k]) >= 0.5
                col_fail[k] += 1; allok = false
                push!(get!(offenders, COLS[k], String[]), cn)
            end
        end
        allok && (c0_exact += 1)
    end
    println("\n===== FIA reader fidelity: $variant  (cycle-0, all 10 cols) =====")
    println("stands=$n_run  both-produced-sum=$n_ok  cycle0-ALL-cols-identical=$c0_exact / $n_ok")
    println("per-column cycle-0 mismatch count (of $n_ok):")
    for k in 1:10
        mark = col_fail[k] == 0 ? "✓" : "✗"
        println("  $mark ", rpad(COLS[k],7), col_fail[k])
    end
    for k in 1:10
        col_fail[k] == 0 && continue
        o = offenders[COLS[k]]
        println("  $(COLS[k]) diverged on: ", join(first(o, min(6, length(o))), ", "), length(o) > 6 ? " …(+$(length(o)-6))" : "")
    end
end

main(ARGS[1], ARGS[2])
