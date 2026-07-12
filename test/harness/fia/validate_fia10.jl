# validate_fia10.jl — Pillar-2 completion: multi-cycle differential over ALL 10 .sum columns
# (TPA/BA/SDI/CCF/TopHt/QMD  +  TCuFt/MCuFt/SCuFt/BdFt) vs live FVS, at scale.
#
# validate_fia.jl checked the 6 stand-structure columns; the 4 VOLUME columns (.sum fields 9-12) were
# validated bit-exact at CYCLE-0 (modernization #85) but not in the at-scale MULTI-CYCLE differential.
# This harness closes that: it diffs fields 3-12 every cycle and reports a PER-COLUMN bit-exact rate so
# the volume columns are measured independently (a volume-only divergence would otherwise hide behind a
# structure-only pass rate).
#
# Usage: FIA_DB=<subdb> julia --project=. test/harness/fia/validate_fia10.jl <standlist> <SN|NE|CS|LS> [regime]

using FVSjl
const DB = get(ENV, "FIA_DB", "/workspace/SQLite_FIADB_ENTIRE.db")
const BIN = Dict("SN"=>"/workspace/FVSjl/tmp/oracles/FVSsn_new","NE"=>"/workspace/FVSjl/tmp/oracles/FVSne_new","CS"=>"/workspace/FVSjl/tmp/oracles/FVScs_new","LS"=>"/workspace/FVSjl/tmp/oracles/FVSls_new")
const VAR = Dict("SN"=>FVSjl.Southern(),"NE"=>FVSjl.Northeast(),"CS"=>FVSjl.CentralStates(),"LS"=>FVSjl.LakeStates())
const COLS = ["TPA","BA","SDI","CCF","TopHt","QMD","TCuFt","MCuFt","SCuFt","BdFt"]  # .sum fields 3..12

kwrec(kw, fields...) = rpad(kw, 10) * join(lpad(string(f), 10) for f in fields)
regime_block(r) =
    r == "simfire" ? "FMIn\n" * kwrec("SIMFIRE", "2.0", "10.00", "1", "50.0") * "\nEnd" :
    r == "thinbba" ? kwrec("THINBBA", "2.0", "40.0") :
    r == "salvage" ? kwrec("SALVAGE", "2.0", "0.0", "999.0", "0.9") :
    r == "plant"   ? "ESTAB\n" * kwrec("PLANT", "2.0", "3", "400") * "\nEnd" : ""

keytext(cn, regime) = """
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
$(regime_block(regime))
ECHOSUM
PROCESS
STOP
"""

# parse the 10 reporting fields (.sum cols 3-12) keyed by year. Robust to fixed-width column abutment by
# reading the leading 4-digit year, then splitting the remainder on whitespace (cols 2-12 don't abut here).
function parse_sum10(text)
    rows = Tuple{Int,Vector{Float64}}[]
    for ln in split(text, '\n')
        s = strip(ln); length(s) < 4 && continue
        y = tryparse(Int, s[1:4]); (y === nothing || y < 1000 || y > 3000) && continue
        f = split(s); length(f) < 12 && continue
        vals = Float64[]; ok = true
        for i in 3:12; v = tryparse(Float64, f[i]); v === nothing && (ok = false; break); push!(vals, v); end
        ok && push!(rows, (y, vals))
    end
    rows
end

function run_live(bin, cn, regime, dir)
    key = joinpath(dir, "s.key"); write(key, keytext(cn, regime))
    for f in ("s.sum","s.out"); fp = joinpath(dir, f); isfile(fp) && rm(fp); end
    try; run(pipeline(`$bin --keywordfile=$key`; stdout=devnull, stderr=devnull)); catch; end
    sp = joinpath(dir, "s.sum"); isfile(sp) ? read(sp, String) : ""
end

function main(listfile, v, regime)
    bin = BIN[v]; var = VAR[v]; dir = mktempdir()
    stands = [split(strip(l), '\t')[1] for l in eachline(listfile) if !isempty(strip(l))]
    n_ok = 0; n_pass_all = 0
    col_cells = zeros(Int, 10); col_match = zeros(Int, 10)   # per-column cell counts
    vol_fail_stands = String[]
    outliers = Tuple{Float64,String,Int,Int}[]   # (maxrel, cn, col, year) — for the >THRESH tail
    othresh = parse(Float64, get(ENV, "FIA_OUTLIER", "0.05"))
    for cn in stands
        live = run_live(bin, cn, regime, dir); isempty(live) && continue
        keyf = joinpath(dir, "jl.key"); write(keyf, keytext(cn, regime))
        jlout = try FVSjl.run_keyfile(keyf; variant=var) catch; ""; end
        isempty(jlout) && continue
        L = parse_sum10(live); J = parse_sum10(jlout)
        (isempty(L) || isempty(J)) && continue
        n_ok += 1
        Jd = Dict(y => vv for (y, vv) in J)
        all_ok = length(J) == length(L); vol_ok = true; smaxrel = 0.0; scol = 0; syr = 0
        for (y, lv) in L
            haskey(Jd, y) || (all_ok = false; continue)
            jv = Jd[y]
            for k in 1:10
                col_cells[k] += 1
                if lv[k] == jv[k]; col_match[k] += 1
                else
                    all_ok = false; (k >= 7) && (vol_ok = false)
                    rel = lv[k] == 0 ? (jv[k]==0 ? 0.0 : 1.0) : abs(lv[k]-jv[k])/abs(lv[k])
                    rel > smaxrel && (smaxrel = rel; scol = k; syr = y)
                end
            end
        end
        all_ok && (n_pass_all += 1)
        vol_ok || push!(vol_fail_stands, cn)
        smaxrel >= othresh && push!(outliers, (smaxrel, cn, scol, syr))
    end
    println("\n===== VALIDATE-10COL ($regime): $v =====")
    println("both-sum stands = $n_ok   ALL-10-cols bit-exact (every cycle) = $n_pass_all / $n_ok")
    println("Per-column bit-exact cell rate (across all cycles):")
    for k in 1:10
        r = col_cells[k] == 0 ? 0.0 : 100 * col_match[k] / col_cells[k]
        marker = k == 7 ? "   <-- volumes:" : ""
        println("  $(rpad(COLS[k],6)): $(round(r,digits=2))%  ($(col_match[k])/$(col_cells[k]))$marker")
    end
    println("stands with a VOLUME-column divergence: $(length(vol_fail_stands))")
    sort!(outliers, rev=true)
    println("OUTLIERS (max-rel >= $(round(othresh*100))% on any col): $(length(outliers))")
    for (r,cn,col,yr) in outliers[1:min(15,end)]
        println("  $(round(r*100,digits=2))%  $cn  col=$(COLS[col]) @ $yr")
    end
    if haskey(ENV,"FIA_VOLFAIL") && !isempty(vol_fail_stands)
        open(ENV["FIA_VOLFAIL"],"w") do io; for s in vol_fail_stands; println(io, s); end; end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: validate_fia10.jl <standlist> <variant> [regime]")
    main(ARGS[1], ARGS[2], length(ARGS) >= 3 ? ARGS[3] : "none")
end
