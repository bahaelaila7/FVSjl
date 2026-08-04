# BC V2 (LV2ATV) regime END-TO-END regression guard. Runs all_BC_essf.key (ESSFdk/01 → V2/IE-form regime)
# through the full FVSjl engine and diffs the produced .sum against the instrumented-FVSbc oracle
# (all_BC_essf.oracle.sum). Validates the whole V2 growth+mortality+density spine end-to-end:
# bc_lv2atv dispatch + bc_v2_dgf! + HCOR calib + Hamilton RIP mortality + the per-TRIPLE ZZRAN small-tree
# draw (commit f030564) + tripling + summary.
#
# BIT-EXACT-OR-CORNERED: the all_BC_essf fixture carries GARBAGE heights (uncorrelated with DBH), so the
# AVHT40 top-height + the ZZRAN multiplicative small-tree error leave a small, stable residual that is the
# accepted cornered class (plus the oracle .sum's plot-area TPA rounding, 2089.12→2087, task #129). We assert
# each column tracks within a cornered tolerance chosen to PASS the converged state yet FAIL a real regression:
#   - reverting the per-triple ZZRAN fix pushed TopHt to +12% / SDI to +9.5% (both blow these tolerances);
#   - the earlier WK1 over-kill pushed TPA to +57%.
# So this guard catches those without demanding bit-exactness the garbage fixture can't deliver.
import Pkg; Pkg.activate("/workspace/FVSjl"; io = devnull)
using FVSjl

const SCEN = "/workspace/FVSjl/test/harness/scenarios/all_BC_essf"

_datarows(txt) = [split(l) for l in split(strip(txt), "\n") if !startswith(strip(l), "-999") && !isempty(strip(l))]

# .sum data cols: 1 Year, 2 Age, 3 TPA, 4 BA, 5 SDI, 6 CCF, 7 TopHt, 8 QMD (metric/ha for BC).
# (col, name, kind, tol) — kind :rel (relative) or :abs (absolute, for coarse small-integer metric cols).
# Tolerances: ~1.5-2× the converged max Δ (giving margin against env/rounding flake) yet well below the
# regression signatures — per-triple-ZZRAN revert = TopHt +12% / SDI +9.5%; WK1 over-kill = TPA +57%.
const COLS = [(3, "TPA", :rel, 0.02), (4, "BA", :abs, 2.0), (5, "SDI", :rel, 0.07),
              (6, "CCF", :abs, 3.0), (7, "TopHt", :rel, 0.08), (8, "QMD", :rel, 0.05)]

function main()
    jl = FVSjl.run_keyfile("$SCEN.key"; variant = FVSjl.BritishColumbia())
    orc = read("$SCEN.oracle.sum", String)
    J = _datarows(jl); O = _datarows(orc)
    nrow = min(length(J), length(O))
    exact = 0; worst = Dict(c[2] => (0.0, "") for c in COLS); fails = String[]
    for r in 1:nrow
        jr = J[r]; orr = O[r]
        all(jr[c] == orr[c] for c in 3:8) && (exact += 1)
        for (c, name, kind, tol) in COLS
            jv = parse(Float64, jr[c]); ov = parse(Float64, orr[c])
            d = kind === :rel ? (ov == 0 ? abs(jv) : abs(jv - ov) / abs(ov)) : abs(jv - ov)
            d > worst[name][1] && (worst[name] = (d, "yr $(orr[1]): jl=$jv or=$ov"))
            d > tol && push!(fails, "$name yr $(orr[1]): jl=$jv or=$ov (Δ$(kind)=$(round(d, digits=4)) > $tol)")
        end
    end
    println("BC V2 e2e regression guard (all_BC_essf, ESSFdk/01): $nrow cycles, $exact/$nrow bit-exact on cols 3-8")
    for (_, name, kind, tol) in COLS
        d, where = worst[name]
        println("  $name  max Δ$(kind)=$(round(d, digits=4)) (tol $tol)  @ $where")
    end
    if isempty(fails)
        println("  ✓ BIT-EXACT-OR-CORNERED — all columns within cornered tolerance")
    else
        println("  ✗ REGRESSION — ", length(fails), " column-cycle(s) exceed cornered tolerance:")
        for f in fails; println("      ", f); end
    end
    return isempty(fails)
end

main() || exit(1)
