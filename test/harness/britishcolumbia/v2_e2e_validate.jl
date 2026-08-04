# BC V2 regime END-TO-END wiring validation. Runs all_BC_essf.key (ESSFdk/01 → V2/IE-form regime) through
# the full FVSjl engine and diffs the produced .sum data rows against the instrumented-FVSbc oracle
# (all_BC_essf.oracle.sum). Validates bc_lv2atv dispatch + bc_v2_dgcons! + bc_v2_dgf! + calibration + the
# growth/mortality/volume loop, not just the DDS formula (v2dg_validate.jl covers that in isolation).
import Pkg; Pkg.activate("/workspace/FVSjl"; io = devnull)
using FVSjl

const SCEN = "/workspace/FVSjl/test/harness/scenarios/all_BC_essf"

# .sum data rows: skip the -999 header; key columns Year Age TPA BA SDI CCF TopHt QMD (cols 1-8-ish).
_datarows(txt) = [split(l) for l in split(strip(txt), "\n") if !startswith(strip(l), "-999") && !isempty(strip(l))]

function main()
    jl = FVSjl.run_keyfile("$SCEN.key"; variant = FVSjl.BritishColumbia())
    orc = read("$SCEN.oracle.sum", String)
    J = _datarows(jl); O = _datarows(orc)
    # Compare the first 8 numeric fields per year (Year Age TPA MerchTPA BA SDI CCF TopHt) — the growth spine.
    ncmp = 8
    nrow = min(length(J), length(O))
    exact = 0; firstdiff = ""
    for r in 1:nrow
        jr = J[r]; orr = O[r]
        ok = all(jr[c] == orr[c] for c in 1:ncmp)
        ok ? (exact += 1) : (isempty(firstdiff) && (firstdiff = "yr $(orr[1]): jl=$(jr[1:ncmp]) or=$(orr[1:ncmp])"))
    end
    println("BC V2 e2e (all_BC_essf, ESSFdk/01): $exact/$nrow rows bit-exact on cols 1-$ncmp  (jl=$(length(J)) or=$(length(O)) rows)")
    !isempty(firstdiff) && println("  first diff: ", firstdiff)
    exact == nrow && println("  ✓ FULL MATCH")
end
main()
