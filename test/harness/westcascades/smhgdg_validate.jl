# =============================================================================
# smhgdg_validate.jl — WC chunk 6 SMHGDG small-tree model, bit-exact vs FVSwc_g16.
#
# Reference: ref_smhgdg_wct01.txt — per-call {sp, h_in, d_in, cr, ptbal, ptba, si_raw, avht} → {hg5, dg5}
# from the FVSwc_g16 DEBUG REGENT (LSTART calibration) dump on wct01 (`DEBUG␠␠1.␠␠1.` / `REGENT`, cycle-1
# scope; the LSTART SMHGDG calls flush before the cycle-0 volume-DEBUG NATCRS segfault). si_raw is the
# SITEAR fed into wc_smhgdg (for DF/sp16 the live dump prints the POST Curtis→King transform, so the ref
# stores the pre-transform raw; wc_smhgdg re-applies the transform). avht = 0.5·AVH (LSTART ATAVH=0).
#
# Groups exercised: WF(2), ES(10), DF(16) — the wct01 small-tree species. Feeds the live per-call inputs
# into the shipped WC_SMH_* coefficients + wc_smhgdg, matching the chunk-3/4/5 "live-inputs" method.
# =============================================================================
using FVSjl
const M = FVSjl

function run_validation()
ref = joinpath(@__DIR__, "ref_smhgdg_wct01.txt")
worst_hg = 0.0f0; worst_dg = 0.0f0; nfail = 0; n = 0
for line in eachline(ref)
    (isempty(strip(line)) || startswith(strip(line), "#")) && continue
    f = split(strip(line))
    sp = parse(Int, f[1]); h = parse(Float32, f[2]); d = parse(Float32, f[3])
    cr = parse(Float32, f[4]); ptbal = parse(Float32, f[5]); ptba = parse(Float32, f[6])
    si = parse(Float32, f[7]); avht = parse(Float32, f[8])
    hg5_ref = parse(Float32, f[9]); dg5_ref = parse(Float32, f[10])
    hg5, dg5 = M.wc_smhgdg(sp, h, d, cr, ptbal, ptba, si, avht)
    dh = abs(hg5 - hg5_ref); dd = abs(dg5 - dg5_ref)
    worst_hg = max(worst_hg, dh); worst_dg = max(worst_dg, dd)
    n += 1
    if dh > 1f-3 || dd > 1f-3
        nfail += 1
        println("MISMATCH sp=$sp h=$h d=$d: hg5 jl=$hg5 ref=$hg5_ref (Δ$dh); dg5 jl=$dg5 ref=$dg5_ref (Δ$dd)")
    end
end
println("SMHGDG: $n calls, $(n-nfail) bit-exact, worst |ΔHG5|=$worst_hg, worst |ΔDG5|=$worst_dg")
nfail == 0 || error("SMHGDG validation FAILED ($nfail mismatches)")
end
run_validation()
