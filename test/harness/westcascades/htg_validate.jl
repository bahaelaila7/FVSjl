# WC large-tree HTG cyc0 validation vs the live FVSwc_clean oracle.
#
# Reference: ref_htg_wct01.txt — per-tree HTGF/FINDAG/HTCALC DEBUG dump from FVSwc_clean on
# wct01 (stand 248112, forest JFOR=6, 10-yr cycle). Feeding the live per-tree inputs
# {SINDX, D, H, ICR, AVH, BA, DG} into the SHIPPED WC height functions isolates the westside
# FINDAG→HTCALC site-curve + HGMDCR/HGMDRH modifier math from the (not-yet-ported) site-index
# and stand-stats chunks (SINDX/AVH are stand/species constants read from the dump, like the
# DGF harness fed CONSPP). MEASURED: HTCON≡0 (no height calibration), XHT≡1, SCALE=1.0 (FINT/YR
# = 10/10). All 24 trees with a printed final HTG are DEFAULT-branch (groups WF,ES,LP,SP,PP,DF);
# the 3 H≥HTMAX trees (4,5,10) print no final HTG and are excluded (source-faithful, unvalidated).
#
# RESULT (2026-08-13): see the printed worst |Δ|.
using FVSjl, Printf
const M = FVSjl

function main()
    ref = joinpath(@__DIR__, "ref_htg_wct01.txt")
    rows = [split(strip(l)) for l in eachline(ref) if !isempty(strip(l)) && !startswith(strip(l), "#")]
    n = 0; worst = 0f0; wl = ""; agefail = 0
    for r in rows
        ispc = parse(Int, r[2]); sindx = parse(Float32, r[3]); d = parse(Float32, r[4])
        h = parse(Float32, r[5]); icr = parse(Float32, r[6]); avh = parse(Float32, r[7])
        ba = parse(Float32, r[8]); dg = parse(Float32, r[9])
        sitage_l = parse(Float32, r[10]); live = parse(Float32, r[12])
        # D2 = D + DG/BARK; bark not loaded in this formula harness. HTMAX2 = HDRAT1·D2+HDRAT2
        # only binds when H+HTG>HTMAX2 (never, for these below-HTMAX trees). d2≈d+dg (bark→1) is a
        # conservative under-estimate of the true (larger) D2 ⇒ if it doesn't bind here it never does.
        d2 = d + dg
        sitage, sitht, agmax, htmax, htmax2 = M.wc_findag(ispc, d, d2, h, sindx)
        sitage != sitage_l && (agefail += 1)
        htg = M.wc_htg_default(ispc, sindx, d, h, icr, avh, ba, dg, d2, sitage, sitht, agmax, htmax2)
        # scale=1, exp(HTCON)=1, MISHGF=1, no SIZCAP breach (tall-tree cap far above these heights)
        e = abs(htg - live); n += 1
        if e > worst; worst = e; wl = @sprintf("i=%s ispc=%d jl=%.6f live=%.6f (Δage %g→%g)",
                                               r[1], ispc, htg, live, sitage, sitage_l); end
    end
    @printf("WC HTG DEFAULT-branch: n=%d  SITAGE mismatches=%d  worst|Δ|=%.6f  %s\n", n, agefail, worst, wl)
    @assert agefail == 0 "WC FINDAG produced a different SITAGE than live"
    @assert worst <= 1f-3 "WC HTG diverged beyond print precision"
    println("PASS — WC large-tree HTG bit-exact vs live (within DEBUG print precision).")
end

main()
