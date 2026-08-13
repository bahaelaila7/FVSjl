# WC crown-ratio cyc0 validation vs the live FVSwc_clean oracle (chunk 5).
#
# Reference: ref_crown_wct01.txt — the cycling-CROWN DEBUG dump from FVSwc_clean on wct01
# (stand 248112, forest JFOR=6). Each row = one tree's per-species 9001 header (SDIAC, SDIDEF,
# live A/B/C) + its 9002 (X, CRNEW). Feeding the live X into the SHIPPED WC crown coefficients +
# Weibull inverse-CDF isolates the WC-specific CRCONS table + the crown math from the shared
# rank/SCALE engine (already validated in TT/CI/NC). Two checks per tree:
#   (a) CONSTRUCTION: reconstruct B,C from WC_CRC0/1 + WEIB* using RELSDI=SDIAC/SDIDEF (printed
#       SDIAC F8.2 ⇒ ~F4 precision) and compare to the live A/B/C.
#   (b) FORMULA: CRNEW = A + B·(−ln(1−X))^(1/C) with the LIVE A/B/C and live X, vs live CRNEW.
# All 81 tripled trees span groups WF(2), ES(11), LP(16), SP(5), PP(6), DF(7) — every wct01 group.
using FVSjl, Printf
const M = FVSjl

function main()
    ref = joinpath(@__DIR__, "ref_crown_wct01.txt")
    rows = [split(strip(l)) for l in eachline(ref) if !isempty(strip(l)) && !startswith(strip(l), "#")]
    cimap = M.coefficients(M.WestCascades()).species[:crown_imap]
    n = 0
    worst_form = 0f0; wlf = ""
    worst_bc = 0f0; wlbc = ""
    groups = Set{Int}()
    for r in rows
        ispc = parse(Int, r[1]); x = parse(Float32, r[2]); crnew_live = parse(Float32, r[3])
        sdiac = parse(Float32, r[4]); sdidef = parse(Float32, r[5])
        A_live = parse(Float32, r[6]); B_live = parse(Float32, r[7]); C_live = parse(Float32, r[8])
        grp = Int(cimap[ispc]); push!(groups, grp)
        # (a) construction from shipped coefficients
        relsdi = sdiac / sdidef; relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = M.WC_CRC0[grp] + M.WC_CRC1[grp] * relsdi * 100f0
        A = M.WC_WEIBA[grp]
        B = M.WC_WEIBB0[grp] + M.WC_WEIBB1[grp] * acrnew; B < 3f0 && (B = 3f0)
        C = M.WC_WEIBC0[grp] + M.WC_WEIBC1[grp] * acrnew; C < 2f0 && (C = 2f0)
        ebc = max(abs(B - B_live), abs(C - C_live), abs(A - A_live))
        if ebc > worst_bc; worst_bc = ebc
            wlbc = @sprintf("ispc=%d grp=%d B %.4f/%.4f C %.4f/%.4f", ispc, grp, B, B_live, C, C_live); end
        # (b) formula with the LIVE A/B/C
        crnew = A_live + B_live * (-log(1f0 - x))^(1f0 / C_live)
        ef = abs(crnew - crnew_live); n += 1
        if ef > worst_form; worst_form = ef
            wlf = @sprintf("ispc=%d x=%.5f jl=%.5f live=%.5f", ispc, x, crnew, crnew_live); end
    end
    @printf("WC CROWN: n=%d groups=%s\n", n, sort(collect(groups)))
    @printf("  (a) B/C construction worst|Δ| = %.5f  (%s)\n", worst_bc, wlbc)
    @printf("  (b) CRNEW formula   worst|Δ| = %.5f  (%s)\n", worst_form, wlf)
    @assert worst_form <= 1f-3 "WC CRNEW formula diverged beyond print precision"
    @assert worst_bc  <= 5f-3 "WC B/C construction diverged beyond SDIAC print precision"
    println("PASS — WC crown ratio bit-exact vs live (within DEBUG print precision).")
end

main()
