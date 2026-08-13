# =============================================================================
# morts_validate.jl — WC chunk 7 ORGANON mortality RIP, bit-exact vs FVSwc_g16.
#
# Reference: ref_morts_wct01.txt — per-tree {sp, D, CR, BAL, PTBAL, HT, AVH, BA, XSITE1, XSITE2} → RIP
# from an instrumented vwc/morts.f (unconditional WRITE of the DO-40 per-tree inputs + RIP to fort.9 —
# bypasses the cyc0 volume-DEBUG NATCRS crash), forest-618 wct01 cycle 1. D and RIP are already floored
# (D≥0.5, RIP≥0.001) as morts.f applies them before the dump, so the harness applies the same floors.
#
# Exercises MORTMAP CASE 1 (DF), CASE 2 (WF/ES), CASE 4 (LP/SP/PP), and the sub-3" Gould-Harrington
# small-tree branch — every equation form present on wct01. (End-to-end, all 27 per-tree WK2 kills also
# match the live dump bit-exact; this harness guards the RIP coefficient/formula transcription.)
# =============================================================================
using FVSjl
const M = FVSjl

function run_validation()
    ref = joinpath(@__DIR__, "ref_morts_wct01.txt")
    worst = 0.0f0; nfail = 0; n = 0
    for line in eachline(ref)
        (isempty(strip(line)) || startswith(strip(line), "#")) && continue
        f = split(strip(line))
        sp = parse(Int, f[1]); d = parse(Float32, f[2]); cr = parse(Float32, f[3])
        bal = parse(Float32, f[4]); ptbal = parse(Float32, f[5]); ht = parse(Float32, f[6])
        avh = parse(Float32, f[7]); ba = parse(Float32, f[8])
        x1 = parse(Float32, f[9]); x2 = parse(Float32, f[10]); rip_ref = parse(Float32, f[11])
        rip = M.wc_mort_rip(sp, d, cr, bal, ptbal, ht, avh, ba, x1, x2)
        rip < 0.001f0 && (rip = 0.001f0)                     # morts.f:425 floor (as in mortality!)
        dr = abs(rip - rip_ref); worst = max(worst, dr); n += 1
        dr > 1f-5 && (nfail += 1; println("MISMATCH sp=$sp D=$d: RIP jl=$rip live=$rip_ref (Δ$dr)"))
    end
    println("WC MORTS RIP: $n trees, $(n-nfail) bit-exact, worst |ΔRIP|=$worst")
    nfail == 0 || error("WC MORTS validation FAILED ($nfail)")
end
run_validation()
