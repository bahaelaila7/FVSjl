# WC large-tree DGF cyc0 validation vs the live FVSwc_clean oracle.
#
# Reference data: ref_dds_wct01_block3.txt — the growth-pass (block 3) per-tree DEBUG DGF dump
# from FVSwc_clean on wct01 (stand 248112, forest JFOR=6), extracted from fort.16:
#   cols: I ISPC CONSPP D BA CR BAL PCCF RELDEN HT AVH  LN(DDS_live)
# CONSPP (=DGCON+COR) is the live per-species constant; feeding it into the SHIPPED WC DEFAULT
# DDS formula + coefficients isolates the westside Wykoff DDS math (calibration COR is baked in).
# Oracle build: /workspace/.wcwork/relink_wc.sh clean  (bin/FVSwc_buildDir/*.o + isoc23 shim).
# Debug keyfile: NUMCYCLE 1 + `DEBUG␠␠␠␠1.␠␠␠␠1.` / `DGF` (cycle-1 scope: dumps DGF before the
# volume-DEBUG segfault). All 27 block-3 trees are DEFAULT-branch (groups 2,5,6,7,11,16).
#
# RESULT (2026-08-13): 27/27 bit-exact vs live LN(DDS); worst |Δ| = 0.00007 = CONSPP F11.4
# print rounding. Proves WC westside-Wykoff DGF reuse.
using FVSjl, Printf
const M = FVSjl

function dds_default(jspc, CONSPP, D, BA, CR, BAL, PCCF, RELHT; jfor::Int = 6)
    dgdsq = M.WC_DGDS[jspc, M.WC_MAPDSQ[jspc, jfor]]
    v = CONSPP + M.WC_DGLD[jspc]*log(D) + CR*(M.WC_DGCR[jspc] + CR*M.WC_DGCRSQ[jspc]) +
        dgdsq*D*D + M.WC_DGDBAL[jspc]*BAL/log(D+1f0) + M.WC_DGPCCF[jspc]*PCCF +
        M.WC_DGHAH[jspc]*RELHT + M.WC_DGLBA[jspc]*log(BA) + M.WC_DGBAL[jspc]*BAL + M.WC_DGBA[jspc]*BA
    v < -9.21f0 && (v = -9.21f0)
    return v
end

function main()
    ref = joinpath(@__DIR__, "ref_dds_wct01_block3.txt")
    rows = [split(strip(l)) for l in eachline(ref) if !isempty(strip(l))]
    n = 0; worst = 0f0; wl = ""
    for r in rows
        ISPC=parse(Int,r[2]); CONSPP=parse(Float32,r[3]); D=parse(Float32,r[4]); BA=parse(Float32,r[5])
        CR=parse(Float32,r[6]); BAL=parse(Float32,r[7]); PCCF=parse(Float32,r[8]); HT=parse(Float32,r[10])
        AVH=parse(Float32,r[11]); live=parse(Float32,r[12])
        RELHT = AVH > 0f0 ? min(HT/AVH, 1.5f0) : 0f0
        my = dds_default(M.WC_MAPSPC[ISPC], CONSPP, D, BA, CR, BAL, PCCF, RELHT)
        e = abs(my - live); n += 1
        e > worst && (worst = e; wl = @sprintf("ISPC=%d jl=%.4f live=%.4f", ISPC, my, live))
    end
    @printf("WC DGF DEFAULT-branch: n=%d  worst|Δ|=%.5f (%s)\n", n, worst, wl)
    @assert worst <= 0.0006f0 "WC DGF diverged beyond CONSPP print precision"
    println("PASS — WC large-tree DDS bit-exact vs live (within F11.4 CONSPP rounding).")
end

main()
