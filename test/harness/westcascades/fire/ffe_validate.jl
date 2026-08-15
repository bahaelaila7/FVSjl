# WC FFE cyc0 (1990) fuel-loading validation vs the instrumented FVSwc_clean oracle (chunks F2/F4).
#
# Reference (ref_ffe_wct01.txt): the fmcba.f WRITE dumps from FVSwc_clean on wct01 (stand 248112,
# forest 618, habitat 52) — COVTYP, PERCOV, FLIVE, STFUEL, per-tree crown width, and the FMCFMD
# fuel-model selection. The Fortran was instrumented (unconditional WRITE to fort.9), measured, then
# restored to pristine (grep-verified 0 markers, relinked clean).
#
# This feeds the 1990 inventory stand into wc_cwcalc + wc_live/dead_fuel_loading and checks:
#   (a) per-tree crown width (WC R6 Crookston CAMAP + forest-618 BF),
#   (b) PERCOV = 100·(1−exp(−ΣcrownArea/43560)),
#   (c) FLIVE (herb/shrub) and STFUEL (11 dead size classes) at the oracle PERCOV.
using FVSjl, Printf
const M = FVSjl

function main()
    key = "/tmp/wct01_ffe.key"
    s = nothing
    for st in M.each_stand(key; variant = M.WestCascades(), faithful = true)
        M.notre!(st); M.setup_growth!(st); s = st; break
    end
    ba = s.plot.basal_area; el = s.plot.elevation
    @printf("stand BA(FFE)=%.5f  EL=%.2f  (oracle BAREA=85.13199 EL=35)\n", ba, el)

    # (a)+(b): per-tree crown width + PERCOV
    totcra = 0f0; worst_cw = 0f0
    for i in 1:s.trees.n
        s.trees.tpa[i] > 0f0 || continue
        sp = Int(s.trees.species[i])
        cw = M.wc_cwcalc(sp, s.trees.dbh[i], s.trees.height[i], Float32(s.trees.crown_pct[i]), ba, el, 0f0)
        totcra += 3.1415927f0 * cw * cw / 4f0 * s.trees.tpa[i]
    end
    percov = (1f0 - exp(-totcra / 43560f0)) * 100f0
    @printf("(b) TOTCRA jl=%.4f oracle=29106.4473 | PERCOV jl=%.5f oracle=48.73655  Δ=%.5f\n",
            totcra, percov, percov - 48.7365494f0)

    # (c) live + dead fuel loading at the ORACLE percov (isolate the interpolation from the tiny PERCOV Δ)
    op = 48.7365494f0
    herb, shrub = M.wc_live_fuel_loading(2, op)
    @printf("(c) FLIVE  herb jl=%.6f oracle=0.183790  shrub jl=%.6f oracle=0.528011\n", herb, shrub)
    stf = M.wc_dead_fuel_loading(2, op)
    ref_stf = Float32[0.654946208, 0.654946208, 2.77473092, 6.05387020, 6.05387020, 0, 0, 0, 0, 0.532419324, 22.0715027]
    wmax = maximum(abs.(stf .- ref_stf))
    @printf("    STFUEL worst|Δ| = %.6f  (jl[1,3,4,10,11]=%.4f %.4f %.4f %.4f %.4f)\n",
            wmax, stf[1], stf[3], stf[4], stf[10], stf[11])
    println(abs(percov - 48.7365494f0) < 0.05f0 && wmax < 1f-4 ? "PASS (fuel loading bit-exact-or-cornered)" : "CHECK")
end

main()
