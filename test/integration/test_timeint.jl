# test_timeint.jl — TIMEINT cycle-length control (the cycle calendar) vs live Fortran.
#
# By default a cycle is 5 years (sn/grinit.f:149, IY=5). TIMEINT sets the cycle length (IY);
# `TIMEINT 0 10` makes every cycle 10 years. The growth models scale to the cycle length:
#   * diameter — DDS · (FINT/YR), FINT = cycle length, YR = 5 (SN measurement base); a 10-yr
#     cycle ⇒ ×2 the squared-diameter change (dgdriv.f:325/715, verified vs a Fortran DEBUG dump
#     where DDS is FINT-independent but the realized DG scales 2×).
#   * height — HTG · (FINT/5);  mortality — rate^FINT;  year/age step by FINT.
#
# STRUCTURAL correctness (the .sum year/age stepping and the ×2 growth scaling) is exact and is
# the asserted contract here. A small residual remains on the volume/BA columns of the 10-yr run
# (≈2%): the calibrated-species (WK3) COR evolution under a non-5 period — the same class as the
# DGSCOR cubic-volume tail (DIVERGENCES.md §1) — pending the full YR-vs-FINT calibration split.

using Test, FVSjl

const _TI_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_ti_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 11 && tryparse(Int, first(split(l))) !== nothing]
_ti_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 11 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2200)]

@testset "TIMEINT cycle calendar vs Fortran" begin
    key = joinpath(_TI_DIR, "timeint10.key"); sav = joinpath(_TI_DIR, "timeint10.sum.save")
    if !isfile(key) || !isfile(sav)
        @test_skip "timeint10 scenario not available"
    else
        jl = _ti_rows(FVSjl.run_keyfile(key; faithful = true))
        ft = _ti_base(sav)
        @test length(jl) == length(ft)
        if length(jl) == length(ft)
            # year + age step by 10 EXACTLY (the calendar) — and must match Fortran's rows.
            for i in 1:length(jl)
                @test jl[i][1] == ft[i][1]                              # year (10-yr steps)
                @test jl[i][2] == ft[i][2]                              # age
            end
            # TPA bit-close (≤8); volume/BA within the ~2% calibrated-species residual.
            # The AUTCOR old-period/PVMLT-carry fix (diameter_growth.jl) + the d10n linear-G
            # recompute (mortality.jl:295) make the FIRST 10-yr cycle BIT-EXACT, and the cycle-2
            # 2010 TPA (360 vs 350) is now ALSO bit-exact after fixing the BAMAX/size-cap kill to
            # use the LINEAR FINT-extrapolated G = (DG/BARK)·(FINT/5) (morts.f:692/714/721) instead
            # of the sqrt fint-year diam_growth — the under-kill was that, NOT a WK3 DG tail.
            # Re-measured post-fix: the old 3% bounds were wildly stale. BA is BIT-EXACT every cycle; TPA
            # drifts ≤2; only cuft accumulates — to 16 (≈0.3%) by 2090. Same PROVEN growth-transcendental
            # class as the CS/SN grown-cycle envelope (docs/TOLERANCE_AUDIT.md, proven via DENSE-DEBUG AVH):
            # the DGF/HTGF Float32 exp/power leaves a sub-render per-cycle residual that is inert in DBH/BA
            # (both bit-exact) but compounds into the nonlinear cuft sum over the long run; the non-native
            # 10-yr cycle (SN calib is 5-yr) just amplifies the per-cycle transcendental step. Irreducible
            # without bit-matching FVS's libm exp/power. Bounds = the observed envelope.
            for i in 1:length(jl)
                @test abs(parse(Float64, jl[i][3]) - parse(Float64, ft[i][3])) <= 2   # TPA — mortality-timing (non-native cycle)
                @test parse(Float64, jl[i][4]) == parse(Float64, ft[i][4])            # BA  — BIT-EXACT
                @test abs(parse(Float64, jl[i][9]) - parse(Float64, ft[i][9])) <= 16  # cuft — deferred non-native DGSCOR tail (≈0.3%)
            end
        end
    end
end
