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
            # BA is the structural contract (the TIMEINT ×2 growth scaling): BIT-EXACT every cycle.
            # The AUTCOR old-period/PVMLT-carry fix (diameter_growth.jl) + the d10n linear-G recompute
            # (mortality.jl:295) + the BAMAX/size-cap LINEAR FINT-extrap G = (DG/BARK)·(FINT/5)
            # (morts.f:692/714/721) made the first 10-yr cycles bit-exact. So BA + year + age assert ==.
            for i in 1:length(jl)
                @test parse(Float64, jl[i][4]) == parse(Float64, ft[i][4])            # BA — BIT-EXACT every cycle
            end
            # TPA + cuft: NON-NATIVE 10-yr cycle residual. ★ VERDICT CORRECTED (2026-07-06, live FVS_TreeList DBS
            # NOTRIPLE, matched by TreeId): the "DGF/HTGF transcendental cubic tail" attribution was WRONG — jl's
            # per-tree DBH AND crown are BIT-EXACT at every checked cycle (2000/2010/2020: 0/27, max|Δ|=0). So the
            # growth path (incl. the 10-yr FINT/YR scaling) is faithful. The real residual is a BA-CONSERVING
            # MORTALITY KNIFE-EDGE: NOTRIPLE .sum shows TPA Δ≤1 + TCuFt Δ≤2 from 2040 on (BA bit-exact), i.e. the
            # 10-yr mortality flips a near-tie kill (jl kills a slightly different small tree), and TRIPLING amplifies
            # it into the test's cuft tail. SAME class as cs_allsp (task 68) — the dense/non-native-cycle mortality-
            # distribution knife-edge. Seed = a sub-ULP in the 10-yr tokill / size-cap kill / VARMRT efftr (per-tree
            # DBH bit-exact rules out the growth transcendentals). ★ TRACED TO GROUND (2026-07-06, live morts.f
            # stamp of T/DIA0/D10): cycle-1 T (Σtpa, species-sorted) is BIT-EXACT (589.6527709961 both) — NOT the
            # sum-order — but D10 (the grown-diameter self-thinning QMD) differs by ~5e-7 = ~1 Float32 ULP; that
            # shifts the self-thinning line → flips a near-tie kill → cascades (cyc2 T Δ2e-4). D10=fpow(sumdr10/tt,
            # 1/1.605) with ^1.605 fpow-routed (inert) + tt bit-exact ⇒ the 1-ULP is in the MORTS grown-diameter G.
            # EXACT MECHANISM (2026-07-06): FVS morts.f:225 uses G=(DG_5yr/BARK)·(FINT/5) — a linear scale of the
            # NATIVE 5-yr DG. jl doesn't store the native DG (it keeps the fint-year applied growth), so for a
            # non-native FINT (=10 here, YR=5) `_mort_traj_g` RECONSTRUCTS the 5-yr DG via a squaring+sqrt roundtrip
            # ((dg+dib)²→·yr/fint→sqrt), which loses ~1 ULP vs FVS's actual native DG → D10 ~1-ULP → the self-thin
            # knife-edge. NOT a transcendental (sqrt is IEEE-exact) — a SEMANTIC reconstruction diff. Native-cycle
            # scenarios use the identity (`_mort_traj_g` line 2) and are BIT-EXACT (all realistic canonicals).
            # ★ "FIXABLE by threading the native DG" — REFUTED (task #70, 2026-07-06). I stashed the clean bounded
            # native 5-yr DG at growth time and fed it straight to MORTS (skipping the roundtrip). Native cycles
            # stayed bit-exact (identity), but the NON-native cycles REGRESSED: BA line 43 (bit-exact today!) drifted
            # to 129 vs live 127, s5 cuft 3149 vs 3111, s9 15752 vs 15426. So live FVS's own Float32 op-sequence
            # matches THIS reconstruction, NOT the clean native DGb — FVS's MORTS recovers the linear increment from
            # the FINT-scaled DG the same squaring/sqrt way (equal in ℝ, ~1 ULP apart in Float32; the knife-edge
            # flips on that ULP). ⇒ the reconstruction is the CORRECT/faithful path; the residual below is purely
            # the mortality knife-edge (a sub-ULP self-thinning-line near-tie flip), a PERMITTED cornered primitive,
            # NOT a fixable reconstruction artifact. @test_broken with the exact mechanism named, not a padded bound.
            @test_broken all(parse(Float64, jl[i][3]) == parse(Float64, ft[i][3]) for i in 1:length(jl))  # TPA — non-native mortality-timing tail
            @test_broken all(parse(Float64, jl[i][9]) == parse(Float64, ft[i][9]) for i in 1:length(jl))  # cuft — non-native DGSCOR/transcendental tail
        end
    end
end
