# test_fixmort.jl — FIXMORT forced-mortality override (morts.f:781) vs live Fortran.
#
# FIXMORT is a one-shot override applied AFTER the BA-check (the last word on the kill). In
# the cycle whose range holds the keyword date, for each tree of the chosen species (0=all)
# with d1≤DBH<d2 it sets the kill by PRM(5): 0 replace (P·rate), 1 add, 2 max, 3 multiply
# (kill·rate). Deterministic — no RNG. Each .sum.save is live-Fortran output.
#   * fixmort_replace — 0.5 replace, all DBH: ~50% of every record killed in the 1995 cycle.
#   * fixmort_mult    — 0.5 multiply (PRM(5)=3): halves the already-predicted mortality.
#   * fixmort_big     — 0.9 replace for DBH≥10" with PRM(6)=10 (KBIG=1, BOTTOM UP): the window's
#       mortality is pooled and re-imposed whole-record on the SMALLEST grown trees first
#       (morts.f:838 size concentration), so the large trees survive — BA/SDI stay high and QMD
#       climbs, unlike a flat per-record kill. Validates the XMORE-pool + ∓grown-DBH RDPSRT sort.
#   * fixmort_kpoint  — same window, PRM(6)=1 (KPOINT): the pooled mortality is concentrated point
#       by point (the 11-point base stand), killing whole records on the first points until XMORE
#       is spent (morts.f:937). A different survivor set than KBIG — a real, non-flat distribution.
#   * fixmort_kpbig   — PRM(6)=11 (KPOINT+KBIG): points have priority, smallest-first within each
#       point (morts.f:978) — the combined traversal.
# These exercise the fix that makes the recovery match: TPAMRT (the self-thinning line-reset
# carried to next cycle) is locked from the BA-check survivors BEFORE FIXMORT, so the forced
# kill doesn't move the self-thinning line (morts.f:772 precedes the FIXMORT block at 781).

using Test, FVSjl

const _FM_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_fm_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 11 && tryparse(Int, first(split(l))) !== nothing]
_fm_base(path) = [split(l) for l in eachline(path)
                  if length(split(l)) >= 11 &&
                     (y = tryparse(Int, first(split(l))); y !== nothing && 1900 < y < 2100)]
_fmcol(r, c) = parse(Float64, r[c])

@testset "FIXMORT forced mortality vs Fortran" begin
    have(nm) = isfile(joinpath(_FM_DIR, nm * ".key")) && isfile(joinpath(_FM_DIR, nm * ".sum.save"))
    runjl(nm) = (_fm_rows(FVSjl.run_keyfile(joinpath(_FM_DIR, nm * ".key"); faithful = true)),
                 _fm_base(joinpath(_FM_DIR, nm * ".sum.save")))

    for nm in ("fixmort_replace", "fixmort_mult", "fixmort_big", "fixmort_kpoint", "fixmort_kpbig")
        if !have(nm); @test_skip "$nm scenario not available"; continue; end
        @testset "$nm" begin
            jl, ft = runjl(nm)
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for i in 1:length(jl)
                    for c in (4, 7, 8)                      # BA / TopHt / QMD — BIT-EXACT (measured 0 all scenarios)
                        @test _fmcol(jl[i], c) == _fmcol(ft[i], c)
                    end
                    # TPA — BIT-EXACT (measured Δ0) for 4 of 5 scenarios; ONLY fixmort_kpoint lands the
                    # per-acre kill·rate TPA on the +0.5 render knife-edge (Δ1 at one cycle). So == for the
                    # bit-exact scenarios; the kpoint residual is exposed below the loop.
                    if nm != "fixmort_kpoint"
                        @test _fmcol(jl[i], 3) == _fmcol(ft[i], 3)
                    end
                end
                # fixmort_kpoint: the TPA 1-step render residual is the non-associative Float32 tree-SUM
                # accumulation order (doctrine #9: exposed as @test_broken, not a passing ≤1 hiding in green).
                if nm == "fixmort_kpoint"
                    @test_broken all(_fmcol(jl[i], 3) == _fmcol(ft[i], 3) for i in 1:length(jl))  # TPA — non-associative tree-SUM order
                end
            end
        end
    end

    # the replace override must visibly act: the 1995-cycle kill roughly halves TPA by 2000.
    if have("fixmort_replace")
        jl, ft = runjl("fixmort_replace")
        r2000 = findfirst(r -> r[1] == "2000", jl)
        @test r2000 !== nothing && _fmcol(jl[r2000], 3) < 0.65 * _fmcol(ft[r2000 - 1], 3)
    end
end
