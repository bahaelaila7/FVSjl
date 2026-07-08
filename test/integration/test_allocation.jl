# test_allocation.jl — pillar-2 guard: the per-cycle hot path (`grow_cycle!` and everything it
# transitively calls) must stay allocation-free down to the documented+justified floor. This LOCKS
# that floor into the suite so a future refactor can't silently regress it (previously the S59
# metric was measured once, ad-hoc, and left unguarded).
#
# Documented floor (docs/MODERNIZATION_AUDIT.md, S59/S99): `@allocated grow_cycle!` on a warmed
# non-fire stand is a small, constant per-cycle cost — ALL of it the justified "Base sort scratch"
# (3× `compute_density!` descending-DBH stat sorts), NOT per-cycle Dict/Vector/comprehension churn.
# Measured net01 NE = 10,656 B/cycle, dead-stable across calls. The CEILING below is generous
# (catches a real ≥1.5× regression — e.g. a reintroduced per-cycle temporary — while staying robust
# to Julia-version / warmup noise). If this fails HIGH, something started allocating per cycle;
# profile with `--track-allocation=user` before raising the ceiling.

@testset "Pillar-2 — grow_cycle! per-cycle allocation floor (net01 NE)" begin
    key = joinpath(@__DIR__, "..", "fixtures", "canonical", "net01.key")
    if !isfile(key)
        @test_skip "net01.key fixture not available"
    else
        s = first(FVSjl.each_stand(key; variant = FVSjl.Northeast()))
        FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
        for _ in 1:4                      # warm (compile + steady-state buffers)
            FVSjl.grow_cycle!(s; fint = 5f0)
        end
        allocs = [ @allocated FVSjl.grow_cycle!(s; fint = 5f0) for _ in 1:5 ]
        floor_bytes = minimum(allocs)     # min = GC-noise-free steady state
        # Guard the documented per-cycle floor (10,656 B). Generous ceiling catches real regressions.
        @test floor_bytes <= 16_000
        # And it must be a small CONSTANT per cycle (no growth-with-cycle churn) — spread stays tight.
        @test (maximum(allocs) - minimum(allocs)) <= 4_000
    end
end
