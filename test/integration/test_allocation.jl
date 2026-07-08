# test_allocation.jl — pillar-2 + pillar-4 hot-path guards, PER VARIANT (SN/NE/CS/LS).
# Locks the two performance-pillar metrics into the suite so a future refactor can't silently regress
# them (previously each was measured once, ad-hoc, and left unguarded).
#
# Pillar 2 (allocation-free): `@allocated grow_cycle!` on a warmed non-fire stand is a small, CONSTANT
# per-cycle cost — ALL of it the justified "Base sort scratch" (3× `compute_density!` descending-DBH
# stat sorts), NOT per-cycle Dict/Vector/comprehension churn. Measured floors (docs/MODERNIZATION_AUDIT
# S59/S100/S102, dead-stable spread=0): SN 15,392 · NE 10,640 · CS 10,048 · LS 9,984 B/cycle. The 20 KB
# ceiling is generous (covers SN with margin; catches a real regression — a reintroduced per-cycle
# temporary is thousands of bytes × tree count, far past 20 KB) while robust to warmup/version noise.
# If this fails HIGH: profile with `--track-allocation=user` before raising the ceiling.
#
# Pillar 4 (type-stable hot path): `@inferred` on the per-cycle ENTRY points — concrete inferred return,
# no `Any`/`Union` leaking out. Robust stdlib assertion (no JET/version flakiness, no new dependency).

const _HOTPATH_CASES = [
    ("snt01", FVSjl.Southern()),
    ("net01", FVSjl.Northeast()),
    ("cst01", FVSjl.CentralStates()),
    ("lst01", FVSjl.LakeStates()),
]

# Build a warmed, cyc0-ready stand and step it a few cycles (compile + reach steady-state buffers).
function _warm_hotpath_stand(key, variant)
    s = first(FVSjl.each_stand(key; variant = variant))
    FVSjl.notre!(s); FVSjl.setup_growth!(s); FVSjl.compute_volumes!(s)
    for _ in 1:4
        FVSjl.grow_cycle!(s; fint = 5f0)
    end
    return s
end

@testset "Pillar-2 — grow_cycle! per-cycle allocation floor (per variant)" begin
    for (nm, variant) in _HOTPATH_CASES
        key = joinpath(@__DIR__, "..", "fixtures", "canonical", "$nm.key")
        if !isfile(key)
            @test_skip "$nm.key fixture not available"
            continue
        end
        s = _warm_hotpath_stand(key, variant)
        allocs = [ @allocated FVSjl.grow_cycle!(s; fint = 5f0) for _ in 1:5 ]
        @test minimum(allocs) <= 20_000                       # generous per-cycle ceiling
        @test (maximum(allocs) - minimum(allocs)) <= 4_000    # constant per cycle (no growth-with-cycle churn)
    end
end

@testset "Pillar-4 — hot-path entry type-stability (per variant)" begin
    for (nm, variant) in _HOTPATH_CASES
        key = joinpath(@__DIR__, "..", "fixtures", "canonical", "$nm.key")
        if !isfile(key)
            @test_skip "$nm.key fixture not available"
            continue
        end
        s = _warm_hotpath_stand(key, variant)
        @test (@inferred FVSjl.grow_cycle!(s; fint = 5f0)) isa
              NamedTuple{(:accretion, :mortality)}            # concrete boundary return
        @test (@inferred FVSjl.compute_density!(s)) isa FVSjl.StandState
    end
end
