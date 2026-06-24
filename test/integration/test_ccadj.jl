# test_ccadj.jl — CCADJ keyword (opt 444, crown-competition-factor adjustment).
#
# CCADJ sets CCCOEF (sstage.f:923 UPDATECCCOEF). CCCOEF is read ONLY by covolp.f (the
# COVER canopy-cover report), sstage.f (the SSTAGE structure-stage code), and evldx.f:430
# (the cover event-monitor variable) — never by the core growth / mortality / density. So
# CCADJ is `.sum`-inert in SN (verified byte-identical vs live Fortran, only the run
# timestamp differs) and is a recognized no-op until COVER/SSTAGE output is ported (C6).
# Checks:
#   1. RECOGNIZED — CCADJ dispatches as a no-op (no "unknown keyword" error / no crash);
#   2. INERT — injecting `CCADJ` leaves every projected `.sum` row unchanged.

using Test, FVSjl

const _CC_DIR = joinpath(@__DIR__, "..", "harness", "scenarios")
_cc_rows(txt) = [split(l) for l in split(txt, "\n")
                 if length(split(l)) >= 12 && (y = tryparse(Int, first(split(l)));
                                               y !== nothing && 1900 < y < 2100)]

@testset "CCADJ → recognized .sum-inert no-op (SN)" begin
    # 1. RECOGNIZED — CCADJ is in the no-op set for SN (universal ∪ variant).
    s = FVSjl.StandState(FVSjl.Southern())
    @test ("CCADJ" in FVSjl.KNOWN_NOOP) ||
          ("CCADJ" in FVSjl.variant_noop_keywords(s.variant))

    key = joinpath(_CC_DIR, "fire_early.key")
    if !isfile(key)
        @test_skip "fire_early scenario not available"
    else
        base = _cc_rows(FVSjl.run_keyfile(key; faithful = true))

        # 2. INERT — inject a CCADJ line after NUMCYCLE; the .sum must be unchanged. The
        # variant key must share the .tre basename, so write it beside the scenario's .tre.
        cckey = joinpath(_CC_DIR, "_ccadj_on.key")
        cp(joinpath(_CC_DIR, "fire_early.tre"), joinpath(_CC_DIR, "_ccadj_on.tre"); force = true)
        open(cckey, "w") do io
            for l in eachline(key)
                println(io, l)
                startswith(strip(l), "NUMCYCLE") && println(io, "CCADJ             0.5")
            end
        end
        try
            on = _cc_rows(FVSjl.run_keyfile(cckey; faithful = true))
            @test length(on) == length(base)
            @test all(on[i] == base[i] for i in eachindex(base))   # byte-identical projection
        finally
            rm(cckey; force = true); rm(joinpath(_CC_DIR, "_ccadj_on.tre"); force = true)
        end
    end
end
