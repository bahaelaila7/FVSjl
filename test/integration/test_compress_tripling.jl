# =============================================================================
# test_compress_tripling.jl — COMPRESS does NOT suppress tripling in its own cycle.
#
# The divergence sweep flagged compress.key (`COMPRESS 1 15 50`) at ~50% Scuft/Bdft. Root cause: jl
# suppressed record tripling in the SAME cycle a COMPRESS fired, but FVS latches LTRIP at cycle start
# (grincr.f:74) BEFORE COMCUP (:391) sets NOTRIP=.TRUE. — so the compress cycle STILL triples the merged
# records (15 → 45), and NOTRIP suppresses tripling only from the NEXT cycle. Without the triple the
# merged records land coarsely relative to the sawtimber DBH threshold, halving Scuft/Bdft.
#
# FIX: apply_compress! sets control.no_tripling (persists to later cycles); the tripling gate reads
# no_tripling CAPTURED BEFORE the compress. Validated: compress.key 1995 (the compress cycle) Scuft/Bdft
# become BIT-EXACT vs live FVSsn (253/1040). Later cycles carry the accepted COMPRESS eigensolver/merge-
# order residual (~3-4%), so only the compress-cycle row is asserted bit-exact here.
# =============================================================================

using Test
using FVSjl

@testset "COMPRESS still triples its own cycle (merch volume) vs live FVSsn" begin
    key = joinpath(@__DIR__, "..", "harness", "scenarios", "compress.key")
    if !isfile(key)
        @test_skip "compress.key not available"
    else
        txt = FVSjl.run_keyfile(key; variant = FVSjl.Southern(), output = :sum)
        row = nothing
        for ln in split(txt, '\n')
            t = split(strip(ln))
            length(t) >= 12 && t[1] == "1995" && (row = t; break)
        end
        @test row !== nothing
        if row !== nothing
            @test parse(Int, row[3])  == 496    # TPA  (live 496) — compress+triple keeps the density right
            @test parse(Int, row[4])  == 104    # BA   (live 104)
            @test parse(Int, row[5])  == 213    # SDI  (live 213)
            @test parse(Int, row[11]) == 253    # Scuft (live 253) — was 125 (~50% low) before the fix
            @test parse(Int, row[12]) == 1040   # Bdft  (live 1040) — was 566 before the fix
        end
    end
end
