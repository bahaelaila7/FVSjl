# =============================================================================
# test_opt01_cyc0.jl — OP (Olympic) chunk-3 cyc0 `.sum` BIT-EXACT vs live FVSop_clean.
#
# Reference: /workspace/.opwork/opdbg.sum (FVSop_clean, stand S248112, opdbg.key). The 1990
# inventory-year row is the deterministic (DGSD=0, ICL4=0) initial-stand summary — it depends
# only on the CRATET missing-height dub (op/htdbh.f MODE=0 for the FVS-native records; the
# stand's big-6 ORGANON trees all carry inventory heights so PREPARE dubs nothing) and the
# BLM Behre volume + op/ccfcal.f CCF. Multi-cycle is NOT validated here: FVSop_clean SIGSEGVs
# oracle-side in the FVS-native fvsvol path on cycle ≥1, and the jl per-cycle DGDRIV coexistence
# driver (ORGANON-NWO / op-Wykoff) is a follow-on.
#
# Oracle 1990 row (opdbg.sum):
#   1990  60  536  77  184  114  63  5.1  1472  972  0  5003  ...
# All volume + density columns are reproduced bit-exact. The FORTYP / size-stocking classes
# (999 / 5 / 5 vs the oracle's 201 / 2 / 3) are the SAME non-volume classification delta declared
# out-of-scope for the sibling ORGANON variant OC (fortypv/strclass), so they are not asserted.
# =============================================================================
using Test
using FVSjl
const F = FVSjl

@testset "OP opt01 cyc0 .sum bit-exact (S248112 vs FVSop_clean)" begin
    key = joinpath(@__DIR__, "..", "harness", "scenarios", "opt01.key")
    if !isfile(key)
        @test_skip "opt01 scenario not available"
    else
        local r = nothing
        for s in F.each_stand(key; variant = F.Olympic(), faithful = true)
            F.notre!(s)
            F.setup_growth!(s)
            F.compute_volumes!(s)
            r = F.summary_row(s)
            break
        end
        @test r !== nothing
        # density + stand-structure columns
        @test r.year  == 1990
        @test r.age   == 60
        @test r.tpa   == 536
        @test r.ba    == 77
        @test r.sdi   == 184
        @test r.ccf   == 114
        @test r.topht == 63
        @test round(r.qmd; digits = 1) == 5.1
        # the four stand volumes (TCuFt / MCuFt / SawCuFt / BdFt)
        @test r.cuft  == 1472
        @test r.mcuft == 972
        @test r.scuft == 0
        @test r.bdft  == 5003
    end
end
