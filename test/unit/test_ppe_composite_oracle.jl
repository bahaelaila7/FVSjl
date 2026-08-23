# PPE COMPOSITE YIELD aggregation — oracle-baked, oracle-free at test time.
#
# src/engine/ppe_landscape.jl:_ppe_aggregate is the CMADDS/CMPRT2 composite-yield
# aggregation (ppbase/cmadds.f + cmprt2.f) — the routine that builds FVSppe's
# "COMPOSITE YIELD STATISTICS" table from the per-stand summary rows. This test bakes
# the historical FVSppe oracle's OWN per-stand yield tables + composite (a 4-stand East
# Cascades 15-cycle run, weights 1/11/11/11, TOTAL SAMPLE WEIGHT 34; ecpp.{key,tre,out}
# rebuilt from FVS rev bc6e2377^) as literals and asserts _ppe_aggregate reproduces the
# oracle composite BIT-EXACT for every area-weighted column.
#
# This isolates the ORCHESTRATION (SPLAEX/CMADDS area-weighting) from per-stand growth:
# it feeds the oracle's per-stand values straight into the aggregator, so a match proves
# the aggregation semantics are exactly the oracle's, independent of FVSjl's EC growth
# model. (End-to-end FVSjl-projected stands are covered by test_ppe_landscape; the EC
# per-stand growth itself is validated/cornered elsewhere.)
#
# The aggregation rules confirmed against the oracle (see cmadds.f loop 220 + 15/16):
#   TREES, TOTAL/MERCH CU FT, MERCH BD FT (IOSUM 3..6): Σ(v·w)/Σw            (÷ PRBSUM)
#   ACCRETION, MORTALITY (IOSUM 15,16): Σ(v·w·prd)/Σ(w·prd)                  (÷ Σ w·prd)
#   PERIOD (IOSUM 14): Σ(prd·w)/Σw ;  TOTAL SAMPLE WEIGHT (IOSUM 17): Σw
# Composite rounding is IFIX(x+.5) (round-half-up on the non-negative totals).

using Test
using FVSjl

include(joinpath(@__DIR__, "..", "fixtures", "ppe", "ppe_composite_oracle.jl"))

@testset "PPE CMADDS/CMPRT2 composite-yield aggregation (FVSppe oracle-baked)" begin
    ifix(x) = floor(Int, x + 0.5)          # FVS IFIX(x+.5) round-half-up (non-negative)

    # Feed the oracle's per-stand rows straight into the aggregator (weight == area).
    stand_rows = [(Float64(w), rows) for (w, rows) in PPE_ORACLE_STANDS]
    agg = FVSjl._ppe_aggregate(stand_rows)

    byyear = Dict(g.year => g for g in agg)
    @test length(agg) == length(PPE_ORACLE_COMPOSITE)          # same set of report years

    for c in PPE_ORACLE_COMPOSITE
        g = byyear[c.year]
        # area-weighted means (÷ Σw): TREES + the three before-thin volume columns
        @test ifix(g.avbtpa)   == c.trees
        @test ifix(g.avbtcuft) == c.cuft
        @test ifix(g.avbmcuft) == c.mcuft
        @test ifix(g.avbbdft)  == c.bdft
        # area·period-weighted means (÷ Σ w·prd): accretion + mortality
        @test ifix(g.avbacc)   == c.acc
        @test ifix(g.avbmort)  == c.mort
        # PERIOD = Σ(prd·w)/Σw ; TOTAL SAMPLE WEIGHT = Σw
        @test g.msperiod       == c.prd
        @test ifix(g.totalwt)  == c.weight
    end

    # Σw over ALL stands = the oracle's TOTAL SAMPLE WEIGHT (34) on every reported year.
    @test all(ifix(g.totalwt) == 34 for g in agg)
    # every stand contributed to every year (4 stands report all 16 years).
    @test all(g.nstands == length(PPE_ORACLE_STANDS) for g in agg)
end
