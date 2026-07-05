# test_hcor_calib.jl — small-tree height-growth calibration (REGENT HCOR_init) vs Fortran.
#
# The REGENT small-tree height calibration regression (regent.f:411-516) computes
# HCOR_init = ln( Σ(measured HTG·P) / Σ(predicted EDH·P) ) for species with ≥ NCALHT(5)
# small (dbh<5) trees that have measured height growth. No standard SN scenario carries
# such data, so this is a purpose-built stand: s29's 30 records plus 6 small SM (sp 22)
# trees with measured HTG, all on an existing plot (hcor_smalltree.key/.tre).
#
# Two things are validated, both bit-exact against live Fortran on this exact stand:
#  1. the regression value HCOR_init = -0.893823 (instrumented regent.f HCORF dump);
#  2. that HCOR then ATTENUATES per cycle (HCOR = WCI + cormlt_h·DIFH) — LDGCAL defaults
#     true for every species (grinit.f:102), so even an un-diameter-calibrated species
#     (WCI=0) decays as cormlt_h·HCOR_init rather than holding the value constant. The
#     end-to-end .sum reproduces Fortran's cycle-1 row exactly (TPA 5511, cuft 3764);
#     getting the attenuation wrong put cycle-1 TPA at 5564.

using Test, FVSjl

const _HCOR_DIR  = joinpath(@__DIR__, "..", "harness", "scenarios")
const _HCOR_KEY  = joinpath(_HCOR_DIR, "hcor_smalltree.key")
const _HCOR_TRE  = joinpath(_HCOR_DIR, "hcor_smalltree.tre")
const _HCOR_BASE = joinpath(_HCOR_DIR, "hcor_smalltree.sum.save")

_hcor_rows(txt) = [split(l) for l in split(txt, "\n")
                   if !occursin("-999", l) && length(split(l)) >= 11]

@testset "REGENT small-tree HCOR calibration vs Fortran" begin
    if !isfile(_HCOR_KEY) || !isfile(_HCOR_TRE)
        @test_skip "hcor_smalltree scenario not available"
    else
        # (1) the calibration regression value, bit-exact to Fortran
        s, _ = FVSjl.initialize(_HCOR_KEY)
        FVSjl.notre!(s); FVSjl.setup_growth!(s)
        @test isapprox(s.calib.htg_cor_init[22], -0.893823f0; atol = 2f-7)   # Float32-ULP: measured Δ=0 vs the 6-dec live stamp (was 1f-4, ~1000× padded)
        @test count(!=(0f0), s.calib.htg_cor_init) == 1   # only SM clears NCALHT(5)

        # (2) end-to-end: HCOR attenuation is applied correctly in small-tree growth.
        # Cycle 1 (1995) is bit-exact to Fortran; later cycles carry the usual ±2 tail.
        if isfile(_HCOR_BASE)
            jl = _hcor_rows(FVSjl.run_keyfile(_HCOR_KEY; faithful = true))
            ft = _hcor_rows(read(_HCOR_BASE, String))
            @test length(jl) == length(ft)
            if length(jl) == length(ft)
                for col in (3, 4, 5, 7, 9)   # TPA/BA/SDI/TopHt/total-cuft — ALL BIT-EXACT (measured Δ=0)
                    @test parse(Float64, jl[2][col]) == parse(Float64, ft[2][col])
                end
            end
        end
    end
end
