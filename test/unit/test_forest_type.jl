# C3 — FIA forest-type classification (FORTYP + STKVAL, engine/forest_type.jl).
# snt01 classifies as 520 (mixed upland hardwoods), confirmed against live Fortran.
# The full 7-type validation (LP→161, SP→162, SA→142, WO→504, YP→803, RM→809, base
# →520) lives in test/harness against the Fortran oracle; here we pin the baseline
# plus the data-table shapes.

using Test
using FVSjl
using FVSjl: compute_forest_type!, coefficients, Southern

const SNT01 = "/workspace/ForestVegetationSimulator/tests/FVSsn/snt01.key"

@testset "FORTYP classifies snt01 as mixed upland hardwoods (520)" begin
    if isfile(SNT01)
        s, _ = initialize(SNT01)
        FVSjl.notre!(s)
        @test compute_forest_type!(s) == Int32(520)
        @test s.plot.forest_type == Int32(520)
    end
    c = coefficients(Southern())
    @test length(c.stock_b0) == 36 && length(c.stock_b1) == 36
    @test length(c.forest_type_codes) == 141
    @test Int32(520) in c.forest_type_codes
    @test haskey(c.fia_group, 802)        # white oak FIA code mapped to a group
end
