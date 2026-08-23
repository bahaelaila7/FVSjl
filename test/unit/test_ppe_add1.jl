# PPE (PPBASE) — ADD1 internal-stand-number digit increment.
#
# Bit-exact regression for ppe_add1 (increments the decimal digit of an 11-char CISN
# at character position icol+4; icol outside 1..7 is a no-op). PPE source recovered
# from FVS git history (deleted by bc6e2377, 2014); oracle is a gfortran-16 driver
# golden baked as literal expected strings.

using Test
using FVSjl
const _PPE = FVSjl

include(joinpath(@__DIR__, "..", "fixtures", "ppe", "ppe_add1_cases.jl"))  # defines ADD1

@testset "PPE ADD1 internal-stand-number digit increment — bit-exact vs gfortran-16 golden" begin
    @test !isempty(ADD1)
    for (cisn_in, icol, cisn_out) in ADD1
        @test _PPE.ppe_add1(cisn_in, icol) == cisn_out
    end
end
