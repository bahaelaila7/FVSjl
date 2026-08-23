# PPE (PPBASE) — HXINDX hexagonal-grid neighbor indexing.
#
# Bit-exact regression for ppe_hxindx, the spatial-neighbor core of PPE's inter-stand
# disturbance spread (called from HVNEDN). PPE source recovered from FVS git history
# (deleted by bc6e2377, 2014); oracle is a gfortran-16 driver golden baked as literal
# expected pointers. Fixture enumerates every neighbor of every stand in small grids
# plus Float32-sqrt perfect-square boundaries (NR=IFIX(SQRT(FLOAT(NSTND))) is REAL*4,
# so the row/col count must truncate identically — the boundary cases lock that).

using Test
using FVSjl
const _PPE = FVSjl

include(joinpath(@__DIR__, "..", "fixtures", "ppe", "ppe_hxindx_cases.jl"))  # defines HXINDX

@testset "PPE HXINDX hexagonal neighbor indexing — bit-exact vs gfortran-16 golden" begin
    @test !isempty(HXINDX)
    for (nei, istnd, nstnd, neiptr) in HXINDX
        @test _PPE.ppe_hxindx(nei, istnd, nstnd) == neiptr
    end
end
