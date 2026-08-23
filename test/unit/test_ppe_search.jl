# PPE (PPBASE) — keyed binary search over a character array via an ascending index.
#
# Bit-exact regression for ppe_c26bsr (C26BSR/CH8BSR — position in A) and ppe_spbsrx
# (SPBSRX — position in IORD). These pair with ppe_index_qsort!: IORD is an ascending
# index over A. PPE source recovered from FVS git history (deleted by bc6e2377, 2014);
# the full FVSppe exe won't link against modern base-FVS, so the oracle is a
# per-routine gfortran-16 driver golden, baked here as literal expected returns.

using Test
using FVSjl
const _PPE = FVSjl

include(joinpath(@__DIR__, "..", "fixtures", "ppe", "ppe_search_cases.jl"))  # defines SEARCH

@testset "PPE C26BSR/SPBSRX keyed binary search — bit-exact vs gfortran-16 golden" begin
    @test !isempty(SEARCH)
    for (keys, iord, f, ip, im) in SEARCH
        @test _PPE.ppe_c26bsr(keys, iord, f) == ip     # subscript in A (0 if absent)
        @test _PPE.ppe_spbsrx(keys, iord, f) == im     # subscript in IORD (0 if absent)
        # cross-check: when found, iord[imid] must equal ip
        im != 0 && @test iord[im] == ip
    end
end
