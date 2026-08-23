# PPE (Parallel Processing Extension) — PPBASE character index QuickerSort.
#
# Bit-exact regression for ppe_index_qsort! (a faithful @goto transliteration of
# PPBASE C11SRT/C26SRT/CH8SRT, Scowen 1965 Algorithm 271), used by PPMAIN/BYGRPS for
# the master all-stand ordering. The PPE source was recovered from the FVS Fortran
# repo's git history (deleted by commit bc6e2377, 2014); the full FVSppe exe won't
# link against modern base-FVS (getstd/putstd COMMON drift), so the oracle is a
# per-routine gfortran-16 driver golden. The fixture bakes that golden as literal
# expected index vectors (oracle-free at test time).
#
# The physical key array is never permuted — only INDEX is rearranged so that
# keys[index[i]] <= keys[index[i+1]]. Ties are placed by the partition's swap
# sequence (NOT stable), so the exact index vector (not just sortedness) is the
# bit-exact contract — the tie-sensitive cases (all_equal, *_revidx) lock that.

using Test
using FVSjl
const _PPE = FVSjl

include(joinpath(@__DIR__, "..", "fixtures", "ppe", "ppe_sort_cases.jl"))  # defines CASES

@testset "PPE C11SRT character index QuickerSort — bit-exact vs gfortran-16 golden" begin
    @test !isempty(CASES)
    for (name, keys, lseq, idx, gold) in CASES
        n = length(keys)
        index = idx === nothing ? zeros(Int, n) : copy(idx)
        _PPE.ppe_index_qsort!(index, keys; lseq = lseq)
        @test index == gold                       # exact index vector (incl. tie order)
        # sanity: the resulting index really orders the keys ascending
        @test issorted([keys[index[i]] for i in 1:n])
    end
end
