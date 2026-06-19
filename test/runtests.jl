# Top-level test entry. Run with: julia --project -e 'using Pkg; Pkg.test()'
#
# Tests are organized as:
#   unit/        — pure-kernel + struct tests (fast, no oracle)
#   integration/ — full-run parity vs Oracle A (snt01/snt02/sndb/sntest)
#
# As each chunk lands, its tests are added here. Integration tests are guarded so
# the suite stays green while the engine is still being built (they activate once
# the corresponding chunk is done).

using Test
using FVSjl

@testset "FVSjl" begin
    include("unit/test_core.jl")           # C0: state, rng, units, variant
    include("unit/test_species.jl")        # C2: SN species tables + blkdat defaults
    include("integration/test_treedata.jl")# C1: .tre parser vs Oracle A
    include("integration/test_keyword.jl") # C1: keyword lexer vs Oracle A
    include("integration/test_io_formats.jl")# C1b: CSV/format-agnostic round-trips
    include("integration/test_init.jl")    # C2: keyword dispatch + tree loading
    # include("integration/test_snt01.jl")  # enabled at C5
    # include("integration/test_sndb.jl")   # enabled at C6
    # include("integration/test_snt02.jl")  # enabled at C8
end
