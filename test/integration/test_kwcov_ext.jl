# test_kwcov_ext.jl — extension-STUB parity for CS/LS/NE/ON vs the LIVE oracles.
#
# test/fixtures/kwcov_ext/<v>/<v>_<case>.{key,tre|stands.db,live.sum,<Table>.csv,tables}: the canonical stand with ONE extension keyword block
# (MPB DFTM WSBW COVER BRUST MISTOE DFB RDIN BGCIN ECON CLIMATE ORGANON, + FERTILIZ), opener alone / with an unknown sub-keyword + END. The eastern
# builds link the ex*.f stubs (FVS11 for MPB DFTM COVER BRUST MISTOE DFB RDIN CLIMATE; ON also FERTILIZ and ECON; WSBW is the silent exbudl.f stub), so the
# run must equal the no-extension run and the FVS_Error rows (FVS11 / FVS01 per record) must match live. Goldens: test/harness/kwcov_ext/make_cases.jl
# (PROVENANCE.toml). Every .sum field and every DBS cell, no tolerance.   KWEXT_ONLY=regex(stem)   restrict the run.
using Test, FVSjl
include(joinpath(@__DIR__, "..", "harness", "kwcov_ffe", "runner.jl"))

const _KWEXT_OPEN = Dict{String,String}()   # all 116 cases bit-exact vs live (2026-10-06, summ.jl ext cs,ls,ne,on)

@testset "extension stubs vs live (CS/LS/NE/ON)" begin
    root = KCASE_ROOTS["ext"]
    pat = Regex(get(ENV, "KWEXT_ONLY", "."))
    for v in ("cs", "ls", "ne", "on"), stem in kcase_stems(root, v)
        occursin(pat, stem) || continue
        @testset "$stem" begin
            ms, crashed, err = run_kcase(root, v, stem)
            if haskey(_KWEXT_OPEN, stem)
                @test_broken isempty(ms)
            else
                isempty(ms) || foreach(m -> println("  MISMATCH ", key_str(m), " gold=", m.gold, " jl=", m.got), first(ms, 12))
                @test isempty(ms)
            end
        end
    end
end
