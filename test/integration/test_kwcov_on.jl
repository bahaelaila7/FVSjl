# test_kwcov_on.jl — Ontario (ON, metric) keyword-isolation coverage vs the LIVE oracle.
#
# test/fixtures/kwcov_on/on/on_<stand>_[ffe_]<case>.{key,live.sum,<Table>.csv,tables} (+ stands.db): two FIA stands from
# test/fixtures/tiered/on/stands.db (DB-read path) × one isolated keyword — base-model keywords (thinnings, multipliers,
# volume/defect, SETSITE, ESTAB, ...) and the FFE keywords the ON build accepts (every FFE report + *_Metric DBS table
# requested). Goldens = /workspace/.onwork/tiered/FVSon_tiered (test/harness/kwcov_on/make_cases.jl, PROVENANCE.toml). Every
# .sum field and every DBS table cell must match at print precision — no tolerance.
#   KWON_ONLY=regex(on stem)   restrict the run (default: all)
# A case that is not yet bit-exact goes in _KWON_OPEN with the measured diagnosis (reported as @test_broken so an unexpected
# pass surfaces); never silently skipped.
using Test, FVSjl
include(joinpath(@__DIR__, "..", "harness", "kwcov_ffe", "runner.jl"))

const _KWON_OPEN = Dict{String,String}(
)

@testset "ON keyword-isolation vs live" begin
    root = KCASE_ROOTS["on"]
    pat = Regex(get(ENV, "KWON_ONLY", "."))
    for stem in kcase_stems(root, "on")
        occursin(pat, stem) || continue
        @testset "$stem" begin
            ms, crashed, err = run_kcase(root, "on", stem)
            if haskey(_KWON_OPEN, stem)
                @test_broken isempty(ms)
            else
                isempty(ms) || foreach(m -> println("  MISMATCH ", key_str(m), " gold=", m.gold, " jl=", m.got), first(ms, 12))
                @test isempty(ms)
            end
        end
    end
end
