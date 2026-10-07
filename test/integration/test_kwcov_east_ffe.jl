# test_kwcov_east_ffe.jl — eastern (CS/LS/NE) FFE keyword-isolation coverage vs the LIVE oracle.
#
# test/fixtures/kwcov_ffe/<v>/<v>_<case>.{key,tre,live.sum,<Table>.csv,tables}: one stand × one isolated FFE keyword
# (FLAMEADJ, MOISTURE, FUELMOVE, PILEBURN, FUELTRET, POTF*, SNAG*, SALVAGE, DROUGHT, FUELINIT, FUELFOTO, CARBCALC,
# FMODLIST, FUELMODL, DEFULMOD, FIRECALC, SOILHEAT, FMORTMLT, CANCALC, …) with every FFE report + DBS table requested.
# Goldens come from the stock-flag LIVE oracle (test/harness/kwcov_ffe/make_cases.jl, PROVENANCE.toml). Every .sum field
# and every DBS table cell must match at print precision — no tolerance.
#   KWFFE_VARIANTS=cs,ls   KWFFE_ONLY=regex(on stem)   restrict the run (default: all)
# A case that is not yet bit-exact goes in _KWFFE_OPEN with the measured diagnosis (reported as @test_broken so an
# unexpected pass surfaces); never silently skipped.
using Test, FVSjl
include(joinpath(@__DIR__, "..", "harness", "kwcov_ffe", "runner.jl"))

# No OPEN entries (2026-10-07): every case is bit-exact vs the live oracle.
const _KWFFE_OPEN = Dict{String,String}(
)

@testset "eastern FFE keyword-isolation vs live" begin
    root = KCASE_ROOTS["ffe"]
    vs = split(get(ENV, "KWFFE_VARIANTS", "cs,ls,ne"), ",")
    pat = Regex(get(ENV, "KWFFE_ONLY", "."))
    for v in vs, stem in kcase_stems(root, v)
        occursin(pat, stem) || continue
        @testset "$stem" begin
            ms, crashed, err = run_kcase(root, v, stem)
            if haskey(_KWFFE_OPEN, stem)
                @test_broken isempty(ms)
            else
                isempty(ms) || foreach(m -> println("  MISMATCH ", key_str(m), " gold=", m.gold, " jl=", m.got), first(ms, 12))
                @test isempty(ms)
            end
        end
    end
end
