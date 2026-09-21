# report.jl — write EVERY fast-tier mismatch (and tally flags) for the selected variants to a TSV, for triage and
# for populating KNOWN_RESIDUALS.toml from measurement (never by hand-waving).
# Usage: TIERED_VARIANTS=BM julia --project=. test/harness/tiered/report.jl <out.tsv>
include(joinpath(@__DIR__, "tiered_runner.jl"))

function report(out)
    all_ms, _, _ = run_fast_tier(; threads = parse(Int, get(ENV, "TIERED_THREADS", "1")))
    open(out, "w") do io
        println(io, "variant\tstand\tregime\tfile\tcol\tyear\tgold\tjl")
        for m in all_ms
            println(io, join((m.variant, m.stand, m.regime, m.file, m.col, m.year,
                              replace(m.gold, r"\s+" => " "), replace(m.got, r"\s+" => " ")), '\t'))
        end
    end
end
report(ARGS[1])
