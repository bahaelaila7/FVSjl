# test_tiered.jl — FAST TIER of the tiered integration suite (test/harness/tiered/README.md).
#
# For every committed fixture variant (test/fixtures/tiered/<v>/), every stand × regime: run jl on the fixture key
# and compare against the committed LIVE-ORACLE goldens at print precision (.sum: every printed field of every row;
# DBS: every table present/absent + every column of every row, numeric-exact at FVS's stored precision, text exact).
# Plus a signed tally per variant × regime (one-directional material bias ⇒ fail).
#
# Rules (no silent passes): an unlisted mismatch FAILS; a KNOWN_RESIDUALS entry that matches nothing in a full run
# FAILS ("unexpected pass — remove it"); an entry without status/proof/issue FAILS; a missing fixture FAILS (never a
# skip); a jl crash FAILS.
#   TIERED=quick            → first 3 stands × QUICK_REGIMES per variant (refactor iteration, a few minutes)
#   TIERED_VARIANTS=BM,SN   → restrict variants
using Test, FVSjl
include(joinpath(@__DIR__, "..", "harness", "tiered", "tiered_runner.jl"))

@testset "tiered fast tier vs live-oracle goldens" begin
    quick = get(ENV, "TIERED", "") == "quick"
    vs = tiered_variants()
    @test !isempty(vs)                                   # fixtures must exist — never a silent skip
    al, alerrs = load_allowlist()
    @test isempty(alerrs)
    foreach(e -> @warn("allowlist: $e"), alerrs)
    all_ms, ran, percase = run_fast_tier(; quick, threads = parse(Int, get(ENV, "TIERED_THREADS", "1")))
    summary = Dict{Tuple{String,String},Vector{Int}}()    # (v,r) => [cases, all-outputs exact, .sum exact, with residual]
    sumbad = Set((m.variant, m.stand, m.regime) for m in all_ms if m.file in ("sum", "CRASH"))
    for ((v, cn, r), ok) in percase
        s = get!(summary, (v, r), [0, 0, 0, 0]); s[1] += 1; ok && (s[2] += 1)
        (v, cn, r) in sumbad || (s[3] += 1); ok || (s[4] += 1)
    end
    unlisted, stale, used, over = apply_allowlist(all_ms, al, ran)
    # compact summary: variant × regime → cases / exact / with residual(s)
    println("\n== tiered fast tier ($(quick ? "quick" : "full")) ==")
    println(rpad("variant/regime", 22), lpad("cases", 6), lpad("all-exact", 10), lpad(".sum-exact", 11), lpad("resid", 7))
    for k in sort(collect(keys(summary)))
        s = summary[k]; println(rpad("$(k[1])/$(k[2])", 22), lpad(s[1], 6), lpad(s[2], 10), lpad(s[3], 11), lpad(s[4], 7))
    end
    ncor = sum(used[i] for (i, e) in enumerate(al) if e.status == "CORNER"; init = 0)
    nopn = sum(used[i] for (i, e) in enumerate(al) if e.status == "OPEN"; init = 0)
    println("mismatch cells: $(length(all_ms))  allowlisted CORNER: $ncor  OPEN: $nopn  UNLISTED: $(length(unlisted))",
            quick ? "" : "  stale allowlist entries: $(length(stale))")
    for m in first(unlisted, 40)
        println("  UNLISTED ", key_str(m), "  gold=", m.gold, "  jl=", m.got)
    end
    length(unlisted) > 40 && println("  … $(length(unlisted) - 40) more")
    @test isempty(unlisted)
    for (e, n) in over
        println("  RESIDUAL SPREAD — allowlist entry #$(e.idx) $(e.variant)/$(e.regime)/$(e.file)/$(e.col): $n cells > max_cells $(e.max_cells)")
    end
    @test isempty(over)
    if !quick
        for (i, e) in enumerate(al)
            0 < used[i] < e.max_cells < typemax(Int) && println("  improved (tighten max_cells): entry #$(e.idx) $(e.variant)/$(e.regime)/$(e.file)/$(e.col) $(used[i]) < $(e.max_cells)")
        end
        for e in stale
            println("  UNEXPECTED PASS — remove allowlist entry #$(e.idx): $(e.variant)/$(e.stand)/$(e.regime)/$(e.file)/$(e.col)/$(e.year)")
        end
        @test isempty(stale)
    end
end
