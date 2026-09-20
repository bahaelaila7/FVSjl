# coverage_gate.jl — COVERAGE GATE of the tiered suite.
#
# Measures which src FUNCTIONS the tiered fixtures exercise, per variant (fast-tier cases run in crash-isolated
# workers under --code-coverage=user; .cov files aggregated then deleted). A function is identified by
# (file, name, ordinal-of-that-name-in-the-file) — stable across edits that shift line numbers. Function extents are
# found heuristically: a top-level (column-1) `function NAME` / `NAME(...) =` / `@inline function NAME` line opens a
# function that runs to the next top-level definition. A function is covered for a variant if any of its lines ran.
#
#   julia --project=. test/harness/tiered/coverage_gate.jl baseline
#       measure + write test/fixtures/tiered/COVERAGE_BASELINE.tsv (covered functions, by variant) and
#       test/fixtures/tiered/COVERAGE_GAPS.tsv (zero-coverage functions — the tracked known gap list) +
#       docs/TIERED_COVERAGE_<date>.md (summary)
#   julia --project=. test/harness/tiered/coverage_gate.jl check
#       measure + FAIL (exit 1) if a baseline-covered function is now uncovered (lost coverage), or a zero-coverage
#       function is not listed in COVERAGE_GAPS.tsv (new untested code); prints functions newly covered (tighten).
# Env: TIERED_VARIANTS, TIERED=quick (quick subset — only for local iteration; the gate uses the full tier).
include(joinpath(@__DIR__, "tiered_runner.jl"))

const SRC = normpath(joinpath(@__DIR__, "..", "..", "..", "src"))
const BASE_P = joinpath(TIERED_ROOT, "COVERAGE_BASELINE.tsv")
const GAPS_P = joinpath(TIERED_ROOT, "COVERAGE_GAPS.tsv")
const _DEF = r"^(?:@inline\s+|@noinline\s+|Base\.@propagate_inbounds\s+)?(?:function\s+([A-Za-z_][\w!.]*(?:\{[^}]*\})?)|([A-Za-z_][\w!.]*)(?:\{[^}]*\})?\([^=]*\)\s*(?:::\s*[\w{},. ]+)?\s*(?:where\s+[^=]+)?=(?!=))"

"(file_rel, name, ordinal, first_line, last_line) for every top-level function in src."
function src_functions()
    out = Tuple{String,String,Int,Int,Int}[]
    for (root, _, files) in walkdir(SRC), f in files
        endswith(f, ".jl") || continue
        p = joinpath(root, f); rel = relpath(p, SRC); lines = readlines(p)
        starts = Tuple{Int,String}[]
        for (i, l) in enumerate(lines)
            m = match(_DEF, l); m === nothing && continue
            nm = something(m.captures[1], m.captures[2]); nm = replace(nm, r"\{.*" => "")
            push!(starts, (i, nm))
        end
        seen = Dict{String,Int}()
        for (k, (i, nm)) in enumerate(starts)
            last = k < length(starts) ? starts[k+1][1] - 1 : length(lines)
            seen[nm] = get(seen, nm, 0) + 1
            push!(out, (rel, nm, seen[nm], i, last))
        end
    end
    out
end

"Line hit counts per src file from all *.cov files (summed over workers), then delete them."
function harvest_cov()
    hits = Dict{String,Vector{Int}}()
    for (root, _, files) in walkdir(SRC), f in files
        occursin(r"\.jl\.\d+\.cov$", f) || continue
        p = joinpath(root, f); rel = relpath(replace(p, r"\.\d+\.cov$" => ""), SRC)
        for (i, l) in enumerate(readlines(p))
            c = tryparse(Int, strip(first(l, 9)))
            v = get!(hits, rel, Int[])
            while length(v) < i; push!(v, 0); end
            c === nothing || (v[i] += c)
        end
        rm(p)
    end
    hits
end

function clean_cov()
    for (root, _, files) in walkdir(SRC), f in files
        occursin(r"\.jl\.\d+\.cov$", f) && rm(joinpath(root, f))
    end
end

fid(t) = "$(t[1])::$(t[2])#$(t[3])"

"Set of function ids covered, per variant."
function measure()
    ENV["TIERED_COVERAGE"] = "1"
    fns = src_functions(); res = Dict{String,Set{String}}()
    for v in tiered_variants()
        clean_cov()
        withenv("TIERED_VARIANTS" => v) do
            run_fast_tier(; threads = 1)
        end
        hits = harvest_cov(); cov = Set{String}()
        for t in fns
            h = get(hits, t[1], Int[])
            any(i -> i <= length(h) && h[i] > 0, t[4]:t[5]) && push!(cov, fid(t))
        end
        res[v] = cov
        println("$v: $(length(cov)) / $(length(fns)) functions covered")
    end
    (fns, res)
end

function baseline()
    fns, res = measure(); vs = sort(collect(keys(res)))
    allc = union(values(res)...)
    open(BASE_P, "w") do io
        println(io, "# function_id\tcovered_by_variants   (coverage_gate.jl baseline, $(Base.Libc.strftime("%Y-%m-%d", Base.time())))")
        for t in fns
            id = fid(t); id in allc || continue
            println(io, id, '\t', join([v for v in vs if id in res[v]], ","))
        end
    end
    gaps = [fid(t) for t in fns if !(fid(t) in allc)]
    open(GAPS_P, "w") do io
        println(io, "# function_id\treason — zero coverage by the tiered fixtures ($(join(vs, ",")) × regime matrix). Each is a")
        println(io, "# KNOWN coverage gap (not silently tolerated): other-variant code, unfixtured keyword paths, tools. Remove a line")
        println(io, "# when a fixture starts covering it; coverage_gate.jl check fails on a zero-coverage function not listed here.")
        for g in gaps; println(io, g, '\t', "not exercised by tiered fixtures ($(join(vs, ",")))"); end
    end
    date = Base.Libc.strftime("%Y-%m-%d", Base.time())
    md = normpath(joinpath(@__DIR__, "..", "..", "..", "docs", "TIERED_COVERAGE_$(date).md"))
    byfile = Dict{String,Vector{Int}}()
    for t in fns
        v = get!(byfile, dirname(t[1]) == "" ? t[1] : t[1], [0, 0]); v[1] += 1; fid(t) in allc && (v[2] += 1)
    end
    open(md, "w") do io
        println(io, "# Tiered-suite function coverage — $date\n")
        println(io, "Measured by `test/harness/tiered/coverage_gate.jl baseline` (fast-tier fixtures, full regime matrix, ",
                    "crash-isolated workers under `--code-coverage=user`). Function = top-level `function`/short-form definition ",
                    "(heuristic extent: to the next top-level definition).\n")
        println(io, "| variant | functions covered | of |\n|---|---|---|")
        for v in vs; println(io, "| $v | $(length(res[v])) | $(length(fns)) |"); end
        println(io, "| **any** | **$(length(allc))** | $(length(fns)) |\n")
        println(io, "Zero-coverage functions: $(length(gaps)) — listed in `test/fixtures/tiered/COVERAGE_GAPS.tsv`.\n")
        println(io, "## Per source file (functions covered / total), sorted by uncovered count\n")
        println(io, "| file | covered | total |\n|---|---|---|")
        for (f, c) in sort(collect(byfile); by = x -> -(x[2][1] - x[2][2]))
            println(io, "| `src/$f` | $(c[2]) | $(c[1]) |")
        end
    end
    println("baseline → $BASE_P ($(length(allc)) covered), gaps → $GAPS_P ($(length(gaps))), report → $md")
end

function check()
    base = Dict{String,Set{String}}()
    for l in eachline(BASE_P)
        startswith(l, "#") && continue
        id, vs = split(l, '\t'); base[id] = Set(split(vs, ','))
    end
    gaps = Set(first(split(l, '\t')) for l in eachline(GAPS_P) if !startswith(l, "#") && !isempty(l))
    fns, res = measure(); bad = 0
    for (id, vs) in base, v in vs
        haskey(res, v) || continue
        id in res[v] || (println("LOST COVERAGE: $id (variant $v)"); bad += 1)
    end
    allc = union(values(res)...)
    for t in fns
        id = fid(t)
        (id in allc || id in gaps || haskey(base, id)) && continue
        println("UNCOVERED, NOT A LISTED GAP: $id"); bad += 1
    end
    for id in gaps
        id in allc && println("newly covered (remove from COVERAGE_GAPS.tsv): $id")
    end
    println(bad == 0 ? "COVERAGE GATE: OK" : "COVERAGE GATE: $bad problem(s)")
    bad
end

if abspath(PROGRAM_FILE) == @__FILE__
    cmd = isempty(ARGS) ? "check" : ARGS[1]
    cmd == "baseline" ? baseline() : cmd == "check" ? exit(check() == 0 ? 0 : 1) :
        error("usage: coverage_gate.jl baseline | check")
end
