# snapshot.jl — BIT-IDENTITY SNAPSHOT TIER (jl vs jl) of the tiered integration suite.
#
# Oracle comparison cannot see a refactor that shifts an allowlisted OPEN/CORNER case, nor sub-print-precision ULP
# drift. This tier records jl's OWN complete outputs at a blessed commit and requires bit-identity:
#   * the .sum data rows (text; the -999 header carries wall-clock date/time and is excluded),
#   * every DBS output table at full stored precision (CaseID dropped),
#   * the full per-record state at every summary cycle: every numeric/Bool vector field of the TreeList and of the
#     Scratch work arrays (WK/IND/…), sliced to the live record count, hashed as raw bits (Float32 via reinterpret).
# Stored compactly: one SHA-1 per (stand×regime, kind) in test/fixtures/tiered/<v>/SNAPSHOT.tsv (tracked). Full dumps
# are regenerated on demand by `detail`, which checks out the blessed commit in a temp worktree to diff per record.
#
# Usage (julia --project=. test/harness/tiered/snapshot.jl …):
#   bless  "<reason>"          record snapshots for every fixture case at HEAD (commit + reason + date in the manifest)
#   check                      recompute and compare; prints the FIRST differing stand×regime / kind; exit 1 on diff
#   detail <V> <cn>_<regime>   first differing cycle / variable / record vs the blessed commit (temp worktree)
#   dump   <V> <cn>_<regime> <out.tsv>   write the full per-record dump (used by `detail`)
# Env: TIERED_VARIANTS=BM,SN (restrict), TIERED=quick (quick subset), TIERED_THREADS=N (cases run on N parallel
# tasks in the crash-isolated worker) — run `check` at TIERED_THREADS=1 and =N to prove thread determinism (both must
# match the manifest bit-for-bit).
include(joinpath(@__DIR__, "tiered_runner.jl"))
using SHA

_vecfields(x) = [f for f in fieldnames(typeof(x)) if getfield(x, f) isa AbstractVector{<:Union{Real,Bool}}]

function _bits(v::AbstractVector, n::Int)
    m = min(n, length(v)); w = view(v, 1:m)
    eltype(v) === Float32 ? collect(reinterpret(UInt32, collect(w))) :
    eltype(v) === Float64 ? collect(reinterpret(UInt64, collect(w))) :
    eltype(v) === Bool    ? UInt8.(collect(w)) : collect(w)
end

"Per-cycle state records: Vector of (cycle, year, field, bits-vector)."
function _collector()
    recs = Tuple{Int,Int,String,Any}[]
    hook = (s, year, c) -> begin
        n = s.trees.n
        push!(recs, (c, year, "n", [n]))
        for f in _vecfields(s.trees);  push!(recs, (c, year, "trees." * String(f),   _bits(getfield(s.trees, f), n))); end
        for f in _vecfields(s.scratch); push!(recs, (c, year, "scratch." * String(f), _bits(getfield(s.scratch, f), n))); end
    end
    (recs, hook)
end

_h(x) = bytes2hex(sha1(x))
function _hbits(recs)
    io = IOBuffer()
    for r in recs
        write(io, r[3], ':'); v = r[4]
        write(io, Int64(length(v)))
        for x in v; write(io, x); end        # raw bits (UInt32 for Float32, Int32/Int64/UInt8 as stored)
        write(io, ';')
    end
    _h(take!(io))
end

function snapshot_case(v, cn, r)
    recs, hook = _collector()
    txt, db, crashed, err = run_case(v, cn, r; snapshot_hook = hook)
    out = Pair{String,String}[]
    crashed && return [("crash" => _h(err))]
    push!(out, "sum" => _h(join(sum_rows(txt; metric = metric_sum(v)), '\n')))
    if !isempty(db)
        for t in db_tables(db)
            hdr, rows = db_table_rows(db, t)
            push!(out, "table:$t" => _h(join(hdr, ',') * "\n" * join((join(rw, ',') for rw in rows), '\n')))
        end
    end
    for c in sort(unique(first.(recs)))
        push!(out, "cycle:$c" => _hbits([x for x in recs if x[1] == c]))
    end
    out
end

manifest_path(v) = joinpath(fixture_dir(v), "SNAPSHOT.tsv")

"Snapshots for every case of the variants, crash-isolated (tiered_worker.jl); threads from ENV TIERED_THREADS."
function compute_all(vs; threads::Int = parse(Int, get(ENV, "TIERED_THREADS", "1")))
    out = Dict{Tuple{String,String},Vector{Pair{String,String}}}()
    for v in vs
        lines, crashed = run_isolated(v, tiered_cases(v); mode = "snapshot", threads)
        for l in lines
            f = split(l, '\t'); (f[1] == "H" && length(f) >= 5) || continue
            push!(get!(out, (v, "$(f[2])_$(f[3])"), Pair{String,String}[]), String(f[4]) => String(f[5]))
        end
        for (cn, r) in crashed; out[(v, "$(cn)_$(r)")] = ["crash" => "process-crash"]; end
    end
    out
end

function bless(reason)
    isempty(strip(reason)) && error("bless requires a reason")
    vs = tiered_variants(); snaps = compute_all(vs)
    commit = _git(normpath(joinpath(@__DIR__, "..", "..", "..")))
    for v in vs
        open(manifest_path(v), "w") do io
            println(io, "# blessed_commit\t$commit")
            println(io, "# blessed_date\t$(Base.Libc.strftime("%Y-%m-%d %H:%M:%S", Base.time()))")
            println(io, "# reason\t$(replace(reason, '\n' => ' '))")
            println(io, "# threads\t$(get(ENV, "TIERED_THREADS", "1"))")
            for k in sort([k for k in keys(snaps) if k[1] == v]; by = last), (kind, h) in snaps[k]
                println(io, k[2], '\t', kind, '\t', h)
            end
        end
        println("blessed $(manifest_path(v))")
    end
end
_git(dir) = try strip(read(`git -C $dir rev-parse HEAD`, String)) catch; "unknown" end
_git_dirty(dir) = try !isempty(strip(read(`git -C $dir status --porcelain --untracked-files=no`, String))) catch; true end

function read_manifest(v)
    p = manifest_path(v); isfile(p) || error("no snapshot manifest $p — run `snapshot.jl bless \"reason\"`")
    meta = Dict{String,String}(); rows = Dict{String,Vector{Pair{String,String}}}()
    for l in eachline(p)
        if startswith(l, "# ")
            k, val = split(l[3:end], '\t'; limit = 2); meta[k] = val
        elseif !isempty(l)
            stem, kind, h = split(l, '\t'); push!(get!(rows, stem, Pair{String,String}[]), kind => h)
        end
    end
    (meta, rows)
end

"Compare a stem's current snapshot list to the manifest; returns the first differing kind or nothing."
function first_diff(cur, man)
    dc = Dict(cur); dm = Dict(man)
    order = sort(unique(vcat(first.(cur), first.(man))); by = k -> (k == "sum" ? 0 : startswith(k, "cycle:") ?
                 1 + parse(Int, k[7:end]) : 10_000, k))
    for k in order
        get(dc, k, "absent") == get(dm, k, "absent") || return (k, get(dm, k, "absent"), get(dc, k, "absent"))
    end
    nothing
end

function check()
    vs = tiered_variants(); snaps = compute_all(vs); ndiff = 0
    for v in vs
        meta, man = read_manifest(v)
        for (k, cur) in sort(collect(snaps); by = x -> x[1])
            k[1] == v || continue
            if !haskey(man, k[2])
                get(ENV, "TIERED", "") == "quick" || (println("NEW CASE (not in manifest): $v $(k[2])"); ndiff += 1)
                continue
            end
            d = first_diff(cur, man[k[2]])
            d === nothing && continue
            ndiff += 1
            ndiff <= 20 && println("DIFF $v $(k[2]) first differing: $(d[1])  (blessed $(first(d[2], 12)) now $(first(d[3], 12)))",
                                   "  → julia --project=. test/harness/tiered/snapshot.jl detail $v $(k[2])")
        end
        println("$v: blessed at $(get(meta, "blessed_commit", "?")) — $(get(meta, "reason", ""))")
    end
    println(ndiff == 0 ? "SNAPSHOT CHECK: bit-identical ($(length(snaps)) cases, threads=$(get(ENV, "TIERED_THREADS", "1")))" :
                         "SNAPSHOT CHECK: $ndiff differing case(s)")
    ndiff
end

function dump_case(v, stem, out)
    cn, r = let i = findfirst('_', stem); (stem[1:i-1], stem[i+1:end]) end
    recs, hook = _collector(); run_case(v, cn, r; snapshot_hook = hook)
    open(out, "w") do io
        for (c, y, f, b) in recs; println(io, c, '\t', y, '\t', f, '\t', join(string.(b), ',')); end
    end
end

"First differing cycle/variable/record between the blessed commit and HEAD for one case."
function detail(v, stem)
    meta, _ = read_manifest(v); commit = get(meta, "blessed_commit", "")
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    wt = mktempdir(); rm(wt)
    run(`git -C $root worktree add --detach $wt $commit`)
    try
        for f in ("deps/libfvsmath.so",)   # gitignored runtime artifact
            isfile(joinpath(root, f)) && cp(joinpath(root, f), joinpath(wt, f); force = true)
        end
        a = joinpath(mktempdir(), "blessed.tsv"); b = joinpath(mktempdir(), "now.tsv")
        run(`$(Base.julia_cmd()) --project=$wt $(joinpath(wt, "test/harness/tiered/snapshot.jl")) dump $v $stem $a`)
        dump_case(v, stem, b)
        la = readlines(a); lb = readlines(b)
        for i in 1:max(length(la), length(lb))
            x = i <= length(la) ? split(la[i], '\t') : ["-","-","(absent)",""]
            y = i <= length(lb) ? split(lb[i], '\t') : ["-","-","(absent)",""]
            x == y && continue
            println("FIRST DIFF: cycle $(x[1]) year $(x[2]) variable $(x[3]) (now: $(y[3]))")
            va = split(x[4], ','); vb = split(y[4], ',')
            j = findfirst(k -> k > length(va) || k > length(vb) || va[k] != vb[k], 1:max(length(va), length(vb)))
            if j !== nothing
                fa = j <= length(va) ? va[j] : "(absent)"; fb = j <= length(vb) ? vb[j] : "(absent)"
                dec(s) = (u = tryparse(UInt32, s); u === nothing ? s : "$(s) = $(reinterpret(Float32, u))")
                println("  record $j: blessed $(dec(fa))  now $(dec(fb))")
            end
            return
        end
        println("no difference in the per-record dump (difference is in .sum/DBS output only)")
    finally
        run(`git -C $root worktree remove --force $wt`)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    cmd = isempty(ARGS) ? "check" : ARGS[1]
    if cmd == "bless"
        _git_dirty(normpath(joinpath(@__DIR__, "..", "..", ".."))) &&
            @warn "working tree is dirty — the manifest records HEAD, commit first so the blessed state is reproducible"
        bless(length(ARGS) >= 2 ? ARGS[2] : "")
    elseif cmd == "check"
        exit(check() == 0 ? 0 : 1)
    elseif cmd == "detail"
        detail(uppercase(ARGS[2]), ARGS[3])
    elseif cmd == "dump"
        dump_case(uppercase(ARGS[2]), ARGS[3], ARGS[4])
    else
        error("usage: snapshot.jl bless \"reason\" | check | detail <V> <stem> | dump <V> <stem> <out>")
    end
end
