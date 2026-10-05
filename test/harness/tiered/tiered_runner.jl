# tiered_runner.jl — run jl on tiered fixtures and compare against the live goldens (fast tier), with the
# KNOWN_RESIDUALS allowlist and the per variant×regime signed tally. Shared by test/integration/test_tiered.jl,
# snapshot.jl and mutation.jl.
include(joinpath(@__DIR__, "tiered_common.jl"))
using TOML

struct Mismatch
    variant::String; stand::String; regime::String
    file::String      # "sum" | DBS table name | "TABLES" (presence) | "CRASH" | "FIXTURE" | "TALLY"
    col::String; year::String
    gold::String; got::String
end
key_str(m::Mismatch) = "$(m.variant)/$(m.stand)/$(m.regime)/$(m.file)/$(m.col)/$(m.year)"

# Variant groups for TIERED_VARIANTS. CORE = the variants under an active regime-close campaign — the DEFAULT set
# (Pkg.test runtime). WEST / EAST = the western / eastern coverage fixtures (measured residual maps, not yet dug): each
# adds ~1.5-4 min at TIERED_THREADS=3, so they run on request (TIERED_VARIANTS=WEST / EAST / ALL / CORE,TT,…), not in
# the default suite.
const TIERED_GROUPS = Dict(
    "CORE" => ["BM", "EM", "IE", "SN"],
    "WEST" => ["TT", "UT", "CI", "CR", "KT", "NC", "WC", "PN", "EC", "SO", "CA", "WS", "AK", "BC", "OC", "OP"],
    # EAST = the eastern coverage fixtures (CS, LS, NE + Ontario, the LS-derived metric variant): measured OPEN residual
    # maps, run on request like WEST (TIERED_VARIANTS=EAST / ALL), never in the default CORE set.
    "EAST" => ["CS", "LS", "NE", "ON"])

"""
Variants to run (sorted, only those with a fixture directory). ENV TIERED_VARIANTS: unset/empty ⇒ CORE; otherwise a
comma list of variant codes and/or group names (CORE, WEST, EAST, ALL = every fixture directory), e.g. `CORE,TT`.
"""
function tiered_variants()
    isdir(TIERED_ROOT) || return String[]
    vs = sort([uppercase(d) for d in readdir(TIERED_ROOT) if isdir(joinpath(TIERED_ROOT, d))])
    sel = uppercase(strip(get(ENV, "TIERED_VARIANTS", "")))
    toks = isempty(sel) ? ["CORE"] : [String(strip(t)) for t in split(sel, ',') if !isempty(strip(t))]
    "ALL" in toks && return vs
    want = Set{String}()
    for t in toks; haskey(TIERED_GROUPS, t) ? union!(want, TIERED_GROUPS[t]) : push!(want, t); end
    [v for v in vs if v in want]
end

fixture_dir(v) = joinpath(TIERED_ROOT, lowercase(v))
fixture_stands(v) = [strip(l) for l in eachline(joinpath(fixture_dir(v), "stands.txt")) if !isempty(strip(l))]

# QUICK tier (TIERED=quick): first 3 stands × these regimes — for refactor iteration (a few minutes).
const QUICK_REGIMES = ["none", "thinbba", "plant_cyc", "simfire", "mistletoe"]

"""
The regimes a variant's fixture covers: PROVENANCE.toml `regimes` (make_fixtures drops a regime whose extension the
variant's oracle build STUBS — FVS11 "requested extension is not part of this program", e.g. CA/AK have no WRD, AK no
Climate-FVS — and lists it under `skipped_regimes`), in REGIMES order. Falls back to REGIMES when absent.
"""
function fixture_regimes(v)
    p = joinpath(fixture_dir(v), "PROVENANCE.toml")
    isfile(p) || return REGIMES
    rs = get(TOML.parsefile(p), "regimes", nothing)
    rs === nothing ? REGIMES : [r for r in REGIMES if r in rs]
end

function tiered_cases(v; quick::Bool = get(ENV, "TIERED", "") == "quick")
    st = fixture_stands(v); rg = fixture_regimes(v)
    if quick; st = st[1:min(3, length(st))]; rg = [r for r in QUICK_REGIMES if r in rg]; end
    [(cn, r) for cn in st for r in rg]
end

"Run jl on one fixture case in a scratch dir. Returns (sum_text, out_db_path_or_empty, crashed, errmsg)."
function run_case(v::AbstractString, cn::AbstractString, r::AbstractString; dir::AbstractString = mktempdir(),
                  snapshot_hook = nothing)
    fx = fixture_dir(v); stem = "$(cn)_$(r)"
    kp = joinpath(fx, stem * ".key")
    isfile(kp) || return ("", "", true, "missing fixture key $kp")
    cp(joinpath(fx, "stands.db"), joinpath(dir, "stands.db"); force = true)
    outdb = joinpath(dir, "out.db"); isfile(outdb) && rm(outdb)
    # absolute DB paths (no cd: the process cwd is shared, so cd would race when stands run on parallel tasks)
    key = join([l == "out.db" ? outdb : l == "stands.db" ? joinpath(dir, "stands.db") : l
                for l in split(read(kp, String), '\n')], '\n')
    skey = joinpath(dir, "s.key"); write(skey, key)
    txt = ""; crashed = false; err = ""
    try
        txt = if snapshot_hook === nothing
            FVSjl.run_keyfile(skey; variant = VAR[uppercase(v)])
        else
            task_local_storage(:fvsjl_snapshot_hook, snapshot_hook) do
                FVSjl.run_keyfile(skey; variant = VAR[uppercase(v)])
            end
        end
    catch e
        crashed = true; err = sprint(showerror, e)
    end
    (txt, isfile(outdb) ? outdb : "", crashed, err)
end

_num(s) = tryparse(Float64, s)
_same(a::AbstractString, b::AbstractString) = a == b || (let x = _num(a), y = _num(b)
    x !== nothing && y !== nothing && x == y end)

# Row-matching key columns for DBS tables (rows compared within equal keys, in sorted order).
const KEY_COLS = ["StandID", "Year", "TreeId", "TreeIndex", "SpeciesFVS", "Species", "PtIndex", "Removal_Code",
                  "RmvCode", "Rank", "Fire_Type", "Size_Class", "SizeClass", "Stratum", "Condition"]

function compare_case(v, cn, r, jl_sum::AbstractString, jl_db::AbstractString)
    fx = fixture_dir(v); stem = "$(cn)_$(r)"; ms = Mismatch[]
    mk(file, col, yr, g, j) = push!(ms, Mismatch(v, cn, r, file, col, string(yr), string(g), string(j)))
    gsp = joinpath(fx, stem * ".live.sum")
    isfile(gsp) || (mk("FIXTURE", "live.sum", "*", "present", "missing"); return ms)
    gtxt = read(gsp, String)
    if startswith(gtxt, "# LIVE_NO_OUTPUT")
        isempty(sum_rows(jl_sum; metric = metric_sum(v))) || mk("sum", "LIVE_NO_OUTPUT", "*", "no output", "jl output")
        return ms
    end
    grows = sum_rows(gtxt; metric = metric_sum(v)); jrows = sum_rows(jl_sum; metric = metric_sum(v))
    ex = sum_extra_lines(jl_sum; metric = metric_sum(v))
    isempty(ex) || mk("sum", "EXTRA_CONTENT", "*", "none", "$(length(ex)) non-summary lines: $(strip(first(ex)))")
    jy = Dict(strip(r[1:4]) => r for r in jrows); gy = Dict(strip(r[1:4]) => r for r in grows)
    for g in grows
        y = strip(g[1:4]); j = get(jy, y, nothing)
        j === nothing && (mk("sum", "ROW", y, "present", "missing"); continue)
        for (nm, a, b) in sum_fields(v)
            gv = sum_field(g, a, b); jv = sum_field(j, a, b)
            gv == jv || mk("sum", nm, y, gv, jv)
        end
    end
    for j in jrows
        y = strip(j[1:4]); haskey(gy, y) || mk("sum", "ROW", y, "missing", "present")
    end
    # DBS table presence + contents
    gtp = joinpath(fx, stem * ".tables")
    isfile(gtp) || (mk("FIXTURE", "tables", "*", "present", "missing"); return ms)
    gt = [strip(l) for l in eachline(gtp) if !isempty(strip(l))]
    jt = isempty(jl_db) ? String[] : db_tables(jl_db)
    for t in setdiff(gt, jt); mk(t, "PRESENCE", "*", "present", "missing"); end
    for t in setdiff(jt, gt); mk(t, "PRESENCE", "*", "absent", "present"); end
    for t in intersect(gt, jt)
        gp = joinpath(fx, "$(stem).$(t).csv")
        isfile(gp) || (mk("FIXTURE", t * ".csv", "*", "present", "missing"); continue)
        ghdr, grws = read_csv(gp); jhdr, jrws = db_table_rows(jl_db, t)
        for c in setdiff(ghdr, jhdr); mk(t, c, "*", "column", "missing"); end
        for c in setdiff(jhdr, ghdr); mk(t, c, "*", "absent", "extra column"); end
        common = intersect(ghdr, jhdr)
        kc = [c for c in KEY_COLS if c in common]
        gi = Dict(c => findfirst(==(c), ghdr) for c in common); ji = Dict(c => findfirst(==(c), jhdr) for c in common)
        keyof(row, idx) = Tuple(row[idx[c]] for c in kc)
        gg = Dict{Any,Vector{Vector{String}}}(); jj = Dict{Any,Vector{Vector{String}}}()
        for rw in grws; push!(get!(gg, keyof(rw, gi), Vector{String}[]), rw); end
        for rw in jrws; push!(get!(jj, keyof(rw, ji), Vector{String}[]), rw); end
        for (k, gl) in gg
            jl = get(jj, k, Vector{String}[])
            yr = "Year" in kc ? k[findfirst(==("Year"), kc)] : "*"
            length(gl) == length(jl) || mk(t, "ROWCOUNT", yr, length(gl), length(jl))
            for i in 1:min(length(gl), length(jl)), c in common
                gv = gl[i][gi[c]]; jv = jl[i][ji[c]]
                _same(gv, jv) || mk(t, c, yr, gv, jv)
            end
        end
        for (k, jl) in jj
            haskey(gg, k) && continue
            yr = "Year" in kc ? k[findfirst(==("Year"), kc)] : "*"
            mk(t, "ROWCOUNT", yr, 0, length(jl))
        end
    end
    ms
end

# ---------------------------------------------------------------------------------------------------------------
# Signed tally per variant × regime (last common year; TPA / BA / TCuFt; material = >1 unit and ≥5%).
function case_signs(v, cn, r, jl_sum)
    fx = fixture_dir(v); gtxt = read(joinpath(fx, "$(cn)_$(r).live.sum"), String)
    startswith(gtxt, "# LIVE_NO_OUTPUT") && return nothing
    g = Dict(strip(x[1:4]) => x for x in sum_rows(gtxt; metric = metric_sum(v))); j = Dict(strip(x[1:4]) => x for x in sum_rows(jl_sum; metric = metric_sum(v)))
    ys = sort([y for y in keys(g) if haskey(j, y)]); isempty(ys) && return nothing
    y = ys[end]; out = Dict{String,Int}()
    for (nm, a, b) in (("TPA", 9, 14), ("BA", 15, 18), ("TCuFt", 37, 42))
        gv = something(tryparse(Float64, sum_field(g[y], a, b)), 0.0)
        jv = something(tryparse(Float64, sum_field(j[y], a, b)), 0.0)
        d = jv - gv
        mat = abs(d) > 1.0 + 1e-9 && (gv == 0 ? true : abs(d) / abs(gv) >= 0.05)
        out[nm] = mat ? (d > 0 ? 1 : -1) : 0
    end
    out
end

"TALLY mismatches for a variant×regime from per-case sign dicts: one-directional = ≥3 material, all one sign."
function tally_mismatches(v, r, signs::Vector)
    ms = Mismatch[]
    for nm in ("TPA", "BA", "TCuFt")
        s = [d[nm] for d in signs if d !== nothing]
        over = count(==(1), s); under = count(==(-1), s)
        if (over >= 3 && under == 0) || (under >= 3 && over == 0)
            push!(ms, Mismatch(v, "*", r, "TALLY", nm, "*", "balanced", "over=$over under=$under"))
        end
    end
    ms
end

# ---------------------------------------------------------------------------------------------------------------
# Allowlist (test/fixtures/tiered/KNOWN_RESIDUALS.toml)
const ALLOWLIST_PATH = joinpath(TIERED_ROOT, "KNOWN_RESIDUALS.toml")

struct Residual
    variant::String; stand::String; regime::String; file::String; col::String; year::String
    status::String; text::String; idx::Int
    max_cells::Int   # measured mismatch-cell ceiling: MORE cells than this FAIL (a grouped entry cannot hide spread)
end
function load_allowlist(path = ALLOWLIST_PATH)
    isfile(path) || return (Residual[], String["allowlist file missing: $path"])
    d = TOML.parsefile(path); out = Residual[]; errs = String[]
    for (i, e) in enumerate(get(d, "residual", Any[]))
        g(k) = string(get(e, k, "*"))
        st = string(get(e, "status", ""))
        txt = st == "CORNER" ? string(get(e, "proof", "")) : st == "OPEN" ? string(get(e, "issue", "")) : ""
        if !(st in ("CORNER", "OPEN"))
            push!(errs, "entry #$i: status must be CORNER or OPEN (got \"$st\")")
        elseif isempty(strip(txt))
            push!(errs, "entry #$i: $(st) requires $(st == "CORNER" ? "proof" : "issue") text")
        end
        push!(out, Residual(uppercase(g("variant")), g("stand"), g("regime"), g("file"), g("col"), g("year"), st, txt, i,
                            Int(get(e, "max_cells", typemax(Int)))))
    end
    (out, errs)
end
_m(pat, val) = pat == "*" || pat == val
covers(e::Residual, m::Mismatch) = _m(e.variant, m.variant) && _m(e.stand, m.stand) && _m(e.regime, m.regime) &&
                                    _m(e.file, m.file) && _m(e.col, m.col) && _m(e.year, m.year)
"An allowlist entry applies to a run only if its variant/regime (and stand) are part of the run."
in_scope(e::Residual, ran::Set{Tuple{String,String,String}}) =
    any(t -> _m(e.variant, t[1]) && _m(e.stand, t[2]) && _m(e.regime, t[3]), ran)

"""
Classify mismatches against the allowlist. Returns (unlisted, stale_entries, used_counts, over_ceiling) where stale
entries are in-scope entries that matched nothing (an unexpected pass ⇒ must be removed) and over_ceiling are entries
whose matched cell count exceeds their measured `max_cells` (the residual spread ⇒ fail).
"""
function apply_allowlist(ms::Vector{Mismatch}, al::Vector{Residual}, ran)
    unlisted = Mismatch[]; used = zeros(Int, length(al))
    for m in ms
        k = findfirst(e -> covers(e, m), al)
        k === nothing ? push!(unlisted, m) : (used[k] += 1)
    end
    stale = [e for (i, e) in enumerate(al) if used[i] == 0 && in_scope(e, ran)]
    over = [(e, used[i]) for (i, e) in enumerate(al) if used[i] > e.max_cells]
    (unlisted, stale, used, over)
end

# ---------------------------------------------------------------------------------------------------------------
# Crash isolation: run cases in tiered_worker.jl subprocesses; a worker that dies (SIGSEGV) leaves an in-flight case
# (START without DONE) — record it as a crash and restart a fresh worker on the remaining cases.
"""
    run_isolated(v, cases; mode = "compare", threads = 1) -> (lines::Vector{String}, crashed::Vector{Tuple})

Runs `cases` (Vector of (cn, regime)) for variant `v` in worker subprocess(es). Returns the worker result lines and
the cases that crashed the worker (each retried once serially before being declared a crash).
"""
function run_isolated(v, cases; mode::AbstractString = "compare", threads::Int = 1)
    root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    worker = joinpath(@__DIR__, "tiered_worker.jl")
    remaining = collect(cases); lines = String[]; crashed = Tuple{String,String}[]
    retried = Set{Tuple{String,String}}()
    nth = threads
    while !isempty(remaining)
        d = mktempdir(); cf = joinpath(d, "cases.tsv"); of = joinpath(d, "out.tsv"); pf = joinpath(d, "prog.log")
        open(cf, "w") do io; foreach(c -> println(io, c[1], '\t', c[2]), remaining); end
        touch(of); touch(pf)
        cov = get(ENV, "TIERED_COVERAGE", "") == "1" ? `--code-coverage=user` : ``   # coverage_gate.jl
        cmd = `$(Base.julia_cmd()) --project=$root --threads=$nth $cov $worker $mode $v $cf $of $pf`
        p = run(pipeline(ignorestatus(cmd); stdout = devnull, stderr = joinpath(d, "err.log")))
        append!(lines, readlines(of))
        started = Set{Tuple{String,String}}(); done = Set{Tuple{String,String}}()
        for l in readlines(pf)
            t = split(l); length(t) == 3 || continue
            (t[1] == "START" ? started : done) |> s -> push!(s, (String(t[2]), String(t[3])))
        end
        inflight = [c for c in remaining if c in started && !(c in done)]
        success(p) && isempty(inflight) && break
        # worker died: in-flight cases retry once serially (a threaded crash cannot name the culprit), then crash
        for c in inflight
            if c in retried || nth == 1
                push!(crashed, c)
                # discard partial output lines of a crashed case
                filter!(l -> !(occursin("\t$(c[1])\t$(c[2])", l) || occursin("\t$(c[1])\t$(c[2])\t", l)), lines)
            else
                push!(retried, c)
            end
        end
        todo = [c for c in remaining if !(c in done) && !(c in crashed)]
        if !isempty(inflight) && nth > 1
            nth = 1                                   # isolate the retry serially
        end
        if isempty(inflight) && !success(p)
            # died before starting anything: surface the worker error instead of looping
            error("tiered worker failed without running a case: " * read(joinpath(d, "err.log"), String)[1:min(end, 2000)])
        end
        remaining = todo
    end
    (lines, crashed)
end

"Parse compare-mode worker lines (+ crashed cases) into (mismatches, signs::Dict{regime=>Vector})."
function collect_compare(v, lines, crashed)
    ms = Mismatch[]; signs = Dict{String,Vector{Any}}()
    for l in lines
        f = split(l, '\t')
        if f[1] == "M" && length(f) >= 8
            push!(ms, Mismatch(v, f[7], f[8], f[2], f[3], f[4], f[5], f[6]))
        elseif f[1] == "S" && length(f) >= 6
            push!(get!(signs, f[3], Any[]), Dict("TPA" => parse(Int, f[4]), "BA" => parse(Int, f[5]), "TCuFt" => parse(Int, f[6])))
        end
    end
    for (cn, r) in crashed
        push!(ms, Mismatch(v, cn, r, "CRASH", "*", "*", "", "process crash (SIGSEGV/abort) — worker died"))
    end
    (ms, signs)
end

"Run all fast-tier cases for the selected variants crash-isolated; returns (mismatches incl tallies, ran, per-case exact map)."
function run_fast_tier(; quick::Bool = get(ENV, "TIERED", "") == "quick", threads::Int = 1)
    all_ms = Mismatch[]; ran = Set{Tuple{String,String,String}}(); percase = Dict{Tuple{String,String,String},Bool}()
    for v in tiered_variants()
        cases = tiered_cases(v; quick)
        lines, crashed = run_isolated(v, cases; mode = "compare", threads)
        ms, signs = collect_compare(v, lines, crashed)
        for (cn, r) in cases; push!(ran, (v, cn, r)); percase[(v, cn, r)] = true; end
        for m in ms; haskey(percase, (v, m.stand, m.regime)) && (percase[(v, m.stand, m.regime)] = false); end
        append!(all_ms, ms)
        for (r, sg) in signs; append!(all_ms, tally_mismatches(v, r, sg)); end
    end
    (all_ms, ran, percase)
end
