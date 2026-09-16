"""
    FiaSim

The FVSjl simulation adapter — the one stable seam `simulate(plt_cns, plan, horizon)`.

A management **plan** (grow / thin / plant actions) is translated to FVS keywords and
injected into a per-plot keyfile that reads the plot from an indexed FIA sub-DB; FVSjl
runs it and we parse per-cycle **standing** (BA/TPA/QMD/volume) + **removed** volume from
the `.sum`, and per-cycle carbon from the FFE Stand Carbon Report. Results cache by
`(plt_cn, plan-hash, cycles, period)`.

Scale: one plot = one run against an indexed sub-DB (built once per AOI plot-set — the
71 GB FIA master is unindexed on the plot tables, so we extract the AOI's plots once).
The master path honors the FIA_DB env var, else the workspace default.
"""
module FiaSim

using FVSjl
import SQLite, DBInterface, SHA

export ManagementPlan, PlanAction, simulate_plots, build_subdb, CycleMetrics, PlotResult,
       LandscapePolicy, simulate_landscape

# Resolved at RUNTIME (honors FIA_DB set when the server starts). The default is derived
# from this module's own location — the workspace root that holds the repo — not hardcoded:
#   <root>/apps/treemap-explorer/src/FiaSim.jl  ->  up 4 = <workspace>  ->  SQLite_FIADB_ENTIRE.db
master_db() = get(ENV, "FIA_DB",
    normpath(joinpath(@__DIR__, "..", "..", "..", "..", "SQLite_FIADB_ENTIRE.db")))

# ---------------------------------------------------------------------------
# Plan model (mirrors the frontend JSON)
# ---------------------------------------------------------------------------
"One scheduled management action. `kind` ∈ (\"grow\",\"thin\",\"plant\").
Timing is RELATIVE: `cycle` is a 1-based projection cycle (fires ~cycle×period years
out), scheduled via FVS cycle-number dating so it lands at the same relative point on
every plot regardless of that plot's inventory year."
struct PlanAction
    kind::String
    cycle::Int                   # 1-based projection cycle the action fires in
    # thin:
    metric::String               # "BA" | "TPA" | "SDI"  (residual target metric)
    target::Float64              # residual value (0 = clearcut)
    direction::String            # "below" | "above"
    dbh_lo::Float64
    dbh_hi::Float64
    # plant:
    species::String              # 2-letter FVS alpha code (e.g. "DF","PP")
    tpa::Float64                 # trees/acre planted
    survival::Float64            # percent 0-100
end
PlanAction(; kind, cycle=1, metric="BA", target=0.0, direction="below",
             dbh_lo=0.0, dbh_hi=999.0, species="", tpa=0.0, survival=100.0) =
    PlanAction(kind, cycle, metric, target, direction, dbh_lo, dbh_hi, species, tpa, survival)

struct ManagementPlan
    name::String
    cycles::Int
    period::Int                  # cycle length (years)
    actions::Vector{PlanAction}
end
ManagementPlan(; name="grow-only", cycles=10, period=10, actions=PlanAction[]) =
    ManagementPlan(name, cycles, period, actions)

"Stable hash of a plan (for caching)."
plan_hash(p::ManagementPlan) = bytes2hex(SHA.sha1(string(
    p.cycles, "|", p.period, "|",
    join(("$(a.kind):$(a.cycle):$(a.metric):$(a.target):$(a.direction):$(a.dbh_lo):$(a.dbh_hi):$(a.species):$(a.tpa):$(a.survival)"
          for a in p.actions), ";"))))[1:16]

# ---------------------------------------------------------------------------
# Plan -> FVS keyword cards (column-aligned: token cols 1-10, fields 10 cols each)
# ---------------------------------------------------------------------------
_kwrec(kw, fields...) = rpad(kw, 10) * join(lpad(_fmt(x), 10) for x in fields)
_fmt(x::Integer) = string(x)
_fmt(x::AbstractFloat) = (x == round(x) ? string(Int(round(x))) : string(round(x; digits=4)))
_fmt(x) = string(x)

const _THIN_TOKEN = Dict(("BA","below")=>"THINBBA", ("BA","above")=>"THINABA",
                         ("TPA","below")=>"THINBTA", ("TPA","above")=>"THINATA",
                         ("SDI","below")=>"THINSDI", ("SDI","above")=>"THINSDI")

"Emit the management keyword block for a plan (goes after DATABASE…END, before PROCESS)."
function plan_keywords(p::ManagementPlan)::String
    io = IOBuffer()
    thins = filter(a -> a.kind == "thin", p.actions)
    plants = filter(a -> a.kind == "plant", p.actions)
    for a in thins
        allsp = isempty(a.species) || lowercase(a.species) == "all"
        if allsp
            # all species → THINBBA/ABA/BTA/… : date(cycle+1), residual target, cut-eff, DBH lo/hi
            tok = get(_THIN_TOKEN, (a.metric, a.direction), "THINBBA")
            println(io, _kwrec(tok, a.cycle + 1, a.target, 1, a.dbh_lo, a.dbh_hi))
        else
            # species-specific → one THINDBH per listed species (comma-separated list allowed):
            # date, DBH lo/hi, cut-eff, SPECIES(alpha), resid CTPA, resid CBA
            ctpa = a.metric == "TPA" ? _fmt(a.target) : ""
            cba  = a.metric == "TPA" ? "" : _fmt(a.target)     # BA (and SDI→BA fallback)
            for sp in split(a.species, ',')
                s = strip(sp)
                isempty(s) && continue
                println(io, _kwrec("THINDBH", a.cycle + 1, a.dbh_lo, a.dbh_hi, 1, s, ctpa, cba))
            end
        end
    end
    if !isempty(plants)
        println(io, "ESTAB")
        for a in plants
            # PLANT: DATE(cycle a.cycle+1), species(alpha), trees/acre, survival%.
            # A comma-separated species list plants each species at the given tpa.
            for sp in split(a.species, ',')
                s = strip(sp)
                isempty(s) && continue
                println(io, _kwrec("PLANT", a.cycle + 1, s, a.tpa, a.survival))
            end
        end
        println(io, "END")
    end
    String(take!(io))
end

# ---------------------------------------------------------------------------
# Sub-DB: extract the AOI's plots from the master into a small indexed DB (cached)
# ---------------------------------------------------------------------------
"Path of the cached sub-DB for a set of PLT_CNs (built if missing). Returns (path, variant_of)."
function build_subdb(plt_cns::Vector{String}; cache_dir::AbstractString)
    mkpath(cache_dir)
    key = bytes2hex(SHA.sha1(join(sort(plt_cns), ",")))[1:16]
    path = joinpath(cache_dir, "aoi_$key.db")
    varmap = Dict{String,String}()
    # reuse a cached sub-DB, but if it's unreadable (a half-written/corrupt build from an
    # interrupted run) drop it and rebuild rather than surfacing "unable to open".
    if isfile(path)
        try
            db = SQLite.DB(path)
            for r in DBInterface.execute(db, "SELECT STAND_CN, VARIANT FROM FVS_STANDINIT_PLOT")
                varmap[string(r.STAND_CN)] = String(strip(String(r.VARIANT)))
            end
            SQLite.close(db)
            return path, varmap
        catch
            @warn "cached sub-DB unreadable, rebuilding" path
            rm(path; force = true); empty!(varmap)
        end
    end

    mdb = master_db()
    isfile(mdb) || error("FIA database not found at \"$mdb\". Set the FIA_DB environment " *
        "variable to your SQLite_FIADB_ENTIRE.db before starting the server.")

    idlist = join(("'" * cn * "'" for cn in plt_cns), ",")
    tmp = path * ".building"                       # atomic: build to a temp, rename on success
    rm(tmp; force = true)
    src = SQLite.DB("file:$mdb?mode=ro&immutable=1")
    dst = SQLite.DB(tmp)
    try
        for tbl in ("FVS_STANDINIT_PLOT", "FVS_TREEINIT_PLOT")
            ddl = first(DBInterface.execute(src,
                "SELECT sql FROM sqlite_master WHERE name='$tbl'")).sql
            DBInterface.execute(dst, ddl)
            DBInterface.execute(dst, "ATTACH DATABASE 'file:$mdb?mode=ro&immutable=1' AS m")
            DBInterface.execute(dst, "INSERT INTO $tbl SELECT * FROM m.$tbl WHERE STAND_CN IN ($idlist)")
            DBInterface.execute(dst, "DETACH DATABASE m")
            DBInterface.execute(dst, "CREATE INDEX idx_$(tbl)_cn ON $tbl(STAND_CN)")
        end
        for r in DBInterface.execute(dst, "SELECT STAND_CN, VARIANT FROM FVS_STANDINIT_PLOT")
            varmap[string(r.STAND_CN)] = String(strip(String(r.VARIANT)))
        end
    catch e
        SQLite.close(src); SQLite.close(dst); rm(tmp; force = true)
        rethrow(e)
    end
    SQLite.close(src); SQLite.close(dst)
    mv(tmp, path; force = true)
    return path, varmap
end

# ---------------------------------------------------------------------------
# Run one plot + parse per-cycle metrics
# ---------------------------------------------------------------------------
"Per-cycle metrics for one plot (per-acre volumes; carbon metric t/ha)."
struct CycleMetrics
    year::Int
    ba::Float64          # ft2/acre standing
    tpa::Float64
    qmd::Float64
    tcuft::Float64       # ft3/acre standing total cubic
    bdft::Float64        # bd ft/acre standing
    rem_tcuft::Float64   # ft3/acre removed this cycle
    rem_bdft::Float64
    carbon_total::Float64      # t C/ha total stand
    carbon_live_ag::Float64    # t C/ha aboveground live
    carbon_removed::Float64    # t C/ha removed
end

# fixed-column ranges in the SUMOUT row (1-indexed char spans)
_col(s, a, b) = (t = strip(s[a:min(b, lastindex(s))]); isempty(t) ? 0.0 : parse(Float64, t))

function _parse_sum_and_carbon(out::AbstractString)::Vector{CycleMetrics}
    sumrows = Dict{Int,NTuple{7,Float64}}()   # year => (ba,tpa,qmd,tcuft,bdft,rem_tcuft,rem_bdft)
    carbon = Dict{Int,NTuple{3,Float64}}()    # year => (total, live_ag, removed)
    incarbon = false
    for ln in split(out, '\n')
        if occursin("STAND CARBON REPORT", ln); incarbon = true; continue; end
        s = String(ln)
        # .sum period row: starts with a 4-digit year at col 1, long fixed-width line
        if !incarbon && length(s) >= 60 && occursin(r"^\d{4}", s)
            y = Int(_col(s, 1, 4))
            sumrows[y] = (_col(s,15,18), _col(s,9,14), _col(s,32,36),
                          _col(s,37,42), _col(s,55,60), _col(s,67,72), _col(s,85,90))
        elseif incarbon && occursin(r"^\s*\d{4}\s", s)
            f = split(strip(s))
            length(f) >= 11 || continue
            y = parse(Int, f[1])
            carbon[y] = (parse(Float64,f[10]), parse(Float64,f[2]), parse(Float64,f[11]))
        end
    end
    out2 = CycleMetrics[]
    for y in sort(collect(keys(sumrows)))
        ba,tpa,qmd,tcuft,bdft,rtc,rbd = sumrows[y]
        ct,cl,cr = get(carbon, y, (NaN,NaN,NaN))
        push!(out2, CycleMetrics(y, ba,tpa,qmd,tcuft,bdft,rtc,rbd, ct,cl,cr))
    end
    out2
end

dbh_class_lo(dia) = min(2 * floor(Int, dia / 2), 40)   # 2-inch DBH class, 40+ cap

"Per-plot result: per-cycle metrics + per-cycle (by year) species×DBH tally + species symbols."
struct PlotResult
    metrics::Vector{CycleMetrics}
    dist::Dict{Int,Dict{Tuple{Int16,Int},NTuple{2,Float64}}}   # year => (spcd,dbh_lo) => (tpa/ac, ba/ac)
    spsym::Dict{Int16,String}
end

"Read the per-cycle FVS_TreeList (start-of-cycle) from a DBS output DB into a per-year species×DBH tally."
function _read_treelist(outdb::AbstractString)
    dist = Dict{Int,Dict{Tuple{Int16,Int},NTuple{2,Float64}}}()
    spsym = Dict{Int16,String}()
    isfile(outdb) || return dist, spsym
    db = SQLite.DB("file:$outdb?mode=ro")
    try
        isempty(collect(DBInterface.execute(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND name='FVS_TreeList'"))) &&
            return dist, spsym
        for r in DBInterface.execute(db, "SELECT Year, SpeciesFIA, SpeciesPLANTS, DBH, TPA FROM FVS_TreeList")
            (r.TPA === missing || r.DBH === missing || r.Year === missing) && continue
            tpa = Float64(r.TPA); tpa > 0 || continue
            d = Float64(r.DBH);   d   > 0 || continue
            spcd = r.SpeciesFIA === missing ? Int16(0) : Int16(parse(Int, string(r.SpeciesFIA)))
            get!(spsym, spcd, r.SpeciesPLANTS === missing ? "" : String(r.SpeciesPLANTS))
            g = get!(dist, Int(r.Year), Dict{Tuple{Int16,Int},NTuple{2,Float64}}())
            k = (spcd, dbh_class_lo(d))
            c, b = get(g, k, (0.0, 0.0))
            g[k] = (c + tpa, b + 0.005454154 * d * d * tpa)
        end
    finally
        SQLite.close(db)
    end
    dist, spsym
end

const _MEM = Dict{String,PlotResult}()   # cache: "plt_cn|planhash" => PlotResult
const _MEM_LOCK = ReentrantLock()

# live progress of the current simulate_plots call (polled by /api/simprogress)
const _PROG = (done = Threads.Atomic{Int}(0), total = Threads.Atomic{Int}(0),
               running = Threads.Atomic{Bool}(false))
progress() = (; done = _PROG.done[], total = _PROG.total[], running = _PROG.running[])
_tick!() = Threads.atomic_add!(_PROG.done, 1)

"Run one plot; returns its `PlotResult` (metrics + per-cycle species×DBH) or `nothing`."
function _run_one(cn::String, vcode::String, subdb::String, plan::ManagementPlan, mgmt::String)
    key = tempname() * ".key"
    outdb = tempname() * ".db"            # per-run DBS output for the per-cycle FVS_TreeList
    write(key, string(
        "STDIDENT\n", cn, "\n",
        "DATABASE\nDSNin\n", subdb, "\n",
        "StandSQL\nSELECT * FROM FVS_STANDINIT_PLOT WHERE STAND_CN = '%StandID%'\nEndSQL\n",
        "TreeSQL\nSELECT * FROM FVS_TREEINIT_PLOT WHERE STAND_CN = '%StandID%'\nEndSQL\n",
        "END\n",
        "DATABASE\nDSNOut\n", outdb, "\nTREELIDB\nEND\n",
        _kwrec("NUMCYCLE", plan.cycles), "\n",
        mgmt,
        "FMIn\nCARBREPT\nEnd\n",
        "ECHOSUM\nPROCESS\nSTOP\n"))
    try
        out = FVSjl.run_keyfile(key; variant=FVSjl.variant_from_code(vcode), period=plan.period)
        metrics = _parse_sum_and_carbon(out)
        dist, spsym = _read_treelist(outdb)
        return PlotResult(metrics, dist, spsym)
    catch e
        @warn "simulate failed" cn vcode exception=e
        return nothing
    finally
        rm(outdb; force=true)
        rm(key; force=true)
    end
end

"""
    simulate_plots(plt_cns, plan; cache_dir) -> Dict{plt_cn => Vector{CycleMetrics}}

Run every plot under the plan (deduped; cached). Runs are parallelized across Julia
threads (start with `--threads=auto`); each variant is warmed up serially first so a
plot's lazy per-variant caches are populated before the threaded fan-out.
"""
function simulate_plots(plt_cns::Vector{String}, plan::ManagementPlan; cache_dir::AbstractString)
    plt_cns = unique(plt_cns)
    _PROG.total[] = length(plt_cns); _PROG.done[] = 0; _PROG.running[] = true
    try
        subdb, varmap = build_subdb(plt_cns; cache_dir=cache_dir)
        ph = plan_hash(plan)
        mgmt = plan_keywords(plan)
        results = Dict{String,PlotResult}()
        torun = String[]
        for cn in plt_cns
            r = lock(_MEM_LOCK) do; get(_MEM, "$cn|$ph", nothing); end
            if r !== nothing
                results[cn] = r; _tick!()                 # cached = instant
            elseif isempty(get(varmap, cn, ""))
                _tick!()                                  # absent from FIA plot table
            else
                push!(torun, cn)
            end
        end
        isempty(torun) && return results

        # warm up lazy per-variant caches with one serial run per variant present
        warmed = Set{String}()
        for cn in torun
            v = varmap[cn]; v in warmed && continue
            push!(warmed, v)
            r = _run_one(cn, v, subdb, plan, mgmt); _tick!()
            r === nothing && continue
            results[cn] = r
            lock(_MEM_LOCK) do; _MEM["$cn|$ph"] = r; end
        end

        rest = filter(cn -> !haskey(results, cn), torun)
        out = Vector{Union{Nothing,PlotResult}}(undef, length(rest))
        Threads.@threads for i in eachindex(rest)
            out[i] = _run_one(rest[i], varmap[rest[i]], subdb, plan, mgmt); _tick!()
        end
        for (i, cn) in enumerate(rest)
            out[i] === nothing && continue
            results[cn] = out[i]
            lock(_MEM_LOCK) do; _MEM["$cn|$ph"] = out[i]; end
        end
        return results
    finally
        _PROG.running[] = false
    end
end

# ---------------------------------------------------------------------------
# PPE cross-stand harvest budget (landscape flow) — FVS Parallel Processing
# Extension MXHRVP. Instead of scheduling each plot independently, the user sets a
# landscape harvest POLICY (a target resource flow + per-stand priority + credit) and
# `ppe_run_landscape_harvest!` picks WHICH plots to cut each master cycle to meet the
# target, in priority order. Policy variables: BBA/AGE (per-stand), SELECTED/AVBBA/
# TOTALWT (landscape). Expressions reuse the Event-Monitor evaluator.
# ---------------------------------------------------------------------------
struct LandscapePolicy
    common_year::Int         # all plots are harmonized to this inventory year so the
    target_expr::String      #   master-cycle grid aligns across plots (FIA invyears vary)
    priority_expr::String    # per-stand cut priority (e.g. "BBA")
    credit_expr::String      # per-stand resource credit toward the target (e.g. "BBA")
    cycles::Int
    period::Int
    mslabel::String
end
LandscapePolicy(; common_year=2020, target_expr="1000", priority_expr="BBA",
                  credit_expr="BBA", cycles=3, period=10, mslabel="ALL") =
    LandscapePolicy(common_year, target_expr, priority_expr, credit_expr, cycles, period, mslabel)

"""
    simulate_landscape(plt_cns, acres, pol; cache_dir) -> NamedTuple

Run the PPE landscape harvest over the AOI's plots (area-weighted by `acres`, a
Dict plt_cn→acres). Because the PPE reader ingests every stand with one variant, v1
runs the AOI's DOMINANT variant and reports minority-variant plots as excluded.
Returns `(; variant, nplots, nexcluded, master_years, cycles)` where `cycles` is a
per-master-year selection table: which plots are cut, resource achieved vs target.
"""
function simulate_landscape(plt_cns::Vector{String}, acres::AbstractDict{String,Float64},
                            pol::LandscapePolicy; cache_dir::AbstractString)
    plt_cns = unique(plt_cns)
    subdb, varmap = build_subdb(plt_cns; cache_dir=cache_dir)

    # dominant variant (the PPE reads all stands with one variant; don't mis-read others)
    counts = Dict{String,Int}()
    for cn in plt_cns
        v = get(varmap, cn, ""); isempty(v) && continue
        counts[v] = get(counts, v, 0) + 1
    end
    isempty(counts) && error("no projectable plots in the AOI")
    vdom = argmax(counts)
    cns = [cn for cn in plt_cns if get(varmap, cn, "") == vdom]
    nexcluded = length(plt_cns) - length(cns)

    # harmonize INV_YEAR into a policy-specific DB copy so all plots share the grid
    harm = joinpath(cache_dir, "landscape_$(vdom)_$(pol.common_year).db")
    cp(subdb, harm; force=true); chmod(harm, 0o644)
    db = SQLite.DB(harm)
    DBInterface.execute(db, "UPDATE FVS_STANDINIT_PLOT SET INV_YEAR = $(pol.common_year)")
    close(db)

    keyfor(cn) = begin
        k = tempname() * ".key"
        write(k, string("STDIDENT\n", cn, "\n",
            "DATABASE\nDSNin\n", harm, "\n",
            "StandSQL\nSELECT * FROM FVS_STANDINIT_PLOT WHERE STAND_CN = '%StandID%'\nEndSQL\n",
            "TreeSQL\nSELECT * FROM FVS_TREEINIT_PLOT WHERE STAND_CN = '%StandID%'\nEndSQL\nEND\n",
            _kwrec("NUMCYCLE", pol.cycles), "\nPROCESS\nSTOP\n"))
        k
    end
    stands = FVSjl.PPEStand[]; ids = String[]
    for cn in cns
        push!(stands, FVSjl.PPEStand(keyfor(cn); area = get(acres, cn, 1.0)))
        push!(ids, cn)
    end
    lbl = fill(pol.mslabel, length(stands))
    master_years = [pol.common_year + pol.period * k for k in 0:pol.cycles]

    res = FVSjl.ppe_run_landscape_harvest!(stands;
            variant = FVSjl.variant_from_code(vdom),
            labels = lbl, mslabel = pol.mslabel,
            target_expr = pol.target_expr, priority_expr = pol.priority_expr,
            credit_expr = pol.credit_expr, master_years = master_years,
            period = pol.period, lprtct = true)

    cycles = map(res) do r
        (; year = r.year, elapsed = r.year - pol.common_year,
           target = round(r.target; digits = 1),
           resource = round(r.selected_resource; digits = 1),
           pct_of_target = round(r.pct_of_target; digits = 1),
           hvpart = round(r.hvpart; digits = 3),
           cut = String[r.stand_ids[i] for i in eachindex(r.selected) if r.selected[i]],
           ncut = count(r.selected), nstands = length(r.selected))
    end
    (; variant = vdom, nplots = length(cns), nexcluded = nexcluded,
       master_years = master_years, cycles = cycles)
end

end # module
