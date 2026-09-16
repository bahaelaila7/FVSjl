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

export ManagementPlan, PlanAction, simulate_plots, build_subdb, CycleMetrics

# resolved at RUNTIME so FIA_DB set when the server starts is honored (not baked at precompile)
master_db() = get(ENV, "FIA_DB", "/workspace/SQLite_FIADB_ENTIRE.db")

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
        tok = get(_THIN_TOKEN, (a.metric, a.direction), "THINBBA")
        # THIN*: field1=DATE (<1000 ⇒ CYCLE number, relative; our 0-based cycle → FVS
        # cycle a.cycle+1, so cycle 0 = at inventory, cycle 1 = +period …),
        # field2=residual target, field3=cut-eff, field4/5=DBH lo/hi
        println(io, _kwrec(tok, a.cycle + 1, a.target, 1, a.dbh_lo, a.dbh_hi))
    end
    if !isempty(plants)
        println(io, "ESTAB")
        for a in plants
            # PLANT: DATE(cycle a.cycle+1), species(alpha), trees/acre, survival%
            println(io, _kwrec("PLANT", a.cycle + 1, a.species, a.tpa, a.survival))
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

const _MEM = Dict{String,Vector{CycleMetrics}}()   # cache: "plt_cn|planhash" => cycles
const _MEM_LOCK = ReentrantLock()

"Run one plot; returns its per-cycle metrics or `nothing`."
function _run_one(cn::String, vcode::String, subdb::String, plan::ManagementPlan, mgmt::String)
    key = tempname() * ".key"
    write(key, string(
        "STDIDENT\n", cn, "\n",
        "DATABASE\nDSNin\n", subdb, "\n",
        "StandSQL\nSELECT * FROM FVS_STANDINIT_PLOT WHERE STAND_CN = '%StandID%'\nEndSQL\n",
        "TreeSQL\nSELECT * FROM FVS_TREEINIT_PLOT WHERE STAND_CN = '%StandID%'\nEndSQL\n",
        "END\n",
        _kwrec("NUMCYCLE", plan.cycles), "\n",
        mgmt,
        "FMIn\nCARBREPT\nEnd\n",
        "ECHOSUM\nPROCESS\nSTOP\n"))
    try
        out = FVSjl.run_keyfile(key; variant=FVSjl.variant_from_code(vcode), period=plan.period)
        return _parse_sum_and_carbon(out)
    catch e
        @warn "simulate failed" cn vcode exception=e
        return nothing
    finally
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
    subdb, varmap = build_subdb(plt_cns; cache_dir=cache_dir)
    ph = plan_hash(plan)
    mgmt = plan_keywords(plan)
    results = Dict{String,Vector{CycleMetrics}}()
    torun = String[]
    for cn in plt_cns
        r = lock(_MEM_LOCK) do; get(_MEM, "$cn|$ph", nothing); end
        r !== nothing ? (results[cn] = r) :
            (isempty(get(varmap, cn, "")) || push!(torun, cn))
    end
    isempty(torun) && return results

    # warm up lazy per-variant caches with one serial run per variant present
    warmed = Set{String}()
    for cn in torun
        v = varmap[cn]; v in warmed && continue
        push!(warmed, v)
        r = _run_one(cn, v, subdb, plan, mgmt)
        r === nothing && continue
        results[cn] = r
        lock(_MEM_LOCK) do; _MEM["$cn|$ph"] = r; end
    end

    rest = filter(cn -> !haskey(results, cn), torun)
    out = Vector{Union{Nothing,Vector{CycleMetrics}}}(undef, length(rest))
    Threads.@threads for i in eachindex(rest)
        out[i] = _run_one(rest[i], varmap[rest[i]], subdb, plan, mgmt)
    end
    for (i, cn) in enumerate(rest)
        out[i] === nothing && continue
        results[cn] = out[i]
        lock(_MEM_LOCK) do; _MEM["$cn|$ph"] = out[i]; end
    end
    results
end

end # module
