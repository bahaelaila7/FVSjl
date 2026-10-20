"""
    PanEngine

Second simulation-engine family for Forest Explorer: **Pan** (Biomass Succession, PnET-Succession, HBAD v1, HBAD v2) behind the
same `simulate` seam as FVSjl. Pan runs in its own Julia process (its own project/depot — the two Manifests are not
compatible), started lazily and kept warm: `tools/pan_worker.jl` of the Pan repo speaks JSON lines over stdin/stdout
(protocol: Pan repo `docs/scenario_campaign_plan.md`). Results are converted to `FiaSim.PlotResult`, so aggregation,
rasters, charts and the compare report are shared with FVSjl.

Config (env): `PAN_REPO` (default `<workspace root>/landis_fe`: the refactor branch with the scenario engine and HBAD v2), `PAN_THREADS` (default `auto`), `PAN_AUTOSTART` (`1`
starts the worker in the background when the server starts), `PAN_WORKER_LOG` (stderr of the worker).
"""
module PanEngine

using JSON3
import ..FiaSim

export pan_status, pan_call, pan_simulate, pan_load_params, pan_default_params, pan_climate_preview

pan_repo() = get(ENV, "PAN_REPO", joinpath(get(ENV, "WORKSPACE", normpath(joinpath(@__DIR__, "..", "..", "..", ".."))), "landis_fe"))   # default: <workspace root>/landis_fe (root = the folder holding fe_pan, or WORKSPACE)
worker_log() = get(ENV, "PAN_WORKER_LOG", joinpath(tempdir(), "pan_worker.log"))

const ENGINE_PLUGIN = Dict("pan-biomass" => "biomass", "pan-pnet" => "pnet", "pan-hbad" => "hbad", "pan-hbad2" => "hbad2")
# `coming = true`: listed in the engine selector but disabled (no Pan plugin yet); the frontend shows the label as is.
const PAN_ENGINES = [
    (; id = "pan-biomass", label = "Pan · Biomass Succession", plugin = "biomass", climate = false, coming = false,
       detail = "Pan's Julia re-implementation of the LANDIS-II Biomass Succession extension (LANDIS-II: Scheller et al. 2007, Ecological Modelling 201:409–419; Biomass Succession: Scheller & Mladenoff 2004, Ecological Modelling 180:211–229). Species-age cohorts of aboveground biomass; growth limited by ANPP and competition for growing space (B_MAX), age- and growth-related mortality, establishment by shade tolerance and seed availability. Parameters calibrated to FIA in Pan."),
    (; id = "pan-pnet", label = "Pan · PnET-Succession", plugin = "pnet", climate = true, coming = false,
       detail = "Pan's re-implementation of the LANDIS-II PnET-Succession extension (de Bruijn et al. 2014, Ecological Modelling 287:44–57; built on PnET-II, Aber & Federer 1992). Cohort photosynthesis driven by light, temperature, water (soil water balance), VPD and CO₂; carbon allocation to foliage/wood/roots; responds to the climate scenario."),
    (; id = "pan-hbad", label = "Pan · HBAD v1", plugin = "hbad", climate = false, coming = false,
       detail = "Pan cohort model of height, basal area and density (N, QMD, height, age) with structural competition (basal area of larger trees, relative density, Reineke self-thinning) and empirical FIA growth/mortality regressions; establishment as in Biomass Succession. Biomass reported via allometry."),
    (; id = "pan-hbad2", label = "Pan · HBAD v2 (v2j)", plugin = "hbad2", climate = false, coming = false,
       detail = "Biomass-free process cohort model: state N (trees/ha), QMD, mean height, age; height-structured light competition (overtopping) + crowding; species growth curves for diameter and height; growth-stress, background and age mortality; recruitment from seed source and light; calibrated jointly on FIA remeasured trees (Sim A) and stand trajectories (Sim B). Default parameter set: v2j."),
    (; id = "pan-hbad3", label = "HBAD3 (process-based) — coming", plugin = "hbad3", climate = false, coming = true,
       detail = "Process-grounded successor of HBAD v2 (stem-volume growth with pipe-model leaf area and hydraulic limitation, trait-based height, carbon-starvation mortality, water balance). Not in Pan yet.")]

mutable struct Worker
    proc::Base.Process
    next::Int
    info::Any
end
const W = Ref{Union{Nothing,Worker}}(nothing)
const LOCK = ReentrantLock()
const STATE = Ref{Symbol}(:stopped)          # :stopped | :starting | :ready | :failed
const LAST_ERROR = Ref{String}("")
const YAML_BY_ID = Dict{String,String}()     # params_id => YAML text (re-sent if the worker restarts)

function _start!()
    repo = pan_repo()
    isfile(joinpath(repo, "tools", "pan_worker.jl")) || error("Pan repo not found at $repo (set PAN_REPO)")
    th = get(ENV, "PAN_THREADS", "auto")
    # depot: PAN_DEPOT_PATH (set by serve_pan.sh) > the checkout's wt_env.sh > julia_gdal.sh's own default (repo depot / ~/.julia)
    script = "cd '$repo' && if [ -n \"\${PAN_DEPOT_PATH:-}\" ]; then export JULIA_DEPOT_PATH=\"\$PAN_DEPOT_PATH\"; " *
             "elif [ -f wt_env.sh ]; then source wt_env.sh; else unset JULIA_DEPOT_PATH; fi && " *
             "exec ./julia_gdal.sh --project=. --threads=$th tools/pan_worker.jl"
    STATE[] = :starting
    proc = open(pipeline(`bash -c $script`; stderr = worker_log()), "r+")
    line = eof(proc) ? "" : readline(proc)    # {"id":null,"ok":true,"result":{"ready":true}}
    isempty(line) && (STATE[] = :failed; LAST_ERROR[] = "Pan worker exited during startup (see $(worker_log()))"; error(LAST_ERROR[]))
    w = Worker(proc, 0, nothing)
    W[] = w
    w.info = _raw_call(w, "ping")
    STATE[] = :ready
    @info "Pan worker ready" repo info = w.info
    w
end

function _raw_call(w::Worker, op::String; on_progress = nothing, kw...)
    w.next += 1; id = w.next
    req = Dict{String,Any}("id" => id, "op" => op)
    for (k, v) in kw; req[String(k)] = v; end
    println(w.proc, JSON3.write(req)); flush(w.proc)
    while true
        if eof(w.proc)              # worker exited (crash / killed): fail this request; the next call restarts it
            W[] = nothing; STATE[] = :failed
            LAST_ERROR[] = "Pan worker exited (see $(worker_log())); it will be restarted on the next request"
            try; kill(w.proc); catch; end
            error(LAST_ERROR[])
        end
        line = readline(w.proc)
        isempty(strip(line)) && continue
        o = JSON3.read(line)
        if haskey(o, :progress)
            on_progress === nothing || on_progress(o.progress.done, o.progress.total)
            continue
        end
        o.id == id || continue
        o.ok || error(String(o.error))
        return o.result
    end
end

"Call the worker (starting it if needed). Serialized: one request at a time."
function pan_call(op::String; on_progress = nothing, kw...)
    lock(LOCK) do
        w = W[]
        (w === nothing || !process_running(w.proc)) && (w = _start!())
        _raw_call(w, op; on_progress = on_progress, kw...)
    end
end

"Start the worker in the background (warm-up: Pan load + compile takes ~1–2 min)."
autostart!() = (get(ENV, "PAN_AUTOSTART", "1") == "1" && isfile(joinpath(pan_repo(), "tools", "pan_worker.jl"))) &&
    Threads.@spawn try
        pan_call("ping")
    catch e
        LAST_ERROR[] = sprint(showerror, e); STATE[] = :failed
        @warn "Pan worker failed to start" exception = e
    end

function pan_status()
    repo_ok = isfile(joinpath(pan_repo(), "tools", "pan_worker.jl"))
    msg = !repo_ok ? "Pan repository not found at $(pan_repo()) (set PAN_REPO)" :
          STATE[] === :failed ? LAST_ERROR[] :
          STATE[] === :starting ? "Pan engine is starting (first start compiles, ~1–2 min)" : ""
    (; available = repo_ok, state = String(STATE[]), engines = PAN_ENGINES, message = msg)   # a failed worker restarts on the next call
end

_desc(r) = (; params_id = String(r.params_id), plugin = String(r.plugin), name = String(r.name),
              description = String(r.description), n_species = Int(r.n_species), n_ecos = Int(r.n_ecos),
              region = r.region, provenance = r.provenance, species = String[String(x) for x in r.species],
              warnings = String[String(x) for x in r.warnings])   # plain Vector{String}: an empty JSON3 array would serialize as {}

function pan_load_params(yaml::AbstractString)
    r = pan_call("load_params"; yaml = String(yaml))
    YAML_BY_ID[String(r.params_id)] = String(yaml)
    _desc(r)
end
pan_default_params(plugin::AbstractString) = _desc(pan_call("default_params"; plugin = String(plugin)))

# re-send an uploaded YAML when the worker no longer knows its params_id (restart)
function _with_params(f, params_id)
    try
        return f()
    catch e
        msg = sprint(showerror, e)
        if params_id !== nothing && occursin("unknown params_id", msg) && haskey(YAML_BY_ID, params_id)
            pan_call("load_params"; yaml = YAML_BY_ID[params_id])
            return f()
        end
        rethrow()
    end
end

"Plot records for the worker: CN, acres and — when known — the AOI location of its pixels (l4, lon, lat)."
function _plot_records(cns, acres, sites)
    map(cns) do cn
        d = Dict{String,Any}("cn" => cn, "acres" => Float64(acres[cn]))
        s = sites === nothing ? nothing : get(sites, cn, nothing)
        s === nothing || (d["l4"] = s.l4; d["lon"] = s.lon; d["lat"] = s.lat)
        d
    end
end

"""
    pan_simulate(cns, acres, plan_json, engine; params_id=nothing, climate=nothing, sites=nothing) -> (Dict{cn => FiaSim.PlotResult}, warnings, engine_info, status_counts)
"""
function pan_simulate(cns::Vector{String}, acres::AbstractDict, plan, engine::AbstractString; params_id = nothing, climate = nothing,
                      sites = nothing)
    plugin = get(ENGINE_PLUGIN, engine, nothing)
    plugin === nothing && error("unknown Pan engine $engine")
    plots = _plot_records(cns, acres, sites)
    prog(d, n) = (FiaSim._PROG.total[] = n; FiaSim._PROG.done[] = d)
    FiaSim._PROG.total[] = 1; FiaSim._PROG.done[] = 0; FiaSim._PROG.running[] = true
    r = try
        _with_params(params_id) do
            pan_call("simulate"; on_progress = prog, plugin = plugin, params_id = params_id, plots = plots, plan = plan,
                     climate = (plugin == "pnet" ? climate : nothing))
        end
    finally
        FiaSim._PROG.running[] = false
    end
    String(r.engine_info.plugin) == plugin ||
        error("the parameter set is for plugin '$(r.engine_info.plugin)', but the engine is $engine")
    out = Dict{String,FiaSim.PlotResult}()
    counts = Dict{String,Int}()
    for (cnk, v) in pairs(r.plots)
        st = String(v.status); counts[st] = get(counts, st, 0) + 1
        st == "ok" || continue
        mets = [FiaSim.CycleMetrics(Int(m.year), Float64(m.ba), Float64(m.tpa), Float64(m.qmd), Float64(m.tcuft), Float64(m.bdft),
                                    Float64(m.rem_tcuft), Float64(m.rem_bdft), Float64(m.carbon_total), Float64(m.carbon_live_ag),
                                    Float64(m.carbon_removed)) for m in v.metrics]
        dist = Dict{Int,Dict{Tuple{Int16,Int},NTuple{2,Float64}}}()
        for (yk, rows) in pairs(v.dist)
            g = Dict{Tuple{Int16,Int},NTuple{2,Float64}}()
            for x in rows; g[(Int16(x[1]), Int(x[2]))] = (Float64(x[3]), Float64(x[4])); end
            dist[parse(Int, String(yk))] = g
        end
        spsym = Dict{Int16,String}(Int16(parse(Int, String(k))) => String(s) for (k, s) in pairs(v.spsym))
        out[String(cnk)] = FiaSim.PlotResult(mets, dist, spsym)
    end
    warnings = String.(r.warnings)
    nbad = sum(v for (k, v) in counts if k != "ok"; init = 0)
    nbad > 0 && pushfirst!(warnings, "$(nbad) of $(length(cns)) plots not simulated by Pan: " *
        join(["$(v) $(k)" for (k, v) in counts if k != "ok"], ", ") * " (out_of_region = the AOI location is outside the parameter set's ecoregions; no_data = no FIA cohorts and no donor)")
    out, warnings, r.engine_info, counts
end

pan_climate_preview(cns, acres, climate; start_year = 2025, years = 100, sites = nothing) =
    pan_call("climate_preview"; plots = _plot_records(cns, acres, sites),
             climate = climate, start_year = start_year, years = years)

end
