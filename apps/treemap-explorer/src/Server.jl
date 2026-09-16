# Included into TreeMapExplorer. HTTP layer (Oxygen): dynamic tiles + AOI API.

using Oxygen, HTTP, JSON3
import PNGFiles

const APP = Ref{App}()
const PUBLIC = normpath(joinpath(@__DIR__, "..", "public"))

# raster geographic bounds (EPSG:4326) for map fit — CONUS Albers extent corners
const CONUS_BOUNDS_4326 = (west = -127.9, south = 22.7, east = -65.3, north = 51.6)

# Extension points: data sources and simulation engines. Today one of each is
# implemented; the UI presents them as selectable so more can be added later.
const SOURCES = [
    (; id = "treemap2022", name = "TreeMap 2022 CONUS", detail = "USFS · 30 m · circa 2022",
       active = true)]
const ENGINES = [
    (; id = "fvsjl", name = "FVSjl", detail = "Forest Vegetation Simulator (Julia port)",
       active = true)]

_png_bytes(img) = (io = IOBuffer(); PNGFiles.save(io, img); take!(io))

function _cmap_name(cm)
    cm === Color.GREENS  && return "greens"
    cm === Color.MAGMA   && return "magma"
    return "viridis"
end

function _tile_handler(attr::String, z::Int, x::Int, y::Int)
    app = APP[]
    spec = get(ATTRS, attr, ATTRS[DEFAULT_ATTR])
    tmids = Raster.warp_tile(app.rasterctx, z, x, y; size = 256)
    valfn = attr_valfn(app, attr)
    img = Color.colorize(tmids, valfn, spec.cmap, spec.lo, spec.hi)
    HTTP.Response(200, ["Content-Type" => "image/png",
                        "Cache-Control" => "public, max-age=86400"];
                  body = _png_bytes(img))
end

function _aoi_response(geom5070)
    app = APP[]
    counts = Raster.aoi_tally(app.rasterctx, geom5070)
    agg = Aggregate.aggregate_aoi(counts, app.tt, app.st)
    cells = [(; spcd = c.spcd,
                symbol = get(agg.species_symbol, c.spcd, ""),
                common = get(agg.species_name, c.spcd, ""),
                dbh_lo = c.dbh_lo,
                count = round(c.count; digits = 2),
                ba = round(c.ba; digits = 2)) for c in agg.cells]
    w, s, e, n = Aoi.bbox4326(geom5070)
    (; nplots = agg.nplots, npixels = agg.npixels,
       bbox = (; west = w, south = s, east = e, north = n),
       geometry = JSON3.read(Aoi.geojson4326(geom5070)),  # 4326 GeoJSON for the map

       acres = round(agg.acres; digits = 1),
       ba_total = round(agg.ba_total; digits = 1),
       tpa_expanded = round(agg.tpa_expanded; digits = 0),
       carbon_l = round(agg.carbon_l; digits = 1),
       carbon_d = round(agg.carbon_d; digits = 1),
       carbon_dwn = round(agg.carbon_dwn; digits = 1),
       drybio_l = round(agg.drybio_l; digits = 1),
       volcfnet_l = round(agg.volcfnet_l; digits = 0),
       cells = cells)
end

# --- Phase 2: run FVSjl on the AOI's plots under a management plan -----------
const MAX_SIM_PLOTS = 200          # cap for a responsive synchronous request
const HA_PER_ACRE = 0.404685642

# the projected layer of the last simulation: tm_id => per-cycle metrics (for /simtiles)
const SIM_LAYER = Ref{Dict{UInt32,Vector{FiaSim.CycleMetrics}}}(Dict{UInt32,Vector{FiaSim.CycleMetrics}}())
# map-colorable per-area sim metrics: name => (field, lo, hi, colormap)
const SIM_METRICS = Dict(
    "carbon" => (:carbon_live_ag, 0.0, 80.0,  Color.GREENS),   # t C/ha live aboveground
    "ba"     => (:ba,             0.0, 250.0, Color.VIRIDIS),  # ft²/ac
    "volume" => (:tcuft,          0.0, 8000.0, Color.VIRIDIS)) # ft³/ac standing

"AOI geometry (5070) -> (Dict{PLT_CN=>acres}, Dict{tm_id=>PLT_CN}), capped to largest plots."
function _resolve_plots(geom5070)
    app = APP[]
    counts = Raster.aoi_tally(app.rasterctx, geom5070)      # tm_id => pixels
    acres = Dict{String,Float64}()
    tm2cn = Dict{UInt32,String}()
    for (tm, px) in counts
        rec = VAT.stand(app.st, tm)
        rec === nothing && continue
        cn = string(rec.plt_cn)
        tm2cn[tm] = cn
        acres[cn] = get(acres, cn, 0.0) + px * Aggregate.ACRES_PER_PIXEL
    end
    acres, tm2cn
end

function _parse_plan(o)
    acts = FiaSim.PlanAction[]
    for a in get(o, :actions, ())
        kind = String(a.kind)
        if kind == "thin"
            push!(acts, FiaSim.PlanAction(; kind = "thin", cycle = Int(get(a, :cycle, 1)),
                metric = String(get(a, :metric, "BA")), target = Float64(get(a, :target, 0)),
                direction = String(get(a, :direction, "below")),
                dbh_lo = Float64(get(a, :dbh_lo, 0)), dbh_hi = Float64(get(a, :dbh_hi, 999))))
        elseif kind == "plant"
            push!(acts, FiaSim.PlanAction(; kind = "plant", cycle = Int(get(a, :cycle, 1)),
                species = String(a.species), tpa = Float64(a.tpa),
                survival = Float64(get(a, :survival, 100))))
        end
    end
    FiaSim.ManagementPlan(; name = String(get(o, :name, "plan")),
        cycles = Int(get(o, :cycles, 10)), period = Int(get(o, :period, 10)), actions = acts)
end

function _simulate_response(geom5070, plan)
    app = APP[]
    acres, tm2cn = _resolve_plots(geom5070)
    total_aoi_acres = sum(values(acres); init = 0.0)
    # cap to the largest-area plots for a bounded synchronous run
    cns = sort(collect(keys(acres)); by = cn -> -acres[cn])
    capped = length(cns) > MAX_SIM_PLOTS
    sim_cns = capped ? cns[1:MAX_SIM_PLOTS] : cns
    cache = joinpath(app.datadir, "derived", "sim_cache")
    res = FiaSim.simulate_plots(sim_cns, plan; cache_dir = cache)   # plt_cn => Vector{CycleMetrics}

    ncyc = maximum((length(v) for v in values(res)); init = 0)
    sim_acres = sum(acres[cn] for cn in keys(res); init = 0.0)
    cycles = map(1:ncyc) do i
        baw = 0.0; w = 0.0; vol = 0.0; rem = 0.0; ct = 0.0; cl = 0.0; el = (i - 1) * plan.period
        for (cn, v) in res
            i <= length(v) || continue
            m = v[i]; a = acres[cn]
            baw += m.ba * a; w += a
            vol += m.tcuft * a
            rem += m.rem_tcuft * a
            isnan(m.carbon_total)   || (ct += m.carbon_total   * a * HA_PER_ACRE)
            isnan(m.carbon_live_ag) || (cl += m.carbon_live_ag * a * HA_PER_ACRE)
        end
        (; cycle = i - 1, elapsed = el,
           ba_ac = w > 0 ? round(baw / w; digits = 1) : 0.0,
           volume = round(vol; digits = 0),
           removed_volume = round(rem; digits = 0),
           carbon_total = round(ct; digits = 1),
           carbon_live_ag = round(cl; digits = 1))
    end
    # build the per-pixel projected layer (tm_id => per-cycle metrics) for the time slider
    layer = Dict{UInt32,Vector{FiaSim.CycleMetrics}}()
    for (tm, cn) in tm2cn
        v = get(res, cn, nothing)
        v === nothing || (layer[tm] = v)
    end
    SIM_LAYER[] = layer

    (; plan = plan.name, nplots_aoi = length(acres), nplots_sim = length(res),
       capped = capped, acres = round(total_aoi_acres; digits = 1),
       sim_acres = round(sim_acres; digits = 1),
       ncycles = ncyc, period = plan.period, cycles = cycles)
end

"Colorized projected tile for the last simulation at `cycle` (0-based), metric per-area."
function _simtile_handler(metric::String, cycle::Int, z::Int, x::Int, y::Int)
    layer = SIM_LAYER[]
    isempty(layer) && return HTTP.Response(204)
    field, lo, hi, cmap = get(SIM_METRICS, metric, SIM_METRICS["carbon"])
    tmids = Raster.warp_tile(APP[].rasterctx, z, x, y; size = 256)
    valfn = tm -> begin
        v = get(layer, tm, nothing)
        (v === nothing || cycle + 1 > length(v)) && return NaN
        Float64(getfield(v[cycle+1], field))
    end
    img = Color.colorize(tmids, valfn, cmap, lo, hi)
    HTTP.Response(200, ["Content-Type" => "image/png", "Cache-Control" => "no-store"];
                  body = _png_bytes(img))
end

function start_server!(app::App = init_app(); host = "127.0.0.1", port = 8080)
    APP[] = app

    @get "/tiles/{attr}/{z}/{x}/{y}" function (req, attr::String, z::String, x::String, y::String)
        try
            yi = parse(Int, first(split(y, '.')))   # strip ".png"
            return _tile_handler(attr, parse(Int, z), parse(Int, x), yi)
        catch e
            @warn "tile error" attr z x y exception = e
            return HTTP.Response(204)   # empty tile
        end
    end

    @get "/simtiles/{metric}/{cycle}/{z}/{x}/{y}" function (req, metric::String, cycle::String, z::String, x::String, y::String)
        try
            yi = parse(Int, first(split(y, '.')))
            return _simtile_handler(metric, parse(Int, cycle), parse(Int, z), parse(Int, x), yi)
        catch e
            @warn "simtile error" exception = e
            return HTTP.Response(204)
        end
    end

    @get "/api/attrs" function (req)
        [(; name = s.name, label = s.label, units = s.units, lo = s.lo, hi = s.hi,
            cmap = _cmap_name(s.cmap))
         for s in sort(collect(values(ATTRS)); by = s -> s.name)]
    end

    @get "/api/meta" function (req)
        (; default_attr = DEFAULT_ATTR, bounds = CONUS_BOUNDS_4326,
           nplots_total = length(app.st.by_tm),
           sources = SOURCES, engines = ENGINES)
    end

    # AOI: JSON body {"geometry": <geojson geometry>} (drawn), or raw file upload
    @post "/api/aoi" function (req)
        ct = HTTP.header(req, "Content-Type", "")
        try
            if occursin("application/json", ct)
                body = JSON3.read(String(req.body))
                geom = Aoi.geom_from_geojson(JSON3.write(body.geometry))
                return _aoi_response(geom)
            else
                fname = HTTP.header(req, "X-Filename", "aoi.geojson")
                geom = Aoi.geom_from_upload(Vector{UInt8}(req.body), fname)
                return _aoi_response(geom)
            end
        catch e
            @warn "aoi error" exception = (e, catch_backtrace())
            return HTTP.Response(400, JSON3.write((; error = string(e))))
        end
    end

    # Phase 2: run a management plan on the AOI's plots -> per-cycle aggregate
    @post "/api/simulate" function (req)
        try
            body = JSON3.read(String(req.body))
            geom = Aoi.geom_from_geojson(JSON3.write(body.geometry))
            plan = _parse_plan(body.plan)
            return _simulate_response(geom, plan)
        catch e
            @warn "simulate error" exception = (e, catch_backtrace())
            return HTTP.Response(400, JSON3.write((; error = string(e))))
        end
    end

    # static frontend
    staticfiles(PUBLIC, "/")

    @info "TreeMap Growth Explorer serving — press Ctrl-C to stop" host port public = PUBLIC

    # Run Oxygen ASYNC and own the interrupt ourselves. Oxygen's own async=false
    # path flips exit_on_sigint off and blocks in wait(); the HTTP accept loop then
    # keeps the process alive so Ctrl-C appears ignored. Here SIGINT is delivered to
    # our sleep loop as an InterruptException, we terminate() and return cleanly.
    serve(; host = host, port = port, async = true)
    Base.exit_on_sigint(false)
    try
        while true
            sleep(0.5)
        end
    catch e
        e isa InterruptException || rethrow()
        @info "shutting down"
    finally
        try
            terminate()
        catch
        end
    end
    return nothing
end
