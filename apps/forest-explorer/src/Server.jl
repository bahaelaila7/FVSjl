# Included into ForestExplorer. HTTP layer (Oxygen): dynamic tiles + AOI API.

using Oxygen, HTTP, JSON3
import PNGFiles
import ArchGDAL

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
# PPE harmonises every stand onto one inventory year so the budget's cycles line up across plots
# (FiaSim.simulate_landscape rewrites INV_YEAR in a policy-specific DB copy).
const _PPE_YEAR = 2020
const HA_PER_ACRE = 0.404685642

# the projected layer of the last simulation: tm_id => per-cycle metrics (source for the
# cached AOI images rendered below)
const SIM_LAYER = Ref{Dict{UInt32,Vector{FiaSim.CycleMetrics}}}(Dict{UInt32,Vector{FiaSim.CycleMetrics}}())
# map-colorable per-area sim metrics: name => (field, lo, hi, colormap). lo/hi are CONUS
# fallbacks; the actual domain is auto-scaled per simulation (SIM_DOMAINS) so a pixel's
# change across cycles spans a visible color range.
const SIM_METRICS = Dict(
    "carbon" => (:carbon_live_ag, 0.0, 80.0,  Color.GREENS),   # t C/ha live aboveground
    "ba"     => (:ba,             0.0, 250.0, Color.VIRIDIS),  # ft²/ac
    "volume" => (:tcuft,          0.0, 8000.0, Color.VIRIDIS)) # ft³/ac standing
const SIM_DOMAINS = Ref{Dict{String,Tuple{Float64,Float64}}}(Dict{String,Tuple{Float64,Float64}}())
# the last simulation's AOI polygon in EPSG:3857 (tile CRS) — clips the rendered image to
# the polygon, since a plot's tm_id recurs on many pixels nationwide (imputation).
const SIM_AOI = Ref{Any}(nothing)
# Pre-rendered, cached AOI raster images: (metric, cycle) => PNG bytes. The whole AOI is
# warped + masked ONCE per simulation, then colorized per (metric, cycle) and cached, so
# the time slider just swaps a ready image (an /simimage lookup) with no recompute.
# Per-token cache: token => ((metric,cycle) => PNG bytes). Each simulation (each SCENARIO)
# gets a fresh token so multiple scenarios' rasters coexist; an LRU cap bounds memory.
const SIM_IMAGES = Ref{Dict{Int,Dict{Tuple{String,Int},Vector{UInt8}}}}(Dict{Int,Dict{Tuple{String,Int},Vector{UInt8}}}())
const SIM_IMG_ORDER = Ref{Vector{Int}}(Int[])   # token insertion order for LRU eviction
const SIM_IMG_LOCK = ReentrantLock()
const SIM_IMG_MAXTOKENS = 12                     # keep the last N scenarios' raster sets
const SIM_IMG_META = Ref{Any}(nothing)      # corners (lon/lat) + size + metrics for the frontend
const SIM_TOKEN = Threads.Atomic{Int}(0)    # bumps per simulation → fresh immutable image URLs
const SIM_IMG_MAXPX = 1200                   # cap the longer image side (bounds memory)

"Store a scenario's rendered images under `tok`, evicting the oldest token past the cap."
function _cache_sim_images!(tok::Int, images::Dict{Tuple{String,Int},Vector{UInt8}})
    lock(SIM_IMG_LOCK) do
        SIM_IMAGES[][tok] = images
        push!(SIM_IMG_ORDER[], tok)
        while length(SIM_IMG_ORDER[]) > SIM_IMG_MAXTOKENS
            old = popfirst!(SIM_IMG_ORDER[])
            delete!(SIM_IMAGES[], old)
        end
    end
end

"EPSG:3857 meters -> (lon, lat) degrees."
function _merc_to_lonlat(x::Float64, y::Float64)
    lon = x / Raster.WEBMERC_HALF * 180.0
    lat = rad2deg(2 * atan(exp(y / 6_378_137.0)) - pi / 2)
    (lon, lat)
end

"""
Render the AOI into one cached PNG per (metric, cycle) and return the frontend metadata
(image corners in lon/lat, pixel size, the metrics that have data, and a fresh token).
The tm_id grid + AOI mask are computed ONCE; each (metric, cycle) is just a colorize.
"""
function _render_sim_images(geom3857, layer, doms, ncyc)
    env = ArchGDAL.envelope(geom3857)                 # 3857 meters: MinX,MaxX,MinY,MaxY
    xmin, xmax, ymin, ymax = env.MinX, env.MaxX, env.MinY, env.MaxY
    wm = max(xmax - xmin, 1.0); hm = max(ymax - ymin, 1.0)
    longer = max(wm, hm)
    # ~30 m native, capped; keep aspect ratio so the overlay isn't distorted
    L = clamp(round(Int, longer / 30), 128, SIM_IMG_MAXPX)
    w = max(1, round(Int, wm / longer * L)); h = max(1, round(Int, hm / longer * L))

    tmids = Raster.warp_extent(APP[].rasterctx, xmin, ymin, xmax, ymax, w, h)   # [row,col]
    mask = Raster.extent_mask(geom3857, xmin, ymin, xmax, ymax, w, h)
    @inbounds for i in eachindex(tmids)
        mask[i] == 0x00 && (tmids[i] = Raster.NODATA)
    end

    images = Dict{Tuple{String,Int},Vector{UInt8}}()
    for (mname, spec) in SIM_METRICS
        haskey(doms, mname) || continue               # skip metrics with no data (all-NaN)
        field, lo0, hi0, cmap = spec
        lo, hi = doms[mname]
        for c in 0:(ncyc - 1)
            valfn = tm -> begin
                v = get(layer, tm, nothing)
                (v === nothing || c + 1 > length(v)) && return NaN
                Float64(getfield(v[c + 1], field))
            end
            img = Color.colorize(tmids, valfn, cmap, lo, hi)
            images[(mname, c)] = _png_bytes(img)
        end
    end
    tok = Threads.atomic_add!(SIM_TOKEN, 1) + 1
    _cache_sim_images!(tok, images)          # per-token store (multi-scenario) + LRU cap
    # MapLibre image-source corners: TL, TR, BR, BL in [lon,lat]
    tl = _merc_to_lonlat(xmin, ymax); tr = _merc_to_lonlat(xmax, ymax)
    br = _merc_to_lonlat(xmax, ymin); bl = _merc_to_lonlat(xmin, ymin)
    meta = (; token = tok, width = w, height = h,
            corners = [collect(tl), collect(tr), collect(br), collect(bl)],
            metrics = [m for m in keys(SIM_METRICS) if haskey(doms, m)])
    SIM_IMG_META[] = meta
    meta
end

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

const ECO_L4_SHP = Ref(get(ENV, "ECO_L4_SHP", "/workspace/EPA_ECOREGIONS/us_eco_l4/us_eco_l4_no_st.shp"))

"Location of each plot inside the AOI (majority EPA L4 of its pixels + their centroid, lon/lat) — Pan's site, see Raster.aoi_sites."
function _plot_sites(geom5070, tm2cn)
    app = APP[]
    sites = try
        Raster.aoi_sites(app.rasterctx, geom5070, ECO_L4_SHP[])
    catch e
        @warn "AOI site lookup failed; Pan falls back to the FIA plots' own locations" exception = e
        return Dict{String,Any}()
    end
    acc = Dict{String,Any}()          # cn => (l4 counts, sx, sy, n) over its TM_IDs
    for (tm, s) in sites
        cn = get(tm2cn, tm, nothing); cn === nothing && continue
        d, sx, sy, n = get(acc, cn, (Dict{String,Int}(), 0.0, 0.0, 0))
        isempty(s.l4) || (d[s.l4] = get(d, s.l4, 0) + s.n)
        acc[cn] = (d, sx + s.x * s.n, sy + s.y * s.n, n + s.n)
    end
    out = Dict{String,Any}()
    src = ArchGDAL.importEPSG(5070); dst = ArchGDAL.importEPSG(4326)
    ArchGDAL.GDAL.osrsetaxismappingstrategy(dst.ptr, ArchGDAL.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    ArchGDAL.createcoordtrans(src, dst) do ct
        for (cn, (d, sx, sy, n)) in acc
            p = ArchGDAL.createpoint(sx / n, sy / n); ArchGDAL.transform!(p, ct)
            out[cn] = (l4 = isempty(d) ? "" : first(argmax(last, collect(d))), lon = ArchGDAL.getx(p, 0), lat = ArchGDAL.gety(p, 0))
        end
    end
    out
end

function _parse_plan(o)
    acts = FiaSim.PlanAction[]
    for a in get(o, :actions, ())
        kind = String(a.kind)
        if kind == "thin"
            push!(acts, FiaSim.PlanAction(; kind = "thin", cycle = Int(get(a, :cycle, 1)),
                metric = String(get(a, :metric, "BA")), target = Float64(get(a, :target, 0)),
                direction = String(get(a, :direction, "below")),
                species = String(get(a, :species, "")),
                dbh_lo = Float64(get(a, :dbh_lo, 0)), dbh_hi = Float64(get(a, :dbh_hi, 999)),
                mode = String(get(a, :mode, "rule")),
                target_expr   = String(get(a, :target_expr, "1000")),
                priority_expr = String(get(a, :priority_expr, "BBA")),
                credit_expr   = String(get(a, :credit_expr, "BBA"))))
        elseif kind == "plant"
            push!(acts, FiaSim.PlanAction(; kind = "plant", cycle = Int(get(a, :cycle, 1)),
                species = String(a.species), tpa = Float64(a.tpa),
                survival = Float64(get(a, :survival, 100))))
        end
    end
    FiaSim.ManagementPlan(; name = String(get(o, :name, "plan")),
        cycles = Int(get(o, :cycles, 10)), period = Int(get(o, :period, 10)), actions = acts)
end

function _parse_policy(o)
    FiaSim.LandscapePolicy(;
        common_year   = Int(get(o, :common_year, 2020)),
        target_expr   = String(get(o, :target_expr, "1000")),
        priority_expr = String(get(o, :priority_expr, "BBA")),
        credit_expr   = String(get(o, :credit_expr, "BBA")),
        cycles        = Int(get(o, :cycles, 3)),
        period        = Int(get(o, :period, 10)),
        mslabel       = String(get(o, :mslabel, "ALL")))
end

"PPE cross-stand harvest budget over the AOI's plots (dominant variant)."
function _landscape_response(geom5070, pol)
    app = APP[]
    acres, _ = _resolve_plots(geom5070)
    cns = collect(keys(acres))
    isempty(cns) && return (; error = "no plots resolved in the AOI")
    capped = length(cns) > MAX_SIM_PLOTS
    if capped
        cns = sort(cns; by = cn -> -acres[cn])[1:MAX_SIM_PLOTS]
    end
    cache = get(ENV, "FE_CACHE_DIR", joinpath(app.datadir, "derived", "sim_cache"))   # FE_CACHE_DIR keeps runs out of the data dir
    r = FiaSim.simulate_landscape(cns, acres, pol; cache_dir = cache)
    (; variant = r.variant, nplots = r.nplots, nexcluded = r.nexcluded, capped = capped,
       master_years = r.master_years, cycles = r.cycles)
end

const MAX_STANDS = 100              # stands per simulation request (delineated-stands vegetation source)

"""
    _simulate_stands_response(stands, plan; engine, raw_plan, params_id, climate)

Vegetation source "delineated stands": the simulation unit is the STAND. `stands` = `[(id, geometry5070)]`. For each stand
the existing AOI pipeline resolves its polygon to TreeMap pixels and FIA plots (`_resolve_plots`: plot acres = pixel
share × pixel area), so a stand's vegetation is the pixel-share-weighted composite of the plots of its pixels. Pan engines
are called once per stand (the site — majority L4 ecoregion and centroid of the stand's pixels — differs per stand); FVSjl
simulates the union of the plots once. Results are aggregated per stand and over all stands; no raster is rendered —
`stands[i].map` holds the per-cycle carbon / ba / volume (area-weighted means) the frontend uses to color the stand polygons.
Each stand is capped to its `MAX_SIM_PLOTS` largest plots.
"""
function _simulate_stands_response(stands, plan; engine::AbstractString = "fvs", raw_plan = nothing,
                                   params_id = nothing, climate = nothing)
    app = APP[]
    cache = get(ENV, "FE_CACHE_DIR", joinpath(app.datadir, "derived", "sim_cache"))
    is_pan = startswith(engine, "pan-")
    info = map(stands) do (id, g)
        acres, tm2cn = _resolve_plots(g)
        cns = sort(collect(keys(acres)); by = cn -> -acres[cn])
        sim_cns = length(cns) > MAX_SIM_PLOTS ? cns[1:MAX_SIM_PLOTS] : cns
        (; id, geom = g, acres, tm2cn, sim_cns)
    end
    warnings = String[]; engine_info = nothing; ppe_out = nothing
    res_for = Vector{Dict{String,FiaSim.PlotResult}}()
    if is_pan
        for st in info
            if isempty(st.sim_cns)
                push!(res_for, Dict{String,FiaSim.PlotResult}()); continue
            end
            r, w, ei, _ = PanEngine.pan_simulate(st.sim_cns, st.acres, raw_plan, engine; params_id = params_id,
                                                 climate = climate, sites = _plot_sites(st.geom, st.tm2cn))
            push!(res_for, r); engine_info = ei
            for x in w                                  # engine-wide notes once, per-stand notes prefixed with the stand
                if startswith(x, "Pan volumes"); x in warnings || push!(warnings, x)
                else push!(warnings, "stand $(st.id): " * x); end
            end
        end
    else
        uacres = Dict{String,Float64}()
        for st in info, (cn, a) in st.acres; uacres[cn] = get(uacres, cn, 0.0) + a; end
        ucns = sort(collect(keys(uacres)); by = cn -> -uacres[cn])
        length(ucns) > MAX_SIM_PLOTS && (ucns = ucns[1:MAX_SIM_PLOTS]; push!(warnings, "FVSjl: the selected stands use more than $(MAX_SIM_PLOTS) distinct plots; the $(MAX_SIM_PLOTS) with the largest area were simulated"))
        plan_for, ppe_out = _ppe_setup(plan, ucns, uacres, cache)
        all = FiaSim.simulate_plots(ucns, plan; cache_dir = cache, plan_for = plan_for)
        for st in info
            push!(res_for, Dict{String,FiaSim.PlotResult}(cn => all[cn] for cn in st.sim_cns if haskey(all, cn)))
        end
    end
    # per stand and overall
    sid = Int[]; sums = Any[]; allunits = Tuple{FiaSim.PlotResult,Float64}[]
    nodata = Int[]
    for (st, res) in zip(info, res_for)
        units = [(v, st.acres[cn]) for (cn, v) in res]
        append!(allunits, units)
        isempty(units) && push!(nodata, st.id)
        ncyc, cycles, _ = _series(units, plan.period)
        mapv = Dict{String,Vector{Float64}}()
        wsum(f, i) = (num = 0.0; den = 0.0;
            for (v, a) in units
                i <= length(v.metrics) || continue
                x = Float64(f(v.metrics[i])); isnan(x) && continue
                num += x * a; den += a
            end; den > 0 ? num / den : NaN)
        for (mname, spec) in SIM_METRICS
            mapv[mname] = [wsum(m -> getfield(m, spec[1]), i) for i in 1:ncyc]
        end
        push!(sums, (; id = st.id, area_ac = round(sum(values(st.acres); init = 0.0); digits = 1), nplots = length(st.acres),
                       nplots_sim = length(res), cycles = cycles, map = mapv))
    end
    ncyc, cycles, dist = _series(allunits, plan.period)
    doms = Dict{String,Tuple{Float64,Float64}}()
    for (mname, _) in SIM_METRICS
        lo = Inf; hi = -Inf
        for st in sums, x in st.map[mname]
            isnan(x) && continue
            lo = min(lo, x); hi = max(hi, x)
        end
        isfinite(lo) && hi > lo && (doms[mname] = (floor(lo), ceil(hi)))
    end
    isempty(nodata) || push!(warnings, "no TreeMap plots simulated for stand(s) " * join(nodata, ", ") *
        " (no TreeMap pixels, e.g. harvested; or plots outside the parameter set's ecoregions)")
    total_ac = sum(st.area_ac for st in sums; init = 0.0)
    (; plan = plan.name, source = "stands", nstands = length(sums), nplots_aoi = sum(st.nplots for st in sums; init = 0),
       nplots_sim = length(allunits), capped = false, acres = round(total_ac; digits = 1),
       sim_acres = round(sum(u[2] for u in allunits; init = 0.0); digits = 1),
       ncycles = ncyc, period = plan.period, cycles = cycles, dist = dist, domains = doms, ppe = ppe_out,
       engine = engine, warnings = warnings, engine_info = engine_info, image = nothing,
       stands = [(; s.id, s.area_ac, s.nplots, s.nplots_sim, s.cycles,
                    map = Dict(k => [isnan(x) ? nothing : round(x; digits = 2) for x in v] for (k, v) in s.map)) for s in sums])
end

"""
Acre-weighted aggregation of simulation units `[(PlotResult, acres)]` into the per-cycle series of a response:
`(ncyc, cycles, dist)` — cycles = BA (ft²/ac, area-weighted), total volume (ft³), removed volume, carbon (t, via HA_PER_ACRE);
dist = per-cycle species × DBH-class tallies (count, BA).
"""
function _series(units, period)
    ncyc = maximum((length(v.metrics) for (v, _) in units); init = 0)
    cycles = map(1:ncyc) do i
        baw = 0.0; w = 0.0; vol = 0.0; rem = 0.0; ct = 0.0; cl = 0.0; el = (i - 1) * period
        for (v, a) in units
            i <= length(v.metrics) || continue
            m = v.metrics[i]
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
    # per-cycle species×DBH distribution (acres-weighted)
    dist = map(1:ncyc) do i
        agg = Dict{Tuple{Int16,Int},NTuple{2,Float64}}()
        syms = Dict{Int16,String}()
        for (v, a) in units
            i <= length(v.metrics) || continue
            g = get(v.dist, v.metrics[i].year, nothing); g === nothing && continue
            merge!(syms, v.spsym)
            for (k, val) in g
                c, b = get(agg, k, (0.0, 0.0))
                agg[k] = (c + val[1] * a, b + val[2] * a)
            end
        end
        [(; spcd = k[1], symbol = get(syms, k[1], ""), dbh_lo = k[2],
            count = round(v[1]; digits = 1), ba = round(v[2]; digits = 1)) for (k, v) in agg]
    end
    ncyc, cycles, dist
end

"""
PPE thinning mode (FVSjl): the cross-stand harvest budget decides WHICH plots are cut in each cycle to meet the resource
flow target; the plan's thin is then applied only to those plots, in those cycles. Every plot therefore carries its own
thin schedule, so the projection runs with a per-plot plan. Returns `(plan_for, ppe_out)`, or `(nothing, nothing)` when
the plan has no PPE thin.
"""
function _ppe_setup(plan, sim_cns, acres, cache)
    ppe_act = findfirst(a -> a.kind == "thin" && a.mode == "ppe", plan.actions)
    ppe_act === nothing && return (nothing, nothing)
    a = plan.actions[ppe_act]
    pol = FiaSim.LandscapePolicy(; common_year = _PPE_YEAR, target_expr = a.target_expr,
            priority_expr = a.priority_expr, credit_expr = a.credit_expr,
            cycles = plan.cycles, period = plan.period, mslabel = "ALL")
    lsc = FiaSim.simulate_landscape(sim_cns, acres, pol; cache_dir = cache)
    # cycle (1-based projection cycle) => set of plot CNs the budget selected
    cutmap = Dict{String,Vector{Int}}()
    for (i, c) in enumerate(lsc.cycles), cn in c.cut
        push!(get!(cutmap, cn, Int[]), i)
    end
    others = [x for x in plan.actions if !(x.kind == "thin" && x.mode == "ppe")]
    plan_for = function (cn)
        cyc = get(cutmap, cn, Int[])
        acts = copy(others)
        for k in cyc
            k >= a.cycle || continue                 # the budget starts at the action's cycle
            push!(acts, FiaSim.PlanAction(; kind = "thin", cycle = k, metric = a.metric,
                    target = a.target, direction = a.direction, species = a.species,
                    dbh_lo = a.dbh_lo, dbh_hi = a.dbh_hi))
        end
        FiaSim.ManagementPlan(; name = plan.name, cycles = plan.cycles, period = plan.period,
                                actions = acts)
    end
    ppe_out = (; variant = lsc.variant, nplots = lsc.nplots, nexcluded = lsc.nexcluded,
                 cycles = [(; elapsed = c.elapsed, target = c.target, resource = c.resource,
                              pct_of_target = c.pct_of_target, ncut = c.ncut,
                              nstands = c.nstands) for c in lsc.cycles])
    plan_for, ppe_out
end

function _simulate_response(geom5070, plan; engine::AbstractString = "fvs", raw_plan = nothing,
                            params_id = nothing, climate = nothing)
    app = APP[]
    acres, tm2cn = _resolve_plots(geom5070)
    total_aoi_acres = sum(values(acres); init = 0.0)
    # cap to the largest-area plots for a bounded synchronous run
    cns = sort(collect(keys(acres)); by = cn -> -acres[cn])
    capped = length(cns) > MAX_SIM_PLOTS
    sim_cns = capped ? cns[1:MAX_SIM_PLOTS] : cns
    cache = get(ENV, "FE_CACHE_DIR", joinpath(app.datadir, "derived", "sim_cache"))   # FE_CACHE_DIR keeps runs out of the data dir

    # PPE thinning mode (FVSjl only): see `_ppe_setup`
    is_pan = startswith(engine, "pan-")
    plan_for, ppe_out = is_pan ? (nothing, nothing) : _ppe_setup(plan, sim_cns, acres, cache)
    warnings = String[]; engine_info = nothing
    res = if is_pan
        # Pan engines take the frontend plan JSON as-is (grow/thin/plant + the LANDIS-style `harvest` kind)
        r, warnings, engine_info, _ = PanEngine.pan_simulate(sim_cns, acres, raw_plan, engine;
                                                             params_id = params_id, climate = climate,
                                                             sites = _plot_sites(geom5070, tm2cn))
        r
    else
        FiaSim.simulate_plots(sim_cns, plan; cache_dir = cache, plan_for = plan_for)
    end

    units = [(v, acres[cn]) for (cn, v) in res]      # iteration order of `res`: keeps the float sums reproducible
    ncyc, cycles, dist = _series(units, plan.period)
    sim_acres = sum(u[2] for u in units; init = 0.0)

    # build the per-pixel projected layer (tm_id => per-cycle metrics) for the time slider
    layer = Dict{UInt32,Vector{FiaSim.CycleMetrics}}()
    for (tm, cn) in tm2cn
        v = get(res, cn, nothing)
        v === nothing || (layer[tm] = v.metrics)
    end
    SIM_LAYER[] = layer
    # auto-scale each map metric's color domain to this simulation's actual value range,
    # so the pixel colors visibly change across cycles (CONUS defaults are far too wide)
    doms = Dict{String,Tuple{Float64,Float64}}()
    for (mname, spec) in SIM_METRICS
        field = spec[1]; lo = Inf; hi = -Inf
        for v in values(layer), m in v
            x = Float64(getfield(m, field)); isnan(x) && continue
            lo = min(lo, x); hi = max(hi, x)
        end
        isfinite(lo) && hi > lo && (doms[mname] = (floor(lo), ceil(hi)))
    end
    SIM_DOMAINS[] = doms
    geom3857 = Aoi.to3857(geom5070)
    SIM_AOI[] = geom3857
    # pre-render + cache every (metric, cycle) AOI image; the slider just swaps them
    imgmeta = ncyc > 0 && !isempty(layer) ? _render_sim_images(geom3857, layer, doms, ncyc) : nothing

    (; plan = plan.name, nplots_aoi = length(acres), nplots_sim = length(res),
       capped = capped, acres = round(total_aoi_acres; digits = 1),
       sim_acres = round(sim_acres; digits = 1),
       ncycles = ncyc, period = plan.period, cycles = cycles, dist = dist,
       domains = doms,        # auto-scaled color domains per map metric (for the legend)
       ppe = ppe_out,         # PPE thinning mode: the per-cycle budget flow + how many plots it cut
       engine = engine, warnings = warnings, engine_info = engine_info,
       image = imgmeta)       # cached AOI raster overlay: token, corners (lon/lat), size, metrics
end

"Serve a scenario's pre-rendered cached AOI raster for (token, metric, cycle). The token
selects the SCENARIO (each simulation caches its own set); browser-cacheable immutably."
function _simimage_handler(token::Int, metric::String, cycle::Int)
    imgs = get(SIM_IMAGES[], token, nothing)
    imgs === nothing && return HTTP.Response(404)
    png = get(imgs, (metric, cycle), nothing)
    png === nothing && return HTTP.Response(404)
    HTTP.Response(200, ["Content-Type" => "image/png",
                        "Cache-Control" => "public, max-age=31536000, immutable"];
                  body = png)
end


# JSON3 objects → plain Dict/Vector (so they can be re-serialized to the Pan worker)
_plain(x) = x === nothing ? nothing : x isa JSON3.Object ? Dict{String,Any}(String(k) => _plain(v) for (k, v) in pairs(x)) :
            x isa JSON3.Array ? Any[_plain(v) for v in x] : x
# FVS plan parser input without the Pan-only `harvest` actions
_strip_harvest(p) = (d = _plain(p); d["actions"] = [a for a in get(d, "actions", Any[]) if get(a, "kind", "") != "harvest"]; JSON3.read(JSON3.write(d)))

const SAMPLES = normpath(joinpath(@__DIR__, "..", "samples"))
function _samples_index()
    f = joinpath(SAMPLES, "index.json")
    isfile(f) ? JSON3.read(read(f, String)) : (; params = [], aois = [], bundles = [], climates = [])
end
function _sample_file(name::AbstractString)
    p = normpath(joinpath(SAMPLES, name))
    (startswith(p, SAMPLES) && isfile(p)) || return HTTP.Response(404, JSON3.write((; error = "no sample $name")))
    ct = endswith(p, ".yaml") || endswith(p, ".yml") ? "text/yaml" : "application/json"
    HTTP.Response(200, ["Content-Type" => ct], read(p))
end

"AOI-weighted climate series for the climate modal (Pan PnET)."
function _climate_preview_response(geom5070, climate, start_year, years)
    acres, tm2cn = _resolve_plots(geom5070)
    cns = sort(collect(keys(acres)); by = cn -> -acres[cn])
    length(cns) > MAX_SIM_PLOTS && (cns = cns[1:MAX_SIM_PLOTS])
    isempty(cns) && return (; error = "no plots resolved in the AOI")
    PanEngine.pan_climate_preview(cns, acres, climate; start_year = start_year, years = years, sites = _plot_sites(geom5070, tm2cn))
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

    # cached AOI raster overlay: /simimage/{token}/{metric}/{cycle}.png . The token
    # (bumped per simulation/scenario) selects which scenario's cached raster set to serve.
    @get "/simimage/{token}/{metric}/{cycle}" function (req, token::String, metric::String, cycle::String)
        try
            return _simimage_handler(parse(Int, token), metric, parse(Int, first(split(cycle, '.'))))
        catch e
            @warn "simimage error" exception = e
            return HTTP.Response(404)
        end
    end

    @get "/api/attrs" function (req)
        [(; name = s.name, label = s.label, units = s.units, lo = s.lo, hi = s.hi,
            cmap = _cmap_name(s.cmap))
         for s in sort(collect(values(ATTRS)); by = s -> s.name)]
    end

    @get "/api/simprogress" function (req)
        FiaSim.progress()
    end

    @get "/api/meta" function (req)
        (; default_attr = DEFAULT_ATTR, bounds = CONUS_BOUNDS_4326,
           nplots_total = length(app.st.by_tm),
           sources = SOURCES, engines = ENGINES, pan = PanEngine.pan_status())
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
            engine = String(get(body, :engine, "fvs"))
            engine == "fvsjl" && (engine = "fvs")
            plan = engine == "fvs" ? _parse_plan(body.plan) :
                   _parse_plan(_strip_harvest(body.plan))    # Pan: parsed only for name/cycles/period
            pid = get(body, :params_id, nothing)
            if haskey(body, :stands) && body.stands !== nothing
                length(body.stands) in 1:MAX_STANDS || error("select 1–$(MAX_STANDS) stands (got $(length(body.stands)))")
                stands = [(Int(st.id), Aoi.geom_from_geojson(JSON3.write(st.geometry))) for st in body.stands]
                return _simulate_stands_response(stands, plan; engine = engine, raw_plan = _plain(body.plan),
                                                 params_id = pid === nothing ? nothing : String(pid),
                                                 climate = _plain(get(body, :climate, nothing)))
            end
            geom = Aoi.geom_from_geojson(JSON3.write(body.geometry))
            return _simulate_response(geom, plan; engine = engine, raw_plan = _plain(body.plan),
                                      params_id = pid === nothing ? nothing : String(pid),
                                      climate = _plain(get(body, :climate, nothing)))
        catch e
            @warn "simulate error" exception = (e, catch_backtrace())
            return HTTP.Response(400, JSON3.write((; error = string(e))))
        end
    end

    # PPE cross-stand harvest budget (landscape flow) over the AOI's plots
    @post "/api/landscape" function (req)
        try
            body = JSON3.read(String(req.body))
            geom = Aoi.geom_from_geojson(JSON3.write(body.geometry))
            pol = _parse_policy(body.policy)
            return _landscape_response(geom, pol)
        catch e
            @warn "landscape error" exception = (e, catch_backtrace())
            return HTTP.Response(400, JSON3.write((; error = string(e))))
        end
    end

    # ---- Pan engine (Biomass Succession / PnET / HBAD) ----
    @get "/api/pan/samples" function (req)
        _samples_index()
    end
    @get "/api/pan/sample" function (req)
        q = HTTP.queryparams(HTTP.URI(req.target))
        _sample_file(get(q, "file", ""))
    end
    @get "/api/pan/default_params" function (req)
        q = HTTP.queryparams(HTTP.URI(req.target))
        try
            return PanEngine.pan_default_params(get(q, "plugin", "biomass"))
        catch e
            return HTTP.Response(400, JSON3.write((; error = sprint(showerror, e))))
        end
    end
    @post "/api/pan/params" function (req)
        try
            return PanEngine.pan_load_params(String(req.body))
        catch e
            return HTTP.Response(400, JSON3.write((; error = sprint(showerror, e))))
        end
    end
    @post "/api/pan/climate_preview" function (req)
        try
            body = JSON3.read(String(req.body))
            geom = Aoi.geom_from_geojson(JSON3.write(body.geometry))
            return _climate_preview_response(geom, _plain(get(body, :climate, nothing)),
                                             Int(get(body, :start_year, 2025)), Int(get(body, :years, 100)))
        catch e
            @warn "climate preview error" exception = (e, catch_backtrace())
            return HTTP.Response(400, JSON3.write((; error = sprint(showerror, e))))
        end
    end

    # optional layer: FIA plots -> tree trajectories (src/FiaTrajectories.jl; read-only DuckDB)
    FiaTrajectories.register!()

    # static frontend
    staticfiles(PUBLIC, "/")

    PanEngine.autostart!()   # warm the Pan worker in the background (PAN_AUTOSTART=0 to disable)
    @info "Forest Growth Explorer serving — press Ctrl-C to stop" host port public = PUBLIC

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
