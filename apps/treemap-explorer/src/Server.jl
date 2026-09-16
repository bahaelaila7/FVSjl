# Included into TreeMapExplorer. HTTP layer (Oxygen): dynamic tiles + AOI API.

using Oxygen, HTTP, JSON3
import PNGFiles

const APP = Ref{App}()
const PUBLIC = normpath(joinpath(@__DIR__, "..", "public"))

# raster geographic bounds (EPSG:4326) for map fit — CONUS Albers extent corners
const CONUS_BOUNDS_4326 = (west = -127.9, south = 22.7, east = -65.3, north = 51.6)

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

    @get "/api/attrs" function (req)
        [(; name = s.name, label = s.label, units = s.units, lo = s.lo, hi = s.hi,
            cmap = _cmap_name(s.cmap))
         for s in sort(collect(values(ATTRS)); by = s -> s.name)]
    end

    @get "/api/meta" function (req)
        (; default_attr = DEFAULT_ATTR, bounds = CONUS_BOUNDS_4326,
           nplots_total = length(app.st.by_tm))
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
