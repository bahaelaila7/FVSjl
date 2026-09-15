"""
    TreeMapExplorer

Localhost web app: display the USFS TreeMap 2022 CONUS raster, resolve an AOI
polygon to its FIA plots, and chart species / DBH-class distributions of the
forest — with FVSjl growth projection behind a single `simulate` boundary.

Data (raster, VAT, tree table, offset index) lives OUTSIDE the repo, under a
data directory resolved from the `init_app` argument, else the `TREEMAP_DATA`
env var, else a default derived from this module's location (see `default_datadir`).
"""
module TreeMapExplorer

include("Raster.jl")
include("TreeList.jl")
include("VAT.jl")
include("Color.jl")
include("Aoi.jl")
include("Aggregate.jl")

using .Raster, .TreeList, .VAT, .Color, .Aoi, .Aggregate

export App, init_app, AttrSpec, ATTRS

# --- display attribute registry ---------------------------------------------
struct AttrSpec
    name::String
    label::String
    units::String
    cmap::Color.Colormap
    lo::Float64
    hi::Float64
    field::Symbol          # StandRecord field to read
end

const ATTRS = Dict{String,AttrSpec}(
    "carbon_l"   => AttrSpec("carbon_l","Live aboveground carbon","tons/ac",Color.GREENS, 0, 60, :carbon_l),
    "drybio_l"   => AttrSpec("drybio_l","Aboveground live biomass","tons/ac",Color.GREENS, 0,120, :drybio_l),
    "volcfnet_l" => AttrSpec("volcfnet_l","Live net volume","ft³/ac",       Color.VIRIDIS,0,8000,:volcfnet_l),
    "balive"     => AttrSpec("balive","Live basal area","ft²/ac",           Color.VIRIDIS,0,250, :balive),
    "qmd"        => AttrSpec("qmd","Quadratic mean diameter","in",          Color.MAGMA,  0, 24, :qmd),
    "tpa_live"   => AttrSpec("tpa_live","Live trees per acre","/ac",        Color.VIRIDIS,0,1500,:tpa_live),
    "standht"    => AttrSpec("standht","Dominant stand height","ft",        Color.VIRIDIS,0,150, :standht),
    "canopypct"  => AttrSpec("canopypct","Live canopy cover","%",           Color.GREENS, 0,100, :canopypct),
)
const DEFAULT_ATTR = "carbon_l"

# --- application state -------------------------------------------------------
struct App
    datadir::String
    rasterctx::Raster.RasterCtx
    st::VAT.StandTable
    tt::TreeList.TreeTable
end

"""
    init_app(datadir=default_datadir()) -> App

Load the VAT (65k plots, ~instant) and bind the tree-table offset index and the
raster. The 4.85 GB raster is memory-mapped by GDAL, not read up front.

The data directory is resolved (in order): the `datadir` argument, the
`TREEMAP_DATA` environment variable, else a default derived from this module's
own location (`<repo-parent>/treemap`) — no absolute path is hard-coded.
"""
function default_datadir()
    haskey(ENV, "TREEMAP_DATA") && return ENV["TREEMAP_DATA"]
    # this file: <root>/apps/treemap-explorer/src/TreeMapExplorer.jl
    # data sits beside the repo root: <root>/../treemap  ==  <repo-parent>/treemap
    normpath(joinpath(@__DIR__, "..", "..", "..", "..", "treemap"))
end

function init_app(datadir::AbstractString = default_datadir())
    dd = joinpath(datadir, "Data")
    tif = joinpath(dd, "TreeMap2022_CONUS.tif")
    vat = joinpath(dd, "TreeMap2022_CONUS.tif.vat.dbf")
    csv = joinpath(dd, "TreeMap2022_CONUS_Tree_Table.csv")
    idx = joinpath(datadir, "derived", "treetable_index.bin")
    @info "loading VAT" vat
    st = VAT.load_vat(vat)
    @info "binding tree table" csv idx
    tt = TreeList.TreeTable(csv, idx)
    App(String(datadir), Raster.RasterCtx(tif), st, tt)
end

"Return a `tm_id -> Float64` value function for a display attribute (NaN = none)."
function attr_valfn(app::App, attr::AbstractString)
    spec = get(ATTRS, attr, ATTRS[DEFAULT_ATTR])
    f = spec.field
    st = app.st
    return tm -> begin
        rec = VAT.stand(st, tm)
        rec === nothing ? NaN : Float64(getfield(rec, f))
    end
end

include("Server.jl")   # defines start_server! into this module

end # module
