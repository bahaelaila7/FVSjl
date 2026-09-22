using ForestExplorer
const TME = ForestExplorer

app = TME.init_app()
println("VAT plots: ", length(app.st.by_tm))
println("tree-table plots: ", length(app.tt.tm_id))

# --- tile warp + colorize + PNG ---
# a low-zoom tile covering the US interior west; z=5 tile x=5 y=11 ~ N Rockies
z, x, y = 5, 5, 11
tm = TME.Raster.warp_tile(app.rasterctx, z, x, y; size = 256)
nvalid = count(!=(TME.Raster.NODATA), tm)
println("tile z$z/$x/$y valid px: ", nvalid, " / ", length(tm))
spec = TME.ATTRS["carbon_l"]
img = TME.Color.colorize(tm, TME.attr_valfn(app, "carbon_l"), spec.cmap, spec.lo, spec.hi)
import PNGFiles
import ColorTypes: alpha
outpng = joinpath(mktempdir(), "test_tile.png")
PNGFiles.save(outpng, img)
println("wrote ", outpng, " (", count(p -> alpha(p) > 0, img), " colored px)")

# --- AOI clip + aggregate: ~6 km box in the Idaho panhandle (dense IE forest) ---
using ArchGDAL
const AG = ArchGDAL
lon, lat, d = -115.7, 47.6, 0.03
gj = """{"type":"Polygon","coordinates":[[
  [$(lon-d),$(lat-d)],[$(lon+d),$(lat-d)],[$(lon+d),$(lat+d)],[$(lon-d),$(lat+d)],[$(lon-d),$(lat-d)]]]}"""
geom = TME.Aoi.geom_from_geojson(gj)
@time counts = TME.Raster.aoi_tally(app.rasterctx, geom)
println("AOI unique plots: ", length(counts), "  total forested pixels: ", sum(values(counts)))
agg = TME.Aggregate.aggregate_aoi(counts, app.tt, app.st)
println("acres=", round(agg.acres; digits=1),
        "  BA_total=", round(agg.ba_total; digits=0), " ft²",
        "  live trees=", round(agg.tpa_expanded; digits=0),
        "  carbon_L=", round(agg.carbon_l; digits=0), " t",
        "  biomass_L=", round(agg.drybio_l; digits=0), " t",
        "  vol_L=", round(agg.volcfnet_l; digits=0), " ft³")
println("top species×dbh cells:")
for c in first(agg.cells, min(8, length(agg.cells)))
    println("  ", rpad(get(agg.species_symbol, c.spcd, "?"), 6),
            " dbh ", lpad(c.dbh_lo, 2), "-", c.dbh_lo+2, "\"",
            "  n=", round(c.count; digits=1), "  BA=", round(c.ba; digits=1), " ft²")
end
println("SMOKE_OK")
