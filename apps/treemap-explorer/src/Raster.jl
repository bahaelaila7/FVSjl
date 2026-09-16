"""
    Raster

TreeMap raster access via ArchGDAL/GDAL.

* Dynamic XYZ tiler: for a web-mercator tile (z,x,y) GDAL warps the EPSG:5070
  TM_ID raster (using the external overviews) into a 256×256 tile of TM_IDs,
  which the caller colors by a VAT attribute. GDAL picks the overview level.
* AOI clip: rasterize an AOI polygon (already in EPSG:5070) over the covering
  window and tally TM_ID → forested-pixel count, so an AOI resolves to its set
  of unique plots and their pixel weights.

The raster CRS is NAD83 / CONUS Albers (EPSG:5070); pixel value == TM_ID;
nodata == 0xFFFFFFFF.
"""
module Raster

using ArchGDAL
const AG = ArchGDAL

export RasterCtx, warp_tile, warp_extent, extent_mask, aoi_tally, aoi_tile_mask, NODATA

const NODATA = typemax(UInt32)                    # 4294967295
const WEBMERC_HALF = 20037508.342789244           # half-extent of EPSG:3857

struct RasterCtx
    path::String
end

# --- web-mercator XYZ tile bounds in EPSG:3857 meters ---
function tile_bounds_3857(z::Integer, x::Integer, y::Integer)
    res = (2 * WEBMERC_HALF) / (2.0^z)
    xmin = -WEBMERC_HALF + x * res
    xmax = xmin + res
    ymax = WEBMERC_HALF - y * res
    ymin = ymax - res
    (xmin, ymin, xmax, ymax)
end

"""
    warp_tile(ctx, z, x, y; size=256) -> Matrix{UInt32}

Warp the TM_ID raster into one XYZ tile (nearest resampling, EPSG:3857).
Returns a `size×size` matrix of TM_IDs in image order `[row, col]`,
`NODATA` outside the raster / non-data.
"""
function warp_tile(ctx::RasterCtx, z::Integer, x::Integer, y::Integer; size::Int=256)
    xmin, ymin, xmax, ymax = tile_bounds_3857(z, x, y)
    warp_extent(ctx, xmin, ymin, xmax, ymax, size, size)
end

"""
    warp_extent(ctx, xmin, ymin, xmax, ymax, w, h) -> Matrix{UInt32}   [row, col]

Warp the TM_ID raster over an arbitrary EPSG:3857 extent into a `w×h` matrix
(nearest resampling). `NODATA` outside the raster. This is the tiler generalized
to any window — used to render a whole AOI into one cached image.
"""
function warp_extent(ctx::RasterCtx, xmin, ymin, xmax, ymax, w::Integer, h::Integer)
    opts = [
        "-of", "MEM",
        "-t_srs", "EPSG:3857",
        "-te", string(xmin), string(ymin), string(xmax), string(ymax),
        "-ts", string(w), string(h),
        "-r", "near",
        "-dstnodata", string(NODATA),
        "-wo", "NUM_THREADS=ALL_CPUS",
    ]
    AG.read(ctx.path) do src
        AG.gdalwarp([src], opts) do dst
            band = AG.getband(dst, 1)
            arr = AG.read(band)             # [x(col), y(row)]
            permutedims(arr, (2, 1))        # -> [row, col]
        end
    end
end

"""
    extent_mask(geom3857, xmin, ymin, xmax, ymax, w, h) -> Matrix{UInt8}   [row, col]

Rasterize the AOI polygon (EPSG:3857) over the same `w×h` 3857 window as
`warp_extent` — 1 inside the polygon, 0 outside — to clip the rendered image.
"""
function extent_mask(geom3857, xmin, ymin, xmax, ymax, w::Integer, h::Integer)
    srs = AG.toWKT(AG.importEPSG(3857))
    permutedims(_rasterize_mask(geom3857, w, h, xmin, ymin, xmax, ymax, srs), (2, 1))
end

# --- geotransform helper on the source (EPSG:5070) ---
"world(5070) -> fractional pixel (col,row)"
function world_to_pixel(gt, X, Y)
    col = (X - gt[1]) / gt[2]
    row = (Y - gt[4]) / gt[6]
    (col, row)
end

"""
    aoi_tally(ctx, geom5070) -> Dict{UInt32,Int}

Rasterize the AOI polygon (EPSG:5070 ArchGDAL geometry) over its covering
full-resolution window and tally TM_ID → forested-pixel count for pixels
inside the polygon. Excludes NODATA.
"""
function aoi_tally(ctx::RasterCtx, geom5070; max_pixels::Int = 200_000_000)
    AG.read(ctx.path) do src
        gt = AG.getgeotransform(src)
        W = AG.width(src); H = AG.height(src)
        env = AG.envelope(geom5070)                       # MinX,MaxX,MinY,MaxY
        c0, r0 = world_to_pixel(gt, env.MinX, env.MaxY)   # top-left
        c1, r1 = world_to_pixel(gt, env.MaxX, env.MinY)   # bottom-right
        col0 = clamp(floor(Int, min(c0, c1)), 0, W - 1)
        col1 = clamp(ceil(Int,  max(c0, c1)), 0, W - 1)
        row0 = clamp(floor(Int, min(r0, r1)), 0, H - 1)
        row1 = clamp(ceil(Int,  max(r0, r1)), 0, H - 1)
        wcols = col1 - col0 + 1
        wrows = row1 - row0 + 1
        counts = Dict{UInt32,Int}()
        (wcols <= 0 || wrows <= 0) && return counts
        wcols * wrows > max_pixels &&
            error("AOI window $(wcols)×$(wrows) exceeds $(max_pixels) px; refine AOI")

        # window bounds in 5070 (edges), aligned to source pixels
        wxmin = gt[1] + col0 * gt[2]
        wxmax = gt[1] + (col1 + 1) * gt[2]
        wymax = gt[4] + row0 * gt[6]
        wymin = gt[4] + (row1 + 1) * gt[6]

        mask = _rasterize_mask(geom5070, wcols, wrows, wxmin, wymin, wxmax, wymax,
                               AG.getproj(src))            # [x,y] UInt8

        band = AG.getband(src, 1)
        tm = AG.read(band, col0 + 1, row0 + 1, wcols, wrows)  # 1-based, [x,y]
        @inbounds for cix in 1:wcols, rix in 1:wrows
            mask[cix, rix] == 0x00 && continue
            v = tm[cix, rix]
            v == NODATA && continue
            counts[v] = get(counts, v, 0) + 1
        end
        return counts
    end
end

"""
    aoi_tile_mask(geom3857, z, x, y; size=256) -> Matrix{UInt8}   [row, col]

Rasterize the AOI polygon (already in EPSG:3857) onto one web-mercator XYZ tile
grid — 1 inside the polygon, 0 outside — aligned pixel-for-pixel with
`warp_tile`. Used to clip a per-plot layer to the AOI (TreeMap imputes each plot
to many pixels nationwide, so a tm_id lookup alone would paint far outside it).
"""
function aoi_tile_mask(geom3857, z::Integer, x::Integer, y::Integer; size::Int=256)
    xmin, ymin, xmax, ymax = tile_bounds_3857(z, x, y)
    srs = AG.toWKT(AG.importEPSG(3857))
    m = _rasterize_mask(geom3857, size, size, xmin, ymin, xmax, ymax, srs)  # [x,y]
    permutedims(m, (2, 1))                                                  # -> [row,col]
end

"Burn geom5070 into a wcols×wrows UInt8 mask ([x,y]); 1 inside, 0 outside."
function _rasterize_mask(geom5070, wcols, wrows, xmin, ymin, xmax, ymax, srs_wkt)
    # in-memory OGR datasource holding the single polygon
    vds = AG.create(""; driver = AG.getdriver("Memory"))
    sref = AG.importWKT(srs_wkt)
    lyr = AG.createlayer(; dataset = vds, geom = AG.wkbPolygon, spatialref = sref)
    AG.createfeature(lyr) do f
        AG.setgeom!(f, AG.clone(geom5070))
    end
    opts = [
        "-of", "MEM", "-ot", "Byte",
        "-burn", "1", "-init", "0",
        "-te", string(xmin), string(ymin), string(xmax), string(ymax),
        "-ts", string(wcols), string(wrows),
    ]
    AG.gdalrasterize(vds, opts) do rds
        AG.read(AG.getband(rds, 1))     # [x,y]
    end
end

end # module
