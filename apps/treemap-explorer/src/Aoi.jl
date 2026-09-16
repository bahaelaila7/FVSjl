"""
    Aoi

Parse an area-of-interest polygon — either a GeoJSON geometry drawn in the
browser (assumed EPSG:4326) or an uploaded GeoJSON/shapefile — and reproject
it to the raster CRS (EPSG:5070) as an ArchGDAL geometry ready for `aoi_tally`.
"""
module Aoi

using ArchGDAL
const AG = ArchGDAL

export geom_from_geojson, geom_from_upload, bbox4326, geojson4326

const EPSG5070_WKT = Ref{String}()

function _wkt5070()
    if !isassigned(EPSG5070_WKT)
        sr = AG.importEPSG(5070)
        EPSG5070_WKT[] = AG.toWKT(sr)
    end
    EPSG5070_WKT[]
end

"Reproject a geometry from `src_epsg` (default 4326) to EPSG:5070; returns a clone."
function to5070(geom, src_epsg::Integer = 4326)
    src = AG.importEPSG(src_epsg)
    dst = AG.importEPSG(5070)
    # honor traditional lon/lat axis order for 4326 input
    AG.GDAL.osrsetaxismappingstrategy(src.ptr, AG.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    g = AG.clone(geom)
    AG.createcoordtrans(src, dst) do ct
        AG.transform!(g, ct)
    end
    g
end

"Bounding box (west,south,east,north) in EPSG:4326 of a 5070 geometry."
function bbox4326(geom5070)
    src = AG.importEPSG(5070)
    dst = AG.importEPSG(4326)
    AG.GDAL.osrsetaxismappingstrategy(dst.ptr, AG.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    g = AG.clone(geom5070)
    AG.createcoordtrans(src, dst) do ct
        AG.transform!(g, ct)
    end
    e = AG.envelope(g)
    (e.MinX, e.MinY, e.MaxX, e.MaxY)
end

"GeoJSON *geometry* string (EPSG:4326 lon/lat) of a 5070 geometry — for the frontend to draw."
function geojson4326(geom5070)
    src = AG.importEPSG(5070)
    dst = AG.importEPSG(4326)
    AG.GDAL.osrsetaxismappingstrategy(dst.ptr, AG.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    g = AG.clone(geom5070)
    AG.createcoordtrans(src, dst) do ct
        AG.transform!(g, ct)
    end
    AG.toJSON(g)
end

"Build a 5070 geometry from a raw GeoJSON *geometry* string (drawn polygon, 4326)."
function geom_from_geojson(geojson::AbstractString)
    g = AG.fromJSON(String(geojson))
    to5070(g, 4326)
end

"""
    geom_from_upload(bytes, filename) -> geom5070

Open an uploaded AOI — a **GeoJSON** file, or a **.zip bundle of a shapefile**
(.shp/.shx/.dbf/.prj) which GDAL reads in place via `/vsizip/`. A bare `.shp`
alone cannot be read (its `.shx`/`.prj` are missing), so upload the zip.
Returns the union of all feature geometries reprojected to EPSG:5070.
"""
function geom_from_upload(bytes::Vector{UInt8}, filename::AbstractString)
    base = basename(filename)
    vpath = "/vsimem/aoi_" * string(hash(bytes); base=16) * "_" * base
    AG.GDAL.vsifilefrommembuffer(vpath, bytes, length(bytes), false)
    # a zipped shapefile bundle is opened through the /vsizip/ virtual filesystem
    open_path = endswith(lowercase(base), ".zip") ? "/vsizip/" * vpath : vpath
    try
        ds = AG.read(open_path)
        lyr = AG.getlayer(ds, 0)
        # union all feature geometries
        acc = nothing
        for f in lyr
            g = AG.getgeom(f)
            g === nothing && continue
            acc = acc === nothing ? AG.clone(g) : AG.union(acc, g)
        end
        acc === nothing && error("no geometry found in AOI upload")
        # Reproject to 5070 using the layer's OWN spatial ref (handles any CRS —
        # a shapefile .prj, or a GeoJSON read as OGC:CRS84 which has no EPSG code).
        # Falls back to WGS84 lon/lat when the layer carries no CRS (bare GeoJSON).
        srcsr = AG.getspatialref(lyr)
        srcsr === nothing && (srcsr = AG.importEPSG(4326))
        AG.GDAL.osrsetaxismappingstrategy(srcsr.ptr, AG.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
        dst = AG.importEPSG(5070)
        g5070 = AG.clone(acc)
        AG.createcoordtrans(srcsr, dst) do ct
            AG.transform!(g5070, ct)
        end
        return g5070
    finally
        AG.GDAL.vsiunlink(vpath)
    end
end

end # module
