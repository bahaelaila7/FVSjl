"""
    Color

Colormaps and tile colorization for the TreeMap display layer.

`colorize` is generic over a `valfn :: tm_id -> Float64` (NaN = no value), so
the identical path serves a static VAT attribute now and a per-cycle projected
metric under the time slider later. TM_ID is thematic, so display always means
"look up an attribute for this plot, then map that scalar to a color".
"""
module Color

using ColorTypes, FixedPointNumbers
using ..Raster: NODATA

export Colormap, colorize, VIRIDIS, GREENS, MAGMA, ramp

struct Colormap
    stops::Vector{NTuple{3,Float64}}   # RGB in 0..1, evenly spaced
end

"Linear-interpolate a colormap at t∈[0,1]."
function ramp(cm::Colormap, t::Float64)
    t = isnan(t) ? 0.0 : clamp(t, 0.0, 1.0)
    n = length(cm.stops)
    n == 1 && return cm.stops[1]
    f = t * (n - 1)
    i = floor(Int, f)
    i >= n - 1 && return cm.stops[n]
    a = cm.stops[i+1]; b = cm.stops[i+2]
    w = f - i
    (a[1] + w*(b[1]-a[1]), a[2] + w*(b[2]-a[2]), a[3] + w*(b[3]-a[3]))
end

# viridis (8 control points)
const VIRIDIS = Colormap([
    (0.267,0.005,0.329),(0.283,0.141,0.458),(0.254,0.265,0.530),
    (0.207,0.372,0.553),(0.164,0.471,0.558),(0.128,0.567,0.551),
    (0.478,0.821,0.318),(0.993,0.906,0.144)])

# sequential greens (forest biomass/carbon feel)
const GREENS = Colormap([
    (0.969,0.988,0.961),(0.780,0.914,0.753),(0.455,0.769,0.462),
    (0.137,0.545,0.271),(0.000,0.267,0.106)])

const MAGMA = Colormap([
    (0.001,0.000,0.014),(0.279,0.059,0.412),(0.639,0.190,0.408),
    (0.945,0.377,0.365),(0.996,0.702,0.384),(0.987,0.991,0.749)])

@inline _u8(x::Float64) = N0f8(clamp(x, 0.0, 1.0))

"""
    colorize(tmids, valfn, cmap, lo, hi) -> Matrix{RGBA{N0f8}}

`tmids` is `[row,col]` from `Raster.warp_tile`. `valfn(tm_id)` returns the
scalar to color (NaN → fully transparent). Values are normalized on `[lo,hi]`.
NODATA and unknown plots are transparent.
"""
function colorize(tmids::AbstractMatrix{UInt32}, valfn, cmap::Colormap,
                  lo::Real, hi::Real)
    rows, cols = size(tmids)
    img = Matrix{RGBA{N0f8}}(undef, rows, cols)
    invspan = hi > lo ? 1.0 / (hi - lo) : 0.0
    transparent = RGBA{N0f8}(0, 0, 0, 0)
    @inbounds for j in 1:cols, i in 1:rows
        tm = tmids[i, j]
        if tm == NODATA
            img[i, j] = transparent; continue
        end
        v = valfn(tm)
        if isnan(v)
            img[i, j] = transparent; continue
        end
        t = (Float64(v) - lo) * invspan
        r, g, b = ramp(cmap, t)
        img[i, j] = RGBA{N0f8}(_u8(r), _u8(g), _u8(b), N0f8(1))
    end
    img
end

end # module
