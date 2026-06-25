# =============================================================================
# site_index.jl (northeast) — NE SITSET site-index fan-out (ne/sitset.f)
#
# NE converts the input site index (for the site species ISISP) to a site index
# for EVERY species via a 28×28 SICOEF conversion matrix indexed by each species'
# site group IPOINT (the `site_group` column). The flag values 0/1 in SICOEF select
# "same", "linear-inverse", or "intercept+slope" conversion (sitset.f:194-200).
# Structurally different from SN's master-group SITSET, so it's its own method.
# =============================================================================

# SICOEF[28,28] loaded once from data/northeast/site_coef.csv (row_i, col1..col28).
const _NE_SICOEF = let
    path = joinpath(NE_DATADIR, "site_coef.csv")
    m = zeros(Float32, 28, 28)
    for (k, line) in enumerate(eachline(path))
        k == 1 && continue
        f = split(strip(line), ',')
        i = parse(Int, f[1])
        for j in 1:28; m[i, j] = parse(Float32, f[j+1]); end
    end
    m
end

"""
    ne_site_index_setup!(s)

SITSET (ne/sitset.f): fill `plot.sp_site_index[sp]` for every species by converting the
site species' site index through the SICOEF matrix (keyed by each species' `site_group`).
Default site species = 27 (sugar maple), default site index = 56.
"""
function ne_site_index_setup!(s::StandState)
    p = s.plot; v = s.variant
    ipoint = s.coef.species[:site_group]
    sea = p.sp_site_index
    isisp = Int(p.site_species)
    isisp <= 0 && (isisp = 27; p.site_species = Int32(27))
    sea[isisp] <= 0f0 && (sea[isisp] = 56f0)
    base = sea[isisp]
    jsite = Int(ipoint[isisp])
    @inbounds for ispc in 1:nspecies(v)
        sea[ispc] > 0.0001f0 && continue
        jspc = Int(ipoint[ispc])
        c = _NE_SICOEF[jsite, jspc]
        if c == 0f0
            sea[ispc] = base
        elseif c == 1f0
            sea[ispc] = (base - _NE_SICOEF[jspc, jsite]) / 1.104f0
        else
            sea[ispc] = c + 1.104f0 * base
        end
    end
    return s
end

# NE forest location codes (ne/forkod.f JFOR) and the default IFOR (grinit.f:189 sets IFOR=2).
const _NE_JFOR = Int32[914, 922, 919, 920, 921, 911, 930]
const _NE_DEFAULT_IFOR = 2
# Per-IFOR default (lat, long, elev) — forkod.f:187-208. CASE(6)/(7) are remapped to 1/4 first.
const _NE_FOR_DEFAULTS = Dict(
    1 => (39.33f0, 82.10f0,  9f0),
    2 => (43.53f0, 71.47f0, 20f0),
    3 => (41.84f0, 79.15f0, 17f0),
    4 => (43.61f0, 72.97f0, 19f0),
    5 => (38.93f0, 79.85f0, 30f0),
)

"""
    ne_forkod_defaults!(s)

FORKOD phase 3 (ne/forkod.f): map the user forest code to an NE forest index (`IFOR`) via
`JFOR`, applying the Wayne-Hoosier (6→1) and Finger-Lakes (7→4) remaps, then set the default
TLAT/TLONG/ELEV for that forest when the stand has none. An unrecognized/absent code keeps the
variant default `IFOR=2` (White Mtn — 43.53/71.47/20), matching grinit.f's `IFOR=2` init. The
Hopkins index (used by the crown-width CCF) needs a real lat/long, so this must run even when a
stand carries no STDINFO/LOCATION keyword (e.g. net01).
"""
function ne_forkod_defaults!(s::StandState)
    p = s.plot
    code = Int32(p.user_forest_code)
    ifor = _NE_DEFAULT_IFOR
    idx = findfirst(==(code), _NE_JFOR)
    idx !== nothing && (ifor = idx)
    ifor == 6 && (ifor = 1)                 # Wayne-Hoosier (911) → Wayne (914)
    ifor == 7 && (ifor = 4)                 # Finger Lakes (930) → Green Mtn (920)
    lat0, long0, elev0 = _NE_FOR_DEFAULTS[ifor]
    p.latitude  == 0f0 && (p.latitude  = lat0)
    p.longitude == 0f0 && (p.longitude = long0)
    p.elevation == 0f0 && (p.elevation = elev0)
    return s
end

function site_setup!(s::StandState, ::Northeast)
    ne_forkod_defaults!(s)
    ne_site_index_setup!(s)
end
