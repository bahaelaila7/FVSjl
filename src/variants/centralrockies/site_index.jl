# =============================================================================
# site_index.jl (centralrockies) — CR per-species site-index bounds (cr/blkdat.f
# SITELO/SITEHI). Used by the site-index between-species conversion (cr/sitset.f)
# and small-tree/regen relative site index (cr/regent.f:213 RELSI). Full sitset
# ISISP mapping + IMODTY model-type defaults + habitat-type groups land in the
# site/habitat chunk; this provides the bound data + the relative-SI helper.
# =============================================================================

"CR relative site index for species `sp` (cr/regent.f:209-213): clamp SI to SITELO+0.5, then normalize."
@inline function cr_relative_si(sd, sp::Integer, si::Real)::Float32
    lo = sd[:site_lo][sp]; hi = sd[:site_hi][sp]
    s = Float32(si)
    s <= lo && (s = lo + 0.5f0)
    return (s - lo) / (hi - lo)
end
