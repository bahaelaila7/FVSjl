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

site_setup!(s::StandState, ::Northeast) = ne_site_index_setup!(s)
