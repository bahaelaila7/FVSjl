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

# cr/sitset.f:479-497 — fill per-species site index (SITEAR) from the site species (ISISP) SI, scaled by the
# SITELO/SITEHI bounds. `sitear[i]<=0` means "not set by SITECODE keyword" => compute the default.
# `imodty` model type (1-5); `isisp` site-species index (0 => default per imodty). Modifies `sitear` in place.
const _CR_ISISP_DEFAULT = (3, 13, 13, 18, 11)     # IMODTY 1..5 -> DF/PP/PP/ES/LP
const _CR_TEM_DEFAULT   = (70.0f0, 70.0f0, 57.0f0, 75.0f0, 65.0f0)

function cr_site_index_defaults!(sitear::AbstractVector{Float32}, sd, imodty::Int, isisp::Int)
    lo = sd[:site_lo]; hi = sd[:site_hi]; n = length(sitear)
    tem = _CR_TEM_DEFAULT[imodty]
    (isisp > 0 && sitear[isisp] > 0.0f0) && (tem = sitear[isisp])
    isisp == 0 && (isisp = _CR_ISISP_DEFAULT[imodty])
    tem < lo[isisp] && (tem = lo[isisp])
    @inbounds for i in 1:n
        sitear[i] <= 0.0f0 && (sitear[i] = lo[i] + (tem - lo[isisp]) / (hi[isisp] - lo[isisp]) * (hi[i] - lo[i]))
    end
    return sitear
end

# cr/sitset.f DEFMT (23 forests -> IMODTY default when MODTYPE not given). IMODTY 1..5 =
# SW-mixed / SW-ponderosa / Black-Hills-ponderosa / spruce-fir / lodgepole.
const _CR_DEFMT = (5, 3, 4, 5, 3, 4, 5, 5, 4, 4, 5, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)

"""
    cr_site_index_setup!(s)

CR sitset.f: resolve the GENGYM model type (IMODTY = MODTYPE keyword override, else DEFMT[forest]),
then fan the site species' SI across all species via the SITELO/SITEHI conversion (`cr_site_index_defaults!`).
"""
function cr_site_index_setup!(s::StandState)
    p = s.plot; sd = s.coef.species
    imodty = Int(s.plot.model_type)
    if !(1 <= imodty <= 5)
        ifor = Int(p.forest_idx)
        imodty = (1 <= ifor <= length(_CR_DEFMT)) ? _CR_DEFMT[ifor] : 5
        s.plot.model_type = Int32(imodty)
    end
    cr_site_index_defaults!(p.sp_site_index, sd, imodty, Int(p.site_species))
    # SDIDEF (cr/sitset.f): per-species max SDI = SDICON unless BAMAX (or SDIMAX kw) overrode it.
    bamax = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi > 0.0f0 ? p.pct_sdimax_mort_hi : 0.85f0
    sdicon = sd[:sdi_max_default]
    @inbounds for i in 1:nspecies(s.variant)
        if p.sp_sdi_def[i] <= 0.0f0
            p.sp_sdi_def[i] = bamax > 0.0f0 ? bamax / (0.5454154f0 * pmsdiu) : sdicon[i]
        end
    end
    return s
end

site_setup!(s::StandState, ::CentralRockies) = cr_site_index_setup!(s)
