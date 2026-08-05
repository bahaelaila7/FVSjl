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

# cr/forkod.f JFOR (national-forest code KODFOR=Region·100+Forest -> forest subscript IFOR = JFOR index).
# The DEFAULT case is `IFOR = findfirst(JFOR .== KODFOR)` (forkod.f:586-592); the first 23 index _CR_DEFMT.
const _CR_JFOR = (202, 203, 204, 206, 207, 209, 210, 211, 212, 213, 214, 215, 301, 302, 303, 304, 305, 306,
                  307, 308, 309, 310, 312, 201, 205, 208, 224, 311, 216)

"cr/forkod.f: map KODFOR (Region·100+Forest, `user_forest_code`) to the forest subscript IFOR via JFOR. Sets
`p.forest_idx` so the MODTYPE default DEFMT[IFOR] resolves — else DB-input stands (no MODTYPE) fall back to
IMODTY=5 (e.g. San Juan NF 213 → IFOR 10 → DEFMT 4=spruce-fir, not 5). Not-found leaves forest_idx unchanged."
# forkod.f:640-678 SECOND pass — consolidate the pseudo/duplicate NF subscripts into their combined units
# (Arapaho 201→Arapaho-Roosevelt, Gunnison 205→GMUG, Pike 208→Pike-San Isabel, Grand Mesa 224→GMUG,
# Sitgreaves 311→Apache-Sitgreaves, Routt 211→Medicine Bow-Routt, McKelvie 216→Nebraska). Keyed on the
# FIRST-pass IFOR (the _CR_JFOR index), then KODFOR=JFOR(IFOR). jl previously stopped after the first pass,
# leaving these 7 forests at their pseudo index (e.g. Pike 208→26) — INERT for growth on stands that supply
# elevation/MODTYPE, but wrong for the BFMIND/SCFMIND merch rule (IFOR<IGFOR=13) and any DEFMT/site default.
const _CR_FORKOD2 = Dict(24 => 7, 25 => 3, 26 => 9, 27 => 3, 28 => 13, 8 => 4, 29 => 5)
function _cr_forkod!(p)
    p.forest_idx > 0 && return p                         # already resolved (e.g. STDINFO keyword)
    kodfor = Int(p.user_forest_code)
    kodfor <= 0 && return p
    idx = findfirst(==(kodfor), _CR_JFOR)
    if idx === nothing
        # forkod.f CASE DEFAULT / .NOT.FORFOUND error trap for location codes not in JFOR (forkod.f:596-626):
        # pick a surrogate forest by the INPUT model type (MODTYPE keyword, 0 if none) so DEFMT[IFOR] resolves.
        # Without this, an out-of-CR-region DB stand (e.g. an R8 FIA stand) left forest_idx=0 ⇒ jl's imodty
        # fell back to 5 (lodgepole); FVS instead maps it to IFOR=15 (Cibola, DEFMT=2) or 10 (San Juan, DEFMT=4).
        input_imodty = Int(p.model_type)
        ifor = if input_imodty <= 2                      # CASE(:2) — includes the no-MODTYPE default 0
            (0 < kodfor < 300) ? 10 : 15                 # SAN JUAN (DEFMT 4) : CIBOLA (DEFMT 2)
        elseif input_imodty == 3
            2                                            # BIGHORN
        else
            4                                            # CASE DEFAULT (imodty 4/5) — GMUG
        end
        ifor = get(_CR_FORKOD2, ifor, ifor)
        p.forest_idx = Int32(ifor)
        p.user_forest_code = Int32(_CR_JFOR[ifor])
        return p
    end
    ifor = get(_CR_FORKOD2, idx, idx)                    # second-pass consolidation
    p.forest_idx = Int32(ifor)
    p.user_forest_code = Int32(_CR_JFOR[ifor])           # KODFOR = JFOR(IFOR) (forkod.f:680)
    return p
end

"""
    cr_site_index_setup!(s)

CR sitset.f: resolve the GENGYM model type (IMODTY = MODTYPE keyword override, else DEFMT[forest]),
then fan the site species' SI across all species via the SITELO/SITEHI conversion (`cr_site_index_defaults!`).
"""
function cr_site_index_setup!(s::StandState)
    p = s.plot; sd = s.coef.species
    _cr_forkod!(p)                       # resolve IFOR from KODFOR (JFOR) so DEFMT[IFOR] works for DB-input stands
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
