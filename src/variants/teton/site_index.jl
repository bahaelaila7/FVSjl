# =============================================================================
# site_index.jl (teton) — TT site index + SDImax + forest code (tt/sitset.f,
# tt/siterange.f, tt/forkod.f). Chunk 2.
#
# TT's site chunk is SIMPLER than EM's (no habtyp.f / habitat-type-group site
# index). SITSET interpolates the SITE SPECIES' site index (TEM) into each
# species' [SITELO,SITEHI] range:
#     pos = (TEM - SLOSSP) / (SHISSP - SLOSSP)      # SLOSSP/SHISSP = site species' range
#     SITEAR(I) = SITELO(I) + pos*(SITEHI(I)-SITELO(I))
# SDIDEF(I) = SDICON(I) (or BAMAX/(0.5454154·PMSDIU/100) when the BAMAX keyword is set).
# forkod maps KODFOR → IFOR (location subscript for DG DGFOR) via JFOR=[403,405,415,416]
# (+ reservation special-cases 7306→Bridger, 8107→Caribou); IGL=KFOR(IFOR) (KFOR uninit ⇒ 0).
#
# VALIDATED vs live ttt01 (forest 415→IFOR=3, site species DF=3, TEM=50 default): all 18
# SITEAR bit-match (WB 43.75→44, DF 50, PM 16.25→16, NC 97.5→98, …) and SDI MAX = SDICON.
# =============================================================================

# tt/sitset.f DATA SDICON (SDImax constant per species, 18).
const TT_SDICON = Float32[621, 409, 570, 358, 620, 562, 679, 620, 602, 446,
                          497, 411, 619, 680, 452, 501, 409, 452]
# tt/siterange.f DATA SITELO / SITEHI (per-species site-index range, 18).
const TT_SITELO = Float32[25, 25, 20,  5, 40, 30, 20, 40, 40, 40,  5,  5,  5,  5, 30,  5, 20,  5]
const TT_SITEHI = Float32[50, 50, 60, 20,100, 70,100,100, 90, 80, 15, 15, 30, 30,120, 15, 50, 20]

# tt/forkod.f: JFOR forest codes (KODFOR→IFOR subscript for DG DGFOR/DGDS). KFOR (geographic
# location IGL) has NO DATA statement in tt/forkod.f ⇒ zeros ⇒ IGL=0 for all TT forests.
const TT_JFOR = Int[403, 405, 415, 416]

# tt/forkod.f — KODFOR → IFOR (location) + IGL. Reservation special-cases map to a host NF.
function tt_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 1; useigl = true
    if kodfor == 7306
        ifor = 1                      # Wind River Reservation → Bridger NF (403)
    elseif kodfor == 8107
        ifor = 2                      # Fort Hall Reservation → Caribou NF (405)
    else
        idx = findfirst(==(kodfor), TT_JFOR)
        idx === nothing ? (useigl = false; ifor = 1) : (ifor = idx)
    end
    p.forest_idx = Int32(ifor)
    useigl && (p.geo_location = Int32(0))   # IGL = KFOR(IFOR); KFOR uninitialized ⇒ 0
    return ifor
end

# tt/sitset.f — fill SITEAR (site index by species) from the site species, + SDImax defaults.
function tt_sitset!(s::StandState)
    p = s.plot
    isisp = Int(p.site_species)
    tem = 50f0                                          # TEM default (tt/sitset.f:56)
    (isisp > 0 && p.sp_site_index[isisp] > 0f0) && (tem = p.sp_site_index[isisp])
    isisp == 0 && (isisp = 3; p.site_species = Int32(3))   # default site species = DF (tt/sitset.f:60)
    slossp = TT_SITELO[isisp]; shissp = TT_SITEHI[isisp]   # SITERANGE for the site species
    @inbounds for i in 1:nspecies(s.variant)
        tem < slossp && (tem = slossp)
        if p.sp_site_index[i] <= 0f0
            p.sp_site_index[i] = TT_SITELO[i] + (tem - slossp) / (shissp - slossp) * (TT_SITEHI[i] - TT_SITELO[i])
        end
    end
    # SDIDEF: SDICON per species, unless the BAMAX keyword is set (Zeide: BAMAX/(0.5454154·PMSDIU/100)).
    bamax_kw = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 85f0   # tt uses PMSDIU as a PERCENT
    @inbounds for i in 1:nspecies(s.variant)
        if p.sp_sdi_def[i] <= 0f0
            p.sp_sdi_def[i] = bamax_kw > 0f0 ? bamax_kw / (0.5454154f0 * (pmsdiu / 100f0)) : TT_SDICON[i]
        end
    end
    return s
end

function tt_site_index_setup!(s::StandState)
    tt_forkod!(s.plot)          # IFOR → p.forest_idx (DG DGFOR/DGDS); IGL → p.geo_location
    tt_sitset!(s)               # SITEAR (p.sp_site_index) + SDIDEF (p.sp_sdi_def)
    return s
end

site_setup!(s::StandState, ::Teton) = tt_site_index_setup!(s)
