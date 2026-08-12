# =============================================================================
# site_index.jl (utah) — UT site index + SDImax + forest code (ut/sitset.f,
# ut/siterange.f, ut/forkod.f, ut/habtyp.f). Chunk 2.
#
# UT's site chunk mirrors TT's (both Region-4): SITSET interpolates the SITE SPECIES'
# site index (TEM) into each species' [SITELO,SITEHI] range —
#     pos = (TEM - SLOSSP) / (SHISSP - SLOSSP)          # SLOSSP/SHISSP = site species' range
#     SITEAR(I) = SITELO(I) + pos*(SITEHI(I)-SITELO(I))
# SDIDEF(I) = SDICON(I) (or BAMAX/(0.5454154·PMSDIU/100) when the BAMAX keyword is set; Zeide).
# TEM default 46, site species (ISISP) default 7 (LP) — ut/sitset.f:56-60.
# HABTYP: KODTYP → ITYPE via CRDECD against R4HABT(363). ★ ut/habtyp.f R4HABT is BYTE-IDENTICAL to
# tt/habtyp.f — REUSE TT_R4HABT_CODE directly (verified by diffing the DATA blocks). ITYPE feeds the
# DG habitat coefficient dimension in later chunks.
# forkod maps KODFOR → IFOR (location subscript for DG DGFOR/DGDS) via JFOR + 16 reservation pseudo-codes,
# then the "forest mapping correction" (Cache 7→6 Wasatch, Humboldt 8→3 Fishlake, Toiyabe 9→3 Fishlake);
# IGL = KFOR(IFOR) = 1.
# =============================================================================

# ut/sitset.f DATA SDICON (SDImax constant per species, 24).
const UT_SDICON = Float32[621, 409, 570, 634, 620, 562, 679, 620, 602, 446,
                          348, 272, 652, 358, 411, 497, 621, 452, 452, 501,
                          619, 344, 409, 652]
# ut/siterange.f DATA SITELO / SITEHI (sitset range, IWHO=1; distinct from the regent SLO/SHI, 24).
const UT_SITELO = Float32[20, 20, 30, 40, 40, 30, 40, 40, 50, 40,  5,  5,  5,  5,  5,  5, 20, 30, 30,  5,  5, 15, 20,  5]
const UT_SITEHI = Float32[50, 50, 70, 80, 90, 70, 85,100, 90, 80, 20, 15, 20, 20, 15, 15, 60,120, 90, 15, 30, 40, 50, 20]

# ut/forkod.f DATA JFOR (KODFOR→IFOR subscript for DG DGFOR/DGDS), KFOR all 1 (IGL=1). NUMFOR=9.
const UT_JFOR = Int[401, 407, 408, 410, 418, 419, 404, 409, 417]
# Reservation pseudo-codes → pre-correction IFOR (ut/forkod.f CASE list).
const UT_FOR_RESERV = Dict{Int,Int}(
    7702=>9, 7711=>8, 7712=>8, 7714=>8, 7715=>9, 7716=>9, 7717=>6, 7718=>1,
    7721=>8, 7722=>8, 7723=>8, 7728=>2, 7729=>2, 7835=>9, 7920=>4, 8001=>4)

# ut/forkod.f — KODFOR → IFOR (location) + IGL. Reservations + JFOR match, then the mapping correction.
function ut_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 1; useigl = true
    r = get(UT_FOR_RESERV, kodfor, 0)
    if r != 0
        ifor = r
    else
        idx = findfirst(==(kodfor), UT_JFOR)
        idx === nothing ? (useigl = false; ifor = 1) : (ifor = idx)   # not found → ERRGRO(3), IFOR stays 1
    end
    # Forest mapping correction (ut/forkod.f:172-188): Cache 7→6, Humboldt 8→3, Toiyabe 9→3.
    ifor == 7 ? (ifor = 6) : ifor == 8 ? (ifor = 3) : ifor == 9 && (ifor = 3)
    p.forest_idx = Int32(ifor)
    useigl && (p.geo_location = Int32(1))   # IGL = KFOR(IFOR) = 1
    return ifor
end

# ut/habtyp.f + crdecd.f — decode KODTYP into ITYPE via the R4HABT(363) table (reused from TT, identical).
# Returns 0 if unmatched (KODTYP≤0 or absent) — CRDECD IHB=-1 leaves ITYPE 0 (default).
function ut_habtyp(kodtyp::Integer)::Int32
    kodtyp <= 0 && return Int32(0)
    @inbounds for i in 1:length(TT_R4HABT_CODE)
        TT_R4HABT_CODE[i] == kodtyp && return Int32(i)
    end
    return Int32(0)
end

# ut/sitset.f — fill SITEAR (site index by species) from the site species, + SDImax defaults.
function ut_sitset!(s::StandState)
    p = s.plot
    isisp = Int(p.site_species)
    tem = 46f0                                          # TEM default (ut/sitset.f:56)
    (isisp > 0 && p.sp_site_index[isisp] > 0f0) && (tem = p.sp_site_index[isisp])
    isisp == 0 && (isisp = 7; p.site_species = Int32(7))   # default site species = LP (ut/sitset.f:60)
    slossp = UT_SITELO[isisp]; shissp = UT_SITEHI[isisp]   # SITERANGE for the site species
    @inbounds for i in 1:nspecies(s.variant)
        tem < slossp && (tem = slossp)
        if p.sp_site_index[i] <= 0f0
            p.sp_site_index[i] = UT_SITELO[i] + (tem - slossp) / (shissp - slossp) * (UT_SITEHI[i] - UT_SITELO[i])
        end
    end
    # SDIDEF: SDICON per species, unless the BAMAX keyword is set (Zeide: BAMAX/(0.5454154·PMSDIU/100)).
    bamax_kw = s.control.ba_max
    pmsdiu = p.pct_sdimax_mort_hi > 0f0 ? p.pct_sdimax_mort_hi : 85f0   # ut uses PMSDIU as a PERCENT
    @inbounds for i in 1:nspecies(s.variant)
        if p.sp_sdi_def[i] <= 0f0
            p.sp_sdi_def[i] = bamax_kw > 0f0 ? bamax_kw / (0.5454154f0 * (pmsdiu / 100f0)) : UT_SDICON[i]
        end
    end
    return s
end

# ut/cratet.f — adjust SITEAR to a 50-YEAR age base for the species whose growth eqns were fit on a
# 50-yr-base site index: WB/LM/LP/OS (1,2,7,23) via Alexander-Tackle-Dahms RM-29; WF/BS/ES/AF (4,5,8,9)
# via Alexander RM-32; PP (10) via Meyer 1961 (TB-630). Uses stand CCF (floored 125, DBH-only open-grown).
# Inert unless one of these species/site-species is present (why utt01 — PJ/woodland — was bit-exact).
function ut_cratet_site_adjust!(s::StandState)
    p, t = s.plot, s.trees
    temccf = 0f0
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0f0 && continue
        temccf += ut_tree_ccf(Int(t.species[i]), t.dbh[i]) * t.tpa[i]
    end
    temccf < 125f0 && (temccf = 125f0)
    @inbounds for sp in 1:nspecies(s.variant)
        si = p.sp_site_index[sp]
        if sp == 1 || sp == 2 || sp == 7 || sp == 23          # Alexander-Tackle-Dahms RM-29
            p.sp_site_index[sp] = 9.89311f0 - 0.19177f0 * 50f0 + 0.00124f0 * 50f0^2 -
                0.00082f0 * (temccf - 125f0) * si + 0.01387f0 * 50f0 * si -
                0.0000455f0 * 50f0^2 * si
        elseif sp == 4 || sp == 5 || sp == 8 || sp == 9       # Alexander RM-32
            p.sp_site_index[sp] = 4.5f0 + (2.75780f0 * si^0.83312f0) *
                (1f0 - exp(-0.015701f0 * 50f0))^(22.71944f0 * si^(-0.63557f0))
        elseif sp == 10                                        # Meyer 1961 (TB-630) — PP
            p.sp_site_index[sp] = (3.635794f0 * si^0.916307f0) /
                (1f0 + exp(6.09478f0 - 0.96483f0 * log(50f0) - 0.277025f0 * log(si)))
        end
    end
    return s
end

function ut_site_index_setup!(s::StandState)
    ut_forkod!(s.plot)          # IFOR → p.forest_idx (DG DGFOR/DGDS); IGL → p.geo_location
    s.plot.habitat_input = ut_habtyp(Int(s.plot.habitat_code))   # KODTYP → ITYPE (ut/habtyp.f)
    ut_sitset!(s)               # SITEAR (p.sp_site_index) + SDIDEF (p.sp_sdi_def)
    ut_cratet_site_adjust!(s)   # CRATET 50-yr-base site adjust (WB/LM/WF/BS/LP/ES/AF/PP/OS)
    return s
end

site_setup!(s::StandState, ::Utah) = ut_site_index_setup!(s)
