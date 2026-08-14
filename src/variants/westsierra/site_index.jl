# =============================================================================
# site_index.jl (westsierra) — WS site index + forest code + SICHG fan. Chunk 2.
# (ws/forkod.f, ws/sichg.f, ws/htcalc.f, ws/sitset.f)
#
# WS's California R5 site chain mirrors SO/EC: forkod maps the FVS location code to IFOR (NEAREST-forest
# MINLOC, unlike SO's exact match); the site species' SITEAR is fanned to every unset species via SICHG
# (per-species reference age SIAGE) + ws_htcalc (the site-species height curve at SIAGE). ws_htcalc is a
# COMPACT 5-group SELECT CASE (all 43 species covered, no no-curve species): Dunning-Levitan (most conifers),
# Powers-oak ×0.80 (the oak/hardwood group), Dolph red-fir (RF sp7), Curtis (MC sp41), Alexander (GB sp21).
# WS is Zeide SDI throughout (ws/grinit.f LZEIDE=.TRUE., NO sitset reset). sp21(GB)/sp41(MC) get a special
# SI(7)-based interpolation instead of their own curve (ws/sitset.f DO30). SDImax (SDICON) fan = follow-on.
# =============================================================================

# ws/forkod.f — location code → IFOR. 5-digit KODFOR (≥40000): KFOR1=KODFOR÷100, MINLOC nearest JFOR.
# 3-digit forest code (wst01 = 511): EXACT match KODFOR==JFOR(i). Then a post-remap SELECT CASE(IFOR)
# folds NFs 7-11→3, 12→1, 13→5. (Tribal reservation pseudo-codes 7712…7860→specific IFOR = a follow-on;
# wst01 is a plain forest code.) ⚠ EARLIER MINLOC-only ws_forkod! was WRONG (gave IFOR=13 for 511; the
# SITEAR 43/43 check MISSED it because ws/htcalc.f ignores IFOR) — this is the exact-match + gate + remap.
const WS_JFOR = Int[503, 511, 513, 515, 516, 517, 501, 502, 504, 507, 512, 519, 417]

# ws/sitset.f DATA SDICON — per-species default SDImax (used by the DO40 fan when no BAMAX/SDIMAX keyword).
const WS_SDICON = Float32[
  561, 570, 800, 1052, 576, 365, 1000, 365, 679, 621,
  272, 358, 790, 679, 409, 365, 409, 365, 214, 365,
  409, 570, 1052, 687, 272, 497, 272, 667, 667, 214,
  406, 440, 667, 785, 785, 562, 406, 515, 406, 629,
  501, 365, 406]

@inline function _ws_forkod_remap(ifor::Int)::Int      # ws/forkod.f final SELECT CASE(IFOR)
    (7 <= ifor <= 11) && return 3
    ifor == 12 && return 1
    ifor == 13 && return 5
    return ifor
end

function ws_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 0
    if kodfor >= 40000                                 # 5-digit forest×100+district → nearest JFOR
        kfor1 = kodfor ÷ 100
        best = typemax(Int)
        for (i, f) in enumerate(WS_JFOR); d = abs(f - kfor1); d < best && (best = d; ifor = i); end
    else                                               # DEFAULT: exact 3-digit forest-code match
        for (i, f) in enumerate(WS_JFOR); kodfor == f && (ifor = i; break); end
    end
    ifor == 0 && (ifor = 1)                            # not-found fallback (ws errgro path)
    ifor = _ws_forkod_remap(ifor)
    p.forest_idx = Int32(ifor)
    return ifor
end

# ws/htcalc.f — predicted total height (ft) at age `ag` for species `ispc`, site index `sindx`.
const WS_DUNL1 = Float32[-88.9, -82.2, -78.3, -82.1, -56.0, -33.8]
const WS_DUNL2 = Float32[49.7067, 44.1147, 39.1441, 35.4160, 26.7173, 18.6400]
const WS_DUNL3 = Float32[2.375, 2.025, 1.650, 1.225, 1.075, 0.875]

function ws_htcalc(ifor::Int, sindx::Float32, ispc::Int, ag::Float32)::Float32
    if (1 <= ispc <= 6) || (8 <= ispc <= 20) || (22 <= ispc <= 27) || ispc == 42
        indx = sindx <= 44f0 ? 6 : sindx <= 52f0 ? 5 : sindx <= 65f0 ? 4 :
               sindx <= 82f0 ? 3 : sindx <= 98f0 ? 2 : 1
        return ag <= 40f0 ? WS_DUNL3[indx] * ag : WS_DUNL1[indx] + WS_DUNL2[indx] * log(ag)
    elseif (28 <= ispc <= 40) || ispc == 43                     # Powers oak
        a = sqrt(ag) - sqrt(50f0)
        return (sindx * (1f0 + 0.322f0 * a) - 6.413f0 * a) * 0.80f0
    elseif ispc == 7                                            # Dolph red-fir
        term = ag * exp(ag * (-0.0440853f0)) * 1.4151f-6
        b = sindx * term - 3.0495f6 * term * term + 5.72474f-4
        term2 = 50f0 * exp(50f0 * (-0.0440853f0)) * 1.4151f-6
        b50 = sindx * term2 - 3.0495f6 * term2 * term2 + 5.72474f-4
        return (sindx - 4.5f0) * (1f0 - exp(-b * fpow(ag, 1.51744f0))) /
               (1f0 - exp(-b50 * fpow(50f0, 1.51744f0))) + 4.5f0
    elseif ispc == 41                                           # Curtis
        return (sindx - 4.5f0) / (0.6192f0 - 5.3394f0 / (sindx - 4.5f0) +
               240.29f0 * fpow(ag, -1.4f0) + (3368.9f0 / (sindx - 4.5f0)) * fpow(ag, -1.4f0)) + 4.5f0
    elseif ispc == 21                                           # Alexander (bristlecone GB)
        return 4.5f0 + (2.75780f0 * fpow(sindx, 0.83312f0)) *
               fpow(1f0 - exp(-0.015701f0 * ag), 22.71944f0 * fpow(sindx, -0.63557f0))
    end
    return 0f0
end

# ws/sichg.f — SIAGE(i) per species (reference age for the site-species curve). refloc 1='T', 0='B'.
function ws_sichg(s::StandState, isisp::Integer, ssite::Float32)
    sd = s.coef.species
    a = sd[:sichg_a]; b = sd[:sichg_b]; refage = sd[:sichg_refage]; refloc = sd[:sichg_refloc]
    maxsp = nspecies(s.variant)
    isiloc = refloc[isisp]
    siage = Vector{Float32}(undef, maxsp)
    @inbounds for i in 1:maxsp
        diff = 0
        (isiloc == 1f0 && refloc[i] == 0f0) && (diff = -1)       # site 'T', species 'B'
        (isiloc == 0f0 && refloc[i] == 1f0) && (diff = 1)        # site 'B', species 'T'
        age2bh = diff != 0 ? a[i] + b[i] * ssite : 0f0
        siage[i] = refage[i] + age2bh * Float32(diff)
    end
    return siage
end

function ws_sitset!(s::StandState)
    p = s.plot; maxsp = nspecies(s.variant)
    ifor = Int(p.forest_idx)
    # WS is Zeide throughout (ws/grinit.f LZEIDE=.TRUE., no IFOR reset). zeide_sdi already true from ws_grinit!.

    isisp = (1 <= Int(p.site_species) <= maxsp) ? Int(p.site_species) : 0
    isisp <= 0 && (isisp = 8)                                    # ws/sitset.f default site species = PP(8)
    p.sp_site_index[isisp] <= 0f0 && (p.sp_site_index[isisp] = 100f0)   # ws/sitset.f default SITEAR

    sindx = p.sp_site_index[isisp]
    siage = ws_sichg(s, isisp, sindx)
    # SI(i) = the site-species curve (ISISP) at species-i's reference age SIAGE(i).
    si = Vector{Float32}(undef, maxsp)
    @inbounds for ispc in 1:maxsp
        si[ispc] = ws_htcalc(ifor, sindx, isisp, siage[ispc])
    end
    # Fan SITEAR (ws/sitset.f DO30). Authoritative override of every non-site species (like so_sitset!).
    @inbounds for i in 1:maxsp
        i == isisp && continue
        v = if i == 21                                          # GB: SI(7)-based interp → [20,60]
            t = si[7]; t < 10f0 && (t = 10f0); t > 134f0 && (t = 134f0)
            20f0 + (t - 10f0) / (134f0 - 10f0) * (60f0 - 20f0)
        elseif i == 41                                          # MC: SI(7)-based interp → [5,23]
            t = si[7]; t < 10f0 && (t = 10f0); t > 134f0 && (t = 134f0)
            5f0 + (t - 10f0) / (134f0 - 10f0) * (23f0 - 5f0)
        else
            si[i]
        end
        v < 0f0 && (v = 0f0)
        p.sp_site_index[i] = v
    end
    # SDIDEF (per-species SDImax) fan — ws/sitset.f DO40: BAMAX>0 ? BAMAX/(0.5454154·PMSDIU/100) : SDICON.
    # Needed by crown/mortality RELSDI (= sdiac/sp_sdi_def). wst01 sets no BAMAX ⇒ SDIDEF = SDICON per species.
    bamax = p.ba_max; pmsdiu = p.pct_sdimax_mort_hi
    @inbounds for i in 1:maxsp
        p.sp_sdi_def[i] > 0f0 && continue                        # already set (SDIMAX keyword) → keep
        p.sp_sdi_def[i] = bamax > 0f0 ? bamax / (0.5454154f0 * (pmsdiu / 100f0)) : WS_SDICON[i]
    end
    p.site_species = Int32(isisp)
    return s
end

function ws_site_index_setup!(s::StandState)
    ws_forkod!(s.plot)
    ws_sitset!(s)
    return s
end

site_setup!(s::StandState, ::WestSierra) = ws_site_index_setup!(s)
