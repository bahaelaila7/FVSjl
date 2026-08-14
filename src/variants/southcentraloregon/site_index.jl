# =============================================================================
# site_index.jl (southcentraloregon) — SO site index + forest code + SICHG fan. Chunk 2.
# (so/forkod.f, so/sichg.f, so/htcalc.f, so/sitset.f)
#
# SO's R6 site chain mirrors EC: forkod maps the FVS location code to IFOR (and, for R6 forests IFOR≤3 or
# =10, resets LZEIDE=.FALSE. ⇒ Reineke); the site species' SITEAR is fanned to every unset species via
# SICHG (per-species reference age SIAGE) + so_htcalc (the site-species height curve inverted to SIAGE).
# so_htcalc is a per-species SELECT CASE of published site curves (Brickell WP / Cochran DF-WF-WL /
# Dolph-red-fir + Dunning-Levitan SH-WO / Dahms LP / Alexander ES / Herman NF / Wiley WH / Curtis misc /
# Harrington RA / Barrett SP-PP). sot01 gives SI=60 via STDINFO (NSISET>0) so the ECOCLS default path is
# inert; the ECOCLS/SDIDEF default-lookup (for no-SI stands) is deferred (sot01 SDImax rides the species
# sdimax column via stand_sdimax, Reineke).
# =============================================================================

# so/forkod.f — location code → IFOR via the JFOR table (11 R6-Oregon/NE-Cal forests).
const SO_JFOR = Int[601, 602, 620, 505, 506, 509, 511, 701, 514, 799, 702]

function so_forkod!(p)
    kodfor = Int(p.user_forest_code)
    ifor = 0
    # reservation pseudo-code crosswalk (so/forkod.f): 7710→602, others → R6 forests
    if kodfor == 7711 || kodfor == 7710
        ifor = 2                       # Fort McDermitt → Fremont(602)=IFOR 2
    else
        for (i, f) in enumerate(SO_JFOR)
            kodfor == f && (ifor = i; break)
        end
        ifor == 0 && (ifor = 1)        # default Deschutes(601)=IFOR 1
    end
    p.forest_idx = Int32(ifor)
    p.user_forest_code = Int32(SO_JFOR[ifor])
    return ifor
end

# so/htcalc.f — predicted total height (ft) at age `ag` for species `ispc`, site index `sindx`, forest `jfor`.
const SO_DUNL1 = Float32[-88.9, -82.2, -78.3, -82.1, -56.0, -33.8]
const SO_DUNL2 = Float32[49.7067, 44.1147, 39.1441, 35.4160, 26.7173, 18.6400]
const SO_DUNL3 = Float32[2.375, 2.025, 1.650, 1.225, 1.075, 0.875]

@inline function _so_dunlevitan(sindx::Float32, ag::Float32)::Float32
    indx = sindx <= 44f0 ? 6 : sindx <= 52f0 ? 5 : sindx <= 65f0 ? 4 :
           sindx <= 82f0 ? 3 : sindx <= 98f0 ? 2 : 1
    return ag <= 40f0 ? SO_DUNL3[indx] * ag : SO_DUNL1[indx] + SO_DUNL2[indx] * log(ag)
end

function so_htcalc(jfor::Int, sindx::Float32, ispc::Int, ag::Float32)::Float32
    if ispc == 1                                                  # WP Brickell
        return sindx / (0.37504453f0 * (1f0 - 0.92503f0 * exp(-0.0207959f0 * ag))^(-2.4881068f0))
    elseif ispc == 3 || ispc == 32                                # DF/OS Cochran-251
        return 4.5f0 + exp(-0.37496f0 + 1.36164f0 * log(ag) - 0.00243434f0 * log(ag)^4) -
               79.97f0 * (-0.2828f0 + 1.87947f0 * (1f0 - exp(-0.022399f0 * ag))^0.966998f0) +
               (sindx - 4.5f0) * (-0.2828f0 + 1.87947f0 * (1f0 - exp(-0.022399f0 * ag)^0.966998f0))
    elseif ispc == 4 || ispc == 6 || ispc == 12 || ispc == 14     # WF/IC/GF/SF Cochran-252
        la = log(ag)
        x2 = -0.30935f0 + 1.2383f0 * la + 0.001762f0 * la^4 - 5.4f-6 * la^9 + 2.046f-7 * la^11 - 4.04f-13 * la^18
        x3 = -6.2056f0 + 2.097f0 * la - 0.09411f0 * la^2 - 0.00004382f0 * la^7 + 2.007f-11 * la^16 - 2.054f-17 * la^24
        return exp(x2) - 84.73f0 * exp(x3) + (sindx - 4.5f0) * exp(x3) + 4.5f0
    elseif ispc == 5                                              # MH interim-means (metric)
        h = (22.8741f0 + 0.950234f0 * sindx) * (1f0 - exp(-0.00206465f0 * sqrt(sindx) * ag))^(1.365566f0 + 2.045963f0 / sindx)
        return (h + 1.37f0) * 3.281f0
    elseif ispc == 7                                              # LP Dahms
        return sindx * (-0.0968f0 + 0.02679f0 * ag - 0.00009309f0 * ag * ag)
    elseif ispc == 8                                              # ES Alexander
        return 4.5f0 + ((2.75780f0 * sindx^0.83312f0) * (1f0 - exp(-0.015701f0 * ag))^(22.71944f0 * sindx^(-0.63557f0)))
    elseif ispc == 9                                              # SH: R6 Dolph red fir else Dunning-Levitan
        if jfor <= 3 || jfor == 10
            term = ag * exp(ag * (-0.0440853f0)) * 1.41512f-6
            b = sindx * term - 3.04951f6 * term * term + 5.72474f-4
            term2 = 50f0 * exp(50f0 * (-0.0440853f0)) * 1.41512f-6
            b50 = sindx * term2 - 3.04951f6 * term2 * term2 + 5.72474f-4
            return (sindx - 4.5f0) * (1f0 - exp(-b * ag^1.51744f0)) / (1f0 - exp(-b50 * 50f0^1.51744f0)) + 4.5f0
        else
            return _so_dunlevitan(sindx, ag)
        end
    elseif ispc == 2 || ispc == 10                                # SP/PP Barrett
        t = -0.7864f0 + 2.49717f0 * (1f0 - exp(-0.0045042f0 * ag))^0.33022f0
        return (128.89522f0 * (1f0 - exp(-0.016959f0 * ag))^1.23114f0) - (t * 100.43f0) + (t * (sindx - 4.5f0)) + 4.5f0
    elseif ispc == 13                                             # AF Johnson/DeMars
        return sindx * (-0.07831f0 + 0.0149f0 * ag - 4.0818f-5 * ag * ag)
    elseif ispc == 15                                             # NF Herman
        x1 = -564.38f0 + 22.250f0 * (sindx - 4.5f0) - 0.04995f0 * (sindx - 4.5f0)^2
        x2 = 6.8f0 + 2843.21f0 / (sindx - 4.5f0) + 34735.54f0 / (sindx - 4.5f0)^2
        return 4.5f0 + (sindx - 4.5f0) / (x1 * (1f0 / ag)^2 + x2 * (1f0 / ag) + 1f0 - 0.0001f0 * x1 - 0.01f0 * x2)
    elseif ispc == 17                                             # WL Cochran-424
        t = -0.12528f0 + 0.039636f0 * ag - 0.0004278f0 * ag * ag + 1.7039f-6 * ag^3
        return 4.5f0 + 1.46897f0 * ag + 0.0092466f0 * ag * ag - 0.00023957f0 * ag^3 + 1.1122f-6 * ag^4 +
               (sindx - 4.5f0) * t - 73.57f0 * t
    elseif ispc == 18                                             # RC Hegyi
        return 1.3283f0 * sindx * ((1f0 - exp(-0.0174f0 * ag))^1.4711f0)
    elseif ispc == 19                                             # WH Wiley
        z = 2500f0 / (sindx - 4.5f0)
        return ag * ag / (-1.73070f0 + 0.1394f0 * z + (-0.0616f0 + 0.0137f0 * z) * ag +
               (0.00192f0 + 0.00007f0 * z) * (ag * ag)) + 4.5f0
    elseif ispc == 20 || ispc == 21 || ispc == 23 || ispc == 25 || ispc == 26 ||
           (28 <= ispc <= 31) || ispc == 33                       # misc Curtis
        return (sindx - 4.5f0) / (0.6192f0 - 5.3394f0 / (sindx - 4.5f0) + 240.29f0 * ag^(-1.4f0) +
               (3368.9f0 / (sindx - 4.5f0)) * ag^(-1.4f0)) + 4.5f0
    elseif ispc == 22                                             # RA Harrington
        c = 59.5864f0 + 0.79530f0 * sindx
        return sindx + c * (1f0 - exp((0.00194f0 - 0.000740f0 * sindx) * ag))^0.9198f0 -
               c * (1f0 - exp((0.00194f0 - 0.000740f0 * sindx) * 20f0))^0.9198f0
    elseif ispc == 27                                             # WO: R6 Powers else Dunning-Levitan
        if jfor <= 3 || jfor == 10
            term = sqrt(ag) - sqrt(50f0)
            return (sindx * (1f0 + 0.322f0 * term)) - 6.413f0 * term
        else
            return _so_dunlevitan(sindx, ag)
        end
    else                                                          # WJ(11)/WB(16)/AS(24) + future = no curve
        return 0f0
    end
end

# so/siterange.f — per-species site-index range (SITELO/SITEHI); used to interpolate SITEAR for the
# no-curve species WJ(11)/WB(16)/AS(24).
const SO_SITELO = Float32[13,27,21,5,5,5,5,12,10,7, 5,9,6,4,7,20,60,29,6,5, 5,56,108,30,10,10,21,20,5,5, 5,5,5]
const SO_SITEHI = Float32[137,178,148,195,133,169,140,227,134,176, 40,173,127,221,210,65,147,152,203,75,
                          100,192,142,66,191,104,85,93,100,75, 75,175,125]

# so/sichg.f — SIAGE(i) per species (reference age for the site-species curve). RF(5) is metric.
function so_sichg(s::StandState, isisp::Integer, ssite::Float32)
    sd = s.coef.species
    a = sd[:sichg_a]; b = sd[:sichg_b]; refage = sd[:sichg_refage]; refloc = sd[:sichg_refloc]
    maxsp = nspecies(s.variant)
    isiloc = refloc[isisp]
    siage = Vector{Float32}(undef, maxsp)
    @inbounds for i in 1:maxsp
        diff = 0
        (isiloc == 1f0 && refloc[i] == 0f0) && (diff = -1)        # site 'T', species 'B'
        (isiloc == 0f0 && refloc[i] == 1f0) && (diff = 1)         # site 'B', species 'T'
        age2bh = 0f0
        if diff != 0
            age2bh = a[i] + b[i] * ssite
            (isisp != 5 && i == 5) && (age2bh = a[i] + b[i] * (ssite / 3.281f0))
            (isisp == 5 && i != 5) && (age2bh = a[i] + b[i] * (ssite * 3.281f0))
        end
        i == 9 && (age2bh = 18f0)                                 # SH fixed
        i == 27 && (age2bh = 2f0)                                 # WO fixed
        siage[i] = refage[i] + age2bh * Float32(diff)
    end
    return siage
end

function so_sitset!(s::StandState)
    p = s.plot; maxsp = nspecies(s.variant)
    ifor = Int(p.forest_idx)
    # so/sitset.f:81 — R6 forests (IFOR≤3 or =10) use Reineke (LZEIDE=.FALSE.); others keep Zeide.
    (ifor <= 3 || ifor == 10) && (s.control.zeide_sdi = false)

    isisp = (1 <= Int(p.site_species) <= maxsp) ? Int(p.site_species) : 0
    isisp <= 0 && (isisp = 10)                                    # R6 default site species = PP(10)
    p.sp_site_index[isisp] <= 0f0 && (p.sp_site_index[isisp] = 70f0)   # so/sitset.f default SITEAR

    sindx = p.sp_site_index[isisp]
    siage = so_sichg(s, isisp, sindx)
    slossp = SO_SITELO[isisp]; shissp = SO_SITEHI[isisp]
    @inbounds for ispc in 1:maxsp
        ispc == isisp && continue
        if ispc == 11 || ispc == 16 || ispc == 24       # WJ/WB/AS: no site curve → SITERANGE interpolation
            tem = sindx < slossp ? slossp : sindx
            v = SO_SITELO[ispc] + (tem - slossp) / (shissp - slossp) * (SO_SITEHI[ispc] - SO_SITELO[ispc])
        else
            v = so_htcalc(ifor, sindx, isisp, siage[ispc])
            ispc == 5 && (v = v / 3.281f0; v > 28f0 && (v = 28f0))   # MH metric conversion, cap 28
        end
        v < 0f0 && (v = 0f0)
        p.sp_site_index[ispc] = v                        # SO sitset is authoritative (overrides the generic reader)
    end
    p.site_species = Int32(isisp)
    return s
end

function so_site_index_setup!(s::StandState)
    so_forkod!(s.plot)
    so_sitset!(s)
    return s
end

site_setup!(s::StandState, ::SouthCentralOregon) = so_site_index_setup!(s)
