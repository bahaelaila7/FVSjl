# =============================================================================
# height_growth.jl (eastcascades) — EC potential-height curves (ec/htcalc.f) + large-tree HTG. Chunks 4a/4b.
#
# ec_htcalc reproduces ec/htcalc.f's per-species SELECT CASE potential-height curves (used by SITSET to
# fan SITEAR, and by findag/htgf). ec_findag (ec/findag.f) + height_growth! (ec/htgf.f) are the large-tree
# height-growth driver: findag inverts the species' OWN site curve to a stand age, htgf grows potential
# height 10 yr and applies the Hoerl (CRA/CRB/CRC) + generalized-Chapman-Richards (RHR/RHYXS/RHM/RHB)
# modifiers. Unlike WC/PN, EC's findag uses per-species AGMAX/AHMAX/BHMAX (not MAPHD/HDRAT). All Float32.
# =============================================================================

"""
    ec_htcalc(sindx, ispc, ag) -> Float32

Potential height (ft) for species `ispc` at site index `sindx` and age `ag` (ec/htcalc.f). CASE mapping:
WP(1) Brickell; WL/LL(2,17) Cochran PNW-424; DF(3) Cochran PNW-251; SF/GF/WF(4,6,16) Cochran PNW-252;
RC(5) Hegyi; LP(7) Alexander RM-29; ES(8) Alexander RM-32; AF(9) Johnson/DeMars PNW-119; PP(10) Barrett
PNW-232; WH(11) Wiley; MH/OS(12,31) Means metric; Curtis(13,14,18:21,23:27,29,30,32); NF(15) Herman
PNW-243; RA(22) Harrington PNW-358; WO(28) King. Default 0.
"""
function ec_htcalc(sindx::Float32, ispc::Int, ag::Float32)::Float32
    S = sindx; A = ag; L = log(A)
    if ispc == 1                                   # WP — Brickell INT-75
        return S / (0.37504453f0 * (1f0 - 0.92503f0 * exp(-0.0207959f0 * A))^(-2.4881068f0))
    elseif ispc == 2 || ispc == 17                 # WL/LL — Cochran PNW-424
        c = -0.12528f0 + 0.039636f0 * A - 0.0004278f0 * A * A + 1.7039f-6 * A^3
        return 4.5f0 + 1.46897f0 * A + 0.0092466f0 * A * A - 0.00023957f0 * A^3 +
               1.1122f-6 * A^4 + (S - 4.5f0) * c - 73.57f0 * c
    elseif ispc == 3                               # DF — Cochran PNW-251 (asymmetric parens verbatim)
        return 4.5f0 + exp(-0.37496f0 + 1.36164f0 * L - 0.00243434f0 * L^4) -
               79.97f0 * (-0.2828f0 + 1.87947f0 * (1f0 - exp(-0.022399f0 * A))^0.966998f0) +
               (S - 4.5f0) * (-0.2828f0 + 1.87947f0 * (1f0 - exp(-0.022399f0 * A)^0.966998f0))
    elseif ispc == 4 || ispc == 6 || ispc == 16    # SF/GF/WF — Cochran PNW-252
        x2 = -0.30935f0 + 1.2383f0 * L + 0.001762f0 * L^4 - 5.4f-6 * L^9 +
             2.046f-7 * L^11 - 4.04f-13 * L^18
        x3 = -6.2056f0 + 2.097f0 * L - 0.09411f0 * L^2 - 0.00004382f0 * L^7 +
             2.007f-11 * L^16 - 2.054f-17 * L^24
        return exp(x2) - 84.73f0 * exp(x3) + (S - 4.5f0) * exp(x3) + 4.5f0
    elseif ispc == 5                               # RC — Hegyi et al. 1981
        return 1.3283f0 * S * ((1f0 - exp(-0.0174f0 * A))^1.4711f0)
    elseif ispc == 7                               # LP — Alexander/Tackle/Dahms RM-29
        return 9.89331f0 - 0.19177f0 * A + 0.00124f0 * A * A - 0.00082f0 * 0f0 * S +
               0.01387f0 * A * S - 0.0000455f0 * A * A * S
    elseif ispc == 8                               # ES — Alexander RM-32
        return 4.5f0 + ((2.75780f0 * S^0.83312f0) *
               (1f0 - exp(-0.015701f0 * A))^(22.71944f0 * S^(-0.63557f0)))
    elseif ispc == 9                               # AF — Johnson/DeMars PNW-119
        return S * (-0.07831f0 + 0.0149f0 * A - 4.0818f-5 * A * A)
    elseif ispc == 10                              # PP — Barrett PNW-232
        return (128.8952205f0 * (1f0 - exp(-0.016959f0 * A))^1.23114f0) -
               ((-0.7864f0 + 2.49717f0 * (1f0 - exp(-0.004504f0 * A))^0.33022f0) * 100.43f0) +
               ((-0.7864f0 + 2.49717f0 * (1f0 - exp(-0.004504f0 * A))^0.33022f0) * (S - 4.5f0)) + 4.5f0
    elseif ispc == 11                              # WH — Wiley 1978
        z = 2500f0 / (S - 4.5f0)
        return (A * A / (-1.7307f0 + 0.1394f0 * z + (-0.0616f0 + 0.0137f0 * z) * A +
               (0.00192f0 + 0.00007f0 * z) * (A * A))) + 4.5f0
    elseif ispc == 12 || ispc == 31                # MH/OS — Means (metric, +1.37 m → ft)
        h = (22.8741f0 + 0.950234f0 * S) *
            (1f0 - exp(-0.00206465f0 * sqrt(S) * A))^(1.365566f0 + 2.045963f0 / S)
        return (h + 1.37f0) * 3.281f0
    elseif ispc == 15                              # NF — Herman PNW-243
        x1 = -564.38f0 + 22.25f0 * (S - 4.5f0) - 0.04995f0 * (S - 4.5f0)^2
        x2 = 6.80f0 + 2843.21f0 * (S - 4.5f0)^(-1) + 34735.54f0 * (S - 4.5f0)^(-2)
        return 4.5f0 + (S - 4.5f0) / (x1 * (1f0 / A)^2 + x2 * (1f0 / A) + 1f0 -
               0.0001f0 * x1 - 0.01f0 * x2)
    elseif ispc == 22                              # RA — Harrington PNW-358
        return S + (59.5864f0 + 0.7953f0 * S) * (1f0 - exp((0.00194f0 - 0.00074f0 * S) * A))^0.9198f0 -
               (59.5864f0 + 0.7953f0 * S) * (1f0 - exp((0.00194f0 - 0.00074f0 * S) * 20f0))^0.9198f0
    elseif ispc == 28                              # WO — King (DF Weyerhaeuser Paper 8)
        z = 2500f0 / (S - 4.5f0)
        return (A * A / (-0.954038f0 + 0.109757f0 * z + (5.58178f-2 + 7.92236f-3 * z) * A +
               (-7.33819f-4 + 1.97693f-4 * z) * (A * A))) + 4.5f0
    elseif ispc in (13, 14, 18, 19, 20, 21, 23, 24, 25, 26, 27, 29, 30, 32)  # Curtis
        h = (S - 4.5f0) / (0.6192f0 - 5.3394f0 / (S - 4.5f0) +
            240.29f0 * A^(-1.4f0) + (3368.9f0 / (S - 4.5f0)) * A^(-1.4f0))
        return h + 4.5f0
    else
        return 0f0
    end
end

# ── ec/findag.f + ec/htgf.f DATA (per-species, index = ISPC).
const EC_HT_AGMAX = Float32[200,110,180,130,250,130,140,150,150,200,200,180,200,200,200,130,200,200,200,200,200,200,200,200,200,200,200,200,200,200,180,200]
const EC_HT_AHMAX = Float32[2.3,12.86,-2.86,21.29,52.27,21.29,2.3,20.0,45.27,-5.00,4.0156013,-2.06,3.2412923,3.2412923,4.3149844,21.29,3.2412923,3.2412923,3.2412923,3.9033821,3.9033821,3.9033821,3.9033821,3.9033821,3.9033821,3.9033821,3.9033821,3.9033821,3.9033821,3.9033821,-2.06,3.9033821]
const EC_HT_BHMAX = Float32[2.39,1.32,1.54,1.24,1.14,1.24,1.75,1.1,1.24,1.30,51.9732476,1.54,62.7139427,62.7139427,39.6317079,1.24,62.7139427,62.7139427,62.7139427,59.3370816,59.3370816,59.3370816,59.3370816,59.3370816,59.3370816,59.3370816,59.3370816,59.3370816,59.3370816,59.3370816,1.54,59.3370816]
const EC_HT_RHR   = Float32[15,12,15,20,20,16,12,16,16,13,20,15,20,15,15,16,12,16,13,20,20,13,13,15,20,12,12,15,13,12,15,12]
const EC_HT_RHYXS = Float32[0.10,0.01,0.10,0.20,0.20,0.15,0.01,0.15,0.15,0.05,0.20,0.10,0.20,0.10,0.10,0.15,0.01,0.15,0.05,0.20,0.20,0.05,0.05,0.10,0.20,0.01,0.01,0.10,0.05,0.01,0.10,0.01]
const EC_HT_RHM   = fill(1.10f0, 32)
const EC_HT_RHB   = Float32[-1.45,-1.60,-1.45,-1.10,-1.10,-1.20,-1.60,-1.20,-1.20,-1.60,-1.10,-1.45,-1.10,0.10,-1.45,-1.20,-1.60,-1.20,-1.60,-1.10,-1.10,-1.60,-1.60,-1.45,-1.10,-1.60,-1.60,-1.45,-1.60,-1.60,-1.45,-1.60]
const EC_HT_CRA = 100.0f0; const EC_HT_CRB = 3.0f0; const EC_HT_CRC = -5.0f0
const EC_HT_RHK = 1.0f0; const EC_HT_RHXS = 0.0f0
const EC_HT_WCFORM = Set{Int}([11,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32])  # WC-form findag/htgf branch

# ec/findag.f — invert the species' OWN site curve (ec_htcalc(SINDX, ISPC, AG)) to a stand age.
@inline function ec_findag(ispc::Int, d1::Float32, d2::Float32, h::Float32, sindx::Float32)
    agmax = EC_HT_AGMAX[ispc]
    htmax1 = 0f0; htmax2 = 0f0; ag = 2.0f0
    if ispc in EC_HT_WCFORM                              # WC-form: htmax from AHMAX·D+BHMAX
        htmax1 = EC_HT_AHMAX[ispc]*d1 + EC_HT_BHMAX[ispc]
        htmax2 = EC_HT_AHMAX[ispc]*d2 + EC_HT_BHMAX[ispc]
        ag = 2.0f0
    elseif ispc == 12 || ispc == 31                     # MH/OS: metric site
        htmax1 = EC_HT_AHMAX[ispc] + EC_HT_BHMAX[ispc]*sindx*3.281f0
        ag = 0.5f0
    else                                                # EC-native {1:10,16}
        htmax1 = EC_HT_AHMAX[ispc] + EC_HT_BHMAX[ispc]*sindx
        ag = 0.5f0
        ispc == 10 && (ag = 98.38f0*exp(sindx*(-0.0422f0)) + 1.0f0; ag < 0.5f0 && (ag = 0.5f0))
        ispc == 3  && (ag = 18.0f0)
    end
    if h >= htmax1
        return (agmax + (h - htmax1)/0.10f0, h, agmax, htmax1, htmax2)
    end
    incrng = false; hguess = 0f0
    while true
        oldhg = hguess
        hguess = ec_htcalc(sindx, ispc, ag)
        if hguess >= 1.0f0
            diff = abs(hguess - h)
            (diff <= 2.0f0 || h < hguess) && return (ag, hguess, agmax, htmax1, htmax2)
            d2h = hguess - oldhg
            (oldhg != 0f0 && d2h >= 0.05f0) && (incrng = true)
            (incrng && d2h < 0.05f0) && return (ag, hguess, agmax, htmax1, htmax2)
        end
        ag += 2.0f0
        ag > agmax && return (agmax, h, agmax, htmax1, htmax2)
    end
end

function height_growth!(s::StandState, ::EastCascades; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    sd = s.coef.species
    avh = p.avg_height; ba = p.basal_area
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0.0f0 && continue
        ispc = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]; dg = t.diam_growth[i]
        sindx = p.sp_site_index[ispc]
        icr_pct = Float32(t.crown_pct[i])
        htcon = c.htg_cor[ispc]
        et = Int(sd[:bark_imap][ispc]); brat = wc_bratio(sd[:bark1][ispc], sd[:bark2][ispc], et, d)
        d1 = d; d2 = d + dg/brat
        sitage, sitht, agmax, htmax, htmax2 = ec_findag(ispc, d1, d2, h, sindx)
        wcform = ispc in EC_HT_WCFORM
        # HTMAX gate (ec/htgf.f SELECT before POTHTG)
        if wcform
            if h > htmax
                local htg::Float32 = 0f0
                if h >= htmax2
                    htg = 0.5f0*dg; htg < 0.1f0 && (htg = 0.1f0)
                    htg = scale*htg*exp(htcon)
                end
                cap = s.control.sp_size_cap[ispc, 4]
                (h + htg > cap) && (htg = cap - h; htg < 0.1f0 && (htg = 0.1f0))
                t.ht_growth[i] = htg; continue
            end
        else
            if h >= htmax
                htg = scale*0.1f0*exp(htcon)
                cap = s.control.sp_size_cap[ispc, 4]
                (h + htg > cap) && (htg = cap - h; htg < 0.1f0 && (htg = 0.1f0))
                t.ht_growth[i] = htg; continue
            end
        end
        # POTHTG
        local pothtg::Float32
        if sitage >= agmax
            pothtg = 0.10f0
            ispc == 10 && (pothtg = -1.31f0 + 0.05f0*sindx; pothtg < 0.1f0 && (pothtg = 0.1f0))
        else
            hguess = ec_htcalc(sindx, ispc, sitage + 10.0f0)
            pothtg = hguess - sitht
            if ispc == 28                                # OWO patch
                maxg = sindx - 18.6024f0/log(2.7f0 + ba)
                dd2 = d + dg; dd2 < 0f0 && (dd2 = 0.1f0)
                hg2 = 4.5f0 + maxg*(1f0-exp(-0.137428f0*dd2))^1.38994f0
                hg1 = 4.5f0 + maxg*(1f0-exp(-0.137428f0*d))^1.38994f0
                pothtg = hg2 - hg1
            end
            pothtg < 0.1f0 && (pothtg = 0.1f0)
        end
        # modifiers
        cr = icr_pct/100f0
        hgmdcr = EC_HT_CRA * cr^EC_HT_CRB * exp(EC_HT_CRC*cr); hgmdcr > 1f0 && (hgmdcr = 1f0)
        relht = avh > 0f0 ? h/avh : 0f0; relht > 1.5f0 && (relht = 1.5f0)
        rhx = relht
        fctrkx = (EC_HT_RHK/EC_HT_RHYXS[ispc])^(EC_HT_RHM[ispc]-1f0) - 1f0
        fctrrb = -1f0*(EC_HT_RHR[ispc]/(1f0-EC_HT_RHB[ispc]))
        fctrxb = rhx^(1f0-EC_HT_RHB[ispc]) - EC_HT_RHXS^(1f0-EC_HT_RHB[ispc])
        fctrm  = -1f0/(EC_HT_RHM[ispc]-1f0)
        hgmdrh = EC_HT_RHK * (1f0 + fctrkx*exp(fctrrb*fctrxb))^fctrm
        htgmod = 0.25f0*hgmdcr + 0.75f0*hgmdrh
        htgmod >= 2f0 && (htgmod = 2f0); htgmod <= 0f0 && (htgmod = 0.1f0)
        htg = pothtg*htgmod
        # cap at site max height
        temph = h + htg
        if wcform
            temph > htmax2 && (htg = htmax2 - h)
        else
            temph > htmax && (htg = htmax - h)
        end
        htg < 0.1f0 && (htg = 0.1f0)
        htg = scale*htg*exp(htcon)
        cap = s.control.sp_size_cap[ispc, 4]
        (h + htg > cap) && (htg = cap - h; htg < 0.1f0 && (htg = 0.1f0))
        t.ht_growth[i] = htg
    end
    return s
end
