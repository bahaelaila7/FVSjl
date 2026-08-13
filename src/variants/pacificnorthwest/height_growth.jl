# =============================================================================
# height_growth.jl (pacificnorthwest) — PN large-tree HTG. Chunk 4.
#
# pn/htgf.f + pn/findag.f are LITERALLY the WC files (WC $Id$) ⇒ the FINDAG age-finder, the RHR/RHYXS/
# RHM/RHB + CRA/CRB/CRC height-growth modifiers, the OWO(28) King patch, and the RW(17) LTHTG special are
# byte-identical — PN reuses the shared WC_HT_* consts. ONLY pn/htcalc.f (the site-species height curve)
# differs: SS(6) & RC(18) use a NEW Farr PNW-326 curve; DF(16) uses King (z=2500, not the misc Curtis).
# So pn_htcalc is the only new function; pn_findag/pn_htg_default/height_growth! mirror WC's but call it.
# (WC HTG is validated 24/24 bit-exact; PN's identical modifiers + the DGCON-validated site fan-out give
#  the same guarantee once site/density run — validated at the pnt01 end-to-end stage.)
# =============================================================================

# pn/htcalc.f — site-species height curve. Identical to wc_htcalc EXCEPT ispc 6/18 (Farr) and 16 (King).
@inline function pn_htcalc(sindx::Float32, ispc::Int, ag::Float32)
    if ispc == 6 || ispc == 18                               # SS, RC — Farr PNW-326 (NEW vs WC)
        la = log(ag)
        az = -0.2050542f0 + 1.449615f0*la - 0.01780992f0*la^3 + 6.519748f-5*la^5 - 1.095593f-23*la^30
        bz = -5.611879f0 + 2.418604f0*la - 0.2593110f0*la^2 + 1.351445f-4*la^5 -
             1.701139f-12*la^16 + 7.964197f-27*la^36
        return 4.5f0 + exp(az) - 86.43f0*exp(bz) + (sindx - 4.5f0)*exp(bz)
    elseif ispc == 16 || ispc == 28                          # DF (King, ≠ WC misc), WO — King z=2500
        z = 2500f0/(sindx-4.5f0)
        return ag*ag/(-0.954038f0 + 0.109757f0*z + (5.58178f-2+7.92236f-3*z)*ag +
               (-7.33819f-4+1.97693f-4*z)*(ag*ag)) + 4.5f0
    elseif ispc == 1                                         # SF — Hoyer PNW-418
        sm45 = sindx - 4.5f0; k = 0.0071839f0 + 0.0000571f0 * sm45
        return sm45 * (1f0 - exp(-k*ag))^1.39005f0 / (1f0 - exp(-k*100f0))^1.39005f0 + 4.5f0
    elseif ispc == 2 || ispc == 3                            # WF, GF — Cochran PNW-252
        la = log(ag)
        x2 = -0.30935f0 + 1.2383f0*la + 0.001762f0*la^4 - 5.4f-6*la^9 + 2.046f-7*la^11 - 4.04f-13*la^18
        x3 = -6.2056f0 + 2.097f0*la - 0.09411f0*la^2 - 4.382f-5*la^7 + 2.007f-11*la^16 - 2.054f-17*la^24
        return exp(x2) - 84.93f0*exp(x3) + (sindx-4.5f0)*exp(x3) + 4.5f0
    elseif ispc == 4 || ispc == 10                           # AF, ES — Alexander RM-32
        return 4.5f0 + (2.75780f0*sindx^0.83312f0) *
               (1f0 - exp(-0.015701f0*ag))^(22.71944f0*sindx^(-0.63557f0))
    elseif ispc == 5                                         # RF — Dolph PSW-206
        term  = ag*exp(ag*(-0.0440853f0))*1.4151f-6
        b     = sindx*term - 3.0495f6*term*term + 5.72474f-4
        term2 = 50f0*exp(50f0*(-0.0440853f0))*1.4151f-6
        b50   = sindx*term2 - 3.0495f6*term2*term2 + 5.72474f-4
        return (sindx-4.5f0)*(1f0-exp(-b*(ag^1.51744f0))) / (1f0-exp(-b50*(50f0^1.51744f0))) + 4.5f0
    elseif ispc == 7                                         # NF — Herman PNW-243
        s45 = sindx - 4.5f0
        x1 = -564.38f0 + 22.25f0*s45 - 0.04995f0*s45^2
        x2 = 6.80f0 + 2843.21f0*s45^(-1) + 34735.54f0*s45^(-2)
        return 4.5f0 + s45 / (x1*(1f0/ag)^2 + x2*(1f0/ag) + 1f0 - 0.0001f0*x1 - 0.01f0*x2)
    elseif ispc == 9 || ispc == 12 || ispc == 15             # IC, JP, PP — Barrett PNW-232
        t = (-0.7864f0 + 2.49717f0*(1f0-exp(-0.0045042f0*ag))^0.33022f0)
        return 128.8952205f0*(1f0-exp(-0.016959f0*ag))^1.23114f0 - t*100.43f0 + t*(sindx-4.5f0) + 4.5f0
    elseif ispc == 11                                        # LP — Dahms PNW-8
        return sindx*(-0.0968f0 + 0.02679f0*ag - 0.00009309f0*ag*ag)
    elseif ispc == 13 || ispc == 14                          # SP, WP — Curtis PNW-423
        num = 1f0 - exp(-exp(-4.62536f0 + 1.346399f0*log(ag) - 135.354483f0/sindx))
        den = 1f0 - exp(-exp(-4.62536f0 + 1.346399f0*log(100f0) - 135.354483f0/sindx))
        return num/den * (sindx-4.5f0) + 4.5f0
    elseif ispc == 19                                        # WH — Wiley 1978
        z = 2500f0/(sindx-4.5f0)
        return ag*ag/(-1.7307f0 + 0.1394f0*z + (-0.0616f0+0.0137f0*z)*ag +
               (0.00192f0+0.00007f0*z)*(ag*ag)) + 4.5f0
    elseif ispc == 20                                        # MH — Means
        hh = (22.8741f0 + 0.950234f0*sindx)*(1f0-exp(-0.00206465f0*sqrt(sindx)*ag))^
             (1.365566f0 + 2.045963f0/sindx)
        return (hh + 1.37f0)*3.281f0
    elseif ispc == 22                                        # RA — Harrington PNW-358
        cc = (59.5864f0 + 0.7953f0*sindx); k = (0.00194f0 - 0.0007403f0*sindx)
        return sindx + cc*(1f0-exp(k*ag))^0.9198f0 - cc*(1f0-exp(k*20f0))^0.9198f0
    elseif ispc == 30                                        # LL — Cochran PNW-424
        t = (-0.12528f0 + 0.039636f0*ag - 4.278f-4*ag*ag + 1.7039f-6*ag^3)
        return 4.5f0 + 1.46897f0*ag + 0.0092466f0*ag*ag - 2.3957f-4*ag^3 + 1.1122f-6*ag^4 +
               (sindx-4.5f0)*t - 73.57f0*t
    elseif ispc == 8 || ispc == 17 || ispc == 21 || (23 <= ispc <= 27) ||
           ispc == 29 || (31 <= ispc <= 37) || ispc == 39    # "misc" — Curtis (16/18 removed vs WC)
        s45 = sindx - 4.5f0
        return s45 / (0.6192f0 - 5.3394f0/s45 + 240.29f0*ag^(-1.4f0) + (3368.9f0/s45)*ag^(-1.4f0)) + 4.5f0
    else                                                     # blank/future (38)
        return 0.0f0
    end
end

# pn/findag.f — byte-identical to wc_findag except it calls pn_htcalc. Reuses the WC_HT_* consts.
@inline function pn_findag(ispc::Int, d::Float32, d2::Float32, h::Float32, sindx::Float32)
    agmax = WC_HT_AGMAX; mp = WC_HT_MAPHD[ispc]
    htmax  = WC_HT_HDRAT1[mp]*d  + WC_HT_HDRAT2[mp]
    htmax2 = WC_HT_HDRAT1[mp]*d2 + WC_HT_HDRAT2[mp]
    h >= htmax && return (agmax + (h - htmax)/0.10f0, h, agmax, htmax, htmax2)
    ag = 2.0f0; hguess = 0.0f0
    while true
        oldhg = hguess; incrng = false
        hguess = pn_htcalc(sindx, ispc, ag)
        if hguess >= 1.0f0
            diff = abs(hguess - h)
            (diff <= 2.0f0 || h < hguess) && return (ag, hguess, agmax, htmax, htmax2)
            d2h = hguess - oldhg
            (oldhg != 0.0f0 && d2h >= 0.05f0) && (incrng = true)
            (incrng && d2h < 0.05f0) && return (ag, hguess, agmax, htmax, htmax2)
        end
        ag += 2.0f0
        ag > agmax && return (agmax, h, agmax, htmax, htmax2)
    end
end

# pn/htgf.f DEFAULT + OWO(28) — byte-identical to wc_htg_default except it calls pn_htcalc.
@inline function pn_htg_default(ispc::Int, sindx::Float32, d::Float32, h::Float32,
                                icr_pct::Float32, avh::Float32, ba::Float32, dg::Float32,
                                d2::Float32, sitage::Float32, sitht::Float32, agmax::Float32, htmax2::Float32)
    local pothtg::Float32
    if sitage >= agmax
        pothtg = (ispc == 2 || ispc == 3) ? (0.2f0 + 0.00264f0*sindx)*10f0 : 0.10f0
    else
        hguess = pn_htcalc(sindx, ispc, sitage + 10.0f0)
        pothtg = hguess - sitht
        if ispc == 28
            maxg = sindx - 18.6024f0/log(2.7f0 + ba)
            dd2 = d + dg; dd2 < 0f0 && (dd2 = 0.1f0)
            hg2 = 4.5f0 + maxg*(1f0-exp(-0.137428f0*dd2))^1.38994f0
            hg1 = 4.5f0 + maxg*(1f0-exp(-0.137428f0*d))^1.38994f0
            pothtg = hg2 - hg1
        end
        pothtg < 0.1f0 && (pothtg = 0.1f0)
    end
    cr = icr_pct/100f0
    hgmdcr = WC_HT_CRA * cr^WC_HT_CRB * exp(WC_HT_CRC*cr); hgmdcr > 1f0 && (hgmdcr = 1f0)
    relht = avh > 0f0 ? h/avh : 0f0; relht > 1.5f0 && (relht = 1.5f0)
    rhx = relht
    fctrkx = (WC_HT_RHK/WC_HT_RHYXS[ispc])^(WC_HT_RHM[ispc]-1f0) - 1f0
    fctrrb = -1f0*(WC_HT_RHR[ispc]/(1f0-WC_HT_RHB[ispc]))
    fctrxb = rhx^(1f0-WC_HT_RHB[ispc]) - WC_HT_RHXS^(1f0-WC_HT_RHB[ispc])
    fctrm  = -1f0/(WC_HT_RHM[ispc]-1f0)
    hgmdrh = WC_HT_RHK * (1f0 + fctrkx*exp(fctrrb*fctrxb))^fctrm
    htgmod = 0.25f0*hgmdcr + 0.75f0*hgmdrh
    htgmod >= 2f0 && (htgmod = 2f0); htgmod <= 0f0 && (htgmod = 0.1f0)
    htg = pothtg*htgmod
    (h + htg > htmax2) && (htg = htmax2 - h)
    htg < 0.1f0 && (htg = 0.1f0)
    return htg
end

function height_growth!(s::StandState, ::PacificNorthwest; scale::Float32 = 1.0f0)
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
        d2 = d + dg/brat
        local htg::Float32
        if ispc == 17                                        # REDWOOD LTHTG (identical to WC)
            dg10 = h < 4.5f0 ? 0.1f0 : dg/brat
            lthtg = exp(1.412947f0 - 0.000204f0*d*d + 0.31971f0*log(d) +
                        0.394005f0*log(sindx) + 0.399888f0*log(dg10) - 0.451708f0*log(h))
            hgbnd = if h >= 217.0f0 && h < 380.0f0
                        max(1.0f0 - (h-217.0f0)/(380.0f0-217.0f0), 0.1f0)
                    elseif h < 217.0f0; 1.0f0 else 0.1f0 end
            htg = lthtg*hgbnd; htg < 0.1f0 && (htg = 0.1f0)
            htg = scale*htg*exp(htcon)
        else
            sitage, sitht, agmax, htmax, htmax2 = pn_findag(ispc, d, d2, h, sindx)
            if h > htmax
                if h >= htmax2
                    htg = 0.5f0*dg; htg < 0.1f0 && (htg = 0.1f0)
                    htg = scale*htg*exp(htcon)
                else
                    htg = 0.0f0
                end
            else
                htg = pn_htg_default(ispc, sindx, d, h, icr_pct, avh, ba, dg, d2, sitage, sitht, agmax, htmax2)
                htg = scale*htg*exp(htcon)
            end
        end
        cap = s.control.sp_size_cap[ispc, 4]
        if h + htg > cap
            htg = cap - h; htg < 0.1f0 && (htg = 0.1f0)
        end
        t.ht_growth[i] = htg
    end
    return s
end
