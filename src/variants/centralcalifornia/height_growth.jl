# =============================================================================
# height_growth.jl (centralcalifornia) — CA large-tree height growth (ca/htgf.f + findag.f + htcalc.f). Chunk 4b.
#
# POTHTG = HGUESS(SITAGE+10) − SITHT, where FINDAG inverts the site curve (ca/htcalc.f) to the age SITAGE at
# which predicted height == the tree's H, then HTG = POTHTG · XMOD (Ritchie-Hann vigor·competition modifier).
#   XMOD = 1.016605·(1−exp(−4.26558·CR))·exp(2.54119·(RELHT^0.250537 − 1)); RELHT = min(H/AVH,1), =1 if PCCF<100.
# HTCALC: R5 forests (JFOR≤5) = Dunning/Levitan; R6 = 8 species-group curves (Hann-Scrivani A/B, Dolph red fir,
# Dahms lodgepole, Powers black oak, Porter-Wiant tanoak/madrone/red-alder). RW/GS(23,50) use a separate
# LN(HI) increment + a 217–380 ft bounding function (no RW/GS in cat01, ported source-faithful).
# =============================================================================

const CA_DUNL1 = Float32[-88.9, -82.2, -78.3, -82.1, -56.0, -33.8]
const CA_DUNL2 = Float32[49.7067, 44.1147, 39.1441, 35.4160, 26.7173, 18.6400]
const CA_DUNL3 = Float32[2.375, 2.025, 1.650, 1.225, 1.075, 0.875]

# ca/htcalc.f — predicted total height (ft) at age `ag` for species `kspec`, site index `sindx`, forest `jfor`.
function ca_htcalc(jfor::Int, sindx::Float32, kspec::Int, ag::Float32)::Float32
    if jfor <= 5                                    # R5 Dunning/Levitan
        indx = sindx <= 44f0 ? 6 : sindx <= 52f0 ? 5 : sindx <= 65f0 ? 4 :
               sindx <= 82f0 ? 3 : sindx <= 98f0 ? 2 : 1
        return ag <= 40f0 ? CA_DUNL3[indx] * ag : CA_DUNL1[indx] + CA_DUNL2[indx] * log(ag)
    end
    # R6 species-group site curves
    if kspec in (1,2,3,4,7,8,16,22,23,25,50)        # Hann-Scrivani A (PC/IC/RC/WF/DF/WH/SP/BR/GS/OS/RW)
        top = 1f0 - exp(-exp(-6.21693f0 + 0.281176f0*log(sindx-4.5f0) + 1.14354f0*log(ag)))
        bot = 1f0 - exp(-exp(-6.21693f0 + 0.281176f0*log(sindx-4.5f0) + 1.14354f0*log(50f0)))
        return (((sindx-4.5f0)*top/bot) + 4.5f0) * 1.05f0
    elseif kspec in (15,17,18,19,20)                # Hann-Scrivani B (JP/WP/PP/MP/GP)
        top = 1f0 - exp(-exp(-6.54707f0 + 0.288169f0*log(sindx-4.5f0) + 1.21297f0*log(ag)))
        bot = 1f0 - exp(-exp(-6.54707f0 + 0.288169f0*log(sindx-4.5f0) + 1.21297f0*log(50f0)))
        return (((sindx-4.5f0)*top/bot) + 4.5f0) * 1.05f0
    elseif kspec in (5,6,9)                          # Dolph red fir (RF/SH/MH)
        term = ag*exp(ag*(-0.0440853f0))*1.41512f-6
        b = sindx*term - 3.04951f6*term*term + 5.72474f-4
        term2 = 50f0*exp(50f0*(-0.0440853f0))*1.41512f-6
        b50 = sindx*term2 - 3.04951f6*term2*term2 + 5.72474f-4
        return ((sindx-4.5f0)*(1f0-exp(-b*(ag^1.51744f0)))) / (1f0-exp(-b50*(50f0^1.51744f0))) + 4.5f0
    elseif kspec in (10,11,12,13,14,21)             # Dahms lodgepole (WB/KP/LP/CP/LM/WJ)
        return sindx*(-0.0968f0 + 0.02679f0*ag - 0.00009309f0*ag*ag) * 1.10f0
    elseif kspec in (24,30,31,32,35,39,40)          # Powers black oak (PY/WO/BO/VO/BU/DG/FL)
        term = sqrt(ag) - sqrt(50f0)
        return ((sindx*(1f0 + 0.322f0*term)) - 6.413f0*term) * 0.70f0
    elseif kspec == 42                              # Porter-Wiant tanoak (TO)
        return (sindx/(0.204f0 + 39.787f0/ag)) * 0.80f0
    elseif kspec in (26,27,28,29,33,37,38)          # Porter-Wiant madrone (LO/CY/BL/EO/IO/MA/GC)
        return (sindx/(0.375f0 + 31.233f0/ag)) * 0.80f0
    elseif kspec in (34,36,41,43,44,45,46,47,48,49) # Porter-Wiant red alder (BM/RA/WN/SY/AS/CW/WI/CN/CL/OH)
        return (sindx/(0.649f0 + 17.556f0/ag)) * 0.80f0
    end
    return 0f0
end

# ca/findag.f — invert the site curve: find SITAGE where ca_htcalc==H. Returns (sitage, sitht, agmax).
function ca_findag(ispc::Int, h::Float32, sindx::Float32, ifor::Int)
    toler = 2f0
    agmax = ifor <= 5 ? 400f0 : 200f0               # AGEMAX=200 (R5=400)
    ag = 2f0; incrng = 0; hguess = 0f0
    while true
        oldhg = hguess
        hguess = ca_htcalc(ifor, sindx, ispc, ag)
        if hguess >= 1f0
            if abs(hguess - h) <= toler || h < hguess
                return (ag, hguess, agmax)
            end
            d = hguess - oldhg
            (oldhg != 0f0 && d >= 0.05f0) && (incrng = 1)
            if incrng == 1 && d < 0.05f0
                return (ag, hguess, agmax)
            end
        end
        ag += 2f0
        ag > agmax && return (agmax, h, agmax)
    end
end

function height_growth!(s::StandState, ::CentralCalifornia; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    sd = s.coef.species
    avh = p.avg_height
    ifor = Int(p.forest_idx)
    pccf = s.density.point_ccf
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        ispc = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]; dglt = t.diam_growth[i]
        sindx = p.sp_site_index[ispc]
        htcon = c.htg_cor[ispc]
        local htg::Float32
        if ispc == 23 || ispc == 50                 # RW / GS (ca/htgf.f LN(HI) + bounding)
            brat = wc_bratio(sd[:bark1][ispc], sd[:bark2][ispc], Int(sd[:bark_imap][ispc]), d)
            dg10 = dglt / brat; h < 4.5f0 && (dg10 = 0.1f0)
            lthtg = exp(1.412947f0 - 0.000204f0*d^2 + 0.31971f0*log(d) + 0.394005f0*log(sindx) +
                        0.399888f0*log(dg10) - 0.451708f0*log(h))
            hgbnd = (h >= 217f0 && h < 380f0) ? max(1f0 - ((h-217f0)/(380f0-217f0)), 0.1f0) :
                    (h < 217f0 ? 1f0 : 0.1f0)
            htg = lthtg * hgbnd
        else                                        # all other species
            sitage, sitht, agmax = ca_findag(ispc, h, sindx, ifor)
            local pothtg::Float32
            if sitage > agmax
                pothtg = 0.10f0
            else
                hguess = ca_htcalc(ifor, sindx, ispc, sitage + 10f0)
                pothtg = hguess - sitht
            end
            pt = Int(t.plot_id[i])
            ppcf = (1 <= pt <= length(pccf)) ? pccf[pt] : 0f0
            relht = avh > 0f0 ? h/avh : 0f0
            relht > 1f0 && (relht = 1f0)
            ppcf < 100f0 && (relht = 1f0)
            cratio = Float32(t.crown_pct[i]) / 100f0
            crmod = 1f0 - exp(-4.26558f0*cratio)
            rhmod = exp(2.54119f0*(relht^0.250537f0 - 1f0))
            xmod = 1.016605f0 * crmod * rhmod
            htg = pothtg * xmod
        end
        htg < 0.1f0 && (htg = 0.1f0)
        htg = scale * htg * exp(htcon)              # XHT=1 (no HTGMULT), MISHGF=1 (no DM)
        cap = s.control.sp_size_cap[ispc, 4]
        (h + htg > cap) && (htg = cap - h; htg < 0.1f0 && (htg = 0.1f0))
        t.ht_growth[i] = htg
    end
    return s
end
