# =============================================================================
# diameter_growth.jl (bluemountains) — large-tree DG (bm/dgf.f). Chunk 3.
#
#   bm_dgcons!(s)        — per-species per-stand DGCON/DGDSQ/DGCCF/SMCON/ATTEN (bm/dgf.f ENTRY DGCONS).
#   dgf!(s, ::BlueMountains) — per-tree WK2 = DDS (bm/dgf.f main body, 3-group SELECT CASE(ISPC)).
#
# BM is a Wykoff DDS variant with one novel feature: the BM-original species (1:5,7:10,17) use an
# MSS SPLINE — a large-tree DDSL (Wykoff, >=10") blended with a small-tree DDSS (SM* coeffs, fit to
# 5yr then adjusted to 10yr) by XWT=(10-D)/7 over D in [3,10]. Groups (6,11,12,15) reuse UT/TT forms
# (WJ DF-projection, AS DGFASP, WB/LM simple Wykoff); (13,14,16,18) use an extended WC Wykoff.
# CONSPP = DGCON + COR + 0.01*DGCCFA*RELDEN.  BAL = (1 - PCT/100)*BA.
# =============================================================================

# bm/dgf.f ENTRY DGCONS — load site/habitat-dependent DGCON, SMCON, DGDSQ, DGCCF, ATTEN, bark.
function bm_dgcons!(s::StandState)
    c = s.calib; p = s.plot; ctl = s.control; sd = s.coef.species
    isisp = Int(p.site_species); (isisp < 1 || isisp > 18) && (isisp = 10)
    ifor  = Int(p.forest_idx);   (ifor < 1 || ifor > 4) && (ifor = 1)
    elev = p.elevation; aspect = p.aspect; slope = p.slope
    # ICL5 = KODTYP habitat (SMMAPH row). A MISSING/OOR habitat defaults to 79 (bm/habtyp.f:67-69 "DEFAULT
    # CONDITIONS — PA = CWG113" → ITYPE=79), NOT 1. #140: jl defaulted to 1, which applied a real DF SMHAB
    # coefficient (SMMAPH(1,2)=1→SMHAB(2,2)=−0.337) instead of the neutral default (SMMAPH(79,·)=0→SMHAB(1,·)=0),
    # under-shooting DF small-tree DG ~6% on habitat-less FIA stands → self-thin under-kill (the 9:2 skew).
    icl5 = Int(p.habitat_code); (icl5 < 1 || icl5 > 92) && (icl5 = 79)
    xsite = p.sp_site_index[isisp]                                      # XSITE = SITEAR(ISISP)
    lsi = trunc(Int, xsite / 10f0)
    isic = lsi < 2 ? 1 : (lsi >= 5 ? 5 : (lsi >= 4 ? 4 : (lsi >= 3 ? 3 : 2)))
    @inbounds for sp in 1:18
        c.dg_ccf[sp] = BM_DGCCFA[sp]
        asptem = (sp == 6 || sp == 11 || sp == 12 || sp == 15) ? aspect - 0.7854f0 : aspect
        tsite = p.sp_site_index[sp]
        if sp == 7
            tsite = -43.78f0 + 2.16f0 * tsite; tsite <= 1f0 && (tsite = 1f0)
        end
        sp == 5 && (tsite *= 3.28f0)
        temel = elev
        (sp == 16 || sp == 18) && temel > 30f0 && (temel = 30f0)
        dgcon = BM_DGSIC[sp] * xsite + BM_DGFOR[ifor, sp] +
                BM_DGEL[sp] * temel + BM_DGEL2[sp] * temel * temel +
                (BM_DGSASP[sp] * sin(asptem) + BM_DGCASP[sp] * cos(asptem) + BM_DGSLOP[sp]) * slope +
                BM_DGSLSQ[sp] * slope * slope + BM_DGSITE[sp] * log(tsite)
        c.dg_dsq[sp] = BM_DGDS[sp]
        # SMCON (small trees <10") for the MSS-spline species, via the habitat-group SMHAB dimension.
        if sp in (1,2,3,4,5,7,8,9,10,17)
            indxs = BM_SMMAPS[sp]
            if indxs >= 1
                indxh = BM_SMMAPH[icl5, indxs] + 1                     # SMMAPH(ICL5,INDXS)+1
                smcon = BM_SMHAB[indxh, indxs] + BM_SMFOR[ifor, indxs] +
                        BM_SMEL[indxs] * elev + BM_SMEL2[indxs] * elev * elev +
                        (BM_SMSASP[indxs] * sin(asptem) + BM_SMCASP[indxs] * cos(asptem) +
                         BM_SMSLOP[indxs]) * slope
                c.sm_const[sp] = smcon
            else
                c.sm_const[sp] = 0f0
            end
        else
            c.sm_const[sp] = 0f0
        end
        # ATTEN (observed-growth weight): WJ/WB/LM/AS use IBSERV[ISIC,·]; else OBSERV.
        c.atten[sp] = sp == 6 ? Float32(BM_IBSERV(isic, 1)) :
                      (sp == 11 || sp == 12) ? Float32(BM_IBSERV(isic, 2)) :
                      sp == 15 ? Float32(BM_IBSERV(isic, 3)) : BM_OBSERV[sp]
        if ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0
            dgcon += log(ctl.dg_cor2[sp])
            sp in (1,2,3,4,5,7,8,9,10,17) && (c.sm_const[sp] += log(ctl.dg_cor2[sp]))
        end
        c.dg_const[sp] = dgcon
        # BM bark (bm/bratio.f) — POWER model; store BARK1/BARK2 in bark_a/bark_b for reference, but
        # the DBH update + UT/TT DDS call bm_bratio directly (dispatched on the variant).
        c.bark_a[sp] = sd[:bark2][sp]; c.bark_b[sp] = sd[:bark1][sp]
    end
    return s
end

# bm/dgf.f main body — per-tree WK2 = DDS (log inside-bark DDS). bmt01 exercises the MSS-spline group.
function dgf!(s::StandState, ::BlueMountains)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    sd = s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    relden = p.relative_density
    ba = p.basal_area
    alba = ba > 0f0 ? log(ba) : 0f0
    rmai = 50f0                                                 # grinit RMAI default; WP(sp1)-only term, maical.f RMAI deferred (no WP in bmt01)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i])
        ald = log(d)
        cr = Float32(t.crown_pct[i]) * 0.01f0
        pct = t.crown_ratio[i]                                    # PCT = BA percentile
        pt = Int(t.plot_id[i])
        pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
        conspp = c.dg_const[sp] + c.dg_cor[sp] + 0.01f0 * c.dg_ccf[sp] * relden
        dds = -9.21f0
        if sp in (1,2,3,4,5,7,8,9,10,17)
            # BM-original: MSS spline (large DDSL + small DDSS blended by XWT over D in 3-10).
            bal = (1f0 - pct / 100f0) * ba
            ddsl = conspp + BM_DGLD[sp] * ald + BM_DGLBA[sp] * alba +
                   cr * (BM_DGCR[sp] + cr * BM_DGCRSQ[sp]) + c.dg_dsq[sp] * d * d +
                   BM_DGDBAL[sp] * bal / log(d + 1f0) + BM_DGPCCF[sp] * pccf
            sp == 1 && (ddsl += 0.00121f0 * bal + 0.00001f0 * 0.01f0 * rmai * relden - 0.0000016f0 * relden)
            sp == 2 && (ddsl -= 0.000695f0 * ba)
            ddss = 0f0
            if d < 10f0
                indxs = BM_SMMAPS[sp]
                ddss = c.sm_const[sp] + c.dg_cor[sp] + BM_SMLD[indxs] * ald + BM_SMLBA[indxs] * alba +
                       cr * (BM_SMCR[indxs] + cr * BM_SMCRSQ[indxs]) + BM_SMDS[indxs] * d * d +
                       BM_SMDBAL[indxs] * bal / log(d + 1f0) + BM_SMPCCF[indxs] * pccf
                dsq = exp(ddss)                                   # 5yr → 10yr adjust
                temdg = (sqrt(d * d + dsq) - d) * 2f0
                dsqnew = (temdg + d)^2 - d * d
                ddss = dsqnew <= 0f0 ? 0f0 : log(dsqnew)
            end
            xwt = 1f0
            d >= 3f0 && (xwt = (10f0 - d) / 7f0)
            xwt < 0f0 && (xwt = 0f0)
            dds = xwt * ddss + (1f0 - xwt) * ddsl
        elseif sp == 6 || sp == 11 || sp == 12 || sp == 15
            bark = bm_bratio(sd, sp, d)
            if sp == 6
                dpp = d < 1f0 ? 1f0 : d
                batem = ba < 1f0 ? 1f0 : ba
                si = p.sp_site_index[sp]
                df = 0.25897f0 + 1.03129f0 * dpp - 0.0002025464f0 * batem + 0.00177f0 * si
                (df - dpp) > 1f0 && (df = dpp + 1f0)
                df < dpp && (df = dpp)
                diagr = (df - dpp) * bark
                dds = diagr <= 0f0 ? -9.21f0 : log(diagr * (2f0 * dpp * bark + diagr)) + conspp
            elseif sp == 15
                cr_raw = Float32(t.crown_pct[i])
                si = p.sp_site_index[sp]
                rmsqd = _TT_CUR_RMSQD[] >= 0f0 ? _TT_CUR_RMSQD[] : stand_qmd(s)   # #195: current RMSQD during DGSCOR calibration
                aspdg = _em_dgfasp(d, cr_raw, bark, si, rmsqd, ba)
                cor2 = (s.control.dg_cor2_on && s.control.dg_cor2[sp] > 0f0) ? s.control.dg_cor2[sp] : 1f0
                dds = aspdg + log(cor2) + c.dg_cor[sp]
            else
                # WB/LM simple Wykoff, BAL uses BA/100.
                bal = (1f0 - pct / 100f0) * ba / 100f0
                dds = conspp + BM_DGLD[sp] * ald + BM_DGBAL[sp] * bal +
                      BM_DGCR[sp] * cr + cr * cr * BM_DGCRSQ[sp] + c.dg_dsq[sp] * d * d +
                      BM_DGPCCF[sp] * pccf
            end
        elseif sp == 13 || sp == 14 || sp == 16 || sp == 18
            # WC extended Wykoff.
            bal = (1f0 - pct / 100f0) * ba
            dds = conspp + BM_DGLD[sp] * ald + cr * (BM_DGCR[sp] + cr * BM_DGCRSQ[sp]) +
                  c.dg_dsq[sp] * d * d + BM_DGDBAL[sp] * bal / log(d + 1f0) +
                  BM_DGPCCF[sp] * pccf + BM_DGLBA[sp] * log(ba) + BM_DGBAL[sp] * bal + BM_DGBA[sp] * ba
        end
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end
