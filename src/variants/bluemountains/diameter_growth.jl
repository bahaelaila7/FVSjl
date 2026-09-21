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
                (BM_DGSASP[sp] * fsin(asptem) + BM_DGCASP[sp] * fcos(asptem) + BM_DGSLOP[sp]) * slope +
                BM_DGSLSQ[sp] * slope * slope + BM_DGSITE[sp] * flog(tsite)   # dgf.f:600-609 SIN/COS/ALOG → glibc sinf/cosf/logf
        c.dg_dsq[sp] = BM_DGDS[sp]
        # SMCON (small trees <10") for the MSS-spline species, via the habitat-group SMHAB dimension.
        if sp in (1,2,3,4,5,7,8,9,10,17)
            indxs = BM_SMMAPS[sp]
            if indxs >= 1
                indxh = BM_SMMAPH[icl5, indxs] + 1                     # SMMAPH(ICL5,INDXS)+1
                smcon = BM_SMHAB[indxh, indxs] + BM_SMFOR[ifor, indxs] +
                        BM_SMEL[indxs] * elev + BM_SMEL2[indxs] * elev * elev +
                        (BM_SMSASP[indxs] * fsin(asptem) + BM_SMCASP[indxs] * fcos(asptem) +
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
            dgcon += flog(ctl.dg_cor2[sp])                               # dgf.f:660 ALOG(COR2)
            sp in (1,2,3,4,5,7,8,9,10,17) && (c.sm_const[sp] += flog(ctl.dg_cor2[sp]))   # dgf.f:664
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
    alba = ba > 0f0 ? flog(ba) : 0f0                              # dgf.f:370 ALOG(BA)
    rmai = bm_rmai(p)                                           # maical.f RMAI=ADJMAI(ISPNUM(ISISP),SITEAR(ISISP),10) ≤128 (CRATET) — WP(sp1) DDSL term (dgf.f:417)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i])
        ald = flog(d)                                             # dgf.f ALOG(D)
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
                   BM_DGDBAL[sp] * bal / flog(d + 1f0) + BM_DGPCCF[sp] * pccf
            sp == 1 && (ddsl = ddsl + 0.00121f0 * bal + 0.00001f0 * 0.01f0 * rmai * relden - 0.0000016f0 * relden)   # dgf.f:412 left-assoc ((DDSL+a)+b)-c
            sp == 2 && (ddsl -= 0.000695f0 * ba)
            ddss = 0f0
            if d < 10f0
                indxs = BM_SMMAPS[sp]
                ddss = c.sm_const[sp] + c.dg_cor[sp] + BM_SMLD[indxs] * ald + BM_SMLBA[indxs] * alba +
                       cr * (BM_SMCR[indxs] + cr * BM_SMCRSQ[indxs]) + BM_SMDS[indxs] * d * d +
                       BM_SMDBAL[indxs] * bal / flog(d + 1f0) + BM_SMPCCF[indxs] * pccf
                dsq = fexp(ddss)                                  # 5yr → 10yr adjust (dgf.f:449)
                temdg = (sqrt(d * d + dsq) - d) * 2f0
                dsqnew = fpow(temdg + d, 2f0) - d * d              # dgf.f:451 (TEMDG+D)**2. = powf, not x*x
                ddss = dsqnew <= 0f0 ? 0f0 : flog(dsqnew)
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
                dds = diagr <= 0f0 ? -9.21f0 : flog(diagr * (2f0 * dpp * bark + diagr)) + conspp   # dgf.f:495
            elseif sp == 15
                cr_raw = Float32(t.crown_pct[i])
                si = p.sp_site_index[sp]
                rmsqd = _TT_CUR_RMSQD[] >= 0f0 ? _TT_CUR_RMSQD[] : stand_qmd(s)   # #195: current RMSQD during DGSCOR calibration
                aspdg = _em_dgfasp(d, cr_raw, bark, si, rmsqd, ba)
                cor2 = (s.control.dg_cor2_on && s.control.dg_cor2[sp] > 0f0) ? s.control.dg_cor2[sp] : 1f0
                dds = aspdg + flog(cor2) + c.dg_cor[sp]                   # dgf.f:504
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
                  c.dg_dsq[sp] * d * d + BM_DGDBAL[sp] * bal / flog(d + 1f0) +
                  BM_DGPCCF[sp] * pccf + BM_DGLBA[sp] * flog(ba) + BM_DGBAL[sp] * bal + BM_DGBA[sp] * ba
        end
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end

"""
    bm_cycle0_dg(s) -> Vector{Float32}

The DG(I) array FVS holds at the inventory row (fvs.f:301, after the LSTART calibration pass): bm/dgdriv.f DO 220
(:746-769) — HT≤4.5 ⇒ 0; a measured increment (DG>0, HT>4.5) is kept, capped at the inside-bark DBH when IDG<2;
every other record gets the calibration dub DG = SQRT(D²+EXP(WK2+OLDRN)·SCALE)−D with D = WK3·BARK (WK3 = the
backdated start-of-period DBH, BARK = BRATIO at the CURRENT DBH), SCALE=FINT/YR, capped at D, then DGBND. WK2 is
the :735 DGF(WK3) prediction (final COR, backdated density) stashed by calibrate_diameter_growth! in
`calib.dub_wk2/dub_wk3`; when no calibration stash matches the record set (no calibration pass), WK3 = DBH and
WK2 is re-predicted on the current stand with COR at 0 elapsed years. COVER's inventory-row foliage biomass
(CVCBMS DDS=(2·D·DG+DG²)/FINT) reads this array; jl's `diam_growth` holds only the input increments there.
Side-effect free: the sort tables, scratch and COR it touches are restored.
"""
function bm_cycle0_dg(s::StandState)
    t, c = s.trees, s.calib
    sd = s.coef.species
    n = t.n
    isct_s = copy(s.control.sp_count_tab); ind1_s = copy(s.scratch.idx1)
    wk2_s = s.scratch.wk[2, :]; cor_s = copy(c.dg_cor)
    @inbounds for sp in eachindex(c.dg_cor)
        c.dg_cor[sp] = c.dg_cor_goal[sp] + c.dg_cor_goal[sp]      # diameter_growth!'s COR at 0 elapsed years (cormlt=1)
    end
    species_sort!(s)
    dgf!(s, s.variant)
    wk2 = view(s.scratch.wk, 2, :)
    dlo_v = haskey(sd, :dg_bound_dbh_lo) ? sd[:dg_bound_dbh_lo] : nothing
    dhi_v = haskey(sd, :dg_bound_dbh_hi) ? sd[:dg_bound_dbh_hi] : nothing
    sc = s.control.growth_fint / 10f0
    out = zeros(Float32, n)
    stash = length(c.dub_wk2) == n && length(c.dub_wk3) == n     # dgdriv.f:735 WK2/WK3 from the calibration
    @inbounds for i in 1:n
        t.height[i] <= 4.5f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        bark = bm_bratio(sd, sp, d)
        dib = d * bark
        if t.diam_growth[i] > 0f0
            dg = t.diam_growth[i]
            (s.control.growth_idg < 2 && dg > dib) && (dg = dib)
            out[i] = dg
        else
            dd = stash ? c.dub_wk3[i] * bark : dib                  # D = WK3(I)*BARK (dgdriv.f:749)
            w2 = stash ? c.dub_wk2[i] : wk2[i]
            dub = sqrt(dd * dd + fexp(w2 + t.old_random[i]) * sc) - dd   # EXP → gfortran expf
            dub > dd && (dub = dd)
            out[i] = dg_bound(dlo_v, dhi_v, sp, d, dub, s.control.sp_size_cap)
        end
    end
    s.control.sp_count_tab .= isct_s; s.scratch.idx1 .= ind1_s
    s.scratch.wk[2, :] .= wk2_s; c.dg_cor .= cor_s
    return out
end
