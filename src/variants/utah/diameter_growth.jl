# =============================================================================
# diameter_growth.jl (utah) — UT large-tree DDS (ut/dgf.f). Chunk 3.
#
# Hooks (EM/KT pattern):
#   ut_dgcons!(s)        — per-species per-stand DGCON/DGDSQ/DGCCF/ATTEN (ut/dgf.f ENTRY DGCONS).
#   dgf!(s, ::Utah)      — per-tree WK2 = DDS (ut/dgf.f main body, 5-branch SELECT CASE(ISPC)).
#
# DGCON = DGSIC[IDGSIM[ISMAP[ISISP],sp],sp]·XSITE + DGFOR[ISPFOR,sp]
#       + (DGSASP·sinθ + DGCASP·cosθ + DGSLOP)·SLOPE + DGSLSQ·SLOPE² + DGEL·EL + DGEL2·EL²
#   ISPFOR=MAPLOC[IFOR,sp]; DGDSQ=DGDS[MAPDSQ[IFOR,sp],sp]; θ=ASPECT−0.7854 (raw for MC/BI); XSITE=SITEAR[ISISP]
#   (ln for MC/BI); ATTEN=IBSERV[ISIC,sp] (OBSERV for MC/BI). CONSPP=DGCON+COR+0.01·DGCCF·RELDEN.
#
# DDS branches (ut/dgf.f DO 20):
#   CASE(1:5,7:10,23) conifers — CONSPP + DGLD·lnD + DGBAL·BAL + CR·DGCR + CR²·DGCRSQ + DGDSQ·D² + DGPCCF·PCCF;
#     CR=crown·0.01, BAL=(1−PCT/100)·BA100 (BA100=BA/100); BS(5) result ×0.95.
#   CASE(6) aspen — DGFASP (Utah): DDS = ASPDG + ln(COR2) + COR.  [reuses the EM/IE/TT _em_dgfasp form]
#   CASE(11:16,24) PJ/woodland — DF-projection: DF=0.25897+1.03129·DPP−0.0002025464·BATEM+0.00177·SI (BATEM≥1,
#     cap DF−DPP≤1), DIAGR=(DF−DPP)·BARK, DDS=ln(DIAGR·(2·DPP·BARK+DIAGR))+CONSPP.
#   CASE(20:21) MC/BI — CONSPP + DGLD·lnD + DGBAL·BAL + CR·(DGCR+CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1)
#     + DGPCCF·PCCF + DGBA·BA; BAL=(1−PCT/100)·BA (NOT BA100).
#   CASE(17:19,22) GB/NC/FC/BE (CR-surrogate) — BAU-based DF-proj + DSTAG.  [needs BAU array — DEFERRED; not in utt01]
# =============================================================================

# ut/dgf.f ENTRY DGCONS — per-species per-stand DG constants. Needs ISISP (site_species), IFOR (forest_idx),
# ITYPE (habitat_input), ELEV/ASPECT/SLOPE, SITEAR. Fills calib dg_const/dg_dsq/dg_ccf/atten + bark.
function ut_dgcons!(s::StandState)
    c = s.calib; p = s.plot; ctl = s.control; sd = s.coef.species
    isisp = Int(p.site_species); (isisp < 1 || isisp > 24) && (isisp = 7)
    ifor  = Int(p.forest_idx);   (ifor < 1 || ifor > 6) && (ifor = 1)
    elev = p.elevation; aspect = p.aspect; slope = p.slope
    site_sp = p.sp_site_index[isisp]                    # SITEAR(ISISP)
    lsi = trunc(Int, site_sp / 10f0)
    isic = lsi < 2 ? 1 : (lsi >= 5 ? 5 : (lsi >= 4 ? 4 : (lsi >= 3 ? 3 : 2)))
    isi = UT_ISMAP[isisp]; (isi < 1 || isi > 7) && (isi = 1)   # ISMAP(ISISP) — site-base-species row of IDGSIM
    @inbounds for sp in 1:24
        temel = elev
        xsite = site_sp
        asptem = aspect - 0.7854f0
        if sp == 20
            temel > 30f0 && (temel = 30f0); xsite = log(xsite); asptem = aspect
        elseif sp == 21
            xsite = log(xsite); asptem = aspect
        end
        ispfor = UT_MAPLOC[ifor, sp]; (ispfor < 1 || ispfor > 5) && (ispfor = 1)
        ispdsq = UT_MAPDSQ[ifor, sp]; (ispdsq < 1 || ispdsq > 4) && (ispdsq = 1)
        ksic = UT_IDGSIM[isi, sp]; (ksic < 1 || ksic > 5) && (ksic = 1)
        dgcon = UT_DGSIC[ksic, sp] * xsite + UT_DGFOR[ispfor, sp] +
                (UT_DGSASP[sp] * sin(asptem) + UT_DGCASP[sp] * cos(asptem) + UT_DGSLOP[sp]) * slope +
                UT_DGSLSQ[sp] * slope * slope + UT_DGEL[sp] * temel + UT_DGEL2[sp] * temel * temel
        c.dg_dsq[sp] = UT_DGDS[ispdsq, sp]
        c.dg_ccf[sp] = UT_DGCCFA[sp]
        c.atten[sp] = (sp == 20 || sp == 21) ? UT_OBSERV[sp] : Float32(UT_IBSERV[isic, sp])
        (ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0) && (dgcon += log(ctl.dg_cor2[sp]))
        c.dg_const[sp] = dgcon
        # UT bark (ut/bratio.f) for the DBH update: IMAP=2→BARK1 const (a=0,b=BARK1); IMAP=3→BARK1+BARK2/D
        # (a=BARK2,b=BARK1); IMAP=1 zero-coef→0.9002−0.3089/D. bark_ratio(a,b)=b+a/d.
        b1 = sd[:bark1][sp]; b2 = sd[:bark2][sp]; imapb = round(Int, sd[:bark_imap][sp])
        if imapb == 1 && b1 == 0f0 && b2 == 0f0
            c.bark_a[sp] = -0.3089f0; c.bark_b[sp] = 0.9002f0
        else
            c.bark_a[sp] = b2; c.bark_b[sp] = b1
        end
    end
    return s
end

# ut/dgf.f main body — per-tree WK2 = DDS (outside-bark). utt01 exercises CASE 1 (conifer) + CASE 6 (aspen).
function dgf!(s::StandState, ::Utah)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    wk2 = view(s.scratch.wk, 2, :)
    relden = p.relative_density
    ba = p.basal_area
    ba100 = ba / 100f0
    isisp = Int(p.site_species); (isisp < 1 || isisp > 24) && (isisp = 7)
    # ISC (site class 1-5) for the ES DGCCF override (ut/dgf.f:539).
    isi_c = round(Int, p.sp_site_index[isisp]); isi_c = clamp(isi_c, 20, 60)
    isc = trunc(Int, isi_c / 10f0) - 1
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i])
        ald = log(d)
        cr  = Float32(t.crown_pct[i]) * 0.01f0
        pct = t.crown_ratio[i]                       # PCT = BA percentile (FVSjl convention)
        dgccf = c.dg_ccf[sp]
        (isc >= 4 && sp == 8) && (dgccf = -0.154870f0)     # ut/dgf.f:539 ES site-class override
        conspp = c.dg_const[sp] + c.dg_cor[sp] + 0.01f0 * dgccf * relden
        if sp in (1,2,3,4,5,7,8,9,10,23)
            # Conifer Wykoff (ut/dgf.f:587-591). BAL uses BA100; DGPTCC·PCCF (point CCF).
            bal = (1f0 - pct / 100f0) * ba100
            pt = Int(t.plot_id[i])
            pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
            dds = conspp + UT_DGLD[sp] * ald + UT_DGBAL[sp] * bal +
                  cr * UT_DGCR[sp] + cr * cr * UT_DGCRSQ[sp] +
                  c.dg_dsq[sp] * d * d + UT_DGPCCF[sp] * pccf
            sp == 5 && (dds *= 0.95f0)                # BS spruce (ut/dgf.f:587)
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        elseif sp == 6
            # Aspen (Utah DGFASP), ut/dgf.f:582-584. Raw crown pct (÷10 inside _em_dgfasp).
            cr_raw = Float32(t.crown_pct[i])
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            si = p.sp_site_index[sp]
            rmsqd = stand_qmd(s)
            aspdg = _em_dgfasp(d, cr_raw, bark, si, rmsqd, ba)
            cor2 = (s.control.dg_cor2_on && s.control.dg_cor2[sp] > 0f0) ? s.control.dg_cor2[sp] : 1f0
            dds = aspdg + log(cor2) + c.dg_cor[sp]
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        elseif sp in (11,12,13,14,15,16,24)
            # PJ/woodland DF-projection (ut/dgf.f:561-578).
            dpp = d < 1f0 ? 1f0 : d
            batem = ba < 1f0 ? 1f0 : ba
            si = p.sp_site_index[sp]
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            df = 0.25897f0 + 1.03129f0 * dpp - 0.0002025464f0 * batem + 0.00177f0 * si
            (df - dpp) > 1f0 && (df = dpp + 1f0)
            df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            dds = diagr <= 0f0 ? -9.21f0 : log(diagr * (2f0 * dpp * bark + diagr)) + conspp
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        elseif sp == 20 || sp == 21
            # MC/BI Wykoff-plus (ut/dgf.f:589-599). BAL uses BA (not BA100); adds DGDBAL·BAL/ln(D+1) + DGBA·BA.
            bal = (1f0 - pct / 100f0) * ba
            pt = Int(t.plot_id[i])
            pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
            dds = conspp + UT_DGLD[sp] * ald + UT_DGBAL[sp] * bal +
                  cr * (UT_DGCR[sp] + cr * UT_DGCRSQ[sp]) + c.dg_dsq[sp] * d * d +
                  UT_DGDBAL[sp] * bal / log(d + 1f0) + UT_DGPCCF[sp] * pccf + UT_DGBA[sp] * ba
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        else
            # CR-surrogate GB/NC/FC/BE (ut/dgf.f:606-624): BAU-based DF-proj + DSTAG. DEFERRED — needs the
            # BADIST BAU (BA-by-diameter-class) array + DSTAG (Zeide RELSDI stagnation); not exercised by utt01.
            # Faithful placeholder: GB(17) uses the case-3 DF form; NC/FC/BE use the BAU form (BAU=0 here → TODO).
            dpp = d < 1f0 ? 1f0 : d
            batem = ba < 5f0 ? 5f0 : ba
            si = p.sp_site_index[sp]
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            df = if sp == 17
                0.25897f0 + 1.03129f0 * dpp - 0.0002025464f0 * batem + 0.00177f0 * si
            else
                (1.55986f0 + 1.01825f0 * dpp - 0.29342f0 * log(batem) + 0.00672f0 * si) * 1.05f0   # BAU term TODO
            end
            sp == 17 && (df - dpp) > 1f0 && (df = dpp + 1f0)
            df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            dds = diagr <= 0f0 ? -9.21f0 : log(diagr * (2f0 * dpp * bark + diagr)) + c.dg_cor[sp] + c.dg_const[sp]
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        end
    end
    return s
end

diameter_growth!(s::StandState, ::Utah) = dgf!(s, Utah())
