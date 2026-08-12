# =============================================================================
# diameter_growth.jl (easternmontana) — EM large-tree DDS (em/dgf.f).
#
# Hooks (LS/KT pattern):
#   em_dgcons!(s)            — per-species per-stand DGCON/DGDSQ/DGCCF/ATTEN (em/dgf.f ENTRY DGCONS).
#   dgf!(s, ::EasternMontana)— per-tree WK2 = DDS (em/dgf.f main body). EM stores outside-bark DDS.
#
# EM DDS is species-dispatched (em/dgf.f DO 20 loop):
#   MAIN (conifers 1-3,7-10,18 = WB/WL/DF/LP/ES/AF/PP/OS): DDS = CONSPP + DGLD·lnD + DGBAL·BAL +
#     CR·(DGCR+CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1) + DGPCCF·RELDEN² + DGLCCF·ln(RELDEN) +
#     (managed?DGPCC1:DGPCC2)·PCCF.  CONSPP = DGCON+COR+0.01·DGCCF·RELDEN. BAL=(1−PCT/100)·BA.
#   NI (sp4 LM, sp5 LL): same minus the PCCF/RELDEN² terms; BAL uses BA100.
#   RM(6)/CO(11,13-16,19): CR-GENGYM DIAGR path.  ASPEN(12,17): DGFASP (Utah).  [deferred — not in emt01]
# DGCON (DGCONS) = DGHAB[ISPHAB]+DGFOR[ISPFOR]+DGEL·ELEV+DGEL2·ELEV²+(DGSASP·sinθ+DGCASP·cosθ+DGSLOP)·XSLOPE
#   +DGSLSQ·XSLOPE²; ISPHAB=MAPHAB[JDTYPE(Wykoff)|ITYPE, sp]; ISPFOR=MAPLOC[IFOR]; DGDSQ=DGDS[MAPDSQ[IFOR]].
# =============================================================================


# em/dgdriv.f:106 PSIGSQ — per-species prior variance for the DGSCOR COR calibration (read by the shared
# calibrate_diameter_growth!). EM DG residual SD dg_resid_sd = SIGMAR (em/dgdriv.f:414 SIGMA=SIGMAR).
const EM_PSIGSQ = Float32[0.0408, 0.0586, 0.1556, 0.0586, 0.0970, 0.07, 0.0636, 0.0970, 0.0970, 0.0636,
                          0.07, 0.1433, 0.07, 0.07, 0.07, 0.07, 0.1433, 0.0858, 0.07]

@inline _em_is_wykoff(sp::Int) = sp <= 3 || (7 <= sp <= 10) || sp == 18   # MAPHAB uses JDTYPE + XSLOPE÷10

# em/dgf.f ENTRY DGCONS — per-species per-stand DG constants (needs IEMTYP=habitat_code, ITYPE=habitat_input,
# IFOR=forest_idx, ELEV/ASPECT/SLOPE, SITEAR). Fills calib.dg_const(DGCON)/dg_dsq(DGDSQ)/dg_ccf(DGCCF)/atten.
function em_dgcons!(s::StandState)
    c = s.calib; p = s.plot; ctl = s.control; sd = s.coef.species
    jdtype = Int(p.habitat_code); jdtype > 117 && (jdtype = 30); jdtype < 1 && (jdtype = 1)
    itype  = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    ifor   = Int(p.forest_idx);    (ifor < 1 || ifor > 7) && (ifor = 1)
    elev = p.elevation; aspect = p.aspect; slope = p.slope
    ispccf = Int(EM_MAPCCF[itype])
    @inbounds for sp in 1:19
        xsite = p.sp_site_index[sp]
        lsi = trunc(Int, xsite / 10f0)
        isic = lsi < 2 ? 1 : (lsi >= 5 ? 5 : (lsi >= 4 ? 4 : (lsi >= 3 ? 3 : 2)))
        wy = _em_is_wykoff(sp)
        isphab = Int(EM_MAPHAB[wy ? jdtype : itype, sp]); (isphab < 1 || isphab > 8) && (isphab = 1)
        ispfor = Int(EM_MAPLOC[ifor, sp]); (ispfor < 1 || ispfor > 6) && (ispfor = 1)
        ispdsq = Int(EM_MAPDSQ[ifor, sp]); (ispdsq < 1 || ispdsq > 4) && (ispdsq = 1)
        tmpasp = aspect
        (sp == 4 || sp == 6 || sp == 12 || sp == 17) && (tmpasp -= 0.7854f0)
        xslope = wy ? slope / 10f0 : slope
        sp == 18 && (xslope = 0f0)
        dgcon = EM_DGHAB[isphab, sp] + EM_DGFOR[ispfor, sp] + EM_DGEL[sp] * elev + EM_DGEL2[sp] * elev * elev +
                (EM_DGSASP[sp] * sin(tmpasp) + EM_DGCASP[sp] * cos(tmpasp) + EM_DGSLOP[sp]) * xslope +
                EM_DGSLSQ[sp] * xslope * xslope
        c.dg_dsq[sp] = EM_DGDS[ispdsq, sp]
        ccf = 0f0
        if sp == 4
            # LM CCF coefficient: the base DGCCF(4)=−0.199592 (LM table) PLUS the em/dgf.f:489 addition
            # `CONSPP += 0.01·(−0.199592)·RELDEN` — live applies −0.199592 TWICE. jl previously had it once, so
            # CONSPP ran 0.01·0.199592·RELDEN too high ⇒ LM DDS/DG ~+17% over. MEASURED via FVSem_g16 (live CONSPP
            # 0.8055 vs jl 0.9689; post-fix jl LM DDS 2.1602==live 2.1602, 1.6371==live 1.6371). Fold both into ccf.
            ccf = -0.199592f0 - 0.199592f0; c.atten[sp] = EM_OBSERV[isic, sp]; dgcon += 0.001766f0 * xsite
        elseif sp == 5
            ccf = EM_DGCCFA[ispccf]; c.atten[sp] = EM_OBSERV[min(isphab, 6), sp]
        elseif sp == 12 || sp == 17
            c.atten[sp] = EM_OBSERV[isic, sp]; dgcon += 0.006460f0 * xsite
        else
            c.atten[sp] = EM_OBSERV[1, sp]
        end
        c.dg_ccf[sp] = ccf
        (ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0) && (dgcon += log(ctl.dg_cor2[sp]))
        c.dg_const[sp] = dgcon
        # EM bark (em/bratio.f) as bark_ratio(bark_a,bark_b)= (bark_a + bark_b·d)/d = bark_b + bark_a/d, for the
        # DBH update (simulate.jl dbh+=DG/bark). IMAP=2→BARK1 const (bark_a=0,bark_b=BARK1); IMAP=3→BARK1+BARK2/D
        # (bark_a=BARK2,bark_b=BARK1); IMAP=1 zero-coef→0.9002−0.3089/D (bark_a=−0.3089,bark_b=0.9002). Without
        # this, c.bark_a/b=0 ⇒ bark_ratio floors to 0.80 ⇒ DG/bark over-applies (0.934/0.80≈17% too much DBH).
        b1 = sd[:bark1][sp]; b2 = sd[:bark2][sp]; imapb = round(Int, sd[:bark_imap][sp])
        if imapb == 1 && b1 == 0f0 && b2 == 0f0
            c.bark_a[sp] = -0.3089f0; c.bark_b[sp] = 0.9002f0
        else
            c.bark_a[sp] = b2; c.bark_b[sp] = b1      # IMAP 2 (b2=0) or 3
        end
    end
    return s
end

# em/dgfasp.f (Utah aspen large-tree DG) — ASPDG = ln(DDS-equiv). Identical to ie_dgfasp/_tt_dgfasp
# (the shared UT form). cr = raw crown pct (÷10 → ASPCR inside). rmsqd = stand QMD.
@inline function _em_dgfasp(d::Float32, cr::Float32, bark::Float32, si::Float32, rmsqd::Float32, ba::Float32)::Float32
    rel = rmsqd > 0f0 ? d / rmsqd : 0f0
    aspcr = cr / 10f0
    pot = (0.4755f0 - 3.8336f-6 * d^4.1488f0) + (4.510f-2 * aspcr * d^0.67266f0)
    pot <= 0f0 && (pot = 0.01f0)
    fofr = 1.07528f0 * (1f0 - exp(-1.89022f0 * rel))
    gofad = 2.1963f-1 * (rmsqd + 1f0)^0.73355f0
    baact = ba >= 310f0 ? 305f0 : ba
    valmod = 1f0 - exp(-fofr * gofad * ((310f0 - baact) / 310f0)^0.5f0)
    predgr = pot * valmod * (0.48630f0 + 0.01258f0 * si)
    return log(2f0 * d * bark * predgr + predgr * predgr)
end

# em/dgf.f main body — per-tree WK2 = DDS (outside-bark). emt01 exercises only the MAIN Wykoff path;
# RM(6)/aspen(12,17)/CO-hardwood(11,13-16,19) are the non-conifer DIAGR + DGFASP forms (this chunk).
function dgf!(s::StandState, ::EasternMontana)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    wk2 = view(s.scratch.wk, 2, :)
    relden = p.relative_density
    alccf = relden > 0f0 ? log(relden) : 0f0
    ba = p.basal_area
    managed = p.managed == Int32(1)
    dum1 = managed ? 1f0 : 0f0; dum2 = managed ? 0f0 : 1f0
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i])
        ald = log(d)
        cr  = Float32(t.crown_pct[i]) * 0.01f0
        bal = (1f0 - t.crown_ratio[i] / 100f0) * ba                 # PCT = BA percentile
        conspp = c.dg_const[sp] + c.dg_cor[sp] + 0.01f0 * c.dg_ccf[sp] * relden
        if sp == 4 || sp == 5
            # NI section (LM/LL, em/dgf.f:563-566): BAL uses BA100=BA/100; no PCCF/RELDEN²/DGLCCF terms. conspp
            # carries the sp4 CCF = base DGCCF(4)+line-489 addition (BOTH −0.199592, see em_dgcons!). VALIDATED
            # 2026-08-12 vs FVSem_g16 on a real LM FIA stand (DDS bit-exact for the CCF-dominated terms; small
            # residual on mid-DBH from the BAL PCT=crown_ratio approximation, shared with the main-conifer path).
            bal100 = (1f0 - t.crown_ratio[i] / 100f0) * (ba / 100f0)
            dds = conspp + EM_DGLD[sp] * ald + EM_DGBAL[sp] * bal100 +
                  cr * (EM_DGCR[sp] + cr * EM_DGCRSQ[sp]) + c.dg_dsq[sp] * d * d +
                  EM_DGDBAL[sp] * bal100 / log(d + 1f0)
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        elseif sp == 6
            # RM juniper DIAGR (from UT), em/dgf.f:517-531. DF linear in DPP/BATEM/SI, capped +1"/period.
            dpp = d < 1f0 ? 1f0 : d
            batem = ba < 1f0 ? 1f0 : ba
            si = p.sp_site_index[sp]
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            df = 0.25897f0 + 1.03129f0 * dpp - 0.0002025464f0 * batem + 0.00177f0 * si
            (df - dpp) > 1f0 && (df = dpp + 1f0)
            df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            dds = diagr <= 0f0 ? -9.21f0 : log(diagr * (2f0 * dpp * bark + diagr)) + conspp
            dds < -9.21f0 && (dds = -9.21f0)                # shared floor (em/dgf.f:593)
            wk2[i] = dds
        elseif sp == 12 || sp == 17
            # Aspen (Utah DGFASP), also used for paper birch, em/dgf.f:535-538. Raw crown pct (÷10 inside).
            cr_raw = Float32(t.crown_pct[i])
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            si = p.sp_site_index[sp]
            rmsqd = _TT_CUR_RMSQD[] >= 0f0 ? _TT_CUR_RMSQD[] : stand_qmd(s)   # #195: current RMSQD during DGSCOR calibration
            aspdg = _em_dgfasp(d, cr_raw, bark, si, rmsqd, ba)
            cor2 = (s.control.dg_cor2_on && s.control.dg_cor2[sp] > 0f0) ? s.control.dg_cor2[sp] : 1f0
            dds = aspdg + log(cor2) + c.dg_cor[sp]
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        elseif sp == 11 || (13 <= sp <= 16) || sp == 19
            # CO/hardwood DIAGR (from CR), em/dgf.f:542-559. DF capped at 36"; floors the log BEFORE + AFTER
            # adding COR+DGCON (NOT conspp — no CCF·RELDEN term).
            dpp = d < 1f0 ? 1f0 : d
            si = p.sp_site_index[sp]
            bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
            df = 0.24506f0 + 1.01291f0 * dpp - 0.00084659f0 * ba + 0.00631f0 * si
            df > 36f0 && (df = 36f0)
            df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            if diagr <= 0f0
                dds = -9.21f0
            else
                dds = log(diagr * (2f0 * dpp * bark + diagr))
                dds < -9.21f0 && (dds = -9.21f0)            # internal floor (em/dgf.f:557)
            end
            dds = dds + c.dg_cor[sp] + c.dg_const[sp]
            dds < -9.21f0 && (dds = -9.21f0)                # shared floor (em/dgf.f:593)
            wk2[i] = dds
        else
            pt = Int(t.plot_id[i])
            pccf = (1 <= pt <= length(dens.point_ccf)) ? dens.point_ccf[pt] : 0f0
            dds = conspp +
                  EM_DGLD[sp]   * ald +
                  EM_DGBAL[sp]  * bal +
                  cr * (EM_DGCR[sp] + cr * EM_DGCRSQ[sp]) +
                  c.dg_dsq[sp]  * d * d +
                  EM_DGDBAL[sp] * bal / log(d + 1f0) +
                  EM_DGPCCF[sp] * relden * relden +
                  EM_DGLCCF[sp] * alccf +
                  (dum1 * EM_DGPCC1[sp] + dum2 * EM_DGPCC2[sp]) * pccf
            dds < -9.21f0 && (dds = -9.21f0)
            wk2[i] = dds
        end
    end
    return s
end
