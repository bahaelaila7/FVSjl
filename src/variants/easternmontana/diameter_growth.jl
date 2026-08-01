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
            ccf = -0.199592f0; c.atten[sp] = EM_OBSERV[isic, sp]; dgcon += 0.001766f0 * xsite
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

# em/dgf.f main body — per-tree WK2 = DDS (outside-bark). emt01 exercises only the MAIN Wykoff path.
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
            # NI section (LM/LL): BAL uses BA100 (TODO verify) — deferred, not in emt01.
            error("EM dgf! NI path (sp $sp) not yet ported (needs BA100)")
        elseif sp == 6 || sp == 11 || (13 <= sp <= 16) || sp == 19 || sp == 12 || sp == 17
            error("EM dgf! DIAGR/aspen path (sp $sp) not yet ported")
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
