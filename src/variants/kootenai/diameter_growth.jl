# =============================================================================
# diameter_growth.jl (kootenai) — KT large-tree DDS (kt/dgf.f, western Wykoff).
#
# Uses coefficients in dg_coefficients.jl + habitat_tables.jl. Two hooks (LS pattern):
#   kt_dgcons!(s)         — per-species DG constant DGCON (kt/dgf.f ENTRY DGCONS)
#   dgf!(s, ::Kootenai)   — per-tree WK2 = DDS (kt/dgf.f main body). KT stores the OUTSIDE-bark DDS
#                           directly (NO OB->IB conversion, unlike LS).
#
# DDS (kt/dgf.f:341): CONSPP + DGLD*ln(D) + CR*(DGCR+CR*DGCRSQ) + DGDBAL*BAL/ln(D+1) + CCFSQ*RELDEN^2
#   + DGDS*D^2 + DGLBA*ln(BA) + (MANAGD? DGPCC1 : DGPCC2)*PCCF1;  clamp >= -9.21.
#   CONSPP = DGCON(sp) + COR(sp) + DGCCFA(sp)*RELDEN  [sp11(OT): 0.01*DGCCFA*RELDEN]
#   DGCON  = DGHAB[ISPHAB] + DGFOR[ISPFOR] + DGEL*ELEV + DGEL2*ELEV^2
#            + (DGSASP*sin(ASP)+DGCASP*cos(ASP)+DGSLOP)*SLOPE + DGSLSQ*SLOPE^2 (+ ln(COR2))
#   ISPHAB = MAPHAB[KKTYPE, sp]  (KKTYPE = plot.habitat_code, set by site_setup!);  ISPFOR = MAPLOC[KOTFOR, sp]
#   BAL = (1 - PCT/100)*BA  [sp11: /100];  RELDEN = plot.relative_density
# =============================================================================

"""KT per-stand DG setup (kt/dgf.f ENTRY DGCONS): fill calib.dg_const (DGCON), atten (OBSERV), bark copy."""
function kt_dgcons!(s::StandState)
    c = s.calib; ctl = s.control; p = s.plot
    kktype = Int(p.habitat_code)            # KKTYPE (site_setup! stored it here)
    kotfor = Int(p.forest_idx)              # KOTFOR forest index (1..10) for MAPLOC
    elev = p.elevation; slope = p.slope; asp = p.aspect
    @inbounds for sp in 1:nspecies(Kootenai())
        isphab = (1 <= kktype <= 175) ? Int(KT_MAPHAB[kktype, sp]) : 1
        ispfor = (1 <= kotfor <= 10)  ? Int(KT_MAPLOC[kotfor, sp]) : 1
        (isphab < 1 || isphab > 9) && (isphab = 1)
        (ispfor < 1 || ispfor > 7) && (ispfor = 1)
        dgcon = KT_DGHAB[isphab, sp] + KT_DGFOR[ispfor, sp] +
                KT_DGEL[sp] * elev + KT_DGEL2[sp] * elev * elev +
                (KT_DGSASP[sp] * sin(asp) + KT_DGCASP[sp] * cos(asp) + KT_DGSLOP[sp]) * slope +
                KT_DGSLSQ[sp] * slope * slope
        (ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0) && (dgcon += log(ctl.dg_cor2[sp]))
        c.dg_const[sp] = dgcon
        c.atten[sp]    = Float32(KT_OBSERV[isphab, sp])   # ATTEN = OBSERV(ISPHAB, ISPC)
        c.bark_a[sp]   = 0f0                              # KT bark: intercept 0, slope BKRAT (bratio chunk)
        c.bark_b[sp]   = KT_BKRAT[sp]
    end
    return s
end

"""KT `dgf!` hook — fill scratch.wk[2,i] with the (outside-bark) DDS per live tree (kt/dgf.f main body)."""
function dgf!(s::StandState, ::Kootenai)
    p, t, c = s.plot, s.trees, s.calib
    wk2 = view(s.scratch.wk, 2, :)
    relden = p.relative_density
    ba = p.basal_area
    lnba = ba > 0f0 ? log(ba) : 0f0
    managed = p.managed == Int32(1)
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d <= 0f0 && continue
        sp = Int(t.species[i])
        rd = sp == 11 ? 0.01f0 * KT_DGCCFA[sp] * relden : KT_DGCCFA[sp] * relden
        conspp = c.dg_const[sp] + c.dg_cor[sp] + rd
        ald = log(d)
        cr  = Float32(t.crown_pct[i])                          # crown ratio (CR)
        bal = (1f0 - t.crown_ratio[i] / 100f0) * ba            # PCT = BA percentile
        sp == 11 && (bal = bal / 100f0)
        pt = Int(t.plot_id[i])
        pccf1 = (1 <= pt <= length(p.point_ccf)) ? p.point_ccf[pt] : 0f0
        dgpc = managed ? KT_DGPCC1[sp] : KT_DGPCC2[sp]
        dds = conspp +
              KT_DGLD[sp]   * ald +
              cr * (KT_DGCR[sp] + cr * KT_DGCRSQ[sp]) +
              KT_DGDBAL[sp] * bal / log(d + 1f0) +
              KT_CCFSQ[sp]  * relden * relden +
              KT_DGDS[sp]   * d * d +
              KT_DGLBA[sp]  * lnba +
              dgpc * pccf1
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end
