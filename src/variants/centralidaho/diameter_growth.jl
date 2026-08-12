# =============================================================================
# diameter_growth.jl (centralidaho) — CI large-tree DDS (ci/dgf.f, western Wykoff). Chunk 3.
#   ci_dgcons!(s)             — per-species DGCON (ci/dgf.f ENTRY DGCONS)
#   dgf!(s, ::CentralIdaho)   — per-tree WK2 = DDS (ci/dgf.f main body)
# DGCON = DGHAB[ICHBCL(ICINDX,sp)+1] + DGFOR[MAPLOC(IFOR,sp)] + DGEL·TEMEL + DGEL2·TEMEL²
#       + (DGSASP·sinAsp + DGCASP·cosAsp + DGSLOP)·SLOPE + DGSLSQ·SLOPE² + site-term(+ln COR2).
# CONSPP = DGCON + COR + 0.01·DGCCFA·RELDEN.
# DDS default (main conifers 1-10,18): CONSPP + DGLD·lnD + DGLBA·ln(BA) + DGDSQ·D²
#       + DGDBAL·PBAL/ln(D+1) + DGBA·BAL/ln(D+1) + DGPCCF·PCCF.
# BAL = (1−PCT/100)·BA  (·BA/100 for sp 11,12,16). PBAL = (1−PCT/100)·PTBAA[point].
# Specials: 13 AS aspen(DGFASP); 14 WJ + 17 CW + 19 OH DIAGR; 15 MC full-Wykoff-with-CR.
# =============================================================================

const CI_NSP = 19

"""CI aspen large-tree DG (ci/dgfasp.f, shared Utah form). Returns ASPDG = ln(DDS-equivalent)."""
@inline function ci_dgfasp(d::Float32, cr::Float32, bark::Float32, si::Float32, rmsqd::Float32, ba::Float32)
    rel = rmsqd > 0f0 ? d / rmsqd : 1f0
    gofad = 1.0f0 - exp(-1.0f0 * max(rel, 0f0))
    dg = (0.22f0 + 0.01f0 * si) * gofad * exp(-0.005f0 * ba) * (cr / 10f0 + 0.5f0)
    dds = dg <= 0f0 ? -9.21f0 : log((dg * bark) * (2f0 * d * bark))
    return dds
end

function ci_dgcons!(s::StandState)
    c = s.calib; ctl = s.control; p = s.plot
    icindx = Int(p.habitat_input); ifor = Int(p.forest_idx)
    (icindx < 1 || icindx > 130) && (icindx = 1)
    (ifor < 1 || ifor > 6) && (ifor = 1)
    elev = p.elevation; slope = p.slope; asp = p.aspect
    @inbounds for sp in 1:CI_NSP
        xsite = p.sp_site_index[sp]
        lsi = trunc(Int, xsite / 10f0)
        isic = lsi < 2 ? 1 : lsi >= 5 ? 5 : lsi     # ci/dgf.f ISIC (LSI clamp 1..5)
        maphab = CI_ICHBCL[icindx, sp] + 1
        (maphab < 1 || maphab > 12) && (maphab = 1)
        ispfor = CI_MAPLOC[ifor, sp]
        (ispfor < 1 || ispfor > 3) && (ispfor = 1)
        tmpasp = (sp == 11 || sp == 12 || sp == 13 || sp == 16) ? asp - 0.7854f0 : asp
        temel = sp == 15 ? min(elev, 30f0) : elev
        dgcon = CI_DGHAB[maphab, sp] + CI_DGFOR[ispfor, sp] +
                CI_DGEL[sp] * temel + CI_DGEL2[sp] * temel * temel +
                (CI_DGSASP[sp] * sin(tmpasp) + CI_DGCASP[sp] * cos(tmpasp) + CI_DGSLOP[sp]) * slope +
                CI_DGSLSQ[sp] * slope * slope
        if sp == 11 || sp == 12 || sp == 16
            dgcon += 0.001766f0 * xsite
        elseif sp == 13
            dgcon += 0.006460f0 * xsite
        elseif sp == 15
            dgcon += 0.227307f0 * log(xsite)
        end
        # ATTEN (calibration observation count)
        if sp == 11 || sp == 12 || sp == 16
            c.atten[sp] = Float32(CI_IBSERV[isic, 1])
        elseif sp == 13
            c.atten[sp] = Float32(CI_IBSERV[isic, 2])
        elseif sp == 14 || sp == 17 || sp == 19
            c.atten[sp] = Float32(CI_IBSERV[isic, 3])
        else
            c.atten[sp] = Float32(CI_OBSERV[sp])
        end
        (ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0) && (dgcon += log(ctl.dg_cor2[sp]))
        c.dg_const[sp] = dgcon
        # linear-bark cache is not exact for CI's POWER model; DIAGR species use ci_bratio directly in dgf!.
        c.bark_a[sp] = 0f0; c.bark_b[sp] = 0.9f0
    end
    return s
end

"""CI `dgf!` hook — fill scratch.wk[2,i] with the DDS per live tree (ci/dgf.f main body)."""
function dgf!(s::StandState, ::CentralIdaho)
    p, t, c, ctl = s.plot, s.trees, s.calib, s.control
    dens = s.density; sd = s.coef.species
    wk2 = view(s.scratch.wk, 2, :)
    relden = p.relative_density
    ba = p.basal_area; ba100 = ba / 100f0
    logba = ba > 0f0 ? log(ba) : 0f0
    rmsqd = _TT_CUR_RMSQD[] >= 0f0 ? _TT_CUR_RMSQD[] : stand_qmd(s)   # #195: current RMSQD during DGSCOR calibration (aspen DGFASP)
    @inbounds for i in 1:t.n
        d = t.dbh[i]; d <= 0f0 && continue
        sp = Int(t.species[i])
        pt_i = Int(t.plot_id[i])
        pccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        ptba = (1 <= pt_i <= length(dens.point_ba)) ? dens.point_ba[pt_i] : 0f0
        conspp = c.dg_const[sp] + c.dg_cor[sp] + 0.01f0 * CI_DGCCFA[sp] * relden
        si = p.sp_site_index[sp]
        ald = log(d)
        balBA = (sp == 11 || sp == 12 || sp == 16) ? ba100 : ba
        bal = (1f0 - t.crown_ratio[i] / 100f0) * balBA
        pbal = (1f0 - t.crown_ratio[i] / 100f0) * ptba
        if sp == 13
            cr = Float32(t.crown_pct[i]); bark = ci_bratio(sd, sp, d)
            aspdg = ci_dgfasp(d, cr, bark, si, rmsqd, ba)
            cor2 = (ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0) ? ctl.dg_cor2[sp] : 1f0
            dds = aspdg + log(cor2) + c.dg_cor[sp]
        elseif sp == 14
            bark = ci_bratio(sd, sp, d)
            dpp = d < 1f0 ? 1f0 : d; batem = ba < 1f0 ? 1f0 : ba
            df = 0.25897f0 + 1.03129f0 * dpp - 0.0002025464f0 * batem + 0.00177f0 * si
            (df - dpp) > 1f0 && (df = dpp + 1f0); df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            dds = diagr <= 0f0 ? -9.21f0 : log(diagr * (2f0 * dpp * bark + diagr)) + conspp
        elseif sp == 15
            cr = Float32(t.crown_pct[i]) * 0.01f0
            dds = conspp + CI_DGLD[sp] * ald + CI_DGBAL[sp] * bal +
                  cr * (CI_DGCR[sp] + cr * CI_DGCRSQ[sp]) + CI_DGDS[sp] * d * d +
                  CI_DGBA[sp] * bal / log(d + 1f0)
        elseif sp == 17 || sp == 19
            bark = ci_bratio(sd, sp, d)
            dpp = d < 1f0 ? 1f0 : d
            df = 0.24506f0 + 1.01291f0 * dpp - 0.00084659f0 * ba + 0.00631f0 * si
            df > 36f0 && (df = 36f0); df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            dds = diagr <= 0f0 ? -9.21f0 : log(diagr * (2f0 * dpp * bark + diagr))
            dds < -9.21f0 && (dds = -9.21f0)
            dds = dds + c.dg_cor[sp] + c.dg_const[sp]
        else                                                 # DEFAULT: main conifers (1-10,16,18)
            cr = Float32(t.crown_pct[i]) * 0.01f0
            dds = conspp + CI_DGLD[sp] * ald + CI_DGLBA[sp] * logba +
                  CI_DGDS[sp] * d * d + CI_DGDBAL[sp] * pbal / log(d + 1f0) +
                  CI_DGBA[sp] * bal / log(d + 1f0) + CI_DGPCCF[sp] * pccf +
                  cr * (CI_DGCR[sp] + cr * CI_DGCRSQ[sp]) + CI_DGBAL[sp] * bal
        end
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end
