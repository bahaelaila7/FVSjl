# =============================================================================
# diameter_growth.jl (inlandempire) — IE large-tree DDS (ie/dgf.f, classic NI Wykoff).
#
# Coefficients live in dg_coefficients.jl. Two hooks (KT/LS pattern):
#   ie_dgcons!(s)            — per-species DGCON setup (ie/dgf.f ENTRY DGCONS)
#   dgf!(s, ::InlandEmpire)  — per-tree WK2 = DDS (ie/dgf.f main body).
#
# DDS forms (ie/dgf.f), dispatched by species:
#   NI (all except the special sp below):  BA100 = BA/100
#     CONSPP = DGCON + COR + 0.01*DGCCF*RELDEN
#     BAL    = (1 - PCT/100)*BA100
#     DDS = CONSPP + DGLD*ln(D) + DGBAL*BAL + CR*(DGCR + CR*DGCRSQ) + DGDSQ*D^2
#           + DGDBAL*BAL/ln(D+1)         [CR = ICR*0.01]
#   sp15,16 (PI/JU, UT logic):  DF = 0.25897+1.03129*DPP-0.0002025464*BATEM+0.00177*SI
#     DIAGR=(DF-DPP)*BARK;  DDS = ln(DIAGR*(2*DPP*BARK+DIAGR)) + CONSPP  (or -9.21)
#   sp18,20,21 (aspen, Utah DGFASP):  DDS = ASPDG + ln(COR2) + COR
#   sp19,22 (CO/OH, CR logic):  DF = 0.24506+1.01291*DPP-0.00084659*BA+0.00631*SI (cap 36)
#     DIAGR=(DF-DPP)*BARK;  DDS = ln(DIAGR*(2*DPP*BARK+DIAGR)) [clamp -9.21] + COR + DGCON
#   all: clamp DDS >= -9.21;  WK2(I) = DDS
#
# DGCON (ENTRY DGCONS), resolved once/stand:  ISPHAB=MAPHAB(ITYPE,sp), ISPFOR=MAPLOC(IFOR,sp),
#   ISPDSQ=MAPDSQ(IFOR,sp), ISPCCF=MAPCCF(ITYPE,sp).  TMPASP=ASPECT-0.7854 for sp{13,15-18,20,21}.
#   DGCON = DGHAB(ISPHAB)+DGFOR(ISPFOR)+DGEL*ELEV+DGEL2*ELEV^2
#           +(DGSASP*sin+DGCASP*cos)(TMPASP)+DGSLOP)*SLOPE+DGSLSQ*SLOPE^2
#   DGDSQ=DGDS(ISPDSQ,sp); DGCCF=DGCCFA(ISPCCF,sp).  ISIC = site class from XSITE=SITEAR (5-class).
#   ATTEN = OBSERV(ISPHAB) for sp{<=12,14,23}, else OBSERV(ISIC); sp13/17 DGCON+=.001766*XSITE,
#   sp18/20/21 DGCON+=.006460*XSITE.  If LDCOR2 & COR2>0: DGCON += ln(COR2).
# =============================================================================

const IE_NSP = 23

"""IE bark ratio (ie/bratio.f): IMAP-dispatched BARK1/BARK2 form. sp15/16 use the 0.9002-0.3089/D
default (D cap 19); most species a constant BARK1; sp19/22 BARK1+BARK2/D. Clamped [0.80, 0.99]."""
@inline function ie_bratio(sp::Int, d::Float32)
    ieqn = IE_BRK_IMAP[sp]
    temd = d < 1f0 ? 1f0 : d
    br = if ieqn == 1
        if IE_BARK1[sp] == 0f0 && IE_BARK2[sp] == 0f0
            temd > 19f0 && (temd = 19f0)
            0.9002f0 - 0.3089f0 * (1f0 / temd)
        else
            IE_BARK1[sp] + IE_BARK2[sp] * (1f0 / temd)
        end
    elseif ieqn == 2
        IE_BARK1[sp]
    else
        IE_BARK1[sp] + IE_BARK2[sp] * (1f0 / temd)
    end
    br > 0.99f0 && (br = 0.99f0)
    br < 0.80f0 && (br = 0.80f0)
    return br
end

"""IE aspen large-tree diameter growth (ie/dgfasp.f, Utah). Returns ASPDG = ln(DDS-equivalent)."""
@inline function ie_dgfasp(d::Float32, cr::Float32, bark::Float32, si::Float32, rmsqd::Float32, ba::Float32)
    xsite = si
    rel = d / rmsqd
    aspcr = cr / 10f0
    pot = (0.4755f0 - 3.8336f-6 * d^4.1488f0) + (4.510f-2 * aspcr * d^0.67266f0)
    pot <= 0f0 && (pot = 0.01f0)
    fofr = 1.07528f0 * (1f0 - exp(-1.89022f0 * rel))
    gofad = 2.1963f-1 * (rmsqd + 1f0)^0.73355f0
    baact = ba >= 310f0 ? 305f0 : ba
    valmod = 1f0 - exp(-fofr * gofad * ((310f0 - baact) / 310f0)^0.5f0)
    predgr = pot * valmod * (0.48630f0 + 0.01258f0 * xsite)
    return log(2f0 * d * bark * predgr + predgr * predgr)
end

"""IE per-stand DG setup (ie/dgf.f ENTRY DGCONS): fill calib.dg_const (DGCON) + calib.atten (ATTEN)."""
function ie_dgcons!(s::StandState)
    c = s.calib; ctl = s.control; p = s.plot
    itype = Int(p.habitat_input); ifor = Int(p.forest_idx)
    (itype < 1 || itype > 30) && (itype = 1)
    (ifor  < 1 || ifor  > 11) && (ifor  = 1)
    elev = p.elevation; slope = p.slope; asp = p.aspect
    @inbounds for sp in 1:IE_NSP
        xsite = p.sp_site_index[sp]
        lsi = trunc(Int, xsite / 10f0)
        isic = lsi < 2 ? 1 : lsi >= 5 ? 5 : lsi == 4 ? 4 : lsi == 3 ? 3 : 2
        isphab = Int(IE_MAPHAB[itype, sp]); ispfor = Int(IE_MAPLOC[ifor, sp])
        (isphab < 1 || isphab > 6) && (isphab = 1)
        (ispfor < 1 || ispfor > 6) && (ispfor = 1)
        tmpasp = asp
        (sp == 13 || (15 <= sp <= 18) || sp == 20 || sp == 21) && (tmpasp -= 0.7854f0)
        dgcon = IE_DGHAB[isphab, sp] + IE_DGFOR[ispfor, sp] +
                IE_DGEL[sp] * elev + IE_DGEL2[sp] * elev * elev +
                (IE_DGSASP[sp] * sin(tmpasp) + IE_DGCASP[sp] * cos(tmpasp) + IE_DGSLOP[sp]) * slope +
                IE_DGSLSQ[sp] * slope * slope
        if sp <= 12 || sp == 14 || sp == 23
            c.atten[sp] = Float32(IE_OBSERV[isphab, sp])
        else
            c.atten[sp] = Float32(IE_OBSERV[isic, sp])
            if sp == 13 || sp == 17
                dgcon += 0.001766f0 * xsite
            elseif sp == 18 || sp == 20 || sp == 21
                dgcon += 0.006460f0 * xsite
            end
        end
        (ctl.dg_cor2_on && ctl.dg_cor2[sp] > 0f0) && (dgcon += log(ctl.dg_cor2[sp]))
        c.dg_const[sp] = dgcon
        # Set the linear bark (bark_a + bark_b*d)/d so the SHARED engine bark_ratio calls (DBH-apply in
        # grow_cycle!, _backdate_dbh!, etc.) reproduce ie_bratio (ie/bratio.f) — without this bark_a/bark_b
        # stay 0 ⇒ bark_ratio floors to 0.80 ⇒ DG over-applied ~9% ⇒ QMD/BA over-growth.
        #   IEQN=2 → a=0,b=BARK1 (constant BARK1); IEQN=1/3 → a=BARK2,b=BARK1 (BARK1+BARK2/d);
        #   sp15/16 (BARK1=BARK2=0) → a=-0.3089,b=0.9002 (the 0.9002-0.3089/d default).
        if IE_BRK_IMAP[sp] == 2
            c.bark_a[sp] = 0f0; c.bark_b[sp] = IE_BARK1[sp]
        elseif IE_BARK1[sp] == 0f0 && IE_BARK2[sp] == 0f0
            c.bark_a[sp] = -0.3089f0; c.bark_b[sp] = 0.9002f0
        else
            c.bark_a[sp] = IE_BARK2[sp]; c.bark_b[sp] = IE_BARK1[sp]
        end
    end
    return s
end

"""IE `dgf!` hook — fill scratch.wk[2,i] with the DDS per live tree (ie/dgf.f main body)."""
function dgf!(s::StandState, ::InlandEmpire)
    p, t, c, ctl = s.plot, s.trees, s.calib, s.control
    wk2 = view(s.scratch.wk, 2, :)
    relden = p.relative_density
    ba = p.basal_area
    ba100 = ba / 100f0
    rmsqd = p.qmd
    itype = Int(p.habitat_input); ifor = Int(p.forest_idx)
    (itype < 1 || itype > 30) && (itype = 1)
    (ifor  < 1 || ifor  > 11) && (ifor  = 1)
    @inbounds for i in 1:t.n
        d = t.dbh[i]
        d <= 0f0 && continue
        sp = Int(t.species[i])
        ispdsq = Int(IE_MAPDSQ[ifor, sp]); ispccf = Int(IE_MAPCCF[itype, sp])
        (ispdsq < 1 || ispdsq > 4) && (ispdsq = 1)
        (ispccf < 1 || ispccf > 5) && (ispccf = 1)
        dgdsq = IE_DGDS[ispdsq, sp]; dgccf = IE_DGCCFA[ispccf, sp]
        conspp = c.dg_const[sp] + c.dg_cor[sp] + 0.01f0 * dgccf * relden
        si = p.sp_site_index[sp]
        ald = log(d)
        if sp == 15 || sp == 16
            # PI/JU logic from UT
            bark = ie_bratio(sp, d)
            dpp = d < 1f0 ? 1f0 : d
            batem = ba < 1f0 ? 1f0 : ba
            df = 0.25897f0 + 1.03129f0 * dpp - 0.0002025464f0 * batem + 0.00177f0 * si
            (df - dpp) > 1f0 && (df = dpp + 1f0)
            df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            dds = diagr <= 0f0 ? -9.21f0 : log(diagr * (2f0 * dpp * bark + diagr)) + conspp
        elseif sp == 18 || sp == 20 || sp == 21
            # aspen logic from Utah (DGFASP)
            bark = ie_bratio(sp, d)
            cr = Float32(t.crown_pct[i])                 # NOTE: ICR (percent), NOT *0.01
            aspdg = ie_dgfasp(d, cr, bark, si, rmsqd, ba)
            dds = aspdg + log(ctl.dg_cor2[sp]) + c.dg_cor[sp]
        elseif sp == 19 || sp == 22
            # CO logic from CR
            bark = ie_bratio(sp, d)
            dpp = d < 1f0 ? 1f0 : d
            df = 0.24506f0 + 1.01291f0 * dpp - 0.00084659f0 * ba + 0.00631f0 * si
            df > 36f0 && (df = 36f0)
            df < dpp && (df = dpp)
            diagr = (df - dpp) * bark
            dds = diagr <= 0f0 ? -9.21f0 : log(diagr * (2f0 * dpp * bark + diagr))
            dds < -9.21f0 && (dds = -9.21f0)
            dds = dds + c.dg_cor[sp] + c.dg_const[sp]
        else
            # ORIGINAL NI SECTION
            cr = Float32(t.crown_pct[i]) * 0.01f0
            bal = (1f0 - t.crown_ratio[i] / 100f0) * ba100
            dds = conspp + IE_DGLD[sp] * ald + IE_DGBAL[sp] * bal +
                  cr * (IE_DGCR[sp] + cr * IE_DGCRSQ[sp]) +
                  dgdsq * d * d + IE_DGDBAL[sp] * bal / log(d + 1f0)
        end
        dds < -9.21f0 && (dds = -9.21f0)
        wk2[i] = dds
    end
    return s
end
