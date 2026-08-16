# =============================================================================
# organon_nwo.jl — OP (Olympic) ORGANON NWO (VERSION=2) growth engine.
#
# The NWO follow-on to OC's SWO engine. The ORGANON *structure* is version-agnostic and the
# version-agnostic helpers (OrganonBuffer, GET_BAL, GET_CCFL, CALTST, QUAD1, SPMIX, DFORTY,
# OLDGROWTH, OLDGRO, PM_FERT) are REUSED from the oregoncoast port. What NWO needs on top is the
# NWO coefficient tables + the handful of NWO structural differences, all ported here verbatim from
# organon/*.f (the DG_NWO diameter core is in organon_diamgro_nwo.jl):
#
#   • DIAMGRO_RUN site: SITE = SI_2 if FIA==263 (WH) else SI_1 (diagro.f:45-50).
#   • Height: HTGRO1 CASE(2) uses B_HG (Bruce 1981, DF/GF) / F_HG (Flewelling WH, group 3) potential
#     height + HG_NWO modifier (htgrowth.f:99-118). Minor species → HTGRO2 via HD_NWO.
#   • Crown: CROWGRO CASE(2) applies CALIB(2,ISPGRP) to PCR1/PCR2 (crngrow.f:70-73,96-101) — the KEY
#     NWO/SWO crown difference (SWO uses no CALIB(2)). HCB_NWO/MAXHCB_NWO/MCW/LCW/HLCW/CW _NWO tables.
#   • Mortality: PM_NWO — a DIFFERENT functional form per group (DF uses √DBH & CR^0.25; GF uses
#     BAL/DBH), MPAR(11,7) (mortality.f:541-616). MORTAL_RUN driver is the VERSION≤3 shared structure.
#   • SUBMAX VERSION=2: A2=0.62305, TEMPA1 default 6.19958, OCMOD default 1.014293245, PWH branch.
#   • PREPARE VERSION=2: IB=3, NSPN=11, SI conversion DF↔WH, HDCALIB via B_H40/HD40_NWO/HD_NWO,
#     CRCALIB via HCB_NWO — produces ACALIB(1,·) (height) + ACALIB(2,·) (crown, consumed by CROWGRO).
#
# The NWO species group (1..11 = DF,GF,WH,RC,PY,MD,BL,WO,RA,PD,WI) comes from op_spgroup_nwo
# (species.jl). Big species IB=3 (DF,GF,WH → HTGRO1); groups 4..11 minor (HTGRO2). Bark ratio =
# op_bratio (op/bratio.f, 39-species). VALIDATED bit-exact vs the live FVSop_clean oracle DGDRIV/
# HTGF/CROWN/MORTS dumps, stand S248112 (27-tree buffer). See docs/OP_VARIANT_PORT_AUDIT.md.
# =============================================================================

# --- op/bratio.f — OP (39-species) bark ratio BARKB(4,16)/JBARK(39) ----------------------------
const OP_BRATIO_JBARK = Int[
    2,2,2,2,2,11,2,5,5,9, 10,4,4,4,3,1,16,12,13,12,
    6,14,14,6,7,14,14,8,12,10, 15,15,15,14,14,14,14,10,10]
# BARKB rows 1..16: fia, a (col2), b (col3), type (col4).
const OP_BRATIO_A = Float32[0.903563,0.904973,0.809427,0.859045,0.837291,0.08360,0.15565,0.8558,0.9,0.9,0.958330,0.949670,0.933710,0.075256,0.933290,0.7012]
const OP_BRATIO_B = Float32[0.989388,1.0,1.016866,1.0,1.0,0.94782,0.90182,1.0213,0.0,0.0,1.0,1.0,1.0,0.94373,1.0,1.04862]
const OP_BRATIO_TYPE = Int[1,1,1,1,1,2,2,1,3,3,1,1,1,2,1,1]

"op/bratio.f BRATIO — OP bark ratio for FVS species index `is` (1..39) at DBH `d`."
@inline function op_bratio(is::Integer, d::Float32)
    d <= 0.0f0 && return 0.99f0
    r = OP_BRATIO_JBARK[is]
    a = OP_BRATIO_A[r]; b = OP_BRATIO_B[r]; ty = OP_BRATIO_TYPE[r]
    if ty == 1
        br = (a*fpow(d, b))/d
    elseif ty == 2
        br = (a + b*d)/d
    else
        br = a
    end
    br > 0.99f0 && (br = 0.99f0)
    br < 0.80f0 && (br = 0.80f0)
    return br
end

"op/dgdriv.f:433-450 — ORGANON DGRO → BARK → DIAGR → DDS (WK2) copy-back."
@inline function op_organon_dds(is::Integer, dbh::Float32, dgro::Float32)
    bark = op_bratio(is, dbh)
    diagr = dgro*bark
    dds = flog(diagr*(2.0f0*dbh*bark + diagr))
    dds < -9.21f0 && (dds = -9.21f0)
    return bark, diagr, dds
end

"Shared FVS DDS→DG (inside-bark 5-yr increment) for an ORGANON tree (op/dgdriv.f copy-back)."
@inline function op_organon_dg(is::Integer, dbh::Float32, dds::Float32)
    bark = op_bratio(is, dbh)
    d = dbh*bark
    return sqrt(d*d + fexp(dds)) - d
end

# =============================================================================
# NWO coefficient leaf functions (organon/*.f, groups 1..11 = DF,GF,WH,RC,PY,MD,BL,WO,RA,PD,WI)
# =============================================================================

# --- crngrow.f MCW_NWO MCWPAR(11,4): B0,B1,B2,PKDBH -------------------------------------------
const OP_MCW_NWO_B0 = Float32[4.6198,6.1880,4.3586,4.0,4.5652,3.4298629,4.0953,3.0785639,8.0,2.9793895,2.9793895]
const OP_MCW_NWO_B1 = Float32[1.8426,1.0069,1.57458,1.65,1.4147,1.3532302,2.3849,1.9242211,1.53,1.5512443,1.5512443]
const OP_MCW_NWO_B2 = Float32[-0.011311,0.0,0.0,0.0,0.0,0.0,-0.0102651,0.0,0.0,-0.01416129,-0.01416129]
const OP_MCW_NWO_PK = Float32[81.45,999.99,76.70,999.99,999.99,999.99,102.53,999.99,999.99,54.77,54.77]

"organon/crngrow.f MCW_NWO — maximum crown width, group `g`, DBH `d`, height `h`."
@inline function op_mcw_nwo(g::Int, d::Float32, h::Float32)
    b0 = OP_MCW_NWO_B0[g]; b1 = OP_MCW_NWO_B1[g]; b2 = OP_MCW_NWO_B2[g]; pk = OP_MCW_NWO_PK[g]
    dbh = d > pk ? pk : d
    if h < 4.501f0
        return h/4.5f0*b0
    else
        return b0 + b1*dbh + b2*dbh*dbh
    end
end

# --- crngrow.f LCW_NWO LCWPAR(11,3): B1,B2,B3 --------------------------------------------------
const OP_LCW_NWO_B1 = Float32[0.0,0.0,0.105590,-0.2513890,0.0,0.118621,0.0,0.3648110,0.3227140,0.0,0.0]
const OP_LCW_NWO_B2 = Float32[0.00436324,0.00308402,0.0035662,0.006925120,0.0,0.00384872,0.0,0.0,0.0,0.0,0.0]
const OP_LCW_NWO_B3 = Float32[0.6020020,0.0,0.0,0.985922,0.0,0.0,1.470180,0.0,0.0,1.61440,1.61440]

"organon/crngrow.f LCW_NWO — largest crown width (SCR=0 branch)."
@inline function op_lcw_nwo(g::Int, mcw::Float32, cr::Float32, scr::Float32, dbh::Float32, ht::Float32)
    b1 = OP_LCW_NWO_B1[g]; b2 = OP_LCW_NWO_B2[g]; b3 = OP_LCW_NWO_B3[g]
    if scr > cr
        cl = scr*ht
        return mcw*fpow(scr, b1 + b2*cl + b3*(dbh/ht))
    else
        cl = cr*ht
        return mcw*fpow(cr, b1 + b2*cl + b3*(dbh/ht))
    end
end

# --- crngrow.f HLCW_NWO DACBPAR(11) -----------------------------------------------------------
const OP_HLCW_NWO_B1 = Float32[0.062000,0.028454,0.355270,0.209806,0.209806,0.0,0.0,0.0,0.0,0.0,0.0]

"organon/crngrow.f HLCW_NWO — height to largest crown width."
@inline function op_hlcw_nwo(g::Int, ht::Float32, cr::Float32, scr::Float32)
    b1 = OP_HLCW_NWO_B1[g]
    cl = (scr > cr ? scr : cr)*ht
    return ht - (1.0f0 - b1)*cl
end

# --- crngrow.f CW_NWO CWAPAR(11,3): B1,B2,B3 --------------------------------------------------
const OP_CW_NWO_B1 = Float32[0.929973,0.999291,0.461782,0.629785,0.629785,0.5,0.5,0.5,0.5,0.5,0.5]
const OP_CW_NWO_B2 = Float32[-0.135212,0.0,0.552011,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const OP_CW_NWO_B3 = Float32[-0.0157579,-0.0314603,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

"organon/crngrow.f CW_NWO — crown width at height XL (above largest-crown-width height)."
@inline function op_cw_nwo(g::Int, hlcw::Float32, lcw::Float32, ht::Float32, dbh::Float32, xl::Float32)
    b1 = OP_CW_NWO_B1[g]; b2 = OP_CW_NWO_B2[g]; b3 = OP_CW_NWO_B3[g]
    rp = (ht - xl)/(ht - hlcw)
    ratio = ht/dbh
    if g == 1
        ratio > 50.0f0 && (ratio = 50.0f0)
    elseif g == 2
        ratio > 31.0f0 && (ratio = 31.0f0)
    end
    return lcw*fpow(rp, b1 + b2*sqrt(rp) + b3*ratio)
end

# --- crngrow.f HCB_NWO HCBPAR(11,7): B0..B6 ---------------------------------------------------
const OP_HCB_NWO_B0 = Float32[1.94093,1.04746,1.92682,4.49102006,0.0,2.955339267,0.9411395642,1.05786632,0.56713781,0.0,0.0]
const OP_HCB_NWO_B1 = Float32[-0.0065029,-0.0066643,-0.00280478,0.0,0.0,0.0,-0.00768402,0.0,-0.010377976,0.0,0.0]
const OP_HCB_NWO_B2 = Float32[-0.0048737,-0.0067129,-0.0011939,-0.00132412,0.0,0.0,-0.005476131,-0.00183283,-0.002066036,-0.005666559,-0.005666559]
const OP_HCB_NWO_B3 = Float32[-0.261573,0.0,-0.513134,-1.01460531,0.0,-0.798610738,0.0,-0.28644547,0.0,-0.745540494,-0.745540494]
const OP_HCB_NWO_B4 = Float32[1.08785,0.0,3.68901,0.0,2.030940382,3.095269471,0.0,0.0,1.39796223,0.0,0.0]
const OP_HCB_NWO_B5 = Float32[0.0,0.0,0.00742219,0.01340624,0.0,0.0,0.0,0.0,0.0,0.038476613,0.038476613]
const OP_HCB_NWO_B6 = Float32[0.0,0.0,0.0,0.0,0.0,0.700465646,0.0,0.0,0.0,0.0,0.0]

"organon/crngrow.f HCB_NWO — height to crown base (group `g`); SI_2 for WH group 3, SI_1 else."
@inline function op_hcb_nwo(g::Int, ht::Float32, dbh::Float32, ccfl::Float32, ba::Float32,
                            si_1::Float32, si_2::Float32, og::Float32)
    b0 = OP_HCB_NWO_B0[g]; b1 = OP_HCB_NWO_B1[g]; b2 = OP_HCB_NWO_B2[g]; b3 = OP_HCB_NWO_B3[g]
    b4 = OP_HCB_NWO_B4[g]; b5 = OP_HCB_NWO_B5[g]; b6 = OP_HCB_NWO_B6[g]
    si = g == 3 ? si_2 : si_1
    return ht/(1.0f0 + fexp(b0 + b1*ht + b2*ccfl + b3*flog(ba) + b4*(dbh/ht) + b5*si + b6*(og*og)))
end

# --- crngrow.f MAXHCB_NWO MAXPAR(11,5): B0,B1,B2,B3,LIMIT -------------------------------------
const OP_MAXHCB_NWO_B0  = Float32[0.96,0.96,1.01,0.96,0.85,0.981,1.0,1.0,0.93,1.0,0.985]
const OP_MAXHCB_NWO_B1  = Float32[0.26,0.31,0.36,0.31,0.35,0.161,0.45,0.3,0.18,0.45,0.285]
const OP_MAXHCB_NWO_B2  = Float32[-0.900721383,-2.450718394,-0.944528054,-1.059636222,-0.922868139,-1.73666044,-1.020016685,-0.95634399,-0.928243505,-1.020016685,-0.969750805]
const OP_MAXHCB_NWO_B3  = Float32[1.0,1.0,0.6,1.0,0.8,1.0,1.0,1.1,1.0,1.0,0.9]
const OP_MAXHCB_NWO_LIM = Float32[0.95,0.95,0.96,0.95,0.80,0.98,0.95,0.98,0.92,0.95,0.98]

"organon/crngrow.f MAXHCB_NWO — maximum height to crown base (group `g`)."
@inline function op_maxhcb_nwo(g::Int, ht::Float32, ccfl::Float32)
    b0 = OP_MAXHCB_NWO_B0[g]; b1 = OP_MAXHCB_NWO_B1[g]; b2 = OP_MAXHCB_NWO_B2[g]
    b3 = OP_MAXHCB_NWO_B3[g]; lim = OP_MAXHCB_NWO_LIM[g]
    maxbr = b0 - b1*fexp(b2*fpow(ccfl/100.0f0, b3))
    maxbr > lim && (maxbr = lim)
    return maxbr*ht
end

# --- htgrowth.f HG_NWO HGPAR(3,8): P1..P8 (big-3 conifers DF,GF,WH only) ----------------------
const OP_HG_NWO_P1 = Float32[0.655258886,1.0,1.0]
const OP_HG_NWO_P2 = Float32[-0.006322913,-0.0328142,-0.0384415]
const OP_HG_NWO_P3 = Float32[-0.039409636,-0.0127851,-0.0144139]
const OP_HG_NWO_P4 = Float32[0.5,1.0,0.5]
const OP_HG_NWO_P5 = Float32[0.597617316,6.19784,1.04409]
const OP_HG_NWO_P6 = Float32[2.0,2.0,2.0]
const OP_HG_NWO_P7 = Float32[0.631643636,0.0,0.0]
const OP_HG_NWO_P8 = Float32[1.010018427,1.01,1.03]

"organon/htgrowth.f HG_NWO — height-increment modifier (big-3 groups 1..3)."
@inline function op_hg_nwo(g::Int, phtgro::Float32, cr::Float32, tcch::Float32)
    p1 = OP_HG_NWO_P1[g]; p2 = OP_HG_NWO_P2[g]; p3 = OP_HG_NWO_P3[g]; p4 = OP_HG_NWO_P4[g]
    p5 = OP_HG_NWO_P5[g]; p6 = OP_HG_NWO_P6[g]; p7 = OP_HG_NWO_P7[g]; p8 = OP_HG_NWO_P8[g]
    fcr = (-p5*fpow(1.0f0 - cr, p6))*fexp(p7*fpow(tcch, 0.5f0))
    b0 = p1*fexp(p2*tcch)
    b1 = fexp(p3*fpow(tcch, p4))
    modifer = fcr < -20.0f0 ? p8*b0 : p8*(b0 + (b1 - b0)*fexp(fcr))
    cradj = 1.0f0
    cr <= 0.17f0 && (cradj = 1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0)))
    return phtgro*modifer*cradj
end

"organon/htgrowth.f B_HG — Bruce (1981) DF/GF potential height growth (5-yr, GP=5)."
@inline function op_b_hg(si::Float32, ht::Float32, gp::Float32)
    x1 = 13.25f0 - si/20.0f0
    x2 = 63.25f0 - si/20.0f0
    b2 = -0.447762f0 - 0.894427f0*si/100.0f0 + 0.793548f0*fpow(si/100.0f0, 2.0f0) -
         0.171666f0*fpow(si/100.0f0, 3.0f0)
    b1 = flog(4.5f0/si)/(fpow(x1, b2) - fpow(x2, b2))
    xx1 = flog(ht/si)/b1 + fpow(x2, b2)
    geage = xx1 > 0.0f0 ? fpow(xx1, 1.0f0/b2) - x1 : 500.0f0
    pht = si*fexp(b1*(fpow(geage + gp + x1, b2) - fpow(x2, b2)))
    return geage, pht - ht
end

"organon/start2.f B_H40 — Bruce (1981) dominant top height at age (for HDCALIB HT40)."
@inline function op_b_h40(age::Float32, si::Float32)
    a2 = -0.447762f0 - 0.894427f0*si/100.0f0 + 0.793548f0*fpow(si/100.0f0, 2.0f0) -
         0.17166f0*fpow(si/100.0f0, 3.0f0)
    a1 = flog(4.5f0/si)/(fpow(13.25f0 - si/20.0f0, a2) - fpow(63.25f0 - si/20.0f0, a2))
    return si*fexp(a1*(fpow(age + 13.25f0 - si/20.0f0, a2) - fpow(63.25f0 - si/20.0f0, a2)))
end

# --- htgrowth.f HD_NWO HDPAR(11,3): B0,B1,B2 (height-diameter, all groups; HTGRO2 + HDCALIB) ---
const OP_HD_NWO_B0 = Float32[7.04524,7.42808,5.93792,6.14817441,9.30172,5.84487,5.21462,4.69753118,5.59759126,4.49727,4.88361]
const OP_HD_NWO_B1 = Float32[-5.16836,-5.80832,-4.43822,-5.40092761,-7.50951,-3.84795,-2.70252,-3.51586969,-3.19942952,-2.07667,-2.47605]
const OP_HD_NWO_B2 = Float32[-0.253869,-0.240317,-0.411373,-0.38922036,-0.100000,-0.289213,-0.354756,-0.57665068,-0.38783403,-0.388650,-0.309050]

"organon/htgrowth.f HD_NWO — predicted height from DBH for group `g`."
@inline op_hd_nwo(g::Int, dbh::Float32) = 4.5f0 + fexp(OP_HD_NWO_B0[g] + OP_HD_NWO_B1[g]*fpow(dbh, OP_HD_NWO_B2[g]))

# --- start2.f HD40_NWO HD40PAR(2,3): DF,WH (IEQ 1,2) ------------------------------------------
const OP_HD40_NWO_B0 = Float32[-2.857232223,-2.790360488]
const OP_HD40_NWO_B1 = Float32[-0.393885195,-0.235470605]
const OP_HD40_NWO_B2 = Float32[-0.000521583,-0.002374673]

"organon/start2.f HD40_NWO — predicted height via 40-largest H40/D40 (IEQ 1=DF,2=WH)."
@inline function op_hd40_nwo(ieq::Int, ht40::Float32, d40::Float32, dbh::Float32)
    b0 = OP_HD40_NWO_B0[ieq]; b1 = OP_HD40_NWO_B1[ieq]; b2 = OP_HD40_NWO_B2[ieq]
    exd   = fexp(b0*fpow(dbh, b1 + b2*(ht40 - 4.5f0)))
    exd40 = fexp(b0*fpow(d40, b1 + b2*(ht40 - 4.5f0)))
    return 4.5f0 + (ht40 - 4.5f0)*(exd/exd40)
end

"organon/htgrowth.f LIMIT (VERSION=2) — cap HG under the predicted height-from-DBH curve. `isp` = FIA."
function op_limit(isp::Int32, dbh::Float32, ht::Float32, dg::Float32, hg::Float32)
    if isp == 202 || isp == 263
        a0 = 19.04942539f0; a1 = -0.04484724f0; a2 = 1.0f0
    elseif isp == 17 || isp == 15
        a0 = 16.26279948f0; a1 = -0.04484724f0; a2 = 1.0f0
    elseif isp == 122
        a0 = 17.11482201f0; a1 = -0.04484724f0; a2 = 1.0f0
    elseif isp == 117
        a0 = 14.29011403f0; a1 = -0.04484724f0; a2 = 1.0f0
    elseif isp == 351
        a0 = 60.619859f0; a1 = -1.59138564f0; a2 = 0.496705997f0
    else
        a0 = 15.80319194f0; a1 = -0.04484724f0; a2 = 1.0f0
    end
    ht1 = ht - 4.5f0; ht2 = ht1 + hg; ht3 = ht2 + hg
    dbh2 = dbh + dg; dbh3 = dbh2 + dg
    pht1 = a0*dbh/(1.0f0 - a1*fpow(dbh, a2))
    pht2 = a0*dbh2/(1.0f0 - a1*fpow(dbh2, a2))
    pht3 = a0*dbh3/(1.0f0 - a1*fpow(dbh3, a2))
    phgr1 = (pht2 - pht1 + hg)/2.0f0
    phgr2 = pht2 - ht1
    if ht2 > pht2
        hg = phgr1 < phgr2 ? phgr1 : phgr2
    elseif ht3 > pht3
        hg = phgr1
    end
    hg < 0.0f0 && (hg = 0.0f0)
    return hg
end

# --- mortality.f PM_NWO MPAR(11,7): B0,B1,B2,B3,B4,B5,POW ------------------------------------
const OP_PM_NWO_B0  = Float32[-4.13142,-7.60159,-0.761609,-0.761609,-4.072781265,-6.089598985,-2.976822456,-6.00031085,-2.0,-3.020345211,-1.386294361]
const OP_PM_NWO_B1  = Float32[-1.13736,-0.200523,-0.529366,-0.529366,-0.176433475,-0.245615070,0.0,-0.10490823,-0.5,0.0,0.0]
const OP_PM_NWO_B2  = Float32[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.015,0.0,0.0]
const OP_PM_NWO_B3  = Float32[-0.823305,0.0,-4.74019,-4.74019,-1.729453975,-3.208265570,-6.223250962,-0.99541909,-3.0,-8.467882343,0.0]
const OP_PM_NWO_B4  = Float32[0.0307749,0.0441333,0.0119587,0.0119587,0.0,0.033348079,0.0,0.00912739,0.015,0.013966388,0.0]
const OP_PM_NWO_B5  = Float32[0.00991005,0.00063849,0.00756365,0.00756365,0.012525642,0.013571319,0.0,0.87115652,0.01,0.009461545,0.0]
const OP_PM_NWO_POW = Float32[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]

"""
    op_pm_nwo(g, dbh, cr, si_1, si_2, bal) -> (pm, pow)

organon/mortality.f PM_NWO — the per-group linear predictor (DIFFERENT form per group: DF uses
√DBH & CR^0.25; GF uses BAL/DBH; WH/RC use SI_2; WO uses log(BAL+5)). `si_1`/`si_2` = SITE−4.5.
"""
@inline function op_pm_nwo(g::Int, dbh::Float32, cr::Float32, si_1::Float32, si_2::Float32, bal::Float32)
    b0 = OP_PM_NWO_B0[g]; b1 = OP_PM_NWO_B1[g]; b2 = OP_PM_NWO_B2[g]; b3 = OP_PM_NWO_B3[g]
    b4 = OP_PM_NWO_B4[g]; b5 = OP_PM_NWO_B5[g]; pow = OP_PM_NWO_POW[g]
    if g == 1
        pm = b0 + b1*fpow(dbh, 0.5f0) + b3*fpow(cr, 0.25f0) + b4*(si_1 + 4.5f0) + b5*bal
    elseif g == 2
        pm = b0 + b1*dbh + b4*(si_1 + 4.5f0) + b5*(bal/dbh)
    elseif g == 3 || g == 4
        pm = b0 + b1*dbh + b2*dbh*dbh + b3*cr + b4*(si_2 + 4.5f0) + b5*bal
    elseif g == 8
        pm = b0 + b1*dbh + b2*dbh*dbh + b3*cr + b4*(si_1 + 4.5f0) + b5*flog(bal + 5.0f0)
    else
        pm = b0 + b1*dbh + b2*dbh*dbh + b3*cr + b4*(si_1 + 4.5f0) + b5*bal
    end
    return pm, pow
end

# =============================================================================
# NWO structural wrappers (SSTATS/SSUM/CRNCLO — call the NWO MCW/crown-width leaves) + SUBMAX
# =============================================================================

"organon/statsorg.f SSTATS (VERSION=2) — SBA/BAL(500)/BALL(51)/CCFL/CCFLL via MCW_NWO."
function op_sstats(spgrp::Vector{Int32}, dbh::Vector{Float32}, ht::Vector{Float32},
                   expan::Vector{Float32}, ntrees::Int)
    bal = zeros(Float32, 500); ball = zeros(Float32, 51)
    ccfl = zeros(Float32, 500); ccfll = zeros(Float32, 51)
    sba = 0.0f0; tpa = 0.0f0; sccf = 0.0f0
    @inbounds for i in 1:ntrees
        expan[i] < 0.0001f0 && continue
        g = Int(spgrp[i]); d = dbh[i]; h = ht[i]; ex = expan[i]
        ba = d*d*ex*0.005454154f0
        sba += ba; tpa += ex
        mcw = op_mcw_nwo(g, d, h)
        ccf = 0.001803f0*(mcw*mcw)*ex
        sccf += ccf
        if d > 50.0f0
            l = trunc(Int, d - 49.0f0); l > 52 && (l = 52)
            for k in 1:500; ccfl[k] += ccf; bal[k] += ba; end
            for k in 1:(l-1); ccfll[k] += ccf; ball[k] += ba; end
        else
            l = trunc(Int, d*10.0f0 + 0.5f0)
            for k in 1:(l-1); ccfl[k] += ccf; bal[k] += ba; end
        end
    end
    return sba, bal, ball, ccfl, ccfll, tpa, sccf
end

"organon/diamcal.f SSUM (N=2 ending stats) — (sba, ccfl, ccfll) via MCW_NWO."
function op_ssum(spgrp::Vector{Int32}, dbh::Vector{Float32}, ht::Vector{Float32},
                 expan::Vector{Float32}, ntrees::Int, npts::Int)
    ccfl = zeros(Float32, 500); ccfll = zeros(Float32, 51)
    sba = 0.0f0
    @inbounds for i in 1:ntrees
        g = Int(spgrp[i]); ex = expan[i]/Float32(npts)
        ex <= 0.0f0 && continue
        ba = (dbh[i]*dbh[i]*ex)*0.005454154f0
        sba += ba
        mcw = op_mcw_nwo(g, dbh[i], ht[i])
        ccf = 0.001803f0*mcw*mcw*expan[i]/Float32(npts)
        if dbh[i] > 50.0f0
            l = trunc(Int, dbh[i] - 49.0f0); l > 52 && (l = 52)
            for k in 1:500; ccfl[k] += ccf; end
            for k in 1:(l-1); ccfll[k] += ccf; end
        else
            l = trunc(Int, dbh[i]*10.0f0 + 0.5f0)
            for k in 1:(l-1); ccfl[k] += ccf; end
        end
    end
    return sba, ccfl, ccfll
end

"organon/crngrow.f CALC_CC — accumulate one tree's crown area into CCH(1..40) via CW_NWO."
function op_calc_cc!(g::Int, hlcw::Float32, lcw::Float32, ht::Float32, dbh::Float32,
                     hcb::Float32, expan::Float32, cch::Vector{Float32})
    if hcb > hlcw
        xhlcw = hcb
        xlcw = op_cw_nwo(g, hlcw, lcw, ht, dbh, xhlcw)
    else
        xhlcw = hlcw
        xlcw = lcw
    end
    @inbounds for ii in 40:-1:1
        l = ii - 1
        xl = Float32(l)*(cch[41]/40.0f0)
        if xl <= xhlcw
            cw = xlcw
        elseif xl < ht
            cw = op_cw_nwo(g, hlcw, lcw, ht, dbh, xl)
        else
            cw = 0.0f0
        end
        cch[ii] += (cw*cw)*(0.001803f0*expan)
    end
    return
end

"organon/crngrow.f CRNCLO (IND=0, CTMUL=0, SCR=0) — crown-closure profile CCH(41) via NWO crown."
function op_crnclo(spgrp::Vector{Int32}, dbh::Vector{Float32}, ht::Vector{Float32},
                   cr::Vector{Float32}, expan::Vector{Float32}, ntrees::Int)
    cch = zeros(Float32, 41)
    hmax = ht[1]
    @inbounds for i in 2:ntrees; ht[i] > hmax && (hmax = ht[i]); end
    cch[41] = hmax
    @inbounds for i in 1:ntrees
        g = Int(spgrp[i]); d = dbh[i]; h = ht[i]; c = cr[i]
        cl = c*h; hcb = h - cl
        mcw = op_mcw_nwo(g, d, h)
        lcw = op_lcw_nwo(g, mcw, c, 0.0f0, d, h)
        hlcw = op_hlcw_nwo(g, h, c, 0.0f0)
        op_calc_cc!(g, hlcw, lcw, h, d, hcb, expan[i], cch)
    end
    return cch
end

"organon/submax.f SUBMAX (VERSION=2) — max size-density line A1/A2."
function op_submax(spgrp::Vector{Int32}, dbh::Vector{Float32}, expan::Vector{Float32},
                   ntrees::Int, msdi_1::Float32, msdi_2::Float32, msdi_3::Float32)
    a2 = 0.62305f0; kb = 0.005454154f0
    tempa1 = msdi_1 > 0.0f0 ? (flog(10.0f0) + a2*flog(msdi_1)) : 6.19958f0
    bagrp = zeros(Float32, 18)
    @inbounds for i in 1:ntrees
        bagrp[Int(spgrp[i])] += kb*(dbh[i]*dbh[i])*expan[i]
    end
    totba = bagrp[1] + bagrp[2] + bagrp[3]
    pdf = ptf = pwh = 0.0f0
    if totba > 0.0f0
        pdf = bagrp[1]/totba; ptf = bagrp[2]/totba; pwh = bagrp[3]/totba
    end
    tfmod = msdi_2 > 0.0f0 ? (flog(10.0f0) + a2*flog(msdi_2))/tempa1 : 1.03481817f0
    ocmod = msdi_3 > 0.0f0 ? (flog(10.0f0) + a2*flog(msdi_3))/tempa1 : 1.014293245f0
    if pdf >= 0.5f0
        a1mod = 1.0f0
    elseif pwh >= 0.5f0
        a1mod = ocmod
    elseif ptf >= 0.6666667f0
        a1mod = tfmod
    else
        a1mod = pdf + ocmod*pwh + tfmod*ptf
    end
    a1mod <= 0.0f0 && (a1mod = 1.0f0)
    return tempa1*a1mod, a2
end

# =============================================================================
# NWO diameter / height / crown / mortality passes (off the /ORGANON/ buffer)
# =============================================================================

"""
    op_dg_nwo_pass(buf, spgrp; si_1, msdi_1, msdi_2, msdi_3, cyclg=0)
        -> (dgro, sba1, bal1, ball1, a1, a2)

organon/grow.f growth-1 (VERSION=2) — the DG_NWO diameter pass off the buffer. SITE = SI_2 if the
ORGANON FIA species is WH (263) else SI_1 (diagro.f:45-50). Fert/thin = 1.0 on cyc0.
"""
function op_dg_nwo_pass(buf::OrganonBuffer, spgrp::Vector{Int32}; si_1::Float32, si_2::Float32,
        msdi_1::Float32=0f0, msdi_2::Float32=0f0, msdi_3::Float32=0f0, cyclg::Int=0)
    n = buf.ntrees
    sba1, bal1, ball1, _, _, _, _ = op_sstats(spgrp, buf.dbh1, buf.ht1or, buf.expan1, n)
    a1, a2 = op_submax(spgrp, buf.dbh1, buf.expan1, n, msdi_1, msdi_2, msdi_3)
    dgro = zeros(Float32, n)
    @inbounds for i in 1:n
        buf.expan1[i] <= 0.0f0 && continue
        g = Int(spgrp[i])
        site = buf.species[i] == 263 ? si_2 : si_1
        sbal1 = oc_get_bal(buf.dbh1[i], ball1, bal1)
        dgro[i] = op_dg_nwo(g, buf.dbh1[i], buf.cr1[i], site, sbal1, sba1)   # ADJ·CRADJ inside; CALIB3=fert=thin=1
    end
    return dgro, sba1, bal1, ball1, a1, a2
end

"""
    op_hg_nwo_pass(buf, dgro, spgrp; si_1, si_2, cyclg=0, calib1) -> hgro

organon/grow.f growth-2 (VERSION=2) — HTGRO1 (B_HG/F_HG + HG_NWO, groups 1..3) / HTGRO2 (HD_NWO
ratio, groups >3) off the buffer + C3 dgro. `calib1` = ACALIB(1,·) (HTGRO2 height calibration).
"""
function op_hg_nwo_pass(buf::OrganonBuffer, dgro::Vector{Float32}, spgrp::Vector{Int32};
        si_1::Float32, si_2::Float32, cyclg::Int=0, calib1::Union{Nothing,Vector{Float32}}=nothing)
    n = buf.ntrees
    cch = op_crnclo(spgrp, buf.dbh1, buf.ht1or, buf.cr1, buf.expan1, n)
    hgro = zeros(Float32, n)
    ib = 3
    @inbounds for i in 1:n
        buf.expan1[i] <= 0.0f0 && continue
        g = Int(spgrp[i]); dbh = buf.dbh1[i]; ht = buf.ht1or[i]; cr = buf.cr1[i]
        if g <= ib
            # TCCH — crown competition at this tree's height (htgrowth.f:67-77)
            xi = 40.0f0*(ht/cch[41])
            ii = trunc(Int, xi) + 2
            xxi = Float32(ii) - 1.0f0
            if ht >= cch[41]
                tcch = 0.0f0
            elseif ii == 41
                tcch = cch[40]*(40.0f0 - xi)
            else
                tcch = cch[ii] + (cch[ii-1] - cch[ii])*(xxi - xi)
            end
            gp = 5.0f0
            if g == 3
                _, phtgro = op_f_hg(si_2 + 4.5f0, ht, gp)   # WH — Flewelling
            else
                _, phtgro = op_b_hg(si_1 + 4.5f0, ht, gp)   # DF/GF — Bruce
            end
            hg = op_hg_nwo(g, phtgro, cr, tcch)
            # HG_FERT/HG_THIN (VERSION≤3) = 1.0 with no fert/thin — reuse OC's (ispgrp==1 branch)
            hg *= oc_hg_thin(cyclg, g, 0f0, zeros(Float32,5), zeros(Float32,5))
            hg *= oc_hg_fert(cyclg, g, si_1, zeros(Float32,5), zeros(Float32,5))
            hgro[i] = op_limit(buf.species[i], dbh, ht, dgro[i], hg)
        else
            c1 = calib1 === nothing ? 1.0f0 : calib1[g]
            hgro[i] = op_htgro2(buf.species[i], g, dbh, dgro[i], ht, c1)
        end
    end
    return hgro
end

"organon/htgrowth.f HTGRO2 (VERSION=2, minor species g>3) — HD_NWO height-ratio (no red alder here)."
@inline function op_htgro2(isp::Int32, g::Int, dbh_start::Float32, dgro::Float32, ht::Float32, calib1::Float32)
    dbh_end = dbh_start + dgro
    prdht1 = 4.5f0 + calib1*(op_hd_nwo(g, dbh_start) - 4.5f0)
    prdht2 = 4.5f0 + calib1*(op_hd_nwo(g, dbh_end) - 4.5f0)
    return (prdht2/prdht1)*ht - ht
end

# F_HG (Flewelling WH) — NOT in the S248112 stand (no WH); raise so a WH stand flags the gap.
function op_f_hg(si::Float32, ht::Float32, gp::Float32)
    error("OP HG_NWO: Flewelling F_HG (western hemlock, group 3) is UNPORTED (no WH in the S248112 " *
          "reference stand). Port SITECV_F (organon/whphg.f) before running a WH-bearing stand.")
end

"""
    op_mortal_nwo(buf, dgro, spgrp, bal1, ball1, a1, a2; si_1, si_2, cyclg=0, mort=true) -> deadexp

organon/mortality.f MORTAL_RUN (VERSION=2) — whole-stand ORGANON mortality via PM_NWO. Same
VERSION≤3 driver structure as SWO (QUAD1, RDCC=0.60, cyc0 init) — only the PM call differs.
"""
function op_mortal_nwo(buf::OrganonBuffer, dgro::Vector{Float32}, spgrp::Vector{Int32},
        bal1::Vector{Float32}, ball1::Vector{Float32}, a1::Float32, a2::Float32;
        si_1::Float32, si_2::Float32, cyclg::Int=0, mort::Bool=true)
    n = buf.ntrees
    ib = 3; a3 = 14.39533971f0; rdcc = 0.60f0; kb = 0.005454154f0
    dbh0 = buf.dbh1; ht0 = buf.ht1or; cr0 = buf.cr1
    expan = copy(buf.expan1)
    deadexp = zeros(Float32, n)
    pmk = zeros(Float32, n); pow = ones(Float32, n)

    stba = 0.0f0; stn = 0.0f0
    @inbounds for i in 1:n
        stba += (dbh0[i]*dbh0[i])*kb*expan[i]; stn += expan[i]
    end
    sqmda = sqrt(stba/(kb*stn))
    rd = stn/fexp(a1/a2 - flog(sqmda)/a2)
    og1 = oc_oldgro(spgrp, dbh0, ht0, expan, dgro, zeros(Float32,n), deadexp, n, ib, 0.0f0)

    @inbounds for i in 1:n
        expan[i] <= 0.0f0 && continue
        g = Int(spgrp[i])
        fertadj = oc_pm_fert(g, cyclg, zeros(Float32,5), zeros(Float32,5))
        sbal1 = oc_get_bal(dbh0[i], ball1, bal1)
        pm, p = op_pm_nwo(g, dbh0[i], cr0[i], si_1, si_2, sbal1)
        pow[i] = p; pmk[i] = pm + fertadj
    end

    na = 0.0f0; baa = 0.0f0
    @inbounds for i in 1:n
        cr = cr0[i]
        cradj = cr <= 0.17f0 ? (1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0))) : 1.0f0
        xpm = 1.0f0/(1.0f0 + fexp(-pmk[i]))
        ps = fpow(1.0f0 - xpm, pow[i])
        pm = 1.0f0 - ps*cradj
        na  += expan[i]*(1.0f0 - pm)
        baa += kb*((dbh0[i]+dgro[i])*(dbh0[i]+dgro[i]))*expan[i]*(1.0f0 - pm)
    end

    apply_base = true
    if mort
        qmda = sqrt(baa/(kb*na))
        ind = 0; a1max = a1; no = 0.0f0
        rda = na/fexp(a1/a2 - flog(qmda)/a2)
        if cyclg == 0
            if rd >= 1.0f0
                a1max = rda > rd ? (flog(sqmda) + a2*flog(stn)) : (flog(qmda) + a2*flog(na))
                ind = 1; a1max < a1 && (a1max = a1)
            else
                ind = 0
                rd > rdcc && (no = stn*fpow(flog(rd)/flog(rdcc), -1.0f0/a3))
            end
        end
        qmdp = (ind == 0 && no > 0.0f0) ? oc_quad1(na, no, rdcc, a1) : fexp(a1max - a2*flog(na))
        if rd <= rdcc || qmdp > qmda
            apply_base = true
        else
            apply_base = false
            kr1 = 0.0f0
            for kk in 1:7
                nk = 10.0f0/fpow(10.0f0, Float32(kk))
                while true
                    kr1 += nk
                    naa = 0.0f0; baaa = 0.0f0
                    @inbounds for i in 1:n
                        expan[i] < 0.001f0 && continue
                        cr = cr0[i]
                        cradj = cr <= 0.17f0 ? (1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0))) : 1.0f0
                        xpm = 1.0f0/(1.0f0 + fexp(-(kr1 + pmk[i])))
                        ps = fpow(1.0f0 - xpm, pow[i]); pm = 1.0f0 - ps*cradj
                        naa  += expan[i]*(1.0f0 - pm)
                        baaa += kb*fpow(dbh0[i]+dgro[i], 2.0f0)*expan[i]*(1.0f0 - pm)
                    end
                    qmda = sqrt(baaa/(kb*naa))
                    qmdp = ind == 0 ? oc_quad1(naa, no, rdcc, a1) : fexp(a1max - a2*flog(naa))
                    if qmdp >= qmda
                        kr1 -= nk; break
                    end
                end
            end
            @inbounds for i in 1:n
                if expan[i] <= 0.0f0
                    deadexp[i] = 0.0f0; expan[i] = 0.0f0
                else
                    cr = cr0[i]
                    cradj = cr <= 0.17f0 ? (1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0))) : 1.0f0
                    xpm = 1.0f0/(1.0f0 + fexp(-(kr1 + pmk[i])))
                    ps = fpow(1.0f0 - xpm, pow[i]); pm = 1.0f0 - ps*cradj
                    deadexp[i] = expan[i]*pm; expan[i] = expan[i]*(1.0f0 - pm)
                end
            end
        end
    end

    if apply_base
        @inbounds for i in 1:n
            cr = cr0[i]
            cradj = cr <= 0.17f0 ? (1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0))) : 1.0f0
            xpm = 1.0f0/(1.0f0 + fexp(-pmk[i]))
            ps = fpow(1.0f0 - xpm, pow[i]); pm = 1.0f0 - ps*cradj
            deadexp[i] = expan[i]*pm; expan[i] = expan[i]*(1.0f0 - pm)
        end
    end
    return deadexp
end

"""
    op_cr_nwo(buf, dgro, hgro, spgrp, deadexp; si_1, si_2, calib2, cyclg=0) -> cr2

organon/crngrow.f CROWGRO (VERSION=2) — crown recession with the NWO CALIB(2) multiplier on
PCR1/PCR2 (crngrow.f:73,96-101). Runs after DG/HG advance DBH/HT and after mortality reduces EXPAN.
"""
function op_cr_nwo(buf::OrganonBuffer, dgro::Vector{Float32}, hgro::Vector{Float32},
        spgrp::Vector{Int32}, deadexp::Vector{Float32}; si_1::Float32, si_2::Float32,
        calib2::Vector{Float32}, cyclg::Int=0)
    n = buf.ntrees; ib = 3
    gdbh = Vector{Float32}(undef, n); ght = Vector{Float32}(undef, n); surv = Vector{Float32}(undef, n)
    @inbounds for i in 1:n
        gdbh[i] = buf.dbh1[i] + dgro[i]
        ght[i]  = buf.ht1or[i] + hgro[i]
        surv[i] = buf.expan1[i] - deadexp[i]
    end
    sba1, _, _, ccfl1, ccfll1, _, _ = op_sstats(spgrp, buf.dbh1, buf.ht1or, buf.expan1, n)
    sba2, _, _, ccfl2, ccfll2, _, _ = op_sstats(spgrp, gdbh, buf.ht1or, surv, n)
    og1 = oc_oldgro(spgrp, gdbh, ght, surv, dgro, hgro, deadexp, n, ib, -1.0f0)
    og2 = oc_oldgro(spgrp, gdbh, ght, surv, dgro, hgro, deadexp, n, ib, 0.0f0)

    cr2 = zeros(Float32, n)
    @inbounds for i in 1:n
        g = Int(spgrp[i])
        pht = buf.ht1or[i]; pdbh = buf.dbh1[i]
        sccfl1 = oc_get_ccfl(pdbh, ccfll1, ccfl1)
        hcb1 = op_hcb_nwo(g, pht, pdbh, sccfl1, sba1, si_1, si_2, og1)
        pcr1 = calib2[g]*(1.0f0 - hcb1/pht)             # NWO: CALIB(2) applied (crngrow.f:73)
        phcb1 = (1.0f0 - pcr1)*pht
        ht = ght[i]; dbh = gdbh[i]
        sccfl2 = oc_get_ccfl(dbh, ccfll2, ccfl2)
        hcb2 = op_hcb_nwo(g, ht, dbh, sccfl2, sba2, si_1, si_2, og2)
        maxhcb = op_maxhcb_nwo(g, ht, sccfl2)
        pcr2 = calib2[g]*(1.0f0 - hcb2/ht)              # NWO: CALIB(2) applied (crngrow.f:100)
        phcb2 = (1.0f0 - pcr2)*ht
        hcbg = phcb2 - phcb1
        hcbg < 0.0f0 && (hcbg = 0.0f0)
        ahcb1 = (1.0f0 - buf.cr1[i])*pht                # actual HCB (start)
        shcb1 = pht                                     # SCR=0 ⇒ SHCB1=(1-0)*pht=pht
        ahcb2 = ahcb1 + hcbg
        shcb2 = shcb1 + hcbg
        if ahcb1 > shcb1                                # false for cr>0 (ahcb1<pht=shcb1) — SCR path, inert
            cr2[i] = ahcb1 > shcb2 ? (1.0f0 - ahcb1/ht) : (1.0f0 - shcb2/ht)
        else
            if ahcb1 >= maxhcb
                cr2[i] = 1.0f0 - ahcb1/ht
            elseif ahcb2 >= maxhcb
                cr2[i] = 1.0f0 - maxhcb/ht
            else
                cr2[i] = 1.0f0 - ahcb2/ht
            end
        end
    end
    return cr2
end

# =============================================================================
# C1 marshalling + C2 PREPARE (VERSION=2) + C7 EXECUTE — off the /ORGANON/ buffer
# =============================================================================

"op/dgdriv.f — the 11 valid ORGANON species (IORG=1): GF,DF,RC,WH,BM,RA,MA,WO,PY,DG,WI."
const OP_ORGANON_VALID = (3, 16, 18, 19, 21, 22, 23, 28, 33, 34, 37)
"op/dgdriv.f — the big-6 stand-gate species (NBIG6 counts these): GF, DF."
const OP_ORGANON_BIG6 = (3, 16)

"""
    op_build_organon_buffer!(s) -> OrganonBuffer

Port of op/dgdriv.f:211-258 — per-tree IORG eligibility (`HT>4.5 AND DBH>=0.1 AND sp∈valid`), the
big-6 stand gate (NBIG6 counts GF/DF), and the /ORGANON/ input buffer fill. Reuses OC's OrganonBuffer.
"""
function op_build_organon_buffer!(s::StandState)
    t = s.trees; n = t.n
    buf = OrganonBuffer(n); buf.ntrees = n
    nbig6 = 0
    @inbounds for i in 1:n
        sp = Int(t.species[i])
        if t.height[i] > 4.5f0 && t.dbh[i] >= 0.1f0
            (sp in OP_ORGANON_BIG6) && (nbig6 += 1)
            buf.iorg[i] = (sp in OP_ORGANON_VALID) ? Int32(1) : Int32(0)
        else
            buf.iorg[i] = Int32(0)
        end
        buf.species[i] = op_organon_fia(sp)         # ORGSPC — surrogate for non-valid
    end
    buf.nbig6 = nbig6
    if nbig6 == 0
        @inbounds for i in 1:n; buf.iorg[i] = Int32(0); buf.mortexp[i] = 0f0; end
        buf.nvalid = 0; buf.runs = false; return buf
    end
    buf.runs = true
    nvalid = 0
    @inbounds for i in 1:n
        buf.treeno[i] = Int32(i); buf.ptno[i] = t.plot_id[i]
        d = t.dbh[i]; d < 0.1f0 && (d = 0.1f0); buf.dbh1[i] = d
        h = t.height[i]; (h > 0f0 && h < 4.6f0) && (h = 4.6f0); buf.ht1or[i] = h
        buf.cr1[i] = Float32(t.crown_pct[i])/100f0
        buf.expan1[i] = t.tpa[i]; buf.user[i] = t.special[i]
        nvalid += Int(buf.iorg[i])
    end
    buf.nvalid = nvalid
    return buf
end

"""
    op_prepare_nwo(species, dbh0, ht0, cr0, expan0, ntrees, npts, stage, bhage, si_1, si_2;
                   acalib0=nothing) -> OrganonCalib

organon/prepare.f PREPARE (VERSION=2) + the op/cratet.f TMPCAL→ACALIB load. Mirrors OC's validated
`organon_prepare_swo` with the NWO leaves: IB=3, NSPN=11, SI conversion DF↔WH, HDCALIB via
B_H40/HD40_NWO/HD_NWO (PDF/PWH dispatch), CRCALIB via HCB_NWO. Produces ACALIB(1,·) (height) and
ACALIB(2,·) (crown — consumed by CROWGRO). `si_1`/`si_2` are the RAW site indices (SITE_1/SITE_2).
"""
function op_prepare_nwo(species::Vector{Int32}, dbh0::Vector{Float32}, ht0::Vector{Float32},
        cr0::Vector{Float32}, expan0::Vector{Float32}, ntrees::Int, npts::Int, stage::Int,
        bhage::Int, si_1::Float32, si_2::Float32; acalib0::Union{Nothing,Matrix{Float32}}=nothing)
    ib = 3; nspn = 11
    dbh = copy(dbh0); ht = copy(ht0); cr = copy(cr0); expan = copy(expan0)
    spgrp = Vector{Int32}(undef, ntrees)
    missht = false; misscr = false
    @inbounds for i in 1:ntrees
        spgrp[i] = op_spgroup_nwo_fia(species[i])
        ht[i] <= 0.0f0 && (missht = true)
        cr[i] <= 0.0f0 && (misscr = true)
    end
    # SI conversion (prepare.f:363-367, VERSION=2): DF↔WH
    if si_1 <= 0.0f0 && si_2 > 0.0f0
        si_1 = 0.480f0 + 1.110f0*si_2
    elseif si_2 <= 0.0f0
        si_2 = -0.432f0 + 0.899f0*si_1
    end
    tmpcal = ones(Float32, 3, 18)

    # --- HDCALIB (VERSION=2): B_H40 HT40, HD40_NWO/HD_NWO predicted heights ---
    entht = zeros(Int, 18); entdbh = zeros(Int, 18)
    yxs = zeros(Float32, 18); xss = zeros(Float32, 18); yss = zeros(Float32, 18)
    ptrht = zeros(Float32, ntrees); d40 = 0.0f0; ht40 = 0.0f0
    pdf = ptf = ppp = pwh = pra = 0.0f0
    even = true                                       # OP even-aged (INDS(4)=1)
    if even
        age1 = Float32(bhage)
        (pdf, ptf, ppp, pwh, pra) = oc_spmix(species, dbh, expan, ntrees, npts)
        d40 = oc_dforty(spgrp, dbh, expan, ntrees, npts, ib)
        if pdf >= 0.80f0
            ht40 = op_b_h40(age1, si_1)
        elseif pwh >= 0.80f0
            error("OP PREPARE: PWH≥0.80 F_H40 (Flewelling WH HT40) UNPORTED (no WH in S248112).")
        end
    end
    @inbounds for i in 1:ntrees
        g = Int(spgrp[i]); entdbh[g] += 1
        ht[i] <= 0.0f0 && continue
        if pdf >= 0.80f0 && g == 1
            ptrht[i] = op_hd40_nwo(1, ht40, d40, dbh[i])
        elseif pwh >= 0.80f0 && g == 3
            ptrht[i] = op_hd40_nwo(2, ht40, d40, dbh[i])
        else
            ptrht[i] = op_hd_nwo(g, dbh[i])
        end
        y = ht[i] - 4.5f0; x = ptrht[i] - 4.5f0; wt = dbh[i]
        yxs[g] += y*x/wt; xss[g] += x*x/wt; yss[g] += y*y/wt; entht[g] += 1
    end
    @inbounds for g in 1:nspn
        (entdbh[g] == 0 || entht[g] < 2) && continue
        tmpcal[1, g] = oc_caltst(yxs[g], xss[g], yss[g], entht[g])
    end
    # --- PRDHT: impute missing heights ---
    if missht
        @inbounds for i in 1:ntrees
            g = Int(spgrp[i])
            (dbh[i] <= 0.0f0 || ht[i] > 0.0f0) && continue
            if pdf >= 0.80f0 && g == 1
                ptrht[i] = op_hd40_nwo(1, ht40, d40, dbh[i])
            elseif pwh >= 0.80f0 && g == 3
                ptrht[i] = op_hd40_nwo(2, ht40, d40, dbh[i])
            else
                ptrht[i] = op_hd_nwo(g, dbh[i])
            end
            newht = 4.5f0 + tmpcal[1, g]*(ptrht[i] - 4.5f0)
            newht < 4.6f0 && (newht = 4.6f0)
            ht[i] = newht
        end
    end
    # --- CRCALIB (VERSION=2): HCB_NWO predicted crown ---
    sba, ccfl, ccfll = op_ssum(spgrp, dbh, ht, expan, ntrees, npts)
    og = oc_oldgrowth(spgrp, dbh, ht, expan, ntrees, npts, ib)
    entcr = zeros(Int, 18)
    cyxs = zeros(Float32, 18); cxss = zeros(Float32, 18); cyss = zeros(Float32, 18)
    xsi_1 = si_1 - 4.5f0; xsi_2 = si_2 - 4.5f0
    @inbounds for i in 1:ntrees
        cr[i] <= 0.0f0 && continue
        g = Int(spgrp[i])
        sccfl = oc_get_ccfl(dbh[i], ccfll, ccfl)
        phcb = op_hcb_nwo(g, ht[i], dbh[i], sccfl, sba, xsi_1, xsi_2, og)
        pcr = 1.0f0 - phcb/ht[i]
        cyxs[g] += pcr*cr[i]; cxss[g] += pcr*pcr; cyss[g] += cr[i]*cr[i]; entcr[g] += 1
    end
    @inbounds for g in 1:nspn
        (entdbh[g] == 0 || entcr[g] < 2) && continue
        tmpcal[2, g] = oc_caltst(cyxs[g], cxss[g], cyss[g], entcr[g])
    end
    # --- PRDCR: impute missing crown ratios ---
    if misscr
        @inbounds for i in 1:ntrees
            cr[i] > 0.0f0 && continue
            g = Int(spgrp[i]); dbh[i] <= 0.0f0 && continue
            sccfl = oc_get_ccfl(dbh[i], ccfll, ccfl)
            phcb = op_hcb_nwo(g, ht[i], dbh[i], sccfl, sba, xsi_1, xsi_2, og)
            newcr = (1.0f0 - phcb/ht[i])*tmpcal[2, g]
            newcr > 1.0f0 && (newcr = 1.0f0); newcr < 0.05f0 && (newcr = 0.05f0)
            cr[i] = newcr
        end
    end
    # --- op/cratet.f TMPCAL → ACALIB load ---
    acalib = acalib0 === nothing ? ones(Float32, 3, 18) : copy(acalib0)
    @inbounds for i in 1:18, j in 1:3
        if acalib[j, i] == 1.0f0 && (tmpcal[j, i] > 0.0f0 && tmpcal[j, i] != 1.0f0)
            acalib[j, i] = tmpcal[j, i]
        end
    end
    return OrganonCalib(acalib, tmpcal, dbh, ht, cr, expan, 0)
end

"organon/prepare.f SPGROUP_EDIT (VERSION=2) — ORGANON FIA code → NWO group 1..11 (or -9999)."
@inline function op_spgroup_nwo_fia(fia::Integer)
    @inbounds for g in 1:length(OP_SCODE2)
        OP_SCODE2[g] == fia && return Int32(g)
    end
    return Int32(-9999)
end

"""
    op_execute_nwo(buf, isp_fvs; si_1, si_2, msdi, calib1, calib2, cyclg=0) -> OrganonGrowth

organon/execute2.f EXECUTE + grow.f GROW (VERSION=2) off the /ORGANON/ buffer, in the faithful order
DG_NWO → HG_NWO → MORTAL_RUN → CROWGRO, then the op/dgdriv.f DGRO→DDS copy-back. `si_1`/`si_2` are
SITE−4.5. Deterministic (DGSD=0). Reuses OC's OrganonGrowth container.
"""
function op_execute_nwo(buf::OrganonBuffer, isp_fvs::AbstractVector{<:Integer}; si_1::Float32,
        si_2::Float32, msdi::Float32=0f0, msdi_1::Float32=NaN32, msdi_2::Float32=NaN32,
        msdi_3::Float32=NaN32, calib1::Vector{Float32}=ones(Float32,11),
        calib2::Vector{Float32}=ones(Float32,11), cyclg::Int=0)
    # MSDI_1/2/3 = RVARS(3/4/5) = SDIDEF(16=DF)/SDIDEF(3=GF)/SDIDEF(19=WH) (op/sitset.f:325-327). Back-compat:
    # a single `msdi` fills all three (the validated 27-tree unit test path); pass msdi_1/2/3 for the real fan.
    isnan(msdi_1) && (msdi_1 = msdi); isnan(msdi_2) && (msdi_2 = msdi); isnan(msdi_3) && (msdi_3 = msdi)
    n = buf.ntrees
    spgrp = Int32[op_spgroup_nwo(Int(isp)) for isp in isp_fvs]
    dgro, sba1, bal1, ball1, a1, a2 = op_dg_nwo_pass(buf, spgrp; si_1=si_1, si_2=si_2,
                                          msdi_1=msdi_1, msdi_2=msdi_2, msdi_3=msdi_3, cyclg=cyclg)
    hgro = op_hg_nwo_pass(buf, dgro, spgrp; si_1=si_1, si_2=si_2, cyclg=cyclg, calib1=calib1)
    deadexp = op_mortal_nwo(buf, dgro, spgrp, bal1, ball1, a1, a2; si_1=si_1, si_2=si_2, cyclg=cyclg)
    cr2 = op_cr_nwo(buf, dgro, hgro, spgrp, deadexp; si_1=si_1, si_2=si_2, calib2=calib2, cyclg=cyclg)
    dds = zeros(Float32, n)
    @inbounds for i in 1:n
        buf.iorg[i] == 1 || continue
        _, _, dds[i] = op_organon_dds(Int(isp_fvs[i]), buf.dbh1[i], dgro[i])
    end
    return OrganonGrowth(dgro, hgro, cr2, deadexp, dds, spgrp, sba1, a1, a2)
end
