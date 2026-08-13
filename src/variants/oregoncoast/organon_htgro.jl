# =============================================================================
# organon_htgro.jl — OC (Oregon Coast) ORGANON SWO height growth (chunk C4).
#
# Ported from (ORGANON edition SWO / VERSION=1 path only):
#   • organon/htgrowth.f — HTGRO1 (big-6 height-growth driver), HS_HG (Hann-Scrivani potential
#                          height growth + growth-effective age), HG_SWO (the height-increment
#                          modifier equation), LIMIT (DG-based HG cap), HG_FERT / HG_THIN.
#   • organon/crngrow.f  — CRNCLO (crown-closure profile CCH), CALC_CC, LCW_SWO (largest crown
#                          width), HLCW_SWO (height to largest crown width), CW_SWO (crown width
#                          above LCW). MCW_SWO is in organon_setup.jl (C2). These build the TCCH
#                          crown-competition value HG_SWO consumes.
#   • oc/htgf.f:95-101   — the FVS-side copy-back: for IORG=1 trees `HTG(I)=SCALE·XHT·HGRO(I)·
#                          EXP(HTCON(ISPC))`. On cyc0 defaults SCALE=FINT/YR=1, XHT=1, HTCON=0, so
#                          HTG(I)=HGRO(I) — the C4 validation target.
#
# SCOPE: the ORGANON SWO height-growth core — the deterministic (DGSD=0) per-tree HGRO that FVS
# copies into HTG for the valid ORGANON trees, bypassing the native OC HTGF. HTGRO1 (big-6, species
# groups 1..5) is the validated path; the 'other'-species HTGRO2 (minor ORGANON species, HD-ratio
# form) is ported separately when a stand exercises it. VERSION=1 (SWO) ONLY; NWO(2)/SMC(3)/RAP(4)
# (OP) raise a clear error. HG_SWO depends on the C3 DGRO (LIMIT caps HG using the 5-yr diameter
# growth), so `organon_hg_swo` takes the DGRO vector from `organon_dg_swo`.
#
# MEASURED bit-exact vs the live FVSoc_clean oracle (scoped `DEBUG 1 / DGDRIV HTGF`, stand S248112 /
# ocmin, 17 valid ORGANON trees) — the `HTGF ORGANON … HGRO …` dump. See docs/OC_VARIANT_PORT_AUDIT.md.
# =============================================================================

# --- HG_SWO HGPAR(5,8): P1..P8 (big-5 conifer groups DF,GW,PP,SP,IC only) ----------------------
const OC_HG_SWO_P1 = Float32[1.0,1.0,1.0,1.0,1.0]
const OC_HG_SWO_P2 = Float32[-0.02457621,-0.01453250,-0.14889850,-0.14889850,-0.01453250]
const OC_HG_SWO_P3 = Float32[-0.00407303,-0.00407303,-0.00322752,-0.00678955,-0.00637434]
const OC_HG_SWO_P4 = Float32[1.0,1.0,1.0,1.0,1.0]
const OC_HG_SWO_P5 = Float32[2.89556338,7.69023575,0.92071847,0.92071847,1.27228638]
const OC_HG_SWO_P6 = Float32[2.0,2.0,2.0,2.0,2.0]
const OC_HG_SWO_P7 = Float32[0.0,0.0,0.0,0.0,0.0]
const OC_HG_SWO_P8 = Float32[1.0,1.0,1.0,1.0,1.0]

# --- crown-width tables (18 groups) ------------------------------------------------------------
# organon/crngrow.f LCW_SWO LCWPAR(18,3): B1,B2,B3.
const OC_LCW_SWO_B1 = Float32[0.0,0.0,0.355532,0.0,-0.251389,0.0,-0.251389,0.0,0.118621,0.0,0.0,0.0,0.0,0.364811,0.0,0.3227140,0.0,0.0]
const OC_LCW_SWO_B2 = Float32[0.00371834,0.00308402,0.0,0.00339675,0.00692512,0.0,0.00692512,0.0,0.00384872,0.0,0.0111972,0.0207676,0.0,0.0,0.0,0.0,0.0,0.0]
const OC_LCW_SWO_B3 = Float32[0.808121,0.0,0.0,0.532418,0.985922,0.0,0.985922,0.0,0.0,1.161440,0.0,0.0,1.47018,0.0,1.27196,0.0,1.161440,1.161440]
# organon/crngrow.f HLCW_SWO DACBPAR(18).
const OC_HLCW_SWO_B1 = Float32[0.062000,0.028454,0.05,0.05,0.20,0.209806,0.20,0.209806,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
# organon/crngrow.f CW_SWO CWAPAR(18,3): B1,B2,B3.
const OC_CW_SWO_B1 = Float32[0.929973,0.999291,0.755583,0.755583,0.629785,0.629785,0.629785,0.629785,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]
const OC_CW_SWO_B2 = Float32[-0.135212,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
const OC_CW_SWO_B3 = Float32[-0.0157579,-0.0314603,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]

"organon/crngrow.f LCW_SWO — largest crown width (SCR=0 ⇒ CL=CR·HT branch)."
@inline function oc_lcw_swo(g::Int, mcw::Float32, cr::Float32, scr::Float32, dbh::Float32, ht::Float32)
    b1 = OC_LCW_SWO_B1[g]; b2 = OC_LCW_SWO_B2[g]; b3 = OC_LCW_SWO_B3[g]
    if scr > cr
        cl = scr*ht
        return mcw*fpow(scr, b1 + b2*cl + b3*(dbh/ht))
    else
        cl = cr*ht
        return mcw*fpow(cr, b1 + b2*cl + b3*(dbh/ht))
    end
end

"organon/crngrow.f HLCW_SWO — height to largest crown width."
@inline function oc_hlcw_swo(g::Int, ht::Float32, cr::Float32, scr::Float32)
    b1 = OC_HLCW_SWO_B1[g]
    cl = (scr > cr ? scr : cr)*ht
    return ht - (1.0f0 - b1)*cl
end

"organon/crngrow.f CW_SWO — crown width at height XL (above the largest-crown-width height)."
@inline function oc_cw_swo(g::Int, hlcw::Float32, lcw::Float32, ht::Float32, dbh::Float32, xl::Float32)
    b1 = OC_CW_SWO_B1[g]; b2 = OC_CW_SWO_B2[g]; b3 = OC_CW_SWO_B3[g]
    rp = (ht - xl)/(ht - hlcw)
    ratio = ht/dbh
    if g == 1
        ratio > 50.0f0 && (ratio = 50.0f0)
    elseif g == 2
        ratio > 31.0f0 && (ratio = 31.0f0)
    end
    return lcw*fpow(rp, b1 + b2*sqrt(rp) + b3*ratio)
end

"organon/crngrow.f CALC_CC — accumulate one tree's crown area into the CCH(1..40) profile."
function oc_calc_cc!(g::Int, hlcw::Float32, lcw::Float32, ht::Float32, dbh::Float32,
                     hcb::Float32, expan::Float32, cch::Vector{Float32})
    if hcb > hlcw
        xhlcw = hcb
        xlcw = oc_cw_swo(g, hlcw, lcw, ht, dbh, xhlcw)
    else
        xhlcw = hlcw
        xlcw = lcw
    end
    @inbounds for ii in 40:-1:1
        l = ii - 1
        xl = Float32(l)*(cch[41]/40.0f0)
        if xl <= xhlcw
            cw = xlcw
        elseif xl > xhlcw && xl < ht
            cw = oc_cw_swo(g, hlcw, lcw, ht, dbh, xl)
        else
            cw = 0.0f0
        end
        ca = (cw*cw)*(0.001803f0*expan)
        cch[ii] += ca
    end
    return
end

"""
    oc_crnclo(spgrp, dbh, ht, cr, expan, ntrees) -> cch::Vector{Float32}(41)

organon/crngrow.f CRNCLO (IND=0, CTMUL=0, SCR=0, MGEXP=0) — the crown-closure-by-height profile
`CCH`, from the CURRENT (start-of-cycle) tree list. `CCH[41]` = max tree height; `CCH[1..40]` =
cumulative crown area at 40 relative-height strata. Consumed by HTGRO1 as `TCCH`.
"""
function oc_crnclo(spgrp::Vector{Int32}, dbh::Vector{Float32}, ht::Vector{Float32},
                   cr::Vector{Float32}, expan::Vector{Float32}, ntrees::Int)
    cch = zeros(Float32, 41)
    hmax = ht[1]
    @inbounds for i in 2:ntrees
        ht[i] > hmax && (hmax = ht[i])
    end
    cch[41] = hmax
    @inbounds for i in 1:ntrees
        g = Int(spgrp[i]); d = dbh[i]; h = ht[i]; c = cr[i]
        scr1 = 0.0f0
        cl = c*h
        hcb = h - cl
        ex = expan[i]                              # CTMUL*MGEXP = 0
        mcw = oc_mcw_swo(g, d, h)
        lcw = oc_lcw_swo(g, mcw, c, scr1, d, h)
        hlcw = oc_hlcw_swo(g, h, c, scr1)
        oc_calc_cc!(g, hlcw, lcw, h, d, hcb, ex, cch)
    end
    return cch
end

"""
    oc_hs_hg(isp, si, ht) -> (geage, phtgro)

organon/htgrowth.f HS_HG — Hann & Scrivani (1987) growth-effective age + 5-yr potential height
growth. `isp` 1 = DF-type, 2 = PP-type; `si` = SITE − 4.5 (species-appropriate).
"""
@inline function oc_hs_hg(isp::Int, si::Float32, ht::Float32)
    if isp == 1
        b0 = -6.21693f0; b1 = 0.281176f0; b2 = 1.14354f0
    else
        b0 = -6.54707f0; b1 = 0.288169f0; b2 = 1.21297f0
    end
    bbc = b0 + b1*flog(si)
    x50 = 1.0f0 - fexp((-1.0f0)*fexp(bbc + b2*3.912023f0))
    a1a = 1.0f0 - (ht - 4.5f0)*(x50/si)
    if a1a <= 0.0f0
        return 500.0f0, 0.0f0
    else
        geage = fpow((-1.0f0*flog(a1a))/(fexp(b0)*fpow(si, b1)), 1.0f0/b2)
        xai  = 1.0f0 - fexp(-1.0f0*fexp(bbc + b2*flog(geage)))
        xai5 = 1.0f0 - fexp((-1.0f0)*fexp(bbc + b2*flog(geage + 5.0f0)))
        phtgro = (4.5f0 + (ht - 4.5f0)*(xai5/xai)) - ht
        return geage, phtgro
    end
end

"organon/htgrowth.f HG_SWO — height-increment modifier (species groups 1..5)."
@inline function oc_hg_swo(g::Int, phtgro::Float32, cr::Float32, tcch::Float32)
    p1 = OC_HG_SWO_P1[g]; p2 = OC_HG_SWO_P2[g]; p3 = OC_HG_SWO_P3[g]; p4 = OC_HG_SWO_P4[g]
    p5 = OC_HG_SWO_P5[g]; p6 = OC_HG_SWO_P6[g]; p7 = OC_HG_SWO_P7[g]; p8 = OC_HG_SWO_P8[g]
    fcr = (-p5*fpow(1.0f0 - cr, p6))*fexp(p7*fpow(tcch, 0.5f0))
    b0 = p1*fexp(p2*tcch)
    b1 = fexp(p3*fpow(tcch, p4))
    if fcr < -20.0f0
        modifer = p8*b0
    else
        modifer = p8*(b0 + (b1 - b0)*fexp(fcr))
    end
    cradj = 1.0f0
    if cr <= 0.17f0
        cradj = 1.0f0 - fexp(-fpow(25.0f0*cr, 2.0f0))
    end
    return phtgro*modifer*cradj
end

"""
    oc_limit(isp, dbh, ht, dg, hg) -> hg

organon/htgrowth.f LIMIT (VERSION=1) — cap HG so the tree's height stays under the predicted
height-from-DBH curve, given the 5-yr diameter growth `dg`. `isp` = ORGANON FIA code.
"""
function oc_limit(isp::Int32, dbh::Float32, ht::Float32, dg::Float32, hg::Float32)
    jsp = isp == 263 ? Int32(15) : isp        # VERSION=1: WH→GW-form
    if jsp == 202 || jsp == 263
        a0 = 19.04942539f0; a1 = -0.04484724f0; a2 = 1.0f0
    elseif jsp == 17 || jsp == 15
        a0 = 16.26279948f0; a1 = -0.04484724f0; a2 = 1.0f0
    elseif jsp == 122
        a0 = 17.11482201f0; a1 = -0.04484724f0; a2 = 1.0f0
    elseif jsp == 117
        a0 = 14.29011403f0; a1 = -0.04484724f0; a2 = 1.0f0
    elseif jsp == 351
        a0 = 60.619859f0; a1 = -1.59138564f0; a2 = 0.496705997f0
    else
        a0 = 15.80319194f0; a1 = -0.04484724f0; a2 = 1.0f0
    end
    ht1 = ht - 4.5f0
    ht2 = ht1 + hg
    ht3 = ht2 + hg
    dbh1 = dbh
    dbh2 = dbh1 + dg
    dbh3 = dbh2 + dg
    pht1 = a0*dbh1/(1.0f0 - a1*fpow(dbh1, a2))
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

"organon/htgrowth.f HG_FERT — height-growth fertilizer adjustment (VERSION<=3). 1.0 with no fert."
function oc_hg_fert(cyclg::Int, ispgrp::Int, si_1::Float32, pn::Vector{Float32}, yf::Vector{Float32})
    if ispgrp == 1
        pf1 = 1.0f0; pf2 = 0.333333333f0; pf3 = -1.107409443f0; pf4 = -2.133334346f0; pf5 = 1.5f0
    else
        pf1 = 0.0f0; pf2 = 1.0f0; pf3 = 0.0f0; pf4 = 0.0f0; pf5 = 1.0f0
    end
    faldwn = 1.0f0
    xtime = Float32(cyclg)*5.0f0
    fertx1 = 0.0f0
    @inbounds for i in 2:5
        fertx1 += (pn[i]/800.0f0)*fexp((pf3/pf2)*(yf[1]-yf[i]))
    end
    tem1 = pf3*(xtime - yf[1]) + pf4*fpow(si_1/100.0f0, pf5)
    tem2 = max(-86.0f0, tem1)
    fertx2 = fexp(tem2)
    return 1.0f0 + (pf1*fpow((pn[1]/800.0f0) + fertx1, pf2)*fertx2)*faldwn
end

"organon/htgrowth.f HG_THIN — height-growth thinning adjustment (VERSION<=3). 1.0 with no thin."
function oc_hg_thin(cyclg::Int, ispgrp::Int, babt::Float32, bart::Vector{Float32}, yt::Vector{Float32})
    if ispgrp == 1
        pt1 = -0.3197415492f0; pt2 = 0.7528887377f0; pt3 = -0.2268800162f0
    else
        pt1 = 0.0f0; pt2 = 1.0f0; pt3 = 0.0f0
    end
    gp = 5.0f0
    xtime = Float32(cyclg)*gp
    thinx1 = 0.0f0
    @inbounds for i in 2:5
        thinx1 += bart[i]*fexp((pt3/pt2)*(yt[1]-yt[i]))
    end
    thinx2 = thinx1 + bart[1]
    thinx3 = thinx1 + babt
    prem = thinx3 <= 0.0f0 ? 0.0f0 : thinx2/thinx3
    prem > 0.75f0 && (prem = 0.75f0)
    return 1.0f0 + pt1*fpow(prem, pt2)*fexp(pt3*(xtime - yt[1]))
end

"""
    oc_htgro1(isp, ispgrp, dbh, ht, cr, dgro, cch, si_1, si_2; cyclg, pn, yf, babt, bart, yt) -> HGRO

organon/htgrowth.f HTGRO1 (VERSION=1, big-6 species groups 1..5) — the full height-growth pass for
one tree: TCCH (from `cch`) → HS_HG potential HG → HG_SWO modifier → HG_FERT/HG_THIN → LIMIT (using
the 5-yr diameter growth `dgro`). `isp` = ORGANON FIA code, `ispgrp` = SWO species group.
"""
function oc_htgro1(isp::Int32, ispgrp::Int, dbh::Float32, ht::Float32, cr::Float32, dgro::Float32,
        cch::Vector{Float32}, si_1::Float32, si_2::Float32; cyclg::Int=0,
        pn::Vector{Float32}=zeros(Float32,5), yf::Vector{Float32}=zeros(Float32,5),
        babt::Float32=0.0f0, bart::Vector{Float32}=zeros(Float32,5), yt::Vector{Float32}=zeros(Float32,5))
    # TCCH — crown competition at this tree's height (htgrowth.f:68-77)
    xi = 40.0f0*(ht/cch[41])
    i = trunc(Int, xi) + 2
    xxi = Float32(i) - 1.0f0
    if ht >= cch[41]
        tcch = 0.0f0
    elseif i == 41
        tcch = cch[40]*(40.0f0 - xi)
    else
        tcch = cch[i] + (cch[i-1] - cch[i])*(xxi - xi)
    end
    # potential height growth (Hann-Scrivani), SWO SITE selection (htgrowth.f:87-96)
    if isp == 122
        site = si_2; isisp = 2
    else
        site = si_1
        isp == 81 && (site = (si_1 + 4.5f0)*0.66f0 - 4.5f0)
        isisp = 1
    end
    _, phtgro = oc_hs_hg(isisp, site, ht)
    hg = oc_hg_swo(ispgrp, phtgro, cr, tcch)
    fertadj = oc_hg_fert(cyclg, ispgrp, si_1, pn, yf)
    thinadj = oc_hg_thin(cyclg, ispgrp, babt, bart, yt)
    hg = hg*thinadj*fertadj
    return oc_limit(isp, dbh, ht, dgro, hg)
end

"""
    organon_hg_swo(buf, dgro, spgrp; si_1, si_2, cyclg=0) -> hgro::Vector{Float32}

The ORGANON SWO height-growth pass for one cycle — the `GROW` "growth-2" HG sequence
(organon/grow.f:133-152) driven off the C1 `/ORGANON/` buffer and the C3 `dgro`/`spgrp`. Builds the
`CRNCLO` crown-closure profile from the start-of-cycle tree list, then runs `HTGRO1` for every
big-6 tree (species group ≤ 5), returning `hgro[i]` — the value `oc/htgf.f:96` copies into `HTG`.

`si_1`/`si_2` are SITE_1−4.5 / SITE_2−4.5 (the ORGANON SI the potential-height model consumes).
Minor ORGANON species (group > 5) use HTGRO2 (not yet ported); `hgro[i]` stays 0 for them here.
"""
function organon_hg_swo(buf::OrganonBuffer, dgro::Vector{Float32}, spgrp::Vector{Int32};
        si_1::Float32, si_2::Float32, cyclg::Int=0)
    n = buf.ntrees
    cch = oc_crnclo(spgrp, buf.dbh1, buf.ht1or, buf.cr1, buf.expan1, n)
    hgro = zeros(Float32, n)
    @inbounds for i in 1:n
        buf.expan1[i] <= 0.0f0 && continue
        g = Int(spgrp[i])
        g <= 5 || continue                         # big-6 (HTGRO1); minor species = HTGRO2 (TODO)
        hgro[i] = oc_htgro1(buf.species[i], g, buf.dbh1[i], buf.ht1or[i], buf.cr1[i], dgro[i],
                            cch, si_1, si_2; cyclg=cyclg)
    end
    return hgro
end
