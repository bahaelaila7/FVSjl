# =============================================================================
# height_growth.jl (olympic) — OP FVS-NATIVE large-tree height growth (op/htgf.f +
# op/findag.f + op/htcalc.f) for the NON-ORGANON species. Chunk 1b.
#
# op/htgf.f is the WC-family site-curve height growth: FINDAG inverts the site curve for
# (SITAGE,SITHT), HTCALC gives the potential height at SITAGE+5 (op is a 5-YEAR step —
# AGP10=SITAGE+5.0, op/htgf.f:305 — unlike WC's 10-yr), POTHTG=HGUESS−SITHT, then the
# HGMDCR crown + HGMDRH generalized-Chapman-Richards relative-height modifiers (WTCR=.25).
# op/findag.f is byte-identical to wc/findag.f; op/htcalc.f is the PN edition (DF/WO use
# King, SS/RC use Farr, and the Curtis "misc" set EXCLUDES 16:18 — differs from WC).
#
#   op_htcalc(sindx,ispc,ag)   — site-index height-at-age curve (op/htcalc.f species branches)
#   op_findag(ispc,d,d2,h,si)  — invert the curve → (sitage,sitht,agmax,htmax,htmax2)
#   op_htg_default(...)        — DEFAULT+OWO potential·HGMDCR/HGMDRH, size-cap at HTMAX2
#   height_growth!(s,::Olympic)— per-tree HTG; ORGANON trees (IORG=1) use HGRO (skipped here)
#
# ORGANON trees take HTG=SCALE·XHT·HGRO·exp(HTCON) (op/htgf.f:167) from the NWO engine; only
# the FVS-native IORG=0 trees are grown here. SCALE=FINT/YR with op YR=5 ⇒ SCALE=1 at the 5-yr
# cycle, HTCON≡0, XHT≡1 (MEASURED, opt01 DEBUG HTGF dump).
#
# VALIDATED per-tree vs live FVSop_clean (stand S248112, DEBUG HTGF, cyc0): WF/ES/LP/SP/PP
# reproduce the oracle HTG(I) (op/htgf.f:427 format-901 dump) to its F9.2 print — see
# test/unit/test_op_native_growth.jl.
# =============================================================================

# op/htgf.f DATA — relative-height (RH) modifier coeffs by species shade tolerance (39 species).
const OP_HT_RHR = Float32[
    16.0,16.0,16.0,20.0,16.0, 16.0,15.0,16.0,20.0,16.0, 12.0,13.0,15.0,15.0,13.0,
    15.0,20.0,20.0,20.0,20.0, 20.0,13.0,13.0,15.0,15.0, 12.0,12.0,15.0,13.0,12.0,
    15.0,12.0,20.0,20.0,12.0, 13.0,12.0,15.0,15.0]
const OP_HT_RHYXS = Float32[
    0.15,0.15,0.15,0.20,0.15, 0.15,0.10,0.15,0.20,0.15, 0.01,0.05,0.10,0.15,0.05,
    0.10,0.20,0.20,0.20,0.20, 0.20,0.05,0.05,0.10,0.10, 0.01,0.01,0.10,0.05,0.01,
    0.10,0.01,0.20,0.20,0.01, 0.05,0.01,0.10,0.10]
const OP_HT_RHM = fill(1.10f0, 39)               # op/htgf.f DATA RHM / MAXSP*1.10 /
const OP_HT_RHB = Float32[
    -1.20,-1.20,-1.20,-1.10,-1.20, -1.20,-1.45,-1.20,-1.10,-1.20, -1.60,-1.60,-1.45,-1.45,-1.60,
    -1.45,-1.10,-1.10,-1.10,-1.10, -1.10,-1.60,-1.60,-1.45,-1.45, -1.60,-1.60,-1.45,-1.60,-1.60,
     0.10,-1.60,-1.10,-1.10,-1.60, -1.60,-1.60,-1.45,-1.45]     # sp31 (WB) = +0.10 per DATA
const OP_HT_CRA = 100.0f0; const OP_HT_CRB = 3.0f0; const OP_HT_CRC = -5.0f0
const OP_HT_RHK = 1.0f0;   const OP_HT_RHXS = 0.0f0

# op/findag.f DATA (byte-identical to wc/findag.f): MAPHD (39→8 HT/DBH classes), HDRAT1/2, AGMAX=200.
const OP_HT_MAPHD = Int[
    1,1,1,2,2, 6,2,3,3,2, 4,5,5,5,5, 6,6,3,7,7, 8,8,8,8,8, 8,8,8,3,3, 3,3,3,8,8, 8,8,6,6]
const OP_HT_HDRAT1 = Float32[4.3396271,4.3149844,3.2412923,2.3475244,5.5324838,6.3657425,4.0156013,3.9033821]
const OP_HT_HDRAT2 = Float32[43.9957174,39.6317079,62.7139427,65.7622908,18.6043842,16.2223589,51.9732476,59.3370816]
const OP_HT_AGMAX = 200.0f0

"""
    op_htcalc(sindx, ispc, ag) -> hguess

op/htcalc.f (PN edition) — potential height at age `ag` for species `ispc` on site `sindx`.
Blank/future species (38) and unlisted cases return 0. NOTE vs WC: DF/WO(16,28) use King,
SS/RC(6,18) use Farr, and the Curtis "misc" set is {8,17,21,23:27,29,31:37,39} (excludes 16:18).
"""
@inline function op_htcalc(sindx::Float32, ispc::Int, ag::Float32)
    if ispc == 1                                             # SF — Hoyer PNW-418
        sm45 = sindx - 4.5f0
        k = 0.0071839f0 + 0.0000571f0 * sm45
        return (1f0 - exp(-k*ag))^1.39005f0 / (1f0 - exp(-k*100f0))^1.39005f0 * sm45 + 4.5f0
    elseif ispc == 2 || ispc == 3                            # WF, GF — Cochran PNW-252
        la = log(ag)
        x2 = -0.30935f0 + 1.2383f0*la + 0.001762f0*la^4 - 5.4f-6*la^9 +
             2.046f-7*la^11 - 4.04f-13*la^18
        x3 = -6.2056f0 + 2.097f0*la - 0.09411f0*la^2 - 4.382f-5*la^7 +
             2.007f-11*la^16 - 2.054f-17*la^24
        return exp(x2) - 84.93f0*exp(x3) + (sindx-4.5f0)*exp(x3) + 4.5f0
    elseif ispc == 4 || ispc == 10                           # AF, ES — Alexander RM-32
        return 4.5f0 + (2.75780f0*sindx^0.83312f0) *
               (1f0 - exp(-0.015701f0*ag))^(22.71944f0*sindx^(-0.63557f0))
    elseif ispc == 5                                         # RF — Dolph PSW-206
        term  = ag*exp(ag*(-0.0440853f0))*1.4151f-6
        b     = sindx*term - 3.0495f6*term*term + 5.72474f-4
        term2 = 50f0*exp(50f0*(-0.0440853f0))*1.4151f-6
        b50   = sindx*term2 - 3.0495f6*term2*term2 + 5.72474f-4
        return (sindx-4.5f0)*(1f0-exp(-b*(ag^1.51744f0))) / (1f0-exp(-b50*(50f0^1.51744f0))) + 4.5f0
    elseif ispc == 7                                         # NF — Herman PNW-243
        s45 = sindx - 4.5f0
        x1 = -564.38f0 + 22.25f0*s45 - 0.04995f0*s45^2
        x2 = 6.80f0 + 2843.21f0*s45^(-1) + 34735.54f0*s45^(-2)
        return 4.5f0 + s45 / (x1*(1f0/ag)^2 + x2*(1f0/ag) + 1f0 - 0.0001f0*x1 - 0.01f0*x2)
    elseif ispc == 9 || ispc == 12 || ispc == 15             # IC, JP, PP — Barrett PNW-232
        t = (-0.7864f0 + 2.49717f0*(1f0-exp(-0.0045042f0*ag))^0.33022f0)
        return 128.8952205f0*(1f0-exp(-0.016959f0*ag))^1.23114f0 - t*100.43f0 +
               t*(sindx-4.5f0) + 4.5f0
    elseif ispc == 11                                        # LP — Dahms PNW-8
        return sindx*(-0.0968f0 + 0.02679f0*ag - 0.00009309f0*ag*ag)
    elseif ispc == 13 || ispc == 14                          # SP, WP — Curtis PNW-423
        num = 1f0 - exp(-exp(-4.62536f0 + 1.346399f0*log(ag) - 135.354483f0/sindx))
        den = 1f0 - exp(-exp(-4.62536f0 + 1.346399f0*log(100f0) - 135.354483f0/sindx))
        return num/den * (sindx-4.5f0) + 4.5f0
    elseif ispc == 19                                        # WH — Wiley 1978
        z = 2500f0/(sindx-4.5f0)
        return ag*ag/(-1.7307f0 + 0.1394f0*z + (-0.0616f0+0.0137f0*z)*ag +
               (0.00192f0+0.00007f0*z)*(ag*ag)) + 4.5f0
    elseif ispc == 20                                        # MH — Means (unpublished)
        h = (22.8741f0 + 0.950234f0*sindx)*(1f0-exp(-0.00206465f0*sqrt(sindx)*ag))^
            (1.365566f0 + 2.045963f0/sindx)
        return (h + 1.37f0)*3.281f0
    elseif ispc == 22                                        # RA — Harrington PNW-358
        c = (59.5864f0 + 0.7953f0*sindx); k = (0.00194f0 - 0.0007403f0*sindx)
        return sindx + c*(1f0-exp(k*ag))^0.9198f0 - c*(1f0-exp(k*20f0))^0.9198f0
    elseif ispc == 30                                        # LL — Cochran PNW-424
        t = (-0.12528f0 + 0.039636f0*ag - 4.278f-4*ag*ag + 1.7039f-6*ag^3)
        return 4.5f0 + 1.46897f0*ag + 0.0092466f0*ag*ag - 2.3957f-4*ag^3 + 1.1122f-6*ag^4 +
               (sindx-4.5f0)*t - 73.57f0*t
    elseif ispc == 16 || ispc == 28                          # DF, WO — King Weyerhaeuser #8
        z = 2500f0/(sindx-4.5f0)
        return ag*ag/(-0.954038f0 + 0.109757f0*z + (5.58178f-2+7.92236f-3*z)*ag +
               (-7.33819f-4+1.97693f-4*z)*(ag*ag)) + 4.5f0
    elseif ispc == 6 || ispc == 18                           # SS, RC — Farr PNW-326
        la = log(ag)
        az = -0.2050542f0 + 1.449615f0*la - 0.01780992f0*la^3 +
              6.519748f-5*la^5 - 1.095593f-23*la^30
        bz = -5.611879f0 + 2.418604f0*la - 0.2593110f0*la^2 +
              1.351445f-4*la^5 - 1.701139f-12*la^16 + 7.964197f-27*la^36
        return 4.5f0 + exp(az) - 86.43f0*exp(bz) + (sindx-4.5f0)*exp(bz)
    elseif ispc == 8 || ispc == 17 || ispc == 21 || (23 <= ispc <= 27) ||
           ispc == 29 || (31 <= ispc <= 37) || ispc == 39   # "misc" — Curtis (EXCLUDES 16:18)
        s45 = sindx - 4.5f0
        return s45 / (0.6192f0 - 5.3394f0/s45 + 240.29f0*ag^(-1.4f0) +
               (3368.9f0/s45)*ag^(-1.4f0)) + 4.5f0
    else                                                     # blank/future (38)
        return 0.0f0
    end
end

"""
    op_findag(ispc, d, d2, h, sindx) -> (sitage, sitht, agmax, htmax, htmax2)

op/findag.f — step AG by 2 up the site curve to solve (SITAGE,SITHT) from tree height `h`
(TOLER=2, INCRNG flatten-detection). HTMAX/HTMAX2 = HDRAT1(MAPHD)·D(+D2)+HDRAT2. If H≥HTMAX,
age is dubbed at 0.10 ft/yr above AGMAX.
"""
@inline function op_findag(ispc::Int, d::Float32, d2::Float32, h::Float32, sindx::Float32)
    agmax = OP_HT_AGMAX
    mp = OP_HT_MAPHD[ispc]
    htmax  = OP_HT_HDRAT1[mp]*d  + OP_HT_HDRAT2[mp]
    htmax2 = OP_HT_HDRAT1[mp]*d2 + OP_HT_HDRAT2[mp]
    if h >= htmax
        return (agmax + (h - htmax)/0.10f0, h, agmax, htmax, htmax2)
    end
    ag = 2.0f0; hguess = 0.0f0
    while true
        oldhg = hguess
        incrng = false
        hguess = op_htcalc(sindx, ispc, ag)
        if hguess >= 1.0f0
            diff = abs(hguess - h)
            (diff <= 2.0f0 || h < hguess) && return (ag, hguess, agmax, htmax, htmax2)
            d2h = hguess - oldhg
            (oldhg != 0.0f0 && d2h >= 0.05f0) && (incrng = true)
            (incrng && d2h < 0.05f0) && return (ag, hguess, agmax, htmax, htmax2)
        end
        ag += 2.0f0
        ag > agmax && return (agmax, h, agmax, htmax, htmax2)
    end
end

"""
    op_htg_default(ispc, sindx, d, h, icr_pct, avh, ba, dg, d2, sitage, sitht, agmax, htmax2) -> htg

op/htgf.f DEFAULT (+OWO 28) height increment for one tree, BEFORE SCALE·XHT·exp(HTCON)·MISHGF.
op is a 5-YEAR step (AGP10=SITAGE+5.0, asymptotic POTHTG=(0.2+0.00264·SI)·5 for WF/GF). Assumes
PROB>0 and H≤HTMAX (the H>HTMAX and RW branches are handled in height_growth!).
"""
@inline function op_htg_default(ispc::Int, sindx::Float32, d::Float32, h::Float32,
                                icr_pct::Float32, avh::Float32, ba::Float32, dg::Float32,
                                d2::Float32, sitage::Float32, sitht::Float32,
                                agmax::Float32, htmax2::Float32)
    local pothtg::Float32
    if sitage >= agmax
        pothtg = (ispc == 2 || ispc == 3) ? (0.2f0 + 0.00264f0*sindx)*5f0 : 0.10f0
    else
        hguess = op_htcalc(sindx, ispc, sitage + 5.0f0)
        pothtg = hguess - sitht
        if ispc == 28                                        # OWO King HT-DBH patch (op/htgf.f:322)
            maxg = sindx - 18.6024f0/log(2.7f0 + ba)
            dd2 = d + dg; dd2 < 0f0 && (dd2 = 0.1f0)
            hg2 = 4.5f0 + maxg*(1f0-exp(-0.137428f0*dd2))^1.38994f0
            hg1 = 4.5f0 + maxg*(1f0-exp(-0.137428f0*d))^1.38994f0
            pothtg = (hg2 - hg1) / 2.0f0                      # DG is 10-yr here ⇒ /2 to 5-yr
        end
        pothtg < 0.1f0 && (pothtg = 0.1f0)
    end
    cr = icr_pct/100f0
    hgmdcr = OP_HT_CRA * cr^OP_HT_CRB * exp(OP_HT_CRC*cr)
    hgmdcr > 1f0 && (hgmdcr = 1f0)
    relht = avh > 0f0 ? h/avh : 0f0
    relht > 1.5f0 && (relht = 1.5f0)
    rhx = relht
    fctrkx = (OP_HT_RHK/OP_HT_RHYXS[ispc])^(OP_HT_RHM[ispc]-1f0) - 1f0
    fctrrb = -1f0*(OP_HT_RHR[ispc]/(1f0-OP_HT_RHB[ispc]))
    fctrxb = rhx^(1f0-OP_HT_RHB[ispc]) - OP_HT_RHXS^(1f0-OP_HT_RHB[ispc])
    fctrm  = -1f0/(OP_HT_RHM[ispc]-1f0)
    hgmdrh = OP_HT_RHK * (1f0 + fctrkx*exp(fctrrb*fctrxb))^fctrm
    htgmod = 0.25f0*hgmdcr + 0.75f0*hgmdrh
    htgmod >= 2f0 && (htgmod = 2f0)
    htgmod <= 0f0 && (htgmod = 0.1f0)
    htg = pothtg*htgmod
    (h + htg > htmax2) && (htg = htmax2 - h)
    htg < 0.1f0 && (htg = 0.1f0)
    return htg
end

"""
    height_growth!(s, ::Olympic; scale=1.0f0)

op/htgf.f per-tree periodic height increment into `trees.ht_growth`. ORGANON trees (IORG=1) get
HTG=SCALE·XHT·HGRO·exp(HTCON) from the NWO engine (not implemented here — skipped, left 0). RW(17)
uses the LTHTG special form. Requires the tree's `diam_growth` (DG) already computed.
"""
function height_growth!(s::StandState, ::Olympic; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    avh = p.avg_height; ba = p.basal_area
    org_ran = length(c.op_iorg) == t.n           # ORGANON ran this cycle ⇒ IORG=1 trees take stashed HGRO
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        ispc = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]; dg = t.diam_growth[i]
        # ORGANON-grown trees (IORG=1) take HTG=SCALE·XHT·HGRO·exp(HTCON) (op/htgf.f:167).
        if org_ran && c.op_iorg[i] == 1
            t.ht_growth[i] = scale * c.htg_mult[ispc] * c.op_hgro[i] * exp(c.htg_cor[ispc])
            continue
        end
        sindx = ispc <= length(p.sp_site_index) ? p.sp_site_index[ispc] : 0f0
        icr_pct = Float32(t.crown_pct[i])
        htcon = c.htg_cor[ispc]
        brat = op_bratio(ispc, d)
        d2 = d + dg/brat
        local htg::Float32
        if ispc == 17                                        # REDWOOD — LTHTG special (op/htgf.f CASE 17)
            dg10 = h < 4.5f0 ? 0.1f0 : (dg*2f0)/brat
            lthtg = exp(1.412947f0 - 0.000204f0*d*d + 0.31971f0*log(d) +
                        0.394005f0*log(sindx) + 0.399888f0*log(dg10) - 0.451708f0*log(h)) * 0.5f0
            hgbnd = if h >= 217.0f0 && h < 380.0f0
                        max(1.0f0 - (h-217.0f0)/(380.0f0-217.0f0), 0.1f0)
                    elseif h < 217.0f0; 1.0f0 else 0.1f0 end
            htg = lthtg*hgbnd
            htg < 0.1f0 && (htg = 0.1f0)
            htg = scale*htg*exp(htcon)
        else                                                 # DEFAULT (+OWO 28)
            sitage, sitht, agmax, htmax, htmax2 = op_findag(ispc, d, d2, h, sindx)
            if h > htmax
                if h >= htmax2
                    htg = 0.5f0*dg; htg < 0.1f0 && (htg = 0.1f0)
                    htg = scale*htg*exp(htcon)
                else
                    htg = 0.0f0
                end
            else
                htg = op_htg_default(ispc, sindx, d, h, icr_pct, avh, ba, dg, d2,
                                     sitage, sitht, agmax, htmax2)
                htg = scale*htg*exp(htcon)
            end
        end
        cap = s.control.sp_size_cap[ispc, 4]
        if cap > 0f0 && h + htg > cap
            htg = cap - h
            htg < 0.1f0 && (htg = 0.1f0)
        end
        t.ht_growth[i] = htg
    end
    return s
end
