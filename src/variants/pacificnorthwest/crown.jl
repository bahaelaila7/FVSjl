# =============================================================================
# crown.jl (pacificnorthwest) — PN crown ratio (pn/crown.f + pn/dubscr.f) + CCF (pn/ccfcal.f). Chunks 5 + density.
#
# Rank-based Weibull crown ratio, IDENTICAL formula to WC (crown_ratio_update! mirrors the WC one; the
# WC version is validated 81/81 bit-exact), only the coefficients differ: pn/crown.f has **17 crown
# groups** (WC 16; SS added) and pn/ccfcal.f has **19 CCF groups** (WC 16). All PN_* consts extracted
# bit-exact from pn/*.f DATA. DUBSCR (sub-1"/dead crown dub) is byte-identical to WC (6 BCR groups + RW).
# RW(sp 17) logistic + DUBSCR ported source-faithful (unexercised on pnt01: no RW, no missing crowns).
# =============================================================================

# pn/crown.f DATA WEIBA/WEIBB0/WEIBB1/WEIBC0/WEIBC1/C0/C1 — 17 crown groups (index = crown_imap[ISPC]).
const PN_WEIBA  = Float32[0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,1.0,0.0,0.0]
const PN_WEIBB0 = Float32[-0.17168,0.130939,-0.981113,-0.135807,0.019948,-0.036696,-0.012061,-0.062693,0.073435,0.162672,0.196054,-0.818809,0.035786,-0.238295,-0.811424,-0.13121,-0.107413]
const PN_WEIBB1 = Float32[1.161549,1.093406,1.092273,1.147712,1.108738,1.132792,1.119712,1.139657,1.107183,1.073404,1.073909,1.054176,1.121389,1.180163,1.05619,1.15976,1.140775]
const PN_WEIBC0 = Float32[2.8263,1.355139,1.326047,3.017494,2.62123,2.876094,3.2126,1.7664,2.6237,3.288501,0.345647,-2.366108,2.0408,3.044134,-3.831124,2.598238,3.0712]
const PN_WEIBC1 = Float32[0.0,0.350472,0.318386,0.0,0.186734,0.0,0.0,0.0,0.0,0.0,0.620145,1.202413,0.0,0.0,1.401938,0.0,0.0]
const PN_CRC0   = Float32[5.073342,5.212394,4.860467,5.568864,4.279655,5.073273,5.666442,4.48133,5.671345,6.484942,5.417431,4.42,4.656659,4.625125,5.20055,4.890318,5.812879]
const PN_CRC1   = Float32[-0.01143,-0.011623,-0.006173,-0.021293,-0.002484,-0.020988,-0.025199,-0.018092,-0.023463,-0.023248,-0.011608,-0.01066,-0.022612,-0.016042,-0.01489,-0.018837,-0.028504]

# pn/dubscr.f — byte-identical to WC (6 BCR groups). IMAP[39]→6 groups.
const PN_DUB_IMAP = Int32[1,1,1,1,1,1,1,3,3,1,4,4,4,4,4,2,2,3,3,3,5,5,5,5,5,5,5,5,6,4,4,4,4,5,5,5,5,5,5]
const PN_DUB_BCR0 = Float32[8.042774,8.477025,7.558538,6.489813,5.0,9.0]
const PN_DUB_BCR1 = Float32[0.007198,-0.018033,-0.015637,-0.029815,0.0,0.0]
const PN_DUB_BCR2 = Float32[-0.016163,-0.01814,-0.009064,-0.009276,0.0,0.0]
const PN_DUB_CRSD = Float32[1.3167,1.3756,1.9658,2.0426,0.5,0.5]

@inline function _pn_dubscr(rng, sp::Integer, d::Real, h::Real, ba::Real, prd::Real, qmdplt::Real)::Float32
    g = Int(PN_DUB_IMAP[sp])
    if sp == 17
        hdr = Float32(h) * 12f0 / Float32(d)
        cr = -1.021064f0 + 0.309296f0 * flog(max(hdr, 1f-6)) + 0.869720f0 * Float32(prd) -
             0.116274f0 * (Float32(d) / Float32(qmdplt))
        sd = 0.15f0; fcr = 0f0
        while true; fcr = bachlo(rng, 0f0, sd); abs(fcr) > sd && continue; break; end
        cr = 1f0 / (1f0 + fexp(cr + fcr))
    else
        cr = PN_DUB_BCR0[g] + PN_DUB_BCR1[g] * Float32(h) + PN_DUB_BCR2[g] * Float32(ba)
        sd = PN_DUB_CRSD[g]; fcr = 0f0
        while true; fcr = bachlo(rng, 0f0, sd); abs(fcr) > sd && continue; break; end
        cr = (((cr + fcr) - 1f0) * 10f0 + 1f0) / 100f0     # pn/dubscr.f: CR=((CR-1.0)*10.0+1.0)/100.
    end
    cr > 0.95f0 && (cr = 0.95f0); cr < 0.05f0 && (cr = 0.05f0)
    return cr
end

function crown_ratio_update!(s::StandState, ::PacificNorthwest; fint::Float32 = 10.0f0, lstart::Bool = false,
                             crown_sdi::Float32 = 0f0, kwargs...)
    p, t = s.plot, s.trees
    n = t.n; (n == 0 && !lstart) && return s
    sd = s.coef.species; cimap = sd[:crown_imap]
    relden = p.relative_density; sdiac = crown_sdi; ba = p.basal_area
    # crown.f ISORT: whole-stand descending-DBH rank on the CURRENT DBH (grown at cycling, as-read at
    # LSTART) — shared crown_isort, see crown_init.jl.
    isort = crown_isort(s; lstart = lstart)
    # pn/crown.f:87-104 (2021 edits): XMAXPT via SDICAL + per-point Zeide ZRD via SDICLS, live trees only,
    # computed once on entry. PRD = ZRD/XMAXPT and QMDPLT = sqrt((PTBAA/PTPA)/0.005454) (floored at 1) are
    # PER INVENTORY POINT for every tree — they feed only the RW (sp 17) logistic and DUBSCR's RW branch.
    xmaxpt, zrd, _ = point_zeide!(s)
    dens = s.density
    _prd(pt) = (1 <= pt <= length(xmaxpt) && xmaxpt[pt] > 0f0) ? zrd[pt] / xmaxpt[pt] : 0f0
    function _qmdplt(pt)
        baplt = (1 <= pt <= length(dens.point_ba)) ? dens.point_ba[pt] : 0f0
        tpaplt = (1 <= pt <= length(dens.point_tpa)) ? dens.point_tpa[pt] : 0f0
        q = tpaplt > 0f0 ? sqrt((baplt / tpaplt) / 0.005454f0) : 1f0
        q <= 1f0 ? 1f0 : q
    end
    # pn/crown.f:121 `IF((ITRN.LE.0).AND.(IREC2.LT.MAXTP1)) GO TO 74`: an all-dead inventory still dubs the dead crowns.
    t.n == 0 && @goto dead79
    # pn/crown.f:150-178 — the CRNMULT keyword (a scheduled activity) overwrites CRNMLT/DLOW/DHI per species.
    cur_year = current_cycle_year(s)
    # crown.f DO 70 ISPC … I=IND1(I3): SPECIES-MAJOR — the DUBSCR/RANN draws follow this order, not storage.
    @inbounds for i in species_major_order(s)
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        (sp < 1 || sp > 39) && continue
        (lstart && t.crown_pct[i] > 0) && continue
        icr = Int(t.crown_pct[i])
        cmult, cdlow, cdhi = crn_mult_band(s.control, sp, cur_year; lstart = lstart)
        inband = cdlow <= d <= cdhi                      # pn/crown.f:344/351/400 `.GE. DLOW .AND. .LE. DHI`
        prd = _prd(Int(t.plot_id[i])); qmdplt = _qmdplt(Int(t.plot_id[i]))
        if d < 1f0 && lstart
            icr != 0 && continue
            cr = _pn_dubscr(s.rng, sp, d, h, ba, prd, qmdplt)
            icri = trunc(Int, cr * 100f0 + 0.5f0)
            inband && (icri = trunc(Int, Float32(icri) * cmult))          # pn/crown.f:400-401
            icri > 95 && (icri = 95); (cmult == 1f0 && icri < 10) && (icri = 10); icri < 1 && (icri = 1)
            t.crown_pct[i] = Int32(icri); continue
        end
        grp = Int(cimap[sp])
        relsdi = p.sp_sdi_def[sp] > 0f0 ? sdiac / p.sp_sdi_def[sp] : 1f0
        relsdi > 1.5f0 && (relsdi = 1.5f0)
        acrnew = PN_CRC0[grp] + PN_CRC1[grp] * relsdi * 100f0
        local crnew::Float32
        if sp == 17
            hdr = d > 0f0 ? h * 12f0 / d : 1f0
            xl = -1.021064f0 + 0.309296f0 * flog(max(hdr, 1f-6)) + 0.869720f0 * prd - 0.116274f0 * (d / qmdplt)
            x = 1f0 / (1f0 + fexp(xl)); x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = x * 10f0
        else
            A = PN_WEIBA[grp]
            B = PN_WEIBB0[grp] + PN_WEIBB1[grp] * acrnew; B < 3f0 && (B = 3f0)
            C = PN_WEIBC0[grp] + PN_WEIBC1[grp] * acrnew; C < 2f0 && (C = 2f0)
            scale = 1f0 - 0.00167f0 * (relden - 100f0)
            scale > 1f0 && (scale = 1f0); scale < 0.30f0 && (scale = 0.30f0)
            x = d > 0f0 ? (Float32(isort[i]) / Float32(n)) * scale : rann!(s.rng) * scale
            x < 0.05f0 && (x = 0.05f0); x > 0.95f0 && (x = 0.95f0)
            crnew = A + B * fpow((-flog(1f0 - x)), 1f0 / C)
        end
        crnew *= 10f0
        if !(lstart || icr == 0)
            chg = crnew - Float32(icr); pdifpy = chg / Float32(icr) / fint
            pdifpy > 0.01f0 && (chg = Float32(icr) * 0.01f0 * fint)
            pdifpy < -0.01f0 && (chg = Float32(icr) * (-0.01f0) * fint)
            crnew = Float32(icr) + (inband ? chg * cmult : chg)   # stmt 41, pn/crown.f:344-348
        end
        icri = trunc(Int, crnew + 0.5f0)
        ((lstart || icr == 0) && inband) && (icri = trunc(Int, Float32(icri) * cmult))   # 9052, :350-353
        if !(lstart || icr == 0)
            htg = t.ht_growth[i]; crln = h * Float32(icr) / 100f0
            crmax = (crln + htg) / (h + htg) * 100f0
            (cmult == 1f0 && icri < 10) && (icri = trunc(Int, crmax + 0.5f0))   # pn/crown.f:375
            Float32(icri) > crmax && (icri = trunc(Int, crmax + 0.5f0))
        end
        lstart && (icri = topkill_icri(t, i, icri))   # crown.f stmt 55 (crown_init.jl)
        icri > 95 && (icri = 95); (cmult == 1f0 && icri < 10) && (icri = 10); icri < 1 && (icri = 1)
        t.crown_pct[i] = Int32(icri)
    end
    # ---- CYCLE-0 DEAD-RECORD CROWN DUB (pn/crown.f:422-475 `DO 79 I=IREC2,MAXTRE`) ----
    # FVS dubs MISSING crowns on the inventory standing-dead records with DUBSCR (a rejection-bounded
    # main-stream BACHLO per record), iterating IREC2→MAXTRE = the REVERSE of jl's dead storage (same layout
    # as bm/ie crown.f DO 79). QMDPLT/PRD are per inventory point (PTBAA/PTPA; ZRD/XMAXPT). Skipping this left
    # the dead crowns 0 AND consumed none of these draws, so every later DGSCOR draw was shifted.
    @label dead79
    lstart && dub_dead_crowns!(s) do i
        pt = Int(t.plot_id[i])
        icri_round(_pn_dubscr(s.rng, Int(t.species[i]), t.dbh[i], t.height[i], ba, _prd(pt), _qmdplt(pt)))
    end
    return s
end

# pn/ccfcal.f MODE=1 — per-tree CCF (19 groups). stand_ccf sums CCFT·P = RELDEN; point_density → PCCF.
const PN_CCF_INDCCF = Int[1,2,2,3,4,8,3,5,5,7,10,19,19,11,19,12,12,13,14,14,15,16,16,17,15,17,18,17,6,19,9,9,6,15,17,15,15,10,10]
const PN_CCF_RD1 = Float32[0.10142,0.0690403,0.0245276,0.0172,0.0194415,0.0318054,0.0288484,0.0761779,0.01925,0.0220871,0.0185728,0.0387616,0.0288484,0.037577,0.0160051,0.115394,0.0170887,0.000450757,0.0219]
const PN_CCF_RD2 = Float32[0.0432725,0.0224682,0.0114741,0.00876,0.0142461,0.0215065,0.0173091,0.0421908,0.01676,0.0252424,0.014621,0.0268821,0.0237999,0.0232893,0.0166659,0.0441381,0.0213617,0.0029209,0.0168]
const PN_CCF_RD3 = Float32[0.00461575,0.00182799,0.0013419,0.00112,0.00260979,0.00363562,0.00259636,0.0058418,0.00365,0.0072121,0.0028775,0.00466086,0.00490874,0.00360853,0.00433848,0.0042207,0.00667579,0.00473186,0.00325]

# ---------------------------------------------------------------------------
# FFE crown-biomass group per species (pn/fmcrow.f:103 DATA ISPMAP). pn/fmcrow.f routes
# CASE(24,26,27,34,35,36,37,39)=PB/AS/CW/DG/HT/CH/WI/OT → FMCROWE (eastern Jenkins TOTABV), all others
# → FMCROWW (western crown-width, shared cr_crownw). Identical to WC's ISPMAP except slot 6 = 18 (SS,
# uses ES; WC blank=0). FMCROWE routing CASE is byte-identical to WC. PN FFE chunk F1.
# ---------------------------------------------------------------------------
const PN_ISPMAP = Int[
   4,  4,  4,  1,  4, 18,  4,  8, 20, 18,
  11, 15, 15, 15, 13,  3, 19,  7,  6, 24,
   5, 23, 10, 43, 17, 41, 17, 17, 16,  1,
  14, 11,  7, 56, 57, 61, 64,  0, 41]
@inline pn_uses_fmcrowe(sp::Integer) = (sp == 24 || sp == 26 || sp == 27 || sp == 34 ||
                                        sp == 35 || sp == 36 || sp == 37 || sp == 39)

@inline function pn_tree_ccf(sp::Integer, d::Real)::Float32
    (sp < 1 || sp > 39) && return 0f0
    ic = PN_CCF_INDCCF[sp]; D = Float32(d)
    return D < 1.0f0 ? D * (PN_CCF_RD1[ic] + PN_CCF_RD2[ic] + PN_CCF_RD3[ic]) :
                       PN_CCF_RD1[ic] + PN_CCF_RD2[ic] * D + PN_CCF_RD3[ic] * fpow(D, 2f0)   # pn/ccfcal.f:109 RD3*D**2.0 (powf)
end
