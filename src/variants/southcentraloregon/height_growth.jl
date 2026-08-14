# =============================================================================
# height_growth.jl (southcentraloregon) — SO large-tree height growth (so/htgf.f + so/findag.f). Chunk 4b.
#
# so/findag.f inverts each species' OWN site curve (so_htcalc, chunk 2) to a stand age SITAGE at the tree's
# current height. so/htgf.f then grows potential height 10 yr and applies one of THREE sub-methods:
#   • DEFAULT + PP(10): potential HGUESS = so_htcalc(SITAGE+10) − SITHT, then the Hoerl CR modifier
#     (CRA/CRB/CRC) + generalized-Chapman-Richards RH modifier (RHR/RHYXS/RHM/RHB), weighted 0.25/0.75.
#   • SH/WO(9,27) — CA method: potential + Ritchie&Hann small-tree modifier (CRMOD/RHMOD, 1.016605).
#   • WB/AS(16,24) — Johnson SBB (Schreuder-Hafley) via COF1(TT)/COF6(UT) tables; NOT modifier-adjusted.
#   • JU(11): height comes from REGENT (skip here).
# All capped at the species/site HTMAX, scaled FINT/YR·XHT·exp(HTCON). MEASURED vs FVSso_g16 htgf per-tree HTG.
# =============================================================================

# so/findag.f DATA AGMAX/HTMAX (per-species asymptotic age / max height), 33 species.
const SO_FINDAG_AGMAX = Float32[500,400,900,450,650,650,350,400,500,900,900,350,350,400,400,400,400,550,550,350,350,100,100,100,100,50,250,50,75,50,50,400,100]
const SO_FINDAG_HTMAX = Float32[165,160,180,180,150,150,130,165,180,175,80,165,120,165,165,85,175,165,165,50,50,100,100,75,125,30,75,30,30,20,25,165,100]

# so/htgf.f DATA — Hoerl CR modifier + generalized-Chapman-Richards RH modifier (shade-tolerance based).
const SO_HTGF_CRA = 100.0f0; const SO_HTGF_CRB = 3.0f0; const SO_HTGF_CRC = -5.0f0
const SO_HTGF_RHK = 1.0f0; const SO_HTGF_RHXS = 0.0f0
const SO_HTGF_RHR   = Float32[15,15,15,20,20,20,12,16,16,13,13,13,20,20,15,15,12,20,20,20,13,13,20,15,12,13,15,12,15,15,15,15,15]
const SO_HTGF_RHYXS = Float32[0.10,0.10,0.10,0.20,0.20,0.20,0.01,0.15,0.15,0.05,0.10,0.20,0.05,0.10,0.10,0.10,0.10,0.10,0.20,0.20,0.05,0.05,0.20,0.10,0.01,0.05,0.10,0.01,0.10,0.10,0.10,0.10,0.10]
const SO_HTGF_RHM   = fill(1.10f0, 33)
const SO_HTGF_RHB   = Float32[-1.45,-1.45,-1.45,-1.10,-1.10,-1.10,-1.60,-1.20,-1.20,-1.60,-1.45,-1.10,-1.20,-1.10,-1.45,-1.60,-1.60,-1.10,-1.10,-1.10,-1.60,-1.60,-1.10,-1.45,-1.60,-1.60,-1.45,-1.60,-1.45,-1.45,-1.45,-1.45,-1.45]

# so/htgf.f DATA COF1 (WB, from TT/HTGF) / COF6 (AS, from UT/HTGF) — Johnson-SBB coefficients.
# Column-major DATA COF(9,3): SO_HTGF_COF1[k, keycr] = COF(k, keycr), k=1..9 coef, keycr=1..3 crown group.
const SO_HTGF_COF1 = Float32[
    37.0     45.0     45.0;
    85.0     100.0    90.0;
     1.77836  1.66674  1.64770;
    -0.51147  0.25626  0.30546;
     1.88795  1.45477  1.35015;
     1.20654  1.11251  0.94823;
     0.57697  0.67375  0.70453;
     3.57635  2.17942  2.46480;
     0.90283  0.88103  1.00316]
const SO_HTGF_COF6 = Float32[
    30.0     30.0     35.0;
    85.0     85.0     85.0;
     2.00995  2.00995  1.80388;
     0.03288  0.03288 -0.07682;
     1.81059  1.81059  1.70032;
     1.28612  1.28612  1.29148;
     0.72051  0.72051  0.72343;
     3.00551  3.00551  2.91519;
     1.01433  1.01433  0.95244]

"""
    so_findag(ispc, ifor, h, sindx) -> (sitage, sitht, agmax1, htmax1)

so/findag.f: invert species `ispc`'s own site curve (so_htcalc) to the stand age matching height `h`.
"""
function so_findag(ispc::Int, ifor::Int, h::Float32, sindx::Float32)
    toler = 2.0f0
    agmax1 = SO_FINDAG_AGMAX[ispc]
    (ifor > 3 && ifor < 10) && (agmax1 = 400f0)          # R5 forests
    htmax1 = SO_FINDAG_HTMAX[ispc]
    ag = 2.0f0
    (ispc == 2 || ispc == 10) && (ag = 98.38f0 * exp(sindx * (-0.0422f0)) + 1.0f0)
    ag < 2.0f0 && (ag = 2.0f0)
    ispc == 3 && (ag = 18.0f0)
    if h >= htmax1                                        # H exceeds site max ⇒ 0.10 ft/yr extrapolation
        return (agmax1 + (h - htmax1) / 0.10f0, h, agmax1, htmax1)
    end
    if ispc == 24                                         # AS: Sheppard age eqn (no iteration)
        return ((h * 2.54f0 * 12.0f0 / 26.9825f0)^(1.0f0 / 1.1752f0), h, agmax1, htmax1)
    elseif ispc == 11 || ispc == 16                       # JU/WB: age irrelevant
        return (0f0, h, agmax1, htmax1)
    end
    incrng = false; hguess = 0f0
    while true
        oldhg = hguess
        hguess = so_htcalc(ifor, sindx, ispc, ag)
        if hguess >= 1.0f0
            diff = abs(hguess - h)
            (diff <= toler || h < hguess) && return (ag, hguess, agmax1, htmax1)
            d2 = hguess - oldhg
            (oldhg != 0f0 && d2 >= 0.05f0) && (incrng = true)
            (incrng && d2 < 0.05f0) && return (ag, hguess, agmax1, htmax1)
        end
        ag += 2.0f0
        ag > agmax1 && return (agmax1, h, agmax1, htmax1)
    end
end

# so/htgf.f Johnson-SBB (Schreuder-Hafley) height growth for WB(16)/AS(24). Returns HTG (ft), 0.1 floor.
@inline function _so_sbb_htg(ispc::Int, d::Float32, h::Float32, dg::Float32, bark::Float32,
                             icr::Int, elev::Float32, iage::Int, icyc::Int, pct::Float32,
                             relden::Float32, fint::Float32, yr::Float32)::Float32
    (d <= 0.1f0 || h <= 4.5f0) && return 0.1f0            # from REGENT; SBB undefined
    cof = ispc == 16 ? SO_HTGF_COF1 : SO_HTGF_COF6
    xi2 = 4.5f0; xi1 = 0.1f0
    iicr = trunc(Int, Float32(icr) / 10.0f0 + 0.5f0)      # crown → 0..9
    iicr > 9 && (iicr = 9)
    keycr = iicr <= 1 ? 1 : iicr <= 7 ? 2 : 3             # GO TO(101,101,102,102,102,102,102,103,103)
    # bounds: SBB undefined if DBH/HT exceed the fitted range ⇒ HTG=0.1
    (h <= 4.5f0) && return 0.1f0
    ((xi1 + cof[1, keycr]) <= d) && return 0.1f0
    ((xi2 + cof[2, keycr]) <= h) && return 0.1f0
    y1 = (d - xi1) / cof[1, keycr]
    y2 = (h - xi2) / cof[2, keycr]
    fby1 = log(y1 / (1.0f0 - y1))
    fby2 = log(y2 / (1.0f0 - y2))
    z = (cof[4, keycr] + cof[6, keycr] * fby2 - cof[7, keycr] * (cof[3, keycr] + cof[5, keycr] * fby1)) *
        (1.0f0 - cof[7, keycr]^2)^(-0.5f0)
    # ZBIAS is 0 for WB/AS (AZBIAS/BZBIAS all 0), except the AS elev-window ZADJ below.
    if ispc == 24
        zadj = 0.1f0 - 0.10273f0 * z + 0.00273f0 * z * z
        zadj < 0f0 && (zadj = 0f0)
        z = z + zadj
    end
    # young small-lodgepole HTG accelerator (Targhee), only cyc1 with known age
    if !(icyc > 1 || iage <= 0)
        ixage = iage                                      # IAGE + IY(ICYC) − IY(1) = IAGE at cyc1
        if ixage < 40 && ixage > 10 && d < 9.0f0 && z <= 2.0f0
            zadj = 0.3564f0 * dg * fint / yr
            closur = pct / 100.0f0
            relden < 100.0f0 && (closur = 1.0f0)
            zadj *= closur
            (iicr == 9 || iicr == 8) && (zadj *= 1.1f0)
            z += zadj
            z > 2.0f0 && (z = 2.0f0)
        end
    end
    dia = d + dg / bark
    ((xi1 + cof[1, keycr]) > dia) || return 0.1f0         # DIA exceeds range ⇒ HTG=0.1 (GO TO 185 inverted)
    psi = cof[8, keycr] * ((dia - xi1) / (xi1 + cof[1, keycr] - dia))^cof[9, keycr] *
          exp(z * ((1.0f0 - cof[7, keycr]^2))^0.5f0 / cof[6, keycr])
    hht = (psi / (1.0f0 + psi)) * cof[2, keycr] + xi2
    hht < h && (hht = h)
    htg = hht - h
    htg < 0.1f0 && (htg = 0.1f0)
    return htg
end

function height_growth!(s::StandState, ::SouthCentralOregon; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    sd = s.coef.species
    avh = p.avg_height; ba = p.basal_area
    relden = p.relative_density
    elev = p.elevation
    ifor = Int(p.forest_idx)
    icyc = Int(s.control.cycle) + 1
    fint = s.control.growth_fint; yr = htg_period(s.variant)
    dens = s.density
    @inbounds for i in 1:t.n
        t.ht_growth[i] = 0f0
        t.tpa[i] <= 0f0 && continue
        ispc = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]; dg = t.diam_growth[i]
        sindx = p.sp_site_index[ispc]
        icr = Int(t.crown_pct[i])
        htcon = c.htg_cor[ispc]
        pt_i = Int(t.plot_id[i])
        pccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        pct = t.crown_ratio[i]
        iage = t.age_known[i] ? round(Int, t.birth_age[i]) : 0
        sitage, sitht, agmax, htmax = so_findag(ispc, ifor, h, sindx)

        htg = 0f0; done_161 = false
        local pothtg::Float32 = 0f0; to_1320 = false; skip_mod = false
        if ispc == 9 || ispc == 27                        # SH/WO — CA Ritchie&Hann method
            if sitage > agmax
                pothtg = 0.10f0; to_1320 = true
            else
                hguess = so_htcalc(ifor, sindx, ispc, sitage + 10.0f0)
                pothtg = hguess - sitht
                cratio = Float32(icr) / 100.0f0
                relht = avh > 0f0 ? h / avh : 0f0
                relht > 1.0f0 && (relht = 1.0f0)
                pccf < 100.0f0 && (relht = 1.0f0)
                crmod = 1.0f0 - exp(-4.26558f0 * cratio)
                rhmod = exp(2.54119f0 * (relht^0.250537f0 - 1.0f0))
                xmod = 1.016605f0 * crmod * rhmod
                htg = pothtg * xmod
                htg < 0.1f0 && (htg = 0.1f0)
                skip_mod = true                            # GO TO 999
            end
        elseif ispc == 11                                 # JU — height from REGENT
            skip_mod = true                                # htg stays 0 → 999
        elseif ispc == 10                                 # PP — no H≥HTMAX gate (CASE(10))
            if sitage > agmax
                pothtg = -1.31f0 + 0.05f0 * sindx
                pothtg < 0.1f0 && (pothtg = 0.1f0)
            else
                hguess = so_htcalc(ifor, sindx, ispc, sitage + 10.0f0)
                pothtg = hguess - sitht
            end
            to_1320 = true
        else                                              # CASE DEFAULT (incl WB/AS 16,24)
            if h >= htmax
                htg = 0.1f0
                htg = scale * htg * exp(htcon)
                t.ht_growth[i] = htg; done_161 = true
            elseif ispc == 16 || ispc == 24               # Johnson-SBB (no 1320 modifier)
                bark = so_bratio(sd, ispc, d)
                htg = _so_sbb_htg(ispc, d, h, dg, bark, icr, elev, iage, icyc, pct, relden, fint, yr)
                skip_mod = true                            # GO TO 999
            elseif sitage > agmax
                pothtg = 0.10f0
                to_1320 = true
            else
                hguess = so_htcalc(ifor, sindx, ispc, sitage + 10.0f0)
                pothtg = hguess - sitht
                to_1320 = true
            end
        end
        done_161 && continue

        if to_1320 && !skip_mod                            # Hoerl CR + Chapman-Richards RH modifiers
            relht = avh > 0f0 ? h / avh : 0f0
            relht > 1.5f0 && (relht = 1.5f0)
            crf = Float32(icr) / 100.0f0
            hgmdcr = SO_HTGF_CRA * crf^SO_HTGF_CRB * exp(SO_HTGF_CRC * crf)
            hgmdcr > 1.0f0 && (hgmdcr = 1.0f0)
            rhx = relht
            fctrkx = (SO_HTGF_RHK / SO_HTGF_RHYXS[ispc])^(SO_HTGF_RHM[ispc] - 1.0f0) - 1.0f0
            fctrrb = -1.0f0 * (SO_HTGF_RHR[ispc] / (1.0f0 - SO_HTGF_RHB[ispc]))
            fctrxb = rhx^(1.0f0 - SO_HTGF_RHB[ispc]) - SO_HTGF_RHXS^(1.0f0 - SO_HTGF_RHB[ispc])
            fctrm = -1.0f0 / (SO_HTGF_RHM[ispc] - 1.0f0)
            hgmdrh = SO_HTGF_RHK * (1.0f0 + fctrkx * exp(fctrrb * fctrxb))^fctrm
            htgmod = 0.25f0 * hgmdcr + 0.75f0 * hgmdrh
            htgmod >= 2.0f0 && (htgmod = 2.0f0)
            htgmod <= 0.0f0 && (htgmod = 0.1f0)
            htg = pothtg * htgmod
        end

        # 999: cap at HTMAX, scale FINT/YR · XHT · exp(HTCON). XHT (XHMULT) = 1 (no GROWTH mult).
        # Applies to EVERY path reaching 999 (9/27, 11, 16/24, and the 1320 fall-through).
        temph = h + htg
        temph > htmax && (htg = htmax - h)
        htg < 0.1f0 && (htg = 0.1f0)
        htg = scale * htg * exp(htcon)
        t.ht_growth[i] = htg
    end
    return s
end
