# =============================================================================
# regent.jl (inlandempire) — IE small-tree growth (ie/regent.f, chunk 6). Coefficients in
# regent_coefficients.jl (dumped from live). Reuses the KT subcycle scaffold (kper/banext/rdnext/DELMAX);
# the equations are IE-specific (measured — DIFFER from KT).
#
#   ie_regcons!(s)            — REGCON: RHCON(sp) = REGCH + RHSC(sp) + RHHAB(MAPHAB(ITYPE,sp),sp),
#                               REGCH = RHGL(IGL) + (RSAB0+RSAB1·cosASP+RSAB2·sinASP)·SLOPE.
#   small_tree_growth!(...)   — NPER subcycles; NIVAR height HTGRL=CON+RHLH·lnH+RHCCF·RDJ+RHBAL·BAL,
#                               H2=H1+exp(HTGRL)·SCALE·XRHGRO·WK4; DG via AX/BX (=HHT1/HHT2) power H-D model.
# NIVAR (sp1:12,14,23) fully ported+validated; TTVAR(13,17)/CRVAR(19,22)/UTVAR(15,16,18,20,21) faithful
# (stochastic BETA/POTHTG/aspen-FINDAG) — not in iet01, validated later.
# =============================================================================

"""IE REGCON (ie/regent.f:1479): per-species small-tree height constant RHCON."""
function ie_regcons!(s::StandState)
    p = s.plot
    itype = Int(p.habitat_input); (itype < 1 || itype > 30) && (itype = 1)
    igl = Int(p.geo_location); (igl < 1 || igl > 3) && (igl = 2)
    asp = p.aspect; slope = p.slope
    regch = IE_RG_RHGL[igl] + (IE_RG_RSAB[1] + IE_RG_RSAB[2]*cos(asp) + IE_RG_RSAB[3]*sin(asp)) * slope
    rhcon = Vector{Float32}(undef, 23)
    @inbounds for sp in 1:23
        # regent.f:1479-1492 — the REGCH+RHSC+RHHAB formula is NIVAR-ONLY; every non-NIVAR (TT/CR/UT)
        # species gets RHCON = 1.0 (their CON = 1.0·EXP(HCOR); the site effect rides in via HCOR calib).
        if sp <= 12 || sp == 14 || sp == 23
            irhhab = clamp(Int(IE_RG_MAPHAB[itype, sp]), 1, 6)
            rhcon[sp] = regch + IE_RG_RHSC[sp] + IE_RG_RHHAB[irhhab, sp]
        else
            rhcon[sp] = 1.0f0
        end
    end
    return rhcon
end

"""IE `small_tree_growth!` (ie/regent.f). Overrides DG/HTG for small trees (D<XMAX). NIVAR path is the
bulk (iet01); special species ported faithfully."""
function small_tree_growth!(s::StandState, stash, ::InlandEmpire; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    sd = s.coef.species                                 # blkdat HT-DBH :ht1/:ht2 for the aspen log-DK
    t.n == 0 && return s
    n = t.n
    rhcon = ie_regcons!(s)
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height
    dgsd = s.control.dg_sd
    regyr = IE_RG_REGYR
    yr = s.control.year
    # subcycle count + lengths (regent.f:186-203) — reuse KT scaffold
    ntyr = Int(round(fint)); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    scale2 = ntyr > 0 ? yr / Float32(ntyr) : 1.0f0
    # per-subcycle stand density from the LARGE trees (regent.f:236-256) — reuse KT scaffold
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    if nper > 1
        @inbounds for i in 1:n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = ie_bratio(sp, d1)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            c1 = ie_tree_ccf(sp, d1); c2 = ie_tree_ccf(sp, d2)
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]
                pn = pr * 0.985f0^k
                # regent.f:269 RDNEXT += K*CI/P*PN where CI uses CCFCAL (= CCFT*P). Our ie_tree_ccf is CCFT
                # (no P), so CI here lacks the P that FVS's /P cancels ⇒ do NOT divide by pr (pn carries it).
                rdnext[j] += k * ci * pn
                banext[j] += k * bi * pn
            end
        end
    end
    # DELMAX (regent.f:333), AH=AVH
    ah = avh
    delmax = (ah / 36.0f0) * (0.01232f0 * relden - 1.75f0); delmax > 0.0f0 && (delmax = 0.0f0)
    # PCTRED (regent.f:338-343) — CR/UT density modifier from CCF·top-height; used by the special-species
    # potential-height model (PI/JU here). Stand-level, computed once (does not vary by subcycle).
    xpr = ah * (relden / 100.0f0); xpr > 300.0f0 && (xpr = 300.0f0)
    pctred = 1.11436f0 + xpr*(-0.011493f0 + xpr*(0.43012f-4 + xpr*(-0.72221f-7 +
             xpr*(0.5607f-10 - xpr*0.1641f-13))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    scale_ut = yr > 0f0 ? Float32(ntyr) / yr : 1.0f0    # CR/UT SCALE = NTYR/YR (regent.f:393/397)
    # per-tree height + diameter accumulators (WK3=H, WK5=D), start at HT/DBH
    wk3 = Float32[t.height[i] for i in 1:n]
    wk5 = Float32[t.dbh[i] for i in 1:n]
    zrand_tt = fill(-999f0, n)          # TTVAR (sp13/17) persistent ZRAND per tree (regent.f:513); drawn fresh
                                        # each call (cross-cycle persistence = accepted ZZRAN-class residual)
    cur_year = current_cycle_year(s)
    # ---- subcycle loop (regent.f:250-620) ----
    @inbounds for j in 1:nper
        baj = banext[j]; rdj = rdnext[j]
        scale = Float32(kper[j]) / regyr
        ky = 0; for jj in 1:j; ky += kper[jj]; end            # KY = cumulative years thru subcycle j (regent.f:353)
        surv = 0.985f0 ^ ky
        for i in 1:n
            sp = Int(t.species[i]); d0 = t.dbh[i]
            d0 >= IE_RG_XMAX[sp] && continue
            t.tpa[i] <= 0.0f0 && continue
            nivar = sp <= 12 || sp == 14 || sp == 23
            if !nivar
                # UTVAR PI/JU (sp15,16): potential-height model, ONE pass only (regent.f:407 J>1 skip).
                if (sp == 15 || sp == 16) && j == 1
                    con = rhcon[sp] * exp(c.htg_cor_small[sp])   # non-NIVAR CON = RHCON·EXP(HCOR) (regent.f:414)
                    h1 = wk3[i]; d = wk5[i]
                    sitear = p.sp_site_index[sp]
                    xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
                    xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
                    sj = sitear                                   # POTHTG uses raw SITEAR (regent.f:459/467); H==H1 at j=1
                    pothtg = ((sj/5f0)*(sj*1.5f0 - h1)/(sj*1.5f0)) * 0.83f0
                    crx = Float32(t.crown_pct[i]) / 100f0
                    vigor = 150f0*crx^3*exp(-6f0*crx) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
                    vigor = 1f0 - (1f0 - vigor)/3f0               # PI/JU cut vigor by 2/3 (regent.f:552)
                    htgrl = pothtg * pctred * vigor * con
                    h2 = h1 + htgrl * scale_ut                    # regent.f:600 (CR/UT else branch)
                    wk3[i] = h2
                    # subcycle diameter for density feedback (regent.f:646-651)
                    d2 = h2 <= 4.5f0 ? d + 0.001f0*h2 : max((h2 - 4.5f0)*10f0/(sitear - 4.5f0), 0.1f0)
                    wk5[i] = d2
                    if j < nper
                        pr = t.tpa[i]
                        c1 = ie_tree_ccf(sp, d); c2 = ie_tree_ccf(sp, d2)
                        rdnext[j+1] += Float32(ky) * pr * (c2 - c1) / 10.0f0 * surv
                        banext[j+1] += (0.005454154f0*d2*d2 - 0.005454154f0*d*d) * pr * surv
                    end
                elseif (sp == 18 || sp == 20 || sp == 21) && j == 1
                    # UTVAR aspen: Sheppard height curve (regent.f:556-582), ONE pass. CON=RHCON(=1)·EXP(HCOR).
                    con = rhcon[sp] * exp(c.htg_cor_small[sp])
                    h1 = wk3[i]
                    si = p.sp_site_index[sp]
                    si > IE_RG_SHI[sp] && (si = IE_RG_SHI[sp])
                    si <= IE_RG_SLO[sp] && (si = IE_RG_SLO[sp] + 0.5f0)
                    relsi = (si - IE_RG_SLO[sp]) / (IE_RG_SHI[sp] - IE_RG_SLO[sp])
                    rsimod = 0.5f0 * (1f0 + relsi)
                    sitage = (h1 * 12f0 * 2.54f0 / 26.9825f0) ^ 0.8509f0    # FINDAG (inverse Sheppard, regent.f:565)
                    hite1 = 26.9825f0 * sitage ^ 1.1752f0
                    hite2 = 26.9825f0 * (sitage + 10f0) ^ 1.1752f0
                    htgrl = (hite2 - hite1) / (2.54f0 * 12f0) * rsimod * con * 0.75f0
                    wk3[i] = h1 + htgrl * scale_ut                          # regent.f:600 (·SCALE=NTYR/YR)
                    # (aspen subcycle DBH/density-feedback omitted — single UT pass; final assembly is authoritative)
                elseif sp == 13 || sp == 17
                    # TTVAR (LM/PY): BETA/ZRAND height — EVERY subcycle (NOT one-pass; regent.f:507-538,598).
                    con = rhcon[sp] * exp(c.htg_cor_small[sp])
                    h1 = wk3[i]
                    cr = Float32(t.crown_pct[i])
                    tpccf = clamp(relden, 25f0, 300f0)                      # PCCF≈stand CCF (single-point)
                    beta1 = exp(1.17527f0 - 0.42124f0*log(tpccf))
                    beta2 = exp(-2.56002f0 - 0.58642f0*log(tpccf))
                    htg1 = beta1 + beta2*cr
                    stddev = htg1*(1.08720f0 - 0.00230f0*cr)
                    if zrand_tt[i] == -999f0                                # draw once (regent.f:513, no DGSD gate)
                        while true
                            z = bachlo(s.rng, 0.0f0, 1.0f0)
                            (z >= -2.0f0 && z <= 2.0f0) && (zrand_tt[i] = z; break)
                        end
                    end
                    htgrl = htg1 + zrand_tt[i]*stddev
                    (htgrl <= 0.1f0) && (htgrl = 0.1f0; zrand_tt[i] = -999f0)   # reset ⇒ redraw next subcycle
                    xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
                    wk3[i] = h1 + htgrl * scale * xrhgro * con              # ·SCALE=kper/regyr ·CON ·WK4(=1)
                    # DLESS3 diameter into wk5 (= DK in final assembly)
                    h2 = wk3[i]
                    if h2 > 4.5f0
                        hl4 = h2 - 4.5f0
                        dless3 = 0.000231f0*hl4*cr - 0.00005f0*hl4*tpccf + 0.001711f0*cr + 0.17023f0*hl4
                        wk5[i] = max(dless3 + 0.3f0, IE_RG_DIAM[sp])
                    end
                elseif (sp == 19 || sp == 22) && j == 1
                    # CRVAR CO (sp19,22): POTHTG height (VIGOR NOT cut, unlike PI/JU); ONE pass. Diameter stays
                    # large-tree dgf (CO CR-logic works; CO TPA already tracks live — no seedling-DG problem).
                    con = rhcon[sp] * exp(c.htg_cor_small[sp])
                    h1 = wk3[i]
                    sj = p.sp_site_index[sp]                       # POTHTG uses raw SITEAR
                    pothtg = ((sj/5f0)*(sj*1.5f0 - h1)/(sj*1.5f0)) * 0.83f0
                    crx = Float32(t.crown_pct[i]) / 100f0
                    vigor = 150f0*crx^3*exp(-6f0*crx) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
                    wk3[i] = h1 + pothtg * pctred * vigor * con * scale_ut
                end
                continue                                          # (all special species handled)
            end
            con = rhcon[sp] + c.htg_cor_small[sp]              # CON = RHCON + HCOR (HCOR=0 until calib)
            h1 = wk3[i]; d = wk5[i]
            pct = t.crown_ratio[i]
            bal = baj * (100.0f0 - pct) * 0.0001f0             # ie/regent.f:442 — NOTE 0.0001 (NOT KT's 0.01)
            xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
            xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
            relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h1 - 4.5f0) / (ah - 4.5f0)
            relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
            dadj = delmax*relh*relh - 2.0f0*delmax*relh + 0.65f0
            htgrl = con + IE_RG_RHLH[sp]*log(h1) + IE_RG_RHCCF[sp]*rdj + IE_RG_RHBAL[sp]*bal
            h2 = h1 + exp(htgrl) * scale * xrhgro            # WK4≈1 (healthy; damage factor omitted)
            wk3[i] = h2
            # NIVAR diameter (regent.f:598-610): skip if last subcycle or D≥3 or H2≤4.5
            d2 = d
            if !(j >= nper || d >= 3.0f0 || h2 <= 4.5f0)
                ax = IE_RG_HHT1[sp]; bx = IE_RG_HHT2[sp]
                d1v = IE_RG_DIAM[sp] + dadj
                h1 > 4.5f0 && (d1v = ax * (h1 - 4.5f0)^bx + dadj)
                d2v = ax * (h2 - 4.5f0)^bx + dadj
                dgj = (d2v - d1v) * xrdgro; dgj < 0.0f0 && (dgj = 0.0f0)
                d2 = d + dgj
            end
            wk5[i] = d2
            # small-tree density feedback into the NEXT subcycle (regent.f:665-666): the growing small tree
            # adds its CCF/BA increase to RDNEXT/BANEXT(J+1). C1/C2 use CCFCAL (= CCFT·P) ⇒ ·pr; no /P here.
            if j < nper
                pr = t.tpa[i]
                c1 = ie_tree_ccf(sp, d); c2 = ie_tree_ccf(sp, d2)
                rdnext[j+1] += Float32(ky) * pr * (c2 - c1) / 10.0f0 * surv
                banext[j+1] += (0.005454154f0*d2*d2 - 0.005454154f0*d*d) * pr * surv
            end
        end
    end
    # ---- final assembly (regent.f:441-600): HTGR1 + ZZRAN + XWT blend toward the large-tree HTG, then
    #      the D<3 diameter dub. XWT blend is CRITICAL — trees D∈[XMIN,XMAX] blend to the (validated) htgf
    #      value, and D≥3 keep their large-tree DG (only D<3 gets the small-tree dub). ----
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= IE_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        if !(sp <= 12 || sp == 14 || sp == 23)
            # UTVAR PI/JU (sp15,16): height increment + ZZRAN + linear height→DBH (regent.f:756-985).
            # regent OVERRIDES the large-tree dgf/htgf for pinyon/juniper (XMAX=99, XMIN=90 ⇒ XWT≡0).
            if sp == 15 || sp == 16
                h = t.height[i]
                sitear = p.sp_site_index[sp]
                xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
                xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
                htgr1 = wk3[i] - h                                # HK−H (UT: NOT floored; regent.f:756)
                zzran = 0f0
                if dgsd >= 1.0f0
                    while true
                        zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                        (zzran <= 0.5f0 && zzran >= -2.0f0) && break   # CR/UT bound (regent.f:813)
                    end
                end
                htgr = (htgr1 + zzran*0.1f0) * xrhgro              # UT: ZZRAN·0.1 (regent.f:814)
                htgr < 0.1f0 && (htgr = 0.1f0)
                xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
                xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)       # ≡0 for PI/JU (xmn=90)
                htg = htgr*(1f0 - xwt) + xwt*t.ht_growth[i]
                cap = s.control.sp_size_cap[sp, 4]
                (h + htg > cap) && (htg = max(cap - h, 0.1f0))
                t.ht_growth[i] = htg
                # diameter: linear height→DBH; DG on the DDS scale (regent.f:879-985)
                hk = h + htg
                if hk < 4.5f0
                    t.diam_growth[i] = 0f0                          # regent.f:883 DG(K)=0
                else
                    dk = (hk - 4.5f0)*10f0/(sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                    dkk = (h - 4.5f0)*10f0/(sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
                    h < 4.5f0 && (dkk = d)                          # regent.f:897 override
                    dgk = (dk - dkk) * ie_bratio(sp, d) * xrdgro    # regent.f:960
                    dgmx = IE_RG_DGMAX[sp]; dgk > dgmx && (dgk = dgmx)
                    dgk < 0f0 && (dgk = 0f0)
                    bark = ie_bratio(sp, d)
                    dds = dgk*(2f0*bark*d + dgk)*scale2             # regent.f:980 (DG(K)=DGK for CR/UT)
                    dgv = sqrt((d*bark)^2 + dds) - bark*d           # regent.f:981
                    (d + dgv) < IE_RG_DIAM[sp] && (dgv = IE_RG_DIAM[sp] - d)
                    t.diam_growth[i] = dgv
                end
            elseif sp == 18 || sp == 20 || sp == 21
                # UTVAR aspen: height increment + ZZRAN + XWT blend, then log-DK diameter (regent.f:756-985).
                # XMAX=4, XMIN=2 ⇒ height blends toward large-tree htgf on D∈[2,4]; diameter only D<3 (else DGFASP).
                h = t.height[i]
                xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
                xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
                htgr1 = wk3[i] - h
                zzran = 0f0
                if dgsd >= 1.0f0
                    while true
                        zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                        (zzran <= 0.5f0 && zzran >= -2.0f0) && break       # CR/UT bound (regent.f:813)
                    end
                end
                htgr = (htgr1 + zzran*0.1f0) * xrhgro
                htgr < 0.1f0 && (htgr = 0.1f0)
                xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
                xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
                htg = htgr*(1f0 - xwt) + xwt*t.ht_growth[i]
                cap = s.control.sp_size_cap[sp, 4]
                (h + htg > cap) && (htg = max(cap - h, 0.1f0))
                t.ht_growth[i] = htg
                # diameter: aspen log-DK, ONLY D<3 (regent.f:860 D≥3 ⇒ GO TO 23, keeps large-tree DGFASP)
                if d < 3f0
                    hk = h + htg
                    bx = sd[:ht2][sp]                                       # blkdat Wykoff HT2
                    ax = c.ht_dbh_iabflg[sp] == 1 ? sd[:ht1][sp] : c.ht_dbh_aa[sp]   # regent.f:900-904
                    if hk < 4.5f0
                        t.diam_growth[i] = 0f0
                    else
                        dk = (bx / (log(hk - 4.5f0) - ax)) - 1f0; dk < 0.1f0 && (dk = 0.1f0)   # regent.f:905-906
                        dkk = h <= 4.5f0 ? d : (bx / (log(h - 4.5f0) - ax)) - 1f0              # regent.f:907-911 (no DKK floor)
                        bark = ie_bratio(sp, d)
                        dgk = (dk - dkk) * bark * xrdgro                                       # regent.f:960
                        dgmx = IE_RG_DGMAX[sp]; dgk > dgmx && (dgk = dgmx)
                        dgk < 0f0 && (dgk = 0f0)
                        dds = dgk*(2f0*bark*d + dgk)*scale2                                    # regent.f:980
                        dgv = sqrt((d*bark)^2 + dds) - bark*d                                  # regent.f:981
                        (d + dgv) < IE_RG_DIAM[sp] && (dgv = IE_RG_DIAM[sp] - d)
                        t.diam_growth[i] = dgv
                    end
                end
            elseif sp == 13 || sp == 17
                # TTVAR (LM/PY): HTGR=HTGR1 (regent.f:800 skips ZZRAN for TT) + XWT blend; DLESS3 DK−DKK diameter.
                h = t.height[i]
                xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
                htgr = wk3[i] - h                                    # HTGR1 (no ZZRAN block for TTVAR)
                xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
                xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
                htg = htgr*(1f0 - xwt) + xwt*t.ht_growth[i]
                cap = s.control.sp_size_cap[sp, 4]
                (h + htg > cap) && (htg = max(cap - h, 0.1f0))
                t.ht_growth[i] = htg
                if d < 3f0
                    hk = h + htg
                    if hk >= 4.5f0                                   # regent.f:921
                        cr = Float32(t.crown_pct[i]); tpccf = clamp(relden, 25f0, 300f0)
                        hl4 = h - 4.5f0                              # DKK from current height HT(K) (regent.f:922)
                        dless3 = 0.000231f0*hl4*cr - 0.00005f0*hl4*tpccf + 0.001711f0*cr + 0.17023f0*hl4
                        dkk = max(dless3 + 0.3f0, IE_RG_DIAM[sp])
                        dgk = (wk5[i] - dkk) * xrdgro                # DK(subcycle wk5) − DKK (regent.f:927)
                        dgmx = IE_RG_DGMAX[sp] * (fint / 10f0)       # TTVAR DGMX = FINT·DGMAX (regent.f:736)
                        dgk > dgmx && (dgk = dgmx); dgk < 0f0 && (dgk = 0f0)
                        bark = ie_bratio(sp, d)
                        dg0 = dgk * bark                            # TTVAR DG(K)=DGK·BARK (regent.f:976)
                        dds = dg0*(2f0*bark*d + dg0)*scale2
                        dgv = sqrt((d*bark)^2 + dds) - bark*d
                        (d + dgv) < IE_RG_DIAM[sp] && (dgv = IE_RG_DIAM[sp] - d)
                        t.diam_growth[i] = dgv
                    else
                        t.diam_growth[i] = 0f0
                    end
                end
            elseif sp == 19 || sp == 22
                # CRVAR CO: HTGR=(HTGR1+ZZRAN·0.2)·XRHGRO (regent.f:815), XWT blend [0.5,2]. Diameter = dgf (not overridden).
                h = t.height[i]
                xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
                htgr1 = wk3[i] - h
                zzran = 0f0
                if dgsd >= 1.0f0
                    while true
                        zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                        (zzran <= 0.5f0 && zzran >= -2.0f0) && break
                    end
                end
                htgr = (htgr1 + zzran*0.2f0) * xrhgro                # CRVAR: ZZRAN·0.2
                htgr < 0.1f0 && (htgr = 0.1f0)
                xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
                xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
                htg = htgr*(1f0 - xwt) + xwt*t.ht_growth[i]
                htg < 0.1f0 && (htg = 0.1f0)                         # CRVAR HTG floor (regent.f:842)
                cap = s.control.sp_size_cap[sp, 4]
                (h + htg > cap) && (htg = max(cap - h, 0.1f0))
                t.ht_growth[i] = htg
            end
            continue                                              # (all special species handled)
        end
        h = t.height[i]
        xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
        ax = IE_RG_HHT1[sp]; bx = IE_RG_HHT2[sp]
        htgr1 = wk3[i] - h; htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 1.0f0 && zzran >= -1.5f0) && break
            end
        end
        htgr = htgr1 * exp(zzran * IE_RG_HSIGMA)                # NIVAR multiplicative randomization (regent.f:533)
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        htg = htgr * (1.0f0 - xwt) + xwt * t.ht_growth[i]        # blend toward the large-tree htgf value
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        # diameter dub for D<3. We set DBH=grown-DK directly + zero diam_growth ⇒ mortality sees the grown DBH.
        # RESIDUAL RESOLVED (2026-08-01, live regent trace): the small-tree base models are CORRECT. For the stand-1
        # sp3/d=1.2/H=11 tree (identical jl↔live population, valid pre-mortality match), jl's pre-ZZRAN height growth
        # htgr1=8.61 MATCHES live HTGR1=8.59 (subcycle HTGRL 1.493/1.420, CON=1.486, HCOR=0); the final HK differs
        # (jl 17.68 vs live 24.54) ONLY because the ZZRAN draw differs (line 142: live drew +0.455, jl negative). That
        # is the ZZRAN RNG stream-order difference — the documented ACCEPTED residual, NOT a bug. Earlier "DG formula
        # under-sizing / tripling-order / 4× granularity" hypotheses were all wrong turns from mismatched comparisons.
        # DBH-old+increment is FVS's true structure (live dump: DBH(K)==D) but restructuring regressed the .sum via
        # the same ZZRAN-driven HK spread; DK-direct is oracle-.sum-closest (matches live BA=115 @2000) so kept.
        if d < 3.0f0
            relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h - 4.5f0) / (ah - 4.5f0)
            relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
            dadj = delmax*relh*relh - 2.0f0*delmax*relh + 0.65f0
            hk = h + htg
            new_dbh = if hk < 4.5f0
                0.1f0 + IE_RG_DIAM[sp] * 0.01f0 + hk * 0.001f0   # regent.f:881
            else
                dk = ax * (hk - 4.5f0)^bx + dadj                 # regent.f:938 (DK<DIAM→DIAM; +HK*.001)
                dk < IE_RG_DIAM[sp] && (dk = IE_RG_DIAM[sp])
                dk + hk * 0.001f0
            end
            new_dbh < d && (new_dbh = d)                         # no shrink
            t.dbh[i] = new_dbh
            t.diam_growth[i] = 0f0                               # DBH already applied ⇒ GRADD-apply is a no-op
        end
    end
    return s
end
