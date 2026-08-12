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

"""
IE regent small-tree HEIGHT calibration (ie/regent.f:1138-1337, the LHTCAL/mode-40 pass). Computes the RAW
regent HCOR into `c.htg_cor_init[sp]` for NIVAR species: for each sub-5" tree with a measured height
increment, accumulate the predicted regent height growth (EDH = HK−H, HK grown over the subcycles by the
NIVAR model with HCOR=0) and the measured TERM=HTG·SCALE3; CORNEW = Σ(TERM·P)/Σ(EDH·P); HCOR_raw =
ln(CORNEW), trapped to [0.0821, 12.1825]. `calibrate_diameter_growth!`'s shared attenuation (dgdriv.f:188-194,
`htg_cor_small = dg_cor_goal + cormlt_h·(htg_cor_init − dg_cor_goal)`) then produces the applied HCOR.

Without this, IE NIVAR species had htg_cor_init=0 ⇒ the diameter COR (dg_cor_goal) leaked into the regent
height CON, over-growing the small-tree cohort on dense stands where the calibration fires (#171). Inert where
fewer than NCALHT(5) sub-5" trees carry a measured height increment (e.g. iet01) ⇒ htg_cor_init stays 0.
"""
function ie_regent_hcor_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector,
                              saved_dbh::AbstractVector)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    rhcon = ie_regcons!(s)                              # raw RHCON (no HCOR)
    # BACKDATED stand density (t.dbh is backdated here): the calibration predicts the PAST growth period, so
    # RDNEXT(1)/BANEXT(1) use the start-of-period BA/RELDEN summed over LIVE + RECENTLY-DEAD records (dense.f:79-86,
    # the notre.f-inflated dead added back at their backdated dbh), not the current live-only stand.
    ba = 0f0; relden = 0f0
    @inbounds for i in 1:(t.n + t.ndead)
        d = t.dbh[i]; pr = t.tpa[i]
        ba += 0.005454154f0 * d * d * pr
        relden += ie_tree_ccf(Int(t.species[i]), d) * pr
    end
    regyr = IE_RG_REGYR
    finth = s.control.growth_finth > 0f0 ? s.control.growth_finth : Float32(htg_period(s.variant))
    scale3 = regyr / finth                              # regent.f:1065 SCALE3 = REGYR/FINTH
    fint = s.control.growth_fint > 0f0 ? s.control.growth_fint : Float32(htg_period(s.variant))
    ntyr = Int(round(fint)); iyr = Int(regyr)
    nper = ntyr ÷ iyr; (ntyr % iyr != 0) && (nper += 1); nper < 1 && (nper = 1)
    kper = zeros(Int, nper); itot = ntyr; nn = nper
    @inbounds for i in 1:nper
        if nn == 1; kper[i] = itot; break; end
        kper[i] = itot ÷ nn; itot -= kper[i]; nn -= 1
    end
    # per-subcycle density from the large trees (inert at calibration: diam_growth≈0 ⇒ flat = ba/relden)
    banext = fill(ba, nper); rdnext = fill(relden, nper)
    if nper > 1
        @inbounds for i in 1:t.n
            d1 = t.dbh[i]; d1 < 3.0f0 && continue          # backdated dbh
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = ie_bratio(sp, d1)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            c1 = ie_tree_ccf(sp, d1); c2 = ie_tree_ccf(sp, d2)
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * 0.985f0^k
                rdnext[j] += k * ci * pn; banext[j] += k * bi * pn
            end
        end
    end
    @inbounds for sp in 1:23
        (sp <= 12 || sp == 14 || sp == 23) || continue    # NIVAR only
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        snx = 0f0; sny = 0f0; snp = 0f0; nh = 0
        for k in i1:i2
            i = ind1[k]
            saved_dbh[i] >= 5.0f0 && continue             # DBH<5 (regent.f:1159)
            hg = t.ht_growth[i]; hg < 0.001f0 && continue # measured HTG required
            hb = t.height[i] - hg; hb < 0.01f0 && continue # backdated H (IHTG<2)
            pct = t.crown_ratio[i]
            hk = hb
            for j in 1:nper
                bal = banext[j] * (100.0f0 - pct) * 0.0001f0
                htgrl = rhcon[sp] + IE_RG_RHLH[sp]*log(hk) + IE_RG_RHCCF[sp]*rdnext[j] +
                        IE_RG_RHBAL[sp]*bal
                hk += exp(htgrl)                          # regent.f:1174-1175 (NO scale in the calib pass)
            end
            edh = hk - hb                                 # regent.f NIVAR EDH = HK−H
            term = hg * scale3
            pr = t.tpa[i]
            snx += edh * pr; sny += term * pr; snp += pr; nh += 1
        end
        nh < 5 && continue                                # NCALHT
        snx /= snp; sny /= snp
        cornew = snx > 0f0 ? sny / snx : 1f0
        cornew <= 0f0 && (cornew = 1f-4)
        (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
        c.htg_cor_init[sp] = log(cornew)
    end
    return s
end

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
    # ---- final assembly (regent.f DO 30 ISPC=1,MAXSP): HTGR1 + ZZRAN + XWT blend toward the large-tree HTG,
    #      then the D<3 diameter dub. XWT blend is CRITICAL — trees D∈[XMIN,XMAX] blend to the (validated) htgf
    #      value, and D≥3 keep their large-tree DG (only D<3 gets the small-tree dub). FVS iterates SPECIES-
    #      SORTED (DO 30 ISPC; DO 25 I3=I1,I2 via IND1) — the per-record ZZRAN (BACHLO) draws MUST happen in
    #      this order or the RNG stream desyncs vs live on multi-species stands (CR proved this). ----
    _sp_order = sortperm(view(t.species, 1:n); alg = Base.Sort.MergeSort)   # stable ⇒ record order within sp
    @inbounds for oi in 1:n
        i = _sp_order[oi]
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
                # CRVAR CO DIAMETER (ie/regent.f:637-640): D<1 seedlings get the sub-breast-height nominal
                # D2=D+0.0001·H2 (H2≤4.5) else D2=D — NOT the large-tree DG, which over-extrapolates a 0.1"
                # seedling ⇒ dense hardwood-seedling BA over-growth (same bug fixed for EM CRVAR, 5116924).
                # D≥1 keeps the large-tree DG (Fortran `IF(D.GE.1.0)GO TO 15`).
                if d < 1.0f0
                    h2 = h + htg
                    d2 = h2 <= 4.5f0 ? d + 0.0001f0*h2 : d
                    dgnew = d2 - d; dgnew < 0f0 && (dgnew = 0f0)
                    t.diam_growth[i] = dgnew
                end
            end
            continue                                              # (all special species handled)
        end
        h = t.height[i]
        xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
        ax = IE_RG_HHT1[sp]; bx = IE_RG_HHT2[sp]
        htgr1 = wk3[i] - h; htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        cap = s.control.sp_size_cap[sp, 4]
        large_htg = t.ht_growth[i]                              # large-tree htgf value for the blend (before l=0 overwrites)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        small_d = d < 3.0f0
        diam = IE_RG_DIAM[sp]; bark = ie_bratio(sp, d)          # BARK at the pre-growth DBH (regent.f:983)
        # deterministic dub-bias terms (no ZZRAN) — computed once. DADJ (regent.f:866-872), D1 (regent.f:875-877).
        relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h - 4.5f0) / (ah - 4.5f0)
        relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
        dadj = delmax*relh*relh - 2.0f0*delmax*relh + 0.65f0
        d1v = diam + dadj; h > 4.5f0 && (d1v = ax * (h - 4.5f0)^bx + dadj)
        # TRIPLING (regent.f:801-810 "IF TRIPLING, EACH TRIPLE GETS A NEW RANDOM"): draw a FRESH ZZRAN per
        # tripled record; central (l=0) → the tree, copies (l=1,2) → the stash (htgU/htgL + is_small + dgU/dgL).
        # The DBH dub is the FAITHFUL non-ESTAB path (regent.f:955-989): DG(K)=(DK−D1)·XRDGRO·BARK on the DDS
        # scale — DBH is UNCHANGED and grows later via GRADD (t.dbh += DG/bark). (The old code did DBH-direct,
        # which is the ESTAB-only branch (regent.f:945-951 LESTB) — wrong for the growth cycle; combined with a
        # 1-draw-per-tree ZZRAN it desynced the RNG and doubled regen BA.) Only HK<4.5 sets DBH directly
        # (regent.f:881), a tiny-tree edge that does not fire on realistic small-tree stands (H+HTG > 4.5).
        nrec = stash !== nothing ? 3 : 1
        central_dbh = d
        for l in 0:(nrec - 1)
            zzran = 0.0f0
            if dgsd >= 1.0f0
                while true
                    zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                    (zzran <= 1.0f0 && zzran >= -1.5f0) && break
                end
            end
            htgr = htgr1 * exp(zzran * IE_RG_HSIGMA)             # NIVAR multiplicative randomization (regent.f:924)
            htg = htgr * (1.0f0 - xwt) + xwt * large_htg        # blend toward the large-tree htgf value
            (h + htg > cap) && (htg = max(cap - h, 0.1f0))
            # diameter dub (D<3 only). dg_inc = inside-bark increment (added to DBH via GRADD); dbh_dir≥0 ⇒ set DBH.
            dg_inc = 0.0f0; dbh_dir = -1.0f0
            if small_d
                hk = h + htg
                if hk < 4.5f0
                    dbh_dir = 0.1f0 + diam * 0.01f0 + hk * 0.001f0          # regent.f:881 (DBH set, DG=0)
                else
                    dk = ax * (hk - 4.5f0)^bx + dadj                        # regent.f:938
                    dk < diam && (dk = diam)
                    dk = dk + hk * 0.001f0
                    dgk = (dk - d1v) * xrdgro; dgk < 0.0f0 && (dgk = 0.0f0) # regent.f:958 (DK−D1)·XRDGRO
                    dg0 = dgk * bark                                        # DG(K)=DGK·BARK (regent.f:981)
                    dds = dg0 * (2.0f0*bark*d + dg0) * scale2               # regent.f:984 (FINT→10yr via SCALE2)
                    dg_inc = sqrt((d*bark)^2 + dds) - bark*d               # regent.f:985
                    (d + dg_inc) < diam && (dg_inc = diam - d)             # regent.f:987 DIAM floor
                    dg_inc = dg_bound(nothing, nothing, sp, d, dg_inc, s.control.sp_size_cap)  # DGBND (SIZCAP)
                end
            end
            if l == 0
                t.ht_growth[i] = htg
                if small_d
                    if dbh_dir >= 0.0f0
                        t.dbh[i] = dbh_dir; t.diam_growth[i] = 0.0f0; central_dbh = dbh_dir
                    else
                        t.diam_growth[i] = dg_inc                           # increment; DBH unchanged (=d)
                    end
                end
            elseif l == 1
                stash.htgU[i] = htg; stash.is_small[i] = true
                small_d && (stash.dgU[i] = dbh_dir >= 0.0f0 ? (dbh_dir - central_dbh)*bark : dg_inc)
            else
                stash.htgL[i] = htg
                small_d && (stash.dgL[i] = dbh_dir >= 0.0f0 ? (dbh_dir - central_dbh)*bark : dg_inc)
            end
        end
    end
    return s
end

# ie/esgent.f (CALL REGENT(.TRUE.,ITRNIN)) — grow the JUST-ESTABLISHED regen IN its birth cycle. IE was OMITTED
# from the esgent dispatch (simulate.jl had CR/TT/EM/UT/CI/BM), so planted/established IE seedlings never got their
# first-cycle height growth (IE BARE-PLANT: TopHt frozen at ~2 ft, BA 50-60% under). Same class as EM #137 /
# UT #184 / CI/BM #185. SINGLE birth-subperiod pass (esgent grows the new records one partial cycle, not the full
# NPER subcycle machinery): mirrors small_tree_growth!'s NIVAR height (HTGRL=CON+RHLH·lnH+RHCCF·RDEN+RHBAL·BAL) +
# AX/BX power H→D dub over the birth fraction (subyr=FINT−GENTIM=5=WK4), applying HT/DBH directly. NIVAR conifers
# (sp≤12,14,23) — the planted-conifer case. Non-NIVAR planted species (PI/JU 15,16 / TT 13,17 / CR 19,22) are rare
# as planting stock and left un-birth-grown here (would need their special-species branches; see #186 follow-up).
function ie_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    sd = s.coef.species
    nstart >= t.n && return s
    rhcon = ie_regcons!(s)
    ba = p.basal_area; relden = p.relative_density; avh = p.avg_height; ah = avh
    dgsd = s.control.dg_sd
    regyr = IE_RG_REGYR; yr = s.control.year
    ntyr = Int(round(fint))
    scale2 = ntyr > 0 ? yr / Float32(ntyr) : 1.0f0
    gentim = max(fint - 5.0f0, 0.0f0)
    bscale = (fint - gentim) / regyr                     # birth-cycle fraction (WK4; =0.5 for fint=10)
    cur_year = current_cycle_year(s)
    delmax = (ah / 36.0f0) * (0.01232f0 * relden - 1.75f0); delmax > 0.0f0 && (delmax = 0.0f0)
    @inbounds for i in (nstart+1):t.n
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= IE_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        (sp <= 12 || sp == 14 || sp == 23) || continue   # NIVAR conifers only (planted-conifer case)
        h = t.height[i]
        con = rhcon[sp] + c.htg_cor_small[sp]             # CON = RHCON + HCOR
        pct = t.crown_ratio[i]
        bal = ba * (100.0f0 - pct) * 0.0001f0
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h - 4.5f0) / (ah - 4.5f0)
        relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
        dadj = delmax*relh*relh - 2.0f0*delmax*relh + 0.65f0
        htgrl = con + IE_RG_RHLH[sp]*log(h) + IE_RG_RHCCF[sp]*relden + IE_RG_RHBAL[sp]*bal
        h2 = h + exp(htgrl) * bscale * xrhgro            # birth-cycle subperiod (was ·SCALE=kper/regyr)
        htgr1 = h2 - h; htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
        ax = IE_RG_HHT1[sp]; bx = IE_RG_HHT2[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        diam = IE_RG_DIAM[sp]; bark = ie_bratio(sp, d)
        d1v = diam + dadj; h > 4.5f0 && (d1v = ax * (h - 4.5f0)^bx + dadj)
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 1.0f0 && zzran >= -1.5f0) && break
            end
        end
        htgr = htgr1 * exp(zzran * IE_RG_HSIGMA)
        htg = htgr * (1.0f0 - xwt)                       # new tree: large-tree HTG(K)=0
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        hk = h + htg
        t.height[i] = hk; t.ht_growth[i] = htg
        if d < 3.0f0                                     # small-tree DBH dub (regent.f:938-987)
            if hk < 4.5f0
                t.dbh[i] = 0.1f0 + diam * 0.01f0 + hk * 0.001f0; t.diam_growth[i] = 0.0f0
            else
                dk = ax * (hk - 4.5f0)^bx + dadj; dk < diam && (dk = diam); dk = dk + hk * 0.001f0
                dgk = (dk - d1v) * xrdgro; dgk < 0.0f0 && (dgk = 0.0f0)
                dg0 = dgk * bark
                dds = dg0 * (2.0f0*bark*d + dg0) * scale2
                dg_inc = sqrt((d*bark)^2 + dds) - bark*d
                (d + dg_inc) < diam && (dg_inc = diam - d)
                dg_inc = dg_bound(nothing, nothing, sp, d, dg_inc, s.control.sp_size_cap)
                dg_inc > 0.0f0 && (t.dbh[i] = d + dg_inc; t.diam_growth[i] = dg_inc)
            end
        end
    end
    return s
end
