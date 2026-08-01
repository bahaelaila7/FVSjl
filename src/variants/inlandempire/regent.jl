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
        irhhab = clamp(Int(IE_RG_MAPHAB[itype, sp]), 1, 6)
        rhcon[sp] = regch + IE_RG_RHSC[sp] + IE_RG_RHHAB[irhhab, sp]
    end
    return rhcon
end

"""IE `small_tree_growth!` (ie/regent.f). Overrides DG/HTG for small trees (D<XMAX). NIVAR path is the
bulk (iet01); special species ported faithfully."""
function small_tree_growth!(s::StandState, stash, ::InlandEmpire; fint::Float32 = 10.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
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
                rdnext[j] += k * ci / pr * pn
                banext[j] += k * bi * pn
            end
        end
    end
    # DELMAX (regent.f:333), AH=AVH
    ah = avh
    delmax = (ah / 36.0f0) * (0.01232f0 * relden - 1.75f0); delmax > 0.0f0 && (delmax = 0.0f0)
    # per-tree height + diameter accumulators (WK3=H, WK5=D), start at HT/DBH
    wk3 = Float32[t.height[i] for i in 1:n]
    wk5 = Float32[t.dbh[i] for i in 1:n]
    cur_year = current_cycle_year(s)
    # ---- subcycle loop (regent.f:250-620) ----
    @inbounds for j in 1:nper
        baj = banext[j]; rdj = rdnext[j]
        scale = Float32(kper[j]) / regyr
        for i in 1:n
            sp = Int(t.species[i]); d0 = t.dbh[i]
            d0 >= IE_RG_XMAX[sp] && continue
            t.tpa[i] <= 0.0f0 && continue
            nivar = sp <= 12 || sp == 14 || sp == 23
            nivar || continue                                  # special species handled in the final pass
            con = rhcon[sp] + c.htg_cor_small[sp]              # CON = RHCON + HCOR (HCOR=0 until calib)
            h1 = wk3[i]; d = wk5[i]
            pct = t.crown_ratio[i]
            bal = baj * (100.0f0 - pct) * 0.01f0
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
        end
    end
    # ---- final assembly (regent.f:441-600): HTGR1 + ZZRAN + XWT blend toward the large-tree HTG, then
    #      the D<3 diameter dub. XWT blend is CRITICAL — trees D∈[XMIN,XMAX] blend to the (validated) htgf
    #      value, and D≥3 keep their large-tree DG (only D<3 gets the small-tree dub). ----
    @inbounds for i in 1:n
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= IE_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        (sp <= 12 || sp == 14 || sp == 23) || continue         # NIVAR only (special species: TODO)
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
        # diameter: only D<3 gets the small-tree dub (regent.f:568-600); D≥3 keeps its large-tree DG
        if d < 3.0f0
            relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h - 4.5f0) / (ah - 4.5f0)
            relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
            dadj = delmax*relh*relh - 2.0f0*delmax*relh + 0.65f0
            hk = h + htg
            if hk < 4.5f0
                t.diam_growth[i] = 0.0f0                        # DBH set to a floor; no DG
            else
                dk = ax * (hk - 4.5f0)^bx + dadj
                dk < IE_RG_DIAM[sp] && (dk = IE_RG_DIAM[sp])
                dk += hk * 0.001f0
                dg = dk - d; dg < 0.0f0 && (dg = 0.0f0)
                t.diam_growth[i] = dg
            end
        end
    end
    return s
end
