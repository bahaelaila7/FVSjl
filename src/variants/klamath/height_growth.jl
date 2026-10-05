# =============================================================================
# height_growth.jl (klamath) — NC large-tree height growth (nc/htgf.f + findag.f + htcalc.f). Chunk 4.
#   nc_htcalc(si,sp,ag)       — site-index height-at-age curve (6 species-group branches, htcalc.f)
#   nc_findag(h,sp,si)        — solve site age from height (findag.f: step AG by 2 to match the curve)
#   height_growth!(::Klamath) — redwood LTHTG special + DEFAULT FINDAG→POTHTG (site-curve 5-yr rise)·XMOD, scaled
# =============================================================================

"nc/htcalc.f — site-index height (HGUESS) at age AG for species sp, site index si."
@inline function nc_htcalc(si::Float32, sp::Int, ag::Float32)
    if sp == 1 || sp == 3 || sp == 12                 # OS/DF/RW — Chapman
        z = 2500.0f0 / (si - 4.5f0)
        a = -0.954038f0 + 0.109757f0 * z
        b = 0.055818f0 + 0.0079224f0 * z
        c = -0.0007338f0 + 0.0001977f0 * z
        return (ag * ag) / (a + b * ag + c * ag * ag) + 4.5f0
    elseif sp == 4 || sp == 6 || sp == 9              # WF/IC/RF
        x1 = 38.0202f0 * fpow(ag, -1.05213f0) * fexp(0.009557f0 * ag)
        x2 = 101.842894f0 * (1.0f0 - fexp(-0.001442f0 * fpow(ag, 1.679259f0)))
        return (si - 69.91f0 + x1 * x2) / x1 + 4.5f0
    elseif sp == 5                                    # MA
        return si / (0.375f0 + 31.233f0 / ag)
    elseif sp == 7                                    # BO
        aa = sqrt(ag) - sqrt(50.0f0)
        return (si * (1.0f0 + 0.322f0 * aa) - 6.413f0 * aa) * 0.80f0
    elseif sp == 8 || sp == 11                        # TO/OH
        return si / (0.204f0 + 39.787f0 / ag) * 0.85f0
    else                                              # SP(2)/PP(10)
        return (1.88f0 * si - 7.178f0) * fpow(1.0f0 - fexp(-0.025f0 * ag), 0.001f0 * si + 1.64f0)
    end
end

"nc/findag.f — solve (SITAGE, SITHT) from tree height H via the site curve (AGMAX=200, HTMAX=300)."
@inline function nc_findag(h::Float32, sp::Int, si::Float32)
    agmax = 200.0f0; htmax = 300.0f0; toler = 2.0f0
    h >= htmax && return (agmax + (h - htmax) / 0.10f0, h)
    ag = 2.0f0; hguess = 0.0f0; incrng = false
    while true
        oldhg = hguess
        hguess = nc_htcalc(si, sp, ag)
        if hguess >= 1.0f0
            (abs(hguess - h) <= toler || h < hguess) && return (ag, hguess)
            diff = hguess - oldhg
            (oldhg != 0.0f0 && diff >= 0.05f0) && (incrng = true)
            (incrng && diff < 0.05f0) && return (ag, hguess)
        end
        ag += 2.0f0
        ag > agmax && return (agmax, h)
    end
end

"NC height_growth! — nc/htgf.f (per-tree HTG; redwood LTHTG special, DEFAULT FINDAG→site-curve POTHTG·XMOD)."
function height_growth!(s::StandState, ::Klamath; scale::Float32 = 1.0f0)
    p, t, c = s.plot, s.trees, s.calib
    dens = s.density; sd = s.coef.species
    avh = p.avg_height; slope = p.slope; asp = p.aspect
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0.0f0 && continue
        sp = Int(t.species[i]); h = t.height[i]; d = t.dbh[i]
        dglt = t.diam_growth[i]
        si = p.sp_site_index[sp]
        pt_i = Int(t.plot_id[i])
        pccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        htcon = c.htg_cor[sp]                          # HTCON calib (0 baseline)
        local htg::Float32
        if sp == 12                                    # REDWOOD
            brat = nc_bratio(sd[:bark1][sp], sd[:bark2][sp], Int(sd[:bark_imap][sp]), d)
            # nc/htgf.f:156-161 — DGLT (the per-cycle large-tree DG) is scaled to a 10-YEAR increment
            # (DGLT·2) BEFORE the outside-bark conversion DG10 = DGLT/BRAT. The ·2 was MISSING ⇒ DG10 half
            # of FVS ⇒ redwood large-tree HTG understated (log(DG10) term), leaving RW stands short even
            # after the site-index fix.
            dg10 = h < 4.5f0 ? 0.1f0 : (dglt * 2f0) / brat
            lthtg = fexp(1.412947f0 - 0.000204f0 * (d * d) + 0.31971f0 * flog(d) +
                        0.394005f0 * flog(si) + 0.399888f0 * flog(dg10) - 0.451708f0 * flog(h)) * 0.5f0
            hgbnd = h < 217.0f0 ? 1.0f0 :
                    h < 380.0f0 ? max(1.0f0 - (h - 217.0f0) / (380.0f0 - 217.0f0), 0.1f0) : 0.1f0
            htg = lthtg * hgbnd
        else                                           # DEFAULT — FINDAG + site-curve POTHTG
            sitage, sitht = nc_findag(h, sp, si)
            htmax = 300.0f0; agmax = 200.0f0
            if h >= htmax
                htg = 0.1f0                              # nc/htgf.f:218-222: GO TO 2600 — skips the HTG≤0.1⇒0.01 trap
            else
                # nc/htgf.f:226-232: SITAGE≥AGMAX ⇒ POTHTG=0.10 and GO TO 140 — the label that applies the RELHT/CR
                # XMOD modifier, exactly as the site-curve path does (jl skipped XMOD there: MEASURED FVSnc_g16
                # 7690091010901 2011 DF H 163 HtG 0.1064 live / 0.0100 jl).
                if sitage >= agmax
                    pothtg = 0.10f0
                else
                    hguess = nc_htcalc(si, sp, sitage + 5.0f0)
                    pothtg = hguess - sitht
                end
                xmod = 1.0f0
                if pccf >= 50.0f0
                    relht = avh > 0f0 ? h / avh : 1f0
                    pccf < 100.0f0 && (relht = (relht + 1.0f0) / 2.0f0)
                    relht > 1.0f0 && (relht = 1.0f0)
                    cr = Float32(t.crown_pct[i]) / 10.0f0
                    xmod = -0.02647f0 + 0.71338f0 * relht * relht + 0.06851f0 * cr
                end
                htg = pothtg * xmod
            end
        end
        (sp != 12 && h >= 300.0f0) || (htg <= 0.1f0 && (htg = 0.01f0))
        htg = scale * htg * fexp(htcon)                 # SCALE·XHT(=1)·EXP(HTCON); MISHGF=1 baseline
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = cap - h; htg < 0.1f0 && (htg = 0.1f0))   # nc/htgf.f:2600 SIZCAP, floored at 0.1
        t.ht_growth[i] = htg
    end
    return s
end
