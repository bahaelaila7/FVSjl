# =============================================================================
# small_tree_growth.jl (centralstates) — CS small-tree (regen) growth (cs/regent.f)
#
# Same SHAPE as NE's REGENT (height-driven model overriding the large-tree DGF/HTGF,
# blended over [XMIN, XMAX]), with the CS specifics from the cs/regent.f vs ne/regent.f
# diff:
#   * XMIN = 3.0 (NE 1.5) — the lower end of the blend range.
#   * the BAL modifier is CS's `cs_balmod` (BAL/BA form), not NE's exp(−b3·BAL).
#   * the HTCALC curve uses MAPCS (IVAR=2) via cs_htcalc_* (shared LTBHEC).
#   * the budwidth DIAM floor = SNDBAL = the htdbh_db column (cs/regent.f:340/384).
# REGYR=10, XMAX=5, DGMAX=5 (same as NE).
# =============================================================================

const CS_REGENT_XMIN = 3.0f0
const CS_REGENT_XMAX = 5f0
const CS_REGENT_REGYR = 10f0
const CS_REGENT_YR = 10f0
const CS_REGENT_DGMAX = 5f0

function small_tree_growth!(s::StandState, stash, ::CentralStates; fint::Float32 = 10f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    regent_diam = sd[:htdbh_db]                     # DIAM budwidth floor (= SNDBAL)
    b1 = sd[:balmod_b1]; b2 = sd[:balmod_b2]; b3 = sd[:balmod_b3]
    sizcap = s.control.sp_size_cap
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    ba = p.basal_area
    avh = p.avg_height
    species_sort!(s)
    trip = stash !== nothing
    nrec = trip ? 3 : 1
    scale  = fint / CS_REGENT_REGYR                 # FNT/REGYR (=1 at FINT=10)
    scale2 = CS_REGENT_YR / fint                    # YR/FNT
    random_on = s.control.dg_stddev_bound >= 1f0
    cur_year = current_cycle_year(s)
    @inbounds for sp in 1:nspecies(s.variant)
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        si = p.sp_site_index[sp]
        htmax = cs_htcalc_htmax(sp, si)
        dgmx = CS_REGENT_DGMAX * scale
        con = exp(c.htg_cor_small[sp])
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        for k3 in i1:i2
            i = ind1[k3]
            d = t.dbh[i]
            (d >= CS_REGENT_XMAX || t.tpa[i] <= 0f0) && continue
            h = t.height[i]
            if htmax - h <= 1f0
                htgr_s = 0.1f0
            else
                aget = cs_htcalc_age(sp, si, h)
                htg1 = cs_htcalc_incr(sp, si, aget)
                htgr_s = htg1 * con * scale * xrhgro
            end
            # BAL + relative-height modifier (cs/regent.f: NO ·0.8, unlike the large-tree HTGF)
            bal = (1f0 - t.crown_ratio[i] / 100f0) * ba
            gmod = cs_balmod(b1[sp], b2[sp], b3[sp], bal, ba, d)
            relht = avh > 0f0 ? min(h / avh, 1f0) : 0f0
            gmod = 1f0 - (1f0 - gmod) * (1f0 - relht)
            htgr_s = max(htgr_s * gmod, 0.1f0)
            xwt = d <= CS_REGENT_XMIN ? 0f0 : (d - CS_REGENT_XMIN) / (CS_REGENT_XMAX - CS_REGENT_XMIN)
            htgr = max(htgr_s * (1f0 - xwt) + xwt * t.ht_growth[i], 0.1f0)
            dgk = t.diam_growth[i]
            dgkU = trip ? stash.dgU[i] : dgk
            dgkL = trip ? stash.dgL[i] : dgk
            for l in 0:(nrec - 1)
                dgk_l = l == 0 ? dgk : (l == 1 ? dgkU : dgkL)
                ran = 0f0
                if random_on
                    while true
                        ran = bachlo(s.rng, 0f0, 1f0)
                        (-1f0 <= ran <= 1f0) && break
                    end
                end
                htg = max(htgr + ran * 0.1f0 * htgr, 0.1f0)
                (h + htg) > sizcap[sp, 4] && (htg = max(sizcap[sp, 4] - h, 0.1f0))
                hk = h + htg
                if hk <= 4.5f0
                    dg = 0.001f0 * hk
                else
                    bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
                    dkk = _htdbh_dbh(sd, sp, hk, Int(p.forest_idx))
                    dk  = h <= 4.5f0 ? d : _htdbh_dbh(sd, sp, h, Int(p.forest_idx))
                    if dk < 0f0 || dkk < 0f0
                        dg = htg * 0.2f0 * bark * xrdgro
                    else
                        dgsm = (dkk - dk) * bark * xrdgro
                        dgsm < 0f0 && (dgsm = 0f0)
                        dds = dgsm * (2f0 * bark * d + dgsm) * scale2
                        dgsm = sqrt((d * bark)^2 + dds) - bark * d
                        dgsm < 0f0 && (dgsm = 0f0)
                        dggr = dgsm * (1f0 - xwt) + xwt * dgk_l
                        dg = max(dggr, 0.1f0)
                        dg > dgmx && (dg = dgmx)
                        (d + dg) < regent_diam[sp] && (dg = regent_diam[sp] - d)
                    end
                end
                dg = dg_bound(nothing, nothing, sp, d, dg, sizcap)
                if fint != CS_REGENT_YR && dg > 0f0
                    bk = bark_ratio(c.bark_a, c.bark_b, sp, d); dib = d * bk
                    dds_e = dg * (2f0 * dib + dg) * (fint / CS_REGENT_YR)
                    dg = sqrt(dib * dib + dds_e) - dib
                end
                if l == 0
                    t.diam_growth[i] = dg; t.ht_growth[i] = htg
                elseif l == 1
                    stash.dgU[i] = dg; stash.htgU[i] = htg; stash.is_small[i] = true
                else
                    stash.dgL[i] = dg; stash.htgL[i] = htg
                end
            end
        end
    end
    return s
end
