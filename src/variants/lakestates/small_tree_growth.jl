# =============================================================================
# small_tree_growth.jl (lakestates) — LS small-tree (regen) growth (ls/regent.f)
#
# Same SHAPE as NE/CS REGENT (height-driven model overriding the large-tree DGF/HTGF,
# blended over [XMIN, XMAX]). ls/regent.f == cs/regent.f EXCEPT the two variant hooks:
#   * the HTCALC curve uses MAPLS (IVAR=1) — the shared LTBHEC via `_ls_htcoef`;
#   * the competition modifier is `ls_balmod(sp, D, BA, RMSQD)` (ls/balmod.f:109), not cs_balmod.
# XMIN=3.0, XMAX=5.0, DGMAX=5.0, REGYR=10, YR=10 (all identical to CS/NE); DIAM budwidth
# floor = SNDBAL (htdbh_db). See centralstates/small_tree_growth.jl for the shared logic notes.
# =============================================================================
const LS_REGENT_XMIN = 3.0f0
const LS_REGENT_XMAX = 5f0
const LS_REGENT_REGYR = 10f0
const LS_REGENT_YR = 10f0
const LS_REGENT_DGMAX = 5f0
function small_tree_growth!(s::StandState, stash, ::LakeStates; fint::Float32 = 10f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    regent_diam = sd[:htdbh_db]                     # DIAM budwidth floor (= SNDBAL)
    mapls = sd[:htcalc_mapls]
    check = sd[:balmod_check]; b1 = sd[:balmod_b1]; b2 = sd[:balmod_b2]; b3 = sd[:balmod_b3]
    b4 = sd[:balmod_b4]; c1 = sd[:balmod_c1]; c2 = sd[:balmod_c2]; bamax1 = sd[:balmod_bamax1]
    sizcap = s.control.sp_size_cap
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    ba = p.basal_area
    avh = p.avg_height; rmsqd = stand_qmd(s)        # RMSQD (dense.f:250) — p.qmd is never stored
    species_sort!(s)
    trip = stash !== nothing
    nrec = trip ? 3 : 1
    scale  = fint / LS_REGENT_REGYR                 # FNT/REGYR (=1 at FINT=10)
    scale2 = LS_REGENT_YR / fint                    # YR/FNT
    random_on = s.control.dg_stddev_bound >= 1f0
    cur_year = current_cycle_year(s)
    @inbounds for sp in 1:nspecies(s.variant)
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        si = p.sp_site_index[sp]
        coef = _ls_htcoef(mapls, sp)
        htmax = _ls_htmax(coef, si)
        dgmx = LS_REGENT_DGMAX * scale
        con = exp(c.htg_cor_small[sp])
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        for k3 in i1:i2
            i = ind1[k3]
            d = t.dbh[i]
            (d >= LS_REGENT_XMAX || t.tpa[i] <= 0f0) && continue
            h = t.height[i]
            if htmax - h <= 1f0
                htgr_s = 0.1f0
            else
                aget = _ls_age(coef, si, h)
                htg1 = _ls_incr(coef, si, aget)
                htgr_s = htg1 * con * scale * xrhgro
            end
            # LS balmod + relative-height modifier (ls/regent.f: NO ·0.8, unlike the large-tree HTGF)
            gmod = ls_balmod(sp, d, ba, rmsqd, check, b1, b2, b3, b4, c1, c2, bamax1)
            relht = avh > 0f0 ? min(h / avh, 1f0) : 0f0
            gmod = 1f0 - (1f0 - gmod) * (1f0 - relht)
            htgr_s = max(htgr_s * gmod, 0.1f0)
            xwt = d <= LS_REGENT_XMIN ? 0f0 : (d - LS_REGENT_XMIN) / (LS_REGENT_XMAX - LS_REGENT_XMIN)
            htgr = max(htgr_s * (1f0 - xwt) + xwt * t.ht_growth[i], 0.1f0)
            dgk = t.diam_growth[i]
            dgkU = trip ? stash.dgU[i] : dgk
            dgkL = trip ? stash.dgL[i] : dgk
            if fint != LS_REGENT_YR
                bk0 = bark_ratio(c.bark_a, c.bark_b, sp, d); dib0 = d * bk0
                un10(g) = g > 0f0 ? sqrt(dib0 * dib0 + g * (2f0 * dib0 + g) * scale2) - dib0 : g
                dgk = un10(dgk); dgkU = un10(dgkU); dgkL = un10(dgkL)
            end
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
                    dkk = _htdbh_dbh(sd, sp, hk, Int(p.forest_idx); db_floor = true)
                    dk  = h <= 4.5f0 ? d : _htdbh_dbh(sd, sp, h, Int(p.forest_idx); db_floor = true)
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
                if fint != LS_REGENT_YR && dg > 0f0
                    bk = bark_ratio(c.bark_a, c.bark_b, sp, d); dib = d * bk
                    dds_e = dg * (2f0 * dib + dg) * (fint / LS_REGENT_YR)
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
