# =============================================================================
# small_tree_growth.jl (northeast) — NE small-tree (regen) growth (ne/regent.f)
#
# Trees below XMAX = 5" DBH grow by a HEIGHT-driven model that OVERRIDES the
# large-tree DGF/HTGF, blended with the large-tree prediction over [XMIN=1.5, 5].
# Runs after height_growth!+diameter_growth!, before mortality (grincr.f order).
#
# Height: NC-128 increment (ne_htcalc_incr, 10-yr) × CON × SCALE × HGADJ × XRHGRO,
#   then the SAME BAL+relative-height modifier as the large-tree HTGF but WITHOUT
#   the ·0.8 (GMOD = 1−(1−bal)(1−relht)), then blended with the large-tree HTG by
#   xwt, ±10% random (BACHLO when DGSD≥1), size-capped (regent.f:184-281).
# Diameter: Wykoff HT→DBH inverse at the new vs old height (DKK−DK), DDS-scaled,
#   blended with the large-tree DG, DIAM-floored, DGBND-capped (regent.f:284-390).
#
# NE constants (regent.f:96-99): REGYR=10, XMIN=1.5, XMAX=5, DGMAX=5, HGADJ=1.
# CON = RHCON·exp(HCOR); RHCON=1 (no RCOR2 keyword) and the small-tree height
# calibration HCOR is 0 for net01 (the .tre carries measured DG, not measured HTG,
# so the regent.f:431-558 ratio estimator never accumulates) ⇒ CON=1.
# =============================================================================

const NE_REGENT_XMIN = 1.5f0
const NE_REGENT_XMAX = 5f0
const NE_REGENT_REGYR = 10f0
const NE_REGENT_YR = 10f0
const NE_REGENT_DGMAX = 5f0

function small_tree_growth!(s::StandState, stash, ::Northeast; fint::Float32 = 10f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    regent_diam = sd[:regent_min_diam]
    b3_dg = sd[:dg_b3]                 # BALMOD uses the DG b3 (same as height_growth!)
    sizcap = s.control.sp_size_cap
    isct = s.control.sp_count_tab; ind1 = s.scratch.idx1
    avh = p.avg_height
    ebau = zeros(Float32, 50); ne_badist!(ebau, s)   # cycle-start BAL basis
    species_sort!(s)
    trip = stash !== nothing
    nrec = trip ? 3 : 1
    scale  = fint / NE_REGENT_REGYR    # FNT/REGYR (=1 at FINT=10)
    scale2 = NE_REGENT_YR / fint       # YR/FNT — convert the FNT-yr small-tree DG to the 10-yr basis
    random_on = s.control.dg_stddev_bound >= 1f0   # DGSD (default 2) ⇒ ±10% random
    cur_year = current_cycle_year(s)
    @inbounds for sp in 1:nspecies(s.variant)
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        si = p.sp_site_index[sp]
        htmax = ne_htcalc_htmax(sp, si)
        dgmx = NE_REGENT_DGMAX * scale
        con = exp(c.htg_cor_small[sp])     # RHCON·exp(HCOR); RHCON=1, htg_cor_small=HCOR (0 for net01)
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        for k3 in i1:i2
            i = ind1[k3]
            d = t.dbh[i]
            (d >= NE_REGENT_XMAX || t.tpa[i] <= 0f0) && continue
            h = t.height[i]
            # small-tree height increment (HGADJ=1)
            if htmax - h <= 1f0
                htgr_s = 0.1f0
            else
                aget = ne_htcalc_age(sp, si, h)
                htg1 = ne_htcalc_incr(sp, si, aget)
                htgr_s = htg1 * con * scale * xrhgro
            end
            # BAL + relative-height modifier (regent.f:232-238; NO ·0.8, unlike large-tree HTGF)
            gmod = ne_balmod(b3_dg[sp], ebau, d)
            relht = avh > 0f0 ? min(h / avh, 1f0) : 0f0
            gmod = 1f0 - (1f0 - gmod) * (1f0 - relht)
            htgr_s = max(htgr_s * gmod, 0.1f0)
            xwt = d <= NE_REGENT_XMIN ? 0f0 : (d - NE_REGENT_XMIN) / (NE_REGENT_XMAX - NE_REGENT_XMIN)
            htgr = max(htgr_s * (1f0 - xwt) + xwt * t.ht_growth[i], 0.1f0)
            dgk = t.diam_growth[i]   # large-tree DG (inside-bark) to blend toward
            # regent.f:373 blends the small-tree DG with the PER-TRIPLE large-tree DG: central uses DG(I),
            # upper uses DG(ITRIPU)=stash.dgU, lower uses DG(ITRIPL)=stash.dgL (the FRU/FRL serial-corr
            # records DGDRIV set). Capture them before the l-loop overwrites the stash slots.
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
                # diameter increment: Wykoff HT→DBH inverse (DKK−DK), DDS-scaled, blended
                hk = h + htg
                if hk <= 4.5f0
                    dg = 0.001f0 * hk
                else
                    bark = bark_ratio(c.bark_a, c.bark_b, sp, d)
                    dkk = _htdbh_dbh(sd, sp, hk, Int(p.forest_idx))
                    dk  = h <= 4.5f0 ? d : _htdbh_dbh(sd, sp, h, Int(p.forest_idx))
                    if dk < 0f0 || dkk < 0f0
                        dg = htg * 0.2f0 * bark * xrdgro      # regent.f:359 degenerate fallback
                        dgsm = dg
                    else
                        dgsm = (dkk - dk) * bark * xrdgro
                        dgsm < 0f0 && (dgsm = 0f0)
                        dds = dgsm * (2f0 * bark * d + dgsm) * scale2
                        dgsm = sqrt((d * bark)^2 + dds) - bark * d
                        dgsm < 0f0 && (dgsm = 0f0)
                    end
                    dggr = dgsm * (1f0 - xwt) + xwt * dgk_l   # blend with PER-TRIPLE large-tree DG (regent.f:373)
                    dg = max(dggr, 0.1f0)
                    dg > dgmx && (dg = dgmx)
                    (d + dg) < regent_diam[sp] && (dg = regent_diam[sp] - d)   # DIAM budwidth floor
                end
                dg = dg_bound(nothing, nothing, sp, d, dg, sizcap)   # DGBND = SIZCAP-only for NE
                # GRADD re-expand the 10-yr basis to the FINT cycle (identity at FINT=10=YR)
                if fint != NE_REGENT_YR && dg > 0f0
                    bk = bark_ratio(c.bark_a, c.bark_b, sp, d); dib = d * bk
                    dds_e = dg * (2f0 * dib + dg) * (fint / NE_REGENT_YR)
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
