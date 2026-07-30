# =============================================================================
# small_tree_growth.jl (centralrockies) — CR small-tree/regen growth (cr/regent.f)
#
# For trees with DBH < XMAX[sp]: a height increment from a site/density/vigor model
# (POTHTG·PCTRED·VIGOR·CON, or the Sheppard aspen/birch curve), weighted (XWT) with
# the large-tree htgf estimate; then a diameter increment from the height change via
# a species-specific inverse height–DBH curve, DDS-scaled to the cycle.
#
# Math: ALOG->flog, EXP->fexp, **realexp->fpow, SQRT->Base. Stochastic ZZRAN
# (DGSD≥1, ASYMMETRIC reject window [-2, 0.5]) injected by the caller.
# =============================================================================

"""
    _cr_regent_tree(...) -> (htg, dg)

Per-tree small-tree height + diameter increment (cr/regent.f growth path). `htg_large` = the htgf
large-tree height increment for the XWT blend (0 for pure seedlings). `pothtg`/`pctred`/`rsimod` are
stand/species terms from the caller. `zzran` is the injected stochastic deviate (0 ⇒ deterministic).
"""
function _cr_regent_tree(sp::Int, d::Float32, h::Float32, icr::Int, abirth::Float32,
                         rsimod::Float32, pothtg::Float32, pctred::Float32, con::Float32,
                         xrhgro::Float32, xrdgro::Float32, scale::Float32, scale2::Float32,
                         wk4::Float32, htg_large::Float32, sizcap4::Float32, sitear::Float32,
                         bark::Float32, ivflag::Bool, lskiph::Bool, zzran::Float32,
                         dgmax::Float32, break_sp::Float32, xmn::Float32, xmx::Float32,
                         diam_sp::Float32, ax::Float32, ht2::Float32)
    dgmx = dgmax * scale
    # ---- HEIGHT increment ----
    local htg::Float32
    if lskiph
        htg = 0.0f0
    else
        if sp == 20 || sp == 28                            # aspen / paper birch (Sheppard curve)
            ag1 = abirth < 5.0f0 ? 5.0f0 : abirth
            hite1 = 26.9825f0 * fpow(ag1, 1.1752f0)
            ag2 = ag1 + 10.0f0
            hite2 = 26.9825f0 * fpow(ag2, 1.1752f0)
            htgr = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con * 0.75f0
        else
            x = Float32(icr) / 100.0f0
            vigor = 150.0f0 * fpow(x, 3.0f0) * fexp(-6.0f0 * x) + 0.3f0
            vigor > 1.0f0 && (vigor = 1.0f0)
            ivflag && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)
            htgr = pothtg * pctred * vigor * con
        end
        htgr = (htgr + zzran * 0.2f0) * xrhgro * scale * wk4
        xwt = (d - xmn) / (xmx - xmn)
        d <= xmn && (xwt = 0.0f0)
        htg = htgr * (1.0f0 - xwt) + xwt * htg_large
        htg < 0.1f0 && (htg = 0.1f0)
        if h + htg > sizcap4
            htg = sizcap4 - h; htg < 0.1f0 && (htg = 0.1f0)
        end
    end
    # ---- DIAMETER increment (only for D < BREAK[sp]) ----
    dg = 0.0f0
    if d < break_sp
        hk = h + htg
        if hk <= 4.5f0
            # regent.f:344-346: sub-breast-height (HK≤4.5) ⇒ DG(K)=0 but DBH(K)=D+0.001·HK (a tiny
            # height-tied diameter bump). jl applies dbh += dg/bark with the SAME cr_bratio(sp,D,imodty),
            # so dg = 0.001·HK·bark lands dbh = D+0.001·HK exactly. WITHOUT this the seedling DBH is pinned
            # at inventory (e.g. Gambel-oak regen stuck at D=0.1) — it never crosses the ccfcal CCF cliff
            # (D>0.1 ⇒ RDA·D^RDB vs D≤0.1 ⇒ 0.001), collapsing stand CCF and driving the dense-regen
            # structure_densephase divergence.
            dg = 0.001f0 * hk * bark
        else
            local dk::Float32, dkk::Float32
            if sp == 13 || sp == 36                        # ponderosa / Chihuahua pine
                dk = (hk - 8.31485f0 + 0.59200f0 * 7.0f0) / 3.03659f0; dk < 0.1f0 && (dk = 0.1f0)
                dkk = (h - 8.31485f0 + 0.59200f0 * 7.0f0) / 3.03659f0; dkk < 0.1f0 && (dkk = 0.1f0)
                h < 4.5f0 && (dkk = d)
            elseif ivflag                                  # pinyon/juniper/oak/bristlecone
                dk = (hk - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                dkk = (h - 4.5f0) * 10.0f0 / (sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
                h < 4.5f0 && (dkk = d)
            else                                           # all other species: AX = cratet-fitted AA or HT1
                bx = ht2
                dk = (bx / (flog(hk - 4.5f0) - ax)) - 1.0f0; dk < 0.1f0 && (dk = 0.1f0)
                dkk = h <= 4.5f0 ? d : (bx / (flog(h - 4.5f0) - ax)) - 1.0f0
            end
            if dk < 0.0f0 || dkk < 0.0f0
                dg = htg * 0.2f0 * bark * xrdgro
                dk = d + dg
            else
                dg = (dk - dkk) * bark * xrdgro
            end
            dg < 0.0f0 && (dg = 0.0f0)
            dg > dgmx && (dg = dgmx)
            dds = dg * (2.0f0 * bark * d + dg) * scale2
            dg = sqrt(fpow(d * bark, 2.0f0) + dds) - bark * d
            # DIAM floor is INSIDE the HK>4.5 branch (regent.f:421-423); DBH(K)≈D for cycling.
            (d + dg) < diam_sp && (dg = diam_sp - d)
        end
    end
    return htg, dg
end

# AB density-modifier polynomial (regent.f DATA AB, 6 terms used).
const _CR_AB = (1.11436f0, -0.011493f0, 0.43012f-4, -0.72221f-7, 0.5607f-10, -0.1641f-13)
# IVFLAG species (vigor cut + pinyon/juniper/oak diameter path): regent.f SELECT CASE.
const _CR_IVFLAG = (9, 12, 16, 23, 24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35)

"""
    small_tree_growth!(s, stash, ::CentralRockies; fint=10)

CR small-tree height + diameter increment (cr/regent.f growth path), overriding htgf's estimate for
trees with DBH < XMAX[sp]. Reads the cratet-fitted HT-DBH intercept (`s.calib.ht_dbh_aa`/`ht_dbh_iabflg`)
for the diameter-from-height inverse. Draws the ZZRAN stochastic deviate (BACHLO, asymmetric [-2,0.5]
reject) per tree — its per-tree VALUES bit-match live only once the cycle RNG sequence is reconciled (ch9).
"""
function small_tree_growth!(s::StandState, stash, ::CentralRockies; fint::Float32 = 10.0f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    t.n == 0 && return s
    dgmax = sd[:st_dgmax]; xmaxv = sd[:st_xmax]; xminv = sd[:st_xmin]; diamv = sd[:st_diam]
    htadj = sd[:st_htadj]; brkv = sd[:st_break]; ht2v = sd[:ht2]; ht1v = sd[:ht1]
    lo = sd[:site_lo]; hi = sd[:site_hi]
    aa = c.ht_dbh_aa; iabflg = c.ht_dbh_iabflg
    imodty = Int(p.model_type)
    regyr = 10.0f0
    fnt = fint
    scale = fnt / regyr
    scale2 = s.control.year / fnt                                   # YR / FNT (p.year = CR YR = 10)
    dgsd = s.control.dg_sd
    # density modifier PCTRED from AVHT * CCF (regent.f:185-190)
    ccf = stand_ccf(s); avht = p.avg_height
    x = avht * (ccf / 100.0f0); x > 300.0f0 && (x = 300.0f0)
    pctred = _CR_AB[1] + x*(_CR_AB[2] + x*(_CR_AB[3] + x*(_CR_AB[4] + x*(_CR_AB[5] + x*_CR_AB[6]))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)

    # REGENT is evaluated for EACH tripled record (regent.f:433-436 L-loop): a fresh ZZRAN per record,
    # central→trees, upper/lower→the tripling stash (dgU/dgL/htgU/htgL/is_small, filled here after the
    # driver's LARGE-tree pass). WITHOUT this the tripled small-tree records inherit the large-tree gemdg
    # DG/HTG — and CR's gemdg is explosive on tiny DBH (limber pine sp10 balloons D 1.3→13). nrec=3 while
    # tripling, else 1 (non-tripling stands draw once, identical to before).
    nrec = stash !== nothing ? 3 : 1
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0.0f0 && continue
        sp = Int(t.species[i])
        d = t.dbh[i]
        d >= xmaxv[sp] && continue                          # regent.f:271 D≥XMAX → skip
        h = t.height[i]
        # per-species site terms
        si = p.sp_site_index[sp]
        si > hi[sp] && (si = hi[sp])
        si <= lo[sp] && (si = lo[sp] + 0.5f0)
        relsi = (si - lo[sp]) / (hi[sp] - lo[sp])
        rsimod = 0.5f0 * (1.0f0 + relsi)
        pothtg = p.sp_site_index[sp] / (15.0f0 - 4.0f0 * relsi) * htadj[sp]
        con = exp(c.htg_cor_small[sp])                       # RHCON·EXP(HCOR) (regent.f:204); HCOR from the CR REGCAL calib
        ivf = sp in _CR_IVFLAG
        ax = iabflg[sp] == 0 ? aa[sp] : ht1v[sp]
        bark = cr_bratio(sd, sp, d, imodty)
        htg_large = t.ht_growth[i]     # large-tree HTG (height_growth!) for the XWT blend — read before l=0 overwrites
        small_d = d < brkv[sp]
        for l in 0:(nrec - 1)
            # ZZRAN: BACHLO(0,1) if DGSD≥1, reject if >0.5 or <-2.0 (regent.f:308-310) — per tripled record
            zzran = 0.0f0
            if dgsd >= 1.0f0
                while true
                    zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                    (zzran <= 0.5f0 && zzran >= -2.0f0) && break
                end
            end
            htg, dg = _cr_regent_tree(sp, d, h, Int(t.crown_pct[i]), t.birth_age[i], rsimod, pothtg,
                pctred, con, 1.0f0, 1.0f0, scale, scale2, 1.0f0, htg_large, s.control.sp_size_cap[sp, 4],
                p.sp_site_index[sp], bark, ivf, false, zzran, dgmax[sp], brkv[sp], xminv[sp], xmaxv[sp],
                diamv[sp], ax, ht2v[sp])
            if l == 0
                t.ht_growth[i] = htg
                small_d && (t.diam_growth[i] = dg)
            elseif l == 1
                stash.htgU[i] = htg; stash.is_small[i] = true
                small_d && (stash.dgU[i] = dg)              # D<BKPT ⇒ regent DG; else keep the driver's gemdg dgU
            else
                stash.htgL[i] = htg
                small_d && (stash.dgL[i] = dg)
            end
        end
    end
    return s
end

# cr_esgent! (esgent.f) — CR-specific: grow the JUST-ESTABLISHED regen records IN their creation cycle via
# REGENT(LESTB=T), then apply HT+=HTG·WK4, derive DBH if WK4<1, cap at HHTMAX. Eastern variants leave birth-cycle
# regen ungrown (GRADD order, bit-exact) — this is CR-only. `nstart` = tree count BEFORE establish! (the new
# records are nstart+1..t.n). Applies the increment directly to t.height/t.dbh (not via the wk2 stash).
function cr_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0)
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    nstart >= t.n && return s
    dgmax = sd[:st_dgmax]; xmaxv = sd[:st_xmax]; xminv = sd[:st_xmin]; diamv = sd[:st_diam]
    htadj = sd[:st_htadj]; brkv = sd[:st_break]; ht2v = sd[:ht2]; ht1v = sd[:ht1]
    lo = sd[:site_lo]; hi = sd[:site_hi]
    aa = c.ht_dbh_aa; iabflg = c.ht_dbh_iabflg
    imodty = Int(p.model_type)
    # Birth-cycle regen grows only FINT−GENTIM years (established mid-cycle at GENTIM=FINT−5, estab.f:448), not the
    # full cycle — so the regent SCALE uses the PARTIAL period (else the birth-cycle HTG ~2× over-shoots).
    gentim = max(fint - 5.0f0, 0.0f0)
    scale = (fint - gentim) / 10.0f0; scale2 = s.control.year / fint; dgsd = s.control.dg_sd
    ccf = stand_ccf(s); avht = p.avg_height
    x = avht * (ccf / 100.0f0); x > 300.0f0 && (x = 300.0f0)
    pctred = _CR_AB[1] + x*(_CR_AB[2] + x*(_CR_AB[3] + x*(_CR_AB[4] + x*(_CR_AB[5] + x*_CR_AB[6]))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    @inbounds for i in (nstart+1):t.n
        t.tpa[i] <= 0.0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= xmaxv[sp] && continue
        h = t.height[i]
        si = p.sp_site_index[sp]; si > hi[sp] && (si = hi[sp]); si <= lo[sp] && (si = lo[sp] + 0.5f0)
        relsi = (si - lo[sp]) / (hi[sp] - lo[sp]); rsimod = 0.5f0 * (1.0f0 + relsi)
        pothtg = p.sp_site_index[sp] / (15.0f0 - 4.0f0 * relsi) * htadj[sp]
        con = exp(c.htg_cor_small[sp]); ivf = sp in _CR_IVFLAG
        ax = iabflg[sp] == 0 ? aa[sp] : ht1v[sp]; bark = cr_bratio(sd, sp, d, imodty)
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 0.5f0 && zzran >= -2.0f0) && break
            end
        end
        htg, _ = _cr_regent_tree(sp, d, h, Int(t.crown_pct[i]), t.birth_age[i], rsimod, pothtg,
            pctred, con, 1.0f0, 1.0f0, scale, scale2, 1.0f0, t.ht_growth[i], s.control.sp_size_cap[sp, 4],
            p.sp_site_index[sp], bark, ivf, false, zzran, dgmax[sp], brkv[sp], xminv[sp], xmaxv[sp],
            diamv[sp], ax, ht2v[sp])
        # esgent.f: HTG*=WK4 (=1, no FIXHTG here); HT+=HTG; WK4<1 ⇒ DBH-derive (skipped at WK4=1); cap HHTMAX.
        nh = h + htg
        nh > _CR_ES_HHTMAX[sp] && (nh = _CR_ES_HHTMAX[sp])
        t.height[i] = nh
        t.ht_growth[i] = htg
    end
    return s
end
