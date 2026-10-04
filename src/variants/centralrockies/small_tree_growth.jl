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
    cw = clim_wk4(s, Float32(current_cycle_year(s)) + fint / 2f0)   # WK4 = CLGMULT (cr/regent.f:313 ·WK4(I)); nothing ⇒ 1
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
    # FVS regent.f:197-239 processes trees SPECIES-SORTED (DO ISPC=1,MAXSP; DO I3=ISCT(1),ISCT(2) via IND1 —
    # SPESRT's chain sort ⇒ record order WITHIN a species). The per-record ZZRAN (BACHLO) draws must happen in
    # THIS order, else the per-tree deviate — and EVERY downstream RNG draw — desyncs vs live on multi-species
    # stands (jl's record-order interleaving ≠ FVS's species grouping). Iterate species-then-record to match.
    # IND1 within a species is the SPESRT lineage order (post-TRIPLE: copy1, original, copy2), not storage order —
    # species_major_order, as TT/IE (a storage-order walk hands each triple member the wrong ZZRAN).
    _sp_order = species_major_order(s)
    @inbounds for i in _sp_order
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
        con = fexp(c.htg_cor_small[sp])                       # RHCON·EXP(HCOR) (regent.f:204); HCOR from the CR REGCAL calib
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
                pctred, con, 1.0f0, 1.0f0, scale, scale2, (cw === nothing ? 1f0 : cw[i]), htg_large, s.control.sp_size_cap[sp, 4],
                p.sp_site_index[sp], bark, ivf, false, zzran, dgmax[sp], brkv[sp], xminv[sp], xmaxv[sp],
                diamv[sp], ax, ht2v[sp])
            # regent.f:343-346: HK=H+HTG(K)≤4.5 ⇒ DG(K)=0 and DBH(K)=D+0.001·HK set IMMEDIATELY — MORTS (SDQ0/
            # SUMDR0 on DBH, G=DG/BARK=0) and every other pre-UPDATE consumer reads the bumped DBH with a zero
            # increment. Carrying the bump as dg=0.001·HK·BARK instead left the start-of-cycle Zeide DR0 low on
            # seedling-heavy stands (46279527020004: SUMDR0 16531.86 vs live 16533.72) ⇒ a lower self-thinning
            # target and a ~1e-3 relative MortPA drift on every record. Tripled copies get their own D+0.001·HK_L
            # through stash.dbhU/dbhL (applied at TRIPLE, as UT/CI).
            direct = small_d && (h + htg) <= 4.5f0
            dbhk = direct ? d + 0.001f0 * (h + htg) : d
            direct && (dg = 0f0)
            if l == 0
                t.ht_growth[i] = htg
                small_d && (t.diam_growth[i] = dg)
                direct && (t.dbh[i] = dbhk)
            elseif l == 1
                stash.htgU[i] = htg; stash.is_small[i] = true
                small_d && (stash.dgU[i] = dg)              # D<BKPT ⇒ regent DG; else keep the driver's gemdg dgU
                small_d && hasproperty(stash, :dbhU) && (stash.dbhU[i] = dbhk)
            else
                stash.htgL[i] = htg
                small_d && (stash.dgL[i] = dg)
                small_d && hasproperty(stash, :dbhL) && (stash.dbhL[i] = dbhk)
            end
        end
    end
    return s
end

# cr_esgent! — cr/esgent.f: SPESRT, REGENT(.TRUE.,ITRNIN) over the records this cycle's ESTAB appended (ITRNIN =
# nstart+1), then the WK4=HTIMLT tail. REGENT(LESTB) (cr/regent.f) differs from the cycle REGENT in:
# - FNT = FINT−5 (LSKIPH, HTG=0, when FINT ≤ 5); SCALE = FNT/REGYR (:155-163).
# - CCF/AVHT blend the gradd.f:192 DENSE (post-growth, PRE-regen RELDEN/AVH) with the start-of-cycle ATCCF/ATAVH:
#   CCF = (5/FINT)·RELDEN + ((FINT−5)/FINT)·ATCCF (:178-181) for PCTRED.
# - each new record (species-major IND1 order, I ≥ ITRNIN) first draws its open-grown crown on the main stream
#   CR = 0.89722 − 0.0000461·PCCF(ITRE) + 0.07985·RAN, RAN ∈ [−1,1], clamp [.20,.90] (:247-261), THEN its ZZRAN (:308-310)
#   — interleaved per record, so CR owns the crown draw (not establish!'s shared phase 2).
# - HTGR·…·WK4(I) with WK4 = HTIMLT (:313), XWT = 0 (:320), no tripling (:433).
# - D < BREAK: HK ≤ 4.5 ⇒ DG = 0, DBH = D + 0.001·HK; else DBH = DK (floored at DIAM) + 0.001·HK and DG = DBH (:343-399).
# esgent.f:49-65: HTEMP = HT+HTG; HTG·=WK4; HT+=HTG; WK4 < 1 ⇒ HT < 4.5: DBH = 0.1+0.001·HT, DG = 0, else DBH·=HT/HTEMP,
# DG = DBH·HT/HTEMP; HT capped at HHTMAX. estab.f:700 then adds GENTIM to ABIRTH.
function cr_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0,
                    atavh::Float32 = -1.0f0, atrelden::Float32 = -1.0f0,
                    relden_pre::Float32 = -1.0f0, avh_pre::Float32 = -1.0f0,
                    pccf_pre::Vector{Float32} = Float32[])
    p, t, c, sd = s.plot, s.trees, s.calib, s.coef.species
    nstart >= t.n && return s
    dgmax = sd[:st_dgmax]; xmaxv = sd[:st_xmax]; diamv = sd[:st_diam]
    htadj = sd[:st_htadj]; brkv = sd[:st_break]; ht2v = sd[:ht2]; ht1v = sd[:ht1]
    lo = sd[:site_lo]; hi = sd[:site_hi]
    aa = c.ht_dbh_aa; iabflg = c.ht_dbh_iabflg
    imodty = Int(p.model_type)
    n = t.n
    pccfv = isempty(pccf_pre) ? s.density.point_ccf : pccf_pre
    @inline _pccf(i) = (pt = Int(t.plot_id[i]); (1 <= pt <= length(pccfv)) ? pccfv[pt] : 0f0)
    cur_year = Int(current_cycle_year(s))
    # regent.f:155-176
    fnt = fint; lskiph = false
    fint <= 5f0 ? (lskiph = true) : (fnt = fnt - 5f0)
    scale = fnt / 10.0f0                                  # REGYR = 10
    relden = relden_pre >= 0f0 ? relden_pre : p.relative_density
    avh = avh_pre >= 0f0 ? avh_pre : p.avg_height
    ccf = relden; avht = avh
    if fnt > 0f0
        ccf = (5.0f0 / fint) * relden + ((fint - 5.0f0) / fint) * (atrelden >= 0f0 ? atrelden : relden)
        avht = (5.0f0 / fint) * avh + ((fint - 5.0f0) / fint) * (atavh >= 0f0 ? atavh : avh)
    end
    x = avht * (ccf / 100.0f0); x > 300.0f0 && (x = 300.0f0)
    pctred = _CR_AB[1] + x*(_CR_AB[2] + x*(_CR_AB[3] + x*(_CR_AB[4] + x*(_CR_AB[5] + x*_CR_AB[6]))))
    pctred > 1.0f0 && (pctred = 1.0f0); pctred < 0.01f0 && (pctred = 0.01f0)
    dgsd = s.control.dg_sd
    @inbounds for i in (nstart + 1):n                      # estab.f: a new record starts with DG = HTG = 0
        t.diam_growth[i] = 0f0; t.ht_growth[i] = 0f0
    end
    order = sort(collect((nstart + 1):n); by = i -> (Int(t.species[i]), i))   # SPESRT (esgent.f:44), I ≥ ITRNIN
    @inbounds for i in order
        sp = Int(t.species[i]); d = t.dbh[i]; h = t.height[i]
        # regent.f:247-261 crown of the new record (before the D ≥ XMAX skip)
        cr0 = 0.89722f0 - 0.0000461f0 * _pccf(i)
        ran = 0f0
        while true; ran = bachlo(s.rng, 0f0, 1f0); (-1f0 <= ran <= 1f0) && break; end
        cr0 = cr0 + 0.07985f0 * ran
        cr0 > 0.90f0 && (cr0 = 0.90f0); cr0 < 0.20f0 && (cr0 = 0.20f0)
        icr = trunc(Int32, cr0 * 100f0 + 0.5f0)
        t.crown_pct[i] = icr; t.crown_ratio[i] = Float32(icr)
        d >= xmaxv[sp] && continue                         # regent.f:265
        bark = cr_bratio(sd, sp, d, imodty)
        local htg::Float32
        if lskiph
            htg = 0f0
        else
            xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
            con = fexp(c.htg_cor_small[sp])                # RHCON·EXP(HCOR)
            si = p.sp_site_index[sp]
            si > hi[sp] && (si = hi[sp]); si <= lo[sp] && (si = lo[sp] + 0.5f0)
            relsi = (si - lo[sp]) / (hi[sp] - lo[sp]); rsimod = 0.5f0 * (1.0f0 + relsi)
            local htgr::Float32
            if sp == 20 || sp == 28                        # aspen / paper birch (Sheppard), SITAGE = ABIRTH
                ag1 = t.birth_age[i]; ag1 < 5f0 && (ag1 = 5f0)
                hite1 = 26.9825f0 * fpow(ag1, 1.1752f0)
                hite2 = 26.9825f0 * fpow(ag1 + 10.0f0, 1.1752f0)
                htgr = (hite2 - hite1) / (2.54f0 * 12.0f0) * rsimod * con
                htgr = htgr * 0.75f0
            else
                xv = Float32(icr) / 100f0
                vigor = 150.0f0 * fpow(xv, 3.0f0) * fexp(-6.0f0 * xv) + 0.3f0
                vigor > 1.0f0 && (vigor = 1.0f0)
                (sp in _CR_IVFLAG) && (vigor = 1.0f0 - (1.0f0 - vigor) / 3.0f0)
                pothtg = p.sp_site_index[sp] / (15.0f0 - 4.0f0 * relsi) * htadj[sp]
                htgr = pothtg * pctred * vigor * con
            end
            zzran = 0f0
            if dgsd >= 1.0f0
                while true
                    zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                    (zzran <= 0.5f0 && zzran >= -2.0f0) && break
                end
            end
            htgr = (htgr + zzran * 0.2f0) * xrhgro * scale * t.htimlt[i]
            htg = htgr * (1.0f0 - 0f0) + 0f0 * t.ht_growth[i]    # XWT = 0 under LESTB
            htg < 0.1f0 && (htg = 0.1f0)
            sizcap4 = s.control.sp_size_cap[sp, 4]
            if h + htg > sizcap4
                htg = sizcap4 - h; htg < 0.1f0 && (htg = 0.1f0)
            end
        end
        t.ht_growth[i] = htg
        d >= brkv[sp] && continue                          # regent.f:342 D ≥ BKPT ⇒ no diameter assignment
        hk = h + htg
        local dbh::Float32, dg::Float32
        if hk <= 4.5f0
            dg = 0f0; dbh = d + 0.001f0 * hk
        else
            local dk::Float32
            if sp == 13 || sp == 36                        # ponderosa / Chihuahua pine
                dk = (hk - 8.31485f0 + 0.59200f0 * 7.0f0) / 3.03659f0
            elseif sp in _CR_IVFLAG                        # pinyon/juniper/oak/bristlecone (IVFLAG = 1)
                dk = (hk - 4.5f0) * 10.0f0 / (p.sp_site_index[sp] - 4.5f0)
            else
                ax = iabflg[sp] == 0 ? aa[sp] : ht1v[sp]
                dk = (ht2v[sp] / (flog(hk - 4.5f0) - ax)) - 1.0f0
            end
            dk < 0.1f0 && (dk = 0.1f0)
            dbh = dk; dbh < diamv[sp] && (dbh = diamv[sp])
            dbh = dbh + 0.001f0 * hk
            dg = dbh
            (dbh + dg) < diamv[sp] && (dg = diamv[sp] - dbh)
        end
        dg = dg_bound(nothing, nothing, sp, dbh, dg, s.control.sp_size_cap)   # DGBND
        t.dbh[i] = dbh; t.diam_growth[i] = dg
    end
    # esgent.f:49-65
    @inbounds for i in (nstart + 1):n
        sp = Int(t.species[i])
        htemp = t.height[i] + t.ht_growth[i]
        t.ht_growth[i] = t.ht_growth[i] * t.htimlt[i]
        t.height[i] = t.height[i] + t.ht_growth[i]
        if t.htimlt[i] < 1f0
            if t.height[i] < 4.5f0
                t.dbh[i] = 0.1f0 + 0.001f0 * t.height[i]; t.diam_growth[i] = 0f0
            else
                t.dbh[i] = t.dbh[i] * (t.height[i] / htemp)
                t.diam_growth[i] = t.dbh[i] * (t.height[i] / htemp)
            end
        end
        t.height[i] > _CR_ES_HHTMAX[sp] && (t.height[i] = _CR_ES_HHTMAX[sp])
    end
    esgent_add_gentim!(s, nstart, fint)                    # estab.f:700 ABIRTH += GENTIM
    return s
end
