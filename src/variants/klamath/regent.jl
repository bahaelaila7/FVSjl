# =============================================================================
# regent.jl (klamath) — NC small-tree growth (nc/regent.f + htgr5.f + htdbh.f). Chunk 6.
#   nc_htgr5(sp,ssite,baa,relht,cr,h) — small-tree height increment (3 methods, htgr5.f)
#   nc_htdbh(sp,d,h,mode)             — Curtis-Arney HT-DBH (SISKIY, htdbh.f; mode1 H→D)
#   small_tree_growth!(::Klamath)     — HTGR5→XWT-blend w/ large-tree HTG + small-tree DG (D<DGMIN), blended
# For NC LHTDRG=.FALSE. ⇒ regent recomputes the small-tree DK/DKK via HTDBH (SISKIY), NOT the HT1/HT2 Wykoff.
# =============================================================================

const NC_REGYR = 5.0f0
const NC_ST_XMAX = Float32[5,5,5,5,5,5,5,5,5,5,5,10]
const NC_ST_XMIN = Float32[2,2,2,2,2,2,2,2,2,2,2,2]
const NC_ST_DGMIN = Float32[3,3,3,3,3,3,3,3,3,3,3,7]
const NC_ST_DIAM = Float32[0.3,0.4,0.3,0.3,0.2,0.2,0.2,0.2,0.3,0.5,0.2,0.3]
# htgr5.f
const NC_HG_IMETH = Int[1,1,1,1,2,1,2,2,1,1,2,3]
const NC_HG_HCON = Float32[-2.193,-2.193,-2.193,-2.193,3.560,-2.193,3.817,3.385,-2.193,-2.193,3.385,-2.193]
const NC_HG_HBA  = Float32[-0.00828,-0.00828,-0.00828,-0.00828,-0.54648,-0.00828,-0.78296,-0.58984,-0.00828,-0.00828,-0.54984,-0.00828]
# htdbh.f SISKIY(sp,1:3) = P2,P3,P4
const NC_HD_P2 = Float32[523.0987,819.8690,523.0987,604.8450,160.6821,1530.3300,48.6795,679.1972,202.8860,1348.0419,679.1972,595.1068]
const NC_HD_P3 = Float32[5.7243,6.4531,5.7243,5.9835,4.1677,7.0811,8.9420,5.5698,8.7469,7.0463,5.5698,5.8103]
const NC_HD_P4 = Float32[-0.4109,-0.3434,-0.4109,-0.3789,-0.4954,-0.2544,-1.4832,-0.3074,-0.8317,-0.3076,-0.3074,-0.3821]

"nc/htgr5.f — small-tree height increment (3 methods by species). CR is the crown ratio on the ICR/10 (0–10)
scale, NOT the 0–1 proportion (regent.f:183 `CR=ICR(I)/10.0`). Final floor `IF(HTGR.LE.0.0)HTGR=0.01`."
@inline function nc_htgr5(sp::Int, ssite::Float32, baa::Float32, relht::Float32, cr::Float32, h::Float32)
    im = NC_HG_IMETH[sp]
    htgr = if im == 2
        b = baa <= 5.0f0 ? 5.0f0 : baa
        fexp(NC_HG_HCON[sp] + NC_HG_HBA[sp] * flog(b))
    elseif im == 1
        NC_HG_HCON[sp] + relht * 4.292f0 + 0.0566f0 * cr * cr +
            0.1699f0 * h + NC_HG_HBA[sp] * baa + 0.00768f0 * ssite
    else                                              # sp12 redwood — site-age curve
        htmax = 2.242202f0 * ssite
        if htmax - h <= 1.0f0
            0.0f0
        else
            age1 = (1.0f0 / -0.010742f0) * flog(1.0f0 - fpow(h / 2.242202f0 / ssite, 1.0f0 / 0.919076f0))
            age2 = age1 + 5.0f0
            h1 = 2.242202f0 * ssite * fpow(1.0f0 - fexp(-0.010742f0 * age1), 0.919076f0)
            h2 = 2.242202f0 * ssite * fpow(1.0f0 - fexp(-0.010742f0 * age2), 0.919076f0)
            h2 - h1
        end
    end
    htgr <= 0.0f0 && (htgr = 0.01f0)                  # htgr5.f: IF(HTGR.LE.0.0)HTGR=0.01
    return htgr
end

"nc/htdbh.f — Curtis-Arney HT-DBH (SISKIY, all forests). mode 1: H→D."
@inline function nc_htdbh_d(sp::Int, h::Float32)
    p2 = NC_HD_P2[sp]; p3 = NC_HD_P3[sp]; p4 = NC_HD_P4[sp]
    hlim = 4.5f0 + p2 * fexp(-p3 * fpow(3.0f0, p4))            # H at D=3 (curve→linear break)
    if h > hlim
        return fexp(flog((flog(h - 4.5f0) - flog(p2)) / (-p3)) / p4)
    else
        return ((h - 4.51f0) * 2.7f0) / (hlim - 4.51f0) + 0.3f0
    end
end

"nc/htdbh.f MODE=0 — Curtis-Arney HT-DBH (SISKIY, all forests): predict total height H from DBH `d`.
D≥3 uses the Curtis-Arney exponential; 0.3<D<3 the linear break; result floored at 4.5 (htdbh.f:79-90).
Used by the CRATET missing-height / broken-top NORMHT dub (dub_missing_heights!) — NC has no :htdbh_p2
blockdata column, so the generic `_htdbh_height` gave a too-short height ⇒ broken-top redwood NORMHT
clamped to the recorded (broken) height ⇒ 16-31% low cubic/board volume on large old redwood."
@inline function nc_htdbh_h(sp::Int, d::Float32)::Float32
    p2 = NC_HD_P2[sp]; p3 = NC_HD_P3[sp]; p4 = NC_HD_P4[sp]
    h = d >= 3.0f0 ? 4.5f0 + p2 * fexp(-p3 * fpow(d, p4)) :
        ((4.5f0 + p2 * fexp(-p3 * fpow(3.0f0, p4)) - 4.51f0) * (d - 0.3f0) / 2.7f0) + 4.51f0
    h <= 4.5f0 && (h = 4.5f0)
    return h
end

"nc/dgbnd.f (IE) — cap the small-tree diameter increment at the TREESZCP size cap (inert at the 999 default)."
@inline function nc_dgbnd(sp::Int, dbh::Float32, ddg::Float32, sizcap1::Float32, sizcap3::Float32)::Float32
    if (dbh + ddg) > sizcap1 && sizcap3 < 1.5f0
        ddg = sizcap1 - dbh
        ddg < 0.01f0 && (ddg = 0.01f0)
    end
    return ddg
end

"NC redwood (sp12) small-tree DG blend + DIAM floor + DGBND for one record (regent.f:305-327): DG =
DGSM·(1−XDWT) + DGLT·XDWT, where DGLT is THIS record's large-tree DG (central/upper/lower differ under tripling)."
@inline function _nc_rw_blend(dgsm2::Float32, xdwt::Float32, dglt::Float32, d::Float32,
                              sp::Int, cap1::Float32, cap3::Float32)::Float32
    dg = dgsm2 * (1f0 - xdwt) + dglt * xdwt
    (d + dg) < NC_ST_DIAM[sp] && (dg = NC_ST_DIAM[sp] - d)     # regent.f:319 DIAM floor
    return nc_dgbnd(sp, d, dg, cap1, cap3)                     # regent.f:327 DGBND
end

"""nc/regent.f LSTART calibration (DO 90 loop): the small-tree HEIGHT-increment CON = RHCON·fexp(HCOR).
HCOR = ln(CORNEW), CORNEW = Σ(observed·P)/Σ(predicted·P) over DBH<5 trees with a measured HTG — observed =
the measured height increment scaled to 5-yr (HTG·SCALE3, SCALE3=REGYR/FINTH), predicted = the raw HTGR5 on
the BACKDATED height/stand (RHCON=1). Trapped to CORNEW∈[0.0821,12.1825] (±2.5 SD of ln), else reset to 1.
Called from calibrate_diameter_growth! where t.dbh is backdated and t.ht_growth holds the measured increment.
The RAW HCOR goes into htg_cor_init; the shared dgdriv per-cycle attenuation (diameter_growth.jl:1133-1136)
produces the applied htg_cor_small = WCI + cormlt_h·(HCOR_init−WCI), WCI=dg_cor_goal (=0 for a species with no
diameter COR), cormlt_h=fexp(−0.02773·elapsed_end). At cyc1 (elapsed 0, +5) ⇒ CON=fexp(0.8705·HCOR_raw). Verified
vs FVSnc_g16 dumps: BO raw HCOR −0.8059 (CORNEW 0.4467) → applied −0.7016 (CON 0.4958) → cyc1 HTG/DG bit-exact."""
function nc_regent_hcor_init!(s::StandState, isct::AbstractMatrix, ind1::AbstractVector,
                              saved_dbh::AbstractVector, avh::Float32)
    p, t, c = s.plot, s.trees, s.calib
    t.n == 0 && return s
    s.control.growth_ifinth == 0 && return s                  # regent.f:359 IF(IFINTH.EQ.0) GOTO 100
    # regent.f:388-390 XBA=BA: the plot BA the nc/cratet.f:172 backdating DENSE left (dead-inclusive at FINT/FINTM-
    # scaled PROB, IMC=9 older dead at D=0) — the crown-init snapshot. RELHT=H/AVH uses AVH from cratet.f:530 AVHT40
    # (current heights), passed in. The former ad-hoc BA (live backdated + unscaled HISTORY-6/7 dead) and the
    # backdated-state AVH put every species' SNX ~2.7% off live on the REGCAL fixture.
    ba = c.cratet_ba
    ba <= 0f0 && (ba = 0.1f0)
    regyr = NC_REGYR
    finth = s.control.growth_finth > 0f0 ? s.control.growth_finth : 10f0   # NC measured HTG period = IFINT=10
    scale3 = regyr / finth                                                 # regent.f:363 SCALE3=REGYR/FINTH
    @inbounds for sp in 1:12
        i1 = isct[sp, 1]; i1 == 0 && continue
        i2 = isct[sp, 2]
        ssite = p.sp_site_index[sp]
        snx = 0f0; sny = 0f0; snp = 0f0; nh = 0
        for k in i1:i2
            i = ind1[k]
            saved_dbh[i] >= 5.0f0 && continue                 # regent.f:378 DBH<5 (current dbh)
            hg = t.ht_growth[i]; hg < 0.001f0 && continue     # regent.f:379 measured HTG≥0.001
            hb = s.control.growth_ihtg < 2 ? t.height[i] - hg : t.height[i]   # nc/regent.f:384 IF(IHTG.LT.2) H=H-HTG
            hb < 0.01f0 && continue
            cr = Float32(t.crown_pct[i]) * 0.1f0              # ICR/10
            relht = avh > 0f0 ? hb / avh : 1f0; relht > 1.5f0 && (relht = 1.5f0)
            xhtgr = nc_htgr5(sp, ssite, ba, relht, cr, hb)    # predicted; RHCON=1 ⇒ EDH=XHTGR
            term = hg * scale3                                # observed, scaled to 5-yr
            pr = t.tpa[i]
            snx += xhtgr * pr; sny += term * pr; snp += pr; nh += 1
        end
        nh < 5 && continue                                    # NCALHT
        snx /= snp; sny /= snp
        cornew = snx > 0f0 ? sny / snx : 1f0
        cornew <= 0f0 && (cornew = 1f-4)
        (cornew < 0.0821f0 || cornew > 12.1825f0) && (cornew = 1f0)
        # RAW HCOR into htg_cor_init; the shared dgdriv attenuation (calibrate/diameter_growth.jl:1133-1136)
        # produces the applied htg_cor_small = WCI + cormlt_h·(HCOR_init − WCI) each cycle (WCI=dg_cor_goal).
        c.htg_cor_init[sp] = flog(cornew)
    end
    return s
end

function small_tree_growth!(s::StandState, stash, ::Klamath; fint::Float32 = 10.0f0)
    p, t = s.plot, s.trees
    dens = s.density; sd = s.coef.species
    ba = p.basal_area; avh = p.avg_height
    scale = fint / NC_REGYR
    scale2 = NC_REGYR / fint         # regent.f:126 SCALE2=YR/FNT (=1 for the native 5-yr NC cycle)
    bark_a = s.calib.bark_a; bark_b = s.calib.bark_b
    hcor = s.calib.htg_cor_small     # regent.f:159 CON = RHCON(=1)·exp(HCOR)
    trip = stash !== nothing         # tripling active: replicate the small-tree DG/HTG onto upper/lower records
    @inbounds for i in 1:t.n
        t.tpa[i] <= 0f0 && continue
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= NC_ST_XMAX[sp] && continue
        h = t.height[i]
        ssite = p.sp_site_index[sp]
        pt_i = Int(t.plot_id[i])
        tpccf = (1 <= pt_i <= length(dens.point_ccf)) ? dens.point_ccf[pt_i] : 0f0
        relht = (h > 0f0 && avh > 0f0) ? h / avh : 1f0
        tpccf <= 75.0f0 && (relht = 1.0f0 - ((relht - 1.0f0) / 75.0f0) * tpccf)
        relht > 1.5f0 && (relht = 1.5f0)
        cr = Float32(t.crown_pct[i]) * 0.1f0                   # regent.f:183 CR=ICR(I)/10.0 (0–10 scale, NOT /100)
        htgr = nc_htgr5(sp, ssite, ba, relht, cr, h) * scale * fexp(hcor[sp])   # ·CON(=exp(HCOR)); XRHMLT=1
        # height: XWT blend with the large-tree HTG (already in t.ht_growth[i])
        xmn = NC_ST_XMIN[sp]; xmx = NC_ST_XMAX[sp]
        xwt = d <= xmn ? 0f0 : (d - xmn) / (xmx - xmn)
        lthg = t.ht_growth[i]
        htg = sp == 12 ? ((htgr + lthg) / 2f0) * (1f0 - xwt) + xwt * lthg :
                         htgr * (1f0 - xwt) + xwt * lthg
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        t.ht_growth[i] = htg
        # small-tree DG (D < DGMIN): HT-DBH (SISKIY) DK/DKK (regent.f:243-320).
        if d < NC_ST_DGMIN[sp]
            cap1 = s.control.sp_size_cap[sp, 1]; cap3 = s.control.sp_size_cap[sp, 3]
            hk = h + htg
            if hk <= 4.5f0
                # regent.f:245-247 — DBH(K)=D+HK·0.001, DG(K)=0 (tiny sub-breast-height nudge).
                t.dbh[i] = d + hk * 0.001f0
                t.diam_growth[i] = 0f0
                if trip                                 # tripled records: same nudge ⇒ DG=0 (regent.f re-entry)
                    stash.dgU[i] = 0f0; stash.dgL[i] = 0f0
                end
            else
                dk = nc_htdbh_d(sp, hk)
                dkk = h <= 4.5f0 ? d : nc_htdbh_d(sp, h)
                bark = variant_bratio(s, sp, d)          # nc/bratio.f
                dgsm = (dk < 0f0 || dkk < 0f0) ? htg * 0.2f0 * bark : (dk - dkk) * bark
                dgsm < 0f0 && (dgsm = 0f0)
                # regent.f:305-318 — DDS-space transform (SCALE2=YR/FNT), then back to DBH increment.
                # NON-redwood uses the pure small-tree DGSM (NO large-tree XWT blend, unlike height);
                # ONLY redwood (sp12) blends small/large DG via XDWT=(D-XMN)/(DGMIN-XMN).
                dds = dgsm * (2f0 * bark * d + dgsm) * scale2
                if sp == 12
                    # REDWOOD: DG = DGSM·(1−XDWT) + DGLT·XDWT where DGLT = the record's OWN large-tree DG
                    # (regent.f:305-311). At D≤XMN=2 XDWT=0 ⇒ pure small-tree (all tripled records equal); at
                    # 2<D<DGMIN=7 XDWT>0 ⇒ each tripled record blends with its own large-tree DG (upper/lower
                    # DIFFER, verified vs FVSnc_g16). The DIAM floor + DGBND (regent.f:319-327) apply per record.
                    dgsm2 = sqrt(fpow(d * bark, 2f0) + dds) - bark * d
                    xdwt = d <= xmn ? 0f0 : (d - xmn) / (NC_ST_DGMIN[sp] - xmn)
                    dg = _nc_rw_blend(dgsm2, xdwt, t.diam_growth[i], d, sp, cap1, cap3)
                    if trip
                        stash.dgU[i] = _nc_rw_blend(dgsm2, xdwt, stash.dgU[i], d, sp, cap1, cap3)
                        stash.dgL[i] = _nc_rw_blend(dgsm2, xdwt, stash.dgL[i], d, sp, cap1, cap3)
                    end
                    t.diam_growth[i] = dg
                else
                    dg = sqrt(fpow(d * bark, 2f0) + dds) - bark * d
                    (d + dg) < NC_ST_DIAM[sp] && (dg = NC_ST_DIAM[sp] - d)   # regent.f:319 DIAM floor
                    dg = nc_dgbnd(sp, d, dg, cap1, cap3)
                    t.diam_growth[i] = dg
                    if trip                             # non-RW small tree: DG is DETERMINISTIC ⇒ tripled = central
                        stash.dgU[i] = dg; stash.dgL[i] = dg
                    end
                end
            end
            # REGENT tripling override (regent.f:339-342 GO TO 18 re-entry): for a small tree (D<DGMIN)
            # the L=1,2 tripled records inherit the small-tree HTG (deterministic) and the DG set above
            # (identical for non-RW / RW below XMN; the RW large-tree-blend spread above 2"). is_small=true
            # routes htgU/htgL through triple_records!. Without this the tripled upper/lower records (0.40 of
            # TPA) kept the LARGE-tree DDS DG stashed by diameter_growth! — far smaller than the small-tree DG
            # for the CA hardwoods (BO/TO/RW, newly mapped by the 442-row crosswalk) ⇒ ~31% BA deficit.
            if trip
                stash.htgU[i] = htg; stash.htgL[i] = htg; stash.is_small[i] = true
            end
        end
    end
    return s
end

# nc/esgent.f: SPESRT, CALL REGENT(.TRUE.,ITRNIN), then HT += HTG·WK4 (WK4=HTIMLT=1) and the HHTMAX cap. jl had NO NC
# birth-cycle growth (planted/natural cohorts sat at their ESTAB height for the cycle). REGENT(LESTB) (nc/regent.f:
# 101-330): FNT=FINT−5 (LSKIPH when FINT≤5), SCALE=FNT/REGYR, SCALE2=YR/FNT; DO 13 (STORAGE order) draws each new
# record's crown; DO 30 ISPC / DO 25 I3 (IND1) grows only I≥ITRNIN with XWT=0: HTGR5·CON·(XRHMLT·SCALE), SIZCAP; below
# DGMIN the LESTB diameter DBH=max(DK,DIAM)+0.001·HK (DG=DBH) or D+0.001·HK under breast height, then DGBND. REGENT
# reads the gradd.f:192 DENSE (post-growth, pre-regen): ba/avh/pccf.
const NC_HHTMAX = Float32[27, 31, 25, 25, 26, 24, 28, 20, 20, 18, 26, 25]   # nc/blkdat.f:80-81
function nc_esgent!(s::StandState, nstart::Int; fint::Float32 = 5.0f0, ba_pre::Float32 = -1f0,
                    avh_pre::Float32 = -1f0, pccf_pre::Vector{Float32} = Float32[])
    p, t = s.plot, s.trees
    n = t.n; nstart >= n && return s
    ba = ba_pre >= 0f0 ? ba_pre : p.basal_area
    avh = avh_pre >= 0f0 ? avh_pre : p.avg_height
    pccfv = isempty(pccf_pre) ? s.density.point_ccf : pccf_pre
    @inline _pccf(i) = (pt = Int(t.plot_id[i]); (1 <= pt <= length(pccfv)) ? pccfv[pt] : 0f0)
    yr = htg_period(s.variant)
    lskiph = fint <= 5.0f0
    fnt = lskiph ? fint : fint - 5.0f0
    scale = fnt / NC_REGYR; scale2 = yr / fnt
    cur_year = current_cycle_year(s)
    @inbounds for i in (nstart + 1):n                         # estab.f: new records enter with DG=HTG=0
        t.diam_growth[i] = 0f0; t.ht_growth[i] = 0f0
    end
    @inbounds for i in (nstart + 1):n                         # regent.f:136-147 DO 13 crown (storage order)
        cr0 = 0.89722f0 - 0.0000461f0 * _pccf(i)
        ran = 0f0
        while true; ran = bachlo(s.rng, 0f0, 1f0); (-1f0 <= ran <= 1f0) && break; end
        cr0 = cr0 + 0.07985f0 * ran
        cr0 > 0.90f0 && (cr0 = 0.90f0); cr0 < 0.20f0 && (cr0 = 0.20f0)
        t.crown_pct[i] = trunc(Int32, cr0 * 100f0 + 0.5f0)
    end
    order = sort(collect((nstart + 1):n); by = i -> (Int(t.species[i]), i))   # SPESRT IND1 (new records only grow)
    hcor = s.calib.htg_cor_small
    @inbounds for i in order
        sp = Int(t.species[i]); (1 <= sp <= 12) || continue
        d = t.dbh[i]; d >= NC_ST_XMAX[sp] && continue
        h = t.height[i]
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year) * scale   # regent.f:154 XRHGRO=XRHMLT·SCALE
        local htg::Float32
        if lskiph
            htg = 0f0
        else
            xba = ba <= 0f0 ? 0.1f0 : ba
            cr = Float32(t.crown_pct[i]) / 10.0f0
            relht = (h > 0f0 && avh > 0f0) ? h / avh : 1f0
            tpccf = _pccf(i)
            tpccf <= 75.0f0 && (relht = 1.0f0 - ((relht - 1.0f0) / 75.0f0) * tpccf)
            relht > 1.5f0 && (relht = 1.5f0)
            htgr = nc_htgr5(sp, p.sp_site_index[sp], xba, relht, cr, h) * fexp(hcor[sp]) * xrhgro   # ·CON·XRHGRO
            htg = htgr                                         # XWT=0 (LESTB); RW HTGR2=HTGR
            cap = s.control.sp_size_cap[sp, 4]
            (h + htg > cap) && (htg = cap - h; htg < 0.1f0 && (htg = 0.1f0))
        end
        t.ht_growth[i] = htg
        if d < NC_ST_DGMIN[sp]                                 # regent.f:243-327
            hk = h + htg
            local dbhk::Float32, dg::Float32
            if hk <= 4.5f0
                dbhk = d + hk * 0.001f0; dg = 0f0
            else
                dk = nc_htdbh_d(sp, hk)                        # LHTDRG off ⇒ HTDBH (SISKIY)
                dbhk = dk < NC_ST_DIAM[sp] ? NC_ST_DIAM[sp] : dk
                dbhk = dbhk + 0.001f0 * hk
                dg = dbhk
                (dbhk + dg) < NC_ST_DIAM[sp] && (dg = NC_ST_DIAM[sp] - dbhk)
            end
            dg = nc_dgbnd(sp, dbhk, dg, s.control.sp_size_cap[sp, 1], s.control.sp_size_cap[sp, 3])
            t.dbh[i] = dbhk; t.diam_growth[i] = dg
        end
    end
    @inbounds for i in (nstart + 1):n                         # esgent.f:58-71
        sp = Int(t.species[i])
        t.height[i] = t.height[i] + t.ht_growth[i]
        (1 <= sp <= 12 && t.height[i] > NC_HHTMAX[sp]) && (t.height[i] = NC_HHTMAX[sp])
    end
    return s
end
