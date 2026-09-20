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
            # calibration DO 49 I=1,ITRN (regent.f:1090-1109) has NO D>=3 gate (unlike the growth pass DO 6 :246):
            # every live record's measured increment feeds the subcycle projection.
            d1 = t.dbh[i]                                   # WK3(I) = backdated dbh
            sp = Int(t.species[i]); pr = t.tpa[i]
            bark = ie_bratio(sp, d1)
            d2 = d1 + t.diam_growth[i] / bark
            b1 = 0.005454154f0 * d1 * d1; b2 = 0.005454154f0 * d2 * d2
            c1 = ie_tree_ccf(sp, d1) * pr; c2 = ie_tree_ccf(sp, d2) * pr     # CCFCAL = CCFT·P (regent.f:1098-1106)
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]; pn = pr * fpowi(0.985f0, k)
                rdnext[j] += Float32(k) * ci / pr * pn; banext[j] += Float32(k) * bi * pn
            end
        end
    end
    let f = get(ENV, "JLRH", ""); isempty(f) || open(f, "a") do io; println(io, "HCI reldm1=", s.plot.relative_density_prev, " relden=", relden, " rdnext1=", rdnext[1], " pccf23=", (23 <= t.n ? s.density.point_ccf[Int(t.plot_id[23])] : -1f0), " plot23=", (23 <= t.n ? t.plot_id[23] : -1)); end; end   # TEMP-DEBUG
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
                htgrl = rhcon[sp] + IE_RG_RHLH[sp]*flog(hk) + IE_RG_RHCCF[sp]*rdnext[j] +
                        IE_RG_RHBAL[sp]*bal
                hk += fexp(htgrl)                          # regent.f:1174-1175 (NO scale in the calib pass)
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
        c.htg_cor_init[sp] = flog(cornew)
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

"""
REGENT DUBSCR crown dub (ie/regent.f:999-1021). A small tree that will reach DBH≥3 THIS cycle and does not
yet carry a crown (`ICRK==0`) gets a crown ratio dubbed from DUBSCR — `IF(DNEW.GE.3.0)THEN; IF(ICRK.EQ.0)
CALL DUBSCR(ISPC,DNEW,HK,...); ICR(K)=ICRK`. DNEW is the projected end-of-cycle DBH: DG(K) (the record's
FINT-basis increment) is put on the full-cycle DDS scale (`DDS2=DG·(2·BARK·D+DG)·SCALE`, SCALE=NTYR/YR) and
`DNEW=D+SQRT((D·BARK)²+DDS2)−BARK·D`.

Without this the crown-less regen — most visibly quaking aspen, whose per-cycle crown update is deferred in
`crown_ratio_update!` (crown.f UTTVAR), so nothing else ever assigns it — keeps `crown_pct==0`. Aspen with
3≤D<XMAX(=4) SKIP the regent diameter dub (regent.f:861 `IF(D.GE.3.0)GO TO 23`) and take their diameter from
the large-tree DGFASP, whose POT carries a crown term (`4.510E-2·ASPCR·D^.67266`); with CR=0 the DDS collapses
(~0.42 vs ~1.45 measured on the oracle) so the cohort never crosses 4″ into the productive regime and the
stand BA/QMD stall (~19% aspen BA deficit). Fires only for `crown_pct==0` records crossing 3″ (the oracle draws
the DUBSCR error here too, so this also re-aligns the RNG for those trees); inert where the crown already exists
(every conifer regen — the NIVAR/CRVAR/LPIJU crown paths in crown_ratio_update! keep them non-zero).
"""
@inline function _ie_regent_dubscr_crown!(s::StandState, i::Int, sp::Int, d::Float32, hk::Float32,
                                          dg_this::Float32, ntyr::Int, yr::Float32, ba::Float32, dgsd::Float32)
    s.trees.crown_pct[i] != 0 && return                      # regent.f:1014 ICRK.EQ.0 guard
    bark = ie_bratio(sp, d)                                   # regent.f:1006 BRATIO(ISPC,DBH(K),HT(K))
    scale = yr > 0f0 ? Float32(ntyr) / yr : 1f0              # DO-30 SCALE = FLOAT(NTYR)/YR
    dds2 = dg_this * (2f0 * bark * d + dg_this) * scale       # regent.f:1007
    dg2 = sqrt((d * bark)^2 + dds2) - bark * d                # regent.f:1008
    dg2 < 0f0 && (dg2 = 0f0)                                  # regent.f:1009
    dnew = d + dg2                                            # regent.f:1010
    dnew < 3f0 && return                                      # regent.f:1013 DNEW.GE.3.0
    s.trees.crown_pct[i] = Int32(ie_dubscr(s.rng, sp, dnew, hk, ba, dgsd))   # regent.f:1015-1018
    return
end

# HTGF per-tripled-copy large-tree height growth (ie/htgf.f:318-347). Each tripled copy's large-tree HTG
# used in the REGENT XWT blend is recomputed from the COPY's OWN diameter growth DG(ITFN)/DG(ITFN+1) —
# NOT the central record's — using the SAME CON (central DBH/HT). NIVAR (recompute) species only; the
# Weibull species (13,17-22) keep the central value (htgf HTG(ITFN)=TEMHTG). `con`,`hdgcof` are the
# central tree's htgf terms; `dg` the copy's large-tree DG; tail = *scale*xht then the HT-based SIZCAP.
@inline function _ie_htgf_copy_large(con::Float32, hdgcof::Float32, dg::Float32,
                                     scale::Float32, xht::Float32, hti::Float32, cap::Float32)
    dg <= 0.0f0 && return 0.1f0 * scale * xht               # ln(DG) undefined ⇒ FVS's pre-tail 0.1 floor
    v = fexp(con + hdgcof * flog(dg)) + IE_HTBIAS
    v < 0.1f0 && (v = 0.1f0)
    v = v * scale * xht
    (hti + v > cap) && (v = max(cap - hti, 0.1f0))          # htgf.f:329-331 SIZCAP on HT(ITFN)=central HT
    return v
end

"""IE `small_tree_growth!` (ie/regent.f). Overrides DG/HTG for small trees (D<XMAX). NIVAR path is the
bulk (iet01); special species ported faithfully."""
function small_tree_growth!(s::StandState, stash, ::InlandEmpire; fint::Float32 = 10.0f0)
    let f = get(ENV, "JLRH", ""); isempty(f) || open(f, "a") do io; println(io, "RGE", lpad(Int(s.control.cycle)+1,3), " ", s.rng.s0); end; end   # TEMP-DEBUG
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    sd = s.coef.species                                 # blkdat HT-DBH :ht1/:ht2 for the aspen log-DK
    t.n == 0 && return s
    n = t.n
    rhcon = ie_regcons!(s)
    # htgf.f HTCONS habitat terms — needed to recompute each tripled COPY's large-tree HTG (its own DG) for
    # the NIVAR XWT blend (ie/htgf.f:318-347). Stand-level; resolved once.
    itype_hg = Int(p.habitat_input); iht_hg = (1 <= itype_hg <= 30) ? IE_HTMAPHAB[itype_hg] : 1
    hghch_hg = IE_HGHC[iht_hg]; h2cof_hg = IE_HGH2[iht_hg]; hdgcof_hg = IE_HGLDD[iht_hg]
    htgf_scale = fint / 10.0f0                              # htgf SCALE=FINT/YR, YR=10 (matches height_growth!)
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
            # CCFCAL returns CCFT·P (ccfcal.f:81) ⇒ C1/C2 carry P, and FVS evaluates RDNEXT += ((K*CI)/P)*PN with
            # PN=P*0.985**K (integer power ⇒ __powisf2) — regent.f:262-273. Cancelling the P algebraically (the old
            # `k*ci*pn` on bare CCFT) rounds differently: IE cycle-4 subcycle-2 RDJ was 2 ULP off (bare-plot fixture).
            c1 = ie_tree_ccf(sp, d1) * pr; c2 = ie_tree_ccf(sp, d2) * pr
            bi = (b2 - b1) / 10.0f0; ci = (c2 - c1) / 10.0f0
            k = 0
            for j in 2:nper
                k += kper[j-1]
                pn = pr * fpowi(0.985f0, k)
                rdnext[j] += Float32(k) * ci / pr * pn
                banext[j] += Float32(k) * bi * pn
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
    # WK4(I) = clgmult's per-tree CLIMATE growth multiplier (regent.f:596/598 H2=H1+…·WK4). In the oracle
    # CLGMULT fills WK4 once (dgdriv.f:153) and BOTH the large-tree DDS (dgdriv.f:217, jl apply_climate_dds!)
    # AND the small-tree regent height read it. jl applied it to diameter but OMITTED it in the height path,
    # so under a suppressing climate score the NIVAR/TTVAR seedling over-grew ~1/WK4 (tripping breast-height a
    # cycle early ⇒ escaping self-thinning — the IE climate-regime dense-phase under-kill). WK4≡1 without a
    # CLIMATE keyword (clgmult.f:55,60 return early) ⇒ IEEE-exact no-op on every non-climate stand.
    wk4 = ones(Float32, n)
    if s.climate !== nothing && s.climate.active
        _c = s.climate; _cd = _c.data; _ix = _c.indices
        if !(_ix[:mtcm]==0 || _ix[:mmin]==0 || _ix[:dd0]==0 || _ix[:d100]==0 || _ix[:dd5]==0 || _ix[:gsp]==0)
            _ty = Float32(cur_year) + fint / 2f0                 # THISYR = IY(ICYC)+FINT/2 (matches apply_climate_dds!)
            _A(sym, yr) = algslp(yr, _cd.years, view(_cd.attrs, :, _ix[sym]))
            _smi(yr) = (g = _A(:gsp, yr); g > 0f0 ? _A(:dd5, yr) / g : 0f0)
            _xgsite = _ix[:pSite] > 0 ? clim_xgsite(_A(:pSite, _ty), _A(:pSite, Float32(_c.inv_year))) : 1f0
            _mtcm_now = _A(:mtcm, _ty); _mmin_now = _A(:mmin, _ty); _smi_now = _smi(_ty)
            _ns = length(_c.plant_symbols)
            _vscore = ones(Float32, _ns)
            @inbounds for sp in 1:_ns; _, _vscore[sp] = species_vscore(_cd, _c.plant_symbols[sp], _ty); end
            @inbounds for i in 1:n
                t.dbh[i] <= 0f0 && continue
                sp = Int(t.species[i]); (sp < 1 || sp > _ns) && continue
                _by = _ty - t.birth_age[i]
                _xdf = leites_xdf(_mtcm_now, _A(:mtcm, _by))
                _xwl = leites_xwl(_mmin_now, _A(:mmin, _by), _A(:dd0, _by))
                _xpp = leites_xpp(_smi_now, _smi(_by), _A(:d100, _by))
                _xr = clim_xrelgr(_c.plant_symbols[sp], _xdf, _xpp, _xwl)
                _, _tm = clim_treemult(_xgsite, _xr, _vscore[sp], _c.growmult[sp])
                wk4[i] = _tm
            end
        end
    end
    # FVS walks every REGENT tree loop species-major over IND1 (DO ISPC=1,MAXSP / DO I3=I1,I2 / I=IND1(I3)):
    # the subcycle density feedback RDNEXT/BANEXT(J+1) (regent.f:665-670) is a Float32 running sum, so its
    # accumulation ORDER is part of the result; record order drifted RDJ/BAJ by ULPs from subcycle 2 on.
    _isct = s.control.sp_count_tab; _ind1 = s.scratch.idx1
    _sp_order = Vector{Int}(undef, n); _no = 0
    @inbounds for sp in 1:MAXSP
        i1 = _isct[sp, 1]; i1 == 0 && continue
        i2 = _isct[sp, 2]
        for k in i1:i2
            (1 <= k <= length(_ind1)) || continue
            ii = Int(_ind1[k]); (1 <= ii <= n) || continue
            _no += 1; _sp_order[_no] = ii
        end
    end
    # ---- subcycle loop (regent.f:250-620) ----
    @inbounds for j in 1:nper
        baj = banext[j]; rdj = rdnext[j]
        scale = Float32(kper[j]) / regyr
        ky = 0; for jj in 1:j; ky += kper[jj]; end            # KY = cumulative years thru subcycle j (regent.f:353)
        surv = fpowi(0.985f0, ky)                              # 0.985**KY: integer power ⇒ libgcc __powisf2 (single precision)
        for _oi in 1:_no
            i = _sp_order[_oi]
            sp = Int(t.species[i]); d0 = t.dbh[i]
            d0 >= IE_RG_XMAX[sp] && continue
            t.tpa[i] <= 0.0f0 && continue
            nivar = sp <= 12 || sp == 14 || sp == 23
            if !nivar
                # UTVAR PI/JU (sp15,16): potential-height model, ONE pass only (regent.f:407 J>1 skip).
                if (sp == 15 || sp == 16) && j == 1
                    con = rhcon[sp] * fexp(c.htg_cor_small[sp])   # non-NIVAR CON = RHCON·EXP(HCOR) (regent.f:414)
                    h1 = wk3[i]; d = wk5[i]
                    sitear = p.sp_site_index[sp]
                    xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
                    xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
                    sj = sitear                                   # POTHTG uses raw SITEAR (regent.f:459/467); H==H1 at j=1
                    pothtg = ((sj/5f0)*(sj*1.5f0 - h1)/(sj*1.5f0)) * 0.83f0
                    crx = Float32(t.crown_pct[i]) / 100f0
                    vigor = 150f0*crx^3*fexp(-6f0*crx) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
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
                    con = rhcon[sp] * fexp(c.htg_cor_small[sp])
                    h1 = wk3[i]
                    si = p.sp_site_index[sp]
                    si > IE_RG_SHI[sp] && (si = IE_RG_SHI[sp])
                    si <= IE_RG_SLO[sp] && (si = IE_RG_SLO[sp] + 0.5f0)
                    relsi = (si - IE_RG_SLO[sp]) / (IE_RG_SHI[sp] - IE_RG_SLO[sp])
                    rsimod = 0.5f0 * (1f0 + relsi)
                    sitage = fpow(h1 * 12f0 * 2.54f0 / 26.9825f0, 0.8509f0)    # FINDAG (inverse Sheppard, regent.f:565)
                    hite1 = 26.9825f0 * fpow(sitage, 1.1752f0)
                    hite2 = 26.9825f0 * fpow(sitage + 10f0, 1.1752f0)
                    htgrl = (hite2 - hite1) / (2.54f0 * 12f0) * rsimod * con * 0.75f0
                    wk3[i] = h1 + htgrl * scale_ut                          # regent.f:600 (·SCALE=NTYR/YR)
                    # (aspen subcycle DBH/density-feedback omitted — single UT pass; final assembly is authoritative)
                elseif sp == 13 || sp == 17
                    # TTVAR (LM/PY): BETA/ZRAND height — EVERY subcycle (NOT one-pass; regent.f:507-538,598).
                    con = rhcon[sp] * fexp(c.htg_cor_small[sp])
                    h1 = wk3[i]
                    cr = Float32(t.crown_pct[i])
                    tpccf = clamp(relden, 25f0, 300f0)                      # PCCF≈stand CCF (single-point)
                    beta1 = fexp(1.17527f0 - 0.42124f0*flog(tpccf))
                    beta2 = fexp(-2.56002f0 - 0.58642f0*flog(tpccf))
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
                    wk3[i] = h1 + htgrl * scale * xrhgro * con * wk4[i]     # regent.f:598 ·SCALE ·CON ·WK4(I)
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
                    con = rhcon[sp] * fexp(c.htg_cor_small[sp])
                    h1 = wk3[i]
                    sj = p.sp_site_index[sp]                       # POTHTG uses raw SITEAR
                    pothtg = ((sj/5f0)*(sj*1.5f0 - h1)/(sj*1.5f0)) * 0.83f0
                    crx = Float32(t.crown_pct[i]) / 100f0
                    vigor = 150f0*crx^3*fexp(-6f0*crx) + 0.3f0; vigor > 1f0 && (vigor = 1f0)
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
            htgrl = con + IE_RG_RHLH[sp]*flog(h1) + IE_RG_RHCCF[sp]*rdj + IE_RG_RHBAL[sp]*bal
            h2 = h1 + fexp(htgrl) * scale * xrhgro * wk4[i]   # regent.f:596 ·WK4(I) = clgmult climate multiplier
            _rhdump("RH1", s, i, sp, (h1, htgrl, fexp(htgrl), scale, xrhgro, wk4[i], h2, con)); _rhdump("RH0", s, i, sp, (rdj, baj, bal))   # TEMP-DEBUG
            wk3[i] = h2
            # NIVAR diameter (regent.f:598-610): skip if last subcycle or D≥3 or H2≤4.5
            d2 = d
            if !(j >= nper || d >= 3.0f0 || h2 <= 4.5f0)
                ax = IE_RG_HHT1[sp]; bx = IE_RG_HHT2[sp]
                d1v = IE_RG_DIAM[sp] + dadj
                h1 > 4.5f0 && (d1v = ax * fpow(h1 - 4.5f0, bx) + dadj)
                d2v = ax * fpow(h2 - 4.5f0, bx) + dadj
                dgj = (d2v - d1v) * xrdgro; dgj < 0.0f0 && (dgj = 0.0f0)
                d2 = d + dgj
            end
            wk5[i] = d2
            # small-tree density feedback into the NEXT subcycle (regent.f:665-666): the growing small tree
            # adds its CCF/BA increase to RDNEXT/BANEXT(J+1). C1/C2 use CCFCAL (= CCFT·P) ⇒ ·pr; no /P here.
            if j < nper
                pr = t.tpa[i]
                # CCFCAL returns CCFT·P (ccfcal.f:81) ⇒ C1/C2 already carry P; FVS evaluates
                # ((KY*(C2-C1))/10.)*(0.985**KY) and ((BACON*D2*D2-B1)*P)*(0.985**KY) (regent.f:443-444, 668-670).
                c1 = ie_tree_ccf(sp, d) * pr; c2 = ie_tree_ccf(sp, d2) * pr
                rdnext[j+1] += Float32(ky) * (c2 - c1) / 10.0f0 * surv
                banext[j+1] += (0.005454154f0*d2*d2 - 0.005454154f0*d*d) * pr * surv
            end
        end
    end
    # ---- final assembly (regent.f DO 30 ISPC=1,MAXSP): HTGR1 + ZZRAN + XWT blend toward the large-tree HTG,
    #      then the D<3 diameter dub. XWT blend is CRITICAL — trees D∈[XMIN,XMAX] blend to the (validated) htgf
    #      value, and D≥3 keep their large-tree DG (only D<3 gets the small-tree dub). FVS iterates SPECIES-
    #      SORTED (DO 30 ISPC; DO 25 I3=I1,I2 via IND1) — the per-record ZZRAN (BACHLO) draws MUST happen in
    #      this order or the RNG stream desyncs vs live on multi-species stands (CR proved this). ----
    # ★ The DO-30 within-species order is the ISCT/IND1 LINKED-LIST order (setup.f: IND1 follows IBEGIN→IND2,
    #   the tree-insertion order), NOT the storage order. A plain `sortperm(species)` is species-major but
    #   ORDERS RECORDS BY STORAGE INDEX within a species — which only coincides with IND1 until establishment/
    #   mortality reshuffle the linked list. Once they diverge (the cycle AFTER the first AUTOES tally on a
    #   dense sp-9 stand), the per-record ZZRAN gets paired to the WRONG record, so a heavy small tree draws a
    #   different height error than the oracle, under-grows, and lingers below REGNBK=2.999 — inflating the
    #   small-tree pool that tips NSTORE (INT(ΣTPA/(prob1·300)+0.5)) and thus the AUTOES ingrowth count. Iterate
    #   the maintained ISCT/IND1 (the SAME order dgdriv.f's DGSCOR draws use — the diameter randomization is
    #   validated bit-exact) so the ZZRAN stream stays RNG-aligned with live FVS's DO-30 (regent.f:696-744).
    # REGENT stale-BARK (ie/regent.f): the small-tree DGK=(DK−DKK)·BARK for a CRVAR/UTVAR species uses a LEFTOVER
    # BARK — the recompute BARK=BRATIO(ISPC,DBH(K),HT(K)) (regent.f:978) happens AFTER, so DGK sees the PREVIOUS
    # tree's bark. FVS's first (subcycle) loop grows UTVAR/CRVAR species in J=1 ONLY (regent.f:407) but NIVAR every
    # subcycle, so the bark entering the second loop is the LAST NIVAR small tree's bark from the final subcycle.
    # Replicate that leftover, then carry it forward (updated by every tree that reaches the diameter dub, in the
    # same species-major order as FVS's DO-30 loop). Without this, a hardwood/aspen seedling that follows a
    # different-bark conifer seedling over-grows (its own bark ~0.95 vs the stale ~0.85–0.87) ⇒ under-mortality.
    prevbark = NaN32
    @inbounds for j in 1:nper, oi in 1:_no
        i = _sp_order[oi]; sp = Int(t.species[i]); d = t.dbh[i]
        (d >= IE_RG_XMAX[sp] || t.tpa[i] <= 0.0f0) && continue
        isniv = sp <= 12 || sp == 14 || sp == 23
        (!isniv && j > 1) && continue                    # regent.f:407 — UTVAR/CRVAR grow in subcycle 1 only
        prevbark = ie_bratio(sp, d)                      # regent.f:445 BARK=BRATIO in the first (subcycle) loop
    end
    @inbounds for oi in 1:_no
        i = _sp_order[oi]
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= IE_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        if !(sp <= 12 || sp == 14 || sp == 23)
            # UTVAR PI/JU (sp15,16): height increment + ZZRAN + linear height→DBH (regent.f:756-985).
            # regent OVERRIDES the large-tree dgf/htgf for pinyon/juniper (XMAX=99, XMIN=90 ⇒ XWT≡0).
            if sp == 15 || sp == 16
                # UTVAR PI/JU tripling (regent.f DO-25 918-loopback, line 1042-1044 GO TO 918): each record —
                # central + 2 copies — draws its OWN ZZRAN and gets its OWN height + linear-DK diameter dub; the
                # copies' DG/HTG MUST be stashed (dgU/dgL, htgU/htgL, is_small) or they FREEZE at the tiny
                # large-tree DGFASP DG (triple_records! always overwrites the copies' DG with dgU/dgL). XWT≡0
                # (XMIN=90); no D<3 gate (regent.f:862 GO TO 1020 for 15/16). DGK uses the STALE leftover BARK.
                h = t.height[i]
                sitear = p.sp_site_index[sp]
                xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
                xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
                htgr1 = wk3[i] - h                                # HK−H (UT: NOT floored; regent.f:756)
                xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
                xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)       # ≡0 for PI/JU (xmn=90)
                cap = s.control.sp_size_cap[sp, 4]
                bark = ie_bratio(sp, d)
                large_htg = t.ht_growth[i]                         # large-tree htgf value for the XWT blend
                nrec = stash !== nothing ? 3 : 1
                for l in 0:(nrec - 1)
                    zzran = 0f0
                    if dgsd >= 1.0f0
                        while true
                            zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                            (zzran <= 0.5f0 && zzran >= -2.0f0) && break   # CR/UT bound (regent.f:813)
                        end
                    end
                    htgr = (htgr1 + zzran*0.1f0) * xrhgro          # UT: ZZRAN·0.1 (regent.f:814)
                    htgr < 0.1f0 && (htgr = 0.1f0)
                    htg = htgr*(1f0 - xwt) + xwt*large_htg
                    (h + htg > cap) && (htg = max(cap - h, 0.1f0))
                    # diameter: linear height→DBH; DG on the DDS scale (regent.f:879-985)
                    hk = h + htg
                    dgv = 0f0; has_dg = false                      # hk<4.5 ⇒ DG=0, DBH unchanged (regent.f:883)
                    if hk >= 4.5f0
                        has_dg = true
                        dk = (hk - 4.5f0)*10f0/(sitear - 4.5f0); dk < 0.1f0 && (dk = 0.1f0)
                        dkk = (h - 4.5f0)*10f0/(sitear - 4.5f0); dkk < 0.1f0 && (dkk = 0.1f0)
                        h < 4.5f0 && (dkk = d)                      # regent.f:897 override
                        dgkbark = isnan(prevbark) ? bark : prevbark # regent.f:961 STALE BARK (see prevbark note)
                        dgk = (dk - dkk) * dgkbark * xrdgro         # regent.f:960
                        dgmx = IE_RG_DGMAX[sp]; dgk > dgmx && (dgk = dgmx)
                        dgk < 0f0 && (dgk = 0f0)
                        dds = dgk*(2f0*bark*d + dgk)*scale2         # regent.f:980 (DG(K)=DGK for CR/UT)
                        dgv = sqrt((d*bark)^2 + dds) - bark*d       # regent.f:981
                        (d + dgv) < IE_RG_DIAM[sp] && (dgv = IE_RG_DIAM[sp] - d)
                    end
                    if l == 0
                        t.ht_growth[i] = htg
                        t.diam_growth[i] = dgv
                    elseif l == 1
                        stash.htgU[i] = htg; stash.is_small[i] = true
                        stash.dgU[i] = dgv
                    else
                        stash.htgL[i] = htg
                        stash.dgL[i] = dgv
                    end
                    has_dg && (prevbark = bark)                    # regent.f:978 recompute (per record)
                end
            elseif sp == 18 || sp == 20 || sp == 21
                # UTVAR aspen (regent.f:756-985): per-RECORD ZZRAN height draw + XWT blend, then log-DK diameter
                # (D<3 only; else keep large-tree DGFASP). XMAX=4, XMIN=2. TRIPLING (regent.f:801-829 DO-25 over
                # IND1): each record — central + 2 copies — draws its OWN ZZRAN and gets its OWN height+diameter
                # dub; the copies' DG MUST be stashed (dgU/dgL) or they FREEZE at the tiny large-tree DGFASP DG and
                # get over-killed in later cycles (the dominant IE none-regime aspen/hardwood over-mortality). The
                # DGK bark is the STALE leftover BARK (prevbark): the recompute at regent.f:978 is only afterwards,
                # so the central sees the previous small tree's bark and each record then updates it for the next.
                h = t.height[i]
                xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
                xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
                htgr1 = wk3[i] - h
                xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
                xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
                cap = s.control.sp_size_cap[sp, 4]
                bx = sd[:ht2][sp]                                          # blkdat Wykoff HT2
                ax = c.ht_dbh_iabflg[sp] == 1 ? sd[:ht1][sp] : c.ht_dbh_aa[sp]   # regent.f:900-904
                bark = ie_bratio(sp, d)
                large_htg = t.ht_growth[i]                                 # large-tree htgf value for the XWT blend
                nrec = stash !== nothing ? 3 : 1
                for l in 0:(nrec - 1)
                    zzran = 0f0
                    if dgsd >= 1.0f0
                        while true
                            zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                            (zzran <= 0.5f0 && zzran >= -2.0f0) && break   # CR/UT bound (regent.f:813)
                        end
                    end
                    htgr = (htgr1 + zzran*0.1f0) * xrhgro
                    htgr < 0.1f0 && (htgr = 0.1f0)
                    htg = htgr*(1f0 - xwt) + xwt*large_htg
                    (h + htg > cap) && (htg = max(cap - h, 0.1f0))
                    dgv = 0f0; has_dg = false
                    if d < 3f0                                            # regent.f:860 D≥3 ⇒ keep large-tree DGFASP
                        hk = h + htg
                        if hk >= 4.5f0
                            has_dg = true
                            dk = (bx / (flog(hk - 4.5f0) - ax)) - 1f0; dk < 0.1f0 && (dk = 0.1f0)   # regent.f:905-906
                            dkk = h <= 4.5f0 ? d : (bx / (flog(h - 4.5f0) - ax)) - 1f0              # regent.f:907-911
                            dgkbark = isnan(prevbark) ? bark : prevbark                            # regent.f:961 STALE BARK
                            dgk = (dk - dkk) * dgkbark * xrdgro                                    # regent.f:960
                            dgmx = IE_RG_DGMAX[sp]; dgk > dgmx && (dgk = dgmx)
                            dgk < 0f0 && (dgk = 0f0)
                            dds = dgk*(2f0*bark*d + dgk)*scale2                                    # regent.f:980
                            dgv = sqrt((d*bark)^2 + dds) - bark*d                                  # regent.f:981
                            (d + dgv) < IE_RG_DIAM[sp] && (dgv = IE_RG_DIAM[sp] - d)
                        end
                    end
                    if l == 0
                        t.ht_growth[i] = htg
                        if d < 3f0
                            t.diam_growth[i] = dgv                        # dgv=0 when hk<4.5 (regent.f:883)
                            # regent.f:999-1021 — aspen D<3 reaches the DUBSCR block; assign a crown when it
                            # crosses 3″ (the 3≤D<4 range then uses DGFASP, which needs the crown).
                            _ie_regent_dubscr_crown!(s, i, sp, d, h + htg, dgv, ntyr, yr, ba, dgsd)
                        end
                    elseif l == 1
                        stash.htgU[i] = htg; stash.is_small[i] = true
                        (d < 3f0) && (stash.dgU[i] = dgv)
                    else
                        stash.htgL[i] = htg
                        (d < 3f0) && (stash.dgL[i] = dgv)
                    end
                    has_dg && (prevbark = bark)                           # regent.f:978 recompute (per record)
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
                        prevbark = bark                             # regent.f:978 recompute (TTVAR reaches the dub)
                    else
                        t.diam_growth[i] = 0f0
                    end
                end
                # TTVAR tripling (regent.f:1029-1041): the 2 copies are EXACT DUPLICATES of the central
                # (DBH,DG,HT,HTG all copied — NO new ZZRAN, NO recompute; TTVAR skips the 918-loopback). Stash the
                # central's HTG/DG so the copies GROW instead of freezing at the large-tree dgU (triple_records!
                # always overwrites the copies' DG). XMAX=3 ⇒ every processed TTVAR record is D<3 (dub always ran).
                if stash !== nothing
                    stash.htgU[i] = t.ht_growth[i]; stash.htgL[i] = t.ht_growth[i]; stash.is_small[i] = true
                    stash.dgU[i] = t.diam_growth[i]; stash.dgL[i] = t.diam_growth[i]
                end
            elseif sp == 19 || sp == 22
                # CRVAR CO/OH tripling (regent.f DO-25 918-loopback): each record — central + 2 copies — draws its
                # OWN ZZRAN (bound [−2,0.5]) and gets its OWN height (HTGR=(HTGR1+ZZRAN·0.2)·XRHGRO + XWT blend
                # [XMIN,XMAX] + HTG floor 0.1) and, for D<1, its OWN log-DK diameter dub; the copies' DG/HTG MUST
                # be stashed or they FREEZE at the tiny large-tree DGFASP DG. D≥1 keeps the large-tree dgf
                # (regent.f:859 GOTO 23) but STILL draws the per-copy ZZRAN and blends the height. DGK uses the
                # STALE leftover BARK (prevbark); hk<4.5 sets DBH directly (regent.f:883 DBH(K)=D+.001·HK), encoded
                # for the copies as a DG vs central_dbh (copy_tree! copies the central DBH, GRADD adds DG/bark).
                h = t.height[i]
                xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
                xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
                htgr1 = wk3[i] - h
                xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
                xwt = d <= xmn ? 0f0 : (d - xmn)/(xmx - xmn)
                cap = s.control.sp_size_cap[sp, 4]
                bark = ie_bratio(sp, d)
                bx = sd[:ht2][sp]                                             # blkdat Wykoff HT2 (regent.f:900)
                ax = c.ht_dbh_iabflg[sp] == 1 ? sd[:ht1][sp] : c.ht_dbh_aa[sp]   # regent.f:901-905
                large_htg = t.ht_growth[i]                                    # large-tree htgf value for the blend
                small_d = d < 1.0f0
                nrec = stash !== nothing ? 3 : 1
                central_dbh = d
                for l in 0:(nrec - 1)
                    zzran = 0f0
                    if dgsd >= 1.0f0
                        while true
                            zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                            (zzran <= 0.5f0 && zzran >= -2.0f0) && break
                        end
                    end
                    htgr = (htgr1 + zzran*0.2f0) * xrhgro                     # CRVAR: ZZRAN·0.2
                    htgr < 0.1f0 && (htgr = 0.1f0)
                    htg = htgr*(1f0 - xwt) + xwt*large_htg
                    htg < 0.1f0 && (htg = 0.1f0)                              # CRVAR HTG floor (regent.f:842)
                    (h + htg > cap) && (htg = max(cap - h, 0.1f0))
                    # CRVAR CO/OH DIAMETER (ie/regent.f:859-987). D≥1 keeps the large-tree DG (regent.f:859
                    # `IF(CRVAR.AND.D.GE.1.0)GOTO 23`). D<1 gets the regent log-DK Wykoff diameter growth on the
                    # DDS scale; hk<4.5 sets DBH directly (regent.f:883, DG=0).
                    dgv = 0f0; dbh_dir = -1.0f0; has_dg = false
                    if small_d
                        hk = h + htg
                        if hk < 4.5f0
                            dbh_dir = d + 0.001f0*hk                          # regent.f:883 DBH(K)=D+.001·HK, DG=0
                        else
                            has_dg = true
                            dk = (bx / (flog(hk - 4.5f0) - ax)) - 1f0; dk < 0.1f0 && (dk = 0.1f0)   # regent.f:906-907
                            dkk = h <= 4.5f0 ? d : (bx / (flog(h - 4.5f0) - ax)) - 1f0              # regent.f:908-912
                            dgkbark = isnan(prevbark) ? bark : prevbark      # regent.f:961 STALE BARK
                            if dk < 0f0 || dkk < 0f0                          # regent.f:957-959
                                dgk = htg*0.2f0*dgkbark*xrdgro
                            else
                                dgk = (dk - dkk)*dgkbark*xrdgro              # regent.f:961
                            end
                            dgmx = IE_RG_DGMAX[sp]; dgk > dgmx && (dgk = dgmx)   # regent.f:963
                            dgk < 0f0 && (dgk = 0f0)                          # regent.f:966
                            dds = dgk*(2f0*bark*d + dgk)*scale2               # regent.f:978-981 (DG(K)=DGK for CRVAR)
                            dgv = sqrt((d*bark)^2 + dds) - bark*d             # regent.f:982
                            (d + dgv) < IE_RG_DIAM[sp] && (dgv = IE_RG_DIAM[sp] - d)   # regent.f:984-986
                            dgv = dg_bound(nothing, nothing, sp, d, dgv, s.control.sp_size_cap)   # DGBND (regent.f:991)
                        end
                    end
                    if l == 0
                        t.ht_growth[i] = htg
                        if small_d
                            if dbh_dir >= 0.0f0
                                t.dbh[i] = dbh_dir; t.diam_growth[i] = 0f0; central_dbh = dbh_dir
                            else
                                t.diam_growth[i] = dgv
                            end
                        end
                    elseif l == 1
                        stash.htgU[i] = htg; stash.is_small[i] = true
                        small_d && (stash.dgU[i] = dbh_dir >= 0.0f0 ? (dbh_dir - central_dbh)*bark : dgv)
                    else
                        stash.htgL[i] = htg
                        small_d && (stash.dgL[i] = dbh_dir >= 0.0f0 ? (dbh_dir - central_dbh)*bark : dgv)
                    end
                    has_dg && (prevbark = bark)                              # regent.f:978 recompute (per record)
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
        # PER-COPY large-tree HTG for the XWT blend (ie/htgf.f:318-347). FVS's HTGF recomputes HTG(ITFN)/
        # HTG(ITFN+1) from the COPY's OWN large-tree DG (dgU=upper/FRU, dgL=lower/FRL) using the central CON,
        # so a tripled NIVAR copy blends toward a DIFFERENT large-tree value than the central. jl formerly used
        # the central `large_htg` for all three copies ⇒ the tripled copies' blended small-tree HEIGHTS were
        # wrong, shifting the <3" pool that feeds NSTORE/AUTOES. The copies' large-tree DG lives in stash.dgU/dgL
        # and is read HERE, before the l-loop overwrites it for small_d records.
        large_htg_u = large_htg; large_htg_l = large_htg
        if stash !== nothing
            htcon_hg = hghch_hg + IE_HGSC[sp]
            (s.control.htg_cor2_on && s.control.htg_cor2[sp] > 0.0f0) && (htcon_hg += flog(s.control.htg_cor2[sp]))
            con_hg = htcon_hg + h2cof_hg * h * h + IE_HGLD[sp] * flog(d) + IE_HGLH * flog(h)
            xht_hg = active_multiplier(s.control, :htg, sp, cur_year)
            large_htg_u = _ie_htgf_copy_large(con_hg, hdgcof_hg, stash.dgU[i], htgf_scale, xht_hg, h, cap)
            large_htg_l = _ie_htgf_copy_large(con_hg, hdgcof_hg, stash.dgL[i], htgf_scale, xht_hg, h, cap)
        end
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        small_d = d < 3.0f0
        diam = IE_RG_DIAM[sp]; bark = ie_bratio(sp, d)          # BARK at the pre-growth DBH (regent.f:983)
        # deterministic dub-bias terms (no ZZRAN) — computed once. DADJ (regent.f:866-872), D1 (regent.f:875-877).
        relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h - 4.5f0) / (ah - 4.5f0)
        relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
        dadj = delmax*relh*relh - 2.0f0*delmax*relh + 0.65f0
        d1v = diam + dadj; h > 4.5f0 && (d1v = ax * fpow(h - 4.5f0, bx) + dadj)
        # TRIPLING (regent.f:801-810 "IF TRIPLING, EACH TRIPLE GETS A NEW RANDOM"): draw a FRESH ZZRAN per
        # tripled record; central (l=0) → the tree, copies (l=1,2) → the stash (htgU/htgL + is_small + dgU/dgL).
        # The DBH dub is the FAITHFUL non-ESTAB path (regent.f:955-989): DG(K)=(DK−D1)·XRDGRO·BARK on the DDS
        # scale — DBH is UNCHANGED and grows later via GRADD (t.dbh += DG/bark). (The old code did DBH-direct,
        # which is the ESTAB-only branch (regent.f:945-951 LESTB) — wrong for the growth cycle; combined with a
        # 1-draw-per-tree ZZRAN it desynced the RNG and doubled regen BA.) Only HK<4.5 sets DBH directly
        # (regent.f:881), a tiny-tree edge that does not fire on realistic small-tree stands (H+HTG > 4.5).
        nrec = stash !== nothing ? 3 : 1
        central_dbh = d
        # regent.f carries the Fortran variable D ACROSS the tripling L-passes: D=DBH(I) once before label 918, then the
        # end-of-pass `D=DBH(K)` (regent.f:1000, reached only via the D<3 small-tree path) — so pass L uses the PREVIOUS
        # record's post-REGENT DBH (a central/copy whose HK<4.5 had DBH(K) set directly at :881) in XWT (:838), the
        # D≥3 gate (:861) and DDS/SQRT (:984-985); DBH(K)-based terms (DIAM floor, DGBND, BRATIO) keep the own DBH.
        dcar = d
        for l in 0:(nrec - 1)
            xwt_l = dcar <= xmn ? 0.0f0 : (dcar - xmn) / (xmx - xmn)
            small_l = dcar < 3.0f0
            zzran = 0.0f0
            if dgsd >= 1.0f0
                while true
                    zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                    (zzran <= 1.0f0 && zzran >= -1.5f0) && break
                end
            end
            htgr = htgr1 * fexp(zzran * IE_RG_HSIGMA)             # NIVAR multiplicative randomization (regent.f:924)
            l == 0 && _rhdump("RH2", s, i, sp, (h, wk3[i], htgr1, zzran, htgr))   # TEMP-DEBUG
            lhtg = l == 0 ? large_htg : (l == 1 ? large_htg_u : large_htg_l)  # per-copy large-tree HTG (htgf.f)
            htg = htgr * (1.0f0 - xwt_l) + xwt_l * lhtg             # blend toward the (per-copy) large-tree htgf value
            (h + htg > cap) && (htg = max(cap - h, 0.1f0))
            # diameter dub (D<3 only). dg_inc = inside-bark increment (added to DBH via GRADD); dbh_dir≥0 ⇒ set DBH.
            dg_inc = 0.0f0; dbh_dir = -1.0f0
            if small_l
                hk = h + htg
                if hk < 4.5f0
                    dbh_dir = 0.1f0 + diam * 0.01f0 + hk * 0.001f0          # regent.f:881 (DBH set, DG=0)
                else
                    dk = ax * fpow(hk - 4.5f0, bx) + dadj                        # regent.f:938
                    dk < diam && (dk = diam)
                    dk = dk + hk * 0.001f0
                    dgk = (dk - d1v) * xrdgro; dgk < 0.0f0 && (dgk = 0.0f0) # regent.f:958 (DK−D1)·XRDGRO
                    dg0 = dgk * bark                                        # DG(K)=DGK·BARK (regent.f:981)
                    dds = dg0 * (2.0f0*bark*dcar + dg0) * scale2            # regent.f:984 (carried D)
                    dg_inc = sqrt((dcar*bark)^2 + dds) - bark*dcar         # regent.f:985 (carried D)
                    (d + dg_inc) < diam && (dg_inc = diam - d)             # regent.f:987 DIAM floor
                    dg_inc = dg_bound(nothing, nothing, sp, d, dg_inc, s.control.sp_size_cap)  # DGBND (SIZCAP)
                end
            end
            if l == 0
                t.ht_growth[i] = htg
                if small_l
                    if dbh_dir >= 0.0f0
                        t.dbh[i] = dbh_dir; t.diam_growth[i] = 0.0f0; central_dbh = dbh_dir
                    else
                        t.diam_growth[i] = dg_inc                           # increment; DBH unchanged (=d)
                        prevbark = bark                                    # regent.f:978 recompute (NIVAR reaches the dub)
                        # regent.f:999-1021 — crown dub for a NIVAR small tree crossing 3″ with no crown.
                        # Inert for the usual FIA conifer regen (crown_ratio_update! keeps NIVAR crowns
                        # non-zero); present so a genuinely crown-less NIVAR seedling is not left at CR=0.
                        _ie_regent_dubscr_crown!(s, i, sp, d, h + htg, dg_inc, ntyr, yr, ba, dgsd)
                    end
                end
            elseif l == 1
                stash.htgU[i] = htg; stash.is_small[i] = true
                if small_l
                    stash.dgU[i] = dbh_dir >= 0.0f0 ? 0.0f0 : dg_inc
                    stash.dbhU[i] = dbh_dir >= 0.0f0 ? dbh_dir : d
                end
            else
                stash.htgL[i] = htg
                if small_l
                    stash.dgL[i] = dbh_dir >= 0.0f0 ? 0.0f0 : dg_inc
                    stash.dbhL[i] = dbh_dir >= 0.0f0 ? dbh_dir : d
                end
            end
            small_l && (dcar = dbh_dir >= 0.0f0 ? dbh_dir : d)   # regent.f:1000 D=DBH(K) (post-REGENT DBH of this record)
            # ↑ Each copy K carries its OWN DBH into UPDATE (TRIPLE copies the CENTRAL, whose DBH REGENT may already have
            # SET directly — HK<4.5 ⇒ regent.f:881): a directly-set copy gets DBH(K)=its value with DG(K)=0 (⇒ next WK1=0);
            # an increment-path copy keeps its PRE-growth DBH D and DG(K)=its increment, so UPDATE forms exactly
            # D+DG/BRATIO(D) (MEASURED live copy 512: 0.10091+0.3169/0.915=0.4472; WK1 of copy 257 = 0).
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
function ie_esgent!(s::StandState, nstart::Int; fint::Float32 = 10.0f0,
                    atavh::Float32 = -1.0f0, atba::Float32 = -1.0f0, atrelden::Float32 = -1.0f0,
                    relden_pre::Float32 = -1.0f0, ba_pre::Float32 = -1.0f0)
    p, t, c, dens = s.plot, s.trees, s.calib, s.density
    sd = s.coef.species
    nstart >= t.n && return s
    rhcon = ie_regcons!(s)
    # REGENT(LESTB) reads the START-of-cycle (post-thin, pre-growth) stand density as TEMBA/TEMCCF/TEMAHT
    # (regent.f:218-222 ← grincr.f ATBA/ATCCF/ATAVH; :326-328 AH=TEMAHT,R=TEMCCF; :1082-1083 BANEXT=TEMBA,
    # RDNEXT=TEMCCF), NOT the current post-growth stand. By establishment time p.avg_height reflects the GROWN
    # stand (top height ≫ 4.5), collapsing the DADJ relh clamp (relh<0⇒0 ⇒ DADJ=0.65 instead of the pre-growth
    # all-seedling relh=1 ⇒ DADJ=0.65−DELMAX≈0.77) so every crossing seedling's dubbed DBH is ~0.12" too small
    # (halving the birth-cohort BA). Use the captured start-of-cycle trio (avh/ba/relden); TEMBA/TEMCCF fall
    # back to BA/RELDEN when ≤0 (regent.f:219,221), TEMAHT=ATAVH is unconditional (0 on a bare stand ⇒ DADJ=0.65
    # both ways). Missing (-1) sentinel ⇒ the legacy post-growth scalars (defensive; the IE call site passes them).
    ba = atba > 0f0 ? atba : p.basal_area
    relden = atrelden > 0f0 ? atrelden : p.relative_density
    avh = atavh >= 0f0 ? atavh : p.avg_height; ah = avh
    # The REGENT height-growth terms (HTGRL BCCF·RDJ + BBAL·BAL, regent.f:500/442) use the CURRENT
    # (post-thin/post-fire) stand RELDEN/BA — RDNEXT/BANEXT = RELDEN/BA for the single LESTB subcycle
    # (NPER=1 at FINT=10) — NOT the start-of-cycle TEMCCF/TEMBA (which feed ONLY DELMAX/AH, regent.f:327-333).
    # Post-fire the start-of-cycle density is the pre-fire overstory (e.g. relden 372 vs current ~9), which
    # collapses HTGRL and stunts the birth cohort ~2.5-3× (the dominant IE simfire post-fire BA deficit).
    # atba/atrelden/atavh stay on delmax/relh (the birth-DBH DADJ dub, #194); only bal/htgrl move to current.
    # regent.f:288-295 (LESTB, NTYR=5 ⇒ one 5-yr period): RDNEXT(1)=RELDEN, BANEXT(1)=BA of the gradd.f:192 DENSE —
    # post-growth/post-disturbance but PRE-ESNUTR (no sprouts, no AUTOES/PLANT cohort). The caller snapshots them
    # before esuckr! (relden_pre/ba_pre); p.* here already include the new cohort (establish! recomputes density),
    # which fed a bare planted stand RDJ≈0.19 instead of 0 ⇒ every birth HTGRL ~0.12% low (live instrumented).
    ba_now = ba_pre >= 0f0 ? ba_pre : p.basal_area
    relden_now = relden_pre >= 0f0 ? relden_pre : p.relative_density
    # regent.f:276-294 (LESTB, label 8): BANEXT(1)=TEMBA+NYR*BAYR with BAYR=(BA-TEMBA)/ITOT, NYR=5 — and likewise
    # RDNEXT(1)=TEMCCF+NYR*CCFYR. Algebraically = BA/RELDEN, but the Float32 divide-then-multiply round trip is not the
    # identity (IE 4727120010690 cyc3: BAJ 417662A9 vs BA 417662A8 ⇒ birth HT 1 ULP). TEMBA/TEMCCF = ATBA/ATCCF
    # (start of cycle), falling back to BA/RELDEN when ≤0 (regent.f:218-221). ITOT = what the regent.f:205-211 KPER
    # loop leaves (NTYR for one subcycle, else the last KPER). LSKIPH ⇒ BAYR=CCFYR=0 (unused: no height growth).
    temba = atba > 0f0 ? atba : ba_now
    temccf = atrelden > 0f0 ? atrelden : relden_now
    ntyr_lestb = Int(round(fint)) - 5
    ba_htg = temba; relden_htg = temccf
    if ntyr_lestb > 0
        nper_l = ntyr_lestb ÷ Int(IE_RG_REGYR); (ntyr_lestb % Int(IE_RG_REGYR) != 0) && (nper_l += 1)
        itot_l = ntyr_lestb; nn_l = nper_l
        for _ in 1:nper_l
            nn_l == 1 && break
            kp = itot_l ÷ nn_l; itot_l -= kp; nn_l -= 1
        end
        bayr = (ba_now - temba) / Float32(itot_l)
        ccfyr = (relden_now - temccf) / Float32(itot_l)
        ba_htg = temba + 5f0 * bayr
        relden_htg = temccf + 5f0 * ccfyr
    end
    dgsd = s.control.dg_sd
    regyr = IE_RG_REGYR; yr = s.control.year
    ntyr = Int(round(fint))
    scale2 = ntyr > 0 ? yr / Float32(ntyr) : 1.0f0
    # REGENT(LESTB) grows the birth cohort over NTYR−5 years (regent.f:200; the 5-yr GENTIM lead is not grown).
    # Single 5-yr subcycle at FINT=10 ⇒ SCALE = KPER(1)/REGYR = 5/5 = 1. NTYR−5≤0 ⇒ LSKIPH (no height growth).
    ntyr_est = fint - 5.0f0
    est_scale = ntyr_est > 0.0f0 ? ntyr_est / regyr : 0.0f0
    lskiph = ntyr_est <= 0.0f0
    cur_year = current_cycle_year(s)
    delmax = (ah / 36.0f0) * (0.01232f0 * relden - 1.75f0); delmax > 0.0f0 && (delmax = 0.0f0)
    # ★ regent.f:301-320 "DUB CROWN RATIO FOR NEWLY ESTABLISHED SEEDLINGS" (the LESTB DO-13 pass). BEFORE the
    # height-growth ZZRAN loop, REGENT(.TRUE.) draws ONE main-stream BACHLO per NEW record (I=ITRNIN..ITRN,
    # ALL species, tree-index/storage order), rejection-bounded to |RAN|≤1, and dubs its crown from PCCF:
    #   CR = 0.89722 − 0.0000461·PCCF(ITRE(I)); 12 RAN=BACHLO(0,1,RANN); IF(|RAN|>1)GOTO 12; CR=CR+0.07985·RAN.
    # jl OMITTED this pass entirely ⇒ on an establishment cycle jl consumed ~5.3 fewer main-stream RANN per new
    # tree than the oracle (64 new trees ⇒ 342 missing draws on stand 1856050746290487 cyc-2), desyncing every
    # downstream DGSCOR/REGENT-ZZRAN draw from the FIRST establishment cycle onward — the dominant cause of the
    # "#206 ZZRAN straddle" (the divergence appears exactly one cycle after establishment fires). The draw is
    # RNG-independent of PCCF (BACHLO(0,1) regardless), so this re-aligns the stream bit-for-bit; the crown VALUE
    # (from point_ccf) is second-order. Unconditional (NOT dgsd-gated — regent.f has no DGSD guard on this draw).
    @inbounds for i in (nstart+1):t.n
        ipccf = Int(t.plot_id[i])
        pccf = (ipccf >= 1 && ipccf <= length(dens.point_ccf)) ? dens.point_ccf[ipccf] : 0f0
        cr = 0.89722f0 - 0.0000461f0 * pccf
        local ran::Float32
        while true
            ran = bachlo(s.rng, 0.0f0, 1.0f0)
            (ran >= -1.0f0 && ran <= 1.0f0) && break        # regent.f:310 IF(RAN<−1 .OR. RAN>1) redraw
        end
        cr += 0.07985f0 * ran
        cr > 0.90f0 && (cr = 0.90f0); cr < 0.20f0 && (cr = 0.20f0)
        t.crown_pct[i] = Int32(trunc(cr * 100.0f0 + 0.5f0))    # regent.f:314 ICR(I)=INT(CR*100+0.5)
    end
    # ★ REGENT(LESTB) DO-30/DO-25 iterates the NEW records SPECIES-MAJOR in IND1 (linked-list) order
    # (regent.f:695 DO 30 ISPC; :743 DO 25 I3=I1,I2 via IND1, gated IF(LESTB.AND.I.LT.ITRNIN)GO TO 25 to skip
    # the pre-existing trees). Each new record then draws its OWN per-record ZZRAN height error (regent.f:812).
    # Iterating the new records in STORAGE order (nstart+1:n) instead mis-PAIRS each drawn ZZRAN to the wrong
    # tree — the SAME wrong-pairing desync Bug P (7f284b29) fixed for the LESTB=F growth path (small_tree_growth!'s
    # _sp_order) but which was left un-fixed here in the birth-cycle (LESTB=T) path. The draw COUNT is unchanged
    # (same NIVAR set), so the main stream stays aligned; only the height↔tree assignment is corrected — which
    # de-scrambles the birth-cohort DBH pool that feeds the next cycle's AUTOES NSTORE tally. New-tree sort_key =
    # storage index (establishment.jl:2148), so species-major-with-increasing-index reproduces the oracle IND1.
    _es_order = Int[]
    @inbounds for sp in 1:MAXSP, i in (nstart+1):t.n
        Int(t.species[i]) == sp && push!(_es_order, i)
    end
    @inbounds for i in _es_order
        sp = Int(t.species[i]); d = t.dbh[i]
        d >= IE_RG_XMAX[sp] && continue
        t.tpa[i] <= 0.0f0 && continue
        (sp <= 12 || sp == 14 || sp == 23) || continue   # NIVAR conifers only (planted-conifer case)
        h = t.height[i]
        wk4 = t.htimlt[i]                                 # per-tree WK4=HTIMLT (PLANT/existing=1.0; AUTOES<1)
        con = rhcon[sp] + c.htg_cor_small[sp]             # CON = RHCON + HCOR
        # PCT(I) of a NEW record is 0 (estab.f:1252/1341) at REGENT(LESTB): FVS runs no DENSE between creating the record
        # and ESGENT (the post-establishment DENSE is gradd.f:244). jl's establish! already re-ran stand_pct! ⇒ read 0 here.
        pct = 0.0f0
        bal = ba_htg * (100.0f0 - pct) * 0.0001f0
        xrhgro = active_multiplier(s.control, :regh, sp, cur_year)
        xrdgro = active_multiplier(s.control, :regd, sp, cur_year)
        relh = abs(ah - 4.5f0) < 0.01f0 ? 0.0f0 : (h - 4.5f0) / (ah - 4.5f0)
        relh > 1.0f0 && (relh = 1.0f0); relh < 0.0f0 && (relh = 0.0f0)
        dadj = delmax*relh*relh - 2.0f0*delmax*relh + 0.65f0
        # REGENT(LESTB) NIVAR height increment (regent.f:596): H2=H1+EXP(HTGRL)·SCALE·XRHGRO·WK4 — WK4 enters HERE
        # (first application). LSKIPH ⇒ no height growth (regent.f:764 HTG=0).
        htgr1 = 0.0f0
        if !lskiph
            htgrl = con + IE_RG_RHLH[sp]*flog(h) + IE_RG_RHCCF[sp]*relden_htg + IE_RG_RHBAL[sp]*bal
            h2 = h + fexp(htgrl) * est_scale * xrhgro * wk4   # regent.f:596 EXP(HTGRL)*SCALE*XRHGRO*WK4 (left-assoc Float32)
            _rhdump("RH1", s, i, sp, (h, htgrl, fexp(htgrl), est_scale, xrhgro, wk4, h2, con)); _rhdump("RH0", s, i, sp, (relden_htg, ba_htg, bal))   # TEMP-DEBUG
            htgr1 = h2 - h; htgr1 < 0.0f0 && (htgr1 = 0.0f0)
        end
        xmn = IE_RG_XMIN[sp]; xmx = IE_RG_XMAX[sp]
        ax = IE_RG_HHT1[sp]; bx = IE_RG_HHT2[sp]
        xwt = d <= xmn ? 0.0f0 : (d - xmn) / (xmx - xmn)
        diam = IE_RG_DIAM[sp]; bark = ie_bratio(sp, d)
        d1v = diam + dadj; h > 4.5f0 && (d1v = ax * fpow(h - 4.5f0, bx) + dadj)
        let f = get(ENV, "JLRH", ""); isempty(f) || open(f, "a") do io; println(io, "RNG", lpad(Int(s.control.cycle)+1,3), lpad(i,5), " ", s.rng.s0); end; end   # TEMP-DEBUG
        zzran = 0.0f0
        if dgsd >= 1.0f0
            while true
                zzran = bachlo(s.rng, 0.0f0, 1.0f0)
                (zzran <= 1.0f0 && zzran >= -1.5f0) && break
            end
        end
        htgr = htgr1 * fexp(zzran * IE_RG_HSIGMA)
        _rhdump("RH2", s, i, sp, (h, h + htgr1, htgr1, zzran, htgr))   # TEMP-DEBUG
        htg_regent = htgr * (1.0f0 - xwt)                # REGENT HTG(I): new tree large-tree HTG(K)=0
        # ★ ESGENT (esgent.f:56-71): scale REGENT's HTG by WK4 AGAIN (second application ⇒ effective WK4²), then
        # for WK4<1 (AUTOES) reset the sub-breast-height DBH to the birth-cycle nominal 0.1+0.001·HT (DG=0) or,
        # once HT≥4.5, shrink REGENT's DBH/DG by the height ratio HT/HTEMP. WK4≥1 (PLANT) keeps REGENT's DBH dub.
        htemp = h + htg_regent
        htg = htg_regent * wk4
        cap = s.control.sp_size_cap[sp, 4]
        (h + htg > cap) && (htg = max(cap - h, 0.1f0))
        hk = h + htg
        t.height[i] = hk; t.ht_growth[i] = htg
        # REGENT(LESTB) NIVAR diameter: the ESTAB birth path sets DBH DIRECTLY to the absolute height→DBH
        # value DK — it does NOT accumulate a DDS increment (regent.f:938-946 `IF(LESTB) DBH(K)=DK; DG(K)=DK`).
        # DK is computed on the REGENT height HK = H + HTG(K), i.e. the ONE-WK4 height (esgent.f applies the
        # SECOND WK4 afterward), which for a new record is `htemp = h + htg_regent` here (XWT≡0 under LESTB).
        # The former code ran the NON-LESTB increment path (DGK=(DK−D1)·BARK on the DDS scale, DBH=d+DG): with
        # the tiny birth DBH d≈0.1 that gives DGK=(DK−D1)<0 (D1=DIAM+DADJ≈1.05 > DK≈0.9), so it FLOORED to
        # DBH=DIAM(=0.3) and the whole PLANT/AUTOES cohort entered the next cycle at ~0.3″ instead of DK≈0.8–1.0″
        # — a one-directional ~9% BA/QMD deficit that compounds every later cycle (the plant-regime residual).
        hk_reg = htemp                                   # REGENT HK (one WK4); == hk when WK4=1 (PLANT)
        if wk4 < 1.0f0
            if hk < 4.5f0                                # esgent.f:60 — HT (the DOUBLE-WK4 height) below breast height
                t.dbh[i] = 0.1f0 + 0.001f0 * hk; t.diam_growth[i] = 0.0f0   # esgent.f:61
            elseif d < 3.0f0
                # REGENT LESTB DBH(K)=DK on the REGENT height, then esgent HT/HTEMP shrink (esgent.f:64-65).
                dk = ax * fpow(hk_reg - 4.5f0, bx) + dadj; dk < diam && (dk = diam); dk = dk + hk_reg * 0.001f0
                ratio = htemp > 0.0f0 ? hk / htemp : 1.0f0
                t.dbh[i] = dk * ratio; t.diam_growth[i] = dk * ratio
            end
        elseif d < 3.0f0                                 # PLANT/existing (WK4≥1)
            if hk_reg < 4.5f0                            # regent.f:881-882 REGENT HK<4.5 ⇒ nominal sub-BH DBH
                t.dbh[i] = 0.1f0 + diam * 0.01f0 + hk_reg * 0.001f0; t.diam_growth[i] = 0.0f0
            else
                dk = ax * fpow(hk_reg - 4.5f0, bx) + dadj; dk < diam && (dk = diam); dk = dk + hk_reg * 0.001f0
                t.dbh[i] = dk; t.diam_growth[i] = dk    # regent.f:939,941 DBH(K)=DK; DG(K)=DK
            end
        end
    end
    return s
end

function _rhdump(tag, s, i, sp, vals)   # TEMP-DEBUG
    f = get(ENV, "JLRH", ""); isempty(f) && return
    h(x) = uppercase(string(reinterpret(UInt32, Float32(x)), base=16, pad=8))
    open(f, "a") do io; println(io, tag, lpad(Int(s.control.cycle)+1,3), lpad(i,5), lpad(Int(sp),3), " ", join(h.(vals), " ")); end
end
