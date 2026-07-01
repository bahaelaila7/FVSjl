# =============================================================================
# establishment.jl — regeneration / establishment (ESNUTR → ESTAB)
#
# Ported from: base/esnutr.f (cycle hook) + base/estab.f (tree creation) +
# base/estab_helpers.f (ESSUBH/ESTIME) + base/esinit.f.
#
# SN's PARTIAL (keyword-driven) establishment model: no auto-ingrowth. When an
# ESTAB packet scheduled PLANT(430)/NATURAL(431) activities are due, ESTAB creates
# the regen trees. A single bare plot (NPTIDS=1) is replicated MINREP=50 times: each
# replicate independently draws an established HEIGHT per species (ESSUBH height-at-
# age + a BACHLO draw on the establishment RNG), and contributes one record per
# species carrying plantedTPA·survival/100 / dupnpt TPA (400/50 = 8) ⇒ 50×2 = 100
# records, 800 TPA. The trees enter AFTER growth+mortality (GRADD order) so they're
# fresh (full TPA) this period; their DBH is derived from the established height.
# =============================================================================

const _ES_MINREP = 50          # MINREP: target plot replication (esinit.f)

# XMIN: per-species establishment min height (blkdat.f) lives in
# data/southern/species_coefficients.csv as the `estab_min_ht` column.
# HHTMAX: per-species max establishment height (blkdat.f).
const _ES_HHTMAX = Float32[23.0,27.0,21.0,21.0,22.0,20.0,24.0,18.0,18.0,17.0,22.0,
    (20.0 for _ in 12:90)...]

# NE HHTMAX: per-species MAX established height (ne/blkdat.f DATA HHTMAX/, 108 values). A HARD cap on the
# grown establishment height (not the soft site-curve HTMAX) — e.g. YB(30)=22, WO(55)=16 are reached exactly
# by live (all trees clamped). The SN _ES_HHTMAX above is wrong for NE (it uses a 20-ft fill for sp≥12).
const _NE_ES_HHTMAX = Float32[
    20,24,18,16,18,16,16,18,20,14, 14,16,16,16,16,16,16,16,14,14,
    16,18,12,20,16,20,16,16,18,22, 20,18,18,18,14,14,14,14,18,14,
    24,24,18,24,28,24,18,20,20,24, 24,20,20,26,16,14,12,12,16,16,
    14,16,14,16,16,12,20,16,16,14, 16,12,12,12,12,12,12,18,20,12,
    20,20,20,20,16,16,16,24,14,24, 32,18,16,16,16,16,12,10,16,18,
    30,20,20,18,16,20,20,30]

# NE ESSUBH per-species reference age CARAGE (essubh.f DATA MAPNE/, 108 values — DISTINCT from the htcalc
# curve-index MAPNE). The planted base height is (NC-128 height at this age / this age) · min(5, TIME−DELAY).
const _NE_ESSUBH_REFAGE = Int[
    20,10,15,20,15,20,20,20, 5,20, 15,20,20,10,20,20,20,20,20,10,
    10,15,20,15,10,20,20,20,20,20, 20,20,20,20,20,20,20,20,20,20,
    20,20,20,35,35,20,10,20,20,20, 10,20,15,20,10,10,10,10,10,10,
    10,30,10,10,30,30,20,10,10,20, 10,10,10,20,10,10,10,20,20,10,
    25,25,10,25,25,10,10,20,10,10, 10,10,20,20,20,20,20,10,10,10,
    10,10,10,10,10,10,10,10]

# CS ESSUBH per-species reference age CARAGE (cs/essubh.f DATA MAPCS/, 96 values — DISTINCT from the
# htcalc curve-index MAPCS). Same NE-style base height = (NC-128 height at CARAGE / CARAGE)·min(5, TIME−DELAY).
const _CS_ESSUBH_REFAGE = Int[
    10,10,10,15,20,10,5,20,20,25, 25,25,25,20,20,20,20,20,20,20,
    20,10,20,20,20,35,20,15,20,10, 15,20,20,20,10,20,20,20,20,20,
    20,20,20,20,20,35,10,20,10,10, 10,10,10,30,10,30,10,10,10,10,
    20,10,35,30,10,10,20,10,10,20, 20,10,10,20,20,20,10,20,20,20,
    20,10,10,10,10,10,10,10,10,20, 10,20,25,20,10,10]

# CS planted/regen height cap HHTMAX (cs/blkdat.f DATA HHTMAX/, 96 species). Clamps the REPORTED
# seedling height (estab.f:496 / esgent.f:64); the DBH is taken from the UNCAPPED grown height.
const _CS_ES_HHTMAX = Float32[
    16,27,14,14,14,16,20,20,18,16, 20,20,16,14,14,14,18,14,14,14,
    14,14,14,14,18,28,20,24,20,16, 18,26,16,14,12,20,16,20,12,20,
    24,16,16,24,24,24,16,20,16,16, 16,20,12,16,14,12,12,20,16,20,
    14,14,20,16,20,14,20,20,18,20, 20,12,20,24,20,20,24,20,24,20,
    18,18,20,32,10,20,20,18,16,20, 12,20,20,20,20,16]

"""
    establish!(state; fint=5f0) -> Bool

Create scheduled PLANT/NATURAL regen for the current cycle (ESNUTR/ESTAB). Runs at
the end of `grow_cycle!` (GRADD order). Idempotent per year. Returns whether any
tree was created. No-op unless an ESTAB packet is active.
"""
function establish!(s::StandState; fint::Float32 = 5f0)::Bool
    s.estab.active || return false
    t = s.trees; sd = s.coef.species
    es_xmin = sd[:estab_min_ht]   # per-species establishment min height (CSV)
    es_hhtmax = s.variant isa Northeast ? _NE_ES_HHTMAX :
                s.variant isa CentralStates ? _CS_ES_HHTMAX : _ES_HHTMAX   # per-variant HHTMAX (base + grown caps)
    per = round(Int, fint)
    yr = Int32(current_cycle_year(s))   # IY schedule; yr+per below = next boundary (fint is per-cycle)
    yr in s.estab.years_done && return false
    due = [a for a in s.control.schedule
           if (a.icflag == Int32(430) || a.icflag == Int32(431)) && yr <= a.year < yr + per]
    isempty(due) && return false

    # NPTIDS = IPTINV − NONSTK (esplt2.f:74): the STOCKABLE inventory points, not the raw
    # plot count. Driving DUPNPT/IDUP and so the regen record count + its per-record RNG draws.
    nptids = max(1, Int(s.plot.points_inv) - Int(s.plot.nonstockable))
    # estab.f:199-207: IDUP = smallest I with NPTIDS·I ≥ MINREP = CEIL(MINREP/NPTIDS) (not floor); the
    # MAXPLT cap doesn't bind for the divergent 1<NPTIDS<MINREP cases. NPTIDS=1 ⇒ ceil=floor=50 (BARE stand).
    idup   = max(1, cld(_ES_MINREP, nptids))
    dupnpt = Float32(nptids * idup)
    # ESSUBH base height from age uses the variant's site-curve: SN Chapman-Richards (ht_curve_b*),
    # NE NC-128 (ne_htcalc_height). bc is SN-only (NE has no ht_curve_b* coefs).
    bc = (s.variant isa Northeast || s.variant isa CentralStates) ? nothing :
         (sd[:ht_curve_b1], sd[:ht_curve_b2], sd[:ht_curve_b3], sd[:ht_curve_b4], sd[:ht_curve_b5])
    montane = !isempty(s.plot.eco_unit) && s.plot.eco_unit[1] == 'M'
    ifor = Int(s.plot.forest_idx)
    # Natural-height random-draw acceptance window (estab.f:482-483 SN/CS vs :489-490 NE). FVS draws
    # RAN~N(0.5,0.25) and REDRAWS until RAN falls in the window; the window is VARIANT-SPECIFIC:
    # NE accepts [-2.5, 2.5]; SN and CS accept [0.0, 1.5]. The narrower SN/CS window truncates the
    # tails (no negative RAN ⇒ fewer trees pinned to the XMIN floor, a longer upper tail) AND rejects
    # more draws, so it also changes how many :estab draws each replicate consumes — a different window
    # desyncs the whole establishment RNG stream (and, downstream, the shared small-tree growth RANN
    # stream), which is the bare_natural sawtimber-tail divergence (D10). jl previously hardcoded the
    # NE window on the shared path.
    ran_lo, ran_hi = s.variant isa Northeast ? (-2.5f0, 2.5f0) : (0f0, 1.5f0)
    # gentim/delay/trage timing (esnutr/estab/essubh): age = FINT − delay − gentim + trage.
    # estab.f:448-449 — GENTIM = FINT−5 (clamped ≥0), depends ONLY on FINT, never IDSDAT/calendar
    # year. (Was `yr − idsdat`, a confirmed bandaid B5; masked today by the es_xmin height floor.)
    gentim = max(per - 5, 0)
    # Each new regen tree's crown ratio uses the per-point CCF computed by DENSE from the EXISTING (pre-regen)
    # overstory: regent.f:178 `CR=0.89722−0.0000461·PCCF(IPCCF)` with `IPCCF=ITRE(I)` (the tree's point). We now
    # carry that exact per-point value (`density.point_ccf`, filled by `point_density!` at start-of-cycle) and
    # index it by each record's point below — replacing the prior whole-stand `stand_ccf` approximation. The
    # coefficient is tiny (4.6e-5), so a bare/sparse stand (CCF≈0) is unchanged to print resolution.
    created = false
    nstart = t.n        # tree count before establishment (phase-2 crown pass starts here)
    # REGENT-LESTB's BALMOD competition uses the PRE-establishment density — the new seedlings do NOT compete
    # in their own creation cycle (live FVSne debug: GMOD=1.0 / AVH=0 for a BARE stand; the DENSE/BAL the cycle
    # uses predates the regen). Snapshot the BAL over the existing overstory (1:nstart) NOW, before any seedling
    # is added; computing it AFTER (over the cohort, the old code) over-counted the seedlings' own BA and
    # under-grew the established cohort ~4% (dbh 1.12 vs live 1.17 ⇒ the cyc-1 SDI/CCF deficit).
    ebau_pre = zeros(Float32, 50)
    s.variant isa Northeast && ne_badist!(ebau_pre, s)
    # estab.f:175-205 — pre-replicate :estab RNG setup, consumed BEFORE any per-replicate height draw.
    # On the FIRST tally (NTALLY==1) FVS draws once to derive ESDRAW = INT(DRAW·1e5+0.5) and SAVEs it;
    # every tally then reseeds the establishment stream with ESRNSD(.TRUE.,ESDRAW) (odd-forced) and
    # consumes IDUP·NPTIDS draws filling the WK6 site-prep vector. Both advance the :estab stream ahead
    # of the height draws, so jl MUST consume them or every replicate's BACHLO height is off — which
    # (D10) shifts the sp3 seedling sizes, hence the cycle each crosses 3" DBH into the large-tree DGF,
    # which desyncs the sp13 DGSCOR serial-correlation stream and spreads the sawtimber tail.
    s.estab.ntally += Int32(1)
    if s.estab.ntally == Int32(1)
        s.estab.es_seed = floor(esrann!(s.rng) * 100000f0 + 0.5f0)   # fresh ESDRAW (NTALLY==1)
    end
    esd = s.estab.es_seed
    (esd % 2f0 == 0f0) && (esd += 1f0)                               # ESRNSD odd-force (esrann.f:56)
    s.rng.es0 = Float64(esd)
    for _ in 1:(nptids * idup)                                       # WK6 site-prep fill (estab.f:202-205)
        esrann!(s.rng)
    end
    # estab.f outer loop: `for nn in 1:NPTIDS` (each inventory point) × `idup` replicates
    # → NPTIDS·idup records total. For a BARE stand every point is identical (BAAA=0,
    # uniform slope/aspect/habitat), so the per-point variables don't vary; only the
    # record count and the ESRANN draw count scale with NPTIDS. (ptree already divides by
    # dupnpt = NPTIDS·idup, so the planted TPA is conserved across all the records.)
    @inbounds for nn in 1:nptids, rep in 1:idup
        # per-replicate establishment RNG draws (estab.f:216-221): two for emsqr
        # (unused on the no-treeht path), one for esdraw (the re-seed value).
        esrann!(s.rng); esrann!(s.rng)
        esdraw = floor(esrann!(s.rng) * 100000f0 + 0.5f0)
        for a in due
            sp = round(Int, a.params[1]); (1 <= sp <= MAXSP) || continue
            ptree = a.params[2] * (a.params[3] / 100f0) / dupnpt
            ptree <= 0f0 && continue
            delay  = Int(a.year) - Int(yr)
            trage  = a.params[4] < 0.5f0 ? 2f0 : a.params[4]; trage > 10f0 && (trage = 10f0)
            age = Float32(per) - Float32(delay) - Float32(gentim) + trage; age < 1f0 && (age = 1f0)
            si  = s.plot.sp_site_index[sp]
            # ESSUBH base height (essubh.f:73-82). NE uses its OWN formula — NOT the site-curve height at the
            # tree age: a per-species reference age CARAGE (essubh.f MAPNE, distinct from the htcalc curve map),
            # H = NC-128 site-curve height at CARAGE, then HHT = (H/CARAGE)·min(5, TIME−DELAY) (avg juvenile rate
            # × available time). The `age` above is FVS's REGENT-start AGE (essubh.f:93), used by growth, not the
            # planted height. SN keeps the Curtis-Arney htcalc_height(age).
            hht = if s.variant isa Northeast
                carage = Float32(_NE_ESSUBH_REFAGE[sp])
                (ne_htcalc_height(sp, si, carage) / carage) * min(5f0, Float32(per) - Float32(delay))
            elseif s.variant isa CentralStates
                # CS ESSUBH (cs/essubh.f:72-81): identical NE-style base height with the CS refage map +
                # the CS NC-128 forward curve — H at CARAGE, then (H/CARAGE)·min(5, TIME−DELAY).
                carage = Float32(_CS_ESSUBH_REFAGE[sp])
                (cs_htcalc_height(sp, si, carage) / carage) * min(5f0, Float32(per) - Float32(delay))
            else
                htcalc_height(bc, sp, si, age, montane)
            end
            treeht = a.params[5]
            if treeht >= 0.1f0                                      # PLANT specified a height
                hht = treeht; xh = log(hht)
                while true
                    xxh = exp(bachlo(s.rng, xh, 0.5f0; stream = :estab))
                    (0.5f0 * hht <= xxh <= 2f0 * hht) && (hht = xxh; break)
                end
                hht < 0.05f0 && (hht = 0.05f0)                      # PLANT floor 0.05 (estab.f:1034), HTADJ=0
            else                                                   # default: RAN~N(0.5,0.25), accept RAN∈[ran_lo,ran_hi]
                while true
                    ran = bachlo(s.rng, 0.5f0, 0.25f0; stream = :estab)
                    (ran_lo <= ran <= ran_hi) && (hht += ran; break)  # estab.f:483/490 (variant-specific window)
                end
                hht < es_xmin[sp] && (hht = es_xmin[sp])           # default/natural floor XMIN (estab.f:1037)
            end
            hht > es_hhtmax[sp] && (hht = es_hhtmax[sp])
            ibrkup = floor(Int, ptree / 10f0 + 1f0); brk = Float32(ibrkup)
            # Establishment DBH from the grown seedling height (esgent.f:55-62). A seedling still BELOW
            # breast height (HT < 4.5 ft) has no real DBH — FVS assigns the nominal `DBH = 0.1 + 0.001·HT`
            # (esgent.f:56), NOT the HTDBH⁻¹ inverse. jl previously ran HTDBH⁻¹ for every seedling, which
            # over-sized sub-breast-height regen (bare_natural: DBH 0.225 vs live 0.10 at HT~3.4 ft),
            # inflating stand BA ~0.26% and biasing large-tree DGF growth (D10). Only HT ≥ 4.5 uses the
            # inverse, floored to the species min DIAM + the height-proportional add.
            if hht < 4.5f0
                dbh = 0.1f0 + 0.001f0 * hht
            else
                dbh = _htdbh_dbh(sd, sp, hht, ifor); dbh < 0.1f0 && (dbh = 0.1f0)
                dbh += 0.001f0 * hht
            end
            for _ in 1:ibrkup
                n = t.n + 1; n > length(t.dbh) && break
                t.n = n
                t.species[n]     = Int32(sp)
                t.dbh[n]         = dbh
                t.height[n]      = hht
                t.tpa[n]         = ptree / brk
                # Records go on inventory point `nn` (estab.f:313 ITRE=IPTIDS[nn]).
                # point_ba scales each point's raw BA by PI/GROSPC with PI=NPTIDS, so with
                # the planted TPA spread evenly over NPTIDS points each point_ba comes back
                # to the full stand BA — matching the oracle's pba=ba_v fallback (PTBAA≤0)
                # for fresh establishment, for any NPTIDS (NPTIDS=1 ⇒ this is point 1).
                t.plot_id[n]     = Int32(nn)
                t.crown_pct[n]   = Int32(0)            # crown set in phase 2 (REGENT lestb)
                t.crown_ratio[n] = 0f0
                t.norm_ht[n]     = Int32(0)
                t.sort_key[n]    = Float64(n)
                created = true
            end
        end
        es = esdraw; (es % 2f0 == 0f0) && (es += 1f0); s.rng.es0 = Float64(es)  # ESRNSD(true,esdraw)
    end
    # PHASE 2 — ESGENT → REGENT(lestb): assign each new tree its open-grown crown in
    # SPESRT (species-then-record) order (regent.f:107-116). cr = 0.89722 −
    # 0.0000461·PCCF + 0.07985·N(0,1)[±1], clamp [0.20,0.90]; the crown draw uses the
    # MAIN RANN stream (separate from the ESRANN heights). The per-cycle CROWN
    # (crown_ratio_update!, run after) then applies its ±1%/yr change limit (~85).
    if created
        newidx = sort(collect((nstart + 1):t.n); by = i -> (Int(t.species[i]), i))
        # NE only: REGENT(LESTB) also GROWS each new seedling its creation cycle (esgent.f:48). SN's
        # essubh assigns the full height-at-age directly, so SN needs no growth here; NE's essubh gives
        # a BASE height that this grows to the cycle-end height (the BARE-stand TopHt fix). XWT=0 for LESTB.
        # CS shares NE's REGENT(LESTB) shape: ESSUBH gives a BASE height that this grows to the cycle-end
        # height, with the CS NC-128 increment + CS balmod (cs/regent.f:118-340, FNT=FINT−5).
        ne_estab = s.variant isa Northeast
        cs_estab = s.variant isa CentralStates
        local ebau_e, b3_e, avh_e, scale_e, rdiam_e, rnd_e
        local cb1_e, cb2_e, cb3_e, ba_e
        if ne_estab
            ebau_e = ebau_pre                              # PRE-establishment BAL (snapshot above), not the cohort's
            b3_e = sd[:dg_b3]; avh_e = s.plot.avg_height
            # REGENT LESTB period: FNT = FINT−5 (regent.f:118-124; LSKIPH ⇒ no ht growth when FINT≤5).
            scale_e = per > 5 ? Float32(per - 5) / NE_REGENT_REGYR : 0f0   # CON=HGADJ=XRHGRO=1
            rdiam_e = sd[:regent_min_diam]
            rnd_e = s.control.dg_stddev_bound >= 1f0        # DGSD random ±10%
        elseif cs_estab
            cb1_e = sd[:balmod_b1]; cb2_e = sd[:balmod_b2]; cb3_e = sd[:balmod_b3]
            avh_e = s.plot.avg_height
            ba_e = s.plot.basal_area                       # PRE-establishment overstory BA (seedlings not yet in)
            scale_e = per > 5 ? Float32(per - 5) / 10f0 : 0f0   # FNT/REGYR, REGYR=10 (CON=HGADJ=XRHGRO=1)
            rdiam_e = sd[:htdbh_db]                         # DIAM floor (= cs/regent.f DIAM == htdbh_db)
            rnd_e = s.control.dg_stddev_bound >= 1f0        # DGSD random ±10%
        end
        @inbounds for i in newidx
            ran_cr = 0f0
            while true
                ran_cr = bachlo(s.rng, 0f0, 1f0)
                -1f0 <= ran_cr <= 1f0 && break
            end
            pccf = s.density.point_ccf[Int(t.plot_id[i])]      # PCCF(IPCCF), IPCCF=ITRE(I) (regent.f:160,178)
            cr = clamp(0.89722f0 - 0.0000461f0 * pccf + 0.07985f0 * ran_cr, 0.20f0, 0.90f0)
            icr0 = floor(Int32, cr * 100f0 + 0.5f0)
            t.crown_pct[i]   = icr0
            t.crown_ratio[i] = Float32(icr0)
            if ne_estab                                        # REGENT(LESTB) height growth + new DBH
                sp = Int(t.species[i]); h = t.height[i]; si = s.plot.sp_site_index[sp]
                if ne_htcalc_htmax(sp, si) - h <= 1f0
                    htgr = 0.1f0
                else
                    htgr = ne_htcalc_incr(sp, si, ne_htcalc_age(sp, si, h)) * scale_e
                end
                gmod = ne_balmod(b3_e[sp], ebau_e, t.dbh[i])
                relht = avh_e > 0f0 ? min(h / avh_e, 1f0) : 0f0
                htgr = max(htgr * (1f0 - (1f0 - gmod) * (1f0 - relht)), 0.1f0)
                if rnd_e
                    rh = 0f0
                    while true; rh = bachlo(s.rng, 0f0, 1f0); -1f0 <= rh <= 1f0 && break; end
                    htgr = max(htgr + rh * 0.1f0 * htgr, 0.1f0)
                end
                hk = h + htgr
                # DBH is derived from the UNCAPPED grown height (the HHTMAX clamp below only bounds the
                # REPORTED height, not the diameter — live YB: dbh from the grown ~23.5 ⇒ 1.8, height clamped
                # to HHTMAX 22). Computing dbh from the clamped height under-sized it (SDI/CCF dropped).
                if hk <= 4.5f0                       # regent.f:290-293: DG=0, DBH=D+0.001·HK (no Wykoff inverse)
                    t.dbh[i] = t.dbh[i] + 0.001f0 * hk
                else
                    dnew = _htdbh_dbh(sd, sp, hk, ifor); dnew < 0.1f0 && (dnew = 0.1f0)
                    dnew < rdiam_e[sp] && (dnew = rdiam_e[sp])
                    t.dbh[i] = dnew + 0.001f0 * hk
                end
                hk > _NE_ES_HHTMAX[sp] && (hk = _NE_ES_HHTMAX[sp])   # HARD HHTMAX clamp on the REPORTED height
                t.height[i] = hk
            elseif cs_estab                                    # CS REGENT(LESTB) height growth + new DBH
                sp = Int(t.species[i]); h = t.height[i]; si = s.plot.sp_site_index[sp]
                if cs_htcalc_htmax(sp, si) - h <= 1f0          # HTMAX−H ≤ 1 ⇒ HTG=0.1 (cs/regent.f:206)
                    htgr = 0.1f0
                else
                    htgr = cs_htcalc_incr(sp, si, cs_htcalc_age(sp, si, h)) * scale_e
                end
                bal = (1f0 - t.crown_ratio[i] / 100f0) * ba_e  # point BAL ≈ overstory BA for the smallest tree
                gmod = cs_balmod(cb1_e[sp], cb2_e[sp], cb3_e[sp], bal, ba_e, t.dbh[i])
                relht = avh_e > 0f0 ? min(h / avh_e, 1f0) : 0f0
                htgr = max(htgr * (1f0 - (1f0 - gmod) * (1f0 - relht)), 0.1f0)
                if rnd_e
                    rh = 0f0
                    while true; rh = bachlo(s.rng, 0f0, 1f0); -1f0 <= rh <= 1f0 && break; end
                    htgr = max(htgr + rh * 0.1f0 * htgr, 0.1f0)
                end
                hk = h + htgr
                # CS LESTB dbh (cs/regent.f:338-341): DBH = htdbh⁻¹(hk), floored to DIAM (or DIAM if hk<4.5),
                # THEN + 0.001·hk. (The htdbh inverse needs hk>4.5; for hk<4.5 the DIAM floor applies first.)
                if hk < 4.5f0
                    dbhk = rdiam_e[sp]
                else
                    dbhk = _htdbh_dbh(sd, sp, hk, ifor)
                    dbhk < rdiam_e[sp] && (dbhk = rdiam_e[sp])
                end
                t.dbh[i] = dbhk + 0.001f0 * hk
                hk > _CS_ES_HHTMAX[sp] && (hk = _CS_ES_HHTMAX[sp])   # HARD HHTMAX clamp on the REPORTED height
                t.height[i] = hk
            end
        end
        # ESGENT calls SPESRT to RE-ESTABLISH the species-order sort after adding
        # regen (esgent.f:41-44). SPESRT/LNKCHN visit records in ascending-record
        # order, so reset the lineage key to the physical record position: otherwise
        # stale TRIPLE lineage keys (3·K+offset) from earlier cycles, which are never
        # reconciled without a thinning compaction, scramble the post-establishment
        # species_sort! order and desync the per-tree DGSCOR RNG stream from FVS.
        @inbounds for i in 1:t.n
            t.sort_key[i] = Float64(i)
        end
        compute_density!(s)
    end
    push!(s.estab.years_done, yr)
    return created
end
