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

const _ES_MINREP = 50          # MINREP: target plot replication (esinit.f) — the DEFAULT for
                               # Establishment.minrep (state.jl); the live value is per-stand via the
                               # MINPLOTS keyword (esin.f opt 20), read at both idup call sites as est.minrep.

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

# LS planted/regen height cap HHTMAX (ls/blkdat.f DATA HHTMAX/, 68 species) + ESSUBH reference age CARAGE
# (ls/essubh.f DATA MAPLS/, 68 species — DISTINCT from the htcalc MAPLS curve map).
const _LS_ES_HHTMAX = Float32[14,20,18,18,20,18,18,20,16,24,16,16,16,16,18,24,24,18,20,26,16,12,20,22,16,16,16,14,24,16,16,14,12,20,16,20,20,14,14,20,20,24,18,20,18,20,20,24,10,16,18,20,20,20,12,18,16,20,16,24,30,20,20,20,32,20,18,20]
const _LS_ESSUBH_REFAGE = Int[20,15,20,20,5,15,15,20,20,10,20,20,10,10,20,35,15,15,20,20,20,20,20,20,20,20,20,20,20,10,30,10,10,20,10,10,20,20,20,20,20,20,20,20,20,20,10,10,10,10,10,10,10,20,10,10,10,10,25,20,10,10,10,10,10,10,10,10]

# Establishment min-height (XMIN) + max seedling height (HHTMAX) per species — from each variant's
# blkdat.f (VERIFIED: IE blkdat.f:62 XMIN == _IE_ES_XMIN). EM/BM/UT/CI had no establishment.jl ⇒ the dispatch
# fell to the missing `:estab_min_ht` coef ⇒ KeyError crash on ESTAB/PLANT-keyword stands (full utt01/emt01/
# bmt01). NOTE: blkdat.f XMIN (establishment), NOT regent.f XMIN (small-tree/regen — a DIFFERENT array).
const _EM_ES_XMIN   = Float32[1,1,1,1,0.5,0.5,1,0.5,0.5,1,3,6,3,3,3,3,6,0.5,3]
const _EM_ES_HHTMAX = Float32[23,27,21,27,18,6,24,18,18,17,16,16,16,16,16,16,16,22,16]
const _BM_ES_XMIN   = Float32[0.9,1.7,1,1,0.5,0.5,1.3,0.5,0.5,1,1,1,1,1,6,1,1,1]
const _BM_ES_HHTMAX = Float32[23,27,21,21,22,6,24,18,18,17,23,9,20,20,16,20,17,20]
const _UT_ES_XMIN   = Float32[1,1,1,0.5,0.5,6,1,0.5,0.5,1,0.5,0.5,0.5,0.5,0.5,0.5,0.5,3,3,0.5,0.5,3,0.5,0.5]
const _UT_ES_HHTMAX = Float32[9,9,10,7,7,16,10,7,7,10,6,6,10,6,6,6,9,16,16,6,6,16,9,10]
const _CI_ES_XMIN   = Float32[1,1,1,0.5,0.5,0.5,1,0.5,0.5,1,1,1,6,0.5,0.5,1,3,0.5,3]
const _CI_ES_HHTMAX = Float32[23,27,21,21,22,20,24,18,18,17,27,27,16,6,6,27,16,22,16]
# EC (East Cascades, MAXSP=32) establishment min height XMIN / max seedling height HHTMAX — ec/blkdat.f:140,153.
const _EC_ES_XMIN   = Float32[1,1,1,0.5,0.5,0.5,1,0.5,0.5,1,1,0.5,1,1,1,0.5,1.5,1,1,1,1,1,1,1,1,1,1,1,1,1,0.5,1]
const _EC_ES_HHTMAX = Float32[23,27,21,21,22,20,24,18,18,17,20,22,20,20,20,20,20,20,20,20,20,50,20,20,20,20,20,20,20,20,22,20]
# NC (Klamath) establishment per-species min height (XMIN) / max sprout height (HHTMAX) — nc/blkdat.f:73,80.
const _NC_ES_XMIN   = Float32[1,1,1,0.5,1,0.5,0.5,1,0.5,1,1,1]
const _NC_ES_HHTMAX = Float32[27,31,25,25,26,24,28,20,20,18,26,25]
# NC subsequent/planted base height (nc/essubh.f): a FIXED per-species table (no age/site/EMSQR),
# clamped [XMIN,HHTMAX] by the shared engine (like UT/TT).
const _NC_ESSUBH_HHT = Float32[1,1,1,1,7,1,7,7,1,0.8,7,2]

"""
    establish!(state; fint=5f0) -> Bool

Create scheduled PLANT/NATURAL regen for the current cycle (ESNUTR/ESTAB). Runs at
the end of `grow_cycle!` (GRADD order). Idempotent per year. Returns whether any
tree was created. No-op unless an ESTAB packet is active.
"""
function establish!(s::StandState; fint::Float32 = 5f0)::Bool
    s.estab.active || return false
    t = s.trees; sd = s.coef.species
    es_xmin = s.variant isa CentralRockies ? _CR_ES_XMIN :
              s.variant isa InlandEmpire ? _IE_ES_XMIN :
              s.variant isa Teton ? _TT_ES_XMIN :
              s.variant isa EasternMontana ? _EM_ES_XMIN :
              s.variant isa BlueMountains ? _BM_ES_XMIN :
              s.variant isa Utah ? _UT_ES_XMIN :
              s.variant isa CentralIdaho ? _CI_ES_XMIN :
              s.variant isa EastCascades ? _EC_ES_XMIN :
              s.variant isa Klamath ? _NC_ES_XMIN :
              sd[:estab_min_ht]   # per-species establishment min height (eastern SN/NE/CS/LS have this column)
    es_hhtmax = s.variant isa Northeast ? _NE_ES_HHTMAX :
                s.variant isa CentralStates ? _CS_ES_HHTMAX :
                s.variant isa LakeStates ? _LS_ES_HHTMAX :
                s.variant isa CentralRockies ? _CR_ES_HHTMAX :
                s.variant isa InlandEmpire ? _IE_ES_HHTMAX :
                s.variant isa Teton ? _TT_ES_HHTMAX :
                s.variant isa EasternMontana ? _EM_ES_HHTMAX :
                s.variant isa BlueMountains ? _BM_ES_HHTMAX :
                s.variant isa Utah ? _UT_ES_HHTMAX :
                s.variant isa CentralIdaho ? _CI_ES_HHTMAX :
                s.variant isa EastCascades ? _EC_ES_HHTMAX :
                s.variant isa Klamath ? _NC_ES_HHTMAX : _ES_HHTMAX   # per-variant HHTMAX (base + grown caps)
    per = round(Int, fint)
    yr = Int32(current_cycle_year(s))   # IY schedule; yr+per below = next boundary (fint is per-cycle)
    yr in s.estab.years_done && return false
    # PLANT/NATURAL dates < 1000 are CYCLE NUMBERS (FVS 1-based), not calendar years — the same OPNEW/OPFIND
    # convention cuts! applies (cuts.jl:203-208). Without the cycle-number clause a `PLANT 2 ...` (cycle 2) never
    # matched `yr <= a.year` (2016 <= 2 is false) ⇒ planting silently never fired (bit vs live only on real FIA
    # stands scheduled by cycle; thin/salvage already resolved this way, ESTAB was the omission).
    fvscyc = Int(s.control.cycle) + 1
    due = [a for a in s.control.schedule
           if (a.icflag == Int32(430) || a.icflag == Int32(431)) &&
              ((yr <= a.year < yr + per) || (0 < Int(a.year) < 1000 && Int(a.year) == fvscyc))]
    isempty(due) && return false

    # NPTIDS = IPTINV − NONSTK (esplt2.f:74): the STOCKABLE inventory points, not the raw
    # plot count. Driving DUPNPT/IDUP and so the regen record count + its per-record RNG draws.
    nptids = max(1, Int(s.plot.points_inv) - Int(s.plot.nonstockable))
    # estab.f:199-207: IDUP = smallest I with NPTIDS·I ≥ MINREP = CEIL(MINREP/NPTIDS) (not floor); the
    # MAXPLT cap doesn't bind for the divergent 1<NPTIDS<MINREP cases. NPTIDS=1 ⇒ ceil=floor=50 (BARE stand).
    idup   = max(1, cld(Int(s.estab.minrep), nptids))   # MINPLOTS keyword (esin.f MINREP; default 50)
    dupnpt = Float32(nptids * idup)
    # ESSUBH base height from age uses the variant's site-curve: SN Chapman-Richards (ht_curve_b*),
    # NE NC-128 (ne_htcalc_height). bc is SN-only (NE has no ht_curve_b* coefs).
    bc = (s.variant isa Northeast || s.variant isa CentralStates || s.variant isa LakeStates ||
          s.variant isa CentralRockies || s.variant isa InlandEmpire || s.variant isa Teton ||
          s.variant isa EasternMontana || s.variant isa BlueMountains || s.variant isa Utah ||
          s.variant isa CentralIdaho || s.variant isa EastCascades ||
          s.variant isa Klamath) ? nothing :   # western variants use a fixed/XMIN base, not the SN ht-curve
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
    # Establishment default-height RAN acceptance window (estab.f:483/490): SN = [0,1.5] (sn/estab.f:486);
    # NE, CS, AND LS all = [-2.5,2.5] (ne/cs/ls estab.f:490). The old `Northeast ? … : (0,1.5)` wrongly gave
    # CS AND LS the SN window [0,1.5], which REJECTS the low tail (RAN<0) ⇒ biased the planted-seedling
    # heights HIGH (esp. the smallest, whose small-RAN draws live accepts) — the BARE-PLANT over-sizing.
    ran_lo, ran_hi = (s.variant isa Southern || s.variant isa CentralRockies || s.variant isa InlandEmpire || s.variant isa Teton || s.variant isa EastCascades) ? (0f0, 1.5f0) : (-2.5f0, 2.5f0)   # CR/IE/TT/EC = SN window (cr/estab.f:486; ec/estab.f:486 RAN∈[0,1.5])
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
    # IPTIDS (esplt2.f:77-131): the STOCKABLE inventory point indices, indexed by the estab outer loop
    # (estab.f:313 `ITRE=IPTIDS[nn]`) — NOT the raw loop counter nn. A nonstockable plot (its `mort_code==8`
    # ".tre" record is skipped in treeinput.jl:91) has NO stored tree, so it is absent from the overstory
    # plot_ids; the stockable points are exactly the distinct plot_ids that DO carry a record. Using nn
    # directly put regen on the nonstockable point and skipped a stockable one, reading the wrong
    # `density.point_ccf[plot_id]` ⇒ the estab_pccf 7-tree crown residual (plant_stocked point 7 nonstockable).
    # FALLBACK to nn when the count doesn't match NPTIDS (bare stands: no overstory ⇒ empty ⇒ identity nn),
    # so every no-nonstockable-point scenario stays bit-exact.
    iptids = sort(unique(Int(t.plot_id[i]) for i in 1:nstart))
    use_iptids = length(iptids) == nptids
    # REGENT-LESTB's BALMOD competition uses the PRE-establishment density — the new seedlings do NOT compete
    # in their own creation cycle (live FVSne debug: GMOD=1.0 / AVH=0 for a BARE stand; the DENSE/BAL the cycle
    # uses predates the regen). Snapshot the BAL over the existing overstory (1:nstart) NOW, before any seedling
    # is added; computing it AFTER (over the cohort, the old code) over-counted the seedlings' own BA and
    # under-grew the established cohort ~4% (dbh 1.12 vs live 1.17 ⇒ the cyc-1 SDI/CCF deficit).
    ebau_pre = zeros(Float32, 50)
    s.variant isa Northeast && ne_badist!(ebau_pre, s)
    # LS ls_balmod reads RMSQD from the DENSE common, which for the establishment cohort is the
    # PRE-establishment stand QMD (live FMEFF stamp: BARE stand → RMSQD=0, so ls_balmod takes the
    # rmsqd≤0 → omega=b4 branch, GM 0.745 for jack pine). Snapshot it BEFORE the seedlings are added;
    # `stand_qmd(s)` recomputed after would include the cohort (0.626) and flip ls_balmod to the else
    # branch (GM 1.0) ⇒ the cohort over-grows (the BARE-PLANT seedling over-sizing).
    rmsqd_pre = stand_qmd(s)
    # CS REGENT-LESTB BALMOD needs the POST-growth OVERSTORY BA/AVH (pre-seedling) — same snapshot pattern as
    # ebau_pre/rmsqd_pre above. plot.basal_area/avg_height are STALE here (set pre-growth at cycle start;
    # compute_density! refreshes them only AFTER establish!, simulate.jl:463), and stand_ba/stand_top_height
    # computed later in the phase-2 loop would over-count the new seedlings' own BA (breaks the BARE-GROUND case).
    # Recompute NOW over the overstory (1:nstart, no seedlings yet): cs_estab overstory 134.1/70.2; bare stand 0/0.
    ov_ba_pre = stand_ba(s)
    ov_avh_pre = stand_top_height(s)
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
        # estab.f:646-650 EMSQR = ±DRAW2 (sign from DRAW1<0.5). These two draws align with the live oracle only
        # for the FIRST plot; from plot 2 on, live's per-plot ESRANN count is inflated by the AUTOES natural-regen
        # tally (STOADJ block + species tally, estab.f:651+) that jl does not model (that is #143). So the CI essubh
        # `disp = EMSQR·DILATE·BNORM` term is #143-entangled and is left at 0 (the deterministic median) — the essubh
        # MEAN (PN) is bit-exact vs live; disp is a stochastic realization straddling 0, .sum-inert on cit01.
        esrann!(s.rng); esrann!(s.rng)
        esdraw = floor(esrann!(s.rng) * 100000f0 + 0.5f0)
        for a in due
            sp = round(Int, a.params[1]); (1 <= sp <= MAXSP) || continue
            ptree = a.params[2] * (a.params[3] / 100f0) / dupnpt
            ptree <= 0f0 && continue
            # a cycle-number date (<1000) resolves to the calendar year at that cycle (cycle_year_at) before the
            # DELAY offset — else `delay = 2 - 2016 = -2014` ⇒ age≈2019 ⇒ grossly over-sized "seedlings". A
            # calendar-year date carries its own sub-cycle offset unchanged.
            pyr    = (0 < Int(a.year) < 1000) ? Int(cycle_year_at(s.control, Int(a.year))) : Int(a.year)
            delay  = pyr - Int(yr)
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
            elseif s.variant isa LakeStates
                # LS ESSUBH (ls/essubh.f:68-77): identical NE/CS-style base height — carage from ls/essubh.f's
                # own MAPLS map, H via the LS NC-128 curve (htcalc IVAR=1), then (H/CARAGE)·min(5, TIME−DELAY).
                carage = Float32(_LS_ESSUBH_REFAGE[sp])
                (ls_htcalc_height(sp, si, carage) / carage) * min(5f0, Float32(per) - Float32(delay))
            elseif s.variant isa CentralRockies
                _CR_ESSUBH_HHT[sp]        # cr/essubh.f: a FIXED per-species base height (not a height-at-age curve)
            elseif s.variant isa InlandEmpire
                # IE NATURAL/PLANT base height — FIRST-CUT placeholder (=XMIN); the DF NATURAL height source
                # is the estb tally path (not essubh, cont.56), to be pinned via esnutr trace + refined vs live.
                _IE_ES_XMIN[sp]
            elseif s.variant isa Teton
                _TT_ESSUBH_HHT[sp]        # tt/essubh.f fixed per-species base height (PP=placeholder); clamped [XMIN,HHTMAX]
            elseif s.variant isa CentralIdaho
                # CI subsequent/planted base height (ci/essubh.f, HHT=EXP(PN + disp·SIG)). IHTSER from the
                # shared estb habitat-bracket chain (em_ihtser == the shared estab MYGRUP→MYHTS map, estab.f:493);
                # IPREP=1/IPHY=3 defaults (esplt2.f:191-192). BAA=overstory competition clamp[1,400]; XCOS/XSIN=
                # cos/sin(aspect)·slope (estab.f:480). disp = EMSQR·DILATE·BNORM (per-stand 2-draw EMSQR × per-species
                # sqrt-shrink DILATE × deterministic BNORML[IAGE]); MEASUREMENT PASS uses disp=0 (deterministic mean),
                # validated .sum-inert on cit01 (planted seedlings stay sub-threshold, never enter the summary TPA).
                let _slo = s.plot.slope
                    ci_essubh(sp, age, clamp(s.plot.basal_area, 1f0, 400f0),
                              em_ihtser(Int(s.plot.habitat_code)), 1, 3,
                              _slo*cos(s.plot.aspect), _slo*sin(s.plot.aspect), _slo, s.plot.elevation, 0f0)
                end
            elseif s.variant isa BlueMountains
                # BM base height (bm/essubh.f): HHT = SMHTGF(sp, MODE=0, DTIME=AGE) = the small-tree height-at-total-
                # age curve. Deterministic — NO EMSQR/DILATE/ELEV (bm/essubh.f discards them). SI = per-species SITEAR.
                bm_essubh_hht(sp, si, age)
            elseif s.variant isa Utah
                # UT base height (ut/essubh.f): a FIXED per-species table (no age/site/EMSQR), clamped [XMIN,HHTMAX]
                # by the shared engine. NOTE: full utt01 validation is gated on the SEPARATE UT sprout crash
                # (esuckr!→essprt_sn dispatch gap: UT uses the ut/esuckr.f Crouch aspen model, not SN essprt) — that
                # is a distinct sprout-subsystem bug, not this establishment-height gap.
                _UT_ESSUBH_HHT[sp]
            elseif s.variant isa EasternMontana
                # EM subsequent/planted base height (em/essubh.f, deterministic EXP(PN)). IHTSER from the habitat
                # code bracket search; IPHY=3 / IPREP=1 defaults (esplt2.f). BAA=overstory competition BA clamp[1,400].
                _slo = s.plot.slope
                em_essubh_hht(sp, log(age), clamp(s.plot.basal_area, 1f0, 400f0),
                              _slo*cos(s.plot.aspect), _slo*sin(s.plot.aspect), _slo, s.plot.elevation,
                              em_ihtser(Int(s.plot.habitat_code)), 3, 1)
            elseif s.variant isa Klamath
                _NC_ESSUBH_HHT[sp]        # nc/essubh.f fixed per-species base height; clamped [XMIN,HHTMAX]
            elseif s.variant isa EastCascades
                # EC base height (ec/essubh.f → ec/smhtgf.f MODE=0): the small-tree height-at-total-age
                # curve HHT = SMHTGF(sp, AGE) with SI = the species' SITEAR. Deterministic (no EMSQR/DILATE/
                # ELEV). The PLANT-no-height branch below then adds the [0,1.5] RAN draw (ec/estab.f:485).
                ec_essubh_hht(sp, si, age)
            else
                htcalc_height(bc, sp, si, age, montane)
            end
            treeht = a.params[5]
            # HTADJ (esin.f opt 15 → esnutr.f 442): per-species height adjustment added to HHT BEFORE the
            # XMIN/0.05 floor and HHTMAX clamp (estab.f:932/1033/1036). Default 0 (empty dict) ⇒ inert.
            hadj = isempty(s.estab.ht_adj) ? 0f0 : get(s.estab.ht_adj, Int32(sp), 0f0)
            if treeht >= 0.1f0                                      # PLANT specified a height
                hht = treeht; xh = log(hht)
                while true
                    xxh = exp(bachlo(s.rng, xh, 0.5f0; stream = :estab))
                    (0.5f0 * hht <= xxh <= 2f0 * hht) && (hht = xxh; break)
                end
                hht += hadj                                        # estab.f:1033 HHT=HHT+HTADJ (before the 0.05 floor)
                hht < 0.05f0 && (hht = 0.05f0)                      # PLANT floor 0.05 (estab.f:1034)
            elseif s.variant isa EasternMontana || s.variant isa CentralIdaho ||
                   s.variant isa BlueMountains || s.variant isa Utah || s.variant isa Klamath
                # Shared estb/estab.f:1035-1037 PLANT (no user height): HHT = essubh + HTADJ(default 0), floor XMIN —
                # NO RAN draw. Only the user-specified-height branch (treeht≥0.1, estab.f:1026-1034) draws the lognormal
                # BACHLO perturbation. jl already consumes the per-replicate EMSQR/ESDRAW draws (line ~218) for stream
                # sync, so skipping this extra RAN keeps the :estab stream aligned vs the live oracle. EM validated;
                # CI added #154 (was wrongly taking the else RAN-branch below → +~0.5 ft spurious height + a stream
                # desync). NOTE (CR/IE/TT): same shared-source no-draw applies, latent behind their essubh branches
                # (their validation used no-PLANT DB stands); fold them in when a PLANT .key is validated per variant.
                hht += hadj                                        # estab.f:1036 HHT=HHT+HTADJ (before the XMIN floor)
                hht < es_xmin[sp] && (hht = es_xmin[sp])
            else                                                   # default: RAN~N(0.5,0.25), accept RAN∈[ran_lo,ran_hi]
                while true
                    ran = bachlo(s.rng, 0.5f0, 0.25f0; stream = :estab)
                    (ran_lo <= ran <= ran_hi) && (hht += ran; break)  # estab.f:483/490 (variant-specific window)
                end
                hht += hadj                                        # estab.f:932 HEIGHT(N)=HHT+HTADJ (before the XMIN floor)
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
            elseif s.variant isa EastCascades
                # ec/estab.f:626 assigns the establishment DBH = 0.1 flat; ec/esgent.f only recomputes DBH
                # when WK4<1 (a partial birth cycle). A full-birth-cycle PLANT/NATURAL tree (WK4=1) keeps
                # DBH=0.1 even after its height exceeds breast height — height grows, DBH does not. (EC also
                # lacks the shared :htdbh_* coef arrays — it has its own ec_htdbh_dbh — so the shared inverse
                # both mis-modeled EC and KeyError-crashed.)
                dbh = 0.1f0
            else
                dbh = _htdbh_dbh(sd, sp, hht, ifor; isne = s.variant isa Northeast); dbh < 0.1f0 && (dbh = 0.1f0)
                dbh += 0.001f0 * hht
            end
            for _ in 1:ibrkup
                n = t.n + 1; n + Int(t.ndead) > length(t.dbh) && break   # leave room for the dead block (t.n+1…t.n+ndead); else the volume loop `1:(t.n+ndead)` overruns the MAXTRE arrays (intermittent SIGSEGV on dense ESTAB stands with inventory dead records)
                t.n = n
                t.species[n]     = Int32(sp)
                t.dbh[n]         = dbh
                t.height[n]      = hht
                t.tpa[n]         = ptree / brk
                t.htimlt[n]      = 1.0f0     # PLANT/NATURAL: full birth-cycle HTG (TRAGE≥GENTIM ⇒ WK4≈1; guards slot reuse). #193
                # ABIRTH = AGEPL + GENTIM (estab.f:628/707) — the REGENT-start `age` already computed above IS
                # FVS's tree age (essubh.f:93). jl left birth_age=0 ⇒ established trees ran ~AGEPL+GENTIM (=7 for a
                # default PLANT) years too YOUNG ⇒ htgf's even-aged site curve (steeper when young) over-predicted
                # height growth as planted stands approached the site asymptote (late-cycle TopHt jl-high). CR-gated:
                # the eastern variants share this latent gap but are separately validated (avoid unvalidated churn).
                (s.variant isa CentralRockies || s.variant isa Teton) && (t.birth_age[n] = age)   # ABIRTH=AGEPL+GENTIM (TT mirrors CR; western even-aged htgf curve)
                # Records go on inventory point `nn` (estab.f:313 ITRE=IPTIDS[nn]).
                # point_ba scales each point's raw BA by PI/GROSPC with PI=NPTIDS, so with
                # the planted TPA spread evenly over NPTIDS points each point_ba comes back
                # to the full stand BA — matching the oracle's pba=ba_v fallback (PTBAA≤0)
                # for fresh establishment, for any NPTIDS (NPTIDS=1 ⇒ this is point 1).
                t.plot_id[n]     = use_iptids ? Int32(iptids[nn]) : Int32(nn)   # IPTIDS[nn] = nn-th STOCKABLE point
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
        ls_estab = s.variant isa LakeStates
        local ebau_e, b3_e, avh_e, scale_e, rdiam_e, rnd_e
        local cb1_e, cb2_e, cb3_e, ba_e
        local lcheck_e, lb1_e, lb2_e, lb3_e, lb4_e, lc1_e, lc2_e, lbamax_e, rmsqd_e
        if ne_estab
            ebau_e = ebau_pre                              # PRE-establishment BAL (snapshot above), not the cohort's
            b3_e = sd[:dg_b3]; avh_e = s.plot.avg_height
            # REGENT LESTB period: FNT = FINT−5 (regent.f:118-124; LSKIPH ⇒ no ht growth when FINT≤5).
            scale_e = per > 5 ? Float32(per - 5) / NE_REGENT_REGYR : 0f0   # CON=HGADJ=XRHGRO=1
            rdiam_e = sd[:regent_min_diam]
            rnd_e = s.control.dg_stddev_bound >= 1f0        # DGSD random ±10%
        elseif cs_estab
            cb1_e = sd[:balmod_b1]; cb2_e = sd[:balmod_b2]; cb3_e = sd[:balmod_b3]
            avh_e = ov_avh_pre                             # POST-growth overstory (pre-seedling snapshot above) —
            ba_e = ov_ba_pre                               # NOT the stale plot.* (regent.f end-of-cycle regen timing)
            scale_e = per > 5 ? Float32(per - 5) / 10f0 : 0f0   # FNT/REGYR, REGYR=10 (CON=HGADJ=XRHGRO=1)
            rdiam_e = sd[:htdbh_db]                         # DIAM floor (= cs/regent.f DIAM == htdbh_db)
            rnd_e = s.control.dg_stddev_bound >= 1f0        # DGSD random ±10%
        elseif ls_estab
            # LS shares NE/CS's REGENT(LESTB) shape: ESSUBH gives a BASE height (5-yr) that this grows to the
            # cycle-end height via the LS NC-128 increment (MAPLS) + ls_balmod (ls/regent.f). BARE stand: BA≈0,
            # RMSQD≈0 ⇒ ls_balmod omega=b4/gm≈1 and AVH=0 ⇒ no competition suppression (live FVSls GMOD=1).
            lcheck_e = sd[:balmod_check]; lb1_e = sd[:balmod_b1]; lb2_e = sd[:balmod_b2]; lb3_e = sd[:balmod_b3]
            lb4_e = sd[:balmod_b4]; lc1_e = sd[:balmod_c1]; lc2_e = sd[:balmod_c2]; lbamax_e = sd[:balmod_bamax1]
            avh_e = s.plot.avg_height; ba_e = s.plot.basal_area; rmsqd_e = rmsqd_pre  # pre-establishment RMSQD (DENSE)
            scale_e = per > 5 ? Float32(per - 5) / 10f0 : 0f0   # FNT/REGYR, REGYR=10 (CON=HGADJ=XRHGRO=1)
            rdiam_e = sd[:htdbh_db]
            rnd_e = s.control.dg_stddev_bound >= 1f0
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
                # XRHGRO = REGHMULT (regent.f HTGR = HTCALC·CON·SCALE·HGADJ·XRHGRO); the LESTB path must apply
                # it too (was hardcoded 1 ⇒ REGHMULT ignored for the establishment cohort, mult_reghmult diverged).
                xrhgro = active_multiplier(s.control, :regh, sp, Int(yr))
                if ne_htcalc_htmax(sp, si) - h <= 1f0
                    htgr = 0.1f0
                else
                    # regent.f:224 HTGR = HTCALC·CON·SCALE·HGADJ·XRHGRO — the LESTB path must apply CON =
                    # exp(htg_cor_small) (= RHCON·exp(HCOR)) too, as NE small_tree_growth.jl:48 does. It was OMITTED
                    # here ⇒ planted seedlings over-grew (WP: CON=0.914, live rawHTGR 7.98 vs jl 8.73; live-stamped).
                    htgr = ne_htcalc_incr(sp, si, ne_htcalc_age(sp, si, h)) *
                           exp(s.calib.htg_cor_small[sp]) * scale_e * xrhgro
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
                    dnew = _htdbh_dbh(sd, sp, hk, ifor; isne = s.variant isa Northeast); dnew < 0.1f0 && (dnew = 0.1f0)
                    dnew < rdiam_e[sp] && (dnew = rdiam_e[sp])
                    t.dbh[i] = dnew + 0.001f0 * hk
                end
                hk > _NE_ES_HHTMAX[sp] && (hk = _NE_ES_HHTMAX[sp])   # HARD HHTMAX clamp on the REPORTED height
                t.height[i] = hk
            elseif cs_estab                                    # CS REGENT(LESTB) height growth + new DBH
                sp = Int(t.species[i]); h = t.height[i]; si = s.plot.sp_site_index[sp]
                xrhgro = active_multiplier(s.control, :regh, sp, Int(yr))   # REGHMULT (was hardcoded 1)
                if cs_htcalc_htmax(sp, si) - h <= 1f0          # HTMAX−H ≤ 1 ⇒ HTG=0.1 (cs/regent.f:206)
                    htgr = 0.1f0
                else
                    htgr = cs_htcalc_incr(sp, si, cs_htcalc_age(sp, si, h)) * scale_e * xrhgro
                end
                # regent.f:156 BAL=(1-PCT/100)·BA, PCT=BA-percentile. A new seedling (D=0.1) is the smallest ⇒
                # PCT≈0 ⇒ BAL≈full overstory BA (live BAL=134.1=BA). crown_ratio was a wrong proxy (gave ~0.5·BA).
                bal = ba_e
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
                    dbhk = _htdbh_dbh(sd, sp, hk, ifor; isne = s.variant isa Northeast)
                    dbhk < rdiam_e[sp] && (dbhk = rdiam_e[sp])
                end
                t.dbh[i] = dbhk + 0.001f0 * hk
                hk > _CS_ES_HHTMAX[sp] && (hk = _CS_ES_HHTMAX[sp])   # HARD HHTMAX clamp on the REPORTED height
                t.height[i] = hk
            elseif ls_estab                                    # LS REGENT(LESTB) height growth + new DBH
                sp = Int(t.species[i]); h = t.height[i]; si = s.plot.sp_site_index[sp]
                xrhgro = active_multiplier(s.control, :regh, sp, Int(yr))
                if ls_htcalc_htmax(sp, si) - h <= 1f0
                    htgr = 0.1f0
                else
                    # regent.f:224 CON = exp(htg_cor_small) (= RHCON·exp(HCOR)), as LS small_tree_growth.jl:40 applies.
                    # Was omitted here (inert for JP where CON≈1, but a latent bug for CON≠1 species — cf. the NE fix).
                    htgr = ls_htcalc_incr(sp, si, ls_htcalc_age(sp, si, h)) *
                           exp(s.calib.htg_cor_small[sp]) * scale_e * xrhgro
                end
                gmod = ls_balmod(sp, t.dbh[i], ba_e, rmsqd_e, lcheck_e, lb1_e, lb2_e, lb3_e, lb4_e, lc1_e, lc2_e, lbamax_e)
                relht = avh_e > 0f0 ? min(h / avh_e, 1f0) : 0f0
                htgr = max(htgr * (1f0 - (1f0 - gmod) * (1f0 - relht)), 0.1f0)
                if rnd_e
                    rh = 0f0
                    while true; rh = bachlo(s.rng, 0f0, 1f0); -1f0 <= rh <= 1f0 && break; end
                    htgr = max(htgr + rh * 0.1f0 * htgr, 0.1f0)
                end
                hk = h + htgr
                if hk < 4.5f0                       # ls/regent.f LESTB dbh: DIAM floor (or htdbh⁻¹), + 0.001·hk
                    dbhk = rdiam_e[sp]
                else
                    dbhk = _htdbh_dbh(sd, sp, hk, ifor; isne = s.variant isa Northeast)
                    dbhk < rdiam_e[sp] && (dbhk = rdiam_e[sp])
                end
                t.dbh[i] = dbhk + 0.001f0 * hk
                hk > _LS_ES_HHTMAX[sp] && (hk = _LS_ES_HHTMAX[sp])   # HARD HHTMAX clamp on the REPORTED height
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
