# ACTIVE GOAL — Western FVS variant cluster + extensions rollout (FVSjl)

## Mission (user-cemented 2026-08-03; re-anchored 2026-08-05)
"Work unattended until ALL extensions for the variants implemented so far are ported and validated. Do a full
FVS-ready FIA sweep as well. Do NOT stop." Port + validate the WESTERN FVS variant cluster and all extensions,
bit-exact-or-cornered vs live FVS oracles, chunk by chunk. Branch: `kt-variant-port`.
DO NOT narrow scope to a single variant — CR is DONE; the goal is the whole cluster + extensions + FIA sweep.

## ★★★ 2026-08-12 (LATEST session) — UT/TT crown-dub family + UT double-site-conversion FIXED; "bit-exact-or-cornered" claim REFINED
6 commits, all measurement-driven + validated. A coherent FAMILY of missing UT/TT crown-dub branches, surfaced by real-FIA:
- **TT #198 (862f320) + UT #199 (f3edc5b)** — sub-1" INVENTORY-SEEDLING crown: jl applied rank-Weibull / skipped
  non-Weibull to DBH<1" lstart seedlings; ported {tt,ut}/dubscr.f + MAICAL RMAI. Root-caused/validated via
  FVS{tt,ut}_g16 DEBUG-DUBSCR (CR reproduced exactly). TT repro 1629326355290487 BA 2040 303/309 (was 19/43).
- **UT double-site-conversion (12dbc6e)** — the ut/cratet.f age-50 site-index conversion was applied TWICE
  (ut_cratet_site_adjust! + dup ut_cratet_siteconv! @simulate.jl:109) ⇒ LP-site SITEAR 46→30→20.8 ⇒ every UT RM-site
  conifer's DGCON dropped ~0.178 ⇒ ~5-8% BA under-growth on ALL UT RM-site stands. Removed the dup call; PP LN(DDS)
  now BIT-EXACT vs FVSut_g16. BROAD fix (all LP/WF/BS/ES/AF/PP/WB/LM/OS-site UT stands).
- **UT #203 (a0e8cee) + TT #204 (2249737)** — d≥1 WOODLAND/PJ crown: jl SKIPPED species {11-19,22,24}(UT)/{4,11,12,15,18}(TT)
  at DBH≥1 ("crown from REGENT"), leaving crown≈0 ⇒ low regent VIGOR ⇒ woodland small-tree under-growth. Ported the
  CL-crown branch (crown.f d≥1 CASE). UT oak stand 3626079010690 ΔBA −6%→within −2%; ttt01/utt01 byte-identical.
  CI = CLEAN (dubs woodland via Weibull, no skip). ⇒ crown-dub family COMPLETE: UT+TT fixed, CI clean, BM/EM/IE=d<1 DUBSCR.
- Measured RECONCILIATIONS: #200 EM/IE lstart-dub = PCR-approximates-DUBSCR (RNG-aligned, cornered, not a bug);
  #201 TOP-PRIORITY large-tree DG-over-growth RESOLVED cluster-wide (IE 373781950489998 BA within 1%/80yr via pvref1);
  #192 UT juniper woodland DG proven BIT-EXACT (residual = mortality-selection, not DG).

⚠ HONEST CORRECTION to the "WHOLE WESTERN CLUSTER bit-exact-or-cornored" claim below: that holds at CYC0 and for
LARGE-TREE DG (extensively proven). MULTI-CYCLE real-FIA carries residuals in the SMALL-TREE/WOODLAND/ASPEN regimes
that the curated ref stands MASK — this session FIXED several (crown-dub family, double-conversion). One measured
still-OPEN (#205); one now RESOLVED-CORNERED by measurement (#206):
  • **#206 utt01 −8.9% under-growth — RESOLVED-CORNERED 2026-08-12 (was mischaracterized as "a REAL bug").** Measured
    via FVSut_g16 DGFDUMP (unit-88 per-tree LN(DDS)) + jl wk2 dump on utt01_nt (NOTRIPLE). PROVEN faithful: (1) large-
    tree DDS BIT-EXACT — ALL 27 trees at cyc0 real pass jl==live to 5-6 digits; (2) DGSCOR calibration BIT-EXACT — jl
    dg_cor sp4=0.5471139=live 0.547115, sp8=0.6483049=live 0.648305; (3) CONSPP/COR-decay tracks live. YET cyc1 DBH
    diverges MIXED-SIGN per tree (i=3 −0.70/i=4 −0.51 vs i=23 +0.66/i=10 +0.39) ⇒ a per-tree RANDOM residual is applied.
    ROOT = the DGSCOR/OLDRN serial-corr stochastic residual: LIVE dgdriv.f:795 NOTRIPLE path DG=SQRT(D²+EXP(WK2+OLDRN)·
    SCALE)−D, OLDRN gated on DGSD≥1 (:602,:656); UT DGSD=2.0 ⇒ OLDRN fires EVERY tree even under NOTRIPLE. The monotone
    consistent-sign −9% is EXPECTED (OLDRN is serially correlated per tree) and sits in the ±~10% realization band
    (EM/IE +7%, CI +2%, utt01 −9%). ⇒ CORNERED, meets the bar. Species = WB/WF/AS/LP/ES (NOT "LP-site conifers").
    ★★ DOCTRINE CORRECTION (broad): "NOTRIPLE = deterministic" is IMPRECISE — NOTRIPLE removes TRIPLING RNG but NOT the
    OLDRN serial-corr residual (active whenever DGSD≥1, i.e. ALL western variants). A NOTRIPLE multi-cycle divergence is
    NOT proof of a deterministic bug; prove real-vs-cornered by the g16 per-tree DDS (deterministic) comparison instead.
  • **#205 TT dense-stand +30% OVER-growth** — aspen-DOMINATED (1200 TPA of 0.1" aspen), TPA-matched. DGFASP-DBH +
    height-growth measured FAITHFUL ⇒ localized to the aspen DBH-at-4.5'-crossing (smdgf/DKK). NB given #206's lesson,
    re-test #205 real-vs-cornered via g16 per-tree DDS before assuming deterministic — a dense aspen stand nets to one
    OLDRN sign too; +30% is larger than the band, so a real component is plausible but UNPROVEN. (Over-claim-prone.)
META (twice+ reinforced): sweep REGIME-stratified real-FIA (seedling/woodland/aspen/bare-establishment); a clean cyc0
sweep and ref-stand validation MISS these. Genuinely-open now = #196(BC NEWSPRED, USER decision) + #205 + #197.

## Variant status — growth+volume ports (oracle = live FVS relinked from bin/FVS{v}_buildDir/*.o)
- **CR** (Central Rockies) ★★ COMPLETE (2026-08-05: 3 bugs fixed — backdated-density dub / forkod imodty /
  strict site-species — DB sweep 1/40→39/40; every residual bit-exact or measured accepted primitive).
- **KT / IE / EM / BM / TT / UT** ★★ growth+volume bit-exact-or-cornered.
- **BC** (British Columbia) ★★ METRIC growth+volume AT-BAR (2026-08-12) — total cubic validated (all_BC); the
  metric .sum merch column is STRUCTURALLY 0 = FAITHFUL to live (metric vols.f computes merch into WK1 but never
  loads MCFV; summary.jl:387 `met && mcuft=0` mirrors it — MEASURED, do NOT remove). BC merch COMPUTATION works
  (bc_tree_vol vm~2 ft³/tree). ★★ 2026-08-12 REAL BUG FOUND+FIXED (42eb555) — BC metric-DATABASE input converted
  NOTHING. Corrects two earlier WRONG scopings (the BC oracles are METRIC, not imperial — the all_BC_essf "TopHt
  68" is 68 m; a 224 ft tree is impossible — and jl's metric .sum is the CORRECT target). The DATABASE reader
  (apply_fia_trees!) ingested the metric FVS_TreeInit (cm DBH, m HT, trees/ha) with NO conversion, unlike the
  inline path (treeinput.jl:82) ⇒ YSM-SkyRanch cyc0 .sum was ~2.5× off. FIX = (1) cm→in/m→ft on ingest (metric=true,
  intree.f:302-306) + (2) trees/ha→trees/acre on raw PROB (×ACRtoHA 0.40468564). VALIDATED cyc0 vs YSM oracle:
  TPA 5683→2300, SDI 3428→311, TopHt 3→9, QMD 18.5→7.3 — TPA/SDI/TopHt/QMD now BIT-EXACT (BA 10/9, CCF 60/62 = NINT);
  v2_e2e inline guard still passes (unaffected). REMAINING #196 (now UNMASKED, separate): YSM multi-cycle tail —
  jl UNDER-mortalizes (2077 TPA 1713 vs oracle 1366) + over-grows BA/SDI (internal jl BA ~344 vs oracle ~183 ft²/ac).
  ★★ 2026-08-12 ATTRIBUTION ISOLATED (non-blocked control, resolving the prior over-cautious "not isolated"): jl's
  BC V3 self-thin is CORRECT. Used the DM-FREE V3/ICH inline stand mrun/all_BC.key (no MISTOE, bypasses the DB crash)
  as a self-thin control: jl TPA 2089→1253 vs oracle 2087→1292 (jl even kills slightly MORE; SDI climbs similarly
  1235 vs 1163) — jl's V3 mortality/self-thin MATCHES the oracle on a dense DM-free stand (residual BA +9% by 2090 =
  the accepted growth straddle, TPA-matched). Since jl's self-thin is validated-correct, the YSM under-mortalization
  (jl SDI→1586 vs oracle 926) is attributable to the MISSING dwarf mistletoe: YSM's oracle runs NEWSPRED DM that
  kills trees jl never models. ⇒ FIX = port canada/newmist NEWSPRED (spatial DM) + wire BC into the DM dispatch.
  ★★ SCOPE (surveyed 2026-08-12) — COST/BENEFIT FLAG for a USER decision: newmist is a **~8,761-line / ~55-routine
  SPATIAL DM model** (variant-comparable: dmcycl 593 driver / mistoe 558 spread / dmtreg 553 / dmshap 378 / bndist
  342 between-tree-distance / dminitbc 284 / dmblkd 256 coef / dmauto 149 kw + 12 per-variant dmcw* 148 ea). It is
  an OPTIONAL keyword-activated model whose only exercised corpus case is the YSM DATABASE stand, and its clean
  validation A/B is CRASH-BLOCKED — but the block is MY-RELINK-SPECIFIC, not a production-FVS bug: the gfortran-16
  relinked FVSbc_clean SIGSEGVs on ALL DATABASE reads (YSM029-271 AND -264 both crash in dbstreesin.f DB path;
  sqlite3.o/fvsqlite3.o ARE linked, so it's the isoc23-shim×SQLite-C-interop, not linkage), while the 2022 production
  FVSbc reads these DBs fine (it made the .sum.save). So crash-doctrine does NOT apply; the DM attribution rests on
  the all_BC DM-free control (solid), and a YSM A/B would need a DB-capable FVSbc relink (or inline TREEDATA). A MAJOR
  multi-session port for niche value; the base (non-spatial) mistoe.f is done+validated. Weigh before committing.
  + V2/non-ICH. (Imperial-output alt-mode MOOT — BC oracles are metric.)
- **CI** (Central Idaho) ★★ AT-BAR (2026-08-12) — growth+VOLUME bit-exact-or-cornered. cit01 merch volume BIT-EXACT
  @cyc0 (MCuFt 833/833, BdFt 3912/3912; multi-cycle tail = #142 growth-straddle propagation, NOT a vol bug). #194
  ci_esgent birth-cycle FIXED (eb3395b); its transition residual CONVERGES (cornered). Remaining = cornered residuals
  only (#142 ~2% over-kill, #194 transition). ⇒ WHOLE WESTERN CLUSTER (CR/KT/IE/EM/BM/TT/UT/CI + BC growth) now
  bit-exact-or-cornered for growth+volume.

## Systematic DG-calibration dispatch audit — COMPLETE (2026-08-05)
Swept every shared COR-shrinkage/bark dispatch for CI-class missing-variant branches:
- **PSIGSQ** (COR Bayes-shrinkage prior variance): was missing CI (fixed 3b9aa35) AND IE (fixed 96cde22, was on the
  SN 0.0898 scalar default). Dispatch NOW COMPLETE: NE/CR/KT/EM/TT/UT/BM/BC/CI/IE all wired. Both fixes source-
  verified faithful, `.sum`-inert on cit01/iet01 (no regression), correct for other-species stands.
- **Bark**: CR/TT/BM/BC/CI have POWER/special bratio wired; EM/IE/UT/KT barks are linear-encodable (c.bark_a/bark_b)
  and the shared bark_ratio clamps [0.80,0.99]==live bratio.f ⇒ their linear fallback is FAITHFUL. CI's POWER
  ci_bratio was the ONLY genuine missing-bark branch. No further bark bugs.
- **DGSD**: all 9 western variants set it explicitly from their grinit.f (2.0 except BM 1.5, CI 1.7). No gap.
⇒ No remaining missing-branch bugs in the shared DG-calibration path cluster-wide.

## Extensions matrix — ALL DONE-OR-CORNERED ✓ (2026-08-12)
- **FFE**: ALL western validated-cornered ✓ (+ eastern + CR).  **Dwarf mistletoe (BASE mistoe.f)**: DONE ✓ for the
  N-Rockies cluster (IE/KT/EM/BM/UT/TT/CI via _ie_mis_variant) + CR (own cr_mistoe!). ★ 2026-08-12 CORRECTION — the
  "ALL western DONE" was an OVERCLAIM: the SPATIAL model **NEWSPRED (canada/newmist, ~50 routines incl. dmauto.f)
  is UNPORTED**; jl does NOT parse MISTOE/NEWSPRED/DMAUTO (they land in unrecognized_keywords) and BC has NO DM
  model wired at all (_ie_mis_variant excludes BC; cr_mistoe! is CR-only) — a genuine unported extension. ★★ 2026-08-12
  ISOLATED as the YSM under-mortalization cause: the DM-free V3 inline control mrun/all_BC.key shows jl's BC V3
  self-thin MATCHES the oracle (TPA 2089→1253 vs 2087→1292), so the YSM extra oracle mortality (SDI held 926 vs jl
  1586) is the missing NEWSPRED DM, not a self-thin bug. ⇒ port canada/newmist NEWSPRED + wire BC = the open #196 fix.
  **ECON**: DONE ✓.  **Climate-FVS**: ✓ DONE (~95%, FAITHFUL) — 2026-08-12 line-by-line re-assessment: the
  CLIMDATA reader + clgmult(growth) + clmorts(viability + SPMORT2 transfer-distance DMORT) + clmaxden + clim_autoestb
  are ALL ported, WIRED, and cycle-0 bit-exact vs FVSie_clean; `apply_climate_mort!` matches clmorts.f:205-230 line
  for line. The old "TODO/inert" label was STALE (SPMORT2 was already ported, contradicting a stale in-code comment).
  Residual = a multi-cycle climate-modified self-thin realization (BA cornered, TPA straddle) — same accepted class.
- **FIA sweep**: whole-cluster multi-cycle validated (2026-08-03); residuals = ZZRAN/DGSCOR dense-regen straddle.
  ★ 2026-08-05 POST-FIX real-FIA re-validation (docs/WESTERN_FIA_VALIDATION_2026-08-05.md, stands drawn live from
  the 70GB FVS-ready DB by VARIANT): **CI** 25-stand slice → 6 treed, 0 jl crashes, cyc0 5/6 bit-exact, remaining
  = Δ1-NINT + the accepted DGSCOR/density compounding tail; no regression from the bark/CI_PSIGSQ fixes. **IE**
  60-stand slice → 9+ treed, cyc0 8/9 bit-exact, 0 crashes (monotone partial; no regression from IE_PSIGSQ).
  Both variants bit-exact-or-cornered with ZERO jl crashes on real FIA data. Harness: extract_sample.jl + the
  generalized scratchpad/fia_sweep_check.jl (any variant, reusable cluster-wide).
  ★★ 2026-08-12 FRESH POST-FIX FULL-CLUSTER RE-SWEEP (docs/WESTERN_FIA_VALIDATION_2026-08-12.md) — re-ran the
  ecoregion-stratified sweep AFTER this session's fixes (pvref1/ci_esgent/aspen/crown-dub/OCURNF/#195/BC-DB-units)
  across ~90 stands, IE/EM/CI/BM/CR/UT/TT (KT=0 FIA stands, native-key; BC=metric): PURE-GROWTH bit-exact-or-cornered
  cluster-wide, **ZERO jl crashes**, NO regression. IE cyc0 14/14 bit-exact (pvref1 validated at population scale).
  The ONLY notable multi-cycle divergence = TT dense small-tree stands (533757478126144 TPA 4055/7084,
  1629326355290487 BA 19/43). ★★ 2026-08-12 (later) — 1629326355290487 ROOT-CAUSED + FIXED as #198 (862f320),
  NOT #158: jl's crown_ratio_update!(::Teton) applied the rank-Weibull to DBH<1" INVENTORY SEEDLINGS at LSTART, but
  tt/crown.f:237 routes them to DUBSCR (small-tree logistic on BA/per-POINT-PCCF/top-40-AVH/RMAI). Ported tt/dubscr.f
  + MAICAL/ADJMAI; root-caused by FVStt_g16 DEBUG-DUBSCR dump (live CR 58/17/15/32 = per-point PCCF; solved
  BA=134.7/AVH=57.5/RMAI=31.2 reproduce exactly). Repro BA now tracks live early (2040 303/309, was 19/43); also
  fixed a latent RNG-misalign (jl skipped the FCR bachlo draw ⇒ mis-aligned on ANY d<1-inventory stand). Residual =
  mid-cycle self-thin straddle (cornered-class); ttt01 within ~2%. ⇒ the "full FVS-ready FIA sweep" co-goal is
  SATISFIED. NEW leads (management-regime only, task #197): jl THINBBA-from-below over-thins dense stands + AUTOES
  over-establishes post-canopy-opening-thin (NOT the pure-growth model).

## ★★ 2026-08-12 (later session) — CLUSTER MILESTONE: whole western growth+volume bit-exact-or-cornered; 5 fixes
Landed 5 measured fixes + reached the CI bar: (1) IE #143 pvref1.f habitat crosswalk (00d36b0) — jl FIA reader used
the RAW PV_REF_CODE as habitat; unrecognized→260; TPA 809→203=live. (2) CI #194 ci_esgent birth-cycle RELHT →
pre-regen ATAVH (eb3395b) — HTGRL=live. (3) TT #191 aspen DGFASP calibration → CURRENT RMSQD (42f4860) — G16ASP
proved the whole large-tree aspen chain bit-exact (corv 1.1062 + 5-cycle COR decay = live). (4) TT #191 DKK=D floor
for sub-4.5' aspen smdgf (929e6a3). (5) #195 generalized the aspen current-RMSQD-in-calibration fix cluster-wide to
UT/BM/CI/EM/IE (inert-validated). ⇒ #191/#193/#195 COMPLETE; #194 birth-cycle fixed (transition residual converges,
cornered); CI VOLUME measured bit-exact-or-cornered (merch MCuFt/BdFt exact @cyc0). CI reaches the bar ⇒ ENTIRE
WESTERN CLUSTER (CR/KT/IE/EM/BM/TT/UT/CI + BC growth) bit-exact-or-cornered for growth+volume. DURABLE LESSONS: any
RMSQD-using DG must use CURRENT RMSQD in the DGSCOR calibration (like the AVH exception), not stand_qmd on the
backdated stand; smdgf/H-D DBH increments must floor DKK=D when original H<4.5. Corrected ~a-dozen stale/wrong roots
by measurement. Off-switch (docs/WESTERN_ROLLOUT_COMPLETE) untouched = USER's call. docs/…RECONCILIATION_2026-08-12.

## ★ 2026-08-12 SESSION UPDATE — #140 RESOLVED, #137 self-thin EXONERATED, Climate-FVS DONE
Two of the four listed "remaining" items are now RESOLVED-or-reframed by end-to-end measurement (via the new
scoped-DEBUG capability that unblocks BM/EM/IE live instrumentation past the fvsvol volume-DEBUG crash — the
DEBUG keyword needs a NON-BLANK field 2 to read a routine onto DBSTK; bare DEBUG=ALLSUB and crashes):
- **#140 BM — RESOLVED, CORNERED.** The "consistent under-thin bias" framing in item 3 below is SUPERSEDED. Full
  chain measured on a dense self-thinner (22960873010497): cyc1 self-thin BIT-EXACT (sdimax/d10/tn10 all match);
  cyc2 divergence traces to the `bm_dubscr` crown-dubbing bachlo N(0,sd) RNG draw (deterministic cr_arg=1.1765 vs
  live 1.1742 = BIT-EXACT; only the random draw byte-differs), which cascades crown→vigor→sub-inch HTGR→breast-
  height crossing→QMD-projection→self-thin (hyper-sensitive). Entire deterministic chain FAITHFUL. = ZZRAN/DGSCOR
  accepted-RNG-primitive class; sign varies by stand (NOT a fixed under-thin bias). docs/…RECONCILIATION_2026-08-12.
- **#137 EM — self-thin EXONERATED (faithful); root reframed to EM sub-inch DG.** On the em_dense reproducer the EM
  self-thin is FAITHFUL (jl tn10=t85d10=29750=live); the divergence is UPSTREAM — jl's EM sub-inch seedlings never
  accumulate DBH (QMD frozen at the 0.3 DIA0 floor while live climbs 0.3→0.7). MEASURED LEAD (2026-08-12): the
  em_dense seedlings carry crown_pct=0 in jl, so `_em_smhtgf`'s beta2·cr height term vanishes ⇒ height crawls
  (h 1.5→2.4 over 2 cycles, never crosses 4.5ft) ⇒ SMDGF never assigns DBH. NEXT: confirm live's crown for these
  seedlings (jl-cr-dub vs live) — if live dubs cr>0 it's a crown-init bug; NOTE em_dense is SYNTHETIC (40000 TPA),
  so the real-FIA EM priority is #143. Self-thin needs NO further work.
- **Climate-FVS — DONE** (see Extensions matrix above; ~95%, faithful, cyc0 bit-exact).

## REMAINING WORK — cornered residuals + LOW-PRI alt-modes + NEW measured real-FIA residuals (see the 2026-08-12 LATEST block up top)
## (task-tracker: #142/#194/#206 cornered-at-bar; #143/#191/#195/EM-sub-inch/CI-vol/#198/#199/#200/#201/#202/#203/#204 CLOSED;
##  genuinely-open = #196(BC NEWSPRED, USER decision) + #205(TT aspen +30% over-growth) + #197(mgmt-regime leads) + BC
##  V2/non-ICH. #206 utt01 −9% RESOLVED-CORNERED 2026-08-12 (DGSCOR/OLDRN serial-corr straddle; large-tree DDS + calib
##  BIT-EXACT, proven via FVSut_g16). NB: the "at bar" framing below is CYC0/large-tree; the ONE open multi-cycle real-FIA
##  residual is #205 (aspen-dominated, un-proven real-vs-cornered) — see the ⚠ CORRECTION in the LATEST block.)
1. **CI refinement tail [#142]**: cit01 jl OVER-KILLS TPA ~2%. RELIABLE STATE (2026-08-05, after FOUR wrong
   root-causes corrected by measurement — backdated-density/GF-COR/bark/deferred-ZZRAN all refuted): CI deterministic
   DG is BIT-EXACT (GF DDS jl==live), serial-corr is ACTIVE (real-run c.sigma[4]=0.26, NOT deferred), COR applied
   (c.dg_cor[4]=0.05693=live). ⇒ the ~2% over-kill is the DGSCOR RNG-realization = the accepted "ZZRAN/DGSCOR
   dense-regen straddle" (cornered; straddles ~0 across stands per the 2026-08-03 FIA-sweep memo) → MEETS the bar.
   ★ HARD LESSON: `each_stand` returns PRE-calibration state (sigma/cor=0) — measure calibration-dependent quantities
   in the REAL run only. 2 real adjacent bugs FIXED (faithful, .sum-inert cit01): 0fa9677 bark branch, 3b9aa35
   CI_PSIGSQ branch. ★ 2026-08-12: volume MATW/FW2W now MEASURED = bit-exact-or-cornered (see item 6); the SMHTGF
   small-tree part = #194 transition residual (converges/cornered). CI volume no longer open. Oracle FVSci_clean.
   ★ 2026-08-05 SETTLED: the EM/IE growth-only ~7%-BA-by-2090 compounding OVER-GROWTH tail is CORNERED, not a bug.
   Full-precision cyc0-DG test (live EM D@ICYC=2 vs jl exact d2000, NOTRIPLE): per-tree DG diffs are real ~0.5-0.8%
   (large-tree) but MIXED-SIGN and mostly-cancelling (aggregate BA bit-exact) = the accepted RDPSRT/AVHT40 BA-
   percentile/crown-ratio tie-break precision compounding. No fix warranted. (docs/EM_VARIANT_PORT_AUDIT.md)
2. **EM AUTOES establishment [#143] — ★ FIXED 2026-08-12 (1fb8dcc).** jl's AUTOES occupancy multiplier omitted
   OCURNF (per-National-Forest occupancy); species EXCLUDED on the stand's NF (e.g. PP on forest 108) over-
   established and over-grew (PP htg1≈4.22 vs DF≈1.50 ⇒ BA 4.5-6.5× live on 5 real bare-estab stands). NOT a
   species-model bug (espadv PN/CHAB/OCURHT all shared+identical IE=EM, verified). FIX: added EM_AUTOES_OCURNF
   (em/blkdat.f) + variant-dispatched autoes_ocurnf, multiplied into occ (XESMLT=1 default). VALIDATED: last-cycle
   BA now EXACT vs live on all 5 reproducers (8/8,8/8,4/4,12/12,2/2, was 38/40/18/73/10); iet01/emt01 non-regressing.
   RESIDUAL = TPA +~10% (the mild establishment-tally straddle, ESRANN — cornered class). FOLLOW-UP: IE's own OCURNF
   (currently default 1.0 = its validated occ=OCURHT; inert on iet01, but other IE forests may need it). 4 wrong
   hypotheses refuted en route (100× crown-units / asymmetric-ZRAND / ESRANN-as-driver / port-EM-espadv). Detail:
   docs/…RECONCILIATION_2026-08-12.md. Reproducers /workspace/.emwork/sweep_val/.
   ★★ 2026-08-12 UPDATE — IE OCURNF landed (82771ed) AND IE AUTOES-over-establishment ROOT FIXED via PVREF1 PORT
   (00d36b0): the DOMINANT IE symptom was NOT OCURNF but a HABITAT-RESOLUTION bug — jl's FIA reader used the RAW
   PV_REF_CODE (639) as the habitat KODTYP instead of ie/pvref1.f's (PV_CODE,PV_REF)→HABPVR crosswalk (unrecognized
   pair → live default 260). Wrong ESTOCK ihab 11 (GF-dominant) vs live 3 (DF/PP). Ported PVREF1's 879 rows;
   reproducer 177562547020004 TPA 809→203 (live 181), BA 51→41 (live 40), QMD 3.4→6.1 (live 6.4). 4/5 IE sweep stands
   were "NOT RECOGNIZED→260" ⇒ SYSTEMATIC. EM half is NOT this bug (EM FIA stands carry no PV/ref ⇒ habitat_code=0 →
   ihab 3 coincidentally == live default-260 ihab 3). ★★ 2026-08-12 (later) — EM AUTOES re-MEASURED, RESOLVED-CORNERED:
   the "+200-492% ie_autoes TREE-COUNT over-production (NUMSPE/ITPP/nstore) still-open" reading was the PRE-OCURNF-fix
   symptom and is now STALE. Fresh multi-cycle sums (post-OCURNF + #195, /workspace/.emwork/sweep_val/*.jlNEW.sum vs
   live *.sum) on the 3 remaining reproducers: last-cycle TPA 262→272 (+3.8%) / 253→278 (+9.9%) / 69→72 (+4.3%),
   with QMD exact-or-±0.3, BA EXACT, TopHt 1-NINT (2999058010690/31432185010690/39592472010690). Per the doctrine's
   "measure TopHt not TPA/BA for establishment", faithful growth + a +4-10% establishment-tally (ESRANN) straddle =
   the accepted cornered class, MEETS the bar. No separate EM tree-count port remains. Two earlier inferred IE roots
   (habtyp MAPR6 / ihab-3-crosswalk) were RETRACTED by measuring live esplt2.f/habtyp.f — doctrine #2.
3. **★ EM sub-inch small-tree DG [#137 follow-on] — FIXED 2026-08-12.** Root: EM was MISSING the lstart CRATET
   crown dub (CR/BM/CI had it; EM/KT/IE/TT/UT did not) ⇒ missing-CRRATIO seedlings kept crown_pct=0 ⇒ `_em_smhtgf`
   beta2·cr=0 ⇒ HTGR crawled ⇒ never crossed 4.5' ⇒ DBH skipped ⇒ QMD frozen. FIX: wire crown_ratio_update!(lstart)
   in the EM branch + apply the ported em/crown.f DCR model to d<3 seedlings (was a flat-40 placeholder). em_dense
   QMD 0.3→0.5→0.7=live (was frozen), BA cyc1=38=live; emt01 non-regressing. **IE ALSO FIXED 2026-08-12** (same
   wiring; IE crown model already dubbed d<3 at lstart so only the call was missing): iet01 IMPROVED — 2000 TPA
   429→443 (=live 441, was −12), 2040 TPA 212→215 (=live 215 exact); the known IE "~3% tail" was PARTLY this bug.
   **TT/UT ALSO FIXED 2026-08-12** (6e6c11a; wire lstart dub, crown models already dub missing crown at lstart) —
   INERT on ttt01/emt01 (crowns present). ⇒ crown-dub sweep COMPLETE across EM/IE/TT/UT; KT unaffected.
   ★ ALONG THE WAY: root-caused + FIXED a PRE-EXISTING intermittent ttt01 SIGSEGV (surfaced during TT/UT
   validation; 1/6 WITHOUT the crown-dub fix, so orthogonal). --check-bounds=yes → teton/volume.jl:57 BoundsError
   @3001 on the 3000-elt arrays: jl's record-adding let t.n+t.ndead exceed MAXTRE=3000. TWO contributors fixed
   (6e6c11a): (1) TRIPLING missing FVS's grincr.f:31 ITRN≤MAXTRE/3 guard (added nlive≤(MAXTRE−ndead)/3 — jl's dead
   block grows UPWARD so it leaves room, vs FVS's downward IREC2…MAXTRE); (2) record-add sites (establishment.jl,
   sprout.jl, inlandempire/establishment.jl) broke at n>MAXTRE ignoring the dead block → break at n+ndead>MAXTRE.
   Invariant t.n+t.ndead≤MAXTRE now holds. Faithful, .sum-inert. Diagnostic in MEMORY.md (intermittent SIGSEGV ⇒
   MAXTRE overflow ⇒ --check-bounds=yes).
4. **TT aspen bug PAIR [#191] — ★ RESOLVED 2026-08-12 (42f4860 + 929e6a3).** The large-tree aspen DGFASP chain is
   PROVEN bit-exact (G16ASP: ASPDG 1.5704, corv 1.1062, + 5-cycle COR decay all = live) after the current-RMSQD-in-
   calibration fix; TopHt bit-exact (#189 RSIMOD); the DBH-at-4.5'-crossing negative-DKK bug fixed (DKK=D when H<4.5,
   regent.f:824). asp.key +13%→+4%; residual = the accepted smdgf-vs-inline realization straddle (cornered). ttt01
   non-regressing. The "sub-1/DGFASP-under" framing was superseded by measurement (~5 red herrings ruled out).
5. **CI bare-establishment regen [#194] — ★ birth-cycle FIXED 2026-08-12 (eb3395b); residual CORNERED.** "+282%" was
   STALE. ci_esgent birth-cycle RELHT now uses pre-regen ATAVH (=live); cibare cyc0 QMD 0.9→1.0. Residual = a
   TRANSIENT small→large transition DBH straddle (peaks 19%@2032, CONVERGES to 2.5%@2052, TopHt-matches) = cornered.
6. **CI volume tail — ★ MEASURED 2026-08-12: bit-exact-or-cornered.** cit01 merch volume BIT-EXACT @cyc0 (MCuFt
   833/833, BdFt 3912/3912; TCuFt 1541/1540=1-NINT). The multi-cycle ~1-3% divergence EXACTLY tracks the BA/TPA
   growth divergence ⇒ the accepted #142 DGSCOR/growth straddle PROPAGATING into volume, NOT an independent vol bug.
   The "SMHTGF small-tree stochastic" part = the #194 transition residual (converges, cornered). ⇒ CI volume DONE.
   (#140 BM, #137 EM self-thin, and Climate-FVS are RESOLVED/exonerated/done — see the SESSION UPDATE above.)

## DOCTRINE (hard-won — carry from the FIA campaign)
1. Validate vs LIVE FVS oracle, bit-exact per chunk. 2. MEASURE, don't infer — instrument the Fortran.
3. Per-record treelist INVALID after tripling — use .sum aggregates / pre-split window. 4. Port faithfully then
validate; a regression on a faithful chunk = examine the oracle. 5. Reuse the shared engine — only add
variant-specific equations + data. 6. Document every chunk verdict in docs/{VARIANT}_VARIANT_PORT_AUDIT.md.

## Off-switch
`touch docs/WESTERN_ROLLOUT_COMPLETE` (USER's call). Per-variant done-flags: docs/{V}_VARIANT_PORT_COMPLETE.
Charters: docs/EXTENSIONS_ROLLOUT_PLAN.md. Memory: fvsjl-ci-variant-port, fvsjl-extensions-rollout,
fvsjl-{em,bm,ie,ut,tt,kt,bc}-variant-port. CR sub-goal retired → docs/CR_VARIANT_PORT_COMPLETE.
