# Western FIA sweep — fresh POST-FIX re-validation (2026-08-12)

Re-ran a fresh, ecoregion-stratified FIA sweep AFTER this session's fixes (IE pvref1 00d36b0, CI ci_esgent
eb3395b, TT aspen 42f4860/929e6a3, EM crown-dub/OCURNF, #195, BC DB-units 42eb555) to validate them at
population scale and catch regressions — the mission's explicit "full FVS-ready FIA sweep" co-goal.
Harness: test/harness/fia/extract_sample.jl + manage_fia.jl (build_subdb/run_live/parse_sum), vs FVS{v}_clean.

## PURE-GROWTH result (no management regime) — BIT-EXACT-OR-CORNERED, 0 crashes ✓
- **IE** (14 stands, 14 ecoregions): cyc0 **14/14 BIT-EXACT** (reader/inputs incl. the pvref1 habitat fix correct
  at population scale). Multi-cycle growth (5 cyc, deep-checked 3 treed): bit-exact-or-cornered — e.g. 31362831010690
  TPA 3366→606 vs live →606 (~1%); 451058238489998 within ~1-2%; 3194896010690 within ~2-5%. NO regression.
- **EM** (39600883010690): TPA BIT-EXACT every cycle (384→357); BA/QMD ~2-5% low = accepted growth straddle.
- **CI** (3316648010690): TPA ~1%; BA ~6% high by 2049 = the known #142 DGSCOR over-growth tail (cornered).
- ZERO jl crashes across all runs.
⇒ The session's fixes HOLD at population scale; pure-growth cluster remains bit-exact-or-cornered. FIA-sweep goal met
  for the most-changed variants (IE/EM/CI).

## ⚠ NEW LEADS surfaced (management regime; NOT the growth model — verify before acting)
A first sweep pass accidentally applied the harness DEFAULT regime (THINBBA to BA 40 @cyc2) and showed large
jl-vs-live TPA divergences. Isolating on stand 31362831010690 (2225 TPA QMD 3.9 pre-thin), the 2028 post-thin row:
  LIVE 138 / QMD 7.4    jl(AUTOES) 652 / QMD 3.4    jl(NOAUTOES) 12 / QMD 25.0
Live sits BETWEEN jl's two modes ⇒ TWO distinct effects, BOTH management-triggered:
1. **THINBBA-from-below OVER-THINS on dense stands:** jl+NOAUTOES keeps only 12 huge trees (QMD 25) where live
   keeps 138 (QMD 7.4) at the same residual BA~40. jl's thin sort direction is correct (cuts.jl:526, −DBH from
   below), so this is a residual-BA-target / stopping-rule difference, not a wrong-end thin. NEEDS VERIFICATION —
   could be a real over-thin bug OR a dense-QMD-3.9 edge / harness artifact (extreme 182→40 BA thin). Check vs a
   clean simple-stand THINBBA and vs the examples/thinba validation.
2. **AUTOES over-establishes AFTER canopy-opening thin:** with AUTOES on, jl 12→652 (adds ~640 regen) where live
   adds far less (→138). This is the #143 AUTOES class, but POST-DISTURBANCE — the no-thin (closed-canopy) sweep
   showed NO over-establishment, so jl's AUTOES over-produces specifically after the canopy opens. Distinct from the
   bare-stand #143 (fixed/cornered); a real candidate follow-up.
Both are under an ARTIFICIAL harness thin, not real FIA management — flagged as leads, not confirmed bugs.

## Full-cluster extension (BM/TT/UT/CR, pure-growth 5cyc) — 2026-08-12
Completed the fresh sweep across the remaining western variants (KT has 0 FIA stands, native-key validated):
- **BM** (23 treed): 16/23 strict bit-exact-or-cornered (TPA/BA/QMD within 3%); the other 7 are all small ~1-8% BA
  straddles (accepted). 0 crashes.
- **CR** (7 treed): 6/7 (lone residual QMD Δ8%). 0 crashes.
- **UT** (7 treed): all within ~4-9% BA = the accepted UT compounding-BA tail (tolerance 3% too tight). 0 crashes.
- **TT** (13 treed): 3/13 strict; TWO LARGE divergences — 533757478126144 (TPA 4055→jl 7084) + 1629326355290487
  (BA 19→jl 43) — both DENSE small-tree stands = the KNOWN #158 TT small-tree-regent dense-stand gap (already
  documented; not new, not a regression). The other TT DIVs are ~3-20% BA straddles. 0 crashes.
⇒ WHOLE WESTERN CLUSTER (IE/EM/CI/BM/CR/UT/TT + KT-native + BC-metric) re-validated FRESH on real FIA data
  post-this-session's-fixes: bit-exact-or-cornered, ZERO jl crashes across ~90 stands, no regression. The only
  notable multi-cycle divergence is the pre-existing, documented TT #158 dense small-tree gap. The mission's
  "full FVS-ready FIA sweep" co-goal is satisfied for the current code.

## TT #158 sweep-divergence CHARACTERIZED (2026-08-12) — confirmed the documented gap, not new
Dumped the 2 TT sweep divergences full multi-cycle (FVStt_clean vs jl, no-thin):
- 533757478126144 (9549 TPA QMD 1.1 = dense sub-1" seedlings): cyc0 bit-exact; jl UNDER-mortalizes the cohort
  (2047 TPA 7084 vs live 4055, ~1.7×) + UNDER-grows QMD (2.3 vs 3.1).
- 1629326355290487 (5325 TPA QMD 0.1 = bare sub-inch seedlings): jl OVER-grows sub-1" DBH EARLY (2030 QMD 1.2 vs
  live 0.8) then retains the excess (2050 TPA 4880 vs 4186).
⇒ EXACTLY the documented TT #158 (sub-1" over-growth via the SMDGF/regent small-tree DG). The fresh full-cluster
  sweep found NO NEW bug — the one divergence is the one documented gap. On NORMAL stands the SMDGF residual is
  small (~7%, TT audit "bit-close" i34 0.637/0.646); on DENSE SUB-1" SEEDLING stands the tiny per-tree DG diff
  COMPOUNDS (self-thin timing is knife-edge — the #140-BM hyper-sensitivity class).
OPEN QUESTION (bug vs cornered): is #158 a real SMDGF-single-step-suppression MODEL gap (fixable port) or a
  hyper-sensitive dense-seedling realization straddle (cornered, like #140 BM)? RESOLVE via FVStt_g16 (exists,
  /workspace/.ttwork/FVStt_g16): instrument live's sub-1" per-tree DG on 1629326355290487 @2030 (the over-growth
  point, NOTRIPLE) vs jl's SMDGF/regent DBH assignment — if the DETERMINISTIC per-tree DG differs = real gap to
  port; if it matches = accepted straddle. This is the priority genuinely-open GROWTH item (real, on-FIA, not gated).

## ★ TT #158 bug-vs-cornered SETTLED (2026-08-12): REAL DETERMINISTIC BUG (via NOTRIPLE)
Re-ran the sub-1" reproducers with NOTRIPLE (tripling OFF ⇒ no per-triple ZZRAN — pure deterministic path):
- **1629326355290487** (sp108/LODGEPOLE, dense sub-inch, cyc0 QMD 0.1 BIT-EXACT jl==live): the divergence PERSISTS
  and GROWS deterministically — 2030 BA jl 9 vs live 27 (~3× UNDER), QMD 0.6 vs 1.0; by 2070 TPA jl 4406 vs live
  1724 (2.5× over-retained), QMD 2.6 vs 5.2. jl UNDER-grows the sub-1" LP seedlings ⇒ they never cross the
  self-thin thresholds ⇒ jl retains 2.5× too many. NOT an RNG/tripling straddle (NOTRIPLE removes tripling and the
  gap remains) ⇒ a REAL deterministic small-tree-DG bug in TT's smhtgf/smdgf/regent for sub-1" seedlings.
- 533757478126144 is closer under NOTRIPLE (jl BA 254 vs 234 @2067) — milder.
⇒ VERDICT: the one non-cornered growth divergence on real FIA data is a CONFIRMED REAL BUG (TT sub-1" small-tree DG
  under-grows dense LP/sp108 seedlings), NOT the accepted straddle. Note the SIGN: jl UNDER-grows here (opp. to the
  memory's "aspen OVER-grows #158" — likely a DIFFERENT species/path than the aspen case, or the aspen framing was
  tripling-flipped). NEXT (fix): instrument FVStt_g16 (/workspace/.ttwork/FVStt_g16) sub-1" per-tree DG for sp108 on
  1629326355290487 @2030 (NOTRIPLE) vs jl's teton/regent.jl small_tree path — find where jl under-assigns DBH/HTGR
  to sub-1" LP seedlings. This is the PRIORITY genuinely-open GROWTH bug (real, on-FIA, deterministic, not gated).

## ★★ TT #198 ROOT-CAUSED (2026-08-12): crown-model density excludes recent-mortality trees
Reproducer 1629326355290487 (TT sp108/LP) is a POST-MORTALITY stand: 4 LIVE seedlings (D0.1) + 89 DEAD/mortality
large trees (D5-9.8, tpa 6 each) — jl loads all (t.n=4 live, t.ndead=89). MEASURED chain:
- jl seedling CR 59/45/94/73 vs live (.trl) 58/17/15/31 — jl OVER-assigns CR to 3 of 4.
- NOT the DBH-rank: the .trl %-TILE (63/10/100/97) matches ranking over the 4 LIVE trees (both jl+live); same
  %-TILE gives different CR (live %-TILE-100 seedling CR 15 vs jl 94).
- ROOT = the Weibull crown DENSITY terms (teton/crown.jl:74-87): jl relative_density=5.3 / crown_sdi~0 (computed
  over ONLY the 4 live seedlings, BA 0.29) ⇒ scale=1−0.00167·(5.3−100)→capped 1.0 (NO suppression) ⇒ HIGH crnew.
  Live includes the 89 dead trees (high SDI/relden) ⇒ scale<1 + higher acrnew ⇒ SUPPRESSED crnew (low CR).
- The inflated CR → _tt_smhtgf BETA2·CR → seedling HEIGHT over-growth (2070 TopHt jl 34 vs 22) → tall/thin/low-BA
  → self-thin never fires → 2.5× TPA over-retention (4406 vs 1724).
⇒ FIX (task #198): include the recent-mortality trees in the TT crown-model density (relative_density + the
  crown_sdi passed to crown_ratio_update!), matching live cratet.f. VERIFY vs cratet.f (does it use mortality-
  inclusive SDIAC/RELDEN?); after fix jl seedling CR should → 58/17/15/31 = live and NOTRIPLE .sum → TopHt 9/QMD
  1.0 @2030. Non-regress ttt01 + the fresh FIA sweep. CAVEAT: mortality-inclusive density may touch other variants'
  crown/density — scope to the crown dub if needed. This is the ONE real growth bug the fresh sweep surfaced,
  now mechanistically root-caused (measure-don't-infer, doctrine #2).
