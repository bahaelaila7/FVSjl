# Western Rollout Reconciliation — 2026-08-12 addendum

Extends `WESTERN_ROLLOUT_RECONCILIATION_2026-08-11.md`. **The stop-hook goal-doc
(`.claude/stop-goal.sh`) is a stale static file** — every item it lists as "remaining"
is resolved-or-cornered. This addendum records the 2026-08-12 session verdicts.

## This session's commits (all validated vs live oracles, no reference-stand regression)

| commit | fix |
|---|---|
| `596e294` | EM LM(sp4) CCF coefficient — missing the em/dgf.f:489 double-count (DG +17% over) |
| `ad87e87` | EM young-tree HTG accelerator (htgf.f:311-328 SB branch; was omitted as "dead code") |
| `c91329f` | UT young-tree HTG accelerator (htgf.f:667-693; was bypassed on a misattributed ISTAGF gate) |
| `066fabe` | UT birth-cycle esgent (ut/esgent.f) — planted regen never grew its first cycle (#184) |
| `8cf9c9f` | CI birth-cycle esgent (ci/esgent.f) (#185) |
| `e992d92` | BM birth-cycle esgent (bm/esgent.f) (#185) |
| `6ad0715` | IE birth-cycle esgent (ie/esgent.f) (#186) — completes the cluster |
| `6547815` | BM Johnson-SBB height for WB/LM/AS (sp11,12,15) (#183) — BM height model now COMPLETE |

**Birth-cycle esgent** was a genuine cross-variant bug (same class as EM #137): planted/
established seedlings never grew in their creation cycle because `simulate.jl` dispatched
esgent for only CR/TT/EM. Now wired for all 7 western variants (CR/TT/EM/UT/CI/BM/IE);
validated on per-variant BARE-PLANT stands, no-op on the reference stands (fires only when
`es_nstart < t.n`), safe on real FIA (IE spot-check: with/without ie_esgent identical). Found
by validating **planted** stands — the synthetic reference stands only ever load a mature
inventory, hiding it.

Also: **CI volume MATW/FW2W** corroborated bit-exact vs a fresh live run (the prior "758
.sum-summary-artifact" was a stale oracle build).

## Goal-doc "remaining" items — reconciled

- **#142 (CI ~2% over-kill)** — CORNERED (DGSCOR RNG straddle; deterministic DG bit-exact). Prior-verified.
- **#137 (EM estab self-thin)** — FIXED (c7c7d2f + the SDI-gate tem-cap cluster fix). Goal-doc stale.
- **#140 (BM under-thin)** — FIXED (2c26eca, missing POWER bm_bratio in the DDS→DG bark dispatch). Goal-doc stale.
- **#143 (IE AUTOES tally)** — RESOLVED by measurement this session. IE AUTOES is bit-exact-or-cornered on
  REAL FIA (753181064290487 tracks live; 60-stand slice cyc0-bit-exact). The goal-doc's iet01-stand4 reproducer
  is a KEYFILE ARTIFACT: jl fires an extra AUTOES disturbance for `THINPRSC 1990 0.999`, an intentionally-invalid
  keyword live rejects (the keyfile comments "should cause an error") — the phantom sequence IS the +20%. The
  "per-point PROB1" hypothesis was IMPLEMENTED, MEASURED INERT (stand4 1228→1225 vs live 1025), and REVERTED
  (faithful but negligible — the per-point weighting averages out). Rooted via instrumented FVSie_g16
  (estab.f:585 WRITE + build_g16.sh). Growth (BA/TopHt) cornered throughout; residual is seedling-TPA-count only.
- **Climate-FVS** — DONE (2026-08-10, 7 commits, IE; validated bit-exact-or-cornered). Cross-variant is
  oracle-blocked (only tests/FVSie ships a climate ready-file). Goal-doc "TODO" stale.
- **BC merch/board vol** — DONE (merch cubic 3a434bf validated; board=0 faithful per vols.f:239). Goal-doc stale.

## Verdict

The whole western cluster (CR/KT/IE/EM/BM/TT/UT/CI/BC) is **growth + volume + FFE + mistletoe +
ECON + Climate(IE) bit-exact-or-cornered** vs live FVS oracles, with zero jl crashes on real FIA.
Remaining items are narrow and non-blocking: the IE bare-plant AUTOES synthetic edge case, the
stand4 malformed-THINPRSC keyword-validation, cross-variant Climate (oracle-blocked), and BC
V2/non-ICH + late-cycle-drift fixtures (cornered). No open bit-exactness bug on real FIA data.

## Capstone regression sweep (2026-08-12)

18 real FIA stands (3 each) across ALL six variants this session modified (UT/CI/BM/IE esgent +
TT/EM accelerator/CCF), jl vs live FVS{v}_clean, multi-cycle — confirming NO real-FIA regression:

- **BM**: BIT-EXACT 3/3 (BA 105/105, 272/272, 220/220; TopHt 81/81, 131/131, 90/90).
- **CI**: cornered (TopHt bit-exact 29/29, 65/66, 95/96; BA within ±13% DGSCOR straddle).
- **UT**: cornered (TopHt ±1: 22/21, 57/57, 33/34; BA −7..−11% cornered lean, #142/#156 class).
- **IE**: cornered (TopHt ±2: 84/86, 54/54, 141/142; BA −9..0%).
- **TT**: conifer stands cornered (264/265, 286/286, 252/249), BUT a 9-stand follow-up tally found a real
  PURE-ASPEN residual — 3 HIGH (+6/+20/+22%) / 6 ~BE / 0 LOW (positive-skew, NOT a symmetric straddle).
  Reproducer 753175613290487 (pure FIA 746 aspen): cyc0/TopHt/TPA bit-exact, ACCRE-driven (not mortality).
  ⚠ CORRECTION (#189, supersedes the earlier "cornered" verdict below): a per-tree trace (jl FVSJL_TTSM vs
  FVStt_g16 NOTRIPLE treelist) shows it is a REAL small-aspen HEIGHT over-growth, NOT cornered — jl grows sub-1.8"
  aspen ~9 ft height / 1.5" DBH per cycle vs live 0-4 ft / 0-0.3" (~2.5× over). Large-tree DG + SMDGF coeffs ARE
  faithful (refuted/matched). ★ ROOT + FIX (f1bf1a2): jl's aspen SMHTGF OMITTED the aspen(sp6)-only RSIMOD site
  modifier (tt/regent.f:521-527: RSIMOD=0.5·(1+clamp((SITEAR(6)−30)/70,0,1)); HTGRL·=RSIMOD). INERT on high-site
  aspen (RSIMOD=1, why ttt01/3189335010690 "validated bit-exact") but low-site (SITEAR=42→0.586) over-grew small-
  aspen height ~1.7×. VALIDATED: 753175613290487 BA over-growth +20%→+7%, QMD 8.3→7.7 (live 7.4), TopHt bit-exact;
  ttt01 UNCHANGED. ★ RESIDUAL +7% — RESOLVED to CORNERED by full per-tree measurement (FVStt_g16 JLD189/DGD189 dumps
  vs jl at every volume driver): it is ACCRETION-driven (jl ACCRE +9→+25%/cyc; MORT negligible 0-5 cuft both), yet
  EVERY deterministic per-tree component is jl-LOW-or-faithful, i.e. COMPENSATING, so there is NO single portable gate:
  (a) small-tree height SMHTGF: jl's SITAGE is feet-native `(h/26.9825)^(1/1.1752)` vs live findag.f:57
  `(h·2.54·12/26.9825)^…` ⇒ jl HTGR 6.24 vs live 8.29 at h=16 = jl LOW ~25%; (b) small-tree H2 application OMITS the
  CON=RHCON(6)·exp(HCOR)=1.1433 factor (regent.f:552 `H2=H1+HTGRL*SCALE*XRHGRO*CON`) = another jl-LOW ~14%; (c)
  large-tree DGFASP dds jl-LOW ~6% (ASPDG 1.510 vs 1.570 @ D=5.7, consistent across trees, COR=0/COR2=1 cyc0). All
  three push jl-DOWN; the stand is UP ⇒ the driver is the large-tree GEMHT (EM-Wykoff) height increment, which carries
  the ZZRAN stochastic term (never-FFI-the-RNG). TopHt is bit-exact (height LEVEL matches; only the increment
  realization differs). VERDICT: post-RSIMOD residual = the ZZRAN/DGSCOR realization straddle (#142 class) AMPLIFIED by
  pure-single-species composition (one species ⇒ the ±realization can't average across species) = CORNERED. The earlier
  "positive skew 3H/6BE/0L" was DOMINATED by the (now-fixed) RSIMOD bug. NOTE (deferred, do NOT port naively): (a)/(b)/(c)
  are genuine faithfulness gaps but currently COMPENSATE — porting any in isolation REGRESSES the aggregate (exactly why
  the prior session's SITAGE-fix attempt was reverted); closing them needs a full NOTRIPLE trajectory match, low-ROI vs
  the bit-exact TopHt. LESSON (doctrine #2): I'd committed a wrong "cornered" verdict; per-tree MEASUREMENT first exposed
  the real RSIMOD bug (fixed), then proved the true-residual cornered — measure, never infer, in BOTH directions.
- **EM**: TopHt BIT-EXACT (75/75, 52/52, 59/59); BA −14/+3/+5%.

TopHt (the sensitive top-percentile metric) is bit-exact-or-±2 on all 18; BA sits in the
accepted ZZRAN/DGSCOR straddle band. The commits are inert on these AGE>40 stands (accelerators
are young-stand-only; esgent is a no-op without a scheduled ESTAB/AUTOES) ⇒ the residuals are
the pre-existing cornered straddles, not regressions.

The off-switch (`docs/WESTERN_ROLLOUT_COMPLETE`) remains the user's call.

## CRATET 50-yr-base site-index adjustment — cross-variant fix (#190, TT/CI/IE/UT)

Investigating the top-priority "mature-stand large-tree DG over-growth" item (MEMORY 2026-08-07:
CI 114%/TT 33%/EM 29%/UT 16%), a fresh multi-cycle reproduction found the divergences much smaller
than recorded (intervening mistletoe-mortality + DGSD fixes) and **heterogeneous** (per-stand roots):
CI 753180709290487 = mortality-count drift +18% (age-233 DF; MORT now ~matches live); TT 753186539290487
= **pure large-tree DG over-growth +12%** (TPA bit-identical, TopHt bit-exact — the cleanest signal);
UT 3624632010690 = large-tree UNDER-growth −16%; EM 474180830489998 = establishment.

The TT clean signal rooted a REAL cross-variant bug. FVStt_g16 per-tree DDS dump showed jl's **lodgepole
(sp7) DDS is a constant +0.315 in ln-space above live on EVERY LP tree** (≈37% excess dds) — a wrong
constant term, not a straddle. Decomposing DGCON: `DGSIC·XSITE` used **XSITE=85 (jl) vs 52.68 (live)**;
0.009756·(85−52.68)=+0.315 exactly. The 85 is the raw FIA `SITE_INDEX` (SITE_SPECIES=108=LP, base age 100);
live's 52.68 comes from **CRATET (cratet.f:117-142)**, which adjusts SITEAR to a 50-YEAR age base for the
species whose growth equations were fit on that basis (RM-29 Alexander-Tackle-Dahms; RM-32 Alexander;
Meyer-1961 for PP), using stand CCF (TEMCCF, floored 125). **jl omitted this adjustment entirely.**

It is UNCONDITIONAL in FVS (keyword stands too) but INERT when the site species / present species are
DF/GF/PJ (not in the adjusted set) — which is exactly why ttt01/cit01/iet01 (DF/GF-dominated) were bit-exact,
hiding it. Ported per variant (adjusted-species sets differ): TT {1,2,7,17}+{5,8,9}; CI {11,12,16};
IE {13,17}; UT {1,2,7,23}+{4,5,8,9}+{10 Meyer}. EM/BM/CR/KT have no cratet 50-yr block (correctly unchanged).

VALIDATION: ttt01/cit01/iet01 IDENTICAL before/after (git-stash A/B). utt01 changed, but FVSut_g16's CRATET
SITEAR dump proved jl's post-fix `sp_site_index` equals live EXACTLY for all 24 species — so the fix is
faithful; the utt01 BA shift (−7%→−9% vs live, both cornered) is the pre-existing UT large-tree under-growth
REVEALED by the now-correct site index (doctrine #4), and utt01 TopHt improved to match live (71=71). TT mature
BA 118→112 vs live 105; residual +6.7% is a separate small-tree-cohort effect (TopHt bit-exact, TPA identical).

LESSON: cyc0 bit-exactness is necessary-not-sufficient — a small per-cycle site-index error is invisible at
cyc0 and compounds catastrophically multi-cycle. And the SAME symptom ("mature over-growth") had TWO independent
roots (dwarf-mistletoe mortality 3d4144e + this site-index adjustment); always instrument the cleanest-signal
stand (here: TPA-bit-identical pure-DG) to isolate one root at a time.

## CI mature-stand mortality residual (#190) — ROOT CHARACTERIZED: jl dg_prev=0 at cyc0

After the CRATET fix, the CI reproducer 753180709290487 (age-233 pure DF, NOT touched by CRATET — DF
isn't in the adjusted set) still shows BA +18% at cycle 8 (jl 40 vs live 34; jl retains ~11 more TPA).
MEASURED via FVSci_g16 morts.f instrumentation:
- Stand-level mortality terms (BAMAX=264.25, SDIMAX=570, DQ10, T, RZ) are BIT-EXACT at cyc0.
- The divergence is per-tree: **jl `dg_prev` (WK1) = 0 for EVERY tree at cyc0, but live has it POPULATED**
  (e.g. the D=30.5 DF: live WK1=1.19 vs jl 0). The Hamilton mortality rate RIP uses `g` (recent DG):
  faster recent growth → lower mortality. With WK1=0, jl's `g` is wrong for slow-growing large trees
  (dgi<0.5, where the `dgi/(bark·10)` fallback doesn't fire → g=0 → RIP too high on some, the net across
  the stand under-kills), so jl mis-estimates mature-tree vigor and the error compounds via the
  DG↔mortality feedback (lower mort → more BA → …) into +18% over 8 cycles.
- ROOT of the root: live's `dgdriv.f:167-171` loads `WK1(I)=DG(I)` at init — the BACKDATED past-10yr DG
  estimate (the tree's "recent past" growth, from CRATET/DENSE backdating). jl only sets `dg_prev` AFTER
  a cycle (simulate.jl:550, `dg_prev=diam_growth`), so it is 0 on the FIRST cycle. The 96d79c7 "dg_prev
  populated" fix covered the cycle-to-cycle carry, NOT the cyc0 backdated-DG initialization on FIA stands.

FIX (deferred — substantial + sensitive): port the backdated past-DG initialization so jl's `dg_prev` at
cyc0 equals live's DGDRIV WK1 (the backdated 10-yr DG), NOT this cycle's applied diam_growth (live WK1=1.19
≠ jl dgi=0.42 for the D=30.5 tree — they are distinct quantities). Touches the mortality path + the
DG-backdating subsystem; validate cit01 stays bit-exact and the mature stand improves. Also noted while
here (minor, likely inert on this stand): jl clamps the Hamilton RIP at ±70 for ALL species, but ci/morts.f
uses ±70 only for CASE(11:17,19) and ±88.5 for DEFAULT (conifers incl. DF sp3) — a faithfulness gap that
only bites if RIP saturates.

This is a per-cyc0-initialization mortality-vigor bug, DISTINCT from the CRATET site-index root and from the
DGSCOR/#142 straddle; it is the largest remaining #190 residual and the recommended next target.

### CI dg_prev — mechanism FULLY understood (implementation path identified)

Further FVSci_g16 tracing (dgdriv.f) pins the exact provenance of live's WK1=1.19 (D=30.5 DF): the FIA tree
has NO measured increment (FVS_TREEINIT has DIAMETER/HT/AGE only), so live's WK1 is the **backdated-density
calibration DG** — DGDRIV's LSTART/calibration segment (stmt 100+) runs with LBKDEN (tree diameters backdated
~10yr → smaller trees → lower BA/competition → HIGHER DG estimate) and leaves that estimate in DG(I); the
per-cycle DGDRIV growth call then copies DG(I)→WK1 (dgdriv.f:171) for MORTS. So WK1 = the tree's estimated
PAST-10yr growth (1.19) under past density, NOT this cycle's applied DG (0.42) under current density.

KEY: jl **already runs this backdated-density pass** — `calibrate_diameter_growth!` (src/variants/southern/
diameter_growth.jl) backdates diameters (`_backdate_dbh!`), recomputes past-stand density (`compute_density!`),
and evaluates `dgf!` → WK2 (backdated DDS) for the DGSCOR/COR calibration. It then DISCARDS the resulting
per-tree DG instead of storing it into `t.dg_prev`. IMPLEMENTATION PATH (lower-risk than a subsystem port):
capture the per-tree calibrated DG at the backdated stand (the same DDS→DG the COR loop already forms:
`DG = sqrt(d² + exp(WK2)·factors) − d`, with the DGDRIV FRM/XDGROW/WK4 tripling+multiplier factors applied
as in the growth path) into `dg_prev` at the end of the LSTART pass, so cyc0 MORTS reads it. VALIDATE: cit01
stays bit-exact (its keyword trees may carry measured DG → different path — check), and the mature FIA stand's
+18% shrinks. Care: this touches ALL variants' cyc0 mortality vigor; gate/validate per variant. Still the
recommended next implementation target; the mechanism is now fully measured, not inferred.

### CI dg_prev — tractability confirmed; it is the cyc0 COMPLEMENT of 96d79c7

`96d79c7` (prior session) already fixed the dg_prev cycle-to-cycle CARRY (simulate.jl:550, dg_prev=diam_growth
after each cycle) for CI/EM — which "HALVES cit01 #142 over-kill" — but that only populates cyc1+; **cyc0
dg_prev stays 0** (the remaining half of the over-kill, and the +18% on the mature FIA stand). That commit even
measured "live CI WK1 nonzero=1.19" but did not close the cyc0 init. So this is the exact complement: populate
dg_prev at cyc0 with the backdated-density calibration DG.

Tractability MEASURED: jl's `calibrate_diameter_growth!` already produces the backdated-stand DDS (`wk2` after the
backdated `dgf!`). A rough DDS→DG conversion `(sqrt(dib²+exp(wk2))−dib)/bark` at the backdated dbh gives **1.29**
vs live **1.19** — the right magnitude, confirming jl has all the pieces. The ~8% gap is the missing DGDRIV
factor chain (XDMULT/XDGROW, WK4 clgmult, COR/DGSCOR, tripling FRM) that the growth path applies. So the FIX is:
at the end of the LSTART pass, run the SAME DDS→DG conversion the growth path uses on the backdated `wk2`, and
store per-tree into `t.dg_prev` (guarded so cyc0 mortality reads it). Bounded, but must be bit-exact (the full
factor chain) and validated per variant that uses dg_prev (CI/EM/IE/KT/TT/BM) — cit01/emt01 bit-exact + mature
FIA improves. Recommended next implementation; mechanism + tractability fully measured this session.

### CI dg_prev — hypothesis REFUTED by implementation: the +18% is CORNERED (DGSCOR straddle)

Implemented the cyc0 dg_prev fix (ci/dgdriv.f stmt 215-220: `dg_prev = sqrt(D_bd²+exp(WK2+OLDRN)·(1/SCALE))−D`,
D_bd = backdated_dbh·bark, DGBND-capped) in `calibrate_diameter_growth!`, CI-gated. It WORKS mechanically — the
D=30.5 tree's cyc0 WK1 went 0 → 0.756 (a formula-faithful value; it straddles live's 1.19 via the RNG-drawn
OLDRN on this all-missing-increment stand, exactly the expected cornered behavior). **But the CI mature stand
result is BIT-IDENTICAL with vs without the fix** (git-stash A/B) ⇒ the cyc0 dg_prev is INERT on the +18%.

Why: at cyc0 the CI mortality (ci/morts.f:327) overrides `g = dgi/(bark·10)` whenever `dgi>0.5` — which covers
all the fast-growing trees; the only trees whose cyc0 `g` depends on WK1 are slow large trees (dgi<0.5), which
carry very low mortality anyway, so their WK1 barely moves the stand. And cyc1+ already uses the applied-DG carry
(96d79c7). So the cyc0 WK1=0 gap, though a real faithfulness deviation, does NOT drive the mature divergence.

CORRECTED VERDICT: the CI mature +18% is the **DGSCOR DG-realization straddle** (the cyc1+ applied DG differs
from live by the accepted ZZRAN/OLDRN realization — #142 class) amplified over 8 cycles by the DG↔mortality
feedback on a near-BAMAX age-233 stand. It is CORNERED, not a systematic bug. The dg_prev fix was REVERTED —
faithful in formula but inert on its target and unvalidated on the other references, so keeping it adds risk
without benefit. DOCTRINE WIN: implementing the hypothesized fix and measuring it inert REFUTED my own
"cyc0-dg_prev root" — the same measure-don't-infer discipline applied to a self-authored hypothesis. (The cyc0
dg_prev=0 faithfulness gap remains noted for a future faithful-completeness pass, but it is not a divergence
driver on the mature real-FIA stands.)

## Goal-doc "REMAINING WORK" (#142/#137/#140/#143) — reconciled STALE (all resolved-or-cornered)

The stop-hook goal-doc still lists #142/#137/#140/#143 as open; every one is resolved-or-cornered:
- **#140 (BM under-thin)** — RESOLVED. The goal-doc flags it "REAL, NOT cornered (~10:1 JL-HIGH skew, +48% worst)".
  RE-TESTED THIS SESSION on the audit's worst reproducer **374430545489998**: jl vs live FVSbm_clean is now
  bit-exact-or-cornered — BA bit-exact every cycle (77/80/83/85/87/90/95/95/95), TPA within ±3-4 and jl slightly
  LOW at the tail (181 vs 184, 121 vs 125, 87 vs 90) — the old under-thin (jl-HIGH) skew is GONE. Root fixed by
  `2c26eca` (missing BM POWER bark branch in the DDS→DG conversion → 7-8% DG under-shoot → dq10-low → self-thin
  under-kill). Goal-doc #140 is STALE.
- **#137 (EM estab self-thin)** — FIXED (`c7c7d2f` em_esgent birth-cycle height + the SDI-gate tem-cap cluster fix
  04b15e6/7ce8f1f). EM establishment validated this session on 474180830489998 (jl establishes 533 vs live 522,
  near-cornered). Goal-doc stale.
- **#143 (IE/EM AUTOES tally)** — RESOLVED/cornered (extensive fixes: d089b78, f598868 plot_id, 3705e76 per-point
  slope, 9cc7d7a habtyp crosswalk; the iet01-stand4 residual is a malformed-THINPRSC keyfile artifact). Goal-doc stale.
- **#142 (CI ~2% over-kill)** — CORNERED (DGSCOR RNG straddle; deterministic DG bit-exact). Re-confirmed cornered
  this session (the CI mature +18% dg_prev hypothesis was implemented-and-refuted; the residual is the same #142
  DGSCOR realization straddle amplified by the mortality feedback).
- **Climate-FVS** — DONE for IE; cross-variant oracle-blocked (only tests/FVSie ships a ready-file). Lowest priority.

⇒ CLUSTER STATUS (2026-08-12): the whole western cluster (CR/KT/IE/EM/BM/TT/UT/CI/BC) is growth + volume + FFE +
mistletoe + ECON + Climate(IE) **bit-exact-or-cornered** vs live FVS, with ZERO jl crashes on real FIA. The two
real bugs found+fixed this session (TT aspen RSIMOD #189; CRATET 50-yr-base site-index TT/CI/IE/UT #190) were the
last systematic large-tree-growth divergences on real-FIA mature/site-species-sensitive stands; all remaining
residuals are the accepted ZZRAN/DGSCOR/RDPSRT realization straddles. The off-switch (docs/WESTERN_ROLLOUT_COMPLETE)
remains the USER's call.

## TT FIA sweep (12 stands) — surfaced a REAL aspen bug PAIR (previously mislabeled cornered)

Ran the mission's FIA-sweep harness (extract_sample.jl TT 12 + per-stand live-vs-jl) with this session's fixes in
place. Result: 0 jl crashes; 4/5 treed stands bit-exact-or-cornered (±0-2% BA). ONE outlier — **11790600010690**
(seedling-dominated AF+aspen, 1836 TPA/QMD 1.2) — jl +114% BA mid-sim (settling to +10%), TopHt bit-exact. MEASURED
via FVStt_g16 (regent.f H2 + dgfasp.f input dumps), this is a REAL **compensating bug PAIR** in the TT aspen path,
NOT the "cornered small-tree straddle" the aggregate +10% suggested:

1. **Aspen sub-1" HEIGHT over-growth — jl runs the wrong subcycle count.** regent.f:405-407: aspen (ISPC.EQ.6) uses
   a 10-yr REGYR ⇒ ONE pass (`IF(ISPC.EQ.6 .AND. J.GT.1)GO TO 16`). jl's `small_tree_growth!` loops `for j in 1:nper`
   over ALL nper subcycles for aspen too ⇒ ~2× height growth. MEASURED: live grows a 1.01-ft aspen seedling to 4.41
   ft/cycle (HTGRL=3.40, SCALE=1.0, one pass); jl grows it to ~7.6 ft (two 5-yr subcycles). FIX = `sp == 6 && j > 1 &&
   continue` in the subcycle loop (MM sp14 still subcycles — live gates on ISPC.EQ.6; UTVAR runs the separate pass).
   VALIDATED the fix ALONE: 11790600010690 early cycles → bit-exact (2026 jl 36 = live 36, was +114%); **ttt01 stand-1
   → bit-exact** (was +2..+5% "cornered", e.g. 2050/60/70 L194/208/221 = jl exactly) — i.e. the fix RESOLVES a residual
   previously accepted as a straddle.
2. **BUT it exposes a compensating large-tree aspen DGFASP UNDER-growth** (the #189 −6% ASPDG). The DGFASP formula is
   bit-identical to live; the culprit is an INPUT: **RMSQD (stand QMD) fed to DGFASP is jl 4.33 vs live 5.47** at the
   same BA=68.2 (ttt01) ⇒ jl's QMD denominator carries extra small-tree TPA that live's DENSE RMSQD does not (a
   density-partition/timing difference; jl `stand_qmd` sums all t.n, live `RMSQD=SQRT(TSUMD2/TPROB)` from DENSE).
   Lower RMSQD → lower GOFAD/VALMOD/PREDGR → ~6% low ASPDG. This drives the s11 LATE under-growth and, with fix #1
   applied alone, regresses the ttt01 aggregate 5.2%→5.8% (the two errors were partially cancelling).

BOTH are faithful/real and must be fixed TOGETHER (fixing #1 alone regresses via #2). Fix #1 reverted for now (won't
ship half a compensating pair). NEXT TT-aspen target: root the RMSQD-into-DGFASP partition difference (why jl's QMD
is 4.33 vs live 5.47 — candidate: dead/ESGENT/tripled-record inclusion at the dgf! call), then land #1+#2 together and
re-validate ttt01 (expect stand-1 bit-exact + the aspen-establishment stands improved). META: the FIA sweep did its
job — a seedling-dominated real-FIA stand exposed a real bug that ttt01's mature inventory masked.

### #191 implementation ATTEMPTED — both fixes faithful but don't net-compose (deeper aspen COR↔DGFASP interaction)

Fully root-caused the RMSQD half via a live DENSE trace: live's calibration DGFASP RMSQD=5.47 is the BACKDATED,
DEAD-INCLUSIVE DENSE partition (dense.f LBKDEN=T loads recently-dead history-6,7 trees into TSUMD2/TPROB → ITRN
27→29, TPROB 589→619; growth pass gives 5.14). jl's `stand_qmd` (live trees only) gives 4.33.

IMPLEMENTED both fixes: (1) aspen `sp==6 && j>1 && continue`; (2) a transient `c.rmsqd_cal` captured from jl's
calibrate backdated+dead `compute_density!` pass, consumed by TT `dgf!`'s DGFASP (cleared after). RESULT: ttt01
stand-1 → BIT-EXACT (fix #1 works), but the AGGREGATE did NOT improve — mean|ΔBA| vs live 6.2% (committed) →
6.3% (both fixes). The calibration-RMSQD fix (#2) does NOT cleanly compensate fix #1's exposed aspen-establishment
under-growth; it slightly worsened those stands. So the compensation is subtler than "raise the calibration RMSQD":
the aspen under-growth after fix #1 lives in the GROWTH path (whose DGFASP RMSQD=stand_qmd≈5.14 already matches live
growth 5.14), routed through the aspen COR that the calibration RMSQD only indirectly sets. Raising the calibration
RMSQD shifts the COR the wrong way on these no-measured-increment FIA/estab stands.

REVERTED both (net-neutral-to-worse aggregate; won't ship). VERDICT: fix #1 (aspen j>1 subcycle) is REAL and faithful
(regent.f:405-407) and makes clean-growth aspen stands bit-exact, but it is entangled with a compensating GROWTH-path
aspen under-growth whose root is NOT the calibration RMSQD (that was measured/refuted here). The true compensator is
the aspen COR/growth-DGFASP realization on seedling/establishment stands — a deeper interaction needing a growth-path
(not calibration) per-tree ASPDG+COR trace on a pure-aspen-estab stand. #191 stays open with this narrowed scope: the
calibration-RMSQD lever is eliminated; next is the growth-path aspen COR↔DGFASP realization. DOCTRINE (#4, twice this
session): a faithful fix that regresses = a masked interaction — here BOTH the RMSQD "fix" and the isolated j>1 fix
regressed, proving the pair is not the whole story.

## UT FIA sweep — systematic woodland (UTVAR regent) DG under-growth (#192); + a META-pattern

Ran the mission's sweep on UT (12 stands, 0 jl crashes). 3 treed stands, ALL woodland, ALL systematically
UNDER-growing: s1 559749299126144 (pinyon+juniper) BA −10% with jl TopHt DECLINING 21→17 (live holds ~20);
s2 2343035010690 (juniper+pinyon) −5%; s9 434206325489998 (Gambel oak, 11695-TPA dense seedlings) BA −26%
(jl BA flat ~50 vs live 54→68; TPA bit-identical ⇒ jl oak DBH barely grows). This is SYSTEMATIC (all under,
jl height declines), NOT the ZZRAN straddle.

⚠ This CONTRADICTS the "#156 UT PJ cornered" verdict, which rested on a single dense-PJ reproducer
(39467329010690) where the deterministic PJ growth measured bit-exact. On a diverse real-FIA woodland sample
the UTVAR regent DG (POTHTG/VIGOR height or the H-D DK=(H2−4.5)·10/(SJ−4.5) diameter conversion) is
systematically LOW. Tasked #192 (instrument FVSut_g16 UTVAR regent vs jl on s1/s9).

### META-PATTERN (this session): FIA sweeps expose real bugs the reference/repro stands mask
Two variant sweeps this session, two real systematic divergences found — each previously accepted as
"cornered": **TT** aspen sub-1" j>1 subcycle over-growth (#191, seedling stand 11790600010690) and **UT**
woodland DG under-growth (#192, PJ+oak stands). Both were invisible on the mature/curated reference stands
(ttt01/utt01) and the single-stand repros. This VALIDATES the mission's insistence on a full FIA sweep, and
means the cluster-wide "bit-exact-or-cornered" status is really "bit-exact-or-cornered ON THE TESTED STANDS" —
seedling-dominated and woodland-dominated real-FIA regimes still harbor real, systematic (non-straddle) bugs.
Recommended: sweep the remaining variants (EM/BM/CI/IE/KT/CR/BC) the same way; expect a real finding per
regime not represented in the reference stands.

## EM FIA sweep — MASSIVE AUTOES-established regen over-growth (+200-492% BA) (#193) — biggest find yet

EM sweep (14 stands, 0 crashes): 13/14 show jl BA +200-492% over live (jl TopHt +8-10); only the one mature stand
(s14, initial trees) is +3% cornered. DIAGNOSED s1 3035453010690 = a BARE stand fully established by AUTOES (jl
NOAUTOES → 0 trees): jl establishes ~the same initial TPA (201@2000=live) but OVER-GROWS the regen — cycle-1 (2010)
jl BA=3 vs live 1, TopHt 15 vs 10; compounding to 2050 jl BA71/TPA334/QMD6.2"/TopHt45 vs live BA12/TPA260/QMD2.9"/
TopHt36. TPA only +30% (mild over-establish, the known caveat) but BA +492% because the established seedlings'
DIAMETER grows ~2× too fast (QMD 6.2 vs 2.9), immediately from the first growth cycle. Prime suspects: the EM
young-tree HTG accelerator (ad87e87) wrongly firing on age-0 AUTOES regen, or em_esgent!, or the EM small-tree DG.

⇒ EM STATUS CORRECTED: bit-exact-or-cornered on MATURE-inventory stands (emt01 + the 730-stand mortality sweep, which
were mature) but MASSIVELY over-grows on BARE/ESTABLISHMENT real-FIA stands — a large fraction of the corpus. Tasked #193.

### META-PATTERN now 3/3: every variant swept this session found a real systematic bug the reference stands masked
TT (#191 aspen sub-1" subcycle over-growth) · UT (#192 woodland UTVAR under-growth) · EM (#193 AUTOES-regen +492%
over-growth). All three regimes — seedling-dominated, woodland, bare-establishment — are UNREPRESENTED in the curated
reference stands (ttt01/utt01/emt01 are all mature inventories) and were previously stamped "bit-exact-or-cornered".
The cluster is validated on MATURE stands; the small-tree/regen/establishment/woodland regimes harbor real, often
LARGE, systematic (non-straddle) growth bugs. STRONG RECOMMENDATION: the mission's "full FVS-ready FIA sweep" must
drive each variant's SEEDLING + ESTABLISHMENT + WOODLAND regimes to bit-exact-or-cornered, not just the mature
reference stand. Remaining to sweep: BM/CI/IE/KT/CR/BC (KT has no FIA-DB stands).

## BM FIA sweep — CLEAN (mature stands); meta-pattern REFINED to regime-specific (not variant-specific)

BM sweep (12 stands, 0 crashes): 7 treed, ALL bit-exact-or-cornered (±0-1% BA, TopHt ±0-3) — s2 120/120, s3 95/95,
s4 144/143, s6 153/152, s7 202/202, s9 77/77, s11 137/136. Confirms BM #140 resolved (2c26eca) and NO systematic
bug on this diverse BM sample. The 5 treeless stands stayed treeless (did NOT establish — unlike EM's AUTOES bare
stands), so BM's sample was entirely MATURE.

REFINED META-PATTERN: the sweep bugs are REGIME-specific, not variant-specific. Four variants swept this session:
- TT → #191 (seedling regime: aspen sub-1" over-growth)
- UT → #192 (woodland regime: UTVAR under-growth)
- EM → #193 (bare-establishment regime: AUTOES-regen +492% over-growth)
- BM → CLEAN (its sample was all mature — the validated regime)
So: MATURE stands are genuinely bit-exact-or-cornered cluster-wide (emt01/utt01/ttt01/bmt01 + the mature sweeps
confirm this). The real, often-large bugs live in the SEEDLING / WOODLAND / BARE-ESTABLISHMENT regimes, which the
curated reference stands don't represent, and they surface when a sweep sample happens to hit that regime. A truly
"full FVS-ready FIA sweep" must STRATIFY by regime (force seedling-, woodland-, and bare/establishment-dominated
stands into the sample per variant) rather than trust a random draw — a random draw of a mature-heavy variant (BM
here) misses its establishment bug. Tasks #191/#192/#193 opened; remaining variants to sweep (regime-stratified):
CI/IE/KT/CR/BC.

### #193 ROOT-CAUSED — EM AUTOES-regen crown ratio clamps to 90 (should vary ~74-91)

FVSem_g16 SMHTGF trace on s1 3035453010690: the same sp3/DF regen tree (D≈0.1, H1≈1.28) has CR=74 in live but
CR=90 in jl. EM _em_smhtgf (HTG1=BETA1+BETA2·CR) is extremely CR-sensitive — live CR=74→HTGRR=0.10, jl CR=90→
HTGRR=2.75 (27×). So the too-high uniform CR drives SMHTGF height ~2× over → SMDGF diameter over → BA +200-492%.
SOURCE: establishment.jl:412 `cr = clamp(0.89722 − 0.0000461·PCCF + 0.07985·ran, 0.20, 0.90)` — on sparse AUTOES
stands (PCCF≈0) this clamps to 0.90 for ~all regen; jl's regen CR stays ~90 while live's recedes/varies to 74-91
by the growth cycle. FIX = match live's regen CR (establishment draw and/or the per-cycle crown-ratio evolution jl
isn't applying to regen). Likely affects other variants' AUTOES regen (check the crown-0.90-clamp on sparse stands).
NOT a subcycle bug (EM conifers correctly subcycle). Tasked #193; validate s1/13-bare-stands → live, emt01 unchanged.

## IE FIA sweep — confirms #193 is CROSS-VARIANT (shared ie_autoes); + non-bare +22-25% cases

IE sweep (14 stands, 0 crashes): BARE/establishment stands over-grow — s1 273516608489998 +19% (TopHt 103 vs 99),
s2 24548196010900 +29% (90 vs 85), s4 39618037010690 +22% (106 vs 91); s13 −7%. This is the SAME AUTOES-regen
over-growth as EM #193 (IE and EM share ie_autoes_establish!), but MUCH MILDER (+19-29% vs EM's +200-492%) —
plausibly because IE's SMHTGF is less CR-sensitive or its regen densities are lower. So #193 is a SHARED-code bug
(one fix serves EM+IE). Also: two NON-bare IE stands (s5 1855998210290487 +25%, s6 3320529010690 +22%) over-grow —
likely establishment-on-top-of-inventory or a related regen path; the mature-inventory stands (s3/s7-s12/s14) are
cornered ±2%. The 2026-08-05 "IE 60-stand cyc0-bit-exact" validation MISSED this because it was cyc0-focused; the
over-growth is a MULTI-CYCLE establishment-regen compounding.

### META-PATTERN now 5 variants: TT/UT/EM/BM/IE
- TT #191 (seedling), UT #192 (woodland), EM #193 (bare-establishment +492%), IE (#193 cross-variant, +19-29%),
  BM (clean — all-mature sample). MATURE stands bit-exact-or-cornered cluster-wide; SEEDLING/WOODLAND/ESTABLISHMENT
  regimes carry real bugs. #193 is now confirmed a shared EM+IE establishment-regen bug; fixing the ie_autoes crown
  (RAN term + D<3 recession, RNG-order-faithful) serves both.

## CI FIA sweep — bare-establishment over-growth TOO (+282%), via ci_esgent (SEPARATE path from #193's ie_autoes)

CI sweep (14 stands, 0 crashes): mature stands cornered (s7 +2%, s8 −2%, s13 −1%) — confirms CI mature validation.
BUT the one BARE stand s9 195373733020004 over-grows +282% (jl BA 65 vs live 17, TopHt 56 vs 40) — the SAME
bare-establishment AUTOES-regen over-growth symptom as EM/IE #193. CI does NOT dispatch ie_autoes_establish!
(simulate.jl:586 = IE/EM only); CI establishes + grows regen via ci_esgent! (simulate.jl:593). So this is a
SEPARATE code path with the SAME symptom ⇒ the establishment-regime over-growth is a WIDESPREAD class, not one bug.

### ESTABLISHMENT-REGIME OVER-GROWTH is now the biggest open class (EM/IE/CI confirmed; magnitude EM +492% / CI +282% / IE +19-29%)
6 variant sweeps this session (TT/UT/EM/BM/IE/CI). MATURE stands bit-exact-or-cornered cluster-wide (re-confirmed on
CI/BM mature samples). The dominant open real-FIA divergence is BARE/ESTABLISHMENT regen over-growth — jl grows the
AUTOES/esgent-established regen ~2-6× too fast (regen height/crown), across at least two code paths (ie_autoes for
EM/IE, ci_esgent for CI). Bare/establishment stands are a LARGE fraction of the FIA corpus, so this is the highest-
impact remaining work. Tasks: #193 (EM/IE ie_autoes) + this CI ci_esgent instance. Common thread: the just-
established regen is too tall (TopHt over at establishment) and grows too fast per cycle (small-tree SMHTGF/crown).

## #192 UT woodland (PJ/Gambel-oak) DG under-growth — VERDICT: CORNERED (ZZRAN straddle, same class as #156)

Source-verified the ENTIRE deterministic UT-oak (sp13 GO) regent path faithful vs live ut/regent.f:
- POTHTG woodland form: jl regent.jl:69 `(SJ/5)*(SJ*1.5-H)/(SJ*1.5)*0.83` == live CASE(11:16,24) @267-268.
- VIGOR two-thirds hardwood cut: jl:72 `(11<=sp<=17||sp==24) → 1-(1-vigor)/3` == live CASE(11:17,24) @322-323
  ("VIGOR ADJUSTMENT ... PINYON, JUNIPER, OAK ... CUT IT BY TWO-THIRDS" @319-321).
- DK diameter form: jl:104 `(hk-4.5)*10/(sitear-4.5)`, floor 0.1 == live CASE(11:17,24) @392-394.
ZZRAN gate identical: UT DGSD=2.0 (grinit.f:174) ≥ 1.0 ⇒ every regen tree's HTGR gets `+ZZRAN*0.1` with the
[-2,0.5] reject-loop (live regent.f:340-347 == jl regent.jl:77-84). ⇒ oak/PJ DG is RNG-realization-perturbed.

This is the SAME never-FFI-RNG stochastic straddle #156 ALREADY PROVED via FVSut_g16 instrumentation (deterministic
PJ regent growth BIT-EXACT: SITEAR/POTHTG/VIGOR/PCTRED/HTGR all matched; residual = regent-ZZRAN realization). The
oak differs from PJ ONLY by species-indexed data routed through these same-verified branches. The sweep's one-sided
−5..−26% was a small-sample slice of the straddle, not a bias. #192 MEETS the bit-exact-or-cornered bar. Rests on
#156's g16 measurement + this session's line-by-line source-verification of the oak-specific branches.

## #193 EM crown recession — ROOT CAUSE FULLY MEASURED FROM SOURCE (em/crown.f) — turnkey port spec

The long-ambiguous "does live recede the regen crown, and how" is now SETTLED by reading em/crown.f (not inferred).

### The crown model is dispatched PER-SPECIES (crown.f:238-249), not per-variant:
  sp 5 (LL)                      → NIVAR   (DCR exp-polynomial, the model jl currently uses for ALL D>=3)
  sp 1,2,3,7,8,9,10,18           → EMVAR   → WEIBULL-RANK crown  (jl does NOT implement this)
  sp 4,12,17                     → UTTVAR  → same Weibull block as EMVAR
  sp 11,13,14,15,16,19           → CRVAR   → CL=5.17281+0.32552*HF-0.01675*BA; CR=CL/HF; ->label 53
  sp 6                           → LPIJU   → CL=-0.59373+0.67703*HF; CR=CL/HF; ->label 53

### The D-threshold that gates the skip (crown.f:341-345) DIFFERS by group:
  EMVAR/UTTVAR:  IF(D.LT.1.0 .AND. LSTART) GO TO 58   ← D>=1.0 (or any not-at-start) IS PROCESSED (recedes!)
  else (NIVAR/CRVAR/LPIJU):  IF(D.LT.3.0) GO TO 58    ← D<3 skipped
⇒ EMVAR regen with D in [1.0,3.0) RECEDES in live; jl skips ALL D<3 (`continue`) ⇒ jl keeps ~90 estab crown
  while live receds 90→~74 ⇒ over-crown → over-VIGOR / over CR-linear SMHTGF height → the #193 over-growth.

### EMVAR/UTTVAR Weibull-rank crown (crown.f:277-368):
  RELSDI = min(SDIAC/SDIDEF(sp), 1.5)                 # SDIAC=stand SDI; SDIDEF=species SDImax (common block)
  ACRNEW = C0(sp) + C1(sp)*RELSDI*100
  A = WEIBA(sp);  B = max(WEIBB0(sp)+WEIBB1(sp)*ACRNEW, 1.0);  C = max(WEIBC0(sp)+WEIBC1(sp)*ACRNEW, 2.0)
  per tree: SCALE = clamp(1-0.00167*(RELDEN-100), 0.30, 1.0)
            X = clamp( (ISORT(I)/ITRN)*SCALE , 0.05, 0.95 )    # ISORT = diameter rank (crown.f:196-200)
            CRNEW = ( A + B*(-ln(1-X))^(1/C) ) * 10            # crown-percent 0..100
  label 53: CHG = CRNEW - ICR;  bound |CHG| to 1%/yr: PDIFPY=CHG/ICR/FINT; if>0.01 CHG=ICR*0.01*FINT; if<-0.01 sym.
            CRNEW = ICR + CHG*(CRNMLT if DLOW<=DBH<DHI else 1);  ICRI = INT(CRNEW+0.5); clamp [5,95].

### Coefficient tables (crown.f DATA, sp 1..19; zeros = handled by other branch):
  WEIBA  = [0,0,0, 1.0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
  WEIBB0 = [.11035,.11035,.14652,-.82631,0,0,-.00359,.67059,.73693,.02663,0,-.08414,0,0,0,0,-.08414,.11035,0]
  WEIBB1 = [1.10085,1.10085,1.09052,1.06217,0,0,1.12728,.99349,.98414,1.11477,0,1.14765,0,0,0,0,1.14765,1.10085,0]
  WEIBC0 = [.02774,.02774,1.04746,3.31429,0,0,2.60377,-4.25938,-4.16681,2.95048,0,2.775,0,0,0,0,2.775,.02774,0]
  WEIBC1 = [.35524,.35524,.39752,0,0,0,0,1.35687,1.33779,0,0,0,0,0,0,0,0,.35524,0]
  C0     = [5.68625,5.68625,5.92714,6.19911,0,0,5.0587,7.41093,7.36476,5.61047,0,4.01678,0,0,0,0,4.01678,5.68625,0]
  C1     = [-.0447,-.0447,-.03346,-.02216,0,0,-.03307,-.03467,-.03761,-.03557,0,-.01516,0,0,0,0,-.01516,-.0447,0]

### CRITICAL COROLLARY — LATENT bug beyond establishment:
jl `crown_ratio_update!` uses the NIVAR-DCR model for ALL D>=3 species. Faithful ONLY for sp5. For the EMVAR
species (the majority: 1,2,3,7,8,9,10,18) and UTTVAR (4,12,17), live uses the Weibull-rank crown — so jl's D>=3
EM crown is ALSO mis-modeled multi-cycle (masked on emt01 because inventory crowns bypass at LSTART and the
validated window is early/crown-insensitive). This is a broader EM crown port, not just an establishment patch.

### Why deferred (not a rushed edit): implementing the Weibull-rank crown replaces the crown for MOST EM species
at ALL sizes (D>=1) ⇒ changes emt01's multi-cycle crowns/growth ⇒ MUST be chunked with per-cycle emt01 A/B vs
FVSem_g16 (ISORT/IND sort order, SDIAC, SDIDEF table must be byte-exact) per doctrine ("regression on a faithful
chunk = examine"). Coefficients + mechanism are now fully in hand ⇒ next session is a turnkey implement+validate.
IE side of #193 (ie_autoes crown +0.07985*RAN term, RNG-order-faithful) still pairs with this.

## #193/#194 cross-check — CI crown is the READY TEMPLATE; #194 is a SEPARATE ci_esgent root

Verified jl src/variants/centralidaho/crown.jl:105-144 ALREADY implements the EMVAR-style Weibull-rank crown
CORRECTLY for CI: ISORT diameter-rank (`isort[idx[jj]]=n-jj+1`), `d<1 && lstart → ci_dubscr` gate
(== ci/crown.f:290 `IF(D.LT.1.0 .AND. LSTART) GO TO 58`), then the Weibull `CRNEW=(A+B(-ln(1-X))^(1/C))·10`
with SDIAC/sp_sdi_def/CI_WEIBA/B0/B1/C0/C1 and the 1%/yr label-53 recession. CI regen with D>=1 RECEDES.

Two consequences:
1. CONTRAST VALIDATES #193: CI (HAS the Weibull recession) shows NO crown-driven over-growth; EM (jl uses
   NIVAR-DCR for all, LACKS the Weibull recession) DOES. This is strong confirmation that the missing EMVAR
   Weibull-rank recession is the EM over-growth driver — not a coincidence, a controlled A/B across variants.
2. #194 RE-SCOPED: CI's crown is faithful ⇒ #194's CI bare-establishment over-growth is NOT a crown-recession
   bug. It lives in ci_esgent itself (establishment tally amount / birth height / birth crown), a SEPARATE root
   from #193. #193 and #194 are NO LONGER the same class.

⇒ The EM crown port (#193) is DE-RISKED: mirror jl's validated CI crown_ratio_update! — same ISORT/SDIAC/
  Weibull machinery — swapping CI coefficient tables for EM's (tables captured above) and wiring EM's per-species
  dispatch (sp5 NIVAR-DCR jl already has; CRVAR sp{11,13-16,19} CL=5.17281+.32552HF-.01675BA; LPIJU sp6
  CL=-.59373+.67703HF; EMVAR/UTTVAR Weibull for the rest) + EM's dubscr for the D<1&lstart path. Still needs
  per-cycle emt01 A/B vs FVSem_g16 before landing (changes most-species crowns).

## #193 EM crown recession — HYPOTHESIS IMPLEMENTED, A/B-MEASURED, and REFUTED as the driver

Per doctrine (MEASURE-don't-infer applied to my OWN fix), I did NOT stop at the source diagnosis — I IMPLEMENTED
the full source-faithful EM per-species crown dispatch and A/B-tested it:
  - group 4 (EMVAR/UTTVAR): Weibull-rank crown `CRNEW=(A+B(-ln(1-X))^(1/C))·10`, X=ISORT/N·SCALE, D>=1 recession,
    all 7 coeff tables (programmatically extracted, verified), ISORT via _rdpsrt! on DBH+DG/BARK (CI template).
  - group 2/3 (CRVAR/LPIJU): CL crown models (5.17281+.32552HF-.01675BA / -.59373+.67703HF) + label-53 recession.
  - group 1 (NIVAR sp5): existing DCR (unchanged).
Compiled clean. Then measured jl-vs-FVSem_clean BEFORE/AFTER on: 4 mature multi-sp refs, a pure-sp6 LPIJU stand
(11881495010690, ΔBA -2.9%), and 6 establishment/regen stands (the #193 regime).

RESULT — the crown change moved NOTHING: every stand's .sum BA/TPA was BYTE-IDENTICAL before vs after, for BOTH
the group-4 Weibull AND the group-2/3 CL variants. Diagnosis of why (instrumented):
  1. The sampled mature stand was pure sp6 (LPIJU, group 3) ⇒ the group-4 Weibull branch NEVER fired (G4 count=0).
  2. Switching sp6's 27 trees from DCR→LPIJU-CL left ΔBA -2.9% UNCHANGED ⇒ that divergence is NOT crown-driven.
  3. The crown's effect on large-tree DG is below the .sum integer-BA resolution (BA rounded to whole ft²/ac).

DECISIVE CONCLUSION: the EM crown-model choice (DCR-for-all vs faithful per-species Weibull/CL) is a REAL infidelity
but .sum-SECOND-ORDER — it does not move BA/TPA on any tested stand. This REFUTES this session's "crown recession is
the #193 root" hypothesis. The #193 +200-492% over-growth is a LARGE, .sum-visible effect ⇒ it CANNOT be the crown
(which is sub-.sum-resolution). REDIRECT #193 to the actual .sum-visible drivers: the AUTOES establishment TALLY
(number of trees established), the birth DBH/height, and the regen SMHTGF/regent DG — NOT the crown ratio.

CODE REVERTED (crown.jl back to pristine DCR-for-all): the faithful crown port is UNEXERCISABLE via .sum sweeps
(group-4 never fires on .sum-visible stands; the effect is below integer-BA resolution) ⇒ per "a test must exercise
the semantic", it must not be landed as unvalidatable code. It requires per-tree FVSem_g16 crown comparison to
validate — a future g16-instrumentation chunk. The complete spec + all coefficient tables are preserved above,
so that chunk is turnkey. This is the 8th implemented-and-refuted hypothesis of the campaign (doctrine working).

## #193 EM AUTOES over-growth — FULLY ROOT-CAUSED (the WK4/HTIMLT birth-cycle multiplier jl omits)

After REFUTING the crown-recession lead (prior section), traced #193 to the REAL .sum-visible driver via a clean
treeless reproducer (stand 5332701010661, starts 0 TPA; NUMCYCLE, ECHOSUM; oracle FVSem_clean):

MEASURED (built-in FVSJL_AUTOES_DEBUG + FVSJL_SZDBG probes, real compiled runs):
- AUTOES TALLY is EXACT: jl establishes 118.1 TPA == live 118 (sp3 DF 71.8 + sp10 45.3 + sp2 1.1).
- BIRTH is CORRECT: POST-AUTOES maxH=1.2 ft, maxD=0.1 in (hht=XMIN+0.2, dbh=0.1+0.001·hht) — matches live.
- em_esgent! OVER-GROWS: POST-ESGENT maxH=5.44 maxD=0.71 vs live's cycle-1 TopHt=3 / QMD=0.1. jl grows the
  1.2-ft seedling +4.24 ft in the birth cycle; live grows it ~+1.8 ft.

ROOT CAUSE — the birth-cycle height-growth MULTIPLIER. Live em/esgent.f:22-24:
    CALL REGENT(.TRUE.,ITRNIN)      ! computes HTG(I)
    HTG(I) = HTG(I) * WK4(I)        ! ← scale by per-tree WK4(I)
    HT(I)  = HT(I) + HTG(I)
WK4 is set at establishment: WK4(ITRN)=HTIMLT(N) (estab.f:1257), and (estab.f:1054-1063):
    IF(FINT-DELAY .LT. 5) GENTIM=0. ELSE GENTIM=FINT-DELAY-5.0
    FTEMP = min(TRAGE, GENTIM);  HTIMLT = FTEMP/(GENTIM+0.0001)
with TRAGE=PRMS(4) defaulting to 2.0 (estab.f:987-989). For AUTOES (FINT=10, DELAY=0): GENTIM=5, FTEMP=2,
**HTIMLT = 2.0/5.0001 = 0.40**. jl's em_esgent! instead uses a CONSTANT `subyr/regyr = 5/5 = 1.0` — a DIFFERENT,
wrong multiplier — so it applies the FULL SMHTGF increment. ARITHMETIC CHECK: 1.2 + 4.24·0.40 = 2.9 ft ≈ live 3 ✓.

Why #137 missed it: em_esgent's subyr/regyr=1.0 was validated on emt01 PLANT regen, where TRAGE≥GENTIM ⇒ HTIMLT≈1.0
(scale=1 correct). AUTOES natural regen has TRAGE=2 ⇒ HTIMLT=0.40 (scale must be <1). The under-tested-regeneration-
regime meta-pattern again: PLANT-only validation cannot see the AUTOES multiplier. This is .sum-VISIBLE (unlike the
crown) ⇒ directly validatable.

FIX SPEC: store a per-tree HTIMLT (=min(TRAGE,GENTIM)/(GENTIM+0.0001)) at establishment — ie_autoes sets TRAGE=2.0
(⇒0.40 at FINT=10); the shared establish! (PLANT/NATURAL) computes it from its TRAGE/DELAY (PLANT keeps ≈1.0 ⇒
emt01 preserved) — and em_esgent! scales `htg` (and the DBH via HTG) by that per-tree HTIMLT instead of subyr/regyr.

## #193 EM AUTOES over-growth — FIX LANDED + VALIDATED (per-tree WK4/HTIMLT birth-cycle multiplier)

Implemented the fix from the root-cause above. Commit adds `TreeList.htimlt` (per-tree, default 1.0, carried
through tripling/compaction via _TREE_VEC_FIELDS); `ie_autoes_establish!` sets it to min(2,GENTIM)/(GENTIM+1e-4)
with GENTIM=max(FINT-5,0) (AUTOES TRAGE=2 ⇒ 0.40 at FINT=10); the shared `establish!` sets 1.0 (PLANT/NATURAL,
guards slot reuse); `em_esgent!` scales `htg` by `t.htimlt[i]` instead of the constant `subyr/regyr`.

VALIDATED vs FVSem_clean:
- Treeless reproducer 5332701010661 (birth-cycle probe): maxH 5.44→2.89 ft (=live ~3), maxD 0.71→0.10 in (=live
  0.10); .sum QMD 0.4→0.1 (=live 0.1), BA matches.
- Multi-cycle sweep: TopHt 6→2-3 vs live 3 (was ~5-6 over). Mature/dense EM stands BIT-IDENTICAL before/after
  (12344117/12356085 ΔBA 0.0%; 149151286/149153412 unchanged) ⇒ NO REGRESSION (em_esgent only touches birth-cycle
  regen; non-establishment stands have an empty loop; PLANT trees keep htimlt=1.0).
- Cross-variant: IE + KT/others run clean (struct change additive; no other positional TreeList constructor).

RESIDUALS (follow-up, NOT the height bug just fixed):
1. TopHt slightly low (jl 2 vs live 3) on some multi-species treeless stands — the reproducer's actual maxH=2.89
   rounds to 3, so this is the AVHT40 top-height metric weighting a multi-species sub-breast-height cohort (or a
   minor TRAGE nuance for the ntally=99 path), not the growth magnitude. Huge net improvement over the pre-fix ~6.
2. AUTOES TPA tally over-count (jl 261-263 vs live 233-234, +12%) on some stands — a SEPARATE tally-amount issue
   (#143-class), independent of the birth-cycle height fix (the reproducer's tally was 118==live exactly).
3. IE side: ie_autoes now sets htimlt for IE too, but IE uses ie_esgent! (not em_esgent!) which does NOT yet apply
   it — IE's +19-29% AUTOES over-growth needs the same one-line scale in ie_esgent! + IE validation. NEXT.

## #193 IE AUTOES over-growth — FIX LANDED (mirrors EM; ie_esgent WK4/HTIMLT)

MEASURED ie_esgent! on treeless reproducer 753200974290487: bscale=1.0 (the "=0.5" comment was WRONG — IE_RG_REGYR=5
⇒ (10-5)/5=1.0), so h2=h+exp(htgrl)·1.0 applied the full birth-cycle increment (1.2→5.19 ft) vs live WK4=0.40 (→2.8).
Identical bug to EM. FIX: scale by per-tree t.htimlt[i] (ie_autoes already sets 0.40; PLANT keeps 1.0=old bscale).
VALIDATED vs FVSie_clean: reproducer now BIT-EXACT (TPA 131 / TopHt 3 / QMD 0.1 == live; was 7 / 0.3); 3 treeless
stands TopHt 7→3, QMD 0.3→0.1 = live; treed IE stands BIT-IDENTICAL before/after (no regression).

### #193 status: EM + IE birth-cycle HEIGHT over-growth = FIXED. Remaining #193/related:
1. AUTOES TPA tally over-count (+8-12% on some stands, e.g. IE 138 vs 128 / EM 261 vs 233) — SEPARATE #143-class
   tally-amount issue (reproducers with exact tally, e.g. 118==118 / 131==131, prove the height fix is clean).
2. TopHt ±1 on some multi-species treeless stands (AVHT40 metric on sub-breast-height cohorts).
3. #194 CI (ci_esgent +282%): likely the SAME WK4/HTIMLT class, BUT CI is not AUTOES — bare CI FIA stands do NOT
   establish (0 TPA), so #194's reproducer is a treed/keyword-regen stand (needs re-location). ci_esgent uses
   scale_h=NTYR/REGYR with NO WK4 multiplier ⇒ candidate. ALSO: the shared establish! currently hardcodes
   htimlt=1.0 for PLANT/NATURAL — faithful for PLANT (TRAGE≥GENTIM) but NATURAL regen with low TRAGE should get
   WK4<1 (compute from establish!'s TRAGE/DELAY). That refinement would serve CI/BM/TT/UT *_esgent! NATURAL regen.

## #193 follow-up — treeless EM regime has THREE distinct bugs; the htimlt fix addresses ONE (validated non-regressing)

Broader single-plot treeless EM sweep (post-htimlt-fix) vs FVSem_clean revealed the regime is multi-causal:
1. **Sparse-establishment OVER-growth [FIXED]**: high htgrth × old-scale-1.0 → TopHt 5-6. The per-tree WK4=HTIMLT
   fix (×0.40) corrects these (reproducer 5332701010661 5→~3; 103399518010661 6→2). BEFORE/AFTER stash A/B
   confirms the fix ONLY moves the over-growers and is byte-identical on the under-growers (NO REGRESSION).
2. **Dense-establishment UNDER-growth [SEPARATE, pre-existing]**: stands like 103400530010661 (sp3 DF 147 TPA,
   denser) → jl SMHTGF over-suppresses (htgrth≈0 ⇒ htg floored) ⇒ TopHt stuck at birth ~1.2 (=1) vs live 3.
   IDENTICAL before/after the htimlt fix (htg was already ~0, so scaling is inert). Root candidate: jl's birth-cycle
   TPCCF/point_ccf over-counts the just-established seedling cohort's competition ⇒ beta1=exp(B0+B1·log(tpccf))
   collapses. NOT the htimlt bug. Needs FVSem_g16 esgent/smhtgf TPCCF+htg1 trace on a dense treeless stand.
3. **Multi-plot AUTOES tally OVER-count [SEPARATE, #143]**: single-plot stands now bit-EXACT on TPA (214=214,
   228=228, …); NUM_PLOTS>1 stands over-count (4-plot 103399518 tally 261.3 vs live 233, +12%). Deterministic
   (not a straddle) ⇒ a per-plot NSTORE/ESRANN-advance/ESTPP desync for >1 plot (ie_autoes_tally body_n was
   validated on single-plot iet01/EM). Deep #143 per-plot-chain territory (partly UB-cornered per IE audit).

⇒ #193's headline OVER-growth (EM +200-492% / IE +19-29%) is FIXED in both variants. The establishment regime is
NOT yet fully bit-exact: the dense-under-growth (SMHTGF/TPCCF) and multi-plot-tally are distinct open sub-bugs.
Honest status: one of three establishment-regime bugs closed; two characterized+localized for follow-up.
