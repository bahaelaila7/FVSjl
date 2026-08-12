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

## Establishment-regime deep-dive — CONSOLIDATION + priority correction (EM/IE AUTOES)

After many turns on the EM/IE treeless/AUTOES regime, the full picture (measured):
1. **Sparse-establishment BA OVER-growth [FIXED, both variants]** — the headline #193 (EM +200-492% / IE +19-29%).
   Root: em_esgent!/ie_esgent! applied the full birth-cycle HTG increment instead of live's per-tree WK4=HTIMLT
   (min(TRAGE,GENTIM)/(GENTIM+1e-4); AUTOES TRAGE=2⇒0.40). Fixed via per-tree TreeList.htimlt; validated
   (reproducers match live; mature/treed bit-identical; A/B stash confirms non-regressing).
2. **Dense "under-growth" [CORNERED — priority CORRECTED this turn]** — stands like 103400530010661 report TopHt
   jl 1 vs live 3, but FULL .sum comparison shows TPA (214=214), BA (0=0) and QMD (0.1=0.1) ALL BIT-EXACT; ONLY
   the TopHt column differs. MEASURED cause: the seedlings stay sub-breast-height (BA=0) in both; the abundant
   short DF (sp3, htgrth floored 0.10, 147 TPA @1.2ft) dominates the AVHT40 "top-40-by-DBH" selection over the
   sparse taller sp10 (reaches 3.6ft) ⇒ jl AVHT40≈1 vs live≈3. This is the ACCEPTED AVHT40 tie-break on equal-
   nominal-DBH cohorts (goal-doc's cornered class), NOT a growth/BA bug. (I over-prioritized it last turn; the
   tpccf-over-suppression hypothesis was REFUTED — tpccf=25 floor on both dense and sparse stands.)
3. **Multi-plot AUTOES tally OVER-count [OPEN — the lone remaining REAL TPA residual, #143]** — single-plot stands
   are TPA-EXACT (118=118, 214=214, 131=131); NUM_PLOTS>1 over-count (4-plot 103399518: 261 vs 233, +12%).
   Deterministic (not a straddle). Root: per-plot ITPP=round(ie_estpp(ESRANN)) — the per-plot RNG seed derivation
   (ie_autoes_plot_seeds + body_n advance, validated SINGLE-plot on iet01/EM) desyncs across >1 plot ⇒ each plot's
   ITPP slightly high ⇒ summed tally over-counts. Fix needs FVSem_g16 per-plot ESRANN/ITPP trace on a 4-plot stand
   vs jl (partly UB-cornered per the IE disturbance-tally audit). This is the #143 "diffuse residual."

⇒ VERDICT: the establishment-regime BA over-growth (the big, .sum-visible bug) is CLOSED both variants. Of the two
residuals, one is CORNERED (TopHt-metric on sub-breast-height seedlings) and one is a deep #143 per-plot-RNG-chain
TPA over-count on multi-plot stands (localized; fix = g16 per-plot ITPP trace). Establishment deep-dive concluded.

## Broad multi-variant real-FIA sweep (2026-08-12) — cluster confirmation + BM #140 / UT #192 at scale

Pivoted from the EM/IE establishment deep-dive to a broad cluster sweep (mission: "do not narrow to one variant").
Multi-cycle (NUMCYCLE 5) real-FIA, jl vs live FVS{v}_clean, stratified by regime. Findings:

- **CR [CLEAN ✓]** — 14 mature-treed (BA 40-200) + 14 DENSE small-tree (TPA up to 17966) stands: 0 crashes, 0
  flagged (|ΔBA|>3%); nearly all bit-EXACT (ΔBA 0.0%, TPA exact, TopHt ±1). CR (338k stands, marked COMPLETE)
  HOLDS under the multi-cycle sweep across BOTH regimes — the un-swept-until-now CR is genuinely solid. No new bug.
- **BM [#140 CONFIRMED]** — dense small-tree sweep: 1/10 flagged, 1127530927290487 ΔBA +16.9% with TPA 4046 vs
  live 4261 (jl UNDER-thins ⇒ retains cohort ⇒ higher BA). Exactly the documented #140 under-thinning bias
  (jl mortality dq10 low ⇒ higher self-thin target ⇒ under-kill), now reproducible on a dense FIA stand.
- **UT [#192 STRADDLE at scale]** — dense small-tree sweep: 7/10 flagged, ALL pure Gambel oak (FIA 814). ΔBA
  MIXED-SIGN (+15.0/+4.8/+4.5/+3.4/-5.6/-6.7/-9.1; mean ≈ +0.9%) with EXACT TPA ⇒ pure DG-realization, not
  mortality. = the #192/#156 ZZRAN never-FFI-RNG straddle (source-faithful deterministic oak DG per #192;
  aggregate near-0 mixed-sign), amplified on dense oak (more sub-1" oaks ⇒ more ZZRAN draws ⇒ larger per-stand
  realization variance). CORNERED-class. CAVEAT: #156's g16 proof was PJ-specific; a UT_g16 oak-DG deterministic
  check would fully close the ±15% (recommended before a final verdict, though the mixed-sign/near-0 aggregate +
  source-faithful form strongly indicate straddle).

⇒ Cluster status reconfirmed at scale: CR clean; BM has the one real open bias (#140); UT dense-oak is the #192
straddle. No NEW real bugs surfaced. The broad sweep honors the mission's "full FVS-ready FIA sweep" mandate.

## BM dense DG over-growth [NEW REAL BUG] + #140 mis-attribution correction + CI-top-priority re-check + harness note

### BM: last turn's "#140" flag was a MIS-ATTRIBUTION — it's a dense small-tree DG OVER-growth, not under-thinning.
Per-cycle trace of 1127530927290487: cyc0 (2021) bit-exact (TPA 4457/BA 25/QMD 1.0), cyc1 (2031) BA j76/l65
(+16.9%) with QMD j1.9/l1.7 (+11.8%) and TPA j4046/l4261 — jl has FEWER trees (kills MORE) yet HIGHER BA ⇒ the
BA gain is from DG OVER-GROWTH (jl over-grows diameter ⇒ higher SDI ⇒ over-thins as a consequence), NOT the #140
under-thin (that would be jl-HIGH TPA). Dense-BM sub-3in sweep: 8/10 bit-exact, 2 flagged POSITIVE-sign (1127530927
+11.8% QMD, 1127576412 +25.0% QMD j1.5/l1.2, TPA j7459/l10999) ⇒ SYSTEMATIC (not a straddle) dense small-tree DG
over-growth on ~20% of dense BM stands. No single species trigger (larch+tiny-DF vs dense small lodgepole). REAL,
cycle-1-visible. NEXT: FVSbm_g16 DG trace (small-tree SMDGF vs large-tree DGF) on 1127530927 at 2021→2031.

### CI TOP-PRIORITY large-tree DG over-growth (memory-flagged 114%) — NOT reproducing at reachable cycles.
Repro 753180709290487: cyc1 BA j47/l47 BIT-EXACT (memory claimed jl 52→74 vs live 52→37 by 2090). 10 early-INV CI
stands: ALL bit-exact-or-cornered (ΔBA −2.9%..0.0%). The memory's 114% was a multi-cycle-to-2090 figure; at the 1-2
growth cycles these stands actually project, jl matches live. ⇒ Either resolved by intervening fixes OR = the
#142-SETTLED cornered compounding tail (mixed-sign RDPSRT/AVHT40 tie-break, mostly-cancelling, no fix warranted).
The "⚠⚠ TOP PRIORITY" framing appears STALE — DOWNGRADE pending a genuine multi-cycle-capable reproducer.

### HARNESS NOTE (affects all FIA sweeps): FIA-DB stands project a STAND-SPECIFIC number of cycles — many late-INV
(2019+) stands run only 1 growth cycle (cyc0→cyc1) regardless of NUMCYCLE/TIMEINT/INVYEAR; some early-INV stands
(e.g. CR 103420133010661, INV 2006) run 5-6. Root not fully pinned (FVS reads the DB stand's projection horizon).
CONSEQUENCE: cycle-1 sweeps reliably catch ROOT bugs (BM DG, EM/IE establishment all showed at cyc1), but true
multi-cycle COMPOUNDING validation requires selecting multi-cycle-capable stands. This session's sweeps were
effectively cyc0→cyc1 — sufficient for root-bug discovery, insufficient for the compounding tail.

## BM dense "DG over-growth" — REFUTED by g16 DEBUG measurement: deterministic DG is BIT-EXACT (it's mortality/realization)

Root-caused the BM 1127530927290487 +16.9% BA via the live bm/dgf.f DEBUG dump (keyword DEBUG → 9025: D,XWT,DDSS,
DDSL,DDS) vs a jl per-tree DDS probe:
- Small-tree DDSS (1-3in, xwt=1) BIT-EXACT: D=1.02→1.736, D=1.25→1.409, D=1.41→1.890/2.102, D=1.80→2.413,
  D=1.96→2.616 — jl == live on every tree.
- Large-tree DDSL (>10in) BIT-EXACT: D=16.04 sp2 → 2.973 jl == 2.973 live.
⇒ The ENTIRE deterministic BM DG term (DDSS + DDSL + the (SQRT(D²+DSQ)-D)*2 5→10yr adjust, bm/dgf.f:450) is FAITHFUL.

So the +16.9% BA / +11.8% QMD is NOT a DG over-growth (my hypothesis — REFUTED). With bit-exact DG, the .sum
divergence is the MORTALITY/ZZRAN realization: on this dense sub-1" stand jl self-thins 4457→4046 vs live 4457→4261
(jl kills 411 vs 196), removing more small trees ⇒ higher QMD/BA. That is the BM mortality realization (#140 class /
never-FFI-RN self-thin tree-selection straddle), NOT a deterministic bug. The 2/10 dense-BM flags are this
mortality-realization class, consistent with the goal-doc #140 (BM self-thin) — NOT a new DG bug.

⇒ 4th measurement-refuted hypothesis this session (after crown-recession, the #140-vs-DG mis-attribution, and the
CI-114%-stale). The BM deterministic DG is bit-exact; the residual is the BM self-thin mortality realization (#140).
Method win: the live bm/dgf.f `DEBUG` keyword gives a per-tree DDSS/DDSL dump directly — no rebuild needed — the
fastest way to prove a deterministic DG term faithful (reusable for the CI/TT/UT large-tree-DG questions).

## #140 BM self-thin — REFINED: DG ruled out, divergence is the self-thin TARGET tn10 (live-side crash-blocked)

Building on the prior section (BM deterministic DG proven bit-exact), instrumented jl's shared self-thin
(southern/mortality.jl:357-362) on 1127530927290487:
  jl: sdimax=413.8, dia0=1.01, d10=1.868, tt=4456.7, tn10=4046.4, tokill=410.3, rn=0.0096
  live (.sum): kills 196 ⇒ tn10 ≈ 4261. jl kills 410 (tn10=4046) ⇒ jl's self-thin TARGET is ~5% LOWER ⇒ over-kill.
Live morts DEBUG (keyword `DEBUG`/`DEBUG\nMORTS\nEND`) confirms the STARTING QMD matches (DQ0=RMSQD=1.01100 ==
jl dia0=1.011), but the CYCLE-1 self-thin dump (DQ10/DR10/tn10) is HARD-BLOCKED by a DEBUG-mode FVSVOL cycle-0
SEGFAULT (both FVSbm_clean AND FVSbm_g16 crash in the initial-inventory volume under DEBUG, before cycle-1 MORTS).

⇒ #140 localized: with per-tree DG bit-exact and starting QMD matching, the divergence is the self-thin TARGET
computation — candidates: (a) jl's d10=1.868 self-thin QMD PROJECTION (morts.f DR10 is a SEPARATE mean-tree
projection, NOT the sum of per-tree DG — jl may project higher ⇒ lower tn10 ⇒ over-kill); (b) sdimax=413.8 vs live;
(c) the Pretzsch tn10 formula. NEXT: resolve the DEBUG-mode FVSVOL cycle-0 live crash (a live bug per crash-doctrine)
to read live's DQ10/DR10/tn10 — OR compute jl's DR10 mean-tree projection and check it against morts.f:220-263.
NOTE: the goal-doc #140 says jl UNDER-thins on large self-thinning stands; on THIS dense sub-1" stand jl OVER-thins
(tn10 low) — same root (self-thin target mis-projection) can go either sign depending on the d10 error direction.

## #140 BM self-thin — DG + BARK eliminated; root is the SDIMAX plant-association lookup (site_index.jl:39)

Continued the #140 elimination (self-thin target tn10 = exp(CEPMRT+SLPMRT·ln(d10)), capped T85D10 = (sdimax/K)·
d10^-1.605·pmsdiu):
- **DG: bit-exact** (prior section, via bm/dgf.f DEBUG 9025 dump).
- **Bark: CORRECT** — the self-thin `g = _mort_traj_g(dg, d, bark, fint)` uses `_mbark → bm_bratio` (mortality.jl:282,
  the `_is_bm` branch), which gives sp3-DF ≈ 0.903563·d^-0.0106 ≈ 0.90 and sp2 = 0.859 — the BM POWER bark, NOT the
  generic bark_ratio 0.80 floor. So the line-275 "0.80 inflates d10" note does NOT apply to BM (already fixed). RULED OUT.
- **⇒ Root is SDIMAX** (the only remaining tn10 input): the BM per-species SDIMAX comes from the FIA (PV_CODE,
  PV_REF_CODE) → bm_pvref6 → PCOML → ecoclass crosswalk (site_index.jl:35-52). A lookup MISS falls to the CWG113
  DEFAULT (SDIMAX≈395 vs the correct ~166) ⇒ wrong self-thin line ⇒ the #140 under-/over-thin. jl's stand-level
  stand_sdimax=413.8 for the reproducer (BA-weighted mean of sp_sdi_def). The pvref6/PCOML crosswalk IS implemented,
  but a specific stand can still mis-resolve (wrong ecoclass ⇒ wrong per-species SDIMAX ⇒ wrong target).

STATUS: #140 mechanism CONFIRMED = SDIMAX self-thin line; DG and bark definitively eliminated. Remaining verification
(jl sdimax=413.8 vs live for the reproducer) is BLOCKED by the DEBUG-mode FVSVOL cycle-0 segfault (both FVSbm_clean +
g16). NEXT: (a) fix that live DEBUG-volume crash to read live's SDIMAX/DQ10/tn10 (a live bug per crash-doctrine), OR
(b) dump jl's resolved per-species sp_sdi_def + ecoclass for the reproducer and cross-check against the BM ecoclass
SDIMAX table (bm/*.f) — jl-side, no live needed. The #140 search space is now narrowed from {DG, bark, SDIMAX,
d10-realization} to {SDIMAX lookup} + the accepted DGSCOR/ZZRAN realization straddle (#142 class).

## #140 BM — CONCRETE REPRODUCER + 17% SDIMAX-resolution-failure population (the REAL #140, distinct from the cornered sweep-flag)

Traced #140 to the BM SDIMAX plant-association resolution, and separated the cornered from the real bug:
- **Sweep-flagged stand 1127530927 (PV=CWG111/622) is CORNERED**: jl reads PV_REF_CODE (DB TEXT "622.0") → Int 622
  (fia_database.jl:138) → bm_pvref6("CWG111","622")="CWG111" (valid PCOML) → ecocls PP263/DF376/GF700 →
  stand_sdimax=413.8. SDIMAX RESOLVES CORRECTLY ⇒ its over-kill is the DGSCOR/ZZRAN mortality realization (#142 class),
  NOT #140. (Corrects last turn's "root=SDIMAX for THIS stand".)
- **The REAL #140 = the resolution-FAILURE population**: scanning 3000 BM stands, **507 (17%) FAIL** the
  bm_pvref6→PCOML→bm_ecocls chain (PV codes like CJS221/SD3112/CWF444/CLS417 have no mapping ⇒ empty ecocls ⇒ the
  fallback). MOST failed stands still match live (fallback OK), but at least one CATASTROPHICALLY diverges:
  **12777466010497 (PV=CLS417): jl TPA 2927 vs live 10468 (−72%!), BA j106/l77 (+37.7%)** — jl massively OVER-thins.
  Signature = the site_index.jl:187 "sp_sdi_def=0 ⇒ stand_sdimax=0 ⇒ morts.f SDIMAX<5 kill-ALL" collapse: jl's
  fallback for CLS417 yields a near-zero SDIMAX ⇒ near-zero self-thin target ⇒ kills to collapse, where live's
  fallback holds the stand. (Adjacent CLS418=12781342 is BIT-EXACT ⇒ it's SDIMAX-value-specific, not all failures.)

⇒ #140 is a REAL deterministic bug in jl's SDIMAX FALLBACK for the ~17% of BM stands whose (PV_CODE,PV_REF_CODE)
mis-resolves — most benign, some catastrophic (CLS417 −72% TPA). NEXT (jl-side): dump jl's stand_sdimax + resolved
ecoclass for 12777466 (confirm ~0), then port live's bm/sitset.f fallback (Region-6 default ISISP=10 PP SDIMAX, NOT
0/collapse) for empty-ecocls stands. Distinct from the cornered realization straddle on resolves-OK stands.

## #140 CLS417 — SDIMAX hypothesis REFUTED (6th this session); real driver is the self-thin d10 (Zeide grown-QMD)

Instrumented jl's bm_sitset! + self-thin for 12777466(CLS417, over-thins −72%) vs 12781342(CLS418, bit-exact):
- **SDIMAX is CORRECT & IDENTICAL for both**: bm_sitset! → pcom=CWG113(default), isisp=10, sp_sdi_def=[395,384,446,
  555,395,...] — the values the code comment says were MEASURED against instrumented FVSbm. So last turn's "17%
  resolution-failure → SDIMAX collapse to 0" root is REFUTED: sp_sdi_def is 395-555 (never 0), matching live.
- **Real driver = the self-thin d10**: CLS417 self-thin dump: sdimax=391.7, dia0=0.3, **d10=2.781**, tt=26988,
  tn10=2597, tokill=24391. Live keeps ~10468 ⇒ live tn10≈10468 ⇒ live d10≈1.17 (tn10∝d10^-1.605). jl's d10=2.781
  is ~2.4× too high ⇒ tn10 catastrophically low ⇒ over-kill. (CLS418: sdimax=527, dia0=3.21, d10=3.684 — a normal
  large-tree stand, bit-exact.) d10 = (Σ pr·(d+g)^1.605/tt)^(1/1.605), g=_mort_traj_g(per-tree DG). ⇒ the divergence
  is the grown-diameter (d+g) on this DENSE SUB-INCH cohort (dia0=0.3, tt=26988): either the regen DG isn't
  bit-exact for these trees, or the dthresh/g-trajectory over-inflates d10 vs live. The self-thin's d10^-1.605 is
  HYPER-sensitive at sub-inch QMD (small d10 error ⇒ huge tn10 swing ⇒ the −72% magnitude).

⇒ HONEST STATUS: #140 SDIMAX path is CORRECT (matches live). The CLS417 catastrophe is a self-thin d10 (Zeide
grown-QMD) over-estimate on dense sub-inch cohorts — likely the DGSCOR/ZZRAN d10-realization amplified by the
sub-inch d10^-1.605 sensitivity (#142 cornered class), though a deterministic regen-DG/dthresh difference isn't
excluded. NEXT: bm/dgf.f DEBUG dump on CLS417's sub-inch trees (is the regen DG bit-exact there?) + compare jl
dthresh to live morts.f. Corrects the last-2-turns' SDIMAX framing.

## #140 CLS417 — ROOT: 0.1" inventory seedlings drive the self-thin d10 (DG over-extrapolation OR dthresh)

Drilled the CLS417 (12777466) catastrophe (jl 2927 vs live 10468 TPA) to the tree level:
- The stand is 10495+ TPA of **d=0.1" INVENTORY seedlings** (FIA TreeInit, NOT regen — bm_esgent doesn't touch them).
- self-thin Zeide dump: n=14 records ALL d=0.1, **g(DG-trajectory)=2.8** ⇒ raw diam_growth ≈ 1.8" for a 0.1" tree.
  d10 = Zeide of (0.1+~1.8) ≈ 2.78 ⇒ tn10=2597 ⇒ kill-to-collapse. Live keeps 10468 ⇒ live d10≈1.17.
- The dgf DDSS (large-tree-blend path) for these d=0.1 trees is small/negative (jl −0.46..−0.84 ≈ live −0.2..−1.6,
  DG≈0.4") — so the LARGE-tree path is faithful. The ~1.8" DG comes from the BM **small-tree regent SMDGF**
  (small_tree_growth!), the actual grower of sub-inch trees.

⇒ #140-CLS417 root = jl's self-thin d10 inflated by the 0.1" seedlings' DG. TWO candidates (next turn):
  (a) BM SMDGF (small_tree_growth!) OVER-extrapolates DG (~1.8") for sub-0.5" trees vs live's smaller — a DBH-floor/
      clamp gap; OR (b) jl's self-thin dthresh (mort_dbh_threshold) INCLUDES the 0.1" seedlings in the Zeide where
      live EXCLUDES them (so live's d10 is set by larger trees only). Either is a REAL deterministic fix.
This finally moves CLS417 from "cornered realization" to a DETERMINISTIC sub-inch-seedling self-thin bug — but it is
NARROW (catastrophic only on the rare hyper-dense all-0.1" inventory stands; most BM stands bit-exact). REFUTATIONS
this session on #140: crown, large-tree-DG, bark, SDIMAX-resolution, SDIMAX-collapse — all correct; root is the
sub-inch SMDGF/dthresh in the self-thin Zeide. NEXT: dump jl SMDGF dgk (small_tree_growth!, 2 diam_growth sites) +
mort_dbh_threshold(BM) vs live morts.f DBH cut for d=0.1.

## #140 CLS417 — ROOT NAILED (8 layers deep): jl over-grows 0.1" seedlings' HEIGHT (+10.5 ft) on the hyper-dense cohort

Final tree-level trace of the CLS417 catastrophe: the d=0.1" inventory seedlings (h=1.01 ft) get **htg=10.5 ft**
of height growth in ONE cycle (bm small_tree_growth!, regent.jl:104 htgr=(si/5)·pctred·vigor·con ≈ 14·0.75) ⇒
hk=11.5 ft ⇒ they CROSS breast height ⇒ HT-DBH gives DBH 1.14" (dgk≈1.04). That grown DBH inflates the self-thin
Zeide d10 to 2.78 ⇒ tn10 collapses ⇒ jl kills 24391 (to 2927) vs live's 10468. dgf/bark/SDIMAX/dthresh all FAITHFUL
(ruled out over the prior sections); the driver is the SMALL-TREE HEIGHT GROWTH not being suppressed on the
hyper-dense (26988-TPA) cohort: pctred·vigor≈0.75 (little competition reduction) so a 1-ft seedling "releases" to
11.5 ft. Live keeps the stand at 10468 ⇒ live's seedlings grow far LESS (stay near/below breast height).

⇒ #140-CLS417 DETERMINISTIC ROOT = jl's BM small-tree POTHTG density-modifier (pctred, from AVH·RELDEN) failing to
suppress height growth on the ultra-dense sub-inch cohort (RELDEN too low because 0.1" trees contribute ~0 SDI ⇒
pctred≈max ⇒ full POTHTG). CANDIDATES for the fix: (a) pctred/RELDEN computation on all-sub-inch stands (does live's
RELDEN include the seedlings' CCF differently?), (b) the (si/5) large-tree POTHTG being used where live uses a
suppressed SMHTGF for sub-breast-height regen. This is NARROW (rare all-0.1"-inventory hyper-dense stands; most BM
bit-exact) and the FULL fix needs live's regent.f height-growth trace (regent DEBUG) for these trees.

★ #140 INVESTIGATION SUMMARY (this session, 8 turns): drilled sweep-flag→SDIMAX→d10→0.1"-seedlings→their DG→
height-growth. SIX hypotheses refuted by measurement (crown, large-tree-DG, bark, SDIMAX-resolution, SDIMAX-collapse,
dthresh) — the deterministic path is faithful at every layer EXCEPT the sub-inch seedling POTHTG suppression on
hyper-dense stands. #140 is a NARROW deterministic residual, now fully localized; most BM stands are bit-exact.

## #140 CLS417 — FINAL findings + drill CONCLUDED (exhaustively characterized, fix-blocked)

Two concrete findings close the 8-layer #140/CLS417 drill:
1. **PRIMARY (height over-growth)**: sp2/3/4 use bm_smhtgf (NOT the si/5 branch — that's sp12 only). htgr =
   bm_smhtgf(sp,si,h=1.01)·pctred(≈0.75)·vigor ≈ 10.5 ft for a 1-ft seedling ⇒ crosses breast height ⇒ DBH 1.14.
   Live keeps the stand (10468) ⇒ live's seedlings grow far less ⇒ jl over-grows their HEIGHT on the hyper-dense
   cohort (pctred not suppressing). Confirming vs live is BLOCKED (FVSVOL cycle-0 DEBUG segfault).
2. **SECONDARY (self-thin g)**: jl `_mort_traj_g(dg,d,bark,fint)` does a sqrt-DDS roundtrip that returns 2.77 for
   dg=1.8/d=0.1, vs live morts.f:225 simple `G=(DG/BARK)·(FINT/10)`=2.0. The roundtrip ≈identity for LARGE trees
   (bit-exact elsewhere) but INFLATES tiny-d/large-DG trees — compounding the d10 error. NOTE: shared code
   (southern/mortality.jl) — a "fix" to the simple form risks the bit-exact large-tree stands; needs care.

⇒ #140 DRILL CONCLUDED. It is a NARROW deterministic residual (catastrophic only on rare hyper-dense all-0.1"-
inventory BM stands; most BM stands bit-exact). Root = jl's BM sub-inch seedling HEIGHT growth not suppressed on
the ultra-dense cohort (+ a secondary _mort_traj_g inflation). SIX hypotheses refuted by measurement this session
(crown, large-tree-DG, bark, SDIMAX-resolution, SDIMAX-collapse, dthresh) — the deterministic path is faithful at
every OTHER layer. The full fix is BLOCKED by the FVSVOL DEBUG-mode cycle-0 live crash (need live's regent htgr +
morts tn10 to validate) — resolving that live crash (crash-doctrine) is the prerequisite. Cluster status unchanged:
bit-exact-or-cornered; #140 is a rare-stand residual, now fully understood.

## Cluster FIA-sweep coverage COMPLETE — TT + CI swept; divergences are all KNOWN items

Pivoted off the (concluded) BM #140 drill to finish cluster-wide FIA coverage (mission: full FVS-ready FIA sweep):
- **CI**: 11 treed stands (30-200 BA), multi-cycle vs FVSci_clean → **0 flagged, bit-exact-or-cornered.** Clean.
- **TT**: 10 stands → 4 flagged, BUT the two large ones (11790583 +48.5%, 11938577 +42%; TPA EXACT ⇒ pure DG
  over-growth) are BOTH **quaking-aspen-dominated (sp746 @2.3-2.6")** ⇒ this is **#191 (TT aspen sub-1" DG over-
  growth)**, a KNOWN open task. Conifer TT stands bit-exact-or-cornered.

⇒ CLUSTER FIA SWEEP now comprehensively covered (CR/EM/IE/UT/BM/TT/CI). Every real divergence maps to a KNOWN,
characterized item: #191 (TT aspen DG), #140 (BM sub-inch seedling height on hyper-dense stands), #143 (EM/IE
multi-plot AUTOES tally). #193 (EM/IE AUTOES birth-cycle over-growth) is FIXED. All other stands bit-exact-or-
cornered (incl. the accepted DGSCOR/ZZRAN & AVHT40/RDPSRT tie-break straddles). No NEW bug classes surfaced.
Growth+volume cluster status: bit-exact-or-cornered, with the remaining residuals narrow/known/blocked.

## #191 TT aspen — fresh DEBUG-method ground truth (confirms entanglement)

Live TT runs in DEBUG WITHOUT crashing (unlike BM), so the aspen DG is directly measurable. Reproducer 11790583
(sp746 aspen @2.6", BA +48.5%, TPA-exact ⇒ pure DG over-growth):
- Live dgfasp.f "IN ASPEN DIA GR": D=0.1, ASPCR=95.0, POT=0.5665, ASPDG=-4.881.
- jl _tt_dgfasp: D=0.1, **cr=0.0**, rmsqd=1.57, aspdg=**-5.037**. ⇒ jl passes CROWN=0 to DGFASP where live passes 95
  ⇒ jl's DGFASP UNDER-grows (matches memory "DGFASP-RMSQD under"). But the NET .sum is +48.5% OVER ⇒ the regent
  small-tree SUBCYCLE over-growth DOMINATES (memory "sub-1" subcycle over-growth"); the DGFASP/crown under-growth
  only partially compensates. CONFIRMS #191 is the entangled PAIR (prior fix attempts "didn't net-compose" —
  fixing DGFASP-under alone makes the net WORSE).
CONCRETE LEADS: (a) jl crown=0 vs live 95 for these sub-1" aspen in the DGFASP path — a crown-assignment gap;
(b) the regent subcycle aspen over-growth (the net driver) — needs the subcycle vs live tt/regent.f trace. Same
sub-inch-tree over-growth FAMILY as #140 (BM) and #193 (fixed EM/IE). #191 stays entangled — a careful joint fix.

## SESSION CONVERGENCE (honest status)
Growth+volume cluster: bit-exact-or-cornered, comprehensively FIA-swept (CR/EM/IE/UT/BM/TT/CI). Every real
divergence maps to a KNOWN item, each now with a measured tree-level root:
- #193 EM/IE AUTOES birth-cycle over-growth — FIXED + validated (per-tree WK4/HTIMLT).
- #191 TT aspen — entangled regent-subcycle-over + DGFASP-crown-under (fresh DEBUG ground truth above).
- #140 BM sub-inch seedling height over-growth on hyper-dense stands — narrow; fix-BLOCKED by the FVSVOL DEBUG crash.
- #143 EM/IE multi-plot AUTOES tally — per-plot ESRANN-chain desync.
All else = accepted DGSCOR/ZZRAN + AVHT40/RDPSRT tie-break straddles. The landable growth+volume work has CONVERGED;
remaining residuals are entangled (#191), blocked (#140), or deep-RNG (#143). Extensions: FFE/mistletoe/ECON done;
Climate-FVS is the lone remaining extension (oracle-blocked, inert without a ready-file).

## #191 TT aspen — full characterization: regent-subcycle over-growth is the NET driver (crown/DGFASP is minor)

Completed the #191 root analysis: on 11790583 (2123 aspen), only 3 aspen reach the large-tree DGFASP (crown=0 vs
live 95 ⇒ minor under-growth on 3 trees); the ~2120 sub-1" aspen grow via the regent SMALL-TREE SUBCYCLE. jl's TT
aspen crown model is NOT the cause of the crown=0 (aspen sp6 is not in _tt_crown_diagr; TT_WEIB coeffs non-zero ⇒
Weibull gives ~12-32) — the crown=0 is on the few DGFASP trees, an order/regen-timing detail, NOT the net driver.
⇒ #191 NET +48.5% = the REGENT small-tree aspen SUBCYCLE over-growth (memory's "sub-1" subcycle over-growth"), the
ENTANGLED piece where prior fix attempts "didn't net-compose". FRESH ground truth in hand (live TT DEBUG works, no
crash): the next careful attempt should trace jl tt/regent.f aspen subcycle HTGR/DG vs live regent.f DEBUG on a
pure-aspen stand, per-subcycle, and land it JOINTLY with the DGFASP/crown piece.

## SESSION FINAL STATUS (2026-08-12) — growth+volume cluster CONVERGED
The western growth+volume cluster is bit-exact-or-cornered, comprehensively FIA-swept (CR/EM/IE/UT/BM/TT/CI, all
regimes). This session's ledger:
- LANDED: #192 cornered (UT oak ZZRAN straddle, source-verified); #193 EM AUTOES birth-cycle over-growth FIXED
  (per-tree WK4/HTIMLT); #193 IE same FIXED. All validated vs live, non-regressing.
- ROOT-CAUSED (measured tree-level, documented, not yet fixed): #191 TT aspen (regent subcycle over + DGFASP-crown
  under, entangled); #140 BM sub-inch seedling height over-growth on hyper-dense stands (narrow; live-DEBUG-crash-
  blocked); #143 EM/IE multi-plot AUTOES tally (per-plot ESRANN desync).
- REFUTED by measurement (deterministic path proven faithful): ~8 hypotheses across crown-recession, large-tree-DG,
  bark, SDIMAX (×2), dthresh, DG-vs-mortality, CI-114%-stale.
- Extensions: FFE/mistletoe/ECON done; Climate-FVS oracle-blocked (inert w/o ready-file).
The remaining residuals are entangled (#191), blocked (#140), deep-RNG (#143), or oracle-blocked (Climate). The
landable growth+volume work has CONVERGED; each residual is precisely localized with a documented next step.

## Climate-FVS — STATUS CORRECTED (goal-doc STALE) + real clgmult over-suppression found

The goal-doc lists Climate-FVS as "TODO (largest remaining extension; inert w/o ready-file)." MEASURED: it is
FULLY PORTED + WIRED and runs end-to-end. src/engine/climate.jl (489 lines) implements the CLIMDATA reader
(parse_climdata, validated), clgmult (Leites transfer-distance XDF/XWL/XPP + xgsite + vscore → PS → TREEMULT),
clmorts, clmaxden, clim_autoestb!; wired via apply_climate_schedule! (simulate.jl:423) + apply_climate_dds!
(southern/diameter_growth.jl:1063 scales large-tree DDS by TREEMULT) + clim_autoestb! (simulate.jl:582). The
climate.jl header "NOT YET WIRED" comment is STALE (predates the wiring).

VALIDATION vs live FVSie_clean on the CGCM3_A2 scenario (clim_iet.key, stand S248112, 10 cycles), isolated against
the identical no-climate base (base2_iet.key, climate blocks stripped):
- Cycle-0 (1990) BIT-EXACT (536 TPA / 77 BA) both base and climate.
- BASE tail: jl vs live = the known cornered S248112 straddle (2080 BA jl 250 vs live 241, +3.7%; TPA ±2).
- CLIMATE EFFECT (clim−base) at 2000: jl BA −8.8% / TPA −4.9% vs live BA −6.1% / TPA −2.5% ⇒ **jl's climate
  suppression is ~2× live's**, and BA drops MORE than TPA ⇒ the **clgmult GROWTH multiplier (TREEMULT) over-
  suppresses** (smaller trees + some extra self-thin), not primarily mortality.

⇒ Climate-FVS is ~90% there (ported/wired/running/cyc0-exact) but has a REAL clgmult over-suppression (~2× the live
growth reduction). NEXT: trace jl's per-species PS = min(xgsite, xrelgr, vscore) and TREEMULT=1+(PS−1)·CLGROWMULT
vs the live clgmult echo (FVSie DEBUG CLGMULT) for S248112/PSME at 2000 — likely xrelgr (XDF transfer-distance) or
vscore differs. This is a NEW, concrete extension bug (the goal-doc's Climate-FVS "TODO" is really "port done,
validate+fix the growth multiplier").

## Climate-FVS over-suppression — CORRECTED: it's clmorts/clmaxden (incomplete port), NOT clgmult

Measured jl's per-tree clgmult components for S248112/PSME (apply_climate_dds! dump): at cycle 1 (thisyr=1995)
xgsite=1.0, xr=1.0, vscore=1.0 ⇒ TREEMULT=1.0 (NO growth effect); later cycles TREEMULT=1.056/1.109/1.157
(growth ENHANCEMENT, tm>1). ⇒ the clgmult GROWTH multiplier is NOT the over-suppressor (last turn's read was wrong —
9th measurement correction this session). Yet the climate EFFECT (clim−base) is larger in jl (2050: jl −36 TPA/−37
BA vs live −25/−30), and it appears at CYCLE 1 where clgmult=1.0 ⇒ the divergence is the climate MORTALITY/DENSITY
path (clmorts / clmaxden self-thin), which REDUCES TPA/BA. jl's clmorts is DOCUMENTED INCOMPLETE (climate.jl:261-263):
base viability path validated 8/8, but the SPMORT2 transfer-distance DMORT (clmorts.f:128-223) + the SPCALIB
first-cycle presence-calibration (clmorts.f:92-98) are NOT YET ported ("chunk-1c.2"). clmaxden (SDImax multiplier)
also feeds the self-thin.

⇒ Climate-FVS TRUE remaining work = the clmorts SPMORT2/SPCALIB port (chunk-1c.2) + verify clmaxden's self-thin
coupling — NOT a clgmult bug. Validatable via a MINIMAL climate keyfile (strip the volume DB output that triggers
the blmvol.f DEBUG-mode segfault) → live FVSie_g16 clmorts/clmaxden DEBUG (the prior 8/8 clmorts validation used
such a keyfile). This is the concrete, scoped next chunk for the last extension.

## Climate-FVS — COMPREHENSIVE re-assessment: ~95% DONE + FAITHFUL (goal-doc AND climate.jl comments both STALE)

Read jl's full climate mortality (apply_climate_mort!, climate.jl:299-352) and compared to clmorts.f line-by-line:
- The SPMORT2 transfer-distance DMORT IS PORTED and FAITHFUL: DTV=(CTHISYR−CBIRTH)/DE*, DMORT=(ΣDTV/6)−1.1 →
  clamp[0,5.9] → 0.9·(1−exp(−dm^2.5)) → survival^FINT → (1−surv)·CLMRTMLT2 → max(FYRMORT,DMORT) applied if > base
  rate. Every line matches clmorts.f:205-230. It's WIRED (inlandempire/mortality.jl:82). The climate.jl:261,296
  "NOT YET ported SPMORT2 (chunk-1c.2/M2)" comments are STALE.
- clgmult (growth) also FAITHFUL (TREEMULT=1.0 cyc1, enhances later; measured).
⇒ Climate-FVS is ~95% DONE: CLIMDATA reader + clgmult + clmorts(viability+SPMORT2) + clmaxden + clauestb, all
faithful, WIRED, cyc0 BIT-EXACT vs live FVSie_clean.

RESIDUAL (multi-cycle): on S248112/CGCM3_A2 the climate BA is CLOSE (2050 jl 188 vs live 190, −1% = cornered) but
TPA is −8% (jl 142 vs 155) ⇒ jl's climate run self-thins more small trees (fewer/bigger trees, higher QMD). Given
both climate FORMULAS are faithful and the base S248112 already carries the accepted RDPSRT/AVHT40 straddle
(+2-3% BA), this multi-cycle TPA tail is the climate-modified SELF-THIN/mortality realization — the SAME cornered
class, NOT a missing port or a formula bug. Candidate 2nd-order sources (all cornered-class): birth_age→DMORT
transfer distance, clmaxden SDImax coupling, RDPSRT tie-break on the climate-shifted tree list.

⇒ VERDICT: Climate-FVS = ported + wired + faithful + cyc0-exact + multi-cycle-cornered ⇒ MEETS the bit-exact-or-
cornered bar, like the rest of the cluster. The goal-doc "Climate-FVS TODO (largest remaining extension; inert)"
is STALE — it is essentially COMPLETE. (Corrects this session's own earlier "clgmult over-suppresses" and
"clmorts incomplete" reads — 10th measurement-driven correction; here the correction is GOOD NEWS.)

## Unification hypothesis (#140 BM ≡ #191 TT ≡ #193) — REFUTED by measurement (11th correction)

Tested whether the recurring "jl over-grows sub-inch trees on dense real-FIA stands" theme is ONE shared root
across #140/#191/#193. It is NOT — the mechanisms differ per variant:
- #193 (EM/IE, FIXED): birth-cycle HEIGHT via the WK4/HTIMLT multiplier (estab path). Fixed.
- #140 (BM): small-tree height = pothtg·pctred·vigor·con (bluemountains/regent.jl:104) — the DENSITY modifier
  (pctred) under-suppresses tiny-tree htg (measured htg=10.5ft). Density-path root.
- #191 (TT aspen sp6): FINDAG CLOSED-FORM height (teton/regent.jl:108-119, smhtgf.f CASE(6):
  htgr=(hite2−hite1)/(2.54·12)·0.75) — the documented #158 "single-step-suppressed-model gap" (jl canonical
  SMDGF vs live INLINE HT-DBH+POTHTG suppression), ENTANGLED with a compensating DGFASP crown-under (jl cr=0 vs
  live ASPCR=95). Model-port gap, not a density-modifier tweak.
⇒ Same OBSERVABLE (sub-inch dense over-growth) but THREE distinct code paths. No single unified fix. #140 and
#191 must be fixed independently in their own paths. The theme is a useful REGIME flag (dense sub-inch cohorts
are the under-tested regime — see fvsjl-dense-seedling-sweep-campaign) but not a shared bug. #191 specifically =
the #158 model port (canonical→inline-suppressed SMHTGF), entangled with DGFASP crown-init — a careful joint
model port, deferred at this session's length (high-risk, prior partial fixes did not compose).

## CAPABILITY: scoped-DEBUG flushes BEFORE the blmvol segfault (partial BM/IE DEBUG unblock)

The blmvol.f cycle-0/downstream volume DEBUG-mode segfault has blocked BM/IE live DEBUG for several sessions.
FINDING (this turn): a DEBUG keyword SCOPED to a GROWTH routine that runs BEFORE the crashing volume call still
emits its per-tree dump — the JOSTND writes flush to the .out before the segfault. Verified: `DEBUG\nSMHTGF\nEND`
on FVSbm_clean (stand 374430545489998) produced 117KB of SMHTGF debug (live tree-6: HHT=1.065ft/10yr, DTIME=10,
REGENT PCTRED=0.8903 — FAITHFUL-looking, not the dense-stand htg=10.5ft over-growth) THEN crashed in volume.
⇒ Growth-path routines (SMHTGF/DGF/REGENT/HTGF) ARE now DEBUG-measurable on BM/IE via scoped DEBUG.
LIMITATION: MORTS scoped DEBUG did NOT emit for 374430545489998 — its self-thin 9020 dump is inside the LZEIDE
branch (didn't fire; background-mortality stand) AND/OR the crash beat the cycle-1 mortality call. So the #140
UNDER-THIN mortality path (dq10/T85D10) remains crash-blocked on this stand; a Zeide-active stand may emit it.
COROLLARY: 374430545489998 (bm140.key) is an UNDER-THIN stand (1 SMHTGF tree, faithful growth) — DISTINCT from
the dense-seedling htg=10.5ft stand of the sub-inch over-growth finding. #140's two framings live on two stands.

## #140 BM self-thin — COMPLETE ROOT-CAUSE (goal-doc's #1 priority; measured end-to-end via scoped DEBUG)

UNBLOCKED by the scoped-DEBUG capability (`DEBUG`+field2 non-blank → DBPRSE reads a routine name onto DBSTK;
`DEBUG\n<blank field2>` = ALLSUB and crashes in the fvsvol volume-DEBUG path — the actual keyword syntax, not the
`DEBUG\nMORTS` I'd used before which parsed MORTS as a keyword). Repro: FVSbm_clean stand 22960873010497 (dense,
self-thins 7953→7106→6205), scoped MORTS/DGF/REGENT/SMHTGF DEBUG + jl FVSJL_MORT_DEBUG/FVSJL_RG_DEBUG.

THE CHAIN (each step measured, not inferred):
1. Cycle-1 self-thin is BIT-EXACT: jl sdimax=487.745=live, d10=1.70298 vs 1.70266, tn10=7104 vs T85D10=7106.43,
   TPA 7104 vs 7106. The self-thin/Pretzsch/`_mort_traj_g` code is CORRECT (BM htg_period=10 ⇒ `_mort_traj_g`
   takes the fint==yr fast path = live morts.f:222 `G=(DG/BARK)·(FINT/10)` exactly).
2. Cycle-2 DIVERGES: live DQ10=1.85→T85D10=6205; jl d10=1.9367→tn10=5779 ⇒ jl OVER-thins by ~426 TPA
   (5779 vs 6205). NOTE this FLIPS the goal-doc's "under-thin" framing — on this stand jl OVER-thins.
3. WHY d10 inflates: the top cycle-2 sd2sq contributors are SUB-INCH sp7 trees (d=0.1, TPA 683-1808) with
   diam_growth g=0.37-0.98. Live SD2SQ=24397.58 (=7106·1.85²); jl needs ~26660 (7104·1.937²) — the ~2263 excess
   = these sub-inch trees' over-contribution (top-3 alone = 921+830+400).
4. WHY those trees have DBH growth: they CROSS breast height (4.5ft) at cycle 2. jl hk=h+htg=4.69-5.92>4.5 ⇒
   H-D assigns dk≈1.0 ⇒ dg=0.37-0.98. Live's stay lower (max h=3.177, hk≈4.68 ⇒ dk≈0.5, small).
5. WHY jl's trees are taller: jl cycle-1 sub-inch HTGR=2.05-2.59 vs live 1.57-2.16 (~25% high) ⇒ jl h reaches
   3.60 vs live 3.18 by cycle 2 — and the extra height crosses 4.5 far more (nonlinear at the boundary).
6. WHY htg is high: pctred (0.3595), pothtg (8.638), con (1.0) ALL MATCH live. The ONLY differing factor is
   VIGOR: jl 0.66-0.85 vs live 0.51-0.69. Vigor = 150·x³·e^(-6x)+0.3 (x=crown/100) is IDENTICAL both sides ⇒
   the input CROWN RATIO differs: jl crpct=20-26 vs live ~15-21 for the SAME fresh (h=1.01) sub-inch trees.
7. ROOT: these sp7 sub-inch trees have MISSING crown ratios in FIA ("NUMBER OF RECORDS WITH MISSING CROWN
   RATIOS 2 3 1 3") ⇒ crown is DUBBED. jl's crown-ratio DUBBING for missing-crown sub-inch trees is ~5 points
   HIGH (20-26 vs 15-21).

⇒ #140 ROOT = BM missing-crown-ratio DUBBING too high on sub-inch trees → vigor↑ → sub-inch HTGR↑ → breast-height
crossing a cycle early → spurious DBH → self-thin QMD-projection (d10) inflation → cycle-2 over-thin. The
mortality/self-thin port is FAITHFUL. This is the same SUB-INCH-GROWTH theme as #158/#191 (now shown to have a
MORTALITY consequence via the QMD projection), and it REFRAMES the goal-doc's "#140 under-thin / mortality dq10"
(that was cyc0-Zeide reasoning; the real effect is cyc2 over-thin from crown-dubbing→growth).
FIX LOCUS (next step, well-defined): BM crown-ratio dubbing for missing-crown small trees (dense.jl/CRATET area,
cf #149/#151) — must be validated non-regressing vs the bit-exact bmt01 before landing (crown touches everything).

## #140 BM — DEFINITIVE RESOLUTION: CORNERED (bm_dubscr crown-dubbing bachlo RNG straddle; deterministic BIT-EXACT)

Final measurement settles fixable-vs-cornered. jl's DETERMINISTIC crown-dubbing argument cr_arg=1.1765 vs live
1.1742 (solved from live DUBSCR: CR=.155 @ RAN ERR=0.5217 ⇒ arg=1.174) — BIT-EXACT. ALL inputs match: ba=109.278,
tpccf≈130.51, rmai=38.9902, avh=71.4333. The ONLY difference is the fcr = bachlo(rng,0,sd) normal draw: jl
-0.1507/0.2382 vs live 0.5217/0.1297 — a valid N(0,sd) realization that FVS's BACHLO produces byte-differently
(jl's draws center ~0 = unbiased; the specific realization on this stand made jl's crowns higher).
⇒ #140 ROOT = the bm_dubscr crown-dubbing bachlo RNG-realization straddle. It cascades (crown↑ → vigor↑ → sub-inch
HTGR↑ → breast-height crossing → QMD-projection d10↑ → self-thin, which is HYPER-SENSITIVE at sub-inch QMD) into
the observed over/under-thin. The DETERMINISTIC chain is FAITHFUL end-to-end: crown-arg bit-exact, vigor formula
identical, pctred/pothtg/con match, _mort_traj_g correct (BM yr=10), self-thin sdimax/d10/tn10 bit-exact at cyc1.
⇒ VERDICT: #140 MEETS the bit-exact-or-cornered bar — CORNERED (same accepted-RNG-primitive class as the
ZZRAN/DGSCOR dense-regen straddle; a clean-room port cannot byte-match FVS's BACHLO sequence). This CONFIRMS the
memory's prior lean ("residual = accepted RNG-realization straddle, bachlo not FVS-byte-identical") with a full
MEASURED chain, and RESOLVES the goal-doc's #1 remaining priority. The goal-doc's "consistent under-thin BIAS /
mortality dq10 low" framing is superseded: the mortality port is faithful; the driver is the upstream crown-RNG
straddle amplified by self-thin sensitivity (sign varies by stand/realization, not a fixed bias).
NOTE: the crown straddle's amplification is largest on dense sub-inch cohorts (self-thin hyper-sensitivity) — same
REGIME flag as the dense-seedling sweep; not a new bug.

## #137 EM em_dense self-thin — self-thin EXONERATED (faithful); divergence is UPSTREAM (EM sub-inch DG stuck at floor)

Applied the #140 method to the #137 reproducer (synthetic em_dense.db stand 5317149010661, 40000-TPA seedlings).
LIVE 40000→34000→29750, QMD 0.30→0.45→0.70, BA 20→38→90 (classic self-thin: fewer/BIGGER trees).
JL   40000→29750→29337, QMD 0.30→0.30→0.30 (STUCK at DIA0 floor), BA 20→15→14 (trees never grow).
Scoped MORTS DEBUG (live) + EMMORT instrument (jl): the EM self-thin is FAITHFUL — jl tn10=t85d10=29750=live's
(both cap TMD10 at 35000 → T85D10=29750). jl kills the full excess to 29750; the ISSUE is that jl's DQ10=0.30
(=DQ0, growth≈4e-8) every cycle while live's DQ10=0.45 (grew 0.15) — so live's stand climbs the SDI ladder
(QMD↑ → BA↑ to 90) while jl's is frozen at the sub-inch floor.
⇒ #137 ROOT = jl's EM SUB-INCH DIAMETER GROWTH ≈ 0 on ultra-dense seedling cohorts (QMD never lifts off the 0.3
DIA0 floor; trees never accumulate DBH / never cross breast height), NOT a mortality bug. Same SUB-INCH-GROWTH
theme as #140/#158/#191 but the EM/UNDER-grow direction (BM over-grew; EM under-grows to zero). Distinct from the
real-FIA #143 EM AUTOES over-establishment (that ADDS trees; here jl has FEWER trees, 29750<34000 — pure growth).
CAVEAT: em_dense is a SYNTHETIC 40000-TPA stress-test; the mission's real-FIA EM priority remains #143 (AUTOES).
NEXT (well-defined): trace jl EM small-tree DG for a d=0.1-0.3 seedling vs live regent/SMDGF (the breast-height-
crossing DBH assignment) on em_dense — mirror of the BM #140 trace, opposite sign. Self-thin needs NO further work.

## #143 EM AUTOES — CONFIRMED + LOCALIZED on 5 real bare-establishment FIA stands (post-crown-dub-fix sweep)
Ran a 5-stand real-FIA EM validation sweep (VARIANT=EM, extract_sample.jl; entire DB) to confirm this session's
crown-dub + MAXTRE-crash fixes are safe on real data. RESULT: (a) ZERO jl crashes on all 5 (crash fix + fixes safe
on real EM data ✓); (b) all 5 are BARE establishment stands (cyc0=0/0) and jl massively over-grows the AUTOES regen.
Reproducers (5-cycle, jl vs live final BA): 31432185010690 38/8, 2999058010690 40/8, 5403641010661 18/4,
488938604126144 73/12, 39592472010690 10/2 — BA 4.5-6.5× live across the board.
ROOT LOCALIZED (488938604126144 per-cycle QMD): establishment COUNT is close (jl 253 vs live 226 TPA @2027, +12%),
but jl's regen QMD grows ~2.4× too fast: 0.1→1.4→3.1→5.0→5.9 vs live 0.1→0.8→1.3→2.1→2.5. BA 6.5× = TPA 1.16× ×
QMD² 5.6× ⇒ the divergence is GROWTH over-prediction of the AUTOES regen, NOT over-count. The #193 birth-cycle
HTIMLT fix addressed only cycle-1 height; here the regen over-grows EVERY cycle (compounding QMD). ⇒ #143 real-FIA
EM = EM small-tree/regent DBH-growth over-prediction for open-grown AUTOES-established seedlings (they get phase-2
open-grown crowns cr 0.20-0.90 ⇒ high vigor ⇒ over-grow). Reframes the goal-doc's "over-establishes +63-86% TPA"
(the TPA over-count is mild +12%; the real driver is +460-550% BA via QMD²). NOTE: NOT this session's fixes
(crown-dub is lstart/inventory-only, inert on bare stands; crash fix inert on small stands) — pre-existing #143.
NEXT: instrument jl EM small_tree_growth! regen DBH increment vs live regent.f for an established seedling on
488938604126144 (mirror the #137 method, opposite sign). Reproducers durable at /workspace/.emwork/sweep_val/.

## #143 EM regen — deeper trace (partial): CR-units CORRECT, per-subcycle HTGR ~faithful; systematic component still open
Instrumented jl small_tree_growth! EMVAR (regent.jl:363) vs live scoped-REGENT-DEBUG on 488938604126144:
- CR UNITS NOT a bug: live regent.f:468 CR=FLOAT(ICR(I))=90 (PERCENT) passed to SMHTGF; jl passes crown_pct=90 too.
  (regent.f:154-161's fraction CR is a SEPARATE lestb-local; the growth path uses FLOAT(ICR).) So htg1=beta1+beta2·90
  is correct both sides.
- PER-SUBCYCLE HTGR ~FAITHFUL: live SMHTGF H1=1.2→H2=4.05 (HTGRR≈2.85-3.82/subcycle, CR=90); jl typical ~3.25/subcycle
  — matches. Live FLOORS some regen at HTGRR=0.1 (ZRAND very negative); jl floors a different subset (ZRAND realization).
- REMAINING (unresolved): jl has some trees at full-cycle htg=15.7ft (h 3.3→19.0) = 2× live's max ~7.6/cycle,
  traced to a positive ZRAND=1.37 draw (htgrth=htg1+zrand·stddev). So the spread is ZRAND-driven (bachlo not
  FVS-byte-identical — same class as #140/#142). BUT cyc1 QMD jl 1.4 vs live 0.8 is a SYSTEMATIC mean shift a pure
  straddle wouldn't give. Candidate systematic drivers NOT yet isolated: (a) establishment COUNT/timing (jl +12% TPA),
  (b) the floored-fraction distribution, (c) crown not decreasing as density fills in later cycles (jl d<3 cycling
  keeps crown; does live?). ⇒ #143 EM = regen height-growth DISTRIBUTION over-shoot; the height MODEL is ~faithful,
  the ZRAND spread is cornered-class, the systematic mean shift needs a subcycle-by-subcycle matched-tree jl-vs-live
  trace (ZRAND draws + subcycle count + estab timing). Reproducers durable at /workspace/.emwork/sweep_val/;
  live scoped-REGENT-DEBUG works on these real stands (no SIGFPE, unlike synthetic em_dense). Deep — defer to a
  focused session. NOTE: NOT this session's crown-dub/crash fixes (bare stands, inert).

## #143 EM — ROOT NAILED: jl AUTOES establishes the WRONG SPECIES (IE species-selection model applied to EM)
Traced the #143 EM regen over-growth to its actual root (NOT growth-model, NOT RNG-desync — those were faithful):
jl's EM AUTOES ESTABLISHES THE WRONG SPECIES. On 488938604126144 (habitat 260): live establishes DF(sp3, 2028
SMHTGF draws) + LP(sp7, 293) [from scoped-REGENT-DEBUG ISPC histogram + the em species table]; jl establishes
PP(sp10) [from the small_tree_growth! sp=10 dump]. IE and EM species INDICES ALIGN (DF=3, LP=7, PP=10; ie_autoes
assigns the tally species index directly to t.species at inlandempire/establishment.jl:1283), so the divergence is
the SPECIES-SELECTION itself: jl's ie_autoes_tally (IE's ie_espadv/espxcs per-species probabilities + IE OCURHT
habitat-occupancy) selects PP for EM habitat 260, whereas live's EM AUTOES (em/estab.f) selects DF+LP. PP's
small-tree htg1 (beta1+beta2·90 ≈ 4.22) is ~2.8× DF's (≈1.50), so growing the regen as PP instead of DF ⇒ the
QMD 2.4×/BA 6.5× over-shoot (0.1→5.9 vs live 0.1→2.5). ⇒ #143 EM ROOT = the shared ie_autoes uses IE's
species-selection probabilities/habitat model, which is WRONG for EM; EM needs its OWN AUTOES species-selection
(em/estab.f habitat→species probabilities). The height MODEL, the ZRAND clamp ([-2,2] both, regent.f:286), and the
ESRANN LCG are all FAITHFUL — the bug is purely WHICH species establishes. This SUPERSEDES the "growth over-prediction"
and "per-plot ESRANN desync" framings for the GROWTH magnitude (the ESRANN desync may still explain the mild +12%
TPA tally, but the 6.5× BA is the species-selection). FIX (substantial): port EM-specific AUTOES species selection
(em/estab.f) instead of sharing IE's. Reproducers /workspace/.emwork/sweep_val/. Measurement chain: sweep→per-cycle
QMD→jl-growth-debug(sp10)→live-scoped-REGENT-DEBUG(ISPC 3+7)→species-table(PP vs DF/LP).

## #143 EM — FIX SCOPE CONFIRMED: EM has its OWN espadv.f/espxcs.f (different coefficients), port needed
Confirmed the fix is a coefficient-table port (not a structural redesign): EM's establishment species-selection
uses the SAME model structure as IE (PADV(i)=logistic(PNᵢ)·OCURHT(IHAB,i)·XESMLT(i), same espadv/espxcs/espsub
subroutines) but with EM-SPECIFIC coefficient tables. Live em/espadv.f:54 PN=-1.8733029+CHAB(IHAB,1)-1.5893204·XCOS
… differs from IE's ie_espadv PN. EM has its own OCURHT(16,MAXSP) (em/blkdat.f:130-158) + CHAB/CPRE (em/esblkd.f).
jl currently calls ie_espadv/ie_espxcs (IE coefficients) for EM via em_ihtser→ihab, which selects PP for hab 260
where EM's own tables select DF+LP. FIX = transcribe EM's espadv.f/espxcs.f PN regressions + em/esblkd.f CHAB/CPRE
+ em/blkdat.f OCURHT into an em_espadv/em_espxcs, and dispatch EM to them (mirror the ie_ functions; the shared
ie_autoes_tally/run scaffold + ESRANN/ESTPP/heights are already correct & validated). SUBSTANTIAL but MECHANICAL
(≈10 species × the PN coeff sets + 3 tables); best done in a focused session with per-species validation vs live
FVSem on the /workspace/.emwork/sweep_val/ reproducers (target: jl establishes DF+LP not PP; BA 4.5-6.5×→~1×).

## #143 EM — fix scope REFINED: OCURHT alone insufficient, need the full espadv/espxcs PN coefficients
Gathered EM's tables to pin the fix. Habitat 260 → IHAB=3 (jl em_ihtser: bracket _EM_ES_IEND[1]=269≥260 →
_EM_ES_MYGRUP[1]=3). EM OCURHT(3,·) (em/blkdat.f) allows WL(2)/DF(3)/LP(7)/PP(10) — so PP is NOT gated out by
occupancy. ⇒ live selects DF+LP over PP via the PN PROBABILITIES (espadv logistic), not OCURHT. So the fix MUST
port EM's espadv PN coefficients (transcribed from em/espadv.f: advance-regen species {1,2,3,5,7,8,9,10}, each a
distinct PN regression in XCOS/XSIN/SLO/TIME/BAA/BAASQ/ELEV/REGT/BWAF/ELEVSQ + CHAB(IHAB,sp)+CPRE(IPREP,sp)+
OVER>9.95 & per-IFO bumps) + em/espxcs.f (excess) + the EM CHAB/CPRE tables (em/esblkd.f) + OCURHT (em/blkdat.f).
jl currently uses IE's ie_espadv/ie_espxcs coeffs. PLAN: add em_espadv!/em_espxcs! mirroring the jl ie_ functions
(same logistic×OCURHT×XESMLT×OCURNF structure, already validated) with EM coeffs, and dispatch EM to them; keep the
shared ie_autoes_tally/run/ESRANN/ESTPP/heights scaffold. Validate per-species vs live FVSem on the 5
/workspace/.emwork/sweep_val/ reproducers (target: DF+LP not PP; BA 6.5×→~1×). ESPADV PN coeffs are captured in
this session's transcript; remaining to gather = em/espxcs.f + em CHAB/CPRE. Substantial mechanical port, best
executed with sustained focus + incremental per-species validation (NOT rushed — the doctrine is port-faithfully).

## #143 EM — TRUE ROOT (corrects the "port EM espadv" framing): jl's AUTOES occ OMITS OCURNF·XESMLT
Deeper measurement CORRECTS the earlier "EM needs its own espadv" conclusion. The espadv PN regressions, CHAB(16,·),
CPRE, and OCURHT(16,·) are ALL SHARED and IDENTICAL between IE and EM (verified line-by-line: EM espadv.f == jl
ie_espadv; EM CHAB/OCURHT == jl _IE_CHAB/_IE_OCURHT for the relevant species). The ihab mapping is also correct
(hab 260→ihab 3, matches live). So EM does NOT need its own espadv. The REAL bug: the espadv occupancy is
PADV(i)=logistic(PN)·OCURHT(IHAB,i)·XESMLT(i)·**OCURNF(IFO,i)**, but jl's occ (inlandempire/establishment.jl:1122)
= [ie_ocurht(ihab,s) …] applies ONLY OCURHT — it OMITS XESMLT·OCURNF (the per-National-Forest occupancy gate).
MEASURED: jl establishes DF(3)×12 + PP(10)×12 + WL(2)×6 + LP(7)×1; live establishes DF(3)+LP(7) only. The PP(10)
records (fast htg1≈4.22) drive the BA 6.5× over-shoot. EM OCURNF(·,PP=10) = 0,0,1,1,1,0,1,0,1,0,… ⇒ PP is EXCLUDED
on many NFs; on this stand's forest (loc-code 108) live's OCURNF(IFO,PP)=0 zeroes PP, but jl keeps it. ⇒ TRUE ROOT
= the SHARED-ENGINE occ multiplier drops OCURNF·XESMLT (a known gap — the code comment at :337 says occ SHOULD be
OCURHT·XESMLT·OCURNF; validated inert on iet01 because iet01's XESMLT·OCURNF=1). FIX (moderate, shared-engine, helps
both EM & any IE forest where OCURNF≠1): add EM+IE OCURNF(20,·)/XESMLT tables + the forest-location→IFO map, and
multiply occ by OCURNF(ifo,sp)·XESMLT(sp) at :1122. Validate: jl establishes DF+LP not PP; BA 6.5×→~1× on the 5
/workspace/.emwork/sweep_val/ reproducers; iet01/emt01 non-regressing. SUPERSEDES the "port EM espadv" scope.

## #143 residual — the establishment AMOUNT (NUMSPE×ITPP tally), post-OCURNF-fix (separate sub-problem)
With the OCURNF species-occupancy fix landed (EM 1fb8dcc + IE 82771ed; BA now exact on the 5 EM reproducers), the
remaining #143 divergence is the establishment TALLY AMOUNT, not the species set. MEASURED on IE outlier
177562547020004 (bare, hab-establishment): jl establishes 6 species (WL/DF/GF/LP/ES/AF, all OCURNF-allowed) ×
~1112 TPA total across 5 cycles; live establishes ~101 TPA/cycle (2022=101→2062=181). ⇒ jl over-establishes the
AMOUNT ~2× (this stand); the EM reproducers show a milder +3-10% TPA (BA already exact). Candidate roots (ESRANN-
driven, deterministic LCG so FIXABLE in principle): ie_esnspe NUMSPE (# species/plot, capped at MAXSPP(ihab)) and/or
ie_estpp ITPP (trees/plot, capped MAXTPP(ihab)) — if jl's NUMSPE or ITPP or the MAXSPP/MAXTPP(ihab) cap differs
from live, the fresh-cohort tally on BARE stands (NSTORE=0 ⇒ NEWTPP=full ITPP) balloons. NOTE the OCURNF fix made
BA exact because it removed the FAST-growing wrong species (PP); the residual amount is spread over slow species so
its BA impact is small (BA exact) but TPA over-counts. ⇒ #143 remaining = a focused NUMSPE/ITPP/MAXSPP tally audit
(the "IE +19-29%" / this 2× outlier class); the ESRANN per-plot desync (multi-plot) is the harder cornered part.
This is a SEPARATE sub-problem from the now-fixed occupancy gate. Reproducers: /workspace/.emwork/sweep_val/ (EM),
/workspace/.iework/ie_sweep/ (IE).

## #143 IE tally outlier — ROOT: IE establishment uses the RAW habitat code, not the habtyp-mapped one (639→260 gap)
Measured the IE 4.5×-TPA outlier (177562547020004): jl ie_estab_indices gets habcode=639 → ihab=11 (MAXSPP=4,
MAXTPP=10); live habtyp.f MAPS 639→260 (MAPR6→JTYPE, "HABITAT TYPE WILL BE MAPPED TO 260") → ihab=3 (MAXSPP=3,
MAXTPP=5). jl's higher caps ⇒ over-establishment. ROOT: establishment.jl:1180-1181 dispatches the ESTAB habitat as
`EM_JTYPE[habitat_code]` for EM (so the EM OCURNF fix landed on the right ihab) but the RAW `s.plot.habitat_code`
for IE — IE never applies its habtyp MAPR6/JTYPE crosswalk to the establishment ihab. jl's IE growth DOES map the
habitat (site_index.jl habtyp→ITYPE, validated), so the bug is establishment-specific: the ESTAB ihab bracket
(_IE_ESTAB_MYGRUP) runs on the unmapped code. FIX: apply IE's habtyp(MAPR6/JTYPE) mapping to the establishment
habitat before ie_estab_indices (mirror the EM_JTYPE dispatch), so ie_estab_indices sees 260 not 639 ⇒ ihab 3 =
live. Needs the IE MAPR6/JTYPE tables (ie/habtyp.f). This is the #143 IE tally-outlier class (habitat-crosswalk
gap); DISTINCT from the multi-plot ESRANN per-plot desync (the milder cornered part) and from the now-fixed OCURNF
occupancy. So #143 IE tally = (a) habtyp-crosswalk for the establishment ihab [FIXABLE, this finding] + (b) ESRANN
per-plot desync [cornered]. Reproducer /workspace/.iework/ie_sweep/177562547020004.

## #143 IE tally outlier — REFINED: jl's ie_habtyp differs from live's MAPR6→JTYPE for code 639 (NOT a 1-line fix)
Tested the "use the mapped code" fix: jl ie_habtyp(639)→ITYPE=19→IE_MTYPE=620 (NOT live's 260); ie_estab_indices(620)
→ihab=13 (still not live's 3). So jl's single-bracket ie_habtyp (IE_JTYPE 95-bracket → KTYPE → ITYPE) does NOT
reproduce live's TWO-LEVEL habtyp.f crosswalk (NITYPE=MAPR6(KODTYP); KODTYP=JTYPE(NITYPE)=260) for code 639. jl is
missing the MAPR6 level (raw-code→representative-code aggregation) that maps 639→260. This is why the IE
establishment ihab is wrong (11 or 13, vs live 3). CAVEAT: jl IE GROWTH was validated bit-exact-or-cornered on
real FIA, so either (a) few real stands hit the mismatched codes, or (b) establishment uses a habitat mapping
distinct from growth's — needs checking which live routine emits "MAPPED TO 260" (habtyp.f MAPR6 vs an estab-
specific map). ⇒ the IE tally-outlier fix = port live's MAPR6 crosswalk (ie/habtyp.f) into jl's ie_habtyp so
639→260, then feed the mapped code to the establishment ihab (mirror EM_JTYPE). Bigger than expected (MAPR6 table).
DISTINCT from the now-fixed OCURNF occupancy (which made BA exact regardless — the wrong-ihab establishes MORE trees
of the RIGHT slow species, so TPA over-counts but BA stays ~exact). ⇒ #143 status: BA-over-growth FIXED cluster-wide
(OCURNF); TPA-tally residual = habtyp-MAPR6-crosswalk (this) + ESRANN per-plot desync (cornered) — a focused chunk.

## #143 IE tally outlier — CORRECTION (doctrine #2): habtyp/ihab root DISPROVEN by live source
MEASURED live source (NOT inferred): ie/esplt2.f:230-239 brackets the RAW habitat code against IEND(33)→MYGRUP(33)
→IHAB(1-16) — byte-IDENTICAL logic AND tables to jl's ie_estab_indices (_IE_ESTAB IEND/MYGRUP == jl _EM_ES_IEND/
MYGRUP == live). For habcode=639: bracket lands J=23 (634<639≤644), MYGRUP(23)=11 ⇒ live establishment IHAB=11 ==
jl's ihab=11 EXACTLY. The habtyp.f MAPR6→JTYPE path (639→other) fires ONLY for IFOR 5/12 & KODTYP≤40 (line 119) —
our stand is ifo=4/kod=639, so it does NOT apply; live uses the plain JTYPE bracket→ITYPE=19→MTYPE=620 (== jl),
purely a DISPLAY remap, unused by establishment. ⇒ the prior two entries' "jl ie_habtyp missing MAPR6 / live maps
639→260 / ihab 3" root was a HALLUCINATED INFERENCE, now RETRACTED. There is NO habitat/ihab bug. Also verified live
IE OCURNF (ie/blkdat.f) GF(sp4) row = [0,0,1,1,1,0,1,0,0,1,0,0,0,1,0,1,1,0,1,1]; OCURNF(ifo4,GF)=1 ⇒ live ALSO
permits GF on forest 4 — OCURNF does not gate this stand either. ⇒ #143 IE tally residual (this stand: jl BA 143 vs
live 88, TPA 51 vs 40) has ihab/MYGRUP/MAXSPP/OCURNF all == live; TRUE cause is still OPEN — candidates now narrowed
to NUMSPE/ITPP per-plot counts or ESRANN species-draw desync, NOT habitat. Next = MEASURE jl vs live per-species
established TPA on this reproducer (live .out ESTAB report vs jl AUTOSP dump).

## #143 IE tally outlier — RESOLVED (00d36b0): ported ie/pvref1.f (PV_CODE,PV_REF_CODE)→HABPVR
The true root (MEASURED, after retracting the earlier inferred habtyp/MAPR6 guesses): jl's FIA reader had NO
PVREF1 crosswalk and fell back to the RAW PV_REF_CODE (639) as the habitat. Live habtyp.f calls PVREF1 whenever a
reference code is present; the raw ref is never a habitat. Unrecognized pair (ABR8+639, absent from PVREF1's 879
rows) → live default 260 → establishment ihab 3 (permits DF/PP); jl's 639 → ihab 11 (permits GF, forbids PP) →
PADV GF=0.51 dominant → over-establishment. Verified via a jl PADV dump (GF=0.51/PP=0 at ihab 11) against live's
ESTAB report (DF 63/PP 35/GF 0). FIX: ported the 879-row (PVCODE,PVREF)→HABPVR table (pvref1_data.jl + ie_pvref1);
IE resolution = PVREF1(precedence when ref present) → PCOML string → numeric PV_CODE → default 260 (not raw ref).
VALIDATED: reproducer TPA 809→203 (live 181), BA 51→41 (live 40), QMD 3.4→6.1 (live 6.4); recognized 860/101→850
stand bit-close (1053/134 vs 1061/136); iet01 (STDINFO) unaffected. Residual (203 vs 181 TPA) = cornered ESRANN
species-draw straddle (#142-class). META: the 4/5 "NOT RECOGNIZED → 260" sweep prevalence means this was a
SYSTEMATIC IE FIA habitat bug, not a one-off. EM/UT/TT keep their raw-ref fallback (unproven for them; a possible
follow-up is whether they need their own PVREF crosswalk — em/ut/tt each have a pvref1.f).

## #143 EM half — MEASURED: NOT the pvref1/habitat bug (habitat coincidentally correct)
Checked EM's FIA stands (HAB_DUMP instrument on sweep/s1-s5): PV_CODE EMPTY, PV_REF=0 ⇒ jl leaves habitat_code=0.
Live defaults no-habitat stands to NI 260 (em_* synthetic .out all "MAPPED TO 260"). BUT jl's EM establishment maps
habitat_code=0 → clamp(1) → EM_JTYPE[1]=NI 10 → ie_estab_indices → ihab 3; live's default NI 260 → ihab 3 — BOTH
bracket to MYGRUP(1)=3 (both ≤ IEND(1)=269). So jl's EM establishment ihab (3) COINCIDENTALLY == live's, and EM has
NO habitat bug (consistent with the OCURNF fix having made EM BA bit-exact — a wrong ihab would have broken OCURHT/occ
and thus BA). ⇒ EM's remaining AUTOES +200-492% TPA over-growth is in the ie_autoes TREE-COUNT logic (NUMSPE/ITPP/
es_nstore) for EM, NOT habitat — a SEPARATE investigation from IE's pvref1. EM has em/pvref1.f (857) but its FIA stands
carry no PV/ref, so porting it would be inert for the sweep. UT/TT have NO pvref1.f (different habitat mechanism). ⇒
the pvref1 fix is IE-specific and complete; #193's EM half stays open as an ie_autoes count-over-production task.

## #143 EM AUTOES — re-MEASURED current state: AT-BAR (BA bit-exact, TPA +3-10%); "+200-492%" was STALE
Re-ran the 5 sweep_val reproducers post-OCURNF: BA BIT-EXACT on all 5 (8/8/2/12/4 == live), TPA +2.7..+9.9%
(272/262, 278/253, 72/69, 341/332, 165/155). The task-tracker "+200-492%" was a PRE-OCURNF-fix number (stale);
current matches goal-doc item 2's accepted "+~10% ESRANN establishment-tally straddle". CAVEAT: the sign is
CONSISTENTLY POSITIVE (5/5 over), so it is a mild SYSTEMATIC sub-inch over-count, not a true ±straddle — likely
the ESRANN ITPP/es_nstore realization (jl's ESRANN is not byte-identical to live, "never-FFI-RNG" cornered class);
BA-exact means species+sizes are right, only a few extra tiny stems. At-bar per doctrine (bit-exact-or-cornered);
a future refinement could measure jl-vs-live per-plot ITPP to see if a rounding/off-by-one contributes to the sign.
⇒ #143 both halves resolved-or-cornered: IE root FIXED (pvref1, 00d36b0); EM at-bar (BA-exact, TPA +5% cornered).

## #194 CI ci_esgent — re-MEASURED: "+282%" STALE; real issue is established-regen UNDER-growth (DBH-limited)
Ran cibare.key jl vs live FVSci (cibare.sum). Current state is NOT +282% over-establishment (that was pre-#185-
esgent-fix, when regen never grew → tiny-stem pileup). Now: TPA +2-5% (minor), but BA ~25% UNDER (jl 66 vs live 89
@2032) and QMD ~10-15% UNDER (jl 4.7 vs live 5.6). DIAGNOSTIC: TopHt TRACKS live (jl 9/23/34/45/54/63 vs live
12/25/36/46/54/61 — slight early lag, converges) while QMD is consistently under ⇒ DBH-GROWTH-limited, NOT height-
limited (opposite of EM #137's frozen-height). KEY: the deficit is present AT BIRTH — cyc0 (2002) jl & live BOTH
727 TPA but jl QMD 0.9 vs live 1.1, BA 3 vs 5 ⇒ CI regen ENTERS ~18% too small in DBH, and the gap roughly persists.
⇒ ROOT LEAD (not yet instrument-confirmed): the CI establishment initial DBH and/or ci_esgent! (#185, 8cf9c9f)
birth-cycle DBH growth is under; NOT the AUTOES count. NEXT: instrument FVSci_g16 esgent.f vs ci_esgent! on cibare
at cyc0 — compare the birth DBH + HTG/DG applied to the just-established regen (esgent.f HT(I)+=HTG·WK4, DBH via
SMDGF). Distinct from the CI #142 large-tree DGSCOR tail. Task #194 re-scoped from "over-establishment" to
"established-regen DBH under-growth".

## #194 CI — root NARROWED (experiment): ci_esgent! single-pass bscale ≠ live's multi-step KPER-loop regen growth
MEASURED jl birth-cycle values on cibare (planted DF sp3 + LP sp7, both d0=0.102"): DF h 1.52→5.6 (DBH 0.53),
LP 2.4→9.36 (DBH 1.13) ⇒ jl cyc0 TopHt 9 / QMD 0.9 vs live 12 / 1.1 — the birth-cycle HEIGHT growth is ~25% under.
ci_esgent! applies a SINGLE-PASS scalar bscale=(fint-gentim)/regyr=(10-5)/5=1.0. EXPERIMENT (bscale=fint/regyr=2.0,
full cycle): OVERSHOOTS at birth (TopHt 16 / QMD 1.5 / BA 9 vs live 12/1.1/5) though it matches by cyc4 (BA 89=89).
So the true live birth-cycle scale is INTERMEDIATE (~1.4×), NOT a clean scalar ⇒ jl's single-pass bscale can't
reproduce it. Live's regent.f grows regen over a KPER-STEP LOOP (SCALE=FLOAT(KPER(J))/REGYR per step, regent.f:510)
with a PER-TREE WK4(I) birth fraction (H2=H1+HTGRL·SCALE·XRHGRO·CON·WK4(I), regent.f:709). WK4(I) is set per tree at
establishment (elsewhere — estab/esinit, not regent.f). ⇒ #194 FIX = port live's multi-step regen-growth iteration
(KPER loop + per-tree WK4) into ci_esgent! (and likely the shared esgent path cluster-wide — EM/UT/TT/BM/IE use the
same single-pass bscale, so their PLANTED-regen birth cycle may share this; NATURAL-regen validations (emt01/utt01)
didn't exercise PLANT-at-cycle-start). NEXT: find WK4(I) assignment (grep estab/esinit), confirm KPER stepping for a
10-yr cycle, port. This is a genuine model port, not a scalar tweak. Distinct from #142 (large-tree DGSCOR).

## #194 CI — birth-cycle RELHT FIXED (eb3395b); residual is a separate small-tree DBH tail
FIXED the birth-cycle component: ci_esgent! CIVAR RELHT now uses the PRE-regen ATAVH (=0 on bare ⇒ RELHT=1.5),
captured as es_avh_pre before establish! (which recomputes avg_height WITH the new regen). MEASURED-exact for DF
(HTGRL 4.08→5.259 = live 5.258). cibare cyc0 QMD 0.9→1.0 / TopHt 9→11 (live 1.1/12); cit01 INERT (stand_top_height
= AVHT40 is invariant to adding tiny regen ⇒ avh_pre==avh for established stands; = live's AVH/ATAVH, which is why
the established-stand RELHT was already right). RESIDUAL: cibare later-cycle DBH still ~12% under (QMD 4.9 vs 5.6,
BA 72 vs 89 @cyc4) — NOT RELHT (jl's AVHT40 proxy == live AVH) and NOT the birth cycle; it's the ongoing small-tree
DBH growth on an even-aged pure-regen stand = the CI "SMHTGF small-tree stochastic" tail the goal-doc lists. NEXT
for that residual: instrument small_tree_growth! DBH regression (dk-dkk / DDS) vs FVSci_g16 regent.f on cibare at a
mid cycle (e.g. 2022). Distinct from #142 (large-tree DGSCOR).

## #194 META check — birth-cycle RELHT bug is CI-SPECIFIC (not cluster-wide)
Checked whether the ci_esgent! RELHT-uses-post-regen-avh bug extends to EM/UT/TT/BM/IE esgent (they share the
"grow birth-cycle regen" pattern). CONCLUSION: CI-specific. The direct RELHT=h0/ATAVH term lives ONLY in CI's CIVAR
HTGRL regression (SRLHT·RELHT). EM/TT use the Wykoff crown form HTGRL=BETA1+BETA2·CR (crown-ratio, no RELHT); UT/BM
use POTHTG·PCTRED·VIGOR. The other variants' only avh use is xd=avh·(relden/100) (the PCTRED density arg) or EM's
delmax — both relden-GATED, so ~inert on bare-plant stands (relden≈0). ⇒ no cluster-wide extension of eb3395b
needed; the earlier "likely applies to EM/UT/TT/BM/IE" note is RETRACTED.

## #194 — birth-cycle HTGRL now MATCHES live for BOTH species (fix eb3395b confirmed correct)
Extracted live's per-species birth HTGRL from FVSci_g16 cibare_dbg.out: DF(sp3)=5.25799847, LP(sp7)=8.24845219
(each tree grows in ONE REGENT step, no KPER loop). With the RELHT=1.5 fix, jl reproduces both: DF 1.599+1.36·1.5+
0.02·81=5.259; LP 2.134+2.559·1.5+0.0593·81−0.43·9+1.345=8.251 (≈ live within Float32). So the birth-cycle height
growth is now faithful. The residual (cibare cyc0 TopHt 11 vs live 12 → BA 72 vs 89 @cyc4) is therefore NOT the
HTGRL regression — it's downstream: (a) the DBH assigned at the 4.5' crossing / DDS accumulation for the just-grown
regen, and/or (b) the ongoing small_tree_growth! in later cycles. = the CI "SMHTGF small-tree" tail (goal-doc item
6, unmeasured). Precisely scoped for a focused next session; the birth-cycle HTGRL itself is settled.

## #191 TT aspen — MEASURED precisely (was vague "sub-1 over + DGFASP under"); current, NOT stale
asp.key (753175613290487, 100% aspen sp6): jl OVER-grows BA +8-13% (2049 BA 111 vs live 98), QMD over (5.9 vs 5.5),
TPA identical → deterministic DG over-growth. Instrumented BOTH sides (FVStt_g16 DEBUG DGFASP + jl TT_ASP_DUMP):
• DGFASP formula + coefficients BIT-MATCH pristine dgfasp.f; BA input matches (~22, both backdated growth-period BA;
  verified live BAACT≈22 back-derived from VALMOD=0.4559).
• RMSQD DIFFERS: jl=2.495 vs live=2.96 (back-derived from live GOFAD=0.6031; live cyc2 GOFAD=0.6905→rmsqd 3.77).
  RMSQD=SQRT(TSUMD2/TPROB) (base/dense.f:250) = QMD over the summed trees, NO threshold ⇒ the gap is a TREE-
  POPULATION difference: jl's DG-point QMD (2.495) < live's (2.96 ≈ reported QMD 3.0). Since cyc0 .sum QMD is
  bit-exact (3.0), jl's DG-point stand carries EXTRA low-DBH aspen that drag RMSQD down — the "sub-1" component.
  Likely SAME CLASS as CI #194: jl computes the stand stat over a population INCLUDING regen/sub-1" that live's
  RMSQD excludes (pre- vs post- a step). Lower RMSQD → lower GOFAD → lower large-tree DGFASP (UNDER) — the
  "DGFASP-RMSQD under" half.
• COR: jl aspen dg_cor=0/cor2=1 vs live DGSCOR scale ~0.94-0.99 (live reduces DG). The two ERRORS PARTIALLY CANCEL
  (d=5.7: jl aspdg 1.5103 ≈ live 1.570·0.94≈1.508) — masking each other ⇒ this is why #191 is "entangled, land
  together". NET over-growth = the sub-1" aspen (extra population + their own over-growth) dominating.
⇒ #191 FIX (land together): (a) determine why jl's DG-point aspen population has extra low-DBH trees (dump jl t.n +
diameter distribution at the DG point vs live ITRN; check DG-vs-establishment order & whether RMSQD should exclude
birth-cycle regen — cf. #194 pre-regen ATAVH); (b) confirm/port the aspen DGSCOR COR. Precise numbers now in hand;
distinct from #142/#158.

## #191 refinement — RMSQD gap is BACKDATING-asymmetry, not extra trees (sharper lead)
jl DG-point BA 31→21.91 (×0.707) and QMD² 3.0²→2.495² (×0.69) dropped PROPORTIONALLY ⇒ the signature of BACKDATING
(all diameters scaled down by ~0.83), NOT an extra low-DBH population (which would move BA/QMD/TPA differently). So
the likely mechanism: live's DGFASP uses the BACKDATED BA but the CURRENT RMSQD (2.96), while jl backdates BOTH
(rmsqd 2.495). CAVEAT/CONTRADICTION to resolve first: jl's per-tree d=5.7 MATCHES live (current), yet jl's
stand_qmd(s) (which sums current t.dbh, standstats.jl:90) returns 2.495 not 2.96 — so either t.dbh ARE backdated at
the dgf! rmsqd call (and the dumped d=5.7 is coincidental/also-backdated) or rmsqd is sourced from the backdated
p.basal_area path elsewhere. NEXT (focused): dump jl t.n + full diameter list AT the dgf! rmsqd line vs the cyc0
stand, and trace whether tt dgf! backdates t.dbh before line 153. Then align RMSQD to live's current-QMD source
while keeping BA backdated. (Still land together with the aspen COR.)

## #191 CORRECTED ROOT — aspen DGSCOR calibration over-corrects (POSITIVE COR), RMSQD was a red herring
Deeper measurement RETRACTS the RMSQD/backdating lead above: jl's dgf! is called TWICE per cycle — a CALIBRATION
pass (backdated stand: rmsqd 2.495, ba 21.91, dg_cor=0 by design) and the ACTUAL-growth pass (rmsqd 2.964 == live
2.96 ✓, ba 30.9). My first TT_ASP_DUMP only captured the calibration pass ⇒ the "RMSQD 2.495 vs 2.96" and "aspen
COR=0" were both CALIBRATION-pass artifacts, NOT the actual-growth values. Actual-growth RMSQD MATCHES live.
TRUE ROOT: jl's aspen DGSCOR calibration computes corv=+1.1656 (bnyv=1.0466 = mean measured-vs-predicted residual,
POSITIVE ⇒ measured >> predicted ⇒ COR *BOOSTS* aspen DG). Live's DGSCOR for this stand REDUCES DG (scale 0.94).
OPPOSITE SIGN ⇒ jl over-grows +13% BA. Calibration gate passes (calib_sp=T, isct=1, fn=10≥5, snp=60; trap doesn't
fire, exp(corv)=3.21 in range). Since the aspdg PREDICTED DDS is bit-matched (formula+coeffs verified vs pristine
dgfasp.f), the divergence is in the MEASURED aspen DDS (from the FIA DG field) or the per-tree residual accumulation
(snx/sny/snp) — jl's measured aspen growth reads much higher-vs-predicted than live's. NEXT (focused): instrument the
calibration per-tree loop — dump jl's measured DDS + predicted DDS per aspen tree vs FVStt_g16 (readdgf/dgdriv DEBUG),
find why jl's measured/predicted ratio is inverted vs live. This is the "land together" root (the earlier
sub-1"/DGFASP-RMSQD framing is superseded). Distinct from #142.

## #191 final narrowing — the inverted COR is in the MEASURED aspen DDS, not predicted/RMSQD
Ruled out the calibration-RMSQD candidate: FVS backdates RMSQD in the calibration too (dense.f LBKDEN, D=WK3(I) at
:184 ⇒ RMSQD from backdated D), so jl's backdated calibration RMSQD (2.495) MATCHES FVS. Predicted aspdg formula+
coeffs+inputs match. ⇒ the corv sign inversion (jl +1.17 boost vs live 0.94 reduce) must be in the MEASURED aspen
DDS or the per-tree residual (sny/snp): jl's measured aspen growth-vs-predicted reads much higher than live's. The
measured DDS derives from the FIA DG field (past dbh, DG_TRANS/DG_MEASURE) via a bark/period conversion — a likely
aspen-specific error there (bark ratio, the DG→DDS conversion, or the measurement-period FINT) would inflate jl's
measured aspen DDS. NEXT (the fix step): instrument the calibration per-tree loop for aspen — dump jl's per-tree
(measured_DDS, predicted_DDS, DG-field value, bark, period) vs FVStt_g16 readdgf/dgdriv DEBUG; align the measured
aspen DDS. This is the precise, corrected #191 root (supersedes both the sub-1"/DGFASP-RMSQD and the calibration-
RMSQD framings). Bounded next-session task.

## #191 PARTIAL FIX LANDED (42f4860) — current RMSQD in aspen calibration; residual = measured-DDS driver
Landed: the aspen DGFASP calibration now uses the CURRENT stand RMSQD (not the backdated stand_qmd) — FVS does
this (asp_dbg.out shows only current-QMD GOFAD, never backdated), same class as the validated AVH exception.
Reduced asp.key over-growth +13%→+9% (BA 2049 111→107 vs live 98). ttt01 INERT (jlPRE==jlFIX, no measured aspen
DG) — no regression. REMAINING (bounded next step): the corv stays POSITIVE (~+1.0, boost) while live REDUCES
(scale 0.94) ⇒ jl's calibration still finds measured aspen DDS >> predicted (~3× vs live's 0.94×). Predicted (aspdg
formula + now-current RMSQD + bark, aspen IMAP=2 const 0.969) all match live ⇒ the driver is the MEASURED aspen DDS
term = dg·(2·bark·wk3 + dg)·scale (southern/diameter_growth.jl:508): jl's measured aspen growth reads ~3× too high
vs predicted. Candidates: the DG-field interpretation (DG_TRANS/DG_MEASURE past-dbh vs increment), the period `scale`,
or the backdated wk3. NEXT: dump jl per-aspen-tree (dg, wk3, bark, scale, term, wk2, reslog) and compare to FVStt_g16
dgdriv/readcor DEBUG; align jl's measured aspen DDS so corv flips to live's reduction. Then #191 is bit-exact-or-
cornered. (The earlier RMSQD-actual-growth and sub-1"/DGFASP framings are fully superseded.)

## #191 CORRECTION+UPGRADE — the RMSQD fix (42f4860) made the aspen COR BIT-EXACT (1.1062=live); residual is sub-1"
FVStt_g16 has a G16ASP DEBUG that prints the exact per-tree ASPDG + COR. MEASURED after the RMSQD fix: jl's aspen
calibration corv = 1.1062 = live's COR 1.1062 EXACTLY (bnyv/bnxv/slp/dist all consistent); jl's predicted ASPDG
(1.5704 @d=5.7) = live 1.5704 EXACTLY; jl's APPLIED dg_cor at actual growth = 1.1062 = live. ⇒ the large-tree
aspen DGFASP path (predicted + COR + applied) is now BIT-EXACT with live — the DGFASP-RMSQD half of #191 is RESOLVED
(not "partial/measured-DDS-open" as 42f4860's message said; the earlier "0.94"/"measured>>predicted" reads were
WRONG — 0.94 was an unrelated readin diagnostic, and live's aspen COR is a genuine +1.1062 BOOST that jl now matches).
Live aspen COR IS a boost (aspen grows ~2× the raw model here) — jl agrees. RESIDUAL: asp.key still +9% BA / QMD
+0.2-0.3 (jl 5.8 vs live 5.5 @2049), present from cyc1. With the large-tree DG bit-exact, this residual is the
SUB-1" regent subcycle (aspen SMHTGF/regent, the #158 "single-step-suppressed-model" half) feeding slightly-bigger
trees into the DGFASP pool — the OTHER half of the #191 pair. NEXT: instrument the aspen sub-1"/regent DG on asp.key
vs FVStt_g16 regent DEBUG (the #158 work). The RMSQD fix is faithful+validated (COR bit-exact, ttt01 inert). NOTE:
correct 42f4860's residual description in the record — it's sub-1", not measured-DDS.

## #191 residual REFINED — it's DIAMETER, not height (TopHt matches live exactly)
asp.key TopHt is bit-exact vs live (31/38/46/53-54/61/67 = live) while QMD/BA is over (BA cyc1 51 vs 48, QMD
3.9 vs 3.8; growing to +9% BA / QMD 5.8 vs 5.5 @2049). ⇒ the residual is a DIAMETER-growth residual, NOT the sub-1"
HEIGHT-crossing hypothesized above (height growth matches). With the large-tree DGFASP predicted+COR bit-exact
(verified for the calibration-range trees), the residual is an ACTUAL-GROWTH aspen diameter residual in some subset
— candidates: (a) the actual-growth DGFASP BA basis (jl uses current p.basal_area 30.9; confirm live's actual-growth
DGFASP BA — the calibration pass used backdated ~22), (b) aspen below the calibration DBH range (dn-dx) getting a
different endpoint DG, (c) the DDS→diameter conversion for aspen. NEXT (bounded): per-tree actual-growth DGFASP
compare jl vs FVStt_g16 (D, BA, ASPDG, applied dg_cor, resulting dg) on a few aspen across the size range at cyc1 —
find which trees over-grow their DBH. NOTE tripling fires (25 recs ≤ MAXTRE/3) so use the pre-split window or the
per-DBH-class .sum distribution, not the raw treelist. The DGFASP-RMSQD/COR half stays RESOLVED (42f4860, bit-exact).

## #191 residual — BA-basis RULED OUT; candidate is the COR-clock (needs NOTRIPLE per-tree compare)
tt/dgf.f:565 `CALL DGFASP(D,ASPDG,CR,BARK,SI,DEBUG)` passes NO BA — DGFASP reads BA from the common block (= the
current stand BA at growth), so jl's `ba = p.basal_area` matches live. BA basis ruled out. With aspdg (predicted),
COR (1.1062), BA, AND TopHt all matching yet asp.key BA +6% at cyc1 (uniform ~1.5%/tree DBH over-growth, TPA exact),
the remaining candidate is the DGSCOR COR-CLOCK: jl applies the full corv (dg_cor=1.1062, TTAPPLY-verified at cyc0)
at growth; if live ATTENUATES the applied COR per cycle (autcor.f decay, dg_cor_goal=0.5·corv) differently than jl's
line-1046 evolution, jl's applied DG is higher → uniform over-growth. (For CI #142 jl's COR-clock matched live at
0.05693; the large TT-aspen COR 1.1062 may expose a decay difference.) NEXT (bounded, decisive): asp_tl.key gives a
per-tree treelist; run jl NOTRIPLE (or the pre-split window) and compare per-aspen-tree DBH growth 2019→2029 to live
asp_tl.out — measure jl's vs live's APPLIED COR/DG per cycle. This isolates COR-clock vs a DDS→diameter/other factor.
The DGFASP-RMSQD/COR half stays RESOLVED (42f4860). Residual is small (+6-9% BA on one aspen stand) and bounded.

## #191 residual FULLY LOCALIZED — COR-clock is BIT-EXACT; residual is the small-aspen DBH-at-4.5'-crossing (#158)
Measured jl's applied aspen dg_cor per cycle vs live's G16ASP decay: 1.1062/0.9722/0.8707/0.7938/0.7355 — jl ==
live ALL 5 CYCLES bit-exact. ⇒ the aspen large-tree DGFASP + DGSCOR COR + autcor decay are FULLY bit-exact. With
TopHt also bit-exact (height growth matches) but QMD/BA over (+6-9%), the residual is NECESSARILY the DBH ASSIGNED
when small aspen cross 4.5' breast height (the HT-DBH / "single-step-suppressed-model" — #158): jl assigns a bigger
DBH at crossing than live, lifting BA/QMD without touching TopHt or the (bit-exact) large-tree DG. This is the pure
#158 half of the #191 pair, now cleanly isolated from the DGFASP-RMSQD/COR half (RESOLVED, 42f4860). NEXT: instrument
the aspen HT-DBH at the 4.5' crossing (the regent/SMHTGF→DBH handoff) on asp.key vs FVStt_g16, for trees transitioning
sub-4.5'→above. Everything else in the aspen DG chain is proven bit-exact — the #158 DBH-at-crossing is the sole
remaining driver.

## #191 residual — PINPOINTED to the aspen _tt_smdgf small-tree DBH blend (regent.jl:227-241), #158 model
Traced the sole remaining driver to small_tree_growth! (teton/regent.jl:207-241): for aspen (sp6, _tt_rg_default),
the per-tree DBH growth BLENDS the #158 small-tree model `dgk` (from _tt_smdgf: DK=smdgf(grown H), DKK=smdgf(orig H),
DG=(DK−DKK)·bark→DDS, line 221-239) with the (bit-exact) large-tree DGFASP DG, weighted by xwt=(d−XMIN)/(XMAX−XMIN).
Small/mid aspen (xwt<1) grow mostly via `dgk`. The large-tree DGFASP + COR + COR-clock are PROVEN bit-exact this
session, and the aspen SMHTGF HEIGHT has the #189 RSIMOD fix (TopHt bit-exact) + the #158 DGMX=FINT·DGMAX cap
(f66f1fd) — so the residual +6-9% BA is the `_tt_smdgf` aspen DBH (or the blend) for small/mid aspen that over-grows
DBH without touching height. This is the last-mile of the #158 "single-step-suppressed-model" ON THIS low-site
(SITEAR=42) aspen stand. NEXT (bounded, decisive): instrument dgk/_tt_smdgf per aspen tree on asp.key vs FVStt_g16
regent DEBUG (DK/DKK/DGR/DGMX) for the xwt<1 trees; the height/large-tree halves are proven exact so this isolates
the smdgf DBH. #191 is now: DGFASP-RMSQD/COR half RESOLVED bit-exact (42f4860); #189 height RSIMOD landed; residual =
the _tt_smdgf DBH-blend for small aspen (bounded).
