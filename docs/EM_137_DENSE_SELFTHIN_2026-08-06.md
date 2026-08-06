# EM #137 dense-cohort self-thin — REPRODUCED + characterized (2026-08-06)

## Reproduction (constructed dense PLANT scenario)
em_plant.key plants only 400 TPA (non-self-thinning; jl TPA BIT-EXACT there, see #137 note). To exercise the
dense-cohort self-thin, cloned it to em_plant_dense.key with `PLANT 1990 3.0 6000.0` (6000 TPA). jl vs FVSem_clean:

  year   liveTPA  jlTPA   liveBA  jlBA
  2000   5455     5455    2       0
  2010   5350     5350    14      9
  2020   5246     5247    30      23
  2030   5145     5145    51      39
  2040   5045     5045    78      58     <- TPA still bit-exact, BA already −26%
  2050   4810     4947    108     81     <- TPA diverges: jl UNDER-thins
  2060   4354     4783    133     106
  2070   4001     4369    158     127
  2080   3699     3973    182     154
  2090   3167     3566    190     184    <- jl TPA +12.6%, BA −3% (converged)

## Characterization
- DIRECTION: jl UNDER-thins (TPA HIGH: 3566 vs 3167 = +12.6% by 2090). This CORRECTS the goal-file's "over-kill"
  label — it is UNDER-kill, the SAME direction as BM #140 (resolved).
- ROOT is UPSTREAM of mortality: TPA is bit-exact 2000-2040 while BA is already one-directional LOW (58 vs 78 at
  2040, −26%). So jl under-grows the dense cohort's BASAL AREA first; the self-thin under-kill (2050+) is the
  DOWNSTREAM consequence (smaller QMD → higher T85D10/BAMAX self-thin target → fewer trees killed → QMD-feedback).
- The BA lag is a TIMING lag (converges to −3% at 2090): the dense planted small-tree cohort grows diameter slower
  in jl early, so it crosses the small→large-tree threshold later. NOT the crown-init gap (planted trees get crowns
  at PLANT, crown≠0) — it is the small-tree DG/HTG under DENSE competition (SMHTGF/SMDGF) or the DGSCOR calibration
  on the dense cohort.

## Next (fresh chunk)
Instrument jl small-tree DG + HTG (em/regent.f SMHTGF/SMDGF) vs live on the dense planted cohort at an early cycle
(2010/2020) — measure per-tree height + diameter growth to see whether jl under-grows HEIGHT (delaying the 4.5'
crossing) or DIAMETER directly, under the high-competition (dense TPCCF/BAL) inputs. Whether it is a real bug or the
accepted dense-cohort small-tree timing tail (cf. the CORNERED EM/IE growth-only ~7% BA tail) depends on that
measurement — it converges by 2090, which leans cornered, but the +12.6% TPA under-thin is a notable persistent
divergence worth the instrument pass. Harness: em_plant_dense.key in /workspace/.emwork (live) + run_keyfile (jl).

## KEY CLUE — the lag starts at cyc0 (planted-tree small-tree HTG timing)
The BA gap is present from the FIRST reported cycle: 2000 (10 yr after the 1990 PLANT) live BA=2 (QMD≈0.26",
5455 TPA) vs jl BA=0. So jl's 6000-TPA planted cohort does NOT reach measurable breast-height DBH in the first
10-yr cycle while live's does — jl's small-tree HEIGHT growth crosses the 4.5' breast-height gate LATER, so
diameter growth (which only starts >4.5') is delayed ⇒ BA lag ⇒ later self-thin ⇒ under-kill. Same 4.5'-threshold
mechanism as BM #149, but here for PLANTED trees under DENSE competition (crown≠0, so NOT the crown-init — it's the
SMHTGF height RATE under high TPCCF/BAL). ⇒ the instrument target is em/regent.f SMHTGF htg1=beta1+beta2·CR under
the dense-cohort competition inputs: does jl's height growth (or the competition attenuation of it) run low vs live?

## Root NARROWED — small-tree HTG RATE under dense competition (initial height RULED OUT)
- jl's EM planted base height (em/essubh.f, HHT=EXP(PN)) is ALREADY validated bit-exact vs FVSem on em_plant.key
  (DF PN=0.156307 HHT=1.18492; establishment.jl:4). And essubh clamps BAA∈[1,400] — at cyc0 both the 400- and
  6000-TPA cohorts have BA≈0 ⇒ BAA clamps to 1 ⇒ IDENTICAL initial HHT. ⇒ the divergence is NOT the planted
  initial height.
- em_plant (400 TPA, low competition) is TPA-bit-exact + BA converges; em_plant_dense (6000 TPA) diverges. The ONLY
  difference is DENSITY. ⇒ jl's small-tree HEIGHT growth (em/regent.f SMHTGF) OVER-SUPPRESSES under high competition
  (dense TPCCF/BAL/RELDEN): over the first 10-yr cycle jl's cohort grows shorter → stays <4.5' → dbh/BA=0 at 2000
  (vs live BA=2) → BA lag compounds → later self-thin → under-kill (+12.6% TPA by 2090).
- DECISIVE NEXT MEASUREMENT: instrument-replay em/regent.f SMHTGF on em_plant_dense.key (live) — dump the per-tree
  small-tree HTG and its competition inputs (TPCCF/BAL/CCF/attenuation) at cyc0 — vs jl's em small-tree HTG for the
  same tree. Whichever competition term jl applies differently (or an attenuation/cap jl over-applies at high density)
  is the root. Then decide real-bug-vs-cornered. (Doctrine #2: MEASURE; #3: use .sum/pre-split, cohort is tripled.)

## PRECISE SOURCE TARGET (2026-08-06 cont.) — em/regent.f:468 sub-BH HTGRL, NOT SMHTGF
Ruled out SMHTGF: em/smhtgf.f:66-70 returns HTGRTH=0 when DBH≤0. The cyc0 planted cohort is dbh=0 (sub-breast-
height), so it grows height via the em/regent.f sub-BH path, NOT SMHTGF. That path:
  regent.f:468  HTGRL = CON + BH·ALOG(H1) + BCCF·RDJ + BBAL·BAL     (H1<4.5' "NI section")
Competition enters via **RDJ = RDNEXT(J)** (the per-period PROJECTED relative-density CCF, built at regent.f:237-278
from RELDEN/TEMCCF/CCFYR) and **BAL** (BA in larger trees). On the 6000-TPA cohort RDJ+BAL are high ⇒ BCCF·RDJ +
BBAL·BAL (both coeffs negative) suppress HTGRL. If jl's RDNEXT trajectory or BAL on the dense cohort differs from
live, HTGRL under-grows ⇒ trees stay <4.5' ⇒ dbh/BA=0 longer ⇒ the observed cyc0 BA=0-vs-2 and the compounding lag.
DECISIVE INSTRUMENT (fresh chunk): patch em/regent.f:468 to WRITE I,H1,RDJ,BAL,CON,BH,BCCF,BBAL,HTGRL → FVSem_trc on
em_plant_dense.key cyc0; dump jl's em sub-BH HTGRL + its RDNEXT/BAL for the same cohort; the diverging competition
input (RDNEXT projection or BAL) is the root. If RDJ/BAL/HTGRL all match and only the post-4.5' SMHTGF ZRAND realization
differs, it's the accepted dense-regen stochastic straddle (cornered) — but the one-directional cyc0 BA=0 lag points at
the DETERMINISTIC HTGRL competition terms, i.e. a real RDNEXT/BAL discrepancy worth fixing.

## CORRECTION (measurement REFUTES regent.f:468) — it's the EM establishment height subsystem
Instrumented em/regent.f:468 (NI-section HTGRL) → dumped to unit 16 on em_plant_dense.key → **0 hits** (even ungated).
So the planted cohort's height growth does NOT go through regent.f:468 at all. My source inference was WRONG; the
measurement caught it (doctrine #2). The PLANTED-tree height trajectory is driven by the EM ESTABLISHMENT height
routines — em/esadvh.f (advance/subsequent height), em/espsub.f, em/essubh.f (base height), em/esdlay.f (delay) — NOT
the regent.f small-tree loop (which is for the REGEN blend / read small trees). So #137's cyc0 BA=0-vs-2 lag is in the
EM ESTAB subsystem's per-cycle height growth under dense competition, and jl's essubh BASE height is validated only at
cyc0/low-density (BAA clamp) — the SUBSEQUENT/advance-height growth (esadvh) under the dense cohort's rising BAA is the
unvalidated path where jl likely diverges. REVISED instrument target: em/esadvh.f + espsub.f (dump the per-cycle
established-tree height + its BAA/density inputs) vs jl's EM establishment height on em_plant_dense.key. This is a
distinct subsystem from the small-tree regent path — a fresh, well-scoped chunk. NET this session: #137 reproduced,
direction corrected (under-thin), and localized to the EM estab height subsystem with regent.f/SMHTGF/initial-height
all RULED OUT by measurement.

## MEASUREMENT CHAIN (2026-08-06 cont.) — per-cycle growth is em/htgf.f, NOT the small-tree/estab paths
Instrument-replay on em_plant_dense.key ruled out, by ZERO hits each, the paths I'd hypothesized:
 - regent.f:468 (NI-section HTGRL) — 0 hits.
 - regent.f:495 (TT-section HTGRL=HTG1+ZRAND) — 0 hits.
 - smhtgf.f:75 (HTGRTH=HTG1+ZRAND·STDDEV) — 0 hits (SMHTGF returns early at dbh≤0).
 - REGENT is called ONLY from cratet.f:102/553 (INIT-time CALIBRATION of the small-tree model), NOT per-cycle; and
   em_plant_dense has no inventory small trees at cratet ⇒ REGENT's growth loop never executes. SMHTGF/SMDGF are
   called ONLY from regent.f ⇒ also init-only.
⇒ The planted cohort's per-cycle HEIGHT growth is em/htgf.f (the MAIN height model), which has a small-tree branch at
htgf.f:225 `HTG=EXP(CON+HDGCOF·ALOG(DG))+0.4809` (height-from-diameter-growth) plus the POTHTG potential-height path
(htgf.f:191-217, crown-modified). So #137's dense-cohort BA lag lives in em/htgf.f (and its DG input from dgf.f) under
the dense cohort — NOT the establishment/regen small-tree subsystem (that only ADDS + calibrates at init). This
also means jl's routing must match: if jl grows the small planted trees via its em small-tree (regent-port) model
while live grows them via htgf.f main, that ROUTING mismatch is the root. NEXT (decisive): instrument htgf.f:217/225
(dump ISPC,DBH,HT,DG,PHTG/HTG) on em_plant_dense cyc1 + check jl's EM per-cycle growth dispatch for a dbh<1 planted
tree (does jl call the small-tree model or the main htgf-equivalent?). The routing/branch discrepancy is the fix.
NET this iteration: #137 reproduced + direction corrected + FIVE candidate paths ruled out by measurement, localized
to em/htgf.f main-height (+ its DG) vs jl's small-tree routing. Doctrine #2 in action (each inference measured, not assumed).

## CORNERED to the establishment essubh model + a likely jl ROUTING MISMATCH (2026-08-06 cont.)
Read em/htgf.f: for EM species (ISPC≤3,7-10,18) line 185 `IF(H.LE.4.5)GO TO 60` sends sub-breast-height trees to
label 60 with HTG(I) STILL 0 (initialized htgf.f:166) ⇒ **htgf gives sub-BH trees ZERO height growth**. Combined with
the measured facts (regent/SMHTGF are cratet-init-only, not per-cycle; htgf=0 for H≤4.5), the ONLY remaining per-cycle
height driver for the sub-BH PLANTED cohort is the EM ESTABLISHMENT model (essubh/esadvh height-BY-AGE, which re-heights
established trees each cycle until they exceed 4.5'). 
★ LIKELY ROOT (jl routing mismatch, to be MEASURED next): jl's EM height growth (height_growth.jl:106) routes ALL
h≤4.5 trees to small_tree_growth!/_em_smhtgf (the stochastic small-tree model). But live grows ESTABLISHED sub-BH
trees via essubh-by-age (deterministic), NOT the small-tree model — that model (regent) only runs at cratet-init.
So jl likely GROWS the planted sub-BH cohort with the wrong model (small-tree stochastic) where live re-heights via
essubh(age); under dense competition these diverge ⇒ the cyc0 BA=0-vs-2 lag ⇒ under-thin. (This wouldn't hit READ
small-tree stands the same way — those ARE inventory trees the small-tree model is meant for; the issue is specifically
ESTABLISHED/planted sub-BH trees.) DECISIVE NEXT: instrument em essubh/esadvh per-cycle on em_plant_dense (does live
re-height the planted cohort via essubh each cycle?) + confirm jl routes them to _em_smhtgf ⇒ the fix is to re-height
jl's established sub-BH trees via em_essubh_hht(age) instead of the small-tree model.
NET (2 iterations): #137 reproduced, direction corrected (under-thin), SIX candidate paths eliminated by measurement,
cornered to the establishment-vs-small-tree ROUTING of sub-BH planted trees. All by instrument-replay (doctrine #2).

## CORRECTION #2 (measurement) — REGENT IS per-cycle (base/grincr.f:449); crux = first-cycle BH crossing
My "REGENT is cratet-only" claim above was WRONG — it came from grepping only em/. Widening to base/: base/grincr.f:449
`CALL REGENT(.FALSE.,1)` runs REGENT PER-CYCLE (grincr.f:445 "CALL REGENT TO COMPUTE HEIGHT AND DIAMETER INCREMENT").
So regent.f IS the per-cycle small-tree growth. Within regent.f the EM species take the ELSE/"EM VARIANT" branch
(regent.f:540-546) → `CALL SMHTGF` → H2=H1+HTGRR·(KPER/REGYR)·XRHGRO·CON. SMHTGF returns HTGRTH=0 for DBH≤0
(smhtgf.f:66-70), so my SMHT dump (after the HTG1+ZRAND line) legitimately never fired for the sub-BH cohort.
NET corrected picture: BOTH htgf (H≤4.5 → HTG=0) AND regent-SMHTGF (dbh≤0 → HTGRTH=0) give the sub-breast-height
planted trees ZERO height growth. So the trees must get from the establishment HHT (1.18') to >4.5' via ANOTHER
mechanism — the establishment model's per-cycle height (estb/estab.f CALL ESSUBH/ESADVH re-heighting the regen cohort
each ESTAB cycle), OR a regent sub-BH path still not instrumented, OR the planted trees are seeded above 4.5'. The
DIVERGENCE is the FIRST cycle (1990→2000): live's cohort reaches dbh>0 (BA=2) by 2000, jl's stays dbh=0 (BA=0). ⇒
the decisive measurement is PER-TREE first-cycle height (not .sum): instrument estb/estab.f ESSUBH/ESADVH height +
the regent first-pass H2 for the planted cohort at 1990→2000, vs jl's established-tree height at 2000. This needs
jl-side instrumentation too (the .sum aggregates can't show the per-tree BH-crossing). A dedicated chunk.
★ HONEST STATUS: #137 reproduced + direction corrected; the per-cycle growth dispatch fully mapped (regent per-cycle
via grincr; EM→SMHTGF; both height paths zero for sub-BH); crux localized to the first-cycle establishment-height
BH-crossing. One self-correction logged (REGENT per-cycle, not cratet-only) — doctrine #2 caught it.
