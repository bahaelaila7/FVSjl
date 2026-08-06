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

## ★★ REFRAMED (per-tree TopHt measurement) — #137 is a HEIGHT-GROWTH bug, general (not density-specific)
Compared TopHt (tripling-safe stand aggregate) at 2000 (age 10, 10 yr after the 1990 PLANT):
  em_plant_dense (6000 TPA): live TopHt=6  jl TopHt=1
  em_plant      (400  TPA): live TopHt=5  jl TopHt=1
So jl UNDER-grows the EM planted cohort's HEIGHT at BOTH densities — it is NOT density-specific. jl's cohort sits at
~1' (≈ the essubh establishment HHT 1.18') while live's grows to 5-6' by 2000. Since jl's essubh base height is
validated bit-exact (establishment.jl:4), live's trees are growing ABOVE the establishment height via the per-cycle
tree model, whereas jl's stay STUCK at it (~0 height growth). ⇒ #137 is fundamentally a HEIGHT-GROWTH bug: jl gives
the EM planted/established SUB-breast-height cohort ~0 per-cycle height growth where live grows ~0.5'/yr (1.18'→5-6'
in 10 yr). This is REAL, one-directional, ~5-6×, NOT the cornered dense-regen straddle. It only CASCADES to the
self-thin TPA under-kill on DENSE stands (where the height/BA lag delays the QMD-driven self-thin); at low density
(em_plant) TPA stays bit-exact and BA merely converges, which is why it read as "bit-exact" before — the .sum TPA
masked a real 5× TopHt height lag (META: measure TopHt, not just TPA/BA).
UNRESOLVED: my regent(NI/TT)/SMHTGF/htgf instrument-replay showed 0 growth for sub-BH — yet LIVE clearly grows them
1.18'→5-6'. So live's sub-BH established-tree height growth is in a path my dumps didn't capture (a regent branch I
mis-identified, or the DIAM/HITE sub-BH path, or an ESTAB per-cycle re-height). The FIX is to make jl grow the EM
sub-BH established cohort's height at live's rate. NEXT: instrument the LIVE per-tree HT trajectory of one planted DF
across 1990→2000 (dump HT(I) at each cycle for a tagged planted tree) to find the exact growing path, then port it.
★ This REFRAMES #137 from "self-thin over-kill (mortality)" to "planted/established sub-BH HEIGHT under-growth
(the self-thin is downstream)". The goal-file's "cyc0 mortality divergence" framing is superseded.

## Establishment height CONFIRMED equal (~1') — #137 is definitively a per-cycle GROWTH bug
Instrumented estb/estab.f:1023 (the PLANT/NATURAL tree-creation ESSUBH, called with TIME=FINT=10): live establishes
the planted DF cohort at HHT ≈ 1.00-1.12' (varies per-tree by BAA). jl's em_essubh_hht gives ~1.18' (validated). So
BOTH establish the cohort at ~1' — the establishment height is NOT the divergence. ⇒ #137 is DEFINITIVELY the
per-cycle GROWTH of the established sub-breast-height cohort: live grows them 1'→5-6' by 2000 (TopHt), jl leaves them
at ~1'. The growth path remains UNCAPTURED by my instrument-replay (regent NI/TT/SMHTGF, htgf all 0 for these trees) —
a genuine contradiction (every path I dump gives 0, yet live grows them). Likely I'm instrumenting the wrong branch,
OR the height jump is a re-assignment not an increment, OR esgent.f:53 REGENT(LESTB=.TRUE.) grows the establishment
cohort via a path/branch my dumps missed. DECISIVE NEXT: a CATCH-ALL per-tree HT trace — dump HT(I)+HTG(I) for the
established cohort at the END of base/grincr.f (after ALL growth routines) across 1990→2000, to see whether HT jumps
via HTG (increment — then trace which routine set HTG) or via a direct HT re-assignment (establishment re-height).
NET (this iteration): #137 reframed to a HEIGHT-GROWTH bug + establishment height RULED OUT (both ~1') ⇒ it's the
per-cycle growth of established sub-BH trees, ~5-6× under-grown, general (both densities), REAL not cornered.

## LIVE growth trajectory CAPTURED (base/update.f:65 HT=HT+HTG catch-all) — jl gives ~0, live 0.2-0.5'/cyc
Instrumented base/update.f:65 (HT(I)=HT(I)+HTG(I), the universal height application — catches HTG from ANY routine):
at ICYC=2 (first growth, →2000) live's small-tree cohort (HT<12') has HTG = 0.2-0.5'/tree and spans HT 1.3'→11.97'
(n=1857) — a real height SPREAD (⇒ TopHt=6 is the top-40 of that spread). Over cycles the spread + growth carries
trees across 4.5' into large-tree (htgf-POTHTG) growth. jl instead keeps the whole cohort flat at ~1' (TopHt=1) ⇒
jl's per-cycle EM small-tree HEIGHT growth returns ~0 for the established/planted sub-BH cohort where live returns
0.2-0.5'+. THAT is the #137 root: jl under-applies (≈0) the EM small-tree height increment to established sub-BH
trees. (The live HTG=0.2-0.5 with per-tree variation ⇒ the ZRAND-perturbed SMHTGF-style increment IS being applied
in live even at small size — so my earlier "SMHTGF returns 0 for dbh≤0" reading was the wrong branch; live's small
established trees DO get a nonzero height increment via update.f, source TBD but now MEASURED to be 0.2-0.5'/cyc.)
DECISIVE NEXT (jl-side, fast — Julia not gfortran): instrument jl's small_tree_growth!/_em_smhtgf for the em_plant_dense
cohort at cyc1 — dump the per-tree htgrth. It will show ~0; then trace why (the dbh≤0 guard, the KPER/REGYR scaling,
or the routing) and make jl apply the 0.2-0.5'/cyc increment live gives. This is the fix.
★ #137 fully MEASURED end-to-end: establishment ~1' (both) → live grows small trees 0.2-0.5'/cyc (spread to 12' by
2000, TopHt 6) → jl gives ~0 (flat at 1', TopHt 1) → jl BA lag → dense self-thin under-kill (+12.6% TPA). REAL bug,
jl-side height-increment under-application. Not cornered.

## jl-side MEASURED — jl COMPUTES healthy htgrth but it doesn't reach cyc-1 output (ordering, not the equation)
Instrumented jl small_tree_growth! _em_smhtgf on em_plant_dense (fast Julia): for the planted DF cohort jl computes
htgrth = 0.1-0.98'/subcycle with crown cr=82% (NOT 0 — RULES OUT any crown-init/VIGOR connection for #137) and h1
starting at 1.169' (=planting HHT). So jl's EM small-tree height INCREMENT model is working (healthy values). YET jl
reports TopHt=1' at 2000 (cyc1) — the computed increment does NOT reach the cycle-1 output. ⇒ #137 is NOT the growth
equation and NOT the crown; it's a CYCLE-1 establishment/growth/report ORDERING issue: jl's newly-PLANTED cohort is
established at cyc1 but its cyc1 height growth isn't applied/reported (likely established AFTER the growth pass, or the
wk3e→t.height application is skipped for freshly-established trees), so the cohort shows flat at ~1' at 2000 then grows
from 2010 (jl TopHt 2010=10 vs live 14) — a ~1-cycle lag. Live grows the planted cohort IN its establishment cycle
(update.f HTG=0.2-0.5' applied at cyc→2000, TopHt=6). DECISIVE NEXT: trace jl's cyc1 order — does establish_regen!/
the PLANT insertion run BEFORE or AFTER small_tree_growth! in the first cycle? If after, the planted trees miss cyc1
growth. FIX = establish planted trees before the cyc1 growth pass (or grow them in their establishment cycle).
★ #137 FULLY LOCALIZED (both sides measured): equation OK, crown OK (82%), establishment height OK (~1' both) ⇒ it's
the cyc1 establish-vs-grow ORDERING (jl planted cohort misses its first-cycle growth → ~1-cycle height lag → BA lag →
dense self-thin under-kill). REAL bug, not cornered. Reproduction: em_plant_dense.key / em_plant.key.

## ★★★ FIXED (commit c7c7d2f) — em_esgent! grows birth-cycle regen; EM was missing from the esgent list
ROOT (simulate.jl:566, confirmed): only CR (cr_esgent!) + TT (tt_esgent!) grew just-established regen in its birth
cycle; EM was omitted (comment literally: "eastern leaves them ungrown"). Added em_esgent! (mirrors tt_esgent! + EM
SMHTGF/SMDGF EMVAR branch, records nstart+1:n over the birth-cycle subperiod), wired after tt_esgent!.
VALIDATION vs FVSem_clean:
  em_plant_dense BA: jl 1,18,33,52,75,101,121,146,178,192  vs live 2,14,30,51,78,108,133,158,182,190
    (was the ~25% persistent-under 0,9,23,39,58,81,106,127,154,184). TopHt 2000: 1→6 = live.
  em_plant (400): TopHt 2000 1→4 (live 5).
  ZERO regression: emt01 bit-identical with/without fix (536-TPA read stand dominates); EM read-tree sweep stands
  (427473386 etc.) unchanged (em_esgent! is inert without PLANT/NATURAL — nstart==n → early return, no RNG draws).
RESIDUAL (refinement, not the root): em_plant_dense BA runs +10-29% in the EARLY cycles (2010 18 vs 14, 2020 33 vs 30)
then slightly under mid (2060 121 vs 133), converging at the ends. ⇒ the birth-cycle GENTIM/subperiod SCALING in
em_esgent! (I used gentim=fint-5, subyr/regyr=1 per tt_esgent!) isn't yet bit-exact vs live's em/esgent.f REGENT(.TRUE.)
period. NEXT (bit-exact refinement): instrument live em/esgent.f (or estb/esgent.f) for the birth-cycle KPER/subcycle +
XRHGRO/CON scaling on em_plant_dense cyc1, match em_esgent!'s subyr/con. The ROOT (missing birth-cycle growth) is
FIXED; this is a scaling-precision tail.

## #152 birth-cycle scaling — MEASURED discrepancy = TPCCF (birth-cycle density), not the period
Instrumented live em/regent.f:545 (EMVAR birth cycle, gated LESTB=.TRUE.) vs jl em_esgent! for DF (sp3) on
em_plant_dense: PERIOD matches (NTYR=fint−5=5 ⇒ NPER=1, KPER=5; XRHGRO=XRHMLT=1 default; CON=1 at establishment).
But the SMHTGF increment differs ~7-10×:
  live:  H1=1.00  HTGRR≈1.04  CON=1  H2≈2.04
  jl:    h=1.169  htgrth≈0.10-0.15  con=1  h2≈1.32
The cause is TPCCF: em_esgent! reads dens.point_ccf AFTER establish! recomputed density WITH the just-established
6000-TPA dense regen ⇒ point_ccf is huge ⇒ clamps to 300 ⇒ _em_smhtgf BETA1=exp(B0ACCF+B1ACCF·ln(300)) is suppressed
(htgrth~0.15). Live's birth-cycle DENSITY (the LESTB path, em/regent.f label 8, GO TO 8 skipping the DO-4 RDNEXT
build) EXCLUDES the just-established regen from its own competition ⇒ a LOW TPCCF ⇒ HTGRR~1.04. ⇒ #152 FIX: em_esgent!
must use a birth-cycle TPCCF that EXCLUDES the newly-established regen (the pre-establishment point_ccf, or the LESTB
RDNEXT that regent.f:label-8 builds), NOT the post-establish dense point_ccf. That will raise jl htgrth 0.15→~1.04 and
should close the em_plant_dense +10-29% early residual. (Note: cr_esgent!/tt_esgent! got away with the post-establish
point_ccf because their test stands lack a 6000-TPA dense-regen event; this is the SAME "dense regime exposes it" META.)

## #152 ROOT (corrected) — it's the CROWN, not TPCCF: em_esgent! runs before the regen crown is dubbed
Tested _em_smhtgf(3=DF, cr, tpccf, 0) directly: cr=0/tpccf=25 → 0.11 (=jl's em_esgent! value!); cr=82/tpccf=25 → 1.38;
cr=82/tpccf=100 → 0.90 (live HTGRR≈1.04 sits at cr~82, tpccf~70). So jl's em_esgent! htgrth=0.1-0.15 is the cr=0 curve.
ROOT: establish! adds the regen with ICR=0 (simulate.jl:562 "adds regen (ICR=0)"), and the crown is dubbed by
crown_ratio_update! at simulate.jl:573 — AFTER em_esgent! at :569. So em_esgent! reads crown_pct=0 ⇒ _em_smhtgf HTG1
= BETA1 + BETA2·0 = BETA1 only ⇒ suppressed (0.1). Live's establishment (estb/estab.f) assigns the regen crown BEFORE
esgent grows it (cr~82 for planted DF) ⇒ HTGRR~1.04. This is a SHARED esgent ordering issue (cr_esgent!/tt_esgent!
also read t.crown_pct before the :573 dub) — masked for CR/TT because their test stands lack a dense-regen event where
cr matters, and because the crown TERM dominates only for the EMVAR SMHTGF form. #152 FIX (direction): dub the
just-established regen's crown BEFORE the birth-cycle esgent (mirror live estab.f's ICR assignment, or run the crown
dub for the new records before em_esgent!), so em_esgent! sees cr~82 not 0. Needs the establishment-crown value/order
matched to live (instrument estab.f ICR for the planted cohort). ⇒ #152 is a CROWN-BEFORE-ESGENT ordering fix, not a
TPCCF fix (supersedes the prior TPCCF hypothesis — corrected by the direct _em_smhtgf test).

## #152 root — em_esgent! sees crown=0; the fix is crown-before-esgent (needs path audit)
Confirmed via _em_smhtgf test: jl em_esgent! htgrth=0.1 IS the cr=0 curve (cr=82 → 0.9-1.4). So for em_plant_dense's
PLANT regen, crown_pct=0 when em_esgent! runs. jl DOES have regen-crown dubs — establish! phase-2 (establishment.jl:379
cr=clamp(0.89722−0.0000461·pccf+0.07985·N(0,1),0.20,0.90)→~85) AND ie_autoes_establish! (:1173) — but for THIS PLANT
path the crown isn't set on the records em_esgent! iterates (nstart+1:n) before em_esgent! runs. FIX (scoped): ensure
the just-established regen crown is dubbed BEFORE em_esgent! for the EM PLANT/AUTOES path — either (a) verify establish!
phase-2 covers the EM PLANT regen and runs before em_esgent!, or (b) dub the crown inside em_esgent! when crown_pct==0
using the phase-2 formula (with the matching MAIN-RANN draw order). Then htgrth 0.1→~1.0 and the em_plant_dense +10-29%
early residual should close. NOTE: cr/tt_esgent! read t.crown_pct too — audit whether their regen crown is set first
(may share this latent ordering, masked without a dense-regen + crown-sensitive SMHTGF test).
★ #137 SUMMARY: ROOT (EM missing from birth-cycle esgent list) FIXED (c7c7d2f) — em_plant_dense BA tracks live (was
~25% under), 0 regression. #152 residual (birth-cycle crown=0 → SMHTGF suppressed ~7×) precisely root-caused; it's a
crown-before-esgent ordering fix. Both are the SAME theme as #149/#150 (regen/read-seedling crown feeding small-tree
growth) — the crown must be dubbed before the growth reads it.

## ★ CORRECTION (2026-08-06) — em_esgent! IS faithful; #152 residual is the accepted straddle, NOT crown/TPCCF
Both prior #152 root-cause commits (657c44a TPCCF, 79e910f crown) were WRONG — artifacts of a `sort -u | head` that
cherry-picked the LOWEST htgrth (0.1) values. DIAGNOSTIC (instrumented phase-2 crown-set + em_esgent! crown-read on
em_plant_dense): em_esgent! reads the CORRECT crown — cr=82/90 (phase-2 establishment.jl:379 sets it, runs before
em_esgent!). And _em_smhtgf(3, cr=82/90, tpccf=100-300) = 0.90-0.98 ≈ live HTGRR 1.04. So em_esgent!'s birth-cycle
increment MATCHES live; the crown is NOT 0 and NOT the issue.
The em_plant_dense BA residual is MIXED-SIGN and converging: diffs jl−live = −1,+4,+3,+1,−3,−7,−12,−12,−4,+2 across
2000-2090 (not one-directional) — the accepted dense-cohort DGSCOR/AVHT40 tie-break + ZRAND-realization STRADDLE, the
SAME cornered primitive as BM #149 (+9%), CI #142 (~2%), the whole-cluster dense-regen straddle. ⇒ #137 is now
bit-exact-OR-CORNERED: the em_esgent! fix (c7c7d2f) resolved the REAL bug (the ~25% one-directional under-growth from
the missing birth-cycle growth); the ±10% mixed-sign remainder is the accepted straddle. #152 is CORNERED, not a real
refinement — em_esgent! is faithful. ★ LESSON (self-inflicted): NEVER characterize a distribution from `sort|head` —
it cherry-picks the tail; sample representatively (mean/median or the actual per-tree with its inputs). Two wrong
root-causes committed then refuted by the direct diagnostic — doctrine #2 (measure) caught it, but only after I
inferred from a biased sample. Measure the RIGHT thing.

## #137 fix — AUTOES-path no-regression CONFIRMED (+ surfaced #143 residual)
Checked em_esgent! against EM AUTOES (auto-established, empty-inventory) stands 103399881010661/103406000010661/
11849257010690 (em_sub.db): jl cyc1 (2016) TPA = 246/252/212 WITH the fix, IDENTICAL 246/252/212 WITHOUT (stashed).
⇒ em_esgent! does NOT regress the AUTOES tally amount (it grows the established regen's height/DBH, not the TPA count;
these are sub-BH BA=0 so the growth is .sum-inert on TPA/BA here). Live = 219/224/189, so jl over-establishes ~+12% —
that is the PRE-EXISTING #143 AUTOES tally-amount residual (ESTOCK/ESPROB tally, a separate subsystem from em_esgent!),
NOT caused by the #137 fix. The earlier "#143 AUTOES amount bit-exact" was on the 4 NULL-elev stands only; these
elev-present stands show the +12% diffuse residual. ⇒ #137 fully validated (0 regression incl. AUTOES); #143 amount
residual remains open (unchanged).

## #143 AUTOES amount +12% — LOCALIZED past the ESTOCK PN (bit-exact) to the downstream tally
On elev-present EM AUTOES stand 103399881010661 (jl cyc1 246 vs live 219, +12%): instrumented live estb/estock.f PN
and jl ie_estock — BOTH give ifo=11, elev=16.8, PN=−0.44419 (BIT-EXACT). So the +12% is NOT the stocking logit / elev
term. ⇒ the residual is DOWNSTREAM in the tally chain: FTEMP = 1/(1+exp(−(PN+ESB−ESB1)))·STOADJ → PROB1 → ITPP =
INT((PLPROB·DUPNPT)/(FTEMP·300)+0.5) [ingrowth] or the TPP draw [disturbance] (estab.f:541/589/679). Candidates for
the +12%: (a) ESB−ESB1 inventory calibration (ESB from inventory, ESB1=ie_estock at inventory — the cyc1/2 correction),
(b) STOADJ, (c) the per-plot ITPP/PLPROB amount or MAXTPP/MAXING cap. NEXT: instrument live estab.f FTEMP/PROB1/ITPP/
NSTORE for this stand vs jl ie_autoes_establish! — the diverging tally factor is the +12% root. (PN ruled out by
measurement.) This is the IE-AUTOES-shared tally (#143); the elev-default fix 29d449d closed the NULL-elev PN path,
this is the amount-draw path.

## #143 +12% — TIGHTLY LOCALIZED to the ITPP/plot-tally amount (FTEMP stocking prob is bit-exact)
Instrumented live estab.f:585 (PROB1=FTEMP) on stand 103399881: PN=−0.444186, ESB=0, ESB1=0, STOADJ=1.0 ⇒
FTEMP=1/(1+exp(0.444186))·1 = 0.390744. jl's PN is bit-exact (−0.4441864) and ESB/ESB1=0, STOADJ=1 for an empty
stand ⇒ jl FTEMP=0.3907 too. So the STOCKING PROBABILITY is bit-exact; the +12% (jl 246 vs live 219) is PURELY in
the amount draw: ITPP (trees-per-plot) × plot count → TPA. The ESTPP disturbance path (estab.f:679, INT(TPP+0.5)) did
NOT fire for this empty stand (my TALLY dump 0 hits), so the amount comes from a different tally path (the "predicted
naturals" / per-plot PLPROB path). NEXT: instrument the ACTUAL amount path for this empty stand — dump ITPP/NEWTPP/
NSTORE/the plot count + final TPACRE per species — vs jl ie_autoes_establish!'s ITPP draw. The +12% is one of: the
ITPP-per-plot draw (ESTPP/ESTPP-equivalent), the plot count (DUPNPT/NPTIDS), or the TPA scaling. FTEMP/PN/ESB/STOADJ
all RULED OUT bit-exact. This is a tight, well-scoped close-out for the IE-shared AUTOES tally (#143).

## #143 +12% — ROOT MEASURED: the INGROWTH-path ITPP raw-draw skew (cap correct, RNG/formula)
Instrumented live estab.f for stand 103399881: INGRO=1, MATCH=0, NTALLY=1, IHAB=3, MAXTPP=5. So this empty
auto-establishing stand is INGROWTH (INGRO=1) in live — jl's is_ingro path (cap=_IE_MAXING[3]=3) is CORRECT (with
INGRO=1, estab.f:682 makes MAXING=3 the binding cap over MAXTPP=5). Both cap at 3. Yet the ITPP DISTRIBUTION differs:
  live ITPP: 22×1 / 15×2 / 15×3  = 97 tree-slots
  jl   ITPP: 18×1 / 11×2 / 23×3  = 109 tree-slots  = +12.4% (= the exact amount residual)
jl's raw ITPP draws skew HIGHER (more 3s, incl. draws that would be 4/5 clamped to MAXING=3). Since the cap + INGRO
classification match, the +12% is the INGROWTH-path ITPP RAW DRAW: jl's ie_estpp(ie_esrann()) gives higher trees-per-
plot than live's ESTPP on the ingrowth path. This is the KNOWN open residual flagged at establishment.jl:676 ("the
ingrowth (is_ingro) path has a separate open residual"). NEXT: instrument live ESTPP DRAW/TPP vs jl ie_estpp draw on
the ingrowth path (the ie_esrann seed-chain position, or the ie_estpp VAL/BB/CC formula) — the ingrowth-path RNG
advance or the ESPROB-driven NEWTPP normalization (estab.f:589 ITPP=INT((PLPROB·DUPNPT)/(FTEMP·300)+0.5)) is the
skew source. ★ #143 converted from "diffuse ~22%" to "ingrowth ITPP raw-draw +12%=109-vs-97 tree-slots; PN/FTEMP/ESB/
STOADJ/cap/INGRO all bit-exact-or-correct". Same residual class as the IE ingrowth audit item.

## #143 +12% DEFINITIVE ROOT — ingrowth per-plot RNG seed-chain DESYNC (not the formula)
Instrumented live estab.f ESTPP (DRAW,TPP) vs jl ie_estpp(val,tpp) on stand 103399881:
  live DRAW: 0.251→1.21, 0.360→1.48, 0.153→1.04, 0.276→1.26, 0.748→3.92, 0.571→2.38, 0.131→1.01, 0.311→1.35
  jl   val : 0.346→1.44, 0.351→1.45, 0.212→1.13, 0.216→1.14, 0.992→19.6, 0.932→8.9, 0.415→1.66, 0.122→1.00
The DRAW/val values DIFFER ⇒ the ie_estpp FORMULA is fine (val→tpp mapping matches: e.g. jl 0.212→1.13 ≈ live
0.153→1.04 shape); jl consumes the WRONG ie_esrann DRAW at the ESTPP point. jl even draws 0.99/0.93 (→TPP 19.6/8.9,
capped to MAXING=3) that live has no counterpart for ⇒ jl's ITPP skews high (+12%). ⇒ #143 root = the EM INGROWTH
per-plot RNG SEED-CHAIN is misaligned: the per-plot body advance (16+3·nsp+2·MAXTPP[ihab]=83 for EM ihab3) or the
WITHIN-plot ie_esrann order (the DO-71..75 ESRANN loops at estab.f:655-669 that jl must replicate exactly, or the
ESB/species draws) is off by some draws, so the ESTPP draw lands at the wrong RNG position. FIX: align jl
ie_autoes_establish!'s ingrowth per-plot ie_esrann sequence to live estab.f exactly (count + order) so the ESTPP
DRAW matches. This is the known ingrowth residual (establishment.jl:676). ★ #143 fully root-caused: PN/FTEMP/ESB/
STOADJ/cap/INGRO/ie_estpp-formula ALL correct; the sole residual is the ingrowth RNG draw ALIGNMENT.

## #143 — RNG offset QUANTIFIED: jl under-advances ~147 ie_esrann draws before the first ingrowth ESTPP
Dumped the full live ESRANN stream (estb/esrann.f, 4369 draws for the 3-cycle run). Live's first ESTPP DRAW
(0.251398891) is at stream position ~205; jl's first ESTPP val (0.346301585) EXACTLY matches live's draw at
position ~58 (jl's ie_esrann LCG stream = live's, so the values align — the desync is purely a COUNT offset). ⇒ jl
consumes its first ingrowth ESTPP draw ~147 draws EARLIER than live: jl's establishment SETUP (before the plot loop
reaches ESTPP) draws ~147 FEWER ie_esrann than live. Live's pre-ESTPP setup (ESTIME, EMSQR ±draws, site-prep, per-
species ESB/ESPROB draws, the DO-71..75 ESRANN loops at estab.f:655-669, the per-plot body 16+3·nsp+2·MAXTPP) has
~147 more draws than jl replicates for the EM ingrowth path. FIX: add the missing setup/per-plot ie_esrann advances
to jl ie_autoes_establish!'s ingrowth path so the ESTPP draw lands at live's position (then the ITPP dist + the +12%
close). NEXT: instrument live ESRANN with a per-CALL-SITE tag (or bisect the 147) to identify WHICH setup phase jl
skips — likely a per-plot or per-species ESRANN loop count. ★ #143 root fully quantified: RNG count offset ~147
draws on the EM ingrowth path; all equations/PN/FTEMP/cap correct.

## ★ CORRECTION (2026-08-06) — #143 is NOT a count offset; the ie_esrann SEQUENCES diverge EM-specifically
The prior "pure ~147-draw count offset, streams align" (commit 63685cd) is REFUTED. Searched all 6 of jl's ingrowth
ESTPP vals in the full live ESRANN stream (4369 draws): only 0.346301 appears (position 54) — COINCIDENTAL; the other
5 (0.350846/0.212444/0.215863/0.991703/0.932161) have ZERO occurrences. If it were a count offset in the SAME LCG
stream, ALL of jl's vals would appear in live's stream at earlier positions. They don't ⇒ jl's ie_esrann is at a
DIFFERENT LCG STATE than any live's ESRANN passes through for this stand ⇒ NOT a fixable count offset. Since jl's
ie_estpp/ie_esrann is VALIDATED bit-exact on iet01 stand-4 (memory), the RNG matches THERE but diverges on the EM
stand ⇒ EM-INGROWTH-SPECIFIC seed/state divergence: jl's per-plot ie_esrann seeding (or a count divergence from an
EARLIER point in the EM establishment setup that puts jl in a different LCG state by ESTPP) differs from live. ⇒ the
+12% root is a per-plot RNG STATE mismatch on the EM ingrowth path, deeper than a draw-count adjustment. NEXT: compare
jl's ie_esrann per-plot SEED derivation vs live's ESRANN ESS1 chaining for the EM stand (is jl re-seeding per plot
where live continues a global chain? is the initial ESS1 seed for the stand different?) — the seeding/chaining is the
divergence. ★ SELF-CORRECTION LOGGED (doctrine #2, twice this session incl #152): the 9-digit coincidental match
misled me into a false "offset"; the FULL-sequence search refuted it. Verify the WHOLE distribution, not one point.

## #143 — ACCURATE ROOT (reconciles all prior findings): constant per-plot body_n=83 DRIFTS over many plots
Read jl ie_autoes_establish! (establishment.jl:687-696): jl derives a per-plot ESTPP seed by advancing a CONSTANT
body_n = 16 + 3·nsp + 2·MAXTPP[ihab] (=83 for EM nsp19/ihab3/MAXTPP5) between plots (ie_autoes_plot_seeds, line 690),
then RESEEDS a fresh IEEstabRNG(sd) per plot (line 693) and draws wk6fill + 2 EMSQR + 1 ESTPP. This reconciles both
prior (partial) findings: (a) the ESTPP draws aren't in live's MAIN ESRANN stream because each plot's ESTPP is a
RESEEDED sub-stream IEEstabRNG(sd); (b) the "count offset" intuition was right in spirit — jl's per-plot SEED (sd)
diverges from live's because the ADVANCE between plots is wrong. The body_n=83 was DERIVED + validated on a 6-plot
sequence (jl [1,2,1,1,3,3] = live bit-exact, comment line 686). But stand 103399881 has 52 PLOTS ⇒ a CONSTANT 83
per-plot advance DRIFTS from live's ACTUAL per-plot draw count (which VARIES by plot outcome: the number of ESRANN
draws a plot consumes depends on its ITPP/species/STOADJ path — e.g. the DO-71..75 loops at estab.f:655-669 only fire
when STOADJ<0.0001, and the species/height draws scale with the plot's realized establishment). Over 52 plots the
small per-plot error accumulates ⇒ jl ITPP 109 vs live 97 tree-slots = +12%. ⇒ #143 FIX: replicate live's VARIABLE
per-plot ESRANN advance (not a constant 83) — count the ACTUAL draws each plot consumes in live estab.f (condition on
STOADJ/ITPP/NOFSPE per plot) and mirror it in jl's seed-chain. This is why NULL-elev stands (bit-exact) and the 6-plot
validation passed but the 52-plot elev-present stand drifts. Known residual (establishment.jl:676). ★ This SUPERSEDES
the "count offset" (63685cd) and "RNG state mismatch, not count" (eacf34a) — BOTH partial; the accurate root is the
CONSTANT-vs-VARIABLE per-plot seed-chain advance drifting over many plots.

## ★★ #143 — SOLID MEASURED FACTS (3rd correction; stop speculating the mechanism, these are measured)
Found live's ESTPP DRAW values (0.251398/0.359962/0.153236/0.275537/0.747601/0.570895/0.131472/0.311342) in the full
live ESRANN main-stream dump at positions: 56, 139, 222, 305, 388, 471, 554, 637 — GAPS = 83,83,83,83,83,83,83 EXACTLY
CONSTANT. ⇒ live's per-plot advance IS a constant 83 (= jl's body_n) — my "variable drift" root (76e8cf3) is WRONG too.
MEASURED (solid, not interpreted):
 1. Live ESTPP draws are the CONTINUOUS MAIN ESRANN stream, one per plot at a constant 83-draw interval (from position 56).
 2. jl's 1st ESTPP val (0.346301) = live main-stream position 54 — i.e. 2 draws EARLY vs live's 56.
 3. jl's 2nd–6th ESTPP vals are ABSENT from the live main stream entirely.
INTERPRETATION (tentative — I've been wrong 3×, treat as hypothesis): jl uses per-plot RESEEDED sub-streams
(ie_autoes_plot_seeds → IEEstabRNG(sd), establishment.jl:690-693) whereas LIVE uses ONE CONTINUOUS main stream for the
ESTPP draw (estab.f XSTORE=ESRANN draw). jl's reseed lands ~2 early on plot 1 and then diverges off-stream for plots
2+. The likely fix is to make jl's ESTPP DRAW the CONTINUOUS-stream draw (advance the ONE stream by 83 per plot and
take the draw), NOT a reseeded IEEstabRNG(sd) sub-stream — OR fix the intra-plot draw offset (2) + ensure the reseed
seed is the EXACT ess0 (not a truncated ESDRAW). ★★ META (important): I have committed THREE different #143 root-causes
this session (count-offset 63685cd / state-mismatch eacf34a / variable-drift 76e8cf3), each refuted by the next
measurement. This is a signal that #143's RNG structure (per-plot reseed vs continuous stream) is subtle and I'm
over-drilling at length — it needs FRESH, careful analysis of ie_autoes_plot_seeds vs estab.f's XSTORE/continuous
stream, comparing the EXACT draw indices, NOT another quick hypothesis. The MEASURED positions (56/139/.../637 gap-83;
jl 1st at 54; jl 2nd+ off-stream) are the reliable handoff.
