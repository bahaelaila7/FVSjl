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
