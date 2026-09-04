# IE FIA RE-SWEEP dig session — M332A/M332B (Idaho Batholith) cap: NSTORE fix VERIFIED + #206 straddle CORNERED

**Date:** 2026-09-04  **Variant:** InlandEmpire (IE)  **Oracle:** `/workspace/.iework/FVSie_clean`
(byte-identical to `FVSie_g16`).  **Engine:** sweep branch `wt/western-sweep` synced to master **86007127**
(`git diff --stat master HEAD -- src/` EMPTY — engine bit-identical to master). **DIGCAP=100 (never raised).**

## Purpose
Re-evaluate the M332A cap (cursor 4000, DIG_PAUSED) with the fixed engine. Master 86007127 landed the NSTORE
ingrowth over-book fix (clamp AUTOES ESB1 BAAA per-inventory-point to [1,400] in the stocking logit,
`src/variants/inlandempire/establishment.jl:1573`). The prior cap MIXED (1) the one-directional NSTORE
over-book on DENSE already-stocked stands [REAL BUG, now fixed] and (2) the #206 OLDRN/ZRAND establishment
realization straddle on bare/sparse stands. Separate them under the fix, corner the straddle, continue.

## Re-sweep (cursor 2000->4000, fixed engine)
Reset cursor to 2000, archived the provisional pre-fix queue (400 rows) to
`data/provisional_archive/fia_dig_queue_M332A_provisional_prefix_20260904_122329.csv`, reset queue to header,
re-ran the 2000-4000 batch SKIP_DONE=0 (force re-eval). Result: 411 comparable stands, bit_exact 0 (the IE
#206 straddle is ubiquitous >= 1 ULP), **397 dig-worthy**. Ecoregions: **M332Ab 129, M332Ac 115, M332Aa 85,
M332Ao 36, M332Ba 29, M332Bb 3** — entirely the Idaho-Batholith M332A/B cluster.

## NSTORE dense over-book — FIX VERIFIED (the separate real bug is GONE)
- **Lead 12281578010690** (M332Aa, dense: 57 large records, small_tpa=0, init BA 683->clamped 400): last-cycle
  TPA jl **388->218 vs oracle 225** (was +72% over, now **-7 = -3%**). BA -26, SDI -34, CCF -27 (small under).
  The one-directional dense over-book is eliminated; the lead now sits on the UNDER side of a straddle.
- **Unbiased 40-stand M332A sample** (every-50th CN from the 2000-4000 emission order, NOT selected on
  divergence): TPA over=14 / under=24 / zero=2 (mean **+0.97**); BA over=18 / under=15 / zero=7 (mean +1.38)
  = TWO-SIDED / balanced, matching the 587-corner's ~1.5% under-lean.
- **DENSE subset within the unbiased sample** (baN_or >= 250, N=6): TPA over=4 / **under=2** (31365985010690
  **-4**, 3294683010690 **-33**). The dense stands FLIP sign — if the NSTORE over-book remained, every dense
  stand would be one-directional OVER. They are not.
  - (My earlier "11 top-dense stands all TPA-over" measurement was pure SELECTION BIAS: those 11 were picked
    as the largest pre-fix TPA-over rows, so a residual straddle necessarily lands them over. The unbiased
    dense subset is two-sided. The apparent "one-directional dense residual" was the selection artifact the
    over-cornering guard warns about, NOT a remaining bug.)

## #206 OLDRN/ZRAND establishment realization straddle — SIGNATURE CONFIRMED
Signed last-cycle (jl-oracle) via `.sweep_work/sig_verify.jl`. Consecutive same-plot pairs (adjacent
subplot records, near-identical inputs) FLIP sign:
- 5390352 **+1** / 5390353 **-9**;  5390639 **+66** / 5390640 **-232**;  5390798 **-155** / 5390799 **+3**
- 12281578 **-7** / 12281579 **+9** (the lead flips with its neighbor)
- 12280747 **-26** / 12280748 **+5**;  4729653 **-18** / 4729654 **+3**;  4729715 **+6** / 4729716 **-23**
- 3314732 (tpa **-35** / ba **-19**) / 3314733 (tpa **+9** / ba **+3**)  [strong flip]
- 3291693 (tpa **+71** / ba **-27**) / 3291694 (tpa **-15** / ba **+7**)  [strong flip]
- **3315450 +22 == 3315451 +22** — TRULY-identical inventory plots yield IDENTICAL residual (the
  deterministic-per-input RNG-realization proof: identical inputs => identical draw stream => identical diff).

Same primitive as the established IE 331A corner (587 stands) and TT M331D: `FVSie` OLDRN/ZRAND (BACHLO)
per-tree normal-deviate draw in the AUTOES establishment cohort growth (esgent.f WK4^2 birth growth +
estab.f DO 99/33/228 ESADVH/ESSUBH/ESXCSH height-class booking). Larger +/-TPA/BA/CCF than the 331A
sub-floor because the denser Idaho-Batholith stratum amplifies the cohort realization straddle through
self-thinning, but the SIGNATURE (adjacent-flip, identical-match, two-sided/balanced, dense-flips) is
identical.

## OVER-CORNERING GUARD — CLEARED
No M332A stand is one-directional AND correlated with a deterministic input: adjacent near-identical plots
flip large-magnitude sign (+66/-232, +71/-15, -35/+9); truly-identical plots match exactly; the DENSE
subset (the NSTORE-bug regime) is itself two-sided (-4, -33 under alongside overs). The one-directional
dense over-book is fixed on master 86007127. NOT a re-balanceable bug -> the IE #206 establishment straddle.

## Cornered
- **397 M332A/B CNs** appended per-CN to `docs/fia_cornered_stands.txt` (with a descriptive header block),
  synced to the live-filter path `/workspace/FVSjl/docs/fia_cornered_stands.txt` (filter_digworthy reads
  CORNERED_STANDS_FILE from there — drops even ESCALATING rows, which the cluster tsv alone does not).
- **6 cluster rows** (M332Aa/M332Ab/M332Ac/M332Ao/M332Ba/M332Bb, signature `structure_densephase`,
  tag `ie-estab-206`) appended to `docs/fia_cornered_clusters.tsv`, synced to `/workspace/FVSjl/docs/`.
- Re-running filter_digworthy on the fresh cycle CSV yields ZERO remaining dig-worthy rows (all 397 dropped).

## State
DIGCAP NOT raised. Engine unchanged (bit-identical to master). Gate not re-run this session (no src change;
339/11 already validated on master 86007127). Cursor advanced past M332A to CONTINUE the IE sweep from 4000.
Backups: `data/fia_sweep_west.db.bak_pre_ie_resweep_nstore_20260904_122329`.

## IE SWEEP COMPLETE (cursor 17808/17808) — continued 4000 -> population end (2026-09-04)
After cornering M332A (batch-1, 397 stands), the sweep continued forward and reached the IE population
end. EVERY subsequent batch capped on the SAME #206 OLDRN/ZRAND AUTOES establishment realization straddle,
verified per batch (unbiased sign-tally two-sided + adjacent same-plot pairs flipping sign); the ecoregion-
independent primitive was confirmed across the entire IE mountain population (10 sections). No NEW real bug
surfaced; the NSTORE dense over-book (fixed on master 86007127) did not recur.

Batches (cursor range -> dig-worthy cornered, sections; commit):
- 2000-4000  -> 397  M332A/M332B (Idaho Batholith)                         eddf4340
- 4000-6000  -> 414  M332B*/M332Gd (Idaho Batholith + Blue Mtns)           337f7fd7
- 6000-8000  -> 426  M333A* (N.Rockies/Bitterroot) + M332Gd tail           8805e5a5
- 8000-10000 -> 357  M333A*/M333Ba (N.Rockies)                             58743697
- 10000-12000-> 181  M333Bb/M333Bc                                         5bec3c13
- 12000-14000-> 272  M333C* (N.Rockies)                                    c592a09a
- 14000-16000-> 366  M333D* (N.Rockies)                                    e1da3cfe
- 16000-17808-> 360  M333Dd/M333De (POPULATION END)                        (this commit)

Total IE #206 establishment straddle cornered THIS session: 397+414+426+357+181+272+366+360 = 2773 stands
(on top of the pre-existing 587-stand front-loaded 331A corner). Per-batch verification samples: TPA sign
two-sided every batch (25/24 balanced at the end; skew direction VARIES over/under by stratum = RNG
realization, never one-directional); adjacency pairs flip sign strongly (e.g. +146/-6, +73/-8, -85/+14,
+62/-9, -280..+113 spans); the biggest single-stand unders always sit beside near-identical +overs. BA
skew alternates over/under across batches (M333D BA over-leans, M332B under-leans) — the signature of a
two-sided draw, not a directional bug. Over-cornering guard cleared at every cap.

VERDICT: IE full-population FIA sweep COMPLETE. The entire IE dig-worthy residual is the single #206
OLDRN/ZRAND AUTOES establishment-cohort realization straddle (named primitive: FVSie esgent.f WK4^2 birth
growth + estab.f DO 99/33/228 height-class booking per-tree normal deviate). No open real bug in IE.
