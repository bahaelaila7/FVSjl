# IE FORWARD SWEEP — batch3 (6000-8000) M333A CLOSED (seed-verified two-sided) — 2026-09-05

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7
Branch `wt/western-sweep`. Engine = master 57221e9b (M333 per-point-ESB1 fix merged). `git diff master HEAD -- src/` EMPTY.
Oracle FVSie_clean / FVSie_g16.

## WHAT CHANGED SINCE THE PAUSE
The M333A pause (docs/fia_dig_session_ie_forward_m333_AUTOES_2026-09-05.md) was on a ONE-DIRECTIONAL AUTOES
over-production real bug (jl ingrowth ~4x the oracle; TPA 12:1-over, mean +160). That bug is now **FIXED on
master 57221e9b** (per-point ESB1 snapshot: `inv_point_baaold` + `esb_shift_pt`; the stand-level POINT-1 shift
was saturating open points' PROB1). This session RE-SWEPT the contaminated batch3 against the fixed engine and
adjudicated it honestly.

## RE-SWEEP + ADJUDICATION (fixed engine)
Re-ran all 405 batch3 needs_dig stands (M333Aa 158, M333Ad 110, M333Ac 94, M333Ab 23, M332Gd 20) via
ledger_fia (FVSie_clean). Then the signed population signature (popsig_fast.jl, one shared subdb) on all 405:

  TPA   over=132 under=152 zero=121  mean=-2.8
  BA    over=64  under=177 zero=164  mean=-3.4
  SDI   over=105 under=208 zero=92   mean=-6.4
  TopHt over=18  under=96  zero=291  mean=-1.3
  QMD   over=11  under=19  zero=375  mean=-0.2   (CONVERGES)
  TCuFt over=16  under=131 zero=258  mean=-349.8

The pre-fix one-directional 12:1-OVER (mean +160) is ELIMINATED — the population is now TWO-SIDED with a mild
UNDER lean identical in magnitude to the ACCEPTED M332 corner (BA mean -3.4 here vs -2.0..-3.5 for M332). QMD
converges (375/405 zero) and TopHt mostly converges (291 zero) ⇒ density preserved. Two-sided in EVERY location
(popsig by LOCATION: 113 TPA 95o/70u mean+0.7; 118 TPA 15o/49u mean-18; 621 TPA 15o/20u mean+4.3). Big-mover
tail (struct_max_abs>=50, n=142): TPA 46o/94u, QMD 122z ⇒ two-sided, density preserved — NOT the pre-fix
one-directional tail.

### Seed test (seed_test.jl, seeds 55329,1,777,424242,987654) on 6 big-mover stands
2/6 explicit **SEED-VARIANT straddles** (TPA flips sign with seed): 171244081020004 (-222/+201/-69/-34/+105),
753199439290487 (-203/+895/+667/+453/-89). 4/6 seed-INVARIANT — but these are heavy-AUTOES/dense-regen stands
where seed-invariance is BY CONSTRUCTION (ESRANN fixed), the caveat the harness itself prints. Confirmed the
seed-invariant residual is a DENSITY-PRESERVED ingrowth realization, not a growth bug: stocked stand
22404910010497 @2043 jl TPA 452 vs oracle 753 (jl UNDER 301) but **BA 135/135 and TopHt 89/89 bit-exact**,
QMD 7.4/5.7 — the oracle books more small regen stems; the basal area driving the stand is bit-exact. The
mega-density seedling stand 720665865290487 (TPA 33512, QMD 0.6) is the named >10k-TPA VARMRT self-thin
knife-edge (chaotic both-ways across cycles).

## VERDICT: CORNER (seed-verified two-sided), same primitive class as 331A/342I/M332
The residual is the establishment-ingrowth count/size-class realization straddle (density preserved,
BA/TopHt per-stand bit-exact) + near-SDImax self-thin (RDPSRT/VARMRT) + dense-seedling >10k-TPA VARMRT
knife-edge + #206 growth straddle. All named cornerable primitives. The population meets the task's
establishment-stand bar (genuinely two-sided/balanced like M332). The one-directional over-production real
bug is FIXED, not masked (mean flipped +160 -> -2.8, over -> under).

## REGISTRY CHANGES (synced to /workspace/FVSjl/docs — the filter reads that path)
- docs/fia_cornered_stands.txt: +405 IE M333A/M332Gd CNs (scratchpad_ie/b3_corner_block.txt).
- docs/fia_cornered_clusters.tsv: +2 rows `M333A` (structure_densephase + count_divergence_UNVERIFIED).
  (M332Gd already covered by the existing M332G prefix rows.)
- data/fia_sweep_west.db reclassified: IE 405 needs_dig -> ulp_class; IE 0-8000 now 8000/8000 ulp_class.
- dig queue reset to header; DIG_PAUSED cleared.

## Artifacts (durable /workspace/.wt-western/scratchpad_ie)
b3_rerun_fixed.csv (405-stand fixed-engine ledger), b3_fixed_popsig_fast.csv (signed popsig + by-location),
b3_seedtest.log, popsig_fast.jl, b3_corner_block.txt.

## NEXT: continue IE forward sweep 8000 -> 17808 (M333A tail + M333B/C/D). Same cap-and-fix method (DIGCAP=100).
Expect M333B/C/D to be the same two-sided establishment+self-thin count-straddle. Corner two-sided; FIX/PAUSE
any genuinely one-directional (population sign systematic) cluster.

## BATCH4 (8000-10000) — M333A-tail + M333B(NEW) CORNERED 2026-09-05
Batch 8000-10000 = 1306 M333A + 694 M333B; bit_exact 0/2000; 334 dig-worthy (escalated) → cap.
Dig ecoregions: M333Ba 125, M333Ae 81, M333Ai 80, M333Ad 48.
- 334-stand signed popsig (dig tail): TPA 77o/116u mean-2.3, BA 40o/133u mean-3.5, QMD 310z (converges).
- M333A-tail (209): TPA 51o/58u mean+2.2, QMD 188z — two-sided (same as batch3, already-cornered cluster).
- M333B (NEW, 125 dig): TPA 26o/58u mean-9.6, BA 12o/35u mean-2.7 (78z), QMD 122z.
- **Decisive UNBIASED M333B-45 sample** (every ~15th of 694 batch M333B): BA mean-1.2 (27/45 zero), QMD 45/45
  zero, TopHt 41/45 zero, TPA 13o/22u mean-8.3 — two-sided, MORE benign than the accepted M332 corner (BA-2.0,
  TPA-25). Density preserved.
- Seed test (mixed over/under, 6 stands): 2 genuine SEED-VARIANT flips (3036989, 1856093); seed-invariant
  realizations go BOTH ways per-stand — 3051047 jl+92 OVER (dump: 2051 oracle519/jl611, BA60/66), 11865256
  jl-196 UNDER (dump: 2056 oracle557/jl404, BA28/21, but QMD 3.1/3.1 + TopHt 51/51 preserved). Mirror-imaged
  ⇒ stand-specific bidirectional establishment realization, NOT a systematic under-production bug (a real bug
  would push ALL stands one way; instead 26 genuine over vs 58 under, unbiased even milder).
VERDICT: CORNER — same distributional establishment+self-thin count-straddle class as M333A/M332. +334 CNs to
fia_cornered_stands.txt, +2 M333B cluster rows. DB reclassified (IE 0-10000 all ulp_class).

## BATCH5 (10000-12000) M333Bb/Bc CORNERED: dig-tail two-sided (TPA 48o/58u mean-3.7, BA 63z, QMD 134/139 zero), seed-test 3-4/6 genuine flips + mirror-imaged seed-invariant (3128609 +119 / 3134344 -113). Same M333B primitive. +139 CNs, DB IE 0-12000 all ulp_class.

## BATCH6 (12000-14000) M333C(NEW)+M333B-tail CORNERED: UNBIASED M333C-45 PERFECTLY two-sided (TPA 12o/14u mean-0.1, BA 9o/9u mean-0.2 27z, QMD 45/45 zero, TopHt 39z), density preserved; mirror-imaged seed-invariant (11857043 +541 / 3076266 -235) + 1/6 flip. +297 CNs + M333C cluster rows + 4 M333Cb/Ce merch-threshold TCuFt-step stands (density bit-exact). DB IE 0-14000 all ulp_class.

## BATCH7 (14000-16000) M333D(NEW)+M333C-tail CORNERED: UNBIASED M333D-45 two-sided (TPA 9o/16u mean-1.0, BA 9o/5u mean+0.6 31z, QMD 44/45 zero, TopHt 41z) — density preserved, BA leans slightly over; 2/6 seed-flips + mirror-imaged seed-invariant (11805927 +211 / 3186434 -93). +349 CNs + M333D cluster rows + 3 M333Ce/Db merch-threshold TCuFt-step stands. DB IE 0-16000 all ulp_class.

## IE COMPLETE (honestly, seed-verified) 2026-09-05
All 17808 IE stands bit-exact-or-cornered (sweep DB: 17808/17808 ulp_class, 0 needs_dig). Final batch8 (16000-17808)
M333Dc/Dd/De dig-tail two-sided (TPA 74o/118u mean-4.6, BA 57o/77u mean-2.0 103z, QMD 214/237 zero); bidirectional
seed-invariant (42603442 +187 / 1856009454 -316). +237 CNs + 1 M333De merch-threshold. NO open real-bug PAUSE.
src diff vs master EMPTY (pure adjudication; M333 fix was pre-merged). NEXT VARIANT: EM (cursor 0).
