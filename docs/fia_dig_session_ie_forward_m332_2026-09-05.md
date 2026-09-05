# IE FORWARD SWEEP (cursor 2000→) — batch 1 (2000→4000): M332 stratum adjudication (2026-09-05)

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7
Branch `wt/western-sweep` (worktree /workspace/.wt-western). Engine = master (merge of b9146275:
IE OH sp22 CRVAR small-tree log-DK diameter fix). `git diff master HEAD -- src/` EMPTY. Oracle FVSie_g16.

## PRE-FLIGHT: the OH-monoculture bug from resweep #2 is FIXED by b9146275
The resweep-#2 PAUSE was the dense-OH-monoculture small-tree under-build (stand 22404926010497). The b9146275
CRVAR log-DK fix RESOLVED it — per-cycle @2054 (FVSie_g16 vs jl):
  PRE-fix : TPA 989/828  BA 75/44  SDI 203/128  CCF 106/62  QMD 3.7/3.1   (jl UNDER-builds, one-directional)
  POST-fix: TPA 989/999  BA 75/79  SDI 203/212  CCF 106/111 QMD 3.7/3.8   (jl tracks; QMD converges)
The one-directional under-build is GONE (QMD now converges 3.7 vs 3.8; residual is a small ±slight-over count).
⇒ the resweep-#2 open bug is closed on master; the forward sweep proceeds.

## BATCH 1 (2000→4000, DIGCAP=100, regime none): 420 dig-worthy → DIG_PAUSE
bit_exact=0/2000. All 420 escalations are in NEW ecoregions (not the cornered 331A/342I):
  M332Ab 137, M332Ac 117, M332Aa 91, M332Ao 39, M332Ba 33, M332Bb 3  (Bitterroot/Belt Mtns, MT — Northern Rockies).
Worst-col: TCuFt 177, TPA 154, CCF 35, BA 30, BdFt 14, MCuFt 6, SDI 2, QMD 2.
The stratum is the GENERAL forested population (ages 0→379, mostly mature 55-200), not a pure-establishment stratum.

## ADJUDICATION — same near-SDImax count-straddle + #206 growth straddle as the cornered 331A/342I
### Population signed signature (jl−oracle, last common cycle), TWO independent unbiased-40 samples:
| metric | sample A (n=40)              | sample B (n=40)             | verdict |
|--------|-----------------------------|-----------------------------|---------|
| BA     | 9o/7u/24z  mean +1.0        | 11o/8u/21z  mean −0.5       | **BALANCED** (density two-sided) |
| QMD    | 0o/0u/40z  mean 0.1         | 2o/1u/37z  mean 0.1         | **CONVERGES** |
| TopHt  | 1o/4u/35z  mean −1.3        | 4o/2u/34z  mean 0.2         | balanced/converges |
| TPA    | 6o/15u/19z mean −2.9        | 13o/20u/7z  mean −13.3      | two-sided, mild UNDER-lean (count knife-edge) |
| TCuFt  | 2o/1u/37z  mean +7.8        | 2o/9u/29z   mean −204       | mostly zero; tail leans under (tracks BA/count) |

BA (the density/growth bug-discriminator) is TWO-SIDED and mean≈0 in BOTH samples — the SAME signature that,
post-crosswalk, licensed cornering 331A/342I. QMD converges, TopHt converges. TPA leans mildly under (jl thins
slightly more at the self-thin knife-edge) but the sign is two-sided — identical to the cornered 331A behaviour.

### Per-cycle mechanism (FVSie_g16 vs jl), mature stand 3296439010690 (age95, M332Ac):
  2000: TPA 2070/2070  BA 153/153  (identical start)
  2010: TPA 1799/1797  BA 198/185  ← **TPA matched, BA already −13** ⇒ a DIAMETER-GROWTH-DEVIATE difference
        (matched mortality, diverging diameter growth) = the #206 OLDRN serial-correlation growth straddle
        hallmark, NOT a mortality bug. Compounds to BA 263/246, TPA 923/834 by 2050; QMD 7.2/7.4 (jl slightly HIGHER).
Direction is stand-specific & seed-invariant (density-triggered self-thin mortality + growth deviate are both
deterministic) but TWO-SIDED across the population: escalated mature stands split 273528477489998 TPA **+98**,
531016134126144 +16, 3281599010690 +18 (OVER) vs 3296439010690 −89, 3285544010690 −97, 3288724010690 −53 (UNDER).

### Seed test (24 escalated stands, seeds 55329/1/777/424242): the residuals are seed-INVARIANT — EXPECTED,
because both drivers (SDImax density-dependent self-thin mortality; the #206 growth deviate assignment) are
deterministic given the stand. Per the seed_test caveat, seed-invariance here is NOT a bug signal; the
population TWO-SIDEDNESS (above) is the decisive test, and it is balanced.

### Volume (177 TCuFt escalations): NOT an independent volume-eq bug and NOT OH-driven.
The TCuFt-escalated sample is CONIFER-dominated (species 108/19/202/321/17; only 1 of 12 has any 768→OH). Its
volume residual tracks BA/count downstream (jl under-projects TPA/BA/TCuFt TOGETHER on the escalated tail) —
i.e. the same self-thin count-straddle showing up in volume, plus small-base % inflation on young stands
(the extreme "100% / 25886-cuft" row 1627644737290487 is a degenerate AGE-0 1-tree/6-TPA establishment stand).

## VERDICT: CORNERED (seed-verified two-sided near-SDImax count-straddle + #206 growth straddle)
Same ecoregion-independent primitive cornered for SN (global), TT (M331D/M331Ja/M332Er), IE (331A/342I).
- 420 M332 CNs → docs/fia_cornered_stands.txt (synced /workspace/FVSjl/docs — the per-CN unconditional-drop list).
- 4 cluster rows → docs/fia_cornered_clusters.tsv: (M332A|M332B) × (structure_densephase|count_divergence_UNVERIFIED),
  tag `ie-sdimax-count-straddle`.
- DB reclassify → the 374 M332 needs_dig rows become ulp_class. Dig queue cleared. DIG_PAUSED cleared. Sweep resumed.

## HONEST RESIDUAL NOTE (sub-material, watched not cornered-away)
There is a consistent MILD TPA/TCuFt under-lean in the unbiased mean (BA sign balanced, but the under-side
magnitudes exceed the over-side). This is within the count-straddle taxonomy (density BA preserved, count at the
knife-edge) and matches the cornered 331A behaviour, so it does NOT meet the one-directional bug bar. It is
flagged to re-examine if it GROWS into a BA-one-directional signal in a later stratum (that would flip to a bug).

## STILL-OPEN forward leads (not in this batch):
- Dense-seedling under-thin (jl ~2x OVER on TPA/BA/SDI): stand 103405608010661 (ecoregion 331Mf, 22801 TPA@0.2"
  → 2056 TPA 7498 vs oracle 3874, BA 306 vs 219). One-directional OVER — a DIFFERENT signature from this
  balanced straddle; will surface when the sweep reaches 331M*. Characterize/fix there.
