# IE FIA RE-SWEEP dig session — master (77b97a90, IE AUTOES complete): establishment #206 straddle CORNERED

**Date:** 2026-09-04  **Variant:** InlandEmpire (IE)  **Oracle:** `/workspace/.iework/FVSie_clean`
(harness default; byte-identical to `FVSie_g16`).  **Engine:** sweep branch `wt/western-sweep` synced to
master (77b97a90 present; `git diff master HEAD -- src/` EMPTY).  **DIGCAP=100 (never raised).**

## Purpose
Re-sweep IE from cursor 0 on the FROZEN master engine (IE AUTOES cohort model complete: b2abab52 +
a3cda0f4 + b6edfd9a + 77b97a90) and CORNER the now-MEASURED establishment residual (#206 OLDRN/ZRAND
realization straddle), or STOP+report if it is a separate real bug. 77b97a90 changed only comments + a
partition-locking test (`test_ie_estab_cohort_wk4.jl`); the projection engine is functionally identical to
b6edfd9a, so the .sum results reproduce b6edfd9a exactly (587 needs_dig over cursor 0–2000).

## Sweep result (batch-1, cursor 0→2000)
Cap TRIGGERED at cursor **2000** (dig-queue ≥ DIGCAP 100; DIG_PAUSED set). Batch: 2000 emitted,
bit_exact **0/2000** (the IE #206 straddle is ubiquitous ⇒ every stand diverges by ≥1 ULP), dig_class
**587 needs_dig / 1413 ulp_class** (bit-identical to the b6edfd9a re-sweep — confirms 77b97a90 is a
no-op on the engine).

- Signature: **585/587 structure_densephase** (+1 count_straddle, +1 threshold_crossing — both TCuFt
  volume-threshold members of the same realization family, struct_abs 9/1 = structure fine).
- worst_col (needs_dig): **BA 303, CCF 139, SDI 69, TPA 39, TCuFt 37** — DENSITY-dominated (the b6edfd9a
  WK4-distribution port moved the residual off merch-volume onto crown/BA, as expected).
- **ALL 587 needs_dig ESCALATE** (bypass a cluster corner: structure_densephase w/ density worst_col &
  max_rel≥15% & struct_abs≥10, or TCuFt w/ vol_abs≥300) ⇒ per-CN cornering required (the TT precedent).
- ecoregions (needs_dig): **331Af 230, 331Ac 165, 331Aa 98, 342Ia 76, 342Ib 18** — the front-loaded IE
  bare/sparse AUTOES-establishment strata in cursor 0–2000.

## SIGNATURE VERIFICATION — CONFIRMED #206 OLDRN/ZRAND establishment realization straddle
Signed last-cycle residual (jl − oracle) via `.sweep_work/sig_verify.jl` (reuses ledger_fia.jl
run_live/keytext/parse_sum10, regime `none`, 5 cycles):

1. **Adjacency triple (77b97a90's example) EXACTLY reproduced:** 3334193010690 **BA +10 / CCF +11 /
   SDI +6**; 3334202010690 **−11 / −11 / −17**; 3334212010690 **−13 / −12 / −27**. Adjacent
   near-identical inventory plots FLIP BA/CCF/SDI sign.

2. **Unbiased 40-stand all-IE sample: BALANCED** — BA over=15 / under=18 / zero=7, **mean −0.17,
   median 0.00**, |mean| 3.58. Matches 77b97a90's 86-stand unbiased (38 over / 44 under, mean −0.85).
   The true population behavior is a two-sided ~50/50 straddle with a ~1.5% under-lean.

3. **Full-run adjacency (all neighbors, escalated + ulp): different survey clusters lean OPPOSITE
   directions, adjacent plots flip sign** — the decisive realization proof:
   - suffix 290487: `+8u −14E +0u +2u +7u +0E +2u −3u +0E −9E +1E +9u` — neighborhood leans OVER
     (ulp 5 over/1 under); a −14 escalated stand sits beside a +8 neighbor. Near-identical inputs → +8
     AND −14 = RNG realization ONLY (a deterministic input→output bug cannot).
   - suffix 489998: `+2u +4u +0u +2u −2E +3E +4u +0u +3u +4u −1u +0u` — neighborhood leans OVER.
   - suffix 010690: leans under (`−1u +1E −14E −13E −13E −13E …`) — but contains +1E and the
     escalated-only run flips (−15/−15/+3/−1/+1/−1).

4. **Escalated needs_dig subset leans under (random-40: 8 over/32 under; escalated-only adjacency:
   2 over/28)** — this is the EXPECTED SELECTION-BIAS artifact: the ±10 struct_abs escalation floor
   preferentially catches the heavier under-tail of a slightly negatively-skewed straddle (|mean| 3.58
   unbiased vs 6.08 escalated). It is NOT one-directional: 20% of the escalated tail is still OVER, and
   its true neighbors straddle/lean over (evidence 3).

## OVER-CORNERING GUARD — CLEARED
No escalated establishment stand is one-directional AND correlated with a deterministic input:
different survey clusters (290487/489998 over-leaning, 010690 under-leaning) share the SAME establishment
regime/ecoregions/model and differ only in the per-stand RNG realization; within a cluster adjacent plots
flip large-magnitude sign (290487: +8/−14/+7/−9). 77b97a90's instrumented measurement (FVSie_g16 estab.f
DO 33/228, per-tree booked WK4+PROB, 86 stands) already showed the jl per-class partition is FAITHFUL
(advance/subsequent/excess within 1-3pp) and the per-stand advance-fraction discrepancy is UNCORRELATED
with BA sign (a stand booking LESS advance yields both +19 and −15). The prior b6edfd9a "one-directional
under" verdict was the selection-biased needs_dig subset + an INFERRED (not measured) partition bias —
refuted here and by 77b97a90. Re-balancing toward the advance class would BREAK faithfulness to the FVS
source without touching the straddle. **NOT a re-balanceable bug; the established IE #206 corner.**

## Named primitive
`FVSie` OLDRN/ZRAND (BACHLO) per-tree normal-deviate draw in the AUTOES establishment-cohort growth
(esgent.f WK4² birth growth + estab.f DO 99/33/228 ESADVH/ESSUBH/ESXCSH height-class booking), same
#206 realization-straddle class as the TT M331D smhtgf-ZRAND corner. jl draws-then-distributes vs FVS's
draw/triple order → the same deviate set reassigned across records → a balanced two-sided .sum straddle on
the bare/sparse establishment cohort's density columns.

## Action taken
- Cornered the **587** batch-1 needs_dig CNs per-CN in `docs/fia_cornered_stands.txt`
  (677 → 1264 lines), tag `# IE establishment #206 OLDRN/ZRAND straddle (<worst_col>|<struct_abs>)`.
- Added **5** cluster rows to `docs/fia_cornered_clusters.tsv` (331Aa/331Ac/331Af/342Ia/342Ib ×
  structure_densephase, label `ie-estab-206`) so future IE-establishment strata below the escalation
  floor drop without re-pausing (the escalation guard still surfaces any genuinely NEW ≥10-unit bug).
- Synced both registries to the live-filter path `/workspace/FVSjl/docs/`.
- Reset the dig-queue (all 584 dig-worthy rows are now cornered); removed DIG_PAUSED; **resumed** the
  sweep past cursor 2000.
- DIGCAP NOT raised. No engine change; gate stays 339/11.

## Backups
- pre-resweep DB (2000-row b6edfd9a-state): `data/fia_sweep_west.db.bak_pre_ie_resweep_master77b_*`.
- pre-resweep dig-queue: `docs/fia_dig_queue_west.csv.bak_pre_ie_resweep_master77b_*`.
