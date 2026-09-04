# IE HONEST RE-SWEEP with the FIXED harness (2026-09-04)

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7
Branch: `wt/western-sweep` (worktree /workspace/.wt-western). Engine = master (`git diff master HEAD -- src/` EMPTY),
i.e. the D2 habitat fix 37c66546 is IN. Oracle = /workspace/.iework/FVSie_g16 (fresh, 18:25). Regime `none`.
Harness fix present: ledger emits `count_divergence_UNVERIFIED` (no auto-corner) + seed_test.jl adjudicator.

## PURPOSE
Re-sweep IE from cursor 0 on the FIXED harness after the 3360-stand "establishment #206 straddle" corner was
adjudicated FALSE and REOPENED (registry 722, no ie-estab entries). Split the escalated cluster into genuine
seed-VARIANT straddle (corner) vs seed-INVARIANT one-directional real bug (do NOT corner), by BOTH the
seed-perturbation test AND — the decisive test for establishment stands (ESRANN fixed ⇒ seed-invariant by
construction) — the POPULATION one-directional signed signature.

## RESET performed
IE file cursor 17808→0; DB progress ie→0; DB `sweep` rows for IE (found 17808 stale rows under variant "IE"
from the old falsely-cornered full sweep — the earlier reset only checked lowercase) DELETED; dig-queue cleared
to header (203 pre-existing **EM** rows saved to .sweep_work/em_digqueue_saved_pre_ie_resweep.csv — the EM sweep
is separately paused on its own real bug); .sweep_work/DIG_PAUSED (EM) removed. Backups:
data/fia_sweep_west.db.bak_pre_ie_honest_resweep_20260904_215141 + the queue .bak of the same ts.

## SWEEP RESULT (cursor 0→2000, DIGCAP=100)
Batch capped after ONE cycle: **bit_exact=0 / 587 comparable stands; 580 needs_dig (escalated), 7 ulp**.
Cluster = the 331A / M332 Idaho-Batholith establishment-heavy stratum. Signature: 585 `structure_densephase`
(worst_col BA 303, CCF 139, SDI 69, TCuFt 38, TPA 31). struct_max_abs median 33, p90 71, max 625.
DB now honestly holds 587 IE rows; cursor paused at 2000/17808.

## SEED-PERTURBATION TEST (12 stands, seeds 55329/1/777/424242)
`seed_test.jl IE …`  →  **9 REAL-DETERMINISTIC, 1 pure-straddle, 2 match.**
- Real-determ (count/Mode-B, seed-INVARIANT): 531011782126144, 672242349126144, 1627683478290487,
  486593012489998, 273527124489998, 195372582020004, 3109970010690, 22404926010497, 2978686010690.
- 3109970010690 (the D2 anchor): TopHt now −2 (was −28 pre-fix) = **D2 fix landed**; residual is pure
  Mode-B TPA −66 seed-invariant.
- 22404926010497 (bitter-cherry): TPA −574/−598/−621/−605, BA −60…−68 seed-invariant (BA does NOT converge).
- 2978686010690: was −1300 in the pre-fix adjudication; post-D2-fix now +15…+40 seed-invariant (habitat fix
  moved it, residual Mode-B remains).
- Pure straddle (TPA FLIPS sign with seed): 12285258010690 (the only one).
- Match/negligible: 3320901010690, 1143317126290487 (the old S1 anchor — now within noise at seed 55329).

## POPULATION ONE-DIRECTIONAL SIGNATURE — the decisive test (40 unbiased escalated stands, default seed)
`.sweep_work/ie_popsig2.jl` — signed jl−oracle at the last common cycle:
| metric | over | under | zero | mean | verdict |
|--------|------|-------|------|------|---------|
| TPA    | 10   | 15    | 15   | −10.0 | mildly under; high-mag seed-invariant tail |
| **BA** | **1**| **23**| 16   | **−5.5** | **ONE-DIRECTIONAL UNDER (23:1)** |
| TopHt  | 2    | 12    | 26   | −1.9 | mild residual under-lean (D2 largely fixed; mostly zero) |
| QMD    | 0    | 0     | 40   | −0.2 | converges |

**BA is one-directionally UNDER 23:1** — categorically NOT a two-sided RNG straddle (a straddle is ~50/50).
DETERMINISM PROOF: five TRULY-IDENTICAL habitat-221 plots (249050692489998, 1287822146290487, 452196010497,
720726861290487, 450543010497 — identical oracle output 19.7/221/14) give the IDENTICAL residual BA −8 /
TopHt −5 / TPA −2 / QMD −0.6. Identical inputs → identical residual = deterministic, not RNG.

## VERDICT: REAL one-directional deterministic bug — NOT a cornerable straddle. IE is PAUSED (honestly).
Per task doctrine, corner an establishment cluster ONLY if its signed count residual is genuinely
two-sided/balanced across an unbiased sample. BA 23:1 under (+ identical-plot identical-residual + the
seed-invariant ±hundreds-TPA tail) is one-directional ⇒ **do NOT corner**. The FALSE 3360-stand corner is
NOT re-created. No genuine-straddle CLUSTER exists here (only the lone seed-flipping 12285258010690 is a true
S1 realization straddle — a swamped minority, left un-cornered pending the fix).

### Mode B (natural-mortality count-partition) — rate/partition bug, NOT an RDPSRT tie-break
Discriminator (measured): an RDPSRT unstable-sort tie-break reorders WHICH equal-key record dies but the SAME
NUMBER die ⇒ TPA is CONSERVED and the sign is BALANCED. Here the TPA total diverges by up to −574/−621 (and
struct_max_abs to 625) SEED-INVARIANTLY (count NOT conserved) and BA is one-directional 23:1 under — neither
is producible by a tie-break. ⇒ **real mortality-RATE / count-partition + density-under-build bug** (the number
of survivors and the BA they carry differ), NOT the RDPSRT primitive. PAUSE + fix, do not corner.

### bitter-cherry (22404926010497-class) — real dense-monoculture establishment/mortality bug
Input = ONE Other-Hardwood record (FIA sp 768 → OH), 1638 TPA @ 0.1" (dense hardwood monoculture). Oracle's
AUTOES (Regen Estab v2.0) fires and the stand builds BA 0→75 (TPA 1638→989 by 2054); jl over-thins the cohort
(TPA −574 seed-invariant) and under-builds BA (13 vs 75 = BA does NOT converge). Real deterministic bug in the
ultra-dense establishment/mortality path. PAUSE + report.

## FIX PLAN (coordinator; EM+IE, likely Wykoff-variant-general)
1. **Mode-B count-partition / density under-build (dominant, BA 23:1 under).** jl retains/builds less BA than
   the oracle on dense IE establishment cohorts, seed-invariantly. Localize the pure-projection density-dependent
   mortality-rate path (morts.f VARMRT rate + which/how-many small trees die) and the small-tree BA accrual on
   the established cohort; A/B per-cycle mortality TPA + BA against FVSie_g16 on a BA-converging Mode-B stand
   (3109970010690: TPA −66, BA ~0 — start there for the count-only path) and on a BA-diverging stand (22404926).
2. **Residual establishment height/density under-build (TopHt 12:2 under, BA −8 identical on habitat-221).**
   The D2 habitat-260 fix removed the worst TopHt deficit (3109970 −28→−2) but a mild cluster-wide one-directional
   BA/TopHt under-lean remains — likely a residual established/small-tree height-or-BA intercept (ESSUBH/ESADVH
   AA-vs-HT1) on resolvable-habitat stands. Trace the habitat-221 identical-plot cohort (BA −8 exactly).
3. **bitter-cherry dense OH monoculture** (22404926010497): jl AUTOES + ultra-dense self-thin under-builds the
   conifer cohort BA. Likely the same Mode-B path amplified at extreme density + a species-availability angle
   in AUTOES for an OH-only stand.
After the fix, re-sweep IE from cursor 0; a genuinely two-sided residual (if any remains, e.g. the sparse
bare-estab S1) is then re-cornered WITH the seed test.

## Artifacts (durable, on /workspace)
- Cluster cycle ledger: .sweep_work/expand/ie_west_cycle.csv (587 rows) + …/ie_west_ledger.csv master.
- Dig queue: docs/fia_dig_queue_west.csv (580 IE escalated rows = the dig worklist).
- Seed test log: .sweep_work/ie_seedtest_honest.log. Population signature: .sweep_work/ie_popsig2_full.log.
- Sample lists: .sweep_work/ie_popsig_cns.txt (40). Scripts: .sweep_work/ie_popsig2.jl.
- DB: data/fia_sweep_west.db (587 IE rows: 580 needs_dig / 7 ulp). Cursor: ie_west.cursor = 2000.
