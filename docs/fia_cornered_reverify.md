# Re-verification of the FIA "cornered" (ulp_class) bucket — per-tree audit (2026-07-18)

**Question (user directive):** are ALL non-bit-exact FIA stands *irreducible* fp-ulp / cornered primitives, or
do real reducible divergences hide in the auto-cornered bucket?

**Answer: NO — real per-tree growth divergences exist in the cornered bucket.** Robustly verified (surviving
every methodology fix below): SN `220315381010854` **17.3%**, NE `1809173020290487` **14.5%** per-tree DBH
divergence on the SAME physical tree with the tree set byte-identical — i.e. genuine diameter-growth
divergence, not the self-thinning/AVHT40 tie-break the classifier assumed.

## How the bucket breaks down (sweep DB `data/fia_sweep.db`, 1.47M stands)
- bit_exact = 1,057,715 ; non-bit-exact = 410,661.
- The non-bit-exact split into STRUCTURAL (TPA/BA/SDI/CCF/TopHt/QMD) vs VOLUME axes (DB columns
  `struct_max_rel_pct` / `vol_max_rel_pct`). Structural tail: **118,183 stands >1%, 51,829 >2%, 9,134 >5%.**
- The huge headline rel% (billions of %) are ALL in the VOLUME columns = the **NE VOLINIT extreme-height
  volume-zeroing FVS bug** (live reports volume 0 at TopHt>200ft, FVSjl computes correct nonzero — FVSjl
  correct, cornered as FVS-bug; verified on NE `1167721956290487`).

## The classifier was never per-tree-verified for the full tail
The full-population sweep auto-cornered the structural tail as "structure_densephase" (self-thin RDPSRT
tie-break). Only a distilled ~29-stand LS candidate subset ever got per-tree verification (dig_verify_treeid).
The remaining ~9k (struct>5%) — and 118k (>1%) — were cornered by heuristic, NOT verified.

## Per-tree verifier (test/harness/fia/verify_treeid_magnitude.jl → TreeIndex version /tmp/verify_ti.jl)
Discriminator: run live+jl with TREELIST, match trees, compare per-tree DBH at cycles where the tree SET is
byte-identical (no tie-break yet) — any DBH divergence there is PURE growth divergence, not tie-break.
**Three methodology traps found & fixed (each inflated false positives — caught by direct dumps):**
1. `(SpeciesFIA,TreeId)` key COLLIDES (TreeId is per-plot, repeats; dead HISTORY=6 records) → mismatched
   pairing → fake 1.4" diffs. Fixed → `(SpeciesFIA,PtIndex,TreeId)`.
2. That STILL collides on multi-count/seedling-split stands (32 trees → 18 keys) → fake 154% diffs. Fixed →
   `TreeIndex` (unique per stand; species+cycle-1-DBH consistency guard confirms correct pairing).
3. Exact `!=` flags sub-0.01" ULP DBH rounding → magnitude threshold (>1% rel, set-identity gated).

**Validation (8 worst, reliable TreeIndex key): 5/8 ESCALATE** incl. SN 17.3%, NE 14.5% — `reliable=true`.
Scaled rate over 160 worst-struct stands: pending (`.sweep_work/cornered_reverify_scale.txt` for the earlier
inflated run; TreeIndex re-run in progress).

## Anchor case — SN `220315381010854`, species 762 (black cherry), tree 1762001
3-tree stand, tree set identical every cycle. DBH bit-exact through 2015 (1.279"), then diverges:
2020 live 2.212/jl 2.202, 2025 3.06/3.00, 2030 **live 3.786 / jl 4.439 (17%)**. Per-cycle DG: 2020 live 0.784
vs jl 0.776 = **1.1% from an identical 2015 state** — a real blend-zone (DBH 1–3", REGENT small-tree) DG
divergence, NOT compounded ULP. Mechanism localized to the SN REGENT small-tree growth (HTCALC height
increment → htdbh inverse; BACHLO ±10% random effect is off at default DGSD=0). Fixable-vs-FVS-correct: TBD
(needs an instrumented both-sides DG trace, like the LS REGENT FIX #7/#8 digs).

## Status
The premise "all cornered = irreducible ULP" is **DISPROVEN**. Scope (reliable rate) + mechanism root-cause
(is FVSjl or live the correct side?) are the open work — a REGENT/small-tree DG dig per variant, analogous to
the LS REGENT fix campaign. This is a real continuation of the FIA-compat campaign, not a closed item.
