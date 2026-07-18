# Re-verification of the FIA "cornered" (ulp_class) bucket — per-tree audit (2026-07-18)

**Question (user directive):** are ALL non-bit-exact FIA stands *irreducible* (fp-ulp / cornered primitive),
or do real reducible growth divergences hide in the auto-cornered bucket?

## FINAL ANSWER: no reducible bug found — the bucket holds up as cornered/irreducible.
A rigorous per-tree audit did **not** surface a reducible FVSjl growth bug. Every divergence examined resolves
to a documented primitive: the **AVHT40 top-height tie-break**, the **self-thinning RDPSRT tie-break**,
small-denominator **rounding** (BA 0↔1 on seedlings), or the **NE volume-zeroing FVS bug** (FVSjl correct).

## IMPORTANT — a mid-audit false alarm, retracted
An intermediate pass reported "real per-tree growth divergences up to 17%." **That was wrong** — it was a
per-record *pairing artifact*. Root cause (found by eyeballing raw treelist rows): once a stand self-thins,
**FVS triples/splits each tree into many sub-records that all share the same `(SpeciesFIA, PtIndex, TreeId)`**
(e.g. SN 220315381010854 sp762 → 9 records at 2030, DBH 3.54–4.50). Any (species,id)-keyed comparison then
pairs mismatched sub-records → fake 17–650% "divergences." Three successive key fixes (TreeId → +PtIndex →
TreeIndex) each still leaked because `TreeIndex` is re-assigned after splitting. The only trustworthy per-tree
window is **before any split/mortality** (full initial count on both sides); there per-tree DBH is **bit-exact**.

## What the divergences actually are (verified stand-level, diff_one)
The sweep's `struct_max_rel_pct` is stand-level and reliable. Re-examined worst cases all show the same shape:
- **Cycle-1 TopHt-only divergence** then a cascade. E.g. NE 1809173020290487 TopHt 13/22, NE 381452868489998
  72/76, CS 225042621010661 40/47 — everything else bit-exact at cycle 1. This is the **AVHT40 top-height
  tie-break**: on dense tied-DBH seedling stands, which trees fall in the largest-40-TPA set is decided by
  FVS's unstable RDPSRT quicksort; the port (stand_top_height double-RDPSRT) matches most stands but is
  stand-dependent (`standstats.jl:127-145` documents the exhaustive prior single-vs-double dig — "no global
  sort choice is bit-exact"). The TopHt seed then cascades into TPA/BA via height-driven competition/mortality.
- **Identical self-thinning with rounding**: SN 220315381010854 — TPA bit-exact every cycle, only BA 0↔1
  (seedling stage), QMD 2.2↔2.3, CCF ±2, TopHt ±1. The DB's "17%" was a small-denominator artifact.
- **Volume tail (billions of %)** = NE VOLINIT extreme-height volume-zeroing (FVS bug; FVSjl correct).

## Method note (for future audits)
Per-record treelist comparison is INVALID after FVS record-splitting/tripling. Valid tests are: (a) stand-level
`.sum` aggregates (tripling-invariant), or (b) per-tree DBH strictly in the pre-split/pre-mortality window
(both sides at full initial count). The verifier `test/harness/fia/verify_treeid_magnitude.jl` and the
no-mortality-window variant encode (b); trust their `worst%` ONLY within that window.

## Residual honesty
The tie-break is a discrete unstable-quicksort permutation on EXACTLY-tied (0.1"-rounded FIA) DBHs — cornered,
not literally "fp-ulp," but no *reducible FVSjl bug*: the trees are bit-identical; only FVS's internal sort
order on ties differs, and the port already matches the majority. Whether a still-more-faithful RDPSRT could
close more of the AVHT40 residual is the only open thread; prior digs concluded it is stand-dependent and no
single sort choice is globally bit-exact. Net: **I could not find a reducible growth bug in the cornered bucket.**
