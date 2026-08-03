# Western-variant FIA-plot validation — results

> **★ RE-VALIDATION ADDENDUM (2026-08-03, later than the sweep below).** The tables further down are **STALE** —
> substantial volume + height-growth work landed AFTER this sweep was recorded. Verified with current code:
> - **UT height growth sp 20/21/22 CRASH — FIXED.** The 3 repro stands (42631952010690 / 198812416020004 /
>   31533289010690) now run clean (CR-surrogate 17:19,22 + MC/BI 20:21 ported). No jl crash across the western set.
> - **"jl=0 missing-volume-species" subclass — largely FIXED.** Repro CI 51052888020004 (juniper sp64): live
>   TCuFt/MCuFt 282/278, **jl now 282/278 (was 0/0)**. EM 2342318010690: **now bit-exact 82/29/120 (was 23-34% low)**.
> - **Current >10% volume-divergence rate is FAR lower than the tables below** (measured on the sub.db treed
>   sample): CI 1/15, EM 0/3, UT 1/8 — vs the recorded CI 24 / EM 8 / UT 25. BM added to multi-cycle (median 0).
> - **Sole remaining volume residual = FW2 merch/board low on large *stunted* DF** (old-growth, low H/D). Per-tree
>   (CI 3200106010690, DBH 34.1″/HT 67ft): TotCu 130.7≈134.4 ✓ but **MchCu 115.8 vs 129.4 (−11%), BdFt 600 vs 740
>   (−19%)** — live merch-fraction ~0.97, jl ~0.87 ⇒ jl places the 6″ merch top too low. Total-cubic matches, so
>   it's a merch top-height / log-segmentation taper-shape issue, NOT missing species. Narrow (≈1/15 CI stands
>   >10%), scoped for a focused fvsvol-taper Fortran trace (do NOT cargo-cult; instrument live merch height first).
> - **Growth core unchanged: density/diameter 100% bit-exact; CCF(±1)/TopHt(AVH-tie) cornered.**
> Bottom line: the western cluster is bit-exact-or-cornered on growth everywhere; volume is now bit-exact-or-
> cornered except the narrow large-stunted-DF FW2 merch/board class. See memory `fvsjl-extensions-rollout`.

Companion to `RESULTS.md` (SN/NE/CS/LS). Seeded stratified random FIA sweep across the six newly-ported
**western** variants — **CR, IE, EM, UT, CI, TT** (BM skipped — mid-port) — validated cycle-0 all-10-`.sum`-
columns vs the freshly-relinked live FVS binaries.

Harness: `sample_stands.jl` (now takes a target-variant list) → `extract_subset.jl` (indexed subset DB) →
`run_sweep_western.jl`. Live oracles: `/workspace/.{cr,ie,em,ut,ci,tt}work/FVS??_clean`. Singletons via
`FVSjl.variant_from_code`. A cell diverges only if the printed `.sum` value differs by **≥0.5**.

## Method note — why a *treed-only* denominator (differs from the eastern harness)
Western live FVS emits an **all-zero `.sum`** for treeless FIA conditions (eastern live emits *none* → "no-sum").
The base `run_sweep.jl` therefore counts western treeless stands as trivially-bit-exact (0 == 0), which inflates
both the denominator and the rate. `run_sweep_western.jl` **excludes treeless** (year-0 TPA & BA both ≈0) and
reports the honest **treed-only** rate, additionally splitting **GROWTH** cols (TPA/BA/SDI/CCF/TopHt/QMD) from
**VOLUME** cols (TCuFt/MCuFt/SCuFt/BdFt) and bucketing the max volume relative-Δ per stand (<2% / 2-10% / >10%).

Treeless conditions are the large majority of random western FIA COND rows (CR/EM/UT/CI ~85-90% treeless), so
the base-harness "both-produced" rate for these variants is dominated by trivial zeros — do **not** read it as
model fidelity. Sample: 140/variant (seed 20260803) + a 500/variant augmentation for the sparse-treed variants
CR/EM/UT/CI (seed 20260804), deduped by STAND_CN.

## Results — treed-only, cycle-0 all-10-col vs freshly-relinked live
| Variant | treed | ALL-10 bit-exact | rate | GROWTH-exact | VOL-exact | live-crash | treeless-excl |
|---|---|---|---|---|---|---|---|
| CR | 119 | 102 | **85.7%** | 117 (98.3%) | 103 (86.6%) | 0 | 521 |
| IE | 101 | 75 | **74.3%** | 96 (95.0%) | 78 (77.2%) | 0 | 39 |
| TT | 57 | 44 | **77.2%** | 56 (98.2%) | 45 (78.9%) | 3 | 80 |
| UT | 106 | 67 | **63.2%** | 105 (99.1%) | 67 (63.2%) | 0 | 507 |
| EM | 47 | 18 | **38.3%** | 44 (93.6%) | 19 (40.4%) | 0 | 592 |
| CI | 182 | 35 | **19.2%** | 138 (75.8%) | 47 (25.8%) | 7 | 445 |
| **TOTAL** | **612** | **341** | **55.7%** | | | 10 | 2044 |

Full per-stand offender lists (with per-stand max volume-Δ%) in `western_sweep_raw.txt`; the treeless-inclusive
140/variant batch-1 in `western_batch1_raw.txt`.

## The growth core generalizes; volume does not (for the woodland-heavy variants)
### GROWTH side — bit-exact-or-cornered everywhere
Across all **612 treed stands, the density/diameter core (TPA / BA / SDI / QMD) is 100% bit-exact — zero
mismatches in any variant.** The only growth-side divergences are the two documented **cornered** classes:
- **CCF** — 43 stands (41 CI, 2 IE). Drilled cases are all **±1 integer-boundary** rounding (live 32 vs jl 31;
  live 1 vs jl 0 on sapling stands). CI hits it often because its sample is rich in low-density/sapling stands
  sitting on the 0.5 crown-width rounding boundary. Cornered, not a model error.
- **TopHt** — 13 stands (all variants), the AVH sort-tie class (±1-2 ft).

So the ported **growth models validate against random FIA inventory** at 93.6–99.1% (the sub-100% being CCF/
TopHt ties). CR/TT/UT growth is 98–99%+; CI's 75.8% is the CCF-tie frequency, not a growth-model gap.

### VOLUME side — the divergence axis (dominant for UT/CI/EM)
Volume drives essentially all real divergence. Max per-stand volume-Δ buckets [<2%, 2-10%, >10%]:
| Variant | <2% (cornered) | 2-10% | >10% (real) |
|---|---|---|---|
| CR | 14 | 1 | 1 |
| IE | 14 | 8 | 1 |
| TT | 4 | 8 | 0 |
| EM | 11 | 9 | 8 |
| UT | 5 | 9 | 25 |
| CI | 47 | 64 | 24 |

- **CR / IE / TT** — volume mostly **cornered** (<2%, board/merch rounding) with a thin 2-10% tail and ≤1 large
  outlier. These three are effectively bit-exact-or-cornered on random FIA, consistent with their COMPLETE status.
- **UT / CI / EM** — a **real, systematic volume divergence class** on species/stands outside each variant's
  reference stand (utt01/cit01/emt01). Two sub-classes, from per-tree drills:
  1. **jl computes ZERO volume where live has volume** (the "Δ=100%" offenders). Repro **CI 51052888020004**:
     live TCuFt/MCuFt = 282/278, **jl = 0/0** (QMD 13.4, real trees) — a species with no ported volume equation.
     Several UT/CI stands show this. This is a genuine **missing-volume-species gap**.
  2. **jl systematically LOWER** than live. Repro **CI 3200106010690** live 659/623/2124 → jl 531/447/1758
     (~20-28% low); **EM 2342318010690** live 82/29 → jl 63/19 (~23-34% low); **UT 198895166020004** small-tree
     TCuFt live 28 → jl 6. A volume-equation / DVEW-woodland-coefficient class for non-reference species.

  These match the "Remaining" notes in the port memories (UT: *non-utt01 PJ-DVEW/hardwood pure stands + TCuFt
  tail*; CI: *vol ±3-10% tail*) — but at population scale the tail is **worse than "±3-10%"**: up to zero-volume
  (100%) and 30-85% on woodland species the reference stand never exercises. SCuFt never diverged (0 mismatches
  anywhere) — the split is entirely in TCuFt/MCuFt/BdFt.

## jl errors (real port gaps) — reported, not suppressed
**UT height growth missing for species 20/21/22** ("CR-surrogate / MC-BI, not yet ported — not in utt01").
3 stands threw; repro STAND_CNs:
- `42631952010690` — UT sp 22
- `198812416020004` — UT sp 20
- `31533289010690` — UT sp 21

`FVSjl.run_keyfile` throws (does not silently mis-model) — a clean, actionable gap: port `height_growth!` for
UT woodland species 20/21/22. **No other variant produced any jl crash/KeyError** across 612 treed stands.

## Live-FVS crashes (excluded, jl runs them)
CI 7, TT 3 (live SIGFPE, no oracle → excluded from the denominator, as in the eastern doctrine). CR/IE/EM/UT: 0.

## Reading
- **Growth cores of all six western variants are validated on random FIA inventory** — density/diameter 100%
  bit-exact, residual growth divergence = only the CCF (±1) and TopHt (AVH-tie) cornered classes.
- **Volume generalizes for CR/IE/TT** (bit-exact-or-cornered) but **not for UT/CI/EM**, which carry a real
  volume class outside their reference stands: (a) missing-volume species (jl=0), (b) systematic under-volume on
  woodland/DVEW species. This is the honest population-scale finding — the memories' "COMPLETE" is anchored on
  the single reference stand (utt01/cit01/emt01); the FIA population exercises far more species.
- **One hard port gap**: UT height growth for sp 20/21/22 (3 repro CNs above).
- CR is the strongest western drop-in on FIA (98.3% growth / 85.7% all-10, volume residual almost entirely <2%).

Reproduce: `sample_stands.jl 140 20260803 <dir> CR,IE,EM,UT,CI,TT` → `extract_subset.jl <dir> <db>` →
`run_sweep_western.jl 60 <dir>=<db> [<augdir>=<augdb>]`.
