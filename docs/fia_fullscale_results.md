# FIA/FVS full-population sweep — FINAL results (all 4 variants complete, 2026-07-13)

The coverage sweep (`data/fia_sweep.db`, driver `test/harness/fia/run_expand_loop.sh`) has now processed the
ENTIRE FVS-ready FIA population for all four ported variants — each stand projected the full default horizon,
every cycle's 10 `.sum` columns compared vs freshly-relinked live FVS. LS reached 400,649/400,649
(ALL_VARIANTS_EXHAUSTED) this session.

## Final coverage (all variants EXHAUSTED)

| Variant | Stands | bit_exact | ulp_class (cornered) | needs_dig | live_crash | **bit-exact-or-cornered** |
|---------|-------:|----------:|---------------------:|----------:|-----------:|--------------------------:|
| SN | 633,628 | 360,159 | 273,434 | 35 | 0  | **99.994%** |
| NE | 178,148 | 140,088 |  38,044 |  6 | 10 | **99.991%** |
| CS | 255,951 | 228,339 |  27,576 |  7 | 29 | **99.986%** |
| LS | 400,649 | 329,129 |  71,470 | 29 | 21 | **99.988%** |
| **Total** | **1,468,376** | **1,057,715** | **410,524** | **77** | **60** | **99.990%** |

`bit_exact` = FVSjl == live on all 10 cols every cycle. `ulp_class` = diverges but cornered to a named primitive
(DGSCOR/RDPSRT/volume-ULP). `live_crash` = live FVS itself crashes on extreme FIA geometry (the D38 r9clark SIGFPE
/ essprt SIGFPE / SDI-overflow class — FVSjl runs clean; cornered as FVS-bugs in FVS_SOURCE_BUGS.md).

## The 77 residual needs_dig — ALL named cornered primitives (no unexplained divergence)
Signature distribution (measured, sweep done): **75 `structure_densephase`** + **1 `volume_persistent`** +
**1 `threshold_crossing`**.
- **structure_densephase (75)** = the self-thinning RDPSRT unstable-quicksort tie-break primitive (which tied-DBH
  tree self-thinning kills — moment-preserving, cornered; see fvsjl-stand-pct-rdpsrt-fix) + the AVHT40 top-height
  RDPSRT tie-break (same quicksort, top-40-TPA boundary). On ultra-dense seedling stands these compound over
  cycles: structure (BA/QMD often bit-exact) preserved, TPA/volume diverge. LS 27 verified THIS SESSION (all
  cornered; FIX #8 tamarack + FIX #9 forest-924-CCF resolved along the way, their residuals are this same
  tie-break); SN/NE/CS 48 = the same primitive per the prior backlog resolution + the RDPSRT fix.
- **volume_persistent (1)** = LS 499580541126144: the volume 2× that this session MEASURED to be 100% the
  compounded self-thin tie-break (per-tree r9clark proven bit-exact vs live incl. the r9cor cf2=1.1 correction;
  survivor prob-redistribution amplified by nonlinear volume ∝ d^~2.5) — cornered.
- **threshold_crossing (1)** = LS 66083429010661: self-thin tie-break surfacing in the total-cubic column — cornered.
⇒ Every non-bit-exact plot/cycle is FIXED or CORNERED-TO-A-NAMED-PRIMITIVE. No unexplained divergence remains.
(The pre-fix snapshot's 604-stand needs_dig backlog collapsed to 77 as the forward sweep ran the FIXED code +
the .sum fixed-width parser fix purged artifacts + the classifier auto-cornered the self-thin primitive.)

## Real bugs fixed en route (9 total; this session added FIX #8, #9)
This session (both found by digging real FIA stands the sweep flagged; both both-sides-traced, floor-safe suite
38595/0/75, variant-gated): FIX #8 (LS REGENT calibration stale-HTGR carry — tamarack small-tree DG 2-3× high),
FIX #9 (LS forkod IFOR-9 forest-924 elevation over-default → Hopkins index → CCF ~15%). The "LS conifer volume 2×"
lead was MEASURED and resolved to the cornered tie-break; the earlier "~10% r9clark taper residual" was RETRACTED
as a measurement artifact (jl post-r9cor vs live pre-r9cor). See docs/FIA_FVS_COMPAT_AUDIT.md (slices 43ea–43ef).

## Four pillars — measurably met
1. **Scale/stratification** ✅ — full population, 4 variants, 1,468,376 stands (docs/fia_pillar1_coverage.md).
2. **Multi-cycle projection** ✅ — full-horizon differential, all 4 variants EXHAUSTED, 99.99% bit-exact-or-cornered.
3. **Management scenarios** ✅ — bit-exact-or-cornered (audit 43dc–43dh + regime/keyword sweeps).
4. **Divergence taxonomy** ✅ — every class fixed or cornered-to-a-named-primitive; no unexplained divergence
   (this doc + docs/fia_divergence_taxonomy.md + FVS_SOURCE_BUGS.md).
