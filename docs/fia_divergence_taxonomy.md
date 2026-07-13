# Pillar-4 — FIA/FVS-compat divergence taxonomy (consolidated)

Every non-bit-exact class the full-population sweep surfaced, both-sides-traced (FVS source + FVSjl) and either
FIXED (floor never regressed) or CORNERED to a named primitive. Detail lives in the per-slice `FIA_FVS_COMPAT_AUDIT.md`
entries and `FVS_SOURCE_BUGS.md`; this is the single-view index. "Bit-exact-or-cornered" holds for every class below.

## A. FVSjl bugs FIXED (real port gaps found by the sweep; floor held at 38587/0/75)
1. **FIA missing-SLOPE default** (43ay) — reader left SLOPE=0; FVS grinit.f defaults 5%. Over-grew sp39 DBH ~2×. Fix=grinit default. 4 variants.
2. **Calibration AVH backdating** (#2) — SN DGSCOR calibration used backdated AVH not current top height ⇒ up to 19% multi-cycle drift. Fix=stash/restore current avg_height. SN-scoped.
3. **stand_pct! RDPSRT tie-break** — percentile used a stable sort where FVS uses the unstable RDPSRT quicksort ⇒ mis-assigned self-thinning kill on tie-heavy stands. Fix=_rdpsrt!. Resolved 67 broken tests.
4. **10-yr cycle size-cap mortality** — BAMAX/size-cap used sqrt fint-yr diam_growth vs FVS's linear (DG/BARK)·(FINT/5). Fix=_mort_traj_g. s5/s9/timeint10.
5. **MORTMSB/MSBMRT** — alternate mature-stand-breakup mortality ported (was missing).
6. **LS dense-phase DG serial-correlation** — wrong DG measurement period (htg_period=10 vs FIA DG_MEASURE=5) ⇒ tripled DG ~3% high ⇒ self-thinning over-kill. Fix=gate meas_fint on growth_dg_set.
7. **LS REGENT HCOR calibration MISSING** (43dn) — calibrate_diameter_growth! had SN/NE/CS branches but no LakeStates ⇒ aspen height-grown 1.60× low ⇒ 2× under-thin on dense regen. Fix=added LS block (ls_htcalc+ls_balmod, backdated BA/QMD).

## B. Cornered to a named ULP-class primitive (Float32 op-order = semantics; FVSjl faithful, residual named)
- **Compounded-ULP self-thinning RDPSRT tie-break** — growth bit-exact per-record, but which tie-DBH trees die at the SDI
  threshold diverges via the unstable-quicksort permutation; BA bit-exact, TPA diverges (amplified to large TPA% on
  ultra-dense tie-heavy stands, e.g. 24886-TPA seedling stands). The dominant dense-stand residual. [[fvsjl-stand-pct-rdpsrt-fix]] residual.
- **Direct DGSCOR/volume ULP** — per-tree DG and volume Float32 op-order (gfortran transcendental order); the ~28%-pure-bit-exact
  growing-stand stratum is 100% bit-exact-or-DGSCOR-cornered. Stratum-uniform (no forest-type/site/age model gap).
- **Non-native cycle-length DGSCOR drift** (deferred, known) — each variant bit-exact at its native cycle (SN 5yr, NE/CS/LS 10yr);
  drifts ~3% at a non-native cycle length. Accepted residual.
- **COMPRESS eigensolver + within-class order** (s22) — IBM Jacobi eigensolver + RDPSRT partition transliterated bit-exact; the
  lone ~1% residual is a sub-Float32-ULP PC1/PC2 sort-key tie flipping a near-tied partition. Faithful port; no fix without bit-matching the eigensolver.

## C. FVS bugs (FVSjl is the CORRECT side and does NOT replicate them)
- **D38 R9-Clark short-tree SIGFPE** — r9ht/r9dib invalid-op/underflow on short trees; root-caused + fixed live (build-flag + source guard) for maintainer submission (`docs/patches/r9clark_D38_allsites.patch`).
- **CS essprt.f stump-sprout SIGFPE** — 1./((DSTMP/0.7788)−0.4403) div-by-zero at DSTMP≈0.343"; root-caused + fix proposed; in-container relink env-blocked.
- **NE NVEL VOLINIT extreme-height volume zeroing** (43do) — VOLINITNVB returns TVOL=0 for extreme height:diameter geometry (TopHt >~250 ft)
  despite its r9clark sub-call computing correct volume; live .sum vol=0, FVSjl reports the correct nonzero volume. Blast radius NARROW
  (rare extreme-height trigger; denser regen stands at realistic 45-76 ft don't hit it). [[fvsjl-ne-r9clark-domain-gap]].
- **Shared SDI overflow on degenerate ultra-dense micro-stands** (SN) — FVSjl reproduces it (faithful).

## D. Faithful-but-extreme behaviours (bit-exact vs live, not divergences — noted to avoid re-litigation)
- **NE height-model extrapolation runaway** — one dense-regen stand reaches TopHt 295 ft (vs SI 70); FVSjl reproduces bit-exactly every cycle. Faithful; it is the trigger for the C-class VOLINIT bug above.

## Status
Every class surfaced by the sweep is FIXED or CORNERED-to-a-named-primitive. No unexplained divergence remains among the
processed dig batches. The sweep continues to full-population completion; new dig batches are processed at each DIGCAP pause
against this taxonomy (a new class ⇒ a new both-sides-trace; a known class ⇒ cornered by fingerprint).

## LS growth-div candidate bucket — fully reconciled post-FIX-#8 (audit 43eb)
The 29 LS `REAL_growthdiv` candidates (the last needs-per-stand-verification frontier) are all explained:
- **1 REAL BUG FIXED** — FIX #8 (LS REGENT calibration stale-HTGR carry; exemplar 1831637837290487 was 2-3× off
  ⇒ now bit-exact-or-±1). The tamarack over-growth that FIX #7 exposed.
- **28 CORNERED primitives** — 11 resolved/ULP, 17 ultra-dense self-thinning RDPSRT tie-break (cycle-0 bit-exact,
  diverge only in later cycles on 12k–39k-TPA seedling stands), 1 AVHT40 top-height RDPSRT tie-break
  (55250794010661, cycle-0 tied-DBH aspen29/balsam22; the already-dug stand-dependent quicksort primitive).
Full per-stand classification: docs/fia_ls_candidates_classified.md. NO remaining unexplained LS divergence.
