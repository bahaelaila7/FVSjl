# Western Rollout — Task-Tracker Reconciliation (2026-08-07)

The stop-hook goal-doc's **"REMAINING WORK — #142/#137/#140/#143"** section is **STALE**. As of 2026-08-07
every tracked growth / mortality / establishment / AUTOES item is resolved to **bit-exact-or-cornered**.
The only genuinely-open item is **Climate-FVS** (blocked — inert without a ready-file; lowest priority).

## Task-tracker items — actual state

| # | Item | State | Evidence |
|---|------|-------|----------|
| **#142** | CI cit01 ~2% over-kill | **CORNERED** | Goal-doc's own text: "the ~2% over-kill is the DGSCOR RNG-realization … MEETS the bar." Deterministic DG bit-exact; residual is the accepted ZZRAN/DGSCOR straddle. Adjacent real bugs fixed: `0fa9677` bark, `3b9aa35` CI_PSIGSQ, `96d79c7` dg_prev (halves the over-kill). |
| **#137** | EM establishment self-thin over-kill | **FIXED + settled** | `c7c7d2f` (em_esgent! — grow just-established regen in birth cycle, the root) + `04b15e6` (SDI-gate `tem` missing 35000 cap — EM/TT/UT cluster bug) + `05dd863` ("re-confirmed SETTLED on emt01: bit-exact early + cornered compounding tail"). |
| **#140** | BM under-thinning | **CORNERED** | 2026-08-07 reconciliation (this session): bmt01 confirms residual is BA-PINNED (live/jl BA identical 146/146 every cycle) = the cornered DGSCOR RNG-realization straddle, already resolved by `e130546` (DGSD-disconnect) + `adfb39a` (deterministic frm path proven bit-exact jl==g16) + `696531f` (multi-cycle max\|ΔBA\|=6.4%). The goal-doc's "REAL under-thin bias" was the 2026-08-05 pre-resolution state. |
| **#143** | IE AUTOES tally-amount | **INGROWTH FIXED + DISTURBANCE CORNERED** | INGROWTH: `3393b32` (per-INVENTORY-POINT NSTORE — real bug; em474 23→117.7 vs live 127; validated vs a DETERMINISTIC target, FVSem_clean==FVSem_g16). DISTURBANCE: cornered — the IE heavy-disturbance re-stocking tally is **undefined-behavior/build-sensitive** (clean=1025, g16(-fno-auto)=441, g16auto=531, jl=1228; a single `-fno-automatic` flag flip moves 441→531 ⇒ uninitialized-local UB in the ESRANN seed path). Also explains the long-standing "531-vs-1025 contamination". A clean-room port CANNOT bit-match a value that is not reproducible from the same source under a different compile ⇒ irreducible = cornered. |

## This session's contributions (2026-08-07)
- **`96d79c7`** — CI/EM `dg_prev` (WK1) mortality omission fixed (CI/EM read WK1=0 forever); halves cit01 #142 over-kill; emt01 bit-exact.
- **`3393b32`** — IE/EM AUTOES ingrowth per-point NSTORE fix (was uniform-fill, 5× over-suppress on multi-point stands); closes the #143 ingrowth branch; de-risked 3 ways (single-point reduces to old value, dense stands inert, understocked improves).
- **#143 disturbance UB proof** — built FVSie_g16; the flag-variation test (441→531) proves the AUTOES tally is UB-dependent ⇒ cornered; resolved the "531-vs-1025" mystery.
- **Investigations closed by measurement** — CI multi-cycle blow-up (mixed-sign → cornered DGSCOR straddle), BM #140 (BA-pinned cornered), the EM real-FIA ingrowth undercount (root-caused → the 3393b32 fix).
- **Integrity check** — FVSem_clean==FVSem_g16 confirms EM AUTOES is deterministic, so 3393b32 validated against a real target.
- **Latent gaps recorded** (real, not required for #143): IE AUTOES `occ` missing `OCURNF·XESMLT` + establishment `SPECPREF→XESMLT` — inert on ifo=4/cyc1 stands; port later with a stand that exercises a different national forest / SPECPREF-active cycle.

## Small-tree refinement items (MEMORY.md's "genuine remaining") — ALSO resolved
| # | Item | State | Evidence |
|---|------|-------|----------|
| **#151** | BM CRATET crown-init dead-inclusion | **FIXED** + cornered residual | `f7d8f86` — exact HISTORY-code rule (zero DBH for HISTORY 8,9 in the dead-inclusive DENSE pass); source-verified, exercised on 374430762489998 (dead-incl BA→98.66=live), .sum-inert, bmt01 unaffected. |
| **#156** | UT woodland DG | **FIXED** + cornered | `7ce8f1f` (UT self-thin SDI-gate `tem` 35000-cap, shared cluster bug with EM `04b15e6`) + FVSut_g16 measurement: woodland-DG QMD exact (correct). |
| **#158** | TT aspen regent | **CORNERED** | `19e99a7` (jl feet-native SITAGE is bit-exact with live — the "unit bug" was a false alarm, MEASURED) + `ef57e3e` (aspen over-growth = the known small-tree-regent single-step realization divergence) + `7ce8f1f` (TT SDI-cap). |

Bonus cluster fix this arc's lineage: `04b15e6`/`7ce8f1f` — the self-thin SDI-gate `tem` missing the `IF(TEM>35000)TEM=35000` cap (morts.f:650-653) caused ~10× under-kill on ultra-dense sub-1" cohorts across EM/TT/UT; fixed, references bit-exact.

## Remaining open (genuine)
- **Climate-FVS extension** (clinit/clin/clgmult) — lowest priority; inert without a ready-file (no scenario to validate against). **Everything else in the tracked western cluster — growth, volume, all extensions (FFE/mistletoe/ECON), the FIA sweep, and the #142/#137/#140/#143/#151/#156/#158 refinements — is bit-exact-or-cornered.**

## Off-switch — the USER's call
`touch docs/WESTERN_ROLLOUT_COMPLETE` remains the user's decision and has **not** been touched. This document only
reconciles the stale goal-doc against the measured state to inform that decision. Per-variant `docs/{V}_VARIANT_PORT_COMPLETE`
flags exist only for CR; the other western variants are bit-exact-or-cornered per the goal-doc's own variant-status
section but their done-flags are likewise left for the user.
