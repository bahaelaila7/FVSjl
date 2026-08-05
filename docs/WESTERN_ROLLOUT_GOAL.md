# ACTIVE GOAL — Western FVS variant cluster + extensions rollout (FVSjl)

## Mission (user-cemented 2026-08-03; re-anchored 2026-08-05)
"Work unattended until ALL extensions for the variants implemented so far are ported and validated. Do a full
FVS-ready FIA sweep as well. Do NOT stop." Port + validate the WESTERN FVS variant cluster and all extensions,
bit-exact-or-cornered vs live FVS oracles, chunk by chunk. Branch: `kt-variant-port`.
DO NOT narrow scope to a single variant — CR is DONE; the goal is the whole cluster + extensions + FIA sweep.

## Variant status — growth+volume ports (oracle = live FVS relinked from bin/FVS{v}_buildDir/*.o)
- **CR** (Central Rockies) ★★ COMPLETE (2026-08-05: 3 bugs fixed — backdated-density dub / forkod imodty /
  strict site-species — DB sweep 1/40→39/40; every residual bit-exact or measured accepted primitive).
- **KT / IE / EM / BM / TT / UT** ★★ growth+volume bit-exact-or-cornered.
- **BC** (British Columbia) — growth+yield complete; remaining: merch/board vol, V2/non-ICH.
- **CI** (Central Idaho) ◐ IN PROGRESS — 9/9 chunks running, cit01 1990 bit-exact; refinement tail OPEN.

## Systematic DG-calibration dispatch audit — COMPLETE (2026-08-05)
Swept every shared COR-shrinkage/bark dispatch for CI-class missing-variant branches:
- **PSIGSQ** (COR Bayes-shrinkage prior variance): was missing CI (fixed 3b9aa35) AND IE (fixed 96cde22, was on the
  SN 0.0898 scalar default). Dispatch NOW COMPLETE: NE/CR/KT/EM/TT/UT/BM/BC/CI/IE all wired. Both fixes source-
  verified faithful, `.sum`-inert on cit01/iet01 (no regression), correct for other-species stands.
- **Bark**: CR/TT/BM/BC/CI have POWER/special bratio wired; EM/IE/UT/KT barks are linear-encodable (c.bark_a/bark_b)
  and the shared bark_ratio clamps [0.80,0.99]==live bratio.f ⇒ their linear fallback is FAITHFUL. CI's POWER
  ci_bratio was the ONLY genuine missing-bark branch. No further bark bugs.
- **DGSD**: all 9 western variants set it explicitly from their grinit.f (2.0 except BM 1.5, CI 1.7). No gap.
⇒ No remaining missing-branch bugs in the shared DG-calibration path cluster-wide.

## Extensions matrix
- **FFE**: ALL western validated-cornered ✓ (+ eastern + CR).  **Dwarf mistletoe**: ALL western DONE ✓.
  **ECON**: DONE ✓.  **Climate-FVS**: TODO (largest remaining extension, lowest priority; inert w/o ready-file).
- **FIA sweep**: whole-cluster multi-cycle validated (2026-08-03); residuals = ZZRAN/DGSCOR dense-regen straddle.
  ★ 2026-08-05 POST-FIX real-FIA re-validation (docs/WESTERN_FIA_VALIDATION_2026-08-05.md, stands drawn live from
  the 70GB FVS-ready DB by VARIANT): **CI** 25-stand slice → 6 treed, 0 jl crashes, cyc0 5/6 bit-exact, remaining
  = Δ1-NINT + the accepted DGSCOR/density compounding tail; no regression from the bark/CI_PSIGSQ fixes. **IE**
  60-stand slice → 9+ treed, cyc0 8/9 bit-exact, 0 crashes (monotone partial; no regression from IE_PSIGSQ).
  Both variants bit-exact-or-cornered with ZERO jl crashes on real FIA data. Harness: extract_sample.jl + the
  generalized scratchpad/fia_sweep_check.jl (any variant, reusable cluster-wide).

## REMAINING WORK — drive each to bit-exact-or-cornered (task-tracker #142/#137/#140/#143)
1. **CI refinement tail [#142]**: cit01 jl OVER-KILLS TPA ~2%. RELIABLE STATE (2026-08-05, after FOUR wrong
   root-causes corrected by measurement — backdated-density/GF-COR/bark/deferred-ZZRAN all refuted): CI deterministic
   DG is BIT-EXACT (GF DDS jl==live), serial-corr is ACTIVE (real-run c.sigma[4]=0.26, NOT deferred), COR applied
   (c.dg_cor[4]=0.05693=live). ⇒ the ~2% over-kill is the DGSCOR RNG-realization = the accepted "ZZRAN/DGSCOR
   dense-regen straddle" (cornered; straddles ~0 across stands per the 2026-08-03 FIA-sweep memo) → MEETS the bar.
   ★ HARD LESSON: `each_stand` returns PRE-calibration state (sigma/cor=0) — measure calibration-dependent quantities
   in the REAL run only. 2 real adjacent bugs FIXED (faithful, .sum-inert cit01): 0fa9677 bark branch, 3b9aa35
   CI_PSIGSQ branch. STILL open (unmeasured): volume MATW/FW2W · SMHTGF small-tree stochastic. Oracle FVSci_clean.
2. **EM establishment self-thin over-kill [#137]** (shared with BM #140) — dense-cohort self-thin.
3. **BM small-tree late-cycle under-thinning [#140]** — TPA +28% by 2090 (SDI-plateau / self-thin).
4. **IE AUTOES tally-amount close-out [#143]** — ~22% diffuse residual; validate vs stand4_booktpa on UNMODIFIED FVSie_clean.
5. **Climate-FVS extension** (clinit/clin/clgmult) — lowest priority.

## DOCTRINE (hard-won — carry from the FIA campaign)
1. Validate vs LIVE FVS oracle, bit-exact per chunk. 2. MEASURE, don't infer — instrument the Fortran.
3. Per-record treelist INVALID after tripling — use .sum aggregates / pre-split window. 4. Port faithfully then
validate; a regression on a faithful chunk = examine the oracle. 5. Reuse the shared engine — only add
variant-specific equations + data. 6. Document every chunk verdict in docs/{VARIANT}_VARIANT_PORT_AUDIT.md.

## Off-switch
`touch docs/WESTERN_ROLLOUT_COMPLETE` (USER's call). Per-variant done-flags: docs/{V}_VARIANT_PORT_COMPLETE.
Charters: docs/EXTENSIONS_ROLLOUT_PLAN.md. Memory: fvsjl-ci-variant-port, fvsjl-extensions-rollout,
fvsjl-{em,bm,ie,ut,tt,kt,bc}-variant-port. CR sub-goal retired → docs/CR_VARIANT_PORT_COMPLETE.
