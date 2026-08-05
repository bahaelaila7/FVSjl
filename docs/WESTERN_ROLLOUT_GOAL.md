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

## Extensions matrix
- **FFE**: ALL western validated-cornered ✓ (+ eastern + CR).  **Dwarf mistletoe**: ALL western DONE ✓.
  **ECON**: DONE ✓.  **Climate-FVS**: TODO (largest remaining extension, lowest priority; inert w/o ready-file).
- **FIA sweep**: whole-cluster multi-cycle validated (2026-08-03); residuals = ZZRAN/DGSCOR dense-regen straddle.

## REMAINING WORK — drive each to bit-exact-or-cornered (task-tracker #142/#137/#140/#143)
1. **CI refinement tail [#142]**: DG under-grows QMD/BA = **backdated-density** (dgf must read growth-START
   BA/RELDEN — live dgf BA=66.80 vs jl current 85.13; SAME class as the CR 8th-bug just fixed) · volume MATW/FW2W
   tuning · un-defer ZZRAN small-tree · SMHTGF stochastic. Oracle /workspace/.ciwork/FVSci_clean. Audit docs/CI_VARIANT_PORT_AUDIT.md.
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
