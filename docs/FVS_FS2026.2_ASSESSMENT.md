# FVS FS2026.1 → FS2026.2 change assessment (impact on FVSjl)

Fetched tag `FS2026.2` (USDAForestService/ForestVegetationSimulator) and diffed against `FS2026.1` (the release
FVSjl is validated against). **Conclusion: FS2026.2 requires NO change to FVSjl within its supported scope
(SN/NE/CS/LS default projection). FVSjl remains a faithful drop-in for FS2026.2 exactly as for FS2026.1.**

## What changed (52 files, +2244/-1888)
- **Test goldens (`tests/FVS{ie,kt,nc,oc,op,pn,so,tt,ut,wc,ws}/*.save`)** — WESTERN variants only. **NONE for
  FVSsn/ne/cs/ls** ⇒ the four FVSjl variants' regression output is byte-identical between FS2026.1 and FS2026.2.
- **#125 (`ci/regent.f`, `em/regent.f`, `ie/regent.f`)** — CR/UT height-growth handling + small-tree ht multiplier.
  WESTERN variants (Central Idaho / E. Montana / Inland Empire). Not SN/NE/CS/LS; not ported by FVSjl.
- **#124 "Fiavbc nfs fix" + #126 "Fiavbc warning adjustment"** (`vbase/setcubicdflts.f`, `vbase/initre.f`,
  `vvolume/fvsvol.f`, `base/errgro.f`, `estb/esout.f`, `vdbsqlite/*`, …) — corrections to the **FIAVBC** (FIA
  National Volume-Biomass Library) merchantability defaults + an ecoregion warning. The commit notes "Impacts LS,
  NE, and PN variants" but the `initre.f` change is gated `IF(LFIANVB)` — i.e. ONLY when the FIAVBC keyword is
  requested. **FVSjl does NOT implement the FIAVBC/NVB path** (southern.jl:28 — "FVSjl has only the R8 Clark
  equations; SN default LFIANVB=.FALSE.; a variant on the NVB path would diverge"). FIAVBC is a recognized-but-
  UNSUPPORTED keyword. So this change is on a path outside FVSjl's implemented scope.
- **`common/INCLUDESVN.F77`** — version string `FS2026.1` → `FS2026.2` (cosmetic). FVSjl hardcodes no FVS version.

## Not changed (verified)
- **`r9clark.f` UNCHANGED** ⇒ the D38 `r9ht` short-tree SIGFPE bug PERSISTS in FS2026.2 (the D38 patch in
  `docs/patches/` still stands).
- **No `sn/ ne/ cs/ ls/` variant source changed**; **no FVSsn/ne/cs/ls test golden changed.**

## Action for FVSjl
- **None required now** — FS2026.2's behavioral changes are FIAVBC (unsupported), western variants (not ported),
  and test/version. The eastern-variant default behavior FVSjl validates against is unchanged.
- **Future (only if FVSjl ever implements FIAVBC/NVB):** port the #124/#126 FIAVBC merch-default corrections
  (DBHMIN/TOPD/STMP/SCF* defaults in initre.f/setcubicdflts.f, LS/NE/PN).
