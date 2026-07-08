# FIA/FVS behaviour-compatibility — working checklist / audit

Goal + doctrine: `docs/FIA_FVS_COMPAT_GOAL.md`. Every slice: plots covered, per-cycle pass rate vs
freshly-relinked live FVS, divergences found → both-sides-traced → fixed or cornered. Never regress the
floor (`julia --project=. test/runtests.jl` = 38527/143/0).

## Infra fixes
- **F1** — `test/harness/fia/validate_fia.jl`: fixed `FVSjl.NorthEast()` → `Northeast()` (would have
  errored EVERY NE FIA stand). Made the harness CLI-arg driven (`julia validate_fia.jl <listfile>
  <SN|NE|CS|LS>`) so it's reusable for the campaign. Committed.

## Slice 1 (Pillar 2 probe) — SN multi-cycle drift is REAL and must be characterized
First multi-cycle differential (3 SN plots from `/tmp/fia_val/sn_feas.txt`, NUMCYCLE 5, 6 stand cols vs
live FVSsn). Result — mean |rel diff| by cycle:

| cyc | TPA | BA | SDI | CCF | TopHt | QMD |
|----:|----:|---:|----:|----:|------:|----:|
| 0 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| 1 | 0.04% | 0.5% | 0.16% | 0.14% | 0.0% | 0.46% |
| 3 | 1.0% | 0.96% | 0.66% | 0.37% | 0.0% | 1.15% |
| 5 | 1.86% | 1.23% | 0.98% | 0.49% | 0.0% | 1.54% |

Worst stand 3.9% (3237541010661). **Cycle-0 bit-exact** (confirms the inventory reader); the drift is
**purely in the projection** and **grows with cycle**. **TopHt = 0.0% every cycle** ⇒ heights match; the
divergence is **DBH-driven** (BA/SDI/QMD directly, TPA via density-dependent mortality).

**Hypothesis (to prove/refute, Pillar 4):** the DGSCOR diameter serial-correlation / grown-Float32
accumulation tail — the SAME accepted class as the `test_multicycle` / `test_dbs_compute` MYBA/MYSDI
`@test_broken`s. BUT the curated snt01/net01 stands are bit-exact multi-cycle, so a ~2-4% drift on real
plots is either (a) that tail hitting harder on these species mixes / larger diverse stands, or (b) a
REAL FIA-specific gap (a species DG coefficient or calibration path these plots trigger that the curated
tests don't). Doctrine #3/#4: root-cause the worst offender BOTH-SIDES before cornering. Next slice:
per-cycle per-species trace of 3237541010661 (are its drivers the WK3-calibrated sp33/65 family = accepted,
or a clean species that SHOULD be bit-exact = real gap?), then scale the differential to all 4 variants.

### Slice 1b — worst-stand species composition (classification clue)
Both worst SN stands are **diverse mixed-species hardwood plots**: 3237541010661 = 5 species
(FIA 68/407/835/541/462), 3196569010661 = **12 species** (FIA 826/837/621/318/931/541/409/701/544/521/
491/404). This matches the documented `test_allspecies` verdict: a diverse many-species stand
"accumulates every species' sub-ULP per-cycle DBH-growth + tripling-spread residual into the nonlinear
density/volume sums = the ACCEPTED aggregate DGSCOR + tripling class." Real FIA plots are inherently
multi-species (unlike the clean curated single/few-species tests), so they AMPLIFY that same cornered
primitive — a strong hypothesis that the ~2-4% drift is the accepted class, not a new gap.

**DECISIVE next diagnostic (do NOT assume — doctrine #3/#4):** the `test_treeszcp` method — compare
PER-TREE DG jl-vs-live (via FVS_TreeList DBS or a debug-FVS dgdriv stamp) on one worst stand. If per-tree
DG is bit-exact to ~1 ULP and only the AGGREGATE sums drift, it IS the accepted grown-Float32/DGSCOR
accumulation class (corner it). If a specific species' per-tree DG diverges beyond ULP, it's a REAL gap
(fix it). Also worth checking: cycle LENGTH used (non-native-cycle DGSCOR is a separate documented
deferred residual — fvsjl-scenario-sweep-findings #2).

## TODO
- [ ] DECISIVE: per-tree DG jl-vs-live on 3196569010661 (accepted aggregate class vs real per-species gap).
- [ ] Scale the multi-cycle differential: larger SN sample + NE/CS/LS (Pillar 1 manifest feeds this).
- [ ] Build the stratified per-variant plot manifest (Pillar 1).
- [ ] Management-scenario differential on real plots (Pillar 3).
