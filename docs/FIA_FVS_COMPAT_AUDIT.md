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

### Slice 1c — DEFAULT (drop-in) trajectory pinpoints a SYSTEMATIC diameter-growth drift (stand 3196569010661)
Ruled out cheap causes first (doctrine #3): **cycle length = 5-yr = SN native** (1998→2023, NOT the
non-native-cycle DGSCOR issue); **NOTRIPLE does NOT converge** them (jl 408 vs live 415 TPA @2003 — a
NOTRIPLE-path confound, not the default behaviour). The apples-to-apples DEFAULT (tripled) comparison:

| year | jl TPA/BA/SDI | live TPA/BA/SDI | BA Δ |
|-----:|---------------|-----------------|-----:|
| 1998 | 416/102/195 | 416/102/195 | 0 (bit-exact cyc0) |
| 2003 | 408/114/214 | 408/115/215 | 1 |
| 2008 | 396/125/231 | 394/128/234 | 3 |
| 2013 | 378/135/244 | 375/139/248 | 4 |
| 2018 | 362/146/259 | 358/151/264 | 5 |
| 2023 | 347/157/273 | 343/163/278 | 6 (3.7%) |

**Signature:** TPA tracks closely (Δ0–4); **BA drifts SYSTEMATICALLY and ONE-DIRECTIONALLY — jl BA always
LOWER**, 0→6 (3.7%) over 5 cycles; TopHt exact all cycles. ⇒ seed = **jl diameter growth marginally but
systematically below live's**, compounding into BA/SDI (and slightly higher TPA via reduced density
mortality).

**Verdict discipline (doctrine #3/#4):** this is **NOT yet cornered**. It is (a) directional (pure Float32
ULP is not) and (b) ~3.7% BA — orders of magnitude larger than the accepted `dbs_compute` MYBA/MYSDI
DGSCOR ULP (~4 ULP / 1.2e-4). A systematic directional diameter-growth bias on diverse real SN plots must
be root-caused per-tree before any corner. It could still resolve to the accepted class (if it's the
sp33/65 WK3-calibrated species' serial-correlation compounding one way on these mixes) OR be a REAL small
DG gap a curated few-species test never exercised. **This is the campaign's first substantive open divergence.**

## TODO
- [ ] DECISIVE next: per-tree DG jl-vs-live on 3196569010661 via FVS_TreeList DBS (live) + jl treelist —
      which species' per-tree DBH increment diverges, and by how much (ULP vs systematic). Classifies
      accepted-DGSCOR vs real-gap and, if real, names the species/coefficient/path.
- [ ] Scale the differential to NE/CS/LS + larger SN sample once the SN driver is understood.
- [ ] Build the stratified per-variant plot manifest (Pillar 1).
- [ ] Scale the multi-cycle differential: larger SN sample + NE/CS/LS (Pillar 1 manifest feeds this).
- [ ] Build the stratified per-variant plot manifest (Pillar 1).
- [ ] Management-scenario differential on real plots (Pillar 3).
