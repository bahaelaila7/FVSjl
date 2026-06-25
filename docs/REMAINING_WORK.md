# FVSjl — remaining work and blockers

Status snapshot (suite **4392 pass + 10 @test_broken**). The natural-process core, the management /
disturbance keywords, AND the full FFE fire/fuels/carbon extension are ported and validated. What remains,
grouped by the **nature of the blocker** (not just "todo") so the path for each is clear.

## ✅ DONE this body of work (was the bulk of "remaining")
- **FFE Stand Carbon Report** — both CARBCALC methods (1 = Jenkins, 0 = FFE), the `.out` CARBREPT report
  integrated into the live `run_keyfile` path, live pools (Aboveground/Merch/Below-Live/Floor/Shrub)
  BIT-EXACT on the bit-exact-growth stand (carbon_snt).
- **Crown-lift** (FMSDIT/FMCADD) — the dominant post-mortality down-wood source — PORTED via per-record
  prev-state carried through DG tripling (`ffe_oldht/olddbh/oldcr` in `_TREE_VEC_FIELDS`); carbon_snt
  DDW final cycle bit-exact.
- **All 9 FFE DBS tables** emit in a live run: per-cycle `FVS_Carbon`/`FVS_Fuels`/`FVS_SnagSum`/
  `FVS_Down_Wood_Vol`/`FVS_Down_Wood_Cov`; fire-event `FVS_BurnReport`/`FVS_Mortality`/`FVS_Consumption`;
  potential-fire `FVS_PotFire` (FMPOFL: canopy bulk density + torching Monte-Carlo + dual-scenario surface
  fire); harvest `FVS_Hrv_Carbon` (FMSCUT/FMCHRVOUT FATE accrual + the FAPROP decay table, extracted to
  `data/southern/fire_hwp_fate.csv`). Each validated by a reconciliation/semantic test (the live-oracle
  diff is binary-blocked — see B — so they are validated by their value-grounding in already-bit-exact
  FFE pools + source-faithful semantics).
- **FFE input-snag seeding** (`ffe_seed_input_snags!`, FMSDIT→FMSADD ITYP=3) — input dead-tree records →
  inventory Stand-Dead/Below-Dead. (Supersedes the old "FMSSEE unported" note.)
- **YARDLOSS** — was a SILENT gap (active in snt01.key/sn.key, unrecognized → ignored). Ported: removed
  merch/saw/board scaled by (1−PRLOST); total cubic + TPA unchanged. Fuel-pool routing of the loss = C7.
- **Unrecognized-keyword surfacing** — `control.unrecognized_keywords`; a test asserts snt01.key/sn.key
  have zero (the guard that caught YARDLOSS).

## A. The one FFE-carbon residual (a timing lag, not a missing extension)
The carbon report's **dead pools** (Below-Dead / Stand-Dead / DDW) are NOT yet bit-exact at the
intermediate cycles — a ~0.5/0.7/1.2 t/ha residual tracked honestly as **10 `@test_broken`** (the
inventory + final cycles ARE bit-exact). ROOT CAUSE: the crown-lift is applied in the NEXT cycle's annual
fuel loop (with decay interleaving, magnitude-correct) while FVS applies it same-cycle. Both loop orders
were tested: the current order is the closer match (the "faithful" reorder overshoots the final cycle).
Closing it fully needs the same-cycle crown-lift + late-bucket mortality-crown scheduling in a single
post-grow annual loop — a coupled refactor of the shared `fmscro!`/death-dating, deferred as the last fine
increment. The input-snag bole also carries a small 3.92-vs-3.8 height-dub residual (heightless dead
records re-estimate height for the R8 Clark stem volume).

## B. Validation-blocked by the stripped ground-truth binary
The rebuilt `/tmp/FVSsn_new` is a **stripped DBS build**: the DATABASE block accepts only
`DSNOUT/SUMMARY/COMPUTDB/TREELIDB`; every FFE/carbon DBS sub-keyword errors. So there is no ground-truth
SQLite to diff the 9 FFE DBS tables against (see `fvsjl-ground-truth-binary-limits`). They are therefore
validated by **value-grounding** (every column is a reuse of an already-bit-exact FFE pool — e.g. DWD
volume = cwd biomass·2000/CWDDEN, fuel lt3+ge3 = the validated DDW) + **source-faithful semantics**, not a
live table diff. A fuller SN binary (compiled with all DBS modules) would let the SQLite emission itself be
diffed; until then the underlying values are the validation surface (and the `.out` CARBREPT text report IS
bit-exact-diffable, and is).

## C. `.sum`-inert + FFE-coupled (low value, hard to validate)
- **FINTM** (GROWTH mortality measurement period) — scales INPUT dead-tree PROB into the FFE snag model
  (cratet.f:486). `.sum`-inert; now reachable via `ffe_seed_input_snags!` but the exact FINTM scaling on
  the seeded snags is validatable only via the FFE snag DBS (B).
- **YARDLOSS fuel-pool routing (C7)** — the loss fraction routes to the FFE down-wood/snag/crown pools
  (DSNG/SSNG/CTCRWN). The `.sum` volume scaling is ported; the fuel-pool side is FFE-coupled (matters only
  with active FFE + a merch-removing thin, which the SN scenarios don't combine).

## D. Out-of-scope for the SN variant (verified, not gaps)
- **8 insect/disease models** (MPB/DFB/DFTM/WSBW/BRUST/MISTOE/RDIN/ANIN/RRIN) — linked into SN as
  **`ex*.f` NO-OP stubs**; the real models are **absent** from the SN build. Nothing to port for SN.
- **FMSNGHT snag height-loss** — `HTX=0` for every SN species (fmvinit.f:1089) ⇒ a no-op; the static snag
  bole is faithful.
- **MORTMSB / MATUREW** extra mortality — keyword-gated, default-inert (QMDMSB=999).
- **COMPRESS** PCA clustering — recognized; the eigensolver partition is not bit-identical (documented),
  the per-cycle aggregate is conserved.

## E. Long-tail accuracy (a residual, not a missing extension)
- **38-species single-species DG-calibration tail** (~1–2% TPA for non-snt01 species). snt01's species
  (22/27/33/65/89) are bit-exact; single-species stands exercise each species' calibration regression
  harder, exposing per-species COR/OLDRN residuals. A long tail — each species needs a per-tree calib
  decomposition vs oracle; partly transcendental-ulp. snt01 is the gate.
- **DGSCOR per-record cubic-volume ±0.03% drift** — documented as **irreducible** (transcendental-ulp
  through the bounded redraw), not a bug. See `DIVERGENCES.md`.
- The multi-cycle volume-sum `.sum` comparisons carry an accumulated Float32-vs-REAL*4 "ULP tail"
  (±1 truncation / ±N cuft on compounding sums) — floating-point in origin, documented in
  `TOLERANCE_AND_COVERAGE_AUDIT.md`.
