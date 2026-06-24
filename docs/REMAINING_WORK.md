# FVSjl — remaining work and blockers

Status snapshot (suite 4223). The natural-process core + most management/disturbance keywords are
ported and validated bit-exact vs live Fortran. What remains, grouped by the **nature of the blocker**
(not just "todo") so the path for each is clear.

## A. FFE surface-fuel dynamics → grown-cycle Stand Carbon Report (7/9 COLUMNS BIT-EXACT)
Inventory cycle bit-exact (all columns). **Grown-cycle: 7 of 9 report columns now reconcile bit-exact**
vs the 4-cycle Fortran baseline (carbon_jenkins) — see test_carbon.jl. Ported + validated this session:
- ✅ `FMCWD` decay (fuel_decay.jl) · ✅ `FMCADD` litterfall + woody breakage (fuel_additions.jl) ·
  ✅ per-cycle driver `ffe_fuel_update!` · ✅ periodic-mortality → snags (mortality.jl) · ✅ dead coarse-
  root pool BIOROOT (Below-Dead column).
- ✅ THE BUG FIX: `crown_biomass` missed the V2T/=2000 rescale (fmvinit.f:1094) → woody cone weights
  ~2000× too large; root-caused via a Fortran XV dump (instrument fmcrow.f/fmcrowe.f, rebuild, dump,
  revert — the `DEBUG` keyword segfaults the stripped binary). One-line fix `sg = v2t * _FM_P2T`.
- ✅ RECONCILE: above/merch/below-live, **below-dead**, forest floor (every cycle), DDW (1990/95/2005),
  shrub — all bit-exact; Stand-Dead populated.
- ⛔ **REMAINING — CWD2B crown-debris scheduling (FMSCRO)**, fully specified in
  FFE_FUEL_DYNAMICS_chunk_plan.md (the `cwd2b[4,6,60]` state, the TFALL table, the at-death crown
  spread, the annual flow → DDW, and the snag = bole-only split). Closes the last two: DDW *timing*
  at the first post-mortality cycle (2000: 2.1 vs 3.8) + the exact Stand-Dead split (6.1 vs 5.2). ⚠ a
  coupled refactor that touches `add_snag!` (so the validated fire path via fmburn! too) — needs its
  own focused effort behind test_fire.jl, not an end-of-session rush.
- ⛔ then: `grow_cycle!` hot-path wiring of the driver + the `.out` Stand-Carbon-Report writer.
- **Blocker:** none external; the most-upstream remaining item is now down to FMSCRO + 2 integration steps.

## B. Validation-blocked by the stripped ground-truth binary
The rebuilt `/tmp/FVSsn_new` is a **stripped DBS build**: the DATABASE block accepts only
`DSNOUT/SUMMARY/COMPUTDB/TREELIDB`; every other DBS sub-keyword errors. So there is no ground-truth
SQLite to diff against (see `fvsjl-ground-truth-binary-limits`).
- **~13–14 DBS output tables** — `FVS_Carbon`, `FVS_StrClass`, `FVS_Fuels`, `FVS_PotFire`,
  `FVS_Mortality`, `FVS_Consumption`, `FVS_SnagSum`, `FVS_DWD`, `FVS_BurnReport`, `FVS_Hrv_Carbon`,
  `FVS_Summary2`, the 2 econ tables, … **Blocker:** need a fuller SN binary compiled with all DBS
  modules. NB the UNDERLYING data of several is validatable via the `.out` **text** reports (carbon via
  `CARBREPT` is done; `FUELREPT`/`POTFIRE`/`MORTREPT` text reports likely work too) — only the SQLite
  *table emission* is unverifiable here.

## C. `.sum`-inert + FFE-coupled (low value, hard to validate)
- **FINTM** (GROWTH mortality measurement period) — its only effect (cratet.f:486) scales INPUT
  dead-tree PROB into the FFE snag model. `.sum`-inert; needs the unported input-mortality→snag path
  (`FMSSEE`; FVSjl's `add_snag!` only fires from fire-kills); validatable only via FFE snag DBS (B).
- **FFE snag-input path (`FMSSEE`)** — seeds the initial snag list from input dead-tree records. Same
  snag-DBS validation blocker.

## D. Out-of-scope for the SN variant
- **8 insect/disease models** (MPB/DFB/DFTM/WSBW/BRUST/MISTOE/RDIN/ANIN/RRIN) — linked into SN as
  **`ex*.f` NO-OP stubs**; the real models (mpbmza/dfbmza/…) are **absent** from the SN build. Keyword-
  recognized but inert without infected/host input. **Blocker:** there is nothing to port for SN; a
  faithful port/validation would need a different variant binary + infected-stand scenarios.

## E. Long-tail accuracy (a residual, not a missing extension)
- **38-species single-species DG-calibration tail** (~1–2% TPA for non-snt01 species). snt01's species
  (22/27/33/65/89) are bit-exact; single-species stands exercise each species' calibration regression
  harder, exposing per-species COR/OLDRN residuals. **Blocker:** a long tail — each species needs a
  per-tree calib decomposition (ran/HTCALC/COR) vs oracle; partly transcendental-ulp. snt01 is the gate.
- **DGSCOR per-record cubic-volume ±0.03% drift** — documented as **irreducible** (transcendental-ulp
  through the bounded redraw), not a bug. See `DIVERGENCES.md`.
