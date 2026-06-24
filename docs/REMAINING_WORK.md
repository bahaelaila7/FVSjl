# FVSjl — remaining work and blockers

Status snapshot (suite 4223). The natural-process core + most management/disturbance keywords are
ported and validated bit-exact vs live Fortran. What remains, grouped by the **nature of the blocker**
(not just "todo") so the path for each is clear.

## A. Portable + validatable, but a multi-routine subsystem (a dedicated effort, not a quick win)
- **FFE surface-fuel dynamics → grown-cycle Stand Carbon Report.** The carbon report is already
  bit-exact at the **inventory cycle** (every column). Grown cycles need the coupled fuel subsystem:
  - `FMCWD` down-wood decay (`cwd *= (1−DKR)^NYRS` per size/decay class, hard→soft transfer) — equations
    + the SN DKR/PRDUFF constants are captured in `FFE_FUEL_DYNAMICS_chunk_plan.md`.
  - `FMCADD` additions — annual **litterfall** = Σ foliage `CROWNW(I,0)/LEAFLF(SP)` per live tree
    (the term that stops litter crashing at its 0.65/yr decay), woody crown breakage (`LIMBRK`), and
    snag-debris falldown (`CWD2B` pool).
  - a per-cycle FFE **driver** (call the fuel update each cycle, not only on a fire event), behind the
    `test_fire.jl` fire-path regression gate.
  - **Blocker:** none external — it is simply a coupled multi-routine port (~3 routines + LEAFLF/LIMBRK
    tables + the CWD2B accumulator + CROWNW wiring). Fully scoped; validatable against the committed
    `carbon_jenkins.report.save` 1995/2000 rows. This is the most-upstream genuinely-remaining item.

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
