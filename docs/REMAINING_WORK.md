# FVSjl — remaining work and blockers

Status snapshot (suite 4223). The natural-process core + most management/disturbance keywords are
ported and validated bit-exact vs live Fortran. What remains, grouped by the **nature of the blocker**
(not just "todo") so the path for each is clear.

## A. FFE surface-fuel dynamics → grown-cycle Stand Carbon Report (IN PROGRESS)
Inventory-cycle report is bit-exact (every column). Grown-cycle progress this session:
- ✅ `FMCWD` down-wood **decay** ported + unit-tested (`fire/fuel_decay.jl`, `fmcwd!`).
- ✅ `FMCADD` **litterfall** ported (`fire/fuel_additions.jl`) — the grown-cycle **FOREST FLOOR
  reconciles BIT-EXACT** (carbon_jenkins 1990→1995: 9.1→6.6 vs Fortran) via the annual fuel loop
  (FMCWD+FMCADD each year, NYRS=1, fmmain.f:226). Implicitly validated crown_biomass FOLIAGE.
- ⛔ **DDW (woody breakage) BLOCKED on FMCROWE woody-component units.** The FMCADD breakage logic is
  right, but `crown_biomass`'s WOODY components (`xv[2..6]`) come out in absurd magnitudes (8″ loblolly:
  8.7/545/14394/20484, not tons — the UMBTW bole-tip `×SG/P2T` cone weights). **Next concrete step:**
  instrument `fmcrowe.f` to dump XV for the carbon_jenkins trees, rebuild, diff vs `crown_biomass`, fix
  the FVSjl woody computation; then the breakage term drops in and DDW (2.1 decay-only vs Fortran 2.5)
  closes. This is the F5/F6 crown-biomass validation.
- ⛔ per-cycle **driver** (wire the annual loop into `grow_cycle!`, behind `test_fire.jl`) — after DDW.
- **Blocker:** none external; the floor half is bit-exact, the remaining sub-steps are the FMCROWE
  woody fix (needs a Fortran XV dump) + the driver. Most-upstream genuinely-remaining item.

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
