# C7 — Fire & Fuels Extension (FFE) — chunk plan

The FFE is the largest SN extension (~100k Fortran LOC across ~50 `fm*.f` files; `fmphotoval.f`
alone is 152 KB). It is the fire-behavior / fuels / snag / carbon model exercised by **snt01 stand 4**
(the "FFE TEST"). Like ESUCKR, it is ported in validated chunks driven by the decision flow, not in one
pass. This plan scopes it from the FFE cycle driver `FMMAIN` (fmmain.f) and stand 4's keyword block.

## What snt01 stand 4 exercises (the validation target)

```
FMIn
  SNAGINIT / SNAGBRK        — initial snags + snag-break
  FLAMEADJ 2003             — flame-length adjustment
  SIMFIRE  2003 …           — a SIMULATED fire in 2003 (the .sum-affecting event)
  SALVAGE  2003             — post-fire salvage
  DEFULMOD / SNAGPSFT       — default fuel model + snag soft/hard fractions
  PotFIRE / POTFTEMP        — potential-fire reports
  SNAGOUT/SOILHEAT/BurnRept/FuelOut/FuelRept/MortRept — reports (DBS/list output, C6 territory)
End
```

The .sum-affecting core is **SIMFIRE → fire behavior (Rothermel) → FMEFF (fire-caused mortality +
top-kill)**; everything else is fuel/snag/carbon accounting + reports.

## FMMAIN call sequence (fmmain.f) → chunks

`FMCBA` (canopy bulk density / crown biomass) → `FMTRET`/`FMFMOV` → `FMUSRFM`/`FMCFMD` (fuel model) →
`FMBURN` (the fire: moisture → Rothermel spread → flame length → `FMEFF` mortality) → `FMSNAG` (snag
dynamics) → `FMCWD` (coarse woody debris) → `FMCADD` (carbon pools).

## Proposed chunks (foundational → fire)

- **F1 — tree biomass (Jenkins):** ✅ **DONE** — `jenkins_biomass` (fmcbio.f), the per-tree
  aboveground/merch/root biomass (tons) from DBH + species Jenkins group. Pure, CSV-driven
  (`fire_biomass.csv` = per-species `BIOGRP`), unit-tested vs hand-computed Fortran values. Foundation
  of the carbon pools (FMCADD) and a building block for crown/surface fuels.
- **F2 — crown biomass by component (FMCROWE/FMCBA):** per-tree foliage / branchwood by size class →
  canopy bulk density + canopy base height (drives crown-fire).
  - **F2-data:** ✅ **DONE** — the two per-species parameter blobs FMCROWE needs are extracted to
    `fire_species_props.csv`: `ls_spi` (ISPMAP, the SN→Lake-States species map that picks the crown
    equation set, fmcrow.f:148) and `v2t` (wood specific gravity lb/cuft = the `SG` arg, fmvinit.f).
    The same file also carries the snag decay/fall classes (`tfall_cls`/`dkr_cls`/`snag_cls`/`leaf_life`)
    that F3/F7 (snag dynamics) will need. 12 unit tests vs the Fortran source values; suite 3060→3072.
  - **F2-fn — PORTED (validation deferred):** `crown_biomass` (`src/engine/fire/crown_biomass.jl`)
    transcribes `FMCROWE` (fmcrowe.f, ~230 Julia ln): Jenkins total-aboveground → foliage/bark/wood/
    branch split → size-class allocation (red-oak / shortleaf-pine / maple / aspen proportion forms) +
    the unmerch bole-tip geometry (`UMBTW`, `LILPCE`) → `XV(0:5)`. Deps RESOLVED: `BRATIO` = FVSjl
    `bark_ratio` (verified identical Clark DIB=a+b·D + [0.80,0.99] clamp, bar Fort-Bragg overrides),
    `HTDBH` = `_htdbh_height`, `FMSVL2` = the SN cubic model via `_R8CLARK_VOL` (`_fm_cuft`), `SG`=`v2t`,
    `P2T`=0.0005. Faithfully reproduces FMCROWE's quirks (cone/frustum ×`SG/P2T`≈×2000, the
    sub-breast-height cylinder added RAW with no SG/P2T) ⇒ FFE-internal units, not literal tons.
    ⚠ **NOT end-to-end validated** — crown biomass has no `.sum` output; its magnitude can only be
    confirmed against live Fortran once F5/F6 (fire behavior) exist. Inert (not yet called); 27
    structural unit tests pin the component split / size-class ordering / species-form selection.
    Suite 3072→3099. Next: FMCBA aggregates XV → canopy bulk density + canopy base height.
- **F3 — FFE state + fuel pools (FMINIT/FMCBA):** the per-stand `FireState` (no globals): surface fuel
  loadings by size/decay class, snag arrays, the down-wood pools; SNAGINIT/DEFULMOD keyword setup.
  - **F3-data:** ✅ **DONE** — the FMCBA initial surface fuel-loading tables extracted to
    `fire_fuel_dead.csv` (FUINI, 9 forest types × 11 size classes: 9 down-wood + litter + duff) and
    `fire_fuel_live.csv` (FULIV, 4 live-fuel types × herb/shrub), loaded into `SpeciesCoefficients`
    (`ffe_fuel_dead`/`ffe_fuel_live`). The two pure forest-type→fuel-category maps (`ffe_dead_fuel_type`
    from FIA forest type, `ffe_live_fuel_type` from IFFEFT) + loading accessors ported in
    `src/engine/fire/fuel_loading.jl`. 17 unit tests vs the Fortran table values; suite 3099→3116.
  - **F3-state — REMAINING:** the per-stand `FireState` struct (STFUEL down-wood pools by size/decay
    class, snag arrays, cover type / percent cover), `FMSNFT` (FIA→IFFEFT FFE-forest-type map, also
    needed by F3's live fuels + F4), the decay-class split (fmcba.f:372+), the FULIV2 understory
    age/SI shrub curve, and the SNAGINIT/DEFULMOD/FUELINIT keyword setup. Wire FMCBA into the cycle.
- **F4 — fuel model classification (FMCFMD):** stand condition → fire-behavior fuel model (static +
  dynamic weighting). The Anderson/SB fuel-model loadings → CSV.
- **F5 — fire behavior (FMBURN core):** fuel moisture → Rothermel surface spread → flame length, with
  FLAMEADJ; the SIMFIRE trigger + fire-type (surface/passive/active crown) logic.
- **F6 — fire effects (FMEFF):** fire-caused mortality + crown/top-kill from flame length & bark
  thickness — the .sum-affecting kill that makes stand 4 diverge today. Wire into the mortality path.
- **F7 — snags + CWD + consumption (FMSNAG/FMCWD/FMCONS):** snag fall/decay, fuel consumption by the
  fire, down-wood transfer.
- **F8 — carbon pools (FMCADD) + reports:** the carbon accounting + the DBS/list reports
  (BurnRept/FuelRept/MortRept/CarbonReport) — overlaps C6 (DBS output).

## Validation

Each chunk: pure pieces unit-tested vs hand-computed Fortran; the integrated fire (F5/F6) validated by
making snt01 stand 4 (or a dedicated SIMFIRE scenario) bit-exact vs live Fortran `/tmp/FVSsn_new`, the
same harness used for ESUCKR (`fortran_baseline.sh` + `.sum` diff). Gated so the non-FFE default path
stays bit-exact until each chunk lands.

## C8 (insects / econ)

Econ (ANNUCST/HRVxxx, sn.key) and the insect models are separate, smaller extension families ported
after the FFE core. Econ is partly scoped already (see the FVS-econ memory); it affects the econ DBS
tables, not the `.sum`.
