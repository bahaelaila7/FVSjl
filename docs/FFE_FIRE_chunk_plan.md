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
  - **F3-classifier:** ✅ **DONE** — `ffe_forest_type` (FMSNFT, fmsnft.f) ported: FIA forest type →
    IFFEFT 1–9 (hardwood / hardwood-pine / pine-hardwood / pine / pine-bluestem / oak-savannah /
    redcedar / St-Francis / nonstocked). The mixed classes split on the stand's pine BA fraction (SN
    species 4–14), the savannah/bluestem classes use top height (`stand_top_height` = AVHT40 = ATAVH)
    and `stocking_class` (ISTCL). 8 unit tests (incl. snt01 ft 520→1, pine-fraction splits); suite
    3116→3124. Drives F3's live fuels and F4's fuel model.
  - **F3-state:** ✅ **DONE (fuel/cover body)** — `FireState` expanded (cover type, percent cover, big
    DBH, live fuels, the `cwd[11 size × 2 dead/soft × 4 decay]` down-wood pools) and `fmcba!`
    (`src/engine/fire/fmcba.jl`) ports FMCBA's deterministic body: cover type = max-BA species, percent
    cover from per-tree crown areas (crown width computed on the fly, CWCALC iwho=0), live fuels by FFE
    forest type, and the first-year dead-fuel load split into decay classes by species BA share
    (`DKRCLS`). Gated on `fire.active` (no-op otherwise). 17 unit tests incl. fuel-pool conservation
    (Σ decay = FUINI) and cover-type/live-fuel correctness; suite 3124→3141.
  - **F3-rest — REMAINING:** the snag arrays + FULIV2 understory age/SI shrub curve, and the FFE
    activation path — parse the **FMIN** keyword block (SNAGINIT/DEFULMOD/FUELINIT/…) to set
    `fire.active` and wire `fmcba!` into the cycle. (Until FMIN parses, FFE stays inactive ⇒ all
    scenarios bit-exact; `fmcba!` is validated standalone.)
- **F4 — fuel model classification (FMCFMD):** stand condition → fire-behavior fuel model (static +
  dynamic weighting). The Anderson/SB fuel-model loadings → CSV.
- **F5 — fire behavior (FMBURN core):** fuel moisture → Rothermel surface spread → flame length, with
  FLAMEADJ; the SIMFIRE trigger + fire-type (surface/passive/active crown) logic.
  - **F5-core:** ✅ **DONE (the Rothermel model)** — `rothermel_surface_fire`
    (`src/engine/fire/rothermel.jl`) transcribes FMFINT's single-fuel-model body: size-class area
    weighting, live moisture of extinction, moisture/mineral damping, optimum reaction velocity → reaction
    intensity, propagating flux, wind & slope factors, spread rate, Byram intensity, and Byram→flame.
    Constants RHOP=32 / TMIN=0.0555 / SILFRE=0.01 per fmfint.f. Produces the Byram intensity the F6 chain
    consumes. 14 tests: physical invariants (wind/slope ↑ spread, moisture ↓ spread, too-moist→0),
    Byram→flame relation, grass-vs-timber, determinism. ⚠ exact bit-validation needs the Fortran oracle
    (no standalone `.sum`); confirmed at the integrated-fire step. Suite 3194→3208.
  - **F5b-inputs:** ✅ **DONE** — `fuel_moisture` (FMMOIS preset dead/live moistures by dryness model
    1–4, fmmois.f) and `fire_wind_reduction` (canopy wind multiplier WMULT = interp of CORFAC over
    CANCLS, fmburn.f:390/fmvinit.f) in `src/engine/fire/fuel_moisture.jl` — the two environmental inputs
    Rothermel needs from the fire weather + stand canopy. 12 tests incl. a full
    `fuel_moisture → wind_reduction → rothermel → scorch → CSV → mortality` chain on realistic inputs;
    suite 3208→3220.
  - **F5b-fuelmodel:** ✅ **DONE** — `build_dynamic_fuel_model` (`src/engine/fire/fuel_model.jl`)
    ports FMCFMD3: the SN custom fuel model assembled from the stand's own fuels — dead loads from the
    `fire.cwd` down-wood pools (1-hr = 0–.25"+litter), live-woody from the understory **crown biomass**
    (`crown_biomass`, foliage+½·fine for trees ≤ CANMHT) + live shrub, live-herb from `fire.flive`, with
    the fmgfmv dead-herb moisture split, USAV/UBD/CANMHT defaults, and the load-weighted depth + moisture
    of extinction. **This makes F2 crown biomass a live input** (no longer inert). 7 tests incl. the full
    `fmcba! → build_dynamic_fuel_model → rothermel_surface_fire` integration producing a fire; suite
    3220→3227.
  - **F5b-driver:** ✅ **DONE (the fire driver)** — `fmburn!` (`src/engine/fire/fmburn.jl`) runs one
    simulated fire: FMCBA → dynamic fuel model → FMMOIS/wind → Rothermel → Van Wagner scorch → per-tree
    RANN draw vs PSBURN → `fire_tree_mortality` + adjust → `CURKIL = PMORT·TPA (+ crown share)` →
    `tpa -= CURKIL` (fmeff.f:542-548). Returns a `FireResult` (killed TPA, flame, Byram, scorch). 14
    tests: kills trees size-dependently (saplings die, big oak survives, ≤1" outright), honors PSBURN /
    mortcode / FFE-off, deterministic, wetter fuel ⇒ fewer kills. Suite 3236→3250. **The whole fire now
    runs and applies the `.sum` kill.**
  - **F5b-keyword:** ✅ **DONE** — `kw_fmin!` parses the FMIN block (SIMFIRE date/wind/FMOIS/temp/
    mortcode/%-burned/season + FLAMEADJ flame-mult/crown), sets `fire.active`, and stores the event in
    `FireState`; `grow_cycle!` fires `fmburn!` once on the SIMFIRE year (before growth). **snt01 stand 4's
    SIMFIRE 2003 now actually burns from the .key file** (468→223 TPA at the right year). The whole FFE
    fire path is live; stands without FMIN are untouched (suite stays 3250/3250 green).
  - **F5b-validate — REMAINING (tuning):** make snt01 stand 4 bit-exact vs live Fortran. Two gaps
    observed: (a) a pre-existing bare-`THINDBH 3.` thinning isn't applied (Fortran thins 245 TPA at 1993,
    FVSjl doesn't) — a thinning-keyword issue, not fire; (b) the fire-mortality precision (the RANN-stream
    order vs growth/mortality, the dynamic fuel-model details, DEFULMOD overrides). Both are now
    observable since the fire runs.
  **All the FFE physics (F1–F6) is now ported**; F5b is the remaining integration/wiring + keyword layer.
  Scoped: `FMFINT` (fmfint.f, ~520 ln) is the Rothermel core — flame `= 0.45·(BYRAMT/60)^0.46`,
  `BYRAMT = XIR·R·384/SIGMA`; it loops the (up to MXFMOD=5) fuel models from `FMCFMD`, each characterized
  by `FMGFMV` (the fuel-model database: `SURFVL` SAV, `FMLOAD` loads, `FMDEP` depth, `MOISEX` moisture of
  extinction — assembled at runtime, not a simple DATA block). The Rothermel computation: size-class
  weighting (area `A`, weighting `F`), reaction intensity `XIR` (optimum reaction velocity `GAMMA`,
  moisture/mineral damping), propagating flux `XIO`, wind `PHIW`/slope `PHIS` factors, heat sink
  `RHOBQIG`, spread rate `R`. This is the hardest single FFE chunk — it needs a dedicated, careful pass
  (the fuel-model DB extraction + the full Rothermel transcription), validated by `BYRAMT`/`FLAME` for
  fm 10 against live Fortran, then feeding F6's now-complete `byram → scorch → CSV → mortality` chain.
- **F6 — fire effects (FMEFF):** fire-caused mortality + crown/top-kill from flame length & bark
  thickness — the .sum-affecting kill that makes stand 4 diverge today.
  - **F6-mort:** ✅ **DONE (the mortality equation)** — `fire_tree_mortality` + `fire_bark_thickness`
    (`src/engine/fire/fire_effects.jl`) port FMEFF's per-tree kill and FMBRKT: SN oaks/hickory/red-
    maple/black-gum use the Regelbrugge-Smith DBH + char-height (0.7·flame) logistic (MORTB0/1/2 by the
    6 mortality groups), other species the Reinhardt bark-thickness + crown-scorch logistic. Bark
    `B1[EQNUM]` (fmbrkt.f) → `bark_eqnum` column in `fire_species_props.csv` (+ the sp-5 Harmon
    quadratic). 41 unit tests vs a hand-transcribed reference (exact per-group values + monotonicity in
    flame/DBH); suite 3141→3182. Self-contained given flame length + crown-volume-scorched.
  - **F6-scorch:** ✅ **DONE** — `scorch_height` (Van Wagner, fmburn.f:470) and
    `crown_volume_scorched` (CSV from scorch height + crown geometry, fmeff.f:170-186) ported and
    tested. The full per-tree fire-effects chain is now closed **given a Byram fireline intensity**:
    `byram → scorch_height → crown_volume_scorched → fire_tree_mortality`. 12 more unit tests (Van
    Wagner exact values, CSV crown-geometry edge cases, the byram→mortality chain); suite 3182→3194.
  - **F6-wire — REMAINING:** the per-tree burn loop (FMEFF body): RANN draw vs %-stand-burned (PSBURN),
    then apply `PMORT` to tree TPA in FMBURN (the actual `.sum` kill). The only missing input is the
    Byram intensity / flame length + PSBURN from **F5** — once F5 lands, the whole chain is bit-exact
    testable and snt01 stand 4 validates end-to-end.
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
