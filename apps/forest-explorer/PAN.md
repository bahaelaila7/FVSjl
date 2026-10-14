# Forest Explorer × Pan (branch `forest-explorer-pan`)

This fork adds **Pan** — LANDIS-II-style succession models (Biomass Succession, PnET-Succession) and the HBAD
height/basal-area/density model, calibrated against FIA — as simulation engines next to **FVSjl**, on the same TreeMap
base. It also adds LANDIS Biomass-Harvest-style harvesting, climate-change scenarios for PnET, uploadable parameter files
and export/import of whole sessions. The FVSjl engine works as before.

## Run
```sh
cd apps/forest-explorer
./serve_pan.sh            # http://127.0.0.1:8088   (HOST/PORT/JULIA_NUM_THREADS override)
```
The script is **relocatable**: every path derives from the workspace root `WORKSPACE` (default: the folder that holds
`fe_pan/`, three levels above the script), so the same command works in the container (`/workspace`) and on a host
checkout (e.g. `~/Documents/projects`). Any single path can still be overridden. On first start (and whenever a
`Manifest.toml` changes) it installs the Julia packages of this app and of the Pan checkout into their private depots
(several minutes; `SETUP=1` forces, `SETUP=0` skips). Missing data paths are reported as warnings at start.
| variable | default | meaning |
|---|---|---|
| `WORKSPACE` | folder containing `fe_pan/` | root the defaults below derive from |
| `JULIA_DEPOT_PATH` | `./.julia_depot_fe` (+ `$WORKSPACE/FVSclaude/.julia_depot` if present) | depot of this app |
| `TREEMAP_DATA` | `$WORKSPACE/FVSclaude/treemap` | TreeMap raster/VAT/tree table (read only) |
| `FIA_DB` | `$WORKSPACE/FVSclaude/SQLite_FIADB_ENTIRE.db` | FVSjl master FIA DB (read only) |
| `FE_CACHE_DIR` | `./cache/sim_cache` | FVSjl sub-DB cache (kept out of the TreeMap data dir) |
| `PAN_REPO` | `$WORKSPACE/landis_scen` | Pan checkout (branch `scenario-engine`) — the Pan worker runs from here |
| `PAN_DEPOT_PATH` | `$PAN_REPO/.julia_depot_wt` (+ `$WORKSPACE/landis_parametrization_julia/.julia_depot` if present) | depot of the Pan worker |
| `PAN_FIA_DB` | `$WORKSPACE/FIA/FIASQLITE2PGSQL/FIADB.duckdb` | curated FIA cohorts for Pan (read only) |
| `PAN_CLIMATE_DB` | `$WORKSPACE/pan_climate/pan_climate.duckdb` | PRISM baseline + forecasts (PnET climate) |
| `PAN_SOIL_CSV`, `PAN_PNET_CLIMATE`, `PAN_CO2_CSV` | `$WORKSPACE/landis_region4/data_r4/…` | PnET plot soils, L4 climate, CO₂ |
| `PAN_THREADS` | `4` | threads of the Pan worker |
| `PAN_AUTOSTART` | `1` | start (and compile) the Pan worker in the background at server start (~1–2 min first time) |
| `PAN_WORKER_LOG` | `./cache/pan_worker.log` | the worker's log |

The Pan worker is a separate Julia process (Pan's own project and depot) talking JSON lines over stdin/stdout; it is
restarted automatically if it dies, and uploaded parameter files are re-sent to it.

## Using it
1. **AOI**: draw or upload polygons as before (samples: `samples/aois/*.geojson`, e.g. Oconee NF GA).
2. **Engine** (plan panel): FVSjl, *Pan · Biomass Succession*, *Pan · PnET-Succession*, *Pan · HBAD*. Each scenario keeps its
   own engine, parameter set and climate.
3. **Parameter set** (Pan engines): the default is the engine's **EPA Region 4** set (AL FL GA KY MS NC SC TN, 14 L3
   ecoregions; in-sample 46–49 of 50 species groups within ±10 %), unless you upload a YAML file
   (*Upload YAML…*) or pick a sample. The file's header (`pan_params: {plugin: …}`) decides which engine it belongs to; a
   mismatching file is refused with an offer to switch engine. Format and full roster: Pan `docs/scenario_engine.md` §1.
   Samples: `r4_{biomass,pnet,hbad}_v1.yaml` (defaults), `gafl_{biomass,pnet,hbad}_v1.yaml` (GA+FL-specific, cross-validated)
   and `gafl_{hbad,biomass,pnet}_alt{1..4}.yaml` — **distinct alternatives** that fit the data about as well as v1 but differ
   in parameters and projections; running a plan with each brackets the projection uncertainty the data cannot resolve
   (see "Alternatives" below).
4. **Plan**: `+ Thin / harvest` (FVS-style residual BA/TPA/SDI, from below/above, DBH range, species — Pan maps it to a
   harvest), `+ Plant`, and for Pan engines `+ Harvest` (LANDIS Biomass Harvest style):
   - cohort selection by species, age range and DBH range;
   - removal: % of each selected cohort's biomass, clearcut, or a residual BA/TPA/SDI target removed from below, from above
     or proportionally;
   - stand selection: all, or ranked (highest biomass, oldest, random) until a % of the AOI area is treated;
   - repeat every N cycles; optional planting right after (species, trees/ac, survival %).
   BA/TPA/SDI/DBH specifications act on Biomass Succession and PnET cohorts through the empirical
   (species, age, biomass) ↔ (TPA, QMD, BA, height) translation; HBAD tracks those variables directly.
5. **Climate…** (PnET only): baseline (PRISM forecast extrapolation, repeat the last T years, or average of the last K
   windows of T years), CO₂ pathway (constant, linear, SSP2-4.5, SSP5-8.5), and any number of changes (variable: mean / min /
   max temperature or precipitation; add °C or % change; step or linear ramp between years; optional months). The modal
   shows AOI-weighted annual charts and a monthly waveform of baseline vs modified; *Reset* restores the default.
6. **Run scenario**, compare scenarios, report — as before. Engine warnings (plots outside the parameter set's ecoregions,
   approximations) appear under the plan and in the status line; Pan volumes are labelled "(approx.)".
7. **Session bundle** (sidebar): *Export* writes a JSON with the AOI patches, every scenario's engine, embedded parameter
   YAML, climate settings and plan; *Import…* restores them (scenarios come back un-run). Sample:
   `samples/bundles/ga_fl_demo_bundle.json` (Oconee AOI, 6 scenarios: FVSjl, HBAD grow/thin+clearcut+replant, Biomass
   hardwood selection, PnET baseline vs hot & dry).

## What Pan reports
Per cycle (start-of-cycle state, before that cycle's management; removals on the row of the cycle they happen in):
BA ft²/ac, TPA, QMD, total cubic volume and board feet (**approximations**: BA × H × 0.42; 4 bd ft per ft³ for cohorts ≥ 9 in),
live aboveground carbon (0.5 × AGB), total carbon (adds roots at 0.25 × AG; no dead pools), removed volume and carbon, and the
species × 2-inch DBH tally (each cohort's QMD stands for its trees). Plots: TreeMap pixel → FIA plot CN → the curated FIA
cohorts of that visit. Coverage: GA+FL ecoregions (8.3.4, 8.3.5, 8.5.3, 8.4.1, 8.4.4, 8.4.9, 15.4.1); other plots are
reported as `out_of_region` (use FVSjl there). With the default Region 4 sets this covers the 8 Region-4 states; e.g. the
NC sample AOI simulates 163 of 208 plots (the rest have no curated FIA cohort data for that visit or lie outside the 14 L3s).

## Design decisions (this fork)
- **DD-P1 Separate Pan process.** FE and Pan have independent Julia projects (different package versions); a worker
  process isolates them, survives crashes (auto-restart) and keeps Pan's compiled code warm.
- **DD-P2 Same result seam.** Pan results are converted to `FiaSim.PlotResult`, so aggregation, maps, charts, compare and
  report are shared and engines are comparable side by side.
- **DD-P3 Plan JSON passes through.** Pan receives the frontend plan as-is; FVSjl ignores the Pan-only `harvest` kind (the UI
  does not offer it for FVSjl). FE `thin` actions are mapped to Pan harvests with the same semantics.
- **DD-P4 Parameter files are YAML with a plugin header** (owner's choice); files stay client-side so bundles can embed them.
- **DD-P5 Out-of-region plots are not silently mapped.** Plots outside the parameter set's ecoregions are skipped and
  reported, rather than borrowing another ecoregion's parameters.
- **DD-P6 Climate only for PnET.** Biomass Succession and HBAD are not climate-driven; the climate button is shown for PnET.
- **DD-P7 Caches outside the data dir.** `FE_CACHE_DIR` keeps FVSjl's sub-DB cache inside the app folder (the TreeMap
  folder is treated as read-only).
- **DD-P8 FVS carbon units fixed.** FFE's Stand Carbon Report is in US tons/acre, but `CycleMetrics` (and the aggregation
  with `HA_PER_ACRE`) treat carbon as t C/ha, so FVSjl carbon was reported 2.24× too low. Verified on the Oconee AOI: FVS
  cycle-0 live aboveground carbon = 35,343 (old) vs TreeMap CARBON_L 86,982 short tons; after the fix FVS and Pan agree with
  TreeMap within ~3 %. `FiaSim._parse_sum_and_carbon` now converts (× 2.2417).

## Alternatives (equifinality)
`gafl_{hbad,biomass}_alt{1..4}` come from Pan's anchor-and-explore search (Pan `docs/cmaes_gdigel.md` §4): lineages start at
the calibrated set and move apart while every species group stays within its bias tolerance on the training plots; among the
behavioural sets found, four are picked whose **100-year projections** differ most (biomass, density, species shares).
- **HBAD**: biomass spread up to ±7 % (CV 6.6 % at year 60) and clearly different stand density (157–295 trees/ac at year 60
  vs 204 for v1) — the biomass data constrain growth × survival, not how it splits between density and tree size.
- **Biomass Succession**: projected biomass stays within CV 1.7 % over 100 years; parameter differences (D, S, B_MAX) compensate
  each other — the data pin the biomass trajectory tightly; only late-century tree counts diverge.
- Held-out (val/test, ±10 %): alternatives pass 25–29 (HBAD) and 14–20 (Biomass Succession) of 46–47 groups, versus 30–32 and
  15–18 for a coordinate calibration fitted on the same training plots — plausible, slightly less supported than the best fit.
- **PnET**: anchored exploration could not move PnET (its influential parameters are per-species traits that affect many
  groups at once, so a block step breaks some group's tolerance); the alternatives are instead independent CMA-ES
  (gdigel optimize-mode) fits from random starts with different species traits. Held-out they match a train-only reference
  fit (val 13–15 vs 16, test 9–15 vs 8 groups at ±10 %), yet their 100-year biomass spans 26–39 kg/m² at year 60
  (CV ≈ 15 %) — PnET's long-run projections are the least constrained by the plot data.
Use v1 as the reference and the alternatives to bracket projection uncertainty.

## Tests
- `python3 test/pan_api_smoke.py` — end-to-end API test (meta, samples, parameter upload incl. a bad file, AOI, Pan HBAD /
  Biomass / PnET (baseline and hot & dry) simulations with harvest/thin/plant, FVSjl simulation, engine/parameter mismatch,
  climate preview).
- `test/ui_walkthrough.py` — headless browser walk-through with screenshots (needs Playwright + Chromium; on this machine
  `/workspace/fe_pan/.uivenv` + libraries in `/workspace/fe_pan/.uitools/env/lib`, see the script header).
- `public/mock_pan.js` (`/?mock=1`) — frontend-only mock of the Pan endpoints.
