# Forest Explorer × Pan (branch `forest-explorer-pan`; UI rework on `fe-ui-2`)

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
| `PAN_REPO` | `$WORKSPACE/landis_fe` | Pan checkout (branch `fe-engine` = the refactor with the scenario engine **and HBAD v2**) — the Pan worker runs from here. `PAN_REPO=$WORKSPACE/landis_scen` gives the old scenario-engine branch (HBAD v1 only) |
| `PAN_DEPOT_PATH` | the `JULIA_DEPOT_PATH` exported by `$PAN_REPO/wt_env.sh` when that file exists, else `$PAN_REPO/.julia_depot_wt` (+ `$WORKSPACE/landis_parametrization_julia/.julia_depot` if present) | depot of the Pan worker |
| `PAN_FIA_DB` | `$WORKSPACE/FIA/FIASQLITE2PGSQL/FIADB.duckdb` | curated FIA cohorts for Pan (read only) |
| `FIA_TRAJ_DB` | `$PAN_FIA_DB` | DuckDB for the optional FIA plots / tree-trajectories layer (read only) |
| `PAN_CLIMATE_DB` | `$WORKSPACE/pan_climate/pan_climate.duckdb` | PRISM baseline + forecasts (PnET climate) |
| `PAN_SOIL_CSV`, `PAN_PNET_CLIMATE`, `PAN_CO2_CSV` | `$WORKSPACE/landis_region4/data_r4/…` | PnET plot soils, L4 climate, CO₂ |
| `PAN_THREADS` | `4` | threads of the Pan worker |
| `PAN_AUTOSTART` | `1` | start (and compile) the Pan worker in the background at server start (~1–2 min first time) |
| `PAN_WORKER_LOG` | `./cache/pan_worker.log` | the worker's log |

The Pan worker is a separate Julia process (Pan's own project and depot) talking JSON lines over stdin/stdout; it is
restarted automatically if it dies, and uploaded parameter files are re-sent to it.

## Using it
1. **AOI**: draw or upload polygons as before (samples: `samples/aois/*.geojson`, e.g. Oconee NF GA).
2. **Engine** (plan panel): FVSjl, *Pan · Biomass Succession*, *Pan · PnET-Succession*, *Pan · HBAD v1*, *Pan · HBAD v2 (v2j)* and a
   disabled placeholder *HBAD3 (process-based) — coming*. Each scenario keeps its own engine, parameter set and climate. The
   *ⓘ About this engine* line under the selector shows the description of the chosen engine (it is also in the report):
   - **FVSjl** — Julia port of the USDA Forest Service Forest Vegetation Simulator (FVS, Southern variant; Dixon 2002, Essential FVS) — tree-list growth and yield.
   - **Biomass Succession** — Pan's Julia re-implementation of the LANDIS-II Biomass Succession extension (LANDIS-II: Scheller et al. 2007, Ecological Modelling 201:409–419; Biomass Succession: Scheller & Mladenoff 2004, Ecological Modelling 180:211–229). Species-age cohorts of aboveground biomass; growth limited by ANPP and competition for growing space (B_MAX), age- and growth-related mortality, establishment by shade tolerance and seed availability. Parameters calibrated to FIA in Pan.
   - **PnET-Succession** — Pan's re-implementation of the LANDIS-II PnET-Succession extension (de Bruijn et al. 2014, Ecological Modelling 287:44–57; built on PnET-II, Aber & Federer 1992). Cohort photosynthesis driven by light, temperature, water (soil water balance), VPD and CO₂; carbon allocation to foliage/wood/roots; responds to the climate scenario.
   - **HBAD v1** — Pan cohort model of height, basal area and density (N, QMD, height, age) with structural competition (basal area of larger trees, relative density, Reineke self-thinning) and empirical FIA growth/mortality regressions; establishment as in Biomass Succession. Biomass reported via allometry.
   - **HBAD v2 (v2j)** — biomass-free process cohort model: state N (trees/ha), QMD, mean height, age; height-structured light competition (overtopping) + crowding; species growth curves for diameter and height; growth-stress, background and age mortality; recruitment from seed source and light; calibrated jointly on FIA remeasured trees (Sim A) and stand trajectories (Sim B). Plugin `hbad2` of the refactor Pan (`docs/hbad2_plugin.md`); no climate.
   - **HBAD3 (coming)** — process-grounded successor of HBAD v2 (stem-volume growth with pipe-model leaf area and hydraulic limitation, trait-based height, carbon-starvation mortality, water balance). Not in Pan yet; the entry is disabled. The texts live in `src/PanEngine.jl` (`PAN_ENGINES[].detail`) and `public/pan.js` (FVSjl).
3. **Parameter set** (Pan engines): the default is the engine's **EPA Region 4** set (AL FL GA KY MS NC SC TN, 14 L3
   ecoregions; in-sample 46–49 of 50 species groups within ±10 %), unless you upload a YAML file
   (*Upload YAML…*) or pick a sample. The file's header (`pan_params: {plugin: …}`) decides which engine it belongs to; a
   mismatching file is refused with an offer to switch engine. Format and full roster: Pan `docs/scenario_engine.md` §1.
   HBAD v2 rosters: `r4_hbad2_v2j.yaml` (**default**, current) and the older `r4_hbad2_v2i/v2h/v2g.yaml`, selectable under *Samples…* when the engine is HBAD v2.
   Samples: `r4_{biomass,pnet,hbad}_v1.yaml` (defaults), `gafl_{biomass,pnet,hbad}_v1.yaml` (GA+FL-specific, cross-validated)
   and `gafl_{hbad,biomass,pnet}_alt{1..4}.yaml` — **distinct alternatives** that fit the data about as well as v1 but differ
   in parameters and projections; running a plan with each brackets the projection uncertainty the data cannot resolve
   (see "Alternatives" below).
4. **Vegetation source** (sidebar, *Data source & display*): *TreeMap 2022 pixels* (the AOI is drawn/uploaded; results are a pixel map) or
   *Delineated stands — TreeMap-derived (prototype)*; the third entry (*generative vegetation, TreeMap-independent*) is a disabled placeholder.
   See "Vegetation sources" below.
5. **Plan**: `+ Thin / harvest` (FVS-style residual BA/TPA/SDI, from below/above, DBH range, species — Pan maps it to a
   harvest), `+ Plant`, and for Pan engines `+ Harvest` (LANDIS Biomass Harvest style):
   - cohort selection by species, age range and DBH range;
   - removal: % of each selected cohort's biomass, clearcut, or a residual BA/TPA/SDI target removed from below, from above
     or proportionally;
   - stand selection: all, or ranked (highest biomass, oldest, random) until a % of the AOI area is treated;
   - repeat every N cycles; optional planting right after (species, trees/ac, survival %).
   BA/TPA/SDI/DBH specifications act on Biomass Succession and PnET cohorts through the empirical
   (species, age, biomass) ↔ (TPA, QMD, BA, height) translation; HBAD tracks those variables directly.
6. **Climate…** (PnET only): baseline (PRISM forecast extrapolation, repeat the last T years, or average of the last K
   windows of T years), CO₂ pathway (constant, linear, SSP2-4.5, SSP5-8.5), and any number of changes (variable: mean / min /
   max temperature or precipitation; add °C or % change; step or linear ramp between years; optional months). The modal
   shows AOI-weighted annual charts and a monthly waveform of baseline vs modified; *Reset* restores the default.
7. **Run scenario**, compare scenarios, report — as before. Engine warnings (plots outside the parameter set's ecoregions,
   donors, method notes) appear under the plan and in the status line.
8. **Session bundle** (sidebar, accordion *Session bundle & samples*): *Export* writes a JSON with the AOI patches, every scenario's engine, embedded parameter
   YAML, climate settings and plan; *Import…* restores them (scenarios come back un-run). Sample:
   `samples/bundles/ga_fl_demo_bundle.json` (Oconee AOI, 6 scenarios: FVSjl, HBAD grow/thin+clearcut+replant, Biomass
   hardwood selection, PnET baseline vs hot & dry).

## What Pan reports
Per cycle (start-of-cycle state, before that cycle's management; removals on the row of the cycle they happen in):
BA ft²/ac, TPA, QMD, total cubic volume and board feet, live aboveground carbon (0.5 × AGB), total carbon (adds roots at
0.25 × AG; no dead pools), removed volume and carbon, and the species × 2-inch DBH tally (each cohort's QMD stands for its trees).
- **Volumes** come from FIA Region-4 per-tree equations by species (→ species group → softwood/hardwood) evaluated at each
  cohort's QMD and height: total-stem cubic feet (FIA VOLTSGRS) and board feet **International ¼-inch** (FIA VOLBFGRS, with a
  fitted sawtimber-qualifying probability; FVSjl reports its own board-foot rule). On held-out FIA plots: cubic R² 0.997, bias
  ≈ 0; board feet R² 0.93, bias +1 % (Pan `docs/volume_equations.md`).
- **Vegetation vs site**: a TreeMap pixel resolves to an FIA plot CN whose curated cohorts give the starting vegetation. The
  site comes from **where the pixels are**: the app sends Pan, per plot, the majority EPA L4 ecoregion of its pixels in the AOI
  (`ECO_L4_SHP`) and their centroid, so the parameter eco cell, PnET climate, soils (SOLUS) and site class are those of the AOI —
  TreeMap often imputes plots measured far away (e.g. Oklahoma plots in a Georgia AOI).
- **Plots without curated cohorts** (13 % of TreeMap plots, mostly FIA visits before 2008) use a TreeMap donor plot of the same
  forest type and similar structure (Pan `docs/treemap_donors.md`); the share of area affected is reported.
- **Coverage** = the parameter set's ecoregions at the AOI location: the default Region 4 sets cover AL FL GA KY MS NC SC TN
  (14 L3s). Sample AOIs: Oconee GA 193/193 plots simulated, Uwharrie NC 200/200 (the app caps a run at 200 plots).

## FIA plots — tree trajectories (optional layer)
A sidebar accordion **FIA plots — tree trajectories** with an on/off switch (off by default; nothing is requested until it is
switched on). On, it draws one dot per FIA plot in the map view (zoom ≥ 6; at most 5,000 dots — an evenly spaced sample is
shown and reported if more plots match). Dot size grows modestly with the number of visits; hovering shows the plot key,
tree and visit counts and the main species. **Click a dot** to open a panel with the plot's info (key, mean coordinates,
visits and years, EPA L3/L4, land use, ownership, cut events) and an SVG chart of the per-tree trajectories. Esc or ✕ closes
the panel; changing any option re-renders the open chart.
- **Coordinates are fuzzed.** Public FIA plot coordinates are perturbed (up to about a mile; some private plots are swapped);
  the dot is the mean LAT/LON of the plot's visits in the FIADB `PLOT` table. The tooltip and the info panel say so.
- **Filters** (also filter the dots): species (multi-select by symbol), ecoregion level (EPA L3 / L4) + ecoregion, land use
  (STDORGCD), cut event (ever), terminal fate, ownership group, soft/hardwood, cycle. As in the marimo notebook a tree is
  kept if any of its rows satisfies all selected filters, and then its *full* trajectory is drawn. The dot filter uses per-plot
  value sets (a plot is shown if each filter has at least one matching row on it) — slightly more permissive than the exact
  per-tree test applied when a plot is opened (the info panel then shows "n match filters").
- **Chart options** (same as `FIA_DATA_PREP/tree_trajectories.py`): X / Y axis (biomass, five age sources, diameter, height),
  colour by (species, land use, terminal fate, cut event, ownership, soft/hardwood, cycle, age source, species changed),
  max trajectories (reservoir sample, seed 42), mark imputed points, growth rate Δy/Δx, points All / First / Last per tree
  (first/last use real measurements only), log X, log Y, linear fit + R² (OLS in the displayed — log10 where log — space).
  Rows with a null on either axis are dropped, so phantom back-fill rows disappear from diameter/height axes. Drawing:
  contiguous lines, dashed phantom-backfill prefix, single-measurement trees as dots, open markers for imputed points,
  terminal-fate marker (○ alive, ▽ dead-natural, ✕ dead-disturbance, ▣ harvest — filled when explicit, ◇ disappeared),
  colours tab10 (≤ 10 levels) or turbo. The fit line and its R² are computed server-side; log axes drop non-positive values.
- **API** (`src/FiaTrajectories.jl`): `GET /api/fiaplots/options` (filter values, axes; builds the plot summary on first use,
  ~15–20 s once per server start, then cached in memory), `POST|GET /api/fiaplots/points` (`west,south,east,north` + filters
  → plots, `n_match`, `truncated`), `POST /api/fiaplots/trajectories` (`{plot, filters, options}` → per-tree series, levels
  with colours, fit). Filter values are validated against the known option lists before they reach SQL.
- **Data**: table `tree_trajectories` + `PLOT` of the curated FIA DuckDB, opened **read-only** (never written). Path:
  `FIA_TRAJ_DB`, else `PAN_FIA_DB`, else `$WORKSPACE/FIA/FIASQLITE2PGSQL/FIADB.duckdb` (`serve_pan.sh` exports it and warns
  if missing). Adds the Julia package `DuckDB` to this app.

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
- `test/ui_walkthrough.py` — headless browser checks (initial state/accordions, engine list + descriptions, one scenario per engine incl. HBAD v2 grow/thin/plant and an FVSjl PPE thin, PPE inside the plan, compare + report with PNG decode, bundle, console errors; PASS/FAIL, exit code). Needs Playwright + Chromium: `LD_LIBRARY_PATH=/workspace/fe_pan/.uitools/env/lib PLAYWRIGHT_BROWSERS_PATH=/workspace/fe_pan/.uitools/pw FONTCONFIG_FILE=/workspace/fe_pan/.uitools/env/etc/fonts/fonts.conf /workspace/fe_pan/.uivenv/bin/python test/ui_walkthrough.py http://127.0.0.1:8089 <outdir>`.
- `test/ui_stands.py` — stands layer + the stands vegetation source (click / polygon selection, HBAD v2 and FVSjl per stand, stand colors per cycle).
- Note: Oxygen serves `public/` from memory — restart the server after editing frontend files.
- `public/mock_pan.js` (`/?mock=1`) — frontend-only mock of the Pan endpoints.

## Stands layer (prototype, 2026-10-19)
`public/stands.js` + `stands.css`: optional sidebar group "Stands — delineated (prototype)". It shows stand polygons from the standmap
pipeline (/workspace/standmap; TreeMap 2022 + SOLUS + AlphaEarth 2022/2025 segmentation, split by roads / water / parcel lines; recent
clear-cuts kept as their own stands).
- Popup: status, area, harvest year, TreeMap structure (BA, height, QMD, TPA, canopy, softwood share), TreeMap coverage, soils.
- **Simulate this stand** sets the stand polygon as the AOI, so the existing AOI pipeline (TreeMap pixels → FIA plots → FVSjl / Pan) runs
  on that stand.
- Data: `public/stands/<area>.geojson` (git-ignored; regenerate with standmap `scripts/24_export_fe.py`). Prototype area: Baker County,
  FL (one of the 5 FL parcel counties); 15,341 stands, 128,121 ha.
- UI check: `test/ui_stands.py` (Playwright; needs `FONTCONFIG_FILE=/workspace/fe_pan/.uitools/env/etc/fonts/fonts.conf` in addition to
  `LD_LIBRARY_PATH`).

## UI layout (fe-ui-2)
Nothing is preallocated for empty content; groups appear when they have something to show.
- **Sidebar = accordions** (`<details>`, keyboard: Tab to a summary, Enter/Space toggles): *Data source & display* (open: data source, vegetation
  source; the display attribute / basemap / opacity live in a nested, closed *Display* section), *Area of interest* (open), *Layers: FIA plots,
  stands* (closed; the two layer groups inside keep their ids `fiaGrp`, `stGrp`), *Session bundle & samples* (closed). Hidden until there is
  content: *Simulated layer*, *AOI summary*, *Species × DBH class* (shown when an AOI/simulation exists).
- **Bottom panel**: one accordion *Scenarios & plan*, closed until an AOI exists (then opened). The *Projection* column (chart, metric, time
  slider) is hidden until a scenario has a result, so the plan uses the full width before that.
- **PPE harvest budget is part of the scenario** (FVSjl): a collapsible *Landscape harvest budget (PPE)* section of the plan panel (ids
  unchanged: `lsTarget`, `lsPriority`, `lsCredit`, `lsYear`, `lsCycles`, `lsPeriod`, `runLandscape`). *Preview* runs `/api/landscape` as before;
  *Add as thin action* writes the same budget into the plan as a Thin action with `mode: ppe` (the existing per-stand scheduling), so it runs
  with the scenario, appears in the plan summary, in the report's "Harvest scheduling" and in saved plans / bundles. As before it is FVSjl-only.
- Ids used by the tests are unchanged (`#engine`, `#runSim`, `#status`, `#cards`, `#file`, `#scenarioChips`, ...); new: `#vegSource`,
  `#standTools`, `#stSelPoly`, `#stSelClear`, `#planDet`, `#planR`, `#engineHelp`/`#engineText`, `#accSource/#accAoi/#accLayers`.

## Vegetation sources and data flow
**TreeMap 2022 pixels** (default, unchanged): AOI polygons → `/api/aoi` (union MultiPolygon → TreeMap pixels → FIA plots) → `/api/simulate`
(`geometry`) simulates the plots once per AOI; the result is rendered as a pixel raster per (metric, cycle) (`/simimage/...`).

**Delineated stands — TreeMap-derived (prototype)**: the simulation unit is the **stand** (`public/stands/<area>.geojson`, Baker County FL).
1. The AOI is a *selection*: click a stand on the map to toggle it (cyan outline), or *Select stands in polygon* (draw a polygon; every stand
   whose centroid lies inside is added). At most 100 stands. The selection is the AOI (patch list shows `Stand <id>`); `/api/aoi` summarises
   it as before. Changing the selection clears scenarios (stale), as any AOI change does.
2. *Run* sends `POST /api/simulate` with `stands: [{id, geometry}]` next to the usual engine/plan/params/climate. For each stand the backend
   resolves its polygon to TreeMap pixels and plots with the existing `_resolve_plots` (plot acres = pixel share × pixel area), so the stand's
   vegetation is the **pixel-share-weighted composite of the TreeMap plots of its pixels**; the Pan site (majority L4, centroid) is that of the
   stand's pixels. Pan engines are called once per stand; FVSjl simulates the union of the plots once. Each stand is capped to its 200 largest plots
   (FVSjl: 200 distinct plots overall; both are reported as warnings). A PPE thin (FVSjl) budgets across the plots of all selected stands.
3. Results are aggregated **per stand** (acre-weighted means of live carbon, BA, volume per cycle → `stands[i].map`, plus the stand's cycle
   series) and over all stands (`cycles`, `dist`, `domains`; `image` is null). Stands without TreeMap pixels (e.g. harvested) are skipped and
   reported.
4. The projection map colors the **stand polygons** (GeoJSON source `simstands`) per cycle and metric with the same colour ramps and
   auto-scaled domain as the raster; the time slider and metric selector work as before; click a stand for its projected values. The
   report snapshot shows the colored stands.
Limitations (prototype): the session bundle stores patches only (a stands session re-imports as TreeMap-source patches, un-run); the
generative (TreeMap-independent) source is not implemented.

## Report
*Compare → ↓ Report* writes a self-contained HTML: AOI map snapshot, metric charts, species × DBH, endpoints, harvest scheduling (PPE), the
scenarios (plan, engine, parameter set, vegetation source) and an **Engines** section with the description of every engine used.
- **Map image fix**: the snapshot was read from the WebGL canvas synchronously at click time, an arbitrary moment at which a drawing buffer can
  be cleared (blank/transparent PNG), and as a full-DPR PNG it could be many MB. It is now read inside a `render` event of a forced repaint,
  downscaled to ≤ 1100 px, checked for content (a blank or cross-origin-tainted canvas is detected) and embedded as `<img id="aoiimg">`; if the
  canvas cannot be read the report draws the AOI outline as an SVG instead of a broken image (with a note). `test/ui_walkthrough.py`
  downloads the report in the headless browser, decodes the PNG and checks signature, `naturalWidth` and pixel variance.

## Pan engine versions (HBAD v2) and parity with the old worker
`PAN_REPO` now defaults to the refactor checkout. Biomass Succession and FVSjl results are identical to the old worker on the test AOI; PnET and
HBAD v1 differ slightly because the refactor fixed order-dependence (Pan commits B32/B37: results used to depend on DuckDB's unordered result
sets). The old worker is *not* reproducible run to run (HBAD v1 grow-only +40 yr BA 151.1–155.9 ft²/ac over 3 runs of the same request); the new one
is deterministic. Comparison on a 62-plot AOI (Oconee, 4 cycles x 10 yr, default sets, same requests to the old worker on 8088 and the new one on 8089): Biomass Succession thin+plant and FVSjl thin+plant identical in every cycle; PnET grow-only within 0.2 % (volume 1,203,483 vs 1,204,241 ft3 at +40 yr); HBAD v1 grow-only +40 yr BA 151.0 (new) vs 155.9 / 151.1 (two old runs); HBAD v1 residual-BA harvest from below +40 yr BA 125.5 (new, reproducible) vs 115.7 / 114.5 / 116.6 (three old runs) - the order in which equal cohorts are removed 'from below' changed with the B32/B37 ordering fixes; the same harvest with proportional removal agrees (139.3 new vs 138.0 / 133.8 old).
