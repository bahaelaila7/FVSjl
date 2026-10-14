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
`serve_pan.sh` sets everything explicitly and writes only inside this app folder:
| variable | default | meaning |
|---|---|---|
| `JULIA_DEPOT_PATH` | `./.julia_depot_fe:/workspace/FVSclaude/.julia_depot` | private depot layered over the workspace FVSjl depot |
| `TREEMAP_DATA` | `/workspace/FVSclaude/treemap` | TreeMap raster/VAT/tree table (read only) |
| `FIA_DB` | `/workspace/FVSclaude/SQLite_FIADB_ENTIRE.db` | FVSjl master FIA DB (read only) |
| `FE_CACHE_DIR` | `./cache/sim_cache` | FVSjl sub-DB cache (kept out of the TreeMap data dir) |
| `PAN_REPO` | `/workspace/landis_scen` | Pan checkout (branch `scenario-engine`) — the Pan worker runs from here |
| `PAN_THREADS` | `4` | threads of the Pan worker |
| `PAN_AUTOSTART` | `1` | start (and compile) the Pan worker in the background at server start (~1–2 min first time) |
| `PAN_WORKER_LOG` | `./cache/pan_worker.log` | the worker's log |

The Pan worker is a separate Julia process (Pan's own project and depot) talking JSON lines over stdin/stdout; it is
restarted automatically if it dies, and uploaded parameter files are re-sent to it.

## Using it
1. **AOI**: draw or upload polygons as before (samples: `samples/aois/*.geojson`, e.g. Oconee NF GA).
2. **Engine** (plan panel): FVSjl, *Pan · Biomass Succession*, *Pan · PnET-Succession*, *Pan · HBAD*. Each scenario keeps its
   own engine, parameter set and climate.
3. **Parameter set** (Pan engines): the default GA+FL set of the engine is used unless you upload a YAML file
   (*Upload YAML…*) or pick a sample. The file's header (`pan_params: {plugin: …}`) decides which engine it belongs to; a
   mismatching file is refused with an offer to switch engine. Format and full roster: Pan `docs/scenario_engine.md` §1.
   Samples: `gafl_{biomass,pnet,hbad}_v1.yaml` (calibrated) and `gafl_hbad_alt{1..4}.yaml` — **distinct alternatives** that fit
   the training data about as well as `gafl_hbad_v1` but differ in growth/mortality/SDImax; running a plan with each brackets
   the projection uncertainty that the data cannot resolve (see "Alternatives" below).
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
reported as `out_of_region` (use FVSjl there).

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
`gafl_hbad_alt{1..4}` come from Pan's anchor-and-explore search (Pan `docs/cmaes_gdigel.md` §4): lineages start at the
calibrated set and move as far apart as possible while every species group stays within its bias tolerance; four
behavioural sets far apart in parameter space are then picked. Over 100 years on 60 GA/FL plots they keep mean biomass
within ±8 % of each other (CV ≤ 6.5 %) but differ in stand structure — e.g. at year 60, 196 vs 357 trees/ac — because
the biomass data constrain growth × survival, not how it splits between density and tree size.

## Tests
- `python3 test/pan_api_smoke.py` — end-to-end API test (meta, samples, parameter upload incl. a bad file, AOI, Pan HBAD /
  Biomass / PnET (baseline and hot & dry) simulations with harvest/thin/plant, FVSjl simulation, engine/parameter mismatch,
  climate preview).
- `test/ui_walkthrough.py` — headless browser walk-through with screenshots (needs Playwright + Chromium; on this machine
  `/workspace/fe_pan/.uivenv` + libraries in `/workspace/fe_pan/.uitools/env/lib`, see the script header).
- `public/mock_pan.js` (`/?mock=1`) — frontend-only mock of the Pan endpoints.
