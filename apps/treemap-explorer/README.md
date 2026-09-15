# TreeMap Growth Explorer

Localhost web app that displays the USFS **TreeMap 2022 CONUS** raster, lets the user
draw/upload an AOI polygon, resolves each forested pixel to its FIA plot, simulates
forest growth with **FVSjl** (as a black box), and shows aggregated species / DBH-class
distributions of biomass, volume, and carbon over a projection — with a time-slider map.

## Run

```sh
./serve.sh                 # http://127.0.0.1:8080
HOST=0.0.0.0 PORT=8087 ./serve.sh
```

`serve.sh` derives every path from its own location — nothing absolute is baked in:

- **Julia depot** — honors an existing `$JULIA_DEPOT_PATH`, else defaults to the
  **persistent** depot beside the workspace (`<workspace>/.julia_depot`). `~/.julia`
  under `$HOME` is wiped on restart; only the workspace volume persists, so the app's
  packages must live in the persistent depot. First run there needs
  `julia --project=. -e 'using Pkg; Pkg.instantiate()'` (with that depot exported).
- **Data dir** — honors `$TREEMAP_DATA`, else resolved by `default_datadir()`
  (env → module-relative default). Override: `TREEMAP_DATA=/path ./serve.sh`.

## Layout
- `src/`      — Julia backend (raster IO, AOI clip, plot resolver, FVSjl adapter, HTTP API)
- `public/`   — static frontend (MapLibre GL + charts)
- `data/`     — placeholder; real data lives OUTSIDE the repo (see below)
- `test/`     — `smoke.jl` end-to-end check

## Data (outside this repo)
TreeMap 2022 CONUS raster + tree table live under the data dir (default a `treemap/`
sibling of the repo; 4.5 GB source zip RDS-2025-0032 from the USFS Research Data
Archive), deliberately **outside** the FVSjl git tree. The tm_id → PLT_CN crosswalk is
the raster attribute table (`.tif.vat.dbf`) shipped with the GeoTIFF — which also carries
precomputed cycle-0 stand state (BA/QMD/TPA/volume/biomass/carbon); per-tree attributes
are in `TreeMap2022_CONUS_Tree_Table.csv` (the FVS tree list).

## Simulate boundary
Single adapter: `simulate(plotInputs, scenario, cycles) -> per-cycle metrics`.
FVSjl sits behind it; the same contract is what a future wasm FVSjl fulfils in-browser.
