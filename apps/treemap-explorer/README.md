# TreeMap Growth Explorer

Localhost web app that displays the USFS **TreeMap 2022 CONUS** raster, lets the user
draw/upload an AOI polygon, resolves each forested pixel to its FIA plot, simulates
forest growth with **FVSjl** (as a black box), and shows aggregated species / DBH-class
distributions of biomass, volume, and carbon over a projection — with a time-slider map.

## Layout
- `src/`      — Julia backend (raster IO, AOI clip, plot resolver, FVSjl adapter, HTTP API)
- `public/`   — static frontend (MapLibre GL + charts)
- `data/`     — small derived/cache artifacts (NOT the raster; raster lives outside the repo)
- `test/`

## Data (outside this repo)
TreeMap 2022 CONUS raster + tree table live at `/workspace/treemap/` (4.5 GB source
zip RDS-2025-0032 from the USFS Research Data Archive), deliberately **outside** the
FVSjl git tree. The tm_id → PLT_CN crosswalk is the raster attribute table shipped with
the GeoTIFF; per-tree attributes are in `TreeMap2022_CONUS_Tree_Table.csv`.

## Simulate boundary
Single adapter: `simulate(plotInputs, scenario, cycles) -> per-cycle metrics`.
FVSjl sits behind it; the same contract is what a future wasm FVSjl fulfils in-browser.
