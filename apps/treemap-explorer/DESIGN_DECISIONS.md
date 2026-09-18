# Forest Explorer — multi-patch AOI, multi-scenario, compare report

Design decisions taken autonomously (2026-09-17). Please validate; each has a rationale + how to change it.

## Feature 1 — Multiple polygons per AOI ("patches")
- **DD1. AOI = an ordered list of polygon "patches."** Each patch is drawn or uploaded
  separately and is individually deletable. The effective AOI sent to the backend is the
  **union as a GeoJSON `MultiPolygon`**. Rationale: matches "several non-contiguous forested
  patches"; lets the owner build the holding incrementally. (Alt: a single freehand
  MultiPolygon — rejected: can't delete one patch.)
- **DD2. Plots resolve over the union** — dedup by PLT_CN, acres summed. No double-count
  if two patches hit the same plot (dedup). Backend `aoi_tally` rasterizes the whole
  MultiPolygon in one pass (GDAL handles it), so no per-patch loop.
- **DD3. Adding/removing a patch re-resolves AOI stats and CLEARS existing scenarios,
  announcing it in the status line** ("AOI changed … N prior scenarios cleared (stale)").
  Rationale: scenarios were simulated on the old AOI; keeping them would be misleading.
  Chose a status-line notice over a blocking modal because a patch is already added/removed
  by the time we'd prompt, so a cancel would need a messy rollback mid-draw. (Alt: keep +
  mark stale and re-simulate on the new AOI — more complex; deferred.)
- **DD4. Upload can contribute multiple patches** — every polygon in an uploaded
  GeoJSON/shapefile becomes a patch (MultiPolygon uploads split into parts too).

## Feature 2 — Multiple scenarios per AOI
- **DD5. A "scenario" = { name, plan, result }.** The plan builder edits a *draft* plan;
  "Run scenario" simulates it and stores the result as a new scenario. Multiple scenarios
  coexist for the same AOI. Rationale: "simulate different scenarios" on one holding.
- **DD6. One ACTIVE scenario drives the map raster, time slider, left species×DBH stats,
  and projection chart** (i.e. the existing single-sim views show the active scenario). A
  scenario selector (chips) switches the active one. **Re-running while a scenario is active
  OVERWRITES it in place** (tweak-and-rerun); **"+ New" forks a fresh draft cloned from the
  current plan** (keeps its actions so you can change one knob and Run to get a new scenario).
- **DD7. Auto-named `Scenario N`, editable inline**, with a plan-summary subtitle
  (e.g. "grow-only", "thin BA→50 @+10yr, plant DF"). Rationale: quick to run, still labelable.
- **DD8. Backend caches sim images PER TOKEN** (was: last-sim-only). Each scenario run gets a
  fresh token and its own cached raster set; the frontend keeps every scenario's
  `{token,corners,domains,…}`. LRU cap = **12 tokens** to bound memory (older scenarios'
  rasters evict first; their numeric results are still kept client-side). Rationale: lets
  every scenario's raster stay viewable without recompute.
- **DD9. Scenario cap = 8 per AOI** (soft; a Run that would add a 9th is refused with a
  status message telling the user to delete one or overwrite the active one — re-running the
  active scenario is always allowed). 8 = the compare palette length. Rationale: keeps compare
  legible + memory bounded (matches the 12-token image LRU with headroom).

## Feature 3 — Compare report (on demand)
- **DD10. A "Compare" button, enabled only when ≥2 scenarios are simulated.** Opens a
  compare panel; nothing is computed until asked. Rationale: user's "separate compare button."
- **DD11. Compare panel = (a) a multi-line projection chart** (one line per scenario) with a
  metric switch (total carbon / live carbon / standing volume / basal area / removed volume),
  **(b) an endpoint table** per scenario (final BA·carbon·volume, cumulative removed volume,
  #plots) with **deltas vs the first scenario as baseline**, and **(c) the final-cycle
  species×DBH comparison is deferred** (kept out of v1 to stay focused).
- **DD12. "Download report" produces a self-contained HTML file** (inline SVG chart + tables +
  AOI summary + each scenario's plan), openable/printable offline. Rationale: "produce a
  report of the comparison." (Alt: PDF — needs a lib; HTML prints to PDF from the browser.)

## Scope / non-goals for this pass
- The PPE cross-stand budget panel stays AOI-wide (not per-scenario) for now.
- Full post-harvest re-projection unchanged.
- Report is comparison-of-simulated-scenarios; not a per-stand ledger.

## Backend changes
- `SIM_IMAGES`: `Dict{token=>Dict{(metric,cycle)=>png}}` + LRU token cap; `/simimage/{token}/…`
  now honors the token.
- `geom_from_geojson` / `_resolve_plots` / `_render_sim_images`: accept `MultiPolygon`
  (verify GDAL path). No new endpoints — compare is client-side over stored results.
