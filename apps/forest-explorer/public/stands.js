'use strict';
/* Optional layer "Stands (delineated)": stand polygons from the standmap pipeline (standmap project: TreeMap + SOLUS + AlphaEarth
   2022/2025 segmentation, cut by roads / water / parcels, recent clear-cuts kept as their own stands). Static GeoJSON per prototype area
   (public/stands/<area>.geojson, EPSG:4326). Click a stand → its attributes; "Simulate this stand" makes the stand polygon the AOI, so the
   existing AOI pipeline (TreeMap pixels → FIA plots → FVSjl / Pan) runs on exactly that stand.
   Vegetation source "Delineated stands" (#vegSource = stands): clicking a stand toggles it into the SELECTION, "Select stands in polygon" adds every
   stand whose centroid lies inside a drawn polygon; the selection is the AOI (app.js setStandAoi) and the simulation unit is the stand
   (POST /api/simulate with `stands`). After a run, showSim() colors the selected stand polygons per cycle (instead of the pixel raster). */
const StandsUI = (() => {
  const SRC = 'stands', L_FILL = 'stands-fill', L_LINE = 'stands-line', L_SEL = 'stands-sel';
  const AREAS = [{ id: 'baker_fl', label: 'Baker County, FL (prototype)', url: '/stands/baker_fl.geojson' }];
  const STATUS = { 1: ['forest', '#2e7d32'], 2: ['harvested (loss since 2017)', '#e65100'], 3: ['regenerating', '#f9a825'], 4: ['forest (not in TreeMap)', '#6a1b9a'] };
  const S = { on: false, data: null, area: AREAS[0], popup: null, color: 'status', sel: new Set(), byId: new Map(), cent: null, timer: null, simOn: false, simPopup: null };
  const MAX_SEL = 100;                                           // server limit (MAX_STANDS)
  const L_PICK = 'stands-pick', L_PICKL = 'stands-pickline', SIM_SRC = 'simstands', L_SIMF = 'simstands-fill', L_SIML = 'simstands-line';
  const RAMPS = { greens: ['#f7fcf5', '#c7e9c0', '#74c476', '#238b45', '#00441b'], viridis: ['#440154', '#414487', '#2a788e', '#22a884', '#7ad151', '#fde725'],
                  magma: ['#000004', '#3b0f70', '#8c2981', '#de4968', '#fe9f6d', '#fcfdbf'] };
  const isStandsMode = () => typeof vegSource !== 'undefined' && vegSource === 'stands';
  const $ = id => document.getElementById(id);
  const esc = t => String(t).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const f1 = v => (v == null || Number.isNaN(+v)) ? '–' : (+v).toLocaleString(undefined, { maximumFractionDigits: 1 });
  const say = t => { const e = $('stStatus'); if (e) e.textContent = t; };

  function colorExpr() {
    if (S.color === 'status') return ['match', ['get', 'status'], 1, STATUS[1][1], 2, STATUS[2][1], 3, STATUS[3][1], 4, STATUS[4][1], '#888'];
    const ramps = { BALIVE: [0, 150], STANDHT: [0, 100], QMD: [0, 16], CANOPYPCT: [0, 100], GEN_BA: [0, 150], GEN_TPA: [0, 600], GEN_QMD: [0, 16], GEN_BA_SD: [0, 100] };
    const [lo, hi] = ramps[S.color];
    return ['case', ['==', ['typeof', ['get', S.color]], 'number'],
      ['interpolate', ['linear'], ['get', S.color], lo, '#f7fcb9', (lo + hi) / 2, '#41ab5d', hi, '#004529'], '#bbbbbb'];
  }
  function legend() {
    const el = $('stLegend'); if (!el) return;
    if (S.color === 'status') { el.innerHTML = Object.values(STATUS).map(([l, c]) => `<span class="stkey"><i style="background:${c}"></i>${esc(l)}</span>`).join(''); return; }
    const u = { BALIVE: 'ft²/ac', STANDHT: 'ft', QMD: 'in', CANOPYPCT: '%', GEN_BA: 'ft²/ac', GEN_TPA: 'trees/ac', GEN_QMD: 'in', GEN_BA_SD: 'ft²/ac' }[S.color];
    el.innerHTML = `<span class="stramp"></span> ${esc(S.color)} (${u}); grey = no TreeMap pixels (e.g. harvested)`;
  }
  async function load() {
    say('loading stands…');
    const r = await fetch(S.area.url); if (!r.ok) throw new Error('stands file not found: ' + S.area.url);
    S.data = await r.json(); S.byId = new Map(S.data.features.map(f => [f.properties.stand_id, f])); S.cent = null; return S.data;
  }
  function addLayers() {
    const before = map.getLayer('aoi-fill') ? 'aoi-fill' : undefined;          // below the AOI so a selected stand-AOI stays visible
    if (!map.getSource(SRC)) map.addSource(SRC, { type: 'geojson', data: S.data, promoteId: 'stand_id' });
    if (!map.getLayer(L_FILL)) map.addLayer({ id: L_FILL, type: 'fill', source: SRC, paint: { 'fill-color': colorExpr(), 'fill-opacity': 0.45 } }, before);
    if (!map.getLayer(L_LINE)) map.addLayer({ id: L_LINE, type: 'line', source: SRC, paint: { 'line-color': '#1b1b1b', 'line-width': ['interpolate', ['linear'], ['zoom'], 9, 0.2, 14, 1.2] } }, before);
    if (!map.getLayer(L_PICK)) map.addLayer({ id: L_PICK, type: 'fill', source: SRC, filter: pickFilter(), paint: { 'fill-color': '#00e5ff', 'fill-opacity': 0.25 } }, before);
    if (!map.getLayer(L_PICKL)) map.addLayer({ id: L_PICKL, type: 'line', source: SRC, filter: pickFilter(), paint: { 'line-color': '#00e5ff', 'line-width': 2.5 } }, before);
    if (!map.getLayer(L_SEL)) map.addLayer({ id: L_SEL, type: 'line', source: SRC, filter: ['==', ['get', 'stand_id'], -1], paint: { 'line-color': '#00e5ff', 'line-width': 3 } }, before);
  }
  function removeLayers() { for (const l of [L_SEL, L_PICKL, L_PICK, L_LINE, L_FILL]) if (map.getLayer(l)) map.removeLayer(l); if (map.getSource(SRC)) map.removeSource(SRC); if (S.popup) S.popup.remove(); }
  function bbox(fc) {
    let w = 180, s = 90, e = -180, n = -90;
    const walk = c => { if (typeof c[0] === 'number') { w = Math.min(w, c[0]); e = Math.max(e, c[0]); s = Math.min(s, c[1]); n = Math.max(n, c[1]); } else c.forEach(walk); };
    fc.features.forEach(f => walk(f.geometry.coordinates)); return [[w, s], [e, n]];
  }
  function popupHtml(p) {
    const st = STATUS[p.status] ? STATUS[p.status][0] : p.status;
    const rows = [['status', st], ['area', f1(p.area_ha) + ' ha'], p.loss_year ? ['harvest (loss) year', p.loss_year] : null,
      ['basal area', f1(p.BALIVE) + ' ft²/ac'], ['stand height', f1(p.STANDHT) + ' ft'], ['QMD', f1(p.QMD) + ' in'], ['live TPA', f1(p.TPA_LIVE)],
      ['canopy cover', f1(p.CANOPYPCT) + ' %'], ['softwood share (BA)', f1(100 * (p.SOFTWOOD_FRAC ?? NaN)) + ' %'],
      ['TreeMap pixels', f1(100 * (p.tm_frac ?? 0)) + ' % of the stand'],
      p.GEN_BA != null ? ['generative (prototype)', `BA ${f1(p.GEN_BA)} ± ${f1(p.GEN_BA_SD)} ft²/ac · ${f1(p.GEN_TPA)} TPA · QMD ${f1(p.GEN_QMD)} in`] : null, ['TreeMap plots', `${p.tm_n ?? 0} (main ${f1(100 * (p.tm_share ?? 0))} %)`],
      ['available water 0–30 cm', p.awc030 == null ? '–' : `${f1(300 * p.awc030)} mm (${(+p.awc030).toFixed(3)} cm/cm)`], ['clay / sand 0–30 cm', `${f1(p.clay030)} / ${f1(p.sand030)} %`]].filter(Boolean);
    return `<div class="stpop"><b>Stand ${p.stand_id}</b><table>${rows.map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(v)}</td></tr>`).join('')}</table>
      <button type="button" class="btn" id="stSim" ${p.status === 2 && !(p.tm_frac > 0) ? 'title="harvested stand: no TreeMap pixels; the AOI pipeline has no vegetation for it yet (donors pending)"' : ''}>Simulate this stand</button></div>`;
  }
  const pickFilter = () => ['in', ['get', 'stand_id'], ['literal', [...S.sel]]];
  function refreshPick() {
    for (const l of [L_PICK, L_PICKL]) if (map.getLayer(l)) map.setFilter(l, pickFilter());
    const c = $('stSelCount'); if (c) c.textContent = `${S.sel.size} stand${S.sel.size === 1 ? '' : 's'} selected`;
  }
  // the selection is the AOI: push it to app.js (debounced — every change re-resolves the AOI and clears stale scenarios)
  function commit() {
    clearTimeout(S.timer);
    S.timer = setTimeout(() => {
      const ids = [...S.sel].filter(id => S.byId.has(id));
      setStandAoi(ids, ids.map(id => S.byId.get(id).geometry));
    }, 500);
  }
  function toggleStand(id) {
    if (S.sel.has(id)) S.sel.delete(id);
    else if (S.sel.size >= MAX_SEL) { setStatus(`at most ${MAX_SEL} stands per simulation — remove one first`); return; }
    else S.sel.add(id);
    refreshPick(); commit();
  }
  function centroid(f) {                                          // area-weighted centroid of the largest polygon's outer ring
    const polys = f.geometry.type === 'Polygon' ? [f.geometry.coordinates] : f.geometry.coordinates;
    let best = null, ba = -1;
    for (const pl of polys) { const r = pl[0]; let a = 0, cx = 0, cy = 0;
      for (let i = 0, n = r.length - 1; i < n; i++) { const k = r[i][0] * r[i + 1][1] - r[i + 1][0] * r[i][1]; a += k; cx += (r[i][0] + r[i + 1][0]) * k; cy += (r[i][1] + r[i + 1][1]) * k; }
      if (Math.abs(a) > ba) { ba = Math.abs(a); best = a ? [cx / (3 * a), cy / (3 * a)] : r[0]; } }
    return best;
  }
  function inRing(pt, ring) { let c = false; for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const [xi, yi] = ring[i], [xj, yj] = ring[j]; if ((yi > pt[1]) !== (yj > pt[1]) && pt[0] < (xj - xi) * (pt[1] - yi) / (yj - yi) + xi) c = !c; } return c; }
  async function selectInPolygon(geom) {
    if (!S.data) { try { await load(); } catch (err) { say('error: ' + err.message); return; } }
    if (!S.cent) S.cent = new Map(S.data.features.map(f => [f.properties.stand_id, centroid(f)]));
    const ring = geom.coordinates[0]; let added = 0, skipped = 0;
    for (const [id, c] of S.cent) {
      if (S.sel.has(id) || !c || !inRing(c, ring)) continue;
      if (S.sel.size >= MAX_SEL) { skipped++; continue; }
      S.sel.add(id); added++;
    }
    refreshPick(); commit();
    setStatus(`${added} stand${added === 1 ? '' : 's'} added inside the polygon (centroid test)` + (skipped ? `; ${skipped} not added (limit ${MAX_SEL})` : ''));
  }
  function clearSelection() { S.sel.clear(); refreshPick(); clearTimeout(S.timer); setStandAoi([], []); }
  function setSelection(ids) { S.sel = new Set(ids); refreshPick(); }                // from app.js (patch list ×): no commit
  function onClick(e) {
    const f = e.features && e.features[0]; if (!f) return;
    if (isStandsMode() && !S.simOn) { toggleStand(f.properties.stand_id); return; }
    if (S.simOn) return;                                          // projected stands: the sim popup answers the click
    map.setFilter(L_SEL, ['==', ['get', 'stand_id'], f.properties.stand_id]);
    if (S.popup) S.popup.remove();
    S.popup = new maplibregl.Popup({ maxWidth: '320px' }).setLngLat(e.lngLat).setHTML(popupHtml(f.properties)).addTo(map);
    const b = document.getElementById('stSim');
    if (b) b.onclick = () => {
      const full = S.data.features.find(x => x.properties.stand_id === f.properties.stand_id);   // full-resolution geometry (not the tile-clipped one)
      aoiParts = [full.geometry]; submitAoiUnion();                                               // the stand becomes the AOI (app.js)
      if (S.popup) S.popup.remove();
      setStatus(`stand ${f.properties.stand_id} (${f1(f.properties.area_ha)} ha) set as the AOI — choose an engine and Simulate`);
    };
  }
  async function toggle(on) {
    S.on = on;
    if (!on) { removeLayers(); say(''); return; }
    try {
      if (!S.data) await load();
      addLayers(); if (!(isStandsMode() && S.sel.size)) map.fitBounds(bbox(S.data), { padding: 30, duration: 600 });
      const n = S.data.features.length, ha = S.data.features.reduce((a, f) => a + (f.properties.area_ha || 0), 0);
      say(`${n.toLocaleString()} stands · ${Math.round(ha).toLocaleString()} ha forest land`);
      legend();
    } catch (err) { say('error: ' + err.message); $('stOn').checked = false; S.on = false; }
  }
  /* ---------------- vegetation source ---------------- */
  async function setVegSource(v) {
    if (v === vegSource) return;
    const had = aoiParts.length;
    S.sel.clear(); clearSim();
    vegSource = v; clearAoi(); refreshPick();
    const st = v === 'stands';
    $('aoiTools').hidden = st; $('standTools').hidden = !st;
    $('aoihint').innerHTML = st
      ? 'Vegetation source <b>delineated stands</b>: the AOI is the set of selected stands. Each stand is simulated on its own (vegetation = the pixel-share-weighted composite of the TreeMap plots of its pixels) and results are shown per stand.'
      : 'Build the holding from one or more <b>patches</b> — draw each (click vertices, double-click to finish) or upload a <b>GeoJSON</b> / <b>zipped shapefile</b>. Patches need not be contiguous; the AOI resolves to the union of their FIA plots.';
    $('vegHint').textContent = st
      ? 'Stands (prototype, TreeMap-derived): click stands on the map to select them, or select all stands inside a drawn polygon. Results color the stand polygons per cycle.'
      : 'TreeMap pixels: the AOI is drawn or uploaded; results are a pixel map.';
    if (st) {
      $('accLayers').open = true; $('stGrp').open = true; $('accAoi').open = true;
      if (!$('stOn').checked) { $('stOn').checked = true; await toggle(true); }
      setStatus('stands source: click stands to select them (or use “Select stands in polygon”), then pick an engine and Run');
    } else setStatus('TreeMap pixels source' + (had ? ' — selection cleared' : ''));
  }

  /* ---------------- projected stands (simulation result) ---------------- */
  function simFeatures(sim, metric, cycle) {
    const by = new Map((sim.stands || []).map(s => [s.id, s]));
    return { type: 'FeatureCollection', features: aoiStandIds.map((id, i) => {
      const s = by.get(id), arr = s && s.map && s.map[metric], v = arr ? arr[cycle] : null;
      return { type: 'Feature', properties: { stand_id: id, v: v == null ? null : v }, geometry: aoiParts[i] }; }) };
  }
  function showSim(sim, metric, cycle, opacity, cmap, dom) {
    const stops = RAMPS[cmap] || RAMPS.viridis, lo = dom[0], hi = dom[1] > dom[0] ? dom[1] : dom[0] + 1;
    const ramp = ['interpolate', ['linear'], ['coalesce', ['get', 'v'], lo]];
    stops.forEach((c, i) => ramp.push(lo + (hi - lo) * i / (stops.length - 1), c));
    const color = ['case', ['==', ['get', 'v'], null], '#7a7a7a', ramp];
    const fc = simFeatures(sim, metric, cycle);
    if (map.getSource(SIM_SRC)) map.getSource(SIM_SRC).setData(fc);
    else map.addSource(SIM_SRC, { type: 'geojson', data: fc });
    const before = map.getLayer('aoi-line') ? 'aoi-line' : undefined;
    if (!map.getLayer(L_SIMF)) map.addLayer({ id: L_SIMF, type: 'fill', source: SIM_SRC, paint: { 'fill-color': color, 'fill-opacity': opacity } }, before);
    else { map.setPaintProperty(L_SIMF, 'fill-color', color); map.setPaintProperty(L_SIMF, 'fill-opacity', opacity); }
    if (!map.getLayer(L_SIML)) map.addLayer({ id: L_SIML, type: 'line', source: SIM_SRC, paint: { 'line-color': '#111', 'line-width': 1 } }, before);
    S.simOn = true; S.sim = sim; S.simCycle = cycle;
  }
  function clearSim() {
    for (const l of [L_SIML, L_SIMF]) if (map && map.getLayer(l)) map.removeLayer(l);
    if (map && map.getSource(SIM_SRC)) map.removeSource(SIM_SRC);
    if (S.simPopup) S.simPopup.remove();
    S.simOn = false; S.sim = null;
  }
  function onSimClick(e) {
    const f = e.features && e.features[0]; if (!f || !S.sim) return;
    const s = (S.sim.stands || []).find(x => x.id === f.properties.stand_id); if (!s) return;
    const c = S.simCycle || 0, cy = s.cycles[c] || {}, m = k => (s.map[k] || [])[c];
    if (S.simPopup) S.simPopup.remove();
    S.simPopup = new maplibregl.Popup({ maxWidth: '300px' }).setLngLat(e.lngLat).setHTML(
      `<div class="stpop"><b>Stand ${esc(s.id)}</b> · cycle ${c}<table>
       <tr><td>area</td><td>${f1(s.area_ac)} ac</td></tr><tr><td>TreeMap plots</td><td>${s.nplots_sim} simulated of ${s.nplots}</td></tr>
       <tr><td>basal area</td><td>${f1(m('ba'))} ft²/ac</td></tr><tr><td>live carbon</td><td>${f1(m('carbon'))} t C/ha</td></tr>
       <tr><td>volume</td><td>${f1(m('volume'))} ft³/ac</td></tr><tr><td>removed volume</td><td>${f1(cy.removed_volume)} ft³ (stand total)</td></tr></table></div>`).addTo(map);
  }

  function init() {
    const host = document.getElementById('fiaGrp'); if (!host) return;
    const g = document.createElement('details'); g.className = 'grp fia'; g.id = 'stGrp';
    g.innerHTML = `<summary><h2>Stands — delineated (prototype)</h2>
      <label class="fswitch" title="Show delineated stand polygons"><input type="checkbox" id="stOn" role="switch" aria-label="Show stands layer"><span></span></label></summary>
      <div class="fbody"><div class="hint">Stand polygons segmented from TreeMap 2022, SOLUS soils and AlphaEarth 2022/2025 embeddings, split by roads,
      water and parcel lines; recent clear-cuts (Hansen loss since 2017) form their own stands. Click a stand for its attributes; <b>Simulate this
      stand</b> makes it the AOI.</div>
      <label class="field" for="stArea">Area</label><select id="stArea">${AREAS.map(a => `<option value="${a.id}">${esc(a.label)}</option>`).join('')}</select>
      <label class="field" for="stColor">Color by</label><select id="stColor"><option value="status">status</option><option value="BALIVE">basal area</option>
      <option value="STANDHT">stand height</option><option value="QMD">QMD</option><option value="CANOPYPCT">canopy cover</option>
      <optgroup label="Generative vegetation (prototype, TreeMap-independent)"><option value="GEN_BA">basal area — generative</option><option value="GEN_TPA">trees per acre — generative</option>
      <option value="GEN_QMD">QMD — generative</option><option value="GEN_BA_SD">basal area — sample spread</option></optgroup></select>
      <div id="stLegend" class="hint" style="margin-top:6px"></div><div class="hint" id="stStatus" role="status" aria-live="polite"></div></div>`;
    host.parentNode.insertBefore(g, host);
    $('stOn').addEventListener('change', e => toggle(e.target.checked));
    $('stColor').addEventListener('change', e => { S.color = e.target.value; if (map.getLayer(L_FILL)) map.setPaintProperty(L_FILL, 'fill-color', colorExpr()); legend(); });
    $('stArea').addEventListener('change', e => { S.area = AREAS.find(a => a.id === e.target.value); S.data = null; if (S.on) { removeLayers(); toggle(true); } });
    map.on('click', L_FILL, onClick);
    map.on('click', L_SIMF, onSimClick);
    $('vegSource').addEventListener('change', e => setVegSource(e.target.value));
    $('stSelPoly').addEventListener('click', () => startDraw(g => selectInPolygon(g)));
    $('stSelClear').addEventListener('click', clearSelection);
    map.on('mouseenter', L_FILL, () => map.getCanvas().style.cursor = 'pointer');
    map.on('mouseleave', L_FILL, () => map.getCanvas().style.cursor = '');
  }
  return { init, setSelection, showSim, clearSim, isSim: () => S.simOn };
})();
