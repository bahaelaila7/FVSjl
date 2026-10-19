'use strict';
/* Optional layer "Stands (delineated)": stand polygons from the standmap pipeline (/workspace/standmap: TreeMap + SOLUS + AlphaEarth
   2022/2025 segmentation, cut by roads / water / parcels, recent clear-cuts kept as their own stands). Static GeoJSON per prototype area
   (public/stands/<area>.geojson, EPSG:4326). Click a stand → its attributes; "Simulate this stand" makes the stand polygon the AOI, so the
   existing AOI pipeline (TreeMap pixels → FIA plots → FVSjl / Pan) runs on exactly that stand. */
const StandsUI = (() => {
  const SRC = 'stands', L_FILL = 'stands-fill', L_LINE = 'stands-line', L_SEL = 'stands-sel';
  const AREAS = [{ id: 'baker_fl', label: 'Baker County, FL (prototype)', url: '/stands/baker_fl.geojson' }];
  const STATUS = { 1: ['forest', '#2e7d32'], 2: ['harvested (loss since 2017)', '#e65100'], 3: ['regenerating', '#f9a825'], 4: ['forest (not in TreeMap)', '#6a1b9a'] };
  const S = { on: false, data: null, area: AREAS[0], popup: null, color: 'status' };
  const $ = id => document.getElementById(id);
  const esc = t => String(t).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const f1 = v => (v == null || Number.isNaN(+v)) ? '–' : (+v).toLocaleString(undefined, { maximumFractionDigits: 1 });
  const say = t => { const e = $('stStatus'); if (e) e.textContent = t; };

  function colorExpr() {
    if (S.color === 'status') return ['match', ['get', 'status'], 1, STATUS[1][1], 2, STATUS[2][1], 3, STATUS[3][1], 4, STATUS[4][1], '#888'];
    const ramps = { BALIVE: [0, 150], STANDHT: [0, 100], QMD: [0, 16], CANOPYPCT: [0, 100] };
    const [lo, hi] = ramps[S.color];
    return ['case', ['==', ['typeof', ['get', S.color]], 'number'],
      ['interpolate', ['linear'], ['get', S.color], lo, '#f7fcb9', (lo + hi) / 2, '#41ab5d', hi, '#004529'], '#bbbbbb'];
  }
  function legend() {
    const el = $('stLegend'); if (!el) return;
    if (S.color === 'status') { el.innerHTML = Object.values(STATUS).map(([l, c]) => `<span class="stkey"><i style="background:${c}"></i>${esc(l)}</span>`).join(''); return; }
    const u = { BALIVE: 'ft²/ac', STANDHT: 'ft', QMD: 'in', CANOPYPCT: '%' }[S.color];
    el.innerHTML = `<span class="stramp"></span> ${esc(S.color)} (${u}); grey = no TreeMap pixels (e.g. harvested)`;
  }
  async function load() {
    say('loading stands…');
    const r = await fetch(S.area.url); if (!r.ok) throw new Error('stands file not found: ' + S.area.url);
    S.data = await r.json(); return S.data;
  }
  function addLayers() {
    const before = map.getLayer('aoi-fill') ? 'aoi-fill' : undefined;          // below the AOI so a selected stand-AOI stays visible
    if (!map.getSource(SRC)) map.addSource(SRC, { type: 'geojson', data: S.data, promoteId: 'stand_id' });
    if (!map.getLayer(L_FILL)) map.addLayer({ id: L_FILL, type: 'fill', source: SRC, paint: { 'fill-color': colorExpr(), 'fill-opacity': 0.45 } }, before);
    if (!map.getLayer(L_LINE)) map.addLayer({ id: L_LINE, type: 'line', source: SRC, paint: { 'line-color': '#1b1b1b', 'line-width': ['interpolate', ['linear'], ['zoom'], 9, 0.2, 14, 1.2] } }, before);
    if (!map.getLayer(L_SEL)) map.addLayer({ id: L_SEL, type: 'line', source: SRC, filter: ['==', ['get', 'stand_id'], -1], paint: { 'line-color': '#00e5ff', 'line-width': 3 } }, before);
  }
  function removeLayers() { for (const l of [L_SEL, L_LINE, L_FILL]) if (map.getLayer(l)) map.removeLayer(l); if (map.getSource(SRC)) map.removeSource(SRC); if (S.popup) S.popup.remove(); }
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
      ['TreeMap pixels', f1(100 * (p.tm_frac ?? 0)) + ' % of the stand'], ['TreeMap plots', `${p.tm_n ?? 0} (main ${f1(100 * (p.tm_share ?? 0))} %)`],
      ['available water 0–30 cm', p.awc030 == null ? '–' : `${f1(300 * p.awc030)} mm (${(+p.awc030).toFixed(3)} cm/cm)`], ['clay / sand 0–30 cm', `${f1(p.clay030)} / ${f1(p.sand030)} %`]].filter(Boolean);
    return `<div class="stpop"><b>Stand ${p.stand_id}</b><table>${rows.map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(v)}</td></tr>`).join('')}</table>
      <button type="button" class="btn" id="stSim" ${p.status === 2 && !(p.tm_frac > 0) ? 'title="harvested stand: no TreeMap pixels; the AOI pipeline has no vegetation for it yet (donors pending)"' : ''}>Simulate this stand</button></div>`;
  }
  function onClick(e) {
    const f = e.features && e.features[0]; if (!f) return;
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
      addLayers(); map.fitBounds(bbox(S.data), { padding: 30, duration: 600 });
      const n = S.data.features.length, ha = S.data.features.reduce((a, f) => a + (f.properties.area_ha || 0), 0);
      say(`${n.toLocaleString()} stands · ${Math.round(ha).toLocaleString()} ha forest land`);
      legend();
    } catch (err) { say('error: ' + err.message); $('stOn').checked = false; S.on = false; }
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
      <option value="STANDHT">stand height</option><option value="QMD">QMD</option><option value="CANOPYPCT">canopy cover</option></select>
      <div id="stLegend" class="hint" style="margin-top:6px"></div><div class="hint" id="stStatus" role="status" aria-live="polite"></div></div>`;
    host.parentNode.insertBefore(g, host);
    $('stOn').addEventListener('change', e => toggle(e.target.checked));
    $('stColor').addEventListener('change', e => { S.color = e.target.value; if (map.getLayer(L_FILL)) map.setPaintProperty(L_FILL, 'fill-color', colorExpr()); legend(); });
    $('stArea').addEventListener('change', e => { S.area = AREAS.find(a => a.id === e.target.value); S.data = null; if (S.on) { removeLayers(); toggle(true); } });
    map.on('click', L_FILL, onClick);
    map.on('mouseenter', L_FILL, () => map.getCanvas().style.cursor = 'pointer');
    map.on('mouseleave', L_FILL, () => map.getCanvas().style.cursor = '');
  }
  return { init };
})();
