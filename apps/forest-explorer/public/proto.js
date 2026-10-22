'use strict';
/* "Prototype layers (research)": every layer listed in /proto/manifest.json, built by the standmap project (scripts/fe_proto_layers.py,
   config/fe_layers.yml): stand-delineation runs and tuning variants, segmentation inputs (land status, superpixels, AlphaEarth / embedding RGB),
   held-out evaluation layers, generative-vegetation stand maps, HBAD3P free-run errors at the public FIA plot coordinates.
   Layer types: fill (polygons, colour by a chosen attribute), circle (points), image (georeferenced PNG). Each layer: on/off, opacity,
   colour-by (vector), legend, click → attributes. Rebuild the layers with the standmap script, then restart the server (static files). */
const ProtoUI = (() => {
  const S = { man: null, state: new Map(), popup: null };
  const VIR = ['#440154', '#414487', '#2a788e', '#22a884', '#7ad151', '#fde725'];
  const DIV = ['#2166ac', '#67a9cf', '#d1e5f0', '#f7f7f7', '#fddbc7', '#ef8a62', '#b2182b'];
  const $ = id => document.getElementById(id);
  const esc = t => String(t).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const sid = L => 'proto-' + L.id, lid = (L, k) => 'proto-' + L.id + '-' + k;
  const fmt = v => (v == null || Number.isNaN(+v)) ? '–' : (typeof v === 'number' ? (+v).toLocaleString(undefined, { maximumFractionDigits: 3 }) : String(v));

  function propOf(L, name) { return (L.props || []).find(p => p.name === name); }
  function colorExpr(L, name) {
    const p = propOf(L, name); if (!p) return '#888';
    if (p.categories) {
      const e = ['match', ['to-string', ['get', p.name]]];
      for (const [k, v] of Object.entries(p.categories)) e.push(k, v.color);
      e.push('#888'); return e;
    }
    const stops = p.diverging ? DIV : VIR, [lo, hi] = p.range, e = ['interpolate', ['linear'], ['get', p.name]];
    stops.forEach((c, i) => e.push(lo + (hi - lo) * i / (stops.length - 1), c));
    return ['case', ['==', ['typeof', ['get', p.name]], 'number'], e, '#bdbdbd'];
  }
  function rampCss(stops) { return `linear-gradient(90deg, ${stops.join(', ')})`; }
  function legendHtml(L, st) {
    if (L.type === 'image') return (L.legend || []).map(g => g.ramp
      ? `<span class="pxramp" style="background:${rampCss(g.ramp)}"></span> ${fmt(g.lo)} … ${fmt(g.hi)} ${esc(g.unit || '')}`
      : g.color ? `<span class="stkey"><i style="background:${g.color}"></i>${esc(g.label)}</span>` : `<span>${esc(g.label)}</span>`).join(' ');
    const p = propOf(L, st.color); if (!p) return '';
    if (p.categories) return Object.values(p.categories).map(v => `<span class="stkey"><i style="background:${v.color}"></i>${esc(v.label)}</span>`).join('');
    return `<span class="pxramp" style="background:${rampCss(p.diverging ? DIV : VIR)}"></span> ${fmt(p.range[0])} … ${fmt(p.range[1])} ${esc(p.unit || '')}; grey = no value`;
  }
  function beforeId(L) {                       // images under every proto vector layer; everything under the AOI outline
    if (L.type === 'image') {
      const v = (map.getStyle().layers || []).find(l => l.id.startsWith('proto-') && !l.id.endsWith('-img'));
      if (v) return v.id;
    }
    return map.getLayer('aoi-fill') ? 'aoi-fill' : undefined;
  }
  function add(L) {
    const st = S.state.get(L.id), b = beforeId(L);
    if (L.type === 'image') {
      map.addSource(sid(L), { type: 'image', url: L.url, coordinates: L.coordinates });
      map.addLayer({ id: lid(L, 'img'), type: 'raster', source: sid(L), paint: { 'raster-opacity': st.opacity, 'raster-resampling': 'nearest' } }, b);
      return;
    }
    map.addSource(sid(L), { type: 'geojson', data: L.url });
    if (L.type === 'fill') {
      map.addLayer({ id: lid(L, 'f'), type: 'fill', source: sid(L), paint: { 'fill-color': colorExpr(L, st.color), 'fill-opacity': st.opacity } }, b);
      map.addLayer({ id: lid(L, 'l'), type: 'line', source: sid(L), paint: { 'line-color': '#1b1b1b', 'line-width': ['interpolate', ['linear'], ['zoom'], 9, 0.15, 14, 1] } }, b);
      map.on('click', lid(L, 'f'), e => onClick(L, e));
    } else {
      map.addLayer({ id: lid(L, 'c'), type: 'circle', source: sid(L), paint: { 'circle-color': colorExpr(L, st.color), 'circle-opacity': st.opacity,
        'circle-radius': ['interpolate', ['linear'], ['zoom'], 5, 3, 10, 6], 'circle-stroke-color': '#222', 'circle-stroke-width': 0.6 } }, b);
      map.on('click', lid(L, 'c'), e => onClick(L, e));
    }
  }
  function remove(L) {
    for (const k of ['img', 'f', 'l', 'c']) if (map.getLayer(lid(L, k))) map.removeLayer(lid(L, k));
    if (map.getSource(sid(L))) map.removeSource(sid(L));
  }
  function restyle(L) {
    const st = S.state.get(L.id);
    if (map.getLayer(lid(L, 'img'))) map.setPaintProperty(lid(L, 'img'), 'raster-opacity', st.opacity);
    if (map.getLayer(lid(L, 'f'))) { map.setPaintProperty(lid(L, 'f'), 'fill-color', colorExpr(L, st.color)); map.setPaintProperty(lid(L, 'f'), 'fill-opacity', st.opacity); }
    if (map.getLayer(lid(L, 'c'))) { map.setPaintProperty(lid(L, 'c'), 'circle-color', colorExpr(L, st.color)); map.setPaintProperty(lid(L, 'c'), 'circle-opacity', st.opacity); }
    const lg = $('pxlg-' + L.id); if (lg) lg.innerHTML = legendHtml(L, st);
  }
  function onClick(L, e) {
    const f = e.features && e.features[0]; if (!f) return;
    if (typeof StandsUI !== 'undefined' && StandsUI.isSim && StandsUI.isSim()) return;
    const rows = Object.entries(f.properties).map(([k, v]) => { const p = propOf(L, k);
      const lab = p ? `${p.label}${p.unit ? ' (' + p.unit + ')' : ''}` : k; return `<tr><td>${esc(lab)}</td><td>${esc(fmt(v))}</td></tr>`; }).join('');
    if (S.popup) S.popup.remove();
    S.popup = new maplibregl.Popup({ maxWidth: '340px' }).setLngLat(e.lngLat).setHTML(`<div class="stpop"><b>${esc(L.label)}</b><table>${rows}</table></div>`).addTo(map);
  }
  function fitTo(L) {
    if (L.type === 'image') { const c = L.coordinates; map.fitBounds([[c[0][0], c[2][1]], [c[1][0], c[0][1]]], { padding: 30, duration: 600 }); return; }
    fetch(L.url).then(r => r.json()).then(fc => { let w = 180, s = 90, e = -180, n = -90;
      const walk = c => { if (typeof c[0] === 'number') { w = Math.min(w, c[0]); e = Math.max(e, c[0]); s = Math.min(s, c[1]); n = Math.max(n, c[1]); } else c.forEach(walk); };
      fc.features.forEach(f => walk(f.geometry.coordinates)); map.fitBounds([[w, s], [e, n]], { padding: 30, duration: 600 }); }).catch(() => {});
  }
  function rowHtml(L) {
    const st = S.state.get(L.id);
    const sel = (L.props && L.props.length) ? `<select id="pxc-${esc(L.id)}" class="pxsel">${L.props.map(p => `<option value="${esc(p.name)}"${p.name === st.color ? ' selected' : ''}>${esc(p.label)}</option>`).join('')}</select>` : '';
    return `<div class="pxrow"><label class="pxlab"><input type="checkbox" id="pxo-${esc(L.id)}"> ${esc(L.label)}</label>
      <div class="pxctl" id="pxctl-${esc(L.id)}" hidden>${sel}<label class="pxop">opacity <input type="range" min="5" max="100" value="${Math.round(st.opacity * 100)}" id="pxa-${esc(L.id)}"></label>
      <button type="button" class="btn pxfit" id="pxf-${esc(L.id)}" title="Zoom to the layer">zoom</button>
      <div class="hint" id="pxlg-${esc(L.id)}">${legendHtml(L, st)}</div>${L.note ? `<div class="hint">${esc(L.note)}</div>` : ''}</div></div>`;
  }
  async function init() {
    const host = document.getElementById('fiaGrp'); if (!host) return;
    let man;
    try { const r = await fetch('/proto/manifest.json', { cache: 'no-cache' }); if (!r.ok) return; man = await r.json(); } catch (err) { return; }
    S.man = man;
    for (const L of man.layers) S.state.set(L.id, { on: false, opacity: L.type === 'image' ? 0.85 : (L.type === 'circle' ? 0.95 : 0.55), color: L.color || null });
    const groups = [...new Set(man.layers.map(L => L.group))];
    const g = document.createElement('details'); g.className = 'grp fia'; g.id = 'pxGrp';
    g.innerHTML = `<summary><h2>Prototype layers (research)</h2></summary><div class="fbody">
      <div class="hint">Results of the standmap / HBAD3P prototypes (built ${esc(man.built || '?')}). Rebuild with <code>standmap/scripts/fe_proto_layers.py</code>, then restart the server.</div>
      ${groups.map(gr => `<h3>${esc(gr)}</h3>` + man.layers.filter(L => L.group === gr).map(rowHtml).join('')).join('')}</div>`;
    host.parentNode.insertBefore(g, host);
    for (const L of man.layers) {
      $('pxo-' + L.id).addEventListener('change', e => { const st = S.state.get(L.id); st.on = e.target.checked; $('pxctl-' + L.id).hidden = !st.on;
        try { if (st.on) add(L); else remove(L); } catch (err) { console.error(err); setStatus('prototype layer error: ' + err.message); } });
      $('pxa-' + L.id).addEventListener('input', e => { S.state.get(L.id).opacity = (+e.target.value) / 100; restyle(L); });
      $('pxf-' + L.id).addEventListener('click', () => fitTo(L));
      const c = (L.props && L.props.length) ? $('pxc-' + L.id) : null; if (c) c.addEventListener('change', e => { S.state.get(L.id).color = e.target.value; restyle(L); });
    }
  }
  return { init };
})();
