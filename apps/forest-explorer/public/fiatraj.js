'use strict';
/* Optional layer "FIA plots — tree trajectories": dots for FIA plot coordinates (filtered), click -> per-tree trajectory
   chart (SVG) for that plot. Backend: /api/fiaplots/{options,points,trajectories} (src/FiaTrajectories.jl).
   Mirrors the options and plotting semantics of FIA_DATA_PREP/tree_trajectories.py. */

const FiaUI = (() => {
  const MIN_ZOOM = 6;                       // below this the dot layer is not requested (too many plots)
  const SRC = 'fia-plots', L_DOTS = 'fia-dots', L_SEL = 'fia-sel';
  const FUZZ_NOTE = 'Public FIA plot coordinates are fuzzed (up to ~0.5–1 mile; some private plots are swapped), so a dot is only approximately where the plot is.';
  const S = { opts: null, on: false, sel: null, seq: 0, pseq: 0, loading: null, popup: null, last: null, ctl: null };
  const $ = id => document.getElementById(id);
  const esc = t => String(t).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const debounce = (fn, ms) => { let t; return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); }; };
  const cssv = (n, d) => (getComputedStyle(document.documentElement).getPropertyValue(n).trim() || d);
  const say = t => { const e = $('fiaStatus'); if (e) e.textContent = t; };
  const post = (url, body, signal) => fetch(url, { method: 'POST', headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(body), signal }).then(async r => { const j = await r.json().catch(() => ({})); if (!r.ok || j.error) throw new Error(j.error || r.statusText); return j; });

  /* ---------------- state read from the accordion ---------------- */
  const checked = name => [...document.querySelectorAll(`#fiaGrp input[data-f="${name}"]:checked`)].map(i => i.value);
  let species = [];                                   // chip multiselect
  function filters() {
    const own = checked('own').map(Number);
    return { species, land_use: checked('land_use'), cut_event: checked('cut_event'), fate: checked('fate'), own,
             swhw: $('fiaSwhw').value, cycle: $('fiaCycle').value, eco_field: $('fiaEcoField').value, eco_val: $('fiaEcoVal').dataset.ok || '' };
  }
  function options() {
    return { x: $('fiaX').value, y: $('fiaY').value, color_by: $('fiaColor').value, max_traj: Math.max(1, Math.min(5000, +$('fiaMax').value || 300)),
             repr: $('fiaRepr').value, imputed: $('fiaImp').checked, derivative: $('fiaDeriv').checked,
             logx: $('fiaLogX').checked, logy: $('fiaLogY').checked, fit: $('fiaFit').checked };
  }

  /* ---------------- build the accordion controls from /options ---------------- */
  function fillSelect(sel, items, val) { sel.innerHTML = ''; items.forEach(([v, l]) => { const o = document.createElement('option'); o.value = v; o.textContent = l; sel.appendChild(o); }); if (val != null) sel.value = val; }
  function checkGroup(id, name, items) {
    $(id).innerHTML = items.map(([v, l]) => `<label class="chk"><input type="checkbox" data-f="${name}" value="${esc(v)}"> ${esc(l)}</label>`).join('');
  }
  function renderChips() {
    $('fiaChips').innerHTML = species.map(s => `<span class="fchip">${esc(s)}<button type="button" data-sp="${esc(s)}" aria-label="remove ${esc(s)}">×</button></span>`).join('');
  }
  function setEcoValues() {
    const f = S.opts.eco_fields.find(e => e.id === $('fiaEcoField').value);
    $('fiaEcoList').innerHTML = f.values.map(v => `<option value="${esc(v)}">`).join('');
    $('fiaEcoVal').value = ''; $('fiaEcoVal').dataset.ok = '';
  }
  async function loadOptions() {
    if (S.opts) return S.opts;
    if (S.loading) return S.loading;
    say('loading plot index (first use takes a few seconds)…');
    S.loading = fetch('/api/fiaplots/options').then(async r => { const j = await r.json(); if (!r.ok || j.error) throw new Error(j.error || r.statusText); return j; })
      .then(o => {
        S.opts = o;
        $('fiaSpList').innerHTML = o.species.map(s => `<option value="${esc(s.id)}">${esc(s.id)} · ${s.plots} plots</option>`).join('');
        checkGroup('fiaLand', 'land_use', o.land_use.map(v => [v, v]));
        checkGroup('fiaCut', 'cut_event', o.cut_event.map(v => [v, v.replace(/_/g, ' ')]));
        checkGroup('fiaFate', 'fate', o.fate.map(v => [v, v.replace(/_/g, ' ')]));
        checkGroup('fiaOwn', 'own', o.own.map(x => [x.id, x.label]));
        fillSelect($('fiaSwhw'), [['', '(all)'], ...o.swhw.map(v => [v, v === 'S' ? 'Softwood (S)' : v === 'H' ? 'Hardwood (H)' : v])], '');
        fillSelect($('fiaCycle'), [['', '(all)'], ...o.cycle.map(c => [c, String(c)])], '');
        fillSelect($('fiaEcoField'), o.eco_fields.filter(e => e.values.length).map(e => [e.id, e.label]), 'epa_l3');
        const ax = o.axes.map(a => [a.id, a.label]);
        fillSelect($('fiaX'), ax, o.defaults.x); fillSelect($('fiaY'), ax, o.defaults.y);
        fillSelect($('fiaColor'), o.color_by.map(c => [c, c]), o.defaults.color_by);
        fillSelect($('fiaRepr'), [['all', 'All points (trajectories)'], ['first', 'First per tree'], ['last', 'Last per tree']], o.defaults.repr);
        $('fiaMax').value = o.defaults.max_traj; $('fiaImp').checked = o.defaults.imputed;
        setEcoValues();
        say(`${o.n_plots.toLocaleString()} FIA plots with tree trajectories. Zoom in to level ${MIN_ZOOM}+ to see them.`);
        return o;
      }).catch(e => { S.loading = null; say('FIA plot layer unavailable: ' + e.message); throw e; });
    return S.loading;
  }

  /* ---------------- map layer ---------------- */
  const emptyFC = () => ({ type: 'FeatureCollection', features: [] });
  function paintColors() {
    if (!map.getLayer(L_DOTS)) return;
    map.setPaintProperty(L_DOTS, 'circle-color', cssv('--fia-dot', '#fbbf24'));
    map.setPaintProperty(L_DOTS, 'circle-stroke-color', cssv('--fia-dot-line', '#1a1200'));
    map.setPaintProperty(L_SEL, 'circle-stroke-color', cssv('--fia-sel', '#ffffff'));
  }
  function ensureLayers() {
    if (map.getSource(SRC)) return;
    map.addSource(SRC, { type: 'geojson', data: emptyFC() });
    map.addLayer({ id: L_DOTS, type: 'circle', source: SRC, paint: {
      'circle-radius': ['interpolate', ['linear'], ['get', 'v'], 1, 3, 6, 6.5],       // modest size by number of visits
      'circle-color': '#fbbf24', 'circle-opacity': 0.9, 'circle-stroke-width': 1, 'circle-stroke-color': '#1a1200' } });
    map.addLayer({ id: L_SEL, type: 'circle', source: SRC, filter: ['==', ['get', 'key'], ''], paint: {
      'circle-radius': 10, 'circle-color': 'rgba(0,0,0,0)', 'circle-stroke-width': 2.5, 'circle-stroke-color': '#fff' } });
    paintColors();
    map.on('mouseenter', L_DOTS, () => { if (!drawing) map.getCanvas().style.cursor = 'pointer'; });
    map.on('mouseleave', L_DOTS, () => { if (!drawing) map.getCanvas().style.cursor = ''; hidePopup(); });
    map.on('mousemove', L_DOTS, e => { const f = e.features && e.features[0]; if (f && !drawing) showPopup(f, e.lngLat); });
    map.on('click', L_DOTS, e => { if (drawing) return; const f = e.features && e.features[0]; if (f) openPlot(f.properties.key); });
  }
  function hidePopup() { if (S.popup) { S.popup.remove(); S.popup = null; } }
  function showPopup(f, ll) {
    const p = f.properties;
    const html = `<b>${esc(p.key)}</b><br>${p.n_trees} trees · ${p.n_visits} visits<br>${esc(p.species)}<br>` +
      `<span style="color:var(--muted)" title="${esc(FUZZ_NOTE)}">approximate location: public FIA coordinates are fuzzed</span><br><i>click for trajectories</i>`;
    if (!S.popup) S.popup = new maplibregl.Popup({ closeButton: false, closeOnClick: false, offset: 10, maxWidth: '260px' });
    S.popup.setLngLat(f.geometry.coordinates).setHTML(html).addTo(map);
  }
  async function refreshDots() {
    if (!S.on || !map.getSource(SRC)) return;
    if (map.getZoom() < MIN_ZOOM) { map.getSource(SRC).setData(emptyFC()); say(`Zoom in to level ${MIN_ZOOM}+ to show plots (now ${map.getZoom().toFixed(1)}).`); return; }
    const b = map.getBounds(), my = ++S.pseq;
    if (S.ctl) S.ctl.abort();
    S.ctl = new AbortController();
    try {
      const r = await post('/api/fiaplots/points', { west: b.getWest(), south: b.getSouth(), east: b.getEast(), north: b.getNorth(), ...filters() }, S.ctl.signal);
      if (my !== S.pseq) return;
      map.getSource(SRC).setData({ type: 'FeatureCollection', features: r.points.map(p => ({ type: 'Feature',
        geometry: { type: 'Point', coordinates: [p.lon, p.lat] },
        properties: { key: p.key, n_trees: p.n_trees, n_visits: p.n_visits, v: p.n_visits, species: p.species.join(', ') } })) });
      say(r.truncated ? `${r.n_match.toLocaleString()} plots match in view; showing an even sample of ${r.shown.toLocaleString()} — zoom in or add filters.`
                      : `${r.n_match.toLocaleString()} plot${r.n_match === 1 ? '' : 's'} match in view.`);
    } catch (e) { if (e.name !== 'AbortError') say('plots: ' + e.message); }
  }
  const refreshDotsSoon = debounce(refreshDots, 350);

  async function setOn(on) {
    S.on = on;
    if (on) {
      $('fiaGrp').open = true;
      try { await loadOptions(); } catch (_) { $('fiaOn').checked = false; S.on = false; return; }
      ensureLayers();
      map.setLayoutProperty(L_DOTS, 'visibility', 'visible'); map.setLayoutProperty(L_SEL, 'visibility', 'visible');
      refreshDots();
    } else {
      hidePopup();
      if (map.getLayer(L_DOTS)) { map.setLayoutProperty(L_DOTS, 'visibility', 'none'); map.setLayoutProperty(L_SEL, 'visibility', 'none'); }
      closePanel(); say('Layer off.');
    }
  }

  /* ---------------- trajectory panel ---------------- */
  function closePanel() {
    $('fiaPanel').hidden = true; S.sel = null; S.last = null; S.seq++;
    if (map.getLayer(L_SEL)) map.setFilter(L_SEL, ['==', ['get', 'key'], '']);
  }
  function openPlot(key) {
    S.sel = key; $('fiaPanel').hidden = false;
    if (map.getLayer(L_SEL)) map.setFilter(L_SEL, ['==', ['get', 'key'], key]);
    $('fiaTitle').textContent = 'Plot ' + key;
    renderChart();
  }
  async function renderChart() {
    if (!S.sel) return;
    const my = ++S.seq;
    $('fiaPlotInfo').innerHTML = '<span class="hint">loading…</span>';
    $('fiaLegend').innerHTML = '';
    try {
      const r = await post('/api/fiaplots/trajectories', { plot: S.sel, filters: filters(), options: options() });
      if (my !== S.seq) return;
      S.last = r; drawInfo(r); drawChart(r);
    } catch (e) {
      if (my !== S.seq) return;
      $('fiaChart').innerHTML = ''; $('fiaPlotInfo').innerHTML = `<span class="err">${esc(e.message)}</span>`;
    }
  }
  const renderChartSoon = debounce(renderChart, 250);

  function drawInfo(r) {
    const p = r.plot, f = n => Number(n).toLocaleString();
    $('fiaPlotInfo').innerHTML =
      `<div class="fkv"><span>Location</span><b title="${esc(FUZZ_NOTE)}">${p.lat.toFixed(4)}, ${p.lon.toFixed(4)} <i class="fz">(fuzzed)</i></b>` +
      `<span>Visits</span><b>${p.n_visits} (${p.first_year}–${p.last_year})</b>` +
      `<span>Trees</span><b>${f(p.n_trees)} on plot · ${f(r.n_match)} match filters</b>` +
      `<span>Ecoregion</span><b>${esc(p.epa_l3)} · ${esc(p.epa_l4)}${p.ecosubcd ? ' · ' + esc(p.ecosubcd) : ''}</b>` +
      `<span>Land use</span><b>${esc(p.land_use.join(', ') || '–')}</b>` +
      `<span>Ownership</span><b>${esc(p.ownership.join(', ') || '–')}</b>` +
      `<span>Cut events</span><b>${esc(p.cut_events.map(v => v.replace(/_/g, ' ')).join(', ') || 'none')}</b></div>`;
  }

  /* ---- SVG chart ---- */
  const NS = 'http://www.w3.org/2000/svg';
  function fmtTick(v) {
    if (v === 0) return '0';
    const a = Math.abs(v);
    if (a >= 1e6) return +(v / 1e6).toPrecision(3) + 'M';
    if (a >= 1e4) return +(v / 1e3).toPrecision(3) + 'k';
    if (a < 0.01) return v.toExponential(0);
    return String(+v.toPrecision(4));
  }
  function niceTicks(lo, hi, n = 6) {
    const span = hi - lo || 1, raw = span / n, p = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / p;
    const step = (m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10) * p, t = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + step * 1e-9; v += step) t.push(+v.toPrecision(12));
    return t;
  }
  function logTicks(lo, hi) {
    const a = Math.floor(Math.log10(lo)), b = Math.ceil(Math.log10(hi)), t = [];
    const mult = (b - a) <= 2 ? [1, 2, 5] : [1];
    for (let e = a; e <= b; e++) for (const m of mult) { const v = m * Math.pow(10, e); if (v >= lo * 0.999 && v <= hi * 1.001) t.push(v); }
    return t;
  }
  function marker(g, shape, cx, cy, color, filled) {
    const mk = (tag, at) => { const e = document.createElementNS(NS, tag); for (const k in at) e.setAttribute(k, at[k]); g.appendChild(e); return e; };
    const fill = filled ? color : 'var(--panel)', st = { stroke: color, 'stroke-width': 1.6, fill };
    const r = 4.5;
    if (shape === 'circle') mk('circle', { cx, cy, r, ...st });
    else if (shape === 'square') mk('rect', { x: cx - r, y: cy - r, width: 2 * r, height: 2 * r, ...st });
    else if (shape === 'diamond') mk('path', { d: `M${cx} ${cy - r - 1}L${cx + r + 1} ${cy}L${cx} ${cy + r + 1}L${cx - r - 1} ${cy}Z`, ...st });
    else if (shape === 'triangle-down') mk('path', { d: `M${cx - r - 1} ${cy - r}L${cx + r + 1} ${cy - r}L${cx} ${cy + r + 1}Z`, ...st });
    else mk('path', { d: `M${cx - r} ${cy - r}L${cx + r} ${cy + r}M${cx - r} ${cy + r}L${cx + r} ${cy - r}`, stroke: color, 'stroke-width': 2.2, fill: 'none' });
  }
  function drawChart(r) {
    const svg = $('fiaChart'); svg.innerHTML = '';
    const W = Math.max(320, svg.parentNode.clientWidth || 600), H = Math.min(440, Math.max(300, Math.round(W * 0.62)));
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`); svg.setAttribute('height', H);
    const M = { l: 64, r: 14, t: 12, b: 46 }, pw = W - M.l - M.r, ph = H - M.t - M.b;
    const el = (tag, at, parent = svg, txt) => { const e = document.createElementNS(NS, tag); for (const k in at) e.setAttribute(k, at[k]); if (txt != null) e.textContent = txt; parent.appendChild(e); return e; };
    const okv = (v, log) => Number.isFinite(v) && (!log || v > 0);
    // collect drawable points (log axes drop non-positive values)
    const segPts = (xs, ys) => { const out = []; for (let i = 0; i < xs.length; i++) if (okv(xs[i], r.logx) && okv(ys[i], r.logy)) out.push([xs[i], ys[i]]); return out; };
    const items = []; const xv = [], yv = [];
    for (const t of r.trees) {
      const it = { t, segs: [], dot: null, imp: [], term: null };
      if (t.dot) { if (okv(t.dotxy[0], r.logx) && okv(t.dotxy[1], r.logy)) it.dot = t.dotxy; }
      else t.segs.forEach(s => { const p = segPts(s.x, s.y); if (p.length) it.segs.push({ dash: s.dash, p }); });
      it.imp = (t.imp || []).filter(p => okv(p[0], r.logx) && okv(p[1], r.logy));
      if (t.term && okv(t.term.x, r.logx) && okv(t.term.y, r.logy)) it.term = t.term;
      [...it.segs.flatMap(s => s.p), ...(it.dot ? [it.dot] : [])].forEach(p => { xv.push(p[0]); yv.push(p[1]); });
      items.push(it);
    }
    if (r.fit) segPts(r.fit.x, r.fit.y).forEach(p => { xv.push(p[0]); yv.push(p[1]); });
    if (!xv.length) {
      el('text', { x: W / 2, y: H / 2, 'text-anchor': 'middle', class: 'ax-txt' }, svg, r.n_match ? 'No plottable points for this axis / log selection.' : 'No trees on this plot match the filters.');
      $('fiaLegend').innerHTML = ''; return;
    }
    const dom = (v, log) => {
      let lo = Math.min(...v), hi = Math.max(...v);
      if (log) { lo = Math.log10(lo); hi = Math.log10(hi); }
      if (hi === lo) { lo -= 0.5; hi += 0.5; }
      const pad = (hi - lo) * 0.04; lo -= pad; hi += pad;
      return log ? [Math.pow(10, lo), Math.pow(10, hi)] : [lo, hi];
    };
    const [x0, x1] = dom(xv, r.logx), [y0, y1] = dom(yv, r.logy);
    const sc = (v, a, b, len, log, flip) => { const u = log ? (Math.log10(v) - Math.log10(a)) / (Math.log10(b) - Math.log10(a)) : (v - a) / (b - a); return flip ? len - u * len : u * len; };
    const X = v => M.l + sc(v, x0, x1, pw, r.logx, false), Y = v => M.t + sc(v, y0, y1, ph, r.logy, true);
    // grid + ticks
    const xt = r.logx ? logTicks(x0, x1) : niceTicks(x0, x1), yt = r.logy ? logTicks(y0, y1) : niceTicks(y0, y1, 5);
    xt.forEach(v => { const x = X(v); el('line', { x1: x, x2: x, y1: M.t, y2: M.t + ph, class: 'grid' }); el('line', { x1: x, x2: x, y1: M.t + ph, y2: M.t + ph + 4, class: 'axl' }); el('text', { x, y: M.t + ph + 17, 'text-anchor': 'middle', class: 'ax-txt' }, svg, fmtTick(v)); });
    yt.forEach(v => { const y = Y(v); el('line', { x1: M.l, x2: M.l + pw, y1: y, y2: y, class: 'grid' }); el('line', { x1: M.l - 4, x2: M.l, y1: y, y2: y, class: 'axl' }); el('text', { x: M.l - 7, y: y + 3.5, 'text-anchor': 'end', class: 'ax-txt' }, svg, fmtTick(v)); });
    el('rect', { x: M.l, y: M.t, width: pw, height: ph, class: 'frame', fill: 'none' });
    el('text', { x: M.l + pw / 2, y: H - 8, 'text-anchor': 'middle', class: 'ax-title' }, svg, r.x_title);
    el('text', { transform: `translate(14 ${M.t + ph / 2}) rotate(-90)`, 'text-anchor': 'middle', class: 'ax-title' }, svg, r.y_title);
    const clip = el('clipPath', { id: 'fiaClip' }); el('rect', { x: M.l, y: M.t, width: pw, height: ph }, clip);
    const gAll = el('g', { 'clip-path': 'url(#fiaClip)' });
    for (const it of items) {
      const c = it.t.color, tip = `tree ${it.t.id}` + (it.t.lv ? ` · ${it.t.lv}` : '');
      const g = el('g', {}, gAll); el('title', {}, g, tip);
      it.segs.forEach(s => {
        if (s.p.length === 1) { el('circle', { cx: X(s.p[0][0]), cy: Y(s.p[0][1]), r: 2.4, fill: c, opacity: 0.85 }, g); return; }
        el('polyline', { points: s.p.map(p => X(p[0]).toFixed(1) + ',' + Y(p[1]).toFixed(1)).join(' '), fill: 'none', stroke: c, 'stroke-width': 1.5,
                         'stroke-opacity': 0.85, 'stroke-dasharray': s.dash ? '5 4' : null, 'stroke-linejoin': 'round' }, g);
      });
      if (it.dot) el('circle', { cx: X(it.dot[0]), cy: Y(it.dot[1]), r: r.repr === 'all' ? 3.4 : 3.2, fill: c, 'fill-opacity': 0.85 }, g);
      it.imp.forEach(p => el('circle', { cx: X(p[0]), cy: Y(p[1]), r: 2.8, fill: 'none', stroke: c, 'stroke-width': 1, 'stroke-opacity': 0.8 }, g));
      if (it.term) marker(g, it.term.shape, X(it.term.x), Y(it.term.y), c, it.term.filled);
    }
    if (r.fit) {
      const p = segPts(r.fit.x, r.fit.y);
      el('polyline', { points: p.map(q => X(q[0]).toFixed(1) + ',' + Y(q[1]).toFixed(1)).join(' '), fill: 'none', class: 'fitline', 'stroke-dasharray': '7 5' }, gAll);
      const lx = r.fit.logx ? 'log₁₀(x)' : 'x', ly = r.fit.logy ? 'log₁₀(y)' : 'y';
      const t1 = `${ly} = ${(+r.fit.slope.toPrecision(3))}·${lx} ${r.fit.intercept < 0 ? '−' : '+'} ${Math.abs(+r.fit.intercept.toPrecision(3))}`;
      const t2 = `R² = ${r.fit.r2 == null ? 'n/a' : r.fit.r2.toFixed(3)}   n = ${r.fit.n.toLocaleString()}`;
      const bw = Math.max(t1.length, t2.length) * 6.4 + 14;
      el('rect', { x: M.l + pw - bw - 6, y: M.t + ph - 42, width: bw, height: 36, rx: 5, class: 'fitbox' });
      el('text', { x: M.l + pw - 13, y: M.t + ph - 26, 'text-anchor': 'end', class: 'fit-txt' }, svg, t1);
      el('text', { x: M.l + pw - 13, y: M.t + ph - 12, 'text-anchor': 'end', class: 'fit-txt' }, svg, t2);
    }
    drawLegend(r);
  }
  function drawLegend(r) {
    const sw = r.levels.map(l => `<span><i class="sw" style="background:${l.color}"></i>${esc(l.level)}</span>`).join('');
    let note = r.repr === 'all'
      ? (r.derivative ? 'Lines: growth rate Δy/Δx between successive visits of each tree.'
        : 'Lines = contiguous trajectories · <b>dashed</b> = phantom back-fill · dots = single-measurement trees' +
          (r.imputed ? ' · open points = imputed (phantom / interpolated)' : '') +
          ' · terminal marker: ○ alive · ▽ dead-natural · ✕ dead-disturbance · ▣ harvest (filled = explicit) · ◇ disappeared')
      : `Dots = each tree's <b>${esc(r.repr)}</b> real measurement (one point per tree).`;
    if (r.fit) note += ' · dashed black line = OLS fit in the displayed axis space.';
    const warn = r.sampled ? `<div class="warnl">${r.n_match.toLocaleString()} trees match — showing a random sample of ${r.n_plotted.toLocaleString()}. Raise “Max trajectories” to see more.</div>` : '';
    $('fiaLegend').innerHTML = `<div class="swatches" title="colour by ${esc(r.color_by)}"><b style="color:var(--muted)">${esc(r.color_by)}:</b>${sw}</div>` +
      `<div class="hint" style="margin-top:6px">${r.n_plotted.toLocaleString()} of ${r.n_match.toLocaleString()} trees · ${note}</div>${warn}`;
  }

  /* ---------------- wiring ---------------- */
  function onFilterChange() { refreshDotsSoon(); if (S.sel) renderChartSoon(); }
  function init() {
    const g = $('fiaGrp'); if (!g) return;
    $('fiaOn').onclick = e => e.stopPropagation();
    $('fiaOn').onchange = e => setOn(e.target.checked);
    $('fiaSwitch').onclick = e => e.stopPropagation();
    $('fiaClose').onclick = closePanel;
    document.addEventListener('keydown', e => { if (e.key === 'Escape' && !$('fiaPanel').hidden && !drawing) closePanel(); });
    // filters (re-filter dots + re-render the open chart)
    g.addEventListener('change', e => {
      const t = e.target;
      if (t.id === 'fiaOn') return;
      if (t.id === 'fiaSp') {
        const vals = t.value.split(',').map(s => s.trim()).filter(Boolean);
        let add = false;
        vals.forEach(v => { if (S.opts && S.opts.species.some(s => s.id === v) && !species.includes(v)) { species.push(v); add = true; } });
        t.value = ''; if (add) { renderChips(); onFilterChange(); } return;
      }
      if (t.id === 'fiaEcoField') { setEcoValues(); onFilterChange(); return; }
      if (t.id === 'fiaEcoVal') {
        const f = S.opts.eco_fields.find(e2 => e2.id === $('fiaEcoField').value);
        const v = t.value.trim();
        if (!v) { t.dataset.ok = ''; t.classList.remove('bad'); } else if (f.values.includes(v)) { t.dataset.ok = v; t.classList.remove('bad'); } else { t.dataset.ok = ''; t.classList.add('bad'); return; }
        onFilterChange(); return;
      }
      if (t.matches('[data-f],#fiaSwhw,#fiaCycle')) { onFilterChange(); return; }
      if (t.matches('.fopt')) { if (S.sel) renderChartSoon(); }        // chart options only affect the open chart
    });
    $('fiaMax').addEventListener('input', () => { if (S.sel) renderChartSoon(); });
    $('fiaChips').addEventListener('click', e => { const b = e.target.closest('button[data-sp]'); if (!b) return; species = species.filter(s => s !== b.dataset.sp); renderChips(); onFilterChange(); });
    $('fiaReset').onclick = () => {
      species = []; renderChips();
      g.querySelectorAll('input[data-f]').forEach(i => i.checked = false);
      $('fiaSwhw').value = ''; $('fiaCycle').value = ''; if (S.opts) setEcoValues(); $('fiaEcoVal').classList.remove('bad');
      onFilterChange();
    };
    window.addEventListener('resize', debounce(() => { if (S.last) drawChart(S.last); }, 150));
    const mq = window.matchMedia && window.matchMedia('(prefers-color-scheme: light)');
    if (mq && mq.addEventListener) mq.addEventListener('change', () => paintColors());
    map.on('moveend', refreshDotsSoon);
    new MutationObserver(() => paintColors()).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme', 'class'] });
  }
  return { init, _state: S };
})();
