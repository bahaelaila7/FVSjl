'use strict';
/* Pan engine family support: engine selector, parameter-set panel, harvest action,
   scenario engine tags, request building, labels.
   Loaded BEFORE app.js; everything here touches app.js globals (el, fmt, setStatus, plan,
   scenarios, renderActions ...) only at call time. */

const PAN = {
  meta: null,                       // /api/meta .pan
  engines: [{ id: 'fvs', label: 'FVSjl', plugin: null, climate: false, enabled: true }],
  engine: 'fvs',
  paramsByPlugin: {},               // plugin -> info (null/absent = default not loaded yet)
  store: {},                        // params_id -> info (incl. yaml text) for bundle export
  climate: null,                    // climate object, null = default
  pending: null,                    // mismatching upload awaiting engine switch
  samples: null,                    // /api/pan/samples (cached promise result)
  samplesP: null,
  fvsSpl: [],                       // FVS species codes for the datalist
};

const PanUI = (() => {
  const esc = s => String(s == null ? '' : s).replace(/[&<>"']/g,
    c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const clone = o => (o == null ? o : JSON.parse(JSON.stringify(o)));
  const FVS_DETAIL = 'Julia port of the USDA Forest Service Forest Vegetation Simulator (FVS, Southern variant; Dixon 2002, Essential FVS) — tree-list growth and yield.';
  const engineById = id => PAN.engines.find(e => e.id === id);
  const engineDetail = id => { const e = engineById(id); return id === 'fvs' || id === 'fvsjl' || !e ? FVS_DETAIL : (e.detail || ''); };
  function updateEngineInfo() {
    const t = el('engineText'), sel = el('engine'); if (!t) return;
    const d = engineDetail(PAN.engine);
    t.textContent = d;
    if (sel) sel.title = (engineLabel() + ': ' + d).slice(0, 600);
    el('engineHelp').hidden = !d;
  }
  const isPanId = id => !!id && id !== 'fvs' && id !== 'fvsjl';
  const isPan = () => isPanId(PAN.engine);
  const plugin = () => { const e = engineById(PAN.engine); return e ? e.plugin : null; };
  const curParams = () => { const pl = plugin(); return pl ? (PAN.paramsByPlugin[pl] || null) : null; };
  const engineLabel = (id = PAN.engine) => { const e = engineById(id); return e ? e.label : (id === 'fvs' ? 'FVSjl' : id); };
  const engineShort = id => ({ fvs: 'FVS', fvsjl: 'FVS', 'pan-biomass': 'BS', 'pan-pnet': 'PnET', 'pan-hbad': 'HBAD1', 'pan-hbad2': 'HBAD2', 'pan-hbad3p': 'HBAD3P' }[id] || id || 'FVS');

  async function getJSON(url) {
    const r = await fetch(url);
    const d = await r.json().catch(() => ({ error: 'bad response (' + r.status + ')' }));
    if (!r.ok || d.error) throw new Error(d.error || ('HTTP ' + r.status));
    return d;
  }

  /* ---------------- init / engine selector ---------------- */
  function init(meta) {
    PAN.meta = meta.pan || { available: false, engines: [], message: 'This backend has no Pan engines.' };
    const avail = !!PAN.meta.available;
    const list = [{ id: 'fvs', label: 'FVSjl', plugin: null, climate: false, enabled: true, detail: FVS_DETAIL }];
    (PAN.meta.engines || []).forEach(e => list.push({ ...e, enabled: avail && !e.coming }));
    PAN.engines = list;
    const sel = el('engine');
    sel.innerHTML = list.map(e => `<option value="${esc(e.id)}"${e.enabled ? '' : ' disabled'}
      title="${esc(e.coming ? 'Not available yet: not in Pan' : e.enabled ? e.label : (PAN.meta.message || 'Pan engines are unavailable on this backend'))}">${esc(e.label)}${e.enabled || e.coming ? '' : ' (unavailable)'}</option>`).join('');
    sel.value = PAN.engine;
    sel.onchange = () => setEngine(sel.value);
    el('addHarvest').onclick = addHarvest;
    el('panUpload').onclick = () => el('panFile').click();
    el('panFile').onchange = onFile;
    el('panDefault').onclick = () => loadDefault(plugin(), true);
    el('panSamples').onchange = onSample;
    if (avail) loadSamples().then(fillSampleSelects).catch(() => { });
    renderPanBar();
    refreshLabels();
    updateEngineInfo();
  }

  function setEngine(id, quiet) {
    const e = engineById(id); if (!e || !e.enabled) return false;
    PAN.engine = id;
    el('engine').value = id;
    PAN.pending = null;
    renderPanBar();
    if (isPan() && !curParams()) loadDefault(plugin(), false);
    if (typeof renderActions === 'function') renderActions();
    refreshSpecies();
    refreshLabels();
    updateSubline(); updateEngineInfo();
    if (!quiet) setStatus('engine: ' + engineLabel());
    return true;
  }

  function renderPanBar() {
    const pan = isPan(), e = engineById(PAN.engine);
    el('planpanel').classList.toggle('pan', pan);
    el('panParams').hidden = !pan;
    el('addHarvest').hidden = !pan;
    el('climateBtn').hidden = !(pan && e && e.climate);
    if (typeof ClimateUI !== 'undefined') ClimateUI.renderSummary();
    if (el('landscapeGrp')) el('landscapeGrp').hidden = pan || !currentAoiGeom;
    if (!pan) return;
    renderParamsInfo();
    fillSampleSelects();
  }

  // compact region summary: states and L3s as short lists; long lists (eco cells) as a count with the list in a tooltip
  function regionText(r) {
    if (r == null || r === '') return '';
    if (typeof r !== 'object') return String(r);
    const lab = { states: 'states', l3: 'EPA L3', eco_cells: 'eco cells' };
    return Object.entries(r).map(([k, v]) => {
      const name = lab[k] || k;
      if (Array.isArray(v)) return v.length > 8 ? `${v.length} ${name}` : `${name} ${v.join(', ')}`;
      return `${name}: ${typeof v === 'object' ? JSON.stringify(v) : v}`;
    }).join(' · ');
  }

  function renderParamsInfo() {
    const host = el('panParamsInfo'), p = curParams(), pl = plugin();
    if (!p) { host.innerHTML = `<span class="dim">Parameter set (${esc(pl)}): loading default…</span>`; return; }
    const warn = (p.warnings || []);
    host.innerHTML =
      `<b title="${esc(p.description || '')}">${esc(p.name || p.params_id)}</b>
       <span class="dim">· ${esc(p.plugin)}${p.isDefault ? ' · default set' : ' · uploaded'} · ${p.n_species != null ? esc(p.n_species) : '?'} species
       · ${p.n_ecos != null ? esc(p.n_ecos) : '?'} ecoregions${p.region ? ' · region ' + esc(regionText(p.region)) : ''}</span>
       ${p.description ? `<div class="dim">${esc(p.description)}</div>` : ''}
       ${(p.species || []).length ? `<details><summary>species (${p.species.length})</summary><div>${p.species.map(esc).join(', ')}</div></details>` : ''}
       ${warn.length ? `<details open><summary style="color:var(--warn)">${warn.length} warning${warn.length === 1 ? '' : 's'}</summary><ul>${warn.map(w => `<li>${esc(w)}</li>`).join('')}</ul></details>` : ''}`;
  }

  function panMsg(kind, text, btn) {
    const m = el('panMsg');
    m.className = 'panmsg ' + (kind || '');
    m.textContent = text || '';
    if (btn) {
      const b = document.createElement('button'); b.className = 'mini'; b.textContent = btn.text; b.onclick = btn.fn;
      m.appendChild(b);
    }
  }

  /* ---------------- parameter sets ---------------- */
  async function postYaml(text) {
    const r = await fetch('/api/pan/params', { method: 'POST', headers: { 'Content-Type': 'text/yaml' }, body: text });
    const d = await r.json().catch(() => ({ error: 'bad response (' + r.status + ')' }));
    if (!r.ok || d.error) throw new Error(d.error || ('HTTP ' + r.status));
    const info = { ...d, yaml: text, isDefault: false };
    PAN.store[d.params_id] = info;
    return info;
  }

  function accept(info) {
    PAN.paramsByPlugin[info.plugin] = info;
    PAN.pending = null;
    renderParamsInfo(); refreshSpecies();
    const w = (info.warnings || []).length;
    panMsg(w ? 'warn' : 'ok', `Loaded “${info.name || info.params_id}”${w ? ` with ${w} warning${w === 1 ? '' : 's'}` : ''}.`);
    setStatus(`parameter set “${info.name || info.params_id}” (${info.plugin}) loaded`);
  }

  async function submitYaml(text) {
    panMsg('', 'uploading parameter set…');
    let info;
    try { info = await postYaml(text); }
    catch (err) { panMsg('err', 'Parameter file rejected: ' + err.message); setStatus('parameter file rejected: ' + err.message); return; }
    if (info.plugin !== plugin()) {
      PAN.pending = info;
      const target = PAN.engines.find(e => e.plugin === info.plugin && e.enabled);
      panMsg('err', `This file is a “${info.plugin}” parameter set but ${engineLabel()} needs “${plugin()}”. It was not applied.`,
        target ? { text: 'Switch engine to ' + target.label, fn: () => { const pend = PAN.pending; if (setEngine(target.id, true) && pend) accept(pend); } } : null);
      return;
    }
    accept(info);
  }

  async function onFile(e) {
    const f = e.target.files[0]; e.target.value = ''; if (!f) return;
    if (f.size > 20e6) { panMsg('err', 'File too large for a parameter set (>20 MB).'); return; }
    submitYaml(await f.text());
  }

  async function loadDefault(pl, announce) {
    if (!pl) return;
    try {
      const d = await getJSON('/api/pan/default_params?plugin=' + encodeURIComponent(pl));
      const info = { ...d, yaml: null, isDefault: true };
      const have = PAN.paramsByPlugin[pl];
      if (!announce && have && !have.isDefault) return;     // an uploaded set arrived meanwhile: keep it
      PAN.paramsByPlugin[pl] = info;
      if (plugin() === pl) { renderParamsInfo(); refreshSpecies(); }
      if (announce) panMsg('ok', `Using the default ${pl} parameter set.`);
    } catch (err) {
      if (plugin() === pl) { PAN.paramsByPlugin[pl] = { isDefault: true, name: 'default (unavailable)', plugin: pl, warnings: [err.message] }; renderParamsInfo(); }
      panMsg('err', 'Could not load default parameters: ' + err.message);
    }
  }

  function loadSamples() {
    if (!PAN.samplesP) PAN.samplesP = getJSON('/api/pan/samples').then(d => (PAN.samples = d)).catch(err => { PAN.samplesP = null; throw err; });
    return PAN.samplesP;
  }
  function fillSampleSelects() {
    const d = PAN.samples; if (!d) return;
    const pl = plugin(), s = el('panSamples');
    if (s) {
      const items = (d.params || []).filter(p => p.plugin === pl);
      s.innerHTML = `<option value="">Samples… (${items.length})</option>` +
        items.map(p => `<option value="${esc(p.file)}" title="${esc(p.description || '')}">${esc(p.name || p.id)}</option>`).join('');
      s.disabled = !items.length;
    }
    if (typeof BundleUI !== 'undefined') BundleUI.fillSamples(d);
  }
  async function onSample() {
    const s = el('panSamples'), file = s.value; if (!file) return;
    s.value = '';
    try {
      const r = await fetch('/api/pan/sample?file=' + encodeURIComponent(file));
      if (!r.ok) throw new Error('HTTP ' + r.status);
      submitYaml(await r.text());
    } catch (err) { panMsg('err', 'Could not fetch sample: ' + err.message); }
  }

  // species suggestions
  function refreshSpecies() {
    const dl = el('splist'); if (!dl) return;
    const p = isPan() ? curParams() : null;
    const list = (p && p.species && p.species.length) ? p.species : PAN.fvsSpl;
    dl.innerHTML = list.map(s => `<option value="${esc(s)}">`).join('');
  }

  /* ---------------- scenario <-> UI state ---------------- */
  function snapshot() {
    const pl = plugin(), p = curParams();
    let params = null;
    if (pl) params = (p && !p.isDefault) ? { id: p.params_id, name: p.name, plugin: p.plugin }
                                         : { id: null, name: (p && p.name) || 'default', plugin: pl, isDefault: true };
    const e = engineById(PAN.engine);
    return { engine: PAN.engine, params, climate: (e && e.climate && PAN.climate) ? clone(PAN.climate) : null };
  }
  // header subline names the ACTIVE engine (it used to stay on the first entry of /api/meta engines)
  function updateSubline() {
    const sub = el('subline');
    if (sub) sub.textContent = sub.textContent.replace(/engine: .*$/, 'engine: ' + engineLabel());
  }

  function restore(s) {
    const eng = s.engine && engineById(s.engine) ? s.engine : 'fvs';
    PAN.engine = eng; el('engine').value = eng; updateSubline(); updateEngineInfo();
    PAN.climate = s.climate ? clone(s.climate) : null;
    const pr = s.params;
    if (pr && pr.plugin) {
      if (pr.id && PAN.store[pr.id]) PAN.paramsByPlugin[pr.plugin] = PAN.store[pr.id];
      else if (pr.isDefault || !pr.id) delete PAN.paramsByPlugin[pr.plugin];
    }
    renderPanBar();
    if (isPan() && !curParams()) loadDefault(plugin(), false);
    refreshSpecies(); refreshLabels();
  }

  function chipTag(s) {
    const id = s.engine || 'fvs';
    return isPanId(id) ? ` <span class="etag" title="${esc(scenMeta(s))}">${esc(engineShort(id))}</span>` : '';
  }
  function scenMeta(s) {
    const id = s.engine || 'fvs';
    if (!isPanId(id)) return 'FVSjl';
    const pn = (s.params && s.params.name) || (s.result && s.result.engine_info && s.result.engine_info.params_name) || 'default params';
    let t = `${engineLabel(id)} · params: ${pn}`;
    if (s.climate && typeof ClimateUI !== 'undefined') t += ' · climate: ' + ClimateUI.summary(s.climate);
    return t;
  }

  function approx(engineId) { return ''; }   // Pan volumes now come from FIA species equations (no longer flagged)
  // relabel volume metrics when the shown/selected engine is a Pan one
  function refreshLabels() {
    const id = (typeof lastSim !== 'undefined' && lastSim && lastSim.engine) || PAN.engine;
    const sfx = approx(id);
    [['simMetric', 'volume', 'Standing volume (ft³)'], ['simMetric', 'removed_volume', 'Removed volume (ft³)'],
     ['simMapMetric', 'volume', 'Standing volume (ft³/ac)']].forEach(([sid, v, base]) => {
      const o = el(sid) && [...el(sid).options].find(x => x.value === v);
      if (o) o.textContent = base + sfx;
    });
  }

  /* ---------------- request building / validation ---------------- */
  function simBody(geom, p) {
    const pp = clone(p);
    const body = { geometry: geom, plan: pp, engine: PAN.engine };
    if (isPan()) {
      pp.actions.forEach(a => { if (a.kind === 'thin') a.mode = 'rule'; });   // PPE is FVS-only
      const pr = curParams();
      if (pr && !pr.isDefault && pr.params_id) body.params_id = pr.params_id;
      const e = engineById(PAN.engine);
      if (e && e.climate) body.climate = PAN.climate ? clone(PAN.climate) : null;
    }
    return body;
  }

  async function recoverFromError(msg) {
    const p = curParams();
    if (!isPan() || !p || p.isDefault || !p.yaml || !/params/i.test(msg || '')) return false;
    try { const info = await postYaml(p.yaml); PAN.paramsByPlugin[info.plugin] = info; return true; }
    catch (_) { return false; }
  }

  const num = v => typeof v === 'number' && isFinite(v);
  function harvestProblems(a, p) {
    const out = [], add = (path, msg) => out.push({ path, msg });
    const s = a.selection || {}, r = a.removal || {}, st = a.stands || {};
    if (!(Number.isInteger(a.cycle) && a.cycle >= 0 && a.cycle <= (p.cycles || 1) - 1)) add('cycle', `cycle must be an integer 0–${(p.cycles || 1) - 1}`);
    if (!(Number.isInteger(a.repeat_every) && a.repeat_every >= 0)) add('repeat_every', 'repeat every must be a whole number ≥ 0 (0 = once)');
    if (!num(s.age_lo) || s.age_lo < 0) add('selection.age_lo', 'min age must be ≥ 0');
    if (!num(s.age_hi) || s.age_hi < (s.age_lo || 0)) add('selection.age_hi', 'max age must be ≥ min age');
    if (!num(s.dbh_lo) || s.dbh_lo < 0) add('selection.dbh_lo', 'min DBH must be ≥ 0');
    if (!num(s.dbh_hi) || s.dbh_hi < (s.dbh_lo || 0)) add('selection.dbh_hi', 'max DBH must be ≥ min DBH');
    if (!String(s.species || '').trim()) add('selection.species', 'species: “all” or comma-separated codes');
    if (r.type !== 'clearcut') {
      if (!num(r.value) || r.value < 0) add('removal.value', 'removal value must be ≥ 0');
      else if (r.type === 'percent_biomass' && r.value > 100) add('removal.value', 'percent must be 0–100');
    }
    if (!num(st.area_pct) || st.area_pct <= 0 || st.area_pct > 100) add('stands.area_pct', 'stand area must be in (0, 100]');
    if (a.plant_after) {
      const pa = a.plant_after;
      if (!String(pa.species || '').trim()) add('plant_after.species', 'plant-after species required');
      if (!num(pa.tpa) || pa.tpa <= 0) add('plant_after.tpa', 'plant-after TPA must be > 0');
      if (!num(pa.survival) || pa.survival < 0 || pa.survival > 100) add('plant_after.survival', 'survival must be 0–100');
    }
    return out;
  }
  function preflight(p) {
    const hv = p.actions.filter(a => a.kind === 'harvest');
    if (!isPan()) return hv.length ? 'FVSjl cannot run harvest actions — remove them or switch to a Pan engine' : '';
    const e = engineById(PAN.engine);
    if (!e || !e.enabled) return 'Pan engines are unavailable on this backend';
    if (!(p.cycles >= 1 && p.period >= 1)) return 'cycles and years per cycle must be ≥ 1';
    for (let i = 0; i < p.actions.length; i++) {
      const a = p.actions[i];
      if (a.kind === 'harvest') { const pr = harvestProblems(a, p); if (pr.length) return `action ${i + 1} (harvest): ${pr[0].msg}`; }
      if (a.kind === 'plant' && (!(a.tpa > 0) || !(a.survival >= 0 && a.survival <= 100))) return `action ${i + 1} (plant): TPA must be > 0 and survival 0–100`;
      if (a.kind === 'thin' && (!isFinite(a.target) || a.target < 0 || a.dbh_hi < a.dbh_lo)) return `action ${i + 1} (thin): check target / DBH range`;
    }
    return '';
  }

  /* ---------------- harvest action ---------------- */
  const REMOVAL = [['percent_biomass', 'remove % of biomass'], ['residual_ba', 'leave BA (ft²/ac)'], ['residual_tpa', 'leave TPA (trees/ac)'],
                   ['residual_sdi', 'leave SDI'], ['clearcut', 'clearcut']];
  const DIRS = [['below', 'from below'], ['above', 'from above'], ['proportional', 'proportional']];
  const RANKS = [['all', 'all stands'], ['max_biomass', 'highest biomass first'], ['oldest', 'oldest first'], ['random', 'random stands']];
  const UNIT = { percent_biomass: '%', residual_ba: 'ft²/ac', residual_tpa: 'tpa', residual_sdi: 'SDI', clearcut: '' };

  function newHarvest() {
    return { kind: 'harvest', cycle: 1, label: 'Harvest', repeat_every: 0,
      selection: { species: 'all', age_lo: 0, age_hi: 999, dbh_lo: 0, dbh_hi: 999 },
      removal: { type: 'percent_biomass', value: 30, direction: 'below' },
      stands: { rank: 'all', area_pct: 100 }, plant_after: null };
  }
  function addHarvest() { plan.actions.push(newHarvest()); renderActions(); }

  function harvestSummary(a, p) {
    const s = a.selection || {}, r = a.removal || {}, st = a.stands || {};
    const rem = r.type === 'clearcut' ? 'clearcut'
      : r.type === 'percent_biomass' ? `remove ${r.value}% biomass ${r.direction || ''}`
      : `${r.type.replace('residual_', 'leave ').toUpperCase().replace('LEAVE ', 'leave ')}→${r.value} ${r.direction || ''}`;
    const sel = [s.species && s.species !== 'all' ? s.species : '',
      (s.age_lo > 0 || s.age_hi < 999) ? `age ${s.age_lo}–${s.age_hi}` : '',
      (s.dbh_lo > 0 || s.dbh_hi < 999) ? `DBH ${s.dbh_lo}–${s.dbh_hi}″` : ''].filter(Boolean).join(' ');
    const when = `@+${(a.cycle || 0) * (p.period || 10)}yr` + (a.repeat_every > 0 ? ` every ${a.repeat_every * (p.period || 10)}yr` : '');
    const stands = (st.rank && st.rank !== 'all' ? ` on ${(RANKS.find(x => x[0] === st.rank) || [0, st.rank])[1]}` : '') +
      (st.area_pct < 100 ? ` (${st.area_pct}% area)` : '');
    const pa = a.plant_after ? `, then plant ${a.plant_after.species} ${a.plant_after.tpa}tpa` : '';
    return `harvest ${rem}${sel ? ' [' + sel + ']' : ''}${stands} ${when}${pa}`;
  }

  const getPath = (o, path) => path.split('.').reduce((x, k) => (x == null ? x : x[k]), o);
  function setPath(o, path, v) { const ks = path.split('.'); const last = ks.pop(); const t = ks.reduce((x, k) => x[k], o); t[last] = v; }

  function harvestRow(a, i, p) {
    const s = a.selection, r = a.removal, st = a.stands, pa = a.plant_after;
    const F = (path, type, extra = '', cls = '') => {
      const v = getPath(a, path);
      const val = (type === 'n' && !num(v)) ? '' : esc(v);
      return `data-i="${i}" data-h="${path}" data-t="${type}" value="${val}" ${extra}`;
    };
    const sp = (path, label) => `<input list="splist" class="wide" ${F(path, 's')} aria-label="${label}" title="${label}: “all” or comma-separated codes">`;
    const n = (path, label, extra = '') => `<input type="number" class="num" ${F(path, 'n', extra)} aria-label="${label}" title="${label}">`;
    const sel = (path, items, label, dis) => `<select data-i="${i}" data-h="${path}" data-t="s" aria-label="${label}" title="${label}"${dis ? ' disabled' : ''}>${_opt(items, getPath(a, path))}</select>`;
    const cc = r.type === 'clearcut';
    return `<div class="hv">
      <div class="hl"><span class="tag harvest">Harvest</span>
        <input class="wide" ${F('label', 's')} aria-label="harvest label" title="label (optional)" placeholder="label">
        <label>at</label>${n('cycle', 'cycle', `min="0" max="${(p.cycles || 1) - 1}" step="1"`)}<span style="color:var(--muted);font-size:11px">${a.cycle === 0 ? 'now' : '+' + (a.cycle * p.period) + 'yr'}</span>
        <label>repeat every</label>${n('repeat_every', 'repeat every N cycles (0 = once)', 'min="0" step="1"')}<span style="color:var(--muted);font-size:11px">cyc (0 = once)</span>
        <button class="rm" data-rm="${i}" title="remove action" aria-label="remove harvest action">×</button></div>
      <div class="hl grpbox"><label>select</label>sp ${sp('selection.species', 'species')}
        <label>age</label>${n('selection.age_lo', 'min age', 'min="0"')}–${n('selection.age_hi', 'max age', 'min="0"')}
        <label>DBH</label>${n('selection.dbh_lo', 'min DBH (in)', 'min="0" step="any"')}–${n('selection.dbh_hi', 'max DBH (in)', 'min="0" step="any"')}″</div>
      <div class="hl grpbox"><label>remove</label>${sel('removal.type', REMOVAL, 'removal type')}
        ${cc ? '' : n('removal.value', 'removal value', 'min="0" step="any"')}<span style="color:var(--muted);font-size:11px">${UNIT[r.type] || ''}</span>
        ${sel('removal.direction', DIRS, 'removal direction', cc)}</div>
      <div class="hl grpbox"><label>stands</label>${sel('stands.rank', RANKS, 'stand ranking')}
        <label>area</label>${n('stands.area_pct', 'percent of stands treated', 'min="1" max="100" step="any"')}%
        <label class="chk" style="margin-left:8px"><input type="checkbox" data-i="${i}" data-h="plant_after.on" data-t="b"${pa ? ' checked' : ''}> plant after</label>
        ${pa ? `${sp('plant_after.species', 'plant-after species')} ${n('plant_after.tpa', 'plant-after trees/ac', 'min="1"')}tpa
          <label>surv</label>${n('plant_after.survival', 'plant-after survival %', 'min="0" max="100"')}%` : ''}</div>
      <div class="hsum" data-hsum="${i}"></div></div>`;
  }

  function refreshHarvestRow(host, i) {
    const a = plan.actions[i], row = host.querySelector(`[data-hsum="${i}"]`); if (!row) return;
    const probs = harvestProblems(a, plan);
    row.textContent = harvestSummary(a, plan) + (probs.length ? '  ⚠ ' + probs.map(x => x.msg).join('; ') : '');
    row.classList.toggle('bad', probs.length > 0);
    host.querySelectorAll(`[data-i="${i}"][data-h]`).forEach(inp => inp.classList.toggle('bad', probs.some(x => x.path === inp.dataset.h)));
    const act = row.closest('.act'); if (act) act.classList.toggle('invalid', probs.length > 0);
  }

  // called by renderActions after the rows are in the DOM
  function bindActions(host) {
    host.querySelectorAll('[data-h]').forEach(inp => {
      const handler = e => {
        const i = +inp.dataset.i, path = inp.dataset.h, t = inp.dataset.t, a = plan.actions[i];
        if (path === 'plant_after.on') {
          a.plant_after = inp.checked ? { species: (curParams() && curParams().species && curParams().species[0]) || 'PITA', tpa: 600, survival: 80 } : null;
          renderActions(); return;
        }
        let v = t === 'n' ? (inp.value === '' ? NaN : +inp.value) : inp.value;
        if (path === 'selection.species' || path === 'plant_after.species') v = String(v).replace(/\s+/g, '');
        setPath(a, path, v);
        if (e.type === 'change' && ['cycle', 'removal.type'].includes(path)) { renderActions(); return; }
        refreshHarvestRow(host, i);
      };
      inp.addEventListener('change', handler);
      if (inp.dataset.t !== 'b' && inp.tagName === 'INPUT') inp.addEventListener('input', handler);
    });
    plan.actions.forEach((a, i) => { if (a.kind === 'harvest') refreshHarvestRow(host, i); });
  }

  /* ---------------- warnings / status ---------------- */
  function showWarnings(d) {
    const box = el('simWarn'), w = (d && d.warnings) || [];
    if (!w.length) { box.hidden = true; box.innerHTML = ''; return; }
    box.hidden = false;
    const det = (d && d.warning_details) || [];
    box.innerHTML = `<b>${w.length} warning${w.length === 1 ? '' : 's'} from the engine</b><ul>${w.map(x => `<li>${esc(x)}</li>`).join('')}</ul>` +
      (det.length ? `<details><summary>per-stand details (${det.length})</summary><ul>${det.map(x => `<li>${esc(x)}</li>`).join('')}</ul></details>` : '');
  }
  function statusWithWarnings(msg, d) {
    const w = (d && d.warnings) || [];
    setStatus(msg + (w.length ? ` — ⚠ ${w.length} warning${w.length === 1 ? '' : 's'}: ${w[0]}${w.length > 1 ? ' …' : ''}` : ''));
    showWarnings(d);
  }

  return { init, setEngine, engineDetail, isPan, plugin, engineLabel, esc, clone, snapshot, restore, chipTag, scenMeta, approx, refreshLabels,
           simBody, recoverFromError, preflight, harvestSummary, harvestRow, bindActions, showWarnings, statusWithWarnings,
           postYaml, loadSamples, fillSampleSelects, refreshSpecies, isPanId, engineById, curParams, accept, newHarvest, harvestProblems,
           renderPanBar, setFvsSpecies: l => { PAN.fvsSpl = l; refreshSpecies(); } };
})();
