'use strict';
/* Export / import of a whole session: AOI patches + scenarios (plan, engine, parameter YAML, climate). */

const BundleUI = (() => {
  const msg = (t, kind) => { const m = el('bundleMsg'); if (m) { m.textContent = t; m.style.color = kind === 'err' ? '#f87171' : kind === 'warn' ? 'var(--warn)' : ''; } setStatus(t); };

  function download(name, text, type) {
    const blob = new Blob([text], { type }), url = URL.createObjectURL(blob), a = document.createElement('a');
    a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 0);
  }

  function buildBundle() {
    const stamp = new Date().toISOString();
    const list = scenarios.length ? scenarios : [{ name: el('scenName').value.trim() || 'Draft', plan, ...PanUI.snapshot() }];
    const out = list.map((s, i) => {
      // when the builder shows this scenario, the live plan edits count
      const pl = (i === activeIdx) ? { ...plan, cycles: +el('pCycles').value || plan.cycles, period: +el('pPeriod').value || plan.period } : s.plan;
      let params = null;
      if (s.params && !s.params.isDefault) {
        const st = PAN.store[s.params.id];
        params = { name: s.params.name, plugin: s.params.plugin, yaml: st ? st.yaml : (s.params.yaml || null) };
      }
      return { name: s.name, engine: s.engine || 'fvs', params, climate: s.climate || null, plan: PanUI.clone(pl) };
    });
    return {
      format: 'forest-explorer-bundle', version: 1, created: stamp,
      aoi: { patches: aoiParts.map((g, i) => ({ type: 'Feature', properties: { name: `Patch ${i + 1}` }, geometry: g })) },
      scenarios: out, active: activeIdx >= 0 ? activeIdx : 0,
    };
  }

  function doExport() {
    if (!aoiParts.length && !scenarios.length) { msg('Nothing to export yet — add an AOI patch and/or a scenario.', 'warn'); return; }
    const b = buildBundle();
    const missing = b.scenarios.filter(s => s.params && !s.params.yaml).length;
    download(`forest-explorer-bundle-${b.created.slice(0, 19).replace(/[:T]/g, '-')}.json`, JSON.stringify(b, null, 2), 'application/json');
    msg(`Exported ${b.aoi.patches.length} patch(es) and ${b.scenarios.length} scenario(s)` +
        (scenarios.some(s => s.result) ? ' (results are not included — Run after importing)' : '') +
        (missing ? `; ${missing} parameter set(s) had no YAML text to embed` : '') + '.', missing ? 'warn' : '');
  }

  function validate(b) {
    if (!b || b.format !== 'forest-explorer-bundle') throw new Error('not a forest-explorer-bundle file');
    if (b.version !== 1) throw new Error('unsupported bundle version ' + b.version);
    if (!Array.isArray(b.scenarios)) throw new Error('bundle has no scenarios array');
  }

  async function importBundle(b) {
    validate(b);
    if ((aoiParts.length || scenarios.length) && !window.confirm('Importing replaces the current AOI patches and scenarios. Continue?')) { msg('import cancelled'); return; }
    const notes = [];
    // 1. AOI patches (same path as an upload: every polygon becomes a patch)
    const patches = ((b.aoi && b.aoi.patches) || []).map(f => (f && f.type === 'Feature') ? f.geometry : f).filter(Boolean);
    if (typeof cancelDraw === 'function') cancelDraw();
    aoiParts = [];
    patches.forEach(g => polysOf(g).forEach(c => aoiParts.push({ type: 'Polygon', coordinates: c })));
    scenarios = []; activeIdx = -1; lastSim = null;
    if (aoiParts.length) await submitAoiUnion(); else clearAoi();
    // 2. scenarios (re-upload embedded YAML -> fresh params_ids; dedupe identical texts)
    const byYaml = new Map(), made = [];
    for (const s of b.scenarios.slice(0, MAX_SCENARIOS)) {
      const eng = (s.engine && PanUI.engineById(s.engine)) ? s.engine : 'fvs';
      if (s.engine && eng !== s.engine) notes.push(`scenario “${s.name}”: engine ${s.engine} is unknown here — using FVSjl`);
      const e = PanUI.engineById(eng);
      if (e && !e.enabled) notes.push(`scenario “${s.name}”: ${e.label} is unavailable on this backend`);
      let params = null;
      if (s.params && s.params.yaml && eng !== 'fvs') {
        try {
          if (!byYaml.has(s.params.yaml)) byYaml.set(s.params.yaml, await PanUI.postYaml(s.params.yaml));
          const info = byYaml.get(s.params.yaml);
          params = { id: info.params_id, name: info.name, plugin: info.plugin };
          if (e && e.plugin && info.plugin !== e.plugin) notes.push(`scenario “${s.name}”: parameter plugin “${info.plugin}” does not match ${e.label}`);
          (info.warnings || []).forEach(w => notes.push(`params “${info.name}”: ${w}`));
        } catch (err) { notes.push(`scenario “${s.name}”: parameter upload failed (${err.message}) — default set will be used`); }
      } else if (s.params && !s.params.yaml && eng !== 'fvs') notes.push(`scenario “${s.name}”: bundle has no YAML for “${s.params.name}” — default set will be used`);
      const p = s.plan || {};
      made.push({
        name: String(s.name || `Scenario ${made.length + 1}`),
        plan: { name: p.name || 'plan', cycles: +p.cycles || 10, period: +p.period || 10, actions: Array.isArray(p.actions) ? p.actions : [] },
        result: null, engine: eng,
        params: params || (e && e.plugin ? { id: null, name: 'default', plugin: e.plugin, isDefault: true } : null),
        climate: s.climate || null,
      });
    }
    if (b.scenarios.length > MAX_SCENARIOS) notes.push(`only the first ${MAX_SCENARIOS} scenarios were imported`);
    scenarios = made; scenSeq = Math.max(scenSeq, made.length);
    if (made.length) activateScenario(Math.max(0, Math.min(made.length - 1, b.active | 0)));
    else renderScenarioChips();
    msg(`Imported ${aoiParts.length} patch(es) and ${made.length} scenario(s) — not simulated, press Run to compute.` + (notes.length ? ' Notes: ' + notes.join(' | ') : ''), notes.length ? 'warn' : '');
  }

  async function onFile(e) {
    const f = e.target.files[0]; e.target.value = ''; if (!f) return;
    try { await importBundle(JSON.parse(await f.text())); }
    catch (err) { msg('import failed: ' + err.message, 'err'); }
  }

  async function onSample() {
    const s = el('bundleSample'), file = s.value; if (!file) return;
    s.value = '';
    try {
      const r = await fetch('/api/pan/sample?file=' + encodeURIComponent(file));
      if (!r.ok) throw new Error('HTTP ' + r.status);
      await importBundle(JSON.parse(await r.text()));
    } catch (err) { msg('sample bundle failed: ' + err.message, 'err'); }
  }

  function fillSamples(d) {
    const s = el('bundleSample'); if (!s) return;
    const items = d.bundles || [];
    s.innerHTML = `<option value="">${items.length ? 'Load sample bundle…' : '(no sample bundles)'}</option>` +
      items.map(b => `<option value="${PanUI.esc(b.file)}" title="${PanUI.esc(b.description || '')}">${PanUI.esc(b.name)}</option>`).join('');
    s.disabled = !items.length;
  }

  function init() {
    el('bundleExport').onclick = doExport;
    el('bundleImport').onclick = () => el('bundleFile').click();
    el('bundleFile').onchange = onFile;
    el('bundleSample').onchange = onSample;
    el('bundleSample').disabled = true;
  }

  return { init, fillSamples, buildBundle, importBundle };
})();
