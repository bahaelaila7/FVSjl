'use strict';
/* DEV ONLY — loaded when the URL contains ?mock=1.
   Intercepts window.fetch for the new /api/pan/* endpoints, the Pan engines of /api/simulate and the
   `pan` block of /api/meta, returning canned responses so the UI can be exercised without the Julia backend.
   Everything else (/api/aoi, tiles, FVS simulate ...) passes through to the real server; if the real
   server is unreachable /api/meta, /api/attrs and /api/aoi fall back to minimal canned replies. */
(function () {
  const realFetch = window.fetch.bind(window);
  const json = (obj, status = 200) => new Response(JSON.stringify(obj), { status, headers: { 'Content-Type': 'application/json' } });
  const text = (t, status = 200, type = 'text/plain') => new Response(t, { status, headers: { 'Content-Type': type } });
  const delay = ms => new Promise(r => setTimeout(r, ms));

  const ENGINES = [
    { id: 'pan-biomass', label: 'Pan · Biomass Succession', plugin: 'biomass', climate: false },
    { id: 'pan-pnet', label: 'Pan · PnET-Succession', plugin: 'pnet', climate: true },
    { id: 'pan-hbad', label: 'Pan · HBAD', plugin: 'hbad', climate: false },
  ];
  const yamlOf = (plugin, name, extra = '') => `pan_params:\n  plugin: ${plugin}\n  name: ${name}\n  description: mock ${plugin} set\nspecies: [PITA, QUAL, LIST, PIEL, QUNI]\n${extra}`;
  const SAMPLE_FILES = {
    'mock_biomass.yml': yamlOf('biomass', 'Mock Biomass v1'),
    'mock_pnet.yml': yamlOf('pnet', 'Mock PnET v1'),
    'mock_pnet_warn.yml': yamlOf('pnet', 'Mock PnET (with warnings)', 'warn: true\n'),
    'mock_hbad.yml': yamlOf('hbad', 'Mock HBAD v1'),
  };
  const BUNDLE = {
    format: 'forest-explorer-bundle', version: 1, created: new Date().toISOString(),
    aoi: { patches: [{ type: 'Feature', properties: {}, geometry: { type: 'Polygon', coordinates: [[[-83.6, 31.0], [-83.4, 31.0], [-83.4, 31.2], [-83.6, 31.2], [-83.6, 31.0]]] } }] },
    scenarios: [
      { name: 'Grow only (PnET)', engine: 'pan-pnet', params: { name: 'Mock PnET v1', plugin: 'pnet', yaml: SAMPLE_FILES['mock_pnet.yml'] },
        climate: { baseline: { mode: 'forecast', T: 30, K: 3 }, co2: { mode: 'linear', rate_ppm_per_yr: 2.0 }, changes: [{ var: 'tmean', kind: 'add', value: 2, shape: 'ramp', start_year: 2030, end_year: 2080, months: null }] },
        plan: { name: 'plan', cycles: 8, period: 10, actions: [] } },
      { name: 'Harvest 40% (BS)', engine: 'pan-biomass', params: null, climate: null,
        plan: { name: 'plan', cycles: 8, period: 10, actions: [{ kind: 'harvest', cycle: 2, label: 'thin from below', repeat_every: 3,
          selection: { species: 'all', age_lo: 0, age_hi: 999, dbh_lo: 0, dbh_hi: 999 }, removal: { type: 'percent_biomass', value: 40, direction: 'below' },
          stands: { rank: 'all', area_pct: 100 }, plant_after: { species: 'PITA', tpa: 600, survival: 80 } }] } },
    ],
    active: 0,
  };
  const SAMPLES = {
    params: [
      { id: 'mock-bs', file: 'mock_biomass.yml', plugin: 'biomass', name: 'Mock Biomass v1', description: 'canned' },
      { id: 'mock-pnet', file: 'mock_pnet.yml', plugin: 'pnet', name: 'Mock PnET v1', description: 'canned' },
      { id: 'mock-pnet-w', file: 'mock_pnet_warn.yml', plugin: 'pnet', name: 'Mock PnET (warnings)', description: 'canned, returns warnings' },
      { id: 'mock-hbad', file: 'mock_hbad.yml', plugin: 'hbad', name: 'Mock HBAD v1', description: 'canned' },
    ],
    aois: [], bundles: [{ name: 'Mock two-scenario bundle', file: 'mock_bundle.json', description: 'PnET grow-only + biomass harvest' }], climates: [],
  };

  let pid = 0;
  const knownParams = {};
  function parseParams(t) {
    const m = /pan_params:\s*[\r\n]+\s+plugin:\s*["']?(\w+)/.exec(t);
    if (!m) return { error: 'missing pan_params header (mock parser)' };
    if (!['biomass', 'pnet', 'hbad'].includes(m[1])) return { error: 'unknown plugin ' + m[1] };
    const nm = /name:\s*(.+)/.exec(t);
    const info = {
      params_id: 'mock-' + (++pid), plugin: m[1], name: nm ? nm[1].trim() : 'uploaded', description: 'mock parse',
      n_species: 5, n_ecos: 3, region: { states: ['GA', 'FL'], l3: ['65', '75'] }, species: ['PITA', 'QUAL', 'LIST', 'PIEL', 'QUNI'],
      warnings: /warn:\s*true/.test(t) ? ['mock warning: 2 species have no shade tolerance, defaulted to 3', 'mock warning: ecoregion 8.5.1 has no soil data'] : [],
    };
    knownParams[info.params_id] = info;
    return info;
  }
  const defaultParams = pl => ({ params_id: 'mock-default-' + pl, plugin: pl, name: 'Default ' + pl, description: 'shipped default (mock)', n_species: 12, n_ecos: 9, region: { name: 'Region 4' }, species: ['PITA', 'QUAL', 'LIST', 'PIEL', 'QUNI', 'ACRU', 'LITU', 'PIVI', 'CAGL', 'CATO', 'FAGR', 'PIST'], warnings: [] });

  /* ---- climate preview ---- */
  function climatePreview(body) {
    const y0 = body.start_year || 2025, n = body.years || 100, c = body.climate || {};
    const defs = [['tmin', 'Minimum temperature', '°C', 9, 7], ['tmax', 'Maximum temperature', '°C', 24, 8], ['tmean', 'Mean temperature', '°C', 16.5, 7.5], ['ppt', 'Precipitation', 'mm/month', 105, 35], ['co2', 'CO₂', 'ppm', 420, 0]];
    const co = c.co2 || { mode: 'constant' };
    const vars = defs.map(([name, label, units, mean, amp]) => {
      const bm = [], mm = [], ba = [], ma = [];
      for (let k = 0; k < n; k++) {
        let sb = 0, sm = 0;
        for (let m = 0; m < 12; m++) {
          let base, mod;
          if (name === 'co2') {
            base = 420;
            mod = base + (co.mode === 'linear' ? (co.rate_ppm_per_yr || 0) * k : co.mode === 'ssp245' ? 1.8 * k : co.mode === 'ssp585' ? 3.2 * k : 0);
          } else {
            const seas = name === 'ppt' ? Math.max(5, mean + amp * 0.6 * Math.sin(2 * Math.PI * (m - 3) / 12) + 10 * Math.sin(k * 0.7 + m)) : mean - amp * Math.cos(2 * Math.PI * (m + 0.5) / 12) + 0.3 * Math.sin(k * 0.9 + m);
            base = seas + (name !== 'ppt' ? 0.02 * k : 0); mod = base;
            (c.changes || []).forEach(ch => {
              if (!(ch.var === name || (ch.var === 'tmean' && (name === 'tmin' || name === 'tmax')))) return;
              if (ch.months && ch.months.length && !ch.months.includes(m + 1)) return;
              const yr = y0 + k; let f;
              if (ch.shape === 'ramp') f = yr < ch.start_year ? 0 : yr >= ch.end_year ? 1 : (yr - ch.start_year) / Math.max(1, ch.end_year - ch.start_year);
              else f = (yr >= ch.start_year && yr <= ch.end_year) ? 1 : 0;
              if (ch.kind === 'add') mod += ch.value * f; else mod *= 1 + (ch.value / 100) * f;
            });
          }
          bm.push(base); mm.push(mod); sb += base; sm += mod;
        }
        ba.push(sb / 12); ma.push(sm / 12);
      }
      return { name, label, units, years: Array.from({ length: n }, (_, i) => y0 + i), baseline_annual: ba, modified_annual: ma, baseline_monthly: bm, modified_monthly: mm };
    });
    return { l4s: ['8.3.5a', '8.5.1b', '8.3.4c'], start_year: y0, vars };
  }

  /* ---- simulate ---- */
  function simulate(body) {
    const p = body.plan || {}, cyc = p.cycles || 10, per = p.period || 10, acts = p.actions || [];
    const warn = ['MOCK: canned response, not a real simulation'];
    if (body.engine === 'pan-pnet') warn.push(body.climate ? 'mock: climate scenario received (' + (body.climate.changes || []).length + ' change(s))' : 'mock: default climate');
    const cycles = []; const dist = []; let ba = 60, rem = 0;
    for (let i = 0; i <= cyc; i++) {
      let removed = 0;
      acts.forEach(a => {
        const rep = a.repeat_every > 0 && i > a.cycle && (i - a.cycle) % a.repeat_every === 0;
        if ((a.kind === 'harvest' || a.kind === 'thin') && i > 0 && (a.cycle === i || rep)) { const f = a.kind === 'harvest' ? (a.removal.type === 'percent_biomass' ? a.removal.value / 100 : 0.4) : 0.35; removed = ba * f * 30; ba *= (1 - f); }
      });
      ba += (130 - ba) * 0.12 * (body.engine === 'pan-pnet' && body.climate ? 0.9 : 1);
      cycles.push({ elapsed: i * per, ba_ac: ba, carbon_total: ba * 11 * 40, carbon_live_ag: ba * 8 * 40, volume: ba * 30 * 40, removed_volume: removed * 40 });
      dist.push([['PITA', 'Loblolly pine', 0.55], ['QUAL', 'Laurel oak', 0.3], ['LIST', 'Sweetgum', 0.15]].flatMap(([s, cm, w]) =>
        [0, 2, 4, 6, 8, 10, 14, 20].map((d, k) => ({ symbol: s, common: cm, dbh_lo: d, ba: ba * w * Math.exp(-Math.pow((d - 8 - i), 2) / 40) / 3 }))));
    }
    return { engine: body.engine, warnings: warn, engine_info: { params_name: body.params_id ? (knownParams[body.params_id] || {}).name || body.params_id : 'default', plugin: (ENGINES.find(e => e.id === body.engine) || {}).plugin },
      cycles, dist, ncycles: cycles.length, period: per, nplots_sim: 12, nplots_aoi: 12, capped: false, image: null, domains: {} };
  }

  const bodyJSON = init => { try { return JSON.parse(init && init.body); } catch (_) { return {}; } };

  window.fetch = async function (input, init) {
    const url = typeof input === 'string' ? input : input.url;
    const path = url.split('?')[0], q = new URLSearchParams(url.split('?')[1] || '');
    const method = ((init && init.method) || 'GET').toUpperCase();
    if (path === '/api/meta') {
      let m;
      try { const r = await realFetch(input, init); m = await r.json(); }
      catch (_) { m = { default_attr: 'carbon_l', bounds: { west: -127.9, south: 22.7, east: -65.3, north: 51.6 }, nplots_total: 0, sources: [{ id: 'mock', name: 'mock source' }], engines: [{ id: 'fvsjl', name: 'FVSjl' }] }; }
      if (!m.pan) m.pan = { available: true, engines: ENGINES, message: '' };
      if (q.get('unavailable') || /[?&]panoff=1/.test(location.search)) m.pan = { available: false, engines: ENGINES, message: 'mock: Pan is not installed on this backend' };
      return json(m);
    }
    if (path === '/api/attrs') {
      try { return await realFetch(input, init); } catch (_) { return json([{ name: 'carbon_l', label: 'Live carbon', units: 't/ha', lo: 0, hi: 100, cmap: 'greens' }]); }
    }
    if (path === '/api/aoi' && method === 'POST') {
      try { return await realFetch(input, init); }
      catch (_) { return json({ nplots: 12, acres: 480, npixels: 2100, carbon_l: 12000, drybio_l: 30000, volcfnet_l: 1.2e6, ba_total: 40000,
        bbox: { west: -83.7, south: 30.9, east: -83.3, north: 31.3 }, cells: [{ symbol: 'PITA', common: 'Loblolly pine', dbh_lo: 6, ba: 40 }, { symbol: 'QUAL', common: 'Laurel oak', dbh_lo: 10, ba: 20 }] }); }
    }
    if (path === '/api/pan/samples') { await delay(80); return json(SAMPLES); }
    if (path === '/api/pan/sample') {
      const f = q.get('file');
      if (f === 'mock_bundle.json') return text(JSON.stringify(BUNDLE), 200, 'application/json');
      return SAMPLE_FILES[f] ? text(SAMPLE_FILES[f], 200, 'text/yaml') : json({ error: 'no such sample' }, 404);
    }
    if (path === '/api/pan/default_params') { await delay(80); return json(defaultParams(q.get('plugin') || 'biomass')); }
    if (path === '/api/pan/params' && method === 'POST') {
      await delay(150);
      const info = parseParams(String(init.body));
      return info.error ? json({ error: info.error }, 400) : json(info);
    }
    if (path === '/api/pan/climate_preview' && method === 'POST') { await delay(250); return json(climatePreview(bodyJSON(init))); }
    if (path === '/api/simulate' && method === 'POST') {
      const b = bodyJSON(init);
      if (b.engine && b.engine !== 'fvs') {
        await delay(500);
        if (b.params_id && !knownParams[b.params_id] && !String(b.params_id).startsWith('mock-default')) return json({ error: 'unknown params_id ' + b.params_id }, 400);
        return json(simulate(b));
      }
    }
    return realFetch(input, init);
  };
  console.info('[mock_pan] Pan endpoints are mocked (?mock=1). Add &panoff=1 to simulate an unavailable Pan install.');
})();
