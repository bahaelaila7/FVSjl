'use strict';
/* Climate scenario editor (PnET-Succession engine only): modal with baseline / CO2 / change
   rows and live preview charts fetched from POST /api/pan/climate_preview. */

const ClimateUI = (() => {
  const esc = s => PanUI.esc(s);
  const VARS = [['tmean', 'Mean temperature'], ['tmin', 'Min temperature'], ['tmax', 'Max temperature'], ['ppt', 'Precipitation']];
  const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];
  const BASE_MODES = [['forecast', 'forecast (use supplied future climate)'], ['repeat', 'repeat last T years'], ['average', 'average of last K windows of T years']];
  const CO2_MODES = [['constant', 'constant'], ['linear', 'linear increase'], ['ssp245', 'SSP2-4.5'], ['ssp585', 'SSP5-8.5']];
  const COLORS = { base: '#94a3b8', mod: '#f59e0b' };

  const defaultClimate = () => ({ baseline: { mode: 'forecast', T: 30, K: 3 }, co2: { mode: 'constant', rate_ppm_per_yr: 2.0 }, changes: [] });
  const isDefault = c => !c || (c.baseline.mode === 'forecast' && c.co2.mode === 'constant' && !(c.changes || []).length);
  const clone = o => JSON.parse(JSON.stringify(o));
  const thisYear = () => new Date().getFullYear();

  function summary(c) {
    if (isDefault(c)) return 'default (forecast baseline, constant CO₂, no changes)';
    const b = c.baseline, parts = [];
    parts.push(b.mode === 'forecast' ? 'forecast baseline' : b.mode === 'repeat' ? `repeat last ${b.T}y` : `avg of last ${b.K}×${b.T}y`);
    const co = c.co2;
    parts.push(co.mode === 'constant' ? 'CO₂ constant' : co.mode === 'linear' ? `CO₂ +${co.rate_ppm_per_yr} ppm/yr` : 'CO₂ ' + co.mode.toUpperCase());
    (c.changes || []).forEach(ch => {
      const val = ch.kind === 'add' ? (ch.value >= 0 ? '+' : '') + ch.value : (ch.value >= 0 ? '+' : '') + ch.value + '%';
      const mo = ch.months && ch.months.length ? ` months ${ch.months.join(',')}` : '';
      parts.push(`${ch.var} ${val} ${ch.shape} ${ch.start_year}–${ch.end_year}${mo}`);
    });
    return parts.join('; ');
  }

  /* ---------------- modal state ---------------- */
  let work = defaultClimate(), modal = null, opener = null;
  const prev = { data: null, seq: 0, timer: null, err: '', loading: false, sel: 'tmean', ws: 0, len: 10 };

  function build() {
    modal = document.createElement('div');
    modal.id = 'climateModal'; modal.setAttribute('role', 'dialog'); modal.setAttribute('aria-modal', 'true'); modal.setAttribute('aria-labelledby', 'climTitle');
    modal.innerHTML = `<div class="box">
      <div class="row" style="align-items:center;gap:8px"><h2 id="climTitle">Climate scenario</h2>
        <span class="hint" style="margin-left:6px">PnET-Succession · applies to the AOI's ecoregions</span><span class="sp" style="flex:1"></span>
        <button class="mini" id="climReset" title="forecast baseline, constant CO₂, no changes">Reset</button>
        <button class="mini go" id="climApply">Apply</button>
        <button class="mini" id="climCancel" aria-label="close without applying">✕ Cancel</button></div>
      <div class="cols">
        <div id="climForm"></div>
        <div><h3>Preview — AOI-area-weighted means</h3><div id="climStatus" class="hint"></div>
          <div class="pgrid" id="climAnnual"></div>
          <h3>Monthly waveform</h3>
          <div class="fld"><label for="climVar">variable</label><select id="climVar"></select>
            <label for="climLen">window</label><select id="climLen"><option value="5">5 yr</option><option value="10" selected>10 yr</option><option value="20">20 yr</option><option value="40">40 yr</option></select>
            <label for="climWs">start</label><input type="range" id="climWs" min="0" max="0" value="0" style="flex:1;min-width:120px"><span id="climWsT" style="font:12px var(--mono)"></span></div>
          <div class="legendrow"><span><i class="sw" style="background:${COLORS.base}"></i> baseline</span><span><i class="sw" style="background:${COLORS.mod}"></i> modified</span><span class="readout" id="climRead" aria-live="off"></span></div>
          <div class="pcard"><svg id="climMonthly" role="img" aria-label="monthly baseline vs modified"></svg></div>
        </div></div></div>`;
    document.body.appendChild(modal);
    modal.addEventListener('click', e => { if (e.target === modal) cancel(); });
    modal.addEventListener('keydown', e => {
      if (e.key === 'Escape') { e.stopPropagation(); cancel(); return; }
      if (e.key === 'Tab') {            // light focus trap
        const f = [...modal.querySelectorAll('button,select,input,[tabindex]')].filter(x => !x.disabled && x.offsetParent !== null);
        if (!f.length) return;
        const first = f[0], last = f[f.length - 1];
        if (e.shiftKey && document.activeElement === first) { e.preventDefault(); last.focus(); }
        else if (!e.shiftKey && document.activeElement === last) { e.preventDefault(); first.focus(); }
      }
    });
    modal.querySelector('#climReset').onclick = () => { work = defaultClimate(); renderForm(); schedule(0); };
    modal.querySelector('#climApply').onclick = apply;
    modal.querySelector('#climCancel').onclick = cancel;
    modal.querySelector('#climVar').innerHTML = VARS.map(v => `<option value="${v[0]}">${v[1]}</option>`).join('');
    modal.querySelector('#climVar').onchange = e => { prev.sel = e.target.value; drawMonthly(); };
    modal.querySelector('#climLen').onchange = e => { prev.len = +e.target.value; clampWin(); drawAll(); };
    modal.querySelector('#climWs').oninput = e => { prev.ws = +e.target.value; drawAll(); };
  }

  function open() {
    if (!modal) build();
    work = PAN.climate ? clone(PAN.climate) : defaultClimate();
    opener = document.activeElement;
    modal.style.display = 'flex';
    renderForm(); drawAll(); schedule(0);
    const f = modal.querySelector('select,input,button'); if (f) f.focus();
  }
  function close() { modal.style.display = 'none'; if (opener && opener.focus) opener.focus(); }
  function cancel() { close(); }
  function apply() {
    const pr = problems(work);
    if (pr.any) { el('climStatus').textContent = 'Fix the highlighted fields before applying.'; return; }
    PAN.climate = isDefault(work) ? null : clone(work);
    renderSummary(); close();
    setStatus('climate: ' + summary(PAN.climate));
  }

  function renderSummary() {
    const s = el('climateSummary'); if (!s) return;
    const e = PanUI.engineById(PAN.engine);
    s.hidden = !(e && e.climate);
    s.textContent = isDefault(PAN.climate) ? 'climate: default' : 'climate: ' + summary(PAN.climate);
    s.title = s.textContent;
    s.style.maxWidth = '48ch'; s.style.overflow = 'hidden'; s.style.textOverflow = 'ellipsis'; s.style.whiteSpace = 'nowrap';
  }

  /* ---------------- validation ---------------- */
  function problems(c) {
    const p = { base: '', co2: '', rows: [], any: false };
    const b = c.baseline;
    if (b.mode !== 'forecast') {
      if (!(Number.isInteger(b.T) && b.T >= 1 && b.T <= 200)) p.base = 'T must be a whole number 1–200';
      else if (b.mode === 'average' && !(Number.isInteger(b.K) && b.K >= 1 && b.K <= 20)) p.base = 'K must be a whole number 1–20';
    }
    if (c.co2.mode === 'linear' && !(isFinite(c.co2.rate_ppm_per_yr) && c.co2.rate_ppm_per_yr >= -10 && c.co2.rate_ppm_per_yr <= 20)) p.co2 = 'rate must be −10…20 ppm/yr';
    c.changes.forEach(ch => {
      const m = {};
      if (!isFinite(ch.value)) m.value = 'value required';
      else if (ch.kind === 'mul' && ch.value <= -100) m.value = 'percent change must be > −100';
      if (!(Number.isInteger(ch.start_year) && ch.start_year >= 1900 && ch.start_year <= 2400)) m.start_year = 'start year 1900–2400';
      if (!(Number.isInteger(ch.end_year) && ch.end_year >= ch.start_year)) m.end_year = 'end year must be ≥ start';
      if (ch.months && !ch.months.length) m.months = 'pick at least one month';
      p.rows.push(m);
    });
    p.any = !!(p.base || p.co2 || p.rows.some(m => Object.keys(m).length));
    return p;
  }

  /* ---------------- form ---------------- */
  const opt = (items, sel) => items.map(([v, t]) => `<option value="${v}"${v === sel ? ' selected' : ''}>${esc(t)}</option>`).join('');

  function renderForm() {
    const host = modal.querySelector('#climForm'), pr = problems(work), b = work.baseline;
    const bad = (m, k) => (m && m[k] ? ' bad' : '');
    host.innerHTML = `<h3>Baseline climate</h3><div class="cc${pr.base ? ' invalid' : ''}">
        <div class="fld"><label for="cBm">mode</label><select id="cBm">${opt(BASE_MODES, b.mode)}</select></div>
        ${b.mode !== 'forecast' ? `<div class="fld"><label for="cT">T (years per window)</label><input type="number" id="cT" min="1" max="200" step="1" value="${esc(b.T)}" class="${pr.base ? 'bad' : ''}">
          ${b.mode === 'average' ? `<label for="cK">K (windows)</label><input type="number" id="cK" min="1" max="20" step="1" value="${esc(b.K)}" class="${pr.base ? 'bad' : ''}">` : ''}</div>` : ''}
        <div class="err">${esc(pr.base)}</div></div>
      <h3>Atmospheric CO₂</h3><div class="cc${pr.co2 ? ' invalid' : ''}">
        <div class="fld"><label for="cCm">pathway</label><select id="cCm">${opt(CO2_MODES, work.co2.mode)}</select>
        ${work.co2.mode === 'linear' ? `<label for="cCr">rate (ppm/yr)</label><input type="number" id="cCr" step="any" value="${esc(work.co2.rate_ppm_per_yr)}" class="${pr.co2 ? 'bad' : ''}">` : ''}</div>
        <div class="err">${esc(pr.co2)}</div></div>
      <h3>Changes <button class="mini" id="cAdd" style="margin-left:8px">+ Add change</button></h3>
      <div id="cRows">${work.changes.length ? '' : '<div class="hint">No changes — the modified series equals the baseline.</div>'}
      ${work.changes.map((ch, i) => rowHTML(ch, i, pr.rows[i] || {})).join('')}</div>
      <div class="hint">Add: absolute offset (°C or mm/month). Mul: percent change (e.g. −20 = 20% less; 0 = none). Step: the change applies from start to end year; ramp: grows linearly from 0 at start to its full value at end. Months: untick to restrict to some months.</div>`;
    bindForm(host);
  }

  function rowHTML(ch, i, m) {
    const u = ch.var === 'ppt' ? (ch.kind === 'add' ? 'mm/month' : '%') : (ch.kind === 'add' ? '°C' : '%');
    const sel = ch.months ? new Set(ch.months) : null;
    return `<div class="cc${Object.keys(m).length ? ' invalid' : ''}" data-row="${i}">
      <div class="fld"><label for="r${i}v">variable</label><select id="r${i}v" data-k="var">${opt(VARS.map(v => [v[0], v[1]]), ch.var)}</select>
        <select data-k="kind" aria-label="add or multiply">${opt([['add', 'add'], ['mul', 'mul (% change)']], ch.kind)}</select>
        <input type="number" step="any" data-k="value" value="${esc(ch.value)}" aria-label="value" class="${m.value ? 'bad' : ''}" title="${esc(m.value || '')}"><span style="color:var(--muted)">${u}</span>
        <select data-k="shape" aria-label="shape">${opt([['step', 'step'], ['ramp', 'ramp']], ch.shape)}</select>
        <span style="flex:1"></span><button class="mini" data-rmrow="${i}" aria-label="remove change ${i + 1}">remove</button></div>
      <div class="fld"><label>from</label><input type="number" step="1" data-k="start_year" value="${esc(ch.start_year)}" aria-label="start year" style="width:80px" class="${m.start_year ? 'bad' : ''}" title="${esc(m.start_year || '')}">
        <label>to</label><input type="number" step="1" data-k="end_year" value="${esc(ch.end_year)}" aria-label="end year" style="width:80px" class="${m.end_year ? 'bad' : ''}" title="${esc(m.end_year || '')}"></div>
      <div class="fld"><label>months</label><div class="months" role="group" aria-label="months affected">${MONTHS.map((n, k) =>
        `<label title="${n}"><input type="checkbox" data-mo="${k + 1}"${!sel || sel.has(k + 1) ? ' checked' : ''}>${n[0]}</label>`).join('')}</div>
        <button class="mini" data-allmo="${i}" type="button">all</button><span class="err">${esc(m.months || '')}</span></div></div>`;
  }

  function bindForm(host) {
    const $ = id => host.querySelector('#' + id);
    $('cBm').onchange = e => { work.baseline.mode = e.target.value; renderForm(); schedule(); };
    if ($('cT')) $('cT').oninput = e => { work.baseline.T = e.target.value === '' ? NaN : +e.target.value; quickValidate(); schedule(); };
    if ($('cK')) $('cK').oninput = e => { work.baseline.K = e.target.value === '' ? NaN : +e.target.value; quickValidate(); schedule(); };
    $('cCm').onchange = e => { work.co2.mode = e.target.value; renderForm(); schedule(); };
    if ($('cCr')) $('cCr').oninput = e => { work.co2.rate_ppm_per_yr = e.target.value === '' ? NaN : +e.target.value; quickValidate(); schedule(); };
    $('cAdd').onclick = () => {
      const y = thisYear();
      work.changes.push({ var: 'tmean', kind: 'add', value: 2.0, shape: 'ramp', start_year: y + 5, end_year: y + 55, months: null });
      renderForm(); schedule();
    };
    host.querySelectorAll('[data-row]').forEach(rowEl => {
      const i = +rowEl.dataset.row, ch = work.changes[i];
      rowEl.querySelectorAll('[data-k]').forEach(inp => {
        const k = inp.dataset.k;
        const h = () => {
          if (['value', 'start_year', 'end_year'].includes(k)) ch[k] = inp.value === '' ? NaN : +inp.value; else ch[k] = inp.value;
          if (['var', 'kind'].includes(k)) { renderForm(); schedule(); } else { quickValidate(); schedule(); }
        };
        inp.addEventListener(inp.tagName === 'SELECT' ? 'change' : 'input', h);
      });
      rowEl.querySelectorAll('[data-mo]').forEach(cb => cb.onchange = () => {
        const on = [...rowEl.querySelectorAll('[data-mo]')].filter(x => x.checked).map(x => +x.dataset.mo);
        ch.months = on.length === 12 ? null : on;
        quickValidate(); schedule();
      });
      rowEl.querySelector('[data-allmo]').onclick = () => { ch.months = null; renderForm(); schedule(); };
      rowEl.querySelector('[data-rmrow]').onclick = () => { work.changes.splice(i, 1); renderForm(); schedule(); };
    });
  }
  // flag invalid fields without re-rendering (keeps focus while typing)
  function quickValidate() {
    const pr = problems(work), host = modal.querySelector('#climForm');
    work.changes.forEach((ch, i) => {
      const row = host.querySelector(`[data-row="${i}"]`); if (!row) return;
      const m = pr.rows[i] || {};
      row.classList.toggle('invalid', Object.keys(m).length > 0);
      row.querySelectorAll('input[data-k]').forEach(inp => { inp.classList.toggle('bad', !!m[inp.dataset.k]); inp.title = m[inp.dataset.k] || ''; });
      const er = row.querySelector('.err'); if (er) er.textContent = m.months || '';
    });
    ['cT', 'cK'].forEach(id => { const x = host.querySelector('#' + id); if (x) x.classList.toggle('bad', !!pr.base); });
    const cr = host.querySelector('#cCr'); if (cr) cr.classList.toggle('bad', !!pr.co2);
  }

  /* ---------------- preview ---------------- */
  function schedule(delay = 450) {
    clearTimeout(prev.timer);
    prev.timer = setTimeout(fetchPreview, delay);
  }
  async function fetchPreview() {
    const st = modal.querySelector('#climStatus');
    if (!currentAoiGeom) { prev.data = null; st.textContent = 'Select an AOI to preview its climate (the preview is the AOI-area-weighted mean over its level-4 ecoregions).'; drawAll(); return; }
    const pr = problems(work);
    if (pr.any) { st.textContent = 'Fix the highlighted fields to refresh the preview.'; return; }
    const seq = ++prev.seq;
    st.textContent = 'loading preview…';
    const years = Math.min(200, Math.max(30, (+el('pCycles').value || 10) * (+el('pPeriod').value || 10)));
    try {
      const r = await fetch('/api/pan/climate_preview', { method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ geometry: currentAoiGeom, climate: work, start_year: thisYear(), years }) });
      const d = await r.json().catch(() => ({ error: 'bad response (' + r.status + ')' }));
      if (seq !== prev.seq) return;                       // superseded
      if (!r.ok || d.error) { prev.data = null; st.textContent = 'Preview error: ' + (d.error || r.status); drawAll(); return; }
      prev.data = d;
      st.textContent = `ecoregions (L4): ${(d.l4s || []).join(', ') || '—'} · ${d.start_year}–${d.start_year + ((d.vars[0] && d.vars[0].years.length) || years) - 1}`;
      if (!d.vars.some(v => v.name === prev.sel) && d.vars[0]) prev.sel = d.vars[0].name;
      modal.querySelector('#climVar').innerHTML = d.vars.map(v => `<option value="${esc(v.name)}"${v.name === prev.sel ? ' selected' : ''}>${esc(v.label || v.name)}</option>`).join('');
      clampWin(); drawAll();
    } catch (err) {
      if (seq === prev.seq) { prev.data = null; st.textContent = 'Preview error: ' + err.message; drawAll(); }
    }
  }
  function nYears() { const v = prev.data && prev.data.vars[0]; return v ? v.years.length : 0; }
  function clampWin() {
    const n = nYears(); prev.len = Math.min(prev.len, Math.max(1, n));
    prev.ws = Math.max(0, Math.min(prev.ws, Math.max(0, n - prev.len)));
  }
  function drawAll() {
    const ws = modal.querySelector('#climWs'), n = nYears();
    ws.max = Math.max(0, n - prev.len); ws.value = prev.ws; ws.disabled = !n;
    const v0 = prev.data && prev.data.vars[0];
    modal.querySelector('#climWsT').textContent = v0 ? `${v0.years[prev.ws]}–${v0.years[Math.min(n - 1, prev.ws + prev.len - 1)]}` : '';
    drawAnnual(); drawMonthly();
  }

  const fmtN = x => (Math.abs(x) >= 100 ? x.toFixed(0) : Math.abs(x) >= 10 ? x.toFixed(1) : x.toFixed(2));
  function ticks(lo, hi, n = 4) { const o = []; for (let i = 0; i <= n; i++) o.push(lo + (hi - lo) * i / n); return o; }

  function drawAnnual() {
    const host = modal.querySelector('#climAnnual');
    if (!prev.data) { host.innerHTML = ''; return; }
    host.innerHTML = prev.data.vars.map(v => {
      const W = 300, H = 140, mL = 40, mR = 6, mT = 6, mB = 20, pw = W - mL - mR, ph = H - mT - mB;
      const a = v.baseline_annual, b = v.modified_annual, n = v.years.length;
      const lo0 = Math.min(...a, ...b), hi0 = Math.max(...a, ...b), pad = (hi0 - lo0) * 0.08 || 1, lo = lo0 - pad, hi = hi0 + pad;
      const X = i => mL + pw * (n > 1 ? i / (n - 1) : 0.5), Y = y => mT + ph - ph * (y - lo) / (hi - lo);
      const path = arr => 'M' + arr.map((y, i) => X(i).toFixed(1) + ',' + Y(y).toFixed(1)).join(' L');
      let g = ticks(lo, hi).map(t => `<line x1="${mL}" x2="${W - mR}" y1="${Y(t)}" y2="${Y(t)}" stroke="#26332b"/><text x="${mL - 4}" y="${Y(t) + 3}" text-anchor="end" font-size="9" fill="#8ba393">${fmtN(t)}</text>`).join('');
      const step = Math.max(1, Math.ceil(n / 5));
      for (let i = 0; i < n; i += step) g += `<text x="${X(i)}" y="${H - 5}" text-anchor="middle" font-size="9" fill="#8ba393">${v.years[i]}</text>`;
      const wx0 = X(prev.ws), wx1 = X(Math.min(n - 1, prev.ws + prev.len - 1));
      g += `<rect x="${wx0}" y="${mT}" width="${Math.max(2, wx1 - wx0)}" height="${ph}" fill="#4ade80" opacity=".12"/>`;
      g += `<path d="${path(a)}" fill="none" stroke="${COLORS.base}" stroke-width="1.6"/><path d="${path(b)}" fill="none" stroke="${COLORS.mod}" stroke-width="1.8"/>`;
      return `<div class="pcard"><div class="t">${esc(v.label || v.name)} (${esc(v.units || '')}) — annual mean · click to move the monthly window</div>
        <svg viewBox="0 0 ${W} ${H}" data-annual="${esc(v.name)}" role="img" aria-label="${esc(v.label || v.name)} annual baseline vs modified" style="cursor:pointer">${g}</svg></div>`;
    }).join('');
    host.querySelectorAll('svg[data-annual]').forEach(svg => svg.onclick = e => {
      const r = svg.getBoundingClientRect(), W = 300, mL = 40, mR = 6;
      const fx = ((e.clientX - r.left) / r.width * W - mL) / (W - mL - mR), n = nYears();
      prev.ws = Math.round(fx * (n - 1) - prev.len / 2); prev.sel = svg.dataset.annual;
      modal.querySelector('#climVar').value = prev.sel;
      clampWin(); drawAll();
    });
  }

  function drawMonthly() {
    const svg = modal.querySelector('#climMonthly'), rd = modal.querySelector('#climRead');
    const v = prev.data && prev.data.vars.find(x => x.name === prev.sel);
    if (!v) { svg.setAttribute('viewBox', '0 0 10 10'); svg.innerHTML = ''; svg.style.height = '0'; rd.textContent = ''; return; }
    svg.style.height = ''; svg.style.width = '100%';
    const W = 720, H = 230, mL = 46, mR = 8, mT = 8, mB = 24, pw = W - mL - mR, ph = H - mT - mB;
    const m0 = prev.ws * 12, m1 = Math.min(v.baseline_monthly.length, (prev.ws + prev.len) * 12), nm = m1 - m0;
    const a = v.baseline_monthly.slice(m0, m1), b = v.modified_monthly.slice(m0, m1);
    const lo0 = Math.min(...a, ...b), hi0 = Math.max(...a, ...b), pad = (hi0 - lo0) * 0.08 || 1, lo = lo0 - pad, hi = hi0 + pad;
    const X = i => mL + pw * (nm > 1 ? i / (nm - 1) : 0.5), Y = y => mT + ph - ph * (y - lo) / (hi - lo);
    const path = arr => 'M' + arr.map((y, i) => X(i).toFixed(1) + ',' + Y(y).toFixed(1)).join(' L');
    let g = ticks(lo, hi).map(t => `<line x1="${mL}" x2="${W - mR}" y1="${Y(t)}" y2="${Y(t)}" stroke="#26332b"/><text x="${mL - 4}" y="${Y(t) + 3}" text-anchor="end" font-size="10" fill="#8ba393">${fmtN(t)}</text>`).join('');
    const ystep = Math.max(1, Math.ceil(prev.len / 10));
    for (let k = 0; k * 12 < nm; k++) {
      const x = X(k * 12);
      g += `<line x1="${x}" x2="${x}" y1="${mT}" y2="${mT + ph}" stroke="#26332b" stroke-dasharray="2 3"/>`;
      if (k % ystep === 0) g += `<text x="${x}" y="${H - 7}" font-size="10" fill="#8ba393" text-anchor="start">${v.years[prev.ws + k]}</text>`;
    }
    g += `<path d="${path(a)}" fill="none" stroke="${COLORS.base}" stroke-width="1.6"/><path d="${path(b)}" fill="none" stroke="${COLORS.mod}" stroke-width="1.8"/>`;
    g += `<line id="climCur" x1="0" x2="0" y1="${mT}" y2="${mT + ph}" stroke="#facc15" stroke-width="1" visibility="hidden"/>`;
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`); svg.innerHTML = g;
    svg.onmousemove = e => {
      const r = svg.getBoundingClientRect(); const fx = ((e.clientX - r.left) / r.width * W - mL) / pw;
      const i = Math.max(0, Math.min(nm - 1, Math.round(fx * (nm - 1))));
      const cur = svg.querySelector('#climCur'); cur.setAttribute('x1', X(i)); cur.setAttribute('x2', X(i)); cur.setAttribute('visibility', 'visible');
      const yr = v.years[prev.ws + Math.floor(i / 12)];
      rd.textContent = `${MONTHS[i % 12]} ${yr}: baseline ${fmtN(a[i])} → modified ${fmtN(b[i])} ${v.units || ''}`;
    };
    svg.onmouseleave = () => { const c = svg.querySelector('#climCur'); if (c) c.setAttribute('visibility', 'hidden'); rd.textContent = ''; };
  }

  function init() {
    const b = el('climateBtn');
    if (b) b.onclick = open;
    renderSummary();
  }

  return { init, open, summary, isDefault, defaultClimate, renderSummary, problems };
})();
