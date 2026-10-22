// Run: (npm i jsdom in a scratch dir) node apps/forest-explorer/test/ui_stands_jsdom.js apps/forest-explorer/public   — offline UI check of stands.js / proto.js (no server).
// jsdom test of the stands / prototype-layer interactions: real index.html markup, stands.js + proto.js, a stub maplibre map, real manifest + files.
const fs = require('fs'), path = require('path'); const { JSDOM } = require('jsdom'); const PUB = process.argv[2];
const html = fs.readFileSync(path.join(PUB, 'index.html'), 'utf8').replace(/<script[\s\S]*?<\/script>/g, '');
const dom = new JSDOM(html, { url: 'http://localhost/', pretendToBeVisual: true, runScripts: 'outside-only' }); const w = dom.window;
const handlers = {}, layers = [], sources = new Map(); let popups = 0;
const container = w.document.createElement('div'); w.document.body.appendChild(container);
w.map = { getContainer: () => container, addSource: (id, s) => sources.set(id, s), removeSource: id => sources.delete(id), getSource: id => sources.get(id),
  addLayer: l => layers.push(l), removeLayer: id => layers.splice(layers.findIndex(l => l.id === id), 1), getLayer: id => layers.find(l => l.id === id),
  getStyle: () => ({ layers }), setPaintProperty() {}, setFilter() {}, fitBounds() {}, getCanvas: () => ({ style: {} }),
  on: (ev, id, fn) => { (handlers[ev + ':' + id] ||= []).push(fn); } };
w.maplibregl = { Popup: class { setLngLat() { return this; } setHTML(h) { this.el = w.document.createElement('div'); this.el.innerHTML = h; return this; } addTo() { popups++; return this; } getElement() { return this.el; } remove() {} } };
w.fetch = async u => { const p = path.join(PUB, u.split('?')[0]); return { ok: fs.existsSync(p), json: async () => JSON.parse(fs.readFileSync(p, 'utf8')) }; };
w.eval('var vegSource = "pixels"; var aoiParts = []; var aoiStandIds = [];');
let aoi = null; w.setStandAoi = (ids, g) => { aoi = ids; }; w.setStatus = () => {}; w.startDraw = () => {}; w.clearAoi = () => {}; w.submitAoiUnion = () => {};
for (const f of ['stands.js', 'proto.js']) w.eval(fs.readFileSync(path.join(PUB, f), 'utf8') + `\n;window.${f === 'stands.js' ? 'StandsUI = StandsUI; window.InfoUI = InfoUI' : 'ProtoUI = ProtoUI'};`);
const $ = id => w.document.getElementById(id); const ok = (c, m) => { if (!c) { console.error('FAIL', m); process.exit(1); } console.log('ok  ', m); };
const click = (layerId, props, shift = false) => (handlers['click:' + layerId] || []).forEach(fn => fn({ features: [{ properties: props }], lngLat: { lng: 0, lat: 0 }, originalEvent: { shiftKey: shift } }));
const sleep = ms => new Promise(r => setTimeout(r, ms));
(async () => {
  w.StandsUI.init(); await w.ProtoUI.init(); await sleep(50);
  const areas = [...$('stArea').options].map(o => o.value); ok(areas.length > 1 && areas.some(a => a.startsWith('proto:')), `stand runs in the Stands area list: ${areas.length} (incl. prototype runs)`);
  ok(w.document.querySelectorAll('.infoModeSel').length === 2, 'an "Info on click" selector in both panels');
  // pixels mode: plain click -> popup info (default mode)
  $('stOn').checked = true; $('stOn').dispatchEvent(new w.Event('change')); await sleep(300);
  const fc = JSON.parse(fs.readFileSync(path.join(PUB, 'stands/baker_fl.geojson'))); const ids = fc.features.slice(0, 3).map(f => f.properties.stand_id);
  click('stands-fill', fc.features[0].properties); ok(popups === 1 && w.vegSource === 'pixels', 'plain click in pixel mode: info popup, no mode change');
  // info off -> no popup; box -> draggable box
  w.InfoUI.setMode('off'); click('stands-fill', fc.features[0].properties); ok(popups === 1, 'info mode off: no popup');
  w.InfoUI.setMode('box'); click('stands-fill', fc.features[1].properties); const box = container.querySelector('.infobox');
  ok(box && !box.hidden && box.textContent.includes(String(ids[1])), 'info box mode: box shows the stand');
  const h = box.querySelector('.infobox-h'); h.dispatchEvent(new w.MouseEvent('mousedown', { clientX: 10, clientY: 10, bubbles: true }));
  w.dispatchEvent(new w.MouseEvent('mousemove', { clientX: 60, clientY: 45 })); w.dispatchEvent(new w.MouseEvent('mouseup', {}));
  ok(box.style.left === '50px' && box.style.top === '35px', `info box dragged to ${box.style.left},${box.style.top}`);
  box.querySelector('.infobox-x').click(); ok(box.hidden, 'info box closes');
  // shift-click: switch to stands mode, build a multi-stand selection
  click('stands-fill', fc.features[0].properties, true); await sleep(700);
  ok(w.vegSource === 'stands' && $('stSelCount').textContent.startsWith('1 stand'), `shift-click: stands mode + selection "${$('stSelCount').textContent}"`);
  click('stands-fill', fc.features[1].properties, true); click('stands-fill', fc.features[2].properties); await sleep(700);
  ok($('stSelCount').textContent.startsWith('3 stands') && aoi && aoi.length === 3, `more stands: "${$('stSelCount').textContent}", AOI ${aoi && aoi.length} stands`);
  click('stands-fill', fc.features[1].properties, true); await sleep(700); ok(aoi.length === 2, 'shift-click again removes a stand');
  // shift-click on a prototype stand layer: Stands panel switches to that run and selects the stand
  const man = JSON.parse(fs.readFileSync(path.join(PUB, 'proto/manifest.json'))); const L = man.layers.find(l => l.id === 'stands_base');
  const pf = JSON.parse(fs.readFileSync(path.join(PUB, L.url))).features[5].properties;
  const cb = $('pxo-' + L.id); cb.checked = true; cb.dispatchEvent(new w.Event('change'));
  click('proto-' + L.id + '-f', pf, true); await sleep(1500);
  ok($('stArea').value === 'proto:' + L.id && aoi && aoi.includes(pf.stand_id), `prototype layer shift-click: Stands run ${$('stArea').value}, stand ${pf.stand_id} selected (${aoi && aoi.length})`);
  console.log('ALL UI CHECKS PASS');
})().catch(e => { console.error('FAIL', e); process.exit(1); });
