'use strict';

const CMAP_CSS = {
  greens:  'linear-gradient(90deg,#f7fcf5,#c7e9c0,#74c476,#238b45,#00441b)',
  viridis: 'linear-gradient(90deg,#440154,#414487,#2a788e,#22a884,#7ad151,#fde725)',
  magma:   'linear-gradient(90deg,#000004,#3b0f70,#8c2981,#de4968,#fe9f6d,#fcfdbf)',
};
// species palette for the chart (top species get a color; rest -> "other")
const SP_PALETTE = ['#4ade80','#38bdf8','#f59e0b','#f472b6','#a78bfa',
                    '#facc15','#fb7185','#2dd4bf','#94a3b8'];

// key-free raster basemaps (Esri/USGS ArcGIS use z/y/x order; OSM uses z/x/y)
const BASEMAPS = {
  satellite: { tiles: ['https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'],
               maxzoom: 19, attribution: 'Esri World Imagery · USFS TreeMap 2022' },
  topo:      { tiles: ['https://basemap.nationalmap.gov/arcgis/rest/services/USGSImageryTopo/MapServer/tile/{z}/{y}/{x}'],
               maxzoom: 16, attribution: 'USGS The National Map · USFS TreeMap 2022' },
  terrain:   { tiles: ['https://server.arcgisonline.com/ArcGIS/rest/services/World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}'],
               maxzoom: 13, attribution: 'Esri World Shaded Relief · USFS TreeMap 2022' },
  osm:       { tiles: ['https://tile.openstreetmap.org/{z}/{x}/{y}.png'],
               maxzoom: 19, attribution: '© OpenStreetMap contributors · USFS TreeMap 2022' },
};
let curBasemap = 'satellite';

const el = id => document.getElementById(id);
const status = el('status');
const fmt = (n, d=0) => Number(n).toLocaleString(undefined,{maximumFractionDigits:d});

// never fail silently — surface any uncaught error to the status bar + console
window.addEventListener('error', e => { try { status.textContent = 'JS error: ' + e.message; } catch(_){} });
window.addEventListener('unhandledrejection', e => { try { status.textContent = 'error: ' + (e.reason && e.reason.message || e.reason); } catch(_){} });

let map, attrs = {}, curAttr = 'carbon_l';
const RASTER_SRC = 'treemap';

async function boot() {
  const meta = await fetch('/api/meta').then(r=>r.json());
  const list = await fetch('/api/attrs').then(r=>r.json());
  list.forEach(a => attrs[a.name] = a);
  curAttr = meta.default_attr;

  const sel = el('attr');
  list.forEach(a => {
    const o = document.createElement('option');
    o.value = a.name; o.textContent = `${a.label}${a.units?` (${a.units})`:''}`;
    if (a.name === curAttr) o.selected = true;
    sel.appendChild(o);
  });

  // source & engine selectors (extension points; one option each for now)
  const fill=(id,items)=>{ const s=el(id); (items||[]).forEach(it=>{ const o=document.createElement('option'); o.value=it.id; o.textContent=it.name; o.title=it.detail||''; s.appendChild(o); }); };
  fill('source', meta.sources); fill('engine', meta.engines);
  const src0=(meta.sources||[])[0], eng0=(meta.engines||[])[0];
  el('subline').textContent = [src0&&src0.detail, eng0&&('engine: '+eng0.name)].filter(Boolean).join(' · ');

  const EMPTY = {type:'FeatureCollection',features:[]};
  map = new maplibregl.Map({
    container: 'map',
    style: {
      version: 8,
      sources: {
        base: { type:'raster', tileSize:256, tiles:BASEMAPS[curBasemap].tiles,
                maxzoom:BASEMAPS[curBasemap].maxzoom, attribution:BASEMAPS[curBasemap].attribution },
        [RASTER_SRC]: { type:'raster', tileSize:256, tiles:[tileUrl(curAttr)] },
        aoi:   { type:'geojson', data:EMPTY },
        draft: { type:'geojson', data:EMPTY }
      },
      // declared here (not via addLayer) so the draw layers can never be missing
      layers: [
        { id:'base', type:'raster', source:'base' },
        { id:'treemap', type:'raster', source:RASTER_SRC, layout:{visibility:'none'}, paint:{'raster-opacity':0.85} },
        { id:'aoi-fill', type:'fill', source:'aoi',
          paint:{'fill-color':'#4ade80','fill-opacity':0.12} },
        { id:'aoi-line', type:'line', source:'aoi',
          paint:{'line-color':'#4ade80','line-width':2} },
        { id:'draft-fill', type:'fill', source:'draft',
          paint:{'fill-color':'#facc15','fill-opacity':0.10} },
        { id:'draft-line', type:'line', source:'draft',
          paint:{'line-color':'#facc15','line-width':2,'line-dasharray':[2,1]} },
        { id:'draft-pt', type:'circle', source:'draft',
          filter:['==',['geometry-type'],'Point'],
          paint:{
            'circle-radius':['case',['boolean',['get','first'],false],7,5],
            'circle-color':['case',['boolean',['get','first'],false],'#22c55e','#facc15'],
            'circle-stroke-width':2,'circle-stroke-color':'#fff'} }
      ]
    },
    center: [(meta.bounds.west+meta.bounds.east)/2, (meta.bounds.south+meta.bounds.north)/2],
    zoom: 3.4
  });
  map.addControl(new maplibregl.NavigationControl({showCompass:false}), 'top-right');
  map.on('error', e => { console.error('map error', e && e.error); setStatus('map error: '+((e&&e.error&&e.error.message)||'see console')); });
  map.on('load', () => {
    map.fitBounds([[meta.bounds.west,meta.bounds.south],[meta.bounds.east,meta.bounds.north]],
                  {padding:20, duration:0});
    setAttrVisible(attrVisible());     // sync: off by default -> no whole-map tiles fetched
    updateLegend();
    setStatus(`${fmt(meta.nplots_total)} plots · draw or upload an AOI to compute stats`);
  });

  sel.onchange = e => { curAttr = e.target.value; swapTiles(); updateLegend();
                        setStatus(attrs[curAttr].label); };
  el('op').oninput = e => {
    const v = +e.target.value; el('opv').textContent = v+'%';
    if (map.getLayer('treemap')) map.setPaintProperty('treemap','raster-opacity', v/100);
  };
  el('showAttr').onchange = e => setAttrVisible(e.target.checked);
  el('basemap').onchange = e => setBasemap(e.target.value);
  el('addThin').onclick = addThin;
  el('addPlant').onclick = addPlant;
  el('savePlan').onclick = savePlan;
  el('loadPlan').onclick = () => el('planFile').click();
  el('planFile').onchange = loadPlan;
  el('runSim').onclick = runSim;
  el('simMetric').onchange = () => { drawSimChart(); if(lastSim) setSimLayer(); };
  el('simSlider').oninput = e => onSlider(e.target.value);
  el('simPlay').onclick = toggleSimPlay;
  el('showSim').onchange = setSimLayer;
  el('simMapMetric').onchange = () => { renderSimLegend(); setSimLayer(); };
  el('simOp').oninput = e => { el('simOpv').textContent=e.target.value+'%';
    if(map.getLayer('simlayer')) map.setPaintProperty('simlayer','raster-opacity',(+e.target.value)/100); };
  el('pCycles').onchange = e => plan.cycles = +e.target.value;
  el('pPeriod').onchange = e => plan.period = +e.target.value;
  renderActions();
  el('draw').onclick = startDraw;
  el('clear').onclick = clearAoi;
  el('download').onclick = downloadAoi;
  el('upload').onclick = () => el('file').click();
  el('file').onchange = onUpload;
  window.addEventListener('keydown', e => {
    if (!drawing) return;
    if (e.key === 'Escape') cancelDraw();
    else if (e.key === 'Enter') finishDraw();
    else if (e.key === 'Backspace') { e.preventDefault(); verts.pop(); renderDraft(); updateDrawHint(); }
  });
}

function setBasemap(key){
  curBasemap = key;
  const b = BASEMAPS[key];
  if (map.getLayer('base')) map.removeLayer('base');
  if (map.getSource('base')) map.removeSource('base');
  map.addSource('base',{type:'raster',tileSize:256,tiles:b.tiles,maxzoom:b.maxzoom,attribution:b.attribution});
  // keep base at the very bottom — insert before the lowest non-base layer (a bare
  // `undefined` beforeId would stack it on TOP, hiding the AOI, sim raster and overlays).
  const first = map.getStyle().layers.find(l => l.id !== 'base');
  map.addLayer({id:'base',type:'raster',source:'base'}, first ? first.id : undefined);
}

const tileUrl = a => `/tiles/${a}/{z}/{x}/{y}.png`;
const attrVisible = () => el('showAttr').checked;
function swapTiles(){
  if (map.getLayer('treemap')) map.removeLayer('treemap');
  if (map.getSource(RASTER_SRC)) map.removeSource(RASTER_SRC);
  if (!attrVisible()) return;                 // don't fetch whole-map tiles when off
  map.addSource(RASTER_SRC,{type:'raster',tileSize:256,tiles:[tileUrl(curAttr)]});
  map.addLayer({id:'treemap',type:'raster',source:RASTER_SRC,
    paint:{'raster-opacity':(+el('op').value)/100}}, aoiBeforeId());
}
function setAttrVisible(on){
  if(on){ if(!map.getLayer('treemap')) swapTiles();
          else map.setLayoutProperty('treemap','visibility','visible'); }
  else if(map.getLayer('treemap')){ map.removeLayer('treemap'); if(map.getSource(RASTER_SRC)) map.removeSource(RASTER_SRC); }
}
function updateLegend(){
  const a = attrs[curAttr];
  el('legbar').style.background = CMAP_CSS[a.cmap] || CMAP_CSS.viridis;
  el('legLo').textContent = fmt(a.lo);
  el('legHi').textContent = fmt(a.hi);
  el('legLabel').textContent = a.units || '';
}
function setStatus(t){ status.textContent = t; }

/* AOI + draft layers are declared in the initial style (see boot). */
const aoiBeforeId = () => (map && map.getLayer('aoi-fill')) ? 'aoi-fill' : undefined;
const emptyFC = () => ({type:'FeatureCollection',features:[]});

/* ---------------- polygon draw ---------------- */
let drawing=false, verts=[], lastClick=0;
function startDraw(){
  clearAoi();                 // clears any prior AOI (also calls cancelDraw)
  drawing=true; verts=[];
  map.getCanvas().style.cursor='crosshair';
  map.doubleClickZoom.disable();
  map.on('click', onDrawClick);
  map.on('dblclick', onDrawDblClick);
  map.on('mousemove', onDrawMove);
  updateDrawHint();
  setStatus('draw mode — click on the map to place points');
  console.log('[draw] mode on; draft layer present:', !!(map.getLayer && map.getLayer('draft-pt')));
}
function onDrawClick(e){
  if (!drawing) return;
  const now = Date.now();
  // swallow the 2nd click of a double-click (dblclick handler finishes instead)
  if (now - lastClick < 300){ lastClick = now; return; }
  lastClick = now;
  const p = [e.lngLat.lng, e.lngLat.lat];
  // click near the first vertex closes the polygon
  if (verts.length >= 3){
    const f = map.project(verts[0]), c = map.project(p);
    if (Math.hypot(f.x-c.x, f.y-c.y) < 12){ finishDraw(); return; }
  }
  verts.push(p); renderDraft(); updateDrawHint();
  setStatus(`drawing — ${verts.length} point${verts.length===1?'':'s'} placed`);
}
function onDrawDblClick(e){ if (e) e.preventDefault(); finishDraw(); }
function onDrawMove(e){ if (drawing && verts.length) renderDraft([e.lngLat.lng,e.lngLat.lat]); }

function renderDraft(hover){
  if (!map.getSource('draft')) return;
  // fold the moving cursor into the preview so the rubber-band stays live at every
  // vertex count (>=3 verts + hover -> the closing polygon tracks the cursor too)
  const chain = hover ? verts.concat([hover]) : verts;
  const feats = verts.map((v,i)=>({type:'Feature',properties:{first:i===0},
                                    geometry:{type:'Point',coordinates:v}}));
  if (chain.length >= 3) feats.push({type:'Feature',geometry:
    {type:'Polygon',coordinates:[chain.concat([chain[0]])]}});
  else if (chain.length >= 2) feats.push({type:'Feature',geometry:
    {type:'LineString',coordinates:chain}});
  map.getSource('draft').setData({type:'FeatureCollection',features:feats});
}
function updateDrawHint(){
  const h = el('drawhint');
  h.style.display = 'block';
  h.textContent = verts.length < 3
    ? `Click to add points — ${verts.length} placed (need 3+). Esc to cancel.`
    : `${verts.length} points · double-click, Enter, or click the first point to finish · Backspace to undo · Esc to cancel`;
}
function finishDraw(){
  if (verts.length < 3){ return; }         // keep drawing until a valid ring
  const geometry = {type:'Polygon',coordinates:[verts.concat([verts[0]])]};
  cancelDraw();
  currentAoiGeom = geometry;               // enable GeoJSON download
  el('download').disabled = false;
  map.getSource('aoi').setData({type:'Feature',geometry});
  submitAoi({body:JSON.stringify({geometry}),
             headers:{'Content-Type':'application/json'}});
}

let currentAoiGeom = null;
function downloadAoi(){
  if (!currentAoiGeom) return;
  const fc = {type:'FeatureCollection', features:[{type:'Feature',
              properties:{name:'AOI', source:'TreeMap Growth Explorer',
                          created:new Date().toISOString()},
              geometry:currentAoiGeom}]};
  const blob = new Blob([JSON.stringify(fc, null, 2)], {type:'application/geo+json'});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `aoi-${new Date().toISOString().slice(0,19).replace(/[:T]/g,'-')}.geojson`;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}
function cancelDraw(){
  drawing=false; verts=[]; lastClick=0;
  map.getCanvas().style.cursor='';
  el('drawhint').style.display='none';
  map.off('click',onDrawClick); map.off('dblclick',onDrawDblClick); map.off('mousemove',onDrawMove);
  map.doubleClickZoom.enable();
  if (map.getSource('draft')) map.getSource('draft').setData(emptyFC());
}
function clearAoi(){
  cancelDraw();
  if (map.getSource('aoi')) map.getSource('aoi').setData(emptyFC());
  currentAoiGeom = null; el('download').disabled = true;
  el('runSim').disabled = true;
  lastSim = null; clearSimLayer(); drawSimChart();
  el('statsGrp').hidden = true; el('chartGrp').hidden = true;
}

/* ---------------- upload ---------------- */
async function onUpload(e){
  const f = e.target.files[0]; if(!f) return;
  const name = f.name.toLowerCase();
  if (name.endsWith('.shp') || name.endsWith('.shx') || name.endsWith('.dbf') || name.endsWith('.prj')){
    setStatus('Upload a .zip of the whole shapefile (.shp+.shx+.dbf+.prj), or a GeoJSON — a single .shp can’t be read.');
    e.target.value=''; return;
  }
  clearAoi();
  const buf = await f.arrayBuffer();
  submitAoi({body:buf, headers:{'Content-Type':'application/octet-stream','X-Filename':f.name}});
  e.target.value='';
}

/* ---------------- submit + render ---------------- */
async function submitAoi(opts){
  setStatus('resolving AOI…');
  try {
    const res = await fetch('/api/aoi',{method:'POST',...opts});
    const data = await res.json();
    if (!res.ok || data.error){ setStatus('AOI error: '+(data.error||res.status)); return; }
    renderResult(data);
  } catch(err){ setStatus('AOI error: '+err.message); }
}

function renderResult(d){
  // draw the AOI the server actually used (covers uploads, which have no client geometry)
  if (d.geometry && map.getSource('aoi')){
    map.getSource('aoi').setData({type:'Feature',geometry:d.geometry});
    currentAoiGeom = d.geometry;            // allow re-download of an uploaded AOI too
    el('download').disabled = false;
  }
  if (d.bbox) map.fitBounds([[d.bbox.west,d.bbox.south],[d.bbox.east,d.bbox.north]],
                            {padding:60,maxZoom:13,duration:600});
  const cards = [
    ['Plots', fmt(d.nplots), ''],
    ['Area', fmt(d.acres), 'acres'],
    ['Live carbon', fmt(d.carbon_l), 'tons AG'],
    ['Live biomass', fmt(d.drybio_l), 'tons AG'],
    ['Live volume', fmt(d.volcfnet_l), 'ft³ net'],
    ['Basal area', fmt(d.ba_total/Math.max(d.acres,1),1), 'ft²/ac'],
  ];
  el('cards').innerHTML = cards.map(([k,v,u])=>
    `<div class="card"><div class="k">${k}</div><div class="v">${v}</div><div class="u">${u}</div></div>`).join('');
  el('statsGrp').hidden = false;
  if(el('cycleBadge')) el('cycleBadge').textContent = 'cycle 0 · current';
  drawChart(d.cells);
  onAoiReady(d.cells);
  setStatus(`${fmt(d.nplots)} plots · ${fmt(d.acres)} ac · ${fmt(d.npixels)} px`);
}

/* species × DBH stacked bar of basal area */
function drawChart(cells){
  const svg = el('chart');
  svg.innerHTML='';
  if(!cells || !cells.length){ el('chartGrp').hidden=true; return; }
  el('chartGrp').hidden=false;

  // top species by total BA
  const spBA={}, spName={};
  cells.forEach(c=>{ spBA[c.symbol]=(spBA[c.symbol]||0)+c.ba; spName[c.symbol]=c.common; });
  const top = Object.keys(spBA).sort((a,b)=>spBA[b]-spBA[a]).slice(0,SP_PALETTE.length-1);
  const color={}; top.forEach((s,i)=>color[s]=SP_PALETTE[i]);
  const OTHER='#64748b';
  const spOf = s => top.includes(s)?s:'other';
  const colOf = s => top.includes(s)?color[s]:OTHER;

  // classes present, ascending
  const classes=[...new Set(cells.map(c=>c.dbh_lo))].sort((a,b)=>a-b);
  const stacks={}; classes.forEach(c=>stacks[c]={});
  cells.forEach(c=>{ const s=spOf(c.symbol); stacks[c.dbh_lo][s]=(stacks[c.dbh_lo][s]||0)+c.ba; });
  const totals=classes.map(c=>Object.values(stacks[c]).reduce((a,b)=>a+b,0));
  const ymax=Math.max(...totals,1);

  const W=svg.clientWidth||300, H=220, mL=40, mB=26, mT=8, mR=6;
  const pw=W-mL-mR, ph=H-mT-mB;
  const bw=pw/classes.length*0.72, gap=pw/classes.length;
  const ns='http://www.w3.org/2000/svg';
  const add=(t,at)=>{const e=document.createElementNS(ns,t);for(const k in at)e.setAttribute(k,at[k]);svg.appendChild(e);return e;};

  // y gridlines + labels
  for(let i=0;i<=4;i++){
    const y=mT+ph-ph*i/4, val=ymax*i/4;
    add('line',{x1:mL,y1:y,x2:W-mR,y2:y,stroke:'#26332b','stroke-width':1});
    const t=add('text',{x:mL-5,y:y+3,'text-anchor':'end',fill:'#8ba393','font-size':10}); t.textContent=fmt(val);
  }
  // bars
  classes.forEach((c,i)=>{
    const x=mL+i*gap+(gap-bw)/2; let yacc=mT+ph;
    const order=[...top,'other'];
    order.forEach(s=>{
      const v=stacks[c][s]; if(!v) return;
      const h=ph*v/ymax; yacc-=h;
      add('rect',{x,y:yacc,width:bw,height:Math.max(h,0.5),fill:colOf(s),
                  rx:1}).append(document.createElementNS(ns,'title'));
      svg.lastChild.lastChild.textContent=`${s} ${c}-${c+2}": ${fmt(v,1)} ft²`;
    });
    const t=add('text',{x:x+bw/2,y:H-8,'text-anchor':'middle',fill:'#8ba393','font-size':9});
    t.textContent = c>=40?'40+':`${c}`;
  });
  add('text',{x:mL,y:H-8,'text-anchor':'start',fill:'#8ba393','font-size':9}).textContent='DBH (in) →';

  // swatches
  el('swatches').innerHTML = top.map(s=>
    `<span><i class="sw" style="background:${color[s]}"></i>${s} · ${spName[s]||''}</span>`).join('')
    + `<span><i class="sw" style="background:${OTHER}"></i>other</span>`;
}

/* ---------------- simulation plan builder ---------------- */
let plan = { name:'plan', cycles:10, period:10, actions:[] };
let lastSim = null;

const _opt = (items, sel) => items.map(it => { const [v,t]=Array.isArray(it)?it:[it,it];
  return `<option value="${v}"${v===sel?' selected':''}>${t}</option>`; }).join('');
const _planYear = () => (new Date().getFullYear()) + (plan.period||10);

function renderActions(){
  const host = el('actions'); if(!host) return; host.innerHTML='';
  if(!plan.actions.length){
    host.innerHTML = '<div class="act"><span class="tag">grow</span>&nbsp;No intervention — grow the stands only.</div>';
    return;
  }
  plan.actions.forEach((a,i)=>{
    const d=document.createElement('div'); d.className='act '+a.kind;
    const cy = a=>`<label>at</label><input type="number" class="num" min="0" max="${plan.cycles-1}" value="${a.cycle}" data-i="${i}" data-f="cycle" title="cycles from now (0 = immediately)"><span style="color:var(--muted);font-size:11px">${a.cycle===0?'now':'+'+a.cycle*plan.period+'yr'}</span>`;
    if(a.kind==='thin'){
      d.innerHTML = `<span class="tag thin">Thin</span>${cy(a)}
        <label>sp</label><input list="splist" value="${a.species||'all'}" data-i="${i}" data-f="species" style="width:80px" title="species to cut: 'all', or one/more FVS alpha codes (comma-separated)">
        <label>to</label><select data-i="${i}" data-f="metric">${_opt(['BA','TPA','SDI'],a.metric)}</select>
        <input type="number" class="num" value="${a.target}" data-i="${i}" data-f="target" title="residual (0 = clearcut)">
        <select data-i="${i}" data-f="direction">${_opt([['below','from below'],['above','from above']],a.direction)}</select>
        <label>DBH</label><input type="number" class="num" value="${a.dbh_lo}" data-i="${i}" data-f="dbh_lo">–<input type="number" class="num" value="${a.dbh_hi}" data-i="${i}" data-f="dbh_hi">″
        <button class="rm" data-rm="${i}" title="remove">×</button>`;
    } else {
      d.innerHTML = `<span class="tag plant">Plant</span>${cy(a)}
        <label>sp</label><input list="splist" value="${a.species}" data-i="${i}" data-f="species" style="width:80px" title="species to plant: one/more FVS alpha codes (comma-separated)">
        <input type="number" class="num" value="${a.tpa}" data-i="${i}" data-f="tpa" title="trees/acre (each species)">tpa
        <label>surv</label><input type="number" class="num" value="${a.survival}" data-i="${i}" data-f="survival">%
        <button class="rm" data-rm="${i}" title="remove">×</button>`;
    }
    host.appendChild(d);
  });
  host.querySelectorAll('[data-f]').forEach(inp => inp.onchange = e => {
    const i=+e.target.dataset.i, f=e.target.dataset.f; let v=e.target.value;
    if(!['species','metric','direction'].includes(f)) v = +v;
    plan.actions[i][f] = v;
    if(f==='cycle') renderActions();     // refresh the "+N yr" label
  });
  host.querySelectorAll('[data-rm]').forEach(b => b.onclick = e => {
    plan.actions.splice(+e.currentTarget.dataset.rm,1); renderActions(); });
}
function addThin(){ plan.actions.push({kind:'thin',cycle:1,metric:'BA',target:80,direction:'below',species:'all',dbh_lo:0,dbh_hi:999}); renderActions(); }
function addPlant(){ plan.actions.push({kind:'plant',cycle:1,species:'DF',tpa:300,survival:85}); renderActions(); }

function savePlan(){
  plan.cycles=+el('pCycles').value; plan.period=+el('pPeriod').value;
  const blob=new Blob([JSON.stringify(plan,null,2)],{type:'application/json'});
  const url=URL.createObjectURL(blob), a=document.createElement('a');
  a.href=url; a.download=`plan-${new Date().toISOString().slice(0,10)}.json`;
  document.body.appendChild(a); a.click(); a.remove(); setTimeout(()=>URL.revokeObjectURL(url),0);
}
async function loadPlan(e){
  const f=e.target.files[0]; if(!f) return;
  try{
    const p=JSON.parse(await f.text());
    plan={name:p.name||'plan',cycles:+p.cycles||10,period:+p.period||10,actions:Array.isArray(p.actions)?p.actions:[]};
    el('pCycles').value=plan.cycles; el('pPeriod').value=plan.period; renderActions();
    setStatus('plan loaded: '+(plan.name||'')+' ('+plan.actions.length+' actions)');
  }catch(err){ setStatus('plan load error: '+err.message); }
  e.target.value='';
}

async function runSim(){
  if(!currentAoiGeom){ setStatus('select an AOI first'); return; }
  plan.cycles=+el('pCycles').value; plan.period=+el('pPeriod').value;
  setStatus('running FVS simulation…');
  el('runSim').disabled=true;
  el('simprog').style.display='block'; el('simprogbar').style.width='4%'; el('simprogtxt').textContent='starting…';
  const poll=setInterval(async()=>{
    try{
      const p=await fetch('/api/simprogress').then(r=>r.json());
      if(p && p.total>0){
        const pct=Math.min(100,Math.max(4,Math.round(100*p.done/p.total)));
        el('simprogbar').style.width=pct+'%';
        el('simprogtxt').textContent = p.done===0 ? `${p.total} plots — first run is compiling FVS…`
                                                   : `${p.done} / ${p.total} plots`;
      }
    }catch(_){}
  }, 400);
  try{
    const res=await fetch('/api/simulate',{method:'POST',headers:{'Content-Type':'application/json'},
      body:JSON.stringify({geometry:currentAoiGeom, plan})});
    const d=await res.json();
    if(!res.ok||d.error){ setStatus('simulate error: '+(d.error||res.status)); return; }
    lastSim=d; plan.period=d.period||plan.period; drawSimChart(); setupSlider(d);
    setStatus(`simulated ${fmt(d.nplots_sim)} of ${fmt(d.nplots_aoi)} plots · ${fmt(d.acres)} ac · ${d.ncycles} cycles${d.capped?' (largest '+d.nplots_sim+' run)':''}`);
  }catch(err){ setStatus('simulate error: '+err.message); }
  finally{ clearInterval(poll); el('runSim').disabled=false; el('simprog').style.display='none'; }
}

/* time slider — recolor the AOI's pixels to the projected state at each cycle */
let simCycle=0, simTimer=null;
const SIM_MAP_META = { carbon:{cmap:'greens',units:'t C/ha'}, ba:{cmap:'viridis',units:'ft²/ac'}, volume:{cmap:'viridis',units:'ft³/ac'} };
// the metric to actually draw: the chosen one if the sim produced data for it, else the
// first metric that DID get a domain (e.g. carbon is empty for unported-FFE variants, so a
// carbon-default raster would be silently all-transparent — fall back to ba/volume).
function effSimMetric(){
  const sel = el('simMapMetric').value, doms = (lastSim && lastSim.domains) || {};
  if(doms[sel]) return sel;
  const first = Object.keys(doms)[0];
  return first || sel;
}
// URL of the pre-rendered, cached AOI image for (metric, cycle); token-versioned so the
// browser caches each immutably and the slider just swaps a ready image.
function simImgUrl(metric, cycle){
  const im = lastSim && lastSim.image; if(!im) return null;
  return `/simimage/${im.token}/${metric}/${cycle}.png`;
}
function renderSimLegend(){
  if(!lastSim || !lastSim.domains) return;
  const m=effSimMetric(), meta=SIM_MAP_META[m]||SIM_MAP_META.carbon, dom=lastSim.domains[m]||[0,1];
  el('simLegbar').style.background = CMAP_CSS[meta.cmap] || CMAP_CSS.viridis;
  el('simLegLo').textContent=fmt(dom[0]); el('simLegHi').textContent=fmt(dom[1]);
  el('simLegLabel').textContent = meta.units + (m!==el('simMapMetric').value ? ' •' : '');
}
function setupSlider(d){
  const sl=el('simSlider'); sl.max=Math.max(0,(d.ncycles||1)-1); sl.value=0; simCycle=0;
  el('timerow').style.display = (d.ncycles>1) ? 'flex' : 'none';
  el('simLayerGrp').hidden = !d.image; el('showSim').checked = true;
  renderSimLegend();
  updateSimTime(); setSimLayer(); drawSimChart(); updateLeftForCycle();
}
// drive the LEFT panel (species×DBH distribution + summary cards) to the projected cycle
function updateLeftForCycle(){
  if(!lastSim || !lastSim.dist) return;
  const c=lastSim.cycles[simCycle]||{}, cells=lastSim.dist[simCycle]||[];
  drawChart(cells);
  const when = simCycle===0 ? 'now' : '+'+simCycle*(lastSim.period||plan.period)+'yr';
  if(el('cycleBadge')) el('cycleBadge').textContent = 'projected · '+when;
  const cards=[
    ['When', when, 'cyc '+simCycle],
    ['Basal area', fmt(c.ba_ac,1), 'ft²/ac'],
    ['Live carbon', fmt(c.carbon_total), 'tons AG'],
    ['Live volume', fmt(c.volume), 'ft³'],
    ['Removed vol', fmt(c.removed_volume), 'ft³'],
    ['Plots', fmt(lastSim.nplots_sim), 'simulated'],
  ];
  el('cards').innerHTML = cards.map(([k,v,u])=>
    `<div class="card"><div class="k">${k}</div><div class="v">${v}</div><div class="u">${u}</div></div>`).join('');
  el('statsGrp').hidden=false; el('chartGrp').hidden=false;
}
function updateSimTime(){
  const t = simCycle===0 ? 'now' : `+${simCycle*(plan.period||10)} yr`;
  el('simTime').textContent = simCycle===0 ? 'now (cycle 0)' : `+${simCycle*(plan.period||10)} yr · cyc ${simCycle}`;
  if(el('simLayerWhen')) el('simLayerWhen').textContent = t;
}
function setSimLayer(){
  if(!map || !lastSim || !lastSim.image) return;
  const im = lastSim.image;
  if(!el('showSim').checked){                        // dedicated toggle, independent of the whole-map layer
    if(map.getLayer('simlayer')) map.removeLayer('simlayer');
    if(map.getSource('simlayer')) map.removeSource('simlayer');
    return;
  }
  const metric = effSimMetric();                     // fall back off an empty metric so it's never blank
  const op = (+el('simOp').value)/100 || 0.92;
  const url = simImgUrl(metric, simCycle);
  const src = map.getSource('simlayer');
  if(src && src.updateImage){
    // cycle/metric changed → point the SAME image source at the cached image for this cycle
    src.updateImage({url, coordinates: im.corners});
    map.setPaintProperty('simlayer','raster-opacity',op);
  } else {
    if(map.getLayer('simlayer')) map.removeLayer('simlayer');
    if(map.getSource('simlayer')) map.removeSource('simlayer');
    // a single georeferenced image over the AOI (pre-rendered + cached server-side) —
    // ABOVE the translucent AOI fill (clearly visible), BELOW the AOI outline
    map.addSource('simlayer', {type:'image', url, coordinates: im.corners});
    const before = map.getLayer('aoi-line') ? 'aoi-line' : undefined;
    map.addLayer({id:'simlayer', type:'raster', source:'simlayer',
      paint:{'raster-opacity':op, 'raster-resampling':'nearest'}}, before);
  }
}
function clearSimLayer(){
  if(simTimer){ clearInterval(simTimer); simTimer=null; el('simPlay').textContent='▶'; }
  el('timerow').style.display='none'; el('simLayerGrp').hidden=true;
  if(map){ if(map.getLayer('simlayer')) map.removeLayer('simlayer'); if(map.getSource('simlayer')) map.removeSource('simlayer'); }
}
function onSlider(v){ simCycle=+v; el('simSlider').value=simCycle; updateSimTime(); setSimLayer(); drawSimChart(); updateLeftForCycle(); }
function toggleSimPlay(){
  if(simTimer){ clearInterval(simTimer); simTimer=null; el('simPlay').textContent='▶'; return; }
  el('simPlay').textContent='❚❚';
  const mx=+el('simSlider').max;
  simTimer=setInterval(()=>{ onSlider(simCycle>=mx ? 0 : simCycle+1); }, 850);
}

function drawSimChart(){
  const svg=el('simChart'); if(!svg) return; svg.innerHTML='';
  el('simEmpty').style.display = (lastSim && lastSim.cycles && lastSim.cycles.length) ? 'none' : 'block';
  if(!lastSim || !lastSim.cycles || !lastSim.cycles.length) return;
  const metric=el('simMetric').value;
  const pts=lastSim.cycles.map(c=>({x:c.elapsed, y:+c[metric]||0}));
  const W=svg.clientWidth||360, H=150, mL=52, mB=22, mT=8, mR=10, pw=W-mL-mR, ph=H-mT-mB;
  const xmax=Math.max(...pts.map(p=>p.x),1), ymax=Math.max(...pts.map(p=>p.y),1)*1.08;
  const X=x=>mL+pw*x/xmax, Y=y=>mT+ph-ph*y/ymax;
  const ns='http://www.w3.org/2000/svg', add=(t,a)=>{const e=document.createElementNS(ns,t);for(const k in a)e.setAttribute(k,a[k]);svg.appendChild(e);return e;};
  for(let i=0;i<=3;i++){ const y=mT+ph-ph*i/3;
    add('line',{x1:mL,y1:y,x2:W-mR,y2:y,stroke:'#26332b','stroke-width':1});
    add('text',{x:mL-5,y:y+3,'text-anchor':'end',fill:'#8ba393','font-size':9}).textContent=fmt(ymax*i/3); }
  const isRem = metric==='removed_volume';
  const path='M'+pts.map(p=>X(p.x)+','+Y(p.y)).join(' L');
  if(isRem){ pts.forEach(p=>{ if(p.y>0) add('rect',{x:X(p.x)-4,y:Y(p.y),width:8,height:mT+ph-Y(p.y),fill:'#f59e0b',rx:1}); }); }
  else { add('path',{d:path,fill:'none',stroke:'#4ade80','stroke-width':2});
         pts.forEach(p=>add('circle',{cx:X(p.x),cy:Y(p.y),r:2.5,fill:'#4ade80'})); }
  const step=Math.ceil(pts.length/6);
  pts.forEach((p,i)=>{ if(i%step===0||i===pts.length-1) add('text',{x:X(p.x),y:H-7,'text-anchor':'middle',fill:'#8ba393','font-size':9}).textContent='+'+p.x; });
  // current-cycle marker (tracks the time slider)
  if(el('timerow').style.display!=='none'){ const mx=X(simCycle*(lastSim.period||plan.period));
    add('line',{x1:mx,y1:mT,x2:mx,y2:mT+ph,stroke:'#facc15','stroke-width':1.5,'stroke-dasharray':'3 2'}); }
}

function onAoiReady(cells){
  el('runSim').disabled = false;
  // seed the plant species picker with the AOI's own species + common western codes
  const dl=el('splist'); if(!dl) return;
  const seen=new Set(), out=[];
  ['DF','PP','LP','ES','AF','WL','GF','WP','WH','RC','WF','JU','PI','AS','LM'].forEach(s=>{seen.add(s);out.push(s);});
  (cells||[]).forEach(c=>{ const s=(c.symbol||'').slice(0,2).toUpperCase(); if(s&&!seen.has(s)){seen.add(s);out.push(s);} });
  dl.innerHTML = out.map(s=>`<option value="${s}">`).join('');
}

boot().catch(err => { console.error(err); try { status.textContent = 'startup error: ' + err.message; } catch(_){} });
