'use strict';

const CMAP_CSS = {
  greens:  'linear-gradient(90deg,#f7fcf5,#c7e9c0,#74c476,#238b45,#00441b)',
  viridis: 'linear-gradient(90deg,#440154,#414487,#2a788e,#22a884,#7ad151,#fde725)',
  magma:   'linear-gradient(90deg,#000004,#3b0f70,#8c2981,#de4968,#fe9f6d,#fcfdbf)',
};
// species palette for the chart (top species get a color; rest -> "other")
const SP_PALETTE = ['#4ade80','#38bdf8','#f59e0b','#f472b6','#a78bfa',
                    '#facc15','#fb7185','#2dd4bf','#94a3b8'];

const el = id => document.getElementById(id);
const status = el('status');
const fmt = (n, d=0) => Number(n).toLocaleString(undefined,{maximumFractionDigits:d});

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

  map = new maplibregl.Map({
    container: 'map',
    style: {
      version: 8,
      sources: {
        base: {
          type: 'raster', tileSize: 256,
          tiles: ['https://a.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png',
                  'https://b.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png',
                  'https://c.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png'],
          attribution: '© OpenStreetMap © CARTO · USFS TreeMap 2022'
        },
        [RASTER_SRC]: { type:'raster', tileSize:256, tiles:[tileUrl(curAttr)] }
      },
      layers: [
        { id:'base', type:'raster', source:'base' },
        { id:'treemap', type:'raster', source:RASTER_SRC, paint:{'raster-opacity':0.85} }
      ]
    },
    center: [(meta.bounds.west+meta.bounds.east)/2, (meta.bounds.south+meta.bounds.north)/2],
    zoom: 3.4
  });
  map.addControl(new maplibregl.NavigationControl({showCompass:false}), 'top-right');
  map.on('load', () => {
    map.fitBounds([[meta.bounds.west,meta.bounds.south],[meta.bounds.east,meta.bounds.north]],
                  {padding:20, duration:0});
    addAoiLayers();
    updateLegend();
    setStatus(`${fmt(meta.nplots_total)} plots · ${attrs[curAttr].label}`);
  });

  sel.onchange = e => { curAttr = e.target.value; swapTiles(); updateLegend();
                        setStatus(attrs[curAttr].label); };
  el('op').oninput = e => {
    const v = +e.target.value; el('opv').textContent = v+'%';
    if (map.getLayer('treemap')) map.setPaintProperty('treemap','raster-opacity', v/100);
  };
  el('draw').onclick = startDraw;
  el('clear').onclick = clearAoi;
  el('upload').onclick = () => el('file').click();
  el('file').onchange = onUpload;
  window.addEventListener('keydown', e => { if (e.key==='Escape') cancelDraw(); });
}

const tileUrl = a => `/tiles/${a}/{z}/{x}/{y}.png`;
function swapTiles(){
  if (map.getLayer('treemap')) map.removeLayer('treemap');
  if (map.getSource(RASTER_SRC)) map.removeSource(RASTER_SRC);
  map.addSource(RASTER_SRC,{type:'raster',tileSize:256,tiles:[tileUrl(curAttr)]});
  map.addLayer({id:'treemap',type:'raster',source:RASTER_SRC,
    paint:{'raster-opacity':(+el('op').value)/100}}, aoiBeforeId());
}
function updateLegend(){
  const a = attrs[curAttr];
  el('legbar').style.background = CMAP_CSS[a.cmap] || CMAP_CSS.viridis;
  el('legLo').textContent = fmt(a.lo);
  el('legHi').textContent = fmt(a.hi);
  el('legLabel').textContent = a.units || '';
}
function setStatus(t){ status.textContent = t; }

/* ---------------- AOI layers ---------------- */
function addAoiLayers(){
  map.addSource('aoi',{type:'geojson',data:emptyFC()});
  map.addLayer({id:'aoi-fill',type:'fill',source:'aoi',
    paint:{'fill-color':'#4ade80','fill-opacity':0.12}});
  map.addLayer({id:'aoi-line',type:'line',source:'aoi',
    paint:{'line-color':'#4ade80','line-width':2}});
  map.addSource('draft',{type:'geojson',data:emptyFC()});
  map.addLayer({id:'draft-line',type:'line',source:'draft',
    paint:{'line-color':'#facc15','line-width':2,'line-dasharray':[2,1]}});
  map.addLayer({id:'draft-pt',type:'circle',source:'draft',
    filter:['==','$type','Point'],
    paint:{'circle-radius':4,'circle-color':'#facc15'}});
}
const aoiBeforeId = () => map.getLayer('aoi-fill') ? 'aoi-fill' : undefined;
const emptyFC = () => ({type:'FeatureCollection',features:[]});

/* ---------------- polygon draw ---------------- */
let drawing=false, verts=[];
function startDraw(){
  clearAoi(); drawing=true; verts=[];
  map.getCanvas().style.cursor='crosshair';
  el('drawhint').style.display='block';
  map.doubleClickZoom.disable();
  map.on('click', onDrawClick);
  map.on('dblclick', onDrawDone);
  map.on('mousemove', onDrawMove);
}
function onDrawClick(e){ verts.push([e.lngLat.lng,e.lngLat.lat]); renderDraft(); }
function onDrawMove(e){ if(drawing && verts.length) renderDraft([e.lngLat.lng,e.lngLat.lat]); }
function renderDraft(hover){
  const pts = hover ? verts.concat([hover]) : verts;
  const feats = verts.map(v=>({type:'Feature',geometry:{type:'Point',coordinates:v}}));
  if (pts.length>=2) feats.push({type:'Feature',geometry:{type:'LineString',coordinates:pts}});
  map.getSource('draft').setData({type:'FeatureCollection',features:feats});
}
function onDrawDone(e){
  if (e) e.preventDefault();
  if (verts.length < 3){ cancelDraw(); return; }
  const ring = verts.concat([verts[0]]);
  const geometry = {type:'Polygon',coordinates:[ring]};
  cancelDraw();
  map.getSource('aoi').setData({type:'Feature',geometry});
  submitAoi({body:JSON.stringify({geometry}),
             headers:{'Content-Type':'application/json'}});
}
function cancelDraw(){
  drawing=false; verts=[];
  map.getCanvas().style.cursor='';
  el('drawhint').style.display='none';
  map.off('click',onDrawClick); map.off('dblclick',onDrawDone); map.off('mousemove',onDrawMove);
  setTimeout(()=>map.doubleClickZoom.enable(),0);
  if (map.getSource('draft')) map.getSource('draft').setData(emptyFC());
}
function clearAoi(){
  cancelDraw();
  if (map.getSource('aoi')) map.getSource('aoi').setData(emptyFC());
  el('statsGrp').hidden = true; el('chartGrp').hidden = true;
}

/* ---------------- upload ---------------- */
async function onUpload(e){
  const f = e.target.files[0]; if(!f) return;
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
  drawChart(d.cells);
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

boot();
