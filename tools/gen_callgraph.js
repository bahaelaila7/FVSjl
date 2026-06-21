#!/usr/bin/env node
// gen_callgraph.js — build an interactive call-flow graph from a Julia source tree.
//
//   node gen_callgraph.js <srcRoot> <outHtml> <rootsCSV> "<title>" "<subtitle>"
//
// Indexes every top-level `function NAME(` / one-liner def, extracts the calls in
// each body (keeping only edges that resolve to a known routine, so math/intrinsics
// fall out as atomic leaves), the routine's leading doc/comment excerpt, and its
// scope from the directory. Embeds the JSON graph into decision_flow.template.html.
// Defaults reproduce the FVS (oracle) graph.
const fs=require('fs'), path=require('path');
const SRC   = process.argv[2] || '/workspace/FVSjulia/src';
const OUT   = process.argv[3] || path.join(__dirname,'..','docs','decision_flow.html');
const ROOTS = (process.argv[4] || 'FVS').split(',');
const TITLE = process.argv[5] || 'FVS-Southern — call graph (Fortran oracle)';
const SUBT  = process.argv[6] || 'Auto-extracted from the faithful oracle. Click a node to break it into the functions it calls (dotted = enter, solid = execution order); hover/click shows the code excerpt. Leaves with no border are atomic.';
const STOP=new Set(('DBCHK DBCHK_FVS DBSCAN GETLUN FVSGETRTNCODE FVSRESTART FVSSTOPPOINT FVSGETRESTARTCODE GETAMSTOPPING CLEARRESTARTCODE CLEARRTNCODE FVSGETRTNCD IF FLOAT MAX MIN ABS INT MOD SQRT EXP LOG LOG10 SIN COS TAN ATAN ASIN ACOS NINT REAL DBLE SIGN AMAX1 AMIN1 MAX0 MIN0 IFIX ANINT AINT SUM PROD ZEROS ONES FILL REF VIEW COPY PUSH TRIM RPAD LPAD TRUNC FLOOR CEIL ROUND LENGTH SPRINT STRING PRINTF FORMAT WRITE READ PI MAXSP MAXTRE MAXTP1 ICYC IY JOSTND ERROR ANY ALL MAP COLLECT SORT SORTPERM FINDFIRST FINDALL FINDNEXT ISNAN ISEMPTY PARSE TRYPARSE TYPEOF CONVERT UNSAFE_TRUNC SIGNBIT HYPOT CBRT EXPM1 LOG1P TANH SINH COSH').split(' '));
let files=[]; (function walk(d){for(const e of fs.readdirSync(d,{withFileTypes:true})){const p=path.join(d,e.name); if(e.isDirectory())walk(p); else if(e.name.endsWith('.jl'))files.push(p);}})(SRC);
const scopeOf=f=> (f.includes('/sn/')||f.includes('/variants/southern'))?'sn'
  : f.includes('/extensions/fire')?'fire' : f.includes('/extensions/econ')?'econ'
  : f.includes('/extensions')?'ext' : f.includes('/common/')?'common' : 'base';
const norm=n=>n.replace(/!$/,'').toUpperCase().replace(/[^A-Z0-9_]/g,'_');

// excerpt = leading Julia docstring ("""...""") and/or contiguous # comments above the def.
function excerptAbove(lines,i){
  let j=i-1, out=[];
  while(j>=0 && lines[j].trim()==='') j--;
  if(j>=0 && /"""\s*$/.test(lines[j]) && !/^\s*""".*"""/.test(lines[j])){
    let k=j-1, ds=[];
    while(k>=0 && !/^\s*"""/.test(lines[k])){ ds.unshift(lines[k]); k--; }
    out = ds; j = k-1;
    while(j>=0 && lines[j].trim()==='') j--;
  }
  let cm=[];
  while(j>=0 && /^\s*#/.test(lines[j])){ cm.unshift(lines[j].replace(/^\s*#+\s?/,'')); j--; }
  let txt=[...cm, ...out].join('\n').replace(/```[a-z]*\n?/g,'').trim();
  if(txt.length>600) txt=txt.slice(0,600).replace(/\s+\S*$/,'')+' …';
  return txt;
}

const DEF={};
for(const f of files){ const base=path.basename(f,'.jl'), scope=scopeOf(f), lines=fs.readFileSync(f,'utf8').split('\n');
  for(let i=0;i<lines.length;i++){
    let m=lines[i].match(/^function\s+([A-Za-z_][\w!]*)/);
    let one=!m && lines[i].match(/^([A-Za-z_][\w!]*)\s*\([^)]*\)\s*=(?!=)/);
    const name=m?m[1]:(one?one[1]:null); if(!name) continue;
    const id=norm(name);
    if(m){ let j=i+1, body=[]; for(;j<lines.length;j++){ if(/^end\b/.test(lines[j]))break; body.push(lines[j]); }
      if(!DEF[id]) DEF[id]={label:name,scope,file:base+'.jl',body,excerpt:excerptAbove(lines,i)};
      else DEF[id].body=DEF[id].body.concat(body);
      i=j; }
    else if(!DEF[id]) DEF[id]={label:name,scope,file:base+'.jl',body:[lines[i].split('=').slice(1).join('=')],excerpt:excerptAbove(lines,i)};
  }
}
const known=new Set(Object.keys(DEF)); const NODES={};
for(const id in DEF){ const d=DEF[id];
  const code=d.body.map(l=>l.replace(/#.*$/,'')).join('\n').replace(/"[^"]*"/g,'""');
  const calls=[]; const re=/\b([A-Za-z_][\w]*!?)\s*\(/g; let m;
  for(;(m=re.exec(code));){ const c=norm(m[1]); if(c===id||STOP.has(c)||!known.has(c))continue; if(!calls.includes(c))calls.push(c); }
  NODES[id]={id,label:d.label,scope:d.scope,file:d.file,calls,excerpt:d.excerpt};
}
const tpl=fs.readFileSync(path.join(__dirname,'decision_flow.template.html'),'utf8');
const html=tpl
  .replace('/*__GRAPH__*/','const GRAPH = '+JSON.stringify(NODES)+';')
  .replace('/*__ROOTS__*/',JSON.stringify(ROOTS)+' || ')
  .replace(/__TITLE__/g,TITLE).replace(/__SUBTITLE__/g,SUBT);
fs.writeFileSync(OUT,html);
const leaves=Object.values(NODES).filter(x=>!x.calls.length).length;
console.log('wrote',OUT,'-',Object.keys(NODES).length,'routines,',leaves,'atomic; roots',ROOTS.join(','));
