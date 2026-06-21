#!/usr/bin/env node
// gen_callgraph.js — regenerate docs/decision_flow.html from the faithful oracle.
// Indexes every top-level `function NAME(` / one-liner def in /workspace/FVSjulia/src,
// extracts the FVS subroutine calls in each body (keeping only edges that resolve to
// a known routine, so math/intrinsics drop out as atomic leaves), tags scope from the
// source directory, and embeds the JSON graph into the HTML viewer template.
const fs=require('fs'), path=require('path');
const ORACLE=process.argv[2]||'/workspace/FVSjulia/src';
const HERE=__dirname, OUT=path.join(HERE,'..','docs','decision_flow.html');
const STOP=new Set(('DBCHK DBCHK_FVS DBSCAN GETLUN FVSGETRTNCODE FVSRESTART FVSSTOPPOINT FVSGETRESTARTCODE GETAMSTOPPING CLEARRESTARTCODE CLEARRTNCODE FVSGETRTNCD IF FLOAT MAX MIN ABS INT MOD SQRT EXP LOG LOG10 SIN COS TAN ATAN ASIN ACOS NINT REAL DBLE SIGN AMAX1 AMIN1 MAX0 MIN0 IFIX ANINT AINT SUM PROD ZEROS ONES FILL REF VIEW COPY PUSH TRIM RPAD LPAD TRUNC FLOOR CEIL ROUND LENGTH SPRINT STRING PRINTF FORMAT WRITE READ PI MAXSP MAXTRE MAXTP1 ICYC IY JOSTND ERROR ANY ALL MAP COLLECT SORT SORTPERM FINDFIRST FINDALL FINDNEXT ISNAN ISEMPTY PARSE TRYPARSE TYPEOF CONVERT UNSAFE_TRUNC SIGNBIT HYPOT CBRT EXPM1 LOG1P TANH SINH COSH').split(' '));
let files=[]; (function walk(d){for(const e of fs.readdirSync(d,{withFileTypes:true})){const p=path.join(d,e.name); if(e.isDirectory())walk(p); else if(e.name.endsWith('.jl'))files.push(p);}})(ORACLE);
const scopeOf=f=>f.includes('/sn/')?'sn':f.includes('/extensions/fire')?'fire':f.includes('/extensions/econ')?'econ':f.includes('/extensions')?'ext':f.includes('/common/')?'common':'base';
const norm=n=>n.replace(/!$/,'').toUpperCase().replace(/[^A-Z0-9_]/g,'_');
const DEFS={};
for(const f of files){ const scope=scopeOf(f); const lines=fs.readFileSync(f,'utf8').split('\n');
  for(let i=0;i<lines.length;i++){ let mm=lines[i].match(/^function\s+([A-Za-z_][\w!]*)/);
    if(mm){ const id=norm(mm[1]); let body=[],j=i+1; for(;j<lines.length;j++){ if(/^end\b/.test(lines[j]))break; body.push(lines[j]); }
      if(!DEFS[id])DEFS[id]={label:mm[1],scope,body}; else DEFS[id].body=DEFS[id].body.concat(body); i=j; continue; }
    let one=lines[i].match(/^([A-Za-z_][\w!]*)\s*\([^)]*\)\s*=(?!=)/);
    if(one){ const id=norm(one[1]); if(!DEFS[id])DEFS[id]={label:one[1],scope,body:[lines[i].split('=').slice(1).join('=')]}; } } }
const known=new Set(Object.keys(DEFS)); const NODES={};
for(const id in DEFS){ const d=DEFS[id]; const code=d.body.map(l=>l.replace(/#.*$/,'')).join('\n').replace(/"[^"]*"/g,'""');
  const calls=[]; const re=/\b([A-Za-z_][\w]*)\s*\(/g; let m;
  for(;(m=re.exec(code));){ const c=norm(m[1]); if(c===id||STOP.has(c)||!known.has(c))continue; if(!calls.includes(c))calls.push(c); }
  NODES[id]={id,label:d.label,scope:d.scope,calls}; }
const tpl=fs.readFileSync(path.join(HERE,'decision_flow.template.html'),'utf8');
fs.writeFileSync(OUT, tpl.replace('/*__GRAPH__*/','const GRAPH = '+JSON.stringify(NODES)+';'));
console.log('wrote',OUT,'-',Object.keys(NODES).length,'routines');
