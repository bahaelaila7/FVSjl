#!/usr/bin/env node
// match_coverage.js — congruence check: every FVS routine reachable from main and
// in-scope (base/sn/common, minus extensions & pure I/O) must be accounted for in
// FVSjl. FVSjl cites the Fortran it implements ("Ported from: x.f"), so a routine
// is "covered" iff FVSjl references its .f file or names it (uppercase token).
const fs=require('fs'), path=require('path');
const ORACLE='/workspace/FVSjulia/src', JL='/workspace/FVSjl/src';
const walk=d=>{let o=[];for(const e of fs.readdirSync(d,{withFileTypes:true})){const p=path.join(d,e.name);
  if(e.isDirectory())o=o.concat(walk(p)); else o.push(p);} return o;};

// ---- FVS graph: routine -> {file, scope, calls} (1:1 oracle .jl ↔ Fortran .f) ----
const STOP=new Set(('DBCHK DBCHK_FVS DBSCAN GETLUN IF FLOAT MAX MIN ABS INT MOD SQRT EXP LOG NINT REAL DBLE SIGN SUM ZEROS ONES FILL REF VIEW COPY PUSH TRIM RPAD LPAD TRUNC FLOOR CEIL ROUND LENGTH SPRINT STRING PI MAXSP MAXTRE MAXTP1 ICYC IY JOSTND ANY ALL MAP SORT FINDFIRST ISNAN ISEMPTY PARSE TRYPARSE').split(' '));
const norm=n=>n.replace(/!$/,'').toUpperCase().replace(/[^A-Z0-9_]/g,'_');
const ojl=walk(ORACLE).filter(f=>f.endsWith('.jl'));
const scopeOf=f=>f.includes('/sn/')?'sn':f.includes('/extensions/fire')?'fire':f.includes('/extensions/econ')?'econ':f.includes('/extensions')?'ext':f.includes('/common/')?'common':'base';
const DEF={};
for(const f of ojl){ const base=path.basename(f,'.jl'), scope=scopeOf(f), L=fs.readFileSync(f,'utf8').split('\n');
  for(let i=0;i<L.length;i++){ let m=L[i].match(/^function\s+([A-Za-z_][\w!]*)/); let one=L[i].match(/^([A-Za-z_][\w!]*)\s*\([^)]*\)\s*=(?!=)/);
    let name=m?m[1]:(one?one[1]:null); if(!name)continue; const id=norm(name);
    let body=[]; if(m){let j=i+1;for(;j<L.length;j++){if(/^end\b/.test(L[j]))break;body.push(L[j]);}i=j;} else body=[L[i]];
    if(!DEF[id])DEF[id]={file:base,scope,body};else DEF[id].body=DEF[id].body.concat(body); } }
const known=new Set(Object.keys(DEF)); const CALLS={};
for(const id in DEF){ const code=DEF[id].body.map(l=>l.replace(/#.*$/,'')).join('\n').replace(/"[^"]*"/g,'""');
  const cs=[]; let m; const re=/\b([A-Za-z_][\w]*)\s*\(/g;
  while(m=re.exec(code)){const c=norm(m[1]); if(c!==id&&!STOP.has(c)&&known.has(c)&&!cs.includes(c))cs.push(c);} CALLS[id]=cs; }

// ---- reachable from FVS, in-scope, excluding extensions & pure I/O/plumbing ----
const seen=new Set(['FVS']),q=['FVS']; while(q.length){const n=q.shift();for(const c of CALLS[n]||[])if(!seen.has(c)){seen.add(c);q.push(c);}}
const DROP=/^(OP[A-Z]|DBS|SV|EC|FM|MIS|PRT|GEN|GHEADS|SUMHED|SUMOUT|DISPLY|EXTREE|FILOPN|FILCLOSE|MYOPEN|UUIDGEN|VERNUM|CMDLINE|ERRGRO|LNK|RD[A-Z]|KEY|UPKEY|FNDKEY|IFWRIT|ISTLNB|UNBLNK|UPCASE|CH2NUM|TVALUE|IAPSRT|TRESOR|LB[A-Z]|_GET_|_WRITE_|PARSE_TREE|STASH|GETID|RESAGE|REVISE|GRDTIM|RCDSET|FIAHEAD|NSPREC|KEYDMP|KEYOPN|KEYRDR|RDIN|RDSTR|RDPRIN|FISHER|DUNN)/;
const scopeOK=id=>['base','sn','common'].includes((DEF[id]||{}).scope);
const inscope=[...seen].filter(id=>scopeOK(id)&&!DROP.test(id)).sort();

// ---- FVSjl references: .f files cited + uppercase tokens (in comments+code) ----
const jl=walk(JL).filter(f=>f.endsWith('.jl')).map(f=>fs.readFileSync(f,'utf8')).join('\n');
const refFiles=new Set([...jl.matchAll(/\b([a-z][a-z0-9_]+)\.f\b/g)].map(m=>m[1]));
const refToks =new Set([...jl.matchAll(/\b([A-Z][A-Z0-9_]{2,})\b/g)].map(m=>m[1]));

const covered=id=>refFiles.has(DEF[id].file)||refToks.has(id);
const miss=inscope.filter(id=>!covered(id));
const cov=inscope.filter(covered);
console.log(`in-scope reachable FVS routines: ${inscope.length}`);
console.log(`  covered (FVSjl cites .f or names it): ${cov.length}`);
console.log(`  NOT referenced anywhere in FVSjl:     ${miss.length}`);
const byscope={}; for(const id of miss)(byscope[DEF[id].scope]=byscope[DEF[id].scope]||[]).push(id+`(${DEF[id].file}.f)`);
for(const s in byscope){ console.log(`\n-- uncovered [${s}] (${byscope[s].length}) --`); console.log('  '+byscope[s].join('  ')); }
