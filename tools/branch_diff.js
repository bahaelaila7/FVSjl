#!/usr/bin/env node
// branch_diff.js — detect missing branches/formulas by magic-number fingerprinting.
// For each FVS (oracle) routine in the C3/C4/C5 core, extract its distinctive numeric
// constants (canonicalised to 4 significant figures by VALUE, so 2.7400 == 2.74 and
// 0.005454154 == 0.0054542); a constant whose value appears NOWHERE in FVSjl
// (src/*.jl or data/*.csv) is a candidate un-ported branch or formula. Strings and
// comments are stripped so printf format specs (%10.4f) don't masquerade as constants.
const fs=require('fs'), path=require('path');
const ORACLE='/workspace/FVSjulia/src', JLSRC='/workspace/FVSjl/src', JLDATA='/workspace/FVSjl/data';
const walk=d=>{let o=[];try{for(const e of fs.readdirSync(d,{withFileTypes:true})){const p=path.join(d,e.name);
  if(e.isDirectory())o=o.concat(walk(p));else o.push(p);}}catch(e){} return o;};
const strip=t=>t.split('\n').map(l=>l.replace(/#.*$/,'')).join('\n').replace(/"(?:[^"\\]|\\.)*"/g,'""').replace(/'(?:[^'\\]|\\.)*'/g,"''");
const TRIV=new Set(['0','1','2','3','4','5','10','100','1000','0.5','0.1','0.25','0.75','0.05','0.01','0.001',
  '1.5','2.5','0.9','0.8','0.95','0.85','0.55','999','998','12','120','1.1','0.2','0.3','0.4','0.6','0.7','20','50','60','40','30']);
function consts(text){ const out=new Map();
  for(const m of strip(text).matchAll(/(?<![\w.])\d+\.\d+(?:[eE]-?\d+)?(?:f0|f)?/g)){
    let raw=m[0].replace(/f0?$/,''); let v=parseFloat(raw); if(!isFinite(v)||v===0)continue;
    let canon=String(parseFloat(v.toPrecision(4)));
    if(TRIV.has(canon))continue; out.set(canon, m[0]); }
  return out; }
let jlText=''; for(const f of [...walk(JLSRC),...walk(JLDATA)]) if(/\.(jl|csv)$/.test(f)) jlText+=fs.readFileSync(f,'utf8')+'\n';
const jlSet=new Set();
for(const m of jlText.matchAll(/(?<![\w.])\d+\.\d+(?:[eE]-?\d+)?(?:f0|f)?/g)){
  let v=parseFloat(m[0].replace(/f0?$/,'')); if(isFinite(v)&&v!==0) jlSet.add(String(parseFloat(v.toPrecision(4)))); }
const CORE={
 'C3 growth':['sn/dgf','sn/dgdriv','base/dgscor','base/autcor','sn/htgf','sn/htcalc','base/regent','sn/dgbnd','sn/crown','base/bachlo'],
 'C4 mort/density':['sn/morts','base/varmrt','base/msbmrt','base/sdical','base/sdichk','base/dense','base/mbacal','base/ccfcal','base/ptbal','base/triple','base/update'],
 'C5 volume':['base/vols','base/cfvol','base/bfvol','base/cftopk','base/bftopk','base/r8clark_vol','base/behre','base/behprm','sn/bratio','base/natcrs'],
};
for(const grp in CORE){ console.log('\n========== '+grp+' ==========');
  for(const r of CORE[grp]){ const f=ORACLE+'/'+r+'.jl'; if(!fs.existsSync(f)){console.log('  '+path.basename(r)+': (no oracle file)');continue;}
    const c=consts(fs.readFileSync(f,'utf8')); const miss=[...c].filter(([k])=>!jlSet.has(k));
    console.log(miss.length? '  '+path.basename(r)+': '+miss.length+'/'+c.size+' MISSING -> '+miss.map(([k,v])=>v).join('  ')
                           : '  '+path.basename(r)+': OK ('+c.size+' consts)'); }
}
