#!/usr/bin/env python3
# Extract WC (West Cascades) per-species coefficient data from the Fortran sources
# and emit data/westcascades/{species_coefficients.csv, ecocls.csv, pcoml.csv}.
import re, os, sys

OUT = os.environ.get("WC_OUT",
        os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "data", "westcascades")))
WC  = os.environ.get("WC_SRC", "/workspace/ForestVegetationSimulator/wc")
os.makedirs(OUT, exist_ok=True)

def expand(spec):
    """Expand a Fortran DATA repeat list like ['2*5.308','5.313'] -> floats."""
    out = []
    for tok in spec:
        tok = tok.strip()
        if not tok: continue
        if '*' in tok:
            n, v = tok.split('*')
            out += [float(v)] * int(n)
        else:
            out.append(float(tok))
    return out

N = 39
alpha = ["SF","WF","GF","AF","RF","__","NF","YC","IC","ES","LP","JP","SP","WP","PP","DF","RW","RC","WH","MH",
         "BM","RA","WA","PB","GC","AS","CW","WO","WJ","LL","WB","KP","PY","DG","HT","CH","WI","__","OT"]
fia   = ["011","015","017","019","020","999","022","042","081","093","108","116","117","119","122","202","211","242","263","264",
         "312","351","352","375","431","746","747","815","064","072","101","103","231","492","500","768","920","998","999"]
plants= ["ABAM","ABCO","ABGR","ABLA","ABMA","","ABPR","CANO9","CADE27","PIEN","PICO","PIJE","PILA","PIMO3","PIPO","PSME","SESE3","THPL","TSHE","TSME",
         "ACMA3","ALRU2","ALRH2","BEPA","CHCHC4","POTR5","POBAT","QUGA4","JUOC","LALY","PIAL","PIAT","TABR2","CONU4","CRATA","PREM","SALIX","","2TREE"]

# --- bratio.f: JBARK[39] -> BARKB[.,14] (fia_ref, b2, b3, eqtype) ---
JBARK = [2,2,2,2,2,2,2,5,5,10, 10,4,4,4,3,1,14,11,12,11, 6,9,9,6,7,9,9,8,11,10, 13,13,13,9,9,9,9,1,10]
BARKB = {  # index -> (b2, b3, eqtype)   eqtype 1=power a*D^b ; 2=linear a+bD
 1:(0.903563,0.989388,1), 2:(0.904973,1.0,1), 3:(0.809427,1.016866,1), 4:(0.859045,1.0,1),
 5:(0.837291,1.0,1), 6:(0.08360,0.94782,2), 7:(0.15565,0.90182,2), 8:(0.8558,1.0213,1),
 9:(0.075256,0.949670,2), 10:(0.9,1.0,1), 11:(0.949670,1.0,1), 12:(0.933710,1.0,1),
 13:(0.933290,1.0,1), 14:(0.7012,1.04862,1) }
bark1 = [BARKB[JBARK[i]][0] for i in range(N)]
bark2 = [BARKB[JBARK[i]][1] for i in range(N)]
bark_imap = [BARKB[JBARK[i]][2] for i in range(N)]   # true Fortran eqtype (1=power,2=linear)

# --- blkdat.f SIGMAR (dg residual SD) ---
sigmar = expand("0.5450 2*0.4390 0.3960 2*0.3102 0.4275 0.3931 2*0.4842 0.3690 0.3222 2*0.5494 0.3222 0.4456 0.6178 0.4442 0.4104 0.3751 0.5107 0.7487 5*0.5357 0.236 0.5357 4*0.4842 6*0.5357".split())
# --- blkdat.f HT1/HT2 (wykoff htcalc height-dub coeffs, >=5" dbh) ---
ht1 = expand("5.288 2*5.308 3*5.313 5.327 5.143 2*5.188 4.865 5.333 2*5.382 5.333 5.288 5.3401 5.271 5.298 5.081 4.700 4.886 7*5.152 4*5.188 6*5.152".split())
ht2 = expand("-14.147 2*-13.624 3*-15.321 -15.450 -13.497 2*-13.801 -9.305 -17.762 2*-15.866 -17.762 -14.147 -15.9354 -14.996 -13.240 -13.430 -6.326 -8.792 7*-13.576 4*-13.801 6*-13.576".split())
# --- sichg.f A/B/REFLOC/REFAGE ---
sichg_a = expand("21.35000 9.01840 9.01840 33.72545 14.81367 14.81367 7.938059 11.56252 8.00000 33.72545 10.65724 8.000000 21.02281 21.02281 8.00000 11.56252 11.56252 11.56252 6.15767 45.24969 11.56252 3.28241 5*11.56252 4.94166 11.56252 8.11668 9*11.56252".split())
sichg_b = expand("-0.10290 -0.05700 -0.05700 -0.27451 -0.11174 -0.11174 -0.02873 -0.05586 -0.04286 -0.27451 -0.10667 -0.04259 -0.15978 -0.15978 -0.04286 -0.05586 -0.05586 -0.05586 -0.03596 -1.28885 -0.05586 -0.02987 5*-0.05586 -0.02419 -0.05586 -0.05661 9*-0.05586".split())
sichg_refage = expand("100 50 50 100 50 50 4*100 50 7*100 50 100 100 20 5*100 50 100 50 9*100".split())
refloc = ['B']*10 + ['T'] + ['B']*10 + ['T'] + ['B']*17    # 10B,T,10B,T,17B
sichg_refloc = [1 if c=='T' else 0 for c in refloc]        # 1=T (total age), 0=B (breast-height)
# --- crown.f IMAP (species -> crown group) ---
crown_imap = [1,2,2,3,3,3,4,15,11,11,16,6,5,5,6,7,11,8,9,10,12,13,14,14,14,14,14,14,14,11,11,11,11,14,14,14,14,14,14]
# --- sitset.f site reduction factors (misc hardwoods); sp20(MH)&sp28(WO) special (in code) ---
site_redux = [1.0]*N
for sp,f in {21:0.75,23:0.65,24:1.5,25:0.70,26:0.75,27:0.85,29:0.23,31:0.70,33:0.25,34:0.60,35:0.25,36:0.50,37:0.50}.items():
    site_redux[sp-1] = f

for name,arr in [("sigmar",sigmar),("ht1",ht1),("ht2",ht2),("sichg_a",sichg_a),("sichg_b",sichg_b),
                 ("sichg_refage",sichg_refage),("crown_imap",crown_imap),("bark1",bark1)]:
    assert len(arr)==N, f"{name} len {len(arr)}!={N}"

hdr = ["species_index","code_alpha","code_fia","code_plants","bark1","bark2","bark_imap",
       "dg_resid_sd","ht1","ht2","sichg_a","sichg_b","sichg_refage","sichg_refloc",
       "site_redux","crown_imap"]
def fmt(x):
    if isinstance(x,int): return str(x)
    return ("%g"%x)
with open(os.path.join(OUT,"species_coefficients.csv"),"w") as f:
    f.write(",".join(hdr)+"\n")
    for i in range(N):
        row = [i+1, alpha[i], fia[i], plants[i], fmt(bark1[i]), fmt(bark2[i]), int(bark_imap[i]),
               fmt(sigmar[i]), fmt(ht1[i]), fmt(ht2[i]), fmt(sichg_a[i]), fmt(sichg_b[i]),
               int(sichg_refage[i]), int(sichg_refloc[i]), fmt(site_redux[i]), int(crown_imap[i])]
        f.write(",".join(str(c) for c in row)+"\n")
print("wrote species_coefficients.csv")

# --- ecocls.f: parse the 139 PA rows: 'PA      ','SCIEN ', SDI,'SPC ', SITE, NUM, FLAG, SEQ, ---
txt = open(os.path.join(WC,"ecocls.f")).read()
# data rows look like: &'CFS551  ','PSME  ', 815.,'DF  ',  73.,   1,   1,  16,
rowre = re.compile(r"'([A-Z0-9 ]{8})'\s*,\s*'[A-Z0-9 ]{6}'\s*,\s*([0-9.]+)\s*,\s*'([A-Z0-9 ]{4})'\s*,\s*([0-9.]+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*[,/]")
ecos=[]
for m in rowre.finditer(txt):
    pa=m.group(1).strip(); sdi=float(m.group(2)); spc=m.group(3).strip()
    site=float(m.group(4)); num=int(m.group(5)); flag=int(m.group(6)); seq=int(m.group(7))
    ecos.append((pa,spc,seq,sdi,site,num,flag))
with open(os.path.join(OUT,"ecocls.csv"),"w") as f:
    f.write("pa,spc,fvsseq,sdimx,site,numbr,iflag\n")
    for r in ecos:
        f.write("%s,%s,%d,%g,%g,%d,%d\n"%r)
print("wrote ecocls.csv rows=",len(ecos))

# --- habtyp.f: PCOML(139) ---
htxt = open(os.path.join(WC,"habtyp.f")).read()
# grab the DATA blocks for PCOML
pcoml=[]
for m in re.finditer(r"'([A-Z0-9]{2,8})\s*'", htxt):
    pass
# simpler: find all quoted 8-char codes inside the PCOML DATA region (between first DATA (PCOML and RETURN)
region = htxt[htxt.index("DATA (PCOML"):htxt.index("LPVREF=")]
for m in re.finditer(r"'([A-Z0-9]{6}) *'", region):
    code=m.group(1).strip()
    pcoml.append(code)
assert len(pcoml)==139, f"pcoml len {len(pcoml)}"
with open(os.path.join(OUT,"pcoml.csv"),"w") as f:
    f.write("index,pa\n")
    for i,c in enumerate(pcoml): f.write("%d,%s\n"%(i+1,c))
print("wrote pcoml.csv rows=",len(pcoml), "PCOML[52]=",pcoml[51])
print("CFS551 eco:", [r for r in ecos if r[0]=="CFS551"])
