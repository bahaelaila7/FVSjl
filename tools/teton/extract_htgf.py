#!/usr/bin/env python3
"""Extract TT height-growth (Schreuder-Hafley SBB) coefficients from tt/htgf.f.
COF(9,33) via COF1..COF11 (each 9x3, EQUIVALENCEd COF(:,1:3),(:,4:6),...). K=(JSPC-1)*3+KEYCR.
AZBIAS/BZBIAS (MAXSP). -> data/teton/htgf_cof.csv (33 rows x 9) + htgf_zbias.csv (18 x 2)."""
import re, csv
txt = open("/workspace/ForestVegetationSimulator/tt/htgf.f").read()
def grab(n):
    m = re.search(r"DATA\s+"+n+r"\s*/(.*?)/", txt, re.S); body = re.sub(r"!.*","",m.group(1)); out=[]
    for t in re.findall(r"-?\d*\.?\d+", body):
        try: out.append(float(t))
        except: pass
    return out
# COF1..COF11 each 9x3 = 27 vals, column-major COF(coeff,crown): DATA fills COF(1,1)..COF(9,1),COF(1,2)...
rows=[]  # 33 rows (K), each 9 coeffs
for g in range(1,12):
    v = grab(f"COF{g}"); assert len(v)==27, f"COF{g}={len(v)}"
    for kc in range(3):                       # 3 crown groups
        rows.append(v[kc*9:(kc+1)*9])
assert len(rows)==33
with open("data/teton/htgf_cof.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["k"]+[f"cof{i+1}" for i in range(9)])
    for k,r in enumerate(rows): w.writerow([k+1]+[repr(x) for x in r])
az=grab("AZBIAS"); bz=grab("BZBIAS"); assert len(az)==18 and len(bz)==18
ALPHA=["WB","LM","DF","PM","BS","AS","LP","ES","AF","PP","UJ","RM","BI","MM","NC","MC","OS","OH"]
with open("data/teton/htgf_zbias.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["species_index","code_alpha","azbias","bzbias"])
    for i in range(18): w.writerow([i+1,ALPHA[i],repr(az[i]),repr(bz[i])])
print("htgf_cof.csv: 33x9 ; htgf_zbias.csv: 18x2")
print("COF row1 (WB crown1):", rows[0])
