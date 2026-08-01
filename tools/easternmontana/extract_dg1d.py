#!/usr/bin/env python3
"""Extract EM 1D per-species large-tree DG coefficients from em/dgf.f (Wykoff DDS).

These are the per-species terms of WK2 = DGLD*ln(D) + DGCR*CR + DGCRSQ*CR^2 + DGBAL*BAL +
DGDBAL*BAL/ln(D+1) + DGLCCF*ln(CCF) + (DGPCC1/DGPCC2/DGPCCF PCCF terms) + slope/aspect/elev
(DGCASP*cos/DGSASP*sin*slope + DGSLOP*slope + DGSLSQ*slope^2 + DGEL*elev + DGEL2*elev^2) + the
2D habitat/forest/dsq terms (DGHAB/DGFOR/DGDS via MAPHAB/MAPLOC/MAPDSQ — extracted separately in
the chunk-3 DG port). Only the ~10 calibrated species are non-zero; others map via MAPHAB.
Emits data/easternmontana/dg_1d_coeffs.csv (species x term). Count-validated (19 each).
"""
import re, csv
txt = open("/workspace/ForestVegetationSimulator/em/dgf.f").read()
def grab(name):
    m = re.search(r"DATA\s+"+name+r"\s*/(.*?)/", txt, re.S)
    body = re.sub(r"!.*", "", m.group(1)); out = []
    for t in re.findall(r"-?[\d.]+\*-?[\d.]+|-?\.?[\d.]+", body):
        if "*" in t: n, v = t.split("*"); out += [float(v)]*int(float(n))
        else:
            try: out.append(float(t))
            except: pass
    return out
TERMS = ["DGLD","DGCR","DGCRSQ","DGBAL","DGDBAL","DGLCCF","DGPCC1","DGPCC2","DGPCCF",
         "DGCASP","DGSASP","DGSLOP","DGSLSQ","DGEL","DGEL2"]
cols = {t: grab(t) for t in TERMS}
for t, v in cols.items():
    assert len(v) == 19, f"{t} has {len(v)} (expect 19)"
ALPHA = ["WB","WL","DF","LM","LL","RM","LP","ES","AF","PP","GA","AS","CW","BA","PW","NC","PB","OS","OH"]
out = "/workspace/FVSjl/data/easternmontana/dg_1d_coeffs.csv"
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["species_index","code_alpha"] + [t.lower() for t in TERMS])
    for i in range(19):
        w.writerow([i+1, ALPHA[i]] + [repr(cols[t][i]) for t in TERMS])
print(f"wrote {out} (19 species x {len(TERMS)} DG terms), all count-validated")
print("calibrated species (DGLD!=0):", [ALPHA[i] for i in range(19) if cols["DGLD"][i] != 0])
