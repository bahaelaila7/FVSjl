#!/usr/bin/env python3
"""Extract TT large-tree DG coefficient tables from tt/dgf.f (Wykoff DDS + DGCONS).

1D per-species (MAXSP=18): DGLD DGCR DGCRSQ DGBAL DGDBAL DGBA DGPCCF DGCASP DGSASP DGSLOP
  DGSLSQ DGEL DGEL2 OBSERV ISMAP  -> data/teton/dg_1d_coeffs.csv
2D (first-dim × MAXSP, Fortran column-major so flat = per-species blocks): DGFOR(5) DGSIC(5)
  DGDS(4) DGCCFC(5) IGCCFM(5) IBSERV(5) IDGSIM(5) MAPLOC(4) MAPDSQ(4) -> data/teton/dg_2d_<name>.csv
Special: DGHAB(7) -> dg_dghab.csv ; ICHBCL(363) -> dg_ichbcl.csv (habitat→DGHAB-index map for PP).
All count-validated. Faithful (programmatic). See docs/TT_VARIANT_PORT_AUDIT.md chunk 3.
"""
import re, csv, os
SRC = "/workspace/ForestVegetationSimulator/tt/dgf.f"
OUT = "/workspace/FVSjl/data/teton"
txt = open(SRC).read()
ALPHA = ["WB","LM","DF","PM","BS","AS","LP","ES","AF","PP","UJ","RM","BI","MM","NC","MC","OS","OH"]
os.makedirs(OUT, exist_ok=True)

def grab(name):
    m = re.search(r"DATA\s+"+name+r"\s*/(.*?)/", txt, re.S)
    assert m, f"DATA {name} not found"
    body = re.sub(r"!.*", "", m.group(1)); out = []
    for t in re.findall(r"-?\d*\.?\d+\*-?\d*\.?\d+|-?\d*\.?\d+", body):
        if "*" in t: n, v = t.split("*"); out += [float(v)]*int(float(n))
        else:
            try: out.append(float(t))
            except: pass
    return out

# --- 1D ---
ONE = ["DGLD","DGCR","DGCRSQ","DGBAL","DGDBAL","DGBA","DGPCCF","DGCASP","DGSASP",
       "DGSLOP","DGSLSQ","DGEL","DGEL2","OBSERV","ISMAP"]
c1 = {t: grab(t) for t in ONE}
for t, v in c1.items(): assert len(v) == 18, f"{t}={len(v)} (expect 18)"
with open(f"{OUT}/dg_1d_coeffs.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["species_index","code_alpha"] + [t.lower() for t in ONE])
    for i in range(18): w.writerow([i+1, ALPHA[i]] + [repr(c1[t][i]) for t in ONE])
print(f"dg_1d_coeffs.csv: 18 species x {len(ONE)} terms")

# --- 2D (first_dim, MAXSP); flat is per-species blocks of first_dim ---
TWO = {"DGFOR":5,"DGSIC":5,"DGDS":4,"DGCCFC":5,"IGCCFM":5,"IBSERV":5,"IDGSIM":5,"MAPLOC":4,"MAPDSQ":4}
for name, nd in TWO.items():
    v = grab(name); assert len(v) == nd*18, f"{name}={len(v)} (expect {nd*18})"
    with open(f"{OUT}/dg_2d_{name.lower()}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["species_index","code_alpha"] + [f"c{k+1}" for k in range(nd)])
        for i in range(18):
            blk = v[i*nd:(i+1)*nd]
            w.writerow([i+1, ALPHA[i]] + [repr(x) for x in blk])
    print(f"dg_2d_{name.lower()}.csv: 18 x {nd}")

# --- specials ---
dghab = grab("DGHAB"); assert len(dghab) == 7, f"DGHAB={len(dghab)}"
with open(f"{OUT}/dg_dghab.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["index","value"])
    for i, x in enumerate(dghab): w.writerow([i+1, repr(x)])
ichbcl = grab("ICHBCL"); assert len(ichbcl) == 363, f"ICHBCL={len(ichbcl)}"
with open(f"{OUT}/dg_ichbcl.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["itype","maphab_base"])
    for i, x in enumerate(ichbcl): w.writerow([i+1, int(x)])
print(f"dg_dghab.csv: 7 ; dg_ichbcl.csv: 363")
print("calibrated (DGLD!=0):", [ALPHA[i] for i in range(18) if c1["DGLD"][i] != 0])
