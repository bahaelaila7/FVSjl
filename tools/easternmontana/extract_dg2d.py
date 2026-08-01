#!/usr/bin/env python3
"""Extract EM 2D large-tree DG arrays from em/dgf.f (chunk-3 DG data).

- DGHAB(8,19)  habitat-class intercepts; MAPHAB(117,19) habitat-type→class index.
- DGFOR(6,19)  forest/location intercepts; MAPLOC(7,19) location-class index.
- DGDS(4,19)   diameter-sq-group intercepts; MAPDSQ(7,19) dsq-group index.
- OBSERV(6,19) DGSCOR calibration observation counts.
All col-major DATA ((ARR(L,K),L=1,n),K=...) (K=species), C-comment separated; OBSERV plain DATA.
Emits data/easternmontana/dg_2d_<name>.csv (rows=species 1..19, cols=1..n). Count-validated.
"""
import re, csv
txt = open("/workspace/ForestVegetationSimulator/em/dgf.f").read()
def nums(body):
    body = re.sub(r"C[^\n]*", "", body); body = re.sub(r"!.*", "", body); out = []
    for t in re.findall(r"-?\d*\.?\d+\*-?\d*\.?\d+|-?\d*\.?\d+", body):
        if "*" in t: c, v = t.split("*"); out += [float(v)]*int(float(c))
        else: out.append(float(t))
    return out
def grab2d(name, n):
    vals = []
    for m in re.finditer(r"DATA\s*\(\(\s*"+name+r"\(L,K\),L=1,\d+\),K=[^)]*\)\s*/(.*?)/", txt, re.S):
        vals += nums(m.group(1))
    return vals
def grab_plain(name):
    m = re.search(r"DATA\s+"+name+r"\s*/(.*?)/", txt, re.S); return nums(m.group(1))
ALPHA = ["WB","WL","DF","LM","LL","RM","LP","ES","AF","PP","GA","AS","CW","BA","PW","NC","PB","OS","OH"]
def save(name, vals, n):
    assert len(vals) == 19*n, f"{name}: {len(vals)} != {19*n}"
    out = f"/workspace/FVSjl/data/easternmontana/dg_2d_{name.lower()}.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["species_index","code_alpha"]+[f"c{j+1}" for j in range(n)])
        for k in range(19):
            # col-major: value(L,K) = vals[K*n + L]
            w.writerow([k+1, ALPHA[k]] + [repr(vals[k*n + l]) for l in range(n)])
    print(f"  {name}: {19}×{n} -> {out}")
for name, n in [("DGHAB",8),("DGFOR",6),("DGDS",4),("MAPHAB",117),("MAPLOC",7),("MAPDSQ",7)]:
    save(name, grab2d(name, n), n)
save("OBSERV", grab_plain("OBSERV"), 6)
print("all 2D DG arrays extracted + count-validated")
