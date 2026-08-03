# Extract ALL CI large-tree DDS coefficient arrays from ci/dgf.f → data/centralidaho/dg_*.csv.
# CI DG (ci/dgf.f) is the standard western Wykoff DDS. Arrays:
#   1-D per-species (19): DGLD DGCR DGCRSQ DGBAL DGDBAL DGBA DGLBA DGPCCF DGEL DGEL2 DGSLOP DGSLSQ
#                         DGCASP DGSASP DGCCFA DGDS(=DGDSQ, dgf.f:602) OBSERV
#   flat 2-D: DGFOR(3,19) MAPLOC(6,19) IBSERV(6,3)  (col-major grab → reshape)
#   per-column 2-D: DGHAB(12,19) [=OCURHT hab-group coeff]  ICHBCL(130,19) [CI-type→hab-group]
import sys, re; sys.path.insert(0, 'tools')
from fortran_data_extract import grab, logical_statements
SRC = "/workspace/ForestVegetationSimulator/bin/FVSci_buildDir/dgf.f"

ONED = ["DGLD","DGCR","DGCRSQ","DGBAL","DGDBAL","DGBA","DGLBA","DGPCCF","DGEL","DGEL2","DGSLOP",
        "DGSLSQ","DGCASP","DGSASP","DGCCFA","DGDS","OBSERV"]

def per_column(name, nrow, ncol):
    """Assemble a 2-D array declared as repeated `DATA (name(I,col),I=1,nrow)/.../` blocks."""
    stmts = logical_statements(SRC)
    cols = {}
    pat = re.compile(rf'\bDATA\s*\(\s*{name}\s*\(\s*I\s*,\s*(\d+)\s*\)\s*,\s*I\s*=\s*1\s*,\s*{nrow}\s*\)\s*/(.*?)/', re.I|re.S)
    for s in stmts:
        for m in pat.finditer(s):
            col = int(m.group(1)); vals = []
            for t in re.split(r'[,\s]+', m.group(2).strip()):
                if not t: continue
                mm = re.match(r"(\d+)\*(.+)", t)
                if mm: vals += [float(mm.group(2))]*int(mm.group(1))
                else:
                    try: vals.append(float(t))
                    except: pass
            cols[col] = vals
    grid = [[0.0]*ncol for _ in range(nrow)]
    for col, vals in cols.items():
        for i, v in enumerate(vals[:nrow]):
            grid[i][col-1] = v
    return grid, sorted(cols.keys())

# 1-D
c1 = {}
for n in ONED:
    v = grab(SRC, n); assert v and len(v) == 19, f"{n}: {len(v) if v else None}"; c1[n] = v
with open("data/centralidaho/dg_coeffs_1d.csv", "w") as f:
    f.write("species_index," + ",".join(ONED) + "\n")
    for i in range(19):
        f.write(str(i+1) + "," + ",".join(repr(c1[n][i]) for n in ONED) + "\n")
print("dg_coeffs_1d.csv: 19 x", len(ONED))

# flat 2-D (col-major → reshape)
def flat2d(name, nrow, ncol, fn):
    v = grab(SRC, name); assert v and len(v) == nrow*ncol, f"{name}: {len(v) if v else None} (want {nrow*ncol})"
    with open(fn, "w") as f:
        for c in range(ncol):
            f.write(",".join(repr(v[c*nrow + r]) for r in range(nrow)) + "\n")
    print(fn.split('/')[-1], f": {nrow}x{ncol}")
flat2d("DGFOR", 3, 19, "data/centralidaho/dg_dgfor.csv")
flat2d("MAPLOC", 6, 19, "data/centralidaho/dg_maploc.csv")
flat2d("IBSERV", 6, 3, "data/centralidaho/dg_ibserv.csv")

# per-column 2-D
for name, nrow, ncol, fn in [("DGHAB", 12, 19, "data/centralidaho/dg_dghab.csv"),
                             ("ICHBCL", 130, 19, "data/centralidaho/dg_ichbcl.csv")]:
    grid, present = per_column(name, nrow, ncol)
    with open(fn, "w") as f:
        for row in grid:
            f.write(",".join(repr(x) for x in row) + "\n")
    print(fn.split('/')[-1], f": {nrow}x{ncol}  columns present={len(present)}")
