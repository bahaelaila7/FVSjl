#!/usr/bin/env python3
"""Extract EM habitat-type mapping tables JTYPE(118) + NIHMAP(118) for em_habtyp.

em/habtyp.f: given input KODTYP, bucket into the sorted JTYPE list → IEMTYP (largest
index with JTYPE(IEMTYP) <= KODTYP), remap KODTYP := JTYPE(IEMTYP), then ITYPE :=
NIHMAP(IEMTYP) (one of the 30 NI habitat types, the index sitset/estab use). JTYPE is
DATA in em/blkdat.f (118 valid codes + `4*0` pad); NIHMAP is DATA in em/habtyp.f.
Emits data/easternmontana/habitat_map.csv (iemtyp, jtype_code, itype).
Validated vs live FVSem emt01 (input hab 260 → "MAPPED TO 260", ITYPE 4).
"""
import re
def grab(path, name):
    m = re.search(r"DATA\s+"+name+r"\s*/(.*?)/", open(path).read(), re.S)
    return [int(x) for x in re.findall(r"-?\d+", m.group(1))]
EM = "/workspace/ForestVegetationSimulator"
jtype  = grab(EM+"/em/blkdat.f", "JTYPE")[:118]     # drop the trailing 4*0 pad
nihmap = grab(EM+"/em/habtyp.f", "NIHMAP")[:118]
assert len(jtype) == 118 and len(nihmap) == 118, (len(jtype), len(nihmap))
assert min(nihmap) == 1 and max(nihmap) == 30, "NIHMAP must index the 30 NI types"

def map_kodtyp(kod):
    J = next((i for i, v in enumerate(jtype, 1) if v > kod), None)
    iem = (J-1) if J else 118
    if iem < 1: iem = 1
    return iem, jtype[iem-1], nihmap[iem-1]

# self-check vs live FVSem emt01 (habitat 260)
iem, code, itype = map_kodtyp(260)
assert code == 260 and itype == 4, f"hab 260 -> code {code}, ITYPE {itype} (expected 260, 4)"

out = "/workspace/FVSjl/data/easternmontana/habitat_map.csv"
with open(out, "w") as f:
    f.write("iemtyp,jtype_code,itype\n")
    for i in range(118):
        f.write(f"{i+1},{jtype[i]},{nihmap[i]}\n")
print(f"wrote 118-row habitat map -> {out}")
print(f"self-check OK: hab 260 -> mapped code {code}, ITYPE {itype} (matches live FVSem 'MAPPED TO 260')")
