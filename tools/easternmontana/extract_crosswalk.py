#!/usr/bin/env python3
"""Extract the EM species crosswalk (ASPT column 10) from bin/FVSem_buildDir/spctrn.f.

The shared western ASPT(442,21) table's column legend is
  ALFA FIA PLNT AK BM CA CI CR EC EM  IE KT NC PN SO TT UT WC WS OC OP
so EM = column 10 (last field of each J=1,10 DATA-block row). Emits the KT 7-col schema
(code_alpha,code_fia,code_plants,target×4 — EM has one target, repeated) to
data/easternmontana/species_translation.csv. Faithful (programmatic, not hand-transcribed).
Self-check: 442 rows, 19 distinct targets (= MAXSP), shared conifers DF/WL/LP/ES/PP map to self.
"""
import re, collections, sys
SRC = "/workspace/ForestVegetationSimulator/bin/FVSem_buildDir/spctrn.f"
OUT = "/workspace/FVSjl/data/easternmontana/species_translation.csv"
rows, in_blockA = [], False
for ln in open(SRC).read().splitlines():
    if re.match(r"\s*DATA\s*\(\(ASPT\(I,J\),J=1,10\)", ln):
        in_blockA = True; continue
    if in_blockA:
        if ln and ln[0] in "Cc*": continue          # Fortran comment — stay in block
        s = ln.lstrip()
        if s.startswith("&"):
            f = re.findall(r"'([^']*)'", ln)
            if len(f) >= 10: rows.append((f[0].strip(), f[1].strip(), f[2].strip(), f[9].strip()))
            if ln.rstrip().endswith("/"): in_blockA = False
        elif s == "": continue
        else: in_blockA = False
with open(OUT, "w") as fo:
    fo.write("code_alpha,code_fia,code_plants,target_em,target_em2,target_em3,target_em4\n")
    for a, fia, p, em in rows: fo.write(f"{a},{fia},{p},{em},{em},{em},{em}\n")
tgt = collections.Counter(em for *_, em in rows)
print(f"{len(rows)} rows, {len(tgt)} distinct targets -> {OUT}")
assert len(rows) == 442, f"expected 442 ASPT rows, got {len(rows)}"
assert len(tgt) == 19, f"expected 19 EM targets (MAXSP), got {len(tgt)}"
for want in ("DF","WL","LP","ES","PP"):
    hit = [em for a,fia,p,em in rows if a == want]
    assert want in hit, f"shared conifer {want} not mapping to self"
print("self-check OK: 442 rows, 19 targets, shared conifers map to self")
