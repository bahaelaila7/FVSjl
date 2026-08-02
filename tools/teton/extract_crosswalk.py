#!/usr/bin/env python3
"""Extract the TT species crosswalk (ASPT column 16) from bin/FVStt_buildDir/spctrn.f.

The shared western ASPT(442,21) table's column legend is
  ALFA FIA PLNT AK BM CA CI CR EC EM  IE KT NC PN SO TT UT WC WS OC OP
   1    2   3    4  5  6  7  8  9  10  11 12 13 14 15 16 17 18 19 20 21
so TT = column 16. The table is stored interleaved: for each group of 10 rows I,
a `DATA ((ASPT(I,J),J=1,10)...` block (cols 1-10: ALFA/FIA/PLNT + 7 targets) then a
`DATA ((ASPT(I,J),J=11,21)...` block (cols 11-21: 11 targets, TT = the 6th field = idx 5).
Join A-rows (alpha/fia/plants) with B-rows (TT target) by row order I=1..442.

Emits the KT 7-col schema (code_alpha,code_fia,code_plants,target×4 — TT has one target,
repeated) to data/teton/species_translation.csv. Faithful (programmatic, not hand-transcribed).
Self-check: 442 rows, 18 distinct targets (= MAXSP), shared conifers DF/LP/ES/PP/AF map to self.
"""
import re, collections
SRC = "/workspace/ForestVegetationSimulator/bin/FVStt_buildDir/spctrn.f"
OUT = "/workspace/FVSjl/data/teton/species_translation.csv"

a_rows, b_rows = [], []          # A: (alpha,fia,plants); B: tt-target
mode = None                      # 'A' | 'B' | None
for ln in open(SRC).read().splitlines():
    if re.search(r"DATA\s*\(\(ASPT\(I,J\),J=1,10\)", ln):   mode = 'A'; continue
    if re.search(r"DATA\s*\(\(ASPT\(I,J\),J=11,21\)", ln):  mode = 'B'; continue
    if mode is None: continue
    if ln and ln[0] in "Cc*": continue                     # Fortran comment — stay in block
    s = ln.lstrip()
    if s.startswith("&") or s.startswith("'"):
        f = [x.strip() for x in re.findall(r"'([^']*)'", ln)]
        if mode == 'A' and len(f) >= 3:  a_rows.append((f[0], f[1], f[2]))
        if mode == 'B' and len(f) >= 6:  b_rows.append(f[5])          # col 16 = TT
        if ln.rstrip().endswith("/"): mode = None
    elif s == "": continue
    else: mode = None

assert len(a_rows) == len(b_rows) == 442, f"A={len(a_rows)} B={len(b_rows)} (expected 442)"
with open(OUT, "w") as fo:
    fo.write("code_alpha,code_fia,code_plants,target_tt,target_tt2,target_tt3,target_tt4\n")
    for (a, fia, p), tt in zip(a_rows, b_rows):
        fo.write(f"{a},{fia},{p},{tt},{tt},{tt},{tt}\n")
tgt = collections.Counter(b_rows)
print(f"{len(a_rows)} rows, {len(tgt)} distinct targets -> {OUT}")
assert len(tgt) == 18, f"expected 18 TT targets (MAXSP), got {len(tgt)}: {sorted(tgt)}"
for want in ("DF", "LP", "ES", "PP", "AF"):
    hit = [tt for (a, fia, p), tt in zip(a_rows, b_rows) if a == want]
    assert want in hit, f"shared conifer {want} not mapping to self (got {hit})"
print("self-check OK: 442 rows, 18 targets, shared conifers map to self")
