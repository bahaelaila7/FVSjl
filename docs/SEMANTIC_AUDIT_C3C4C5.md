# C3/C4/C5 semantic-congruence audit

Branch-by-branch audit of each ported C3 (growth) / C4 (mortality+density+regen) /
C5 (volume+biomass+stats) routine against its FVS Fortran source — to verify the
LOGIC maps over, not just that the test cases pass. Method (the §MORTS template):
for each routine, enumerate the Fortran decision branches, then mark
**[I]** implemented in FVSjl? and **[T]** exercised by a test scenario? A branch that
is unimplemented OR untested is a gap; untested branches get a designed scenario.

Oracle source: `/workspace/FVSjulia/src` (1:1 Fortran transliteration).

Legend: ✅ implemented+tested · 🟡 implemented, UNtested · ⛔ NOT implemented · 🔬 accuracy

---

## GAP LOG (findings; fix + scenario tracked here)

### G1 — Volume: CULL / DEFECT / DECAYCD / WDLDSTEM not applied  ⛔  (C5)
- **Fortran (`vols.f`):** `CULL` and `DEFECT` (decoded to `NCFDEF`/`NBFDEF` cubic/board
  defect %) reduce merch/board volume inside `NATCRS`/`CFVOL`; `DECAYCD>0` selects
  dead-tree decay-class carbon factors (vols.f:154-164); `WDLDSTEM` drives woodland
  multi-stem volume.
- **FVSjl:** `volume.jl` references **none** of `cull/defect/decay_code/woodland_stems`
  → a tree carrying cull/defect/decay/woodland gets the wrong (too-high) volume.
- **Why untested:** snt01's TREEFMT has no cull/defect/decay/woodland columns → all
  default 0 → the branch never fires → snt01 stays bit-exact and hides the gap.
- **Plan:** port the defect/cull/woodland/decay handling in the volume path; add a
  `.tre` scenario with cull%, defect codes, woodland stems, and a dead+decay tree;
  validate the cuft/bdft/carbon vs live Fortran.
- Status: ⛔ OPEN

### G2 — Biomass / carbon not computed at all  ⛔  (C5)
- **Fortran (`vols.f` → `calcbiomass.f`):** `JENKINS` (Jenkins 2003 above-ground
  biomass) + `WOODDEN` (wood/bark density) fill `ABVGRD_BIO/MERCH_BIO/CUBSAW_BIO/
  FOLI_BIO` and the carbon arrays (×carbon fraction, decay-class adjusted) every cycle.
- **FVSjl:** the `TreeList` biomass/carbon fields exist but are **never written** —
  `volume.jl` computes no biomass/carbon. JENKINS/WOODDEN are unported.
- **Why untested:** the `.sum` has no biomass/carbon columns (they feed the DBS Carbon
  tables, C6), so nothing surfaces the gap.
- **Plan:** port JENKINS + WOODDEN (CSV-drive WDBKWT like the other coef tables);
  validate per-tree biomass/carbon vs Oracle A (dump, not `.sum`). Lower priority than
  G1 (not in `.sum`; consumed by C6).
- Status: ⛔ OPEN

---

## Routine-by-routine audit (coverage confirmations)

The 90-species `all_*` sweep + dense_long(30cy) + bare_multipoint + snt01's own trees
exercise the SPECIES- and DENSITY-dependent branches; snt01 carries dead trees
(history 6-8 ×10), 3 dubbed heights, and broken-top/topkill records, so those branches
are covered. Confirmed branch-complete + covered (bit-exact ±2 ulp):

| routine (Fortran → FVSjl) | branches covered by | status |
|---|---|---|
| `DGF` → `diameter_growth!` | forest-type(8 grps)/physiography/Fort-Bragg/size-cap/calib — 90-species sweep + s30/s31 + snt01 | ✅ |
| `MORTS`+`VARMRT` → `mortality!` | §MORTS checklist; BAMAX/line-reset/Hamilton/SDImax — dense_long + all_* | ✅ |
| `DGSCOR`/`AUTCOR`/`TRIPLE` → serial_correlation + triple_records! | tripling cyc1-2 + stochastic cyc3+ — dense_long 30cy | ✅ |
| `CROWN`(MCREQN 1-5) → crown_ratio_update! | all 5 eqn types across 90 species; ±change cap | ✅ |
| `CRATET` → dub_missing_heights! | missing-height dub, broken-top, dead-tree dub, topkill | ✅ (snt01) |
| `HTGF`/`HTCALC`/`HTDBH` → height_growth! | height-age + inverse — 90 species | ✅ |
| `REGENT` → small_tree_growth! + establish! | small-tree blend, ESSUBH, NPTIDS>1 — bare_* | ✅ |
| `DENSE`/`SDICAL`/`PTBAL`/`SDICHK` → compute_density! | BA/SDI/SDImax/point-BA/PCT — 90 species + backdated calib | ✅ |
| `VOLS`/`CFVOL`/`NATCRS`/R8-Clark → compute_volumes! | METHC methods per species; topkill/CFTOPK — 90 species + s30/s31 | ✅ except G1 attributes |
| `COMCUP` → comcup! | zero-PROB delete — dense_long | ✅ |

**Net:** C3/C4 and the species/density side of C5 are branch-covered and bit-exact.
The two genuine semantic gaps are **G1** (tree-attribute volume: cull/defect/decay/
woodland) and **G2** (biomass/carbon) — both invisible because no scenario sets those
attributes and the `.sum` omits biomass. Fix order: G1 (affects `.sum`) then G2.
