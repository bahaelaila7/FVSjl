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

## ⚑ RE-SCOPING (after reading intree.f field layout)

The `.sum` is **semantically congruent** — G1/G2 do NOT silently corrupt it. Verified:
`intree.f` reads exactly 25 fields (same as FVSjl) and **CULL/DEFECT are never read
from a text `.tre`** — they're zeroed for text input and only set via the DBS database
path. For text records: CULL=DEFECT=0 (no volume reduction), DECAYCD=3 dead / 0 live
(carbon only, not `.sum`), WDLDSTEM=1 woodland (single-stem ⇒ no change). So every
attribute branch is a no-op for the text-input `.sum`, which is why the 90-species
sweep is bit-exact. **G1 and G2 are therefore C6-coupled** (DBS database input + the
DBS Carbon/biomass tables), not natural-process `.sum` gaps. Still ported here for
completeness + tested per-tree (not via `.sum`), per the user's request.

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
- **Deeper (C1 input):** FVSjl's `parse_tree_record` reads only through field 25
  (birth_age) — it never reads CULL/DEFECT/DECAYCD at all. And `intree.f:123-126`
  DEFAULTS `DECAYCD=3` for dead trees (ith>5) / 0 for live; FVSjl does neither. So the
  gap spans **C1 (parser doesn't read them) → C2 (no dead-tree DECAYCD default) → C5
  (volume/biomass don't apply them)**.
- **Plan (multi-layer):** (1) extend TreeRecord + parser to read cull/defect/decay
  (and the woodland default); (2) intree DECAYCD=3-for-dead default; (3) apply
  defect%/cull in the volume path + woodland multi-stem; (4) a `.tre` scenario with
  cull/defect/woodland + dead/decay trees, validated vs live Fortran.
- Status: ⏭️ DEFERRED TO C6 (see docs/C6_DBS_TODO.md) — DBS-input + Carbon-table coupled

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

## Designed test coverage for G1/G2

- **G2 biomass/carbon** — extend the 90-species harness to dump **per-tree**
  `ABVGRD_BIO/MERCH_BIO/CUBSAW_BIO/FOLI_BIO` + carbon and diff vs Oracle A (the `.sum`
  can't show it). The 90-species relabel already spans the Jenkins species groups, so
  this gives broad coverage. New test: `test_biomass.jl` (per-tree, ±tol).
- **DECAYCD-dead / WDLDSTEM defaults** — testable from text `.tre` now: assert mine's
  dead trees (history 6-9) get `decay_code=3` and woodland species get
  `woodland_stems=1`, matching the oracle (per-tree state, not `.sum`).
- **G1 cull/defect volume** — CANNOT be set from a text `.tre` (Fortran reads 25 fields;
  cull/defect come from DBS database input). Genuine coverage needs **DBS input = C6**.
  Port the volume-side application now (ready); add the cull/defect scenario as a C6
  DBS-input test. Tracked.

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
