# Bandaid audit — Stump sprouting (ESUCKR / SPROUT)

File audited: `src/engine/sprout.jl`
FVS sources checked: `bin/FVSsn_buildDir/esuckr.f` (STRP), `bin/FVSsn_buildDir/essprt.f`
(ESSPRT/NSPREC/SPRTHT/ESASID entry points), `bin/FVSsn_buildDir/cwcalc.f`,
`bin/FVSsn_buildDir/forkod.f`, `bin/FVSsn_buildDir/grinit.f`, plus
`data/southern/sprout_essprt.csv`.

## Summary

The three pure helpers (`nsprec_sn`, `sprtht_sn`, `essprt_sn`) and the per-species
ESSPRT coefficient CSV are bit-faithful to the SN `SELECT CASE` blocks — I verified
every SN species row in `sprout_essprt.csv` against `essprt.f:514-591`, the NSPREC SN
block (`essprt.f:1102-1120`), and the SPRTHT SN curve set (`essprt.f:1387-1393`). The
aspen omission is correct (ESASID has no SN case → INDXAS=9999, `essprt.f:751-752`, so
ASSPTN never fires). `nint` uses `RoundNearestTiesAway` = Fortran NINT.

Five non-trivial decisions are flagged below: one source-contradicting value (crown
width), and four GAPs on untested paths.

---

## FLAG 1 — BANDAID: crown-width call passes ICR(70) where FVS passes a dummy CR of 1.0

- jl symbol: `esuckr!`, `src/engine/sprout.jl:171-172`
  ```julia
  cw = crown_width(coef, sp2, dbh, ht, 70f0, 0, ...)   # 5th arg = crown ratio
  ```
- FVS source: `esuckr.f:313-315`
  ```
  CRDUM=1.
  CALL CWCALC(ISSP,PROB(ITRN),DBH(ITRN),HT(ITRN),CRDUM,ICR(ITRN),CW,0,JOSTND)
  ```
  CWCALC's 5th arg is `CR = CROWN RATIO IN PERCENT` (`cwcalc.f:1`, doc line 19); the 6th
  arg `IICR` is explicitly discarded (`cwcalc.f:82  IDANUW = IICR`). So FVS feeds the
  equation **CR = CRDUM = 1.0**, not ICR=70.
- The SN forest-grown Bechtold equations use CR directly, e.g. `cwcalc.f:667`
  `CW = 0.6564 + 0.8403*D + 0.0792*CR`. `_cw_eval` (`crown_width.jl:43`) likewise adds
  `e.cr_coef * cr`. With cr=70 vs 1.0 the crown width is inflated by
  `cr_coef*(70-1)` (≈ 5.5 ft for cr_coef≈0.0792) per sprout record.
- Severity: **BANDAID** (source mandates 1.0; jl conflated the CR slot with ICR and
  dropped the unused IICR arg).
- Faithfulness impact: every newly-created sprout gets a too-large CRWDTH at birth.
  Propagation is likely partial (CRWDTH is recomputed by CWIDTH the next growth cycle),
  but it diverges in the regeneration cycle's CCF/density/structure outputs. Fix: pass
  `1f0` in the CR slot.

## FLAG 2 — GAP: special-forest ESSPRT gate keyed on KODFOR, but FVS keys on ISEFOR

- jl symbol: `esuckr!` line 145 `isefor = Int(s.plot.user_forest_code)`, consumed by
  `essprt_sn` / `_es_special_forest` (`sprout.jl:84-85,101`).
- `s.plot.user_forest_code` is **KODFOR** (the 5-digit region/forest/district code, e.g.
  80907) — `state.jl:365` comment "(KODFOR)", set raw from the keyword in
  `keyword_dispatch.jl:430` and treated as KODFOR at `:447` (`div(...,100)==701`).
- FVS ESSPRT tests **ISEFOR** (`essprt.f:545,552,559,569,576`), which `forkod.f:513`
  sets to `JFOR(IFOR)` — the bare 3-digit forest code (809, 810). KODFOR is a different
  value: `forkod.f:514  KODFOR=(JFOR(IFOR)*100)+KODIST` (= 80907). So
  `Int(KODFOR)==809` is essentially never true.
- Consequence: the special-forest survival equations for species 64/66/70/75/77 on NFs
  809/810 (e.g. `essprt.f:547` `(57.3-0.0032*DSTMP**3)/100`) **never fire** in the jl;
  it always falls to the common logistic (the ELSE branches). Also note 905/908 in
  `_es_special_forest` are dead even in FVS (`forkod.f:440` sets `ISEFOR=0` for region-9
  forests), so mirroring them is harmless but the KODFOR mismatch makes the whole branch
  unreachable.
- Severity: **GAP** (faithful for the SN test stands, which are region-8 standard
  forests like Talladega 801 → not in the special set; silently wrong for stands on NF
  809/810). Fix: store/compare the post-`forkod` ISEFOR (JFOR), not KODFOR.

## FLAG 3 — GAP: SPROUT keyword per-species / species-group / DBH-range table not honored

- jl symbol: `esuckr!` uses one stand-level pair `s.control.sprout_smult` /
  `sprout_hmult` for every sprouting species (`sprout.jl:141-143,158,163`); the docstring
  (`sprout.jl:135-137`) admits the OPGET table is "a later refinement."
- FVS `esuckr.f:96-150` builds per-keyword-instance `SPRMLT/HTMSPR/DMIN/DMAX` arrays
  indexed by species (single species, species **group** `J<0`, or all), and
  `esuckr.f:199-205` selects SMULT/HMULT **per tree record** only when
  `DSTMP ∈ [DMIN,DMAX)`. The jl ignores the species selector entirely:
  `keyword_dispatch.jl:1154-1157` applies the multipliers globally whenever field-2 is
  present, so e.g. `SPROUT 33 0.5` (species 33 only) would in jl scale **all** species,
  and any DBH-range form is dropped.
- Severity: **GAP** (faithful only for the common single-pair-for-all SPROUT form;
  diverges for species-targeted or DBH-range keywords — none exercised by the SN test
  suite). Cite `esuckr.f:108-148` (species/group decode) and `:200-204` (DBH-range).

## FLAG 4 — GAP: no ESCPRS compression when the tree list fills; sprouts silently dropped

- jl symbol: `esuckr!` line 160
  `n = t.n + 1; n > length(t.dbh) && break    # no ESCPRS compression — list-overflow guard`
- FVS `esuckr.f:252-256`: when `ITRN >= MXRR` it calls `ESCPRS(ITRGT,DEBUG)` to compress
  the live tree list (down toward 70% of MXRR) and then **continues** creating the
  remaining sprout records. The jl instead `break`s, abandoning the rest of that record's
  sprouts (and, with the outer loop, subsequent records can't grow the list either).
- Severity: **GAP** (untested — SN test stands never approach MAXTRE during regen). On a
  list-saturated stand FVS keeps the sprouts via compaction; jl loses TPA. Note: this is
  the establishment-specific ESCPRS, distinct from the thinning-only compaction noted in
  prior memory, so the "compaction is thinning-only" generalization does not cover it.

## FLAG 5 — GAP (minor): `sprout_dbh` ignores the per-stand AA height–diameter refit

- jl symbol: `sprout_dbh` (`sprout.jl:75-81`) always uses `wykoff_ht1` as AX.
- The docstring's cited defaults are correct: `grinit.f:104 LHTDRG=.FALSE.`,
  `grinit.f:105 IABFLG=1`, so `esuckr.f:298-302` takes `AX=HT1` by default. But when a
  HT-DBH keyword flips `LHTDRG` on, CRATET sets `IABFLG=0` and `AX=AA(ISSP)`
  (`esuckr.f:300-301`), which the jl never applies.
- Severity: **GAP** (faithful for every default/standard SPROUT stand; silently wrong
  only for a stand that both sprouts and enables the species HT-DBH regression — an
  untested combination, and the docstring discloses it).
