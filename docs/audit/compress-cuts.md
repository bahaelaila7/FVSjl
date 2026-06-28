# Bandaid audit — COMPRESS + thinning/cuts

Files audited:
- `/workspace/FVSjl/src/engine/compress.jl`
- `/workspace/FVSjl/src/engine/cuts.jl`

FVS source consulted (SN build dir `bin/FVSsn_buildDir/`): `comprs.f`, `comcup.f`,
`cuts.f`, `cutstk.f`, `sdical.f`, `rdpsrt.f`, plus `sn/grinit.f`, `intree.f`.

The COMPRESS eigensolver replacement (LAPACK vs 1966 IBM-SSP `EIGEN`) is a KNOWN, accepted
divergence and is **out of scope**; the partition / PC scores are explicitly not bit-exact.
The flags below concern the **non-eigensolver** parts only.

Overall the two modules are largely faithful: `_rdpsrt!` is a label-exact port of `rdpsrt.f`
(incl. the tie-break that the comment correctly flags as load-bearing); the partial-last-record
formula `PREM=((XLEFT+CUT)/CUT)*PREM` (cuts.f:1366) matches; `_clsstk` uses the correct
`0.005454154` BA factor (cutstk.f:104); the THINDBH CUTEFF re-derivation matches the
label-355 / label-325 split; YARDLOSS scales merch/saw/board but not total-cubic/BA
(cuts.f:1394-1398); the MINHARV cancel-all gate maps to cuts.f:1567; `_sdi_zeide`/`_rd_curtis`
match SDICLS/RDCLS with SN's `LZEIDE=.TRUE.`, `DBHZEIDE=DBHSTAGE=0` (grinit.f:129,262-263);
`_thin_auto!` matches the AUTSTK/label-150 gate; the IMC clamp to 1..3 in `_cut_pref_wt`
faithfully mirrors intree.f:621-622 (verified — not a bandaid). The class-variable order
(HT, ICR, IMC, lnDBH, DG) and STDDEV floors match comprs.f:170-187, and the eigenvector
sign-fixes on EIVECT(4)/EIVECT(7) are correct.

---

## FLAG 1 — GAP: COMPRESS merge samples decay/defect/woodland-stems instead of averaging them

- jl symbol: `_merge_one!`, compress.jl:286-289
- Claim/comment (line 285): "nominal attributes from the PROB-sampled record (comprs.f:733-741)".
  The code assigns `decay_code`, `defect`, and `woodland_stems` from the RANN-sampled record `sel`.
- FVS source checked: comprs.f:733-741 is the genuinely-nominal block (ITRE, ISP, KUTKOD, ISPECL,
  IMC, IESTAT, IDTREE, NCFDEF, NBFDEF). But `DECAYCD`, `WDLDSTEM`, and `DEFECT` are **weighted-averaged**
  later in the same routine:
  - comprs.f:937 `DECAYCD(IREC1)=DECAYI/TXP`  (DECAYI accumulated `DECAYCD(IREC)*XP`, lines 644/700)
  - comprs.f:939 `WDLDSTEM(IREC1)=WDLDSTEMI/TXP`
  - comprs.f:963 `DEFECT(IREC1)=IDF11*1000000+...` (per-digit weighted mean via DF11..DF44, 673-680/963)
- Severity: GAP. The cited line range (733-741) does **not** cover decay/defect/woodland; those
  three were mis-classified as nominal. Source mandates a PROB-weighted average.
- Faithfulness impact: merged records get a single sampled decay/defect/woodland value rather than the
  class mean. Only exercised under COMPRESS (suite scenario s22, already accepted-divergent), so no
  impact on the validated suite, but it is a genuine source contradiction independent of the
  eigensolver caveat.

## FLAG 2 — GAP: TCONDMLT point-density weights (PBAWT/PCCFWT/PTPAWT) not ported

- jl symbol: `_cut_pref_wt` (cuts.jl:408-413) and the TCONDMLT handler (cuts.jl:163-164).
- Claim/comment (line 404, 163): the removal weight is "IORDER + TCWT·IMC + SPCLWT·ISPECL". The handler
  reads only `act.params[1]` (→ total_wt) and `act.params[2]` (→ special_wt).
- FVS source checked: cuts.f:1424-1428 sets `TCWT=PRMS(1); SPCLWT=PRMS(2); PBAWT=PRMS(3);
  PCCFWT=PRMS(4); PTPAWT=PRMS(5)`, and cuts.f:1071-1072 forms the RDPSRT key as
  `WK2 = XSZ + IORDER + TCWT*IMC + SPCLWT*ISPECL + PBAWT*PTBAA(IP) + PCCFWT*PCCF(IP) + PTPAWT*PTPA(IP)`.
- Severity: GAP. Params 3-5 of TCONDMLT are silently dropped, and the point BA/CCF/TPA terms are absent
  from the priority key.
- Faithfulness impact: when a user sets PBAWT/PCCFWT/PTPAWT, the RDPSRT removal order (and therefore which
  records are thinned and the post-thin RNG/DGSCOR traversal) diverges. Defaults are 0, so the common
  path is unaffected and the omission is untested.

## FLAG 3 — GAP (low): AUTSTK uses the wrong BA constant (0.005454154 vs FVS 0.0054542)

- jl symbol: `_autstk`, cuts.jl:691 — `tba = _BA_PER_TREE * d * d * wk4[i]` with `_BA_PER_TREE = 0.005454154`.
- FVS source checked: cutstk.f:43 (the AUTSTK entry) uses `TBA = 0.0054542 * DBH*DBH*PROB` — note this is
  the *only* place FVS uses 0.0054542 for a cut; CLSSTK (cutstk.f:104) and the cut loop (cuts.f:1209) use
  0.005454154. FVSjl applies the single `_BA_PER_TREE` constant everywhere.
- Severity: GAP (low). The constant cancels in `TMPMAX = TEMBA/TOTBA`, so the SDImax weighting is
  unaffected. It survives only in the gate `TOTBA.LE.1 .OR. TEMBA.LE.1` (cutstk.f:46), where the ~8e-6
  relative difference could flip the result only at a knife-edge TOTBA≈1.
- Faithfulness impact: effectively nil; flagged because it is a literal source mismatch, not an
  output-matched choice.

## FLAG 4 — GAP (low): sorted SDI/RDEN/CC partial-tree + efficiency differs from cuts.f

- jl symbols: `_thin_sdi!`/`_thin_rden!`/`_thin_cc!` sorted branches (icut/dir > 0), e.g. cuts.jl:550,
  565-569 (and the RDEN/CC analogues).
- Claim/comment: "last record is partial (prem = remaining / weight)" using a direct formula and the
  user efficiency `ce = cuteff_p > 0 ? cuteff_p : 1`.
- FVS source checked:
  - cuts.f:890-894 sets the per-record `CUTEFF` for the sorted density paths to effectively
    `max(REMOVE/SDIC, CUTEF1)` (user efficiency only wins when it is the larger), not the raw user value.
  - The last-tree partial is solved **iteratively** (cuts.f:955-1363: ±0.05 for THINSDI, ±0.2 for
    THINRDEN, converging |DIFF|≤0.5 / ≤0.01) rather than by a single closed-form division. For Zeide SDI
    the closed form is exact (weight is fixed `(D/10)^1.605`), so FVSjl matches there; for Curtis RD/CC the
    factors are stand-QMD-dependent and FVS's iteration can land up to ~0.2 TPA off the direct value.
- Severity: GAP (low). Only the sorted (from-below/above) density-target thins with an explicit cut
  efficiency are affected; the common SN form is the proportional throughout path (icut=0), which is
  faithful (`CUTEFF=REMOVE/SDIC`).
- Faithfulness impact: small per-cut redistribution on a single boundary record in an untested edge.

## FLAG 5 — GAP (low): truncated-tree merge rounds NORMHT/ITRUNC instead of truncating (IFIX)

- jl symbol: `_merge_one!`, compress.jl:262 — `normnew = round(Int32, xnr/txp*100); truncnew = round(Int32, xit/txp*100)`.
- FVS source checked: comprs.f:805-806 `NORMHT(IREC1)=IFIX(XNR/TXP*100.)`, `ITRUNC(IREC1)=IFIX(XIT/TXP*100.)`
  — `IFIX` truncates toward zero; ICR uses `NINT` (comprs.f:957), which the port correctly rounds.
- Severity: GAP (low). Off-by-≤1 on merged truncated-tree normal-height / truncation point, only under
  COMPRESS of a class containing truncated trees (accepted-divergent module).
- Faithfulness impact: negligible; noted for completeness.

---

### Minor note (not separately flagged)
MINHARV vs YARDLOSS ordering: `cuts!` evaluates the MINHARV minimum-harvest gate (cuts.jl:210-221) on the
**un-scaled** merch/saw/board totals, then applies the `(1−PRLOST)` yardloss scaling afterward
(cuts.jl:229-234). FVS accumulates CMCUT/BFCUT/SCCUT already multiplied by `(1−PRLOST)` inside the loop
(cuts.f:1395-1397) and tests those scaled values (cuts.f:1567). When BOTH MINHARV and YARDLOSS are active
the gate is marginally easier to pass in FVSjl. Both default off; untested combination.
