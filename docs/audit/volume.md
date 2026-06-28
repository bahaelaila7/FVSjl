# Bandaid Audit — Volume (Clark R8/R9, CFVOL/NATCRS, board-ft, topkill/Behre)

Module files audited:
- `/workspace/FVSjl/src/engine/volume.jl` (VOLS/CFVOL driver + CFTOPK/BFTOPK + CRATET height dubbing + VOLUME/BFVOLUME overrides + defect)
- `/workspace/FVSjl/src/engine/r8clark_vol.jl` (`_R8CLARK_VOL`, the R8 New-Clark taper, Scribner board feet)
- `/workspace/FVSjl/src/engine/r9clark_vol.jl` (NE R9 Clark, WIP — not on the SN path)
- `/workspace/FVSjl/src/engine/volume_equations.jl` (VOLEQDEF / `_r8_ceqn` geoa + SNFIA→SNSP crosswalk)

FVS sources read: `bin/FVSsn_buildDir/behprm.f`, `behre.f`, `cftopk.f`, `bftopk.f`, `fvsvol.f`, `vols.f`, `r8prep.f`, `r9clark.f`, `r9logs.f`, `mrules.f`, `sn/cratet.f`, `voleqdef.f`.

## Verified faithful (not flagged)
Spot-confirmed against source and judged faithful: BEHPRM hyperbola constants and BEHRE
integral (behprm.f:27-43 / behre.f:28-31, byte-for-byte); CFTOPK mcf+scf truncation
branches and the `MCF>TCF`/`SCF>MCF` clamps (cftopk.f:70-132); BFTOPK cone/non-cone
ratios (bftopk.f:46-67); `vmax = vol[1]` passed to BOTH cftopk and bftopk (fvsvol.f:495
`BFMAX=TVOL(1)`, :510 `VMAX=TCF`, :511 `IF(BFMAX.EQ.0)BFMAX=VMAX`); MCF = vol[4]+vol[7]
gated on D≥DBHMIN and SCF = vol[4] gated on D≥SCFMIND (fvsvol.f:512-517); the Region-8
"≥10 ft of product" rule on the primary cubic call (fvsvol.f:343-346 — the jl's use of
the *un-zeroed* `rawSawHt` as HT1PRD is provably equivalent because r9clark.f:320 zeros
sawHt only below 9.5 ft, which is < 10 and fires the rule regardless); MRULES R8 defaults
maxLen=8/minLen=2/merchL=8 (`08`→12)/trim=0.5/stump 1.0|0.5/mTopp 7|9/mTops 4
(mrules.f:337-369); short-tree topHt=17.4 + ×shrtHt/17.3 scaling (r8prep.f:328-338,
r9clark.f:226); `_r8_fcmin_adj` form-class minima and the THT<47.5 floor
(r8prep.f:346-365); `_r8_remap_spec` species folding (r8prep.f:99-116); the `_r8clark_lookup`
binary search + geoa=9 retry + spgrp∈{100,300,500} guard (r8prep.f:142-192); the cubic
defect ICDF/PULPV path (vols.f:294-332, eastern branch MCFV=PULPV+SCFV); R9LOGS/R9LOGLEN
sawtimber-only log segmentation incl. the even-foot resegmentation (r9logs.f:61-102,
246-279; board uses NOLOGP only, r9clark.f:379); R9BDFT Scribner scaling at top-of-log
INT(DIB+0.499) with per-log NINT (r9clark.f:1447-1488, scribner table matches DATA);
`_vol_geoa` 12-forest map + SNFIA/SNSP crosswalk + the spec<300?22:110 fallback
(voleqdef.f R8_CEQN:1929-2089); KODFOR decode in `setup_volume_equations!`. The extra
`(X-Y)>1e-10` guard in `r8clark_vol._r9cuft`/`_r9ht` (absent from r9clark.f:1020) and the
`max(dib17,0.1)` floor are inert defensive additions: AFI is negative for every SN species
(checked `r8clark_cf.csv`), so `(dib17-AFI)/BFI > dib17` always and the omitted
r8prep.f:507 `COEFFSO%DIB17=max(.,COEFFS%DIB17)` clamp can never fire; X>Y always for real
trees. None change tested output.

---

## FLAG 1 — Input board-defect not applied to sawtimber cubic when board feet is zero
- **Symbol / line:** `compute_volumes!`, volume.jl:430-446 (the `if bf > 0f0 … end` board-defect block, and `mcf = pulpv + scf` at :446)
- **Claim:** The comment (volume.jl:410-413) says board feet AND the sawtimber cubic are
  cut by IBDF, "applied only where board feet exist."
- **FVS source checked — vols.f:393-421:** the IBDF *curve/log-linear* augmentation is
  gated on `IF(BFV(I).GT.0.0 .AND. LBVOLS)` (line 393), but the **application** of IBDF is
  NOT: lines 415-421 run unconditionally after label 100 —
  `IF(IBDF.LT.99) THEN BFV=BFV*(1-IBDF/100); IF(.NOT.LFIANVB) SCFV=SCFV*(1-IBDF/100)`.
  So when a tree carries an *input* per-tree board defect (IBDF>0 from the DEFECT field)
  but its board feet are zero — e.g. `D<BFMIND` with `SCFV>0` after a BFVOLUME override
  raises BFMIND above SCFMIND (vols.f:354 `GO TO 100` keeps the defect step) — FVS still
  reduces the reported sawtimber cubic SCFV by IBDF. The jl wraps the entire IBDF
  application inside `if bf > 0f0`, so it leaves `scf` untouched and reports
  `mcf = pulpv + scf` with the full sawtimber cubic.
- **Severity:** GAP
- **Faithfulness impact:** Reported SCFV / merch-cubic too high by the input board-defect %
  for trees with zero board feet but a per-tree board-defect code. Unreachable on snt01
  (no defect input, BFMIND=SCFMIND); requires a BFVOLUME override + per-tree board defect.

## FLAG 2 — `round(Int32, …)` used to emulate Fortran `INT(…)` in height dubbing
- **Symbol / line:** `dub_missing_heights!`, volume.jl:213 (`norm_ht = round(Int32, h_v*100+0.5)`),
  :216 / :218 (`trunc = round(Int32, 80*ht+0.5)`), :230 (`norm_ht = round(Int32, ht*100)`)
- **Claim:** Comment cites cratet.f:212-265 for the CRATET height resolution (the actual
  dubbing math is cratet.f:376-397).
- **FVS source checked — cratet.f:381,384,386,397:** `NORMHT=INT(H*100.0+0.5)`,
  `ITRUNC=INT(80.0*HT+0.5)`, `NORMHT=INT(HT*100.0)`. Fortran `INT` truncates toward zero;
  `INT(x+0.5)` is round-half-**up**, and `INT(x*100)` is a pure truncation. Julia
  `round(Int32, y)` is round-to-**nearest, ties-to-even**. Thus `round(Int32, x+0.5)`
  double-applies the bias and `round(Int32, x*100)` rounds instead of truncating — the
  results differ by 1 unit for roughly half of the fractional inputs.
- **Severity:** GAP
- **Faithfulness impact:** NORMHT/ITRUNC are stored in **centi-feet**, so the deviation is
  ≤0.01 ft of normal/break height, feeding the Behre topkill redistribution (H=NORMHT/100).
  Magnitude is negligible, but it is a systematic, non-bit-faithful divergence from the
  cited `INT` semantics (would matter only as sub-cuft noise, never structurally).

## FLAG 3 — Final volume rounding uses ties-to-even instead of NINT (ties-away)
- **Symbol / line:** `_R8CLARK_VOL`, r8clark_vol.jl:433-436
  (`vol[1]=Float32(round(vol[1]*10)/10)`, same for vol[4], vol[7]; `vol[10]=round(vol[10])`)
- **Claim:** Comment cites "r9clark.f:456-461" / "Round to match Fortran nint()".
- **FVS source checked — r9clark.f:456-462:** `vol(1)=nint(vol(1)*10.0)/10.0`,
  `vol(4)=nint(vol(4)*10.0)/10.0`, `vol(7)=nint(...)`, `vol(2)/vol(10)=nint(...)`. Fortran
  `NINT` rounds halves **away from zero**; Julia `round` defaults to `RoundNearest`
  (ties-to-**even**). On an exact `vol*10 == k+0.5` tie the two differ by 0.1 cuft.
  (vol[10] board feet is already an integer sum of per-log `RoundNearestTiesAway` values,
  so its final `round` is a no-op and is unaffected.)
- **Severity:** GAP
- **Faithfulness impact:** Effectively ULP-class — a tie requires `cfVol*10` to land
  exactly on a representable `.5`, which essentially never occurs for the transcendental
  taper integral. Listed for completeness; the cubic columns could in principle differ by
  0.1 cuft on a measure-zero tie. Replacing with `RoundNearestTiesAway` would make it
  bit-faithful to NINT.

---

### Net
The SN R8-Clark volume spine (taper integral, CFTOPK/BFTOPK Behre topkill, BFMAX/VMAX,
DBHMIN/SCFMIND gating, Region-8 ≥10-ft rule, MRULES defaults, defect, Scribner board feet,
species crosswalk/geoa) is faithful to source. The only genuine logic gap is FLAG 1 (board
defect on zero-board sawtimber trees), reachable only via a non-default BFVOLUME override.
FLAGS 2-3 are source-grounded rounding-mode deviations of negligible (≤0.01 ft / ≤0.1 cuft,
tie-only) magnitude. No output-matching bandaids were found.
