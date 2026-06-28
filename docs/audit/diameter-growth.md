# Bandaid audit — Large-tree diameter growth + calibration + tripling

Module: `src/variants/southern/diameter_growth.jl`, `bark_and_bounds.jl`,
`serial_correlation.jl`
FVS sources checked: `sn/dgf.f`, `sn/dgdriv.f`, `sn/dgbnd.f`, `sn/bratio.f`,
`base/autcor.f`, `base/dgscor.f`, `base/grincr.f`, `base/gradd.f`, `base/dense.f`,
`sn/grinit.f`, `sn/blkdat.f`.

## Summary

The core large-tree DG spine is **faithful**. I traced and confirmed against FVS
source: the DGF DDS terms and the folding of site/slope/aspect/physiographic terms
into DGCON (dgf.f:252,290-296, debug block 305-310 shows the site/slope terms live in
CONSPP, not the DDS sum); the −9.21 / RELHT-1.5 / BA-25 clamps (dgf.f:266,270,368);
the Fort Bragg DG5 longleaf/loblolly equations and bark overrides (dgf.f:347-361,
bratio.f:151-164); BRATIO 0.80/0.99 clamps (bratio.f:173-174); DGBND taper 0.048 +
0.90 interpolation + size-cap `SIZCAP(,3)<1.5` (dgbnd.f:131-146); AUTCOR variance/
covariance loops (autcor.f:44-92, YR=5 from blkdat.f:61); DGSCOR bounded-normal taper
at DDS>4/>5 (dgscor.f:35-49); the tripling constants FU/FM/FL=1.271/−0.14228/−1.549
(dgdriv.f:638-640); PSIGSQ=0.089827273 (dgdriv.f:96); the empirical-Bayes COR fit,
DIST-weighted ratio/regression blend, exp(COR)∈(0.0821,12.1825) trap, and SIGMA/VARDG
(dgdriv.f:486-634); the backdating BAGR / WK3=√(D²·R) and IDG=1 bark skip
(dense.f:100-128); OLDRN seeding (regression endpoints vs BACHLO redraw) and the
±DGSD·SIGMA clamp (dgdriv.f:544-606,649-653); and the species-sorted RNG order /
TRIPLE interleaved append (triple.f layout). All bit-exact-consistent with source.

I also verified two items that *look* suspicious but are faithful:

- **COR attenuation clock (diameter_growth.jl:563-577).** The jl uses
  `cormlt = exp(-0.02773·elapsed)` with `elapsed = IY(ICYC)−IY(1)` for COR but
  `cormlt_h = exp(-0.02773·(elapsed+sfint))` (= `IY(ICYC+1)−IY(1)`) for HCOR. The
  comment mis-cites "FVS SFINT = IY(icyc)−IY(1)"; the literal Fortran is
  `SFINT=IY(ICYC+1)-IY(1)` (dgdriv.f:147) used for **both** COR and HCOR
  (dgdriv.f:192-194). It is nevertheless faithful: FVS calls `DGF` at dgdriv.f:136
  **before** attenuating COR at line 193, so the COR that DGF actually consumes in
  cycle ICYC was attenuated at the *end of the previous cycle* with
  `SFINT=IY(ICYC)−IY(1)`. HCOR has no such lag (consumed by the height model after
  the attenuation), so it correctly uses `IY(ICYC+1)−IY(1)`. The jl reproduces both
  clocks exactly. Not a bandaid.

- **`t.crown_ratio` in the DGF competition term (dgf.jl:165-168).** Despite the name
  and the loose comment "PCT is the crown-modeled crown ratio", the value stored is
  the BA percentile PCT (filled by `stand_pct!`/the calibrate percentile loop), which
  is exactly FVS `PCT(I)` in `BAL=(1-PCT/100)·BA` / `PBAL=PBA·(1-PCT/100)`
  (dgf.f:272,284). The crown-ratio used in the `ln(ICR)` term is the separate
  `crown_pct` field. Faithful (only the comment's terminology is wrong).

Two real concerns are flagged below.

---

## FLAG 1 — DGBND/size-cap applied AFTER the FINT cycle-length scaling (GAP)

- jl symbol/line: `diameter_growth!`, diameter_growth.jl:622-631 (and the tripling
  bounds 627/629/631), using `dg_bound` (bark_and_bounds.jl:41-57).
- What the jl does: it folds the cycle-length scaling into the DDS —
  `dds = exp(wk2[i]) * xbai * (sfint/5)` — then computes the **FINT-year** increment
  `sqrt(d_ib² + dds·frm) − d_ib` and applies `bnd(...)` (the full taper + 0.048 floor
  + size cap) to that already-scaled increment.
- FVS source checked: `dgdriv.f:206,214,223,267-269` computes DDS=EXP(WK2+XDGROW) on
  the **YR(=5)-year** basis, builds the 5-year DG, and calls `DGBND` on that 5-year
  DG. The cycle-length scaling happens **later** in `gradd.f:77-90`:
  ```
  IF (ITRN.GT.0 .AND. FINT.NE.YR) THEN
     SCALE=FINT/YR
     ...
     DDS=(DG(I)*(2.0*BARK*D+DG(I)))*SCALE
     DG(I)=SQRT((D*BARK)**2+DDS)-BARK*D
  ```
  — and `gradd.f` does **not** re-call `DGBND` after scaling. (grincr.f:432-433 even
  documents "DG WILL RETURN [from DGDRIV] WITH DIAMETER GROWTHS SCALED TO YR-YEAR
  BASIS"; gradd then rescales to FINT.) The jl comment cites "dgdriv.f:325/715" for
  this scaling, but those lines are the *calibration* `SCALE=YR/FINT` and the
  missing-increment dub-in `SCALE=1./SCALE` — neither is the growth-mode scaling,
  which is `gradd.f:80`.
- Why it matters / faithfulness impact: the DDS-space scaling algebra is identical, so
  for the standard SN 5-year cycle (`sfint=5`, all of snt01) the two orders coincide
  and results are bit-exact. They **diverge whenever FINT≠YR**, because the bound is
  applied on different-magnitude increments and is itself not scaled:
  - A tree above its DLODHI upper bound gets `DDG=0.048` in FVS *on the 5-year DG*,
    which gradd then scales up (≈0.096 for a 10-yr cycle); the jl caps the 10-yr DG at
    a flat `0.048` → ~½ the large-tree increment.
  - The size cap (`DBH+DDG>SIZCAP`) is enforced by FVS on the 5-year DG and then
    scaled past the cap; the jl enforces `DBH+DG_fint≤SIZCAP` (tighter).
  Severity GAP: faithful on the tested 5-yr path; for the 10-yr scenarios (s5/s9/
  timeint10) it could diverge if any tree sits above DLODHI/SIZCAP. I could not
  confirm whether those scenarios exercise a tree above the bound, so the practical
  impact on the existing suite is unverified, but the ordering is a genuine departure
  from `gradd.f:79-90`.

## FLAG 2 — Calibration overrides PTBAA to CURRENT diameters (UNVERIFIED)

- jl symbol/line: `calibrate_diameter_growth!`, diameter_growth.jl:267-281 and 320
  (the `cur_point_ba` capture at current/dead-exposed dbh, then
  `s.density.point_ba .= cur_point_ba` after the backdated `compute_density!`).
- The claim: an elaborately-commented block asserts that the calibration `dgf!` must
  read PTBAA built at **current** diameters (with every dead tree, history 6 and 8,
  contributing its fixed point-BA), while stand BA/AVH/PCT stay backdated — comment
  "dgf.f:493 reads the live PTBAA the last DENSE pass filled at current diameters".
- FVS source checked: `dgf.f` only has 1172 lines (no line 493 in DGF; the PTBAA read
  is `dgf.f:282 PBA=PTBAA(ITRE(I))`). The DENSE that fills PTBAA does so via
  `dense.f:280 CALL PTBAL`, which consumes `WK5`. In the backdating (LREDO) pass
  `dense.f:184` sets `D=WK3(I)` (the **backdated** dbh) and `WK5(I)=D*DP`
  (dense.f:187). So the PTBAL/PTBAA produced by the backdating DENSE pass is built
  from **backdated**, not current, diameters — the opposite of the jl override. Long-
  dead trees (`IMC=9`) are zeroed in WK3 (dense.f:86), again contradicting the jl's
  "PTBAA exposes long-dead trees at current dbh".
- Why I cannot fully resolve it: the calibration LSTART DENSE/DGF call sequence
  (whether a *current-dbh* DENSE pass runs after the LREDO backdating pass and leaves
  its PTBAA in the common array before `DGDRIV(LSTART)`'s `CALL DGF(WK3)`) is in the
  LSTART driver I did not trace (grincr.f:305-311 is the per-cycle post-thin DENSE,
  not the calibration path). The jl produces snt01 bit-exact, which argues it is
  faithful — but the source I did read (dense.f:184,280,86) points to *backdated*
  PTBAA, so the "current dbh + all-dead-exposed" override reads as a match-tuned
  choice until the LSTART DENSE ordering is confirmed.
- Faithfulness impact: if FVS in fact presents backdated PTBAA to the calibration DGF,
  this override mis-predicts the calibration DDS (the PBAL term) for any stand with
  measured growth and a different live/dead mix than snt01, shifting COR/OLDRN.
  Need: the LSTART calibration driver's DENSE/DGF sequence (LBKDEN/LREDO passes and
  any trailing current-dbh DENSE) to settle BANDAID-vs-faithful.
