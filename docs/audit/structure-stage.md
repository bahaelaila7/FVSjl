# Faithfulness audit — Stand structure stage (SSTAGE)

Module: `src/engine/structure_stage.jl`
FVS sources checked: `base/sstage.f` (SSTAGE / SSTGHP / SSTGHTPA / UPDATECCCOEF),
`base/covolp.f`, `base/pctile.f`, `bin/FVSsn_buildDir/evtstv.f`, `bin/FVSsn_buildDir/grincr.f`,
`bin/FVSsn_buildDir/sdical.f`, `bin/FVSsn_buildDir/PLOT.F77`.

The core machinery is a faithful port: the height-gap stratification (two-largest-gap
tracking, SUMPRB<2 ladder-tree skipping, 10-ft min gap, upper-stratum swap), the
cover-includes-gap `I2=MAX(IS1I2,IS2I1-1)` window, COVOLP's `(1-exp(-cccoef·ΣA/43560))·100`
with raw PROB crown area, the 41382 (=0.95 ac) cohort cutoff, the RDPSRT(.FALSE.)
descending single-crown sort, the PCTILE 70th-percentile ±4-tree PROB-weighted DBHNOM/height,
the ICRB crown-base average, the species top-2-by-cover pick, the `NSTR==0 & TPROB>=TPAMIN`
reform-one-stratum path, and the FORMAT 85/90 report layout all match the Fortran. Those are
counted in `decisionsReviewed`, not listed. Three source-grounded concerns below.

---

## FLAG 1 — BANDAID: PCTSMX (SE→SI) demotion uses the user BAMAX keyword, not the computed stand SDImax (BTSDIX)

- jl symbol/line: `structure_class`, lines 169-170:
  ```julia
  xbamax = Float64(s.control.ba_max)
  (xbamax > 0 && _event_bsdi(s) < 0.01 * pctsmx * xbamax) && (cls = 1)   # PCTSMX demotion
  ```
- Claim: this reproduces the FVS class-2 (SE) → class-1 (SI) demotion when current SDI is
  below PCTSMX% of the stand max.
- FVS source checked:
  - `sstage.f:154` (non-FFE branch): `XBAMAX = BTSDIX`.
  - `sstage.f:544-550`: `IF(IBA.EQ.1) XSDI=SDIBC ELSE XSDI=SDIAC; IF (XSDI .LT. .01*TMPPCT*XBAMAX) TMPSCL = 1`.
  - `grincr.f:240`: `CALL SDICAL(0,BTSDIX)` — BTSDIX is recomputed every cycle.
  - `sdical.f:6-7`: SDICAL "COMPUTES THE MAXIMUM SDI IN EFFECT FOR A STAND … a weighted
    average of the SDI maximums by species" — i.e. BTSDIX is a Reineke **SDImax** (carrying
    capacity), always available, in SDI units.
  - `PLOT.F77:70`: `BTSDIX -- MAXIMUM SDI BEFORE TREATMENT`.
- What FVS actually does vs jl: `XBAMAX` is `BTSDIX`, the per-cycle computed stand **SDImax**
  (Reineke units). The jl substitutes `s.control.ba_max`, which `keyword_dispatch.jl:210-217`
  documents as the user **BAMAX keyword** (a basal area, LBAMAX, **0 when the keyword is
  absent**). Two consequences:
  1. When no BAMAX keyword is supplied (the normal case, e.g. snt01), `ba_max == 0`, so the
     `xbamax > 0` guard makes the demotion **never fire** — yet Fortran's BTSDIX is the
     always-computed SDImax, so `SDIBC < 0.30·SDImax` *can* demote SE→SI.
  2. Even with a BAMAX keyword present, the test compares an SDI (`_event_bsdi`=SDIBC) against
     `0.30·basal-area`, a dimensionally different quantity than `0.30·SDImax`.
  There is no source basis for using the user BAMAX here; FVS mandates the SDICAL SDImax.
- Severity: BANDAID.
- Faithfulness impact: any young, understocked stand whose dominant stratum DBH is in
  [SSDBH, SAWDBH) and whose SDI is below 30% of the stand SDImax is classified SE (2) by jl
  but SI (1) by FVS, unless a BAMAX keyword happens to be present with a coincidentally-matching
  magnitude. Passing tests likely never exercise the demotion (class lands at 2 either way),
  which is exactly why the wrong variable survives.

---

## FLAG 2 — GAP: the single-tree (NTREES≤1) classification path is not ported

- jl symbol/line: `_ss_strata` / `structure_class` — n==1 falls through the **general**
  stratification path (no special branch).
- FVS source checked: `sstage.f:235-267` (`IF (NTREES.LE.1) THEN … GOTO 80`). For a single
  tree FVS uses a distinct classifier:
  - cover test is **linear**: `IF (WK6(I) .LT. 435.60*TMPCCM)` (crown-area fraction < CCMIN%),
    not COVOLP's exponential;
  - below that cover it sets `TMPSCL=0`, then `TMPSCL=1` only if `PROB ≥ TPAMIN`;
  - it can yield class 1/2/5 and uses `SDIAC` for the PCTSMX test (`sstage.f:252`); it **never**
    produces class 6 (no DMIND<3 branch).
- What jl does instead: routes n==1 through `_ss_strata`+`structure_class`, using the
  exponential COVOLP cover for the OK test and the multi-stratum classifier (which *can* reach
  class 6 via `dmind<3`). Near the CCMIN boundary the exponential vs linear cover diverges
  (e.g. pccu=0.05 → 4.88% exp vs 5.0% linear), and the below-cover "class=1 if TPA≥TPAMIN
  regardless of DBH" rule is lost.
- Severity: GAP (degenerate single-record stand; untested).
- Faithfulness impact: a 1-record stand can get a different OK/class result; immaterial for the
  multi-tree scenarios actually validated.

Sub-note (same flag, filtering): jl `_ss_strata:90` keeps only `height>0 && tpa>0` trees,
whereas `sstage.f:222` filters solely on `PROB > 0.00001` (HT-zero records are kept and sink to
the bottom stratum via RDPSRT). Immaterial when all live records have HT>0, but it is a literal
divergence from the source loop.

---

## FLAG 3 — GAP: after-thin (Rm=1) row reuses before-thin SDI / cover; CCCOEF2 + SDIAC paths unported

- jl symbol/line: `write_structure_report:266` emits the Rm=1 row via
  `structure_report_row(stand, yr, 1)`, which calls `structure_report → structure_class(s)`
  with the default `iba=1`; comment at line 266 says "same w/o a thin". `structure_class`
  always uses `_event_bsdi(s)` (=SDIBC) for the PCTSMX SDI and never branches on `iba`.
- FVS source checked:
  - `sstage.f:546-548`: after thinning (`IBA≠1`) the demotion uses `XSDI = SDIAC` (SDI **after**
    cutting), not SDIBC.
  - `sstage.f:154-156`: `XBAMAX = BTSDIX`, but `IF(IBA.NE.1 .AND. ONTREM(7).GT.0.) XBAMAX=ATSDIX`
    (after-thin SDImax when there were removals).
  - `sstage.f:98` + `UPDATECCCOEF` (`sstage.f:904-941`): a CCADJ (activity 444) sets `CCCOEF2`,
    and `IF(INICYCLE.GT.1.AND.INBA.EQ.2.AND.CCCOEF2.NE.1) CCCOEF=CCCOEF2` — the after-thin cover
    coefficient. jl has no CCCOEF2/CCADJ handling (`cc_coef` is a single value).
- What jl does instead: the Rm=1 row is byte-identical to Rm=0 for an unthinned cycle (correct
  there), but for a cycle with an actual removal jl would still report before-thin SDI/cover and
  never apply the after-thin SDIAC/ATSDIX/CCCOEF2 substitutions, so the after-thin structural
  class/cover can be wrong.
- Severity: GAP (the validated stands report on unthinned/no-CCADJ cycles, so this never fires).
- Faithfulness impact: thinned-cycle after-thin (Rm=1) structural class/cover may diverge from
  FVS; no impact on the snt01-stand-1 / fire_early validation runs.

---

## Verdict
Core SSTAGE classification + per-stratum report machinery is faithful. One real BANDAID
(the PCTSMX SE→SI demotion keyed to the user BAMAX keyword instead of the SDICAL stand
SDImax/BTSDIX, guarded so it cannot fire on a plain run) and two untested-path GAPs
(single-tree classifier; after-thin SDIAC/ATSDIX/CCCOEF2).
