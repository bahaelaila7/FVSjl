# Bandaid Audit — Crown ratio (CROWN) + init crown (CRATET)

Module: `src/variants/southern/crown_ratio.jl`
FVS sources checked: `sn/crown.f`, `sn/cratet.f`, `base/dense.f`, `base/fvs.f`, `base/grincr.f`.

## Summary

The core CROWN model is faithful. I traced and confirmed against FVS source:
the ISORT rank direction (`crown.f:156-159`/`261`), the percentile X formula and
`[0.05,0.95]` clamp (`crown.f:288-298`), SCALE = `1-0.00167*(RELDEN-100)` clamped
`[0.30,1.0]` (`crown.f:285-287`), RELSDI = `SDIAC/SDIDEF*10` clamped `[1,12]` with
default 6 (`crown.f:176-182`), all five MCR equations and their MCREQN(2/3/4/5/6)
coefficient roles (`crown.f:190-242`), the Weibull `B=b0+b1*acrnew` / `B≥3` / `C≥2`
clamps and the `A+B*(-ln(1-x))^(1/C)` draw (`crown.f:248-299`), the negative-ICR
sign-flip bypass (`crown.f:270-271`), the ±1%/yr change limit (`crown.f:310-314`),
the CRNMULT-on-change window (`crown.f:318`), the crown-length cap geometry
(`crown.f:333-337`), and the final 95/10/1 clamps (`crown.f:366-368`).

I also confirmed the timing is right: per-cycle CROWN gets the **pre-growth** Reineke
SDIAC (`fvs.f:196` SDICLS before CRATET for init; `grincr.f:323` for the cycle, set
before growth and not recomputed in GRADD) — passed via `crown_sdi` — and the
**post-growth** RELDEN (GRADD's DENSE-before-CROWN) via the default `stand_ccf`. The
init path's backdated-CCF override (`bd_relden`) correctly reproduces DENSE's LBKDEN
two-pass result (`dense.f:91-132,257-264`, RELDEN ends = backdated RELDM1). The
IDG missing-exclusion (`g<0` for IDG 1/3, `g≤0` for IDG 0/2) is consistent with
`cratet.f:174-179` (all `DG≤0` inputs become `-1`) + `dense.f:100`, and the bark-divide
skip for IDG==1 matches `dense.f:101,122`. These are all faithful.

The flags below are all narrow GAPs that are correct for the default/tested path but
silently diverge from FVS on untested keyword/inventory edge cases.

---

## GAP 1 — `dbh ≤ 0` uses a fixed 0.5 instead of a RANN draw

- jl: `crown_ratio_update!`, line 125:
  `x = d > 0f0 ? Float32(isort[i]) / Float32(n) * scale : 0.5f0 * scale`
- Claim: silently substitutes a constant 0.5 percentile for a zero-diameter tree.
- FVS `sn/crown.f:288-293`:
  ```
  IF(DBH(I) .GT. 0.0) THEN
     X = (REAL(ISORT(I)) / REAL(ITRN)) * SCALE
  ELSE
     CALL RANN(RNUMB)
     X = RNUMB * SCALE
  ENDIF
  ```
  FVS draws a **random** percentile (and consumes one RANN value from the stream).
- Severity: GAP. Faithful only because every live record at CROWN time has dbh>0 in
  the tested SN scenarios. If a dbh≤0 live record ever reaches CROWN, FVSjl diverges
  both in value (0.5 vs random) AND in the RNG stream (FVSjl skips the RANN consume),
  which would shift every subsequent stochastic draw that cycle.

## GAP 2 — CRNMULT keyword ignored on dubbed crowns and on the <10 floor

- jl: `crown_ratio_update!`, lines 138, 143, 146-147. The `icr_old == 0` (missing/dubbed)
  branch sets `icri = trunc(crnew+0.5)` with **no** CRNMULT scaling, and the `<10`
  floor / crmax-refloor are applied unconditionally.
- FVS `sn/crown.f:323-328` applies CRNMULT to the dubbed crown itself:
  ```
  9052 ICRI = INT(CRNEW(I)+0.5)
       IF(LSTART .OR. ICR(I).EQ.0)THEN
         IF(DBH(I).GE.DLOW(ISPC) .AND. DBH(I).LE.DHI(ISPC))THEN
           ICRI = INT(REAL(ICRI) * CRNMLT(ISPC))
  ```
  and the crmax/`<10` floors are gated on a unit multiplier:
  `crown.f:346` `IF(ICRI.LT.10 .AND. CRNMLT(ISPC).EQ.1.0)ICRI=INT(CRMAX+0.5)` and
  `crown.f:367` `IF (ICRI .LT. 10 .AND. CRNMLT(ISPC).EQ.1) ICRI=10`. FVSjl drops the
  `CRNMLT==1.0` guard.
- Severity: GAP. `active_crn_mult` correctly handles the change-on-existing-crown path
  (`crown.f:318`), so with no CRNMULT keyword (CRNMLT=1.0, the common case) behavior is
  identical. But a CRNMULT keyword with a non-unit multiplier would (a) not scale newly
  dubbed crowns of missing-crown / regen trees, and (b) still apply the <10→10 floor that
  FVS suppresses when CRNMLT≠1. Untested CRNMULT-keyword path only.

## GAP 3 — Inventory top-killed crown reduction (ITRUNC) not ported

- jl: `init_crown_ratios!` / `crown_ratio_update!` — no ITRUNC/NORMHT-based crown
  reduction anywhere in `crown_ratio.jl`.
- FVS `sn/crown.f:350-354` (LSTART branch, live tree loop):
  ```
  55 IF (.NOT.LSTART .OR. ITRUNC(I).EQ.0) GO TO 59
     HN=REAL(NORMHT(I))/100.0
     HD=HN-REAL(ITRUNC(I))/100.0
     CL=(REAL(ICRI)/100.)*HN-HD
     ICRI=INT((CL*100./HN)+.5)
  ```
  At inventory, a live top-killed tree (ITRUNC≠0) with a dubbed crown has that crown
  reduced for the lost top. Only missing-crown records reach here (`crown.f:263` skips
  ICR>0 at LSTART), which matches FVSjl's "only crown_pct==0" init scope — so the case is
  reachable, just unhandled.
- Severity: GAP. Affects only inventory trees that are BOTH flagged top-killed (broken/
  dead top, NORMHT<0 input → ITRUNC set in `cratet.f:380-396`) AND have no input crown.
  Their cycle-0 crown (and downstream DGF/mortality/volume) would be too high vs FVS.

## UNVERIFIED — dead/snag-record crown dubbing (DUBSCR) for the inventory snag list

- jl: `init_crown_ratios!` operates only on the live tree list (`s.trees`). FVS `CROWN`
  also dubs crowns on cycle-0 dead/snag records `IREC2..MAXTRE` via DUBSCR with the same
  ITRUNC reduction: `crown.f:74,380-396`.
- I could not locate a DUBSCR-equivalent for inventory snag-record crowns in FVSjl
  (`grep dubscr` → none; FFE carbon/snag code derives snag biomass from bole volume, not
  an ICR). This may be intentional (SN FFE snag pools may not need a per-record ICR) or a
  real omission feeding the FFE StandDead/CWD2B crown pools.
- Severity: UNVERIFIED. To settle: trace whether any SN FFE inventory-snag carbon path
  consumes a per-snag crown ratio; if so, the DUBSCR dub at `crown.f:380-396` is missing.
