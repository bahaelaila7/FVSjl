# Confirm B4 — SSTAGE PCTSMX demotion uses user BAMAX keyword instead of computed stand SDImax (BTSDIX)

**Verdict: CONFIRMED**

## The flag
In `structure_class` (single-stratum, sawtimber-DBH branch), the class-2 → class-1
"percent-of-MaxSDI" demotion uses the user `BAMAX` keyword value (`s.control.ba_max`,
which is `0` when the keyword is absent) where FVS uses the **computed** stand
maximum SDI `BTSDIX`.

## Julia logic (`src/engine/structure_stage.jl:164-173`)
```julia
if st.nstr == 1
    if tmpdbh < ssdbh
        cls = 1
    elseif tmpdbh < sawdbh
        cls = 2
        xbamax = Float64(s.control.ba_max)
        (xbamax > 0 && _event_bsdi(s) < 0.01 * pctsmx * xbamax) && (cls = 1)   # PCTSMX demotion
    else
        cls = dmind < 3.0 ? 6 : 5
    end
```
- `s.control.ba_max` is the user **BAMAX** keyword, field 1 (a *basal area* max, sq ft/ac),
  set only in `keyword_dispatch.jl:217` (`rec.present[1] && rec.values[1] > 0 ⇒ ba_max`)
  and `:527`. **It is `0` when the user does not key BAMAX** (the common case, e.g. snt01).
- The `xbamax > 0` guard therefore makes the entire demotion a **no-op whenever BAMAX is
  not supplied**.
- `_event_bsdi(s)` = SDIBC (Reineke before-cut SDI) — correct for the IBA=1 path.
- The field intended to hold BTSDIX, `StandState.before_max_sdi` (`state.jl:388`, commented
  `(BTSDIX)`), is declared but **never assigned anywhere in src/** (verified by grep).

## FVS Fortran logic (`bin/FVSsn_buildDir/sstage.f`)
Driver: `grincr.f:240  CALL SDICAL(0,BTSDIX)` then `IBA=1; CALL SSTAGE(IBA,ICYC,.FALSE.)`.
Inside SSTAGE:
```
   (sstage.f:149-156, non-FFE branch)
      ELSE
        XBAMAX = BTSDIX
        IF(IBA.NE.1 .AND. ONTREM(7).GT.0.) XBAMAX = ATSDIX
      ENDIF
   ...
   (sstage.f:544-552, NSTR.EQ.1, sawtimber branch)
        TMPSCL = 2
        IF(IBA.EQ.1)THEN
           XSDI = SDIBC
        ELSE
           XSDI = SDIAC
        ENDIF
        IF (XSDI .LT. .01*TMPPCT*XBAMAX) TMPSCL = 1
```
- `XBAMAX = BTSDIX` = **"MAXIMUM SDI BEFORE TREATMENT"** (`PLOT.F77:70`), a *stand density
  index*, produced by `SDICAL` as the basal-area-weighted average of the per-species SDImax
  defaults (`sdical.f:6-12`). It is always populated/positive for a non-empty stand — it is
  **not** the user BAMAX keyword.
- `TMPPCT = PCTSMX` (default 30 ⇒ `0.01*30 = 0.30`). The test is "single-stratum sawtimber
  stand whose current SDI is below 30 % of stand MaxSDI ⇒ demote to stand-initiation (class 1)".
  Both sides of the comparison are SDI quantities (SDIBC vs 0.30·SDImax), dimensionally
  consistent.

## Precise semantic divergence
| | FVS | FVSjl |
|---|---|---|
| XBAMAX | `BTSDIX` = SDICAL weighted stand **MaxSDI** (always > 0) | `ba_max` = user **BAMAX** keyword = basal-area max, **0 when absent** |
| Effect, no BAMAX keyword | demotion CAN fire (compares SDIBC vs 0.30·MaxSDI) | demotion **never fires** (`xbamax>0` guard fails) → class stays 2 |
| Effect, BAMAX keyword set | n/a (FVS never uses the BA keyword here) | compares SDIBC (an SDI) vs `0.30·BAMAX` (a **basal area**) → wrong units/scale |

Both arms are wrong: the absent-keyword arm silently disables a real FVS demotion; the
present-keyword arm compares an SDI against a fraction of a basal area.

## Faithful fix (do NOT apply)
Replace the user-keyword read with the computed pre-treatment stand MaxSDI. FVSjl already
has the exact SDICAL weighting in `variants/southern/mortality.jl:35` `stand_sdimax(s)`
(BA-weighted `sp_sdi_def`, == `SDICAL(0,BTSDIX)`):

```julia
    elseif tmpdbh < sawdbh
        cls = 2
        xbamax = stand_sdimax(s)            # BTSDIX = SDICAL pre-treatment MaxSDI (sstage.f:150)
        (_event_bsdi(s) < 0.01 * pctsmx * xbamax) && (cls = 1)   # PCTSMX demotion
```
(The `xbamax > 0` guard becomes redundant — `stand_sdimax` returns ≥ 1 — and may be dropped.)
Optionally populate `StandState.before_max_sdi` from `stand_sdimax` at GRINCR-time and read
that, to mirror the Fortran BTSDIX/ATSDIX (before- vs after-treatment) selection exactly; for
the default `iba=1` call `stand_sdimax(s)` is the faithful value.

## Upstream rank: LEAF
`structure_class.class` feeds only the SSTAGE "Structural statistics" report
(`structure_report`), the `BSCLASS` event variable, and the (binary-blocked) `FVS_StrClass`
DBS table. It does **not** feed growth, mortality, density, regen, or other trees. It becomes
MID only if a user IF/THEN keys a treatment off `BSCLASS` (uncommon, not in the default suite).

## Reachability
The demotion branch requires: single valid stratum (`nstr==1`) **and** 70th-pctile DBHNOM in
`[SSDBH, SAWDBH) = [5, 25)` in **and** SDIBC `< 0.30·MaxSDI` — i.e. a young/sparse single-cohort
stand. The SSTAGE port was validated bit-exact vs Fortran on fire_early and snt01 stand-1, but
those mature stands sit well above 30 % of MaxSDI and/or have multiple strata, so this branch is
**almost certainly never exercised** → the bug is currently **silent**. Fixing it will NOT move
existing goldens unless a scenario actually enters this branch; if one does, the corrected value
flips class 2 → 1 to match Fortran.

## Masked-bug watch
No masked bug. The faithful fix only enables/realigns a demotion the FVS source already
performs; any golden that currently shows class 2 in the qualifying condition is itself wrong
(the Fortran-generated reference would be class 1), so a "regression" there would actually be a
correction. No semantically-correct downstream value is at risk.
