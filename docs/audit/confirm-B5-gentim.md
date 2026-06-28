# B5 — establishment GENTIM: `yr − idsdat` vs FVS `max(FINT−5, 0)`

**Verdict: CONFIRMED** (semantic mismatch; masked in the current suite by the XMIN floor and by the TREEHT-specified branch).

## The two implementations

### Julia (src/engine/establishment.jl)
Line 61:
```julia
gentim = (Int(yr) + per - idsdat) - per; gentim < 0 && (gentim = 0)
```
The `+per` and `−per` cancel, so this is algebraically:
```
gentim = yr − idsdat            (floored at 0)
```
where `yr = current_cycle_year(s)` (calendar year, line 45) and `idsdat = _es_idsdat(s)` (line 51) = the TALLY(427) disturbance date, else `cycle_year − 20` (lines 26–31).

`gentim` then feeds the establishment AGE (line 81):
```julia
age = Float32(per) - Float32(delay) - Float32(gentim) + trage; age < 1f0 && (age = 1f0)
```
and `age` drives the established height `hht = htcalc_height(bc, sp, si, age, montane)` (line 83), used only on the **default-height (natural / no-TREEHT)** branch (lines 91–96).

### FVS Fortran (estb/estab.f)
The GENTIM that is passed into ESSUBH (and therefore governs the planted/natural tree HEIGHT) is computed once, at lines 448–449:
```fortran
      GENTIM=FINT-5.0
      IF(GENTIM.LT.0.0) GENTIM=0.0
```
i.e. `GENTIM = max(FINT − 5, 0)` — a pure function of the cycle length FINT. It is handed to the PLANT/NATURAL height call unchanged:
```fortran
1023  CALL ESSUBH (IPNSPE,HHT,EMSQR,DILATE,DELAY,ELEV,IHTSER,GENTIM,TRAGE)
```
Inside ESSUBH (estb/essubh.f), with `TIME` set to `FINT` (estab.f:1020 `TIME=FINT`):
```fortran
89    AGE=TIME-DELAY-GENTIM      ! = FINT − DELAY − GENTIM
92    AGE=AGE+TRAGE              ! AGE = FINT − DELAY − GENTIM + TRAGE
```
So FVS's AGE formula = `FINT − DELAY − GENTIM + TRAGE`, with `GENTIM = max(FINT−5, 0)`. The Julia mirrors the *shape* of this formula (`per − delay − gentim + trage`) but supplies the wrong `gentim`.

### The second GENTIM (estab.f:1055–1059) is NOT the height GENTIM
After ESSUBH returns, FVS re-derives GENTIM for the planted record:
```fortran
1055  IF(FINT-DELAY .LT. 5)THEN
1056    GENTIM = 0.
1058    GENTIM = FINT-DELAY-5.0
```
This second value is used only for the stocking multiplier `HTIMLT`/`STOMLT` (lines 1061–1063), which the Julia port does not compute. It is **also** purely a function of `FINT` and `DELAY`. Neither GENTIM form anywhere in estab.f depends on `IDSDAT` or any calendar year — confirming the flag.

## Precise semantic difference
- FVS height-path GENTIM = `max(FINT − 5, 0)` → for SN 5-yr cycles = **0**; for 10-yr cycles = **5**. Independent of calendar year / disturbance date.
- Julia gentim = `yr − idsdat` (floored 0). With the `_es_idsdat` fallback of `cycle_year − 20`, this is ≈ **20** (or whatever calendar gap to the TALLY-427 date), i.e. it scales with the calendar offset, not the cycle length.

These produce different `age` values and hence different `htcalc_height` outputs on the natural/default-height branch.

## Why it is currently masked
1. **TREEHT branch bypass (lines 85–90):** when the PLANT keyword specifies a height (`a.params[5] = TREEHT ≥ 0.1`), the height comes from the lognormal draw around TREEHT and `gentim`/`age` are never used. (FVS: estab.f:1026–1034 same bypass.)
2. **XMIN floor (line 97 / FVS estab.f:1037):**
   ```julia
   hht < es_xmin[sp] && (hht = es_xmin[sp])
   ```
   On the default-height branch, the established saplings are short (TRAGE 2–10 yr) and the predicted height typically lands at or below the per-species establishment minimum XMIN, so both the wrong-age Julia height and the FVS height collapse to the same XMIN value. With the Julia age driven down to its `age < 1 → 1` floor (because `yr−idsdat ≈ 20` makes `age = 5 − delay − 20 + trage` strongly negative), the result is an even shorter raw height that is then floored to XMIN — identical to FVS once floored.

This is exactly why snt01 BARE establishment is bit-exact on TPA (TPA does not depend on height at all) and on DBH whenever the heights floor to XMIN.

## Faithful fix (do NOT apply)
Replace line 61 with the cycle-length form (height-path GENTIM):
```julia
gentim = round(Int, fint) - 5; gentim < 0 && (gentim = 0)
```
(equivalently `gentim = max(per - 5, 0)`). Keep `idsdat` for `_es_idsdat`'s other legitimate uses; just remove it from `gentim`. Note FVS uses real `FINT` (not the rounded `per`) inside ESSUBH, so for non-integer cycle lengths use `gentim = max(fint - 5f0, 0f0)` as a Float to be fully faithful; for the integer SN cycles the two agree.

## Upstream rank: UPSTREAM (currently masked)
`gentim → age → htcalc_height → hht`, and `hht` feeds `_htdbh_dbh` (line 102) → the established tree's DBH and crown → the tree list → `compute_density!` (BA/SDI, line 153) → next-cycle growth and mortality and every later cycle. So the computation is genuinely upstream of growth/density. Consumers: `htcalc_height`, `_htdbh_dbh`, then density/growth. In practice it is masked today (XMIN floor / TREEHT bypass), so the live downstream effect is presently nil for the tested scenarios.

## Reachability
ESTAB is reached by BARE-stand regen (snt01-like NATURAL/PLANT tail). The path that would expose the bug is the **default-height (no TREEHT) branch whose predicted height exceeds XMIN** — not exercised in a way that escapes the XMIN floor in the current keyword-coverage suite. Fixing `gentim` is therefore expected to be a **silent change** (no suite movement) unless/until a scenario plants a default-height natural tree tall enough to clear XMIN. The fix is correctness-restoring regardless and should be made most-upstream-first.

## Masked-bug watch
Because the golden oracle is live FVS (which uses `max(FINT−5,0)`), the faithful fix can only move Julia *toward* FVS. If any future test regresses after this fix, the regression would be from one of the still-unverified neighbours in the same AGE formula — chiefly the `delay` term (Julia line 79 `delay = a.year − yr` vs FVS estab.f:991 `DELAY = IPYEAR − (KDT+1−FINT)`) and the `TIME=FINT` vs `per` integerization — not a reason to reject the gentim fix. Flag those as the next things to differential-test if a 10-yr-cycle PLANT scenario shows a height divergence.
