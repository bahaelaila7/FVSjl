# B2 — sprout crown-width CR arg: `70f0` vs FVS `CRDUM=1.0`  →  CONFIRMED

## Flag
`src/engine/sprout.jl:171-172` calls `crown_width(coef, sp2, dbh, ht, 70f0, 0, …)`,
passing **70** in the crown-ratio (`cr`) slot. FVS `estb/esuckr.f` passes a dummy
`CRDUM = 1.0` in the corresponding `CWCALC` argument and passes `ICR=70` only in a
*discarded* slot.

## FVS Fortran logic
`estb/esuckr.f` (SN build: `bin/FVSsn_buildDir/esuckr.f`), lines 313-320:
```fortran
      ICR(ITRN)=70
C  CALCULATE A CROWN WIDTH FOR SPROUTS
      CRDUM=1.
      CALL CWCALC(ISSP,PROB(ITRN),DBH(ITRN),HT(ITRN),CRDUM,
     &            ICR(ITRN),CW,0,JOSTND)
      CRWDTH(ITRN)=CW
```
`sn/cwcalc.f` (`bin/FVSsn_buildDir/cwcalc.f`) signature, line 1:
```fortran
      SUBROUTINE CWCALC(ISPC,P,D,H,CR,IICR,CW,IWHO,JOSTND)
```
- 5th arg `CR` = "CROWN RATIO IN PERCENT" (cwcalc.f:19) — **used** in the equations.
- 6th arg `IICR` = "CROWN RATIO - FLOAT" (cwcalc.f:20) — **discarded**:
  cwcalc.f:82 `IDANUW = IICR` is the dummy-arg-not-used suppression; `IICR`
  never appears again.

So in FVS the sprout crown-width equation receives **CR = 1.0**. The `ICR=70`
the sprout sets is consumed elsewhere (it becomes the tree-record crown percent),
but it is NOT what feeds the CW equation. The Bechtold family equations use `CR`
linearly, e.g. cwcalc.f:667 (balsam fir, eqn 01201):
```fortran
           CW = 0.6564 + 0.8403*D + 0.0792*CR
```

## FVSjl logic
`src/engine/crown_width.jl:60` signature:
```julia
crown_width(coef, sp2, d, h, cr, iwho, lat, long, elev)
```
The 5th positional arg is `cr`. For the Bechtold family it is used directly
(`crown_width.jl:39-41`):
```julia
    if fam === :bechtold
        xc = min(x, e.dbh_cap)
        v = e.a + e.b * xc + e.c * xc * xc + e.cr_coef * cr + e.hi_coef * hi
```
`sprout.jl:171` passes `70f0` into that `cr` slot. So the Julia CW equation
receives **cr = 70**, not 1.0.

## Semantic difference
For any sprout whose species uses a **Bechtold-family** CW equation (the only
family with a `cr_coef` term), the crown width differs by
`cr_coef * (70 − 1) = 69 * cr_coef`.
Balsam-fir example: FVS `0.0792*1 = 0.0792`; Julia `0.0792*70 = 5.544` → CW is
~5.46 ft larger in Julia. The other families (`:bragg`, `:braggm`, `:ek`,
`:smith`) ignore `cr`, so those species are unaffected — but any Bechtold-coded
sprout species is materially wrong.

This is the classic "discarded arg vs used arg" transcription slip: the Julia
port mapped FVS's *discarded* `IICR=70` onto the *used* `cr` parameter, and
dropped FVS's actual `CRDUM=1.0`.

## Faithful fix (do NOT apply here)
`src/engine/sprout.jl:171` — change the `cr` argument from `70f0` to `1f0`:
```julia
            cw = crown_width(coef, sp2, dbh, ht, 1f0, 0,
                             s.plot.latitude, s.plot.longitude, s.plot.elevation)
```
Leave `t.crown_pct[n] = Int32(70)` / `t.crown_ratio[n] = 70f0` (sprout.jl:180-181)
untouched — those correctly carry the FVS `ICR=70`. Only the CW-equation `cr`
argument must become the FVS dummy `1.0`.

## Upstream rank: LEAF
The CW value is stored only in `t.crown_width[n]` (sprout.jl:182). Tracing
consumers of the *stored* field across `src/`: the only reader is
`src/io/dbs_output.jl:466` (TreeList DBS `CrWidth` column). `compute_density!`
(simulate.jl:110-116) does NOT read `t.crown_width` — CCF/canopy paths
(standstats.jl:189, cuts.jl, structure_stage.jl) **recompute** crown width on the
fly from DBH/HT and their own CR, so the stored sprout CW never feeds density,
growth, mortality, or later cycles. Consumer = TreeList report column only ⇒ LEAF.

## Reachability
Requires `LSPRUT && cut_log` (sprout.jl:140) — does not fire for the standard
snt01 stands (LSPRUT forced false in esinit). It IS exercised by the dedicated
SPROUT/ESUCKR suite (`sprout.key`, the 113 sprout tests). Fixing it should move
only the TreeList `CrWidth` column for sprout records of Bechtold-family species;
it cannot move any density/growth/summary number. If those sprout tests currently
validate the CW column against live Fortran and pass, that is suspicious (either
the column is not compared, or the test sprout species are non-Bechtold) — worth a
live `CWCALC` differential to confirm the reported column.

## Masked-bug watch
None expected: the fix is isolated to a report column for a narrow path. If a
sprout golden for a Bechtold species currently "passes" with cr=70, the fix will
flip that column; per the validate-against-live-FVS principle, a flip there is the
*correct* value (FVS uses 1.0) and indicates a stale oracle, not a regression.
