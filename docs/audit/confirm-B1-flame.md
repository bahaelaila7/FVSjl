# B1 — fire-behavior flame: weighted flame summed per-model vs recomputed from weighted Byram

**Verdict: CONFIRMED** (semantic mismatch, two related divergences)

## The flag
`fmburn.jl` accumulates a *weighted sum of per-model flame lengths* and uses that as the
stand flame, while FVS `fmfint.f` accumulates that same weighted sum **only to throw it
away** and recompute flame from the weighted-total Byram intensity.

## What the FVS Fortran computes

`fire/base/fmfint.f`, inside the `DO 800 INB = 1,MXFMOD` fuel-model loop (line 109):

```fortran
507  IF ((FTYP .NE. 2) .OR. (ICALL .EQ. 2)) THEN
508    FLAME = FLAME + (0.45 * (BYRAMT / 60.0) ** 0.46) * FWT(INB)   ! weighted Σ of per-model flame
509    BYRAM = BYRAM + BYRAMT * FWT(INB)                            ! weighted Σ of per-model byram
...
535  800 CONTINUE
...
C     COMPUTE FLAME LENGTH AS A FUNCTION OF BYRAM RATHER THAN USING
C     THE WEIGHTED AVERAGE AS COMPUTED IN THE BLOCK ABOVE (NLC 21 AUG 2003)
541  FLAME = (0.45 * (BYRAM/ 60.0) ** 0.46)                        ! OVERWRITES the Σ above
```

So the per-model FLAME accumulation at :508 is dead — line **:541 unconditionally overwrites**
`FLAME` with `0.45*(BYRAM/60)^0.46`, where `BYRAM` is the weighted total of per-model Byram.
The returned flame is a function of the *weighted-total Byram*, never the weighted sum of flames.
(The 2003 NLC comment states this explicitly.)

Downstream in `fire/vbase/fmburn.f`:
```fortran
436  CALL FMFINT(IYR, BYRAM, FLAME, 1, HPA, 1)   ! FLAME := 0.45*(BYRAM/60)^0.46
439  OLDFL = FLAME
443  IF (FLMULT .NE. 1.0) FLAME = OLDFL * FLMULT  ! flame multiplier
459  IF (FLAME .NE. OLDFL) THEN
464    BYRAM = 60.0 * ((FLAME / 0.45) ** (1.0 / 0.46))  ! BACK-COMPUTE byram from multiplied flame
465  ENDIF
470  BYRAM = BYRAM / 60.0
471  SCH = (63.0/(140.0-ATEMP)) * (BYRAM**(7/6) / (BYRAM + FWIND**3)**0.5)  ! scorch from that byram
```

## What the Julia computes

`src/engine/fire/fmburn.jl:58-67` (the live fire, `fmburn!`):
```julia
flame_raw = 0f0; byram = 0f0
for (fm, w) in models
    ...
    flame_raw += r.flame * w      # the weighted Σ that Fortran DISCARDS
    byram     += r.byram * w
end
flame = flame_raw * flmult         # uses the discarded sum, NOT 0.45*(byram/60)^0.46
sch = byram > 0f0 ? scorch_height(byram, atemp, fwind) : 0f0
```
Same pattern in `potential_fire` scenario, `fmburn.jl:299-305` (no flmult there).

## Precise semantic divergence

Let weights `wᵢ` (Σ=1) over models with per-model Byram `bᵢ`.

1. **Flame (the flag).** Fortran: `flame = 0.45·((Σ wᵢ·bᵢ)/60)^0.46`.
   Julia: `flame = Σ wᵢ·(0.45·(bᵢ/60)^0.46)`.
   Because `x^0.46` is **concave**, by Jensen `Σ wᵢ·f(bᵢ) ≤ f(Σ wᵢ·bᵢ)`, i.e. Julia's flame is
   a **systematic under-estimate** vs FVS whenever the stand blends ≥2 fuel models with
   fractional weights. (Identical only when a single model has weight 1.) This is exactly the
   kind of residual seen in the F4 note: snt01 stand-4 flame **3.84 (Julia) vs 3.9 (Fortran)**.

2. **Scorch byram when flmult≠1 (secondary).** Fortran back-computes
   `byram = 60·(flame/0.45)^(1/0.46)` from the *multiplied* flame (fmburn.f:464) and feeds that
   to scorch. Julia always feeds the raw weighted `byram` to `scorch_height` regardless of
   `flmult`. So with FLAMEADJ active the scorch heights diverge too. (With flmult=1 the scorch
   byram matches, since both use the weighted byram and Fortran's back-compute is a no-op.)

## Faithful fix (do NOT apply)

In `fmburn.jl:58-67` drop `flame_raw`; compute flame from the weighted byram and back-compute
byram for scorch when a multiplier is applied:
```julia
byram = 0f0
for (fm, w) in models
    load, sav, depth, mext = standard_fuel_model(coef, fm)
    r = rothermel_surface_fire(load, sav, depth, mext, mois; wind = fwind, slope_tan = s.plot.slope)
    byram += r.byram * w
end
flame  = 0.45f0 * (byram / 60f0) ^ 0.46f0     # fmfint.f:541
oldfl  = flame
flame  = flame * flmult                        # fmburn.f:443
if flame != oldfl                              # fmburn.f:459-464
    byram = 60f0 * (flame / 0.45f0) ^ (1f0 / 0.46f0)
end
sch = byram > 0f0 ? scorch_height(byram, atemp, fwind) : 0f0
```
In `potential_fire` (`fmburn.jl:299-305`) make the same primary change (no flmult/back-compute):
```julia
byram = 0f0
for (fm, w) in models
    ...
    byram += r.byram * w
end
flame = byram > 0f0 ? 0.45f0 * (byram / 60f0) ^ 0.46f0 : 0f0
```
(Guard `0^0.46 = 0` is fine in Float32; the explicit `byram>0` keeps it tidy.)

## Upstream rank: MID
`flame` feeds `fire_tree_mortality` → `pmort` → `curkil` → `t.tpa` (and snag pool), i.e. it
changes killed TPA which propagates into later cycles — so it is upstream of stand trajectory.
`byram`→scorch→`crown_volume_scorched`→`pmort` too. It is downstream of the fuel-model weighting.
Consumers: `fire_tree_mortality`, `scorch_height`, the `burn_reports` flame/scorch columns,
`FireResult.flame/byram/scorch`. Rank MID (event-internal but it does move downstream TPA).

## Reachability
Exercised: the FFE weighted-fuel-model path is the SN active surface-fire path
(FMCFMD/FMDYN), reached by snt01 stand-4 and any FLAMEADJ/burn scenario in the suite. The
divergence is nonzero precisely when ≥2 standard models are blended, which is the normal SN
case. The 3.84-vs-3.9 flame and the BA-81-vs-78 kill residuals recorded in the F4/fuel-model
notes are consistent with this under-estimate. A faithful fix should **move** the suite
(fire-event flame, scorch, killed TPA/BA, and post-fire cycles for fire stands).

## Masked-bug watch
The fix raises flame (Jensen) → higher `pmort` → more fire kill → higher killed BA. The
recorded snt01 stand-4 residual is Julia killing **less** than Fortran (BA 78 vs 81), so this
fix moves toward FVS, not away — no masked regression is expected here; rather it should
shrink an existing residual. If after the fix some report value regressed instead, that would
signal a separate masked bug (e.g. in the per-tree kill distribution), not a reason to reject
this change.
