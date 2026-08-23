# AUTOES ingrowth PROB1 residual — #143 follow-on (NSTORE/TPA gap) — REAL-BUG-FIXED

## VERDICT: REAL-BUG-FIXED (3 measured establishment bugs). Gate 339/11 byte-identical.

The post-#143 "cornered NSTORE under-stocked seedling" residual (jl TPA 7–15% low vs the oracle) was
**re-measured against the now-bit-exact draw stream and found to be THREE real, fixable bugs — NOT a numeric
straddle.** NSTORE itself was already correct (3 vs 3 at ic=1); the residual was the ingrowth stocking
**PROB1** (jl 0.4631 vs live 0.6174), which scales every established tree (ESPROB = PROB1·NEWTPP/ITPP).

Fixture: under-stocked IE `753189105290487`, `ie_understocked.db`, plain ESTAB, NUMCYCLE 5,
`InlandEmpire()`. Oracle: `FVSie_estab2` (FVSie_g16 relinked with instrumented estab.f/esplt2.f/initre.f
dumping ESTOCK inputs+PLPROB+NSTORE+ESB/ESB1+PSLO/PASP). All numbers MEASURED, not inferred.

## Root cause — PROB1 = logistic(PN_endcycle + ESB − ESB1), decomposed at ic=1

| quantity            | oracle (FVSie_g16)        | jl BEFORE            | jl AFTER fix        |
|---------------------|---------------------------|----------------------|---------------------|
| small-tree TPA      | 477.73 (PLPROB·DUPNPT)    | 477.73 ✓             | 477.73 ✓            |
| NSTORE (ingrowth)   | 3                         | 3 ✓                  | 3 ✓                 |
| **PSLO / PASP**     | **0.33 / 2.79 rad (160°)**| **0 / 0**  ✗         | 0.33 / 2.79 ✓       |
| PN (end-cycle BAA)  | 0.0425 (BAA=77.51)        | −0.148 (BAA=62.51) ✗ | ~0.046 (BAA=78.4)   |
| ESB1 (inv BAAOLD)   | −0.3583 (BAAOLD=62.48)    | −0.6627 ✗            | −0.358 ✓            |
| ESB                 | 0.0778                    | 0.0778 ✓             | 0.0778 ✓            |
| esb_shift (ESB−ESB1)| +0.4361 (INADV=0)         | 0 (skipped) ✗        | +0.4358 ✓           |
| **PROB1 (FTEMP)**   | **0.6174**                | **0.4631** ✗         | **0.6188** (≈exact) |

### Bug 1 (DOMINANT) — establishment used the STAND slope/aspect, not the per-plot PSLO/PASP
estab.f runs the plot loop on `SLO=PSLO(NNID)`, `XCOSAS=cos(PASP(NNID))` — the **per-plot** slope/aspect,
which FVS reads from the **TREE records** (esplt1.f:69-70, IPINFO=2), NOT from the stand card. This stand's
`FVS_TREEINIT_COND` carries SLOPE=33/ASPECT=160 while `FVS_STANDINIT_COND` has SLOPE=0/ASPECT=0 (the stand
value drives DGF — correctly kept at 0 in jl, hence near-exact BA). At TIME=0 (ESB1) all the SQSQ/SQREGT
terms vanish and ESTOCK's `−0.017901·XCOSAS·SLO·XBAA − 0.006001·XSINAS·SLO·XBAA` interaction is all that
survives BA/ELEV; with SLO=0 jl dropped it (+0.30 to ESB1). jl already reads these per-plot values into
`p.point_slope`/`p.point_aspect` (used for ESTPP) — the fix feeds them into the stocking PN + ESB1 too.

### Bug 2 — ESB inventory calibration skipped on the first-cycle ingrowth
estab.f applies the ESB−ESB1 actual-vs-predicted correction when `INADV=0` (KDT+1−IY(1) ≤ 20) AND NTALLY=1.
The first-cycle AUTOES ingrowth IS a NTALLY=1/INADV=0 tally (dump: `ntally=1 inadv=0 esb=0.0778` at ic=1;
`inadv=1 esb=0` at ic=3/5). jl's `idsdat==inv_year` gate missed it (ingrowth idsdat=−1) ⇒ esb_shift=0. Fix:
gate on `(next_year−inv_year) ≤ 20 && ntally∈{1,99}` (= INADV=0 on a fresh tally).

### Bug 3 — stocking PN used the PRE-growth (stale) per-point BA
estab.f:572 computes the end-of-cycle PN with `BAA=BAAA(NNID)` = the POST-growth per-point BA. jl's density
was last refreshed pre-growth (simulate.jl:548, before diameter/height growth), so `point_ba` lagged one
cycle (jl 62.51 = live INVENTORY BAAOLD 62.48, not the live end-cycle 77.51). Fix: `compute_density!` before
the ingrowth stocking PN. (ESB1 correctly keeps the pre-growth baaa as its inventory BAAOLD.)

## .sum impact (jl vs live FVSie_g16), MEASURED

| year | live TPA | jl BEFORE | jl AFTER |   | live BA | jl AFTER |
|------|----------|-----------|----------|---|---------|----------|
| 2019 | 763      | 763       | 763      |   | 167     | 167      |
| 2029 | 815      | 755 (−7%) | **815** ✓|   | 201     | 203      |
| 2039 | 603      | 552 (−8%) | 599      |   | 222     | 224      |
| 2049 | 728      | 653 (−10%)| 726      |   | 265     | 263      |
| 2059 | 572      | 501 (−12%)| 563      |   | 288     | 287      |
| 2069 | 793      | 676 (−15%)| 782      |   | 339     | 339      |

Was 7–15% low; now within ~1–2% (2029 bit-exact). **Residual = cornered**: jl post-growth point BA 78.4 vs
live end-cycle BAAA 77.5 (~1.2%) — the pre-existing per-point BA scaling/attribution Δ in the density
machinery (establishment.jl:1392), NOT an establishment-logic bug. PROB1 0.6188 vs 0.6174 is this ~1%.

## Validation
- **Gate `test/integration/test_multicycle.jl` = 339 pass / 11 broken, byte-identical** (SN/BARE; IE-only fix).
- IE establishment unit suites unregressed: test_ie_estock, test_estab_minplots 12/12, test_estab_mechprep
  11/11 (disturbance path — the oracle uses per-plot slope there too, so es_slope converges TO the oracle),
  test_estab_specmult_htadj 8/8.
- iet01 AUTOES anchor inert: inline/TREEDATA ⇒ p.point_slope empty ⇒ es_slope falls back to p.slope; the
  ie_autoes_run compose test (12/12, includes iet01 stand-4) passes.

## Fix (src/variants/inlandempire/establishment.jl, ie_autoes_establish!)
1. `es_slope/es_aspect` = per-plot point_slope[1]/point_aspect[1] (fallback stand) → ESB1 + ie_autoes_run.
2. esb_shift gate `idsdat==inv_year` → `(next_year−inv_year)≤20 && ntally∈{1,99}`.
3. `compute_density!` before the ingrowth stocking PN (baaa_pn = post-growth point_ba[1]).
