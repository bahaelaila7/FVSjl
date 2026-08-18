# LPMPB LPOPDY (population-dynamics) — port handoff (measurement foundation READY)

## What LPOPDY is
The alternate MPB path: POPDYN keyword (mpbin.f opt 26 → LPOPDY=.TRUE.) routes MPBCUP → MPBDRV
(233 ln) → MPBMOD (808 ln) + GARBEL classifier / SURFCE / phloem model / beetle brood dynamics
(~1050 ln total, stochastic). The DEFAULT (Cole rate-of-loss, COLDRV) is DONE + bit-exact; jl's
mpb_apply! early-returns when `lpopdy` (byte-identical to not wiring). LPOPDY is the LAST major
LPMPB payload — large, multi-turn.

## Measurement foundation (BUILT this turn)
- Oracle: `/workspace/.iework/lpmpb/FVSie_lpmpb` — REBUILT to a persistent path
  (scratchpad/lpmpb/build_ie_lpmpb_persist.sh; IE buildDir .o − exmpb stub + real lpmpb/*.o + shim).
  Invoke: `FVSie_lpmpb --keywordfile=<x>.key` (companion `<x>.tre` auto-opened, MUST match key base).
- POPDYN stand: `scratchpad/lpmpb/lp_popdy.key` (= lp_on.key + POPDYN in the MPB block) + lp.tre
  (pure-lodgepole IE stand). Copy lp.tre→lp_popdy.tre in the run dir.
- **TARGET measured**: LPOPDY collapses the stand **1990→2000 TPA 89→0** (catastrophic outbreak) vs
  the Cole default's gradual 89→34→31→28. jl currently gives the Cole trajectory (LPOPDY early-return),
  so the port must reproduce the 89→0 collapse. (Verify it's a clean STOP, not a crash, before trusting.)

## Port plan (dep-ordered, multi-turn — same insect recipe as DFB/DFTM/WPBR)
1. Instrument MPBDRV/MPBMOD (single-.o swap into the ieobj set, instrumented .sum byte-identical first)
   → dump the per-class/per-tree population state (SURFCE, phloem, brood, emergence, GARBEL classes).
2. Port MPBDRV driver + MPBMOD core + GARBEL/SURFCE/phloem/brood (glibc transcendentals; own MINSTD
   LCG seed 55329 for the stochastic draws, never FFI). Wire the LPOPDY branch in mpb_apply! (drop the
   early-return). MEASURE BOTH SIDES per cycle (isolated dump-replay misses glue/feeder bugs).
3. Validate end-to-end lp_popdy .sum-DELTA vs FVSie_lpmpb (cornered per #206 where stochastic).
   multicycle 339/11 must hold (LPOPDY only fires under a POPDYN keyword ⇒ additive/inert otherwise).

## PROGRESS 2026-08-18 (sub-chunk 1 of the port)
- **Phloem model — BIT-EXACT 18/18** (scratchpad/lpmpb/validate_phloem.jl vs FVSie_lpmpb_ph fort.779).
  MPBDRV lines 106-116: DDS5=(2·DBH·DG+DG²)/2; BAI5=DDS5·0.7853982; XPT=exp(-3.17152+.12591·ln(BAI5)+
  .50932·ln(DBH)-.0077·HT) for BAI5>1e-4 (else 0). glibc expf/logf.
- MPBDRV FLOW (the LPOPDY driver): (1) MPGR if MPBYR==0 (growth ratio); (2) phloem XPT per LP tree [DONE];
  (3) GARBEL classify LP trees into NACLAS classes by PROB (garbel.f); (4) SURFCE (surface model); (5) MPBMOD
  (808 ln stochastic brood-dynamics core — the big one); (6) per-class SURVIV=CLASS(I,IMPROB)/SURVIV_pre;
  (7) WK2(I)=max(WK2(I), PROB(I)·(1-SURVIV(class))), cap PROB-1e-6 = the mortality output (DFB-style max-combine).
- Instrument recipe: `bash build_ie_lpmpb_persist.sh <instr-dir> <out>` single-.o swaps into the ieobj set;
  instrumented .sum VERIFIED byte-identical to clean (fort.779 phloem dump inert). Run with
  `FVSie_lpmpb --keywordfile=lp_popdy.key` in a dir with lp_popdy.tre.
## PROGRESS 2026-08-18 (sub-chunk 2 — GARBEL classifier)
- **GARBEL — BIT-EXACT** (scratchpad/lpmpb/validate_garbel.jl vs FVSie_lpmpb_g fort.780/781/782):
  NACLAS=10 + all-10 class membership + class-PROB hex reproduce the oracle exactly on the 18-LP-tree stand.
- Spec (lpmpb/garbel.f + mpbdrv.f call): KEYMPB=`2,3,6*0,1` ⇒ active attrs A1=DBH (BETTER wt 1.0), A2=XPT/phloem
  (wt 4.0); IMP=KEY(9)=1 ⇒ PROB → CLAS col 1. NACLAS=MIN(ILP,NCLASS)=min(18,10)=10. PN1=0.5.
  Chain: GRPSUM standardize each attr (ML mean/stdev, P += wt·(x-ave)/stdv) → RDPSRT descending sort of IPT by P
  → method-1 NCL1=5 max-diff class boundaries (bubble-keep the NC1=4 largest diffs, IQRSRT ascending) → method-2
  split the NCL2=5 largest classes → sector pointers MPISC → CLAS(I,1)=Σ PROB over class members.
- **TWO faithful-port bugs found+fixed**: (1) `NCL1 = NCLAS*PN1+.5` is a Fortran REAL→INTEGER assignment that
  TRUNCATES toward zero (5.5→5); Julia `round(Int,5.5)`=6 (banker's) created a spurious 6th class ⇒ use
  `trunc(Int,·)`. (2) GARBEL's PROB arg is the REAL tree PROB array, NOT the WK3 work array passed as `P`.
- RDPSRT sort order was bit-exact on FIRST try (distinct phloem-driven P ⇒ tie-break moot); a plain descending
  sort matches. GRPSUM sqrt/divides in Float32 match gfortran (IEEE correctly-rounded sqrt).
## PROGRESS 2026-08-18 (sub-chunk 3 — SURFCE surface areas)
- **SURFCE — BIT-EXACT** (extends validate_garbel.jl, vs fort.783): all-10 class avg-DBH + LP surface + SNOHST.
- Spec (surfce.f/surflp.f/grclas.f): SNOHST = Σ over non-LP host trees (WP/WL/DF/PP by MPBSPM) of SUR5·PROB
  (SUR5≥0), per-species SUR/SUR5 = fns of ln(DBH),ln(HT),CFV,HT — this stand is all-LP ⇒ SNOHST=0.
  SURF(I)=SURFLP(CLASS(I,2)) where CLASS(I,2)=GRCLAS PROB-weighted class-avg DBH = Σ(DBH·PROB)/Σ PROB;
  SURFLP(d)= d≤5 ? d·0.672 : 8.835·d−40.82. IE map (mpblkdie.f): IDXWP1/WL2/DF3/LP7/PP10.
## MPBMOD SCOPING 2026-08-18 (the centerpiece — 808 ln, DETERMINISTIC ⇒ bit-exact target, NOT cornered)
- **KEY FINDING: MPBMOD calls NO RNG** (grep ran1/random/iseed = none; tafit.f "ran" is "RANGE" in a comment).
  ⇒ the whole LPOPDY population dynamics is a DETERMINISTIC year-by-year epidemic simulation → validate BIT-EXACT
  via dump-replay (no #206 cornering). The only "beta-dist" is BETIN = the incomplete-beta CDF (deterministic).
- **Helper leaves** (all in lpmpb/, port bottom-up, dump-replay each):
  - `EMERG(BY,T,INCRS)` (28 ln) emergence increment: T=0 → C=1/2^INCRS; else C=C·(INCRS−T+1)/T; EMERG=BY·C.
    ⚠ C is a DOUBLE local with NO SAVE but the build is `-fno-automatic` ⇒ C is STATIC/persists across calls in
    the emergence loop. In jl carry C as explicit state threaded through the INC loop (do NOT re-init per call).
  - `PERCNT(v,base)` = |base|≥1e-30 ? 100·v/base : 0.   `PMSLP(xx,x,y,n)` = piecewise-linear interp, flat-extrapolate.
  - `BETIN(a,b,x)` (78 ln, DOUBLE PRECISION) = regularized incomplete beta I_x(a,b) via continued fraction — the
    numerically sensitive one; port in Float64, match glibc. Used as AGG(I,INC)=TREES(I)·BETIN(DTA,DSMTA,XX).
  - `TAFIT` (90 ln) threshold-of-aggregation fit (solves a curve, picks root in range). `GENO`,`EXLOSS`,`AMP` = the
    genotype/flight-loss/attack-mult arrays (likely block-data or simple fns — read next).
  - `PTSYM/PTGRP/EVSET4` = graphics/event output → NO-OP in jl (MPBGRF path; verify they don't mutate WK2).
- **MPBMOD flow** (main driver): year loop over epidemic; per year: RESIST=PMSLP(...); TREES(I)=CLASS(I,IMPROB);
  emergence BNEW=EMERG; genotype brood B0(IG)=GENO·BNEW+BOLD; aggregation attack AGG=TREES·BETIN(DTA,DSMTA,XX);
  TREES−=AGG (attacked) then re-add survivors; TRKILL/TKYR/TM accumulate kill; brood B3; percentages; loop until
  epidemic ends; CLASS(I,IMPROB)=TREES(I) (survivors). Then MPBDRV: SURVIV(I)=CLASS(I,IMPROB)/SURVIV_pre;
  WK2(I)=max(WK2(I), PROB(I)·(1−SURVIV(class))) — the mortality output (DFB-style max-combine, cap PROB−1e-6).
## ★★ MORTALITY-TARGET FINDING 2026-08-18 (fort.784 dump — MASSIVELY simplifies the validation)
- For lp_popdy: MPBYR=0, NEPIYR=0 (single initial epidemic call, NACLAS=10). The post-MPBMOD class survivors
  CLASS(I,IMPROB)_after are ALL denormal-tiny (~1e-11 … 1e-16) ⇒ **the epidemic wipes out ~ALL lodgepole pine.**
- ⇒ SURVIV(class)=after/before ≈ 0 for every class ⇒ DEAD≈1 ⇒ per-tree the mortality-cap ALWAYS fires:
  **WK2(I) = PROB(I) − 1e-6** (exact; verified tree18 PROB 403758B4 → WK2 403758B0). This IS the "89→0" collapse.
- mpbdrv.f:198-208 mortality loop: X=PROB(I)·DEAD; WK2(I)=max(WK2(I),X); if PROB(I)−WK2(I)<1e-6 → WK2(I)=PROB(I)−1e-6.
  Since SURVIV<3.6e-7 (=1e-6/PROB) for all classes, PROB−X<1e-6 ⇒ the cap dominates ⇒ WK2=PROB−1e-6 REGARDLESS of the
  exact tiny survivor. ⇒ **end-to-end .sum is bit-exact as long as the ported MPBMOD drives each class survival <~3.6e-7**
  (robust to ULP-level survivor differences — do NOT need BETIN's continued fraction matched to the last bit).
- The full per-tree WK2 target is in /workspace/.iework/lpmpb/run/fort.784 (EPI/CLPOP/SURV/WK2 records).
- **NEXT**: port MPBMOD faithfully (year-loop epidemic → total kill) enough to drive SURVIV→~0; the leaf helpers
  (EMERG/PERCNT/PMSLP/TAFIT/BETIN) + GENO/EXLOSS/AMP/coeff block-data support it. Then wire mpb_apply! LPOPDY branch
  (drop early-return) applying WK2(I)=PROB(I)−1e-6 via the mpbdrv loop, validate end-to-end lp_popdy .sum (target 89→0).
  Attributes into MPBMOD: SNOHST=0, SURF(class), CLASS(I,1)=ΣPROB, CLASS(I,2)=avgDBH. Oracle FVSie_lpmpb_g
  @/workspace/.iework/lpmpb (instr/ dumps fort.780 GARBEL/781 inputs/782 sortP+M1/783 SURFCE/784 mortality target;
  .sum byte-identical to clean verified). BETIN pulls FORW/BACK/PQSML/DGAMMA (incomplete-beta pkg, DOUBLE).
