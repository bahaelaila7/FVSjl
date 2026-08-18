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
## ★★★ MPBMOD FULLY READ + ALL PARAMETERS CAPTURED 2026-08-18 (fort.785) — port is now fully specified
Entire LPOPDY chain read: MPBDRV, phloem, GARBEL, GRCLAS, SURFCE, SURFLP, MPBMOD (all 808 ln), EMERG, PERCNT, PMSLP,
TAFIT, BETIN, FORW, BACK, PQSML, MPBGAM. All params dumped from the oracle (fort.785, lp_popdy, IE):
- **Scalars**: NG=2, INCRS=10, IB=1, MPMXYR=10, NACLAS=10. CE=1.0, EXCON=640.0 (acres), SQFTPA=43560.0, STRBUG=500.0
  (init beetles), STRP=0.95 (init P), SEXRAT=0.66, HS=1.0, CF1=CF2=0.01, CF3=0.5, TAFAC=1.7, TAMIN=1.7, TAMAX=3.0,
  **TA=2.099609** (=RESIST via PMSLP(PGR,PGRX,TAY,5) clamped [TAMIN,TAMAX]; LCRES=T; PGR from MPGR — capture TA directly),
  AMP1=1200.0, AMP2=600.0, CRITAD=1.5, ELEV=34.0 (hundreds-ft), FORLAT=44.0, EFELEV=1.0, EFLAT=1.00019(=4.667−.08333·44),
  EPS=1.0D-6, BMIN=1.0D-10 (mpbint.f), BETTER=[1.0,4.0], KEYMPB=[2,3,0,0,0,0,0,0,1], IMPROB=1.
- **Genotype arrays** (NG=2): DST=[3000.0, 500.0] (flight dist), EXODUS=[0.27588, 0.008926] (=1−exp(−CE·DST²/(EXCON·SQFTPA))).
- **Switches**: LGO=T (actual sim from start since NEPIYR≤0 ⇒ skip 3-try TA calibration), LCRES=T, **LAGG=LREP=LPS=LDC=F**
  ⇒ ALL pheromone/spray/direct-control paths INERT (AGGPH/REPL/PSPK/PSDL/PSPF/PSE1/PSE3/DCPF/DCPK unused). Big simplification.
- **Per-class inputs** (fort.785 'C' recs, all match GARBEL/SURFCE): TREES(I)=CLASS(I,1)=ΣPROB, SURF(I), DIAM(I)=CLASS(I,2)
  =avgDBH, PHLOEM(I)=CLASS(I,3)=avgXPT. Derived: SEXDBH(I)=0.918−0.0168·DIAM(I); EFPHLM(I)=max(0,16.67·PHLOEM(I)−0.667).
- **MPBMOD year-loop** (LGO=T path, no partial-epidemic branch): BY=STRBUG, P=STRP. Repeat years (MPBYR++) until BY<1
  or MPBYR≥MPMXYR: GENO=[P,1−P]; GTEB=BY; TEB=BY·GENO. EMERGENCE loop INC=1..NINC(=INCRS+1=11): E1=Σ SURF·TREES;
  E3/EF3/TAGG from prior-inc AGG·AMP (AMP(KK)=max(0,AMP1−AMP2·AD(KK)·(1−SEXRAT))); RHO1=ΣTREES; E2=OS−(E1+E3)+SNOHST
  (OS=Σ SURF·TREES [+SADLPP if LGO]); EFFS=E1+E2+EF3; RHO3=TAGG/SQFTPA, RHO2=max(0,TPROB−TAGG−RHO1)/SQFTPA,
  RHO1/=SQFTPA; BNEW=EMERG(BY,INC−1,INCRS); GENOTYPE loop IG: B0=GENO·BNEW+BOLD; EXLOSS=B0·EXODUS; B0−=EXLOSS;
  B1=B0·E1/EFFS,B2=B0·E2/EFFS,B3=B0·EF3/EFFS; FMi=exp(−CFi·2·√RHOi·DST²/DST(1)) (0 if arg≤−80); Bi−=Bi·FMi;
  BOLD=B1+B2; sums B1INC,B3INC,B3SUM. PIODEN=B1INC/E1; XX=1−exp(−PIODEN); AGG(I,INC)=TREES(I)·BETIN(TA,SURF(I)−TA+1,XX)
  [if XX≠0 & DSMTA>0]; TREES−=AGG; AD(KK)+=AMP(KK)·B3INC/EF3. PRODUCTIVITY loop INC: EGGS=630·(1−exp(−0.117·AD));
  PSURV=1−exp(−AD·.04328) [·6.812·exp(−1.191·√AD) if AD≥2.595]; YOUNG=EGGS·PSURV·EFELEV·EFLAT·EFPHLM(I);
  if AD<CRITAD → strip-kill (TREES+=AGG, AGG=0) else TM(I)+=AGG; SKILL=SURF·AGG; BY+=YOUNG·SKILL·HS·SEXDBH(I).
  Update P from B3SUM (NG≠3 ⇒ P=(B3SUM(1)/B3ALL)²). CLASS(I,IMPROB)=TREES(I). Terminate.
- **BETIN pkg** (DOUBLE): BETIN(a,b,x)=regularized incomplete beta via BACK (backward continued fraction, EPS/BMIN
  convergence) + FORW (forward recurrence) + PQSML (series, ·MPBGAM(p+q)/(MPBGAM(p)·MPBGAM(q))); MPBGAM=exp(Stirling
  log-gamma, shift arg to ≥18 by TERM-multiply). Restrictions a,b>0, 0≤x≤1. Route DEXP/DLOG via glibc for bit-exact.
- **NEXT (the coding)**: write src/engine/lpmpb LPOPDY module: BETIN pkg + EMERG(state C)/PERCNT/PMSLP + MPBMOD year-loop;
  dump-replay a mid-epidemic TREES(I)/AGG/BNEW/BY snapshot from the oracle (add a per-year dump to instr/mpbmod.f) to
  validate the loop; then wire mpb_apply! LPOPDY branch (drop early-return) → WK2(I)=PROB−1e-6 via the mpbdrv loop →
  end-to-end lp_popdy .sum (target 89→0). All params above are exact (fort.785); mortality target in fort.784.

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


## ★★★ MPBMOD YEAR-LOOP PORTED + EPIDEMIC REPRODUCED 2026-08-18 (validate_mpbmod.jl)
The full deterministic MPBMOD year-loop is ported in Julia (LGO path; LAGG/LREP/LPS/LDC=F) and RUNS the epidemic:
- yr1 per-increment AD builds 0→5.097 MATCHING the oracle exactly; yr1 BY=1912.432 = oracle 1912.432 (exact at print prec).
- The epidemic collapses LP: ΣTREES 97→96→92→81→53→14→0 (TOTAL KILL by yr6), BY 1912→6858→22358→47690→40294→10016.
- 2 bugs fixed while debugging: (1) fort.785 'G' parse off-by-one (DST/EXODUS swapped → B0 blew up); (2) top-level
  soft-scope on B1INC accumulator → wrapped the sim in run_mpbmod() (function scope).
- RESIDUAL (the ONE remaining item): the oracle retains tiny survivors in the 3 smallest-DBH classes (8,9,10:
  ~3e-7 / 1.3e-5 / 9.4e-4 TPA); my port over-kills them to 0 due to a small multi-year BY drift (yr2 6858 vs 6853
  ~0.08%, compounding). Candidate causes: (a) gfortran computes PIODEN=B1INC/E1 in Float32 THEN widens to DOUBLE
  (my port divides in Float64); (b) AGG=TREES*BETIN done in DOUBLE then rounded to REAL (my port Float32*Float32);
  (c) BETIN 1.56e-8. Fix = tighten those to gfortran's exact Float32/Float64 operation order.
- ★ BUT the over-kill is ~0.001 TPA TOTAL (classes 8+9+10 survivors) out of ~89 LP TPA — almost certainly BELOW the
  .sum's integer-TPA rounding. ⇒ **the DEFINITIVE test is the end-to-end .sum, not the tail-class survivors.**
- **NEXT**: wire the LPOPDY branch into FVSjl (mpb_apply! — drop the early-return; run phloem→GARBEL→SURFCE→MPBMOD→
  WK2(I)=PROB(I)-1e-6 or PROB(I)*(1-SURVIV) via the mpbdrv loop) and run the end-to-end lp_popdy .sum vs FVSie_lpmpb
  (target 89→0). If the .sum is bit-exact, LPOPDY is DONE (tail survivors round away). If not, tighten the Float32/
  Float64 op-order above. All the ported+validated pieces live in scratchpad/lpmpb/: validate_phloem.jl, validate_garbel.jl
  (GARBEL+SURFCE), validate_betin.jl + betin_pkg.jl, validate_mpbmod.jl.
