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
- NEXT sub-chunks (dep-order): GARBEL classifier → SURFCE → MPBMOD core (+ betin/forw/back/gamma beta-dist,
  emerg emergence) → wire mpb_apply! LPOPDY branch (drop early-return) → end-to-end lp_popdy .sum (target 89→0).
