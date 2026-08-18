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
