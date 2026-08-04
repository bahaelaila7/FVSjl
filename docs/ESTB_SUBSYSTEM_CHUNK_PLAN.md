# ESTB subsystem chunk plan — western establishment (essubh/esadvh + AUTOES tally)

Consolidated scoping from the 2026-08-04 western-regen-sweep session (measured, doctrine #2).
This turns the vague "auto-establishment gap" into an executable, reproduction-backed chunk plan.
Related tasks: #137 (EM essubh), #143 (AUTOES tally). Related fix commits this session: IE #139, KT #141, CI #142.

## What already exists (do NOT re-port)
- The **explicit PLANT/NATURAL establishment engine** works in jl (`src/engine/establishment.jl`): record
  creation, dup/replicate counting (NPTIDS·IDUP), the `:estab` RNG stream (incl. the two per-replicate `emsqr`
  draws consumed at line ~218 but currently discarded), the post-base `bachlo` height perturbation (line ~262),
  and per-variant ESSUBH base-height branches for NE/CS/LS/CR/IE/TT (lines ~237-260).
- `establish!` fires ONLY on scheduled PLANT(icflag 430)/NATURAL(431) actions — see "AUTOES" below.

## Gap A — EM ESSUBH base height (task #137). Fixes a real CRASH.
**Symptom:** an EM stand with an explicit PLANT/NATURAL keyword crashes at `establishment.jl:259`
(`htcalc_height(bc=Nothing)` → `MethodError getindex(::Nothing,::Int64)`) because EM has no ESSUBH branch.
**Reproduction (READY):** `/workspace/.emwork/em_plant.key` (bare + `PLANT DF(sp3) 400`), oracle
`em_plant.sum` runs clean (establishes ~364 DF @ TopHt 5 by 2000 → 60 ft by 2090). `FVSem_clean` builds it.

**Model (em/essubh.f):** 19-species `HHT = EXP(PN + EMSQR·DILATE·BNORML(IAGE)·sigma_sp)`, where per species
`PN = intercept + b1·ALOG(AGE) + b2·BAA + UHAB(IHTSER,sp) + UPRE(IPREP,sp) + UPHY(IPHY,sp) + [per-sp aspect
XCOS/XSIN, slope SLO, elev ELEV/ELEV², habitat flags BWAF/BWB4]`. DF(sp3):
`PN = -2.16416 + 1.28151·ALOG(AGE) - 0.0031363·BAA + UHAB(IHTSER,3) + UPRE(IPREP,3) + UPHY(IPHY,3)
      - 0.09626·XCOS - 0.23946·XSIN - 0.14589·SLO`, sigma=0.55942. LM(4)/RM(6) = fixed HHT=0.5.
Advance-regen counterpart = em/esadvh.f (THAB/TPHY/TPRE tables). Tables UHAB(5,19)/UPRE(4,19)/UPHY(5,19),
BNORML(table) in em/blkdat.f.

**Inputs (em/estab.f:470-493, per-plot):** IPREP=ITYPEP (siteprep; default NONE=1 ⇒ UPRE(1,sp)=0);
IPHY=IPHYS(NNID) (physiographic 1-5); BAA=BAAA(NNID) competition BA clamp[1,400] (BARE ⇒ 1);
IHAB=IPHAB(NNID) habitat GROUP 1-16; IHTSER=MYHTS(IHAB), `MYHTS=/1,3*2,4*3,2*4,.../`,
`MYHABG=/4*1,4*2,3,4,6*5/`. XCOS/XSIN from aspect, SLO=slope fraction, from STDINFO.

**The stochastic wrinkle:** EMSQR=±DRAW (em/estab.f:647-650, DRAW random) and DILATE=FIRST(2,sp) a per-record
running order-statistic product (`FIRST(2)=SQRT(DILATE)` each record, line 837/1039) → the "tallest of N
subsequent" spread. Since ±DRAW is symmetric (mean≈0), **HHT ≈ EXP(PN)** is the deterministic expectation.

**Recommended sub-chunks:**
1. Instrument-replay: patch em/estab.f (buildDir) to WRITE `sp, AGE, BAA, IHTSER, IPREP, IPHY, XCOS, XSIN, SLO,
   ELEV, EMSQR, DILATE, HHT` for each planted record on em_plant.key; relink `relink_em.sh trc estab.o`
   (uses plain `gfortran`; if absent use gfortran-16 + `/workspace/.crwork/isoc23_shim.o`); run → dump the
   ground-truth inputs+HHT. (Pins IHAB→IHTSER and IPHY, which jl's EM estab context may not yet derive.)
2. Port DF ESSUBH deterministically (`HHT=EXP(PN)`, EMSQR variance DEFERRED — exactly like CI's intentional
   ZZRAN deferral, RNG stream stays synced because the emsqr draws are already consumed). Add the EM branch at
   `establishment.jl:258`. Validate BA/QMD/SDI cornered on em_plant.key; TopHt = deferred-variance tail.
3. Extend to the other EM conifers (WB/WL/LP/ES/AF/PP), then the IE-borrowed species (LM/RM/AS/CW/… use IE forms).
4. (Later refinement) faithful stochastic EMSQR/DILATE order-statistic → bit-exact TopHt.

## Gap B — AUTOES automatic-establishment tally (task #143). LARGE.
**Symptom:** stands relying on default automatic natural regen after disturbance collapse in jl (iet01 stand-4
"SHELTERWOOD WITH AUTO REGENERATION": oracle TPA 1025→1788, jl 224→28). jl's `establish!` has NO AUTOES path —
only explicit PLANT/NATURAL. This is the estb/ natural-regen tally/stocking model (~6280 lines; the memory's
`estab_chunk_plan.md` covered only the BARE PLANT path).
**Validation anchor:** iet01.key stand-4, oracle `.sum` rows 40-50 (`/workspace/.iework/ierun/iet01.sum`).
Needs its own chunk plan (stocking trigger + per-species seed/sprout tallies + the essubh heights from Gap A).

## Doctrine reminders for this subsystem
Validate bit-exact-or-cornered vs live per sub-chunk; MEASURE (instrument-replay) don't infer; per-record
treelist diff INVALID after tripling (use .sum aggregates); a deferred stochastic component (EMSQR, like CI
ZZRAN) is acceptable-cornered when documented, but port the deterministic mean faithfully.
