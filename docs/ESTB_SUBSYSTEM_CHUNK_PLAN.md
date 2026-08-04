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
IHAB=IPHAB(NNID) habitat GROUP 1-16; IHTSER=MYHTS(IHAB). XCOS/XSIN from aspect, SLO=slope fraction.
**IHTSER derivation (traced 2026-08-04):** EM has TWO habitat classifications — ITYPE(1-30, growth model,
which jl HAS as habitat_input) and **IHAB(1-16, ESTAB-specific**, indexed into 16 HABTYP strings at em/estab.f:112,
which jl LACKS). Tables (em/estab.f:109-111 DATA): `MYHTS(16)=/1,2,2,2,3,3,3,3,4,4,6*5/`,
`MYHABG(16)=/4*1,4*2,3,4,6*5/`, `MYTYPE(30)=/9*1,2*5,2*2,3*3,4,12*5,1/` (ITYPE→ISER), `MAXTPP(16)`,
`MAXSPP(16)=/4,3*3,5,4,6,4,2*6,4,2*5,4,6,4/`, `MAXING(16)`. So the port needs the EM estab habitat-decode (habitat
code → IHAB 1-16), an estab-setup piece jl doesn't have. For em_plant: IHTSER=2 ⇒ IHAB∈{2,3,4}. First-cut: hardcode
IHTSER=2/IPHY=3 to validate the formula+integration on em_plant, then port the IHAB decode to generalize.
DF-verified coeffs from the dump: UHAB(2,3)=-0.03354, UPRE(1,3)=0.0, UPHY(3,3)=-0.12317.

**The stochastic wrinkle:** EMSQR=±DRAW (em/estab.f:647-650, DRAW random) and DILATE=FIRST(2,sp) a per-record
running order-statistic product (`FIRST(2)=SQRT(DILATE)` each record, line 837/1039) → the "tallest of N
subsequent" spread. Since ±DRAW is symmetric (mean≈0), **HHT ≈ EXP(PN)** is the deterministic expectation.

**Recommended sub-chunks:**
1. ✅ DONE 2026-08-04 — Instrument-replay: patched buildDir/essubh.f DF branch (label 30) to WRITE
   `I,IHTSER,IPREP,IPHY,AGELN,BAA,XCOS,XSIN,SLO,UHAB,UPRE,UPHY,EMSQR,DILATE,BNORM,PN,HHT` to unit 16 (→ .out),
   compiled gfortran-16, relinked manually (`gfortran-16 -o FVSem_trc $(ls buildDir/*.o) crwork/isoc23_shim.o`
   — relink_em.sh's plain `gfortran` is absent in non-interactive bash), ran em_plant.key. **VALIDATED the DF
   formula bit-exact:** for em_plant IHTSER=2, IPREP∈{1,2,3} (per-record, from the site-prep vector!), IPHY=3,
   AGE=7 (AGELN=1.9459), BAA=1.0. Hand-check PN=0.156307 ✓, HHT=EXP(PN+EMSQR·DILATE·BNORM·0.55942)=1.18492 ✓.
   Derived: **XCOS=SLO·cos(aspect), XSIN=SLO·sin(aspect)** (0.3·cos315°=0.2121); IAGE=INT(age−TRAGE+0.5)=5 (age
   BEFORE trage) → BNORML(5)=1.093. mean(HHT over 50 recs)=1.188 ≈ EXP(PN_none)=1.169 (Jensen +1.6%). buildDir
   RESTORED pristine. Ground truth: `/workspace/.emwork/em_essubh_df_groundtruth.txt`.
2. Port DF ESSUBH deterministically (`HHT=EXP(PN)`, EMSQR variance DEFERRED — like CI's intentional ZZRAN
   deferral; RNG stays synced, the emsqr draws are already consumed at establishment.jl:218). Add the EM branch at
   `establishment.jl:258`. **jl-input scouting DONE 2026-08-04** (src/core/state.jl): AVAILABLE — `slope`(0..1),
   `aspect`(rad) ⇒ XCOS=slope·cos(aspect)/XSIN=slope·sin(aspect); `elevation`(100s ft); `basal_area` ⇒ BAA=clamp[1,400];
   `physio_region`(IPHREG, line 379) — CONFIRM == IPHY(1-5 position; em_plant IPHY=3); `habitat_input`(ITYPE)/
   `habitat_code`(KODTYP). MISSING — IHTSER derivation: FVS IHTSER=MYHTS(IHAB), IHAB=IPHAB(NNID) a 1-16 habitat
   GROUP; need the EM habitat→IHAB→MYHTS(=/1,3*2,4*3,2*4,.../) chain (trace IPHAB setup in em/estab.f/esdlay.f;
   em_plant IHTSER=2). Also verify jl's estab `age`==7 (AGELN 1.9459) for em_plant. IPREP: simplest =1(NONE)⇒HHT≈1.169
   cornered; faithful = per-record WK6 site-prep vector (jl fills at line 207). Validate BA/QMD/SDI cornered on
   em_plant.key vs em_plant.sum; TopHt = deferred-variance tail. Coeff tables UHAB(5,19)/UPRE(4,19)/UPHY(5,19) +
   per-sp PN intercepts/slopes: dump from em/essubh.f + em/blkdat.f (same table-dump technique as other chunks).
3. Extend to the other EM conifers (WB/WL/LP/ES/AF/PP — em/essubh.f labels 10/20/70/80/90/100), then the
   IE-borrowed species (LM/RM/AS/CW/… use IE forms). Same instrument-replay per species if uncertain.
4. (Later refinement) faithful stochastic EMSQR/DILATE order-statistic + per-record IPREP → bit-exact TopHt.

**Complete em/essubh.f spec (all 19 species, extracted 2026-08-04) — deterministic HHT=EXP(PN):**
```
sp1  WB: PN=-1.51302+1.24537·L-0.003052·B+UPHY(iphy,1)                                    σ.46010
sp2  WL: PN=-1.36257+1.21548·L-0.003797·B+UHAB(ih,2)+UPRE(ip,2)                            σ.52668
sp3  DF: PN=-2.16416+1.28151·L-0.0031363·B+UHAB(ih,3)+UPRE(ip,3)+UPHY(iphy,3)
         -0.09626·XC-0.23946·XS-0.14589·SLO                                               σ.55942 [VALIDATED]
sp4  LM: HHT=0.5                                     sp6 RM: HHT=0.5
sp5  LL & sp9 AF: PN=-2.06377+1.18184·L-0.0044465·B+0.06615·XC+0.03085·XS-0.37402·SLO      σ.56740
sp7  LP: PN=-0.27105+1.32027·L-0.008208·B+UPRE(ip,7)+UPHY(iphy,7)+UHAB(ih,7)-0.15385·XC
         +0.04156·XS-0.49186·SLO-0.04744·E+0.0003511·E²+0.01105·BWAF+0.02588·BWB4          σ.47557
sp8  ES: PN=-2.93213+1.43503·L-0.002504·B+UPRE(ip,8)+UPHY(iphy,8)+UHAB(ih,8)               σ.48951
sp10 PP: PN=-1.99480+1.53946·L-0.00402·B+UHAB(ih,10)+UPRE(ip,10)-0.01155·E                 σ.49076
sp11-17 GA/AS/CW/BA/PW/NC/PB: HHT=5.0    sp18 OS: PN=-2.42379+1.52366·L-0.003256·B σ.54116   sp19 OH: HHT=5.0
```
L=ALOG(AGE), B=BAA clamp[1,400], XC=SLO·cos(asp), XS=SLO·sin(asp), E=ELEV(100s ft), ih=IHTSER, ip=IPREP, iphy=IPHY.
Tables (rows = the species with nonzero coefs; cols = the index):
```
UHAB(5,·): WL/2=[-.01541,-.03814,.11409,.35334,0] DF/3=[-.21858,-.03354,.22756,.51988,0]
           LP/7=[-.29969,-.15449,.04545,-.00601,0] ES/8=[0,0,.18740,.26511,0] PP/10=[-.02287,-.14710,.19278,.13817,0]
UPRE(4,·): WL/2=[0,-.11310,-.06246,.009632] DF/3=[0,.06961,.19508,.17952] LP/7=[0,.11502,.02486,.13080]
           ES/8=[0,.10587,.27072,.16240] PP/10=[0,.20729,.18491,.11864]
UPHY(5,·): WB/1=[-.18731,-.48682,-.32160,-.16113,0] DF/3=[-.27801,-.20433,-.12317,-.26736,0]
           LP/7=[.32401,.14743,.22165,.24559,0] ES/8=[.41120,.01164,.22217,.15834,0]
```
**INTEGRATION FINDING — RESOLVED (2026-08-04):** the EM no-treeht PLANT path genuinely differs from CR's.
CR estab.f:486-489 (matches jl's shared establishment.jl:269-274): `RAN=BACHLO(0.5,0.25) reject∉[0,1.5]; HHT+=RAN;
HHT+=HTADJ; floor XMIN`. **EM estab.f:1035-1038 OMITS the RAN**: just `HHT+=HTADJ(sp); floor XMIN` — no BACHLO draw.
⇒ jl's EM branch must (a) NOT add the RAN perturbation, and (b) NOT consume the RAN `bachlo` draw at
establishment.jl:271 (else the `:estab` stream desyncs vs FVSem). i.e. the EM path is `HHT = EXP(PN) + HTADJ(sp)`
(HTADJ default 0), floor XMIN, cap HHTMAX — deterministic. Matters for the mid/late .sum (which cycle a seedling
crosses 4.5 ft → large-tree DGF), not the initial DBH (~0.1 for any HT<4.5). BWAF/BWB4 (LP habitat flags) still to
source. Remaining sub-steps: (1) branch establishment.jl on EM to skip the RAN draw + use HHT=EXP(PN)+HTADJ;
(2) derive IHTSER(IHAB decode)/IPHY (em_plant: 2/3); (3) validate em_plant.key (DF) cornered vs em_plant.sum.

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
