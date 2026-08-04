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
2. ✅ DONE (commit d219ddb, 2026-08-04) — src/variants/easternmontana/establishment.jl (em_essubh_hht 19-sp +
   UHAB/UPRE/UPHY + em_ihtser habitat-code bracket search) + EM branch at establishment.jl:258 + RAN-skip.
   em_plant.key: NO CRASH, TPA BIT-EXACT every cycle, density cornered-late (2090 BA 42/46 SDI 101/109 QMD 5.0/5.3).
   TopHt low (52/60) = deferred-EMSQR-variance tail. Guard green. (Original plan text:)
   Port DF ESSUBH deterministically (`HHT=EXP(PN)`, EMSQR variance DEFERRED — like CI's intentional ZZRAN
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
4. ✗ NOT NEEDED — the establishment-height variance is DEFINITIVELY irrelevant to the residual. CONTROLLED
   EXPERIMENT (2026-08-04): injecting the oracle's EXACT 50 per-record HHTs (from em_essubh_df_groundtruth.txt)
   into jl leaves the .sum COMPLETELY UNCHANGED (2090 BA 42/46, TopHt 52/60, QMD 5.0/5.3 — identical to
   deterministic EXP(PN)). So even with perfect heights the establishment cohort grows to TopHt 52 vs live 60 over
   100 yr ⇒ **the ~8% density + TopHt gap is EM REGENT/DG seedling GROWTH on the establishment cohort (D~0.1
   origin), NOT the essubh height.** The EM essubh height port (deterministic EXP(PN)) is therefore COMPLETE +
   CORRECT — the EMSQR/DILATE stochastic reconstruction (below) is confirmed NOT worth pursuing (would not move the
   .sum). REMAINING EM-estab lead = the seedling-growth trajectory (regent/DG on tiny establishment-origin trees),
   a separate low-pri cornered tail needing multi-cycle per-tree instrument-replay. (Earlier reconstruction attempt,
   for reference:) the EMSQR/DILATE reconstruction is **.sum-INERT** so it was NOT kept.
   Implemented HHT=EXP(PN+EMSQR·DILATE·BNORML·σ) using the two emsqr draws jl already consumes + a per-species
   running DILATE. Findings: **DILATE reconstructs BIT-EXACT** (0.1,0.3162,0.5623,… matched live); emsqr record-1
   matched (0.2186) but **desyncs record-2+** (jl's per-replicate draw count diverges from FVS after record 1 —
   likely FVS also draws emsqr for the advance/subsequent essubh paths, not just the plant path); IAGE uses Fortran
   INT (trunc), not round (jl gave 6, FVS 5). CRITICAL: even with the variance the .sum barely moved (2090 TopHt
   52→53, density unchanged) ⇒ **the ~8% density gap is NOT the establishment-height variance — it's seedling
   GROWTH (EM regent on D~0.1 seedlings) or another downstream factor.** So step-4 is LOW-VALUE; the real remaining
   EM-establishment lead is the seedling growth trajectory, not essubh. Deterministic HHT=EXP(PN) kept (cleaner,
   equally cornered, no wrong-draw noise). Reconstruction MECHANISM for reference (em/estab.f:646-650, 800-840):
   `EMSQR = (esrann<0.5 ? -1 : +1) · esrann` — TWO ESRANN draws PER REPLICATE, which jl ALREADY consumes+discards
   at establishment.jl:218 (just use them). `DILATE = FIRST(2,sp)`, a PER-SPECIES running order-statistic:
   init 0.1, then `FIRST(2,sp)=SQRT(DILATE)` after each record (0.1→0.316→0.562→0.75→…→1 — the tallest-of-N
   dilation shrinks as N grows). Per-record IPREP from the WK6 site-prep vector (`DRAW=WK6(NDRAW)`, jl fills WK6 at
   line 207). Then HHT = EXP(PN + EMSQR·DILATE·BNORML(IAGE)·σ_sp), IAGE=INT(age−TRAGE+0.5). Validate TopHt on
   em_plant.key (oracle 2000 TopHt=5 vs current jl=1). Since TPA is bit-exact the stream is already aligned to the
   record level, so this is reconstruction (use the consumed draws), not new RNG plumbing.
5. em/esadvh.f advance-regen (THAB/TPHY/TPRE tables) — the second establishment-height source (natural advance regen).

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

## Gap A-residual — EM establishment seedling-growth ROOT CAUSE PINNED (2026-08-04)
The ~8% density + TopHt lag on EM establishment stands (em_plant + emt01 ESTAB) is NOT the essubh height
(controlled experiment proved injecting exact heights is .sum-inert) — it is **establishment-cycle growth timing**.
jl's `establish!` runs at simulate.jl:558, AFTER the growth/UPDATE (line 533-534) — the SN "fresh this period"
model (establishment.jl:7-14): trees enter at their essubh height with NO growth in the establishment cycle. But
the EM oracle's establishment trees appear GROWN (emt01 ESTAB 2002 TopHt=10 while essubh height is ~3; em_plant
2000 TopHt=5 while essubh ~1.2) ⇒ **FVS grows EM establishment trees IN their establishment cycle; jl does not.**
So jl's establishment cohort is ~1 growth-cycle behind → persistent TopHt lag (density aggregates converge/cross
because mortality compensates). **VERIFIED EXACT MECHANISM (estb/esgent.f, shared):** line 46-58 —
`'GROW' TREES TO THE END OF THE CYCLE` via `CALL REGENT(.TRUE.,ITRNIN)` then `HT(I)=HT(I)+HTG(I)`, DBH re-derived
from the grown height. So FVS grows the freshly-established trees ONE cycle via REGENT before summarizing (SN's
buildDir esgent.f has the SAME call). jl's establish! OMITS this REGENT grow-to-cycle-end. jl handles SN OK because
SN's essubh uses the age-at-cycle-end site curve (htcalc_height at the grown age = already the cycle-end height),
but **EM's essubh gives the establishment-age (age-7) height and RELIES on the esgent REGENT grow** — which jl
skips. FIX (EM-gated, safe from SN): after creating the EM establishment trees, grow them one cycle via
small_tree_growth!(::EasternMontana) (mirror esgent.f: HT+=HTG, DBH from grown HT, HHTMAX cap). Intricate (needs
the regent setup on the fresh cohort); LOW-PRIORITY (density cornered, TopHt = the visible tail). This upgrades the
#137 residual from hypothesis to a SOURCE-VERIFIED mechanism (estb/esgent.f CALL REGENT) with a specific EM-gated fix.

## Gap A-residual — em_esgent! IMPLEMENTED + REVERTED (2026-08-04): exposes a WL+PP density over-growth
Implemented em_esgent! (mirror cr_esgent!: grow the establishment cohort nstart+1..t.n via SMHTGF height + SMDGF
diameter, partial period scale=(fint−gentim)/regyr, HT+=HTG/DBH from grown height/HHTMAX cap; wired at
simulate.jl after cr/tt_esgent!). RESULT on em_plant (DF): TopHt MAXΔ 0.7→0.2 (2000 jl 1→4 vs or 5; 2090 52→58 vs
60) — FIXES the lag; but BA flipped 8%-under → 13%-over. On emt01 (WL+PP): TopHt better (69/71) but BA WORSE
(2090 267/218 = 22% over, up from essubh-only's 13% over). So em_esgent UNIFORMLY adds growth: fixes DF, worsens
WL+PP. Guard stayed green (no establishment there); stand-1 growth-only unchanged. REVERTED (doctrine #4: a
faithful chunk exposing a regression ⇒ examine the oracle, don't ship). The esgent grow IS faithful (esgent.f does
it), so this EXPOSES a pre-existing WL+PP establishment DENSITY over-growth (essubh-only already 13% BA over for
WL+PP while DF was 8% UNDER) — likely the WL/PP essubh height too tall, or the SMDGF/mortality on the WL/PP
establishment cohort. EXAMINED (2026-08-04, doctrine #4): instrument-replayed the FVSem oracle WL(sp2)/PP(sp10)
essubh on emt01 — oracle WL HHT-mean 2.614, PP 1.253 (IHTSER=2/IPREP=1/IPHY=3/AGE=7/BAA=1/ELEV=54). jl matches BOTH
(WL 2.6135, PP 1.2536; jl stores elevation=54 correctly, keyword_dispatch.jl:550). So WL/PP ESSUBH HEIGHTS ARE
CORRECT — elevation-bug hypothesis REFUTED. ⇒ the WL+PP over-growth is NOT essubh; it's the EM REGENT's SPECIES-
DEPENDENT seedling growth on the D~0.1 cohort (DF under-grows, WL/PP over-grow — same class, opposite sign). The
em_esgent grow (faithful) amplifies this pre-existing regent bias. ROOT = EM regent tuning on tiny establishment
seedlings (below the guard's emt01_smallr D0.4-3.8 range) — a deep per-cycle regent examination (cornered tail).
NARROWED (2026-08-04): the .sum deltas split the residual into TWO parts — (1) HTG under for BOTH DF and WL/PP
(the TopHt lag = the missing esgent in-cycle grow), and (2) DG species-split: DF DBH UNDER (em_plant QMD 5.0/5.3),
WL/PP DBH OVER (emt01 QMD/BA over). The DG split MATCHES the SMDGF FORM split: DF=sp3 uses the LINEAR SMDGF
(sp{3,7,8,9}); WL=sp2/PP=sp10 use the HLESS4 form (h−4.5, sp{1,2,10,18}). For the establishment cohort at small h
(~1.2-3 ⇒ hl=h−4.5 < 0) the HLESS4-form EXTRAPOLATES below its validated range (guard's emt01_smallr is D0.4-3.8, h
mostly >4.5) ⇒ over-grows WL/PP DBH. FIXED (commit 876ea0c): confirmed jl _em_smdgf goes NEGATIVE below h~2.5 (WL@1.2=-0.198, PP@1.2=-0.074) while
em/regent.f:578 `IF(H2.LE.4.5)GO TO 14` SKIPS the SMDGF/DBH below 4.5. Added the h>4.5 guard to small_tree_growth!
(only update the SMDGF DBH for h2>4.5; start DBH = actual d, not SMDGF, below 4.5). A real latent bug — guard's
emt01_smallr (h>4.5) never triggered it, no regression. IMPROVES em_plant DF establishment BA 42→44/46 (8.7→4.3%
under). BUT emt01 (WL+PP) UNCHANGED ⇒ the WL/PP BA over-growth is NOT the SMDGF-h<4.5 (hypothesis REFUTED by
measurement). So the WL/PP over-growth root remains OPEN (two hypotheses refuted: elevation, SMDGF-h<4.5) — a
separate species-dependent establishment-cohort density issue, deferred to a per-cycle instrument-replay session.
The SMDGF guard is faithful+shipped regardless. Then re-visit em_esgent (git history) once WL/PP is understood.
NARROWED to MORTALITY (2026-08-04, 3rd hypothesis refuted): instrument-replayed the FVSem SMDGF call on emt01 —
jl _em_smdgf matches FVS BIT-EXACTLY at the dumped inputs (PP: jl 0.4672 vs FVS 0.46724, 0.6947/0.69471,
0.8793/0.87934, 1.0874/1.08735). So the DBH-growth FORMULA is correct; ALL three EM establishment growth models
(essubh, SMHTGF, SMDGF) are now validated/matching. The WL/PP BA/QMD over-growth is therefore NOT growth — the
emt01 ESTAB trajectory shows QMD MATCHES early (2002 jl 0.7/or 0.7) and diverges only LATE (2092 8.7/8.0) as TPA
over-kills (2092 jl 590 vs or 630). ⇒ ROOT = the LATE SELF-THIN MORTALITY on the dense WL/PP establishment cohort
(over-kill leaves fewer, thicker trees → QMD/BA over). SAME CLASS as BM #140 (self-thin on dense stands). NEXT:
per-cycle mortality/self-thin comparison on the dense EM establishment cohort (shared with #140); NOT a growth fix.

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
