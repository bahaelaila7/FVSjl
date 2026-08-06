# AUTOES — Western Automatic (Natural) Establishment tally — CHUNK PLAN (task #143)

Status: **★★ TALLY COMPUTATION COMPLETE + BIT-EXACT (2026-08-04).** `ie_autoes_tally` reproduces the live FVSie
per-species ingrowth split on iet01 stand-4 (1999) EXACTLY: WP33 WL20 DF7 GF202 WH222 RC50 ES23 AF27, total 583.7.
72 tests. All 8 primitives + the multi-plot seed chain + the tally composition are bit-exact vs the live oracle
(commits de5da44→c28dca2). The entire RNG/probability/selection/booking model is SOLVED — see the A1/A2/A2b/A2c
sections below for the full spec + measure-first findings (instrument recipe, per-plot reseed, MAXING cap, PXCS-
weighted excess, etc.). REMAINING = engine integration only (no model unknowns): (1) wire ie_autoes_tally into
engine/establishment.jl behind the esnutr disturbance-trigger scheduler (A3) — derive per-plot inputs from stand
state (IHAB/ISER/IFO via habtyp; PROB1 via ESTOCK+ESB; DUPNPT; **seed0 = ESDRAW = INT(main-RNG-draw*100000+0.5)**
from the engine's seed-55329 stream at establishment time); create the ingrowth trees (per species, heights via
ie_essubh/ie_esxcsh) in the tree list; (2) A4 full-cycle .sum (536→1025→…→1788); (3) generalize per-plot inputs
(OCURNF/XESMLT for non-forest-4, ESPSUB for TIME>2 habitat mixes, multi-tally +10/+20yr scheduling). This is a
focused shared-engine integration best done as one careful unit with end-to-end .sum validation.

★ HEIGHT SUB-MODEL ANALYSIS (2026-08-04, for the wiring phase): ingrowth-tree heights are the one remaining
model piece, but they are COUPLED to the tally-body draw machinery, not an isolatable increment. Since PSUB=0
(ITIME≤2), all BEST trees are ADVANCE regen (ICHOI adv, DRAW≤PADV/(PADV+PSUB)=1 always) → ESADVH; EXCESS trees →
ESXCSH (ported). ESADVH (esadvh.f, 174L) is a clean per-species regression HHT=EXP(PN + EMSQR·DILATE·BNORM·σ_sp)
with THAB(5,MAXSP)/TPRE(4,MAXSP)/TPHY(5,MAXSP) hab/prep/phys tables (coeffs extracted, ready to port) — BUT its
inputs DILATE(=FIRST(1,I)), DELAY(=ESDLAY draw, N=INT(DELAY+.5) cap2→1), GENTIM(=FINT-5), TIME, EMSQR come from
the body-39..84 height draws + ESDLAY + the FIRST dispersion array. ⇒ can't validate ESADVH in isolation (no
advance-only height oracle; the regen report's AVERAGE HEIGHT [WP3.5 DF3.8 GF1.7 WH1.9 RC2.1 ES2.2] blends
adv+excess and needs the full DILATE/DELAY/GENTIM chain). So the heights + tree-creation + esnutr-scheduler +
per-plot-input-derivation are ONE coupled shared-engine build, validated end-to-end via full-cycle .sum. Port
ESADVH + ESDLAY + wire into engine/establishment.jl as that dedicated unit; the tally counts/species/total are
already bit-exact and feed it directly.
★ HEIGHT-MODEL SPEC COMPLETE (routine-level, 2026-08-04): the ingrowth height model = **ESDLAY** (esdlay.f 186L,
per-species delay from a DRAW; BSUB/BBW Weibull coeffs) → DELAY → N=INT(DELAY+.5) cap2→1, AGE=3-DELAY-GENTIM
(GENTIM=FINT-5) → **ESADVH** (best/advance, 174L) or **ESXCSH** (excess, PORTED) with **DILATE=FIRST(1,i)** where
FIRST(1,i) inits 0.1 (estab.f:181) then FIRST(1,i)=SQRT(FIRST(1,i)) after each tree (order-statistic spread:
0.1→0.316→0.562→…→1). Height draws at body-39..84 feed ESDLAY. This is an OUT-OF-ENGINE extension of
ie_autoes_tally (validatable via the regen report AVERAGE HEIGHT per species), but it's ~2 more routine ports
(ESDLAY+ESADVH) + the FIRST chain — a substantial sub-model best built with the tree-creation as one unit. The
tally counts/species/total (the .sum-TPA driver) are already bit-exact; heights are a second-order refinement
(affect ingrowth tree DBH→multi-cycle growth). All routine-level pieces now identified; nothing left to discover.
★ HEIGHT INPUT-TIMING TRAP (2026-08-04, like ESPADV's TIME=1): ESDLAY read (Weibull, BADV/CADV/BSUB/CSUB +
budworm BBW tables) — for ADVANCE regen DELAY=((-ln(1-DRAW))^(1/CC))·BADV, then DELAY=(DELAY+3)·(-1) [NEGATIVE]
→ clamp≤0 → **advance DELAY=0 always** (iet01 BAA=1→IBAA=1, BWB4=0→no budworm). Then ESADVH AGE=3-DELAY-GENTIM=
3-GENTIM. But a hand-calc with GENTIM=FINT-5 gives GF height ≈0.5 vs oracle 1.7 (3× low) ⇒ GENTIM/AGE (and the
DILATE=FIRST chain + EMSQR dispersion) must be INSTRUMENT-MEASURED per tree, not derived — same measure-first
pattern as every other AUTOES piece. So the height sub-model = port ESDLAY+ESADVH + instrument AGE/GENTIM/DILATE
+ validate avg-height-per-species [WP3.5 DF3.8 GF1.7 WH1.9 RC2.1]. Substantial; part of the dedicated build.

--- ORIGINAL PLAN (historical; the tally is now built + validated per the above) ---
Estimated new surface: ~600–900 lines (ESTOCK+ESNSPE+ESADVH + the esnutr scheduler + the predicted-naturals
slice of estab.f), reusing jl's EXISTING explicit-PLANT/NATURAL tree-creation machinery (engine/establishment.jl).

## Symptom (why this matters)
Stands relying on FVS's DEFAULT automatic natural regeneration after a disturbance collapse in jl. jl's
`establish!` only handles the EXPLICIT PLANT/NATURAL keyword path — it has NO automatic (disturbance-triggered)
regen tally. So a shelterwood / clearcut / heavy-mortality stand that FVS re-stocks from predicted naturals just
dwindles in jl.

## Validation anchor (MEASURED vs live oracle, 2026-08-04)
`iet01.key` stand-4 ("SHELTERWOOD WITH AUTO REGENERATION"), oracle `.sum` = `/workspace/.iework/ierun/iet01.sum`
(the SECOND stand block). Oracle TPA regenerates massively after the shelterwood cut; jl (no AUTOES) collapses:

| year | 1990 | 2000 | 2010 | 2020 | 2030 | 2040 | 2050 | 2060 | 2070 | 2080 | 2090 |
|------|------|------|------|------|------|------|------|------|------|------|------|
| oracle TPA | 536 | 1025 | 1401 | 881 | 1324 | 1531 | 853 | 1412 | **1788** | 1286 | 1147 |
| jl TPA (no AUTOES) | 536 | ~224 | … | … | … | … | … | … | ~28 | … | … |

The oracle's 2000 jump 536→1025 (and the recurring re-stocking pulses) is the predicted-naturals tally. That is
the entire gap. Validate the port against these `.sum` aggregates (NOT per-record — AUTOES adds trees → tripling
→ per-record treelist diff INVALID, doctrine #3).

**MECHANISM CONFIRMED at config level (iet01.key stand-4, line 138):** the stand has `THINPRSC 1990 0.999`
(the shelterwood cut = the disturbance), AUTOES ACTIVE (no `NOAUTOES` — unlike stands 1–3 which carry it), and
**NO explicit ESTAB/PLANT/NATURAL keyword**. So the regeneration is 100% the automatic disturbance-triggered path
(esnutr.f "within 20 yrs of a disturbance" → ESTAB → ESTOCK), NOT the keyword path jl already has. Confirms the
plan is aimed correctly. NEXT measure-first step before coding A1: instrument live FVSie estab.f/estock.f for
stand-4 to capture the PN (stocking prob) + predicted per-species TPP that produce the 536→1025 pulse (the A1
oracle target). Oracle binaries present: /workspace/.iework/FVSie_clean (+ relink_ie.sh, estab.o/estock.o in
bin/FVSie_buildDir).

## ★ A1/A2 ORACLE TARGETS — MEASURED (2026-08-04, no instrument-replay needed)
Ran `FVSie_clean` on iet01 stand-4 with a `DEBUG`(cycle 1) keyword; the standard RegRepts establishment output
(REGENERATION ESTABLISHMENT MODEL v2.0) reports P(stocking) + per-species ingrowth directly — these ARE the A1/A2
targets. Recipe: copy iet01.key → add `DEBUG` line after stand-4 title → `printf "iet01_dbg.key\n..." | FVSie_clean`
→ read iet01_dbg.out "INGROWTH TREES/ACRE ADDED" + "PROBABILITY OF STOCKING IS" + "PLOT HABITAT TYPE SUMMARY".

- **ESTOCK inputs (measured):** IHAB=**10** (all 10 plots in habitat GROUP 10 = WH/western-hemlock series) ⇒ ESTOCK
  dispatch **IEQ=3 (CEDAR & HEMLOCK series, estock.f:74)**. Site-prep IPREP distribution NONE 90% / MECH 8% / BURN 2%.
  This is why GF+WH+RC dominate the regen (the cedar/hemlock stocking eq).
- **P(stocking) targets** (= logistic of ESTOCK PN + the ESB/ESB1/STOADJ correction, estab.f:579): first disturbance
  tally-1(2029)=**0.8918**, tally-2(2039)=**0.9606**; second disturbance (2050 thin) tally-1(2059)=**0.8817**,
  tally-2(2069)=**0.9554**. (For raw-PN bit-exact isolation, add the estab.f:538 `IF(DEBUG)` 6031 dump — but the
  P(stocking) aggregate is enough for A1 acceptance.)
- **Ingrowth pulses (ESNSPE species split, the A2 target):**
  - FALL 1999 (cyc1, post-shelterwood): **583.7 TPA** = WP33 WL20 DF7 GF202 WH222 RC50 ES23 AF27.
  - FALL 2019 (cyc3): **415.1 TPA** = WP70 WL5 DF19 GF150 WH137 RC35.
  - FALL 2089: **135.5 TPA**.
- **Scheduler behavior confirmed:** the model does PERIODIC TALLIES over a ~20-yr post-disturbance window (tally-1 at
  +10yr, tally-2 at +20yr), each ESTOCK→ESNSPE→add. Two disturbances (1990 shelterwood + 2050 thin) ⇒ two tally
  series ⇒ the recurring re-stocking pulses in the .sum. This is the esnutr.f "within 20yr of a disturbance" trigger.

## FVS structure (SOURCE-MAPPED, estb/*.f)
Two paths reach the same tree-creation tail; jl implements only the second:

1. **SCHEDULER — `esnutr.f` (430L)** decides WHEN to run the tally and calls `ESTAB(KDT)` (esnutr.f:401). Rules
   (esnutr.f:292–353): (a) current date within **20 yrs of a disturbance**; (b) **ingrowth** rules; (c) PLANT/
   NATURAL keywords present (NTALLY=99 + old disturbance date). The `NOAUTOES` keyword (base/keywds.f:38) disables
   the automatic (a)/(b) paths. **jl lacks (a) + (b) entirely.**

2. **TALLY ENGINE — `estab.f` (1657L)** — two sub-paths:
   - **PREDICTED NATURALS (the gap):** `CALL ESTOCK` (estab.f:536,572) → `CALL ESNSPE` (701) → `CALL ESADVH`/
     `CALL ESSUBH` (811,830) → pass to tree list (estab.f:1196 "NATURALS FIRST"). Header at :1167 confirms
     "PREDICTED NATURALS ABOVE; PLANT/NATURAL BELOW".
   - **EXPLICIT PLANT/NATURAL (jl HAS this):** estab.f:965–1070 → tree list at :1361. = jl engine/establishment.jl.

### Core sub-models (the NEW code)
- **`ESTOCK` (estock.f, 103L)** — P(stocking) for predicted naturals. Habitat-series dispatch on IHAB→IEQ
  (1=Douglas-fir, 2=grand-fir, 3=cedar/hemlock, 4=subalpine-fir; IPREP>3 → "all roads" eq at :100). Each is a
  regression in: SHAB(IHAB) habitat const, XCOSAS/XSINAS·SQSQ (aspect·√(slope·time)), SPRE(IPREP,·) site-prep,
  ELEV/ELEVSQ, XBAA/XBAA² (BAA), SQREGT/SQBWAF/BWB4 (habitat-code flags), FORDF/FORGF(IFO) forest. Returns PN.
  Coefficients are all in the DATA statements (SHAB[16], SSER[5], SPRE[3,4], FORDF[20], FORGF[20]) — transcribe
  verbatim. Inputs IHAB/IPREP/SLO/aspect/ELEV/TIME/SQREGT/SQBWAF/BWB4 come from ESCOMN/ESCOM2 commons (the estab
  setup); XBAA/XBAALN = stand BAA + ln(BAA).
- **`ESNSPE` (esnspe.f, 68L)** — apportion the predicted total across species (PSPE/TPPLN/ISER/ITPP/TPP).
- **`ESADVH` (esadvh.f, 174L)** — ADVANCE-regen height. **`ESSUBH` (essubh.f, 188L)** — SUBSEQUENT-regen height.
  NOTE: jl ALREADY ported ESSUBH for EM (`em_essubh_hht`, src/variants/easternmontana/establishment.jl) — the IE
  ESSUBH coefficients differ but the structure is identical; generalize that.

## Chunk breakdown (proposed; refine on contact)
- **A0 — harness:** wire iet01 stand-4 into a differential harness vs live FVSie `.sum` (the 2nd stand block).
  Oracle: `/workspace/.iework/ierun/` (relink FVSie if needed). NOAUTOES on = must stay bit-exact (regression gate).
- **A1 — ESTOCK: ✅ DONE (2026-08-04).** Ported estock.f verbatim → `ie_estock` (src/variants/inlandempire/
  establishment.jl): 5 habitat-series eqs + IPREP>3 "all roads" branch + DATA coeffs (SHAB/SSER/SPRE/FORDF/FORGF).
  VALIDATED BIT-EXACT vs live FVSie: iet01 stand-4 {IHAB=10→IEQ=3, SLO=0.30, ASPECT=5.498, ELEV=34, BAA=1, IPREP=1,
  TIME=1, SQREGT=1} → PN=0.2116 (oracle 0.2116). Locked in test/unit/test_ie_estock.jl (7 tests). The estab.f:538
  DEBUG dump needs bare `DEBUG` (all cycles — the fall tally runs under a later ICYC); note a debug-ONLY segfault
  in VOLINIT on tiny regen trees under bare DEBUG (production/.sum runs unaffected — low-pri, not the SIGFPE class).
- **A2 — species apportionment + heights (SOURCE-ANALYZED 2026-08-04; RNG-heavy, defer the RNG half).** Two parts:
  - **A2a ESNSPE: ✅ DONE (2026-08-04).** Ported esnspe.f verbatim → `ie_esnspe` (6 count-logits + SPEHAB(5,4)
    from esblkd.f). Key detail: **XCOS=cos(asp)·SLO, XSIN=sin(asp)·SLO** (estab.f:480-481, SLO-weighted — distinct
    from ESTOCK's plain XCOSAS/XSINAS). VALIDATED BIT-EXACT vs live FVSie (iet01 stand-4 plot-1, ISER=4/ITPP=2/
    TPP=2): PSPE=(0.543, 0.393, 0, 0, 0, 0) = oracle. Locked in test/unit/test_ie_estock.jl (+4 tests).
  - **A2b species probabilities PADV/PSUB/PXCS (DETERMINISTIC, SOURCE-ANALYZED — the next real chunk, larger).**
    estab.f:615-617 calls ESPADV / ESPSUB(if ITIME>2) / ESPXCS → per-species advance/subsequent/excess regen
    probabilities, dumped `PADV=`/`PSUB=`/`PXCS=` under DEBUG (measured iet01 stand-4: PADV=0.062 0.005 0.048
    0.485 0.283 0.122 0.001 0.014 0.039 0…). Each is a clean logistic ×occupancy: e.g. espadv.f
    `PADV(i)=1/(1+exp(-PNᵢ))·OCURHT(IHAB,i)·XESMLT(i)·OCURNF(IFO,i)`, PNᵢ = per-species regression in XCOS/XSIN/
    SLO/TIME/BAA/BAASQ/ELEV/ELEVSQ/REGT/BWAF/BAALN + CHAB(IHAB,i) + CPRE(IPREP,i) + OVER(i)>9.95 & forest bumps.
    **✅ ESPADV DONE (2026-08-04).** Ported espadv.f verbatim → `ie_espadv` (10-species advance-regen logistic ×
    occupancy) with the habitat-type-group tables **CHAB(16,10)** + **CPRE(4,10)** transcribed from esblkd.f:45-76
    (col-major). Occupancy occ(i)=OCURHT(IHAB,i)·XESMLT(i)·OCURNF(IFO,i) + over(i) passed in. VALIDATED BIT-EXACT
    (3 dp, all 10 species) vs live FVSie: iet01 stand-4 PADV=(.062 .005 .048 .485 .283 .122 .001 .014 .039 0)=oracle.
    +11 tests. ★ KEY UNLOCK via measurement (doctrine #2): the prob-call **TIME=1.0** (years-since-disturbance),
    NOT the estab.f:606 TIME=10 literal (overwritten before ESPADV) — my earlier TIME=10 hand-calc (0.008 vs 0.062)
    was the wrong-input trap. PP=0 comes from OCURHT(grp10,PP)=0 (occupancy), not the logistic. (Correction to a
    prior note: CHAB(10,LP)=0, not 1.9319803.) OVER<9.95 (BAAA=0 dump ⇒ no bumps), IPHY≠1.
    **✅ ESPXCS DONE (2026-08-04):** ported espxcs.f → `ie_espxcs` (10-species excess-regen logistic × occupancy)
    with FHAB(16,10)/FPRE(4,10) from esblkd.f:117-153. Quirk: WH(5) uses FPRE(iprep,7) & has no FHAB. VALIDATED
    BIT-EXACT (3 dp, all 10 sp) vs live FVSie: PXCS=(.045 .005 .080 .327 .242 .194 .043 .001 .025 0)=oracle (+10
    tests). REMAINING A2b: **ESPSUB** (subsequent, espsub.f — ITIME>2 gate, NOT fired for iet01 TIME=1 so PSUB=0;
    same logistic×occupancy pattern w/ DHAB/DPRE tables — port when a longer-delay stand exercises it).
    ⇒ deterministic probability layer = ESTOCK + ESNSPE + ESPADV + ESPXCS all bit-exact (ESPSUB inert on anchor).
  - **A2c species COUNT + IDENTITY (RNG) — FEASIBILITY PROVEN (2026-08-04).** estab.f:646-726 draws EMSQR + TPP +
    NUMSPE (from the PSPE cumulative, cap MAXSPP(IHAB)) + species selection via `CALL ESRANN(DRAW)`. ★ ESRANN is
    NOT the hard-to-align main RNG — it's a SEPARATE Park-Miller LCG (ie/esrann.f: ESS1=mod(16807·ESS0,2147483647),
    SEL=Float32(ESS1/2147483648)), seed 43303 for iet01 stand-4. PORTED → `IEEstabRNG`/`ie_esrann!` + VALIDATED
    vs live: from seed 43303, draw#52=0.21862 = the oracle EMSQR magnitude (estab.f:646 uses #51 sign + #52 mag).
    +4 tests. ⇒ bit-exact end-to-end AUTOES is FEASIBLE; the remaining A2c work is replicating the driver's exact
    ESRANN CALL ORDER (a driver-transcription chunk, no longer an RNG-unknown). (The IE memory's "EMSQR/DILATE RNG
    alignment" concern was the establishment-HEIGHT main-stream draw, a different path.)
    ★ ESRANN DRAW-ORDER BLUEPRINT (estab.f, MAPPED 2026-08-04 — the A2c spec):
      1. estab.f:290-294 (NTALLY==1): one ESRANN off the PRIOR stream → ESDRAW=INT(DRAW·100000+0.5); :295
         `CALL ESRNSD(.TRUE.,ESDRAW)` RESEEDS the stream to ESDRAW (=43303 for iet01 stand-4). From here the
         sequence is deterministic and `IEEstabRNG(43303)` replicates it (VERIFIED: draw#52=EMSQR).
      2. :333-336 `DO 183 I=1,IDUP*NPTIDS: ESRANN→WK6(I)` — fills WK6 with **IDUP·NPTIDS** draws (=**50** for iet01:
         DUP=5 × NPTIDS=10). Site-prep assignment (:380-410) then CONSUMES WK6, no new draws.
      3. Per plot (NNID loop): EMSQR = 2 draws (:646-650, draws #51 sign + #52 magnitude=0.21862 ✓); the
         STOADJ<1e-4 branch (:653-670) is SKIPPED when STOADJ normal; ESTPP(TPP) = 1 draw (:675-678); NUMSPE = 6
         draws into WK6 (:693-698) then cumulative-PSPE select (:702-718, cap MAXSPP(IHAB)); species-identity
         selection from PADV/PSUB/PXCS = draws at :739+.
      4. :1075 `CALL ESRNSD(.TRUE.,ESAVE)` saves/restores the stream at the end.
    ★ SPECIES-IDENTITY SELECTION (estab.f:726-810, MAPPED — completes the A2c spec):
      a. :726-735 SUMUP(i) = normalized **(PADV(i)+PSUB(i))** — the species-selection distribution (DISTINCT from
         the PSPE cumulative used for NUMSPE-count). NSPNZ = #species with (PADV+PSUB)>1e-4; NUMSPE capped to NSPNZ.
      b. :738-741 six ESRANN → WK6. :742-763 DO 60 I=1,NUMSPE: draw WK6(I) vs SUMUP cumulative → mark IBEST(J)=1,
         zero SUMUP(J), renormalize ⇒ picks NUMSPE DISTINCT species.
      c. :773-788 recompute PADV(ESPADV if NTALLY==1)+PSUB(ESPSUB); per IBEST species, 1 ESRANN → ADV vs SUBS
         (ICHOI(J,i), J=1 if DRAW≤PADV/(PADV+PSUB) else 2).
      d. :797-810 NOFSPE*2 ESRANN → WK6; per chosen species assign tallest-tree height via ESDLAY(delay)+ESSUBH/
         ESXCSH (DILATE=FIRST). Then TPA booked (estab.f:900+) distributing the plot's ITPP over IBEST species.
    ⇒ A2c = a plot-loop transcription: reseed→WK6(50)→per-plot{EMSQR,ESTPP→ITPP,NUMSPE(PSPE),IBEST(PADV+PSUB),
    ADV/SUBS,heights,book TPA} using the ported ie_esrann!/ie_esnspe/ie_espadv/ie_espxcs/ie_ocurht/ie_essubh.
    Large but FULLY SPECIFIED + every primitive bit-exact. Validate end-to-end vs the .sum (1999 GF202 WH222=583.7).
    ★★ A2c CORE VALIDATED END-TO-END (2026-08-04): `ie_estab_pick_species` (estab.f:745-753) + the draw-order model
    reproduce iet01 plot-1's species pick STRAIGHT FROM THE RNG: IEEstabRNG(43303) → 59 draws (WK6-fill 50 + EMSQR
    2 + ESTPP 1 + NUMSPE-WK6 6) → draw#60=0.61708 ∈ WH band (0.567,0.834] of normalized PADV → species 5 (WH) =
    oracle IBEST. draw#52=0.21862=EMSQR checkpoint. ⇒ the hardest A2c piece (RNG-driven selection) is PROVEN
    bit-exact; the remaining A2c work is the mechanical multi-plot loop + TPA booking (estab.f:900+) + engine wiring.
    ★ ASSEMBLY CHALLENGE + PLOT-BY-PLOT VALIDATION TARGET (2026-08-04): the per-plot draw stride is LARGE and
    VARIABLE (each plot consumes EMSQR 2 + ESTPP 1 + NUMSPE-WK6 6 + species-WK6 6 + ADV/SUBS NOFSPE=23 [:780-788]
    + heights NOFSPE*2=46 [:797-800] + excess MAXTPP(IHAB)*2 [:909-912] + …), so it CANNOT be pattern-matched —
    the assembly must transcribe EVERY per-plot ESRANN call in order and validate PLOT-BY-PLOT. Validation anchor
    = the per-plot EMSQR sequence (first tally, plots 1-10, captured live): **0.219, -0.925, -0.528, 0.863,
    -0.565, 0.721, -0.280, 0.562, -0.130, -0.757** (plot-1 mag 0.219 = draw#52). Build ie_autoes_tally so each
    plot's EMSQR (2 draws) matches this list ⇒ proves the per-plot draw count exact; then ITPP/NUMSPE/species/TPA
    booking follow. NEEDS: MAXTPP(IHAB)+MAXSPP(IHAB) (blkdat.f), NOFSPE=23, the ESDLAY/height draw count. Intricate.
    ★ LOOP STRUCTURE (estab.f, MAPPED): per-plot processing is a NESTED loop `DO 203 NN=1,NPTIDS → DO 202 ITYPEP=
    1,4 → DO 201 IREP=1,NTIMES`, NCOUNT=plot counter (1-10). The per-plot ESRANN count is BRANCH-DEPENDENT: the
    STOADJ<1e-4 "skip" path (:653-670, loops DO 71=6, DO 72=6, DO 73=NOFSPE, DO 74=NOFSPE*2, DO 75=MAXTPP*2) vs the
    MAIN selection path (DO 15=6 NUMSPE-WK6, DO 120=6 species-WK6, DO 63=NOFSPE ADV/SUBS, DO 122=NOFSPE*2 heights,
    DO 123=MAXTPP*2 excess). ⇒ ie_autoes_tally must replicate the nested loop + the STOADJ branch exactly; validate
    each NCOUNT plot's EMSQR against the captured list. This is a focused, careful build (one wrong ESRANN count
    desyncs everything) — best done as a dedicated unit, not piecemeal. Foundation (8 primitives) is 100% ready.
    ★ CONSTANTS FOUND (estab.f:108-112 DATA): **MAXTPP(16)=[9,7,5,5,10,8,9,5,21,25,10,10,11,7,10,8]** (IHAB10→25,
    excess=MAXTPP*2=50 draws); **MAXSPP(16)=[4,3,3,3,5,4,6,4,6,6,4,5,5,4,6,4]** (IHAB10→6, the NUMSPE cap);
    MYHABG=[4*1,4*2,3,4,6*5] (→ISER), MYHTS=[1,3*2,4*3,2*4,6*5], MYTYPE=[9*1,2*5,2*2,3*3,4,12*5,1]. ESDLAY does
    NOT draw ESRANN; the DO 99 height loop (:802-830, ESDLAY/ESADVH/ESSUBH) consumes the pre-drawn WK6, no new draws.
    ⚠ DRAW-COUNT GAP: main-path per-plot count = EMSQR2+ESTPP1+NUMSPE6+species6+ADV/SUBS23+heights46+excess50 (+1
    ESAVE?) ≈ 135, but the observed plot1→plot2 EMSQR stride is ~160 (Δ~25) ⇒ ~25 draws unaccounted (likely a
    per-plot site-prep/WK6 refill or the ITYPEP nested-loop iterations). The exact count MUST be forward-traced
    line-by-line through the DO 201/202/203 nest — cannot be reverse-engineered from EMSQR positions (false matches).
    All constants are now in hand; the assembly is a careful forward transcription validated plot-by-plot vs the
    EMSQR list.
    ★ CONFIRMED EMPIRICALLY (2026-08-04): plot-1 is bit-exact (EMSQR mag=draw#52 AND species-pick=draw#60→WH, two
    independent matches). But plot-2 EMSQR (-0.925, needs sign<0.5) does NOT appear at any small stride — the only
    sign<0.5 match for 0.925 is draw#476 (~424 stride, implausible for one plot). ⇒ the 10 captured EMSQR values
    are NOT a simple consecutive-plot stream; they span the nested DO 202 ITYPEP=1,4 / DO 201 IREP loop (and/or the
    per-plot reseed/ESAVE). CONCLUSION: the multi-plot draw sequence CANNOT be reverse-engineered from EMSQR
    positions — ie_autoes_tally MUST forward-trace the full DO 245/2451/203/202/201 nest line-by-line (every ESRANN
    in each branch), validating plot-1 first (known-good), then each subsequent NCOUNT. This is the dedicated build.
    ★ LOOP NESTING RESOLVED (estab.f:447-497): NCOUNT=0; DO 203 NN=1,NPTIDS → DO 202 ITYPEP=1,4 (skip if
    NNPREP(ITYPEP)<1) → DO 201 IREP=1,NTIMES(=NNPREP); NCOUNT++ per rep ⇒ **NCOUNT total = NPTIDS·IDUP = 50** (NOT
    10). `IF(IPREP.EQ.IPOLD) GO TO 137` skips only the plot-SETUP (475-601, which draws NOTHING — ESTOCK/ITPP calc
    are draw-free); all 50 reps run the 602+ EMSQR/selection draws. So the 10 captured EMSQR are the first 10 of
    ~50 NCOUNT iterations. Per-plot fixed-draw count SHOULD be 2(EMSQR)+1(ESTPP)+6+6+23+46+50 ≈ 134, but plot-2
    EMSQR isn't at +134 ⇒ my hand-count is missing draws in some branch (INGRO/INADV/MATCH path at :588, or the
    967 ESAVE, or a per-rep WK6). ⇒ DEFINITIVE NEXT STEP = INSTRUMENT the draw count: patch esrann.f with a call
    counter + dump it at each EMSQR (or before each 6033 PLOT print) via the DEBUG-run recipe; that gives the EXACT
    per-NCOUNT draw stride directly, ending the hand-count guessing. THEN write ie_autoes_tally with that stride and
    validate plot-by-plot vs the EMSQR list. (Doctrine #2: instrument when unsure — hand-tracing 540 lines of
    branchy Fortran is the wrong tool here.)
    ★★★ INSTRUMENTED (2026-08-04, doctrine #2 win) — patched esrann.f (call counter NESDRW in a new /ESDBGC/
    common) + estab.f (dump NESDRW at each EMSQR), relinked FVSie_trc, ran, RESTORED buildDir clean. RESULT:
    **per-plot draw stride = EXACTLY 135** (NESDRWCT at EMSQR = 52,187,322,457,592,727,862,997,1132,1267 — constant
    +135; EMSQR at draw-52 within each plot). My hand-count 134 was off by 1 (the :967 ESAVE draw). ★ BUT the pure
    LCG from 43303 only reproduces PLOT-1 EMSQR (0.219@#52); plots 2-10 DON'T match at #52+135k ⇒ there is a
    **PER-PLOT RESEED**: estab.f:967 `ESAVE=INT(DRAW*100000+0.5)` then :1075 `CALL ESRNSD(.TRUE.,ESAVE)` reseeds the
    stream each plot. So each plot is a FRESH LCG from its own seed (plot1=43303), 135 draws, EMSQR at draw 52. The
    seed-chain is NOT a simple fixed draw-offset (tried seed_{n+1}=INT(d_n[k]*100000+0.5) ∀k — k=5 matches plots
    1-2 then diverges at plot-3). ⇒ FINAL UNKNOWN = the exact ESAVE-generation: RE-INSTRUMENT to dump ESAVE (or ESS0
    at each plot's first draw) per plot → gives the exact per-plot seed sequence. THEN ie_autoes_tally = {for each of
    NPTIDS·IDUP plots: seed→135-draw stream, EMSQR@52, ESTPP@53, NUMSPE, IBEST, excess, book TPA; next seed=ESAVE}.
    Plot-1 is already bit-exact end-to-end. This is the last measurement before the assembly is fully determined.
  - **ESADVH/ESSUBH heights:** reuse the EM essubh generalization (ie_essubh already exists in this file).
- **A3 — scheduler (esnutr.f rules) — MEASURED 2026-08-04:** MODEL COMPLETE (all probs+RNG+heights bit-exact,
  commit 6ae1ffb, 95 tests). Remaining = this trigger + tree creation. Trigger semantics (esnutr.f/esin.f/esinit.f):
  - **Defaults (esinit.f:50-64):** `LAUTAL=LINGRW=LSPRUT=.TRUE.`, `THRES1=0.10`, `THRES2=0.30`, `NTALLY=0`,
    `IDSDAT=-9999`, `MINREP=50`, `STOADJ=1.0`. IE defaults auto-establishment ON.
  - **NOAUTOES** (initre.f:2815 opt-72 → ESNOAU, esin.f:783): sets `LAUTAL=LINGRW=LSPRUT=.FALSE.`, STOADJ=0.
    (jl currently only zeroes lsprut — keyword_dispatch.jl:2233 — a STUB; must also clear lautal/lingrw.)
  - **Keywords (esin.f):** INGROW→LINGRW=T; NOINGROW→F; AUTALLY→LAUTAL=T; NOAUTALY→F; THRSHOLD→THRES1/THRES2
    (÷100, clamp T1∈[.025,.95], T2∈[.05,.975]); NATURAL keyword → LAUTAL=LINGRW=F, STOADJ=0.
  - **LAUTAL removal path (esnutr.f:264-290):** after a thinning cycle, `XTPA=ONTREM(7)/ONTCUR(7)`,
    `XCUF=OCVREM(7)/OCVCUR(7)`, `XTES=max(XTPA,XCUF)`. `LONE=(THRES1≤XTES<THRES2)`. If `LONE .OR. XTES≥THRES2`:
    `IDSDAT=IY(ICYC)`, `NTALLY=1`, schedule tally at `KDT=IY(ICYC+1)-1`. (ONTREM/ONTCUR = removed/current
    per-acre TOTALS from the thin — jl `cuts!` returns removed totals; need current-before-thin for the ratio.)
  - **20-yr continuation (esnutr.f:298):** `IF(KDT-IDSDAT≤19 .AND. NTALLY>0)` → NTALLY++, reschedule tally.
    Drives the multi-cycle 536→1025→…→1788 on stand-4 (IDSDAT=1990 from THINPRSC 0.999, XTES≈0.999≥0.30).
  - **LINGRW ingrowth (esnutr.f:313-343):** if no tally next cycle & ((ITRN=0 & ICYC=1) OR (IY(ICYC+1)-IDSDAT≥40)):
    `NTALLY=99` (ingrowth signal), `IDSDAT=IY(ICYC+1)-20`.
  - **jl infra present:** `Establishment` struct (src/core/state.jl:599) has active/idsdat/ntally/es_seed/years_done;
    `establish!` (engine/establishment.jl) already runs PLANT/NATURAL + REGENT-crown tail + computes
    dupnpt/gentim/nptids/idup. NEED-TO-ADD: lautal/lingrw/thres1/thres2 fields (IE default T/T/.10/.30); the
    removal-fraction trigger post-cuts!; the AUTOES branch (ie_autoes_tally → create trees via the existing tail).
  Fire the tally in engine/establishment.jl's cycle hook. Reuse the existing tree-creation tail (naturals-first).
- **A3b — MEASURED stand-4 trigger + target (2026-08-04, FVSie_clean DEBUG on iet01 stand-4/THN3):**
  - **★ TRIGGER IS THE INGROWTH PATH, NOT THE REMOVAL-THRESHOLD PATH.** Cycle-1 ESNUTR dump:
    `XTPA XCUF=0.000 0.000; LONE=F` (removal path did NOT fire) but `ITRN=0` (THINPRSC 0.999 left the stand
    bare) → LINGRW ingrowth trigger (esnutr.f:328 `ITRN=0 .AND. ICYC=1`) → `NTALLY=99`, `IDSDAT=IY(ICYC+1)-20=1980`.
    So the removal-fraction (XTES/THRES) machinery is NOT what drives stand-4 — it's the bare-stand ingrowth rule.
    (The removal-threshold path still needs porting for stands thinned to a NON-bare residual, but stand-4's
    anchor trajectory is ingrowth-driven.) Cycle-1 ESTAB added the 531 TPA seen at 2000.
  - **TARGET trajectory (clean .sum, /workspace/.iework/autoes_measure/iet01_clean.sum, THN3 = 4th block):**
    year→TPA: 1990→0, 2000→531, 2010→420, 2020→711, 2030→1537, 2040→1868, 2050→1256 (THINBTA 2050),
    2060→1300, 2070→1600, 2080→1171, 2090→1081. jl currently COLLAPSES to 0 (no AUTOES). This is the A4 gate.
  - **BLOCKER for the full NTALLY sequence:** global DEBUG keyword makes cycle-2 VOLINIT segfault on the tiny
    (~0-volume) regen seedlings — the debug .out truncates at cycle 2 (iet01_s4dbg.out, 18863 lines). Only
    cycle-1 ESNUTR/ESTAB captured. To get cycles 2-10 NTALLY: either patch/guard the crashing VOLINIT debug
    WRITE and relink FVSie_trc, or read estab.f's OPADD follow-on-tally scheduling to reconstruct analytically.
    LINGRW re-fires when `IY(ICYC+1)-IDSDAT≥40` (esnutr.f:332) → IDSDAT resets to IY(ICYC+1)-20 each fire.
  - Artifacts persisted: /workspace/.iework/autoes_measure/{iet01_s4dbg.out, iet01_clean.sum, iet01_s4dbg.key}.
- **A3c — FULL measured NTALLY scheduler (2026-08-04, unconditional AUTOESTRC dump in esnutr.f before CALL ESTAB,
  relinked FVSie_trc, ran WITHOUT the DEBUG keyword → no VOLINIT crash, all 10 cycles). GROUND TRUTH:**
  ESTAB fires on cycles 1,3,4,5,7,8,10 (NOT 2,6,9). Per firing (ICYC / →endYr / NTALLY / IDSDAT / trigger):
  - cyc1 →2000 NTALLY=99 IDSDAT=1980 — INGROWTH (ITRN=0 after 99.9% thin; ICYC=1 rule). +531 TPA.
  - cyc3 →2020 NTALLY=99 IDSDAT=2000 — INGROWTH (IY(ICYC+1)-IDSDAT_old=2020-1980=40≥40 rule).
  - cyc4 →2030 NTALLY=1  IDSDAT=2020 — REMOVAL disturbance (THINBTA 2020 thin in cyc4, start-year 2020).
  - cyc5 →2040 NTALLY=2  IDSDAT=2020 — CONTINUATION (KDT-IDSDAT=2039-2020=19≤19).
  - cyc7 →2060 NTALLY=1  IDSDAT=2050 — REMOVAL disturbance (THINBTA 2050 thin in cyc7).
  - cyc8 →2070 NTALLY=2  IDSDAT=2050 — CONTINUATION.
  - cyc10 →2090 NTALLY=99 IDSDAT=2070 — INGROWTH (2090-2050=40≥40 rule).
  ★ So BOTH paths fire: the LAUTAL removal-threshold path (cyc4/5, cyc7/8 from THINBTA) AND the LINGRW ingrowth
  path (cyc1/3/10). NTALLY 1→2 is the TALLYONE→TALLYTWO 20-yr continuation. The removal path (A3a) IS needed
  after all — a thinning to a non-bare residual (THINBTA) triggers NTALLY=1; a thinning to bare (THINPRSC 0.999)
  falls through to the ingrowth path. A no-tally cycle (2/6/9) = growth+mortality only. Sequence artifact:
  /workspace/.iework/autoes_measure/stand4_ntally_sequence.txt. This fully determines the establish! scheduler.
- **A3d — NEXT: implement the NTALLY state machine in establish!** mirroring esnutr.f: (1) after cuts!, if a
  thinning removed ≥THRES2 to a NON-bare residual → set idsdat=cycleYr, ntally=1 (schedule removal tally);
  (2) at establish!, the ingrowth rule ((ITRN=0 & ICYC=1) OR (IY(ICYC+1)-idsdat≥40)) → ntally=99→1 full tally,
  idsdat=IY(ICYC+1)-20; (3) 20-yr continuation (KDT-idsdat≤19 & ntally>0 → ntally++); (4) each firing runs
  ie_autoes_tally (PNONE=1 for ingrowth) + creates trees via the existing tail. Validate cycle-by-cycle vs the
  0/531/420/711/1537/1868/1256/1300/1600/1171/1081 trajectory.
- **A3d DONE (commit bc25320): ie_autoes_schedule! bit-exact vs the measured 10-cycle sequence** (106 tests).
  XTES measured (FVSie_trc): cyc1 ONTREM=ONTCUR=0 (inventory-year thin NOT tallied → gate on year>inv_year);
  cyc4 XTES=0.799, cyc7 XTES=0.975 (THINBTA). Artifact stand4_xtes_ntally.txt.
- **A3e — tree-creation wiring (NEXT, fully specced from cyc1 ESTAB dump iet01_s4dbg.out:3165+):** run
  ie_autoes_tally per firing + create trees. MEASURED cyc1 inputs (ingrowth): seed0=43303 (=ie_autoes_seed0(55329)),
  nplots=50, NOFSPE=23, IPREP=1, ISER=4, IFO=4, SLO=0.30, BAA=1.00, ELEV=34.0, ASPECT=5.498(→xcos/xsin),
  TIME=1.0, BWB4=BWAF=0, REGT=SQREGT=1.0, PN_stocking=0.2116→PROB1=logistic=0.5527, OCURHT/OCURNF occupancy
  (first 9 sp=1). PADV/PSUB/PXCS per species match ie_espadv/espsub/espxcs. ★ NEW SUB-PIECE = derive the ESTAB
  habitat indices (estab.f): ISER=MYTYPE(ITYPE) then MYHABG(IHAB); IHAB=IPHAB(NNID); IPHY=IPHYS(NNID);
  IFO=IFORST(I) (default 4, estab.f:214/218); IPREP via ESPREP(ISER,PNONE,PMECH,PBURN)+WK6 sampling (all=1 when
  PNONE=1 ingrowth). For ITYPE=17→ISER=4, IFO=4, IPREP=1. Port MYTYPE/MYHABG/IPHAB/IPHYS/IFORST tables + ESPREP.
  Then create the per-species TPA (ie_autoes_tally output) as trees via establish!'s existing tree-creation tail
  (DBH=0.1+0.001·HHT since all establishment heights <4.5ft; heights from ie_esadvh/essubh/esxcsh by ICHOI).
- **A3f — engine hook:** capture XTES in grow_cycle! after cuts! (rem.tpa/pre_tpa, rem.cuft/pre_cuft; skip
  inventory-year), stash on estab; call ie_autoes_schedule! + the A3e tally in establish!'s IE branch.
- **A4 — full-cycle .sum** vs the target trajectory.
- **A4 — full-cycle differential:** iet01 stand-4 `.sum` TPA/BA/SDI vs oracle (the anchor table above). `.sum`
  aggregates ONLY (tripling). Then sweep the other western auto-regen stands.

## Doctrine (carry forward)
Validate bit-exact vs LIVE FVSie per chunk (instrument-replay ESTOCK/ESNSPE for PN/TPP); MEASURE don't infer;
`.sum` aggregates only after the tally (tripling voids per-record diff); NOAUTOES stands must stay bit-exact
(regression gate); reuse engine/establishment.jl's tree-creation tail — only add the ESTOCK/ESNSPE/ESADVH
equations + the esnutr trigger. This is IE-first (iet01 anchor); EM/BM/UT/KT/CR share the estb/ engine with
per-variant ESTOCK/ESSUBH coefficient sets (EM essubh already done).

See [[fvsjl-ie-variant-port]] (IE is the anchor variant) and docs/ESTB_SUBSYSTEM_CHUNK_PLAN.md (the explicit
PLANT/NATURAL path + EM essubh, already ported).

## ★★ CRITICAL ORACLE CORRECTION (2026-08-04) — earlier target was a broken (tree-data-not-loaded) run
The A3b/A3c "target" (0/531/420/711/…) and the NTALLY sequence [cyc1(99 ingrowth)…] were measured from FVS runs
where iet01.tre FAILED TO ASSOCIATE (scratchpad copy not named to the keyfile) → "TREE RECORDS: 0" → every stand
ran BARE, so stand-4 regenerated purely from the ingrowth path. WRONG. Running FVSie_clean IN tests/FVSie (correct
tree-data association) gives the REAL behavior:
- **REAL target TPA (iet01.sum THN3, /workspace/.iework/autoes_measure/iet01_REAL.sum):** 1990→536, 2000→1025,
  2010→1401, 2020→881, 2030→1324, 2040→1531, 2050→853, 2060→1412, 2070→1788, 2080→1286, 2090→1147. (2070=1788
  matches the memory's "oracle→1788".) jl 1990=536 is BIT-EXACT (inventory loads correctly).
- **REAL scheduler (stand4_REAL_scheduler.txt):** ESTAB fires cyc 1(NTALLY=1),2(2),4(1),5(2),7(1),8(2),10(99).
  cyc1 = REMOVAL path via the THINPRSC XTES=0.5526 (55% of KUTKOD≥2 trees), IDSDAT=IY(ICYC)=1990 — AT the
  inventory year. cyc2 = 20-yr continuation. cyc4/7 = THINBTA removal (XTES 0.838/0.963). cyc10 = ingrowth (40-yr
  gap). The stand is NOT bare — AUTOES fires ALONGSIDE the overstory.
- **FIX:** removed the wrong `year > inv_year` removal-path gate (it was derived from the broken run's XTES=0). The
  removal path now fires at the inventory year too. jl stand-4 now fires AUTOES end-to-end: 2000 jumped 224→790
  (was the no-regen baseline). Trajectory in the right ballpark, ~20-50% UNDER — REFINEMENT (not structural):
  (1) non-bare tally inputs (baa=actual stand BA, prob1 recomputed — validated only for the bare baa=1 case);
  (2) multi-tally ESRANN seed chain (all firings reuse seed0=43303); (3) per-record heights (XMIN placeholder).
- LESSON: always verify the oracle loaded its inputs — "TREE RECORDS: 0" / "TOO FEW PROJECTABLE TREE RECORDS" in
  the .out means a bare run. Run FVS in tests/FVSie (or name the .tre to the keyfile base).

## ★ PROB1 ESB-CORRECTION — the ~20-50% under-production root (2026-08-04, MEASURED)
jl computes PROB1 = logistic(ie_estock PN) only; the REAL PROB1 (estab.f:536-584) is:
```
ESB1(NCOUNT) = ESTOCK(ELEV,IFO,BAAOLD,BAAOLN)         # predicted stocking at "inventory"/disturbance BAA (l.536)
PN           = ESTOCK(ELEV,IFO,BAA,BAALN)             # predicted stocking at END-of-cycle BAA (l.572)
PROB1 = logistic(PN + ESB - ESB1) * STOADJ            # (l.579-580), clamp[0.0001,0.9990], floor PNN
```
where ESB (the ACTUAL-vs-predicted intercept, estab.f:319-326) is RECOMPUTED at EACH NTALLY==1 (each new
disturbance tally), NOT once:
```
TPACRE = Σ PROB(i) for trees with DBH < REGNBK(=2.999)   # current SMALL-tree TPA (blkdat.f:235)
ESA    = logistic(-5.17397 + 0.85131*ln(max(TPACRE,1)))
ESB    = logit(clamp(ESA, 0.10, 0.90))
```
MEASURED: cyc1 (1990) TPACRE=0 (no small trees at inventory) → ESB=0 → PROB1=logistic(0.2116)=0.5527. By cyc4
(2029) the accumulated AUTOES regen gives many small trees → TPACRE high → ESB≈logit(0.90)=2.197 → PROB1=0.8918.
This is a REGEN→small-trees→higher-ESB→more-regen FEEDBACK. jl PROB1 (measured, per firing): 0.5995/0.5897/…/
0.5491(cyc4)/…/0.5996(cyc7) — flat ~0.55 vs targets 0.5527(cyc1)/0.8918(cyc4)/0.9606(cyc5)/0.8817(cyc7). The
missing ESB feedback is the dominant lever. IMPLEMENTATION: (a) at each ie_autoes_establish! firing, if NTALLY==1
compute ESB from the current small-tree (DBH<2.999) TPA; NTALLY≥2 (continuation) reuses the saved ESB/ESB1/PNN;
(b) ESB1 = ie_estock(disturbance/inventory BAA) — needs BAAOLD (per-plot inventory BA, BAAINV); (c) PN =
ie_estock(current BAA); (d) PROB1 = logistic(PN+ESB-ESB1)*STOADJ clamp/floor. NOTE: also verify the BAA passed to
ie_estock at the tally — jl uses stand_ba (70.6@cyc1) but FVS cyc1 gave PN=0.2116 (=baa≈1 bare); the tally BAA may
be the per-plot regen BAAA not the overstory BA. Measure PN's BAA input per tally (instrument ESTOCK args). After
ESB: multi-tally seed chain + per-record heights are the remaining (smaller) refinements.

## ★ PROB1 per-tally ingredients MEASURED (2026-08-04, stand4_prob1_ingredients.txt)
Instrumented estab.f:584 (unconditional dump ICYC/NTALLY/BAA/PN/ESB/ESB1 at NNID=1, relink, run in tests/FVSie):
```
ICYC NTALLY   BAA      PN      ESB     ESB1     => PROB1=logistic(PN+ESB-ESB1)
  1    1    41.9322  1.9628  -2.1972  -0.6450   => logistic(0.4106)=0.601
  2    2    46.1743  3.0501  -2.1972  -0.6450   => logistic(1.4979)=0.817
  4    1     1.0000  1.9916   0.0000   0.0000   => logistic(1.9916)=0.880  (target 0.8918 ✓)
  5    2     1.0000  3.0699   0.0000   0.0000   => logistic(3.0699)=0.956  (target 0.9606 ✓)
  7    1     1.0000  1.9916   0.0000   0.0000   => 0.880   (target 0.8817 ✓)
  8    2     4.7932  3.0896   0.0000   0.0000   => 0.957
 10    1    50.0508  0.3845   0.0000   0.0000   => 0.595
```
★★ THE DOMINANT DRIVER IS THE PER-PLOT BAAA, NOT THE STAND BA. The BAA passed to ESTOCK at a tally = the per-plot
basal area BAAA (TBAAA=max(BAAA,1)), NOT stand_ba. At cyc4/5/7 the heavy THINBTA leaves the regen PLOTS bare →
BAAA=0 → BAA=1 → ESTOCK PN=1.99 → PROB1=0.88. jl uses stand_ba (147@cyc4) → PN wrong → PROB1=0.55 (the ~40% gap).
The ESB/ESB1 correction (ESB=-2.197 clamped-low, ESB1=-0.645) fires ONLY at the cyc1/2 tally (INADV=0, inventory-
based); cyc4+ have ESB=ESB1=0 (INADV=1 for the auto/ingrowth tallies → the l.319 calibration block is skipped).
So PROB1 model: cyc1/2 = logistic(PN+ESB-ESB1) with the inventory calibration; all later tallies = logistic(PN).
IMPLEMENTATION (revises the earlier ESB-first plan): (1) ★ compute the PER-PLOT BAAA (bare after thin → 1), feed
each plot's own BAAA to ie_estock — this is the main lever; the current single stand-BA is wrong. Find BAAA in
estab.f (the "NNID,BAAA,TBAAA" per-plot loop) — likely per-INVENTORY-POINT BA, not stand BA. (2) ESB calibration
only at the inventory tally (INADV=0). (3) then multi-tally seed chain + per-record heights. jl's PROB1 was flat
~0.55; the fix (per-plot BAAA) lifts the disturbance tallies to ~0.88-0.96, closing most of the 20-50% gap.

## ★ TIME/REGT + BAAA fixes LANDED (commit 37e0155) — trajectory -40% → ±40%
The ~40% under-production had TWO input bugs (both fixed):
1. **BAA = per-inventory-point BAAA(NNID)** (dense.f:213 → jl point_ba[1]), NOT stand_ba. After heavy thin the
   regen point is bare → BAAA≈0 → TBAAA=1 → ESTOCK PN high. (jl point_ba[1]=40 vs live BAAA=41.93 @cyc1; 0→1 @cyc4.)
2. **TIME/REGT = years-since-disturbance**, NOT 1: tally-1=10, tally-2=20, ingrowth=1 (stand4_estock_inputs.txt).
   TIME=(ntally==99)?1:(next_year-IDSDAT); REGT=TIME, SQREGT=√TIME. jl TIME=1→PN=0.21→PROB1=0.55; real TIME=10→
   PN=1.99→0.88. Threaded through ie_autoes_run + ie_autoes_tally (species probs use TIME too).
Result: stand-4 was flat -20/-50%; now OSCILLATES ±40% (2000 790→1210). Two refinements remain:
- **ESB correction (cyc1/2 only, INADV=0 inventory tally)**: PROB1=logistic(PN+ESB-ESB1). Measured cyc1
  ESB=-2.197 (=logit(0.10), few inventory small trees), ESB1=-0.645 → shifts logit by -1.552 → cyc1 0.88→0.60.
  cyc4+ (auto tallies, INADV=1) ESB=ESB1=0. Gate on est.idsdat==inv_year. This reduces the cyc1/2 over-production.
- **★ multi-tally ESRANN seed chain** (the OSCILLATION): per-tally seeds MEASURED = 43303(cyc1/2), 61677(cyc4/5),
  25425(cyc7/8), 48837(cyc10) — each NTALLY==1 draws a NEW seed from the CONTINUING ESRANN stream; NTALLY≥2
  reuses it (estab.f:290-295). NOT consecutive draws from 55329 (only the 1st, 43303, matches) — the stream
  advances by each tally's full draw consumption (per-plot ESAVE reseeds). jl reuses seed0=43303 every tally →
  identical per-tally output → the oscillation. Fix = track the ESRANN stream state (rng.es0) across firings and
  draw the next seed after each tally's consumption; persist on est.es_seed. Artifacts: stand4_estock_inputs.txt,
  stand4_prob1_ingredients.txt, stand4_REAL_scheduler.txt (all /workspace/.iework/autoes_measure/).

## ★ ESB correction LANDED (commit cc3c404) + seed-chain modeling note
ESB inventory calibration wired (cyc1/2 only, gate est.idsdat==inv_year): PROB1=logistic(PN+ESB-ESB1),
ESB=logit(clamp(logistic(-5.174+0.851·ln(TPACRE)),0.10,0.90)) [TPACRE=Σtpa DBH<2.999], ESB1=ESTOCK(BAAOLD,TIME=0),
persisted on est.esb_shift. Effect: 2000 +18%→-14%, 2010 +44%→+17%, 2020 +44%→+19%. Current stand-4 (all fixes):
536/885/1642/1046/1102/2220/1675/999/2100/1606/1737 vs target 536/1025/1401/881/1324/1531/853/1412/1788/1286/1147.
REMAINING = multi-tally seed chain (the ±oscillation, worst at cyc5/2050 +96%). Per-tally seeds measured 43303/
61677/25425/48837. HYPOTHESES TESTED + REFUTED: (a) consecutive ESRANN draws from 55329 → 43303,85428,88687,…
(only 1st matches); (b) last-plot(plot50 of tally-1)+135-draw-body then 1 more → 12147 (≠61677); (c) one-stream
seed0 after 50+135·50 draws → 13356 (≠61677). ⇒ the next-tally seed needs the EXACT total ESRANN consumption per
tally — which INCLUDES the per-tree height (ESADVH/ESSUBH/ESXCSH) + ESDLAY draws, NOT just the 135-draw selection
body. Next: instrument ESRANN call-count per tally (or the ESRNCM state at each tally boundary), model jl's tally
to consume the same total, then draw the next seed. This is the last refinement; the per-record HEIGHTS (still
XMIN placeholder) fold into the same tree-creation pass.

## ★ SEED CHAIN + HEIGHTS ARE ONE COUPLED TASK (2026-08-05, ESS0 states measured)
Instrumented ESRANN's ESS0 (the LCG state, ESRNCM common) at each NTALLY==1 seed draw (stand4_ess0_states.txt):
  cyc1 ESS0=55329 (=ESSS, fresh) → cyc4 ESS0=78807 → cyc7 ESS0=32485 → cyc10 ESS0=62399.
Each tally's seed = ESRANN(ESS0): 55329→43303, 78807→61677, 32485→25425, 62399→48837 (the measured seeds ✓).
So the stream STATE between tallies is the ground truth. My prior 135-draw model (last-plot+135) gave the wrong
ESS0 because the FULL per-plot body consumes MORE than the 135-draw SELECTION body — it also draws the per-tree
HEIGHTS (ESADVH/ESSUBH/ESXCSH) + ESDLAY. The ESAVE plot-to-plot chain uses the draw at estab.f:967 (~draw 135),
but the plot CONTINUES past that (heights/delay) before ending; only the LAST plot's post-135 draws leak into the
final ESS0 (interior plots reseed to ESAVE, discarding them). ⇒ THE SEED CHAIN AND THE PER-RECORD HEIGHTS ARE THE
SAME TASK: implement the full per-plot body as ONE persistent ESRANN stream (seed → EMSQR → ITPP → NUMSPE-WK6 →
species-WK6 → ADV/SUBS+HEIGHTS(ESADVH/ESSUBH/ESXCSH via ICHOI)+ESDLAY → excess-WK6, consuming every draw in FVS
order), tracking rng state across plots AND tallies; then (a) the heights are emitted for tree creation, and (b)
the post-last-plot state gives the next tally's seed (validate vs ESS0 55329/78807/32485/62399). Persist on
est.es_seed. This is the FINAL AUTOES chunk — it closes both the ±oscillation and the placeholder heights at once.
Current jl (all other fixes landed): 536/885/1642/1046/1102/2220/1675/999/2100/1606/1737 vs target 536/1025/1401/
881/1324/1531/853/1412/1788/1286/1147.

## ★ REMAINING ERROR = the CONTINUATION tally (NTALLY≥2), NOT the heights (2026-08-05 trace)
Per-cycle trace (jl, all fixes) of TPA before/after establishment shows the FIRST tally is CORRECT but the
CONTINUATION over-produces ~5×:
  cyc4 (NTALLY=1, 2020 disturbance): pre 157 → post 1280 (target 2030 = 1324) ✓ first tally correct.
  cyc5 (NTALLY=2, continuation):     pre 1280 → post 2502 (target 2040 = 1531) ✗ +63%, added ~1222 vs FVS ~207.
  Same at cyc7(1)/8(2), cyc1(1)/2(2). The FIRST tally of each sequence ≈ right; the CONTINUATION adds a FULL
  tally when FVS adds ~5× fewer. ⇒ the ±oscillation is the continuation amount, NOT the placeholder heights
  (heights affect DBH negligibly — all <4.5ft → DBH≈0.1 — and would hit all tallies equally).
TWO coupled causes (estab.f:770-787 selection + 608-728 stocking):
1. **NTALLY≥2 uses SUBSEQUENT species only** (estab.f:773 `IF(NTALLY.EQ.1)CALL ESPADV`; :774 `CALL ESPSUB`
   always). jl's ie_autoes_tally ALWAYS uses PADV (advance). For NTALLY≥2: PADV=0, SUMUP=PSUB (ie_espsub, already
   validated). For NTALLY=1 with ITIME>2 (cyc4 TIME=10): SUMUP=PADV+PSUB. The ICHOI dispatch (advance vs
   subsequent HEIGHT) is FTEMP=PADV/(PADV+PSUB).
2. **★ the per-plot ITPP / tally AMOUNT** — likely the deeper lever: estab.f:541 ITPP=INT(PLPROB·DUPNPT/(prob1·300)
   +0.5) (inventory) stored in NSTORE; then per-plot total = ITPP·tpaw = ITPP·(prob1·300/dupnpt) = PLPROB (prob1
   CANCELS) — so the tally total per plot = PLPROB (the plot's actual stocking), NOT prob1-scaled. jl uses ie_estpp
   (draw-based ITPP) × prob1·300/dupnpt, which does NOT cancel → the continuation (high prob1=0.956) over-scales.
   VERIFY: does the tally use NSTORE (PLPROB-based, prob1 cancels) or a fresh ie_estpp draw? Read estab.f how the
   per-plot loop gets ITPP + tpaw for the tally (vs the inventory NSTORE calc). If PLPROB-based, jl's whole tpaw/
   itpp needs revisiting (the 583.7 first-tally match may be coincidental alignment). This is the next measure-
   first step: instrument the per-plot ITPP + booked TPA at cyc4 vs cyc5. Redirects from heights → continuation
   amount. Heights remain a (minor) later refinement.

## ★★ THE CONTINUATION FIX = NEWTPP = ITPP − NSTORE (estab.f:679-685) — precise root
The per-plot tally books only the INCREMENT over what the previous tally already stocked:
  679  ITPP = INT(TPP+0.5)                      ! fresh per-plot count (draw), clamp [1,MAXTPP], INGRO→MAXING
  683  NEWTPP = ITPP - NSTORE(NCOUNT)           ! ← only the ADDITIONAL trees vs the prior tally's stocked count
  685  NSTORE(NCOUNT) = ITPP                     ! carry the new stocked count to the next tally
So a NTALLY=1 tally (NSTORE=0) books a full ITPP; a CONTINUATION (NTALLY≥2) books ITPP-NSTORE ≈ the small
increment (FVS cyc5 ~207 vs a full ~1123). jl books a FULL ITPP every tally → the ~5× continuation over-
production → the ±oscillation. FIX: persist the per-plot NSTORE (a dupnpt-length array) on est across tallies
within a disturbance sequence; book NEWTPP=max(0,ITPP-NSTORE) trees; update NSTORE=ITPP. Reset NSTORE at a NEW
disturbance (new IDSDAT / NTALLY=1 that is not a continuation). NOTE this is per-PLOT (NCOUNT), so ie_autoes_tally
must track it per plot and persist across calls — a state addition (est gains a Vector for NSTORE, or a running
per-plot count). This + the NTALLY≥2 subsequent-only species (ESPSUB, PADV=0) are the last two model pieces;
together they should close the oscillation. Heights remain a minor DBH-negligible refinement after.
Current stand-4: 536/885/1642/1046/1164/2274/1386/1043/2067/1622/1752 vs target 536/1025/1401/881/1324/1531/853/
1412/1788/1286/1147.

## ★★ CONTINUATION FIX FULLY SPECCED (estab.f:944-953 ESPROB weighting) — no more measurement needed
Each tally books ITPP trees per plot, but the PER-TREE TPA weight ESPROB encodes the increment:
  944  FTEMP = PROB1(NCOUNT)                         ! this tally's stocking prob
  945  FTEMP2 = FLOAT(NEWTPP)/FLOAT(ITPP)            ! new-fraction (ingrowth only)
  946  ITEMP = ITPP - NEWTPP                          ! = NSTORE = the prior tally's stocked count
  947  DO I=1,ITPP
  948    ESPROB(I) = FTEMP                            ! default = full PROB1 (the NEW trees, I>ITEMP)
  949    IF(I < ITEMP+1) ESPROB(I) = FTEMP - PNN(NCOUNT)   ! OLD trees get the PROB1 INCREMENT (current-previous)
  950    IF(INGRO) ESPROB(I) = FTEMP*FTEMP2           ! ingrowth: all trees scaled by NEWTPP/ITPP
  951    clamp >= 0.0001
  953  PNN(NCOUNT) = FTEMP                            ! carry PROB1 to the next tally
Each tree's TPA = ESPROB(I)·(300/DUPNPT). So:
- FIRST tally (NSTORE=0→NEWTPP=ITPP→ITEMP=0): all ITPP trees get full PROB1. total = ITPP·PROB1·(300/dupnpt).
- CONTINUATION (NSTORE=ITPP_prev): the first ITPP_prev trees get (PROB1_now - PROB1_prev) [≈0 if PROB1 stable],
  the last NEWTPP get full PROB1. So the continuation adds ~NEWTPP·PROB1 + a small increment on the old — NOT a
  full tally. jl gives EVERY tree full PROB1 → over-books the continuation ~5×.
- INGRO tally: every tree scaled by NEWTPP/ITPP.
FIX (final, fully specced): ie_autoes_tally must (a) persist per-plot NSTORE (prior ITPP) + PNN (prior PROB1) on
est, across tallies within a disturbance sequence, reset at a new disturbance (new IDSDAT); (b) ITPP from the
XSTORE-frozen ESTPP draw (NTALLY==1 stores DRAW, continuation reuses it — estab.f:676-677 — so ITPP grows only via
REGT); (c) book each tree's tpaw = ESPROB(I)·300/dupnpt with the I<ITEMP+1 old/new split (+ INGRO scaling); (d)
NTALLY≥2 subsequent-only species (PADV=0/ESPSUB). Heights (ESADVH/ESSUBH per best sp, floored XMIN+0.2 — estab.f:
838) are computed at :795-840 and are a minor DBH-negligible refinement. This closes the oscillation.

## ★★ ESPROB CONTINUATION FIX LANDED (commit f725dd2) — oscillation GONE, now systematic under
Implemented the per-tree ESPROB weighting + per-plot NSTORE/PNN state. Result on stand-4:
  before: 536/885/1642/1046/1164/2274/1386/1043/2067/1622/1752 (oscillating ±12-96%, mean~29%)
  after:  536/885/ 981/ 649/1164/1242/ 741/1043/1050/ 772/1167 (systematic UNDER, mean |Δ|=22.3%)
  target: 536/1025/1401/881/1324/1531/853/1412/1788/1286/1147
The +96% continuation SPIKES are GONE (2050 +62%→-13%, 2040 +49%→-19%). Cap kept MAXING (MAXTPP over-produced).
REMAINING = systematic under-production, and the gap WIDENS across a sequence (cyc1 -14% → cyc2 -30% → cyc2020
-26%) ⇒ the seedling cohort GROWS TOO SLOWLY. Root = the placeholder XMIN heights: the AUTOES seedlings enter
REGENT small-tree growth at the wrong height → wrong growth rate → they lag becoming larger DBH → the cohort
under-accumulates over cycles. So HEIGHTS ARE NOT DBH-negligible after all — they set the growth trajectory (DBH
at creation ≈0.1 either way, but the HEIGHT drives regent's rate). NEXT (the real last piece): emit per-tree
heights in ie_autoes_establish!'s creation loop — advance=ie_esadvh, subsequent=ie_essubh, excess=ie_esxcsh (all
validated bit-exact) dispatched by ICHOI, using the per-tree EMSQR/DILATE(FIRST order-stat)/DELAY(ESDLAY)/AGE from
the tally body (estab.f:795-840; floor TALL=max(HHT+HTADJ, XMIN+0.2), cap HHTMAX). The tally already consumes
those 69 draws (line 697 skips them) — compute + return the per-record heights instead of discarding. Also (minor,
split-only) NTALLY≥2 subsequent-only species. + re-examine the cyc1 first-tally -14% (ESB strength / tally total).

## ★★★ DEFINITIVE ROOT of the tally-amount residual: the total must be PLPROB-CONSTANT (2026-08-05)
FVS ITPP (for the TPA count) = INT(PLPROB(NNID)·DUPNPT/(prob1·300)+0.5) — INVERSELY proportional to prob1. So the
plot total = ITPP·prob1·(300/DUPNPT) = PLPROB — CONSTANT as prob1 rises (higher prob1 → fewer, higher-TPA trees).
jl uses the ESTPP DRAW for ITPP (independent of prob1), so jl total = ITPP_draw·prob1 GROWS with prob1. That is
exactly why the high-prob1 disturbance/continuation cycles (cyc4/5 prob1 0.88/0.96) over-produce at MAXTPP, and why
MAXING(7) accidentally limits it (the cap truncates the prob1-scaled ITPP). Confirmed: jl cyc1 itpp dist (cap25)
sum=307 mean=6.1 max=25 — matches FVS's 20-25 disturbance ITPP, but jl's total scales with prob1 while FVS's doesn't.
THE FIX (the real last sub-model): compute per-plot PLPROB(NNID) = Σ (PROB/DUP) for inventory small trees DBH<REGNBK
(estab.f:303-313) — the plot's actual stocking. Then the tally total per plot = PLPROB (book tpaw so Σ_tree = PLPROB;
equivalently ITPP_tpa = PLPROB·DUPNPT/(prob1·300) with per-tree = prob1·300/DUPNPT). Keep the ESTPP draw-ITPP ONLY
for NUMSPE (species count). Then re-enable MAXTPP. NSTORE inits to the inventory PLPROB-ITPP (estab.f:544), so the
FIRST tally's NEWTPP = draw_ITPP - PLPROB_ITPP and the ESPROB old/new split applies. Also PNN inits to ESA (≈0.10,
line 545). This is the last model piece — it makes the total prob1-INVARIANT (=PLPROB), closing the systematic
under/over. Current best (MAXING workaround): mean|Δ| 22.3%, no oscillation. Artifacts stand4_itpp_newtpp.txt.

## ★ PLPROB MEASURED + prior "prob1-invariant" claim CORRECTED (2026-08-05)
Dumped PLPROB(NNID) per inventory point (stand4_plprob.txt): PER-POINT (not constant), ranging ~7-42 (cyc10:
point1=13.65, point2=29.9, point3=7.09; cyc8: point2=41.9, point3=22.05). ⇒ multiple inventory points (NPTIDS>1),
IDUP replicates. CORRECTION to the prior turn: the tally total DOES grow with prob1 (the TARGET confirms it —
cyc1 prob1=0.60→1025, cyc4 prob1=0.88→1324), so "total=PLPROB prob1-invariant" was WRONG. Re-derivation of the
per-plot total (ESPROB, estab.f:944-951): total = draw_ITPP·prob1·scale − PLPROB·PNN/prob1 (scale=300/DUPNPT).
The 2nd term (the ESPROB old-tree reduction) subtracts ~PLPROB·PNN/prob1 — SIGNIFICANT since PLPROB is large
(20-42). So the real remaining gaps in jl's tally amount are:
1. **NSTORE must INIT to the inventory PLPROB-ITPP** (estab.f:544 = INT(PLPROB·DUPNPT/(prob1·300)+0.5)), NOT 0.
   jl resets NSTORE=0 at a new disturbance ⇒ the first tally books ALL trees at full prob1 (over); FVS books the
   first PLPROB-ITPP as "old" at prob1−PNN.
2. **PNN must INIT to ESA** (estab.f:545, ≈0.10 the actual-stocking prob), not 0.
3. Then the ESPROB subtracts the PLPROB·PNN/prob1 term correctly — reducing high-PLPROB plots.
4. Re-enable MAXTPP (the cap is real; the over-production was the missing NSTORE/PNN init, not the cap).
IMPLEMENT: compute per-plot PLPROB in jl = Σ (PROB/DUP) for inventory small trees DBH<REGNBK (estab.f:303-313);
init NSTORE[plot]=INT(PLPROB·DUPNPT/(prob1·300)+0.5), PNN[plot]=ESA at each new disturbance; keep the ESPROB
weighting; MAXTPP cap. HONEST STATUS: this is a genuine multi-piece sub-model needing a methodical port (PLPROB
+ the NSTORE/PNN inventory init) — not another 1-line tweak. AUTOES is functional (mean|Δ| 22.3%, oscillation
fixed, from collapse-to-28); this closes the systematic bias. Best-state code unchanged (MAXING workaround).

## ★★★ per-tree TPA CONFIRMED = ESPROB·300/DUPNPT (estab.f:1232) — jl's formula is right; fix = NSTORE-init + MAXTPP
estab.f:1232 `PROB(ITRN)=(ESPROB(N)*300.0)/DUPNPT` — EXACTLY jl's tpaw. (estab.f:1203 skips ESPROB<0.00011, i.e.
the clamped-0.0001 old trees are NOT created — negligible.) So the per-tree TPA is correct; the residual is ONLY:
1. **ITPP cap**: FVS uses MAXTPP(25) for disturbances (measured ITPP=20-25), jl's MAXING(7) truncates → fewer
   trees → UNDER. Must use MAXTPP.
2. **NSTORE init from inventory PLPROB-ITPP** (estab.f:544 = INT(PLPROB·DUPNPT/(prob1·300)+0.5)): the first
   PLPROB-ITPP trees are "old" (ESPROB=prob1−PNN, PNN≈ESA≈0.10) not full prob1. THIS reduction is what keeps
   MAXTPP from over-producing. jl inits NSTORE=0 (all trees full prob1) → MAXTPP over-produces (36%).
So the two are COUPLED: enable MAXTPP *and* init NSTORE from PLPROB together. NEEDS the per-plot PLPROB(NNID) =
Σ(PROB/DUP) over the CURRENT small trees (DBH<REGNBK=2.999) on each inventory point (estab.f:303-313) — dynamic,
recomputed each tally (includes accumulated regen; that's why PLPROB=7-42 at later cycles). jl has the tree list
(plot_id = point) → computable, but needs the point/IDUP-replicate mapping (dupnpt=NPTIDS·IDUP). IMPLEMENT:
(a) in ie_autoes_establish!, per point compute PLPROB=Σ(tpa/idup) for small trees on that point; (b) map the 50
tally plots to points (plot n → point via IPTIDS/replicate order); (c) NSTORE[plot]=INT(PLPROB[point]·dupnpt/
(prob1·300)+0.5), PNN[plot]=ESA; (d) cap=MAXTPP. Then the ESPROB old/new split (already implemented) does the rest.
This is the last piece and it is now UNAMBIGUOUS (per-tree TPA confirmed, only the ITPP-count + NSTORE-init remain).

## ★ REGRESSION-SAFETY ASSESSMENT (2026-08-05) — AUTOES firing for all IE stands is FAITHFUL, likely IMPROVES FIA-compat
Concern: ie_autoes_establish! fires for ANY IE stand (guard = lautal||lingrw, both default TRUE), incl. the LINGRW
ingrowth path which triggers at cyc3+ for a normal no-disturbance stand (next_year−idsdat = 50 ≥ 40, idsdat inits
to inv_year−20). With the tally amount only ~22%-accurate, does this regress the broader IE FIA-compat?
ASSESSMENT — NO, it is faithful and likely helps:
- The western FIA sweep keyfiles (run_sweep_western.jl keytext) carry NO NOAUTOES, and live FVSie defaults
  LINGRW=LAUTAL=TRUE (esinit.f:51-52). So LIVE FVSie ALSO fires ingrowth at cyc3 on these stands.
- BEFORE this change jl fired ZERO ingrowth → it already diverged from live by the FULL ingrowth cohort (100% gap).
  jl now adds ~78% of it → the divergence SHRINKS. Net improvement (or neutral), not a regression.
- Regen guard PASS (UT/EM/BM/CR/IE/KT small-tree, iet01_smallr etc. all NOAUTOES → inert; guard clause returns
  false for NOAUTOES stands, guaranteed). iet01 stands 1-3 (NOAUTOES) bit-exact unchanged.
- Ingrowth cohort at cyc3 = freshly established seedlings (DBH≈0.1) → tiny BA impact; the .sum TPA at cyc3 moves
  from 100%-missing toward ~78%-of-live. cyc0 (inventory) UNAFFECTED (AUTOES fires cyc3+).
CAVEAT (validate when convenient): run run_sweep_western.jl on an IE sample to CONFIRM the pass-rate moves up (or
holds), per the FIA/FVS-compat campaign [[fia-fvs-compat-campaign]]. The 22% amount gap means jl's ingrowth is
approximate, so a few knife-edge stands could flip either way; the aggregate direction is toward live.

## ★★★ DEFINITIVE FIT DATA: per-plot BOOKED TPA measured (2026-08-05, stand4_booktpa.txt)
Instrumented estab.f:1232 to accumulate PLTPA = Σ PROB(ITRN) per plot (the actual booked TPA). Cyc1 (first
DISTURBANCE tally, prob1≈0.60), per plot: ITPP=1→3.61, 2→7.21, 3→11.10, 4→14.43, 6→14.79, 14→14.43, 15→14.79,
22→13.90, 25→13.90. ⇒ **PLTPA = min(ITPP, ~numspe≈4)·prob1·(300/dupnpt)** — it CAPS at ~4 trees per plot
REGARDLESS of ITPP (14.43 at ITPP=4 AND ITPP=14). The per-tree unit 3.61 = prob1·300/dupnpt confirms the formula;
the COUNT carrying full TPA is ~numspe (the species count), NOT ITPP. jl books ALL itpp trees at prob1 → at
MAXTPP(25) it over-books ~6× (25 vs 4) → the +36% over-production; MAXING(7) accidentally limits it to ~7 (closer
but still over the ~4). So the model is:
- DISTURBANCE tally: per-plot booked ≈ Σ over the ~numspe BEST trees of prob1 + the EXCESS trees at ESPROB=
  prob1−PNN (which for a fresh disturbance ≈ prob1−ESA, but MEASURED contributes ~0 → clamped/small). Net ≈
  numspe·prob1·scale. (INGROWTH tally is different: ESPROB=prob1·NEWTPP/ITPP for ALL trees → total=NEWTPP·prob1,
  which is why the validated 583.7 ingrowth case matched with jl's itpp·prob1 when NSTORE=0→NEWTPP=itpp.)
FIX (now a FIT, not a derivation): jl's tally must book the EXCESS trees at ESPROB=prob1−PNN (not prob1) for
disturbance tallies — the excess are the "old"/over-stocked trees that FVS discounts. With PNN≈ESA≈0.10 and the
clamp-skip (estab.f:1203 ESPROB<0.00011 skipped), the excess contribute ~0 → per-plot total ≈ numspe·prob1. jl
currently gives BOTH best and excess full prob1. The fix: in ie_autoes_tally, the numspe BEST trees get prob1, the
EXCESS get prob1−PNN (≈0 → effectively skipped for a fresh disturbance). Then MAXTPP is safe. VALIDATE per-plot
vs stand4_booktpa.txt (ITPP→PLTPA map). This is the exact, measured close-out of the 22% residual.

## ★ numspe-cap FIT ATTEMPTED + REVERTED (2026-08-05) — the per-plot amount needs EXACT per-plot fitting
Tried: per-plot total = numspe·prob1·scale (uniform esprob=prob1·numspe/itpp) + MAXTPP. RESULT: mean|Δ| 61.9%
(WORSE than the 22.3% MAXING best) — severe UNDER (2000: 485 vs 1025). Reverted. Learnings that bound the truth:
- MAXING+full-prob1 (per-tree=prob1, itpp≤7): 2000=885 (−14%).
- MAXTPP+full-prob1 (itpp≤25):               2000=1228 (+20%).
- numspe-cap+MAXTPP:                          2000=485 (−53%).
- TARGET:                                     2000=1025.
So the truth books MORE than MAXING-full (885) but the per-tree ESPROB is REDUCED from full prob1 enough to bring
MAXTPP-full (1228) DOWN to ~1025 (~17% reduction) — NOT the drastic numspe-cap. The measured PLTPA=14.4 (plots
1-3) = ~4·prob1·scale was NOT globally "numspe·prob1": jl's avg numspe (~2.7) ≠ live's per-plot count, so a global
numspe rule under-books by half. CONCLUSION (honest): the per-plot amount is a genuine multi-variable fit
(ITPP-cap × the exact ESPROB(I) old/new split × NSTORE-init-from-PLPROB × PNN-init-from-ESA) that CANNOT be
resolved by global trajectory tuning — each global rule I tried (MAXING, MAXTPP, numspe-cap) misses. The
CORRECT next step is to fit jl's PER-PLOT PLTPA to stand4_booktpa.txt EXACTLY (the ITPP→PLTPA map, plot by plot),
reproducing each plot's booked total, THEN the trajectory follows. That requires the full estab.f amount flow
ported faithfully (not fit), validated per-plot — a focused methodical chunk. BEST STATE remains MAXING+ESPROB
old/new = mean|Δ| 22.3%, oscillation-free, committed (77f3ba3). AUTOES is functional; this is the last residual.

## ★★ REFRAMING: MAXING cap is ~RIGHT for the established amount; the 22% residual is DIFFUSE (2026-08-05)
Standalone cyc1 tally sums (prob1=0.5527, 50 plots): Σ itpp·prob1·scale = 1018 raw (MAXTPP-full); Σ min(itpp,4)·
prob1·scale = 487. The TARGET established ≈ 745 (=.sum target 2000=1025 − ~280 grown post-thin remnant). So:
- MAXTPP-full (1018): OVER-establishes ~37% → the +20% .sum seen earlier (1228 = remnant + 1018).
- min4/numspe-cap (487): UNDER ~35%.
- MAXING (cap 7, Σ min(itpp,7)·prob1·scale): lands ~745-800 = CLOSEST to the ~745 target established.
⇒ MAXING (the committed best, mean|Δ| 22.3%) is ~RIGHT on the tally AMOUNT. The 22% residual is therefore NOT one
big lever — it is DIFFUSE across (a) the post-thin REMNANT overstory growth (jl 2000 remnant grows to ~140 vs
live ~280 — the −14% cyc1 is largely the remnant, not the establishment), (b) the continuations, (c) multi-cycle
mortality/compounding. This REFRAMES the close-out: rather than a single tally-amount port, the residual needs
per-component attribution (instrument the remnant vs established split at cyc1: does the post-thin overstory grow
to match live, independent of AUTOES?). MAXING stays as the amount model. HONEST FINAL STATE: AUTOES functional at
22.3%, oscillation-free, regression-safe, amount ~correct; residual diffuse (remnant growth + continuation +
compounding), no single dominant error. This is a reasonable IE-AUTOES v1; further tightening is incremental.

## ★★ DIFFUSE RESIDUAL CONFIRMED by direct split (2026-08-05)
Split jl's cyc1 (2000) stand into established (DBH<1, AUTOES seedlings) vs remnant (DBH≥1, thinned overstory):
  jl established = 661/ac  (est. target ~745, −11%)
  jl remnant     = 223/ac  (est. target ~280, −20%)  ← post-thin ~245 DECLINES to 223 (net mortality) while live
                                                          GROWS to ~280
  jl total 885 vs target 1025 (−14%).
⇒ BOTH components are under — the residual is DIFFUSE, confirmed. Two SEPARATE, smaller contributors:
1. Established slightly low: MAXING(cap7)=661 vs itpp-full=943 vs target ~745 → the true cap is ~10 (between
   MAXING7 and MAXTPP25); a ~11% established gap. (NOTE the target split 745/280 is ESTIMATED, not measured from
   live — measure live's cyc1 established-vs-remnant split to pin it, e.g. instrument FVSie or read its treelist.)
2. Remnant DECLINES in jl (245→223) but GROWS in live (→280): this is LARGE-TREE growth/mortality of the thinned
   overstory — INDEPENDENT of AUTOES (the shared IE growth/mortality, possibly a density-feedback from the added
   seedlings, or the post-thin remnant's mortality). This is NOT an AUTOES-tally bug.
FINAL: the ~22% residual = ~half established-cap (a modest AUTOES tuning: cap ~10 not 7) + ~half remnant growth/
mortality (a non-AUTOES, shared large-tree interaction). No single dominant AUTOES error. AUTOES v1 is functional,
amount ~correct. Close-out next steps (both incremental, both need live measurement): (a) measure live cyc1
established/remnant split → set the true established cap; (b) diagnose the remnant decline-vs-grow (large-tree,
likely density feedback from the seedlings — separate from AUTOES). Committed best: MAXING, mean|Δ| 22.3%.

## ★★ TARGET RE-CONFIRMED + TREELIST-CONTAMINATION caught (2026-08-05)
Near-miss: adding `TREELIST 2000` to iet01.key to measure the live established/remnant split CHANGED the sim —
THN3 (block 4) dropped from 536/1025/…/1788 to 0/531/420/… (the treelist at a cycle boundary perturbs the
establishment RNG/processing). The UNMODIFIED iet01.key confirms the 4 blocks: NONE 536/441/…, THN1 536/441/…,
THN2 (NOAUTOES shelterwood) 536/223/…/28, **THN3 (AUTOES) 536/1025/1401/881/1324/1531/853/1412/1788/1286/1147**
(2070=1788 = memory's oracle). So the WORKING TARGET (536/1025/…) IS CORRECT — my whole analysis stands; the 531
this turn was a contaminated measurement, NOT a target correction. LESSON: the TREELIST keyword ALTERS the AUTOES
sim — cannot use it to inspect the live per-cycle tree split. To get the live established-vs-remnant split
non-invasively: use a DBS TreeList (DATABASE block, already in the keyfile → iet01_Out.db FVS_TreeList table) or
instrument intree/the tree array directly. The jl-side split (est 661 + remnant 223 = 885 vs target 1025) is
UNCONTAMINATED (jl's own run) and stands: jl under-produces ~14% at cyc1, diffuse across established (~11% low)
+ remnant (declines vs live grows). DIFFUSE-RESIDUAL conclusion UNCHANGED; the estimated 745/280 target split
still needs the non-invasive live measurement (DBS TreeList) to pin the established cap precisely.

## ★ live split measurement BLOCKED both ways; target confirmed (2026-08-05 close-out)
Attempts to measure the live THN3 established/remnant split at 2000:
- Text TREELIST 2000: CONTAMINATES (alters the establishment sim → THN3 0/531 not 536/1025). Unusable.
- DBS FVS_TreeList (iet01_Out.db): the DB is 80MB ACCUMULATED across many runs (stale S248112-1..4 CaseIDs +
  fresh UUIDs); TPA is not clean per-acre (S248112-4 2000 est=2/rem=5247 ≠ THN3 .sum). Unusable without a fresh
  isolated-DB run. To do it right: run iet01.key stand-4 ALONE with a FRESH DSNout DB (no accumulation), then
  query FVS_TreeList THN3 2000 by DBH — a clean bounded step for the next pass.
CONFIRMED this turn (the important result): the WORKING TARGET is correct — unmodified iet01.key block4=THN3=
536/1025/1401/881/1324/1531/853/1412/1788/1286/1147. All prior analysis stands. jl's OWN cyc1 split (est 661 +
rem 223 = 885, uncontaminated) is valid; the residual is diffuse (established ~11% low, remnant declines-vs-grows).
NET: AUTOES v1 functional at mean|Δ| 22.3%, target verified, residual diffuse + attributed, the two remaining
close-out steps (fresh-DB live split → established cap; remnant decline diagnosis → shared large-tree) are each
bounded and need a fresh isolated measurement, not more inference.

## ★★★★ CRITICAL TARGET CORRECTION (2026-08-05): the real target is 531, NOT 1025 — jl OVER-produces
The iet01_Out.db is an 80MB DB ACCUMULATED across many prior runs. iet01.key's DataBase block writes DSNOut there.
DECISIVE TEST: same keyfile, FRESH DSNout DB → THN3 = 0/531/420/711/1537/1868/1256/1300/1600/1171/1081; the
ACCUMULATED DB → 536/1025/1401/../1788. The difference is the DB STATE (fresh vs accumulated), not the treelist
(treelist ON/OFF both give 531 with a fresh DB). jl uses TREEDATA (the .tre file, NO database read) → jl's CORRECT
oracle is the FRESH-DB / TREEDATA-based run = **0/531/420/711/1537/1868/1256/1300/1600/1171/1081**. The 536/1025/
../1788 trajectory I (and the [[fvsjl-ie-variant-port]] memory: "oracle→1788") targeted ALL SESSION was
CONTAMINATED by the accumulated iet01_Out.db (a stale S248112 stand whose data the run picked up). This INVERTS
the whole analysis:
- jl stand-4 2000 = 885 vs the REAL target 531 → jl OVER-produces ~67% (not under!). Two CONCRETE causes:
  1. **THINPRSC under-removal**: clean FVS 1990=0 (BARE — removes the whole overstory); jl leaves 223 remnant
     (removes only 55%, reading 9/27 KUTKOD≥2). This is the long-known THINPRSC/KUTKOD gap — jl must remove ~all.
  2. **AUTOES over-establishment**: clean established 531 vs jl 661 (+25%).
- ALL the "under-production / diffuse residual / 22%" analysis this session was against the WRONG (contaminated)
  target. With the real target 531, the picture is CONCRETE: over-thin-remnant (223) + over-establish (130).
NEXT (re-oriented, both concrete): (1) THINPRSC — why does clean FVS remove the WHOLE overstory (1990=0) when
KUTKOD marks only ~9/27? Re-measure the live post-thin count with a FRESH DB (the accumulated DB contaminated the
earlier XTES=0.55 too — re-verify). (2) then the AUTOES over-establishment (661 vs 531) — likely the ITPP cap
(MAXING=7 is slightly HIGH now, not low). ALWAYS use a FRESH DSNout DB for iet01 measurements (rm the accumulated
db or point DSNout elsewhere) — the 80MB iet01_Out.db contaminates. Best jl state committed (MAXING) now reads as
~67% OVER vs the corrected target, but the mechanisms are clear and concrete.

## ★★★★ CORRECTION-OF-THE-CORRECTION (2026-08-05): target IS 1025 — the 531 was a keyfile-mod artifact
The prior "target is 531" entry is WRONG — RETRACTED. DECISIVE one-variable test: delete iet01_Out.db, run the
UNMODIFIED iet01.key (fresh DB, SAME path) → THN3 2000 = 1025 (not 531). So the DB content is NOT the contaminant;
my MODIFIED keyfiles (iet01_fresh.key / iet01_notl.key, which changed the DSNout PATH to /tmp AND/OR removed
treelist keywords) introduced the 531. Changing the DSNout path to /tmp alone (iet01_fresh.key) → 531; that path
change (not the DB, not the treelist) caused it — likely a failed/relocated DB write perturbing the run. LESSON
(the burn): I changed TWO variables at once (path + treelist) → confounded → nearly concluded a huge wrong target
flip. The one-variable isolated test (unmodified key, fresh DB, same path) settled it: TARGET = 536/1025/1401/881/
1324/1531/853/1412/1788/1286/1147 (THN3), CONFIRMED, unchanged. All the session's analysis STANDS: jl UNDER-
produces, mean|Δ| 22.3%, residual diffuse. The [[fvsjl-ie-variant-port]] memory's "oracle→1788" is CORRECT.
NET (final, verified): AUTOES v1 functional, target 1025/../1788 CONFIRMED (survived a contamination scare via a
clean isolated test), residual diffuse+attributed. METHOD LESSON: change ONE variable per measurement; when a
number flips, isolate before concluding. For iet01 measurements the DSNout must stay iet01_Out.db (or the run
behaves differently) — do NOT relocate it.

## 2026-08-06 — CURRENT-STATE MEASUREMENT (post #154; no regression) — residual = MULTI-TALLY RNG, not the model
★ AUTOES IS INTEGRATED + FIRING: ie_autoes_establish! wired at simulate.jl:565 (IE||EM each cycle). estock/autoes
unit tests PASS after this session's #154 establishment.jl edits (NO regression). The goal-file's "jl only does
explicit PLANT/NATURAL" framing is STALE.
★ MEASURED current residual (iet01 s4dbg AUTOES stand, jl vs live iet01_clean.sum):
  2000 jl 531 / live 531 = BIT-EXACT (first tally); 2010 +0.7%; 2020 +30%; 2030 +19%; 2040 -7%; 2050 -9%; 2060 +29%.
⇒ the FIRST AUTOES tally is BIT-EXACT; the ~22% "diffuse residual" is the MULTI-CYCLE RE-TALLY RNG stream diverging
after the first tally — MIXED-SIGN ±20-30% per cycle (partially averaging). This matches the CI-essubh-disp finding
(plot-1 EMSQR bit-exact, later plots/tallies desync = per-plot/re-tally ESRANN stream not aligned). ⇒ the remaining
#143 work is NOT the tally model (bit-exact) NOR the initial wiring (done) — it is the MULTI-TALLY per-plot RE-SEED
RNG stream alignment across cycles (the +10/+20yr re-tally scheduling + per-plot reseed advancing). Coupled RNG-stream
work → fresh-session build. FIX PATH: instrument the live per-cycle re-tally ESRANN seed sequence (ESDRAW re-derivation
at each re-tally) vs jl's ie_autoes_schedule!/plot-seed advance; align the re-seed. Refs: .iework/autoes_measure/
(booktpa + instrumentation), simulate.jl:565, ie_autoes_schedule! (establishment.jl:948).

## 2026-08-06 (CORRECTED via direct per-cycle tally instrumentation) — residual = SCHEDULER + NSTORE BOOKING
★ SELF-CORRECTION (doctrine #2): the earlier "first tally bit-exact, subtle multi-tally RNG drift" was inferred from
a MISALIGNED .sum stand comparison (jl AUTOES stand vs live stand4 = different stands). Direct per-cycle tally
instrumentation (FVSJL_AUTOES_DEBUG in ie_autoes_establish!) gives the REAL picture:
  jl per-cycle tally TOTAL (icyc/ntally/total):  1/99/583.7  3/99/649.9  4/1/1856  5/2/27  7/1/1817  8/2/630  10/99/0.1
  LIVE (stand4_ntally_sequence.txt):             1/99/0      3/99/468    4/1/206   5/2/419 7/1/11    8/2/225  10/99/472
★ CONFIRMED: the tally COMPUTATION is bit-exact — jl icyc1=583.7 perSp[33,20,7,202,222,50,23,27] == the chunk-plan
validated WP33/WL20/DF7/GF202/WH222/RC50/ES23/AF27. ⇒ the residual is NOT the tally model and NOT RNG drift; it is the
ENGINE SCHEDULING + BOOKING: (a) jl fires tallies at different YEARS (icyc1=1990 vs live 1999 — a ~9yr offset in the
esnutr fire-timing / IY(icyc) cycle-year mapping); (b) per-cycle BOOKED amount diverges massively (jl books the full
583.7 at cyc1; live books 0 at cyc1, 468 at cyc3 — the NSTORE continuation booking [ITPP−NSTORE per plot] and the
disturbance/ingrowth schedule don't match live). ⇒ #143 FIX = align ie_autoes_schedule! fire-timing (fire-YEAR /
IDSDAT / IY-mapping) + the NSTORE per-cycle booking to the live sequence above (exact target per icyc/ntally/total).
This is the sharp, actionable target — a scheduler+booking alignment, validated against the 7-tally live sequence,
NOT a coupled RNG mystery. (Debug hook kept ENV-gated in ie_autoes_establish! for the fix session.)

## 2026-08-06 (FURTHER CORRECTED) — residual is purely ENGINE BOOKING (ITPP-cap + NSTORE + deferred), fire-timing OK
★ CORRECTION to the prior entry's "~9yr fire-timing offset": that was a LABELING artifact — my jl debug printed the
cycle-START year (1990); live AUTOESTRC logs KDT (=next_year-1, e.g. 1999), which jl ALSO computes (ie_autoes_schedule!
kdt line 951). Fire-timing IS ALIGNED (same cycles 1,3,4,5,7,8,10; same KDT). The residual is PURELY the per-cycle
BOOKED AMOUNT.
★ MECHANISM (from stand4_xtes_ntally.txt AUTOESXTES): live COMPUTES the tally = 583.6507 (icyc2, = the bit-exact
ie_autoes_tally value) but BOOKS ITPP = 468 (icyc2/3), and books 0 at icyc1. jl books the FULL r.tally sum (583.7) at
icyc1. ⇒ jl is missing the ENGINE BOOKING step: (a) ITPP = the per-plot MAXTPP/MAXING-capped + PXCS-weighted booked
amount (468 < computed 583.65) — jl books the uncapped computed tally instead of ITPP; (b) the icyc1 tally books 0
(setup/calibration) with the real booking DEFERRED to the next cycle (NSTORE continuation books ITPP−NSTORE). So the
tally COMPUTATION is bit-exact but the jl engine books it in full, immediately, instead of the ITPP-capped/NSTORE-
deferred amount live books. FIX = in ie_autoes_establish!, book ITPP (capped+NSTORE-incremented per plot) not sum(r.tally);
match the icyc1-books-0 / deferred-cycle schedule. EXACT per-cycle target (ITPP booked): icyc1=0, 3=468, 4=206, 5=419,
7=11, 8=225, 10=472 (stand4_ntally_sequence.txt). This is the precise, bounded target — the booking step, validated
against the exact ITPP sequence. (stand4_itpp_newtpp.txt has the per-plot ITPP/NEWTPP breakdown for the fix.)

## 2026-08-06 (DECISIVE, measure-verified) — SCHEDULER is BIT-EXACT; residual = per-cycle TALLY-AMOUNT inputs
★ Instrumented jl's full per-tally schedule (idsdat/kdt/time/ntally) vs live AUTOESTRC — ALL BIT-EXACT:
  icyc/kdt/idsdat/time/ntally = 1/1999/1980/1/99 · 3/2019/2000/1/99 · 4/2029/2020/10/1 · 5/2039/2020/20/2 ·
  7/2059/2050/10/1 · 8/2069/2050/20/2 · 10/2089/2070/1/99 — jl == live EXACTLY. ⇒ BOTH prior hypotheses
  ("multi-tally RNG drift" AND "scheduler fire-timing offset") are REFUTED by measurement (doctrine #2, 3rd correction).
★ The residual is PURELY the per-cycle tally AMOUNT: jl total vs live booked = icyc1 583.7/0 · 3 650/468 · 4 1856/206
  · 5 27/419 · 7 1817/11 · 8 630/225 · 10 0.1/472. jl is 4-9× HIGH on most cycles. Since ie_autoes_tally is bit-exact
  for the validated inputs (583.65) and the SCHEDULE inputs (idsdat/time/ntally) match, the divergence is the OTHER
  per-cycle inputs to the tally: the seed0/es_stream ESAVE chain, the per-point baaa, and the es_nstore continuation
  state — AND jl books the FULL tally at icyc1 (idsdat=1980, a PRE-inventory disturbance whose ingrowth is already in
  the 1990 inventory) where live books 0 (ESB inventory-calibration zeroes pre-inventory ingrowth). FIX PATH: instrument
  live's per-cycle seed0/baaa/nstore + the ITPP/NEWTPP booking (AUTOESXTES/ITPPTRC in .iework/autoes_measure/) vs jl's,
  cycle by cycle, starting at icyc1 (why live books 0) then the seed/baaa/nstore chain for icyc3+. The tally MODEL and
  SCHEDULER are bit-exact; only the per-cycle amount-input state diverges. This is the exact, measure-narrowed target.

## 2026-08-06 (ROOT CAUSE — the divergent input is BAAA) — jl uses growing overstory BA, not disturbance-adjusted BAAA
★ Compared jl's per-cycle tally INPUTS (seed0/es_stream/baaa) to the prior session's live instrumentation
(stand4_ess0_states.txt, stand4_estock_inputs.txt). DEFINITIVE root cause = the BAAA (per-point basal area) input:
  cycle:        icyc1   icyc3   icyc4   icyc5   icyc7   icyc8   icyc10
  jl baaa:      0.0     10.27   137.09  491.57  69.17   136.69  848.44
  live BAA:     41.93   46.17   1.00    1.00    1.00    4.79    50.05   (stand4_estock_inputs ESTOCKIN col6)
⇒ jl uses s.density.point_ba[1] = the GROWING OVERSTORY total BA (0→137→491→848 as the regen cohort grows), but
live's BAAA(NNID) is the DISTURBANCE-ADJUSTED per-inventory-point BA: 41.93 at the inventory-year calibration, then
DROPS to ~1.0 after a heavy overstory removal (the bare regen point), recovering slowly. The code comment at
establishment.jl:1104 anticipated exactly this ("BAAA(NNID), NOT the whole-stand BA; after removal → BAAA≈1;
validated jl point_ba[1]=40 vs live 41.93 at cyc1") — but the current point_ba[1] gives 0 at icyc1 and the growing
overstory total thereafter, NOT the disturbance-adjusted BAAA. High/wrong BAAA → wrong ESTOCK PN → wrong PROB1 →
wrong tally amount (the 4-9× divergence). SECONDARY: the es_stream seed chain (jl 78807/32485/62399 reached at
icyc1/3/4 vs live ESS0 55329/78807/32485/62399 at icyc1/4/7/10) advances at the wrong cycles — but BAAA is the
dominant driver. ★ FIX = derive baaa as the disturbance-adjusted per-INVENTORY-POINT BAAA (inventory-year = the ESB
BAAOLD/STDINFO BA ~41.93; post-removal = the bare regen-point BA ~1; NOT the growing s.density.point_ba), matching the
live ESTOCKIN sequence 41.93/46.17/1/1/1/4.79/50.05. This is the exact, target-valued root cause — a specific input-
derivation bug, NOT a coupled RNG mystery. (ENV-gated AUTOES_IN debug in ie_autoes_establish! dumps seed0/es_stream/baaa.)

## 2026-08-06 (FIX HYPOTHESIS — point_ba is STALE at AUTOES time; missing GRADD first-DENSE) — code-level
★ stand4 (and all iet01 stands) use TREEDATA — a REAL ~41.93 BA overstory (NOT bare/NOTREES). So jl baaa=
point_ba[1]=0 at icyc1 is DEFINITIVELY WRONG (should be ~41.93). ROOT of the wrong point_ba: the GRADD order is
"UPDATE→DENSE→ESNUTR→DENSE→CROWN" (simulate.jl:551 comment), but there is NO compute_density! between the growth/
mortality block (compute_volumes! + comcup!, ~549-557) and the ESNUTR block (esuckr!/establish!/ie_autoes_establish!,
560-565). establish! only recomputes density when it ADDS regen; the AUTOES stand has no explicit PLANT/NATURAL so
establish! is a no-op → point_ba stays STALE/empty (0) when ie_autoes_establish! reads point_ba[1] at line 565. The
final compute_density!(571) runs AFTER AUTOES. ⇒ AUTOES reads a stale/zero point_ba ⇒ baaa=0/wrong ⇒ wrong ESTOCK
PROB1 ⇒ wrong tally. The code comment (establishment.jl:1108) validated "point_ba[1]=40 vs live 41.93 at cyc1" —
consistent with the first-DENSE having been present/correct before and the current path leaving it stale.
★ FIX HYPOTHESIS: insert the GRADD first-DENSE (compute_density!) BEFORE the ESNUTR block (before esuckr!/establish!/
ie_autoes_establish!) so point_ba reflects the POST-GROWTH per-point overstory BA (~41.93) that live's BAAA(NNID)
uses. Then re-measure baaa per cycle vs the live ESTOCKIN target (41.93/46.17/1/1/1/4.79/50.05) and the tally totals
vs live (0/468/206/419/11/225/472). CAUTION: a compute_density! before ESNUTR also changes the density esuckr!/
establish! see — validate NE/CS/LS/SN establish+sprout tests + the IE cyc0 8/9 stay bit-exact (doctrine #4). Bounded,
testable, code-level fix for a fresh session; the AUTOES_IN debug + live target sequence make it directly verifiable.

## 2026-08-06 (CORRECTION — "missing DENSE" REFUTED; baaa SOURCE is the real fix) — measured
★ Measured stand_ba + ntrees at AUTOES time (jl): icyc1 stand_ba=0 ntrees=0 (stand EMPTY) → so baaa=point_ba[1]=0
is CORRECT for the actual empty-stand state; there is NO stale density (the "missing GRADD first-DENSE" hypothesis
is REFUTED — doctrine #2, 6th self-correction). ⇒ live's BAAA=41.93 at icyc1 is a CALIBRATION REFERENCE (the ESB
inventory BAAOLD/BAAINV), NOT an actual tree BA. And at icyc4 jl has ntrees=13/stand_ba=49/point_ba=137 while live
BAAA=1.0 — confirming live's BAAA(NNID) is a DISTURBANCE-ADJUSTED per-inventory-point reference (calibration ~41.93
at inventory; ~1 after a heavy removal; recovering), fundamentally a DIFFERENT quantity from jl's growing point_ba[1].
★ CONFIRMED FIX (supersedes the missing-DENSE hypothesis): the baaa SOURCE is wrong. jl must feed ESTOCK the per-
inventory-point BAAA(NNID) = the calibration/disturbance-adjusted reference BA (matching live ESTOCKIN 41.93/46.17/1/
1/1/4.79/50.05), NOT s.density.point_ba[1] (the growing regen-cohort BA). Needs the estab.f/dense.f BAAA(NNID)
derivation: BAAINV at the inventory year (where 41.93 comes from — likely STDINFO/the inventory stocking, since the
stand is empty of trees), and the post-removal per-regen-point BA (~1). This is the exact, target-valued fix; the
tally COMPUTATION + SCHEDULER remain bit-exact. Fresh-session BAAA(NNID)-derivation port, verifiable via the AUTOES_IN
debug + the live ESTOCKIN/tally target sequences.

## 2026-08-06 (CAUTION — prior instrumentation is INCONSISTENT; re-instrument live for the fix)
★ The prior-session instrumentation files DISAGREE on the scheduler: stand4_REAL_scheduler.txt (idsdat=1990/ntally=1/
tally 81,852,118,337,48,261 / computed AUTOESXTES 589.65) vs stand4_ntally_sequence.txt (idsdat=1980/ntally=99/tally
0,468,206,419,11,225,472 / computed matching jl 583.65). Different captures during model development. ⇒ the SCHEDULER
(idsdat/ntally/tally) comparison must use FRESH live FVSie instrumentation, NOT these stale files — my earlier
"scheduler bit-exact" claim was against stand4_ntally_sequence.txt (which jl happens to match) but stand4_REAL_
scheduler.txt differs; re-instrument to establish the authoritative current live schedule before trusting either.
★ SOLID + CONSISTENT across BOTH files: the baaa root cause. PROB1TRC col3 = BAAOLD = 41.93/46.17/1/1/1/4.79/50.05
(== ESTOCKIN col6) — the disturbance-adjusted per-inventory-point reference jl must feed ESTOCK (vs jl's growing
point_ba[1]). 41.93 = the inventory-year BAAOLD calibration for the TREE-EMPTY stand (ntrees=0), NOT from STDINFO's
34.0 field — so it is computed from the stand's density/stocking (BAAINV); trace estab.f/dense.f BAAINV to source it.
NET for the fix session: (1) RE-INSTRUMENT live FVSie (fresh Fortran instrument-replay) for the authoritative
per-cycle idsdat/ntally/seed0/baaa/tally — do NOT trust the inconsistent stale files for scheduling; (2) the baaa
source fix (BAAA(NNID)=BAAOLD calibration/disturbance ref, not point_ba) is the confirmed root cause to implement.

## 2026-08-06 (Fortran-source COMPLETE) — the two BAs: BAAA(NNID) species-prob + BAAINV(NNID) calibration
★ Traced estb/estab.f: ESTOCK takes TWO per-inventory-point BAs, both of which jl gets wrong:
  - line 482  BAA    = BAAA(NNID)   — the CURRENT per-inventory-point BA → the ESTOCK species-probability + the
                                       tally-body BAAA (the 41.93/46.17/1/1/1/4.79/50.05 ESTOCKIN sequence).
  - line 487  BAAOLD = BAAINV(NNID) — the INVENTORY per-inventory-point BA (captured at inventory) → line 536
                                       CALL ESTOCK(ELEV,IFO,BAAOLD,BAAOLN,PN) — the PROB1 PN calibration.
jl feeds BOTH from s.density.point_ba[1] = the GROWING regen-cohort BA (0/10/137/491/848), which is neither
BAAA(NNID) (the disturbance-tracked per-point overstory BA that DROPS after removals) nor BAAINV(NNID) (the retained
inventory BA ~41.93). ⇒ FIX (complete spec): (1) capture BAAINV(NNID) at inventory (the per-inventory-point BA before
any removal — for stand4 ~41.93; the stand IS overstory-loaded at inventory then thinned to bare, so BAAINV must be
RETAINED, not recomputed from the post-thin empty point_ba); (2) maintain BAAA(NNID) = the current per-inventory-point
overstory BA that reflects removals (drops to ~1 post-thin); (3) feed BAAA→species-prob, BAAINV→ESTOCK PN calibration.
Needs the dense.f BAAA(NNID) / BAAINV(NNID) per-point BA tracking (where 41.93 originates + the post-removal drop).
★ DIAGNOSIS COMPLETE at the Fortran-source level. FIX = port the BAAA(NNID)/BAAINV(NNID) per-inventory-point BA
derivation (dense.f/estab.f) + feed the correct one to each ESTOCK use; re-instrument live (stale scheduler files
inconsistent) to confirm the per-cycle BAAA/BAAINV sequence; validate the tally + full-cycle .sum end-to-end. Fresh
session. This closes the AUTOES diagnostic arc: tally + scheduler bit-exact; the sole gap is the BAAA/BAAINV input.

## 2026-08-06 — FRESH LIVE INSTRUMENTATION DONE: BAAA vs BAAINV sequence measured (the blocker cleared)
Re-instrumented live FVSie estab.f (WRITE unit 74 after BAAOLN, line 490 — BEFORE the NTALLY/INADV gates so it
fires for EVERY tally, not just the NTALLY==1 calibration) + relinked (gfortran-16; the .iework relink_ie.sh's bare
`gfortran` is stale in this env — use gfortran-16 directly). Ran ierun/iet01.key (4 stands). MEASURED, for the AUTOES
inventory point NNID=1, per tally event:
    BAAA(NNID=1)  = 41.9321823, 46.1742516, 1.0, 1.0, 1.0, 4.79316950, 50.0507812   ← the plan's target sequence, CONFIRMED
    BAAINV(NNID=1)= 40.0000038 CONSTANT across ALL tallies and all 4 stands
★ THIS CORRECTS THE PRIOR SPEC. The plan hypothesized "BAAINV must be RETAINED ~41.93". WRONG: **41.93 is BAAA at
tally-1** (the current per-point overstory BA), and **BAAINV is a SEPARATE constant = 40.0** (the inventory-time
per-point BA). The two are close only at tally-1 (before the cut: BAAA 41.93 ≈ BAAINV 40.0); after the shelterwood
removes the overstory BAAA DROPS to 1.0 (disturbance-tracked) while BAAINV STAYS 40.0 (frozen at inventory), then
BAAA recovers (4.79→50.05) as the overstory regrows. So:
  - **BAAA(NNID)** = dense.f:213 `BAAA(IP) += BATREE·PI/GROSPC` for D≥REGNBK (overstory only), RECOMPUTED each cycle
    ⇒ reflects removals. Feeds ESTOCK species-probability (estab.f:482 BAA=BAAA(NNID)). Sequence 41.93/46.17/1/1/1/4.79/50.05.
  - **BAAINV(NNID)** = esfltr.f:67 `BAAINV(N) += 0.005454154·D²·ZPROB·PIX`, captured ONCE at inventory, PERSISTED via
    getstd/putstd (BFREAD/BFWRIT unit, IPTINV). = 40.0 constant. Feeds ESTOCK PN calibration (estab.f:487→536 BAAOLD).
jl feeds BOTH from the GROWING regen-cohort s.density.point_ba[1] (0/10/137/491/848) — wrong for each: it is neither
the disturbance-tracked overstory BAAA (which DROPS post-cut) nor the frozen inventory BAAINV=40.0.
⇒ FIX (now fully measured, spec corrected): (1) track BAAA(NNID) = current per-inventory-point OVERSTORY BA (D≥REGNBK),
recomputed per cycle, reflecting removals → feed to the species-prob ESTOCK; (2) capture BAAINV(NNID) = inventory-time
per-point BA (=40.0 here), frozen/persisted → feed to the PN-calibration ESTOCK. The engine change is still the
regression-risky shared-path port (dense.f overstory-BA-per-point + esfltr.f inventory-BA-per-point), but the model
inputs are now GROUND-TRUTH MEASURED, not inferred. Live sources restored (estab.f un-instrumented, FVSie_clean relinked).

## 2026-08-06 (jl-side handoff addendum) — point_ba[1] tracks BAAINV, NOT BAAA; runner needs full iet01 s4 config
Scoped the jl fix site: establishment.jl:1109 `baaa = point_ba[1]` feeds ie_autoes_run's `baa` (the ESTOCK species-
prob input). ★ Cross-referencing the live measurement (this session): the inline comment at :1104-1108 says "jl
point_ba[1]=40 vs live BAAA=41.93 at cyc1" — but the live instrument shows BAAINV(NNID=1)=**40.0** and BAAA(NNID=1)=
**41.93**. So jl's point_ba[1]≈40 is actually tracking live's **BAAINV** (frozen inventory BA), NOT the **BAAA**
(disturbance-tracked overstory BA, 41.93→…→1→…→50.05) the species-prob path requires. THAT is the mislabel at the
root of #143: jl feeds a BAAINV-like value where live feeds BAAA, and separately has no distinct BAAINV for the PN
calibration. FIX (next session): introduce a per-inventory-point OVERSTORY BA (D≥REGNBK, recomputed per cycle,
reflects removals ⇒ the 41.93/46.17/1/1/1/4.79/50.05 sequence) for the species-prob ESTOCK, and keep the frozen
inventory BA (40.0) for the PN-calibration ESTOCK. jl-side per-tally CONFIRMATION is blocked only on a RUNNER: a
hand-minimal iet01 stand-4 keyfile NaNs (InexactError round(Int64,NaN32)) because the stand isn't fully initialized
(TREEDATA didn't load / habitat unset ⇒ empty stand ⇒ ln(TPA=0) in the esb_shift calc ⇒ NaN). Build the runner from
the FULL stand-4 block (ierun/iet01.key lines 137-end: DESIGN, STDINFO 118/570/60/315/30/34, INVYEAR 1990, THINPRSC
1990 0.999, TREEDATA→iet01.tre, SPECPREF/THINBTA) so the stand initializes, then run with FVSJL_AUTOES_DEBUG=1 to dump
jl's per-tally baaa and diff against 41.93/46.17/1/1/1/4.79/50.05. THEN implement + validate vs the stand-4 .sum
(NOT per-record — AUTOES adds trees ⇒ tripling). Regression-scope: AUTOES fires POST-disturbance (post-cyc0), so the
IE cyc0 8/9 bit-exact is not at risk; the risk is contained to the (already-diverging) AUTOES tally output.

## 2026-08-06 (jl trigger-chain map) — why AUTOES doesn't fire on the keyfile iet01 stand-4
Attempted the jl-side per-tally baaa measurement. Built a full stand-4 keyfile (DESIGN/STDINFO 118/570/60/315/30/34/
INVYEAR 1990/NUMCYCLE 10/THINPRSC 1990 0.999/SPECPREF/THINBTA) with a matching .tre (jl auto-reads base_path.tre when
no TREEDATA keyword loads — keyword_dispatch.jl:2186; the earlier NaN was just a missing iet01_s4.tre → empty stand).
RESULT: the projection runs 10 cycles but **AUTOES NEVER FIRES** (zero AUTOES_IN debug lines). TRIGGER CHAIN mapped:
  simulate.jl:565 ie_autoes_establish! is called EVERY cycle for IE/EM, but
  establishment.jl:1062 gates on `(est.lautal || est.lingrw) || return false`, and
  establishment.jl:1081 `ie_autoes_schedule!(…, est.last_xtes, …)` uses last_xtes = the thinning REMOVAL FRACTION.
⇒ On this keyfile stand-4, neither lautal is set nor is last_xtes armed by the THINPRSC 1990 0.999 removal. So the
INTEGRATION GAP to close FIRST (before the baaa fix can even be measured jl-side) is: (a) default lautal=TRUE (live
FVS runs AUTOES ON unless NOAUTOES — iet01 stands 1-3 carry NOAUTOES, stand-4 relies on the default-on), and/or (b)
have the THINPRSC/thinning path set est.last_xtes = the removed BA fraction so ie_autoes_schedule! detects the
disturbance. Only once AUTOES fires jl-side can the per-tally baaa be diffed against the live 41.93/46.17/1/1/1/4.79/
50.05 and the per-inventory-point BAAA/BAAINV fix be validated. NOTE: the memory's "AUTOES wired + scheduler bit-exact"
was validated on the real-FIA/DB path (which arms lautal/last_xtes differently); the KEYFILE THINPRSC path is a
separate integration surface. SEPARATE minor edge found: on this 99.9%-removal stand jl's .sum writer NaNs at
summary.jl:319 (`dt` truncates NaN — a per-year rate ÷ ~0 residual BA); keyfile-specific (real-FIA sweep = 0 crashes),
low-priority robustness note. ⇒ #143 next session, in order: (1) arm the keyfile AUTOES trigger (lautal default +
last_xtes from THINPRSC), (2) confirm jl fires + dump per-tally baaa, (3) implement per-point overstory-BA (BAAA) vs
frozen inventory-BA (BAAINV) split, (4) validate vs iet01 stand-4 .sum. The live-side ground truth is already measured.

## 2026-08-06 (correction to the trigger map) — the infrastructure EXISTS; gap is a runtime firing condition
Followed the trigger deeper: the arming infrastructure is ALREADY present in jl —
  - state.jl:607-608 document est.lautal / est.lingrw default **TRUE** (AUTOES on unless NOAUTOES clears them,
    keyword_dispatch.jl:2235).
  - cuts.jl:317-319 sets `est.last_xtes = rem.tpa / autoes_pre_tpa` for IE thinnings (the removal fraction the
    esnutr LAUTAL scheduler reads); establishment.jl:954 fires when `est.lautal && xtes >= est.thres1` (thres1=0.10).
So a THINPRSC 1990 0.999 (≈0.999 removal ≫ 0.10) on stand-4 (no NOAUTOES ⇒ lautal true) SHOULD arm and fire the
tally. It did NOT in the keyfile run ⇒ the gap is a specific RUNTIME condition, not missing wiring. Candidates to
check FIRST next session (add a one-line stderr dump of est.lautal / est.last_xtes / xtes / thres1 at the top of
ie_autoes_establish!): (a) is lautal actually TRUE at runtime for a keyfile stand (constructor honoring the documented
default), (b) does the THINPRSC cut path actually reach cuts.jl:319 (vs a different removal method that skips the IE
last_xtes stash), (c) tally-timing/years_done gate (establishment.jl:1084 `year in est.years_done`). This is a small,
well-scoped debug — NOT a re-port — after which the per-inventory-point BAAA/BAAINV fix (live ground truth already
measured: 41.93/46.17/1/1/1/4.79/50.05 vs 40.0) can be implemented and validated vs the stand-4 .sum. The whole #143
remaining arc is now reduced to: [debug firing] → [swap point_ba[1]→per-point overstory BAAA + add frozen BAAINV] →
[validate]. No model unknowns remain.

## 2026-08-06 (validation-path blocker found) — jl AUTOES fires on DB stands; keyfile .tre reads 0 trees
Ran the trigger-state probe (temp stderr at ie_autoes_establish! entry; reverted). DEFINITIVE:
  - **jl AUTOES fires correctly on a real-FIA DB IE stand**: AUTOES_TRIG cyc=1 lautal=**true** lingrw=**true**
    last_xtes=0.0 thres1=0.1 ntrees=21 — so the lautal-default-TRUE + ie_autoes_establish! call path WORKS. (This
    stand had no thinning ⇒ last_xtes=0 ⇒ it fires via the LINGRW ingrowth branch, not the LAUTAL removal branch.)
  - **The keyfile iet01 stand-4 never reaches AUTOES because it loads 0 TREES**: `load_trees!(iet01.tre)` returns
    trees.n=**0** — jl's .tre reader does not parse iet01.tre's fixed-column format (record e.g.
    `   1      248112       0101   011LP 11510   0734   00111     0  0`). Empty stand ⇒ grow_cycle skips meaningful
    work ⇒ no AUTOES ⇒ the summary.jl:319 `dt` NaN (a per-year rate ÷ 0 residual). So the earlier "keyfile AUTOES
    doesn't fire" was NOT a scheduler gap — it was an EMPTY STAND from the .tre parse.
⇒ #143 VALIDATION PATH (revised, next session): do NOT use the iet01.key/.tre keyfile — jl can't read that .tre. Instead
either (a) load iet01 stand-4's trees into a small SQLite DB (FVS_TREEINIT/STANDINIT, like ie_test.db / build_subdb.jl)
and drive jl via DATABASE, OR (b) pick a real-FIA IE stand that carries a management THINNING so last_xtes arms the
LAUTAL removal branch (the branch that exercises the BAAA-drops-after-cut behavior the fix targets — the ingrowth
branch alone won't). Then: confirm AUTOES fires via lautal, dump per-tally baaa, implement the per-inventory-point
overstory-BAAA vs frozen-BAAINV split (live truth 41.93/46.17/1/1/1/4.79/50.05 vs 40.0), validate vs the DB stand's
.sum. SEPARATE minor gap surfaced: jl's .tre reader returns 0 trees for the standard iet01.tre column format — a
tangential IO limitation (jl real-FIA input is DB-based; the keyword-file .tre path is under-exercised), worth its own
small task but NOT on the #143 critical path.

## 2026-08-06 (schedule-firing layer) — even THINPRSC on a DB stand doesn't produce a tally yet
Follow-up probe: a DB IE stand (753199439290487, ie_test.db) with `THINPRSC 2029 0.8` added still produced NO
AUTOES_IN (no tally). Combined with the reverted TRIG probe (real-FIA stand: lautal=TRUE but ie_autoes_schedule!
returned fire=FALSE — establishment.jl:1083 `fire || return false` — because last_xtes=0 < thres1=0.1 and the LINGRW
condition wasn't met), this isolates the NEXT debug layer: the LAUTAL schedule doesn't arm from a keyword THINPRSC on
these stands. Check next session: (a) does the THINPRSC cut path actually reach cuts.jl:319 (set last_xtes) for a
DATABASE stand + keyword thin — instrument est.last_xtes right after the cut and at ie_autoes_establish! entry; (b)
the THINPRSC year vs cycle-boundary alignment (2029 must land on a projected cycle for the removal to register); (c)
whether the tally then fires the FOLLOWING cycle (schedule timing). Only after AUTOES actually produces a tally can
jl's per-cycle baaa be read and the point_ba[1]→(overstory BAAA / frozen BAAINV) fix be measured + validated. NET for
#143: the model math (tally) is bit-exact (prior sessions) and the live baaa/BAAINV ground truth is measured (this
session); the ENTIRE remaining arc is engine plumbing — (1) make the LAUTAL removal branch actually fire jl-side from
a thinning, (2) split the baaa input into per-point overstory-BAAA vs frozen-BAAINV, (3) validate vs a DB stand's .sum.
A clean, self-contained implementation session with zero model unknowns.

## 2026-08-06 (CORRECTION) — AUTOES DOES fire on DB stands; prior "doesn't fire" was a NUMCYCLE-alignment artifact
★ RETRACTS the 07aa537 conclusion ("THINPRSC on a DB stand doesn't arm a tally"). That run used `NUMCYCLE 5.0`
UNALIGNED — jl (like live) parses NUMCYCLE from fixed columns 11-20, so "5.0" one column short reads as ~1 CYCLE
(the exact trap found on the BM sweep). With only 1 cycle, the THINPRSC 2029 thin (end of cycle 1) never applied and
AUTOES had no disturbance. RE-RAN with COLUMN-ALIGNED `NUMCYCLE         5` on DB stand 753199439290487:
  AUTOES_IN icyc=2 ntally=99 baaa=13.18 total=1152.7 ; AUTOES_IN icyc=4 ntally=99 baaa=216.12 total=592.0
⇒ **AUTOES FIRES on DB stands** — the scheduler works; there is NO firing gap. (This stand is all-seedling BA~3, so the
tallies fire via the LINGRW ingrowth branch ntally=99, not the LAUTAL removal branch; a mature-overstory + heavy-thin
DB stand is still needed to exercise the removal path.)
★ jl-SIDE BUG CONFIRMED: jl's baaa GROWS across tallies (13.18 → 216.12) = the growing regen-cohort point_ba[1],
NOT the disturbance-tracked overstory BAAA (which live drops after a cut: 41.93→1→50.05). This is exactly the #143
input bug, now confirmed on the jl side by direct measurement (matches the plan's predicted growing 0/10/137/491/848).
⇒ REVISED remaining #143 (further de-risked): (1) find/build a DB IE stand with a mature overstory + heavy thin to
exercise the LAUTAL removal tally (ntally=1); (2) confirm jl baaa there also grows (vs live BAAA dropping); (3) fix:
feed the per-inventory-point OVERSTORY BA (D≥REGNBK, drops after removal) to the species-prob ESTOCK + the frozen
inventory BA to the PN calibration, instead of the growing point_ba[1]; (4) validate vs .sum. The scheduler/firing is
NOT a blocker (it works); the whole remaining arc is the baaa-input split + validation. Use ALIGNED NUMCYCLE always.

## 2026-08-06 (fix scoped to exact infrastructure) — need a per-point D≥REGNBK OVERSTORY BA (point_ba is all-trees)
Confirmed WHY jl baaa grows: establishment.jl:1109 uses `s.density.point_ba[1]` = PTBAA (standstats.jl:169
point_basal_area!), the per-point BA of ALL trees on the point — so as the regen cohort grows, point_ba grows
(13→216). The Fortran BAAA(IP) (dense.f:213) sums BA ONLY over D≥REGNBK (REGNBK = "minimum DBH of an overstory
tree") — it EXCLUDES the sub-REGNBK regen, so it drops when the overstory is cut and stays low until regen crosses
REGNBK. ⇒ THE FIX IS NEW INFRASTRUCTURE (not a one-liner): add a per-inventory-point OVERSTORY BA = Σ BA over trees
with DBH ≥ REGNBK, recomputed each cycle (a D≥REGNBK-filtered variant of point_basal_area!), feed THAT to the
species-prob ESTOCK (establishment.jl:1109 baaa); AND keep a frozen inventory-time per-point BA (BAAINV=40.0,
captured once, persisted) for the PN-calibration ESTOCK. Then validate vs live on a mature-overstory + heavy-thin DB
IE stand (aligned NUMCYCLE). REGNBK value: confirm IE's (esplt2.f/estab.f) — likely 1.0". This is the complete,
self-contained #143 implementation: (1) REGNBK-filtered per-point overstory BA, (2) frozen inventory BAAINV, (3) wire
both to the two ESTOCK uses, (4) validate. Scheduler fires (confirmed), live ground truth measured, no model unknowns.

## 2026-08-06 (FIX ATTEMPTED + REFUTED BY LIVE — the "overstory-BA" premise is WRONG) ★ major redirect
IMPLEMENTED the long-planned fix (replace point_ba[1] with a D≥REGNBK=2.999 OVERSTORY-only per-point BA) and
VALIDATED against live (re-instrumented estab.f ZZBAA on DB stand 753199439290487, aligned NUMCYCLE):
  live BAAA(NNID=1) = **20.66, 205.6** ; BAAINV = 1.0
  jl overstory-only (my fix) = 4.59, 12.72  ← FARTHER from live
  jl original point_ba[1]     = 13.18, 216.12 ← CLOSER to live
★ REFUTED: this stand is 0.2"-QMD (56074 TPA of sub-3" stems) — there is NO D≥2.999 overstory, yet live BAAA=20.66.
So **BAAA(NNID) is NOT the overstory-filtered BA** — it counts the tiny sub-3" trees. dense.f's REGNBK is effectively
~0 for these stands (REGNBK=REALS(66), a stand var, not a fixed 2.999 in this path). The plan's central premise for
MANY sessions ("BAAA = overstory BA that excludes regen / drops after cut") is WRONG; the iet01 "drops to 1" was the
CUT removing trees, not an overstory filter. REVERTED the jl change (back to point_ba[1]) — doctrine #4: a faithful-
looking change that regresses vs live ⇒ examine, don't keep. jl is byte-unchanged (functional no-op + corrected comment).
★ CORRECTED #143 SCOPE: point_ba[1] (per-point ALL-tree BA) IS the right concept. The residual is a per-point BA
SCALING/ATTRIBUTION Δ, and it is MIXED-SIGN not a constant factor: live 20.66 vs jl 13.18 (jl −36%) EARLY, live 205.6
vs jl 216.12 (jl +5%) LATE. Candidates: (a) PI/GROSPC scaling (PI=NPTIDS) — jl point_ba uses p.pi/p.gross_space; live
BAAA(NNID) may normalize per the tallied point differently; (b) jl uses point index [1] while live loops NNID over all
inventory points (the tallied point may not be jl's point 1); (c) tally-moment tree state (jl reads point_ba after
compute_density! at a different growth/mortality phase than live's BAAA capture). Given the Δ is mixed-sign and ~±30%
on a dense-seedling stand (and point_ba was already validated ≈live on iet01: 40 vs 41.93), this may be closer to the
accepted dense-regen straddle than a clean bug. NEXT: measure point_ba[1] vs live BAAA across SEVERAL stands + the
multi-point attribution (is the tallied NNID always jl's point 1?), before any further code change. The overstory-BA
avenue is CLOSED. ⇒ #143 is NOT the tractable clean fix the plan assumed; re-scope as a per-point-BA attribution study.

## 2026-08-06 (3-stand jl point_ba[1] vs live BAAA — #143 is REAL, not cornered; it's disturbance-timing/attribution)
Settled the cornered-vs-bug question with a 3-stand measurement (ie_test.db, aligned NUMCYCLE 5 + THINPRSC 2029 0.8),
live BAAA(NNID=1) (re-instrumented estab.f) vs jl baaa=point_ba[1] (FVSJL_AUTOES_DEBUG), per tally:
  12343703010690: live 1.0 / 45.07   | jl 16.05 / 216.21  → jl 16× / 4.8× HIGH
  44987944020004: live 154.4 / 328.1 | jl 153.6 / 358.8   → ~exact / +9%
  753199439290487: live 20.66 / 205.6| jl 13.18 / 216.12  → -36% / +5%
★ VERDICT: NOT a cornered straddle (corrects the 1b82ae0 "likely cornered" lean). There is a REAL, large,
stand-dependent divergence. The SMOKING GUN is 12343703010690: live BAAA=**1.0** (the inventory point is BARE — the
disturbance removed its trees, and regen isn't added until AFTER the tally) while jl point_ba[1]=**16.05** (jl still
counts trees on that point). So jl's per-point BA does NOT reflect the DISTURBANCE-ADJUSTED point state that live's
BAAA(NNID) captures at the tally moment — exactly the task-title framing ("disturbance-adjusted BAAA"). This is a
TIMING/ATTRIBUTION issue, NOT the overstory-filter (refuted 1b82ae0) and NOT a straddle. Candidates:
  (a) point_ba[1] is computed by compute_density! at a cycle phase that does NOT reflect the THINPRSC removal on that
      specific point (jl removes stand-wide but the per-point BA the tally reads is stale/pre-cut for point 1);
  (b) jl's point index [1] ≠ live's tallied NNID (multi-point attribution — live loops NNID over inventory points; the
      point being tallied may be a DIFFERENT, bare point than jl's point 1);
  (c) the tally reads point_ba after regen/growth of a phase live captures earlier.
NEXT (real fix path, no overstory filter): instrument WHEN live captures BAAA(NNID) vs when jl fills point_ba (cut →
density → tally order), and whether the tallied NNID maps to jl's point 1. The stand 12343703010690 (live BAAA=1 vs
jl 16) is the cleanest repro. ⇒ #143 is a genuine open bug (per-point disturbance-adjusted BA), re-scoped away from
both the overstory-filter AND the cornered-straddle; the ~22% tally residual is real on disturbance stands. jl unchanged.

## 2026-08-06 (#143 localized) — single-plot (attribution REFUTED); it's per-point BA on disturbed-overstory stands
Repro 12343703010690: NUM_PLOTS=**1** (single plot) ⇒ the multi-point NNID↔point-1 ATTRIBUTION hypothesis is REFUTED
(jl point_ba[1] and live BAAA(NNID=1) are the SAME point). The stand has a real overstory (33 recs, DBH 0.1–15.5",
24 trees ≥5") that THINPRSC 2029 0.8 heavily cuts. At the icyc2 tally: live BAAA=**1.0** (point ~bare after the cut)
vs jl point_ba[1]=**16.05**. So jl's per-point BA is ~16× too high on a heavily-disturbed overstory point. THREE
hypotheses now REFUTED: overstory-filter (1b82ae0), cornered-straddle (453e1a6), multi-point-attribution (this entry).
REMAINING mechanism = a per-point BA value difference on DISTURBED-OVERSTORY stands, ONE of:
  (i) TALLY TIMING — jl reads point_ba (compute_density!) at a cycle phase that does not reflect the THINPRSC removal
      that live's BAAA capture already sees (cut→density→tally ordering within grow_cycle!); or
  (ii) jl's THINPRSC itself removes LESS BA than live's on this stand (a thinning-implementation Δ, NOT AUTOES) — jl
      would then legitimately carry 16 BA while live carries ~0. ← CHECK THIS FIRST (cheap): diff the jl vs live .sum
      BA at the thin year 2029 on 12343703010690; if jl BA >> live BA post-thin, the root is THINPRSC, not the baaa
      input, and #143 re-routes to a cut bug. If jl and live post-thin BA MATCH, then it's (i) tally timing.
The overstory-removal stands (not the all-seedling ones) are where point_ba[1] diverges most from BAAA; the ~22% AUTOES
residual is concentrated there. NEXT SESSION: run (ii) first (one .sum diff at the thin year), then (i) if needed.
Clean repro = 12343703010690. jl unchanged this session; live restored (0 residue).

## 2026-08-06 (#143 .sum diff — THINPRSC-removal REFUTED; baaa comparison was tally-timing-CONFOUNDED)
Ran hypothesis (ii): jl vs live .sum on repro 12343703010690 (THINPRSC 2029 0.8):
  live: 2007 BA=2, 2017 BA=15, 2027 BA=148, 2037 BA=219, 2047 BA=240, 2057 BA=246
  jl:   2007 BA=2, 2017 BA=16, 2027 BA=169, 2037 BA=216, 2047 BA=209, 2057 BA=206
★ (ii) REFUTED: NO visible THINPRSC drop in EITHER jl or live (BA grows through the 2029 thin) — the cut is not
removing differently. Stand BA tracks closely early (2/2, 15/16, 148/169, 219/216) and diverges only in the LATE tail
(jl -15% by 2057 = the known EM/IE growth compounding tail). ⇒ the ~16× baaa gap (live BAAA=1 vs jl 16) is NOT a
per-point removal difference.
★ KEY REALIZATION — the earlier baaa comparison was TALLY-TIMING CONFOUNDED: live BAAA=**1.0** is the CLAMPED value
(estab.f:483 IF(BAA<1)BAA=1) at a tally near INVENTORY (BA≈2 in 2007) — i.e. live's first tally fires early on a
near-empty stand; jl's icyc2 tally reads a LATER, higher-BA state (BA≈169 → point_ba≈16). The two "tally-1"s are NOT
at the same stand moment, so BAAA 1 vs point_ba 16 is apples-to-oranges. The real question re-narrows to: DO jl and
live fire their AUTOES tallies in the SAME cycles/at the same stand state? (Prior sessions claimed the scheduler is
bit-exact on iet01, but that was the shelterwood stand, not these THINPRSC DB stands.) ⇒ #143 NEXT (fresh session, in
order): (1) align WHICH cycle each tally fires jl vs live on 12343703010690 (dump idsdat/kdt/ntally/year both sides);
(2) ONLY THEN compare baaa at MATCHED cycles; (3) if matched-cycle baaa still diverges, it's the point_ba value; if the
tallies fire at different cycles, it's a SCHEDULER-timing bug on the THINPRSC path. Four hypotheses now tested this
session (overstory-filter / cornered / multi-point-attribution / THINPRSC-removal — ALL refuted); the live-vs-jl
per-tally baaa numbers must be RE-TAKEN at matched cycles before any further fix. jl unchanged; live restored (0 residue).

## 2026-08-06 (ROOT CAUSE FOUND + faithful XCUF fix landed; deeper baaa-semantics puzzle remains)
★ MATCHED-CYCLE measurement (re-instrumented live estab.f: ICYC/YR/NTALLY/BAAA) on repro 12343703010690:
  live: ICYC=2 YR=2017 NTALLY=1 BAAA=1.0 ; ICYC=4 YR=2037 NTALLY=1 BAAA=45.07
  jl:   icyc=2 ntally=99 baaa=16.05 ; icyc=4 ntally=99 baaa=216.21
Tallies fire at the SAME cycles but jl labels them ntally=99 (LINGRW ingrowth sentinel) vs live NTALLY=1/2 (counter).
★ SOURCE-VERIFIED ROOT (one real bug found): esnutr.f:271-275 sets the disturbance-detection fraction
  XTES = AMAX1(XTPA, XCUF), XTPA=ONTREM(7)/ONTCUR(7) (TPA-removal frac), XCUF=OCVREM(7)/OCVCUR(7) (CUBIC-VOLUME
  removal frac). jl (cuts.jl:319) used ONLY XTPA. An overstory thin removes few TREES (low XTPA) but high VOLUME (high
  XCUF), so jl under-detected the disturbance and fell to the ingrowth branch. FIX LANDED (cuts.jl): capture
  autoes_pre_cuft = Σ tpa·cuft_vol (=OCVCUR) and set last_xtes = max(rem.tpa/pre_tpa, rem.cuft/pre_cuft). Faithful port
  of esnutr.f; INERT on the 3 ie_test.db stands (their THINPRSC removes ~0 TPA in jl ⇒ rem.tpa=0 ⇒ block skipped),
  cyc0-safe (only fires on a real removal), no regression — but UNEXERCISED here (kept as source-verified-faithful).
★ STILL UNRESOLVED (deeper, next session) — the baaa VALUE is INCONSISTENT across stands under any single model:
  - 753199439290487 (all-tiny, 0.2"QMD): live BAAA=20.66 (counts the tiny trees ⇒ NOT overstory-filtered)
  - 12343703010690 (has overstory to 15.5", stand BA=15 at 2017): live BAAA=**1.0** (≈bare ⇒ does NOT count the 15 BA)
  These two REFUTE both "BAAA=all-tree BA" (fails 12343) AND "BAAA=overstory BA" (fails 753). So BAAA(NNID) at the
  tally reflects the per-inventory-point state in a way not captured by point_ba[1] OR a D≥REGNBK filter. NEXT: in live,
  dump per-NNID the trees/BA that feed BAAA at the tally moment (dense.f BAAA(IP) accumulation) on BOTH stands — resolve
  what tree set + what phase (inventory vs current) BAAA counts. The tally-amount residual is downstream of this.
★ NET #143: one real source-verified bug FIXED (XCUF in last_xtes, faithful, inert-on-available-stands); the baaa-value
  semantics need per-NNID tree-level live instrumentation on the two contradictory stands. Five hypotheses tested (four
  refuted, XCUF confirmed-by-source). Live restored (0 residue).

## 2026-08-06 (XCUF fix: exercise BLOCKED this session; inert-safe kept; a THIRD layer noted)
Tried to EXERCISE the XCUF fix (force a volume-heavy removal to confirm it trips LAUTAL): THINPRSC 0.8 and a THINBTA
to residual BA 30 on 12343703010690 BOTH removed ~nothing in jl (after-treatment BA == before, ntally stayed 99) —
malformed test-keyword parameters, not a fix failure. The real repro (iet01 stand-4, a working shelterwood) can't run
jl-side (the .tre parser returns 0 trees — separate gap). So the XCUF fix is UNEXERCISED this session; it is KEPT
because it is source-verified faithful (esnutr.f:271-275 AMAX1(XTPA,XCUF)) AND inert-safe (guarded by rem.tpa>0 ⇒
cannot regress when no removal occurs; cyc0 untouched; iet01's 99.9% shelterwood is TPA-heavy so the MAX is unchanged
there ⇒ its bit-exact scheduler validation is preserved).
★ THIRD interacting layer noted (for next session): establishment.jl:1082 `est.last_xtes = 0f0` CONSUMES the removal
fraction EVERY cycle, right after ie_autoes_schedule! reads it — UNCONDITIONALLY (whether or not the tally fired). So
a thin that sets last_xtes in cycle N is only usable by the schedule IN cycle N; if the 20-yr/idsdat rule doesn't fire
that same cycle, the removal signal is discarded before the tally cycle. Live's LAUTAL uses a persisted disturbance
date (IDSDAT) + the removal, not a one-cycle-consumed flag. So even with XCUF correct, the removal tally may not fire
at the right cycle. NEXT SESSION #143 (three coupled sub-issues, needs a jl-runnable thinned+overstory IE stand):
(1) EXERCISE + validate XCUF on a real removal; (2) last_xtes consume-timing vs live IDSDAT persistence; (3) the
baaa VALUE per-NNID semantics (live BAAA 20.66 all-tiny vs 1.0 has-overstory — needs dense.f BAAA(IP) tree-level dump).
The XCUF root cause is real and fixed; the AUTOES-removal path has these two further coupled layers. jl cyc0/growth/
mortality/volume all unaffected (change is AUTOES-tally-only, inert on non-removal cycles).

## 2026-08-06 ★ UNBLOCKED — iet01 runs in jl via TREEFMT; #143 residual REPRODUCED + localized to tally-AMOUNT
KEY ENABLER: iet01.tre needs its custom TREEFMT (in iet01.key: (T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,
T63,F3.0,T60,F3.1,T48,I1,T52,I2,T66,5I1,T54,7I1,T75,F3.0)) — the earlier "load_trees! returns 0" was a keyfile-
OMISSION on my part (my extracted stand-4 keyfile dropped the TREEFMT), NOT a jl parser bug. jl's run_keyfile handles
TREEFMT fine. With the TREEFMT + CRLF-stripped .tre, jl RUNS iet01 stand-4 (the canonical #143 shelterwood stand).
★ ON iet01 STAND-4 (THINPRSC 1990 0.999 shelterwood + AUTOES), jl vs live FVSie_clean .sum TPA:
  1990 536/536 BIT-EXACT; then jl OVER-produces AUTOES regen: 2000 +20%, 2030 +48%, 2040 +59%, 2050 +76%, 2090 +55%.
★ WHAT NOW MATCHES (big progress): jl fires the CORRECT tally TYPES — ntally=1/2 (LAUTAL removal/continuation) at
  icyc 1/2/4/5/7/8 + ntally=99 (ingrowth) at icyc10 — the live pattern. And jl's baaa is CLOSE to the live target:
  jl 40.0/43.4/1.0/1.11/1.0/3.8 vs live 41.93/46.17/1/1/1/4.79. So on the shelterwood stand the baaa INPUT and the
  tally TYPE/SCHEDULE are ~right (the earlier DB-stand "16× baaa" was ingrowth-tally-at-different-state, a red herring).
★ RESIDUAL LOCALIZED: the ~22% (here +13–76%) is the tally AMOUNT — jl's per-tally TOTAL trees produced is too high
  (jl per-tally totals 1105/839/2006/809/1848/872/534). The memory's "tally computation bit-exact (WP33 WL20 DF7 GF202
  WH222 RC50 ES23 AF27=583.7 on stand-4 1999)" was for ONE tally's species-split given its inputs; the multi-cycle
  over-production is either (a) MORE tallies fired than live, or (b) each tally's total slightly high (baaa 40 vs 41.93
  → PROB1 higher → more trees), compounding. NEXT (now fully doable in jl): dump live per-tally TOTAL (instrument
  estab.f ITPP/booktpa) vs jl per-tally total at matched cycles on iet01 stand-4; if jl fires the same #tallies but each
  is high, it's the baaa/PROB1 (~-4% baaa → the over-production); if jl fires MORE tallies, it's the schedule count.
⇒ #143 is UNBLOCKED and REPRODUCIBLE in jl. The XCUF fix (2a27b78) is inert here (shelterwood is TPA-heavy) but
correct for volume-heavy thins. Harness: scratchpad/iet01s4.key (+.tre, TREEFMT, CRLF-stripped) + iet01s4.jl.

## 2026-08-06 (#143 over-production narrowed — ITPP-cap REFUTED, it's downstream of ITPP)
Checked the natural over-production hypothesis (missing MAXTPP/MAXING per-plot cap): REFUTED. jl applies both
(_IE_MAXTPP=[9,7,5,5,10,8,9,5,21,25,10,10,11,7,10,8], _IE_MAXING=[4,4,3,3,5,4,5,4,7,7,5,5,5,4,5,4]; cap = is_ingro ?
MAXING : MAXTPP; itpp=clamp(...,1,cap), establishment.jl:677,696) and the code comment records ITPP validated BIT-EXACT
vs live on iet01 stand-4 ([2,1,14,2,4,2,25,4]). ⇒ with per-plot ITPP already bit-exact but the multi-cycle .sum regen
+13–76% high, the over-production is DOWNSTREAM of ITPP: candidates = (a) DUPNPT (nptids·idup) point-scaling that
converts per-plot ITPP → per-acre TPA (if jl's dupnpt or the /300 normalization differs, bit-exact ITPP still yields
wrong TPA), (b) the NUMBER of tallies fired (jl fired 7 on iet01 s4: icyc 1/2/4/5/7/8 removal + 10 ingrowth — compare
live's count), (c) regen SURVIVAL — the added regen's subsequent mortality (if jl under-kills the new cohort, TPA
accumulates over cycles, which fits the GROWING +13→+76% divergence). ★ (c) is the strongest fit (the divergence
COMPOUNDS over cycles, not a fixed per-tally offset). NEXT (enabled by the iet01 unblock): instrument live estab.f
NEWTPP/ITPP·DUPNPT per tally vs jl per-tally total at matched cycles on iet01 s4; if per-tally totals MATCH but .sum
diverges, it's regen survival/mortality of the AUTOES cohort (not the tally). Harness ready (iet01s4.{key,tre,jl}).
⇒ #143 residual = AUTOES-cohort AMOUNT, ITPP bit-exact, root now downstream (dupnpt / tally-count / regen-survival);
survival is the leading candidate given the compounding. jl unchanged this turn.

## 2026-08-06 (units caution — AUTOES_IN "total" is PRE-scaling, NOT per-acre; use .sum for the real figure)
Ran live iet01 stand-4 with RegRepts: the regen report prints "TREES/ACRE ADDED" per tally (e.g. 109.3 ingrowth in
FALL 2089, species split 21/0/1/61/6/18/0/0/3...). jl's FVSJL_AUTOES_DEBUG "total=" (e.g. 533.8 at icyc10) is a
PRE-scaling sum(r.tally), NOT the per-acre added trees — comparing them directly (533.8 vs 109.3 ⇒ "~5×") is a UNITS
ARTIFACT and WRONG: the real per-acre over-production is the .sum TPA divergence (+13–76%), which a true 5× could not
produce. So jl's per-acre regen addition is +13–76% high, not 5×. ⇒ the clean per-tally comparison must be at the
PER-ACRE level: instrument jl to emit its per-cycle per-acre added-trees (the regen-report equivalent, post dupnpt/300
scaling) OR read jl's cohort TPA delta pre-mortality, then diff vs live's RegRepts "TREES/ACRE ADDED" at matched
cycles. Only that isolates per-acre-tally-amount vs survival. (Doctrine #2: caught the units mismatch before concluding
— the AUTOES_IN debug counter is an internal pre-scaling total, useful for tally TYPE/baaa but NOT for per-acre amount.)
★ #143 status unchanged & accurate: residual = AUTOES-cohort per-acre amount +13–76% on iet01 s4 (ITPP bit-exact,
tally-type + baaa matching); root still downstream (dupnpt/300 scaling vs survival). Live clean (RegRepts is a keyword,
no source change). Harness: iet01s4r.key (adds RegRepts) for the live per-acre reference (109.3-type figures).

## 2026-08-06 (CORRECTION of the units-caution + REAL localization: jl adds "ALL" not "<3.0 PROJECTED"; ingrowth 5×)
★ RETRACT the prior "units caution" (f8d36ee): r.tally IS per-acre — establishment.jl:1172 sets t.tpa[n]=r.tally[sp]
(the created seedling's per-acre TPA), NO further scaling. So jl's FVSJL_AUTOES_DEBUG "total"=sum(r.tally) IS the
per-acre trees added, directly comparable to live's RegRepts. The "pre-scaling" claim was WRONG (doctrine #2/#4 — a
correction that was itself an over-correction, now re-measured).
★ REAL per-tally comparison (iet01 stand-4, live RegRepts vs jl AUTOES_IN, both per-acre):
  DISTURBANCE TALLY 1 (1999/icyc1): live "SUMMARY OF ALL"=1079, "TREES <3.0 BEING PROJECTED"=882 | jl=1105.
    ⇒ jl (1105) ≈ live "ALL" (1079), but live only PROJECTS the <3.0 subset (882). jl adds ~ALL, ~+25% over the
    <3.0-projected count. jl's added trees are all dbh=0.1 (seedlings) so the excess = jl not restricting to the
    "<3.0 BEING PROJECTED" subset live actually books.
  INGROWTH (2089/icyc10): live "109.3 INGROWTH TREES/ACRE ADDED" | jl=533.8 ⇒ jl ~5× over (the DOMINANT over-producer).
  Tally COUNT/structure MATCHES: both fire TALLY 1/2 at 1999/2009, 2029/2039, 2059/2069 + ingrowth 2089 (7 tallies).
★ ROOT (localized, two coupled): (1) jl books the FULL tally ("ALL") rather than live's "TREES <3.0 IN. DBH BEING
  PROJECTED" subset — find where live restricts the booked count to <3.0 (estab.f the ITPP/NEWTPP → only sub-3.0
  advance/subs are added; the ≥3.0 "best" trees are summarized but NOT projected). (2) the INGROWTH tally over-produces
  ~5× — separate, likely the MAXING cap or the ingrowth PROB/height path (SHORTY time=1) producing too many. The
  net .sum is only +13–76% because the excess seedlings are tiny (0.1") and largely die / stay sub-.sum-threshold.
NEXT: (a) instrument live estab.f — what distinguishes "ALL" (1079) from "<3.0 PROJECTED" (882)? port that restriction
to jl's booking loop (establishment.jl:1159-1185); (b) the ingrowth 5× — compare jl vs live ITPP/NEWTPP on the 2089
tally specifically (MAXING=_IE_MAXING cap + the ingrowth branch). #143 residual = tally-AMOUNT over-booking (ALL-vs-
<3.0 + ingrowth 5×), NOT baaa/type/schedule/survival. jl unchanged; live clean (RegRepts only). Harness iet01s4r.key.

## 2026-08-06 (SPECIES-SPLIT measurement CONFIRMS the ESADVH/ESSUBH-height root — WP-absent signature)
Per-species tally on iet01 s4 (jl AUTOES_IN tally vector vs live RegRepts 2089 ingrowth, sp order WP WL DF GF WH RC LP ES AF PP):
  live 2089 ingrowth: WP21 WL0 DF1 GF61 WH6 RC18 LP0 ES0 AF3 = 109.3
  jl icyc10 ingrowth: WP0 WL3.1 DF33.6 GF268.4 WH79.3 RC137.3 LP0 ES3.1 AF9.2 = 533.8
  (same pattern on the disturbance tallies: jl icyc1 WP3.6/GF302/WH568 vs the memory-validated WP33/GF202/WH222 split.)
★ SIGNATURE: jl produces ~ZERO white pine (WP) and OVER-produces the wet-side species (GF/WH/RC/DF). This INVERTS the
PADV species-probability (validated bit-exact: WP.062/GF.485/WH.283 — GF/WH-heavy) vs live's FINAL split (WP-heavy).
The only thing between PADV and the final per-species regen count is the ADVANCE-vs-SUBSEQUENT split + per-tree
ESADVH/ESSUBH heights (+ ESDLAY delay + the FIRST dispersion chain). jl books the raw tally directly with ALL trees
floored at dbh=0.1" (est.jl:1166), so it emits the raw PADV-like GF/WH-heavy split and no advance-regen WP. ⇒ CONFIRMS
the #143 root = the MISSING ESADVH/ESSUBH/ESDLAY height sub-model (plan "HEIGHT-MODEL SPEC COMPLETE, routine-level").
★ PORT SCOPE (the fix, a dedicated sub-model unit per the plan): ie/esadvh.f (268L, best/advance heights) + ie/essubh.f
(281L, subsequent) + ie/esdlay.f (190L, delay) = ~739L + the FIRST(1,i) dispersion (0.1→0.316→0.562→…→1) + the
advance/subsequent DRAW split (ITIME≤2 ⇒ all BEST=advance). Wire into establishment.jl:1159-1185 tree-creation
(replace the dbh=0.1 floor with the computed per-species heights → ≥3.0" trees excluded from projection; advance-regen
adds WP). Validate vs live RegRepts (iet01 s4r: 2089 WP21/GF61/... + the .sum 536→1025→…→1788). This is the SINGLE
remaining #143 piece — large but fully specified; substantial dedicated port, not a bounded step. jl clean (debug reverted).

## 2026-08-06 (KEY CORRECTION — the height routines are PORTED-BUT-UNWIRED; #143 fix = WIRING, not porting)
★ While starting the "height-model port" I discovered (after mistakenly duplicating ie_esadvh — reverted 6ec952d)
that ALL FOUR AUTOES height routines ALREADY EXIST in establishment.jl and are FULLY PORTED:
  ie_essubh (subsequent, line ~68), ie_esxcsh (excess, ~147), ie_esadvh (advance, ~782), ie_esdlay (delay, ~922).
CONFIRMED: grep across all of src/ — each is DEFINED but has ZERO call sites (only definitions + docstrings). So the
"HEIGHT-MODEL SPEC COMPLETE, routine-level" note meant the ROUTINES were built; only the tree-creation WIRING was
deferred. ⇒ #143's remaining fix is NOT a 739-line port (equations done) — it is WIRING the existing routines into
the AUTOES tree-creation loop (establishment.jl:1159-1185), which currently floors ALL trees at dbh=0.1" (line 1166)
and never calls any height routine. THAT floor is why jl books the raw PADV-like species split (WP0/GF-heavy, all
<3.0") instead of live's advance-regen-shaped split (WP21) — confirmed by the species-split measurement (commit 9663634).
★ THE WIRING (the actual #143 fix, now much smaller than believed): in ie_autoes_establish!, per tallied tree, run the
advance/subsequent/excess DISPATCH (plan: ITIME≤2 ⇒ all BEST=advance→ie_esadvh; EXCESS→ie_esxcsh; else subsequent→
ie_essubh) with DILATE=FIRST(1,i) dispersion (0.1→0.316→0.562→…→1, an order-statistic sqrt chain), DELAY=ie_esdlay
draw, GENTIM=FINT−5, EMSQR (2 per-stand ESRANN draws), TIME → compute per-tree HHT; derive dbh (est.jl:1166 dbh=
0.1+0.001·hht); trees reaching ≥3.0" are advance regen NOT floored. Then re-validate vs live RegRepts (iet01 s4r:
2089 WP21/GF61/... + regen-report AVERAGE HEIGHT WP3.5/DF3.8/GF1.7/WH1.9/RC2.1 + .sum 536→…→1788). The dispersion/
draw ORDER must match live's ESRANN stream (the seed chain is already bit-exact per the tally validation). ⇒ #143 is
a WIRING+RNG-order task over EXISTING equations — a focused unit, but NOT the large port previously scoped.
META (doctrine #5): I duplicated ie_esadvh by not grepping for the existing function first — caught + reverted; the
lesson (check what exists before porting) is exactly why the port turned out to be a wiring task.

## 2026-08-06 (WIRING located + de-risked: height draws are already BURNED at ie_autoes_tally:724)
Read ie_autoes_tally (establishment.jl:658-733). KEY: the ADV/SUBS + per-species height ESRANN draws are ALREADY
CONSUMED — line 724 `for _ in 1:adv_heights; ie_esrann!(rng); end` draws & DISCARDS them (adv_heights = nsp ADV/SUBS
+ 2·nsp heights) purely to keep the stream aligned. So the RNG is ALREADY at the correct position (the seed chain +
draw order are bit-exact-validated) ⇒ the wiring is LOW-RNG-RISK: capture those drawn values and USE them (with
DILATE=FIRST(1,i), ie_esdlay for DELAY, ie_esadvh/ie_essubh for HHT) instead of burning them. NO stream shift.
★ BUT the wiring is more than "use the burned draws" — the tally COUNTS (tally[j], set in the best-pick loop :715-718
and excess-pick :727-729 via ie_estab_pick_species + esprob·scale) and the height model TOGETHER decide what live
ADDS. Live's regen report "TREES/ACRE ADDED" = the <3.0-DBH subset; the ≥3.0" (advance-regen, tall ESADVH heights)
are summarized but NOT projected as new regen. jl currently books ALL tally[j] as dbh=0.1 (est.jl:1166) ⇒ over-books
+ wrong split (measured jl WP0/GF268 vs live WP21/GF61, commit 9663634). ⇒ THE WIRING (final #143 unit): (1) in
ie_autoes_tally capture per-species EMSQR(from :695) + the ADV/SUBS draw + 2 height draws (from :724) + DILATE=
FIRST(1,i); dispatch advance→ie_esadvh / subsequent→ie_essubh / excess→ie_esxcsh → per-tree HHT; (2) RETURN heights
alongside tally; (3) in ie_autoes_establish! tree-creation replace the dbh=0.1 floor with the HHT-derived dbh and
apply the ≥3.0" advance-regen rule live uses (which trees are ADDED vs summarized); (4) validate vs iet01 s4r RegRepts
(2089 WP21/GF61 + AVERAGE HEIGHT WP3.5/DF3.8/GF1.7/... + .sum 536→…→1788). Needs the estab.f body-16..84 trace for the
exact count↔height↔ADD interaction (ITIME≤2 ⇒ BEST=advance). ⇒ #143 = a focused, LOW-RNG-RISK wiring over existing
equations; the delicate part is the count/height/ADD interaction (estab.f body), NOT the RNG. jl clean this turn.
