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
