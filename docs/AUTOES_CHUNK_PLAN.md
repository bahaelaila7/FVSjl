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
- **A3 — scheduler (esnutr.f rules):** the 20-yr-disturbance + ingrowth triggers → fire the tally in
  engine/establishment.jl's cycle hook. Reuse the existing tree-creation tail (naturals-first).
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
