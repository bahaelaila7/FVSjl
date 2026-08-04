# AUTOES — Western Automatic (Natural) Establishment tally — CHUNK PLAN (task #143)

Status: **SCOPED (not started).** Structure mapped + validation anchor measured 2026-08-04. No jl code yet.
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
- **A1 — ESTOCK:** port estock.f verbatim (5 habitat-series eqs + DATA coeffs) + the IHAB/IPREP/aspect/BAA inputs.
  Validate PN vs a live ESTOCK instrument-replay (dump PN at estab.f:536/572 for iet01 stand-4).
- **A2 — ESNSPE + heights:** species apportionment + ESADVH/ESSUBH (reuse the EM essubh generalization). Validate
  the predicted per-species TPP + heights vs instrument.
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
