# IE (Inland Empire) variant port — audit

## AUTOES disturbance-path amount [#143] — MEASURED entry point (2026-08-06)

The IE AUTOES tally-amount residual (iet01 stand4, live establishment ~590) is the DISTURBANCE path
(NTALLY=1, INGRO=0) — a DIFFERENT mechanism from the EM ingrowth path (NTALLY=99) that was driven bit-exact
(EM_VARIANT_PORT_AUDIT §AUTOES amount, 29d449d).

**Live IE disturbance tally (instrument-replay, FVSie_buildDir/estab.f:950 ESPROB dump, iet01):**
```
IEPRB nc=1 PROB1=0.60122 NEWTPP=2 ITPP=2  DUPNPT=50 ESPROB1=0.60122 INGRO=0 MATCH=0
IEPRB nc=2 PROB1=0.60122 NEWTPP=1 ITPP=1  ESPROB1=0.60122
IEPRB nc=3 PROB1=0.60122 NEWTPP=14 ITPP=14 ESPROB1=0.60122   ← ITPP=14 > MAXING(7): live uses MAXTPP(25)
IEPRB nc=4 PROB1=0.60122 NEWTPP=2 ITPP=2  ESPROB1=0.60122
IEPRB nc=5 PROB1=0.60122 NEWTPP=4 ITPP=4  ESPROB1=0.60122
(2nd tally) nc=1 PROB1=0.81726 NEWTPP=1 ITPP=3 ESPROB1=0.21604  ← continuation: ESPROB=FTEMP−PNN, books increment
```

**Confirmed mechanics** (now fully understood from the EM measurement):
- Disturbance ESPROB = FTEMP = PROB1 for NEW trees (I>ITEMP), or FTEMP−PNN(NCOUNT) for the continuation
  increment (I≤ITEMP). NOT the ingrowth FTEMP·(NEWTPP/ITPP).
- Per-tree TPA booked = ESPROB·300/DUPNPT (estab.f:1459). Live: 0.60122·300/50 = 3.607/tree.
- **Cap = MAXTPP (25), NOT MAXING (7)** — live ITPP reaches 14. jl's `cap = _IE_MAXING[ihab]` (establishment.jl:676)
  is the workaround from the AUTOES-v1 approximation (comment lines 671-675).
- Multiple disturbance tallies re-stock over cycles at rising PROB1 (0.601 → 0.817), NSTORE carrying between
  them so each books only ITPP−NSTORE.

**jl measurement (iet01, AUTOES_DBG):** prob1 = 0.60012 (≈ live 0.60122 — elev 34, matches), is_ingro=false,
but **total = 727.3 vs live ~590 = jl OVER-produces ~23%.** So the residual is NOT under-production from the
MAXING cap (jl over-produces even capped LOWER than live) — prob1/per-tree match, so the gap is **sum(ITPP): jl's
ESTPP-drawn ITPP sequence is too high vs live for the disturbance path.** This is an ESTPP-draw/ITPP divergence
(the EM ingrowth body=135→83 fix is INERT for IE: nsp=23/ihab=10 → 16+69+50=135 unchanged), OR a plot-count /
NSTORE-continuation handling difference.

### RESOLVED piece — disturbance cap MAXTPP (2026-08-06, committed)
Instrument-replay (FVSie iet01 ESTPP DRAW dump) settled the disturbance ITPP: **jl's ESTPP draws + TPP are
BIT-IDENTICAL to live** (0.346302→2.486, 0.184747→1.443, 0.835307→14.036, 0.243140→1.755, 0.462182→3.656,
0.230992→1.684, 0.967520→34.495, 0.513879→4.337). The RNG/seed chain is fully synced (my earlier "727 vs 590"
was a stand mismatch). The SOLE divergence: jl capped the disturbance path at MAXING(7) too — truncating ITPP
14→7 and 25→7 ⇒ UNDER-production. estab.f:682 gates MAXING on `INGRO.EQ.1` only; disturbance uses MAXTPP(25).
FIXED (`cap = is_ingro ? MAXING : MAXTPP`) ⇒ jl disturbance ITPP now = live bit-exact [2,1,14,2,4,2,25,4].
Live tally-1 total TTOTTP=1079.33 (was jl 727 at MAXING(7); now matches at MAXTPP). EM ingrowth unchanged (118).

### STILL OPEN (IE establishment residuals, measured 2026-08-06)
1. **iet01 stand-1 ~2.7%** (2000 jl 429 vs live 441): stand-1's establishment does NOT hit the high-ITPP capped
   plots (unchanged by the cap fix) ⇒ a SEPARATE residual — likely the multi-tally NSTORE re-stock (live fires
   NTALLY=1 then NTALLY=2 at rising PROB1 0.601→0.817, each booking ITPP−NSTORE) or the established-cohort
   subsequent DG/mortality. Needs stand-level tally alignment (which of the 4 iet01 stands owns which tally).
2. **Ingrowth path over-produces** (jl tally 629 vs live TTOTTP 109 on the iet01 INGRO=1 tally): the is_ingro
   per-tree TPA (prob1·newtpp/itpp) or the NSTORE continuation runs high — analogous to the EM ingrowth amount
   but with body=135 already correct (IE nsp=23/ihab=10). ROOT likely the ingrowth prob1 (PLPROB-driven
   ITPP=INT(PLPROB·DUPNPT/(FTEMP·300)) at line 589, NOT the ESTPP draw — the ingrowth path OVERRIDES ITPP with
   the PLPROB formula) or the FTEMP·NEWTPP/ITPP weighting. NEXT: instrument live PLPROB(NNID) + the line-589
   ITPP override for the INGRO tally, compare to jl's ingrowth ITPP/newtpp. Also the prob1 0.60012 vs 0.60122 hair.
