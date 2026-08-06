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

**NEXT** (the #143 IE close-out): instrument live IE estab.f to dump the per-plot ESTPP DRAW + ITPP for the
disturbance tally (like the EM EMDRAW measurement), align to jl's per-plot draws, and root why jl's ITPP runs
high. Candidates: (a) jl's disturbance cap MAXING truncates the wrong plots vs MAXTPP; (b) jl fires the tally on
more plots / a different NSTORE-continuation than live's multi-tally re-stock; (c) the ESTPP draw itself diverges
for the disturbance (NTALLY=1) seed state vs the validated ingrowth path. Also validate the prob1 0.60012 vs
0.60122 hair (a tiny ESTOCK input diff — likely the disturbance TIME/REGT or elevation-hundredths).
