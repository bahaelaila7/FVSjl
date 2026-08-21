# ESTAB modifier keywords (STOCKADJ / TALLY / NOINGROW) — CORRECTED 2026-08-21: PRIOR VALIDATION WAS FLAWED

## What the code does (FAITHFUL transcription — this part is correct)
- STOCKADJ (esnutr.f IACTK 440): est.stoadj = PRMS(1); ie_autoes_run applies estab.f:578-580 —
  clamp STOADJ>=0.001, then PROB1 = logistic(PN+ESB-ESB1)*STOADJ. (committed eddf8dc2)
- TALLY/TALLYONE/TALLYTWO (esin.f 16/11/12): kw_estab! pushes ScheduledActivity(date,427/428/429);
  ie_autoes_establish! honors a scheduled tally in-window. (committed 1b9ea7d9)
- NOINGROW (esin.f opt 22): est.lingrw=false; ie_autoes_schedule! path-3 gated on est.lingrw.
  NOTE also handled only at TOP-LEVEL (keyword_dispatch.jl:2422), NOT inside kw_estab! — so
  NOINGROW-inside-the-ESTAB-packet (FVS-correct placement) is silently skipped by jl.

## THE FLAW (measured 2026-08-21, definitive clean A/B)
On the ONLY fixture (dense stand 12343703010690, ~28000 TPA, THINPRSC 2029), a SAME-CODE jl A/B
(ie_sa_base vs ie_sa_050) is BYTE-IDENTICAL every cycle: Delta(STOCKADJ 0.5 - base) = 0/0/0/0/0/0.
STOCKADJ is INERT in jl here. The oracle (FVSie_clean) A/B shows Delta = 0/0/-88/-75/-205/-172.
NOINGROW: jl trace confirms it mechanically disables LINGRW (all paths NONE, itrn frozen) but the
resulting .sum Delta (0/+35/+140/+227) DISAGREES IN SIGN with the oracle (-88/-75/-205/-172).

ROOT CAUSE: on this already-stocked dense stand jl's #143 ingrowth clamp books ~0 ingrowth
(NEWTPP = max(0, ITPP - NSTORE) = 0 because NSTORE >= ITPP), so scaling/disabling the ingrowth
tally changes nothing in jl. The ORACLE books ~88 ingrowth trees that these keywords modulate.
=> jl and the oracle DIVERGE on the underlying establishment-ingrowth booking on dense IE stands
(jl UNDER-books vs oracle). The baseline itself is NOT bit-identical (jl 2027 BA 147 vs oracle 168,
mort 95 vs 121) — a dense-establishment-regime divergence, NOT the clean baseline the prior doc claimed.

## Status of the prior "bit-exact" claim
FLAWED. The prior session reported "Delta(sa050-base) 2027 -88/-88, jl matches oracle every cycle"
and "jl baseline also bit-identical" — BOTH are FALSE (jl Delta is 0; jl baseline diverges). The -88
was the ORACLE's own delta, erroneously attributed to jl (the jl A/B was evidently never run cleanly).

## Disposition
- The KEYWORD CODE (est.stoadj, prob1*sa, TALLY scheduling, est.lingrw gate) is a faithful esnutr.f/
  estab.f transcription, gate-safe (339/11, inert when keyword absent) — KEPT, would be correct on a
  stand where jl books ingrowth. NOT independently validated on any available fixture.
- REAL OPEN ITEM (the actual porting work): jl's IE AUTOES establishment-ingrowth booking diverges
  from the oracle on dense stocked stands (#143 NSTORE clamp interaction) — a #143-class deep dive
  needing a dense-IE fixture + oracle per-plot NSTORE/NEWTPP dump. Until that converges, NONE of the
  ingrowth-modulating ESTAB keywords (STOCKADJ/TALLY/NOINGROW) are validatable on this stand.
