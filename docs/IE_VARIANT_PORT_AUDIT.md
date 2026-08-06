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
2. **Ingrowth path over-produces** (jl tally 629 vs live TTOTTP 109 on the iet01 INGRO=1 tally). ROOT MEASURED
   (estab.f:589-590,683 instrument-replay): the ingrowth pre-loads NSTORE = INT(PLPROB·DUPNPT/(FTEMP·300)+0.5)
   per plot, so NEWTPP = ITPP−NSTORE books only the increment above the standing regen. Measured iet01 ingrowth:
   PLPROB=13.651, FTEMP=0.5949, DUPNPT=50 → NSTORE=4 (NEWTPP per plot 0,3,0,3…); disturbance NSTORE=0 (PLPROB
   0.01). PLPROB = Σ(small-tree TPA<REGNBK=2.999)/DUP, DUP=5 (nptids). jl reset nstore=0 ⇒ booked full ITPP.
   ★ FIX ATTEMPTED + REVERTED (doctrine #4): pre-loading nstore = the PLPROB formula OVER-corrected (jl tally
   →0.1 vs live 109) because **jl reads iet01 as nptids=1 (live DUP=5) and its small-tree TPA is ~2× live's**
   (jl PLPROB 135 vs live 13.651 = 10×: 5× from nptids, 2× from small_tpa). This is a deeper PLOT-STRUCTURE
   discrepancy (jl collapses iet01's 5 inventory points to 1; dupnpt=50 still matches via idup=50 vs 10, so the
   disturbance path + growth were unaffected — only the per-point PLPROB exposes it). The 2× small_tpa is likely
   a COMPOUNDING artifact (the ingrowth tally fires late, after prior establishment adds regen). NEXT: resolve
   why jl reads iet01 nptids=1 vs live 5 (points_inv from the .tre plot IDs / DESIGN card), then the PLPROB
   pre-load formula (confirmed correct) drops in. Blocked on the plot-count read, not the establishment formula.

## Real-FIA IE dense-seedling sweep (2026-08-06) — CRVAR cross-fix + WH under-growth lead
Built scratchpad/ie_test.db (3 dense-seedling IE stands from the 70GB DB, VARIANT='IE'). vs FVSie_clean, 2 cyc:
- **12343703010690** (conifer seedlings PICO/spruce/larch, 35484 TPA): jl cyc1 28915/16 vs live 28930/15 ✓ (SMDGF works).
- **44987944020004** (37528 TPA): jl 32203/154 vs live 32104/164 — cornered (−6% BA).
- **753199439290487** (sp263=WH, 50377 TPA): jl 45030/13 vs live 45880/74 — ★ jl UNDER-grows ~5× (QMD 0.23 vs 0.54).
IE CRVAR CO diameter fix (219d04f) is INERT on these (none are sp19/22 CRVAR) ⇒ non-regressive but unconfirmed-on-target;
it is source-faithful (ie/regent.f:637-640 = the EM-validated block) + iet01-clean. The WH under-growth (#146) is a
SEPARATE new lead (conifer, opposite direction, dense) — NOT the CRVAR fix (jl identical with/without).

### WH under-growth [#146] — subcycle-1 HTGRL + coefficients MATCH (2026-08-06, paradox localized)
Instrument-replay (ie/regent.f:500 vs jl regent.jl:194) on stand 753199439290487 WH seedling (D=0.1, H1=1.01):
jl HTGRL=0.68317 vs live 0.67374 (jl slightly HIGHER, con=0.7259 both, h1 match); rdj jl 11.79 vs live 13.92
(jl LOWER density ⇒ LESS suppression). HHT1[5]/HHT2[5] (H-D power) BIT-MATCH (0.0729/1.1988). ⇒ PARADOX: every
measured subcycle-1 component says jl should grow ≥ live, yet jl under-grows 5× (QMD 0.23 vs 0.54). The divergence
is therefore NOT height-suppression/coefficients — it is downstream: the MULTI-SUBCYCLE density-feedback trajectory
(jl's rdnext/banext accumulation across subcycles diverging so later subcycles over-suppress) OR the final DDS→DG
conversion (regent.jl:414-428 dk→dds→dg_inc). NEXT: dump the FULL per-subcycle (h1,h2,rdj,d2) trajectory for the
WH seedling jl vs live — the height must reach ~9.7' in live (for d=0.54) vs jl's ~6.5' (d=0.23), so the height
DIVERGES over subcycles 2+ despite matching at subcycle 1 ⇒ the rdnext density-feedback compounding is the prime suspect.
