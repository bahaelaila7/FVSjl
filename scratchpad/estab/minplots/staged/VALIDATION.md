# MINPLOTS / PASSALL / PLOTINFO — ESTAB-packet model-effect keywords (unblocked by AUTOES #143)

These three were consumed-INERT ("DEFERRED pending AUTOES #143") in kw_estab!. With the draw stream now
bit-exact (081e0e6c) each was re-examined against the LIVE oracle on the AUTOES-exercising under-stocked IE
fixture (753189105290487, ihab 9, NPTIDS 50, ie_understocked.db, plain ESTAB, NUMCYCLE 5).

Oracles: FVSie_g16 (clean .sum A/B) and FVSie_estab2 (instrumented per-tally ESDRAW seed dump).

## MINPLOTS — WIRED BIT-EXACT

esin.f opt 20 → MINREP (esinit.f default 50; esin.f:587-588 IF(MINREP<20) MINREP=20). Drives the
estab.f:199-207 IDUP loop: IDUP = smallest I with NPTIDS·I ≥ MINREP = ceil(MINREP/NPTIDS) ⇒
DUPNPT = NPTIDS·IDUP = the number of establishment plots looped ⇒ the ESRANN draw-stream length.

Model-affecting on the .sum (FVSie_g16, TPA @2029): default 815 vs MINPLOTS 100 → 761. And bit-exact on the
per-tally ESDRAW seed chain (FVSie_estab2), which the base .sum masks behind the pre-existing NSTORE
under-stocked-seedling corner:

| MINPLOTS | idup | dupnpt | tally-1 | tally-2 | tally-3 |
|----------|------|--------|---------|---------|---------|
| default 50 | 1 | 50 | 43303 | 61997 | 49053 |
| 100        | 2 | 100 | 43303 | **31334** | **6152** |
| 150        | 3 | 150 | 43303 | **53059** | **15474** |

jl reproduces ALL three seed chains byte-for-byte (FVSJL_AUTOES_DEBUG=1 on each keyfile). tally-1 is invariant
(the stand-level ESDRAW before the plot loop); tally-2/3 shift because the ESAVE chain advances through DUPNPT
plots not 50.

**Wiring (4 hunks):**
- `state.jl` Establishment: new `minrep::Int32` field (default 50).
- `keyword_dispatch.jl` kw_estab!: MINPLOTS branch parses `MINREP=IFIX(ARRAY(1))`, clamps `<20→20`, stores.
- `variants/inlandempire/establishment.jl` + `engine/establishment.jl`: `idup = cld(est.minrep, nptids)`
  (was the hardcoded const `_ES_MINREP`).
- **The load-bearing fix:** pass `wk6 = DUPNPT` (was the stale default 50) to `ie_autoes_plot_seeds` at BOTH
  call sites (the es_stream advance + `ie_autoes_tally`). The one-time WK6 site-prep fill is
  `DO 183 I=1,IDUP*NPTIDS` = DUPNPT draws; with wk6 stuck at 50 the idup>1 chain desynced (jl 4681/17715 vs
  oracle 31334/6152). With wk6=DUPNPT it's bit-exact. Inert for the default (dupnpt=50 ⇒ wk6=50).

Test: `test/unit/test_estab_minplots.jl` (12/12) — parse+clamp + the DUPNPT-scaled ESAVE chain + next-tally seed.

## PASSALL — MEASURED-INERT on the reachable regime

esin.f opt 18 → IBLK=1, CONFID=max(ARRAY(1),1.0) → PASMAX (esinit.f default CONFID 5). PASMAX caps ONLY the
per-plot EXCESS-tree probability (estab.f:1318-1321 XCSMAX=min(EXCESS/BRKUP, PASMAX/BRKUP)) — a DETERMINISTIC
post-draw PROB scaling, no RNG effect — and only when a plot overflows MAXTPP so EXCESS(I)≥0.5. IBLK is set but
referenced NOWHERE in estab.f (dead serialization flag, getstd/putstd only).

A/B (FVSie_g16, under-stocked fixture): PASSALL 1 == PASSALL 100 == default → .sum **BYTE-IDENTICAL** (all 6
cycles). The under-stocked stand never overflows a plot, so there are zero excess trees to cap. Consumed-inert;
exercising it would need a dense over-regenerating fixture that also fires AUTOES (untried, not the AUTOES stand).

## PLOTINFO — STRUCTURALLY INAPPLICABLE (legacy plot-card input)

esin.f opt 10 → IPINFO=1 + READs per-plot site CARDS (ID,slope,aspect,habitat,physio,prep) from the keyfile via
esplt1.f/esplt2.f — the LEGACY TREEDATA plot-card input format. jl's DATABASE/FIA pipeline takes per-plot
SLOPE/ASPECT from the DB reader (point_slope/point_aspect); IPINFO stays 0 and a real jl keyfile never carries
PLOTINFO plot cards. Validating it bit-exact would require porting the entire legacy plot-card reader plus a
non-DB keyfile — outside the reachable DB regime and off the AUTOES critical path. Consumed-inert (documented).

## Gate

`test/integration/test_multicycle.jl` = **339 pass / 11 broken** byte-identical (all gate scenarios SN/BARE; the
change is IE-AUTOES only and inert at the default MINREP=50). `test_ie_estock.jl` + `test_estab_specmult_htadj.jl`
still pass (dupnpt=50 compose bit-exact 583.7).

## Reproduce
```
# oracle .sum A/B (MINPLOTS moves it; PASSALL does not):
for k in base minplots100 passall1 passall100; do echo $PWD/$k.key | /workspace/.iework/FVSie_g16; done
# oracle per-tally seed chain:
echo $PWD/minplots100.key | /workspace/.iework/FVSie_estab2; grep ZESDRAW /workspace/.iework/{estabdump.txt,fort.91}
# jl seed chain:
FVSJL_AUTOES_DEBUG=1 julia --project -e 'using FVSjl; FVSjl.run_keyfile("…/minplots100.key"; variant=FVSjl.InlandEmpire())'
```
