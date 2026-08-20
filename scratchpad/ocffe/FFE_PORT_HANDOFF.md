# OC/OP FFE port — turnkey handoff (updated 2026-08-20)

USER "do all user-gated" ([[feedback-do-all-user-gated]]). OC/OP FFE is TRACTABLE (live oracle exists,
produces a full FFE .sum + ALL FUELS REPORT). Multi-chunk subsystem, ported chunk-by-chunk bit-exact.

## ORACLE (validatable — DOES produce FFE output despite the vol-path crash)
- `.sum`: `scratchpad/ocffe/ocffe_oracle.sum` — 13-cycle FFE stand; 2003 SIMFIRE kills TPA 261→72.
- ALL FUELS REPORT: `scratchpad/ocffe/ocffe_oracle_fuels.txt` (the chunk-2 target).
- Reproduce: `cd /workspace/.ocwork && ./FVSoc_clean --keywordfile=oct01_full.key` (oct01_full.{key,tre}
  are the tests/FVSoc/oct01.{key,tre} copies). FFE stand = the `-999 "FFE"` .sum block; fuel report at
  `oct01_full.out` line ~3100.
- jl harness: `scratchpad/ocffe/ocffe_full.key` + `ocffe_full.tre` = a SELF-CONTAINED FFE stand (the real
  TREEFMT inlined — the extracted stand-4-only key lacked it, which mis-read the .tre as d=0/sp=49 garbage
  → the "forest-711 CCF NaN" was a keyfile-extraction artifact, NOT an OC bug).

## CHUNK 1 — crown-biomass + bark + snag species-props  ✅ COMMITTED 342b1afc (gate 339/11)
- `data/oregoncoast/fire_species_props.csv` (50 sp): v2t/leaf_life/dkr_cls/snag_fallx/snag_alldwn from
  oc/fmvinit.f SELECT CASE; snag_decayx=999 (fmvinit.f:484); biogrp/ls_spi = OC_FFE_ISPMAP. Auto-loaded by
  load_species_coefficients (fire_species_props.csv is in its merge list).
- `_OC_FM_BARK_B1` (oc/fmbrkt.f) + OC branch in `fire_bark_thickness` (OP=39sp separate) (fire_effects.jl).
- crown_biomass.jl: OC bark via `oc_bratio` (ORGANON BRATIO); ls_spi/biogrp via OC_FFE_ISPMAP.
- KNOWN CORNER: tfall_cls=1 placeholder. OC TFALL(I,3)=15/20yr EXCEEDS the SN class table (_FM_TFALL3 max
  10, fuel_additions.jl). OC sets TFALL per-species directly in fmvinit.f (TFALL(I,0)=1/3, (I,1)=10,
  (I,2)=15, (I,3)=15/20, (I,4)=(I,5)=(I,3)). ⇒ a variant-aware `_fm_tfall` w/ an OC TFALL matrix is a
  sub-chunk (affects CWD2B multi-year crown-debris fall, not cyc0).

## CHUNK 2 — fmcba initial fuel loading  ⏳ SCOPED, tables extracted (next)
Tables extracted bit-exact (parser-verified vs oc/fmcba.f DATA), committed as data CSVs:
- `data/oregoncoast/fire_fuel_covtype_dead_est.csv`  = FUINIE (50 sp × 11 size classes; established, ~60% cover)
- `data/oregoncoast/fire_fuel_covtype_dead_init.csv` = FUINII (50 sp × 11; initiating, ~10% cover)
- `data/oregoncoast/fire_fuel_covtype_live.csv` = FULIVI(herb,shrub) + FULIVE(herb,shrub) per species
Size classes (11): <.25, .25-1, 1-3, 3-6, 6-12, 12-20, 20-35, 35-50, >50, LITTER, DUFF.

**OC fmcba is NOT the generic ffe_dead_fuel_loading(ifortp) path** (that indexes by FIA forest type). OC
indexes by COVTYP = dominant SPECIES, blended est↔init by canopy cover. Needs an OC fmcba dispatch branch
(like CS/NE have their own). Algorithm (oc/fmcba.f:425-770), MEASURED:
  1. FMTBA(ksp) = Σ FMPROB·DBH²·0.0054542 (basal area per species).
  2. RDPSRT(FMTBA) → ICT; COVTYP=ICT(1); COVCA(1..2)=top-2 species; COVCAWT(j)=FMTBA(ICT(j))/Σtop2.
  3. TOTCRA = Σ (π·CW²/4)·FMPROB;  PERCOV = 100·(1−exp(−TOTCRA/43560)).   [CW = CRWDTH crown width]
  4. LIVE:  FLIVE(i=herb,shrub) = Σ_{j=1,2} ALGSLP(PERCOV, XCOV=[10,60],
              YLOAD=[FULIVI(i,COVCA j)·COVCAWT j, FULIVE(i,COVCA j)·COVCAWT j], 2).
  5. DEAD:  STFUEL(sz,2) = Σ_{j=1,2} ALGSLP(PERCOV, [10,60],
              [FUINII(sz,COVCA j)·wt, FUINIE(sz,COVCA j)·wt], 2)   for sz=1..11.
     (ALGSLP = piecewise-linear interp w/ flat extrapolation; XCOV increasing.)
  6. No-tree cyc0 (COVTYP=0): IFOR≥6 → COVINI6(ITYPE) [R6 habitat]; else COVINI5(ITYPE) [R5 plant assoc];
     fallback COVTYP=7 (DF). Post-cyc0: OLDCOVTYP. (COVINI5/6 maps still to extract, fmcba.f:355+.)
  7. DKRT decay rates (first year; multi-cycle only, NOT cyc0 loading): OC forest 711 = KODFOR≥600 = R6
     Oregon branch (fmcba.f:610-655): base DKRT(size,dkrclass) then ×DKRADJ(TEMP=CAHMC(ITYPE),
     MOIST=CAWMD(ITYPE), K∈{1,2,3 by size}), clamp≤1, monotonic-bump. Litter loss 0.5/yr, duff 0.002/yr,
     PRDUFFT=0.02. (R5 KODFOR<600 adds a Dunning-SI DCYMLT — NOT needed for OC 711.) ⇒ chunk 2b.

**cyc0 target (ocffe_oracle_fuels.txt, 1993):** LITT 0.60, DUFF 15.7, 0-3"(sz1-3) 3.7, >3"(sz4-9) 9.2,
HERB 0.21, SHRUB 0.57, SURF TOTAL 29.9. Validate by instrumenting jl to dump STFUEL/FLIVE at cyc0 (no .sum
until fire behavior ports). Needs: ALGSLP helper + OC fmcba branch reading the 3 new CSVs + PERCOV.

## CHUNK 3 — RUNS END-TO-END ✅ (fbfba49f). cyc0 .sum BIT-EXACT; fire-mortality gap = refinement.
Three fixes cleared the fmburn path: fire_fuel_models.csv (standard Anderson-13, variant-independent);
is_sprouting (blkdat.f ISPSPE {24,26-48,50}); bio_group/biogrp = Jenkins BIOGRP (fmcblk.f), NOT ISPMAP.
FIRST END-TO-END .sum (ocffe_jl.sum vs ocffe_oracle.sum):
  • **1993 (cyc0) BIT-EXACT** — every column. • 2003 SIMFIRE mortality **206 vs oracle 266** (jl
    under-kills) = the open FFE gap. • pre-fire CCF 76 vs 79 = the KNOWN OC ORGANON multi-cycle growth
    drift (documented in [[fvsjl-oc-organon-blm-volume]], NOT FFE). • cyc1 accretion col 105 vs 110.
REMAINING (the 206-vs-266) — **FULLY TRACED 2026-08-20 to SURFACE-fire intensity (NOT crown fire).**
Debug at the 2003 event: jl byram=1189.9 → flame=1.78ft → scorch=4.7ft, at wind=10/temp=50/fmois=1.
  • **The fire is SURFACE, both sides.** Oracle POTENTIAL FIRE REPORT 2003 = type "S" (torch index 56 >
    wind 20; crown index 47.5 > 20). My earlier "CRWNG=1 = crown fire" read was WRONG. jl's crown-fire
    path (now wired, 69141cdb) correctly returns crb=0 (surface). Crown fire is a RED HERRING.
  • **jl's WEATHER is CORRECT.** FVS `SIMFIRE 2003 10.00 1 50.0` explicitly sets wind=10/fmois=1/temp=50
    (fmin.f:326-336 defaults 20/1/70, overridden by fields 2/3/4). The report's SEVERE col (wind 20/temp
    70) is a what-if diagnostic, not the actual fire. jl uses the right weather.
  • **The DKRT decay is .sum-INERT here** (9cce8d0b): with the correct OC decay the 2003 mortality is still
    206. So the surface-fuel decay is not the lever either (still correct OC physics for fuel-limited stands).
  ⇒ **The under-kill = the SURFACE byram (fuel intensity):** jl's byram=1189.9 gives flame 1.78/scorch 4.7,
    slightly too weak to kill the 266. Two candidate levers, both chunk-2 fuel precision:
    (a) **PERCOV precision** — jl 49.8 vs oracle ~47.6 ⇒ the ~10% LITT/shrub gap ⇒ lower surface load ⇒
        lower byram. Pin exact PERCOV by INSTRUMENTING oc/fmcba.f (WRITE PERCOV; g16 single-.o — FMDEBUG
        keyword CRASHES FVSoc). Likely a crown_pct(CR)/tree-state detail in oc_cwcalc or BAREA=1-vs-actual.
    (b) **fmcfmd fuel-model selection** — the oracle weights models 6(62%)/10(38%) at 2003; if jl's
        select_fuel_models picks weaker models, byram drops hard (discrete, high-leverage). DUMP jl's
        selected (fm,weight) at 2003 and A/B vs the oracle's "6 62 10 38" (POTENTIAL FIRE REPORT tail).
    Get the oracle's actual-fire flame/byram (instrument fmburn/fmfint at wind=10) as the exact target.
## CHUNK 4 — snag dynamics (SNAGINIT/SNAGOUT) → STANDING WOOD columns (cyc0 snag already bit-exact in .sum).

## DOCTRINE
Bit-exact-or-cornered per chunk by RUNNING FVSoc_clean; glibc libm ccall for Float32; gate = multicycle
339/11 byte-identical. Wire Olympic (OP) analogously — same 39-vs-50 species list? CHECK op/fmvinit.f,
op/fmbrkt.f, op/fmcba.f (OP MAXSP=39; the tables above are OC's 50 — OP needs its own extraction).
