# FVS Extensions Rollout Plan — FFE / ECON / Mistletoe / Climate across the 12 ported variants

User directive (after CI done, BC re-scoped): "port all of the extensions to those two and the previous
ones" → chosen sequencing "Extensions first (11 done + CI), defer BC's full BEC port." This plan drives
the FFE (+ ECON/mistletoe/climate) rollout across the 12 variants that have a validated growth core:
**sn, ne, cs, ls (eastern) + cr, ie, kt, em, bm, tt, ut, ci (western)**. Validate each per doctrine:
bit-exact-or-cornered vs the LIVE relinked binary, one chunk/variant at a time, instrument don't infer.

## Current extension coverage (measured, src/engine/fire + econ.jl + variants/*)
| Extension | Have | Missing |
|---|---|---|
| FFE (fire/fuels/carbon) | **CR** (western template, 23 dispatch pts) + **SN/CS/LS/NE** (eastern) | **ie, kt, em, bm, ut, tt, ci** (7 western) |
| ECON | econ.jl (variant-general: log-graded HRVRVN, discount) — validated SN/eastern | confirm western wiring (should be ~drop-in; ECON is variant-agnostic) |
| Dwarf mistletoe | CR (variants/centralrockies/dwarf_mistletoe_model.jl); IE MISTOE done (per memory) | per-variant DMR tables for other western (situational) |
| Climate-FVS | none in tree | all (lowest priority — climate is an optional keyword layer) |

## FFE rollout — the dominant phase. Order set by measured table reuse (bin/FVS*_buildDir/fm*.f diffs):
- **fmmois.f (moisture): ALL western identical (ie=kt=em=bm=ci, 0 diff).** Port the western moisture table ONCE.
- **fmcfmd.f (cover-type→fuel-model, the big table ~1250 lines): ie/kt = 0 (IDENTICAL).** em=205, ci=299,
  bm=622 differ from ie. ⇒ **IE port auto-covers KT** for the fuel-model half.
- **fmbrkt.f (bark→fire mortality): per-species DBH·B1, small tables.** ie=23 sp, kt=11 (=first 11 of ie).
  em/bm/ci differ by 36-83 lines (species-count). fmsnag.f: CR=IE identical (shared snag base).

### Sequence
1. **IE** (N-Rockies template) — full western FFE: fmmois(western-shared) + fmcfmd + fmbrkt + crown biomass +
   snag + fuel loading + dispatch `s.variant isa InlandEmpire`. Validate vs live FVSie SIMFIRE on iet01.
2. **KT** — nearly free after IE: fmcfmd identical, fmbrkt = IE[1:11], fmmois shared. Just add dispatch + 11-sp table.
3. **EM, CI, BM, UT, TT** — each: own fmcfmd (fuel-model), own fmbrkt, shared fmmois, crown-biomass species map.
4. Then ECON confirm-western + mistletoe (per-variant DMR, situational) + climate (optional, last).

### Per-variant FFE data needed (the extension points, mirrored from CR in src/engine/fire/*.jl)
- fire_effects.jl: bark B1 (fmbrkt) → mortality.  ✔ IE+KT bark DONE (this commit).
- fuel_moisture.jl: `_FM_MOIS_<v>` (fmmois) — western-shared, port once.
- fuel_model.jl: `_FMD_XPTS_<v>` breakpoints + fmcfmd cover-type→fuel-model map + fuel loading (FUINIE/FUINII).
- fmcba.jl: cover-type default fuel model, live/dead fuel loading, large-fuel photo-series species default.
- crown_biomass.jl: species→crown-biomass-eqn map (ISPMAP) + bratio (already have per-variant bark_ratio).
- fuel_decay.jl: `_FM_DKR_<v>` decay rates.  snag.jl: snag fall/decay (fmsnag shared with CR).
- fmburn.jl: conifer-species mortality gate + flame/byram kill.

## IE fuel-model (fmcfmd) chunk — fully mapped (ie/fmcfmd.f is only 105 lines; the 1250-diff was structural)
IE's FMCFMD = candidate standard fuel models → shared FMDYN interpolation. Pieces:
1. ✔ **XPTS** `_FMD_XPTS_IE` (ie/fmcfmd.f:22-36, ICLSS=14, model10=(15,30), model14=(30,60)) — DONE; KT identical.
2. **IDRY class** from `MAPDRY(2,·)` habitat→dryness table (fmcba.f:423-446): 1=dry-grassy, 2=dry-shrubby, 0=other.
   (`NIFMHAB(IDRY)` just reads IDRYB set there.) → extract MAPDRY for IE (habitat ITYPE → IDRY). TODO.
3. **EQWT candidate weights** (ie/fmcfmd.f:57-99): CASE(IDRY): 1→models{1,9}; 2→models{2,9} weighted by PERCOV
   via ALGSLP(PERCOV,[30,50]); DEFAULT→model 8 (=1.0). Then post-harvest activity fuels:
   AFWT=max(0,1-(IYR-HARVYR)/5); if SLCHNG≥SLCRIT or LATFUEL → EQWT(11)=EQWT(14)=AFWT; EQWT(10)=EQWT(12)=1-AFWT;
   EQWT(13)=1. → build the IE EQWT vector, hand to the ported `_fmdyn` with `_FMD_XPTS_IE`. TODO.
4. **fmcba fuel loading** (live/dead by cover type) + **crown_biomass** species→eqn map for IE. TODO.
5. **Validate**: iet01 SIMFIRE keyword .sum (fire-year TPA/mortality/surface-fuel) bit-exact vs live FVSie
   (relink FVSie oracle in /workspace/.iework). Then KT drops in (fmcfmd identical, bark=[1:11], mois shared).

## Progress log
- **2026-08-03** Assessed architecture; set order (IE→KT free, western-shared moisture, then EM/CI/BM/UT/TT).
  Committed the full IE/KT surface-fire path (7 commits): (1) bark-thickness `_IE_FM_BARK_B1`; (2) moisture
  `_FM_MOIS_IE` (covers ie/kt/em/bm/ci); (3) fuel-model XPTS `_FMD_XPTS_IE`; (4) candidate selection
  `ie_select_fuel_models` (MAPDRY→IDRY + PERCOV weights + natural fuels); (5) fmcba fuel loading
  (`data/inlandempire/fire/ffe_fuel.jl` FULIVE/FUINIE + COVINI); (6) snag/biomass species props
  (`fire_species_props.csv` — v2t/decay/fall from fmvinit, ls_spi=ISPMAP, biogrp=BIOGRP; note BIOGRP grab
  contaminated by "! NN" comments, corrected from source) + standard Anderson-13/Jenkins tables + IE bark in
  crown_biomass; (7) **fuel-decay fix** — IE/KT fell through to the SN `_FM_DKR`; wired to `_FM_DKR_CR`
  (ie/fmcwd.f == cr/fmcwd.f).
  - **IE SIMFIRE VALIDATED-CORNERED vs live FVSie** (iet01_fire.key in /workspace/.iework/ierun). Growth
    bit-close (2020 jl 302 trees/BA 178 vs live 305/178). Fire mortality (2020) driven from 45% → 83% → **99%**
    of the stand by two root-caused fixes:
    1. **decay-table fallthrough** (IE/KT used SN `_FM_DKR` → wired to `_FM_DKR_CR`): MOR 145 → 322.
    2. **crown-fire flame boost** (fmburn crown-fire adjustment was gated CR/NE-only; ie/fmburn.f is
       byte-identical to cr → admitted IE): MOR 322 → **500** (live ≈512), 2030 survivors 51 → **3** (live 0).
    jl now kills 99% matching live's total kill. **Residual = 3 trees/ac survive (MOR 500 vs 512, 2.3%)** — the
    cornered tail (the last few large trees at the scorch boundary; likely a bark/scorch-height ULP or RNG-tie
    at psburn=100). This is bit-exact-or-cornered per doctrine. The surface-fire + crown-fire + mortality path
    for IE is COMPLETE. Optional refinement: instrument the 3 survivors' scorch fraction vs live if exactness
    is wanted. **Then KT drop-in**: needs KT fmbrkt[1:11] bark + kt_bratio in crown_biomass + admit KT to the
    crown-fire gate/Union; fmcfmd/moisture/decay/fuel-loading already shared.
- **Watch-list carried to EM/CI/BM/UT/TT**: each western variant needs its own fmcfmd MAPDRY + fmcba
  FULIVE/FUINIE + fmvinit fire_species_props + verify its `_FM_DKR`/moisture vs CR (don't let them fall
  through to SN defaults — the decay-fallthrough bug class).
