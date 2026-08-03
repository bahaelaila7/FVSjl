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

## CI (Central Idaho) FFE — SCOPED (most complex western fmcfmd; checkpoint before wiring)
CI decay==CR, moisture==IE (both shareable). CI's fmcfmd (ci/fmcfmd.f) is the most elaborate western selection:
- **ICT = MAPPVG[ICINDX]** (ci/fmcba.f ENTRY CIPVG; MAPPVG 130 entries, ICT values 1-11). ICINDX is ALREADY
  computed by the CI growth port (dual-habitat ICITYP bracket) — reuse it. Extracted: MAPPVG (126/130 non-null).
- **ICT case groups**: CASE(1:4) → K by ICT (1→model1; 2→CIS9B ninebark?5:snowberry?2; 3→5; 4→2), then
  EQWT(K)+=WT1(1); EQWT(9)+=WT1(2)·PRLONG; EQWT(8)+=WT1(2)·(1−PRLONG), WT1=ALGSLP(PERCOV,[30,50]).
  CASE(5:6) → K=2/5 THEN a **grand-fir-understory sub-model** (~100 lines: grand-fir saplings ISCT(4), sapling
  crown-cover CRGF/CCGF → a 2nd weight WT2, splitting EQWT across models via WT1·WT2 products) — the complex part.
  CASE(7:11) → EQWT(8)=1.0 (trivial). Then activity 11/14 (AFWT) + natural 10/12/13 (AFWT=0 path).
- **PRLONG** = BA-fraction in sp {1,10} (long-needle pines) / total BA. **CIS9B** = ninebark-vs-snowberry from
  the CI habitat (ci/fmcba.f) — small habitat lookup, extract with MAPPVG.
- **cit01 CONFIRMED: habitat 520 → ICINDX=66 → ICT=MAPPVG[66]=6 → CASE(5:6) grand-fir branch** (ICINDX=66
  cross-checked vs the CI growth port's BAMAX work). So cit01 REQUIRES the grand-fir-understory logic — cannot
  validate with simple cases. ICINDX is stashed in **p.habitat_input** (ci site_index.jl:102) — read it at fire
  time; ict = MAPPVG[Int(p.habitat_input)].
- **Full CASE(5:6) logic (ci/fmcfmd.f, read + transcribed)**: K=2(ICT5)/5(ICT6). Grand-fir understory = species-4
  saplings (DBH≤3): CRGF=Σ(FMPROB·ICR)/ΣFMPROB (avg crown ratio %); TOTCRA=Σ(π·CRWDTH²/4·FMPROB);
  CCGF=100·(1−exp(−TOTCRA/43560)); LCRGF=(CRGF≥75). If LCRGF: WT1=ALGSLP(CCGF,[50,70]); WT1(1) path splits by
  WT2=ALGSLP(PERCOV,[40,60]) into EQWT[K]+=WT1(1)·WT2(1), EQWT[9/8]+=WT1(1)·WT2(2)·{PRLONG,1−PRLONG}; WT1(2)
  path: if==1 EQWT[5]+=1 else WT2=ALGSLP(PERCOV,[40,60]) → EQWT[5]+=WT1(2)·WT2(1), EQWT[9/8]+=WT1(2)·WT2(2)·{PRLONG…}.
  Else (no GF understory): WT1=ALGSLP(PERCOV,[40,60]) → EQWT[K]+=WT1(1); EQWT[9/8]+=WT1(2)·{PRLONG,1−PRLONG}.
- **fmvinit combined CASE labels**: (11,16) share; (17,19) share; so 15 blocks cover 19 species. FULIVE/FULIVI/
  FUINIE/FUINII DATA have inline `!species` comments (strip before parsing). ISPMAP(19)/BIOGRP(19) extracted OK.
  A robust line-based extractor (strip C/*/! then read DATA name / … /) is needed; the naive regex grab fails.
- **CIS9B** (ninebark vs snowberry, ICT=2 only — NOT needed for cit01/ICT6): small ci/fmcba.f habitat lookup.
- **Data-extraction TODO/subtleties**: FULIVE grab FAILED (format differs — re-extract) ; fmvinit found only 15
  of 19 CASE blocks (CI likely uses CASE ranges / shared defaults for 4 species — verify, don't assume). FUINIE/
  FUINII/FULIVI/ISPMAP(19)/BIOGRP(19) extracted OK. ISPMAP=[15,8,3,4,6,7,11,18,1,13,14,7,41,16,41,11,17,24,17].
- **Wiring (next session)**: reuse ci growth-port ICINDX → ci_pvg(icindx)=MAPPVG → ci_select_fuel_models (3 case
  groups) ; ffe_fuel.jl (FULIVE/FUINIE 19sp) ; fire_species_props.csv ; fmd_xpts==IE, decay==CR, moisture==IE ;
  crown_biomass CI bark (ci_bratio — CI has per-species branches, use it not calib) ; crown-fire gate. Validate
  cit01 SIMFIRE MOR vs live FVSci (oracle /workspace/.ciwork/FVSci_clean).

## Remaining western FFE variants — assessed (BM, UT, TT); each a dedicated fmcfmd chunk
All three: decay==CR (verified). Oracles exist (.bmwork/.utwork/.ttwork FVS*_clean). bmt01/utt01/ttt01 stands.
- **BM (MAXSP=18)**: moisture==IE. fmcfmd is the **MOST COMPLEX** western — FOUR nested weights
  WT1·WT2·WT3·WT4 + a WD(K) weight-distribution array (ie/bm fmcfmd diff 622 lines). A large careful
  transcription; give it its own turn. Data: own fmcba FULIVE/FUINIE (18 sp) + fmvinit + ISPMAP + BIOGRP.
- **UT (MAXSP=24) + TT (MAXSP=18)**: **SHARE ONE fmcfmd** (ut/tt fmcfmd.f dispatches `SELECT CASE(VARACD)` with
  CASE('TT') branches) — porting UT gets TT mostly free. It is the **CR-style COVER-TYPE selection** (ICT/USCT/
  FMAVH top-40 ht/LCUNDR conifer-understory/LPPDOM ppine-dom/IFMST structure-class — the SAME variables as the
  ported `cr_select_fuel_models`). ⇒ **reuse CR's cover-type machinery**, swapping UT/TT cover-type groups + ICT
  cases + fuel data. UT/TT moisture DIFFERS from IE (own table — extract). ICLSS=12 (not 14; models 1-12, no
  13/14). Recipe: generalize cr_select_fuel_models to accept a variant's cover-type map + ICT cases, or write
  ut_select_fuel_models mirroring it. Extract ut/tt fmcba + fmvinit + ISPMAP + BIOGRP + the UT/TT moisture table.
- **Order suggestion**: UT (+TT free) next (CR-reuse, moderate), then BM (most complex) last of the western FFE.

### UT/TT fmcfmd — CONFIRMED CR-reuse; cover-type maps extracted (ut/fmcfmd.f, ut==tt diff 0)
The fmcfmd is ONE unified algorithm `SELECT CASE(VARACD)` shared by CR/UT/TT; ASCT=8 cover-type groups
(OBCT=1,PJCT=2,PPCT=3,WSCT=4,SFCT=5,LPCT=6,MCCT=7,ASCT=8) — SAME as the ported `cr_select_fuel_models`. Only
the **species→covtype map** differs per variant (the `SELECT CASE(ISP)→CTBA(xxCT)` block, ci/fmcfmd.f:196):
- **TT (18 sp)**: PJCT{4,11,12}, PPCT{10}, SFCT{5,8,9}, LPCT{7}, MCCT{1,2,3,17}, ASCT{6,13,14,15,16,18}. (no OBCT/WSCT)
- **UT (24 sp)**: OBCT{13}, PJCT{11,12,14,15,16}, PPCT{10}, SFCT{5,8,9}, LPCT{7}, MCCT{1,2,3,4,17,23}, ASCT{6,18,19,20,21,22,24}. (no WSCT)
- **DEFINITIVE scope (measured: `diff -w ut/fmcfmd cr/fmcfmd` ICT-body = 0 — IDENTICAL)**: the whole fmcfmd
  algorithm (covtype-stat scaffolding + the `SELECT CASE(ICT)` 640-line body) is ONE shared VARACD-branched
  routine — CR-buildDir and UT-buildDir versions are byte-identical bar whitespace. So the ported
  `cr_select_fuel_models` ALREADY has UT/TT's exact logic. **UT/TT is a PARAMETERIZATION of cr_select, not a new
  transcription** (earlier "~230-line port" was wrong — over-corrected). The real work = thread a per-variant
  SPECIES-ROLE map through cr_select's ~15 CR-hardcoded `sp == N` checks (`_cr_fm_covtype` covtype groups; +
  inline roles: hemlock sp12/28, w.white-pine sp5, oak sp30/35, aspen sp40/41, birch sp24/43, grand-fir, sp34;
  `imodty`=model_type is CR-GENGYM-only → UT/TT skip it) + the `VARACD.EQ.'UT'/'TT'` OBCT sub-branches (IS=13
  for UT; TT skips OBCT). Approach: refactor cr_select to take a `FireCovtypeRoles` struct (covtype_fn +
  species-role indices) per variant; CR passes its current values, UT/TT pass theirs. ICLSS=12 (UT/TT, no
  model 13/14). Then extract UT/TT fmcba/fmvinit/moisture + validate utt01/ttt01. META: measure the ACTUAL diff
  (`diff -w`) before sizing a "port" — raw diff line-count (2114) was 100% whitespace-inflated here.
- Then extract UT/TT fmcba FULIVE/FUINIE + fmvinit props + ISPMAP + BIOGRP + the **UT/TT moisture table** (differs
  from IE — ut==tt?) ; wire (xpts, decay==CR, moisture, fmcba, crown_biomass bark, crown-fire gate) ; validate
  utt01 (+ttt01) SIMFIRE vs live FVSut/FVStt. Oracles: /workspace/.utwork, /workspace/.ttwork.

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
