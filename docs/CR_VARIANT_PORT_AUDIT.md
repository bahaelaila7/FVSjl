# CR (Central Rockies) variant port — audit log

Charter/doctrine: docs/CR_VARIANT_PORT_GOAL.md. Branch: cr-variant-port. Oracle: tmp/oracles/FVScr_new
(relinked from bin/FVScr_buildDir/*.o, 12.7 MB; runs — prints the FVS banner). Test keyfiles: tests/FVScr/*.key.

## Chunk 0 — scaffold [DONE]
- `CentralRockies <: AbstractVariant` singleton + registration (src/variants/centralrockies/centralrockies.jl):
  variant_code="CR", nspecies=38, htg_period=10 (YR=10, like NE/CS/LS), RNG seed 55329 (shared with eastern).
- Registered: include in src/FVSjl.jl; variant_from_code("CR") → CentralRockies() (variant.jl).
- Verified: `using FVSjl` loads; CentralRockies() resolves (code=CR, nsp=38). Unimplemented hooks error loudly.
- Oracle relinked (tmp/oracles/FVScr_new) + confirmed runnable.
- data/centralrockies/ created (coefficients land per chunk).

## Chunk 1 — infra/species [DONE]
- data/centralrockies/species_coefficients.csv: 38 species (index, code_alpha, code_fia, code_plants) from
  cr/blkdat.f (comment table + PLNJSP). Verified bit-exact vs Fortran: AF/CB/DF/GF/WF (FIA 19/18/202/17/15,
  PLANTS ABLA/ABLAA/PSME/ABGR/ABCO) ... OS/OH (299/998).
- data/centralrockies/species_translation.csv: 38-row self-map PLACEHOLDER (real cross-variant FIA crosswalk +
  a target_cr column is a later chunk; the eastern translation table hardcodes 4 targets).
- coefficients(::CentralRockies) = cached_coefficients(load_species_coefficients(CR_DATADIR), "CR"). Loads;
  code arrays padded to MAXSP=108, CR fills the first 38.
- init_blockdata!(::CentralRockies) (species.jl): copies codes, year=10, zeide_sdi=true (cr/grinit.f:134
  LZEIDE=.TRUE.), RNG both streams 55329. load_species_coefficients!(::CentralRockies)=init_blockdata!.
- Includes placed AFTER core/state.jl (species.jl signature needs StandState). `using FVSjl` loads clean.

## Chunk 2 — site/habitat [PARTIAL: site-index bounds done]
- site_lo/site_hi columns added to species_coefficients.csv (cr/blkdat.f SITELO/SITEHI, 38 species; verified
  bit-exact: 40/30/40… lo, 105/100/120… hi).
- src/variants/centralrockies/site_index.jl: cr_relative_si (cr/regent.f:209-213 RELSI) — verified
  cr_relative_si(DF,90)=0.625=(90-40)/(120-40).
- cr_site_index_defaults! (site_index.jl): the sitset ISISP between-species SI conversion (cr/sitset.f:479-497)
  — TEM default by IMODTY (70/70/57/75/65), ISISP default (DF/PP/PP/ES/LP), per-species scale by SITELO/SITEHI.
  VERIFIED bit-exact: SITEAR[AF]=64.375, SITEAR[DF]=70.0 (imodty=1). This supplies dgf's SSITE=SITEAR(ISPC).
- STILL TODO in this chunk: IMODTY determination (MODTYPE kw / habitat default) + SDIDEF SDI-max defaults +
  cr/habtyp.f habitat-type groups (OCURHT (16,MAXSP)). CR-specific; feed DG/mortality.

## Road ahead (the substantial growth-model chunks)
3. **Large-tree diameter growth = GENGYM (cr/gemdg.f), NOT Wykoff DDS** [SCOPED, correction logged].
   dgf.f sets stand/tree terms then CALL GEMDG(DDS,IMODTY,IS,BAUTBA,SPBA,SI,DP,BAT,BARK,CR,SLOPE,ASPECT,ELEV,
   PBAL,PCCFI,RELDEN,BAL); GEMDG dispatches IMODTY(1-5 model types) x SELECT CASE(species) => per-species
   regression DDS (hardcoded coefficients). WK2=DDS+COR+DGCON. The gem* trio (gemdg/gemht/gemcr) = the CR model.
   This is the largest chunk: faithful transcription of every IMODTY x species branch + calibration (COR2), then
   per-tree DG diff vs live FVScr in the pre-tripling window.
   [PARTIAL DONE] cr_gemdg (diameter_growth.jl): full GENGYM transcription — all ~15 species groups x IMODTY
   branches + the two output paths (DF-forecast / direct ln-DDS) + common tail (DSTAG inert by default, BH-PP
   x0.80, ln-DDS floor -9.21). Math: flog/fexp/fpow + Base sqrt/sin/cos. DBHMAX (cr/sitset.f) -> dbh_max column.
   ISTAGF all-0 default (DSTAG inert). Compiles; DF-path verified vs hand-computed DDS. ★ VALIDATED vs LIVE FVScr (doctrine #1): instrumented
   cr/gemdg.f to dump (IMODTY,IS,inputs,DDS) on crt01.key -> 9982 real per-tree tuples; cr_gemdg reproduces
   **9981/9982 bit-exact (99.99%)**, the 1 residual = 1 ULP (limber pine is=10, rel 5e-7, sin/cos transcendental
   rounding — the accepted ULP floor, not a model error). The GENGYM DDS transcription is FAITHFUL.
   STILL TODO: the dgf.f WRAPPER + wire diameter_growth!(::CentralRockies) + per-tree diff.

### dgf.f wrapper — per-tree term mapping (cr/dgf.f:166-197), for the next step
Stand-level (once/cycle): BA (stand BA); DGQMD/STDSDI/RELSDI => RELDEN; DSTAG = 3.33333*(1-RELSDI) if
RELSDI>0.7 (else 1). IMODTY (model type, from habitat/MODTYPE kw — chunk-2 remainder). SLOPE/ASPECT/ELEV.
Per species ISPC: SSITE=SITEAR(ISPC) (site index per species, sitset — chunk-2 remainder); SPBA=TBA(ISPC)
(species BA). Per tree I: D=DIAM; BARK=BRATIO(ISPC,D,HT); CR=ICR*0.01; ICLS=min(int(D+1),41);
BAUTBA=BAU(ICLS)/BA (BA-above-by-dbh-class / stand BA); IPCCF=ITRE(I) (point idx); PBAL=(1-PCT/100)*PTBAA(IPCCF)
(point BAL); BAL=(1-PCT/100)*BA (stand BAL); PCCFI=PCCF(IPCCF) (point CCF). Then
CALL GEMDG(...) => WK2(I)=DDS+COR(ISPC)+DGCON(ISPC). Most terms (BA, BAL, point-BA/CCF, PCT percentile) are
SHARED engine state the eastern DG already builds; CR-specific are BAU/BAUTBA (BA-above-by-class), RELDEN,
IMODTY, SITEAR. COR/DGCON are the DGSCOR calibration + per-species constant.

### dgf wrapper — engine dependency map (checked against src/engine + core/state)
SHARED (reuse, already built): BA = p.basal_area; PTBAA = density.point_ba; PTBALT = density.point_bal
(point_basal_area!); PCT (BA percentile) as the eastern DG uses (bal=(1-PCT/100)*BA). CR-SPECIFIC to compute
in a CR stand-stats pre-pass: BAU (BA-above-by-dbh-class, ICLS=min(int(D+1),41)) -> BAUTBA=BAU(ICLS)/BA;
RELDEN = STAND CCF (PLOT.F77: "current crown competition factor for the stand" — NOT relative SDI; measure>infer)
=> reuse the shared stand_ccf. RELSDI/DGQMD/STDSDI -> DSTAG=3.33333*(1-RELSDI) if RELSDI>0.7 (DSTAG INERT since
ISTAGF all-0). SPBA=TBA(ISPC)
(species BA sum); PCCF (point CCF, density.point_ccf); RELDEN=stand_ccf (shared). IMODTY = DEFMT(IFOR) (cr/sitset.f:73,94 — 23-forest table [5,3,4,5,3,4,5,5,4,4,5,4,2x11] -> data/
centralrockies/model_type_by_forest.csv; IFOR from forkod; MODTYPE keyword overrides). SITEAR done. =>
ALL GEMDG INPUTS NOW SCOPED. Wrapper = CR stand-stats pre-pass (BAU-by-class only new; rest shared) +
IMODTY=DEFMT[forest] + per-tree loop calling cr_gemdg, then WK2=DDS+COR+DGCON. Ready to implement + diff vs
live FVScr — the next chunk-3 unit is now turnkey (all inputs have a known source).
4-6. Height (htgf), crown (crown/cratet/ccfcal), small-tree (regent). 7. Mortality reconcile. 8. Volume (NVEL).
9. Full-cycle differential vs live FVScr. Best tackled as focused sessions, each port-and-diff per doctrine #1.

### dgf!(::CentralRockies) — COMPLETE implementation spec (turnkey; all data sources resolved)
Engine contract (like ls dgf!): loop trees, write `wk2 = view(s.scratch.wk,2,:)`; the shared
diameter_growth!(::AbstractVariant) driver handles calibration(COR)+tripling. Per-cycle setup:
  ba_v = p.basal_area; relden = stand_ccf(s) (RELDEN=stand CCF); slope=p.slope, aspect=p.aspect;
  imodty = s.control.model_type (= DEFMT[forest] default, MODTYPE kw override — data/centralrockies/
    model_type_by_forest.csv); spba[sp] = species BA sum (0.005454154*sum(tpa*d^2) per species);
  BAU pre-pass: for ICLS=1..41, BAU(ICLS)=BA in trees with class > ICLS (ICLS=min(trunc(d+1),41)),
    i.e. cumulative BA above each dbh class.
Per tree i (d=t.dbh[i]>0, sp=species):
  bark=bark_ratio(c.bark_a,c.bark_b,sp,d); cr = t.crown_pct[i]*0.01 (ICR*0.01, 0-1);
  icls=min(trunc(Int,d+1f0),41); bautba = BAU[icls]/ba_v; pct = t.crown_ratio[i] (BA percentile, shared);
  pbal = (1-pct/100)*point_ba[pt]; bal = (1-pct/100)*ba_v; pccfi = point_ccf[pt]; ssite = p.sp_site_index[sp];
  dds = cr_gemdg(imodty, sp, bautba, spba[sp], ssite, d, ba_v, bark, cr, slope, aspect, pbal, pccfi, relden,
                 bal; dbhmax=sd[:dbh_max][sp]);
  wk2[i] = dds + c.dg_cor[sp] + c.dg_const[sp]   # DGCON=0 default (cr/dgf.f:222); COR=shared calibration.
DATA still to add: dg_const column = 0 (all species); confirm engine has point_ccf + control.model_type +
p.slope/p.aspect accessors (verify; add if missing). VALIDATION: cr_gemdg already 9981/9982 bit-exact vs live;
validate the wrapper's TERM computation by running the engine's density pre-pass on crt01 and diffing the
GEMDG inputs (bautba/spba/pbal/bal/relden/pccfi) vs the captured live GEMTRC, then wk2 vs live WK2.

### cr_bratio + dgf!(::CentralRockies) — WRITTEN; both cores LIVE-VALIDATED
- **cr_bratio** (cr/bratio.f): 3 equation forms dispatched by IMAP (data/centralrockies bark1/bark2/bark_imap
  columns, extracted from cr/sitset.f); species {12,13,16,29..37} zero coefs under IMODTY 3/4/5 (=> default
  IEQN-1 formula 0.9002-0.3089/TEMD, TEMD capped 19); H arg dead (RDANUW=H); clamp [0.80,0.99].
  ★ VALIDATED 9982/9982 BIT-EXACT vs the live BARK column in /tmp/cr_gemtrc.txt (species 2,5,10,13,15,18,20 —
  covers IEQN 1/2/3 incl. D-dependent). crt01 is IMODTY=2 (the 3/4/5 override path is faithful transcription,
  untested by crt01).
- **dgf!(::CentralRockies)** WRITTEN (diameter_growth.jl): BADIST pre-pass (BAU=BA-above-class via the
  TOTBA-BAU(1) then BAU(i)=BAU(i-1)-BAU(i) recurrence; TBA=species BA; coef 0.0054542; HT<4.5 seedlings
  excluded) + AGERNG (age range over ABIRTH>1 & HT>4.5, else 1000 — inventory trees have ABIRTH=0 => 1000,
  only aspen suckers set it via esuckr.f) + per-tree loop: bark=cr_bratio, cr=ICR*0.01, ICLS=min(int(D+1),41),
  BAUTBA=BAU[icls]/BA, SPBA=TBA[sp], PBAL=(1-PCT/100)*point_ba[pt], BAL=(1-PCT/100)*BA, PCCFI=point_ccf[pt],
  SSITE=sp_site_index[sp], RELDEN=stand_ccf(s); wk2[i]=cr_gemdg(...)+dg_cor[sp]+dg_const[sp] (DGCON=0).
  Module LOADS clean. Two cores (cr_gemdg 9981/9982, cr_bratio 9982/9982) live-validated.
- NEXT: validate the wrapper's STAND-STAT pre-pass (BAU/TBA/PBAL/BAL/RELDEN/PCCFI + WK2 assembly) — needs
  either the CR engine stand-load integration (tree input + RCON + density pre-pass) OR a BADIST/dgf.f
  instrumentation dump on crt01 to diff the intermediate terms. dg_const CSV column (=0 all species) still TODO.
