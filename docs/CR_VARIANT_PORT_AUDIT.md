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
   ISTAGF all-0 default (DSTAG inert). Compiles; DF-path (DF spp) verified bit-exact vs hand-computed DDS.
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
4-6. Height (htgf), crown (crown/cratet/ccfcal), small-tree (regent). 7. Mortality reconcile. 8. Volume (NVEL).
9. Full-cycle differential vs live FVScr. Best tackled as focused sessions, each port-and-diff per doctrine #1.
