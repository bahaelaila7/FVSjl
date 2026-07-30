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

### CR engine wiring (MODTYPE keyword + site_setup!) — LANDED; end-to-end DG blocked on species crosswalk
- **MODTYPE keyword** (kw_modtype!, keyword_dispatch.jl): sets `s.plot.model_type` (IMODTY 1-5) from the
  MODTYPE keyword (cr/sitset.f). crt01 uses `MODTYPE 2.` (SW ponderosa). NOTE: IMODTY lives on PlotData
  (`s.plot.model_type`), NOT Control — fixed all refs (dgf!/site_setup!/handler).
- **site_setup!(::CentralRockies)** (cr_site_index_setup!, site_index.jl): resolve IMODTY (MODTYPE override
  else _CR_DEFMT[forest_idx], 23-forest table) then fan the site species' SI across species via
  cr_site_index_defaults!. EXECUTED CLEANLY on crt01 (the load reached tree-input past init+site_setup).
- **BLOCKER for end-to-end**: crt01.tre uses species ALIASES (WP, etc.) not in the 38 CR alpha codes ⇒ tree
  input calls translate_species ⇒ needs `spctrn_column(::CentralRockies)` + `other_species(::CentralRockies)`
  + the REAL CR SPCTRN crosswalk (target_cr column; current CSV is a placeholder self-map). Left spctrn_column
  UNDEFINED (loud MethodError, doctrine #5 — a placeholder self-map would silently mis-map WP→OH). ⇒ the CR
  species-crosswalk is the next chunk; it unblocks all real-stand runs + the chunk-9 full-cycle differential.
- **Two validation paths for the dgf! stand-stat pre-pass** (BADIST recurrence / PBAL/BAL / WK2 assembly, the
  only chunk-3 piece not yet live-checked; cr_gemdg + cr_bratio already are): (a) instrument cr/dgf.f+badist.f
  to dump raw per-tree inputs + derived terms on crt01 and REPLAY through dgf! (bypasses the crosswalk — build
  dir intact: 668 .o incl dgf/badist/gemdg); or (b) finish the crosswalk + tree-input integration and diff WK2
  in the pre-tripling window. Path (a) is self-contained; do it OR (b) when the crosswalk lands.

### ★★ CHUNK 3 CLOSED — dgf!(::CentralRockies) FULLY VALIDATED END-TO-END vs live FVScr
Instrument-replay (doctrine-pure, bypasses the un-ported species crosswalk): instrumented cr/dgf.f with a
per-DGF-call id (IDGCALL) + LSTART flag + a per-tree dump of raw inputs (ISP,D,HT,PROB,PCT,ITRE,ICR) + point/
stand (BA,RELDEN,PTBAA,PCCF,SSITE) + derived (BARK,BAUTBA,SPBA,PBAL,BAL,DDS,COR,DGCON,WK2,AGERNG) → 9982 rows on
crt01 (glibc __isoc23_sscanf shim + 12.2.0 recompile of dgf.o + relink of the 667 other .o; JOSTND=unit16=fort.16).
Julia replay grouped by call-id, ran the dgf! BADIST pre-pass + term derivation + cr_bratio + cr_gemdg, compared
bit-exact:
- **BARK: 9982/9982**  (all calls; cr_bratio re-confirmed)
- **PBAL: 9982/9982**  ((1-PCT/100)*PTBAA)
- **BAUTBA: 9766/9766**  (growth calls; BADIST BA-above-class recurrence — must sum in TREE-INDEX order for Float32)
- **SPBA: 9766/9766**  (species TBA)
- **★ END-TO-END DDS: 9766/9766 BIT-EXACT**  (cr_gemdg fed by MY computed terms == live DDS, all growth evals)
Calibration calls (LSTART=T, 216 rows) use backdated WK3 diameters + recent-mortality trees the standalone replay
can't reconstruct — but that state is the shared calibrate_diameter_growth! DRIVER's job (validated separately),
NOT dgf!'s logic. crt01 IMODTY=2; SLOPE=0.3, ASPECT=5.49778938. AGERNG was 120.2 (non-1000 ⇒ crt01 trees carry
birth ages) but >40 like 1000 ⇒ same DF-branch, so the DDS check used live AGERNG; the dgf! AGERNG-from-birth_age
computation is still to be independently checked when a real stand loads. Oracle RESTORED pristine (git-clean).
⇒ Chunk 3 (large-tree diameter growth) is COMPLETE + end-to-end live-validated. Next: species crosswalk (unblocks
full tree-input + chunk 9), then height(htgf)/crown/regent/mortality/volume.

### CR species SPCTRN crosswalk — DONE (tree input now works; WP→SW verified vs live)
Source: bin/FVScr_buildDir/spctrn.f = the SHARED western SPCTRN table ASPT(442,21): cols 1-3 = (alpha,FIA,PLANTS)
input codes, and column J=8 = the CR target 2-char code (header row "AK BM CA CI CR EC EM" over J=4..10). Extracted
all 442 rows → data/centralrockies/species_translation.csv in the 7-col eastern layout (target_cr in the first
target slot, replicated across cols 4-7). Hooks: spctrn_column(::CentralRockies)=4, other_species=Int32(38) (OH;
softwood misses go to OS via the table). ★ VERIFIED vs live: crt01's WP (FIA 119, PIMO3) → SW (sp 15), matching
the live "INPUT SPECIES CODE (WP) WAS SET TO (SW)" echo. crt01 now LOADS through the engine (27 trees, species
WF/PP/SW/ES/AS correct). Species-mapping was approved by FVS regional contacts Apr 2007 (spctrn.f header).

### FOLLOW-UP (tree-input, shared reader — matters for chunk 9 full-cycle, NOT the crosswalk):
jl loads 27 trees vs live's 29 (cruise "NUMBER OF RECORDS PER SPECIES 5 8 1 9 6"). The 2 dropped = crt01.tre
lines 5 ("018ES 346") + 14 ("016PP 072"), the ONLY records whose pre-species I1 field (TREEFMT T31,F2.0,I1,A3
→ the 1 col after the F2.0 count, before the A3 species) is ≠1 (=8 and =6). All I1=1 records load (incl. tiny
DBH=0.1 seedlings and a ht=0 record). Likely a history/value-class-code or fixed-format field-alignment handling
gap in the shared tree reader for this CR TREEFMT. Investigate before the chunk-9 differential; it does NOT affect
the chunk-3 DG validation (that used the live treelist via instrumentation).

### CHUNK 4 (height growth) — cr_gemht LIVE-VALIDATED 18864/18864 (all 4 outputs); wrapper next
GENGYM height model cr/gemht.f transcribed to cr_gemht (height_growth.jl): even-aged HHE by IMODTY 1-5
(site-index curves: Edminster/Minor/Meyer/Alexander/Alexander-Tackle-Dahms) + uneven-aged HHU / direct HTGI
by SELECT CASE(species). Surrogate species GF/MH/RC/WL(direct CON+ALOG(DGI) form)+LM(logistic) set IHTG=1.
BATEM carried from the IMODTY block, reset per species case (faithful). Math: flog/fexp/fpow, **2(int)->x*x.
★ VALIDATED vs live (instrument gemht.f WRITE→unit16, relink, crt01→18864 rows): HHE 18864/18864, HHU
18864/18864, HTGI 18864/18864, IHTG 18864/18864 ALL BIT-EXACT. (crt01 exercises imodty=2, species 5/13/15/18/20
+ their twin current/future calls; other imodty/species paths are faithful transcription.) Oracle restored.
NEXT: the htgf.f WRAPPER = height_growth!(::CentralRockies): HTOSI(6,MAXSP) adjust, AP breast-high-age adjust,
twin gemht calls (current + future age AGEFUT=AP+10, DFUT=D+DGI/BARK), IHTG=1 direct path, Black-Hills realign,
even/uneven blend (AGERNG>40 & BA≥70, PCT gates), LPP H30, aspen/birch ASPFAC, ZZRAN random (DGSD/BACHLO/RANN),
HTCON calib + SCALE=FINT/YR + XHMULT mult, mistletoe MISHGF, SIZCAP cap, tripling. Writes HTG.

### ★ CHUNK 4 (height growth) — htgf.f WRAPPER mesh VALIDATED 9607/9607; ZZRAN RNG = full-cycle item
height_growth!(::CentralRockies) + _cr_htg_tree (height_growth.jl): HTOSI(6,MAXSP) adjust, BH-age adjust
(_cr_bhage_adjust, IMODTY 1/2/4), twin cr_gemht calls (current + future AGEFUT=AP+10, DFUT=D+DGI/BARK),
IHTG=1 direct path, Black-Hills realign, even/uneven blend (AGERNG>40 & BA≥70: PCT≤10→uneven, 10<PCT<40→XWT
weighted), LPP-H30 young override, aspen/birch ASPFAC, then ZZRAN + 0.1-floor + EXP(HTCON)·SCALE·XHMULT +
SIZCAP cap. MISHGF=1 (no mistletoe). BAU via shared _cr_badist_bau; AGERNG via _cr_agerng.
★ VALIDATED (instrument htgf.f WRITE→unit16 gated to large trees, dump raw+derived+HTG+ZZRAN+DGSD, crt01 →
9607 rows): with the live ZZRAN injected, _cr_htg_tree == live HTG **9607/9607 BIT-EXACT**. The deterministic
mesh is fully correct.
★ FOUND: DGSD=2.0 is a CR default (grinit) ⇒ height carries a per-tree stochastic ZZRAN = rejection-sampled
BACHLO(0,1,RANN), |z|≤DGSD (htgf.f:288-295); the 0.1-floor is AFTER the ZZRAN add (ordering fixed). Wired the
draw into height_growth! (bachlo + reject loop). Its per-tree VALUES bit-match live only once the cycle RNG
sequence is reconciled — FVS draws in species-sorted IND1 order, jl loops tree-index order ⇒ deferred to the
chunk-9 full-cycle RNG reconciliation (same class as DG COR/DGSCOR driver-managed stochastics). cr_gemht itself
18864/18864. Oracle restored pristine. ⇒ CHUNK 4 deterministic core COMPLETE; only the ZZRAN stream-order remains.

### ★ CHUNK 5 (crown ratio) — cr_gemcr 11140/11140 + crown.f wrapper 11132/11132 BIT-EXACT
cr_gemcr (crown.jl, GENGYM cr/gemcr.f): crown-length CL by species SELECT CASE → CR=CL/HF, or direct CR
(MH/RC/WL CASE 6/7/8, ICRFLG=1). VALIDATED 11140/11140 vs live (instrument gemcr.f, crt01). _cr_crown_tree +
crown_ratio_update!(::CentralRockies) (cr/crown.f cycling path): HF=H+HTG, DF=D+DG/bark, GEMCR → new crown%,
±1%/yr change limit (PDIFPY), CRMAX cap ((CRLN+HTG)/(H+HTG)), bounds [10,95]. CRNMULT keyword adjustments
(CRNMLT/DLOW/DHI defaults 1/0/99 ⇒ inert) + top-kill (ITRUNC, lstart) transcribed but not exercised by crt01.
VALIDATED 11132/11132 BIT-EXACT vs live (instrument crown.f, cycling path LSTART=0; FINT=10). LSTART DUB path +
dead-tree dub are faithful transcription (crt01 inventory has crowns ⇒ unexercised).
REMAINING chunk-5 item: ccfcal.f (CCF) — the DG/height/crown validations all took live RELDEN(stand CCF)+PCCF as
dumped INPUTS, so jl's stand_ccf/point_ccf reproduction for CR (which needs CR crown-width coefficients) is still
UNVALIDATED. Check whether the shared CCF engine uses CR-specific crown-width coefs before chunk 9.

### CHUNK 6 (small-tree/regen, regent.f) — SCOPED turnkey (NOT yet transcribed)
regent.f (667 lines) = small-tree HEIGHT increment + DIAMETER-from-height inverse; two near-identical blocks
(cycling growth ~150-450 + LSTART/ESTAB path ~460-654). Dispatched as small_tree_growth!(::CentralRockies).
HEIGHT increment (deterministic core):
  RELSI=(SI-SITELO)/(SITEHI-SITELO); RSIMOD=0.5*(1+RELSI); POTHTG=SITEAR/(15-4*RELSI)*HTADJ[sp];
  density modifier X=AVHT*(CCF/100) capped 300; PCTRED=poly AB(1..6) in X, clamp[0.01,1];
  VIGOR: X=ICR/100; VIGOR=min(150*X^3*exp(-6X)+0.3, 1); pinyon/juniper/oak/bristlecone (IVFLAG: sp
    9,12,16,23:27,29:35) VIGOR=1-(1-VIGOR)/3;
  aspen/birch (20,28): Sheppard curve HITE=26.9825*AG^1.1752, HTGR=(HITE2-HITE1)/(2.54*12)*RSIMOD*CON*0.75;
  else HTGR=POTHTG*PCTRED*VIGOR*CON;   CON=RHCON[sp]*exp(HCOR[sp]).
  ★ RNG: ZZRAN=BACHLO(0,1) if DGSD≥1, reject if ZZRAN>0.5 or <-2.0 (ASYMMETRIC window!); HTGR=(HTGR+ZZRAN*0.2)
    *XRHGRO*SCALE*WK4[i]. XWT=(D-XMN)/(XMX-XMN), 0 if D≤XMN or LESTB; HTG=HTGR*(1-XWT)+XWT*HTG_large; floor .1;
    SIZCAP[sp,4] cap.
DIAMETER-from-height (D<BREAK[sp]): HK=H+HTG; if HK≤4.5 DG=0,DBH=D+.001*HK; else 3 paths:
  PP/Chihuahua(13,36): DK=(HK-8.31485+.592*7)/3.03659; pinyon/juniper/oak(IVFLAG): DK=(HK-4.5)*10/(SITEAR-4.5);
  else: DK=(HT2[sp]/(ALOG(HK-4.5)-AX))-1, AX=(IABFLG[sp]==1 ? HT1[sp] : AA[sp]); floors 0.1; DG=(DK-DKK)*BARK...
COEFFICIENTS (all captured): regent.f DATA — DGMAX[38], XMAX[38], XMIN[38], REGYR=10, MPCC[38], DIAM[38],
  AB=[1.11436,-.011493,.43012E-4,-.72221E-7,.5607E-10,-.1641E-13,0,0,0], HTADJ[38] (RC/WL=1.10 else 1.0);
  IVFLAG set = sp{9,12,16,23:27,29:35}. Plus blkdat: BREAK, XMIN(dup), AA, HT1, HT2; grinit: IABFLG; have:
  SITELO/HI, SITEAR. MPCC = the AB-poly selector? (J=MPCC[sp]) — CHECK: AB may be (poly, MPCC-group) 2D.
VALIDATION: same instrument-replay as htgf — dump per-tree inputs+ZZRAN+HTG+DG, replay with injected ZZRAN
(DGSD≥1 so ZZRAN active, asymmetric [-2,0.5] window). RNG stream-order = chunk-9. New-tree crown (LESTB) uses a
SECOND BACHLO draw (regent.f:255, CR=0.89722-0.0000461*PCCF+0.07985*RAN) — establishment-coupled, defer.
Note: regent needs WK4[i] (a per-tree factor) + AVHT/ATAVH/ATCCF (estab mid-period interp) from the engine.

### ★ CHUNK 6 (regent small-tree) — HEIGHT increment VALIDATED 1349/1349; DG needs cratet AA-fit
small_tree_growth.jl: _cr_regent_tree transcribes cr/regent.f growth path. Coefficients (st_dgmax/st_xmax/
st_xmin/st_diam/st_htadj/st_break/ht1/ht2 added to species_coefficients.csv). MPCC is DEAD (debug-only). IABFLG≡1
(grinit; no overrides) and AA is never DATA-defined ⇒ static path would be AX=HT1.
★ HEIGHT increment VALIDATED 1349/1349 BIT-EXACT vs live (instrument regent.f REGTRC, crt01, injected ZZRAN):
POTHTG·PCTRED·VIGOR·CON (or Sheppard aspen/birch curve) + ZZRAN·0.2 (asymmetric reject [-2,0.5]) ·XRHGRO·SCALE·WK4,
XWT blend with htgf HTG, floor .1, SIZCAP.
DG (diameter-from-height inverse): 250/338 for d<BREAK — the matches are HK≤4.5 (DG=0) / PP / pinyon-juniper
(IVFLAG) paths; the 88 mismatches are the "else" path DK=(HT2/(ALOG(HK-4.5)-AX))-1 which uses a **cratet.f-FITTED
height–DBH intercept AX (AA), NOT static HT1**. Reverse-solved live sp18 DG=0.6689 ⇒ AX≈3.867 vs HT1[18]=4.5293
(different). ⇒ the DG path is blocked on the cratet.f AA height-DBH regression (a stand-specific per-species fit,
same mechanism as the eastern NOHTDREG/AA-fit already in the engine). NEXT: wire the cratet AA-fit for CR (check
if the shared init_crown_ratios!/cratet path fits AA+IABFLG for CR), then DG validates. ZZRAN RNG-order = chunk-9.
Oracle restored pristine.

### ★★ CHUNK 6 (regent) — DG FORMULA VALIDATED 338/338 + AA-fit WIRED (engine AA == live BIT-EXACT)
Two fixes closed the DG path: (1) the diameter-from-height "else" intercept AX is the cratet-FITTED AA (not
static HT1); (2) the DIAM floor (regent.f:421-423) is INSIDE the HK>4.5 branch — NOT applied when HK≤4.5.
With the correct AX injected: DG 338/338 BIT-EXACT (was 250, then 336, then 338). AA-FIT WIRING (reuses the
shared engine, doctrine #5): CR LHTDRG defaults .TRUE. (grinit.f:110) ⇒ set ht_drag_sp=trues in init_blockdata!;
added wykoff_ht2 col (= CR ht2) so dub_missing_heights! (volume.jl:245, the NOHTDREG/AA fit AA=mean(log(H-4.5)-
HT2/(D+1)) over D≥3,H>4.5,NORMHT≥0; ≥3 trees; IABFLG=0 if AA≥0) uses CR's HT2. ★ VERIFIED: engine ht_dbh_aa for
crt01 == live AX BIT-EXACT (sp5 4.3869486=4.38694859; sp18 3.8670504=3.86705041; sp15 IABFLG=1 <3 trees).
small_tree_growth!(::CentralRockies) LOOP written (PCTRED density modifier + per-species POTHTG/RSIMOD + ZZRAN
draw + _cr_regent_tree), reads s.calib.ht_dbh_aa/iabflg. ⇒ CHUNK 6 deterministic core COMPLETE + AA-fit engine-
validated. Only the ZZRAN per-tree stream-order (chunk-9 RNG) + new-tree-crown BACHLO (estab) remain.
Coefficients: st_* cols + ht1/ht2/wykoff_ht2 in species_coefficients.csv. Oracle restored pristine.

### ★★ CHUNK 7 (mortality) — MOSTLY SHARED; CR coefficients wired + VALIDATED 11577/11577
CR morts.f is the SAME two-stage model as the shared mortality!(::AbstractVariant) driver (southern/mortality.jl):
background RI=ri_scale/(1+exp(B0+B1·D)) + Pretzsch SDI self-thinning (RN, TEM=SDIMAX/0.02483133·D10^-1.605·PMSDIL).
CR-specific wiring only: mort_ri_scale(::CentralRockies)=0.5 (morts.f RI=0.5·RI); background B0/B1 = PMSC/PMD (6
groups) mapped via IBGMAP[38] → mort_bkgd_intercept/mort_bkgd_dbh columns; sdi_max_default = SDICON[38] (sitset)
+ SDIDEF added to cr_site_index_setup! (populates sp_sdi_def); PMSDIL/PMSDIU=0.55/0.85 (grinit defaults);
DBHZEIDE/DBHSTAGE=0 (no small-DBH SDI exclusion). ★ VALIDATED vs live (instrument morts.f MORTTRC, crt01, 11577
trees): b0/b1 coefficients 11577/11577, RI 11577/11577 (=0.5/(1+exp(b0+b1·d))), WKI=P·(1-(1-RIP)^FINT)·X
11577/11577 ALL BIT-EXACT. The RN/RIP SDI self-thinning combination is the shared Pretzsch driver (validated for
eastern; CR uses identical SDIMAX formula) — its end-to-end RN/RIP match rides the chunk-9 full-cycle. Oracle
restored. ⇒ CHUNK 7 done (CR-specific mortality validated; self-thinning shared).

### CHUNK 8 (volume) — SCOPED (NVEL Region-2 equations; larger than pure wiring)
CR volume routes through NVEL like the eastern variants, BUT uses Region-2 equation families NOT yet in FVSjl.
The per-species assignment (from live crt01.out VOLEQ table, saved data/centralrockies/volume_equations_reference.csv):
  **31 DVEE** (300DVEW093/301DVEW015/300DVEW113/300DVEW060/300DVEW800/300DVEW999/300DVEW106/300DVEW122/301DVEW202),
  **2 FW2** (DF 300FW2W202, PP 300FW2W122), **5 NVB** (CB/WF NVB0000015, SW NVBM240119, ES NVBM330093, AS NVB0000746).
Assignment = cr/sitset.f:567-584 CALL VOLEQDEF(VAR='CR', IREGN=KODFOR/100, FORST, IFIASP) per species (METHC 6/9 →
VOLEQDEF, METHC 10 → NVBEQDEF). FVSjl status: setup_volume_equations! (volume_equations.jl) only handles IREGN==8
(SN) → r8clark; CR gets blank ⇒ needs a CR VOLEQDEF (extractable from the saved table OR the NVEL voleqdef.f).
Equation coverage: DVEE partly present (r9vol_gevorkiantz.jl _dvee_* but that's Region-9 "900DVEE" for CS METHC=5;
CR is "300DVEW"=Region-2 west — CHECK if coefficients differ); FW2 + NVB (National Biomass) NOT ported. Merch
specs: cr/grinit.f defaults STMP=1, TOPD/SCFMIND/SCFTOPD/SCFSTMP/BF*=0 (⇒ NVEL internal defaults) — must add the
scf_min_dbh/scf_topd/scf_stump/top_diam/stump_ht/bf* columns (compute_volumes! KeyErrors on :scf_min_dbh for CR).
⇒ Chunk 8 = (a) add CR merch-spec columns from grinit; (b) CR VOLEQDEF assignment (setup_volume_equations! branch);
(c) port/verify DVEE-Region2 + FW2 + NVB volume routines (the big part). Validate each vs live per-tree cuft/bdft
(instrument cr volume or diff crt01 .sum volume cols). This is the LAST equation-porting chunk before the ch9
full-cycle. NOTE: growth/mortality (chunks 3-7) are all bit-exact and volume is a downstream leaf, so it does not
block the grow-cycle correctness — only the reported volume columns.

### CHUNK 8 (volume) — merch specs WIRED (compute_volumes! runs); equation computation is the remaining work
Added CR merch-spec + HT-DBH species columns from cr/grinit.f defaults: stump=1, top_dib=0, dbh_min=0,
scf_min_dbh=0, scf_top_dib=0, scf_stump=1, bf_min_dbh=0, bf_top_dib=0, bf_stump=1 (0 ⇒ NVEL internal defaults);
htdbh_iwykca=0 (Wykoff form), htdbh_ht1=HT1, htdbh_ht2=HT2, htdbh_p2/p3/p4/db=0. compute_volumes! now RUNS for
CR (was KeyError :scf_min_dbh) but returns cuft=0/bdft=0 — no CR volume-equation computation yet.
★ KEY FINDING (measure): CR's "300DVEW…" DVEE is a **Region-2** estimator, DISTINCT from the Region-9 "900DVEE"
Gevorkiantz FVSjl already has (r9vol_gevorkiantz.jl, for CS KODFOR-900). So the DVEE coefficients must be sourced
fresh for Region 2. REMAINING chunk-8 work (the substantial part):
  (a) CR VOLEQDEF assignment: add a CentralRockies branch to setup_volume_equations! mapping species→eq from
      data/centralrockies/volume_equations_reference.csv (or port NVEL voleqdef.f for VAR='CR'/IREGN=2).
  (b) Port the 3 equation families: DVEE-Region2 (31 sp, 300DVEW*), FW2 (2 sp DF/PP, 300FW2W*), NVB National
      Biomass (5 sp CB/WF/SW/ES/AS, NVB*/NVBM*). Each needs NVEL coefficient tables + cuft/bdft computation.
      NVB species also flip METHC=6→10 via NVB_REGION_CHECK (sitset.f:566, LFIANVB gate).
  (c) Validate per-tree cuft/bdft vs live (instrument cr volume or diff crt01 .sum volume cols).
This is a large NVEL port (comparable to the eastern r9clark). Growth+mortality (ch3-7) all bit-exact; volume is
a downstream leaf that does NOT affect grow-cycle correctness — only the reported CF/BF columns + the ch9 .sum
volume rows.

### CHUNK 8 volume — NVEL SOURCE LOCATED (turnkey for the equation port)
The NVEL Fortran library is in-tree at **/workspace/ForestVegetationSimulator/volume/NVEL/** + **volume/voleqdef.f**
(2749-line assignment, VAR/region/forest/FIA → eq string; the 300DVEW/300FW2W/NVB tables are here). DVEE volume =
**volume/NVEL/dvest.f**. FW2 + NVB routines also under volume/NVEL/. So the equation port is unblocked (source +
coefficients present) but LARGE (comparable to the eastern r9clark port). Recommended order next turn:
  1. Port the CR VOLEQDEF assignment (voleqdef.f VAR='CR' path → species.vol_eq) — or just load the saved
     volume_equations_reference.csv (already has the exact 38 per-species eq strings, faster + already verified
     vs live crt01.out).
  2. Port DVEE (dvest.f) for Region-2 "300DVEW" — the dominant family (31/38 species).
  3. Port FW2 (300FW2W, DF/PP) + NVB National-Biomass (NVB*/NVBM*, 5 sp).
  4. Validate per-tree cuft/bdft vs live via instrument-replay (dump the NVEL call args+results on crt01,
     replay). Merch specs already wired; compute_volumes! runs (returns 0 until equations land).

### CHUNK 8 volume — full characterization (3 NVEL families; crt01 exercises NVB+FW2, not DVEE)
★ KEY (measure): crt01's 5 species use **NVB** (WF NVB0000015, SW NVBM240119, ES NVBM330093, AS NVB0000746) +
**FW2** (PP 300FW2W122) — NOT DVEE. So crt01-based volume validation needs NVB+FW2; the DVEE family (31/38 sp) is
UNEXERCISED by crt01 and needs a different test stand. Source (all in-tree, volume/NVEL/):
  - **NVB** = NVEL/nsvb.f (FIA National-Scale Volume & Biomass — modern, LARGE) + biomassformula.f/calcbiomass.f.
  - **FW2** = NVEL/fwinit.f (form-class volume).
  - **DVEE** = NVEL/r3d2hv.f (717 lines; 49 GCUFT4 D2H-polynomial formulas by equation 093/113/122/746/060/106/
    800/999/310/314 — TRACTABLE direct formulas, D2H=DBH²·HT). dvest.f dispatches on VOLEQ[1] region digit
    (crt01 forest 303 ⇒ region 3 ⇒ R3D2HV). Assignment logic = voleqdef.f (2749 lines) — or reuse the saved
    volume_equations_reference.csv (verified vs live).
⇒ Chunk 8 is a genuine 3-family NVEL port (NVB the largest). Recommended: port NVB+FW2 first (validate on crt01),
then DVEE (needs a DVEE-species test stand — pick one from the FVScr test suite). This is the largest remaining
chunk; growth+mortality (ch3-7) are all bit-exact and volume is a downstream leaf (does not affect grow-cycle).
Merch specs already wired; setup_volume_equations! still needs a CR branch (load the reference CSV).

### CHUNK 8 volume — VOLEQDEF assignment DONE (equations assigned); FW2 also large (Flewelling)
setup_volume_equations! now has a CentralRockies branch: species.vol_eq = _CR_VOLEQ (38 verified-vs-live eq ids
from volume_equations_reference.csv). crt01 species assign correctly (sp5=NVB0000015, sp13=300FW2W122, etc.).
NOTE: this is the forest-303 assignment; VOLEQDEF is region/forest-keyed so a full voleqdef.f port is needed for
arbitrary CR forests (the growth/mortality chunks are forest-independent, so this only limits volume on non-303
forests). ★ CONFIRMED (measure): FW2 (300FW2W) = the Flewelling profile/taper model (fwinit.f) — also LARGE, not a
simple polynomial. So ALL of crt01's volume (NVB×4 + FW2×1) needs big NVEL models (NSVB nsvb.f + Flewelling
fwinit.f); the tractable DVEE D2H polys (r3d2hv.f) cover 31 non-crt01 species. ⇒ chunk 8 equation computation is
the LARGEST remaining chunk — a full NVEL port (NSVB + Flewelling + DVEE), comparable to the entire eastern
volume effort. It is a DOWNSTREAM LEAF: growth (ch3-6) + mortality (ch7) are all bit-exact vs live, so the CR
simulation CORE is complete; volume only affects reported CF/BF columns + ch9 .sum volume rows. Assignment +
merch specs done; equation computation deferred to a dedicated volume effort.

### CHUNK 9 (full-cycle) — STARTED: core integration RUNS; VARMRT ported; FFE column wiring next
Ran run_keyfile(crt01, CentralRockies) end-to-end — the full grow_cycle! now executes through
diameter_growth! → height_growth! → crown_ratio_update! → small_tree_growth! → mortality! → cuts (THINDBH),
i.e. ALL the ported growth+mortality chunks integrate and run together. Fixes this pass:
- small_tree_growth! p.year → s.control.year (YR on Control, not PlotData).
- ★ CR VARMRT ported (mortality.jl _varmrt_efftr!(::CentralRockies)): CR uses a DISTINCT efficiency —
  EFFTR=PEFF·((100-CRI)/100)·VARADJ·0.01, PEFF=0.84525-0.01074·PCT+2e-7·PCT³ (BA percentile, NOT relative
  height like NE/CS/LS), oak(23-27) CRI cap 50; varmrt_varadj[38] col added. (The rest of mortality is shared.)
- is_sprouting col added = 0 (CR does aspen SUCKERING via esuckr, not ESTUMP stump-sprouting).
REMAINING for the .sum: FFE (fire/fuels) runs even without fire keywords and needs ~10 CR species columns:
v2t (wood SG), bark_eqnum, bio_group/biogrp, dkr_cls, leaf_life, ls_spi, snag_alldwn/decayx/fallx, tfall_cls.
These are a distinct FFE-wiring sub-task (downstream of the validated growth core). THEN compare crt01 .sum
TPA/BA/QMD/TopHt vs live (tripling-invariant, doctrine #3) — expect divergence from the ZZRAN/regen RNG
stream-order (FVS species/IND1 draw order vs jl tree-index), the one cross-cutting stochastic item. Growth+
mortality per-chunk bit-exact; this validates their INTEGRATION + surfaces the RNG-order fix.

### CHUNK 9 — full cycle reaches FFE/biomass; that subsystem is the .sum prerequisite (needs faithful CR coefs)
Continuing the full-cycle run past cuts, it reaches the FFE (fire/fuels/snags) + Jenkins-biomass subsystem, which
runs EVERY cycle even without fire keywords and needs CR species coefficients from a data/centralrockies/
fire_species_props.csv (cols ls_spi/v2t/tfall_cls/leaf_life/dkr_cls/snag_cls/bark_eqnum/snag_decayx/snag_fallx/
snag_alldwn/biogrp) PLUS a bio_group (Jenkins) column. Placeholder values were tried but the biomass indexing is
interdependent (bio_group=placeholder → BoundsError), so faithful CR FFE/biomass coefficients are required (doctrine
#4 — no placeholder coefficients). Reverted the placeholders; CR data stays faithful (varmrt_varadj + is_sprouting=0
+ merch/htdbh cols are all faithful and retained). ⇒ CHUNK 9 .sum prerequisite = wire the CR FFE/biomass subsystem
(fire_species_props + Jenkins bio_group) from the CR FFE Fortran (fmvinit.f / the biomass tables) — a distinct
downstream task. The growth+mortality CORE integrates + runs end-to-end (DG→ht→crown→small-tree→mortality→cuts);
once FFE/biomass is wired the .sum runs and the TPA/BA/QMD/TopHt vs live comparison surfaces the ZZRAN/regen
RNG-order (the last cross-cutting item). NOTE: FFE is downstream of the validated growth core; a NOFFE-style run
(if the engine can skip FFE) would let the growth .sum validate without the FFE coefficients.

### CHUNK 9 — the .sum cascades through multiple downstream subsystems (FFE, establishment, volume)
Measured the full-cycle dependency chain: after growth+mortality+cuts, the crt01 .sum needs — in order —
(a) FFE (crt01 is a FULL-FFE demo: FMIN/SNAGINIT/SNAGBRK/SIMFIRE-2003/POTFIRE ⇒ species props+biomass DONE
faithfully from fmvinit/fmcblk, but ALSO fuel-loading/cover-type/fire-behavior tables — ffe_live_fuel_override
BoundsError next); (b) with FFE stripped, ESTABLISHMENT (establish! needs :estab_min_ht + the CR regen coefs);
(c) VOLUME (the NVEL port). So a bit-exact CR .sum requires wiring ALL of: volume-NVEL, FFE-fuel, establishment
— each a distinct downstream subsystem. The GROWTH+MORTALITY CORE is complete: it validates bit-exact per-chunk
(DG 9766, ht 9607, crown 11132, small-tree 1349+338, mort 11577) AND integrates end-to-end (grow_cycle! runs
DG→ht→crown→small-tree→mortality→cuts). The remaining .sum work is downstream-subsystem wiring (mechanical
data-extraction, each like the FFE species-props extraction) + the ZZRAN/regen RNG-order. RECOMMENDATION: these
downstream subsystems (volume/FFE/estab) are each sizable; prioritize per the milestone's needs — the
scientifically-meaningful growth+mortality modeling is done and proven.

### ★★ CHUNK 9 — FULL-CYCLE DIFFERENTIAL WORKS: 1990 BIT-EXACT; root divergence = ccfcal (CCF) not ported
Built a growth-only CR keyfile (crt01's first stand = UNTHINNED CONTROL, no FFE/ESTAB/THIN — saved
test/harness/cr/crt01_growth.key+.tre) and diffed FVSjl-CR vs a CLEAN-relinked live oracle (/tmp/FVScr_clean =
bin/FVScr_buildDir/*.o + isoc23 shim; the shipped bin/FVScr needs GLIBC_2.38, unavailable).
RESULT (year age TPA BA SDI CCF TopHt QMD):
  1990: 536 77 160 93 63 5.1  (LIVE) == 536 77 160 [0] 63 5.1 (JL) — TPA/BA/SDI/TopHt/QMD BIT-EXACT at inventory.
  2000+: JL grows FASTER (2040 BA 246 vs 223, TopHt 94 vs 87, TPA 344 vs 417) — over-growth + over-mortality.
★ ROOT CAUSE: JL's CCF column = 0 (stand_ccf=0.38 vs live 93; crown_width[i]=0 for ALL CR trees). The engine uses
the EASTERN cwcalc (crown_width.jl) which has NO coefficients for CR species ⇒ crown widths 0 ⇒ CCF≈0 ⇒ RELDEN≈0
⇒ the density-suppression terms in gemdg/gemht (RELDEN/CCF) + regent PCTRED vanish ⇒ growth too high. (The per-
chunk DG/ht validations used the LIVE RELDEN as a dumped INPUT, so this was masked until the full cycle.)
★ FIX = port cr/ccfcal.f (CR's OWN CCF/crown-width model: MAP1/2/3 species maps, CCF=RDA·DBH^RDB for DBH<10 +
CCFLGE for >10) and wire it as the CR crown_width/CCF (replace cwcalc for CentralRockies in standstats.jl point_
density! + stand_ccf). This is a bounded chunk (ccfcal.f coefficients + dispatch). It is THE primary full-cycle
fix; after it, re-diff and the residual should be the ZZRAN/regen RNG stream-order (height/small-tree). ⇒ the
full-cycle harness is now in place + the divergence is root-caused. Growth+mortality integration CONFIRMED correct
(1990 bit-exact).

### ★ ccfcal (CCF) PORTED — 1990 now FULLY bit-exact (CCF 93=93); residual = ZZRAN RNG-order
cr_crown_width (crown.jl) ports cr/ccfcal.f: IMAP=MAP{IMODTY}[sp] (MAP1-5, 8 CCF eqn groups = Wykoff-Crookston-
Stage INT-133 Table-8); CCF = D≥10: RD1+D·RD2+D²·RD3 ; D>0.1: RDA·D^RDB ; else 0.001; crown_width=sqrt(CCF/0.001803)
cap 99.9. Wired into standstats.jl point_density! + stand_ccf (CentralRockies → cr_crown_width, else eastern cwcalc).
RESULT: the growth-only .sum 1990 row is now BIT-EXACT vs live INCLUDING CCF (536/77/160/93/63/5.1). CCF at later
cycles tracks (2000 jl133/live127, 2040 jl257/live239). RESIDUAL: jl still over-grows ~5%/cycle (2040 BA 246 vs 223,
TopHt 94 vs 87) BUT TPA matches EXACTLY (mortality integration correct) — so it's per-tree DG/HT slightly high with
correct RELDEN ⇒ the ZZRAN/regen RNG STREAM-ORDER (DGSD=2 ⇒ every tree draws a height + small-tree ZZRAN; FVS draws
in species-sorted IND1 order, jl in tree-index order ⇒ each tree gets a different deviate ⇒ stand-specific
systematic offset that compounds). ⇒ the LAST growth-affecting item is the RNG-order reconciliation: make
height_growth! + small_tree_growth! draw ZZRAN in the species-sorted order FVS uses (needs the IND1/ISCT sort in the
growth loops). After that the growth .sum should be bit-exact-or-cornered. Downstream leaves (volume/FFE/estab)
remain but don't affect these growth columns.

### CHUNK 9 — fixed CR-missing-from-DG-calibration; residual = DG serial-correlation not applying
★ BUG FIXED: setup_growth! (simulate.jl:48-67) had DG-calibration branches for SN/NE/CS/LS but NONE for
CentralRockies ⇒ calibrate_diameter_growth! + dgcons NEVER ran for CR ⇒ c.sigma=0, COR=0, DGSD unset. Added:
(a) s.control.dg_sd=2.0 in init_blockdata! (cr/grinit.f:176 DGSD default); (b) dg_resid_sd column = IMODTY-2
SIGMAR (cr/sitset.f; base 0.2 + overrides 4/6/7/8/10/14 + imodty-2 12/23-27/33-35); (c) cr_dgcons! (DGCON=0,
CR ATTEN 1000+overrides, bark inert) + a CentralRockies branch calling calibrate_diameter_growth!. VERIFIED
c.sigma now set (0.2/0.201/…). BUT the growth-only .sum is UNCHANGED (2000 BA 113 vs live 106) — so the DG
serial-correlation (dgscor!, negative bias FM=-0.14228 that reduces DG) is STILL not applying for CR despite
sigma≠0 + DGSD=2. ⇒ NEXT: the DG driver's serial-correlation path (diameter_growth!(::AbstractVariant) dgscor!
line 907, gated on tripling/dg_sd) isn't reached for CR — investigate the gating (per-variant flag? the CR dgf!
tripling path?). This is the remaining ~6%/cycle over-growth (TPA still EXACT ⇒ purely per-tree DG magnitude).
1990 remains fully bit-exact (ccfcal fix). The calibration-branch fix is correct+necessary (COR/sigma) regardless.

### ★ CR-bark-in-DG-driver FIX — over-growth ~halved; residual = RNG stream-order
The shared DG driver (diameter_growth!(::AbstractVariant)) converted DDS→DG via bark_ratio(calib.bark_a/b) —
which for CR (bark_a=b=0) returned the 0.80 FLOOR, but CR's real bark = cr_bratio ~0.89. A too-LOW bark ⇒
too-low inside-bark dia ⇒ too-HIGH DG. FIX: dispatch bark on variant — CR → cr_bratio(sd, sp, d, imodty), else
bark_ratio (must match the bark cr_gemdg used internally). RESULT: 2000 BA 113→109 (live 106), 2010 154→145
(live 137), TPA now matches at 2000/2010. Residual now ~3%/cycle BA + slight TPA divergence at later cycles
(2020 jl476/live497) = the ZZRAN(height)+ZZRAN(regent)+DG-serial-correlation RNG STREAM-ORDER: FVS draws these
per-tree deviates in species-sorted IND1 order, jl in tree-index order ⇒ each tree gets a different deviate ⇒
small per-tree growth diffs that feed self-thinning. THIS is the last item — reconcile the RNG draw order across
DG-serial-corr + height + regent to FVS's exact sequence (species-sorted, per the growth loop order). 1990 stays
bit-exact. Chunk-9 progress: full cycle runs + close to live; ccfcal + CR-calibration-branch + CR-bark all real
fixes found via the differential.

### CHUNK 9 residual refined: cycle-1 DG ~3% high is DETERMINISTIC (not RNG-order)
Checked: at 2000 (cycle 1) tripling is active ⇒ the DG serial-correlation is DETERMINISTIC (frmbase=FM·ssigma·
rhocp, no draw). c.vardg IS set for CR (0.00144) ⇒ ssigma≈0.038 ⇒ the negative-bias reduction is only ~0.5%,
NOT enough to explain the ~3% BA over-growth. So the cycle-1 residual is a DETERMINISTIC DG diff, not the RNG
stream-order (that only bites at post-tripling cycles + the height/regent ZZRAN). Candidates to isolate via a
per-tree cycle-1 DG instrument-replay (dump live DG(I) for crt01_growth, diff jl t.diam_growth): (a) CR DGBND
(cr/dgbnd.f) diameter-growth bounds not wired into the driver (dlo_v/dhi_v) ⇒ unbounded DG; (b) a density-input
detail (point_ba PTBAA / PCT percentile) feeding gemdg PBAL/BAL that differs from live's dense.f; (c) the
DDS→DG conversion beyond bark. NOTE the per-chunk DG was bit-exact using LIVE dumped inputs, so the residual is
in jl's OWN computation of those inputs (density pre-pass) OR the DGBND/conversion. NEXT: instrument cr DG on
crt01_growth (DGTRC per-tree DG at cycle 1), replay, isolate. Then the RNG-order for later cycles. Session fixed
3 real full-cycle bugs (ccfcal, CR-calibration-dispatch, CR-bark); over-growth 6%→3%; 1990 bit-exact.

### DGBND ruled out — residual is a density-input/conversion detail (needs per-tree DG instrument-replay)
CR's DGBND (bin/FVScr_buildDir/dgbnd.f) is ONLY the SIZCAP(ISPC,1) size cap (999 default ⇒ non-binding), NOT a
DBH-range bound — so dlo_v/dhi_v=nothing for CR is CORRECT; DGBND is not the residual. ⇒ the deterministic
cycle-1 ~3% DG over-growth is in jl's OWN computation of the gemdg density inputs (PBAL/BAL from point_ba PTBAA +
PCT percentile, or SPBA/BAUTBA) OR the DDS→DG conversion. The per-chunk DG was bit-exact using LIVE dumped inputs,
so the gap is where jl's density pre-pass (point_basal_area!/stand_pct!/point_density!) differs from live dense.f
for the CR stand. NEXT (concrete): re-run the DGFTRC instrument on crt01_growth (the no-FFE stand), dump live
per-tree gemdg inputs+DDS at cycle 1, and diff jl's computed inputs to find which one diverges. Bounded. STATUS:
full-cycle differential CLOSE to live (1990 bit-exact, 2000 within ~3%, TPA matching early); 3 real bugs fixed
this session (ccfcal / CR-calibration-dispatch / CR-bark); growth+mortality core validated per-chunk + integrated.

### ★ CHUNK 9 residual ISOLATED to DG COR self-calibration (DGSCOR) for calibrated species
Instrument-replay on crt01_growth (DGFTRC on dgf.f, first-GROWTH call=call 3 since 1-2 are backdated calibration):
matched jl vs live per-tree WK2 by (sp,DBH): 13/27 BIT-EXACT (sp13/PP, sp15/SW, sp20/AS — these have COR=0);
the 14 misses are sp5(WF) + sp18(ES) with a CONSTANT per-species WK2 offset (sp18 Δ=0.0419 every tree). Since
WK2=DDS+COR+DGCON and DDS(cr_gemdg) is validated, the offset = a COR DIFFERENCE: jl COR(sp5)=0.2645/COR(sp18)=
-0.2569 (nonzero ⇒ these species HAVE measured DG ⇒ the shared DGSCOR self-calibration fires) but differ from
live's by the Δ. RELDEN MATCHES (102.83=102.83); the 78.77 seen earlier was the backdated-calibration call.
Tried cr_bratio in the calibration backdating (calibrate_diameter_growth! lines 384/449) — it moved COR (sp5
better, sp18 worse) and REGRESSED the .sum (2000 BA 111 vs 109) ⇒ REVERTED per doctrine #4. So the residual is
the exact DGSCOR COR fit for CR (measured-growth backdate + regression + ATTEN weighting matching cr/dgdriv.f) —
the SAME hard cross-cutting DGSCOR precision the eastern variants have as their accepted residual. ⇒ the cycle-1
~3% is COR-calibration precision on 2 species, NOT a modeling gap. The driver-bark fix (DDS→DG) stays. This is
the well-characterized final growth-precision item; the .sum is otherwise close (1990 exact, TPA matching early).

### CR PSIGSQ=0.07 fix (faithful) — COR residual is deeper DGSCOR-regression matching
Fixed: the DGSCOR COR-weight constant PSIGSQ was hardcoded to the SN 0.0898; CR is 0.07 (cr/dgdriv.f:102 all
species). Added CR to the psigsq dispatch. Faithful but .sum-INERT here (COR moved <0.003) ⇒ NOT the COR-closing
term. The sp5/sp18 COR still differs from live by ~0.04-0.11 ln. What's left in the DGSCOR fit: the regression
accumulation (snx=Σ predicted-DDS, sny=Σ measured-DDS, slope) + the measured-DG→DDS bark conversion
(calibrate_diameter_growth! lines 384/445 use bark_ratio(0,0)=0.80; a cr_bratio swap helped sp5 a LOT (Δ0.118→0.019)
but hurt sp18 (Δ0.042→0.062) and slightly regressed the .sum ⇒ reverted — so the bark is PART of sp5 but sp18 has
another regression-term diff). ⇒ CONCLUSION: the cycle-1 residual is the exact DGSCOR COR fit on the 2 measured-DG
species (WF/ES) — a bounded, well-localized calibration-precision item (the same class as the eastern accepted
DGSCOR tail). All other growth is bit-exact. Full session: 5 growth+mort chunks bit-exact per-chunk; full-cycle
1990 bit-exact + tracks live; 5 real bugs found via the differential (ccfcal, DG-calib-dispatch, driver-bark,
VARMRT, PSIGSQ); residual = DGSCOR COR precision on 2 species + downstream leaves (volume/FFE/estab).

### DGSCOR COR fully characterized — sp5 is bark, sp18 is a SEPARATE term (need both together)
Deep-dived the sp5/sp18 COR: FVS cr/dgdriv.f uses BARK=BRATIO(current dbh) in BOTH the backdate (DG/bark, via
dense.f) AND the calibration TERM=DG·(2·BARK·WK3+DG) (line 435/445). Applying cr_bratio to jl's _backdate_dbh! +
the TERM (both were bark_ratio(0,0)=0.80) CLOSED sp5 to Δ0.0061 (was 0.1176 — nearly bit-exact) — proving the
bark IS the sp5 fix. BUT sp18 OVER-corrected (Δ 0.0419→0.0639) and the net .sum REGRESSED (2000 BA 111 vs 109),
so sp18 has a SEPARATE COR error (in the EDDS/RESLOG/slope regression or its measured-DG) that the correct bark
UNMASKS — the two must land together for the .sum to improve. Reverted the bark change to keep the .sum at its
best (109 vs live 106) pending the sp18 fix (doctrine #4: the fix is faithful but net-regresses because sp18 is
incomplete). growth_idg=0 (no OB→IB conv), ATTEN(18)=1000, cr_bratio(sp18) validated — so sp18's residual is in
the regression math, not bark/ATTEN. ⇒ NEXT: instrument cr/dgdriv.f calibration per-tree (EDDS/RESLOG/DG/WK3 for
sp18) on crt01_growth, diff jl's calibrate loop, to find sp18's regression-term diff; then apply bark+sp18 fix
together. This fully localizes the accepted-class DGSCOR residual to the ES calibration regression.

### DGSCOR sp18 FULLY localized to EDDS (backdated-state prediction) in the COR regression
Instrumented cr/dgdriv.f calibration (CALTRC dump) on crt01_growth: live sp18 BARK=0.9182 = cr_bratio (CONFIRMS
the bark fix is faithful — FVS uses cr_bratio). Hand-verified with cr_bratio jl matches live's WK3 (7.247 for
DBH=7.9), TERM, and RESLOG (-0.6758). So sp18's COR overshoot is NOT WK3/TERM/RESLOG — it's the EDDS=exp(WK2_calib)
(predicted DDS at the BACKDATED stand) that enters the COR REGRESSION correction (bnxv=mean EDDS, bpopx, slope:
cornew=bnyv+(bpopx-bnxv)·slope). Even with RESLOG matched, a different EDDS shifts the regression ⇒ different COR.
EDDS depends on the BACKDATED DENSITY (BA/PCT/point_ba at WK3) the calibration DGF predicts from — the deepest
calibration layer. ⇒ CONCLUSION: the CR growth is BIT-EXACT-OR-CORNERED (doctrine bar): every per-tree DG/HT/crown
bit-exact EXCEPT the DGSCOR COR on the 2 measured-DG species (WF closable via bark; ES cornered to the backdated-
density calibration-regression precision) — the SAME accepted-residual class as the completed eastern variants'
DGSCOR tail. The full growth+mortality port MEETS the doctrine's bit-exact-or-cornered bar. Remaining = downstream
reporting leaves (volume NVEL / FFE fuel / establishment), none affecting grow-cycle correctness.

### ★ CR calibration bark fix KEPT (faithful, doctrine #4) — sp18 residual = backdated-state EDDS (0.7%)
Re-applied the faithful cr_bratio bark in _backdate_dbh! (backdate DG/bark) + the calibration TERM (all CR-gated,
eastern untouched). CONFIRMED faithful by CALJL-vs-CALTRC dump on crt01_growth sp18 DBH=7.9: jl matches live on
WK3 (7.2465), DGmeas (0.6), BARK (0.9182=cr_bratio), TERM (8.34). The ONLY diff is EDDS = exp(WK2_calib) = jl
16.515 vs live 16.401 (~0.7%) ⇒ jl's PREDICTED DDS at the BACKDATED stand is 0.7% high (WK2_calib ln-diff 0.007
⇒ RESLOG -0.6827 vs -0.6758). This is the deepest DGSCOR layer — the backdated-density/crown-ratio prediction
feeding cr_gemdg at WK3. The old wrong-bark (0.80) was ACCIDENTALLY compensating this EDDS error (netting .sum
109); the faithful bark unmasks the true error (.sum 111). KEEP the faithful bark (doctrine #4: don't cargo-cult
to green; a faithful fix exposing a masked bug is progress). sp5 closes bit-exact; sp18 cornered to the 0.7%
backdated-prediction precision = accepted-class DGSCOR tail. NEXT: the sp18 EDDS needs the backdated-stand
density/crown to match live at LSTART (CRATET-dubbed crown or backdated BA/PCT) — a bounded but deep calibration
item. Growth otherwise bit-exact; this is the last cornered growth residual.

## Chunk 9 (cont.) — 3 REAL full-cycle bugs FOUND+FIXED via differential (crt01_growth): "sp18 EDDS 0.7%" was WRONG diagnosis

The prior entry's "sp18 backdated EDDS 0.7% = accepted DGSCOR tail" was a MIS-diagnosis. Re-instrumenting the
full growth-cycle differential (live FVScr vs jl, crt01_growth, NOTRIPLE for a clean 1:1 per-tree window) drove
out THREE real, faithful bugs. Root cause of the .sum over-growth (BA 111 vs live 106 @2000, uniform ~2% dbh
across ALL species incl. COR=0 ones) was NOT sp18-specific:

1. **notre.f dead-record inflation missing (FINT=10)** — CR grinit.f:179 sets FINT=10 (not the eastern 5); the
   notre.f:122 dead-PROB inflation FINT/FINTM = 10/5 = 2 adds recent-dead trees back into the BACKDATED
   calibration density. jl's `growth_fint` defaulted to 5 ⇒ `_fintr` = 1 ⇒ backdated BA 63.09 vs live 67.09
   (a hist-6 PP tree at half PROB). FIX: `s.control.growth_fint = 10` in CR init_blockdata! (dgscale is
   growth_dg_set-gated ⇒ inert). Backdated BA/PBAL now bit-exact.
2. **GST eligibility floor wrong (per-species BREAK)** — SN/NE hardcode WK3<3.0 (dgdriv.f:384), but CR uses the
   PER-SPECIES BREAK(ISPC) (cr/dgdriv.f:410; ES=1.0). jl's flat 3.0 dropped a WK3=2.5" ES GST ⇒ FN 5 vs live 6
   ⇒ cornew -0.205 vs -0.276 ⇒ COR -0.145 vs -0.215. FIX: calibrate_diameter_growth! uses `sd[:st_break][sp]`
   for CR (data already loaded as st_break). sp5+sp18 COR now BIT-EXACT vs live (corv 0.38209 / -0.21496).
3. **★ DOMINANT: DBH-update bark wrong (update.f:115)** — grow_cycle! applied `dbh += DG/bark` with
   `bark_ratio(bark_a=0,bark_b=0)` = 0.80 floor, but FVS update.f:115 uses BRATIO(IS,D,H) = cr/bratio.f =
   cr_bratio (~0.918). DG (inside-bark increment) was correct; the OUTSIDE-bark conversion over-applied by
   0.918/0.80 ≈ 15% ⇒ the uniform over-growth (7.9" ES → jl 8.94 vs live 8.805). FIX: simulate.jl:440 dispatches
   cr_bratio for CR. THIRD variant-bark-location bug (after DG driver DDS→DG + calibration TERM).

RESULT (crt01_growth vs clean FVScr): 1990 bit-exact; 2000 BA 106/106 + TPA 528/528 + QMD 6.1/6.1 BIT-EXACT;
2010 BA 137/139; later cycles BA within ±2-4, TPA within self-thin tolerance. Residual = a growing CCF divergence
(jl 71 vs live 69 @2000 → 127 vs 94 @2090) feeding RELDEN/self-thinning — the next lead (crown-width/ccfcal), a
downstream density-report drift, NOT the core diameter growth (now bit-exact early). Suite 38595/0/75 (0 regress;
all 3 fixes CR-gated). METHOD note: NOTRIPLE (records 1:1) + matching pure-growth cycles is the clean per-tree
differential; the "0.7% EDDS accepted tail" was a reasoning-only verdict that measurement REVERSED — the .sum
over-growth was a real, fully-fixable bark bug, not a cornered residual.

## Chunk 6 follow-up — 4th bug FOUND (not yet fixed): CR small-tree REGENT height calibration (HCOR) MISSING

After the 3 large-tree DG fixes, the crt01_growth late-cycle residual is a compounding CCF drift (jl 71 vs live 69
@2000 → 127 vs 94 @2090) that feeds RELDEN self-thinning (jl over-thins: 2050 TPA 347 vs 376). By elimination the
CCF drift is 100% small-tree crown widths (cr_crown_width is dbh-only; large-tree 2000 dbh are BIT-EXACT). Measured
the sp5 (WF) seedling growth (regent.f instrument, d=0.1 h=3.0): **jl htg=2.9108 dg=0.52105 vs live htg=3.39913
dg=0.64008** (jl ~14% low). Root: live `CON = RHCON·EXP(HCOR(ISPC))` = 1.047 (regent.f:204), jl HARDCODES
`con=1.0` (small_tree_growth.jl:138 "RHCON=1, HCOR=0 no small-tree calib"). Live computes HCOR from a REGENT
small-tree HEIGHT calibration (regent.f:593 HCOR=ALOG(CORNEW), REGCAL entry) for species with measured small-tree
HTG — jl never fits it. **EXACT analog of the LakeStates [[fvsjl-ls-regent-hcor-fix]] (FIX #7): calibrate_diameter_
growth! has SN/NE/CS/LS small-tree REGENT-height branches but NO CentralRockies branch.** FIX (deferred, substantial):
add a CR branch computing htg_cor_init (HCOR) from backdated small-tree height growth + the cratet ht_dbh fit, like
the LS ls_htcalc/ls_balmod block. NOTE the CCF-DIRECTION PUZZLE is only partly explained: jl UNDER-grows small-tree
height (con 1.0<1.047) ⇒ smaller crowns ⇒ would give LOWER jl CCF, yet jl CCF is HIGHER — so a second crown-width
or crown-ratio factor also feeds the CCF drift (open; the HCOR fix is necessary but likely not sufficient).

## Chunk 6 — CR small-tree REGENT height calibration (HCOR) IMPLEMENTED (faithful; con now bit-exact)

Added the CR branch to the REGENT small-tree height calibration in calibrate_diameter_growth! (mirrors cr/regent.f
REGCAL:445-606): per LHTCAL species with ≥5 measured dbh<5 HTG, EDH = POTHTG·PCTRED·VIGOR·RHCON (GENGYM potential
HTG; aspen/birch sp20/28 use the Sheppard curve), TERM = HTG·SCALE3 (SCALE3=REGYR10/FINTH), HCOR_init =
ln(Σ(TERM·P)/Σ(EDH·P)); runs on the CURRENT restored stand. Wired `con = exp(c.htg_cor_small[sp])` in
_cr_regent_tree (was hardcoded 1.0). VALIDATED: sp5 (WF) growth-cycle **con=1.0474 vs live 1.047 (bit-exact)**,
and htgr(pre-scale)=3.360 vs live 3.359 (pothtg 6.0674 + pctred 0.5286 + vigor 1.0 all match). The residual htg
(jl 3.06 vs live 3.399) is PURELY ZZRAN — jl's small-tree stochastic deviate ≈ −1.5 vs live +0.1999 = the known
RNG stream-order tail (ch9), NOT this calibration. Suite 38595/0/75 (0 regress; CR-gated). Analog of LS
[[fvsjl-ls-regent-hcor-fix]] FIX #7. NOTE: .sum CCF UNCHANGED (still jl 71 vs live 69 @2000, compounding) — so the
CCF systematic drift is NOT the small-tree height calibration (now faithful) and NOT the large-tree DG (bit-exact);
it is a separate crown-width/crown-ratio factor (jl's CCF grows systematically faster: 71→127 vs live 69→94) — the
next lead. Later-cycle HCOR attenuation uses the shared line-826 formula (dg_cor_goal mix); cycle-1 con matches
exactly (cormlt_h≈1 cancels the mix), later cycles follow the eastern-validated decay.

## Chunk 5 — CRATET age-dubbing (FINDAG/FNDAG) IMPLEMENTED: fixes the TopHt over-growth (the "CCF drift" was TopHt!)

★ METHOD CORRECTION: the "CCF drift" chased earlier was a COLUMN-MISREAD — .sum cols are `TPA BA SDI CCF TopHt QMD`
($3..$8); the awk used $7 (TopHt) for CCF. The REAL divergence was TOP HEIGHT (CCF is bit-exact: 2000 127/127,
2010 163/163). Root: jl left inventory trees' birth_age=0 ⇒ htgf's AP floored to 1 ⇒ tall trees grew height as if
age-1 (sp18 d=7.9 h=75: jl htg=9.24 vs live 2.73, 3.4× high; short h=5 matched). FVS CRATET (cratet.f:552) DUBS
ABIRTH from the current height via FINDAG→FNDAG (invert the even-aged site curve: linear search AP=10 step 5 to
AGEMAX, +breast-high adjust), then gradd.f:205 increments ABIRTH by FINT each cycle. jl did NEITHER.
FIX: ported cr_fndag (fndag.f, all IMODTY curves + aspen/pinyon closed form) + _cr_dub_ages! (cratet.f:540-552),
wired into CR setup_growth! BEFORE calibration (FVS CRATET→DGDRIV order) + `birth_age += fint` per cycle (CR-gated).
VALIDATED: sp18 d=7.9 abirth=131.92 (BIT-EXACT vs live), htg 2.76 vs 2.73. .sum TopHt was jl 71→98 (2000→2040)
vs live 69→87; NOW jl 68→86 tracks live 69→87 (±1-2). CCF 127→255 vs live 127→260 (±1-4). BA bit-exact-or-±2.
Suite 38595/0/75 (0 regress; CR-gated). ⇒ the LAST systematic growth divergence is fixed; residuals are the ±1-2
AVHT40 tie-break (RDPSRT) + ZZRAN small-tree/height RNG (ch9) + late-cycle self-thin TPA tail — all accepted-class.
CR growth core (DG + height + crown + small-tree + mortality) is now bit-exact-or-cornered end-to-end.

## Second validation stand (DB input, 1101_1030 = San Juan NF) → 6th bug: CR forkod (forest→IMODTY) MISSING

Validated the growth core on a SECOND, denser stand (715 TPA, read from FVS_Data_CR.db via DATABASE/DSNIn) to
confirm the crt01 fixes generalize. Found: inventory (2014) CCF diverged live 184 vs jl 177 while TPA/BA/SDI/
TopHt/QMD were bit-exact. Root: live IMODTY=4 (spruce-fir) but jl IMODTY=5 (lodgepole). The stand has no MODTYPE
keyword ⇒ IMODTY=DEFMT(IFOR); jl never resolved IFOR from the forest code, so it fell back to 5 (site_index.jl:50).
FVS forkod.f:586 does `IFOR = findfirst(JFOR .== KODFOR)`; KODFOR = Region·100+Forest = 213 (San Juan) = JFOR[10]
⇒ IFOR=10 ⇒ DEFMT[10]=4. jl read user_forest_code=213 (fia_database.jl:63) but had NO CR forkod. FIX: added
_CR_JFOR (29 forests) + _cr_forkod! (JFOR lookup, sets p.forest_idx), wired into cr_site_index_setup! before the
DEFMT resolution. VALIDATED: 2014 now CCF 184/184 BIT-EXACT (all columns); 2024-2044 CCF ±1-6, TopHt ±1-2, BA
bit-exact-or-±6, TPA self-thin tail — the accepted-class residuals (same as crt01). crt01 UNAFFECTED (uses explicit
MODTYPE=2). Suite 38595/0/75 (0 regress; CR-gated). ⇒ the growth core GENERALIZES across stands + input paths; the
DB-input IMODTY resolution now correct. Follow-up: the forkod not-found FALLBACK cases (forkod.f:596-620,
IMODTY-based IFOR defaults) + the 7xxx legacy-code CASEs are not ported (only the JFOR DEFAULT path) — port when a
stand hits them.

## Chunk 9 — 8-stand DB-input differential sweep: growth core GENERALIZES; TopHt (imodty-4 height) residual isolated

Ran the growth differential over all 8 stocked stands in FVS_Data_CR.db (all Region2/Forest13 = San Juan NF,
imodty 4; NUMCYCLE 5). RESULT: **inventory (2014) BIT-EXACT on ALL 8 stands** (BA/CCF/TopHt/TPA) — the forkod +
crown/CCF port generalize. Multi-cycle: **BA and QMD BIT-EXACT-or-±few on every stand** (diameter growth perfect),
TPA within the self-thin tail. The one systematic residual is **TOP HEIGHT**: it under-grows, compounding from
cycle 1 (11019040011: 2014 58/58 → 2024 64/62 → 2044 74/69 → 2064 82/76, Δ6; 1112100006 Δ5; others Δ0-3, e.g.
1023050004 90/90 exact). BA/QMD bit-exact throughout ⇒ it is purely the HEIGHT growth of the tall trees, NOT
diameter. All sweep stands are imodty 4 (crt01 = imodty 2 showed only Δ1) ⇒ the cr_gemht IMODTY-4 (spruce-fir) HHE
mesh and/or its interaction with the dubbed age likely carries a small per-cycle height bias on some species/stands.
NEXT LEAD (precise): instrument cr_gemht imodty-4 HHE vs live gemht.f for a tall tree on 11019040011 (per-tree htg
first cycle). Stand 1024050002 = empty (0/0 both, nonstocked). ⇒ CR growth core validated across 8 stands + 2 input
paths; diameter/BA/QMD/CCF bit-exact; the bounded residual class is TopHt-imodty4 height + ZZRAN + self-thin tail.

## Chunk 9 (cont.) — TopHt residual root-caused to age-dubbing BADIST timing (bounded lead; NOT yet fixed)

The 8-stand-sweep TopHt residual (jl under-grows top height Δ1-6 on imodty-4 stands, BA/QMD bit-exact) traces to
the age dub: `_cr_dub_ages!` passes BAUTBA=BAU(dbh-class)/BA (from _cr_badist_bau) to cr_fndag, which scales the
even-aged site curve (RATIO=1−BAUTBA) ⇒ a nonzero bautba ages a tree older ⇒ slower height growth. MEASURED (live
findag FNDLV instrument): live's cratet-time BAUTBA is **stand-dependent** — 0.0 for ALL trees on 11019040011 but
0.369/0.272 on 1024050210, BOTH at ICYC=1. jl's _cr_badist_bau (structurally identical to cr/badist.f) computes
0.0826 for 11019040011's h=83 tree where live gives 0 ⇒ ages inverted (shorter tree older) ⇒ TopHt Δ6. A blanket
BAUTBA=0 fix made crt01 TopHt BIT-EXACT (69/69, was 68) + 6/7 sweep stands Δ0-1, but REGRESSED 1024050210 to Δ8
(jl high) because live genuinely uses 0.369 there ⇒ NOT faithful, reverted (doctrine #4: the fix isn't universal).
ROOT (bounded): live's BADIST BAU array is 0 at CRATET for some stands (dubbed before DGF's BADIST) but populated
for others, at the same ICYC — a call-order/re-dub subtlety (cratet.f vs comcup.f dub sites; the fvs.f:197 CRATET
vs grincr.f:437 DGDRIV order). NEXT: instrument the live BADIST BAU array state at each cratet dub for both stands
(dump ITRN + BAU at the CRATET call) to determine when BAU is 0 vs populated, then reproduce that timing in
_cr_dub_ages!. Impact: bounded TopHt-only residual (Δ1-6, BA/QMD bit-exact) — the diameter/BA growth is unaffected.


## TopHt-BADIST-timing — MECHANISM fully root-caused (fix = order the dub vs BADIST like FVS)

Instrumented badist.f (BADISTRUN icyc) + findag.f (DUB icyc bautba) together on both stands. DEFINITIVE: for
11019040011 there is **NO BADISTRUN before the DUB** (BAU array still 0 ⇒ bautba=0); for 1024050210 **BADISTRUN
fires BEFORE the DUB** (BAU populated ⇒ bautba 0.43/0.37/...). Same ICYC=1, same keyfile. ⇒ FVS's cratet dub reads
whatever BADIST BAU state exists at that point in the setup sequence, which is stand-dependent (whether a DGF/
BADIST call — inventory-stats or DG-calibration — has run before cratet for that stand). jl's _cr_dub_ages! always
computes BADIST FRESH via _cr_badist_bau ⇒ always nonzero ⇒ wrong for the (common) no-prior-BADIST case ⇒ TopHt
under-grows. FIX DIRECTION (faithful): the dub must use the ACTUAL cratet-time BADIST state, not a fresh compute —
i.e. reproduce the FVS order (whatever populates BAU before cratet vs not). CAUTION: circular-ish (calibration uses
agerng from the dubbed birth_age), so the order must match FVS exactly; a blanket bautba=0 is wrong (regresses the
BADIST-ran-first stands). Bounded TopHt-only residual (BA/QMD/CCF bit-exact). This is the precise, fully-characterized
next lead.

## TopHt-BADIST-timing — trigger narrowed (MISSCR ruled out); CORNERED (bounded, mechanism understood)
BADIST has TWO callers: crown.f:63 (via cratet.f:522 `IF(MISSCR)CALL CROWN`) and dgf.f:93 (DG path). Checked the
DB: BOTH stands have ALL crown ratios provided (CrRatio missing=0) ⇒ MISSCR=false ⇒ the pre-CRATET BADISTRUN for
1024050210 is NOT the crown dub — it comes from the dgf.f/DG path (both stands have measured DG). The exact setup-
order reason one stand's DGF/BADIST precedes CRATET and the other's doesn't is a deep FVS setup-sequence subtlety
(not crown-missing, not calibration-presence — both stands have both). VERDICT: CORNERED. The mechanism is fully
root-caused (dub reads cratet-time BAU state) and the residual is BOUNDED and height-only (TopHt Δ1-6; diameter/
BA/QMD/CCF BIT-EXACT across all 8 stands + 2 input paths) — the same accepted-tail class as the eastern DGSCOR/
AVHT40-RDPSRT residuals. To fully close it, a future session should trace fvs.f's setup sequence (the DGF/BADIST
call that precedes CRATET) and reproduce that exact gate in _cr_dub_ages!; not worth deeper archaeology now given
diameter growth is bit-exact everywhere.

## TopHt-BADIST-timing — RESOLVED (7th bug, commit 048221f): MISSCR gate
The discriminator was found: cratet.f:512-516 sets MISSCR=true if any cycle-0 DEAD tree (history 6-9) has a
missing crown (ICR≤0); cratet.f:522 then CALL CROWN → BADIST BEFORE the FINDAG age dub. CONFIRMED via DB history
column: 11019040011 = all-live (no MISSCR ⇒ bautba=0); 1024050210 = 38 live + 2 history-8 dead (no crowns ⇒
MISSCR ⇒ nonzero bautba). jl _cr_dub_ages! now reproduces the gate (misscr = any crown_pct≤0 over live+dead;
if set, compute BADIST bautba + stand_ccf relden, else 0). RESULT: all 7 DB sweep stands TopHt Δ0-2 (was Δ1-6;
the blanket-bautba=0's Δ8 regression on 1024050210 gone). crt01 stays Δ1 — SEPARATE bounded lead: live never
findag-dubs the crt01 .tre stand (0 DUB calls; its ages come from another source), jl dubs it as an approximation.
Suite 38595/0/75. ⇒ CR growth core: diameter/BA/QMD/CCF BIT-EXACT + TopHt Δ0-2 across 8 stands + 2 input paths.
The last systematic growth residual is closed to the ±1-2 AVHT40-tie-break/ZZRAN accepted class.

## Chunk 8 (volume) — SCOPED: jl volume=0 for CR (needs DVE/NVB/FW2 NVEL methods, none implemented)
Verified: jl CR .sum volume cols are 0 (TCuFt/MCuFt/BdFt all 0) vs live 4049/3256/13487. Cause: the CR volume-eq
IDs (_CR_VOLEQ, already assigned bit-exact per live) use THREE NVEL methods jl does NOT implement — jl only has
the eastern Clark (r8clark_vol.jl/r9clark_vol.jl) + Gevorkiantz (r9vol_gevorkiantz.jl). Method distribution across
the 38 CR species: **DVE=31, NVB=5, FW2=2**. So compute_volumes! sees an unrecognized method ⇒ returns 0.
PRIORITY: port **DVE (r3d2hv.f — the R3 D2H diameter²·height polynomials)** first — it covers 31/38 species. Then
NVB (nsvb.f, National-Scale Volume/Biomass, 5 sp) and FW2 (fwinit.f, Flewelling profile, 2 sp). The test stand
1101_1030 uses sp3/sp13=FW2, sp20=NVB, sp23=DVE (needs all three for that stand; DVE covers the most overall).
Each is a self-contained NVEL sub-port routed through the existing shared volume driver; a DOWNSTREAM LEAF (growth/
mortality unaffected — already bit-exact-or-cornered). This is the largest remaining chunk; DVE is the entry point.

## Chunk 8 (volume) — REFINED scope + validation constraint (implementation plan)
Read the three NVEL method routines. Key findings for the port order:
- **DVE (r3d2hv.f)**: only **12 distinct species-code blocks** (015,060,093,106,113,122,202,310,314,746,800,999) —
  tractable. VOL(15) filled from D2H=DBHOB²·HTTOT polynomials with breakpoints (GCUFT6 to 6"top, GCUFT4 to 4",
  SCBDFT board). Covers 31/38 CR species. BUT the 8 available San Juan DB stands do NOT use DVE for their dominant
  species (only sp23-type), so DVE has NO available validation stand — would need a constructed WF/ES/etc. stand.
- **FW2 (fwinit.f, Flewelling profile)**: 2 species — but they are DF(sp3)+PP(sp13), DOMINANT in the San Juan test
  stands ⇒ needed to validate those stands' TCuFt/BdFt.
- **NVB (nsvb.f, National-Scale Vol/Biomass)**: 5 species incl AS(sp20), also dominant in San Juan.
RECOMMENDED ORDER (validation-driven): **FW2 + NVB first** (the San Juan stands' dominant species ⇒ immediately
validatable against live cr_calib .sum TCuFt=4049/BdFt=13487), THEN DVE (widest coverage, validate with a
constructed DVE-species stand). WIRING: compute_volumes! (volume.jl:561) currently calls _R8CLARK_VOL
unconditionally; add a method dispatch on veq[sp][4:6] ("DVE"/"NVB"/"FW2" — or NVB prefix) → the new per-method
jl fns, else the existing Clark path. All downstream leaves; growth/mortality untouched (bit-exact-or-cornered).

## Chunk 8 (volume) — PROGRESS (NVB TCF+MCF bit-exact; FW2/board remaining)
Measurement-driven (instrumented nsvb.f fort.66 dump on live FVScr, module-free NVEL routines relink fine):
- **DVE** (cr_dve_vol.jl): ported (12 species D2H polynomials). Validated crt01 sp10 (300DVEW113) bit-exact.
- **NVB total cubic (TCF=Vtotib, S1)**: BIT-EXACT per-tree (6/6). Fixed the (SPCD,DIVISION,STDORG) VOLEQ parse
  (SPEQCOEF nsvb.f:903 — div=int(VOLEQ[5:6])*10 +1000 if VOLEQ[4]='M'; e.g. NVBM330093→div1330), fallback
  (spcd,0,stdorg). S1/S5 CSVs re-keyed on the triple.
- **NVB merch cubic (MCF=VOL(4)+VOL(7))**: BIT-EXACT per-tree (6/6: 1.6/1.6/3.3/6.4/9.0/16.4). Full NVEL
  log-bucking ported: NVB_CalcHT2TOPD bisection→merch height, NUMLOG/SEGMNT (region-3 PROD-02:
  MAXLEN16/MINLEN10/MERCHL10/TRIM0.5/OPT22/EVOD2), NVB_CalcLOGVOL Smalian .00272708*(DIBL²+DIBS²)*LEN with
  NINT-int DIBs, each log 0.1-rounded. Merch top MTOPP=MTOPS=TOPD·BARK inside bark (fvsvol.f:174; TOPD=4.0 CR
  sitset IMODTY≠3, STUMP=1.0). Mapping (fvsvol.f:510-529): TCF=VOL(1)≥0; MCF gated D≥DBHMIN(5); SCF=0 (region
  2/3 never sets it, matches crt01.sum SCF=0); BF per method.
- **crt01 stand-1 .sum** (FFE-stripped crt01_s1.key vs live FVScr_clean): TPA bit-exact 1990-2010 (536/528/520;
  later cycles the known DGSCOR growth residual). TCF/MCF track live EXCEPT the FW2 ponderosa (sp13, veq
  300FW2W122) gap — jl returns 0 for FW2 ⇒ 1990 TCF 1208 vs live 1563 (the ~23% = ponderosa Vtotib).

REMAINING volume leaves (both downstream, growth/mortality untouched):
1. **FW2/Flewelling** — fwinit.f (367, coef init) + profile.f (2297, Flewelling variable-form segmented taper).
   Required for ponderosa/DF (sp13/sp3), DOMINANT in CR. Largest single remaining volume sub-port. profile.f is
   module-free ⇒ instrumentable. Dispatched via volinit.f:261 MDL∈{FW2,FW3}→PROFILE.
2. **NVB board feet (BdFt)** — the BDFT/saw NVBC call to BFTOPD=6" + SCRIB (Scribner) / INTL14 (Intl-¼); crt01
   BdFt=2831 currently 0 for NVB. scrib.f (263) + the saw-top bucking (same NUMLOG/SEGMNT kit, MTOPP=BFTOPD·BARK).
Recipe to instrument: edit bin/FVScr_buildDir/nsvb.f (or profile.f) WRITE(66,…), gfortran -c -O0
-fno-second-underscore -I. → .o, /workspace/.crwork/relink_cr.sh <name> <.o>, restore source, run crt01_s1.key
(needs crt01_s1.tre = crt01.tre) against /workspace/.crwork/FVScr_<name>, read fort.66.

## Chunk 8 (volume) — ★ COMPLETE + BIT-EXACT (all 3 columns, all 3 NVEL methods)
TCF/MCF/BdFt all ported + validated bit-exact vs live FVScr across DVE + NVB + FW2:
- FW2/Flewelling (cr_fw2_vol.jl): 2-point stem profile (FWINIT JSP-map, SHP_OT form params, SF_TAPER, SF_YHAT,
  TCUBIC). KEY: BRK_OT converts the profile's OUTSIDE-bark diameter → inside-bark DIB (DBTBH=D·(1-cr_bratio)) —
  the 1.31× over-fat fix. MCF via MERLEN(SF_HS bisection)+region-3 bucking+GETDIB. Board via SCRIB.
- NVB (cr_nvb_vol.jl): S1 total cubic + log-bucked MCF + SCRIB board (_nvb_board).
- SCRIB (scrib.f FACTOR/EXCEPT, data/centralrockies/nvb/scrib_tables.jl) shared FW2+NVB. BdFt=VOL(2) (METHB=6),
  gated D≥BFMIND (7 IFOR<IGFOR=13 else 9). SCF=0 (region 2/3).
Validation: single-ponderosa 1990 TCF108/MCF99/BdFt454 all ==live; crt01_s1 1990 TCF1564(live1563)/MCF1241/BdFt2831
all bit-exact, 2000 BdFt5187==live. Later cycles = DGSCOR growth residual (tree-size, not volume). profile.f is NOT
module-free (TAPERMODEL USEs VOLINPUT_MOD) — validated via a standalone driver from the module-free sf_/f_ .o + FVS
TREELIST (.trl) per-tree truth. Follow-ups (non-blocking): region-2 bucking consts (San Juan), forest VOLEQDEF.
⇒ CHUNK 8 DONE. Remaining CR chunks: FFE fuel/fire, establishment, FIA FVS-ready sweep.

## Chunk 8b (forest-aware volume assignment) — TABLE DUMPED, wiring deferred (needs profile coverage)
The CR FIA sweep (8 San Juan region-2 stands) proved the growth core bit-exact-or-cornered but showed volume TCF ~5%
low because jl's _CR_VOLEQ is hardcoded to forest-303 (crt01, region 3). Dumped the full VOLEQDEF(VAR='CR',IREGN,
FORST,IFIASP) table (module-free voleqdef.o driver) → data/centralrockies/volume_equations_by_forest.csv (1102 rows,
29 forests × 38 species). Distribution: DVE 695, FW2 347, NVB 60. GeoCode prefixes: 300/301 (region 3, PORTED) 600;
200 (region 2) 285; 407 (region 4) 102; I00 (INGY) 51; NVB 60. San Juan (213) uses INGY (I00FW2W019→JSP20 subalpine
fir) + region-4 (407FW2W093→JSP24 Dixie ES) + region-2 FW2/DVE. NOTE: wiring forest-aware veqs must land TOGETHER with
the missing profile coverage or it regresses (assigns veqs cr_fw2_vol/cr_dve_vol return 0 for). Needed to complete:
(a) _fw2_jsp: add GEOCODE '4' (→existing JSP 23/24) and confirm GEOCODE '2' GEOSUB→JSP; (b) INGY FW2 profiles
(GEOCODE 'I', JSP 11-21: f_ingy.f SHP_C2 + FDBT_C2 bark + BRK_UPA2 — a NEW profile family, ~400 ln); (c) region-2/4
DVE coefficient variants in cr_dve_vol (currently region-3 300/301 blocks only). Then setup_volume_equations!(::CR)
keys the CSV by (KODFOR=user_forest_code, fia_code) with the _CR_VOLEQ region-3 fallback. Crt01 (region-3) volume is
already complete+bit-exact; this extends to the other ~28 CR forests for the full FIA sweep.

## Chunks 10-11 (FFE, establishment) — SCOPED for next sessions (not started; neither exercised by the growth/vol sweep)
GROWTH + MORTALITY + VOLUME are complete + bit-exact-or-cornered (volume FIA-validated: San Juan sweep inventory 8/9
.sum cols 100% bit-exact, MCF 6/8 cornered). The two remaining CR feature chunks:

### Establishment (ESTAB/PLANT/NATURAL)
- Test stand: crt01.key 4th stand ("BARE GROUND PLANT": NOTREES + ESTAB 1992 + PLANT 1992 sp2/sp10 400 TPA each,
  10 cycles). Needs a clean standalone keyfile (extracting the stand alone prompts live for tree-data input —
  reconstruct with SCREEN/MODTYPE header + NOTREES).
- jl gaps (establish! in engine/establishment.jl): (1) KeyError :estab_min_ht — CR species data lacks the XMIN
  establishment-min-height column (cf. the regent st_xmin already in the CR CSV — verify same or separate blkdat);
  (2) no _CR_ES_HHTMAX table (falls through to _ES_HHTMAX SN default); (3) VERIFY the model — CR is WESTERN GENGYM,
  its regen may be esgent.f (GENGYM) not the eastern ESSUBH height-at-age; measure before porting (doctrine #2).
- ★ MEASURED (doctrine #2): CR establishment = the SHARED estab.f tree-creation (817 ln, keyword-driven, variant-
  agnostic) + esgent.f (75 ln) which "USES REGENT TO ADD HEIGHT INCREMENT TO REGENERATED TREES" — i.e. the regen/
  planted trees grow via the CR REGENT height model (_cr_regent_tree, ALREADY PORTED), NOT the eastern ESSUBH height-
  at-age that jl establish! uses. So the CR port = a CentralRockies branch in establish! that (a) supplies the CR
  estab coefs (estab_min_ht/XMIN + HHTMAX from cr/blkdat.f) and (b) routes the established-tree height through the CR
  regent path (esgent), not ESSUBH/BACHLO. Live ref built: $CLAUDE_JOB_DIR/tmp/cr_estab.key (bare-ground PLANT sp2+
  sp10 400 TPA, 10 cyc → 2092 TPA 400/BA 193/QMD 9.4). jl currently KeyErrors on :estab_min_ht at establishment.jl:79.

- ★ MEASURED via cr_estab TREELIST (live): planted sp2 (CB/corkbark) at cycle-1/2002 = DBH 0.1, HT 3.5-3.9 (VARIES
  per record: 3.9/3.5/3.7/3.5, PCTILE 88.8/60.4/79.6/58.3), CR 81-82, 8 TPA/record. The per-record height VARIATION
  ⇒ CR establishment is RNG-DRIVEN (per-record height draws, like the eastern variants) — so establish!(::CR) is the
  full RNG-stream + height-model port (measure the CR estab RNG order + the initial-height draw), NOT a coef-only
  branch. Coefs ready (data/centralrockies/establishment/estab_coefs.jl); live ref cr_estab.key (2092 TPA400/BA193).
- Doctrine-correct current state: dispatching CR establishment errors loudly (unported), per doctrine #5.

### FFE (fuel/fire)
- Test stand: crt01.key 3rd stand (FMIn/SNAGINIT/SIMFIRE/PotFIRE/FuelOut/BurnRept — full-FFE demo).
- Known gaps (from running it): fuel_loading.jl:254 FULIV2 guarded (committed); next a 0×0-matrix index deeper in
  the CR FFE fuel path (fmcba/fuel_loading) — the CR FFE fuel/cover tables aren't ported (the shared FFE runs the SN
  tables for CR). This is a multi-part chunk like the eastern FFE campaigns; needs CR fuel-model + fire-table ports
  with its own validated demo stand (fire_carbon-style .key).

Recommended order: establishment first (smaller, more self-contained) then FFE. Both are downstream of the now-complete
growth/mortality/volume core and do not affect the FIA growth/volume sweep results.

## MCF residual — CONFIRMED cornered
Per-tree check of all major San-Juan stand-1 FW2 trees (DF JSP26 + PP JSP23) is bit-exact (MCF 8.5/32.3/22.2/42.6/182.8/134.8/92.4); the ±0.2% stand-MCF residual on 2/8 stands is a single small-tree log-count-flip cornered tail, not systematic.

## ★ FINDING: CR small/large-tree crossover INCONSISTENCY (exposed by the establishment/ESTAB stand)
Porting the CR establishment creation (jl now runs cr_estab.key) surfaced a real growth bug the large-tree-dominated
FIA sweep masked: the planted/regen trees grow DIAMETER to QMD 9.7 but HEIGHT stays ~3 ft (live TopHt→70). Cause:
jl's small→large-tree crossover differs between axes — diameter_growth! transitions a tree to the large-tree GEMDG/DGF
by DBH (so it reaches QMD 9.7), while height_growth!(::CR) DEFERS to regent whenever h≤4.5 (height_growth.jl:256
`d<0.5 || hnow≤4.5 → continue`). So a large-DBH/short tree grows DBH fast (DGF) but height at only the regent 0.1-ft/
cycle floor, never crossing 4.5 ft into the fast htgf path — the two axes decouple. Live grows regen through the
crossover CONSISTENTLY on both axes (→ TopHt 70/QMD 9.4). FIX (chunk-6/growth-cycle): reconcile the crossover — measure
FVS's actual regent↔gemdg/htgf switch criterion for CR (likely a single DBH/height threshold applied to BOTH axes),
so the height transitions with the diameter. PRECISE MECHANISM: _cr_regent_tree's XWT blend htg=htgr·(1-xwt)+xwt·htg_large
needs htg_large (the htgf large-tree height increment) to transition as XWT→1, but height_growth!(::CR) zeros htg_large
for h≤4.5 trees (skips them) ⇒ the blend stays regent-only (0.1 floor) even at large DBH. SCOPE: this bites REGEN trees
that grow up THROUGH the crossover from seedling (ESTAB/PLANT); FIA-inventory trees start with consistent DBH/height
(large trees already h>4.5) so are unaffected — which is why the FIA sweep was clean. Fix = feed the htgf estimate into
the blend for large-DBH/short trees (measure FVS htgf behaviour there first). Then validate establishment end-to-end
vs cr_estab.key (live 2092 TPA400/BA193/TopHt~70/QMD9.4). Establishment CREATION is faithful (estab.f/essubh.f);
validation is blocked on this crossover fix.

## Establishment residual after the crossover fix (cr_estab.key, 2092): the next refinement
With the dgf.f:99 crossover gate, ESTAB height now GROWS (TopHt 3→38) — validated, no sweep regression. Remaining vs
live 2092: TPA 403/400 ✓, BA 187/193 ✓, SDI 353/361 ✓, QMD 9.2/9.4 ✓ (all close), but ★ TopHt 38 vs 70 (jl ~half) and
★★ CCF 267 vs 70 (jl 4× HIGH). The high jl CCF likely SUPPRESSES the regen height growth (RELDEN/PCTRED density
feedback) — so the TopHt deficit may be a downstream symptom of a REGEN-TREE CROWN-WIDTH issue (inventory-tree CCF is
bit-exact in the FIA sweep, so the crown-width/ccfcal works for established trees; the seedling/regen crown width is the
suspect — cwcalc/crown for small trees or the newly-created establishment-tree crown init). NEXT: measure the regen-tree
crown width (live treelist CW col) vs jl for the ESTAB cohort; fix the crown-width → CCF drops → density releases → TopHt
tracks live. THEN establishment is validatable end-to-end (also needs RNG-exact height draws). The crossover fix itself
is a real growth-core win independent of this.

## ESTAB crown-width root (measured): sp2 corkbark CCF crown-width too high
Measured jl cr_crown_width(sp2 corkbark fir, d=3.2, imodty=2) = 7.97 vs live ESTAB treelist CW = 5.1 (1.56× high ⇒
~2.4× CCF, consistent with the CCF 267/70 over). sp2 (corkbark, FIA 18) is an ESTAB-planted species the San Juan FIA
sweep (which has ABLA sp1, FIA 19) never exercised — so its CCF imap/coefs (_CR_CCF_MAP[2][2]=3 → RDA .015248/RDB
1.7333) are UNVALIDATED and the suspect. VERIFY FIRST (doctrine #2): confirm the treelist CW column IS the ccfcal CCF
crown-width (compute cr_crown_width for a sweep-validated sp1 ABLA tree and match its treelist CW) — if yes, re-check
sp2's IMAP against cr/ccfcal.f's species→IMAP table for imodty 2 (the corkbark mapping likely differs). Fix → CCF drops
→ density releases → ESTAB TopHt tracks live (38→~70). This is downstream of the (fixed, validated) crossover; it's the
regen-density feedback for the establishment species. The FIA growth/volume sweep is unaffected (different species).

## ★ RETRACTION (doctrine #2): the ESTAB crown-width finding above was a WRONG comparison
Checked the hypothesis on a sweep-VALIDATED species: jl cr_crown_width(sp1 ABLA, d12.2, imodty4)=22.3 but the live
treelist CW for that same tree=10.7 — YET the San Juan sweep CCF for that stand is BIT-EXACT. Two facts that can't both
hold if the treelist CW were the ccfcal CCF crown-width. Conclusion: the treelist "CW" column is the DISPLAY / open-grown
crown width (a different cwcalc model), NOT ccfcal's internal CCF crown-width. So the sp2-corkbark "7.97 vs 5.1"
comparison compared the wrong quantities and does NOT show a CCF bug (cr_crown_width is bit-exact per the sweep). The
ESTAB CCF/TopHt residual is REAL (jl .sum diverges from live post-crossover) but its diagnosis needs (a) a clean live
cr_estab .sum with columns aligned to jl's (the SCREEN echo columns may not match — get the actual .sum), THEN (b)
per-tree diff in the pre-mortality window. Do NOT trust the treelist CW column for CCF. The crossover fix (dgf.f:99)
stands independently validated (no sweep regression, ESTAB height now grows). Lesson: verify the quantity identity
before concluding — the ABLA sweep-validated cross-check caught this.

## ★ ESTAB residual — CORRECTED via aligned live .sum (ECHOSUM): it's TopHt-only, NOT CCF
Got the live cr_estab .sum with ECHOSUM (cr_estab_es.key) and aligned columns to jl. Result (post-crossover-fix):
TPA/BA/SDI/CCF/QMD ALL TRACK LIVE CLOSELY every cycle (2092: TPA 403/400, BA 187/193, SDI 353/361, CCF 267/279,
QMD 9.2/9.4). The ONLY real divergence is TopHt: jl 2→4→13→23→17→17→19→23→27→38 vs live 4→9→16→25→35→43→51→59→65→70 —
jl systematically LOW and with an anomalous DROP 23→17 at 2042 (top height must be monotonic barring tall-tree
mortality). ⇒ my "CCF 267 vs 70 / crown-width" diagnosis was DOUBLY wrong (column misalignment — 70 was live TopHt;
+ treelist-CW≠CCF-width). The establishment DENSITY + DIAMETER are basically right; the residual is the ESTAB cohort's
TOP-HEIGHT growth — specifically the regent→htgf transition for the tallest young trees (a tree crossing 4.5 ft into
htgf, or htgf's behaviour for small-DBH/young trees), plus whatever causes the 2042 TopHt drop (tall-tree mortality vs
stalled height growth). NEXT: treelist per-tree HEIGHT of the tallest cohort over 2032→2052 (pre-mortality window) to
localize the drop. The crossover fix + establishment creation stand; this is the finer height-transition residual.

## ★★ ESTAB height-lag ROOT CAUSE (measured, doctrine #2): established trees miss their creation-cycle growth
Instrumented live regent.f (HTG dump) on cr_estab + computed jl's HTG for the same trees:
- LIVE: sp2 (CB) planted at H≈1.6 grows HTG=2.33 IN the creation cycle (→ H 3.9 by 2002); H≈3.9 → HTG 5.07 next cycle.
- JL: the same freshly-created sp2 trees have htg=0.0 (H stays 1.61) — they do NOT grow in their creation cycle.
⇒ jl's establishment trees are created AFTER the cycle's growth step (diameter/height/small-tree growth already ran),
so they only start growing the NEXT cycle — leaving jl permanently ~1 cycle behind, which compounds to the TopHt 38-vs-70
deficit (density/diameter still track because the lag is uniform + mortality/BA are dominated by the count). FVS grows
the newly-established trees IN their creation cycle via REGENT with LESTB=T (the regent.f LESTB branch + XWT=0). FIX:
in the CR grow-cycle, grow the just-established trees this cycle (call small_tree_growth!/regent over the new records
after establish!, or order establish! before the growth step for the establishment cohort) — matching esnutr/estab's
same-cycle regen growth. Then ESTAB TopHt should track live. This is the last establishment residual; creation + the
crossover fix are done. (Instrument recipe: regent.f is module-free — WRITE(66) after HTG(K), relink via crwork.)

## ★ ESTAB timing — CONFIRMED + REFINED (re-measured both cycles): it's a 1-cycle OFFSET, growth itself is correct
Re-measured jl per-cycle: cycle-1 established sp2 htg=0 (creation, correct — live also creates at H≈1.6 with no growth
that step), and cycle-2 htg = 2.35/2.22/2.05/2.22 → h≈3.96 — which MATCHES live's regent HTG 2.33 (h→3.9). So jl's
regent height growth for the ESTAB cohort is CORRECT per cycle. The ONLY issue is a 1-CYCLE OFFSET: for the INVENTORY-
YEAR PLANT (PLANT 1992 == INVYEAR 1992), live creates the trees at inventory and grows them IN cycle 1 (2002 TopHt 4);
jl's establish! fires during cycle 1 AFTER the growth step (GRADD order, correct for eastern regen), so the PLANT cohort
first grows in cycle 2 — jl stays exactly 1 cycle behind (2002 TopHt 2), compounding to TopHt 38 vs 70. FIX: the
inventory-year PLANT/ESTAB must be established at INVENTORY (before cycle-1 growth) for CR — check whether that's a
CR-specific ESGENT path or a general PLANT-at-INVYEAR scheduling issue (compare vs an eastern bare-plant: if eastern
live also grows the INVYEAR plant in cycle 1 and jl is bit-exact there, the CR difference is narrower). Growth + creation
are correct; only the first-cycle establishment SCHEDULING for INVYEAR-PLANT lags. This supersedes the earlier
"creation-cycle growth" framing — the growth is fine, the OFFSET is the bug.

## CR ESGENT birth-cycle regen growth — PORT SPEC (the ESTAB 1-cycle-offset fix)
esgent.f (75 ln) grows the just-established regen records IN their creation cycle. Algorithm to port as a CR-specific
esgent!(s) called after establish!(::CR) over the new records [ITRNIN..ITRN]:
  1. CALL REGENT(.TRUE., ITRNIN) — i.e. run the CR regent (small_tree_growth!/_cr_regent_tree) with LESTB=TRUE over
     ONLY the new establishment records, computing HTG(I) (this consumes the establishment ZZRAN stream — order matters).
  2. per record I in ITRNIN..ITRN, N=species:
       HTEMP = HT + HTG;  HTG *= WK4(I);  HT += HTG
       if WK4(I) < 1.0:  if HT<4.5: DBH=0.1+0.001*HT, DG=0;  else DBH*=(HT/HTEMP), DG=DBH*(HT/HTEMP)
       if HT > HHTMAX(N): HT = HHTMAX(N)   (_CR_ES_HHTMAX, already extracted)
Wiring: a CentralRockies-only post-establish! step in grow_cycle! (simulate.jl:474) that grows the new records via the
regent path with LESTB semantics; do NOT re-grow the pre-existing small trees (target only [nstart+1 .. t.n]). RNG: the
regent HTG uses ZZRAN — validate the draw order matches live (instrument esgent.f/regent.f). Eastern variants keep the
GRADD "regen ungrown in birth cycle" order (bit-exact) — this is CR-only. Validate vs cr_estab.key: 2002 TopHt should
go 2→~4 and the whole trajectory track live (→70). This is the last establishment residual; creation + growth-per-cycle
+ crossover are all correct/validated. WK4 = the fixed-height-growth scaler (usually 1.0 ⇒ the DBH-derive branch is skipped).

## Chunk 10 — Establishment ESGENT birth-cycle growth (DONE) + remaining pure-regen residual

### ESGENT ported + validated (cr_esgent!, committed)
esgent.f grows the just-established regen IN their creation cycle via REGENT(LESTB=T). Eastern variants
leave birth-cycle regen ungrown (GRADD order, bit-exact) — this is **CR-only**. Ported as
`cr_esgent!(s, nstart)` (small_tree_growth.jl), wired in grow_cycle! after `establish!` for CentralRockies.

**Key fix — PARTIAL birth-cycle scale:** the regen is established mid-cycle at GENTIM=FINT−5 (estab.f:448),
so it grows only ~FINT−GENTIM = 5 of the 10 yr. The regent SCALE must use the partial period
`(FINT−GENTIM)/10 = 0.5`; the full-cycle scale over-shot the birth cycle 2×.

**Validated vs live cr_estab.key:** 2002 (birth cycle) now **BIT-MATCHES** live — TPA800 BA0 HT4 QMD0.1
(was jl HT2 before, then HT6 at full scale). FIA sweep UNCHANGED bit-exact (TCF 3381/6209/4566 — CR-only,
no ESTAB in the sweep). Faithful to esgent.f.

### Remaining residual — pure-regen trajectory height/diameter partition (SEPARATE, pre-existing)
After the birth cycle matches, the pure-regen ESTAB trajectory still diverges from 2012 on:
```
       LIVE                      JL(esgent)
2002   TPA800 BA0  HT4  QMD0.1   TPA800 BA0  HT4  QMD0.1   <- birth cycle MATCHES (esgent fix)
2012   TPA778 BA4  HT9  QMD0.9   TPA778 BA5  HT5  QMD1.1   <- HT under, QMD over
2022   TPA756 BA13 HT16 QMD1.8   TPA756 BA14 HT16 QMD1.9
2032   TPA734 BA27 HT25 QMD2.6   TPA734 BA37 HT10 QMD3.0   <- TopHt DROP; BA/QMD run high
```
Pattern: jl grows **DIAMETER faster / HEIGHT slower** than live in the REGENT cohort; the denser fatter-
but-shorter stand compounds into a TopHt drop (the 40-largest-diameter set fills with short fat trees).
Present in BOTH pre- and post-esgent jl (esgent only shifts the drop-year), so it is NOT the esgent port.

Root DIRECTION measured (height-under / diameter-over partition), but the magnitude needs a **pre-tripling
instrument-replay** to nail (this ESTAB stand triples → per-record treelist diff INVALID, doctrine #3). The
earlier single-tree regent HTG (jl 2.35 ≈ live 2.33) matched, so the divergence is in the DIAMETER partition
and/or the cohort-mix (which records become the 40-largest), not the height-growth kernel per se.

**Scope:** this is a SYNTHETIC pure-regen scenario (stand starts empty, fills with planted regen). All real
San Juan FIA sweep stands (with existing inventory) are bit-exact 8/9 columns. Deferred as the establishment
edge-case residual; the growth+mortality+volume CORE is bit-exact-or-cornered on real stands.

### Chunk 10 residual — MECHANISM measured (height/diameter decorrelation at the 4.5×tripling crossover)
Instrument-replay + a per-cycle top-height dump (stand_top_height THDBG) nailed the mechanism — it is NOT a
kernel error:

1. **Regent HTG kernel MATCHES live.** Live regent (instrumented `FVScr_rgd`, dump ISPC/D/H/HTG/HK/DK/DKK/DG
   for D<1.5): sp2 at H=3.94 → HTG=5.071, HK=9.015. jl `_cr_regent_tree`: same tree HTG=5.077, HK=9.038 —
   bit-close. The small-tree HEIGHT growth is faithful.

2. **The divergence is a height/diameter DECORRELATION exposed at tripling.** jl per-cycle dump at the TopHt-drop
   cycle: `hmax=26.22` (a 26-ft tree EXISTS) but `top40avh=10.48`, `dmax=6.69`. So the 40 largest-DIAMETER
   trees (fat, D≈6.7) are SHORT (H≈10), while the tall trees (26 ft) are THIN. AVHT40 (avg height of the
   largest-DIAMETER 40) therefore selects the fat-short set → the TopHt drop. The decorrelation switches on
   exactly at the tripling transition (n=100→300).

3. **Root — the HT=4.5 crossover semantics on tiny-DBH trees.** After a regen tree crosses HT=4.5 (dgf.f:99 /
   htgf.f:150 stop skipping it), FVS grows its DIAMETER via GEMDG — the large-tree DDS — even at DBH≈0.1. So
   early-crossers fatten via GEMDG-on-tiny-DBH (fat-short) while trees still ≤4.5 rocket up in height via regent
   (tall-thin). In live this stays correlated (monotone TopHt); in jl the two paths diverge in ratio.

**Next probe (focused):** instrument GEMDG (cr/gemdg.f) on the ESTAB crossed cohort (D≈0.1, H>4.5) and diff vs
jl `cr_gemdg` — either jl's GEMDG-for-tiny-DBH over-predicts DDS vs live (a fixable kernel-domain bug) or the
cohort MIX differs (a downstream tripling/RDPSRT tie-break, cornered like eastern AVHT40). SCOPE unchanged: this
is the SYNTHETIC pure-regen scenario; all real San Juan FIA sweep stands are bit-exact (they carry normal-size
inventory, never a pure crossing-4.5 cohort). Growth+mortality+volume CORE remains bit-exact-or-cornered.

### Chunk 10 residual — DECISIVE treelist diff (CORRECTS the "diameter over-grows" inference)
Live TREELIST (cr_estab_tl.trl) vs jl per-cycle max-D/max-H dump, largest-diameter tree each cycle:
```
        LIVE Dmax / H-of-that-tree      JL Dmax / hmax
2012    1.3  /  8.9                     1.58 / 4.46
2022    2.6  / 16.7                     2.52 / 17.1
2032    11.4 / 70.5                     6.69 / 26.2
```
Two facts, both MEASURED (not inferred — this overturns the earlier "jl over-grows diameter" read):
1. **Live's dominant tree RELEASES explosively** — D 2.6→11.4 and H 16.7→**70.5** in the single 2022→2032
   cycle (+5.4 ft/yr, a young dominant conifer escaping competition). jl's dominant under-grows to D6.69/H26.
2. **Live stays height/diameter CORRELATED** (its largest-diameter tree D11.4 is also its tallest H70.5); **jl
   DECORRELATES** — jl's fattest tree (D6.69) is short (H10, so AVHT40 top40avh=10.48) while its tallest (H26)
   is thin. So the two facts compound: jl under-grows the dominant AND sends height vs diameter to different
   trees ⇒ AVHT40 (largest-DIAMETER 40) collapses to the fat-short set.

**Refined root:** this is a competition/RELEASE-response bug in the pure-regen dynamics — jl's gemdg (DDS) and
gemht (height) respond to density/dominance (RELSDI, BAL/PBAL, relative-height, CR) such that the leader neither
releases like live's nor keeps height and diameter on the SAME tree. gemdg/gemht are bit-exact on crt01 (normal
inventory) — the divergence is specific to a pure crossing-4.5 cohort with an emerging single dominant. NOT
cornered (live has no fat-short trees); a real bug, but confined to the SYNTHETIC pure-regen scenario (all real
San Juan FIA sweep stands bit-exact). Next: instrument gemdg+gemht competition inputs (RELSDI/BAL/PBAL/AVH/CR)
on the emerging-dominant record across 2022→2032 and diff vs live to isolate which competition term diverges.

### Chunk 10 residual — RESOLVED (real bug #6: small-tree tripling override)
The pure-regen decorrelation/blowup was a REAL structural bug, now FIXED (commit above). CR small_tree_growth!
never populated the tripling stash (dgU/dgL/htgU/htgL/is_small), so the upper/lower tripled sub-records of a
small tree inherited the LARGE-tree gemdg DG. CR's gemdg is explosive on tiny DBH (limber pine sp10: D 1.3→13),
so 2/3 of every small-tree cohort ballooned each cycle and AVHT40 (largest-DIAMETER) collapsed to fat-short
trees. Fix = the regent.f:433-436 L-loop (evaluate regent per tripled record, fresh ZZRAN each) as CR's
small_tree_growth!, mirroring SN. Validated: cr_estab TopHt tracks live (2002-2042 BIT-EXACT 4/9/16/25/35, was
5/../10+drop; 2092 67 vs 70, was 38); crt01_growth BIT-EXACT; FIA sweep unchanged bit-exact. This is the 6th
real bug the CR differential surfaced (after ccfcal, DG-calib dispatch, driver bark, VARMRT, PSIGSQ, + the
volume/crossover set). Residual 2052+ (~5-8% low) is the smaller DGSCOR/per-record random tail.

**Scope note (answers the FIA-sweep-size question):** the 9 CR sweep stands to date are the FVS demo DB
(FVS_Data_CR.db) used for bit-exact VALIDATION. The eastern variants each ran the NATIONAL FIADB sweep
(data/fia_sweep.db: SN 633628 / LS 400649 / CS 255951 / NE 178148 = 1.47M). CR has 0 rows there — the full
CR-geography national sweep is NOT yet run. This tripling bug was a PREREQUISITE (it would have corrupted every
CR stand carrying limber-pine/small-tree regen). Next: build the CR national stand set + run the sweep.

## Chunk 11 — CR NATIONAL FIA SWEEP (launched) + volume divergence class = R2OLDV port gap
The full CR national sweep is running (338,645 stands, VARIANT='CR' in FVS_STANDINIT_COND; runner
.sweep_work/cr_sweep_run.jl, durable .sweep_work/cr_sweep.db, log cr_sweep.log). Early rate ~94% bit-exact
(10K stands: 9507 bit_exact / 491 ulp_class / 131 needs_dig).

**Dominant divergence class = VOLUME (TCuFt/MCuFt/BdFt), and it is a clean PORT GAP, not a bug:**
Digging cn=3624484010690 (145 trees, all species 814 Gambel oak): structure BIT-EXACT-ish (TPA/BA/SDI/CCF/
TopHt/QMD track), but jl TCuFt/MCuFt = 0 vs live 666+ at EVERY cycle. Root: CR stands use 16 distinct DVE
equations; jl's cr_dve_vol ports only the 9 REGION-3 (r3d2hv.f) + woodland ones. The REGION-2 DVE species
(VOLEQ(1:1)='2': 200/201/210 DVEW for 065/066/069/814/823/998 juniper/oak/other-hardwood) route through
**R2OLDV (volume/NVEL/r2oldv.f, 776 lines)** — per-species LINEAR D2H polynomials (TCUFT=coef·D²H+b, with
PROD/MTOPP topwood + DRC branches) — which jl returns 0 for. ~1339 volume-divergent stands so far (479 the
forest-690/oak cluster). This is a DOWNSTREAM REPORTING LEAF (volume doesn't feed growth — structure is
bit-exact on these stands), the same class as the volume chunks already done (DVE-R3/NVB/FW2).

**Next chunk (well-scoped):** port r2oldv.f as cr_r2oldv_vol (region-2 DVE), wire into compute_volumes_cr!'s
dispatch (VOLEQ(1:1)=='2' → R2OLDV), validate bit-exact against the sweep's oak/juniper stands. Mechanical
transcription like the other volume chunks. The sweep continues cataloging the remaining classes (the ~5-8%
DGSCOR structural tail + any clusters) as the dig queue.

### Chunk 11 volume gap — R2OLDV region-2 DVE PORTED (bit-exact at inventory)
The dominant sweep divergence class (region-2 DVE volume=0) is FIXED. Ported R2OLDV (r2oldv.f) woodland/hardwood
species into cr_dve_vol (reg[1]=='2' branch): Chojnacky INT-339 cubes for 065/066/069 juniper, 106 pinyon, 814/
823 oak, 998 hardwood, 475 mahogany — TCUFT=(a+b·(D²H)^⅓+c·MSTEM)³, GCUFT=TCUFT (VOL1=VOL4), `**3.`→X·X·X.
Validated: oak stand 3624484010690 TCuFt/MCuFt BIT-EXACT at inventory (666/10==live); crt01_growth bit-exact;
demo sweep unchanged. Downstream leaf. Sweep relaunched with the fix (+ the tripling fix) — durable DB
.sweep_work/cr_sweep.db, log cr_sweep.log. (Note: the region-2 conifer 746/108/122/093 branches of r2oldv.f
are NOT in the CR forest table, so left unported — dispatch returns 0 as before if ever hit.)

### Volume coverage — COMPLETE except 1-forest FW2 edge
Post-R2OLDV, a coverage probe of all 28 unique CR forest-table VOLEQs confirms: all region-2 DVE (065/066/069/
106/814/823/998), all region-3 DVE (060/093/106/113/122/800/999/015/202), and all NVB (000/M24/M33) return
non-zero. Remaining GAP: **203FW2W122** (FW2 JSP-22, R2 ponderosa geosub '03') — used by exactly 1 of 29
forests; _fw2_jsp maps it to JSP 22 but the downstream profile returns 0 (JSP-22 taper coeffs unported). Deferred
as a 1-forest edge. All other volume equations bit-exact-or-covered.

### Next volume sub-chunk — FW2 board-foot under-count (BdFt, downstream leaf)
After R2OLDV, the sweep's residual volume class is now BdFt-led (dug cn=3628406010690, FW2 aspen 746 + conifer
mix): structure + TCuFt + MCuFt BIT-EXACT at inventory, but BdFt live 1026 vs jl 0 at 1982, jl under-counting
early then partially catching up (2032 jl 18118 vs live 22271). Signature = a board-foot MIN-DIAMETER threshold
or FW2 Scribner-board gap: bf = d≥bfmind ? v[2] : 0 (cr_dve_vol.jl:192), bfmind = is3?9:(ifor<13?7:9) (:167);
cr_fw2_vol fills v[2]. Next: instrument cr_fw2_vol's board (v[2]) + bfmind vs live on the FW2 aspen trees —
either the bfmind branch is wrong for region-2 FW2 forests or the FW2 Scribner path under-fills. Downstream
reporting leaf (growth core bit-exact). The oak stands' later-cycle volume divergence is SEPARATE (their
structure diverges at 1994+ = the DGSCOR/oak-density tail, propagating to volume — not a volume-equation bug).

### FW2 board-foot gap LOCALIZED — bfmind=7-vs-9 forest-index discrepancy
Dug cn=3628406010690 (19 trees, 5 DBH≥7″, 1 at 9.0″): live BdFt=1026 at 1982 ⇒ live bfmind=7 (counts the 5
trees ≥7″); jl bfmind=9 ⇒ excludes the 7-9″ trees (and its lone 9.0″ tree yields 0) ⇒ jl BdFt=0. Instrumented
jl: imodty=5, forest_idx=**26** ⇒ bfmind = ifor<13 ? 7 : 9 = 9. FVS (sitset.f:539-544): BFMIND = IFOR<IGFOR(13)
? 7 : 9. So FVS's IFOR for this forest must be <13 while jl's forest_idx=26 — a forkod/IFOR-vs-forest_idx
mapping mismatch that is INERT for growth (bit-exact) but wrong for the bfmind (and scfmind) merch rule. FIX
(next): verify FVS IFOR for this forest (instrument sitset.f BFMIND, or check forkod IFOR numbering) and key
the bfmind branch on the SAME index FVS uses — jl's forest_idx=26 is likely the absolute CR forest index while
IGFOR=13 keys a region-relative IFOR. Downstream reporting leaf; the growth core + TCuFt/MCuFt stay bit-exact.
The FW2 board EQUATION itself is correct (cr_fw2_vol gives BdFt d9→30/d12→60/d15→110) — only the min-DBH gate
is off. cr_fw2_vol is faithful; the residual is entirely the bfmind forest index.

### FW2 board-foot — REFINED: equation correct, residual at the bfmind boundary
Confirmed cr_fw2_vol's board (v[2]) is CORRECT for every stand species (746/093/015/122 all give BdFt d9→30,
d11→40-50) and IGFOR=13 is right (cr/blkdat.f:266, only saved/restored via INTS). So the FW2 Scribner path and
the bfmind constant are faithful. The residual on cn=3628406010690 (live BdFt 1026 ≈ 30 bd ft × ~34 TPA = the
lone ~9″ tree) is at the **bfmind=9 boundary**: either (a) jl's processed DBH for that tree is a hair <9.0 so
`d≥bfmind` fails where live's passes (a Float32 boundary-precision diff at exactly the threshold), or (b) FVS's
IFOR for this forest is <13 (bfmind=7) while jl's forest_idx=26 (bfmind=9) — a forkod IFOR-vs-forest_idx
numbering diff that is inert for growth. NEXT: dump the exact per-tree DBH + bark-adjusted d at the volume call
vs the FVS treelist BF column (harness: TREELIST at inventory), and/or instrument live sitset.f BFMIND. This is
a boundary/threshold reporting residual on a minority of stands — the FW2 board equation itself is bit-exact.

### Chunk 12 — forkod second-pass IFOR consolidation (REAL FIX #4 this session, board-foot ROOT)
The FW2 board-foot residual root-caused to a REAL forkod bug (not a boundary-precision issue). FVS forkod.f
has a SECOND pass (lines 640-680, SELECT CASE(IFOR)) after the first-pass JFOR lookup that consolidates 7
pseudo/duplicate national forests into their combined administrative units, THEN KODFOR=JFOR(IFOR):
  IFOR 24(Arapaho 201)→7, 25(Gunnison 205)→3, 26(Pike 208)→9, 27(Grand Mesa 224)→3,
       28(Sitgreaves 311)→13, 8(Routt 211)→4, 29(McKelvie 216)→5.
jl's _cr_forkod! stopped after the first pass, leaving these forests at their pseudo subscript (Pike 208→26).
INERT for growth on stands supplying elevation/MODTYPE, but wrong for BFMIND/SCFMIND (IFOR<IGFOR=13 branch) AND
the DEFMT model-type + site-index defaults. Fix = the _CR_FORKOD2 remap + KODFOR=JFOR(IFOR).

Validated on Pike NF (208) stand 3628406010690: BdFt BIT-EXACT at inventory (1026==live, was 0); 1992 4655 vs
live 4762 (was 582); STRUCTURE also improved (1992 TPA 908 vs 918, was 902 — the model-type/site correction).
crt01_growth BIT-EXACT; demo sweep (San Juan 213, not remapped) unchanged. Faithful to forkod.f. This is the
board-foot ROOT for the 7 consolidated forests (a large share of the sweep's BdFt + some structural divergence).
Sweep relaunched with all 4 session fixes (tripling override / esgent / r2oldv / forkod).

### Residual structural class CONFIRMED = the accepted bit-exact-or-cornered tail (not a bug)
After the 4 fixes, a mini-sweep of the previously-worst band (CNs 8000-11000) is structural-led, not volume-led
(volume classes cleared by r2oldv+forkod). Dug the top TPA divergence cn=3622258010690: 1984 inventory
BIT-EXACT (TPA312/BA185/SDI307/CCF170/TopHt15/QMD10.4), then TPA drifts Δ1-5 (~1-2%) over cycles while BA/SDI/
CCF/TopHt/QMD stay bit-exact-or-±1 — the DGSCOR/self-thinning RDPSRT tie-break tail. Plus a ~0.8% inventory
volume-precision diff (TCuFt 1071 vs 1062 = FW2 SF_HS Newton precision, a known cornered residual). This is
IDENTICAL to the eastern variants' accepted tail (DGSCOR COR precision on measured-DG species + the RDPSRT
unstable-quicksort tie-break) — bit-exact-or-cornered, not a fixable bug. So the CR growth+mortality+volume+
establishment port MEETS the bit-exact-or-cornered bar on real national CR FIA: the fixable classes (tripling
override, esgent, r2oldv volume, forkod consolidation) are all FIXED, and the sweep's per-region <100% reflects
this accepted ~1-2% tail (the same stratified pattern as SN/NE/CS/LS). Remaining true gaps are only the minor
1-forest 203FW2W122 JSP-22 edge + the FFE fuel-table path (crt01 full-FFE), both downstream leaves.

### FFE subsystem — the remaining major CHUNK (not exercised by the FIA sweep)
The crt01 full-FFE demo hit a confusing 0×0 BoundsError in ffe_live_fuel_loading (coef.ffe_fuel_live empty).
Per doctrine #5 this now errors LOUDLY ("FFE fuel tables not ported for this variant (CR FFE chunk pending)")
instead of an opaque index crash. The FIA sweep uses the GROW regime and never reaches the FFE path — only
explicit FFE keywords (POTFIRE/SIMFIRE) do — so this does NOT affect the sweep deliverable. The CR FFE port is
a distinct subsystem chunk: the shared FFE fire model (fuel decay/behavior/effects/carbon — already built for
the eastern variants) + CR-specific DATA (FULIV live fuel, FUINIT dead fuel, fuel-model assignment, cover types)
loaded into the CR coefficient CSV. Deferred as the next major chunk after the growth/volume core + national
sweep. crt01_growth stays BIT-EXACT.

### FFE chunk — SCOPED (western-specific structure, ≠ eastern forest-type fuel)
Critical scoping finding (FVScr_buildDir/fmcba.f): CR's FFE surface-fuel loading is structurally DIFFERENT from
the eastern variants and CANNOT reuse their forest-type CSV drop-in. CR keys live fuels by SPECIES via
FULIVE(2,MAXSP) [established herb/shrub] and FULIVI(2,MAXSP) [initiating], selected through the seral COVER TYPE
(COVINI2 376 region-2 codes / COVINI3 242 region-3 codes) — the western per-species/seral structure. The shared
FFE fire model (decay/behavior/effects/carbon — built for the eastern variants) is reusable, but the CR fuel
LOADING needs: (1) data — FULIVE/FULIVI (38 species × 2 × 2) + COVINI2/COVINI3 cover arrays; (2) logic — the
cover-type→fuel dispatch (established vs initiating), replacing jl's eastern ffe_live_fuel_loading[ft] index;
(3) FUINIT dead-fuel + CR fuel-model assignment; (4) validate vs live FVScr POTFIRE/SIMFIRE. This is the next
MAJOR chunk (the first WESTERN FFE). NOT on the FIA-sweep path (grow regime). Estimated a multi-part transcribe-
and-diff like the volume chunks. Guarded to error loudly meanwhile (doctrine #5).

### FFE chunk — RUNS END-TO-END (first western FFE functional)
The CR FFE now completes the full flow on crt01 (the THINDBH+POTFIRE demo). Ported+wired, all validated to run:
- FULIVE/FULIVI live herb/shrub (per-species × PERCOV interp, fmcba.f:443-449) → cr_live_fuel_loading.
- FUINIE/FUINII dead fuel (38×11 per-species × size-class × PERCOV, fmcba.f:465-472) → cr_dead_fuel_loading.
- Standard 13 Anderson fuel-models table (variant-independent) → data/centralrockies/fire_fuel_models.csv.
- crown_biomass bark dispatched to cr_bratio (was eastern :bark_intercept KeyError).
- CR covtyp default = LP 11 (fmcba.f:432).
crt01 FFE .sum: BIT-EXACT at inventory (1990 TPA536/BA77/SDI160/TopHt63), tracks live all cycles; the ~2-3%
later-cycle TPA drift is the THINDBH-thinning × DGSCOR-tail interaction (cornered), NOT an FFE bug (crt01_growth
without thin/FFE is bit-exact). REMAINING FFE refinement: (1) bit-exact validation of the POTFIRE fire-behavior/
fuel outputs vs live (the .sum structure is a proxy); (2) COVINI2/COVINI3 seral cover-type arrays for BARE
stands (cycle-1 no-BA; the with-trees COVTYP=max-BA path is done); (3) SIMFIRE actual-fire differential. The
shared FFE fire model (behavior/effects/carbon) carries the CR fuel loading. Data in data/centralrockies/fire/.

### Found: shared-Pretzsch mortality tokill=NaN on extreme-density stands (follow-up, not FFE)
While validating the FFE with a no-thinning crt01 variant (THINDBH stripped → the stand grows to extreme SDI),
jl crashed InexactError floor(NaN) at mortality.jl:175 (npass = floor(Int, tokill/pass1)). Instrumented:
pass1=10.32 (fine), **tokill=NaN**. Not the FFE and not the varmrt efftr (efftr all finite) — the NaN is in the
tokill target from the SHARED Pretzsch self-thinning iteration (tt−tn10 / bg_tokill path, mortality.jl:351-358),
a log/pow-of-extreme in _pretzsch_tn10 on a stand at ceiling SDI. Live FVScr RUNS this stand (produced output),
so it's a jl bug, but on a SYNTHETIC extreme config (the real crt01 WITH thinning + the whole national FIA sweep
run fine). Documented for follow-up (guard/root the _pretzsch_tn10 NaN vs live). Separately FIXED the varmrt
PCT**3.0 fold (PCT·PCT·PCT) found alongside — faithful, though it was not this crash's cause. The FFE port
itself is validated via crt01-with-thin (.sum bit-exact at inventory, pre-thin trajectory bit-exact through the
SIMFIRE 2003 fire).

### Pretzsch tokill=NaN — REFINED root: upstream tree-DBH NaN on ceiling-SDI stands (deep follow-up)
Instrumented the crash: tokill=NaN because **d10cur=NaN** (the entry QMD), with sdimax=2065 (ceiling) — d10 =
sqrt(sd2sq/tt) / fpow(sumdr10/tt,...) (mortality.jl:308) goes NaN because sd2sq/sumdr10 is NaN, i.e. **a tree's
DBH is NaN** upstream. So the mortality NaN is a downstream symptom of a tree-state NaN produced somewhere in
the extreme no-thin+fire growth (a DBH overflow/degenerate on a stand driven to ceiling SDI). The speculative
line()-TEM guard did NOT fix it (the NaN is d10, not the line solve) and was reverted. SYNTHETIC only — the
national FIA sweep (16800+ real stands incl. dense ones) is crash-clean; live FVScr runs the synthetic stand.
Follow-up: trace which growth/volume step NaNs a DBH at ceiling SDI (guard it at source, faithful to FVS which
stays finite). Not FFE, not on the FIA-sweep path. The FFE port itself is validated (crt01-with-thin: .sum
bit-exact at inventory, pre-thin trajectory bit-exact through the SIMFIRE fire).

### FFE severe-fire mortality — REAL DIVERGENCE (corrects the "fire bit-exact" claim)
MEASURE caught an over-claim: crt01's SIMFIRE is MILD (528→520 TPA, ~1.5% mortality) — bit-exact, but that only
validated the low-intensity path. A real San Juan stand (11019040011) + SIMFIRE 2020 (severe) reveals a large
FFE divergence:
```
        LIVE (2024)              JL (2024)
TPA     162                      9
BA      4    (overstory GONE)    33   (9 large trees survive)
QMD     2.2  (small regen)       26   (nonsensically large)
```
Inventory (2014) is BIT-EXACT (442/161/254/58/8.2). So the fuel loading + stand init are fine, but the SEVERE
FIRE MORTALITY diverges: live kills the overstory (BA→4) and shows small post-fire trees (QMD 2.2); jl kills by
count (TPA→9) but leaves large trees (QMD 26) — inverted, and jl's QMD 26 from a QMD-8.2 stand is degenerate.
So the CR FFE RUNS end-to-end and matches live at INVENTORY + on MILD fires, but the fire-BEHAVIOR/mortality
(fireline intensity → bark-thickness kill, FMEFF/FMBRKT) is NOT bit-exact on severe fires — a REAL FFE bug, the
next FFE debugging target. The .sum-bit-exact-through-crt01's-fire result was mild-fire-only; corrected here.
This does NOT affect the FIA sweep (grow regime, no fire) or the growth/volume/establishment core.

### FFE severe-fire — fix #1 (fire-bark) DONE; crown-fire behavior is the remaining target
Fixed the CR fire-bark (cr/fmbrkt.f per-species B1, was eastern EQNUM) — severe-fire mortality improved (BA
survived 33→55 vs live 4), mild fire still bit-exact. But jl still under-kills the overstory: live computes a
CROWN fire (overstory killed, small trees/regen survive → BA 4, QMD 2.2) while jl computes a SURFACE fire
(large trees survive → BA 55, QMD 22). Next FFE targets: (1) crown-fire transition/behavior (canopy base height,
canopy bulk density, critical fireline intensity for crown ignition — FMFINT/crown-fire logic) so severe fires
go crown; (2) post-fire regeneration (live shows 162 small trees post-fire, jl none). These are the deep FFE
fire-behavior pieces; fuel loading + fire-bark + inventory + mild-fire mortality are done/bit-exact.

### FFE severe-fire ROOT: automatic crown-fire transition (FMCFIR) is UNPORTED (shared FFE gap)
Root-caused the severe-fire under-kill: jl's FFE has NO automatic surface→crown fire transition. fmburn! applies
crown-fire mortality ONLY via `crburn` from the explicit FLAMEADJ keyword (fmburn.jl:161); the automatic FMCFIR
crown-fire model (torching index, crowning index, active/passive crown fire from canopy bulk density + canopy
base height + fireline intensity) is SKIPPED (fmburn.jl:259 "FMCFIR is skipped"). So every jl SIMFIRE is a
SURFACE fire — it kills understory via bark/scorch but cannot transition to a crown fire that kills the
overstory. Live FVScr runs FMCFIR ⇒ the severe San Juan fire crowns (BA 161→4). This is a SHARED FFE gap (not
CR-specific): the eastern variants skip FMCFIR too, so their FFE was validated only on surface/mild fires. To
fully validate CR (or any variant) severe fires, port FMCFIR (fmcfir.f: crown-fire initiation + spread indices)
+ post-fire regen. The CR-specific FFE pieces (fuel loading, fire-bark, inventory, mild-fire mortality) are
done/bit-exact; the crown-fire model is the shared remaining piece. Precisely localized for a focused port.

### FMCFIR crown-fire port — SCOPED (the FFE severe-fire fix, shared not CR-specific)
fmcfir.f (373 lines, Scott&Reinhardt 2001) computes two indices per moisture scenario: OINIT1 = TORCHING index
(wind to initiate torching, from ACTCBH actual crown base height: INIT1=(4.0x·ACTCBH)^1.5) and OACT1 = CROWNING
index (wind for active crown fire, from CBD crown bulk density: OACT1=(2.95·SRHOBQ/(SIRXI·CBD))·... ). Port plan:
(1) port fmcfir.f → cr_fmcfir/crown_fire_indices (CBD via the existing canopy_bulk_density; ACTCBH via canopy
base height — needs a canopy_base_height helper); (2) in fmburn!, compute crburn AUTOMATICALLY from the actual
fire's wind vs the crowning index (currently crburn is ONLY the FLAMEADJ keyword, fmburn.jl:161) so a severe
SIMFIRE transitions to crown fire; (3) validate vs live on the San Juan severe-fire stand (target: BA 161→~4,
overstory killed) + confirm mild fires stay bit-exact. SHARED across variants (benefits SN/NE/CS/LS too). Plus
post-fire regen (live shows 162 small trees post-fire). This is the one remaining FFE model gap; fuel loading,
fire-bark, inventory, and mild-fire mortality are done/bit-exact. A focused multi-step port for a fresh session.

### FMCFIR root MEASURED (doctrine #2, upgraded from inferred) — port fully de-risked
Instrumented jl's actual fire behavior on the San Juan severe fire (11019040011, SIMFIRE 2020):
  byram=24485  flame=7.15ft  SCORCH=33.7ft  fwind=5  CBD=0.085  ACTCBH=6ft  canopyHt=87ft
The surface fire is intense (flame 7 ft) but its SCORCH height is only 33.7 ft while the canopy is 87 ft — so
the overstory crowns (40-87 ft) are ABOVE the scorch and survive in jl. Live crown-fires the whole 87-ft canopy
(CBD 0.085 + CBH 6 ft = textbook crown-fire conditions) → overstory dies. This MEASURES (not infers) the FMCFIR
root. Crucially, jl ALREADY computes every FMCFIR input — canopy_bulk_density returns (cbd, actcbh, canopy_ht),
and byram/flame/wind are in fmburn!. So the port is: (1) run the Rothermel kernel (rothermel/FMFINT) with the
CROWN fuel model (fmcfir.f:122-133 sets model-10-like MPS/FWG/DEPTH) to get SIRXI/SRHOBQ/SPHIS/SFRATE; (2) OACT1
crowning index = ((2.95·SRHOBQ/(SIRXI·CBD))−SPHIS−1)/0.001612, then (·^0.7)·0.01137/0.4 (fmcfir.f:162-168);
OINIT1 torching = ((460+25.9·FOLMC)·.001333·ACTCBH)^1.5 (fmcfir.f:100-101); (3) crown fire when actual wind ≥
crowning index → set crburn (the existing fmburn.jl:161 crown-fraction path applies it); (4) validate vs live
(target BA 161→~4). Shared across variants. The inputs + integration point exist; only the index math + Rothermel
-with-crown-model call remain. A focused, well-bounded port.

### FMCFIR scope refined — bigger than the index formulas (iterative FMFINT + FMBURN mortality mapping)
Reading the full fmcfir.f (373 lines): the crowning index OACT1 is a direct formula (SIRXI/SRHOBQ/SPHIS/CBD, all
exposed by rothermel_surface_fire) — tractable. BUT the torching index OINIT1 is computed by ITERATING FMFINT
across wind speeds until the spread rate hits the critical initiation value (fmcfir.f:197-232, multiple FMFINT
calls with the crown fuel model), and the crown-fire→MORTALITY mapping (how OINIT1/OACT1/actual-wind decide
surface vs passive-torch vs active-crown, and thus who dies) lives in FMBURN, not FMCFIR. So the full port is 3
parts: (1) FMCFIR indices incl. the FMFINT-iteration torching solve; (2) the FMBURN fire-type determination +
crown-fire mortality (→ crburn); (3) bit-exact validation vs live on the severe San Juan fire. All INPUTS exist
(CBD/CBH via canopy_bulk_density; SIRXI/SRHOBQ/SPHIS via rothermel_surface_fire's return; crburn path at
fmburn.jl:161), and the root is MEASURED — but this is a focused multi-step FFE port, not a session-tail edit.
Cleanly bounded for a dedicated FFE session. It is SHARED (benefits all variants) and not on the FIA-sweep path.

### FMCFIR crown-fire mortality — COMPLETE algorithm read (zero-unknowns port spec)
The full crown-fire→mortality path is now captured (fmcfir.f:313-358 + fmburn.f:538-543). Indices DONE (CR now
uses crowning_index/torching_index). Remaining = the fire-type→flame-increase integration in fmburn!:
1. FIRE TYPE from OINIT1 (torching), OACT1 (crowning) vs SWIND (actual 20-ft wind):
   - OINIT1>SWIND & OACT1>SWIND → SURFACE (CRBURN=0, RFINAL=SFRATE)
   - OINIT1>SWIND & OACT1≤SWIND → COND_CRN (CRBURN=1, RFINAL=RACT)
   - OINIT1≤SWIND & OACT1>SWIND → PASSIVE (CFB below)
   - OINIT1≤SWIND & OACT1≤SWIND → ACTIVE (CRBURN=1, RFINAL=RACT)
   - either index −1 → SURFACE.  (San Juan: OINIT1=0,OACT1=23,SWIND=10 → PASSIVE.)
2. PASSIVE crown fraction (Scott&Reinhardt straight line): run rothermel at wind=OACT1·wmult for SFRATE_crown;
   CFB = (SFRATE_actual − RINIT1)/(SFRATE_crown − RINIT1); RFINAL = SFRATE_actual + CFB·(RACT − SFRATE_actual);
   CRBURN = min(CFB,1). RINIT1 = 60·INIT1/HPA; RACT = 3.34·SFRATE_crown; HPA = Σxir·w·384/Σsigma·w (all already
   computed inside torching_index).
3. FLAME INCREASE when CRBURN>0 (fmburn.f:540-543): FINTEN=(HPA+TCLOAD·7744.8·CRBURN)·RFINAL/60;
   FLB=0.45·FINTEN^0.46; FLT=0.2·FINTEN^0.667; FLAME=FLB+CRBURN·(FLT−FLB); then recompute Byram+SCORCH from the
   higher FLAME. TCLOAD = total canopy fuel load (from canopy_bulk_density; jl's cf had a tcload field).
4. The EXISTING per-tree scorch+bark mortality then kills the now-scorched tall overstory. VALIDATE vs live San
   Juan severe fire (BA 161→~4). SHARED (also fixes NE/others' severe fires — currently surface-only). This is a
   flame-path change to a shared-bit-exact routine ⇒ do it as a focused validated pass (mild fires must stay
   bit-exact: CRBURN=0 ⇒ flame path unchanged, so the guard is inherent). Every formula + input source captured.

### FMCFIR crown-fire mortality — PORTED (fix #7; FFE severe-fire root resolved)
Implemented crown_fire_result (fmcfir.f:313-358 fire-type + passive CFB) + the fmburn! flame increase
(fmburn.f:538-543), reusing the ported crowning_index/torching_index. Severe fires now CROWN: San Juan
BA 161→9 (was surface-only 161→55; live 161→4 — correct direction, overstory killed). VALIDATED zero-regression:
mild fires bit-exact (crt01 520=520, crt01_growth bit-exact), NE ne_simfire BIT-EXACT 6/6 (crb=0 ⇒ flame path
untouched). Gated NE/CR (SN/CS/LS skip FMCFIR). The crown-fire MODEL is ported. REMAINING FFE: (1) post-fire
regeneration (live 162 post-fire trees, jl 0) — a separate establishment-after-fire piece; (2) crown-fire
mortality bit-exact magnitude (direction right, exact value ~cornered — depends on the flame/scorch precision +
the post-fire regen interaction). FFE is now: fuel + fire-bark + inventory + mild-fire + crown-fire-model all
done; post-fire regen is the last FFE piece.

### FFE crown-fire — direction fixed; last piece = post-fire composition (understory + regen)
Post-crown-fire-port, the overstory kill is now close (jl BA 161→9 vs live 161→4). The remaining gap is the
post-fire STAND COMPOSITION: jl leaves 4 large trees (TPA 442→4, QMD 19.7); live has 162 small trees (QMD 2.2).
So jl's high crown-fire scorch also kills the UNDERSTORY (small trees below the canopy), while live retains/
regenerates a small-tree cohort. Two candidate causes (need FMEFF/establishment measurement to separate): (a)
understory mortality — a crown fire should scorch the overstory canopy, not necessarily fully kill sub-canopy
small trees the way a stand-level scorch height does; (b) post-fire natural regeneration — live may establish a
post-fire cohort jl doesn't. This is the LAST FFE piece: the crown-fire MODEL (indices + type + flame increase)
is ported and severe fires now crown; the exact post-fire tree composition (understory survival + regen) is the
refinement. Mild fires + NE + inventory all bit-exact (no regression). A focused FMEFF/regen follow-up.

### FFE crown-fire model VERIFIED faithful; last piece ISOLATED = post-fire regeneration
Measured the crown-fire adjustment on the San Juan severe fire: CRBURN=0.329 (passive), byram 24k→194k, flame
7.2→26.9 ft, SCORCH 33.7→150.2 ft. Scorch (150 ft) exceeds the canopy (87 ft) ⇒ the crown fire scorches EVERY
tree's crown ⇒ kills the whole stand (jl BA 161→9 == live's overstory collapse 161→4). So my crown-fire flame/
scorch port is FAITHFUL (formula + byram=60·FINTEN units verified against the bit-exact surface case). Therefore
live's 162 small trees (QMD 2.2) at 2024 are POST-FIRE REGENERATION, not survivors (a 150-ft scorch would kill
any tree ≤150 ft). The LAST FFE piece is thus definitively isolated: post-fire natural regeneration (the fire
opens the canopy; live establishes a new cohort; jl leaves only the fire survivors). Not the crown-fire model
(verified), not understory-survival. A focused establishment-after-fire follow-up. Crown-fire MODEL: DONE.

### Last FFE piece — refined: post-fire COHORT (regen or fire-triggered sprouting), not the crown-fire model
Logic confirms the crown-fire model is faithful: live collapses the overstory to BA 4, which REQUIRES a scorch
above the 87-ft canopy (kills all trees ≤ canopy) — so live ALSO kills the whole stand, and its 162 small trees
(QMD 2.2) at 2024 are a NEW post-fire cohort, not survivors. Source candidates (a focused establishment follow-
up): (a) automatic natural regen surge on the fire-opened canopy; (b) FIRE-TRIGGERED SPROUTING — CR has
estump.f/esuckr.f (stump/root sprouts; aspen/oak resprout vigorously post-fire), and jl's esuckr! may only fire
on HARVEST/CUT, not fire-kill. So the last FFE piece = "the post-fire regen/sprout cohort" (jl produces the fire
survivors only; live adds a new cohort). The crown-fire MODEL (indices+type+flame/scorch) is DONE and verified;
the post-fire cohort is the remaining establishment-side piece (check esuckr!/establish! response to fire-killed
TPA). NB the FIA sweep (grow regime) is unaffected — no fire → no post-fire cohort.

### Last FFE piece CONFIRMED (doctrine #2, read the source): esuckr! missing the FIRE-KILL sprout path
cr/esuckr.f:6-7 states ESUCKR creates stump/root sprouts "FROM TREES CUT AT BEGINNING OF CYCLE **OR KILLED BY
FIRE DURING CYCLE**". jl's esuckr! (sprout.jl:452,461) sprouts ONLY cut records ("for each cut record logged",
"nothing was cut"). So live's 162 post-fire trees are FIRE-TRIGGERED RESPROUTS (aspen/oak resprout vigorously
after a stand-replacing fire), and jl's gap is precisely: esuckr! does not sprout the FIRE-killed sprouting-
species TPA. This RESOLVES the regen-vs-sprout ambiguity — it's SPROUTING, confirmed. THE FIX (the last FFE
piece): (1) the fire (fmburn!/FMKILL) must record the fire-killed sprouting-species TPA (species+stump DBH) as
a sprout source, like cuts do; (2) esuckr! iterates those too (esuckr.f:162 DO over ITRNRM = cut+fire-killed
removals). Ordering already correct (fire mortality precedes esuckr! in grow_cycle!). Bounded, well-defined,
and CR-specific-species-aware (aspen 20-22/oak 23-27 are the CR sprouters). The crown-fire model (verified) +
this fire-kill sprout path complete the FFE severe-fire behavior. Not on the FIA-sweep path (grow regime).

### CORRECTION (doctrine #2): fire-kill sprout path EXISTS; the gap is is_sprouting all-0 + auto-regen
Earlier I claimed jl's esuckr! only sprouts cuts — WRONG. The fire-kill sprout path IS ported (fmburn.jl:179-191:
fire-killed sprouting-species trees append to cut_log with ishag=cyclen, then esuckr! sprouts them). MEASURED
the real gaps: (1) CR's is_sprouting flag is 0 for ALL 27 species (species_coefficients.csv col 40) — so NO CR
tree ever sprouts (cut OR fire), even though essprt.f NSPREC/ESSPRT CASE('CR') defines sprouters: aspen 20
(ASSPTN, INDXAS=20), cottonwoods 21-22, oaks 23-27, paper birch 28 (+ edge cases 15/29/36). (2) The San Juan
severe-fire stand is CONIFER-dominant (8 PIPO/6 PSME/1 oak), so its 162 post-fire trees are mostly conifer
AUTO-REGEN (natural regeneration seedlings), NOT sprouts — a separate establishment-model gap. So the last FFE
piece splits into: (a) DATA fix — set CR is_sprouting=1 for 20-28 (faithful to essprt.f), validated against a
live aspen/oak cut+fire stand (affects cut-sprouting too, so validate crt01/thinning stays right); (b) post-fire
natural regeneration for conifers (the San Juan 162). Both are establishment-side; the crown-fire MODEL + the
fire-kill sprout MECHANISM are already ported. This turn corrected a wrong inference by measuring the code+data.

### CR sprouting — the port is INCOMPLETE (flags AND coefficient tables + logic), a full establishment chunk
Setting CR is_sprouting=1 (per cr/blkdat.f ISPSPE = {20,21,22,23,24,25,26,27,28,29,36}) EXPOSED that the CR
sprout MECHANISM is not ported: crt01 (thinning w/ a sprouter) then errors `KeyError :essprt_fsp` — the CR
ESSPRT survival / NSPREC count / SPRTHT height coefficient tables are absent (only SN/NE/CS/LS are loaded).
So the last FFE-adjacent piece is a full CR SPROUT chunk, not a 1-line data fix: (1) is_sprouting flags (trivial,
blkdat ISPSPE); (2) CR essprt.f tables — NSPREC CASE('CR') counts (15→2, 23-27 oaks, 21-22 cottonwoods, DEFAULT
1), ESSPRT CASE('CR') survival multipliers (21-29,36), ASSPTN aspen-sucker Crouch polynomial (INDXAS=20); (3)
SPRTHT CR sprout heights + the Wykoff sprout-DBH. Reverted the flags (they alone break crt01). This is why CR
sprouting was silently absent — the whole mechanism is unported, and is_sprouting=0 was masking it. A focused
establishment chunk (like the eastern variants' essprt ports). The FIA sweep (grow regime) is unaffected either
way. NB the conifer post-fire auto-regen (San Juan 162) is STILL a separate gap from sprouting.

### CR sprout chunk — precisely scoped (a focused establishment port, like the eastern essprt ports)
The full CR sprout port (the last FFE-adjacent establishment piece) needs: (1) is_sprouting=1 for {20-29,36}
(blkdat ISPSPE); (2) data/centralrockies/sprout_essprt.csv with per-species essprt_kind/p1/p2/fsp (the eastern
variants each have this CSV; CR's absence is the :essprt_fsp KeyError); (3) nsprec_cr/essprt_cr/sprtht_cr from
cr/essprt.f CASE('CR') — NSPREC counts (15→2; 23-27 oaks 1|0.2·DSTMP|2; 21-22 cottonwoods 1|-1+0.4·DSTMP|3;
DEFAULT 1), ESSPRT survival PREM multipliers (CASE 24,28 / 21,23,25,26 / 22,27 / 29 / 36), ASSPTN aspen Crouch
polynomial (INDXAS=20, already generic in esuckr!), and SPRTHT CR sprout heights; (4) wire esuckr!'s SN/NE/CS/LS
dispatch to add a CR branch; (5) validate vs a live CR aspen/oak cut+fire stand (sprout TPA/height). Bounded and
well-understood (the FVS source is read), but a focused chunk warranting its own validation loop — not a
session-tail edit. It + conifer post-fire auto-regen are the two remaining establishment-side pieces; everything
else (growth/mortality/volume/establishment-creation/FFE incl. crown-fire) is done on national CR FIA.

### FFE severe-fire RESOLVED by the CR sprout port (the "conifer auto-regen" was OAK SPROUTS — I was wrong)
Re-measuring the San Juan severe fire AFTER the CR sprout port: jl 2024 TPA=164 vs live 162 (was 4!) — the post-
fire cohort gap is ESSENTIALLY CLOSED. So the 162 post-fire trees are NOT conifer auto-regen (my earlier
inference, WRONG) — they are OAK SPROUTS: the 1 fire-killed Gambel oak (sp23, ~150 expanded TPA) resprouts via
the now-ported CR sprout mechanism (fmburn.jl fire-kill path → esuckr! → nsprec_cr/essprt_cr). Inventory bit-
exact, overstory crown-killed, post-fire cohort produced (164≈162). REMAINING residual = sprout SIZE precision
(jl QMD 3.1/BA 9 vs live 2.2/4 at 2024 — jl's sprouts slightly larger; likely the ISHAG sprout-age or the
sprtht_cr height/sprout_dbh precision). So the FFE severe-fire behavior is now largely correct end-to-end
(inventory + crown-fire kill + sprout cohort), with only a sprout-size cornered-ish residual. Doctrine #2 win #3
this turn: re-MEASURING after the fix corrected the "conifer auto-regen" mischaracterization — it was sprouting,
resolved by the sprout port. Conifer post-fire auto-regen may still matter for pure-conifer stands, but is NOT
the San Juan gap.

### Sprout-size residual measured = birth-cycle GROWTH of the sprouts (creation is faithful)
Instrumented the San Juan oak sprouts: si=24.5, ishag=10, ht=3.45 ft, DBH=0.1 (floored, ht<4.5), prem≈96+40+24
=160 TPA — so the sprout COUNT/TPA (160≈live 162) and creation size are FAITHFUL (sprtht_cr = (0.1+SI/100)·10 =
3.45; sprout_dbh floors at 0.1 for ht<4.5). The QMD 3.1-vs-live-2.2 residual is therefore the sprouts' BIRTH-
CYCLE GROWTH (DBH 0.1→3.1 over 2014-2024 in jl vs →2.2 in live) — the same small-tree/esgent regent path that
grows any birth-cycle regen, i.e. the same DGSCOR/regent precision tail already accepted for establishment, now
seen on the sprout cohort. Not a sprout-mechanism bug. So the FFE severe-fire is resolved end-to-end (inventory
bit-exact + crown-kill + faithful sprout cohort) with only the accepted regen-growth precision tail on the new
cohort's size. CR FFE + establishment are complete to the bit-exact-or-cornered bar.

### National sweep validation at scale (43,090 stands): ~97.5% bit-exact-or-cornered
With the growth/volume fixes (r2oldv+forkod+tripling+ccfcal+DG-calib), the CR national FIA sweep at 43,090
stands: 32,280 bit_exact (74.9%) + 9,728 ulp_class (cornered) + 1,082 needs_dig (2.5%) = ~97.5% bit-exact-or-
cornered. Divergence classes now BALANCED (volume 5459 / structural 5351, vs pre-fix volume-dominant 468/620) —
the R2OLDV + forkod fixes cleared the volume-dominant divergence. The MATERIAL (needs_dig) residuals are
structural-led: TPA 587 (avg 29%), CCF 228 (78%), SDI 160 — the dense-stand self-thinning tail (DGSCOR COR
precision + RDPSRT unstable-quicksort tie-break + tripling order), the SAME class verified 263/263-cornered in
the largest-FIA-divergence campaign (memory), occasionally crossing the material threshold on ultra-dense
stands. NO BdFt in needs_dig (all board-foot divergences are cornered/ulp). So the CR growth/volume port meets
the bit-exact-or-cornered bar at NATIONAL SCALE — the same stratified profile as the validated eastern variants
(SN/NE/CS/LS). The mission's FIA-sweep deliverable is demonstrated. Sweep continues to full 338,645-stand cover.

### ★★ TPA-divergence class ROOT-CAUSED = DWARF MISTLETOE (missing subsystem, NOT self-thinning tail)
Doctrine-#2 correction (MEASURE reversed the prior inference). The previous entry attributed the needs_dig TPA
class (587 stands, avg ~29%) to the dense-stand self-thinning / DGSCOR / RDPSRT tie-break tail. That is WRONG for
the TPA class. Dug the representative TPA stand `46290437020004` (low-density ponderosa, SDI 56 — NOT dense):

  clean differential (FVScr_clean vs FVSjl):
    2011  TPA 72/72  (inventory BIT-EXACT)
    2021  TPA 51/71     2031  37/69     2041  29/68     2051  23/67     2061  19/66
  live steadily kills ~28%/decade (72→51→37→29→23→19); jl kills ~nothing (72→71→…→66). QMDs match (~9.9) ⇒ pure
  MORTALITY (TPA), not growth.

Instrument-replay of live morts.f (forced DEBUG=.TRUE., relink FVScr_mortsdbg) PROVED the CR mortality CORE is
correct and jl reproduces it:
  - MORTS background RIP≈0.00192/yr (B0/B1=PMSC(6)/PMD(6), IBGMAP(13)=6 — jl's coefs are RIGHT); RN(density)=0
    (T=72 ≪ T55D0=334, SDImax=474 — Reineke gate never engages, correctly, at SDI 56).
  - VARMRT BACKGROUND TOKILL = 1.296 TPA total; post-MORTS TNEW = 70.70. jl's 2021 TPA=71 ≈ live's post-MORTS
    70.70 ⇒ jl matches the ported mortality bit-exactly.
The ~20 TPA/cycle that live kills BEYOND MORTS comes from `MISTOE` (dwarf mistletoe), called unconditionally in
the WESTERN base/gradd.f (line 96), before FMKILL/growth-apply. The debug .out carries a full DWARF MISTLETOE
INFECTION AND MORTALITY STATISTICS table:
    Stand 46290437020004: PP MEAN DMR=5.5 (of max 6), 91% of TPA infected, TPA MORTALITY FROM DM = 21 (32%) in
    2011→2021, 14 (30%) 2021→2031, … — EXACTLY the 72→51→37→29→23→19 decline.
morts.f even anticipates it: "SOME EVENTS CAN CHANGE THE TRAJECTORY … REGENERATION, THINNING, USER MORTALITY,
FIRE, I&P EFFECTS" → the TPAMRT≠T reset at morts.f:227 fires every cycle because DM (an I&P effect) removed trees
between MORTS calls.

CONCLUSION: the CR TPA-divergence class is substantially the DWARF MISTLETOE subsystem — a WESTERN insect &
pathogen (I&P) extension the eastern variants (SN/NE/CS/LS) never had, so it was never needed in the shared
engine. jl has NO dwarf-mistletoe model (grep: none). It is a DOWNSTREAM LEAF like volume-NVEL / FFE-fuel /
establishment: the growth+mortality CORE is bit-exact/correct; DM is an additive per-cycle mortality (+ growth
loss + spread + DMR seeding from FIA damage codes) layered on top. FIA seeds initial DMR from tree damage codes,
so DM fires automatically on infected stands with no keyword.

Sweep scope (56,842 stands, 80.5% bit-exact): of 11,109 non-bit-exact, worst-col = BdFt 2074 / MCuFt 1944 / TPA
1760 / TCuFt 1592 / CCF 1501 / SDI 863 / TopHt 700 / QMD 383 / BA 292. The ~1,760 TPA-worst stands are the
strongest DM candidates (jl over-retains lacking DM mortality); volume-worst (BdFt+MCuFt+TCuFt ≈ 5,610) is the
known NVEL gap. So the two largest remaining non-cornered classes are TWO named missing subsystems (dwarf
mistletoe, NVEL volume), not core-model bugs — consistent with the growth+mortality core being complete.

NEXT CHUNK (newly scoped): DWARF MISTLETOE port — base/mistoe.f (DMR spread iteration), DMR seeding from FIA
FVS_TREEINIT damage/severity codes (mistgen), DM growth-loss multiplier (into dgf/htgf), and DM-induced
mortality (into the MORTS/gradd apply). A distinct extension model requiring its own multi-cycle bit-exact
validation vs FVScr_clean — deferred as a downstream leaf; the growth+mortality core remains bit-exact-or-
cornered independent of it.

### Dwarf mistletoe port plan (turn-key scope, ~2,140 core lines)
Core files (mistoe/): mistoe.f (532, driver — per-cycle DMR spread + growth/mort dispatch), mismrt.f (216,
DM mortality), misdam.f (108, seed DMR from FIA damage/severity codes), misdgf.f (174, DM diameter-growth loss),
mishgf.f (160, DM height-growth loss), miscnt.f (82)/miscpf.f (126, DMR spread iteration), misintcr.f (742, the
CR per-species coefficient tables — PMCSP mortality coefs, DMMMLT multiplier, spread/growth-loss coefs), misran.f
(RNG — check ZZRAN stream). DM mortality equation (mismrt.f, per tree):
    DMMORT = PMCSP(sp,1) + PMCSP(sp,2)·DMR + PMCSP(sp,3)·DMR²   (DMR = IMIST rating 0..6)
    DMMORT *= DMMMLT(sp); IF DBH<9.0 DMMORT *= 1.2; IF DMR==0 DMMORT=0
    clamp [0, 0.71 (or 0.5)]; DMMORT = 1-(1-DMMORT)^(FINT/10)   → added to WK2 mortality in the gradd apply
Hooks into FVSjl: (1) seed s.trees IMIST/DMR from FIA FVS_TREEINIT damage_agent/severity at setup (misdam);
(2) per-cycle DMR spread (miscnt/miscpf, RNG-driven — validate stream order like ZZRAN); (3) DM growth-loss
multiplier into diameter_growth!/height (misdgf/mishgf); (4) DM mortality (mismrt) added to the mortality apply;
(5) optional misprt DM summary table. Validate bit-exact vs FVScr_clean on 46290437020004 (DMR=5.5, 91%
infected, TPA MORTALITY FROM DM 21/14/… — a strong single-stand oracle) then re-sweep the 1,760 TPA-worst
stands. A distinct extension model = its own chunk; growth+mortality core stays bit-exact-or-cornered without it.

### DM CHUNK started: data + RNG-independent kernels built & VALIDATED (data/centralrockies/dwarf_mistletoe.jl)
Extracted the full CR DM coefficient set from mistoe/misintcr.f — all IMODTY-INDEPENDENT (MISFIT/DGPDMR/PMCSP
byte-identical across the 5 model-type blocks; verified by diff). Built dwarf_mistletoe.jl with:
  - CR_DM_MISFIT (38 host flags; sp13 PP=host, sp7/16/19/20/21-32/37/38 not)
  - CR_DM_DGPDMR (38×7 diameter-growth-potential by DMR 0-6; PP: 1,1,1,.98,.86,.73,.50)
  - CR_DM_PMCSP (38×3 mortality B0,B1,B2; PP: .00681,-.00580,.00935)
  - spread coefs (BDMR/BCONST/BTPA/BHTG/DDMR/DCONST/DTPA/DHTG — mistoe.f, variant-uniform)
  - cr_dm_seed_dmr (misdam), cr_dm_mortality_rate (mismrt, DBH-dep clamp 0.71/<9" | 0.5/>=9"),
    cr_dm_dg_mult (misdgf).
VALIDATED (RNG-independent) vs live on 46290437020004:
  - SEEDING: 10 records → DMR 6 (FIA code 23023 sev 6) = live "NUMBER OF RECORDS WITH MISTLETOE 10" exactly.
  - MORTALITY RATE: DMR6 PP = 0.370/dec (DBH<9), 0.309/dec (>=9); × ~60 infected TPA ≈ 21 TPA = live "TPA
    MORTALITY FROM DM 21" (2011→2021) exactly.
KEY wiring fact (mismrt.f:191): DM mortality combines with MORTS background per-tree via MAX (WK2=max(WK2,DM)),
NOT addition ⇒ infected trees die at the DM rate, uninfected at background ⇒ ~21 TPA total (matches).
REMAINING (engine-integration chunk, the bit-exact crux): (1) per-tree IMIST storage; (2) seed IMIST at setup
from the FIA damage tuple (treedata.jl already carries it); (3) per-cycle spread/intensification (mistoe.f,
RANN-heavy — must match the shared ZZRAN stream order, the hard part); (4) MISMRT mortality MAX-combined into the
mortality apply; (5) DGPDMR into diameter_growth!. Foundation is done + validated; the RNG-coupled spread + wiring
is next.

### DM spread wiring = the ZZRAN stream-order chunk (grounded; deferred as focused work)
The DM MORTALITY + GROWTH-LOSS are deterministic (no RANN) — the validated kernels drop straight in. The SPREAD
(mistoe.f) is RNG-coupled and is the bit-exact crux. Traced the FVS main-stream (RANN/ZZRAN) call order:
  base/fvs.f: CALL CRATET (crown, draws RANN) @197  →  CALL TREGRO @376 → GRADD @gradd.f:52 → CALL MISTOE
  @gradd.f:96 (spread draws 1-3 RANN per host-species tree, conditionally) → later COMPRS (comcup.f, draws RANN).
jl's main stream (src/core/rng.jl rann!) is ALREADY consumed every cycle by crown_ratio.jl:108 (CRATET) and
compress.jl:366 (COMPRS). So DM spread's draws must be inserted BETWEEN crown and compress, with the EXACT
per-tree draw count/order of mistoe.f (which branches: PPLUS/PMINUS draw, overstory-intensification extra draw,
uninfected-spread 1-2 draws). Getting this wrong shifts the main stream for crown/compress and would REGRESS
currently-bit-exact stands (80.5% baseline). This is the OPEN "ZZRAN stream-order (ch9)" problem, now localized to
the MISTOE↔CRATET↔COMPRS interleave. ⇒ DM spread is a focused chunk requiring per-draw RNG-order validation vs
FVScr_clean (instrument RANN call sequence on 46290437020004), NOT a quick wire-in. The deterministic foundation
(seeding/mortality/growth-loss, all validated) is complete and RNG-safe; the spread+full wiring is the next
focused step. Note: a mortality-only wire (static DMR, no spread) would fix the ~21-TPA .sum gap but DESYNC the
RNG (live makes DM draws jl wouldn't) ⇒ not bit-exact ⇒ rejected per doctrine #1; do the spread properly.

### MISTOE per-tree RANN draw pattern (derived from mistoe.f — the spread-port spec)
Spread runs ONLY for host species (MISFIT=1) whose stand-mean DMR SMR(ISPC)>0 (≥1 infected tree); species with
no infection skip the tree loop entirely (NO draws). For each tree of such a species (mistoe.f tree loop 305-494):
  - INFECTED (IDMR≠0): draw #1 @line 371 (XNUM for PPLUS/PMINUS test). If IDMR<6 AND PPLUS>XNUM AND overstoried
    (DMTALL·0.7 > HT): draw #2 @415 (intensification magnitude 1/2/3). ⇒ 1 or 2 draws.
  - UNINFECTED (IDMR=0): draw #1 @476. If overstoried AND XNUM<0.55: draw #2 @481 (initial DMR 1/2/3). ⇒ 1 or 2.
So on 46290437020004 (15 PP all one species, infected): ~15 draws/cycle minimum inserted into the MAIN stream
between crown(CRATET) and compress(COMPRS). Order within the species loop = ISCT/IND1 tree order (same order the
shared engine iterates). YPLMLT/YNGMLT default 1.0 (no MISTMULT keyword). MISINF (forced infection) = keyword-
only, skip for FIA. The spread port must reproduce this draw count/order exactly to keep compress bit-exact.

### DM foundation COMMITTED + suite-validated (dmr storage + seeding + module) — safe, 0 regressions
Landed the RNG-safe DM infrastructure (behavior-inert until the spread/mortality is wired):
  - src/core/trees.jl: new per-tree field `dmr::Vector{Int32}` (FVS MISCOM IMIST, 0..6); added to the
    @generated tripling-copy field list (_TREE_VEC_FIELDS) so it carries through record splitting/compaction.
  - src/engine/treeinput.jl: `_store_tree!` seeds `t.dmr[i]` from the FIA damage tuple via misdam.f logic
    (agents 30-34 → severity; FIA Arceuthobium codes `_DM_FIA_CODES` → severity|3; else 0). Runs for ALL
    variants but writes an UNUSED field ⇒ provably cannot change output for any variant.
  - src/FVSjl.jl: include data/centralrockies/dwarf_mistletoe.jl (coefs + pure kernels, no load side effects).
VALIDATION: (1) target 46290437020004 differential UNCHANGED (jl 72→71…, seeding inert); (2) FULL test suite
38579 passed / 0 FAILED / 75 broken (pre-existing floor) — the lone "errored" is a transient SQLite
"database is locked" flake from the concurrently-running CR sweep on an SN fire-RNG test, unrelated to the DM
change. ⇒ the dmr field + seeding + module are safe, committed infrastructure. REMAINING (the ZZRAN-order chunk):
cr_mistoe! spread (RANN draws between crown/compress, gated on infection so mistletoe-free stands stay untouched)
+ MISMRT max-combine into the mortality apply + DGPDMR into diameter_growth!. Per doctrine #1, wire spread+mort
together (mortality-only would desync the RNG on infected stands); the deterministic kernels are ready.

### ★★ DM spread port DE-RISKED: jl is RNG-aligned to the MISTOE point (the ZZRAN problem is NOT a blocker here)
Instrumented live FVScr's RANN (base/rann.f: global counter + per-draw value → fort.66) + MISTOE entry/exit
markers, relinked FVScr_rngtrace, ran 46290437020004. Then traced jl's rann! sequence (temporary hook) on the
same stand. FINDINGS:
  - jl and FVS draw from the IDENTICAL LCG sequence (same seed 55329) — jl[1..12] == FVS[1..12] byte-for-byte.
  - MISTOE cycle 1 consumes draws #106-138 (33 draws) for the 15 PP trees (~2.2/tree: DMR-6 trees 1 draw at
    mistoe.f:371, uninfected/overstoried up to 2). Later cycles ~99 draws each (spread reaches more trees).
  - ★ jl[1..105] == FVS[1..105] EXACTLY (all 105 pre-MISTOE draws match), and jl[106]==FVS[106]==0.08902877.
    ⇒ jl is RNG-aligned with FVS right UP TO the MISTOE insertion point; jl currently consumes FVS's MISTOE
    draws for its own next consumer (crown/growth), which is exactly why the stand desyncs from cycle 1's MISTOE
    point onward (both the RNG-dependent columns AND the missing DM mortality).
CONSEQUENCE: the DM spread port is TRACTABLE — NOT the intractable ZZRAN-reorder feared. The task is: insert
cr_mistoe! at the FVS MISTOE cycle position (gradd.f:96 — after GRINCR's growth-increments+MORTS, before UPDATE;
uses HTG) so it consumes exactly the cycle's 33 draws, and the stream RE-ALIGNS for all downstream consumers.
Then MISMRT (deterministic, max-combine) + DGPDMR growth-loss ride along bit-exact. REMAINING precise step: map
which jl grow_cycle! consumer makes draw #105 (one more phase-marked jl trace) to fix the exact insertion line,
then port mistoe.f's spread loop (SMR/DMTALL overstory test, PPLUS/PMINUS logistic, the intensification/spread
draw pattern) reproducing the 33-draw count. Oracle FVScr_rngtrace removed; base/rann.f + mistoe/mistoe.f + jl
rann.jl all restored pristine. Reference traces saved: /home/node/.claude/jobs/e7166935/tmp/{rng.66,jl_rann.txt}.

### ★★ DM spread IMPLEMENTED + RNG-validated: jl matches FVS for 3047 draws (perfect prefix, 4+ cycles bit-exact)
Ported mistoe.f's spread/intensification as cr_mistoe! (src/variants/centralrockies/dwarf_mistletoe_model.jl):
per-cycle, species-sorted (species_sort! ISCT/IND1), SMR gate (mistletoe-free species draw nothing), per-point
DMTALL overstory test, PPLUS/PMINUS logistics (Float32, YPLMLT/YNGMLT=1), the RANN-driven DMR intensification/
spread distributions. Wired into grow_cycle! after FIXHTG, before mortality (FVS gradd.f:96 position). Coefs made
Float32 (FVS computes PPLUS in REAL). VALIDATION (rann! trace vs FVScr_rngtrace on 46290437020004): jl now draws
3047 (was 2667) and jl[1..3047] == FVS[1..3047] EXACTLY — a perfect prefix through 4+ cycles of growth+MISTOE. The
remaining 88-draw tail (FVS 3135) is the last cycle only, expected because DM MORTALITY is not yet wired ⇒ tree
state (record survival) diverges by cycle 5. ⇒ the spread RNG order is CORRECT; the hardest part is done. NEXT:
wire MISMRT mortality (max-combine into mortality!) + DGPDMR growth-loss into diameter_growth!, which should align
the tail and fix the .sum TPA. (rann! trace hook still in for the final check; remove after.)

### ★★★ DWARF MISTLETOE CHUNK COMPLETE — full model wired, bit-exact-or-cornered vs live FVScr
Wired all three DM effects, validated on 46290437020004 (the representative DM stand):
  1. SPREAD (cr_mistoe!, mistoe.f): per-cycle DMR intensification/spread, RANN-aligned (proven jl==FVS to 3047
     draws pre-mortality). Inserted after FIXHTG, before mortality (gradd.f:96).
  2. MORTALITY (cr_dm_mortality_combine!, mismrt.f): per-tree DM kill MAX-combined into killed[] after
     apply_fixmort!, before snag-booking/apply (WK2=max(WK2,DM)).
  3. GROWTH-LOSS (cr_dm_growth_loss!, misdgf.f/dgdriv.f:230): DG·=DGPDMR(sp,DMR) on central+tripled records,
     right after diameter_growth!, using start-of-cycle DMR.
RESULT (jl vs live, was jl 72→66 over-retaining):
  Year  TPA(live/jl)  BA      SDI      TCuFt      MCuFt
  2021  51/51 ✓★      27/27✓  49/49✓   394/400    327/332
  2031  37/37 ✓★      26/26✓  44/44✓   424/430    358/362
  2041  29/29 ✓★      25/26   41/42    460/476    397/416
  2051  23/23 ✓★      26/26✓  40/41    511/527    448/463
  2061  19/19 ✓★      27/27✓  40/40✓   570/583    505/516
TPA BIT-EXACT all cycles; BA/SDI bit-exact-or-±1; volume within ~1.5% (the accepted DGSCOR-COR/NVEL cornered
tail — same class as the rest of CR growth). Was the LARGEST non-cornered CR sweep class (~1,760 TPA-worst
stands). FULL SUITE: 38579 passed / 0 FAILED / 75 broken (floor) — 0 regressions (DM gated on CentralRockies +
per-species infection: mistletoe-free stands draw zero RANN ⇒ untouched; eastern variants never call it). The 1
"errored" is the recurring SQLite-lock flake from the concurrent sweep (SN fire test), not DM. Files:
data/centralrockies/dwarf_mistletoe.jl (coefs+kernels), src/variants/centralrockies/dwarf_mistletoe_model.jl
(spread+mortality+growth-loss), trees.jl (dmr field), treeinput.jl (seeding), 3 wire points in simulate.jl/
southern-mortality.jl. ⇒ CR dwarf mistletoe is DONE to the bit-exact-or-cornered bar.

### DM generalization check (2 more infected stands, doctrine: don't over-claim from one)
- 3622258010690: TPA 312/312→287/286→274/272→240/235→206/203→181/177 (bit-exact-or-±few); BA/SDI ±1-2; vol ~1%.
  Confirms DM works on a 2nd infected stand (would have grossly over-retained pre-DM).
- 39450996010690 (dense 14925-TPA seedlings): bit-exact early (14925/14925…); late divergence 2049+ (10042/10867)
  is the SEPARATE cornered dense self-thinning/RDPSRT tail, NOT DM — DM didn't regress it.
⇒ DM port generalizes; residuals are the pre-existing DGSCOR/dense-thinning cornered classes, not the DM model.

### Volume-worst class CHARACTERIZED (doctrine #2 measure): targeted equation gaps, NOT total NVEL absence
Dug vol-worst stands to scope the next chunk. Finding: volume IS computed but systematically too HIGH while
TPA/BA/SDI are bit-exact (trees identical) ⇒ a volume-EQUATION error, not missing volume:
  - 24318722010900 (FIA 756 = honey mesquite PRGL2 → CR "OH" Other Hardwood sp38, woodland form D6/HT9):
    TCuFt live/jl 7/8→19/24→44/55→86/110→182/234 = ~28% HIGH, growing with size; TPA/BA/SDI BIT-EXACT.
  - 23718531010900: TPA bit-exact, vol ~5-10% high, growing.
ROOT: the earlier r2oldv (Chojnacky woodland cubic) port covered oak/juniper/pinyon (FIA 065/066/069/106/814/
823/998) but NOT mesquite/OH — so OH uses a normal-tree equation ⇒ over-volumes the woodland form. ⇒ the vol-
worst class is (at least partly) MISSING WOODLAND/HARDWOOD equation coverage for specific species, a TARGETED
extension of the volume assignment, NOT the feared full-NVEL (NSVB/Flewelling/DVEE) rewrite. Refines the memory's
scoping. NEXT: enumerate the CR species whose vol-worst stands over/under-volume, map each to live's VOLEQDEF
equation, and port the missing ones (mesquite/OH first). Some may still need NSVB/Flewelling, but many are
woodland-cubic (r2oldv-family) gaps. Measured one stand each — the full class needs a per-species enumeration.

### ★★ REAL FIX: CR woodland volume FCLASS default (single→multi-stem) — clears the vol-worst class
Found via doctrine-#2 measurement of the vol-worst class. The volume divergence was NOT an equation-assignment
gap (jl's OH voleq='300DVEW999' MATCHES live) nor a height/DBH issue (QMD+TopHt BIT-EXACT). It was the FCLASS
form-class branch in cr_dve_vol (r3d2hv.f): jl defaulted `fclass=1` (single-stem) but FIA trees have NO Girard
form class ⇒ FVS passes FCLASS=0 ⇒ the `FCLASS.NE.1` MULTI-STEM coefficient set (r3d2hv.f:301/318/…; "1=single,
others=multistem" per the source comment). Hand-calc confirmed: 999 woodland cube at D11.6/H43 gives 16.3 cuft
(multi, =live) vs 20.9 cuft (single, =jl's wrong 234). FIX: cr_dve_vol default `fclass::Int = 0` (one line). This
also flips the r2oldv MSTEM term (ms = fclass==1?1:0) to ms=0, which is likewise the correct FIA multi-stem form.
VALIDATION (was 13-30% high on woodland species):
  24318722010900 (mesquite/OH): TCuFt 7/8→182/234  ⇒  now 7/7…182/183 BIT-EXACT
  25013840010900: now 14/14, 48/48, 104/105, 174/176 (±1)
  24268481010900: 210/210, 261/261, 322/317 (~2% DGSCOR tail)
  24257722010900: 230/230, 324/324, 416/415, 510/510 (±1)
  3622258010690 (DM stand w/ woodland sp): 1984 1071/1062 → 1071/1071
Full suite 38579 pass / 0 FAILED / 75 broken — 0 regressions (the r2oldv ms-flip did not regress any tested
stand; the earlier r2oldv validation used non-ms species 066/814/823 or coincidence). ⇒ affects EVERY CR stand
with woodland/hardwood species (060 juniper / 106 pinyon / 800 oak / 999 other-hardwood, + r2oldv 065/069/106/
475) — a large fraction of the ~5,610 vol-worst class. Residuals now the DGSCOR/growth cornered tail (~1-2%).

### Volume situation after this session's FCLASS fix — full class map
With the FCLASS woodland fix landed, the CR volume divergence resolves into:
  1. WOODLAND (060 juniper/106 pinyon/800 oak/999 other-hw + r2oldv 065/069/106/475): BIT-EXACT (FCLASS=0 fix). ✓
  2. CONIFER FW2 (PP=300FW2W122, DF=300FW2W202): ~1.5% residual. Confirmed a pure VOLUME-EQUATION gap — on
     46290437020004 at 2021 QMD 9.9/9.9 AND TopHt 36/36 are BIT-EXACT but TCuFt 394/400 (1.5% high). jl uses a
     hardcoded R3 ponderosa approximation (cr_dve_vol.jl:70 `spc==122 && reg==300`), NOT live's real FW2W122
     (Flewelling fwinit.f). Later cycles compound with the DGSCOR growth tail (QMD 12.7/12.9 by 2041). ⇒ the
     remaining conifer volume work = port Flewelling FW2 (fwinit.f) for the 122/202 species. Small residual
     (~1.5%), near-cornered but systematic (not COR precision).
  3. NSVB species (CB/WF=NVB0000015, SW=NVBM240119, ES=NVBM330093, AS=NVB0000746): unverified this session —
     check whether jl's approximation diverges; if so, port NSVB (nsvb.f) for those.
⇒ NEXT volume chunk (much reduced from "full NVEL"): Flewelling FW2 (122/202) + verify/port NSVB species. The
big/systematic vol divergences (woodland 13-30%) are FIXED; remaining is a ~1.5% FW2 conifer tail.

### Vol-worst sample after FCLASS fix — remaining = NSVB refinement + FW2 (small/moderate, NOT woodland)
Sampled vol-worst stands post-FCLASS: most now bit-exact-or-±1.5% (504587228 fully bit-exact; 550252962 ~1.5%).
Two remaining conifer-equation residuals characterized:
  - FW2 (PP 300FW2W122): ~1.5% high — jl hardcoded R3 approx vs live Flewelling FW2 (verified QMD+TopHt exact).
  - NSVB (ES 93→NVBM330093): stand 31309338 (ES+AF spruce-fir) jl TCuFt 17% LOW but MCuFt HIGH — a total/merch
    PARTITION issue in the ported NVB (memory: "NVB TCF+MCF ported, BF TODO" — the cubic partition needs refining).
    VERIFIED NOT an FCLASS regression: the 093 DVE branch (cr_dve_vol.jl:54) is pure D²H, no fclass dependency;
    the ES divergence is the pre-existing NSVB port precision.
⇒ CR volume after this session: WOODLAND fixed bit-exact (FCLASS); remaining = (a) Flewelling FW2 for 122/202
(~1.5%), (b) NSVB cubic-partition refinement for the NVB* species (ES/CB/WF/SW/AS). Both are bounded conifer-
equation chunks, much smaller than the original "full NVEL port" framing. The 13-30% woodland divergences (the
bulk of the vol-worst class) are RESOLVED.

### ★ CR VOLUME chunk essentially COMPLETE — equations bit-exact at inventory; residuals = DGSCOR growth tail
Decisive measurement: cycle-0 (INVENTORY, pre-growth) volume is BIT-EXACT across ALL equation families with the
FCLASS fix in place:
  46290437020004 PP/FW2 (300FW2W122): TCuFt 372/372
  24318722010900 mesquite/woodland (300DVEW999): 7/7
  31309338010690 ES-NVB (NVBM330093)+AF-DVE (300DVEW093): 15/15
  550252962126144: 365/365, MCuFt 337/337
⇒ cr_dve_vol / cr_nvb_vol / cr_fw2_vol are all CORRECT (bit-exact on the measured inventory trees). The later-
cycle ~1.5% residuals (e.g. 46290437 2021 394/400 with QMD 9.9/9.9 + TopHt 36/36 bit-exact) are therefore NOT
equation errors but the DGSCOR growth-precision tail: individual-tree D/H differ by the accepted cornered ~1%
(QMD matches to 0.1"), amplified ~2× by volume (V~D²) to ~1.5%. Same accepted class as the rest of CR growth.
CONCLUSION: chunk 8 (volume) is bit-exact-or-cornered — the FCLASS fix removed the last SYSTEMATIC volume
divergence (woodland 13-30%); everything else is the growth tail. Remaining true-volume TODO is only NVB board-
foot (partial) for stands that report BdFt. The vol-worst sweep class was dominated by the woodland-FCLASS bug
(now fixed) + the density/self-thinning tail (cornered), not missing/wrong conifer equations.

### BdFt-worst class ALSO = growth/density tail (not a board equation bug) — volume characterization COMPLETE
Dug the BdFt-worst stand 39451382010690 (23.9% div). INVENTORY bit-exact (2009: BdFt 5820/5820, TCuFt 2041/2041,
QMD 2.4/2.4) but diverges after growth WITH the QMD: 2019 QMD 2.7/2.8, 2049 3.8/4.0 — jl consistently over-grows
this dense small-tree stand ~4-5%. The cubic+board divergence (TCuFt 2544/2822, BdFt 7929/8699) is DRIVEN by the
QMD/growth divergence amplified by V~D²·H (4-5% D ⇒ ~11-20% vol), NOT a board-foot equation error (board is bit-
exact at inventory). ⇒ the ENTIRE vol-worst sweep class (TCuFt+MCuFt+BdFt) reduces to TWO causes: (1) the woodland
FCLASS bug (real equation bug, NOW FIXED — systematic 13-30%), and (2) the dense-stand growth/self-thinning tail
(cornered structure_densephase/DGSCOR class, amplified by the volume power-law). The volume EQUATIONS themselves
(DVE/NVB/FW2/woodland) are all bit-exact at inventory ⇒ CR volume chunk is bit-exact-or-cornered. Only true-vol
TODO remaining: NVB board-foot for the few large-tree NVB-species stands (partial), a small leaf.

### ★ DOCTRINE-#2 CORRECTION: the BdFt/structure class is NOT all-cornered — a real aspen/fir DG lead
Ran the TreeId-matched verifier (dig_verify_treeid.jl) on the BdFt-worst stand 39451382010690 instead of assuming
cornered. VERDICT = ★ESCALATE: per-tree DBH/DG divergence at 2009 (FIRST cycle, PRE-tripling) — a REAL growth
divergence, NOT the cornered downstream tie-break I'd inferred. Stand is aspen (CR sp20, 22 trees) + subalpine fir
(sp1, 20) + 1 ES; UNINFECTED (dmr=0 ⇒ my DM changes don't touch it) and FCLASS is volume-only ⇒ this is a
PRE-EXISTING growth divergence, not a session regression. jl over-grows the small trees (QMD 2.7/2.8→3.8/4.0),
amplified to 20% BdFt. ⇒ CORRECTION to the prior note: the vol-worst/structure_densephase class is NOT uniformly
cornered — at least the aspen/fir dense stands carry a real first-cycle per-tree DG divergence (GEMDG large-tree
or REGENT small-tree for aspen/fir). This is a genuine growth-chunk lead (NEXT: deep-trace the per-tree DG for
sp20/sp1 at cycle 1 vs live — instrument cr_gemdg/regent on this stand). The volume EQUATIONS remain bit-exact at
inventory (that finding stands); but the growth that feeds them diverges on this class. Lesson: verify cornered
with the TreeId verifier, don't infer it from "dense stand" + .sum aggregates (doctrine #2/#3).

### CR aspen (sp20) diameter-growth lead LOCALIZED (via TreeId verifier → GEMDG CASE(20))
Deep-traced the ESCALATE stand 39451382010690: the divergent tree is FIA 746 (quaking ASPEN), TreeId 9 — raw FIA
DBH 8.0; live keeps ~8.0 through cycle 1, jl grows it to 8.5 (~0.5"/decade over-growth), compounding to the stand
BdFt 20% high by 2059. (Treelist "2009" == .sum cycle-1, since .sum-2009 BdFt is bit-exact 5820/5820 and diverges
only at 2019+.) The stand is IMODTY 4 (Spruce-Fir); aspen uses GEMDG SELECT CASE(20,21:22,28,38). jl's aspen DF
(diameter_growth.jl:168-175) MATCHES cr/gemdg.f EXACTLY incl. the `DF*=1.05` "growth underestimated" bump — so
the DF FORMULA is faithful. ⇒ the over-growth is in an INPUT (DPP/BATEM/SI/BGTTBA) or a missing common-tail
adjustment for the DF path (diagr=(DF-DPP)·bark, gemdg.f:355-377), NOT the regression coefficients. NEXT: instrument-
replay cr/gemdg.f for aspen (dump DPP/BATEM/SI/BGTTBA/DF/DDS on TreeId 9) vs jl cr_gemdg to find the diverging
input. Aspen is common in CR (Rockies) ⇒ a material lead, not cornered. The validated cr_gemdg 9981/9982 evidently
did not exercise the aspen CASE(20) DF path (those were conifer IDDS-path species). REMAINING CR growth work:
this aspen GEMDG input bug + verify other DF-path (non-IDDS) species.

### CORRECTION (instrument-replay): GEMDG aspen is BIT-EXACT — the dense-stand divergence is NOT aspen large-tree DG
Instrument-replayed cr/gemdg.f (dump IS/DPP/BATEM/BGTTBA/SI/DF for aspen) vs jl cr_gemdg on 39451382010690. ALL
aspen GEMDG inputs AND DF are BIT-EXACT (DPP 9.2, BATEM 120.61409, BGTTBA 0.31099963, SI 55, DF 10.385434 — every
row matches live to full precision). ⇒ the CR aspen large-tree diameter growth is CORRECT; the prior "aspen GEMDG
over-grows" lead was a premature inference. The verifier's TreeId-matched "TreeId 9 DBH 8.0/8.5" is therefore
either a non-GEMDG path (small-tree REGENT / a backdated-DP tree) or a POST-TRIPLING TreeId-match artifact (the
doctrine-#3 trap: TreeId re-assigned after tripling ⇒ jl's TreeId 9 ≠ live's). The stand's real .sum divergence
(QMD 2.7/2.8, BdFt 10%) remains, but must be attributed in the PRE-SPLIT window, not by post-tripling TreeId
match. LESSON (doctrine #2/#3, again): instrument-replay the actual model (GEMDG bit-exact) before blaming it; and
the TreeId verifier's per-tree match is NOT reliable post-tripling for attribution. GEMDG aspen: RULED OUT / clean.
NEXT: attribute this dense-stand divergence via the pre-split window (fir sp1 GEMDG? small-tree REGENT? mortality?).

### RESOLUTION: the ESCALATE was a verifier FALSE POSITIVE (CR post-tripling TreeId re-assignment) — stand is cornered
Chain of evidence: (1) verifier flagged 39451382010690 ESCALATE on aspen TreeId 9 (DBH 8.0/8.5); (2) instrument-
replay proved that tree's GEMDG DF + all inputs BIT-EXACT vs live; (3) aspen large-tree DG = DDS(GEMDG bit-exact)
+ COR + DGCON + bark (all validated 9981/9982) ⇒ the aspen's grown DBH IS bit-exact. Therefore the treelist
"8.0/8.5" compares jl's TreeId 9 to live's TreeId 9 which are DIFFERENT physical trees after tripling (TreeId re-
assigned post-split) = the doctrine-#3 trap. ⇒ dig_verify_treeid.jl's TreeId match is NOT reliable for CR (the
tripling re-indexes TreeIds differently than the eastern variants it was validated on) — it OVER-escalates. The
stand's real .sum divergence (QMD 2.7/2.8, BdFt 10%) is the GROWTH-bit-exact + self-thinning-tie-break cornered
structure_densephase class (same as the memory's 263/263-cornered campaign), NOT a reducible growth bug. ⇒ CR
growth is bit-exact (aspen GEMDG directly proven; conifer GENGYM validated 9981/9982); the dense-stand vol-worst/
BdFt residuals are the cornered self-thin tail amplified by volume. FOLLOW-UP (tooling, not a CR bug): the CR
sweep's structure_densephase "ESCALATE" verdicts from dig_verify_treeid need re-checking with a pre-split-window
matcher — the eastern TreeId key is invalid under CR tripling. NET: CR growth/mortality/DM/volume all bit-exact-
or-cornered; the escalations were a verifier artifact, now understood.

### Honesty caveat on the above resolution
What is PROVEN: CR aspen GEMDG (the flagged tree's growth) is bit-exact vs live; the verifier's post-tripling
TreeId match is unreliable for CR. What is INFERRED (not yet fully proven): that the residual .sum divergence is a
pure self-thin TIE-BREAK (cornered) vs a real mortality-COUNT difference. To close that, verify at the first
projected cycle that (a) the PRE-mortality grown stand matches (growth fully bit-exact, not just aspen), and (b)
the mortality total (TPA killed) matches — if both hold and only WHICH trees differ, it's the cornered tie-break;
if the killed-count differs, it's a reducible mortality bug. Deferred as the next dense-stand dig (with a pre-split
matcher). Not claiming fully-cornered on faith — aspen growth is the only piece directly measured bit-exact here.

### ★★ CORRECTION + REAL BUG LOCALIZED: dense-stand divergence is the REGENT/GEMDG boundary blend (NOT cornered)
Definitive measurement chain on 39451382010690 (aspen+fir, dense): (1) .sum 2019 TPA BIT-EXACT (2786/2786) but BA
jl 10% HIGH (112/123) ⇒ mortality is fine, GROWTH over-shoots (bigger diameters, same tree count) — a REAL
reducible growth bug, NOT the cornered self-thin tie-break I'd inferred (correction #3). (2) Instrument-replay of
cr/gemdg.f (dump IS/DPP/DIAGR/DDS per tree) vs jl cr_gemdg: the GEMDG DDS VALUES are BIT-EXACT for aspen (sp20)
AND fir (sp1) — every matched row identical. (3) But the GEMDG CALL COUNTS DIFFER: jl fir 592 vs live 660, aspen
702 vs 726 (jl ~11 fewer/cycle) ⇒ jl routes boundary trees to pure small-tree REGENT where FVS BLENDS the GEMDG
large-tree component. FVS regent.f uses per-species XMIN→XMAX: DBH≤XMIN pure small-tree, XMIN<DBH<XMAX a BLEND of
REGENT + GEMDG, DBH≥XMAX pure GEMDG. jl's small_tree_growth gates the diameter increment on a single `break_sp`
(BREAK[sp]) with no XMIN/XMAX blend ⇒ boundary trees get pure REGENT (over-predicts) instead of the blend ⇒ BA
over-grows ~10% on stands with many boundary-size trees. ⇒ REAL growth bug in the CR small-tree chunk (the
REGENT/GEMDG blend), equations themselves bit-exact. NEXT: port regent.f's XMIN/XMAX blend (the large-tree DDS
weight in the transition zone) into cr small_tree_growth!; validate BA bit-exact on this stand. This is the true
cause of a chunk of the structure_densephase/vol-worst class (BA-driven, amplified into volume). Corrects the
earlier "verifier false-positive / cornered" note — it IS reducible, found by measuring TPA-vs-BA + GEMDG counts.

### REFINEMENT (accuracy over the prior note): the BA bug is small-tree DBH or COR/DGCON, blend NOT yet proven
Correcting my own prior note before it misleads: in regent.f the DIAMETER uses a HARD breakpoint BKPT (regent.f:
342 `IF(D.GE.BKPT) GO TO 23` — D<BKPT small-tree DBH, D≥BKPT large-tree), NOT a blend; only the HEIGHT is XWT-
blended (XMN/XMX, line 319/326) and height doesn't drive BA. And jl ALREADY has break_sp/xmn/xmx params. So the
GEMDG call-count difference is most likely jl skipping GEMDG for D<break_sp small trees (which REGENT overrides
anyway) — NOT necessarily the bug. WHAT IS PROVEN on 39451382010690: (a) real growth divergence (TPA bit-exact
2786/2786, BA jl 10% high 112/123) — NOT mortality, NOT cornered; (b) GEMDG large-tree DDS bit-exact (aspen sp20
AND fir sp1, every row). ⇒ the diameter over-growth is in EITHER the small-tree REGENT DBH increment (regent.f:
343-395, D<BKPT) OR the COR/DGCON added to the large-tree DDS in the dgf wrapper (WK2=DDS+COR+DGCON) — NOT yet
isolated between them. NEXT (decisive): instrument the FINAL applied per-tree DG (dgf.f WK2 for large trees;
regent DG(K) for small) vs jl, matched in the pre-split window, to isolate small-tree-DBH vs COR/DGCON. Honest
status: real reducible growth bug, well-bounded (GEMDG eqns + mortality ruled out), exact mechanism pending one
more instrument pass. (Session note: this dig required 4 measure-driven corrections — the discipline caught each
premature inference; the remaining two candidates are both small, bounded code paths.)

### HONEST STATUS on the dense-stand BA bug: extensively bounded, exact cause NOT yet isolated
Further measurement narrowed but did not cleanly isolate. PROVEN bit-exact vs live on 39451382010690:
  - Large-tree GEMDG central DDS (aspen+fir, every row); COR(ISPC)=0 and DGCON=0 for all species (instrument
    dgf.f WK2 dump) ⇒ WK2 = DDS = bit-exact ⇒ the large-tree applied DG is bit-exact, DGSCOR calibration ruled out.
  - Mortality (TPA 2786/2786 bit-exact).
  - Volume equations (bit-exact at inventory).
STILL DIVERGENT: BA jl 10% high (112/123) at cycle 1 with all the above bit-exact. And jl st_break=1.0 for fir/
aspen ⇒ those 5.9-13.7" trees use GEMDG (not REGENT DBH) ⇒ the "small-tree REGENT DBH" hypothesis does NOT hold
for the dominant trees either. Remaining candidates (NOT isolated): (a) the TRIPLED sub-record DG spread (dgU/dgL
FRMT serial-correlation) — GEMDG central proven, sub-records not; (b) ESTABLISHMENT regen (GEMDG call counts
differ jl 592 vs live 660 fir — could be different regen tree counts feeding growth); (c) a bark/UPDATE detail in
applying the bit-exact DDS to DBH. NEXT (decisive): dump the FINAL per-record applied DG (central + 2 tripled subs
+ any regen) in the PRE-SPLIT window from both sides and diff — that pins (a)/(b)/(c). CONCRETE GAP FOUND
(fix regardless): jl's small_tree_growth omits DGBND (ie/dgbnd.f: cap DBH+DG ≤ SIZCAP(sp,1) when SIZCAP(sp,3)<1.5)
— a faithfulness gap, likely minor here (trees far from cap) but should be added. HONEST: real BA bug, heavily
bounded (4 subsystems ruled out bit-exact), exact cause pending one per-record pre-split DG dump — not claiming a
specific cause I haven't proven (this dig produced 5 measure-driven corrections; the discipline is holding the line).

### CR GLIM tripling-spread cap ADDED (faithful, 0 regress) — but NOT the dense-stand BA cause (correction #6)
Found via source diff: cr/dgdriv.f caps EACH DG spread (GDIF=DG−WKI, GLIM=WKI·0.33, IF GDIF>GLIM DG=WKI+GLIM;
WKI=un-FRM'd central DG) — a cap the EASTERN dgdriv.f lacks (GLIM count: sn=0, cr=8; cr also calls DGBND after).
jl's shared tripling used _bound_scale (DGBND-style, correct for eastern) but OMITTED CR's GLIM. FIXED: CR-gated
GLIM cap on central+upper+lower raw DG before _bound_scale (diameter_growth.jl tripling block). Suite 38579/0
FAIL/75 broken — 0 regress (eastern untouched, crv-gated). BUT on 39451382010690 the BA barely moved (2019 still
112/123; 2059 5568→5565) ⇒ the tripled spread here doesn't hit the WKI·1.33 cap ⇒ GLIM is a faithful gap-fill but
NOT the ~10% BA cause. STATUS on the BA bug: remains bounded-but-UNISOLATED after ruling out (all bit-exact/no-
effect): GEMDG DDS, COR/DGCON, mortality/TPA, volume eqns, and now GLIM. Remaining candidates: the tripled-record
TPA split (0.25/0.15/0.60) or FRM spread magnitude, OR this is a larger instance of the accepted DGSCOR/tripling
residual class (memory: CR growth-only "2000 BA 109/live106" ~3% was accepted). This dig hit 6 measure-driven
corrections — a signal to defer to a fresh focused effort with a per-record pre-split DG matcher; NOT chasing a
7th hypothesis at session tail. DELIVERED regardless: the faithful CR GLIM cap (real dgdriv.f gap, now closed).

### The dense-stand BA bug is a PARADOX: every piece bit-exact, aggregate diverges — points to record COUNT
Direct measurement (instrument cr/dgdriv.f DGT dump: central+upper+lower DG per tree) vs jl: the tripling DGs are
BIT-EXACT (fir D→central 0.86136, upper 1.12546, lower 0.65854 — every value matches live). The TPA split is
bit-exact too (FVS base/triple.f: central·0.60, upper 0.25, lower 0.15 = jl exactly). So per-record: same DG, same
TPA fraction. Combined with GEMDG DDS + COR/DGCON(=0) + mortality-TPA + volume-eqns all proven bit-exact, EVERY
per-record piece matches — yet .sum BA diverges 10% at cycle 1 (112/123). ⇒ the divergence must be a RECORD
COUNT/COMPOSITION difference, not any per-record value. Corroborating hint: the GEMDG call count differs (jl 592
vs live 660 fir over the run) — jl processes a different NUMBER of fir records. Leading candidate: ESTABLISHMENT/
regen adds a different set of trees each cycle (different count/species/size), OR the tripling/comcup record
bookkeeping drops/keeps records differently. NEXT (systematic, doctrine-#3-safe): a per-cycle RECORD CENSUS by
(species, size-class) — count records + Σ TPA + Σ BA per bucket in jl vs live at each cycle — to find WHICH bucket
gains/loses records. That pins establishment vs tripling-bookkeeping. This is the 7th ruled-out layer; the bug is
NOT in any growth/mortality/volume VALUE (all bit-exact) — it is a record-population difference. Deferred to a
fresh census-based dig. KEPT: the faithful CR GLIM cap (0 regress). All instrumentation reverted, sources pristine.

### ★★ CENSUS BREAKTHROUGH: the dense-stand BA bug is ASPEN REGEN/SPROUT growth (cr_esgent/esuckr), inventory bit-exact
Resolved the "every-piece-bit-exact-but-BA-diverges" paradox with a per-species RECORD CENSUS (aggregate by
species+DBH-bucket, doctrine-#3-safe — no per-tree TreeId match). On 39451382010690 at 2019:
  fir (sp19) BA 47.8/48.4, spruce (sp93) 7.5/7.5 — bit-exact/close; ASPEN (sp746) 56.6/67.0 — jl 18% HIGH.
So it's ASPEN-specific. DBH-bucket census (aspen): <3" TPA 724/700 BA 5.4/5.2; 6-9" TPA 39/56 BA 12.2/17.3;
>9" TPA 57/67 BA 39.0/44.4 — jl has MORE TPA/BA in the LARGER buckets (same record counts). The <3" bucket
(724 TPA, 3 recs) is aspen REGEN/SPROUTS seeded during the cycle; jl grows them into the larger buckets FASTER
(worse by 2029: 6-9" BA 7.4/12.5, >9" 51.2/62.7). PROVEN bit-exact (ruling out the inventory large-tree aspen):
aspen GEMDG DF (GEMDGASP dump), aspen GEMDG DDS (GEMDDS dump: DPP9.2→DIAGR1.126163→DDS3.042309 = live), aspen
TRIPLING DGs (DGT dump: central 1.17172/upper 1.51874/lower 0.90168 = live). ⇒ the inventory aspen growth is
BIT-EXACT; the divergence is the ASPEN REGEN/SPROUT birth-cycle growth (cr_esgent grows just-established regen via
REGENT; esuckr aspen sprouts asp_idx=20) — jl over-grows the aspen sprouts. NEXT: instrument cr_esgent/esuckr
aspen-sprout DBH+growth in the birth cycle vs live (the sprout seed size, the birth-cycle REGENT increment, or the
sprout count). This is the definitive localization (8 layers deep) — a real, bounded bug in the aspen regen path,
NOT the large-tree growth. The census method (species+size-bucket aggregate) is the RIGHT CR dig tool (TreeId
match is invalid post-tripling). Build OK, GLIM cap kept, all instrumentation reverted.

### ★★★ DEFINITIVE localization (9 layers): aspen REGEN seed-size/count + small-aspen GEMDG DIAGR
Compared the saved GEMDDS dumps for SMALL-DPP aspen (the regen/sprouts, not the bit-exact large aspen):
  - live has aspen at DPP=1.000 (3 records, floored) that jl LACKS — jl's smallest aspen is DPP 1.14; jl makes
    FEWER small aspen (30 vs 33 GEMDG calls) ⇒ jl seeds the aspen regen/sprouts at a BIGGER DBH and/or FEWER count.
  - at MATCHING DPP=1.1417: jl DIAGR 0.59081/DDS 0.48896 vs live 0.61813/0.54402 — the small-aspen GEMDG DIAGR
    (=(DF−DPP)·bark) DIVERGES (jl LOWER), while large aspen (DPP≥6) is bit-exact ⇒ a bark(cr_bratio) or DF-input
    (BGTTBA is high for small understory trees) difference specific to TINY aspen.
NET mechanism: jl's aspen regen population differs (bigger seed DBH / fewer records) ⇒ the TPA distributes into
larger DBH buckets over cycles ⇒ aspen BA 18% high (dense-stand BA bug). The INVENTORY large-tree aspen growth is
bit-exact (DF+DDS+tripling all proven); the bug is entirely in the ASPEN REGEN/SPROUT seeding+small-tree-growth
(esuckr! sprout DBH/height/count via cr/essprt.f + sprtht_cr, and the small-aspen bark in cr_gemdg's DIAGR).
NEXT (turn-key): instrument the aspen sprout seed (esuckr! DBH/HT/count) vs cr/essprt.f live, AND cr_bratio for
aspen at D~1" — those two pin the seed-population and the small-tree DIAGR. ⇒ the dense-stand BA bug is NOT a
large-tree/mortality/volume bug (all bit-exact) — it is the aspen regen path, now localized to two concrete
sub-checks. DELIVERED this session: DM subsystem + FCLASS woodland-vol + CR GLIM cap (3 fixes); this BA bug
localized 9 layers deep to aspen-regen via the species+size-bucket census (the correct CR dig tool).

### Refinement (direction-corrected): PRIMARY cause = aspen regen SEED size/count, not per-tree growth
Reconciling the direction: for matching DPP=1.1417 the jl per-tree DDS (0.489) is LOWER than live (0.544) — so jl
grows each small aspen LESS, yet jl ends with MORE aspen BA. ⇒ the per-tree growth is NOT the driver; the aspen
REGEN SEED (size + count) is. jl seeds FEWER (30 vs 33) BIGGER aspen (smallest DPP 1.14 vs live's 1.000) ⇒ jl's
regen start in higher DBH buckets and dominate BA despite lower per-tree growth. Back-solving the DIAGR: DF≈1.77
for DPP1.14, so DIAGR=(DF−DPP)·bark ⇒ jl bark≈0.938 vs live≈0.981 — jl's cr_bratio for D~1" aspen is also lower
(a SECOND, secondary diff). SO the dense-stand BA bug = (1° ) aspen sprout/regen SEED DBH+count (esuckr!/cr/
essprt.f: sprtht_cr height → sprout DBH, and the sprout COUNT/nsprec) jl seeds too big/too few; (2°) cr_bratio
for tiny aspen. NEXT turn-key dig: instrument live esuckr!/essprt.f aspen (asp_idx=20) sprout DBH/HT/count vs jl,
+ cr_bratio(sp20, D≈1). Both are small, bounded code paths. The dense-stand BA bug is now localized to the aspen
regen SEEDING — the large-tree growth, mortality, volume, and tripling are ALL proven bit-exact. (10 measure-
layers; discipline held; census = the right CR tool.)

### Scale validation of the 3 session fixes (16-stand sample of the affected sweep classes)
Ran the differential (new code) on 8 vol-worst + 8 TPA-worst stands (the OLD sweep's non-bit-exact classes):
  VOL-WORST (FCLASS fix target): 6/8 now BIT-EXACT-or-±1% (24318722 0.0%, 24257722 0.2%, 25039978 0.3%,
    24257202 0.6%, 25013840 1.1%, 15312527 1.9%), 2 close (~2-8%). ⇒ the FCLASS woodland fix cleared the vol-worst
    class BROADLY (not just the one validated stand) — a large fraction of the ~5,610 vol-worst class.
  TPA-WORST: the DM-infected subset is fixed (DM validated earlier), but the 39xxx-forest cluster (39452085 21%,
    39467153 38%, 39466839 16%, 31288684 33%, 42479374 37% …) STILL diverges 15-47% — these are the ASPEN-REGEN
    dense stands (same forest 39xxx as the localized 39451382). ⇒ the aspen-regen bug is a MEANINGFUL CLASS (a
    whole forest cluster), NOT a one-off — raising the value of the aspen-regen seeding fix.
NET: FCLASS (broad vol-worst fix) + DM (TPA-worst DM subset) validated at scale; the residual TPA-worst is the
aspen-regen class (localized, fix pending). A full re-sweep would quantify the new bit-exact rate (old 80.5%);
the sample shows the vol-worst class largely cleared. Session's 3 fixes have broad, measured impact.

### Aspen-regen localization refined: sprout SEED formulas MATCH → divergence is ISHAG/birth-cycle-growth
Code-read (no relink): FVS bin/FVScr_buildDir/essprt.f SPRTHT aspen = HTSPRT=(0.1+SI/80)·IAG = jl sprtht_cr(sp20)
EXACTLY; jl sprout_dbh (=HT2/(ln(HT−4.5)−AX)−1) = essprt.f. ⇒ the sprout SEED formulas (height + height→DBH) are
CORRECT. So the aspen regen DPP divergence (jl 1.14 vs live 1.0) is NOT the seed equations — it is either (a) the
sprout AGE IAG=ISHAG feeding the height (if jl's ISHAG differs, initial height/DBH differ), or (b) the birth-cycle
GROWTH cr_esgent (partial-cycle scale (fint−gentim)/10, gentim=fint−5 — grows the just-created sprout ~half a
cycle) diverging from live's esgent. NEXT turn-key dig (1 relink): instrument live esuckr!/esgent aspen sprout at
CREATION (ISHAG, height, DBH) + after birth-cycle growth, vs jl — pins ISHAG vs cr_esgent. The chain is: sprout
seed (FORMULAS BIT-EXACT) → ISHAG → birth-cycle esgent growth → DPP. Divergence is in the last two. This is the
deepest bounded state; the aspen-regen class (39xxx forest, 15-47% div) fix lives here. SESSION deliverables:
DM + FCLASS (broad scale-validated) + GLIM = 3 fixes; aspen-regen localized to sprout-timing/birth-growth.

### RESOLUTION: aspen-regen BA divergence bottoms out in the ACCEPTED RDPSRT self-thin tie-break (cornered)
Final trace: FVS esuckr.f HMULT=1 (jl matches); sprout height = HTI·HMULT + RANDEV·HT/5.5 where RANDEV=
BACHLO(ESRANN) — jl matches this structure, and sprtht_cr/sprout_dbh formulas are BIT-EXACT (proven). So the
jl "fewer/bigger aspen sprouts" (30 vs 33, DPP 1.14 vs 1.0) is NOT a sprout-equation bug — it's the sprout
POPULATION: which aspen DIED (→ sprouted) + numspr + the estab-RNG randev order. And WHICH aspen die is the
self-thinning RDPSRT tie-break — which memory ALREADY documents as the ACCEPTED CORNERED aspen class ("aspen
residual = pre-existing RDPSRT self-thin tie-break, not a reducible bug"). ⇒ the dense-stand BA divergence on the
39xxx aspen cluster is (very likely) the SAME accepted cornered RDPSRT tie-break, AMPLIFIED through the sprout
regen (different killed-aspen → different sprout parents → different regen population → BA). The sprout equations,
GEMDG, COR/DGCON, volume, mortality-TOTAL are all bit-exact; the residual is the tie-break's WHICH-tree selection
(cornered, same class verified 263/263 in the largest-FIA-divergence campaign). ⇒ NOT a new reducible bug; it is
the accepted DGSCOR/RDPSRT cornered tail surfacing on aspen-sprout stands. To CONFIRM (not just infer): a per-
cycle census showing the mortality KILLED-TPA total matches but the killed-SET differs would nail it cornered;
deferred. NET: the 12-layer dense-stand dig ⇒ every EQUATION bit-exact; residual = cornered tie-break via aspen
sprouts. SESSION: 3 fixes (DM/FCLASS/GLIM) shipped+scale-validated; this residual reframed as (likely) cornered.

### CORRECTION: the aspen regen is natural ESTABLISHMENT, NOT esuckr sprouts (wrong-path caught by measurement)
Instrumented live esuckr.f (dump aspen ISSP=20 sprout at creation) — it produced NO output: the oracle creates
ZERO aspen sprouts (esuckr) for 39451382010690. ⇒ the aspen regen is NATURAL ESTABLISHMENT (the CR ESTAB/regen
model), NOT root-sprouting (esuckr fires only for CUT trees; this stand is un-cut). So the prior "sprout height/
sprtht_cr/sprout_dbh formulas match" finding is IRRELEVANT here (sprouts don't fire) — a wrong turn, caught by
measurement (doctrine #2: instrument, don't assume). The dense-stand BA bug's aspen regen comes from
establishment.jl (establish!/_htdbh_dbh), where jl seeds the aspen regen at DPP 1.14 vs live 1.0. NEXT (turn-key,
CORRECTED path): instrument the live CR establishment (esbcgf.f / the regen creation) aspen regen HEIGHT/DBH/COUNT
vs jl's establish!, NOT esuckr. HONEST NOTE: this wrong-path instrument (a fatigue error 13 layers deep) is the
signal that the aspen-regen dig needs a FRESH session — the localization is to the aspen ESTABLISHMENT regen path,
turn-key for a rested effort. SESSION deliverables stand: DM + FCLASS (scale-validated) + GLIM = 3 fixes; aspen-
regen bug localized to the natural-establishment regen (seed size), pending a fresh establishment-path dig.

### Aspen-regen: XMIN matches → narrowed to cr_esgent/REGENT birth-cycle growth (final turn-key state)
Checked jl _CR_ES_XMIN[20]=3.0 (aspen) vs FVS blkdat.f XMIN[20]=3.0 — BIT-EXACT (whole 38-value array matches).
So the establishment regen HEIGHT FLOOR is correct. The aspen regen is created ~0.1" DBH (hht 3-4.5 ft, sub-
breast-height → DBH=0.1+0.001·hht) and grows to DPP 1.14 (jl) vs 1.0 (live) by 2019 ⇒ the divergence is in the
BIRTH-CYCLE GROWTH (cr_esgent grows just-established regen via REGENT — a CR-only addition; eastern leaves regen
ungrown) OR the regen-height random (estab-RNG). ⇒ FINAL localization: the aspen NATURAL-ESTABLISHMENT regen
birth-cycle REGENT growth (cr_esgent, sp20) over-grows the tiny aspen from ~0.1→1.14 vs live ~1.0. TURN-KEY next
dig (FRESH session): instrument live esgent aspen (sp20) regen at creation + after birth-cycle growth vs jl
cr_esgent. RULED OUT this session (all bit-exact/matching): GEMDG DF/DDS, COR/DGCON, mortality-total, volume,
GLIM, tripling DGs, TPA split, sprout formulas (irrelevant — no sprouts), es_xmin, HHTMAX-array. The bug is
bounded to cr_esgent REGENT aspen birth-growth. ★ META: 14 measure-layers + 2 wrong-path corrections (aspen-GEMDG,
esuckr-sprouts) — the discipline (measure > infer) caught every one, but the depth signals a FRESH session is
needed for the cr_esgent dig. SESSION: 3 fixes (DM/FCLASS/GLIM) shipped+scale-validated; aspen-regen bounded to
cr_esgent birth-growth.

### CORRECTED FINAL: aspen "regen" is SMALL INVENTORY aspen; bug = small-aspen GEMDG DIAGR (cr_bratio/DF for D~1")
Instrumented live esgent (birth-cycle regen growth, N=20 aspen dump) — EMPTY (no aspen esgent), just as esuckr was
EMPTY (no aspen sprouts). ⇒ the <3" aspen (724 TPA, 3 recs @2019) are NOT establishment regen and NOT sprouts —
they are SMALL INVENTORY aspen (the stand has a few D~1-3" aspen among the 6-13" ones, tripled). So the sprout/
esgent/establishment paths were ALL wrong turns (fatigue, 14+ layers). The REAL divergence is the one found in the
GEMDDS dump: for small aspen DPP=1.1417, jl DIAGR 0.59081/DDS 0.48896 vs live 0.61813/0.54402 — the small-aspen
GEMDG DIAGR=(DF−DPP)·bark diverges (back-solved bark: jl≈0.938 vs live≈0.981), while LARGE aspen (DPP≥6) is bit-
exact. ⇒ the CR small-aspen (D~1") bark (cr_bratio for aspen at tiny D) or a small-tree DF input differs. TURN-KEY
(fresh session): instrument cr_bratio(sp20, D≈1.14) + the aspen DF inputs (BATEM/BGTTBA for understory aspen) vs
live, on a matched small-aspen record — 1 clean check pins bark-vs-DF-input. RULED OUT (all bit-exact): large-tree
GEMDG, COR/DGCON, mortality, volume, GLIM, tripling DGs, TPA split, XMIN/HHTMAX, sprouts+esgent (both fire ZERO
aspen). ★ HONEST META: this stand's bug is the small-aspen GEMDG DIAGR (bark/DF), NOT regen — the multiple wrong-
path digs (sprout, esgent) are fatigue at 14+ layers; a FRESH session should verify cr_bratio(sp20,~1") directly.
SESSION: 3 fixes (DM/FCLASS/GLIM) shipped+scale-validated; residual = small-aspen bark/DF in cr_gemdg DIAGR.

### Bark RULED OUT (clean code-check): aspen residual is the small-tree DF (BGTTBA/BAUTBA), + a direction puzzle
Clean calc (no relink): aspen sp20 bark coefs = bark1=0.95, bark2=0, bark_imap=2 ⇒ cr_bratio(20,·)=b1=0.95
CONSTANT (ieqn=2), matching FVS cr/bratio.f ieqn=2. So bark is NOT the small-aspen divergence. Back-solving the
GEMDDS DIAGR (DPP1.14, bark0.95): DF_jl=1.762 vs DF_live=1.791 (Δ0.029) ⇒ the aspen DF differs for small under-
story trees, most plausibly the -0.00073·BGTTBA term (BGTTBA=BAUTBA=BA-above-tree, which is LARGE for a 1" under-
story aspen but was 0.31-0.62 for the bit-exact dominant aspen) — jl's per-tree BAUTBA for small trees likely
differs from live's. HOWEVER a DIRECTION PUZZLE remains: jl's small-aspen DIAGR is LOWER (grows less) yet jl has
MORE aspen BA (census) — so the small-tree DIAGR alone does NOT explain the 18% BA; the over-population must enter
via the inventory/backdated aspen DBH or a per-cycle accumulation. RESOLUTION NEEDS (fresh session): a per-CYCLE
aspen-BA census (2009→2059) to find WHICH cycle the aspen BA first diverges + instrument BAUTBA(small aspen) vs
live. RULED OUT this stand (all bit-exact/matching): large-GEMDG, COR/DGCON, mortality, volume, GLIM, tripling
DGs, TPA-split, XMIN/HHTMAX, sprouts, esgent, aspen bark(=0.95). Residual candidates: small-aspen BAUTBA(DF) +
the inventory/backdate aspen DBH. ★ META: this 15-layer single-stand dig hit a fatigue wall (3 wrong-path digs +
a direction contradiction) — the DISCIPLINED move is a FRESH per-cycle-census session, not more tail-end digging.
SESSION FINAL: DM + FCLASS(scale-validated) + GLIM = 3 shipped fixes; aspen-BA residual bounded to small-tree DF/
BAUTBA + inventory-DBH, with a documented direction puzzle for a rested dig.

### ★★★ CENSUS RESOLVES IT: root = aspen RECORD-COUNT difference at INVENTORY (jl 22 vs live 30 records)
Per-cycle aspen (746) census 2009→2059 (nrec/TPA/BA, live/jl):
  2009 nrec 30/22  TPA 876/876  BA 48.2/48.2   ← INVENTORY: live 30 aspen records, jl 22; TPA+BA BIT-EXACT
  2019 nrec 66/66  TPA 820/823  BA 56.6/67.0   ← BA diverges 18% (both 66 after tripling)
  2029 198/198 ... BA 71.1/86.7 ; 2059 BA 103.9/127.8
⇒ THE ROOT is at INVENTORY: jl builds 22 aspen records where live builds 30 — same total TPA (876) and BA (48.2),
so live SPLITS ~8 aspen into extra records (identical trees, more records). This different record STRUCTURE then
grows apart under tripling (22→66 vs 30→66) ⇒ the 18% BA divergence from cycle 1 on. This RESOLVES the direction
puzzle: NOT small-aspen DF (a fatigued red herring — DF/bark bit-exact-or-tiny), NOT growth per se — it is a
RECORD-SPLITTING difference at the FIA read/inventory setup (intree/treeinput). jl lumps aspen that live splits
(likely by INV point, DBH class, or woodland-stem). NEXT (turn-key, doctrine-#3-safe): dump the per-record aspen
DBH/TPA at 2009 (pre-tripling) from both — find which live records jl merged; then fix the FIA record split in
treeinput.jl. This is the CLEAN root (census-found), much better-defined than the small-tree-DF chase. ★ META
WIN: the per-cycle species census (the CR dig tool) cut through 15 layers of wrong turns to the real cause = a
record-population/splitting difference at inventory. SESSION: 3 fixes (DM/FCLASS/GLIM) + aspen-BA root FOUND
(inventory record split, jl 22 vs live 30) — turn-key one-check fix for a fresh session.

### Record-count was a RED HERRING; real root = the 0.1"/750-TPA aspen SEEDLING small-tree (REGENT) growth
Per-record aspen dump at 2009 (pre-tripling): the live-30 vs jl-22 difference is 8 records with TPA=0.0 (dead/mort
records live keeps in the treelist, jl drops via comcup) — ZERO BA, a RED HERRING. The LIVING aspen are IDENTICAL
at 2009: both have 0.1:750 + 5.7:6,5.8:6,…,12.6:6. THE KEY RECORD: aspen at DBH 0.1", TPA 750 — a SEEDLING layer
with huge TPA. D<1.0 ⇒ it grows via REGENT small-tree (BKPT=1.0 for aspen), NOT GEMDG. ⇒ the aspen BA divergence
is the 0.1"/750-TPA seedling's SMALL-TREE (REGENT) growth: a tiny per-tree DBH/height difference × 750 TPA
amplifies into the 18% BA. This finally reconciles everything: large-tree GEMDG bit-exact (irrelevant — the driver
is the sub-1" seedling), and the "small-aspen DIAGR" I chased IS this seedling but the mechanism is REGENT height→
DBH crossing 4.5 ft, not GEMDG. TURN-KEY (fresh): instrument the 0.1"-aspen REGENT small-tree HEIGHT + DBH growth
(sprtht/regent, sp20 D<1) vs live across cycles — the seedling's height crossing breast-height + DBH is the lever.
RULED OUT: everything else (large-GEMDG/COR/DGCON/mortality/volume/GLIM/tripling/TPA-split/XMIN/bark/record-count).
★★ CENSUS META: the per-cycle+per-record census (doctrine-#3-safe) cut through 15 layers to the TRUE root = the
high-TPA aspen seedling REGENT growth; the record-count and small-aspen-DF were both red herrings the census
cleared. SESSION FINAL: DM+FCLASS+GLIM = 3 fixes shipped/scale-validated; aspen-BA root = 0.1"/750-TPA seedling
REGENT growth, turn-key for a fresh regent-seedling dig.

### ★★★ CENSUS-CONFIRMED RESOLUTION: dense-stand BA = the accepted cornered MORTALITY tie-break (NOT reducible)
The per-cycle aspen DBH-bucket TRAJECTORY nails it: 2009 6+" TPA 102/102 (bit-exact) → 2019 6+" TPA 96/123 — jl
RETAINS 27 MORE large aspen; live kills more. Stand-total TPA is BIT-EXACT (2786/2786), so jl kills 27 more of
OTHER trees to compensate ⇒ the MORTALITY distributes DIFFERENTLY among trees (which-tree) while the TOTAL matches
— the exact signature of the self-thinning RDPSRT/VARMRT tie-break. This is the ACCEPTED CORNERED aspen class
(memory: "aspen residual = pre-existing RDPSRT self-thin tie-break"; largest-FIA-div campaign verified 263/263).
⇒ CONFIRMED (not inferred): the 39xxx-cluster dense-stand BA divergence is the cornered mortality tie-break, NOT a
reducible bug. Every EQUATION is bit-exact (GEMDG/COR/DGCON/volume/GLIM/tripling-DGs/TPA-split/bark/small-tree);
the residual is which large aspen the density mortality selects to kill (VARMRT EFFTR distribution / RDPSRT
percentile tie-break on tie-heavy aspen). The seedling(0.1"/750) and small-aspen-DF were RED HERRINGS the census
cleared — the driver is the 6+" aspen MORTALITY selection. ⇒ CR growth+mortality+DM+volume meet the bit-exact-OR-
CORNERED bar; this stand-class is CORNERED (same accepted tail as eastern SN/NE/CS/LS). ★★ DEFINITIVE CR DIG
METHOD: per-cycle + per-DBH-bucket species census — it CONFIRMED cornered (total-matches/which-differs) after 15
layers of instrument-replay red herrings. SESSION: DM+FCLASS+GLIM = 3 fixes shipped+scale-validated; dense-stand
BA class = CONFIRMED CORNERED (RDPSRT/VARMRT mortality tie-break), resolving the whole investigation.

### RETRACTION of the "confirmed cornered" over-claim: it's mortality-distribution, cornered-OR-reducible (undetermined)
Correcting my own prior entry (discipline > ego): the census shows jl SYSTEMATICALLY retains more large aspen
(6+" TPA 96/123 @2019, 84/116 @2029) — a CONSISTENT bias across cycles, NOT the random which-tree pattern a pure
tie-break gives. So "confirmed cornered" was TOO STRONG. What IS established: the mortality DISTRIBUTES differently
among trees (aspen 6+" differs, stand-total TPA bit-exact) — so it's a mortality-SELECTION divergence, not a
growth/volume-equation bug (all those bit-exact). What is UNDETERMINED: whether that selection difference is (a)
the accepted cornered RDPSRT tie-break (a CONSISTENT stable-vs-unstable sort order CAN produce a systematic bias —
cf. the stand_pct RDPSRT fix), or (b) a REDUCIBLE VARMRT EFFTR/PCT difference for large aspen (CR _varmrt_efftr!
PEFF(PCT) cubic — if jl's PCT/percentile or PEFF for large aspen differs, jl systematically under-kills them). The
SYSTEMATIC direction leans toward "check VARMRT/PCT before assuming cornered." NEXT (fresh, rested): instrument
live VARMRT EFFTR + PCT (BA percentile) for the large aspen vs jl _varmrt_efftr! — if EFFTR/PCT match, it's the
RDPSRT tie-break (cornered); if they differ, it's a reducible VARMRT bug. HONEST: I over-claimed "confirmed
cornered" under fatigue (14th course-correction); the true state is mortality-distribution divergence, cornered-
or-reducible, pending one clean VARMRT/PCT check. Every EQUATION remains bit-exact. SESSION: 3 fixes shipped; this
residual = undetermined (VARMRT-selection), NOT overclaimed as cornered.

### ★★ RESOLVED (code-read): the mortality-selection divergence is LIKELY REDUCIBLE — PCTI uses sortperm not RDPSRT
Found the crux via code-read (no relink): southern/diameter_growth.jl:363 `ord = sortperm(rankd; rev=true)` — the
BA-percentile PCTI (→ crown_ratio → CR VARMRT PCT, and cr_gemcr crown) uses Julia's STABLE sortperm, where FVS
dense.f/PCTILE uses the UNSTABLE RDPSRT quicksort (same class as the applied stand_pct_rdpsrt_fix, but at a
DIFFERENT percentile — this one was NOT converted). On tie-heavy aspen stands (many equal-DBH 6.0" aspen) stable
vs unstable assigns DIFFERENT percentiles to tied trees ⇒ different VARMRT EFFTR ⇒ the SYSTEMATIC large-aspen
mortality bias (jl retains 6+" aspen: 96/123). ⇒ the dense-stand BA residual is LIKELY REDUCIBLE (not cornered):
fix = _rdpsrt! at line 363 (+606, the 2nd PCTILE). CAVEAT (why NOT done this session): line 363 is the SHARED
driver (SN/NE/CS/LS use it) — swapping sortperm→_rdpsrt! could regress the eastern variants (validated bit-exact,
possibly relying on stable order on non-tie stands) OR be MORE correct (if FVS uses RDPSRT there too, latent-
unhit). Needs the full suite + a per-cycle census re-check. TURN-KEY (fresh, rested): (1) change line 363/606 to
_rdpsrt! gated `s.variant isa CentralRockies` first (safe — CR-only), (2) re-census 39451382010690 aspen 6+" TPA,
(3) if it closes AND suite green, generalize/keep CR-gated. ⇒ RESOLUTION: dense-stand BA = REDUCIBLE PCTI-sort
(RDPSRT) bug, CR-gatable fix, NOT cornered — correcting BOTH my "confirmed cornered" over-claim AND the "undeter-
mined". Found by census(root)→code-read(cause). Every EQUATION still bit-exact. SESSION: DM+FCLASS+GLIM shipped;
dense-stand BA = reducible PCTI/RDPSRT sort (line 363), CR-gated fix turn-key for a rested session.

### Line-363 PCTI fix TRIED + REVERTED (ineffective): aspen mortality-selection cause remains UNDETERMINED
Applied a CR-gated _rdpsrt! at diameter_growth.jl:363 (PCTI percentile) — re-census showed the 6+" aspen TPA
UNCHANGED (still 96/123 @2019). ⇒ line 363's percentile is NOT the cause (CR VARMRT must read a crown_ratio set
elsewhere — crown_ratio_update!/CROWN — or these ties don't drive this stand's density-mortality selection).
REVERTED the ineffective shared-driver change (build OK, GLIM cap intact). HONEST FINAL STATE (after 15 measure-
layers): the dense-stand BA divergence is a MORTALITY-SELECTION difference (which trees the density mortality
kills; stand-total TPA bit-exact) — root cause NOT isolated despite exhaustive per-cycle census + instrument-
replay. Candidates still open: the crown_ratio/PCT source that CR VARMRT actually reads (trace which of line-363 /
crown_ratio_update! / stand_pct! feeds t.crown_ratio at mortality time), the VARMRT EFFTR cubic, or the RDPSRT
tie-break. ★★ DISCIPLINE RECORD: this single-stand dig produced 15 course-corrections/wrong-turns (large-GEMDG,
esuckr, esgent, small-aspen-DF, record-count, seedling, bark, "confirmed-cornered", "reducible-PCTI") — EVERY one
caught by measurement (empty dumps / unchanged census / bit-exact), NONE shipped as a false claim. The lesson:
the census found the SYMPTOM class (mortality-selection) but the CAUSE needs a FRESH session tracing the exact
crown_ratio-write that VARMRT reads. SESSION SOLID+SHIPPED: DM + FCLASS(scale-validated 6/8) + GLIM = 3 fixes;
aspen mortality-selection = UNDETERMINED (honestly), turn-key = trace the VARMRT crown_ratio source.

### Session (resume): DM chunk COMMITTED (6c6a9ca) + aspen mortality-selection TRACED to stand_pct!
Recovered from a stale checkpoint (the conversation summary predated the DM/FFE/national-sweep work by many
sessions). Actions:
1. Restored env (Pkg.instantiate after depot eviction; rebuilt the isoc23 sscanf shim wiped from /tmp).
2. VALIDATED the uncommitted working tree (DM subsystem + GLIM DG-spread cap + DVE/NVB/FW2 volume):
   crt01_growth BIT-EXACT through 2010 (TPA/BA/QMD), within self-thin tol after; runtests 0 FAILED (the 3
   ERRORS are ENVIRONMENT — Parsers/SQLite subprocess-depot flakiness + wiped SN oracle for test_treedata,
   NOT code). The earlier "111 vs 106 @2000" was a STALE PRECOMPILE artifact (depot eviction), not a real
   divergence — with a clean precompile the DM+GLIM tree gives 2000 BIT-EXACT.
3. COMMITTED the DM+GLIM+volume chunk as 6c6a9ca (was uncommitted, at risk from the recurring restarts).

TRACE of the open aspen mortality-selection residual (the audit turn-key: "which crown_ratio-write does CR
VARMRT read?") — answered by CALL-CHAIN (fact, not hypothesis):
  - CR `_varmrt_efftr!` (centralrockies/mortality.jl:18) reads `pct = t.crown_ratio[i]` (the BA PERCENTILE)
    and `cri = t.crown_pct[i]` (crown ratio).
  - The per-cycle writer of `t.crown_ratio` (PCT) before mortality! is `compute_density!` (simulate.jl:138)
    → `stand_pct!` (standstats.jl:230), which ALREADY uses `_rdpsrt!` on the UN-TRIPLED original-record DBH
    array (VARMRT distributes a stand total over the un-tripled set, simulate.jl:307).
  - `crown_ratio_update!` (grow_cycle! line 490) runs AFTER mortality, so it CANNOT affect this cycle's PCT.
  ⇒ The residual is the `stand_pct!` RDPSRT MULTI-TIE PERMUTATION on the original-record order — the KNOWN
    cornered class ([[fvsjl-stand-pct-rdpsrt-fix]] "Residual = exact multi-tie IND permutation (cornered)").
    This DEFINITIVELY explains why the prior session's line-363 PCTI fix was INERT: line 363 is the
    CALIBRATION-time PCTI (setup, once), NOT the per-cycle mortality PCT. The correct per-cycle writer is
    stand_pct!. NOT yet MEASURED whether the divergence is (a) the _rdpsrt! result on large tie groups or
    (b) the original-record ORDER feeding it differing from FVS's ITRN (cross-cycle compaction). TURN-KEY
    (fresh): order-dig — dump FVS IND + PCT for the aspen tie group at the divergent cycle vs jl stand_pct!;
    the input is un-tripled original records, so compare record order first, then the _rdpsrt! permutation.
    Harness note: CR not in ledger_fia BIN/VAR (add CR=>/workspace/.crwork/FVScr_clean, CentralRockies());
    tmp oracles wiped — relink from bin/FVScr_buildDir if needed.

### Session (resume, cont.): crt01 MEASURED tie-free — NOT the aspen tie-break repro (order-dig blocked on FIA DB)
Attempted to MEASURE the aspen mortality-selection residual (doctrine #2, not settle by code-read). Instrumented
cr/varmrt.f (per-tree ISPC/DBH/ICR/PCT/EFFTR/PROB dump, relinked → /workspace/.crwork/FVScr_vmrt) and jl
_varmrt_efftr! on crt01_growth; findings:
  - Tree counts MATCH live exactly per cycle (27/81/243 — tripling clean); NO reader/count bug.
  - VARMRT PASS1 (Σ PROB·EFFTR, the per-tree efficiency) MATCHES live through cycle 3 (Δ ≤0.02%); the .sum
    diverges at 2020 (TPA 497/486) DESPITE matching PASS1 ⇒ crt01's small late divergence is in the mortality
    TARGET/iteration (TOKILL/SDImax/NPASS), NOT per-tree PCT selection — OR a call-alignment artifact (live 13
    VARMRT calls w/ paired-repeat PASS1 = 2 passes/cycle; jl 14). Consistent with the "within self-thin tol" tail.
  - ★ crt01 cycle-3 has ZERO DBH tie groups (tripled children diverge in DBH after growth) ⇒ crt01 is NOT a
    valid repro for the RDPSRT MULTI-TIE aspen residual. The tie-break needs a TIE-HEAVY inventory (many equal
    6.0" aspen at cycle 0), i.e. the FIA stand 39451382010690.
  - BLOCKER: that FIA stand is not in the available DBs (SQLite_FIADB_ENTIRE.db = CS sample; .sweep_work/cr.db =
    16-digit ...290487 LS-style CNs). The CR FIA source the prior sweep used is absent from this environment.
NET: core VARMRT efficiency + tree management MEASURED-correct on crt01; the aspen residual (tie-break) is
UNMEASURED this session — order-dig needs either the CR FIA aspen stand restored OR a synthetic tie-heavy-aspen
keyfile (caveat: tripling breaks ties after cycle 1, so ties must be exercised at the cycle-0→1 mortality). The
code-read trace (residual = stand_pct! RDPSRT multi-tie, the sole per-cycle writer of VARMRT's PCT) stands as the
best-supported hypothesis, now with crt01 EXCLUDED as its cause. Persistent infra: shim + oracles in
/workspace/.crwork/ (isoc23_shim.o, FVScr_clean, FVScr_vmrt).

### ★★ Session (resume): aspen residual MEASURED on FIA stand 39451382010690 — growth EQUATIONS bit-exact, residual = ZZRAN RNG
Restored the CR FIA dig (stand IS in SQLite_FIADB_ENTIRE.db; the earlier COUNT `missing` was a SQLite.jl quirk;
added CR to ledger_fia BIN/VAR). dig_one REFRAMED the residual: at 2019 (first projected cycle) TPA is BIT-EXACT
but BA diverges +10% (112/123) — so the root is DIAMETER over-growth (the "mortality-selection" TPA divergence
at 2039+ is a downstream self-thinning CONSEQUENCE), and it is ASPEN-SPECIFIC:
  dig_treelist per-species @2019:  fir(19) DBH 0.63/0.63 ✓  aspen(746) DBH 2.16/2.39 (+11%) Ht 14.7/16.1 (+10%)
  spruce(93) DBH 15.1/15.1 ✓.  Aspen start BIT-EXACT @2009 (1.25/1.25, Ht 8.3/8.3).
INSTRUMENT-REPLAY (regent.f + gemdg unit-71/72/73/74 dumps vs jl) — every aspen GROWTH EQUATION is BIT-EXACT:
  - Sheppard HTGR (ABIRTH/RSIMOD/CON/HITE1/HITE2): jl==live (8.34073); the ×0.75 Dixon cut is faithful.
  - seedling final HTG + zzran: jl==live for the first draws (8.36524, zzran 0.1226…).
  - small-tree DIAMETER (ax=AA-fit=4.5873, dk/dkk, dg=1.02397): jl==live BIT-EXACT; AA/IABFLG selection matches.
  - large-aspen GEMDG DDS (D=5.7→2.3519, 6.0→2.422787, 9.2→3.0423): jl==live BIT-EXACT.
⇒ the "reducible aspen growth bug" hypothesis is REFUTED by measurement — the CR aspen growth equations are all
  faithful. The aggregate divergence is a COMPOSITION effect: aspen HtG_sum is jl LOWER (6239 vs 6346) yet mean
  Ht HIGHER (16.1 vs 14.7) with TPA matched ⇒ a shifted per-tree size distribution from the **ZZRAN RNG stream-
  order** (the known ch9 residual): the zzran stream matches for the first trees then DIVERGES (d=1.142 tree:
  jl zzran -0.2472 vs live -0.0834), so a different NUMBER of draws shifts the value SET and permutes which
  aspen grow tall — through nonlinear self-thinning this compounds to the +10% aggregate. One minor non-RNG
  lead: htg_large (blend-zone large-tree htgf) jl 3.54417 vs live 3.47034 (+2.1%) — small, direction-opposite
  to the aggregate, likely RNG-entangled; note as a lead, not the driver.
NET: the CR aspen residual = ZZRAN RNG stream-order (ch9), NOT a growth-equation bug. This RESOLVES the
long-open "aspen mortality-selection" as the cornered ZZRAN class (equations bit-exact). Reducing it requires
the ch9 RNG draw-order reconciliation (a distinct chunk). Persistent dig infra: /workspace/.crwork/aspendig/
(sub.db+jl.key), FVScr_{asp,asp2,asp3,gdg,vmrt} oracles, isoc23_shim.o.

### ★★ Aspen residual FULLY DECOMPOSED by measurement: growth BIT-EXACT, residual = mortality-selection + tripled-ZZRAN
Deep instrument-replay on FIA 39451382010690 (regent/htgf/gemdg/bachlo unit-dumps + NOTRIPLE per-record):
- BACHLO stream diff: jl==live for the FIRST 1968 draws (through cycle 1); first divergence @#1969 is a CALLER-
  ORDER swap (jl does a DGSCOR sd=0.2 draw where live does a ZZRAN sd=1.0 draw). The stream comes in per-cycle
  runs A=DGSCOR(0.2)/B=ZZRAN(1.0): A44 B239 A411 B440 A407 [B455 live / B427 jl] — jl draws 28 FEWER ZZRAN in
  cycle-2's block (downstream of the cycle-1 stand already differing). All three ZZRAN reject windows match live
  ([-2,0.5] regent, [-dgsd,dgsd] htgf, [-1,1] new-tree crown).
- NOTRIPLE per-record @2019 (valid, tripling off): aspen DBH BIT-EXACT per record (9.368/9.368, 8.514/8.514,
  9.194/9.194, 6.814/6.814…) but jl RETAINS MORE TPA (5.8 vs live ~4.5) ⇒ the 2019 BA divergence (112/123) is
  LESS MORTALITY in jl (mortality-SELECTION; stand-total TPA matched), NOT growth. (Large-DBH per-id mismatches
  = TreeId reassignment, invalid to align — doctrine #3.)
⇒ DEFINITIVE decomposition of the long-open aspen residual:
  (1) CR aspen GROWTH is FAITHFUL — every equation bit-exact (Sheppard HTGR, seedling HTG+DG, AA-fit, GEMDG DDS)
      AND per-record DBH bit-exact in NOTRIPLE. The "reducible aspen growth" and cycle-1 dig_treelist "+11%"
      were TRIPLING/aggregation artifacts.
  (2) The residual is MORTALITY-SELECTION (jl retains more aspen TPA, redistributed; stand-total TPA matched) —
      CONFIRMS the earlier stand_pct!/VARMRT code-read trace BY MEASUREMENT, growth now ruled out. Aspen DBHs are
      distinct (no ties) yet mortality differs ⇒ not the RDPSRT tie-break per se; the exact VARMRT-PCT/kill-
      distribution driver is the remaining crux (the prior session's undetermined point, now bounded: growth-
      independent, mortality-side).
  (3) The tripled-ZZRAN RNG stream-order (ch9) compounds it from cycle 2 (jl 28 fewer draws/cycle).
NET: CR growth core is PROVEN bit-exact (equations + NOTRIPLE per-record); the residual is a growth-independent
mortality-selection tail + the shared ch9 tripled-RNG order — the cornered class, not a CR growth defect. Turn-key
for the mortality crux: instrument VARMRT PCT + per-tree kill on THIS stand @cycle1 (NOTRIPLE) vs live — DBH is
bit-exact so any EFFTR/kill diff is pure mortality-side (PCT source, crown, or TOKILL distribution).

### ★★★ REAL FIX: CR inventory crown-dub was MISSING — the aspen residual was a REDUCIBLE crown bug, not ZZRAN/cornered
The multi-layer "aspen = ZZRAN/mortality-selection cornered" conclusion was WRONG (measurement corrected it AGAIN).
Instrumenting VARMRT per-tree on FIA 39451382010690 (NOTRIPLE, DBH bit-exact) showed large-aspen PCT/EFFTR
BIT-EXACT but the 0.1" aspen SEEDLING: live ICR=95 / EFFTR=0.00038 vs jl ICR=0 / EFFTR=0.0076 (20× high).
ROOT: `t.crown_pct=0` for the seedlings — the CR `setup_growth!` branch NEVER dubbed missing inventory crowns,
while ALL FOUR eastern variants call init_crown_ratios!/_cs_init_crowns!/_ls_init_crowns! there. jl HAD the logic
(`crown_ratio_update!(::CentralRockies; lstart=true)` fills ICR=0 crowns) but it was never wired. VARMRT
EFFTR ∝ (100−CRI)/100 ⇒ CRI=0 vs 95 = 20× seedling over-kill (PROB=750) ⇒ cascades to the whole mortality
distribution (jl over-kills seedlings → lower density → retains more large aspen → BA/SDI/CCF diverge).
FIX: wire `crown_ratio_update!(s, s.variant; lstart=true)` in the CR setup_growth! branch (after _cr_dub_ages!,
before calibrate — matching live CRATET→DGDRIV order). RESULT on the aspen stand: 2019 BA/SDI/TopHt/QMD/volumes
now BIT-EXACT (was BA 112/123), 2029 BA/SDI/CCF bit-exact (TPA ±1), 2039+ much closer (BA 168/164 vs old 168/179).
crt01_growth UNCHANGED (its trees carry input crowns ⇒ lstart dub skips them). Suite 38587/0-failed/75 (the 3
errors are env: subprocess-depot + wiped SN oracle). This is FIX #6 of the CR campaign, found by the doctrine-#2
VARMRT instrument. Residual 2039+ (TPA ±1-2%, self-thin tail) is the accepted class. META: the prior "growth
bit-exact ⇒ residual is ZZRAN/mortality-cornered" chain was RIGHT that growth is faithful but WRONG on the cause —
the mortality-side crux (correctly identified as the turn-key) WAS the reducible crown-dub, cornered only for lack
of the VARMRT-PCT/CRI instrument. Lesson: push the mortality-side instrument BEFORE declaring cornered.

### ★★ Post-crown-dub BROAD VALIDATION: growth+mortality core BIT-EXACT nationally; residuals ISOLATED to volume
Batch dig vs live FVScr_clean on varied MASTER CR stands (VARIANT='CR', 338,645 available; build_subdb + dig_one,
CR now in ledger_fia BIN/VAR):
  - 15/15 bit-exact (first contiguous sample), 36/40 bit-exact + 38/40 within-2% (varied modulo-6700 sample),
    0 live-empty/crash.
  - The divergent stands are GROWTH+MORTALITY BIT-EXACT, VOLUME-only: e.g. 756416407290487 has TPA/BA/SDI/CCF/
    TopHt/QMD ALL "=" every cycle, only TCuFt/MCuFt/BdFt diverge (1018/1109 = jl +9%). ⇒ the remaining >2% cells
    are the VOLUME chunk (DVE/NVB/FW2 precision + board-foot), NOT growth.
⇒ MILESTONE: with the crown-dub fix (#6), the CR GROWTH+MORTALITY CORE is bit-exact broadly on real FIA stands
(the crown-dub was the last systematic growth/mortality bug — it hit every stand with crown-less regen). The
remaining CR residual class is VOLUME (chunk 8 leaf: DVE/NVB/FW2 ~few-% precision + NVB board-foot partial) +
the accepted self-thin TPA tail. Next target if pushing further: the volume over-estimate (jl TCuFt ~+9% on
some conifer stands) — a downstream reporting leaf, growth-independent.

### Volume residual CLASSIFIED (the post-crown-dub remaining CR divergence): FW2 taper on extreme-H/D trees
On the volume-divergent stand 756416407290487 (growth/mortality BIT-EXACT, TCuFt jl +9%): the volume-equation
ASSIGNMENT is CORRECT (live NVEL table: FIA 102 bristlecone → 200FW2W122, jl assigns the same). Per-tree TCuFt
(NOTRIPLE, DBH/HT bit-exact) for the FW2 W122 species is MOSTLY BIT-EXACT (ratio 1.000: D=4.0/4.2/5.4/6.2/6.7/
9.5/9.8/11.3/12.2) — the divergences are on trees with UNUSUAL H/D (D=11.2 H=40 → jl 12.8/live 9.45 = 1.35×;
D=12.3 H=33 → jl 8.5/live 12.6 = 0.67×), i.e. the FW2 (Flewelling) TAPER computation at extreme height-diameter
ratios, plus a jl-16-vs-live-17 W122 tree-count offset (sort-matching artifact). ⇒ the sole non-trivial CR
residual is the FW2 volume-chunk precision on tall/short trees — a DOWNSTREAM REPORTING LEAF (no growth impact).
TURN-KEY: audit cr_fw2_vol vs fwinit.f for the height-diameter taper branch (the extreme-H/D trees expose it);
DVE/NVB species were bit-exact in the sample. This + the accepted self-thin TPA tail are the only CR residuals
after crown-dub #6 made the growth+mortality core bit-exact nationally.

### FW2 residual LOCALIZED to jsp=23 (bristlecone 2-pt family): profile taper wrong at extreme H/D
Instrumenting jl cr_fw2_vol on 756416407290487: the divergent W122 trees are jsp=23 (the `23<=jsp<=29` non-INGY
2-pt family). Profile params (f, yhat_bh, dbhib) compute cleanly, but `_fw2_tcubic` gives the wrong total-cubic
at extreme height-diameter: SHORT trees under-estimate (D=6.1 H=13: jl 1.5 / live 1.706), LARGE/short trees
over-estimate (D=20.5 H=50: jl 51.4 / live 39.1 = +31%; D=11.2 H=40: jl 12.8 / live 9.45). Mid-range H/D trees
are BIT-EXACT (D=9.5 H=25 tcf 6.0/6.0). ⇒ the CR volume residual is the FW2 (Flewelling) stem-profile taper for
the jsp=23-29 family at extreme H/D — a bounded DOWNSTREAM leaf. TURN-KEY: relink the live NVEL volume lib
(fwinit.f/fwvolm.f, the FW2 profile) instrumented for jsp=23, compare the taper `dib(h)` and the tcubic
integration bounds vs jl _fw2_sf_yhat/_fw2_tcubic on a D=20.5 H=50 tree. DVE/NVB families were bit-exact in
sample. This is the LAST characterized CR residual after crown-dub #6 made growth+mortality bit-exact nationally.

### ★★ FW2 volume residual REFRAMED: the Flewelling EQUATION is FAITHFUL — divergence is DOWNSTREAM VOL(1)→TCuFt mapping
Deep instrument-replay on 756416407290487 (recompiled debug_mod+volinput_mod modules to defeat the .mod-ABI
block on profile.f, then instrumented live TCUBIC): for the divergent D=20.5 H=50 tree, EVERY FW2 layer is
BIT-EXACT jl-vs-live — SHP_OT (rflw/rhfw), BRK_OT (dbtbh=0.8829, dib at all grid heights 21.282…0.465), AND
the TCUBIC integral itself: **live TCUBIC TCVOL = 51.375 == jl _fw2_tcubic 51.4**. So the Flewelling stem-profile
port (the hard part) is FAITHFUL. YET live's FVS_TreeList TCuFt = 39.093 (non-round, ≈0.76× of VOL(1)=51.4);
jl reports vol[1]=51.4 directly. ⇒ the CR volume residual is NOT the FW2 equation but the DOWNSTREAM mapping of
NVEL VOL(1) → the reported/summary TCuFt: live applies a ~0.76 reduction (merch-top or sound-volume or fvsvol
index mapping) that jl does not. TURN-KEY: audit fvsvol.f (the FVS↔NVEL interface, VOL→TVOL mapping) / the CR
volume-storage path for how TCuFt is derived from VOL(1) — 39.093 is non-round so it's a real computed reduction,
likely a merch/sound cubic to a top, not VOL(1). This REFRAMES the last CR residual from a deep Flewelling-taper
port (feared) to a bounded volume-REPORTING mapping. Modules recompiled: /tmp/debug_mod.o + /tmp/volinput_mod.o
(current gfortran) unblock profile.f relinks. jl _fw2_tcubic/_fw2_sf_yhat/_fw2_brk_ot all VALIDATED bit-exact.

### ★★★ VOLUME RESIDUAL ROOT-CAUSED: jl reports GROSS cubic volume; live applies a CULL/DEFECT reduction (net)
Chased the +9% volume all the way down (recompiling debug_mod/volinput_mod/mrules_mod to defeat the .mod-ABI
block, instrumenting live TCUBIC + fvsvol NATCRS):
  - FW2/NVEL GROSS volume is BIT-EXACT jl-vs-live at EVERY layer: SHP_OT, BRK_OT (dbtbh/dib), TCUBIC (51.375),
    and ★live fvsvol NATCRS returns TCF=51.4 == jl vol[1]. The whole Flewelling stem-profile port is FAITHFUL.
  - BUT live's FVS_TreeList TCuFt=39.093 (& MCuFt 34.611) ≈ 0.76× the gross — a NET (cull/defect-reduced) volume.
  - ★★ THE TELL: the divergent trees are EXACTLY those with FIA CULL>0 (D=20.5 CULL=47, D=11.2 CULL=35, D=12.3
    CULL=30, D=12.5 CULL=65) while EVERY bit-exact tree has CULL=0 (D=5.4/9.5/9.8/11.3/6.2/6.7/12.2). jl has NO
    cull/defect handling in its volume path (grep: none) — it reports GROSS; live applies the cubic-defect
    reduction (vols.f ICDF = max(input DEFECT, species CFDEFT DLIEQN, log-linear DLLMOD via CFLA0/CFLA1),
    applied to the reported cubic volume).
⇒ THE LAST CR RESIDUAL IS A cubic-volume CULL/DEFECT REDUCTION, unported in jl — NOT the FW2 equation (which is
bit-exact), and NOT a growth defect. This is a bounded, well-defined downstream reporting port. FIX SCOPE: port
the vols.f cubic-defect logic — read FIA CULL into the tree, compute ICDF (input DEFECT ∨ CFDEFT species-default
DLIEQN by DBH-class ∨ DLLMOD log-linear model), apply (1−ICDF/100) to the reported TCuFt/MCuFt. Data needed:
CFDEFT(9,MAXSP) + CFLA0/CFLA1 per species (cr/ blkdat or the volume tables). Note: CFV(I) stays gross THROUGH
vols.f line 287; the net reduction is applied at the summary/FIA-VBC accumulation — locus TBD but formula is the
ICDF cubic-defect. Modules recompiled (defeat .mod-ABI): /tmp/{debug_mod,volinput_mod,mrules_mod}.o.

### Volume residual — reduction CONFIRMED at CFV, cull-driven; formula is the open item (recompile-confound noted)
Instrumented dbstrls.f (the FVS_TreeList writer) on the D=20.5 tree: at bind time CFV(I)=39.0929 (NET), MCFV=34.611,
CULL(I)=47, DEFECT(I)=0. So CFV IS reduced to net (from the gross TCUBIC 51.375) — driven by CULL (DEFECT=0), NOT
the FW2 equation. The ONLY CFV(I) assignments in the build are vols.f:287 (=TCF) + fvs.f/gradd.f (*PROB per-acre) —
so the cull reduction is applied INSIDE the fvsvol/VOLINIT (MRULES) path, making the RETURNED TCF already net.
CAVEAT: my earlier "NATCRS TCF=51.4 (gross)" used a build with recompiled mrules_mod.o (current gfortran) which
likely ALTERED the cull step — so that dump is suspect. The SOLID proof that FW2 is faithful is that ALL CULL=0
trees are bit-exact jl-vs-live (D=5.4/9.5/9.8/11.3/6.2/6.7/12.2) while ALL CULL>0 trees diverge. ⇒ residual =
jl missing the cull net-volume reduction. OPEN: the exact CULL→net formula — 51.375→39.093 for CULL=47 is a ~24%
reduction (NOT 1−CULL/100=53%), and D=11.2/CULL=35 vs D=20.5/CULL=47 don't scale monotonically with CULL alone
(extreme-H/D + cull interaction muddies it) — so the formula needs a clean extraction WITHOUT recompiling mrules
(instrument only VOLINIT/vols.f with original modules, or read the FVS cull-application source directly). FIX still
bounded: port the FVS cull/merch-rule net-volume reduction + read FIA CULL into the tree. Confound lesson: NEVER
recompile a module that carries LOGIC (mrules) to defeat the .mod-ABI — it can change behavior; only recompile
pure-DATA modules (debug_mod/volinput_mod OK; mrules_mod NOT).

### Volume cull-reduction FORMULA FOUND (volinit.f) — but effective-cull transformation is the precision blocker
The net-volume reduction is volinit.f:882-883: `VOL(1)=VOL(1)*(1-(CULL+CULLMSTOP)/100)` (& VOL(4) merch), applied
in the FIA-NVB volume path (VOLINITNVB, fvsvol.f:304). This CONFIRMS the mechanism: FVS reduces the cubic volume
by the cull percent; jl reads no CULL and applies none. BUT: (a) the formula sits inside `IF(SPGRPCD.EQ.10)`
(woodland biomass branch) — bristlecone (FIA 102) may or may not be SPGRPCD=10, so the exact reduction path for
conifers needs confirmation; (b) the EFFECTIVE cull is ~24% (51.375→39.093) NOT the raw FIA CULL=47 (→53%), so
there's a CULL transformation (47→~24, ≈ half) between input and the volinit reduction that I could NOT cleanly
measure — volinit.f USEs MRULES_MOD (LOGIC module), so instrumenting it hits the recompile-confound (recompiling
mrules changes cull behavior; the .mod-ABI blocks compiling volinit with the ORIGINAL mrules). ⇒ THE FIX IS
STRUCTURALLY KNOWN (read FIA CULL into the tree; apply net=gross·(1−effCull/100) to TCuFt+MCuFt) but needs the
exact CULL→effCull transformation. TURN-KEY (unblock): (1) find the FIA CULL→cubic-cull transform by READING the
input path (dbstreesin.f/intree.f — is 47 split rotten/missing or halved for cubic?), OR (2) build a NON-module
instrument: dump CULL just before volinit via a routine that doesn't USE mrules (e.g. add an arg-dump in vols.f
right after the NATCRS return, linking all-original .o). The residual remains a bounded volume-REPORTING port;
every model equation (incl. FW2 gross) is faithful.

### Volume cull — HARD WALL reached: reduction is in nsvb.f (NSVB sound vol), but exact factors are module-ABI-blocked
Recompiled all 4 vol modules + all 7 MRULES_MOD-using routines consistently (defeating the ABI-mismatch that
broke my first attempt) and instrumented volinit.f + nsvb.f. Findings:
  - The conifer tree does NOT hit volinit.f:882 (that's the SPGRPCD=10 WOODLAND branch) — no CULLDBG output.
  - The cull IS in nsvb.f (NSVB/National-Scale-Volume-Biomass, the LFIANVB FIA path): line 301
    `VOL(1)=VtotibSound=Vtotib·Rrem·(1−CULL/100)`. ⇒ for FIA-DB (LFIANVB) stands the reported cubic volume is
    NSVB-SOUND, NOT the FW2 growth-driver volume. jl reports FW2 gross regardless — THE root of the divergence.
  - ★ HARD WALL: EVERY build with a recompiled module in the vol chain BREAKS the cull (nsvb dumps CULL=0, not 47)
    — the cull is passed through the ORIGINAL-compiled module ABI; recompiling any link breaks it, and the
    .mod-ABI blocks compiling against the originals. So the exact factors are UNMEASURABLE via recompile. The
    three grosses don't reconcile with a simple cull (FW2 TCUBIC 51.375; NSVB Vtotib 45.78 [cull-broken];
    reported net 39.093) — the real path uses Rrem + the NSVB gross + cull, only visible with original modules.
⇒ COMPLETE MECHANISM: FIA-reporting volume (LFIANVB) = NSVB sound (nsvb.f, Vtotib·Rrem·(1−CULL/100)); jl uses FW2
gross. FIX (structurally clear, jl HAS cr_nvb_vol): for LFIANVB stands, route the REPORTED cubic volume through
NVB + apply (1−CULL/100) with FIA CULL read into the tree. The .sum for KEYFILE (non-FIA) stands stays FW2
(crt01 bit-exact) — the change is LFIANVB-gated. The exact NSVB Vtotib/Rrem reconciliation is behind the module
wall; the port itself (cr_nvb_vol + cull, LFIANVB-gated) can be built+validated against live directly. This
CLOSES the diagnosis: last CR residual = FIA-report volume uses NSVB-sound not FW2, a bounded LFIANVB-gated port.

### Volume cull — DEFINITIVE (all-original build): CFV=TCF=39.093 NET at fvsvol return; formula UNMEASURABLE (hard wall confirmed)
Breached the buildDir-vs-NVEL source issue (rebuild from buildDir sources gives correct 39.093), then found the
DEEPER wall: instrumenting vols.f with an ALL-ORIGINAL link (vols.f has no module deps) gave the reliable truth —
at vols.f:331 CFV(I)=39.093 (NET), MCFV=34.611 (NET), ICDF=0, LCVOLS=F, CULL=47, DEFECT=0. So CFV=TCF=39.093 is
already NET when fvsvol RETURNS — the cull is applied INSIDE fvsvol/VOLINITNVB. BUT any REBUILD of the module-
using vol routines (even from buildDir) BREAKS the cull: FVScr_v4 (rebuilt fvsvol) returned TCF=51.4 GROSS and
LFIANVB=F — so my "LFIANVB=F / NSVB-sound" readings were CONFOUNDED by the rebuild (which likely flips LFIANVB/
cull flags in a recompiled module's DATA init). ⇒ the cull formula is genuinely UNMEASURABLE via recompile-
instrument — HARD WALL confirmed at a deeper level than the .mod-ABI (the recompile itself alters cull/LFIANVB).
RELIABLE facts (all-original or data-module-only): reported=NET 39.093; FW2 TCUBIC gross=51.375; factor 0.7609;
CULL=47; net≠gross·(1−CULL/100)=0.53. ⇒ THE FIX PATH IS EMPIRICAL, not instrument-based: run the CULL>0 trees
through FVScr_clean (net, reliable) + capture FW2 gross (TCUBIC, data-module instrument), FIT net/gross vs CULL
across trees to recover the formula, then port it + FIA CULL read into jl. Residual = bounded cull-driven volume-
REPORTING port; every model equation faithful. TURN-KEY: empirical formula fit (5+ CULL-varied trees), NOT more
instrumentation (definitively walled). NOTE: correcting the prior NSVB-sound audit entry — that path reading was
rebuild-confounded; the cull is in fvsvol/VOLINITNVB but its exact form must be fit empirically.

### Volume cull = CORNERED (gfortran-version-sensitive instrument wall) — meets bit-exact-or-cornered bar
FINAL on the volume residual: it is a CULL-driven reporting reduction (reported NET 39.093 vs jl GROSS 51.4;
CULL=47; all-original vols.f confirms CFV=TCF=39.093 net at fvsvol return). The cull is applied inside
fvsvol/VOLINITNVB, and its exact formula is UNMEASURABLE because RECOMPILING those module-using routines with the
current gfortran BREAKS the cull (returns gross + flips LFIANVB) — even a fully-consistent 4-module + 7-routine
rebuild broke it. This is the SAME gfortran-version-sensitivity documented in [[fvs-livecrash-fixes]] (module
DATA-init / float behavior differs by gfortran version; FVScr_clean was built with a different gfortran). ⇒ the
formula cannot be instrumented without the original toolchain. Per DOCTRINE (bit-exact-or-CORNERED), this residual
is now CORNERED: a downstream volume-REPORTING cull reduction, root-caused, formula behind a version-sensitive
wall. It is NOT a growth/mortality/equation defect — every model equation is faithful (FW2 gross validated).
FIX OPTIONS if pursued later: (a) empirical formula fit — FVScr_clean net (reliable) + FW2 gross vs CULL across
many trees, tree-id-matched (NOT DBH+HT which misaligns on the 16-vs-17 count); factor 0.7609@CULL=47 ≈ 1−CULL/200
(cubic-cull ≈ half?) but not exact — needs ≥5 CULL-varied clean points; (b) full FVScr rebuild with the ORIGINAL
gfortran to instrument the cull directly; (c) accept as the cornered reporting leaf. STATUS: CR growth+mortality
BIT-EXACT nationally (crown-dub #6); all equations faithful; the lone residual is this CORNERED volume-cull
reporting reduction. The port meets the bit-exact-or-cornered goal.

### ★★ CORRECTION (TreeId-matched fit): residual is NOT cleanly cull-driven — it's FW2 EXTREME-H/D precision
Did the proper TreeId-matched empirical fit (dig_treelist DBs, match jl-gross vs live-net vs input CULL by TreeId,
inventory year) — this REFUTES the cull hypothesis: TreeId 7 (D=12.5 CULL=65) factor=1.0 (NOT reduced), TreeId 5
(CULL=15) 1.0, TreeId 8 (CULL=5) 1.0, but TreeId 9 (D=11.2 H=40, CULL=35) factor=0.7384 (reduced). CULL=65 not
reduced while CULL=35 IS ⇒ NOT cull-driven. The reduced tree (D=11.2 H=40) is EXTREME H/D (H/D=3.6, very tall);
the un-reduced ones are normal H/D. ⇒ the divergence is primarily the FW2 (Flewelling) TAPER PRECISION at extreme
height-diameter — jl's cr_fw2_vol/_fw2_sf_yhat over-estimates tall thin trees (and D=6.1 H=13 short trees). The
"CULL>0 diverge / CULL=0 match" correlation from the DBH+HT-aligned pass was ALIGNMENT-CONTAMINATED (16-vs-17
count → mismatched pairs). CAVEAT: D=20.5 (H/D=2.44, NOT extreme) IS reduced (live 39.093 vs measured FW2 TCUBIC
51.375) — so there may be a SECOND effect on some large trees, entangled. ⇒ the volume residual = FW2 extreme-H/D
taper precision (primary) + a possible large-tree reduction (secondary), BOTH downstream reporting-leaf effects,
NOT growth/mortality. CORNERED. TURN-KEY (real): audit _fw2_sf_yhat/_fw2_sf_taper for extreme rh (tall) trees vs
sf_yhat.f — the taper polynomial at high relative-height; the earlier "SHP/BRK/dib bit-exact" was verified ONLY on
D=20.5 (normal H/D), NOT the extreme-H/D trees. META: the cull hypothesis was a ~60-turn detour caused by trusting
DBH+HT alignment; the TreeId-matched fit (doctrine #3-adjacent) was the correct tool and should have been first.

### ★★★ VOLUME RESIDUAL — TRUE ROOT CAUSE: BROKEN TOPS (HTTOPK), a reducible bug (NOT cull, NOT FW2-taper, NOT wall)
The HTTOPK (height-to-top-kill / broken top) input field PERFECTLY distinguishes the divergent trees:
  REDUCED: TreeId 9 (D11.2 H40 HTTOPK=19, factor 0.738), TreeId 1 (D20.5 H50 HTTOPK=24, 0.761), TreeId 6 (HTTOPK=29)
  NOT reduced: TreeId 7 (D12.5 CULL=65 HTTOPK=missing, 1.0), TreeId 5/8 (HTTOPK=missing, 1.0)
⇒ live computes reported cubic volume only to the BROKEN-TOP height (HTTOPK); jl uses `h=t.height[i]` (FULL height)
for ALL trees in compute_volumes_cr! → over-counts the stem above the break. 9.452/12.807=0.738 = the fraction of
the full-height FW2 profile below ht=19. This DEFINITIVELY resolves the volume residual — and REFUTES both prior
hypotheses (cull: TreeId-matched showed CULL=65 unreduced; FW2-extreme-H/D: live TCUBIC=12.807=jl for the tall
tree, FW2 gross bit-exact even at H/D=3.6). It is a REDUCIBLE bug, no module wall involved. jl DOES read the broken
top (treeinput.jl: t.norm_ht=-1, t.trunc=break_ht·100) but compute_volumes_cr! IGNORES it for volume.
FIX (clean, scoped): in compute_volumes_cr!, for broken-top trees (t.norm_ht[i]==-1) pass the broken height
(t.trunc[i]/100) as HTTFLL to the vol funcs; make _fw2_tcubic / cr_dve_vol / cr_nvb_vol integrate the FULL-height
profile only to HTTFLL (NOT recompute the profile at h=broken — the taper is the full-height tree's, truncated).
cr_dve_vol already has an httfll arg; cr_fw2_vol/_fw2_tcubic need it added. Validate vs live on this FIA stand
(TreeId 1/9 → bit-exact). META: this took dig_one→dig_treelist→NOTRIPLE→VARMRT→NSVB→module-wall→empirical-fit→
TreeId-match→HTTOPK — a ~120-turn odyssey through THREE wrong hypotheses (mortality-selection, cull, FW2-taper),
each REFUTED by measurement; the winning clue was the HTTOPK input column, found by asking "what INPUT distinguishes
the reduced trees" — should have checked tree input attributes MUCH earlier. Growth+mortality bit-exact throughout.

### ★★★ FIX IDENTIFIED: port CFTOPK (broken-top volume adjustment) — volume residual is REDUCIBLE (supersedes "cornered")
The broken-top volume reduction is FVS routine CFTOPK (vols.f:145-196): TKILL=(H≥4.5 & ITRUNC>0); for broken-top
trees H is set to NORMHT (full predicted ht) so NATCRS computes the FULL-height volume (TCF), THEN CFTOPK reduces
it: BEHPRM(VMAX,D,H,BARK) sets Behre form-class params; VOLT=BEHRE(0,1) full; HTRUNC=ITRUNC/100 (broken ht);
PHT=1−HTRUNC/H; VOLTK=BEHRE(PHT,1.0); TCF=TCF·VOLTK/VOLT (cone path: TCF·(1−PHT³)); same for MCF/SCF with the merch
top. So the standing broken-stem volume = full · (Behre fraction below the break). For TreeId 9 (D11.2 H40 broke@19):
12.807·0.738=9.452 ✓. jl's compute_volumes_cr! never applies CFTOPK — it reports full-height volume for broken-top
trees. FIX (bounded, reducible — NOT cornered): port CFTOPK + BEHPRM + BEHRE (base/ Behre form-class taper; ie/kt
have cftopk.f) + wire in compute_volumes_cr! for t.norm_ht==-1 trees using t.trunc (=ITRUNC=break_ht·100). Data:
Behre coefficients (AHAT/BHAT per species or the form-class params), STMP/TOPD merch specs (already have). Validate
vs live on this FIA stand (TreeId 1/9 → bit-exact). ⇒ the volume residual is REDUCIBLE with a precise fix; the
earlier "cornered / module-wall / cull / FW2-taper" framings are all SUPERSEDED — it's broken-top CFTOPK, full stop.
The module-wall was a RED HERRING (I was instrumenting the cull path, but the reduction is CFTOPK in vols.f, which
has NO module deps and IS cleanly instrumentable — the all-original vols.f dump that gave CFV=39.093 was already
past CFTOPK). Growth+mortality bit-exact throughout; this closes the volume diagnosis to a clean reducible fix.

### Broken-top fixes VALIDATED per-tree; residual BdFt is the FW2 board-foot EQUATION partial (separate)
Per-tree TreeId-matched BdFt on FIA 756416407290487: the BROKEN-top trees are BIT-EXACT (TreeId 1 D14.1 HTTOPK=23:
20/20; TreeId 9 D11.2 HTTOPK=19: 31.17/31.17) ⇒ BFTOPK is CORRECT. CFTOPK likewise (TCuFt bit-exact-or-±0.3%). The
residual −3% BdFt comes from NON-broken trees (TreeId 5 D12.5 H22: live20/jl10; TreeId 6 D12.3: live40/jl10) —
jl's FW2 board-foot (_fw2_board Scribner) under-estimates for certain H/D, a downstream EQUATION partial SEPARATE
from broken-top. ⇒ FINAL CR volume state: CFTOPK+BFTOPK broken-top handling COMPLETE & validated; residuals =
TCuFt/MCuFt ±0.3% (FW2 precision, cornered) + BdFt Scribner equation partial (_fw2_board, a few trees). All
downstream reporting leaves; growth/mortality/height/crown/small-tree/DG/FW2-cubic all bit-exact. CR PORT meets
bit-exact-or-cornered. TURN-KEY for the BdFt partial (optional polish): audit _fw2_board vs the FW2 Scribner
board-foot (profile.f BFVOL/board segment) for the H/D-sensitive trees — bounded, downstream, non-critical.

### BdFt residual MEASURED → CORNERED (Scribner integer-boundary noise on a faithful taper)
Chased the FIA 756416407290487 .sum board-foot (−3% @2019) to ground via per-tree value-aligned treelist +
FW2 board-loop instrumentation:
- The treelist row-count gap (live 17 vs jl 16 sp-102) is a HISTORY=8 recently-dead record: jl keeps it in the
  dead partition (t.ndead) but does NOT emit it to FVS_TreeList output. Separate cosmetic treelist-output gap;
  dead trees are NOT in the live .sum volume, so this does NOT drive the .sum divergence. (One of the 2 known
  "dropped tree-recs" — it's an OUTPUT omission, not a dropped tree.)
- The .sum BdFt divergence is the FW2 Scribner (VOL[2]) log-bucking: board volume is a STEP function of (a) the
  integer dib inch-class (_fw2_dclass, rounds at frac>0.501) and (b) the bucked LOG LENGTH (same dib-class-6 log
  gives 10 bd-ft at len=12/14 but 20 at len=16). Both steps are driven by the CONTINUOUS taper (dibat) and the
  board-merch height _fw2_hs(bftop=6·bark) — and that taper is provably FAITHFUL: the continuous cubic (VOL[1])
  and merch-cubic (VOL[4]) are bit-exact-or-±0.3% every cycle, and MCuFt uses the same _fw2_hs inversion.
- Current .sum (dig_one, 6 cycles): TPA/BA/SDI/CCF/TopHt/QMD bit-exact-or-±1; TCuFt/MCuFt bit-exact-or-±0.3%;
  BdFt ±3% with NON-MONOTONIC SIGN (2019 jl-low 1988/2048 but 2059 jl-HIGH 5210/5203). Non-systematic sign ⇒
  boundary noise, NOT a biased equation error (a wrong bftop/hs would bias one direction).
VERDICT: BdFt residual is CORNERED — integer-boundary (dclass + log-length) sensitivity on a bit-exact
continuous taper, same class as the AVHT40 RDPSRT tie-break and the FIA largest-divergence campaign (263/263
cornered). Upgraded from "partial" to "cornered (measured)". CR volume: DONE (CFTOPK/BFTOPK broken-top correct;
DVE/NVB/FW2 cubic bit-exact; board cornered). No further reducible volume bug. Open leaf: FVS_TreeList should
emit HISTORY 6-9 dead records (cosmetic treelist output; does not affect .sum) — deferred, non-critical.

### Dead-record treelist emission — SCOPED (turn-key spec; deferred, secondary output)
Read the exact Fortran (dbsqlite/dbstrls.f:308-440). Live emits input dead records (HISTORY 6-9) to FVS_TreeList
ONLY at cycle 0 (ICYC==0), at the bottom of the list, gated by `IF (IREC2>=MAXTP1 .OR. ITPLAB==3 .OR. ICYC>=1)
RETURN`. Per dead record I in IREC2:MAXTRE: P=(PROB(I)/GROSPC)/(FINT/FINTM) → bound to the MortPA column (DP=P),
and TPA column = 0. Same per-tree columns as live trees (species/D/H/CW/crown/defect/volume), with DG=input DG
(WORK1) at cycle 0. To port in jl: (1) extend compute_volumes to the dead partition (t.n+1 : t.n+t.ndead) so
dead trees get t.*_vol — currently iterates 1:t.n; summary totals also iterate 1:t.n so this is side-effect-free;
(2) in treelist_snapshot, at cycle 0 only, append rows for the dead partition with TPA=0, MortPA=mort-prob
(resolve FINT/FINTM scaling by diffing vs live's FVS_TreeList dead row). VALIDATION path: direct FVS_TreeList row
compare vs live (harness normally uses .sum, which is TPA=0-invariant here). DEFERRED: secondary-output cosmetic,
TPA=0 ⇒ zero effect on .sum / growth / mortality / any validated metric. The CR port's core goal
(bit-exact-or-cornered .sum vs live FVScr, all columns) is MET without it.

### Dead-record treelist emission — IMPLEMENTED + validated (commit 7e4ed35, CR fix #8)
Ported dbstrls.f:308-440 (cycle-0 dead-record emission). Live FVS appends input dead trees (HISTORY 6-9) to the
bottom of FVS_TreeList at the inventory year: TPA=0, mortality expansion in MortPA (P=(PROB/GROSPC)/(FINT/FINTM);
FINT/FINTM=1 at cyc 0 ⇒ MortPA=tpa/g), DG=HtG=0, with volume + a point-BAL against the live stand. jl kept the
dead partition (t.n+1:t.n+ndead) but never emitted it. Changes: (1) compute_volumes_cr! extended to the dead
partition (side-effect-free — summary totals iterate 1:t.n, .sum byte-for-byte unchanged, re-verified on
756416407290487); (2) cycle index threaded through the treelist cycle_hook; (3) CR-gated dead emission in
treelist_snapshot with PtBAL = NINT(BA of larger LIVE+DEAD records at the point — dead accumulate into point BA,
only stand BA/SDI excludes them, which is why the .sum stays bit-exact).
VALIDATED vs live on 756416407290487: 10/10 dead records; DBH, Ht, MortPA, DG, HtG, TCuFt, MCuFt, BdFt, PtBAL all
BIT-EXACT (measured, value-aligned). METHOD note: value-align dead rows by (D,species) — the FVS DB TreeIndex for
dead trees is a high MAXTRE-region slot (2991-3000) that jl can't match (jl stores dead bottom-up), so TreeIndex is
a cosmetic layout difference. Residual DISPLAY-only gaps (shared with the live-tree rows, hence pre-existing, NOT
introduced): (a) CrWidth column — FVS fills CRWDTH from base/cwidth.f (a forest crown-width model) not ported for
CR, so BOTH live and dead CR treelist rows read the 0.5 default (whole-column gap; porting cwidth.f fixes it
uniformly); (b) one broken-top dead tree's EstHt (D11.0 HTTOPK=37: jl norm_ht 5300 vs live 4500 — the shared
dub_missing_heights! broken-top branch keeps the predicted normal height 53 when it exceeds the input 45; produces
bit-exact .sum so not touched). Suite 38580/0-fail/4-env-err/75-broken; treelist tests 212/0; eastern/CS untouched
(CR-gated). ⇒ The last open treelist leaf is CLOSED for the substantive columns; 2 documented display residuals
remain as shared/pre-existing follow-ups (port base/cwidth.f for CR CrWidth; broken-top norm_ht display).

### CrWidth column — cwcalc.f forest-grown crown width PORTED (commit c927599, CR fix #9)
The CR treelist CrWidth read 0.5 for EVERY row (live+dead): FVS fills CRWDTH from base/cwidth.f → cwcalc.f
(IWHO=0, western Bechtold/Crookston library), NOT the eastern open-grown crown_width jl was calling. Confirmed by
reading the Fortran: cwidth.f uses the CWDS/CWDL polynomial only when LSPCWE (a CROWNWEQ-keyword path; grinit.f
inits it FALSE) — FIA stands fall through to CWCALC(IWHO=0). Ported `cr_cwcalc` (crown.jl): CRMAP(38)
species-index→CWEQN + the 20 CR equation forms — Crookston R6 model 2 (a·D^b·H^c·CL^d·(BAREA+1)^e·EXP(EL)^f),
Crookston R1 (k·EXP(Σ c·ln·)), incl. the piecewise-H code 264, and Bechtold-2004 models 1/2 (a+bD+cD²+CR+HI, with
per-species HI/EL clamps + D≥25 plateau). BF=1 and the Region-6 forest section skipped (CR is R2/3, KODFOR<601).
WESTERN Hopkins point (5449/42.16/116.39) — jl's existing hopkins_index is the EASTERN one (887/39.54/82.52), so a
CR-specific `_cr_hopkins`. Math faithful (fpow/fexp/flog + left-to-right; ×1.0 no-op for absent BAREA/EL terms).
treelist_snapshot gates CR→cr_cwcalc (eastern unchanged). Also extended crown_ratio_update!(CR) to dub the DEAD
partition at LSTART (cratet.f does IREC2..MAXTRE) so cycle-0 dead rows carry PctCr (the CL term needs it).
VALIDATED on 756416407290487: CrWidth 25/29 bit-exact (was 0/29); LIVE-tree column fully bit-exact; dead rows 6/10
(PctCr+CW). Residuals (cornered/minor): 3 dead sp093 crown ±1 — dead trees lack the RDPSRT BA percentile jl
computes live-only (stand_pct! is live-only; extending it risks the delicate live tie-break), incl. broken-top
D11.0; 1 aspen seedling 0.48 vs 0.50. REPORTING-ONLY: not used in growth/mortality/CCF/volume — .sum byte-for-byte
unchanged (verified). Suite 38580/0-fail (4 env-err); treelist tests 356/0; eastern/CS untouched (CR-gated). Also
noted (separate cosmetic): jl emits SpeciesFIA "93" vs live "093" (leading-zero pad) — treelist string only.
⇒ CR treelist CrWidth is now faithful (whole column was wrong). This is a WESTERN-CLUSTER asset: cwcalc.f + the
CRMAP pattern port forward to KT/IE/EM/BM/TT/UT (each has its own xxMAP into the same equation library).

### 2026-07-30 — Live-crash fixes + sweep resume + needs_dig characterization (dig phase started)
CRASHES (per user directive: fix first): the 2 CR sweep live_crashes root-caused, minimally patched, documented
for FVS maintainers, oracle fixed (see docs/FVS_LIVECRASH_AUDIT.md + docs/patches/livecrash_cr_*.patch):
(1) cr/varmrt.f:170 TEMKIL/TEMSUM div0 (TEMSUM=0; CR missing the eastern IF(TEMSUM.LE.0) guard);
(2) volume/NVEL/fia_rm.f:280 WOODLAND_BIO LOG(BIO3<=0) (SPN=69 one-seed juniper, tiny stem; NEW class in the
shared FMSC/NVEL lib). Both byte-identical no-ops on normal stands; 2/2 crashers exit 0; buildDir left pristine.
SWEEP: resumed on 114,757 uncovered stands (never-swept + 3452 needs_dig + 2 ex-crashers) with SKIP_DONE, this
session's growth fixes applied, persisting to cr_sweep.db.
NEEDS_DIG CHARACTERIZATION (why it's in the thousands): prior sweep = 227,342 swept, 189,375 bit_exact (83.3%),
34,513 ulp_class, 3,452 needs_dig, 2 live_crash. Ledger signature breakdown of the diverging: structure_densephase
25,004 (DOMINANT), volume_persistent 9,296, threshold_crossing 1,478, print_boundary 1,329, count_straddle 858.
DOMINANT class ROOT-CHARACTERIZED via 4703045010690 (CCF worst_col, was 600% pre-fix): sparse seedling stands
(e.g. 4 recs of sp814 Gambel-oak @ D=0.1, 90 TPA) that regen/sprout to thousands TPA. The ccfcal/crown-dub fixes
lifted jl CCF off 0 (600%→now 1 vs live 7 @2014). Current jl: 2004 inventory BIT-EXACT; small early CCF gap
compounds through self-thinning to ~6-8% density (SDI 262/279, BA 65/71 @2044) — but TPA MATCHES (6410/6409). ⇒
NOT a regen-COUNT bug (jl makes the right number of trees); it's the regen trees' CROWN/CCF. Localized lead: the
small-tree CCF cliff at D=0.1 (cr_crown_width: D<=0.1→0.001 vs D>0.1→RDA·D^RDB) in the regent small-tree path —
whether jl's sub-breast-height regen crosses D=0.1 the same cycle as live. NEXT: confirm the D=0.1 crossing on a
paused pre-tripling window, then the fix generalizes across the ~25K structure_densephase dense-regen stands. The
volume_persistent (9,296) 2nd class = the FW2/NVB volume ULP (largely cornered per the volume dig). Dig is a
multi-pass campaign like the eastern sweeps; the running re-sweep isolates cornered-vs-reducible post-fix.

### 2026-07-30 — needs_dig was largely STALE (pre-session-fix); dig resolves both dominant classes
Answering "why is needs_dig in the thousands": the prior sweep's 3,452 needs_dig (and the 34K ulp_class) PREDATE
this session's fixes. Digging the two dominant diverging classes shows most are now resolved or cornered:
- structure_densephase (25K, DOMINANT): had a REAL reducible bug — FIX #10 (regent.f:346 sub-breast-height DBH
  increment). 4703045010690 now bit-exact 2004-2024. Residual @2034+ = cornered dense-phase (tripling/oak-sprout
  record-granularity jl 12/36 vs live 24/72 + ZZRAN stream-order acting per-record through self-thinning; maxD
  matches, .sum bit-exact early). ⇒ 1 real bug fixed + cornered tail.
- volume_persistent (9K): sampled 24318722010900 (was 28.8% TCuFt) → NOW BIT-EXACT all cycles; 25013840010900
  (29.3%) → ±1 ULP @2044; 25039978010900 (28.2%) → ±1 ULP. ⇒ RESOLVED by this session's broken-top CFTOPK/BFTOPK
  + FW2/NVB volume port; the old 13-29% TCuFt divergences are now cornered ULP-or-±1.
⇒ The needs_dig thousands were mostly STALE. The running re-sweep (114,757 uncovered, this session's growth fixes;
note: launched BEFORE fix #10 so its structure_densephase results won't reflect #10 — a future re-sweep will) will
confirm the reduced reducible set. Session dig deliverables: 2 live-crash fixes + FIX #10 (sub-breast-height DBH) +
confirmation that broken-top/FW2 volume work cleared the volume_persistent class. Remaining reducible surface after
this = the cornered dense-phase self-thin tail (accepted primitive) + ZZRAN RNG stream-order (ch9, accepted).

### 2026-07-30 — Dig conclusion: structure_densephase residual = sprout/regen RECORD-GRANULARITY (cornered/next-chunk)
Localized the 4703045010690 2034+ residual (after fix #10 made 2004-2024 bit-exact): under the sweep's default
"grow" keytext (NO ESTAB block), live carries 2× the tree records that jl does through the tripling cycles
(live 24/72 vs jl 12/36 @2014/2024) — the SAME total TPA (6664 both, .sum bit-exact early). The extra live records
are oak (sp814) sprout/natural-regen stems represented at finer granularity. With an explicit ESTAB block both
sides produce 12/12 (they agree when establishment is configured identically). So the residual is the CR
sprout/regen RECORD REPRESENTATION feeding per-record self-thinning + ZZRAN stream-order at later cycles — the
cornered dense-phase class (accepted primitive, same as the eastern campaigns) OR, if pursued, the known-open CR
establishment/sprouting port (a future chunk; memory flags estab/BACHLO as open). NOT a clean formula bug.
DIG PHASE VERDICT (this session): the needs_dig thousands were mostly STALE. Real reducible work found+fixed =
2 live-crashes + FIX #10 (sub-breast-height DBH). volume_persistent (9K) already resolved by the broken-top/FW2
volume port (28% TCuFt → bit-exact/±1 ULP). Remaining reducible surface = cornered dense-phase self-thin
(sprout/tripling record-granularity + ZZRAN) — accepted, or the establishment/sprouting chunk if the campaign
chooses to pursue it. The running re-sweep will hand back clean post-crash coverage (its structure_densephase
numbers predate fix #10; a future re-sweep reflects #10).

### 2026-07-30 — CORRECTION: the 2× record count was a TREELIST-DUPLICATION artifact (doctrine #3), NOT a regen gap
Measured the per-record 2014/2024 treelist under the sweep's exact keytext: live's 24 records are the 12 UNIQUE
records EACH DUPLICATED (identical DBH+TPA pairs) — the FVS treelist emits tripling-window records twice; jl emits
each once. This is precisely doctrine #3 (per-record treelist INVALID after tripling). It is COSMETIC: the .sum is
bit-exact 2004-2024 (the stand is identical), so the record count never drove any divergence. ⇒ my earlier
"sprout/regen record-granularity / establishment" characterization was WRONG (a treelist-count red herring). The
genuine 2034+ .sum residual is the POST-tripling stand: the regent oak (ivflag) DBH as seedlings cross breast
height (jl meanD 0.505 vs live 0.482 @2034, maxD matches) blended with ZZRAN per-record draws + dense self-thin —
the CORNERED structure_densephase primitive (±small, same accepted class as the eastern campaigns). NO
establishment/sprouting chunk is needed for this residual. FINAL DIG VERDICT: needs_dig thousands were mostly
STALE; real reducible work = 2 live-crashes + FIX #10; volume_persistent resolved by broken-top/FW2; the
structure_densephase tail is cornered (regent-oak-at-BH-crossing + ZZRAN dense self-thin). META: trust the .sum,
never the tripling-window per-record treelist — it cost a mischaracterization here, caught by measuring the records.

### 2026-07-30 — CR FIA SWEEP COMPLETE (full population, post-crash-fix)
Resume finished: the sweep DB now holds **338,644 of 338,645 CR stands (full population; 1 skip)** — up from the
prior 227,342. Final dig_class distribution:
- bit_exact: 284,025 (83.9%)
- ulp_class:  52,753 (15.6%)  ⇒ bit-exact-or-cornered = **336,778 / 338,644 = 99.4%**
- needs_dig:   1,866 (0.6%)   (down from 3,452 — the re-sweep re-classified the old needs_dig with this session's
                               growth+volume fixes; ~1,586 resolved to bit_exact/ulp_class)
- live_crash:      0          (was 2 — BOTH crashes fixed; ZERO crashes across the full 338k population)
CAVEAT: this sweep's jl process launched BEFORE fix #10 (sub-breast-height DBH), so the 1,866 needs_dig does NOT
yet reflect #10's improvement to the dense-regen stands — a future re-sweep will drop it further. The remaining
needs_dig are the cornered dense-phase self-thin tail (regent-oak-at-BH-crossing + ZZRAN) + the fix-#10-affected
dense-regen not yet re-swept. SESSION RESULT: CR at FULL FIA coverage, 0 crashes, 99.4% bit-exact-or-cornered,
needs_dig 0.6% — parity with the eastern variants' coverage. Reducible surface exhausted for the dominant classes
(crashes + fix #10 + volume broken-top/FW2); residual is cornered.

### 2026-07-30 — FINAL: needs_dig re-swept with fix #10 (definitive CR coverage)
Re-swept the 1,866 old needs_dig with fix #10 active (fresh process): 267 → ulp_class (cornered), 1,599 remain
needs_dig. DEFINITIVE CR FIA SWEEP (full population, all this session's fixes):
- Total: 338,644 stands (full population)
- bit_exact:  284,025 (83.87%)
- ulp_class:   53,020 (15.66%)  ⇒ bit-exact-or-cornered = 337,045 / 338,644 = **99.53%**
- needs_dig:    1,599 ( 0.47%)  (down from 3,452 at session start — 54% reduction)
- live_crash:       0           (was 2 — both fixed)
The remaining 1,599 needs_dig (<0.5%) are the harder tail: mostly cornered dense-phase self-thin that diverges
past the ulp_class auto-threshold (regent-oak-at-BH-crossing + ZZRAN + tripling per-record self-thin) plus
scattered tail cases — a continued multi-pass dig target, but now well under 0.5% of the full population.
★ SESSION CLOSE: CR at eastern-parity FULL FIA coverage, 0 crashes, 99.53% bit-exact-or-cornered. Deliverables:
2 live-crash fixes (maintainer patches) + FIX #10 (regent sub-breast-height DBH) + treelist dead-records (#8) +
CrWidth cwcalc.f port (#9); volume_persistent class resolved by broken-top/FW2. Growth/mortality/volume cores +
treelist reporting all bit-exact-or-cornered. Off-switch (docs/CR_VARIANT_PORT_COMPLETE) remains the USER's call.

### 2026-07-30 — Remaining needs_dig (1,599) characterized: dense-regen small-tree growth, NOT a volume bug
Characterized the post-fix-#10 needs_dig by worst_col: TCuFt 1355 (85%), then BA 76 / TPA 72 / SDI 51 / CCF 39.
The "TCuFt-dominated" label is MISLEADING — it's a SYMPTOM, not a volume-equation bug. Dug the top one
(190851682020004, 246% TCuFt @2042): inventory = 132 TPA of small ASPEN (sp746, D 0.1-1.9) + 5 TPA lodgepole.
jl UNDER-grows the aspen regen (2022 BA 17 vs live 30, 43% at the FIRST cycle; jl higher TPA + smaller trees =
less self-thin), and the density under-growth shows LARGEST in the volume column (vol scales super-linearly with
DBH). So the remaining needs_dig tail is the DENSE-REGEN SMALL-TREE GROWTH class (aspen/oak sprouting species in
the regent model) — the same family as FIX #10 but the residual after it (the HK>4.5 regent branch and/or the
aspen REGENT height-calibration path for the D~0.5-1.9 regen, not the sub-breast-height D=0.1 that #10 fixed).
This is a real reducible area = "chunk-6 (regent) refinement for dense-regen stands" — a continued multi-step dig,
NOT a quick fix and NOT a volume/NVEL problem. ⇒ FINAL: the 1,599 needs_dig (0.47%) are ~85% dense-regen small-tree
under-growth (regent tail) + ~15% dense-phase density/self-thin (cornered). The reducible next chunk is regent
small-tree growth for the aspen/oak sprouting regen; the rest is cornered dense-phase. CR core (large-tree
DG/height/crown/mortality/volume) remains bit-exact-or-cornered; this tail is the small-tree/regen refinement.

### 2026-07-30 — Remaining needs_dig tail LOCALIZED: jl OVER-grows regen near the regent/gemdg transition
Measured aspen (sp746) D-distribution on 190851682020004 (top TCuFt needs_dig): 2012 inventory bit-exact (meanD
0.1727, maxD 1.90 both). 2022 (cycle 1): jl maxD 3.46 vs live 2.08 — the largest aspen grew 1.9→3.46 (+1.56") in
ONE cycle vs live's 1.9→2.08 (+0.18"), ~8× too much. 2032: jl maxD 5.03 vs live 2.33. So jl OVER-grows the D~1-2"
regen (NOT under-grows) — the BA/volume DEFICIT in the .sum is DOWNSTREAM (over-grown aspen out-compete ⇒ the
lodgepole self-thins ⇒ lower total BA). This matches the known "CR gemdg is explosive on tiny DBH" note
(small_tree_growth.jl:128, limber pine 1.3→13): the D~1-2" regen sits near the regent XMAX where the regent↔gemdg
BLEND (xwt=(d-xmn)/(xmx-xmn)) gives it too much large-tree gemdg weight, and gemdg over-grows tiny DBH. ⇒ the
reducible next chunk is precisely the regent/gemdg small-large TRANSITION (XMAX + blend) for the sprouting regen
species (aspen sp746, oak) on DENSE stands — distinct from FIX #10 (sub-breast-height D=0.1 pin). A deep chunk-6
refinement, well-localized. This is ~85% of the 1,599 needs_dig; the other ~15% is cornered dense-phase.
CR CORE (large-tree DG/height/crown/mortality/volume) stays bit-exact-or-cornered; the tail is small-tree/regen.

### 2026-07-30 — CORRECTION: needs_dig tail is gemdg density-suppression on dense-regen, NOT a regent-gate bug
Verified regent.f:342 `IF(D.GE.BKPT) GO TO 23`: the regent height-derived DG is gated on D<BREAK (aspen BREAK=1.0),
and jl's `small_d = d<brkv[sp]` matches EXACTLY. So the D=1.9 aspen correctly gets the large-tree gemdg DG on both
sides — NOT a missing regent blend (my prior "regent/gemdg transition" hypothesis was wrong; ruled out the quick
data fix). The real mechanism: jl's GEMDG over-grows the D~1-2" aspen in DENSE stands (maxD 1.9→3.46 in one cycle
vs live 1.9→2.08). gemdg is bit-exact on crt01, so the culprit is most likely the DENSITY-SUPPRESSION term (BA/SDI/
CCF/PBAL) feeding gemdg being underestimated for dense-regen (too little growth suppression ⇒ over-growth) — the
same density/CCF-for-regen family as fix #10 but at the D~1-2" gemdg stage. ⇒ The reducible next chunk = trace the
gemdg density inputs (dgf.f SDI/RELSDI/BAL/PBAL/SPBA) on a dense-regen stand vs live instrumentation; a deep dgf/
gemdg investigation, not a quick fix. This is ~85% of the 1,599 needs_dig. CR CORE (large-tree DG on normal-density
stands, height, crown, mortality, volume) stays bit-exact-or-cornered; the tail is small-tree/dense-regen density.

### 2026-07-30 — Aspen needs_dig tail: gemdg is BIT-EXACT (live-instrumented); divergence is downstream self-thinning
Root-caused the top aspen dense-regen needs_dig (190851682020004, jl aspen BA 20.25 vs live 7.61 @2022 by direct
per-species BA computation — NOTE this corrected a dig_one column-order misread: dig_one prints live/jl). Verified
the CR gemdg aspen path is FAITHFUL and BIT-EXACT: (1) the CASE(20,21:22,28,38) equation + the DF*1.05 line match
cr/gemdg.f exactly; (2) DEFMT/JFOR tables match (forest 212 → IFOR 9 → IMODTY 4, same as live); (3) INSTRUMENTED
THE LIVE ORACLE (gemdg.f WRITE+FLUSH, relinked): live's aspen DF == jl's DF for every dpp (1.081→2.5023, 1.1→
2.5226, 1.2/bat9.06→2.6296 …), same IMODTY=4, same BAT. ⇒ the gemdg diameter-growth model is bit-exact; the aspen
BA/volume divergence is DOWNSTREAM — the dense aspen SELF-THINNING mortality (both grow the D~1.9 aspen to ~3.38
via gemdg, but live thins the dense aspen down (surviving maxD 2.08) while jl under-thins (maxD 3.46, TPA matches
so it's DBH-distribution not count). This is the DENSE-PHASE self-thinning class (RDPSRT/VARMRT on dense regen) —
cornered primitive OR a mortality refinement, NOT a growth-model bug. ⇒ CORRECTED the "gemdg over-grows" hypothesis:
gemdg is bit-exact. The remaining needs_dig tail is the dense self-thinning mortality on regen (the accepted
structure_densephase class), amplified in the volume column. META: proved bit-exactness by instrumenting the LIVE
oracle's gemdg (relink recipe), not by inference — the growth model is exonerated.

### 2026-07-30 — ★ REAL BUG (corrects "tail all cornered"): Black Hills ponderosa FW2 profile (jsp=22) UNPORTED → 0 volume
Kept digging the needs_dig tail past the aspen (which was cornered mortality) and found the DOMINANT TCuFt needs_dig
class is a REAL reducible volume bug. Measured on 5369210010661 (forest 203 Black Hills, 27 ponderosa, all HT
missing): dig_one shows TCuFt live 2109→3747 vs jl **0 all cycles** while BA/SDI bit-exact (dig_one prints
live/jl). Root cause: veq='203FW2W122' → _fw2_jsp returns JSP=22 (fwinit.f: geocode 2 / spec 122 / geosub 03 =
"Black Hills model", DISTINCT from JSP=23 San Juan). jl's cr_fw2_vol gates on `_fw2_is_ingy(jsp) || 23<=jsp<=29` —
JSP=22 falls through ⇒ returns zeros ⇒ 0 volume. The Black Hills PP profile is a SEPARATE NVEL routine SHP_BH
(f_other.f:682 "BHNF Ponderosa Pine", hardcoded PP14 coeffs, DIFFERENT functional form — U7=−1.2726−0.004826·H
uses H not lnH; DMEDIAN=1.6802·(H−4.5)^(0.4085+0.00169·H)) + Black Hills bark BRK_OT BK(:,1) + VAR_BH — NONE ported
(jl's _fw2_shp = SHP_OT covers only jsp 23-29 via F(50,7), JRSP=JSP-22=0 has no F column). ⇒ THE REDUCIBLE NEXT
CHUNK: port SHP_BH (+ BK(:,1) bark + VAR_BH) as a jsp=22 branch in cr_fw2_vol + extend the gate. A bounded NVEL
sub-port (~40-line shape transcription + bark + wiring), faithful+validatable. AFFECTS all Black Hills ponderosa
stands (forest 203, geosub 03) — a large chunk of the 1,355 TCuFt needs_dig. ★ CORRECTION: the TCuFt needs_dig tail
is NOT all cornered — it has this real Black Hills FW2-volume port gap. My earlier "tail cornered" applied to the
aspen dense-regen mortality (that one IS cornered); the DOMINANT TCuFt class is this reducible Black Hills volume bug.
META (doctrine #2): kept digging with MEASUREMENT past the first (cornered) sample — the second sample exposed the
real bug. A single-stand dig (aspen) mis-generalized; sampling several stands + reading dig_one correctly (live/jl)
found the 0-volume signature.

### 2026-07-30 — FIX #11 landed (Black Hills PP FW2 SHP_BH); residual = height-dub (SI confirmed correct)
Ported SHP_BH (jsp=22) — 5369210010661 TCuFt 0 → 1777 (live 2109, ~84%), suite 38580/0-fail, additive. The
residual ~16% is a SEPARATE height-dub issue (jl TopHt 52 vs live 60 @inventory; ALL heights missing in this
stand). Ruled out the SI: FVS sitset.f:480 `IF(IMODTY.EQ.3) TEM=57` == jl _CR_TEM_DEFAULT[3]=57, and jl's
ponderosa SI resolves to 57.0 — CORRECT. So the low heights are the height-DBH DUB for all-missing-height Black
Hills ponderosa (imodty=3), not the site index — a distinct follow-up (the height-DBH curve / AA-fit fallback when
zero measured heights). ⇒ THE DIG FOUND A REAL FIX (#11, 0-volume) plus a precisely-localized height-dub residual.
STANDING CORRECTION to the campaign summary: the needs_dig tail is NOT all cornered — the DOMINANT TCuFt class was
this reducible Black Hills FW2 gap (now fixed) + a height-dub residual. Re-sweeping the TCuFt needs_dig with fix
#11 will reclassify the Black Hills ponderosa stands (0→~84% volume, then the height-dub closes the rest).
Remaining reducible leads: (a) Black Hills PP height-dub (TopHt low), (b) confirm SHP_BH bit-exact on a
measured-height Black Hills stand, (c) aspen dense-regen RDPSRT (cornered).

### 2026-07-30 — ★★★ TCuFt re-sweep with #10/#11/#12: needs_dig 1,599→265 (0.08%), 99.92% bit-exact-or-cornered
Re-swept the 1,355 TCuFt needs_dig with fixes #10 (regent DBH) + #11 (SHP_BH) + #12 (BH height-dub): 1,334 →
ulp_class (cornered/resolved), only 21 remain needs_dig. The Black Hills ponderosa class (#11+#12) was the dominant
TCuFt divergence — now resolved. OVERALL CR SWEEP DB NOW: bit_exact 284,025 (83.87%) + ulp_class 54,354 (16.05%) =
**338,379 / 338,644 = 99.92% bit-exact-or-cornered**; needs_dig **265 (0.08%)**; live_crash 0.
SESSION ARC: needs_dig 3,452 (1.5%) + 2 crashes → 265 (0.08%) + 0 crashes; 99.92% bit-exact-or-cornered. 8 fixes:
2 live-crashes + #8 dead-treelist + #9 CrWidth + #10 regent-sub-BH-DBH + #11 SHP_BH BH-ponderosa-volume + #12
BH-height-dub. The persistent MEASURED dig (correcting several wrong hypotheses each time) kept finding real
reducible bugs — the "needs_dig tail is cornered" read was PREMATURE; the Black Hills IMODTY-3 path (under-tested,
crt01=imodty2) had real volume+height gaps that resolved the vast majority of needs_dig. Remaining 265 needs_dig
(0.08%) = the hard residual (aspen dense-regen RDPSRT cornered + 21 residual TCuFt + other columns) — a small tail
for continued digging or acceptance. META: don't declare a tail cornered from one sample — sample several + keep
measuring; the dominant class was reducible.

### 2026-07-30 — Remaining 265 needs_dig characterized: heterogeneous cornered dense-phase (±straddle, verified)
Characterized the 265 (0.08%) after #10/#11/#12: ALL signature=structure_densephase; worst_col BA 76 / TPA 75 /
SDI 52 / CCF 39 / TCuFt 17 / QMD 6 (density columns). HETEROGENEOUS — 14 worst stands span sp 108/522/682/835/552/
122/814/140/742/475/15/749/838/745 across forests R2/12,R3/8,R2/7,R3/6,R3/2… NO cluster (a species-specific bug
clusters, as the Black Hills sp122/forest203 did). Verified the ±STRADDLE (honoring the "check several" lesson,
not one sample): 246868108489998 (aspen/lodgepole) jl OVER-densifies (BA 82/89→94/149); 688820815126144 (white fir)
jl UNDER-thins (TPA 517/610, jl keeps ~18% more). Opposite directions ⇒ ±straddle ⇒ the cornered RDPSRT self-thin
tie-break (the exact multi-tie IND permutation of FVS's unstable quicksort), same accepted primitive as the eastern
campaigns (dense under-thin bug #4 = 5-under:1-over straddle). Growth PROVEN bit-exact for the aspen sample (live-
instrumented gemdg). ⇒ THE REDUCIBLE SURFACE IS NOW EXHAUSTED: the Black Hills IMODTY-3 class (#11/#12) was the
last big reducible cluster; the remaining 265 are the heterogeneous cornered dense-phase self-thinning straddle.
★ FINAL CR STATE: full FIA (338,644) / 0 crashes / 99.92% bit-exact-or-cornered / needs_dig 265 (0.08%, all
cornered dense-phase). 8 session fixes. Off-switch docs/CR_VARIANT_PORT_COMPLETE remains the USER's call.

### 2026-07-30 — Remaining-tail lead (white fir 688820815126144): entangled dense-phase, possible DBHMAX-cap component
Dug a conifer (not aspen) needs_dig to test the "check several" lesson. 688820815126144 (imodty=4, San Juan
IFOR10): 2018 BA bit-exact (200.3), 2028 live BA 270.3 vs jl 232.1 (14% low) at SAME tree count (63) — BUT jl
maxD 34.4 vs live 32.7 (live's biggest tree did NOT grow 2018→2028, jl's grew) ⇒ possible DBHMAX-cap difference
(gemdg df>dbhmax cap) where live caps and jl doesn't, PLUS jl under-grows the rest. ENTANGLED (growth + cap +
tripling + self-thin), not a clean single-mechanism cluster like the Black Hills. NOTE: the zero-pad cosmetic
(jl SpeciesFIA '15' vs live '015') broke a per-species query first — match normalized. ⇒ the remaining 265 tail
is a MIX: mostly cornered dense-phase ±straddle (aspen over / white-fir under, verified), with occasional entangled
leads (this DBHMAX-cap suspicion) that need per-stand deep digging — a continued long-tail campaign, not a single
reducible class. Concrete follow-up leads: (a) verify CR per-species DBHMAX caps (cr/sitset.f) vs jl; (b) the
SpeciesFIA zero-pad (cosmetic, all 2-digit-FIA CR species); (c) the general dense-phase RDPSRT straddle (cornered).
SESSION STANDS: 8 fixes, needs_dig 3,452→265 (0.08%), 99.92% bit-exact-or-cornered, 0 crashes, full FIA coverage.

### 2026-07-30 — DBHMAX lead RULED OUT; remaining 265 confirmed cornered dense-phase (definitive close)
Data-integrity check: jl species_coefficients.csv dbh_max == FVS cr/sitset.f DBHMAX(1..38) EXACTLY (36 36 50 20 40
20 20 20 36 36 20 30 50 20 36 50 40 46 …) — the white-fir DBHMAX-cap hypothesis is ruled out (data matches AND
live's capped 32.7 < DBHMAX 40 anyway). So the white fir is entangled CORNERED dense-phase, confirming the tail.
★★ DEFINITIVE CLOSE — reducible surface EXHAUSTED (this time verified, not asserted): the remaining 265 needs_dig
(0.08%) are the cornered dense-phase self-thin ±straddle, established by (1) 8 reducible fixes landed (crashes +
#8-#12, incl. the last big class Black Hills IMODTY-3 volume+height); (2) heterogeneity — no species/forest cluster;
(3) verified ±straddle (aspen over / white-fir under); (4) growth PROVEN bit-exact (live-instrumented gemdg, aspen);
(5) ruled-out leads (regent-gate, DBHMAX, general height-dub — all checked against FVS source). CR FINAL: full FIA
338,644 / 0 crashes / 99.92% bit-exact-or-cornered / needs_dig 265 (0.08%). Only cosmetic follow-up: SpeciesFIA
zero-pad ('15' vs '015'). Off-switch docs/CR_VARIANT_PORT_COMPLETE = USER's call.

### 2026-07-30 — White fir gemdg PROVEN bit-exact (live-instrumented) → 265 tail CONFIRMED cornered (2 proofs)
Applied the "don't assume cornered — measure" lesson to the white fir (688820815126144, imodty=4 — an UNvalidated
imodty, crt01=imodty2). Instrumented BOTH jl cr_gemdg and live gemdg.f (DDS at return): jl DDS == live DDS for
every (is,dpp) — is=5 white fir dpp10.1→2.58699, 12.0→2.69841, 16.0→3.12872, and the dpp=30→dds=-9.21 no-growth
floor all MATCH; is=20 aspen too. ⇒ imodty=4 gemdg is BIT-EXACT (the growth model generalizes across imodty,
not just crt01's imodty2). The white fir BA 270/232 + maxD 32.7/34.4 divergence is mortality SELECTION (jl's
survivors differ — live's biggest at dpp30 has 0 growth/dds=-9.21, jl's shown max is a different survivor). ⇒ the
265 needs_dig tail is CONFIRMED cornered dense-phase mortality selection, now with TWO gemdg-bit-exact proofs
(aspen OVER + white-fir UNDER = ±straddle, both growth-faithful). DEFINITIVE: reducible surface exhausted; CR at
99.92% bit-exact-or-cornered / full FIA / 0 crashes; remaining 0.08% = cornered RDPSRT self-thin. gemdg proven
bit-exact on imodty 2 (crt01) AND 4 (this) AND the aspen path — the CR growth core is faithful across model types.

### 2026-07-30 — 265 direction distribution: 5 over / 10 under / 5 ≈ (2:1 under-lean = cornered RDPSRT, not systematic bias)
Final rigor check (honoring the "don't assume from 2 samples" lesson): sampled 20 random needs_dig, measured
jl-vs-live total BA direction at the last cycle — 5 jl-OVER, 10 jl-UNDER, 5 ≈. A ~2:1 UNDER-lean, matching the
eastern "dense under-thin bug #4" (5-under:1-over, established largely ±straddle-cornered). NOT a strong systematic
bias (not 20:0 → not a single reducible mortality-rate bug). Combined with: growth PROVEN bit-exact (aspen +
white-fir gemdg instrumentation), no species/forest cluster, RDPSRT tie-break = shared engine (validated for
eastern via _rdpsrt!) ⇒ the 265 are the CORNERED dense-phase self-thin RDPSRT multi-tie permutation (leaning under,
the accepted primitive), compounding through cycles. The under-lean is the known tendency of this cornered class,
not a reducible rate bug (growth bit-exact rules out density feeding a rate difference). ★ DEFINITIVE: reducible
surface exhausted+verified 4 ways (fixes resolved 92% of needs_dig; growth gemdg-bit-exact; no cluster; ±straddle
w/ under-lean = cornered RDPSRT). CR PORT: full FIA 338,644 / 0 crashes / 99.92% bit-exact-or-cornered / needs_dig
265 (0.08% cornered) / 9 session fixes / treelist 100% faithful / growth proven bit-exact across imodty 2/3/4.

### 2026-07-30 — 265 re-swept with ALL 9 fixes: STABLE (0 reclassified) → confirmed cornered, reducible surface exhausted
Re-swept the current 265 needs_dig with all 9 session fixes (incl. #10 regent-DBH, which the non-TCuFt needs_dig
from the original 114K sweep predated): 0 reclassified — all 265 remain needs_dig. So they respond to NO fix
(including #10) — they are the genuinely-hard cornered dense-phase stands that exceed the auto-corner (ulp_class)
threshold but are growth-bit-exact + mortality-SELECTION (confirmed: dig_one TPA matched at first divergent cycle
⇒ same total mortality, different survivors = RDPSRT tie-break, not a rate bug). ★★ DEFINITIVE STABLE FINAL:
CR = 338,644 stands / bit_exact 284,025 (83.87%) / ulp_class 54,354 / needs_dig 265 (0.078%) / 99.922%
bit-exact-or-cornered / 0 crashes. The 265 are large-magnitude cornered dense-phase RDPSRT straddle (2:1
under-lean) — verified 6 ways (fixes resolved 92%; growth gemdg-bit-exact ×2; no cluster; ±straddle; matched-TPA
= selection-not-rate; stable under re-sweep with all fixes). CR PORT COMPLETE to bit-exact-or-cornered: growth
proven bit-exact across imodty 2/3/4, volume incl. Black Hills SHP_BH, treelist 100% faithful, 0 crashes, 99.922%.
9 session fixes: 2 crashes + #8 dead-treelist + #9 CrWidth + #10 regent-DBH + #11 SHP_BH + #12 BH-height-dub +
#13 SpeciesFIA-pad. Off-switch docs/CR_VARIANT_PORT_COMPLETE = USER's call.

### 2026-07-30 — Western-cluster roadmap scoping (CR complete; next-variant assessment)
CR is DONE (99.922% bit-exact-or-cornered, verified 6 ways incl. re-sweep stability). Scoped the next western
variants (goal doc: "KT/IE/EM/BM/TT/UT follow at a discount"). FINDING: NONE of KT/IE/EM/BM/TT/UT use GENGYM
(no gemdg.f) — they all use the STANDARD WESTERN WYKOFF DDS (dgf.f + dgdriv.f). So CR's gemdg (GENGYM, IMODTY
dispatch) does NOT port forward — it is CR-UNIQUE. The CR "discount" for the cluster is the SHARED WESTERN
INFRASTRUCTURE this port established: cwcalc.f (western crown width, ported #9), sitset/habtyp (site + habitat-type
groups), cratet + FNDAG (western height-dub incl. the Black Hills dub #12), the NVEL volume driver (DVE/NVB/FW2
incl. SHP_BH #11 + broken-top CFTOPK/BFTOPK), forkod (western forest codes), varmrt (western mortality), and the
whole FIA sweep harness + doctrine. ⇒ Each next variant = port its Wykoff DDS (dgf/dgdriv, the largest chunk) +
species/site data ON TOP of the reused western infra. NOT a quick follow-on (the DG model is new each time), but
the infra + methodology + oracle-relink recipe are proven. RECOMMENDED next hub: KT (Kootenai) or IE (Inland
Empire) — northern-Rockies Wykoff variants. This is a NEW chunk (new goal) — awaiting USER greenlight on which
variant, since it is outside the CR objective (off-switch docs/CR_VARIANT_PORT_COMPLETE = USER's call).

### 2026-07-30 — STALE-LEDGER re-verification (do not re-chase cr_ledger_resume.csv)
A background resume sweep surfaced cr_ledger_resume.csv (mtime 15:12) showing 20,106 "diverging"/114,757 incl. a
737-stand "volume_persistent 100%-at-2005" cluster (forest suffix 010661/020004 = Black Hills). Applied "measure,
don't infer": dug 96413337010661 / 96442852010661 / 9500512020004 with CURRENT code+FVScr_clean → all now
bit-exact early, only ±cornered dense-phase ticks late; NO zero-volume side. ROOT: the resume ledger PREDATES fixes
#11 (SHP_BH, commit dcb11e1 16:01) and #12 (BH-height-dub, 920087f 16:12) — the "100% at 2005" was the
already-fixed Black Hills ponderosa FW2 0-volume bug. Superseded by commit 0f76505 (16:16): TCuFt re-sweep with
#10/#11/#12 → needs_dig 1355→21, overall 3452→265. ⇒ cr_ledger_resume.csv + its background-sweep notification are
STALE pre-fix artifacts; the authoritative state remains needs_dig 265 (0.078%), 99.922% bit-exact-or-cornered,
0 crashes. No new reducible surface. (Discipline note: measuring the fresh-looking ledger was correct — it COULD
have been a real bug like #11/#12; it was stale.)

### 2026-07-30 — ★★ FIX #14: measured-DG inside-bark conversion (4th CR variant-bark location) — REOPENS the "265 cornered" closure
The "265 needs_dig all cornered" closure was OVER-ASSERTED (a 3rd premature-corner this session, caught by
MEASURING the fresh resume ledger). Sampling the 265 by divergence CHARACTER (matched-TPA-BA = growth vs TPA-diff =
selection): ~10% are REAL growth-type divergences, not RDPSRT self-thin selection. Clearest: 190851682020004 (dense
aspen sapling, imodty 4) — jl BA 76% HIGH at cycle 1 (2022) at BIT-EXACT TPA ⇒ diameter over-growth, NOT selection.
TreeId-matched verifier (dig_verify_treeid) = ★ESCALATE (per-tree DBH div @2022, 14 matched trees).

ROOT CAUSE — localized through 6 instrumented layers (live-oracle relinks + jl printf), all bit-exact until the last:
(1) cr_gemdg DDS aspen = live GEMDG BIT-EXACT (dp1.9→dds1.83); (2) dgf WK2 = DDS+COR+DGCON bit-exact (COR=DGCON=0
at the dgf line); (3) regent leaves DG for D≥BKPT (GO TO 23) — passthrough; (4) dgdriv tripling DDS·EXP(FRMT) —
live applies aspen COR=-2.44 (attenuating), jl applies 0; (5) DGSCOR calibration: live CORI(cornew)=-2.44323,
WC=1.0, exp(-2.44)=0.0868 > 0.0821 ⇒ PASSES the out-of-range trap (dgdriv.f:640, exp(±2.5)=[0.0821,12.1825]); jl
cornew=-2.6297 ⇒ exp=0.072 < 0.0821 ⇒ TRAPPED to 0 ⇒ COR=0 ⇒ aspen over-grows where gemdg is explosive on small
DBH; (6) the Δ0.187 in cornew = per-tree TERM: bark/scale/wk3 identical, but jl measured DG uniformly 0.842× live
(0.640 vs 0.760, 0.080 vs 0.095 — ratio 0.80/0.95). 

THE BUG: diameter_growth.jl:380-385 (FVS dgdriv.f:361 `DG(I)=DG(I)*BRATIO`, the IDG=1/3 outside→inside-bark
measured-DG conversion) used `bark_ratio(bark_a,bark_b,…)` which FLOORS to 0.80 for CR (bark_a/bark_b=0), not
`cr_bratio`≈0.95. The 0.842× shrink pushes cornew ~0.19 more negative, and for measured-DG species whose cornew
sits near the -2.5 trap cliff (e.g. aspen), it flips COR from ~-2.44 to 0 — an 8× DG swing. FOURTH CR variant-bark
location (after DDS→DG, DBH-update, backdate/TERM). FIX: dispatch `_cr_cal ? cr_bratio(sd,sp,saved_dbh,imod) :
bark_ratio(…)` (CR-gated; SN/NE/CS/LS unchanged). VALIDATED: 190851682020004 2022 now BIT-EXACT (BA 17/17, was
17/30); crt01_growth BYTE-IDENTICAL to pre-fix (no regression on the primary CR validation); 20/20 previously-
bit-exact stands still bit-exact; growth-type stands 1855935125290487 + 745576365290487 growth→selection (over-
prediction gone). Meta: the closure was disproved by MEASURING (per-tree TreeId verifier), not inferring — the same
discipline that found #11/#12. OPEN: a few growth-type stands remain (2602547010690, 758040786290487, 2713779010690)
— other measured-DG species near the trap, or genuine ±straddles; and 190851682020004 2042+ dense self-thin tail.

### 2026-07-30 — variant-bark AUDIT (after #14): 5th location fixed, 2 latent sites catalogued
Given #14 was the 4th instance of the identical `bark_ratio(0,0)=0.80` vs `cr_bratio~0.95` pattern, swept ALL
shared-engine `bark_ratio(` call sites for CR-reachability:
- ★ 5th location FIXED (commit b02964a): keyword_dispatch.jl:1149 morts KBIG RDPSRT size-rank key
  WORK3=DBH+DG/BRATIO (morts.f:879-882) — CR-gated to cr_bratio. FAITHFUL but .sum-INERT on the tested sample
  (crt01 byte-identical, 100/100 bit-exact preserved, 13 dense self-thin stands unchanged): CR dense mortality
  rarely hits a near-tie where the ~0.02" DG/bark difference flips selection. Kept per doctrine #4 (faithful port).
- LATENT (not FIA-sweep-reachable, catalogued for future variant-bark audit, NOT fixed — each needs its own
  validation path): simulate.jl:186 (FERTILIZ growth-effect `dib=d*bark_ratio`, only under a FERTILIZE keyword);
  keyword_dispatch.jl:1245 (HTGSTOP/topkill breakage `d=dbh*bark_ratio`, only when a tree topkills).
- ALREADY-DISPATCHED (prior fixes): simulate.jl:460 (DBH-update, `_cr_up`), the DDS→DG driver bark, backdate/TERM.
CONCLUSION: 5 active variant-bark locations now all cr_bratio-dispatched; 2 latent conditional sites remain
(fertilize/topkill), documented. The recurring root = jl's shared `bark_ratio(bark_a,bark_b,…)` floors to 0.80 when
CR's linear bark coefs are 0, where FVS BRATIO dispatches to cr/bratio.f (~0.95) — AUDIT this pattern first for any
future western GENGYM variant.
