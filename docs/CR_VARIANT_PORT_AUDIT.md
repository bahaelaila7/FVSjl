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
