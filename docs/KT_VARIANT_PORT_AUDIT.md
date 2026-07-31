# KT (Kootenai/Kaniksu/Tally Lake) variant port — audit log

Branch `kt-variant-port`. First of the western Rockies cluster (KT/IE/EM/BM/TT/UT). KT uses the STANDARD
WESTERN WYKOFF DDS (kt/dgf.f: DGLD/DGCR/DGCRSQ/DGDBAL/DGHAB(9)/DGFOR(7)/DGDS/DGEL/DGEL2/DGSASP, no CALL GEMDG)
— NOT CR's GENGYM. Each western-cluster variant = discount of BOTH the eastern SN/NE/CS/LS Wykoff-DDS engine
AND CR's western framework (crown/ccfcal/regent/htgf), plus western terms (elevation, slope-aspect, 9-group
habitat). Oracle: /workspace/.ktwork/FVSkt_clean (relink_kt.sh). Doctrine identical to CR (validate bit-exact
vs LIVE FVSkt per chunk; MEASURE don't infer; reuse shared engine).

KT infra: VARACD='KT', MAXSP=11, IFINT=10 (10-yr cycle), IFINTH=5, seed 55329. Species JSP (11):
WP WL DF GF WH RC LP ES AF PP OT (western Montana conifers).

## Chunk 0 — scaffold: DONE
- Branch kt-variant-port off cr-variant-port.
- Live oracle relinked + smoke-tested (/workspace/.ktwork/FVSkt_clean, banner "KOOTENAI, KANIKSU, TALLY LAKE").
- Kootenai <: AbstractVariant singleton (src/variants/kootenai/kootenai.jl): variant_code="KT", nspecies=11,
  htg_period=10. Registered in variant_from_code (KT|KOOTENAI). Included in FVSjl.jl. Exported.
- Validated: package precompiles clean; variant_from_code("KT")===Kootenai(); nspecies=11; htg_period=10;
  unimplemented hooks have NO Kootenai method (error loudly, doctrine #5). data/kootenai/ created (empty).
- NEXT: chunk 1 species (kt/blkdat.f JSP/FIAJSP/PLNJSP → data/kootenai/species_*.csv + load), then site/habitat,
  then the western-Wykoff dgf port (reuse eastern DDS structure + add DGEL/DGSASP/DGHAB terms + KT coefs).

## Chunk 1 — species + grinit infra: DONE
- 11-species identity CSV data/kootenai/species_coefficients.csv: WP/WL/DF/GF/WH/RC/LP/ES/AF/PP/OT
  (FIA 119/073/202/017/263/242/108/093/019/122/999; PLANTS PIMO3/LAOC/PSME/ABGR/TSHE/THPL/PICO/PIEN/ABLA/PIPO/2TREE).
- SPCTRN crosswalk data/kootenai/species_translation.csv (442 rows): extracted the shared western ASPT(442,21)
  table from kt/spctrn.f — KT target = COLUMN 12 (measured: kt/spctrn.f SELECT CASE(VAR) CASE('KT') ASPT(I,12);
  CR=col 8). Extractor .sweep_work/extract_kt_spctrn.jl. All 442 targets are valid KT species (11 distinct = all
  used); native rows verified (WP/119/PIMO3→WP, WL/073/LAOC→WL, DF/202/PSME→DF incl PSMEG/PSMEM subspecies).
- src/variants/kootenai/species.jl: init_blockdata! (mirrors CR) with KT grinit defaults (MEASURED, kt/grinit.f):
  YR=10, FINT=10 (:173), LHTDRG=.TRUE. (:101), DGSD=2.0 (:170), seed 55329 (blkdat:196); ★ LZEIDE=.FALSE. (:125)
  ⇒ zeide_sdi=false (KT uses STAGE SDI, UNLIKE CR's Zeide — a real KT/CR divergence). spctrn_column=4 (col-12
  target lands in the 7-col CSV's slot 4). other_species=11 (KT "OT", fia 999).
- Validated: coefficients(Kootenai()) loads 11 species w/ correct alpha+FIA; 442 translation rows; package
  precompiles clean. NEXT: chunk 2 site/habitat (kt/sitset.f + kt/habtyp.f — habitat-type groups feed DG), then
  the western-Wykoff dgf port (chunk 3).

## Chunk 2/3 SCOPING (measured from kt/habtyp.f, kt/sitset.f, kt/dgf.f)
KT habitat/site/DG structure, precisely traced:
- kt/habtyp.f: habitat-type CODE -> KKTYPE index. KTYPE(95) maps code-slot -> ITYPE(1..30) (e.g. 5*1,4*2,...);
  MTYPE(30) = representative habitat codes per ITYPE (130,170,250,...,999). (30 habitat-type groups.)
- kt/sitset.f: SDImax/BAMAX defaults. BAMAXA(30) = per-ITYPE max BA; SDIDEF(I)=BAMAX/(0.5454154*(PMSDIU/100)).
- kt/dgf.f DG-habitat flow (CHUNK 3, not 2): ISPHAB=MAPHAB(KKTYPE,ISPC) (MAPHAB(175,MAXSP) maps the raw
  habitat-type index -> a per-species 1..9 class); DGCON(ISPC)=DGHAB(ISPHAB,ISPC) (DGHAB(9,MAXSP) intercepts).
  DG coefficient arrays: DGLD/DGCR/DGCRSQ/DGDBAL/DGDS/DGEL/DGEL2/DGSASP + DGFOR(7,MAXSP) + MAPHAB + DGHAB.
=> CHUNK 2 = port habtyp (extract KTYPE(95)+MTYPE(30) -> code->KKTYPE) + sitset (extract BAMAXA(30) + SDIDEF
   formula). CHUNK 3 = western-Wykoff dgf (reuse eastern DDS structure + KT coefficient DATA + the
   KKTYPE->MAPHAB->ISPHAB->DGHAB habitat-intercept path + DGEL/DGEL2/DGSASP elevation/slope-aspect terms).
NOTE: validating chunks 2+ against live needs a KT harness (BIN["KT"]=/workspace/.ktwork/FVSkt_clean, VAR["KT"]=
Kootenai()) + a KT stand list (check the FIA DB for northern-ID/MT stands in KT geography) — a small harness
add mirroring the CR dig_dir wiring, do at the start of chunk 2's validation.

## KT-stands for validation: FIA DB has NO KT-labeled stands (use IE/EM geography through KT)
Checked FVS_STANDINIT_COND.VARIANT: ~20 variant codes present (SN/LS/CR/CS/NE/EM/UT/AK/SO/IE/WS/CI/EC/CA/WC/
PN/BM/TT/NC + missing) — KT is ABSENT (KT's N-Idaho/Montana geography is labeled IE (17808 stands) / EM (106962)
in the DB). NOT a blocker: for a bit-exact differential BOTH live FVSkt and jl Kootenai() get IDENTICAL stand
input + the KT variant, so IE/EM-geography stands (realistic KT-species trees: WP/WL/DF/GF/WH/RC/LP/ES/AF/PP)
run through KT are a valid differential. Chunk-2 harness plan: pull an IE (or EM) stand list from the FIA DB,
BIN["KT"]=/workspace/.ktwork/FVSkt_clean + VAR["KT"]=Kootenai(), reuse dig_dir/dig_batch/dig_vol verbatim
(they already take a stand list + variant). Species that fall outside KT's 11 map via the SPCTRN crosswalk
(already ported, col 12) to KT's OT/OS-equivalents.

## Chunk 2 — COMPLETE DATA MAP (all tables located; 2-level habtyp flow traced)
KT habtyp.f uses a TWO-INDEX habitat system (measured):
  code -> KODTYP; search KOTHAB(175): K where KODTYP<KOTHAB(K) -> KKTYPE=K-1  (the "Kootenai habitat type",
    1..175 — feeds DG: MAPHAB(KKTYPE,ISPC)->ISPHAB(1..9)->DGHAB, chunk 3);
  KODTYP=KOTHAB(KKTYPE); search JTYPE(95): K where KODTYP<JTYPE(K) -> ITYPE=KTYPE(K-1) (the "Inland Empire
    habitat type", 1..30); KODTYP=MTYPE(ITYPE) — ITYPE feeds sitset BAMAXA(ITYPE).
Tables + locations (all extractable):
  - KOTHAB(175): kt/blkdat.f:55  (Kootenai habitat codes, sorted; 10,65,70,74,79,91,92,93,95,100,110,...)
  - JTYPE(95):   kt/blkdat.f:160 (IE habitat codes: 10,100,110,130,140,160,170,180,190,200,...)
  - KTYPE(95):   kt/habtyp.f DATA KTYPE (code-slot -> ITYPE 1..30: 5*1,4*2,4,3*1,3,4,5,6,7,8,9,8,8,9,7,3,...)
  - MTYPE(30):   kt/habtyp.f DATA MTYPE (ITYPE -> IE rep code: 130,170,250,260,280,290,310,320,330,420,470,
                 510,520,530,540,550,570,610,620,640,660,670,680,690,710,720,730,830,850,999)
  - BAMAXA(30):  kt/sitset.f:32 (ITYPE -> maxBA: 140,220,250,310,240,270,310,310,200,310,290,330,380,440,500,
                 500,390,390,440,180,290,400,350,390,260,300,220,220,160,300)
  - SDIDEF(I) = BAMAX/(0.5454154*(PMSDIU/100.))  [BAMAX defaults to BAMAXA(ITYPE)]
CHUNK 2 TODO: extract KOTHAB/JTYPE (blkdat) + wire site_index.jl (Kootenai) computing KKTYPE (for DG) + ITYPE +
BAMAX/SDIDEF; add site_setup!(::Kootenai) hook (mirror CR site_index.jl); build KT harness (IE stand list +
BIN["KT"]/VAR["KT"]); validate KKTYPE/ITYPE/SDImax vs live FVSkt (instrument like CR). Then chunk 3 dgf.

## Chunk 2 — GROUND TRUTH captured + tables HAND-VERIFIED vs live FVSkt
KT harness: 40 IE-geography stands (.sweep_work/kt_ie_stands.txt), run through KT via /workspace/.ktwork/
FVSkt_clean (IE stands are valid for KT differential — identical input + KT variant both sides). Live probe
.sweep_work/kt_live_probe.jl.
Stand 753200841290487 live FVSkt .out:
  KOOTENAI HABITAT TYPE IS 531   (= KOTHAB(KKTYPE))
  INLAND EMPIRE HABITAT TYPE IS 530  (= MTYPE(ITYPE))
  SDI MAX = 949 for all 11 species
HAND-VERIFICATION against the extracted tables (doctrine #2 — the data reproduces live):
  MTYPE(14)=530  => ITYPE=14
  BAMAXA(14)=440
  SDIDEF = BAMAX/(0.5454154*(PMSDIU/100)) = 440/(0.5454154*0.85) = 949.1 -> 949  (PMSDIU=85 default) ✓
=> The chunk-2 MTYPE/BAMAXA tables + SDIDEF formula are CONFIRMED correct against live output. Remaining
   chunk-2 execution: extract KOTHAB(175)/JTYPE(95) from kt/blkdat.f, write site_index.jl(Kootenai) implementing
   code->KKTYPE (KOTHAB search) -> KOTHAB(KKTYPE) -> JTYPE search -> ITYPE=KTYPE(K-1) -> MTYPE(ITYPE) + BAMAX/
   SDIDEF, wire site_setup!(::Kootenai), and diff jl's KKTYPE/ITYPE/SDImax vs these live values across the 40
   IE stands. (Also: locate each stand's raw INPUT habitat code — the FIA stand record field feeding habtyp.)

## Chunk 2 — tables VALIDATED (habitat_tables.jl); site_index.jl wiring is next
habitat_tables.jl committed + programmatically validated vs live (KODTYP=531 -> ITYPE=14 -> IE=530 -> SDI=949,
all matching FVSkt). Plot fields identified (src/core/state.jl): habitat_code=KODTYP (input), habitat_input=ITYPE,
valid_habitat=JTYPE(122), sp_sdi_def=SDIDEF, sp_site_index=SITEAR, model_type=IMODTY, forest_idx=IFOR.
OPEN DESIGN Q (couples chunk 2<->3): KT's DG needs KKTYPE (1..175, the Kootenai habitat index feeding
MAPHAB(KKTYPE,ISPC)->DGHAB), but the shared plot has no KKTYPE field (only habitat_input=ITYPE 1..30). Options:
(a) store KKTYPE in habitat_input for KT (repurpose; check no shared code reads it as ITYPE), (b) add a KT field,
(c) recompute KKTYPE in the DG chunk from habitat_code. Decide when writing chunk 3 dgf. site_index.jl(Kootenai)
TODO: site_setup! = KOTHAB search(habitat_code)->KKTYPE; KOTHAB[KKTYPE]->JTYPE search->ITYPE; BAMAX default
BAMAXA[ITYPE]; SDIDEF=BAMAX/(0.5454154*pmsdiu) into sp_sdi_def; + KT site-index (SITEAR) per kt/sitset.f site
curves (may defer to height chunk). Validate jl KKTYPE/ITYPE/sp_sdi_def vs live across the 40 IE stands
(.sweep_work/kt_ie_stands.txt) — need each stand's raw input habitat code + a setup-only entry point.

## Chunk 2 — INPUT confirmed (PV_CODE) + full flow validated end-to-end with real input
FIA FVS_STANDINIT_COND habitat column = PV_CODE (also PV_FIA_HABTYPCD1). Stand 753200841290487: PV_CODE=531,
PV_REF_CODE=110. Live FVSkt read it from the DB (build_subdb copies the stand; FVS reads FVS_STANDINIT directly)
=> KOOTENAI habitat 531. FULL FLOW now validated with the REAL input:
  PV_CODE=531 -> plot.habitat_code=531 -> KOTHAB search (531 in KOTHAB) -> KKTYPE (KOTHAB[KKTYPE]=531, =live 531)
  -> JTYPE search(531) -> ITYPE=14 -> MTYPE[14]=530 (=live IE 530) -> BAMAXA[14]=440 -> SDIDEF=949 (=live).
REMAINING chunk-2 wiring (precise):
  (1) jl FIA loader must read PV_CODE -> plot.habitat_code (CR didn't need it — habitat-input was a documented
      gap, keyword_dispatch.jl:630 "non-zero habitat ignored"; the LIVE side reads PV_CODE from the DB
      automatically, so jl must match). Also handle PV_REF_CODE (the CPVREF/PVREF1 path in habtyp for
      reference-code stands) — 110 here; check if it changes KODTYP.
  (2) src/variants/kootenai/site_index.jl: site_setup!(::Kootenai) = KOTHAB search -> KKTYPE (store for DG),
      -> ITYPE -> BAMAX/SDIDEF into sp_sdi_def; include habitat_tables.jl + site_index.jl in FVSjl.jl.
  (3) validate jl KKTYPE/ITYPE/sp_sdi_def vs live (531/ITYPE14/949) across the 40 IE stands.

## Chunk 2 — site_index.jl IMPLEMENTED + unit-validated
src/variants/kootenai/site_index.jl: kt_habtyp(kodtyp) does the 2-level KOTHAB->KKTYPE + JTYPE/KTYPE->ITYPE
mapping; kt_site_index_setup! stores KKTYPE in plot.habitat_code (chunk-3 DG), ITYPE in plot.habitat_input,
computes SDIDEF=BAMAX/(0.5454154*pmsdiu) into sp_sdi_def; site_setup!(::Kootenai) hook. Included in FVSjl.jl.
UNIT-VALIDATED vs live: kt_habtyp(531) -> KKTYPE=87 (KOTHAB[87]=531 = live KOOTENAI 531), ITYPE=14 (MTYPE=530 =
live IE), BAMAXA[14]=440, SDIDEF=949 (= live SDI MAX). Package precompiles; site_setup! method present.
REMAINING for full stand-level integration (bundles with chunk 3, since a full jl stand run needs DG):
  jl FIA loader must read PV_CODE -> plot.habitat_code so site_setup! sees the right input (CR left this a gap).
  Then end-to-end harness diff (grow regime) once DG/height/crown/mort are ported. Site-index SITEAR curves
  (kt/sitset.f site-index-by-species) deferred to the height/DG chunk where they're consumed.

## Chunk 3 (dgf — western Wykoff DDS) — FULLY MEASURED (doctrine #2), ready to execute
The KT large-tree DDS (kt/dgf.f:341-344), a standard western Wykoff regression:
  DDS = CONSPP + DGLDS*ALD + CR*(DGCRS + CR*DGCRS2)
      + DGDBLS*BAL/ln(D+1) + DGCCF2*CCF2 + DGD2*D*D
      + DGLBAS*ln(BA) + DGPC1*PCCF1 + DGPC2*PCCF1
  DDS = max(DDS, -9.21);  WK2 = DDS
  where CONSPP = DGCON(ISPC) + COR(ISPC) + DGCCF(ISPC)*RELDEN   [sp11(OT): 0.01*DGCCF]
        ALD = ln(D);  BAL = (1-PCT/100)*BA   [sp11: /100]
        DGCON (per-species, from DGCONS entry): DGHAB(ISPHAB) + DGFOR(ISPFOR) + DGEL*ELEV + DGEL2*ELEV^2
              + (DGSASP*sin(ASPECT)+DGCASP*cos(ASPECT)+DGSLOP)*SLOPE + DGSLSQ*SLOPE^2 + ln(COR2)
              [ISPHAB=MAPHAB(KKTYPE,ISPC), ISPFOR=MAPLOC(KOTFOR,ISPC)]
Coefficient arrays to extract from kt/dgf.f DATA blocks (all length 11=MAXSP unless noted):
  DGLD, DGCR, DGCRSQ, DGDBAL, DGDS(->DGD2 per-spc), DGCCFA(->DGCCF), + the CCF2/lnBA/PCCF1 coefs
  (DGCCF2/DGLBAS/DGPC1/DGPC2 — locate their per-species source arrays), DGEL, DGEL2, DGSASP, DGCASP,
  DGSLOP, DGSLSQ; 2D: DGHAB(9,11), DGFOR(7,11), MAPHAB(175,11), MAPLOC(10,11), OBSERV(9,11).
EXEC PLAN: (1) extract all arrays (compact — 11 species); (2) write diameter_growth!(::Kootenai) reusing the
shared stand-stat computation (BAL/CCF/RELDEN/PCCF1/CR — same as eastern engine) + the validated
KKTYPE->MAPHAB->DGHAB path (chunk-2 site_index stored KKTYPE in habitat_code); (3) wire loader PV_CODE->
habitat_code; (4) per-tree WK2 instrument-replay vs live FVSkt (relink_kt.sh instrumented dgf.o) on the 40 IE
stands — same recipe as CR's DGFTRC. This is the LARGE core chunk (careful transcription + per-tree diff).

## Chunk 3 — DATA block inventory + layouts CONFIRMED (extraction is mechanical)
All kt/dgf.f DATA blocks located (line): DGLD(93) DGCR(96) DGCRSQ(99) DGDBAL(102) CCFSQ(105) DGPCC1(108)
DGPCC2(111) DGLBA(114) DBHCH(117) ICRLIM(119) OBSERV(121) MAPHAB(138) DGHAB(180) DGCCFA(206) MAPLOC(213)
DGFOR(225) DGDS(241) DGCASP(257) DGSASP(260) DGSLOP(263) DGSLSQ(266) DGEL(269) DGEL2(272).
Scalar per-species sources (dgf.f:294-321): DGCCF2=CCFSQ(ISPC); DGLBAS=DGLBA(ISPC); DGD2=DGDS(ISPC);
DGPC1=DGPCC1(ISPC)*DUM1; DGPC2=DGPCC2(ISPC)*DUM2; DGCCF(ISPC)=DGCCFA(ISPC); CONSPP includes DGCCF*RELDEN.
LAYOUTS confirmed: 1D arrays = 11 values (e.g. DGLD = 0.89068,0.71363,...,0.89778). 2D arrays fill column-major
= grouped by SPECIES (each species gets its inner-dim values in order): DGHAB(9,11)=11 groups of 9 (habitat-
class intercepts); DGFOR(7,11)=11 groups of 7; MAPHAB(175,11)=11 groups of 175 (INTEGER KKTYPE->ISPHAB map);
MAPLOC(10,11)=11 groups of 10; OBSERV(9,11)=11 groups of 9. E.g. DGHAB species1(WP) classes1-9 =
0,0.46877,0.36827,0.25380,0,0,0,0,0. DGFOR species1 = 1.76061,0,0,0,0,0,0.
=> Extraction now purely mechanical (parse DATA -> tokens -> reshape by species). Then diameter_growth!(::Kootenai)
implements the DDS eqn (dgf.f:341) reusing shared stand-stats + KKTYPE->MAPHAB->ISPHAB->DGHAB (chunk-2 stored
KKTYPE in habitat_code) + KOTFOR->MAPLOC->ISPFOR->DGFOR + elev/slope-aspect DGCON terms; validate per-tree WK2
vs live via relink_kt.sh instrumented dgf.o (CR DGFTRC recipe) on the 40 IE stands.

## Chunk 3 — coefficients EXTRACTED+VERIFIED; diameter_growth! implementation scoped (LS template)
src/variants/kootenai/dg_coefficients.jl: all 23 arrays parsed+spot-verified vs source (DGLD/DBHCH(repeat-
syntax)/ICRLIM/DGHAB[,1]/DGFOR[,1] exact). The transcription risk (the big chunk-3 hazard) is ELIMINATED.
IMPLEMENTATION TEMPLATE = LS (src/variants/lakestates/diameter_growth.jl) — KT is Wykoff like LS/CS, NOT CR
GENGYM. LS provides: dgf!(s, ::LakeStates) fills scratch.wk[2,i]=ln(inside-bark DDS) per tree; ls_dgcons!(s)
computes per-species dg_const. KT dgf!(::Kootenai) TODO (mirror LS):
  kt_dgcons!: dg_const[sp] = DGHAB[MAPHAB[KKTYPE,sp], sp] + DGFOR[MAPLOC[KOTFOR,sp], sp] + DGEL[sp]*ELEV +
    DGEL2[sp]*ELEV^2 + (DGSASP[sp]*sin(ASP)+DGCASP[sp]*cos(ASP)+DGSLOP[sp])*SLOPE + DGSLSQ[sp]*SLOPE^2
    (+ ln(COR2) calib). KKTYPE from plot.habitat_code (chunk-2 stored it there).
  dgf! per tree: DDS = dg_const[sp] + COR[sp] + DGCCFA[sp]*RELDEN(sp11:0.01*) + DGLD[sp]*ln(D) +
    CR*(DGCR[sp]+CR*DGCRSQ[sp]) + DGDBAL[sp]*BAL/ln(D+1) + CCFSQ[sp]*CCF2 + DGDS[sp]*D^2 + DGLBA[sp]*ln(BA) +
    (DGPCC1[sp]*DUM1+DGPCC2[sp]*DUM2)*PCCF1; BAL=(1-PCT/100)*BA (sp11:/100); clamp>=-9.21;
    inside-bark: diagro=sqrt(D^2+exp(DDS))-D, then WK2=ln((D*BR+diagro)^2-(D*BR)^2) using BR=KT_BKRAT[sp].
  Need from shared stand-stats: RELDEN, CCF2, PCCF1, DUM1/DUM2, BA, PCT, CR, D — confirm availability/definitions
  vs LS dgf! (BAGE5/BALC etc. may differ; KT uses raw BA + PCCF1). Then wire loader PV_CODE->habitat_code +
  per-tree WK2 instrument-replay vs live FVSkt (relink_kt.sh instrumented dgf.o) on 40 IE stands.

## Chunk 3 — dgf! stand-stats traced; KEY: WK2 = OB DDS (no bark conversion, unlike LS)
Traced kt/dgf.f per-tree stand-stats (doctrine #2):
  - WK2(I)=DDS DIRECTLY (clamp >= -9.21), NO OB->IB bark conversion in dgf.f (UNLIKE LS which converts).
    => KT dgf! is simpler; WK2 = the outside-bark DDS as computed. (BKRAT is for bratio/other chunks, not dgf.)
  - CONSPP = DGCON(ISPC) + COR(ISPC) + DGCCF(ISPC)*RELDEN   [sp11: 0.01*DGCCF*RELDEN]  (RELDEN=stand rel density)
  - CCF2 = RELDEN*RELDEN (dgf.f:324 — NOT crown-competition²; the CCFSQ term is CCFSQ(ISPC)*RELDEN²)
  - PCCF1 = PCCF(IPCCF) (point CCF at the tree's point; IPCCF = point index)
  - DUM1/DUM2 (dgf.f:314-318): default DUM1=0,DUM2=1; under a condition (~:316) DUM1=1,DUM2=0 -> a SWITCH between
    DGPCC1 and DGPCC2 coefficients (DGPC1=DGPCC1*DUM1, DGPC2=DGPCC2*DUM2; both multiply the SAME PCCF1).
  - ALD=ALOG(D); BAL=(1-PCT/100)*BA [sp11:/100]; CR=crown ratio; BA=stand basal area.
REMAINING TRACE before coding dgf!(::Kootenai): (1) RELDEN exact formula (computed before the tree loop in
kt/dgf.f — likely RELSDI/relative-SDI or a density ratio; find its assignment); (2) PCCF point array + IPCCF
index (shared point-CCF — check the engine's per-point CCF); (3) the DUM1/DUM2 condition at :316. Then implement
kt_dgcons! (DGCON from DGHAB/DGFOR/elev/slope-aspect via KKTYPE=habitat_code + KOTFOR) + dgf! (DDS eqn, WK2=DDS,
no bark) + loader PV_CODE->habitat_code + per-tree WK2 instrument-replay vs live on 40 IE stands.

## Chunk 3 — MEASUREMENT COMPLETE (all stand-stats sourced); dgf! is now a pure coding task
Final stand-stat sources traced:
  - RELDEN = shared dense.f (RELDEN=RELDM1 dense.f:261, or RELDT :247) — the stand RELATIVE DENSITY. NOT
    KT-specific; the FVSjl engine already computes it (CR used dense.f RELDM1 for crown/REGCAL). Find the jl
    field (compute_density! output — the relative-density scalar CR reads).
  - PCCF1 = PCCF(ITRE(I)) — POINT CCF at the tree's inventory point (ITRE(I)=tree's point). Shared point-CCF
    (the engine's per-point CCF, as CR estab used density.point_ccf[plot_id]).
  - DUM switch = MANAGD (managed-stand flag): MANAGD==1 -> DGPCC1 else DGPCC2 (both x PCCF1).
  - BA=p.basal_area; PCT=crown percentile (t.crown_ratio); CR=crown ratio (t.crown_pct); D=t.dbh; ELEV/SLOPE/
    ASPECT from plot.
FULL dgf!(::Kootenai) SPEC (all inputs now sourced):
  kt_dgcons!: dg_const[sp] = DGHAB[MAPHAB[KKTYPE,sp],sp] + DGFOR[MAPLOC[KOTFOR,sp],sp] + DGEL[sp]*ELEV
    + DGEL2[sp]*ELEV^2 + (DGSASP[sp]*sin(ASP)+DGCASP[sp]*cos(ASP)+DGSLOP[sp])*SLOPE + DGSLSQ[sp]*SLOPE^2
    (+ ln(COR2)); KKTYPE=plot.habitat_code (ch2), KOTFOR=forest index.
  dgf! per tree: DDS = dg_const[sp] + COR[sp] + DGCCFA[sp]*RELDEN (sp11: 0.01*) + DGLD[sp]*ln(D)
    + CR*(DGCR[sp]+CR*DGCRSQ[sp]) + DGDBAL[sp]*BAL/ln(D+1) + CCFSQ[sp]*RELDEN^2 + DGDS[sp]*D^2 + DGLBA[sp]*ln(BA)
    + (MANAGD? DGPCC1[sp] : DGPCC2[sp])*PCCF(point); BAL=(1-PCT/100)*BA (sp11:/100); clamp>=-9.21; WK2=DDS (no bark).
=> CHUNK-3 MEASUREMENT COMPLETE. Remaining = pure coding (map jl field names for RELDEN/point-CCF/MANAGD/ELEV/
SLOPE/ASPECT) + kt_dgcons!/dgf! + loader PV_CODE->habitat_code + per-tree WK2 instrument-replay vs live.

## Chunk 3 — jl field mapping (dgf! inputs -> StandState); RELDEN=RELDM1 traced to dense.f
jl fields for dgf!(::Kootenai) (from src/core/state.jl):
  PCCF1  -> p.point_ccf[ITRE] (per-point PCCF, dense.f:210)     MANAGD -> p.managed (0/1)
  ELEV   -> p.elevation (hundreds ft)   SLOPE -> p.slope (0..1)   ASPECT -> p.aspect (radians)
  BA     -> p.basal_area   D -> t.dbh   CR -> t.crown_pct   PCT -> t.crown_ratio
  KKTYPE -> p.habitat_code (chunk-2 stored)   KOTFOR -> forest index (p.forest_idx or KOTFOR map)
  RELDEN -> RELDM1: dense.f RELDT=Σ RELDSP(ISPC) (per-species relative density), RELDM1=RELDT, then backdate
            adjust (LBKDEN: TEMP1=(RELDEN-RELDM1)*FINTH/FINT+RELDM1). The stand relative-density scalar from
            compute_density!. NOTE: CR crown reads relden=stand_ccf(s) — CONFIRM whether the jl engine's
            relative-density (RELDM1) == stand_ccf or a distinct helper; KT dgf uses RELDM1 specifically (and
            CCF2=RELDM1^2). This is THE last field to pin (identify/expose RELDM1 in the jl density engine).
REMAINING chunk-3 = (1) confirm/expose RELDM1 jl scalar; (2) code kt_dgcons! + dgf!(::Kootenai) per the full spec;
(3) loader PV_CODE->plot.habitat_code; (4) per-tree WK2 instrument-replay vs live FVSkt (relink_kt.sh dgf.o) on
the 40 IE stands (.sweep_work/kt_ie_stands.txt) — CR DGFTRC recipe.
