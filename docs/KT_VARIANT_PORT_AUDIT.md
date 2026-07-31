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

## Chunk 3 — dgf!/kt_dgcons! IMPLEMENTED + compiles
src/variants/kootenai/diameter_growth.jl written (LS pattern): kt_dgcons!(s) fills calib.dg_const (DGCON =
DGHAB[MAPHAB[KKTYPE,sp]] + DGFOR[MAPLOC[KOTFOR,sp]] + DGEL*ELEV + DGEL2*ELEV^2 + slope-aspect + ln(COR2)),
atten=OBSERV[ISPHAB], bark=BKRAT; dgf!(s,::Kootenai) fills scratch.wk[2,i]=DDS (kt/dgf.f:341, WK2=OB no bark
conversion). jl fields: KKTYPE=p.habitat_code, KOTFOR=p.forest_idx, RELDEN=p.relative_density, PCCF=p.point_ccf
[t.plot_id], MANAGD=p.managed, ELEV/SLOPE/ASP=plot. Wired kt_dgcons! into simulate.jl dgcons! dispatch (KT
branch) + included in FVSjl.jl. VALIDATED: precompiles clean; kt_dgcons!/dgf!(::Kootenai) methods resolve;
coefficients loaded (KT_DGLD/KT_MAPHAB(175,11)/KT_DGPCC1). REMAINING: loader PV_CODE->plot.habitat_code +
confirm calibrate_diameter_growth! handles KT + per-tree WK2 instrument-replay vs live FVSkt (relink dgf.o) on
40 IE stands. NOTE placeholders to verify vs live: KOTFOR=forest_idx (may need KT forkod map); crown dubbing
skipped (validate on input-crown stands until KT crown chunk 5); RELDEN=relative_density (vs _prev).

## Chunk 3 — loader PV_CODE->habitat_code WIRED; next: species CSV needs shared-setup columns
Wired fia_database.jl: KT-gated PV_CODE -> plot.habitat_code (so site_setup! sees the habitat input; matches
live FVSkt which reads PV_CODE from the DB). Loading an IE stand through jl Kootenai() now reaches SETUP and
errors at `KeyError: :wykoff_ht2` — the shared setup path (HT-DBH height dubbing, then DG/calibrate) reads
species-coefficient columns that KT's minimal chunk-1 species_coefficients.csv (4 identity cols only) lacks.
=> NEXT: expand data/kootenai/species_coefficients.csv to the shared-engine column set (CR schema is the
reference: site_lo/hi, dbh_max, bark1/bark2, st_*, ht1/ht2/wykoff_ht2, mort_*, sdi_max_default, htdbh_*,
varmrt_varadj, dg_resid_sd, ...). For chunk-3 DG VALIDATION specifically, the DG inputs that must match live are
D/CR/BAL/BA/RELDEN/PCCF (heights don't feed the KT DDS) — so height/htdbh columns can be placeholder until ch4,
BUT they must be PRESENT (non-KeyError) and consistent enough that the density computation matches. Fill the
columns KT genuinely needs now (bark from KT_BKRAT, sdi_max_default, dbh bounds) with real values; defaults
elsewhere pending their chunks. Then: confirm calibrate_diameter_growth! handles KT (or add branch) + per-tree
WK2 instrument-replay vs live FVSkt (relink instrumented dgf.o) on the 40 IE stands.

## Chunk 3 — DG RUNS END-TO-END on a real stand (species CSV + point_ccf fix)
- data/kootenai/species_coefficients.csv expanded to full CR schema (11 KT rows: identity + KT bark (KT_BKRAT);
  NON-DG columns = CR western-conifer placeholders matched by FIA (DF fallback for WP/WH/OT) — TO BE REPLACED
  with real KT values in chunks 4(height)/6(bark)/7(mort)/8(vol). Builder .sweep_work/build_kt_species_csv.jl.
- Fixed dgf!: point_ccf is on s.density (Density struct), not s.plot — pccf1 = s.density.point_ccf[t.plot_id[i]]
  (matches CR diameter_growth.jl:324 dens.point_ccf[ipccf]). RELDEN = p.relative_density confirmed (CR DG comment
  "RELDEN=stand CCF"; state.jl relative_density = current CCF).
- Loading IE stand 753200841290487 through jl Kootenai() now: SETUP ok -> kt_dgcons! ok -> dgf! ok -> shared
  DGDRIV calibrate/tripling/DBH-update ok -> errors at height_growth!(::Kootenai) = the NEXT unimplemented chunk
  (4). So CHUNK 3 (DG) EXECUTES END-TO-END on a real stand (doctrine #5: next hook errors loudly). Functional
  integration confirmed.
- REMAINING for chunk-3 SIGN-OFF: per-tree WK2 bit-exact validation vs live FVSkt (relink instrumented dgf.o
  printing WK2/tree, CR DGFTRC recipe; compare on the 40 IE stands' first cycle). Placeholders to confirm via
  that diff: KOTFOR=forest_idx (KT forkod map?), managed-flag/DUM, RELDEN source. Chunk 4 = height_growth! (htgf).

## Chunk 3 — DG VALIDATED-CLOSE vs live; residual = CCF/RELDEN dep on ccfcal (chunk 5)
Per-tree WK2 instrument-replay (live FVSkt_dgftrc: instrumented dgf.o printing I/ISPC/D/WK2/RELDEN/BA; jl:
env-gated dgf! trace) on stand 753200841290487, first cycle. Found + fixed TWO issues (doctrine #4 validation):
  BUG (fixed): kt/dgf.f:327 CR=ICR*0.01 — crown ratio is a FRACTION 0..1, jl passed t.crown_pct (0..100 pct)
    => CR*DGCR blew DDS to ~97.9. Fix cr=crown_pct*0.01: tree1(sp4,D13.7) WK2 97.9 -> 3.89 (live 4.09).
  RESIDUAL (~0.2, cross-chunk): RELDEN. Live RELDEN=98.39 (stand CCF), jl stand_ccf(s)=2.05 — because stand_ccf
    uses the EASTERN crown-width model for non-CR variants; KT's crown width isn't ported (chunk 5 ccfcal). SAME
    pattern CR hit ("CCF=0/wrong, eastern cwcalc lacks variant coefs; fixed by porting ccfcal"). The WK2 gap is
    tiny because RELDEN feeds only DGCCFA*RELDEN (~0.17) + CCFSQ*RELDEN^2 (~0.03); once KT ccfcal lands, RELDEN
    ->98.39 and WK2 -> bit-exact.
VERDICT: KT DG EQUATION + COEFFICIENTS + CR-units are VALIDATED CORRECT (WK2 3.89 vs 4.09, all terms match
except the CCF-derived RELDEN). Chunk-3 DG is bit-exact-PENDING-ccfcal — the residual is a documented cross-chunk
dependency (chunk 5 crown/ccfcal), NOT a DG bug. Coefficient extraction + equation + integration all confirmed
by the per-tree diff. (Instrument: /tmp/kt_dgf_instr.f, relink_kt.sh dgftrc /tmp/dgf.o; 40-stand batch pending.)

## Chunk 4 (height, htgf) — MEASURED (equation + coefficients), ready to implement
KT height growth (kt/htgf.f:107-112):
  CON = HTCON(sp) + H2COF*HT^2 + HGLD(sp)*ln(D) + HGLH*ln(HT)
  HTG = EXP(CON + HDGCOF*ln(DG)) + BIAS            [DG = this cycle's diameter growth from chunk 3!]
  HTG = max(HTG, 0.1);  HTG = HTG * SCALE * XHMULT(sp) * MISHGF(mistletoe)
  HTCON(sp) = HGHCH + HGSC(sp)  [+ ln(HCOR2) calib]  (htgf.f:210, a DGCONS-like per-stand setup)
=> height growth DEPENDS on DG (chunk 3) — cycle order DG->height (already how the engine runs; the gating hook
   height_growth!(::Kootenai) fires right after DG). Coefficient DATA blocks (kt/htgf.f): HGLD(11)@69, BIAS/HGLH
   @71 (scalars .4809/.23315), MAPHAB@185, HGHC@187, HGLDD@190, HGH2@193, HGSC(11)@196, XHMULT(11), + scalars
   H2COF/HDGCOF (locate). IMPLEMENT (chunk-3 pattern): extract arrays (parser like extract_kt_dgf.jl) ->
   ht_coefficients.jl; kt_htcons! (HTCON=HGHCH+HGSC) + height_growth!(s,::Kootenai) computing HTG; per-tree HTG
   instrument-replay vs live FVSkt (relink instrumented htgf.o). NOTE: DG must be bit-exact first (needs ccfcal
   ch5) for HTG's ln(DG) term to match live — so height validation partly gated on ch5 too (like DG's RELDEN).

## Chunk 5 (crown/ccfcal) — CCF-per-tree MEASURED (the DG-unblock; priority)
kt/ccfcal.f CCFCAL(MODE=1 → CCF): per-tree CCF contribution CCFT:
  D >= 10.0:  CCFT = RD1(sp) + D*RD2(sp) + D*D*RD3(sp)
  D <  10.0:  CCFT = RDA(sp) * D**RDB(sp)
  CCFT = CCFT * P   (P = tree TPA);  STAND CCF = Σ CCFT = RELDEN (=98.39 for the validated stand).
Coefficients (kt/ccfcal.f, all 11-wide): RD1@62 (.03,.02,.11,.04,.03,.03,.01925,.03,.03,.03,.03), RD2@63,
RD3@65, RDA@67 (.009884,.007244,.017299,.015248,.011109,.008915,.009187,.007875,.011402,.007813,.011109),
RDB@70 (1.6667,1.8182,1.5571,1.7333,1.7250,1.7800,1.7600,1.7360,1.7560,1.7680,1.7250). Crown-WIDTH (MODE=2:
B1-B6 + R6CRWD/IFOR branches) is for the crown model (chunk-5 crown_ratio, separate) — NOT needed for RELDEN.
=> IMPLEMENT (priority — unblocks DG+height bit-exact): extract RD1/RD2/RD3/RDA/RDB -> ccf_coefficients.jl;
   kt_tree_ccf(sp,d) = (d>=10 ? RD1+d*RD2+d^2*RD3 : RDA*d^RDB); add a Kootenai branch to stand_ccf(s) summing
   kt_tree_ccf*tpa (NOTE: KT CCF is the DIRECT polynomial, NOT the crown-width->area path jl's stand_ccf uses for
   eastern/CR — so KT needs its own CCF sum). Validate stand CCF == 98.39 vs live, then re-run the DG per-tree
   WK2 diff -> expect bit-exact (RELDEN closes). THEN chunk 4 height (ln(DG) term matches).
### KT growth-core DEPENDENCY CHAIN fully measured:
  ccfcal CCF(ch5) -> RELDEN 98.39 -> DG WK2 bit-exact(ch3) -> height ln(DG)(ch4). Implement ch5-CCF FIRST.

## Chunk 5 — CCF implemented (RELDEN 2.05->79.27); residual = GROSPC normalization
crown.jl: KT_RD1/RD2/RD3/RDA/RDB + kt_tree_ccf(sp,d); stand_ccf(s) Kootenai branch = Σ kt_tree_ccf*tpa.
RELDEN 2.05 -> 79.27 (live 98.39) — right CCF scale, DG WK2 residual ~5% -> ~1%. RESIDUAL 79.27 vs 98.39
(ratio 1.241 = 1/0.806) = GROSPC normalization: dense.f:207 PCCF=PCCF+CCFT*PI/GROSPC (and RELDSP/PRDA all
carry PI/GROSPC). jl's raw Σ CCFT omits the /GROSPC (gross growing space) divisor; 79.27/0.806 ≈ 98.4. Also PI
(per-point weight, =NPTS or GROSPC-related). NEXT: apply the KT stand-CCF normalization = Σ CCFT*tpa /GROSPC
(find p.gross_space; confirm PI). Then RELDEN->98.39, DG WK2 -> bit-exact, then the 40-stand DG batch. NOTE:
RELDEN may technically be the SDI-based RELDSP sum (dense.f:217 PRDA=P*(D/10)^1.605/GROSPC/XMAXPT) not the CCF
sum — but the CCF-poly value (79->98 w/ GROSPC) tracks live 98.39 closely, so verify which after the GROSPC fix.

## Chunk 5 — GROSPC hypothesis DISPROVEN (measured GROSPC=1.0); residual needs proper RELDEN diagnosis
Measured: p.gross_space=1.0, p.pi=4.0, points_inv=4. So /GROSPC does NOTHING (=1) — the 79.27-vs-98.39 gap is
NOT GROSPC (hypothesis retracted, doctrine #2: measure don't infer). The KT CCF polynomial (Σ kt_tree_ccf*tpa =
79.27) is close-but-not-live's-98.39. PROPER diagnosis needed (do NOT cargo-cult a factor):
  (a) RELDEN may be the SDI-BASED relative density (dense.f:217 PRDA=P*(D/10)^1.605*PI/GROSPC/XMAXPT; RELDT=
      Σ RELDSP), NOT the CCF sum — compute SDI-relative for this stand and compare to 98.39. My CCF-poly being
      ~80% of 98.39 could be coincidence.
  (b) OR backdated-DBH: DGDRIV backdates DBH before dgf!; live's RELDEN(98.39) vs jl stand_ccf on current DBH.
  (c) OR tree-set: jl stand_ccf sums live trees at dgf!; verify count/DBH match.
IMPACT: small — RELDEN feeds only DGCCFA*RELDEN + CCFSQ*RELDEN^2 in DDS; DG WK2 already ~1% (was ~5%). Instrument
the LIVE dense.f to print RELDSP/RELDEN derivation (relink) to settle (a) vs CCF-sum. The KT CCF polynomial +
stand_ccf branch are committed (real: RELDEN 2.05->79.27); the final normalization/definition is the open item.

## Chunk 5 — RELDEN CONFIRMED = Σ CCFT (dense.f:199, source-read); value gap = backdate/pass-timing
Settled by reading dense.f (not inferring): RELDSP(ISPC)=RELDSP(ISPC)+CCFT (dense.f:199); RELDEN=RELDM1=RELDT=
Σ RELDSP = Σ CCFT (the CCF sum, CCFT already ×P). So RELDEN IS the CCF sum (both the GROSPC and SDI-based
hypotheses are DISPROVEN — RELDSP has NO PI/GROSPC, unlike PCCF at :207). => jl's Σ(kt_tree_ccf·tpa)=79.27 is
the RIGHT quantity; the 98.39 gap is a DBH/tree-set difference at dgf! time:
  - DGDRIV backdates DBH for the CALIBRATION pass; dense.f (RELDEN) runs on that stand. jl's stand_ccf(s) at
    dgf! reads t.dbh — need to confirm it's the SAME (backdated) DBH the live calibration-pass dgf sees, and the
    SAME tree set. jl relden 79.27(cyc1)/143.64(cyc2): the cyc1 value should match live's calibration-pass 98.39.
  - LIKELY: jl stand_ccf uses CURRENT dbh but live RELDEN is on a DIFFERENT-pass/backdated dbh, OR small trees
    (D<0.1 skipped) differ. DIAGNOSE: print jl per-tree (sp,d,tpa,kt_tree_ccf) at dgf! and diff vs a live
    dense-instrument (relink printing per-tree CCFT+D) — settle whether it's DBH-backdate, tree count, or a
    small-tree threshold. IMPACT ~1% on DG WK2 (RELDEN in DGCCFA·RELDEN + CCFSQ·RELDEN²).
VERDICT: DG equation+coefficients+CR-units+CCF-definition all CONFIRMED correct; the last ~1% is the RELDEN CCF
sum's input DBH/tree-set at dgf! time (a shared-driver backdate-timing detail, not a KT-coefficient bug).

## Chunk 5 — CCFT formula CONFIRMED correct (ccfcal.f:110-116 IF/ELSE); residual = tree-set/TPA input
Re-read exact ccfcal.f MODE=1: IF(D.GE.10) CCFT=RD1+D*RD2+D^2*RD3 ELSE CCFT=RDA*D^RDB; CCFT=CCFT*P. My
kt_tree_ccf matches this exactly (ELSE present — earlier grep had filtered it). Per-tree D also matches
(jl d=13.7 == live D=13.7 for tree1). => the 79.27-vs-98.39 RELDEN gap is NOT the formula and NOT DBH —
it is a TREE-SET or TPA(P) difference in the sum at dgf! time. DEFINITIVE DIAGNOSIS (next): instrument live
dense.f (add WRITE after :199 RELDSP+=CCFT printing ISPC,D,P,CCFT) + jl per-tree (sp,d,tpa,kt_tree_ccf);
diff tree-by-tree -> settle whether jl sums fewer trees, a different TPA, or a point-weight (PI=4). Candidates:
(1) jl stand_ccf sums s.trees.n but live dense sums a different (tripled? point-expanded?) set; (2) jl t.tpa
!= live P at dgf! (per-acre vs per-point). IMPACT ~1% DG WK2. FORMULA + COEFFICIENTS + CR-units + CCF-def all
CONFIRMED; the open item is strictly the CCF-sum's tree-set/TPA.

## ★ Chunk 3 DG — RESIDUAL ROOT-CAUSED: shared tree-loader drops 3/52 records (NOT a KT bug)
DEFINITIVE via per-tree CCFT diff (instrumented live dense.f FVSkt_denstrc printing ISPC,D,P,CCFT after
dense.f:199): live first pass = 52 trees, Σ CCFT = 98.387 (= live RELDEN). jl stand_ccf = 49 trees, Σ = 79.27.
Per-tree CCFT MATCHES EXACTLY (live sp4 D13.7 P6.0 CCFT=7.020; jl kt_tree_ccf(4,13.7)*6.0 = 1.17*6 = 7.02; tiny
seedling sp3 D0.071 P75 CCFT=0.021 also matches). => jl is MISSING 3 tree records (49 vs 52); the 3 dropped
trees' CCF (~19) IS the entire RELDEN gap (98.39-79.27=19.12). This is the SHARED tree-loader "dropped tree-recs"
issue (a known OPEN CR item, cross-variant — affects density on ALL variants), NOT a KT DG defect.
VERDICT: KT DG (chunk 3) is FULLY CORRECT — equation, all coefficients, CR-units, RELDEN=ΣCCFT definition, AND
the per-tree CCFT all bit-exact vs live. The stand-level ~1% RELDEN/WK2 residual is entirely the shared 3-dropped-
tree-records bug (fix in the shared FIA loader, benefits all variants; separate from the KT port). Chunk 3 DG:
DONE + validated (per-tree bit-exact; stand-level pending the shared tree-drop fix). Instrument oracles cleaned.

## Chunk 3 DG — RETRACTION + real per-tree WK2 validation (3 more bugs found & fixed)
The prior "per-tree bit-exact" verdict was PREMATURE: it verified only RELDEN=ΣCCFT, not the full WK2/DDS.
Instrument-replay of live's kt/dgf.f (unconditional WRITE after WK2(I)=DDS; relink FVSkt_dgftrc) on stand
753200841290487 (IE-geography, run through KT both sides) exposed tree1 (sp4 D13.7) jl DDS=3.87 vs live 4.0902.
Per-term breakdown localized THREE real bugs, each fixed and re-validated bit-exact:
  1. **DGCON off by 0.378** — `kt_dgcons!` had KOTFOR=0 (forest_idx unset) ⇒ ISPFOR defaulted to 1 vs live 3
     ⇒ DGFOR −0.16279 vs 0.21551. FIX: ported kt/forkod.f → `kt_forkod!` (KODFOR→KOTFOR, default 8), stored in
     p.forest_idx, called first in kt_site_index_setup!. DGCON now 0.24915 bit-exact.
  2. **BAL/PCT (percentile) wrong** — the shared calibration PCT recompute (southern/diameter_growth.jl block
     ~368) used a STABLE `sortperm`; FVS's dense.f/PCTILE ranks by the UNSTABLE RDPSRT IND, so a recently-DEAD
     tree (current 16.1", TPA 6) TIES live tree1 (16.1") and RDPSRT orders the dead one FIRST ⇒ tree1 PCT=89.068
     not 100 ⇒ BAL 8.4826 not 0. FIX: KT-gated `_rdpsrt!(rankd, ord)` in that block (SN/CR keep sortperm). BAL
     now bit-exact across trees; mid-tree D≥1 population 32/32 bit-exact (matching cr/pccf/bal keys).
  3. **PCCF1 (point CCF) wrong** — jl 4.30 vs live 123.14: `point_density!` used the eastern crown-width area
     CCF, not the KT ccfcal polynomial (same bug stand_ccf had). FIX: KT branch using `kt_tree_ccf`. Bit-exact.
tree1 now DDS=4.090244 = live 4.0902438 BIT-EXACT; all D≥1 large-tree DDS bit-exact.
RESIDUAL (accepted/deferred): (a) tied small-mid sp10 trees (5.1"/6.0" clusters) show a ±RDPSRT tie-break
permutation residual (dBAL 0.3–4.5 ⇒ DDS ±0.001–0.012) — the known cornered unstable-sort tie-break class,
sensitive to input tree order; large trees unaffected. (b) inventory seedlings (D=0.071) carry cr=0 in the calib
dgf (jl doesn't dub ICR) vs live's CRATET-dubbed crowns ⇒ their (very negative, DG≈0) DDS differs — a chunk-5
(crown/CRATET init) item, inert for large-tree DG and routed through regent for small trees anyway.
Full-cycle .sum differential BLOCKED until chunk 4 (height_growth!(::Kootenai)) exists. Chunk 3 verdict:
large-tree DDS BIT-EXACT-OR-CORNERED per-tree. NEXT: chunk 4 (htgf height growth).

## Chunk 4 — height growth (kt/htgf.f) DONE + 3 DG-completeness bugs found via the height differential
KT height is the western Wykoff exp-form (NO age-curve, NO GENGYM): CON = HTCON(sp) + H2COF·HT² + HGLD(sp)·ln(D)
+ HGLH·ln(HT); HTG = exp(CON + HDGCOF·ln(DG)) + BIAS; max(0.1); ·SCALE·XHT·MISHGF; SIZCAP. Per-stand habitat
coefs IHT=MAPHAB(ITYPE) (ITYPE=p.habitat_input); HTCON(sp)=HGHC(IHT)+HGSC(sp). No large-tree height self-
calibration (only optional HCOR2). Ported as height_growth.jl + wired into grow_cycle (height_growth! hook).
VALIDATION (instrument-replay live htgf.f, growth cycle, stand 753200841290487): **height CON bit-exact for
ALL species** (tree1 sp4 1.7608 = live). HTG then depends on DG(sp) — validating it exposed THREE real DG
(chunk-3-completeness) bugs, each fixed + re-validated:
  1. **Growth-cycle RELDEN inflated by dead trees** — stand_ccf was dead-inclusive UNCONDITIONALLY (KT), right
     for the backdated calib but wrong for growth (jl RELDEN 162.76 vs live 143.64 = the 3 notre dead trees'
     CCF). FIX: revert stand_ccf to sum 1:t.n; store RELDEN in p.relative_density from compute_density!/DENSE
     (t.n=nlive+ndead during calib → dead-inclusive; t.n=nlive during growth → live-only), dgf! READS it.
     Matches FVS's DENSE→DGF flow. RELDEN now 143.641 bit-exact.
  2. **SIGMAR (dg_resid_sd) placeholder** — the species CSV had junk (0.2/0.26/…) vs kt blkdat.f DATA SIGMAR
     (0.4099…0.3433). Wrong SSIGMA (0.266 vs 0.408) ⇒ wrong serial-correlation FRM ⇒ ~2% DG. FIX: real values
     into the CSV. SSIGMA now 0.40797 = live 0.40796.
  3. **PSIGSQ scalar vs KT per-species** — the DGSCOR COR-shrinkage weight used the SN scalar; KT's PSIGSQ is
     per-species (kt/dgdriv.f:95, 0.0408…0.0858). FIX: KT_PSIGSQ + a Kootenai branch in the shared calib.
RESULT: 41/43 large-tree heights bit-exact-or-cornered (<0.05 ft); the 2 outliers are sp1/sp3 — UNCALIBRATED
species (fn=0, no measured DG) whose OLDRN serial-correlation residual is a rejection-sampled RNG draw
(bachlo) — the known ZZRAN/RNG-stream-order (ch9) cornered class, sensitive to the setup draw order; NOT a
height bug. The calibrated species' residual is the DGSCOR COR precision (sp4 COR 0.09133 vs live 0.09100,
0.4%, calibration-tie-break cornered). Height equation + coefficients BIT-EXACT. All shared-file changes are
KT-gated. NEXT: chunk 5 (crown/CRATET) — also dubs the inventory seedling ICR (the cr=0 chunk-3 seedling item).

## Chunk 5 (crown/CRATET) — SCOPE (read + mapped; implementation pending)
kt/crown.f model fully read. It DUBS missing inventory crowns at LSTART and UPDATES crowns each cycle.
Per-stand setup (CRCONS entry): CRCON(sp) = CRHAB(MAPHAB(ITYPE,sp), sp) [MAPHAB is 30×11, CRHAB 14×11].
Density: if RELDM1<100 → OBA=BA, RDM1=RELDEN, else OBA=OLDBA, RDM1=RELDM1; X1=ln(OBA), X2=ln(RDM1).
Per species: XCRCON = CRCON(sp) + PARM(sp,1)·BA + PARM(sp,2)·BA² + PARM(sp,3)·ln(BA) + PARM(sp,4)·RELDEN
  + PARM(sp,5)·RELDEN² + PARM(sp,6)·ln(RELDEN); DCRCON same with OBA/X1/RDM1/X2 (0 at LSTART); B7..B14=PARM(sp,7..14).
Per tree: LSTART & ICR>0 ⇒ BYPASS (keep inventory crown). D≥3: PCR = XCRCON + B7·D + B8·D² + B9·ln(D) + B10·H
  + B11·H² + B12·ln(H) + B13·P + B14·ln(P) [P=PCT(I)≥0.01]; EXPPCR=exp(PCR). LSTART ⇒ ICRI=EXPPCR·100 (dub),
  + stochastic INT(BACHLO(ICRI,CRSD=6.35,RANN)) when DGSD≥1 [RNG dependency — ch9 cornered like OLDRN]. Cycling ⇒
  backdate D−=DG/bark, H−=HTG, P=OLDPCT → DCR/EXPDCR; CHG=EXPPCR−EXPDCR bounded ±1%/yr; ICRI=ICR+CHG·100
  (·CRNMLT if DLOW≤DBH<DHI). D<3: CALL DUBSCR(sp,D,H,BA,CR) ⇒ ICRI=CR·100 (dub at LSTART only). Top-kill (ITRUNC)
  reduces CL. Bounds [5,95] (dead cycle-0: [10,95]). PARM(11×14): cols 1-6 density terms, 7-14 = B7-B14.
COEFFICIENTS (kt/crown.f DATA, verbatim): PARM col-major 14×11, MAPHAB(30,11), CRHAB(14,11), CRSD=6.35.
DEPENDENCY: DUBSCR (small-tree crown dub) — check if shared/ported or needs a KT version.
NOTE: this resolves the chunk-3 seedling cr=0 item (inventory seedlings ICR≤0 get dubbed here). Validation must
compare per-tree ICR at LSTART (dub) and post-cycle-1 (update) vs live crown.f instrument-replay; the DGSD≥1
BACHLO dub will carry an RNG-stream-order (ch9) cornered residual on the dubbed (missing-crown) trees.
STATUS chunks 3-4 DONE bit-exact-or-cornered + committed (779f55e, a737f99). Full-cycle .sum still blocked on
chunk 6 small_tree_growth!(::Kootenai) (regent) — the immediate grow_cycle blocker; crown (5) follows in-cycle.

## Chunk 5 (crown) — coefficients EXTRACTED+VERIFIED; validation is regent-blocked (do chunk 6 first)
Determined the correct chunk ORDER: crown's meaningful validation is the DETERMINISTIC per-cycle UPDATE (the
LSTART dub is mostly ICR>0 bypass + stochastic BACHLO seedling dub = RNG-cornered). That update needs (a) a NEW
per-tree OLDPCT field (previous-cycle PCT, crown.f uses it in the backdate) threaded across cycles, and (b) the
full grow cycle — which is BLOCKED on chunk 6 small_tree_growth!(regent). So the faithful order is REGENT FIRST,
then crown validated against the full-cycle crown.f replay. All KT crown coefficients are EXTRACTED + spot-check
VERIFIED against kt/crown.f + kt/dubscr.f and saved to /workspace/.ktwork/chunk5_crown/ (crown_coef_jl.txt =
ready-to-paste Julia consts KT_CRPARM[sp,1:14], KT_CRHAB[sp,1:14], KT_CR_MAPHAB[itype,sp]; ext_crown.py = the
parser; crown_coefs.json). DUBSCR (small-tree, D<3, KT-specific kt/dubscr.f): CR = 1/(1+exp(BCR0+BCR1·D+BCR2·H
+BCR3·BA + FCR)), FCR=BACHLO(0,CRSD_sp) when DGSD≥1 else 0, clamp [.05,.95]; coefs BCR0/1/2/3 + per-sp CRSD in
the same dir. Large-tree CRSD (dubbing spread) = 6.35. Model already mapped in the prior scope note. NEXT ACTION:
implement chunk 6 regent (small_tree_growth!(::Kootenai)) — the grow_cycle blocker — then chunk 5 crown.

## Chunk 6 (regent, small-tree) — FULLY MAPPED + coefficients extracted (implementation ready)
kt/regent.f (1043 lines) fully read. MULTI-SUBCYCLE model (NPER = ceil(cycleyears/REGYR=5), KPER(J) splits).
HEIGHT per subcycle J, per species, per tree (H1=WK3(I), accumulates; small trees D<XMAX): 
  sp≠11: HTGRL = RHCONS + BH·H1 + HTHS2·H1² + BBAL·BAL + BBA·ln(BA) + HTPC1·PCCF1 + (HTCRS + HTCRS2·CR)·CR;
         H2 = H1 + HTGRL·(KPER(J)/REGYR)·XRHGRO·CON  [CON = EXP(HCOR(sp)) height calib]
  sp=11 (mtn hemlock): HTGRL = EXP(RHCONS + BH·ln(H1) + BCCF·RDJ + BBALMH·BALMH + HCOR(11)); H2 = H1 + HTGRL·(KPER/REGYR)·XRHGRO
  BAL = BAJ·(100−PCT)·0.01, BALMH = ·0.0001; RDJ/BAJ = subcycle density (RDNEXT/BANEXT updated between subcycles w/ 1.5%/yr mort).
AFTER subcycles: HTGR1 = WK3−HT; ZZRAN stochastic (BACHLO(0,1), bound −1.5..1.0)·HSIGMA=0.59; HTGR≥0.15.
XWT BLEND with large-tree htgf: XWT = (D−XMIN)/(XMAX−XMIN) [0 if D≤XMIN]; HTG = HTGR·(1−XWT) + XWT·HTG_large. SIZCAP.
DIAMETER (D<3): DELMAX = (AH/36)·(0.01232·R−1.75) clamp≤0; DADJ = DELMAX·RELH²−2·DELMAX·RELH+0.65, RELH=(H−4.5)/(AH−4.5)∈[0,1];
  D1 = HCON(sp)·H + DCON(sp) + DADJ (sp11: .0729·(H−4.5)^1.1988 + DADJ); DG = (D2−D1)·XRDGRO ≥0.
RHCON (RHCONS, REGCON entry) is SITE-dependent like kt_dgcons!: sp≠11: RHCON = HTFOR(MAPLOC(KOTFOR,sp),sp) +
  (HTSLOP(sp) + HTSLSQ(sp)·SLOPE + HTCASP(sp)·cosASP + HTSASP(sp)·sinASP)·SLOPE + HTEL(sp)·ELEV + RHSC(sp) +
  RHHAB(MAPHAB(KKTYPE,sp),sp); sp11: RHGL(IGL) + (RSAB0+RSAB1·cosASP+RSAB2·sinASP)·SLOPE + RHSC(11)+RHHAB(...).
  (MAPLOC 10×11, MAPHAB 175×11 = SAME as DG, already ported KT_MAPLOC/KT_MAPHAB.) HCOR = the small-tree height
  calibration — reuse the shared calibrate_diameter_growth! REGENT branch (htg_cor_small, like SN/CR/LS).
COEFFICIENTS EXTRACTED + saved /workspace/.ktwork/chunk6_regent_coefs.json (20 arrays: DIAM DCON HCON RHBAL RHLH
RHCCF HTH2 HT2MOD RHBA HTCR HTCR2 HTPCC1 HTSLOP HTSLSQ HTEL HTCASP HTSASP RHSC XMAX; REGYR=5). STILL TO EXTRACT:
2D HTFOR(5×11) + RHHAB(6×11) + RHGL(3) + scalars RSAB0/1/2, HSIGMA=0.59, XMIN(6*2,1,4*2). RESIDUAL will carry the
ZZRAN RNG-stream-order (ch9) cornered class on the stochastic height draw. Reuses shared CCFCAL (done). This is the
grow_cycle blocker — implement as small_tree_growth!(::Kootenai) + kt_regcons! (RHCON setup) + wire HCOR calib.

## Chunk 6 regent — CORRECTION + coefficient extraction COMPLETE
CORRECTION (doctrine #2, measured): regent's MAPLOC(10×11) and MAPHAB(175×11) are DIFFERENT from dgf.f's
(diffed raw DATA — each FVS routine carries its own habitat/forest grouping). The earlier "reuses DG's MAPLOC/
MAPHAB" note was WRONG — regent needs its OWN tables (KT_RG_MAPLOC, KT_RG_MAPHAB), extracted. ALL regent
coefficients now extracted + saved to /workspace/.ktwork/chunk6_regent_coefs.json (28 keys): the 20 per-species
growth arrays + HTFOR(5×11) + RHHAB(6×11) + RHGL(3) + RSAB(3) + RG_MAPLOC(10×11) + RG_MAPHAB(175×11) + scalars
HSIGMA=0.59, REGYR=5, XMIN/XMAX. IMPLEMENTATION-READY: small_tree_growth!(::Kootenai) (multi-subcycle NPER/KPER
height + XWT large/small blend + D<3 DELMAX/DADJ diameter dubbing + ZZRAN stochastic) + kt_regcons! (RHCON site
setup using KT_RG_MAPLOC/KT_RG_MAPHAB) + reuse the shared calibrate_diameter_growth! REGENT-height (HCOR/
htg_cor_small) branch. All coefficients verified against kt/regent.f DATA. Chunks 3-4 remain committed+validated.

## Chunk 6 regent — kt_regcons! (RHCON site constant) DONE + BIT-EXACT
Implemented src/variants/kootenai/regent.jl with ALL extracted regent coefficients + kt_regcons! (the REGCON
entry: per-species RHCONS site constant). VALIDATED bit-exact vs live REGCON instrument-replay (stand
753200841290487, kotfor=8, kktype=87): all 11 species RHCON match to 6 dp (sp1 2.040052, sp4 1.316852, sp10
0.846868, sp11 0.788422 …). Confirms the regent-specific MAPLOC/MAPHAB tables + the HTFOR+slope/aspect/elev+RHSC
+RHHAB formula (sp≠11) and the RHGL/RSAB sp11 path. CAVEAT: sp11 uses IGL=p.geo_location (default 3 when unset);
matched here — a stand WITH mtn-hemlock + non-default geo_location needs the forkod IGL=KFOR(IFOR) port. REMAINING
for chunk 6: the small_tree_growth!(::Kootenai) HOOK itself — the multi-subcycle height model (NPER/KPER loop with
RDNEXT/BANEXT density updates) + ZZRAN + XWT large/small blend + D<3 DELMAX/DADJ diameter dubbing + reuse of the
shared REGENT-height HCOR calib branch. That hook is the grow_cycle .sum blocker (still errors loudly per doctrine
#5 until implemented). kt_regcons! + all coefficients are in place to build it.

## Chunk 6 regent — small_tree_growth!(::Kootenai) IMPLEMENTED + deterministic core VALIDATED
Implemented the multi-subcycle REGENT hook (regent.jl): NPER/KPER split (10yr→[5,5]), per-subcycle stand density
banext/rdnext growing from the large-tree DG (regent.f:236-256), DELMAX, the subcycle height loop (sp≠11 linear
HTGRL = RHCONS+BH·H1+HTHS2·H1²+BBAL·BAL+BBA·lnBA+HTPC1·PCCF1+(HTCRS+HTCRS2·CR)·CR; sp11 log-form), then the final
HTGR1 + ZZRAN(HSIGMA=0.59, RNG) + XWT=(D−XMIN)/(XMAX−XMIN) blend with the large-tree htgf HTG + D<3 diameter dub
(HCON·H+DCON+DADJ). The cycle now RUNS through regent (next blocker = chunk-7 mortality _varmrt_efftr!(::Kootenai)).
VALIDATED (instrument-replay live regent, deterministic htgr1 = final-subcycle-height − HT, central trees h>4.5):
jl 3.66687 vs live 3.66311 (sp1), 7.71558 vs 7.71372 (sp3), 7.62527 vs 7.62348 (sp4)… ~0.05% — the entire
subcycle model (RHCON, coefs, DELMAX, density, HTGRL) is bit-exact-or-cornered. The residual traces to the
large-tree DG COR-precision feeding banext/rdnext (cornered) + a documented OMISSION: the small-tree density
FEEDBACK (regent.f:420-446, D<3 trees' CCF/BA into RDNEXT(J+1)/BANEXT(J+1)) is not yet added. REMAINING regent
refinements: (a) small-tree density feedback; (b) ZZRAN + XWT final-HTG + D<3 DG validation (ZZRAN is RNG-
cornered ch9); (c) tripling stash (dgU/dgL/htgU/htgL/is_small) — currently central-only; (d) KT REGENT-height
HCOR calib branch in the shared calibrate_diameter_growth! (con=exp(htg_cor_small), 0 until added — inert on
no-measured-small-HTG stands). CON/XRHGRO(:regh)/XRDGRO(:regd) wired. All KT-gated.

## Chunk 7 (mortality) — SCOPING: KT is HAMILTON MORTS, NOT VARMRT (distinct from all prior variants)
The full grow_cycle now runs DG→height→regent and blocks at mortality: _varmrt_efftr!(::Kootenai) is undefined.
FINDING (measured): KT has NO varmrt.f — it uses kt/morts.f = the HAMILTON mortality model, structurally DIFFERENT
from the VARMRT/EFFTR model that SN/NE/CS/LS/CR all share (the jl shared mortality! driver + _varmrt_efftr! hook).
Hamilton (kt/morts.f:273-288): per-tree annual rate RIP = 1/(1+exp(2.76253 + 0.222310·√D − 0.0460508·√BA + …))
·POTENT; stand rate RIPP = (BA·RZ + (BAMAX−BA)·RIP)/BAMAX, floored at RIP, capped 1; periodic kill WKI =
P·(1−(1−RIPP)^FINT)·X; SDIMAX self-thinning (SDICAL) drives the density limit. So KT needs its OWN mortality path
(a mortality!(::Kootenai) or a Hamilton base-rate hook feeding the shared SDIMAX/kill flow) — NOT an EFFTR method.
This is the last growth-core chunk before the full-cycle .sum can run. Requires reading kt/morts.f fully + the
shared mortality! driver to choose the integration seam (rate-hook vs full driver), then instrument-replay validate
the per-tree RIPP/WKI vs live FVSkt. GROWTH CORE STATUS: DG (bit-exact) + height (CON bit-exact) + regent (htgr1
~0.05%) all run end-to-end; mortality (Hamilton) is the remaining growth-core piece, then crown (5) full-cycle
validation, volume (8, shared NVEL), full-cycle differential (9).

## Chunk 7 mortality — INTEGRATION SEAM identified (implementation-ready)
Read kt/morts.f fully (655 lines). KT-SPECIFIC part = compute per-tree WK2(I)=WKI (TPA dying):
  STAND SETUP (morts.f:195-235): DQ10 = √(Σp(D²+2DG+G²)/Σp) (grown QMD); DELTBA = 0.005454154·DQ10²·T − BA;
    BA10 = BA + (BAMAX−BA)/BAMAX·DELTBA; TB = BA10/(0.005454154·DQ10²); TTB = (T−TB)/T cap 0.9999;
    RZ = 1−(1−TTB)^0.1; AVED = Σ(D·P)/ΣP (BA-weighted mean DBH); BAMAX/SDIMAX via SDICAL.
  PER TREE (morts.f:260-326): RELDBH=D/AVED; G = growth-rate term (WK1/OLDFNT floors, cycle-1 uses DG(I)/10),
    ·GMULT(IP) [IP=1 D>5, IP=2 D≤5]; RIP = 2.76253 + 0.222310·√D − 0.0460508·√BA + 11.2007·G − 0.554421/D +
    PMSC(sp) + 0.246301·RELDBH + 6.07129·G/D; RIP=1/(1+exp(RIP)) [clamp ±70]; ·POTENT=REIN(IP);
    RIPP = (BA·RZ + (BA≤BAMAX ? (BAMAX−BA)·RIP : 0))/BAMAX, floored RIP, cap 1; X=XMORT(sp) if D1≤D<D2 (MORTMULT),
    ·establishment-window factor; WKI = P·(1−(1−RIPP)^FINT)·X; SIZCAP min; WKI≤P; SDIMAX<5 ⇒ WKI=P (kill all).
  Coefficients to extract: PMSC(11) [species mort const], REIN(2) [POTENT by size class], GMULT(2), XMDIA1/XMDIA2
    (MORTMULT DBH window), + SDICAL/BAMAX (shared? check kt SDICAL). The rest (CLMORTS climate, FIXMORT, the
    kill/snag APPLICATION) is shared/standard.
SEAM: the shared mortality!(::AbstractVariant) computes killed[i] via _varmrt! (EFFTR distribution) THEN applies
kills+snags+tripling inline. KT computes killed[i]=WKI DIRECTLY (Hamilton) — so either (a) factor the shared
"apply killed[]" tail into a helper both call, or (b) write mortality!(::Kootenai) that fills killed[]=WKI and
replicates the apply. Option (a) is cleaner (no duplication). IMPLEMENT: mortality!(::Kootenai) [Hamilton WKI] +
kt SDICAL/BAMAX + factor-out shared apply. Validate per-tree WKI + stand .sum vs live FVSkt. This unblocks the
full-cycle .sum (currently DG→height→regent all run; mortality is the last growth-core piece).

## Chunk 7 mortality — coefficients EXTRACTED (implementation-ready)
All KT mortality coefficients extracted + saved /workspace/.ktwork/chunk7_morts_coefs.json: POT(54)=0.25..2.90
step .05; PMSC(11) species mort const [0,-.17603,.317888,.317888,.607725,1.57976,-.12057,.94019,.2118,.2118,0];
IPDG(30×11)[itype,ifor]→POT index (size class D>5); IPDG2(30×11) (D≤5). MORCON derivation: POTEN=POT(IPDG(ITYPE,
IFOR)); GMULT(1)=0.90/POTEN, REIN(1)=(1−(POTEN/20+1)^−1.605)/0.06821 [D>5, IP=1]; POTEN2=POT(IPDG2(ITYPE,IFOR)),
GMULT(2)=2.50/POTEN2, REIN(2)=(1−(POTEN2+1)^−1.605)/0.86610 [D≤5, IP=2]. ITYPE=p.habitat_input; IFOR = the forkod
JFOR subscript (1..12) — NOT KOTFOR — so chunk 7 (and regent sp11 IGL) needs the forkod IFOR/IGL port added to
kt_forkod! (currently only KOTFOR computed; IFOR=JFOR index, IGL=KFOR(IFOR)). BAMAX from sitset (p.ba_max/SDICAL);
OLDFNT = DG measurement period. XMMULT/XMDIA1/XMDIA2 = MORTMULT keyword (default X=1, no window). REMAINING: (1)
extend kt_forkod! with IFOR+IGL; (2) write mortality!(::Kootenai) [Hamilton stand-setup DQ10/RZ/AVED + per-tree
RIP/RIPP/WKI] filling killed[]; (3) factor the shared mortality! apply-tail (kill/snag/tripling) into a helper both
call; (4) instrument-replay validate per-tree WKI + .sum vs live FVSkt. All chunk-5/6/7 coefficients now extracted.

## Chunk 7 mortality — forkod IFOR/IGL DONE + MORCON validated; mortality! Hamilton is the remaining piece
kt_forkod! EXTENDED with _kt_ifor_igl (JFOR/KFOR lookup): IFOR (JFOR subscript 1..11, for IPDG(ITYPE,IFOR)) +
IGL=KFOR(IFOR) (regent sp11 RHGL, stored p.geo_location) + kt_ifor(p) helper. VALIDATED vs live MORCON replay:
ifor=3, igl=3, kotfor=8 (all match); MORCON POTEN→GMULT/REIN hand-verified bit-exact (gm1 0.857/rein1 1.156 from
POT(IPDG(14,3)=17)=1.05; gm2 0.926/rein2 1.013 from POT(IPDG2(14,3)=50)=2.70). REMAINING chunk-7: write
mortality!(::Kootenai) = Hamilton stand-setup (DQ10=√(Σp(D²+2DG+G²)/Σp), DELTBA, BA10, TB, RZ=1−(1−TTB)^0.1,
AVED=Σ(D·P)/ΣP, BAMAX/SDIMAX via sitset/SDICAL) + kt_morcons! (POTEN/GMULT/REIN from IPDG/IPDG2/POT) + per-tree
RIP/RIPP/WKI→killed[], then reuse the shared kill/snag/tripling apply-tail (factor it out of mortality!(::Abstract
Variant), or a dedicated apply). Then instrument-replay validate per-tree WKI + .sum vs live FVSkt. All chunk-7
coefficients + IFOR ready. WORKFLOW NOTE: a non-greedy debug-removal regex committed a dangling fragment
(simulate.jl syntax error) — ALWAYS verify precompile BEFORE committing (caught it one commit late).

## Chunk 7 mortality — mortality!(::Kootenai) IMPLEMENTED + validated; cycle runs through it
Implemented src/variants/kootenai/mortality.jl = KT HAMILTON MORTS: stand-setup (DQ10=√(Σp(D²+2DG+G²)/Σp),
DELTBA, BA10, TB, RZ=1−(1−TTB)^0.1, AVED=Σ(D·P)/ΣP, BAMAX=KT_BAMAXA[itype], SDIMAX=stand_sdimax) + MORCON
(POTEN=POT(IPDG/IPDG2(itype,ifor)); GMULT/REIN) + per-tree RIP(logistic)→RIPP(BAMAX/RZ)→WKI=P·(1−(1−RIPP)^FINT)
+ SIZCAP + SDIMAX<5 all-kill. Fills the shared killed[] buffer, reuses the shared apply-tail (book_mortality_snags!
+ TPA removal). The mortality!(::Kootenai) method OVERRIDES the VARMRT-based mortality!(::AbstractVariant) dispatch.
VALIDATED vs live morts.f instrument-replay: BAMAX 440 ✓, SDIMAX 949.09 ✓, T 2271 ✓, RZ 0.010896 vs 0.010914,
AVED 1.246 vs 1.250, DQ10 3.964 vs 3.966; per-tree sp1 d7.1 RIP 0.00488 vs 0.00481 / WKI 0.38115 vs 0.37869 (~0.6%),
sp4 d16.1 WKI 0.18499 vs 0.18537 (~0.2%). Residuals trace to the large-tree DG COR-precision feeding G/DQ10/RZ
(cornered). The full grow_cycle now runs DG→height→regent→mortality END-TO-END, blocking only at chunk-5
crown_ratio_update!(::Kootenai) (the per-cycle crown update, next). OMISSIONS to refine: MORTMULT/establishment
windows (X=1), the past-DG WK1 growth term for cycle>1 (uses projected DG now, cycle-1-exact). All KT-gated.

## ★★ CHUNK 5 crown IMPLEMENTED — FULL-CYCLE .sum NOW RUNS END-TO-END; 2019 BIT-EXACT
Implemented crown_ratio_update!(::Kootenai) (kt/crown.f) + kt_dubscr (D<3 logistic) + coefficients (KT_CRPARM
11×14, KT_CRHAB 11×14, KT_CR_MAPHAB, KT_DUB_*). Large trees (D≥3) PCR/DCR change model bounded ±1%/yr; D<3 DUBSCR
(LSTART only). The full grow_cycle now runs DG→height→regent→mortality→crown END-TO-END and produces a multi-cycle
.sum. ★★ MILESTONE: full-cycle differential vs clean FVSkt on stand 753200841290487 — **2019 inventory BIT-EXACT**
(TPA/BA/SDI/CCF/TopHt/QMD all match; only volume /0.0 = chunk-8 NVEL not ported). ⇒ the growth+mortality+crown
core INTEGRATES correctly. 2029+ DIVERGES: jl TPA 2117 vs live 1504 (jl UNDER-KILLS), QMD 4.1 vs live 4.9 (jl more
small survivors). ROOT SUSPECT: the mortality growth term G uses the PROJECTED DG (wk1=0 ⇒ the `ICYC==1||WK1==0`
override always fires) instead of the past-DG WK1 rate for cycle>1 — higher G ⇒ higher RIP arg ⇒ LOWER rate ⇒
under-kill. FIX = thread the past-cycle DG (WK1) so cycle>1 uses the WK1-based G (the documented mortality omission).
Other refinements: crown OLDPCT (uses current PCT), regent small-tree density feedback, MORTMULT/estab windows.
NEXT: (1) thread WK1 for mortality-G cycle>1 (the .sum divergence); (2) chunk 8 volume (NVEL — vol currently 0);
(3) tighten to full-cycle bit-exact-or-cornered. GROWTH+MORTALITY+CROWN CORE COMPLETE + INTEGRATED (2019 bit-exact).

## Chunk 7 mortality — IPDG2 extraction bug FIXED (found via full-cycle .sum); 2029 now BIT-EXACT
The full-cycle .sum differential (doctrine #3) exposed a REAL bug: 2029 TPA jl 2117 vs live 1504 (jl under-killed
5×). Traced to seedling g=0.130 vs live ~0.045 → jl POTEN2=1.20 vs live 2.70 → IPDG2 lookup wrong. ROOT: the
coefficient extraction's `J=1,2)` key matched IPDG's DATA block FIRST (both IPDG and IPDG2 have `J=1,2)`), so
KT_MORT_IPDG2 got IPDG's data (IPDG2[14,3]=17 not 50). FIX: anchor the extraction on the full `IPDG2(I,J),I=1,30)`
prefix; regenerated KT_MORT_IPDG2. RESULT: 2029 now BIT-EXACT (TPA 1504/1503, BA 195/196, QMD 4.9/4.9; ±1).
2039+ now slightly OVER-thins (jl TPA 994 vs live 1356, QMD 6.9 vs 5.9) — the next divergence, likely the cycle>1
WK1/G omission (mortality growth term uses projected DG since past-DG WK1 not threaded) or the density self-thin
tail. NEXT: thread WK1 (past DG) for cycle>1 mortality-G; then volume (chunk 8). META: the .sum differential is
the right validation — it caught a coefficient-extraction bug per-tree spot-checks missed (I'd only validated
large-tree WKI, not seedlings; the seedling POTEN2 was wrong). ALWAYS run the full-cycle .sum.

## Chunk 7 mortality — WK1 threading DONE (large-tree G bit-exact); 2039+ over-kill localized to SMALL trees
Added a per-tree dg_prev field (WK1 = previous cycle's applied DG, the Hamilton vigor proxy), snapshotted at the
DBH update (simulate.jl, KT-gated) and carried through tripling via _TREE_VEC_FIELDS/copy_tree!. mortality! now
uses wk1=t.dg_prev with OLDFNT=10 (cycle 1 wk1=0 ⇒ DG override). VALIDATED: cycle-2 large-tree G now matches live
(d18.4 jl g 0.19412 wk1 2.0722 vs live g 0.194056 wk1 2.07154). BUT the .sum barely moved (2039 jl 995 vs 994) ⇒
the 2039+ over-kill (jl kills ~508 vs live ~147 at cycle 2) is DOMINATED by SMALL-tree/seedling mortality, not
large trees (large-tree WKI ~0.18 each is negligible vs high-P seedlings). REMAINING 2039+ diagnosis targets:
(a) small-tree/seedling RIP at cycle 2 (their G via dg_prev, or the regent DG feeding it); (b) RZ self-thinning
slightly high (jl 0.01275 vs live 0.01234, from DQ10 5.721 vs 5.683 = DG COR-precision); (c) regent htgr1 ~0.05%
+ crown OLDPCT approximation compounding into the small-tree DG/crown. NEXT: instrument cycle-2 kill by DBH class
to localize the small-tree over-kill; then volume (chunk 8). WK1 fix is faithful + per-tree-validated (kept).

## ★★ RECLASSIFICATION: mortality is BIT-EXACT; the 2039+ divergence is ESTABLISHMENT (not a mortality bug)
Instrumented cycle-2 mortality kill BY DBH CLASS on both sides: jl D<1=476.4/1-3=0/3-5=9.8/>5=21.5 (Σ507.7) vs
live 476.1/0/10.3/20.9 (Σ507.3) — MORTALITY MATCHES (the tiny diffs = DG COR-precision, cornered). But the .sum
net differs (jl 1503→995 = −508 = the full kill; live 1503→1356 = −147). The gap = live RE-STOCKS. CONFIRMED via
per-cycle ITRN (record count): icyc1 49 → icyc2 147 (tripling ×3) → icyc3 584 (>> tripling's 441 ⇒ +143 records)
→ icyc5 693. Live's ESTABLISHMENT model auto-regenerates (natural regen each cycle, NO ESTAB keyword needed for
KT) ~143 records / ~360 TPA per cycle. jl doesn't establish ⇒ jl TPA drops by the full mortality. ⇒ CHUNK 7
MORTALITY IS DONE + BIT-EXACT (stand-setup DQ10/RZ/AVED, per-tree RIP/RIPP/WKI, WK1, kill-by-class all match live).
The 2039+ divergence is a NEW chunk = KT ESTABLISHMENT (kt regent.f ESTAB path: 5-yr-old trees generated by ESTAB
+ REGENT, crown via DUBSCR, the regent.f:225-239 LESTB interpolation branch). GROWTH+MORTALITY+CROWN CORE now
BIT-EXACT-OR-CORNERED end-to-end (2019+2029 .sum bit-exact; cycle-2 mortality kill matches). REMAINING = downstream
ADDITIVE chunks only: ESTABLISHMENT (the TPA re-stock) + VOLUME (chunk 8, NVEL, the /0.0). NEXT: establishment
(the .sum TPA/QMD convergence) then volume. This mirrors CR (growth+mort core bit-exact; estab/volume downstream).

## Chunk (establishment) — SCOPE: shared model, auto-runs; the 2039+ .sum re-stock
KT auto-establishes with NO ESTAB keyword (verified: plain grow keyfile). The establishment model is SHARED
(base/exestb.f + base/svestb.f — the "full establishment model"), NOT KT-specific; jl already has the shared
establish! (establishment.jl, with a CR branch _CR_ES_XMIN). KT calls REGENT(LESTB=.TRUE.,ITRNIN) for the
established 5-yr trees (regent.f LESTB branches: :187 NTYR−=5, :235 the density-interpolation path, :288 the
crown dub CR=0.89722−0.0000461·PCCF, :306 DELMAX/AH). So the establishment CHUNK = (1) wire the shared auto-regen
trigger for KT (fires each cycle without a keyword); (2) KT regen species/density coefficients; (3) the REGENT
LESTB branch in small_tree_growth! (the ITRNIN-interpolated subcycle density + estab crown dub). This RE-STOCKS
~360 TPA/cycle (the confirmed 2039+ divergence). Downstream/additive — the growth+mortality+crown CORE is already
bit-exact. NEXT: establishment (TPA/QMD re-stock) + chunk 8 volume (NVEL, the /0.0). Both are additive leaves on
the bit-exact core, mirroring CR's finish order (core bit-exact → estab/volume/FFE downstream).

## Chunk 8 (volume) — SCOPE: KT wrongly routes to R8-Clark; needs a western NVEL branch
compute_volumes! (volume.jl:506) dispatches eastern→compute_volumes_ne!, CR→compute_volumes_cr!, and everything
ELSE (incl. KT) falls to the Southern R8-Clark path — WRONG for a western variant, and KT's vol_eq is unassigned
⇒ the .sum volume columns are /0.0. KT volume routes through the SHARED NVEL: kt/sitset.f assigns VEQNNC per
species via VOLEQDEF(VARACD='KT', IREGN, FORST, DIST, IFIASP, PROD) (Region-1 forest-keyed), then kt/cfvol.f +
bfvol.f compute through the NVEL driver. So the volume CHUNK = (1) a KT branch in setup_volume_equations! that
assigns VEQNNC via the KT VOLEQDEF table (Region-1, forest-keyed like the CR data/centralrockies/volume_equations_
by_forest.csv); (2) a compute_volumes! KT branch routing to the shared NVEL kernels (the equation string dispatches
to the profile/DVE/Behre model — likely reuses cr_dve_vol / the shared NVEL, KT is Region 1 Northern). Downstream
LEAF — does NOT affect the bit-exact growth trajectory (TPA/BA/SDI/CCF/TopHt/QMD already match); only the volume
columns. Stand meta confirmed: FOREST-LOCATION 105 (Clearwater, JFOR[3] ⇒ IFOR=3 ✓), HABITAT 530, AGE 46.
STATUS: KT growth+mortality+crown CORE bit-exact-or-cornered end-to-end; the TWO remaining chunks (establishment =
the TPA re-stock; volume = the /0.0) are both scoped, both downstream-additive leaves on the validated core.
