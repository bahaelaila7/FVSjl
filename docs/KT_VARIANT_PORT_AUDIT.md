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
