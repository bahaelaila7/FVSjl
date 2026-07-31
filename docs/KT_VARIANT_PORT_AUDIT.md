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
