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
