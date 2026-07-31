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
