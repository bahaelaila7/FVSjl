# NC (Klamath Mountains) variant port — chunk audit (doctrine #6)

Branch kt-variant-port. Oracle /workspace/.ncwork/FVSnc_clean. Canonical stand nct01
(tests/FVSnc/nct01.key + nct01.tre + nct01.sum.save). MAXSP=12, imperial, Zeide SDI, DGSD=2.0.
Coefficient extraction (ground-truth): docs/NC_CHUNK1_EXTRACTION.md.

## Chunk status
- **ch0 scaffold** ✅ (0d2437c): struct Klamath, variant_from_code("NC"), FVSjl.jl includes. Package loads, dispatches.
- **ch1 species** ✅ (42a3825 data + f83779b species.jl): species_coefficients.csv (41 cols, measured flat values +
  placeholders for structural cols) + species_translation.csv + ecocls_sdimax.csv (90-row PA→SDIMX). species.jl
  block-data init (Zeide, 10-yr, dg_sd=dg_stddev_bound=2.0, seed 55329). coefficients(Klamath()) LOADS correct.
- **ch2 site_index** ◐ IMPLEMENTED (63c0ead), NOT yet validated: site_setup! = nc_forkod (505/510/514/611/705/800/
  712→IFOR1-7) + nc_sitset (SITEAR site-range interp + provisional SDImax). Pipeline runs through init + into
  grow_cycle!, advancing to the ch3 gate dgf!(::Klamath).
  ⚠ **cyc0 VALIDATION BLOCKER (found)**: jl cyc0 summary is FAR off — jl TPA 28 / BA 7 / SDI 14 / TopHt 32 vs live
  536 / 77 / 160 / 63. ROOT = TREE-INPUT EXPANSION, not a growth chunk: jl reads 27 trees (file has 30) at tpa≈1.0
  each (sum 31, gross_space 1.1 ⇒ 28/ac) but live expands each ~18× to 536/ac. nct01.key `DESIGN 11.0 1.0` (11-plot
  design) + the tree-record count/prob field are not applying the plot expansion. NEXT: trace DEFAULT_TREE_FORMAT
  vs nct01.tre column layout (`011SP 11510 0734 00111`) + the DESIGN 11-plot expansion in the shared tree reader —
  is jl dropping the 3 records + not multiplying by the per-plot expansion? Likely a shared reader/DESIGN path NC's
  stand exercises differently (check how CI/other stands' DESIGN + tree count expand). This gates ALL cyc0 .sum
  validation — resolve before ch3.
- **ch3 DG** ▶ NEXT gate: dgf!(::Klamath). NC-specific (ISCT sections + DGHAH + 2 sets + sp12 redwood + nc_bratio
  POWER). Coefficients measured (extraction doc).
- ch4 height · ch5 crown · ch6 regent · ch7 mortality (reuses EM/UT) · ch8 volume (vollib): pending.

## Validation gate
cyc0 nct01 (live 1990): TPA 536 BA 77 SDI 160 CCF 87 TopHt 63 QMD 5.1 TCuFt 1308 MCuFt 449 (forest 371... wait
STDINFO forest=505 Klamath; .sum ForType col=371). Validate TPA/BA/SDI/TopHt/QMD FIRST (no growth models needed) —
currently blocked by the tree-expansion bug above.
