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
  ✅ **cyc0 "blocker" was a HARNESS ARTIFACT, NOT a bug** (corrected): the jl TPA 28 came from measuring
  `each_stand`'s PRE-EXPANSION state — the documented "each_stand returns pre-calibration state, measure in the REAL
  run only" trap (CI #142 hard lesson). PROOF: CI cit01 (same 248112 stand) gives the IDENTICAL pre-expansion state
  via each_stand (n=27, gross_space 1.1, TPA 28, tpa[1:3]=[1,3,1]) yet its REAL-run cyc0 is bit-exact 536. So NC's
  tree reading is fine (identical to CI); the PROB×plot-design expansion to 536/ac happens in the real run, after
  each_stand's yield. ⇒ cyc0 must be validated from run_keyfile's .sum, NOT each_stand. That currently needs the ch3
  dgf! gate to not crash (a no-op stub suffices for the cyc0 row, which is pre-growth). NEXT: stub/implement dgf!
  → get run_keyfile's 1990 .sum row → validate TPA/BA/SDI/TopHt/QMD vs live 536/77/160/63/5.1.
- **ch3 DG** ▶ NEXT gate: dgf!(::Klamath). NC-specific (ISCT sections + DGHAH + 2 sets + sp12 redwood + nc_bratio
  POWER). Coefficients measured (extraction doc).
- ch4 height · ch5 crown · ch6 regent · ch7 mortality (reuses EM/UT) · ch8 volume (vollib): pending.

## Validation gate
cyc0 nct01 (live 1990): TPA 536 BA 77 SDI 160 CCF 87 TopHt 63 QMD 5.1 TCuFt 1308 MCuFt 449 (forest 371... wait
STDINFO forest=505 Klamath; .sum ForType col=371). Validate TPA/BA/SDI/TopHt/QMD FIRST (no growth models needed) —
currently blocked by the tree-expansion bug above.
