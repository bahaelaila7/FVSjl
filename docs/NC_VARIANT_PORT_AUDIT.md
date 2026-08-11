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

## Chunk-3 DG — data COMPLETE + integration scoped (2026-08-11), ready to code
All DDS coefficients measured (docs/NC_CHUNK1_EXTRACTION.md chunk-3 section). Integration with the shared engine:
- jl PROVIDES the point-level density NC's dgf needs (density.point_ba=PTBAA, point_bal=PTBALT, point_ccf=PCCF,
  point_tpa=PTPA) — the shared standstats even computes PCCF specifically for the dgf DGPCCF term (doctrine #5 win).
- dgf! stores per-tree DDS in scratch.wk[2,i]=wk2 (CI pattern); the shared DDS→DG + serial-corr engine converts.
- ONE gap: PRD = ZRD(pt)/XMAXPT(pt) (point ZEIDE relative density) used ONLY by sp2/6/9 + redwood branches — jl
  has no point-Zeide-RD yet; add it (point Σ(D)^1.605 / point XMAXPT=Σ SDIDEF·point-BA-frac) OR compute inline in
  nc dgf!. nct01 has SP(sp2) trees ⇒ needed for nct01 validation.
- NC redwood DDS is a DIB²-diff form (like CI's diagr species) + POWER bark (nc_bratio) — the shared engine handles
  the diagr convention (CI proves it).
TO CODE (next): src/variants/klamath/diameter_growth.jl = nc_bratio(a,b,eqtype,d) + nc_dgcons!(s) [DGCON/DGDSQ per
sp, 3 branches: DGFOR/MAPLOC default, DGLAT2+site sp2/6/9, redwood ln(SITEAR)] + dgf!(s,::Klamath) [3-branch DDS →
wk2, 5-yr TDDS/2] + point-Zeide-RD. Then validate DDS per-tree vs FVSnc_clean dgf DEBUG dump on nct01.

## Chunk-3 DG — IMPLEMENTED (ce5111c) + first validation vs FVSnc_dbg (2026-08-11): NOT yet bit-exact
- ★ BUG FOUND + FIXED by reading nc/dgf.f: the 5-yr TDDS/2 conversion applies ONLY to the redwood + sp2/6/9
  branches, NOT the DEFAULT branch (nc/dgf.f has TDDS/2 inside CASE(12) and CASE(2,6,9), none in CASE DEFAULT).
  jl had applied /2 to ALL branches ⇒ DEFAULT species off by ln(2)≈0.693. Moved /2 into the two special branches.
- VALIDATION SETUP: instrumented FVSnc dgf.f (WRITE 'ZNC',I,ISPC,D,DDS after WK2(I)=DDS), recompiled dgf.o
  (gfortran-16 -std=legacy -w -fno-automatic -O0 -fPIC -I../../common), relinked → /workspace/.ncwork/FVSnc_dbg
  (ZNC dump in fort.16). Oracle dgf.f/dgf.o restored pristine after. jl side: env NC_DDS_DBG print in dgf!.
- FIRST COMPARISON (not yet aligned): jl sp2 D=10.39 DDS=2.50 vs live sp2 D~10.41 DDS~2.17; jl sp9 D=3.94 DDS=1.675
  vs live 1.171. jl HIGHER on the set-2 (sp2/6/9) species by ~0.3-0.5. ⚠ CAVEAT: both dumps mix CALIBRATION passes
  (backdated/modified DIAM) + the growth pass, and tree indices don't align 1:1 (jl re-sorts) ⇒ the comparison is
  NOT pass-aligned yet (doctrine #3 class). NEXT: align by (ISPC, exact DBH) on the GROWTH pass only (gate the dump
  on a growth flag or match DBH exactly), then diagnose the set-2 residual (candidates: the sp2/6/9 DGCON site form
  DGLAT2/DGSITE, the CRID/PBAL terms, or the -8.52-then-/2 order). DEFAULT-species DDS (post-/2-fix) still to compare.
⇒ Chunk-3 status: IMPLEMENTED + running to ch4; /2 bug fixed; per-tree DDS validation IN PROGRESS (not bit-exact
yet — set-2 species show a residual to diagnose). FVSnc_dbg durable for the aligned comparison.
