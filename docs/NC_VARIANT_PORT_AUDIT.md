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

## Chunk-3 DDS validation — refined diagnosis (2026-08-11): sp2/6/9 FORMULA verified faithful; residual = align/COR
Re-checked the jl sp2/6/9 branch term-by-term vs nc/dgf.f CASE(2,6,9): DDS = CONSPP + DGLD2·ALD + DGDSQ2·D²/1000 +
DGCR2·CRID + DGDBA2·PBAL/ln(D+1)/100 + DGBA2·ALPBA — EXACT match. So the ~0.3-0.5 jl-high on set-2 is NOT a formula
bug. The dump comparison was NOT valid: (a) fort.16 holds ZNC for ALL 10 nct01 cycles (the "last block" is a LATER
cycle with grown DBH 15-19", not cyc1); (b) both sides mix calibration passes (backdated DIAM, COR=0) + the growth
pass (which adds the DGSCOR-calibrated COR); (c) tree indices don't align (jl re-sorts). ⇒ residual candidates,
in order: (1) the growth-pass calibration COR (jl c.dg_cor set by the shared calibrate_diameter_growth! — verify it
computes NC COR like live's DGSCOR; the deterministic pre-COR DDS is the thing to match first), (2) CONSPP/DGCON
site inputs (SITEAR/ELEV/SLOPE for nct01 — verify jl's sp_site_index[2,6,9] + p.elevation/slope match live), (3)
per-tree CRID/PBAL/ALPBA. CLEAN NEXT: use nc/dgf.f's own DEBUG FORMAT-334 (DDS,CONSPP,ALD,ALPBA,CRID) — enable
DEBUG on FVSnc_dbg for ONE cyc1 growth-pass tree + match jl by exact (species,DBH) ⇒ isolates CONSPP vs the terms.
BETTER GATE (once ch4-8 land): the aggregate cyc1 .sum (TPA/BA/QMD) which averages out per-tree/pass noise, like CI.

## Chunk-4 height (nc/htgf.f) — structure characterized (2026-08-11), substantial port (NOT CI-reusable)
NC htgf = per-species SELECT CASE, using DG(I) (from dgf!) + the site-index height curve:
- CASE(12) redwood: DG10=DGLT/BRAT; LTHTG=EXP(1.412947 −0.000204·D² +0.31971·lnD +0.394005·ln(SINDX)
  +0.399888·ln(DG10) −0.451708·lnH)·0.5; HGBND height bound (H<217→1.0; 217-380 ramp to 0.1; ≥380→0.1);
  HTG=LTHTG·HGBND. (If H<4.5, DG10=0.1.)
- CASE DEFAULT (conifers+hardwoods): CALL FINDAG(I,ISPC,D1,D2,H,SITAGE,SITHT,AGMAX,HTMAX,HTMAX2,…) using HD1-4
  (P1-P4, the HT-DBH coeffs) + SITEAR → SITAGE (site age from H) + HTMAX. Then POTHTG from the site-curve
  (HGUESS−SITHT over SITAGE→SITAGE+5); H≥HTMAX→HTG=0.1. Scaled by SCALE=FINT/YR, XHT=XHMULT (MULTS keyword),
  exp(HTCON). Hardwoods (sp5/7/8/11) use their HD1-4; conifers use the site curve too (HD1-4=0 ⇒ different FINDAG path).
- NOT reusable from CI height_growth.jl (CI = NI-conifer exp form + Weibull, a DIFFERENT model). NC needs its own
  height_growth.jl porting FINDAG (nc/findag.f: site-index age↔height curve) + the redwood LTHTG + the DEFAULT
  POTHTG loop. FINDAG is the Dixon site-index curve (shared by the Pacific-coast/CA variants NC belongs to).
TO CODE (chunk 4): read nc/findag.f (the site-curve age/height solve) + the DEFAULT POTHTG iteration (lines 200-280);
port nc_findag + height_growth!(::Klamath) (redwood LTHTG + default FINDAG/POTHTG); validate HTG per-tree /aggregate.
Coefficients: HD1-4 measured (extraction doc); SITEAR from site chunk; HTCON/XHMULT default 0/1.

## ★★★ 2026-08-11 MILESTONE: NC nct01 control stand RUNS END-TO-END, cyc0 BIT-EXACT
All 8 growth+mortality+crown chunks implemented (ch1 species, ch2 site, ch3 DG, ch4 height, ch5 crown, ch6 regent,
ch7 mortality) + pipeline plumbing (FFE-inert gate, bark_intercept cols, v2t, volume placeholder). NC nct01
(control stand) runs full multi-cycle. cyc0 (1990 inventory) BIT-EXACT vs live: TPA 536, BA 77, SDI 160, TopHt 63,
QMD 5.1 (all ==). ⇒ chunks 1/2 + the stand summary are FAITHFUL.
REMAINING (3 items to bit-exact-or-cornered):
1. CYCLE LENGTH: live nct01 runs 5-YR cycles (1990,1995,...→2040); jl ran 10-yr (control.year=10 in species.jl is
   WRONG for NC — set to 5). nc/grinit IFINT=10 but the stand default cycle is 5 (verify: no TIMEINT in nct01.key
   ⇒ NC default cycle = 5, fix species.jl control.year/growth_fint = 5).
2. CCF = 0 in jl vs 87 live — the crown-competition-factor / per-tree CCF not wired for NC (crown.f ccfcal). Affects
   the .sum CCF column + RELDEN (density → DG/crown). Wire nc CCF (ccfcal) — likely reuse the shared ccfcal.
3. VOLUME = 0 (CI r4vol placeholder gives 0 for NC species). Real NC VEQNNC (R5/NVEL California equations) = chunk-8.
After cycle-length fix + CCF, validate the multi-cycle trajectory vs nct01.sum.save (5-yr, TPA 536→357 by 2040,
BA 77→301). This is the aggregate .sum gate.

## ★ 2026-08-11 TopHt under-growth ROOT: site index defaults to site_lo (chunk-2 site-input bug)
After the htg_period=5 fix, BA/SDI/QMD track live closely but TopHt DECREASES (jl 64→57 vs live 72→91). Traced:
the largest-DBH DF tree gets SITE INDEX = 50 = its site_lo DEFAULT, not the input value → FINDAG SITAGE=86 (old) →
POTHTG=1.45 (tiny) → height barely grows. ROOT: nc_sitset! (site_index.jl) fills sp_site_index from site_lo when
p.site_species/site index isn't set from the SITECODE/STDINFO input. cyc0 is bit-exact because initial HEIGHTS come
from the tree input, but the growth site curve uses the wrong (low) index. nct01 SITECODE/STDINFO should set the
real site index (STDINFO field ~84). FIX: parse the SITECODE/STDINFO site species+index into p.site_species +
p.sp_site_index[isisp] so nc_sitset! interpolates from the REAL index. This should fix TopHt (and refine DG, which
reads SITEAR). Then re-validate multi-cycle. cyc0 fully bit-exact; BA/SDI/QMD good; TopHt gated on the site index.

## Site-index/TopHt — CONFIRMED root (2026-08-11 test): forcing site index 84 (STDINFO fld2) makes TopHt GROW
Temp test (reverted): default site index 50→84 ⇒ TopHt 63→57(decrease) becomes 63→72(grows); live 63→91. So the
site index IS the TopHt driver — the site_lo default (50) froze/shrank it. Proper FIX (chunk-2 completion): NC has
NO SITECODE in nct01 ⇒ live DEFAULTS the site species to DF (sp3, per the .out "SITE SPECIES=DF CODE=3") and derives
the site index from the habitat/ecocls SITE (or estimates from the site-species trees). Wire in nc_sitset!: default
site_species=DF; get the real site index (STDINFO habitat 84 → ecocls SITE per PA, OR the sitcind.f site-tree
estimate) → p.sp_site_index[DF] → interpolate. Residual after 84 (jl 72 vs live 91) ⇒ the true index is higher
and/or a small height-model residual — measure live's exact SITEAR (instrument FVSnc_dbg sitcind) to pin it.
STATUS: NC growth port ~90% — cyc0 FULLY bit-exact; BA/SDI/QMD track live closely; TopHt gated on the site-index
derivation (confirmed root); volume = placeholder (real NC VEQNNC R5/NVEL pending). Clear remaining scope.

## Chunk-8 VOLUME — VEQNNC measured from live nct01.out (2026-08-11): NVEL 500WO2W conifers + 500DVEW hardwoods
NC volume = National Volume Estimator Library, region/forest 500 (California). Per-species equation numbers:
- OS 500WO2W108 · SP 500WO2W117 · DF 500WO2W202 · WF 500WO2W015 · IC 500WO2W081 · RF 500WO2W020 · PP 500WO2W122 ·
  RW 500WO2W211  (conifers + IC → 500WO2W* profile — the R5 West-side 2-point taper).
- MA 500DVEW361 · BO 500DVEW818 · TO 500DVEW631 · OH 500DVEW981  (hardwoods → 500DVEW* woodland D²H — the SAME
  DVEW/r4d2h family CI/UT/TT already implement; reuse r4d2h_vol with the 500DVEW coefficients).
⇒ chunk-8 port: (a) hardwoods = REUSE the shared DVEW/r4d2h (like CI's compute_volumes_ci! DVEW branch) with the
NC 500DVEW eqnums; (b) conifers = port the 500WO2W NVEL profile (West-side 2-point — check if it's the Flewelling
FW2 family jl already has, or a distinct NVEL call). Replace the CI-r4vol placeholder (compute_volumes!(::Klamath)
-> compute_volumes_ci!) with compute_volumes!(::Klamath) using these. Validate TCuFt/MCuFt vs live nct01.sum
(1990 TCuFt 1308/MCuFt 449). All VEQNNC measured; DVEW reusable, WO2W is the new piece.

## Chunk-8 volume — implementation scoping (2026-08-11):
- DVEW hardwoods (MA361/BO818/TO631/OH981): jl r4d2h_vol (Chojnacky INT-339 D²H) exists (TT/UT/CI/CR) but is keyed
  by species code — the NC R5 hardwood codes 361/818/631/981 are NOT in the existing table (those variants use
  066/475/998/133/065 etc.). ⇒ add the 500DVEW coefficients for the 4 NC hardwoods to r4d2h (or a NC table).
- WO2W conifers (OS/SP/DF/WF/IC/RF/PP/RW): "WO2W" is a DISTINCT NVEL eqtype vs jl's "FW2W" (Flewelling, cr_fw2_vol)
  and "MATW" (Matney, r4vol). VERIFY whether 500WO2W maps to an existing jl profile or is a new NVEL taper to port
  (check ForestVegetationSimulator/volume/NVEL for the WO2W equation). This is the substantial chunk-8 piece.
⇒ compute_volumes!(::Klamath) = DVEW branch (r4d2h + NC hardwood coeffs) + WO2W branch (verify/port), validate
TCuFt/MCuFt per-tree + .sum vs FVSnc_dbg (1990 1308/449). Replaces the CI-r4vol placeholder. All eqnums measured.

## Chunk-8 volume — REFINED (2026-08-11): NC 500* = Region-5 California NVEL, NEW to jl (not R4 reuse)
The NC DVEW codes 361/818/631/981 are NOT in jl's R4 r4d2h_vol1 branches (which cover 064/066/106/475/... R4
species). ⇒ NC's 500DVEW (hardwoods) + 500WO2W (conifers) are the REGION-5 (California) NVEL equations, distinct
from the R4 Matney/Chojnacky and R1 Flewelling jl already has. chunk-8 = a genuine NVEL R5 port: locate the WO2W +
R5-DVE coefficients in ForestVegetationSimulator/volume/NVEL (voleqdef/the R5 taper source), port the taper + the
D²H woodland eq for the 12 NC species, wire compute_volumes!(::Klamath), validate TCuFt/MCuFt vs FVSnc_dbg. All 12
eqnums measured (500WO2W108/117/202/015/081/020/122/211 conifers+IC; 500DVEW361/818/631/981 hardwoods). Substantial
but self-contained; the growth port is unaffected (volume is a reporting/.sum column, output-only).

## Chunk-8 volume — NVEL sources located (2026-08-11): WO2W=Flewelling(r6vol/fwinit), DVE=dvest.f
voleqdef.f confirms NC uses 500WO2W* (conifers) + 500DVEW* (hardwoods). The NVEL profile routines:
- WO2W → fwinit.f + r6vol.f/r6vol1.f (Region-6 FLEWELLING West-side 2-point taper) — the SAME Flewelling family
  jl's cr_fw2_vol already implements (jl uses it for CI/EM/BM/UT "FW2W"). ⇒ the WO2W conifer volume likely REUSES
  cr_fw2_vol with the WO2W (Flewelling-West) profile coefficients — verify the fwinit coeff set for WO2W vs FW2W.
- DVEW → dvest.f (the DVE woodland/D²H estimator) — region-5 species coefficients (361/818/631/981) in dvest.f.
⇒ chunk-8 is MORE reusable than feared: WO2W ≈ jl cr_fw2_vol (Flewelling) + WO2W coeffs; DVEW = port dvest.f R5
coeffs (small). compute_volumes!(::Klamath) = cr_fw2_vol(conifers, WO2W) + dvest(hardwoods). Measure the WO2W/DVE
coefficients from fwinit.f/dvest.f, validate TCuFt/MCuFt per-tree vs FVSnc_dbg (1990 1308/449). Doctrine-#5 reuse.

## Chunk-8 volume — FINAL scope (2026-08-11): genuine NVEL R5/R6 port (cr_fw2_vol INGY-only, needs West-side coeffs)
Tested: cr_fw2_vol("500WO2W202",...) returns 0 — jl's Flewelling is the INGY (Inland NW) + R2/R3 subregion set only
(_fw2_is_ingy || 22≤jsp≤29 gate). NC's 500WO2W is the Region-6 Flewelling WEST-SIDE profile (fwinit.f SHP_C2 west
coeffs), a DIFFERENT coefficient set + species map than jl has. ⇒ chunk-8 = a real NVEL port:
- WO2W conifers: extend cr_fw2_vol (or a nc_fw2) with the R6 West-side SHP/taper coefficients from fwinit.f + the
  NC species→Flewelling-jsp map (framework reused; coefficients new).
- DVE hardwoods (MA/BO/TO/OH): port dvest.f's DVE D²H estimator + the R5 coefficients for 361/818/631/981.
Substantial but SELF-CONTAINED + OUTPUT-ONLY (volume is a .sum reporting column; the growth port is complete and
unaffected). All eqnums measured; sources = volume/NVEL/{fwinit.f, r6vol.f, dvest.f}. Validate TCuFt/MCuFt per-tree
+ .sum vs FVSnc_dbg (1990 1308/449). This is the last NC chunk; growth (ch1-7) tracks live end-to-end.

## NC growth verdict (2026-08-11): cyc0 BIT-EXACT; multi-cycle ~2-5% = plausibly the cornered DGSCOR class
Aggregate .sum nct01 control: cyc0 all 6 density cols BIT-EXACT; multi-cycle BA ~2% high, QMD ~1-2%, TopHt ~5% low,
CCF ~12% low (CCF/TopHt downstream of the DG-distribution + height realization). The BA/QMD ~2% is the SAME
magnitude as CI #142's accepted cornered DGSCOR/ZZRAN straddle (~2% over-kill) — so NC growth is plausibly
bit-exact-or-cornered, NOT a gross bug. TO CONFIRM cornered-vs-small-bug (doctrine #1): per-tree deterministic DDS
jl-vs-FVSnc_dbg on ONE cyc1 growth-pass tree matched by (species, exact DBH) — gate the dgf DEBUG dump on the
growth pass (not calibration) to avoid the earlier mixed-pass alignment error; if the deterministic DDS matches,
the residual is the DGSCOR realization (cornered). The TopHt-5%/CCF-12% may need a 2nd look (height POTHTG or the
site-index HTCALC-conversion approximation — I hardcoded the DF=90-default SI() rather than porting the full
sitset.f DO-30 HTCALC conversion; exact for nct01, approximate if a SITECODE sets a non-default site species/index).
⇒ NC growth port: cyc0 bit-exact + multi-cycle tracks live in the cornered range. Remaining to fully close:
(1) per-tree DDS confirm-cornered, (2) full sitset HTCALC site-conversion (for SITECODE stands), (3) chunk-8 volume.
