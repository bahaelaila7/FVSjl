# NC (Klamath Mountains) chunk-1 coefficient extraction — WIP (2026-08-11)

Ground-truth coefficient values MEASURED (read directly from literal Fortran DATA blocks in
`bin/FVSnc_buildDir/`, doctrine #2) for building `data/klamath/species_coefficients.csv` (41 cols × 12 sp) +
`dg_coefficients.jl`. Species order = the NC JSP order (1..12). Validate each against `FVSnc_clean` on `nct01`.

## Species identity (blkdat.f:142-154) — MEASURED
| idx | JSP | FIA | PLANTS  |
|-----|-----|-----|---------|
| 1 | OS | 299 | 2TN    |
| 2 | SP | 117 | PILA   |
| 3 | DF | 202 | PSME   |
| 4 | WF | 015 | ABCO   |
| 5 | MA | 361 | ARME   |
| 6 | IC | 081 | CADE27 |
| 7 | BO | 818 | QUKE   |
| 8 | TO | 631 | LIDE3  |
| 9 | RF | 020 | ABMA   |
| 10| PP | 122 | PIPO   |
| 11| OH | 998 | 2TB    |
| 12| RW | 211 | SESE3  |

## Site index (sichg.f:14-24) — MEASURED
- site_lo (SIMIN): 50, 40, 50, 30, 50, 30, 30, 50, 30, 40, 50, 50
- site_hi (SIMAX): 150, 120, 150, 130, 100, 130, 70, 90, 130, 120, 90, 150
- SICHG A (ht-adj): 10, 12, 10, 10, 3, 10, 6, 4, 10, 12, 4, 10
- SICHG B: -0.08, -0.05, -0.08, -0.07, -0.02, -0.05, -0.05, -0.03, -0.06, -0.05, -0.03, -0.08
- REFLOC: B,T,B,B,B,B,B,B,B,T,B,B ; IREFAG: 12*50 (site reference age 50)
- SITEAR default when unset = 90. (sitcind.f:120,129)

## Large-tree DG DDS (dgf.f:88-134,231) — MEASURED. ⚠ NC-SPECIFIC dgf structure (NOT CI drop-in):
per-species diameter SECTIONS via ISCT(ISPC,1/2); DGHAH habitat term; TWO coefficient sets; sp12 REDWOOD special CONSPP.
- **Set 1** (species 1..12):
  - DGLD:   0.88425, 0.0, 0.86990, 1.01718, 1.14082, 0.0, 1.23911, 0.99531, 0.0, 0.96865, 0.99531, 0.0
  - DGCR:   2.83271, 0.0, 2.96040, 3.01884, 2.82796, 0.0, -1.20841, 2.08524, 0.0, 1.5466, 2.08524, 0.0
  - DGCRSQ: -0.84141, 0.0, -1.08219, -1.12464, -2.14739, 0.0, 2.31782, -0.98396, 0.0, 0.07152, -0.98396, 0.0
  - DGDBAL: -0.00358, 0.0, -0.00443, -0.00257, -0.00126, 0.0, -0.00199, -0.00147, 0.0, -0.00408, -0.00147, 0.0
  - DGBA:   0.0, 0.0, -0.01744, -0.16596, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  - DGHAH:  0.0, 0.0, 0.0, 0.0, 0.56348, 0.0, 0.0, 0.50155, 0.0, 0.0, 0.50155, 0.0
  - DGPCCF: 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.00018, 0.0, -0.00002, -0.0018, 0.0
- **Set 2** (used where set-1 is 0.0 / by section — sp2/6/9/12 SP/IC/RF/RW etc.):
  - DGLD2:  1.52284, 1.26883, 1.92426, 1.53339, 1.26883, 1.41389, 1.41389, 1.52284, 1.53339, 1.63568, 1.41389, 0.0
  - DGCR2:  0.51033, 0.27986, 0.40047, 0.35739, 0.27986, 0.32660, 0.32660, 0.51033, 0.35739, 0.27516, 0.32660, 0.0
  - DGDSQ2: -0.26590, -0.35325, -0.44612, -0.47442, -0.35325, -0.48938, -0.48938, -0.26590, -0.47442, -0.54162, -0.48938, 0.0
  - DGBA2:  -0.35579, 0.0, -0.24596, -0.12359, 0.0, -0.25287, -0.25287, -0.35579, -0.12359, -0.20902, -0.25287, 0.0
  - DGDBA2: 0.0, -0.79922, 0.0, -0.44256, -0.79922, -0.16000, -0.16000, 0.0, -0.44256, -0.32497, -0.16000, 0.0
- DGDS (dgf.f:231) + DGFOR, ISCT section indices, CONSPP per-species: STILL TO READ (dgf.f:135-290).

## Bark (bratio.f BRKRAT(4,12)) — MEASURED. ⚠ MIXED TYPES; sp12 (RW) is POWER ⇒ needs a special branch (BM#140/CI class)
BRKRAT row = (idx, a, b, eqtype). Three forms: type1 DBT=a+b·D ⇒ BRATIO=(D−DBT)/D=(1−b)−a/D ; type2 DIB=a+b·D ⇒
BRATIO=b+a/D ; type3 DIB=a·D^b ⇒ BRATIO=a·D^(b−1) [POWER]. Shared bark_ratio(bark_a,bark_b,sp,d)=bark_b+bark_a/d.
| sp | a | b | type | ⇒ bark_b | bark_a | note |
|----|-----|-----|----|-----|-----|-----|
| 1 OS | 0.1429 | 0.1137 | 1 | 0.8863 | −0.1429 | linear |
| 2 SP | 0.1429 | 0.1137 | 1 | 0.8863 | −0.1429 | linear |
| 3 DF | 0.1045 | 0.1661 | 1 | 0.8339 | −0.1045 | linear |
| 4 WF | 0.1593 | 0.1089 | 1 | 0.8911 | −0.1593 | linear |
| 5 MA | −0.01348 | 0.98155 | 2 | 0.98155 | −0.01348 | linear |
| 6 IC | −0.0549 | 0.1626 | 1 | 0.8374 | 0.0549 | linear |
| 7 BO | −0.26824 | 0.95767 | 2 | 0.95767 | −0.26824 | linear |
| 8 TO | −0.26824 | 0.95354 | 2 | 0.95354 | −0.26824 | linear |
| 9 RF | 0.1593 | 0.1089 | 1 | 0.8911 | −0.1593 | linear |
| 10 PP | 0.4448 | 0.1033 | 1 | 0.8967 | −0.4448 | linear |
| 11 OH | −0.26824 | 0.95767 | 2 | 0.95767 | −0.26824 | linear |
| 12 RW | 0.70120 | 1.04862 | **3 POWER** | — | — | **nc_bratio: 0.70120·D^0.04862 (≈0.78-0.81, near 0.80 clamp) — NEEDS a POWER branch in the shared DDS→DG apply-loop like bm_bratio/ci_bratio** |
⇒ sp1-11 encode into c.bark_a/c.bark_b (bark_imap=linear); sp12 RW needs a POWER special-function branch at EVERY
shared bark site (bark dispatch + backdating), else DG bias on redwood (the BM#140 lesson). Verify vs FVSnc_clean.

★ CSV REPRESENTATION RESOLVED (2026-08-11): store the RAW BRKRAT values + eqtype (cleanest, all 3 forms preserved):
- bark1 (=BRKRAT(2,sp), raw a): 0.1429, 0.1429, 0.1045, 0.1593, -0.01348, -0.0549, -0.26824, -0.26824, 0.1593, 0.4448, -0.26824, 0.70120
- bark2 (=BRKRAT(3,sp), raw b): 0.1137, 0.1137, 0.1661, 0.1089, 0.98155, 0.1626, 0.95767, 0.95354, 0.1089, 0.1033, 0.95767, 1.04862
- bark_imap (=BRKRAT(4,sp), eqtype): 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 3
NC needs its OWN `nc_bratio(a,b,eqtype,d)` in chunk-3 diameter_growth (the shared bark_ratio only does type-2 b+a/D):
type1 → (D−(a+b·D))/D = 1−b−a/D ; type2 → b+a/D ; type3 → a·D^(b−1). Dispatch `_nc_dg = s.variant isa Klamath`
at EVERY shared DDS→DG bark site (apply-loop + backdating + mortality bark), like _ci_dg/_bm_dg. CI used bark_imap=4
for its all-POWER set; NC's per-species eqtype 1/2/3 selects the form within nc_bratio.

## Small-tree regent (regent.f:94-97) — MEASURED
- st_diam (DIAM): 0.3, 0.4, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.3, 0.5, 0.2, 0.3
- st_xmax (XMAX): 5,5,5,5,5,5,5,5,5,5,5, 10  (sp12 RW = 10; blend cap)
- st_xmin (XMIN): 12*2.0
- DGMIN: 3,3,3,3,3,3,3,3,3,3,3, 7  (sp12 RW = 7 — NC uses DGMIN, note name vs CI's st_dgmax; verify role in regent.f)
- REGYR = 5.0 (small-tree growth period; SCALE=FNT/REGYR, SCALE2=YR/FNT). Small-tree HTG model = blended, NOT simple
  per-sp DATA — port the regent.f logic (chunk 6). ht-inc model based on 5-yr growth data.

## Large-tree height (htgf.f:58-68) — MEASURED. Hardwoods (sp5/7/8/11) use HD1-4 curve; CONIFERS 0.0 ⇒ site-index path
- HD1: 0,0,0,0, 4.4666, 0, 4.80758, 4.9684, 0,0, 4.9684, 0
- HD2: 0,0,0,0, -0.00179, 0, -0.00224, -0.004057, 0,0, -0.004057, 0
- HD3: 0,0,0,0, 0.002048, 0, -0.000513, 0.000924, 0,0, 0.000924, 0
- HD4: 0,0,0,0, -7.9428, 0, -7.729644, -10.45158, 0,0, -10.45158, 0
⇒ NC htgf chunk-4 = shared conifer site-index potential-height (sp1-4,6,9,10,12) + variant hardwood HD1-4 curve
(sp5 MA, sp7 BO, sp8 TO, sp11 OH). Structured port, not coefficient-only.
NOTE: these HD1-4 are height-GROWTH coeffs; the species_coefficients.csv ht1/ht2/htdbh_* (height-DIAMETER curve) are
SEPARATE — still to read from htdbh.f / the HT-DBH DATA.

## Background mortality (morts.f:99-104) — MEASURED. ⇒ REUSES the EM/UT mortality path (doctrine #5)
NC background mortality = EM/UT form RI=0.5·(1/(1+exp(PMSC+PMD·D))) + Zeide self-thin (LZEIDE=T). sp1-11 PMSC/PMD
are IDENTICAL to the EM/UT coefficients ⇒ chunk-7 mortality reuses the shared easternmontana/utah mortality! path.
- PMSC: 6.5112, 6.5112, 7.2985, 5.1677, 9.6943, 5.1677, 5.9617, 9.6943, 5.1677, 5.5877, 5.1677, **2.59680**
- PMD: -.0052485, -.0052485, -.0129121, -.0077681, -.0127328, -.0077681, -.0340128, -.0127328, -.0077681,
  -.005348, -.0077681, **+0.51261**
- ⚠ sp12 RW (redwood) is SPECIAL: PMSC=2.5968, PMD=+0.51261 (POSITIVE slope — redwood-specific; verify the RI
  form still holds or if redwood takes a separate branch in nc/morts.f). Zeide self-thin: reuse _em_tn10_iter etc.

## HT-DBH curve (htdbh.f) — STRUCTURE NOTED (chunk-4). ⚠ FOREST-SPECIFIC, not simple per-species ht1/ht2
NC htdbh depends on IFOR forest code (1=Klamath/505, 2=Six Rivers/510, 3=Trinity/514, 4=Siskiyou/611, 5=Hoopa/705,
6=Simpson/800, 7=BLM Coos Bay/712) — e.g. a SISKIY coefficient array (htdbh.f:48). So the species_coefficients.csv
ht1/ht2/htdbh_* columns need the FOREST-appropriate curve (nct01 forest=371... verify which IFOR). Port htdbh.f
logic (CRATET height-dubbing + REGENT diameter est, MODE 0/1) rather than a flat per-species ht1/ht2. Still to read
the full SISKIY + per-forest coefficient arrays.

## SDImax (sdical.f) — computed, not DATA: BAMAX=XMAX·0.5454154·PMSDIU (sdical.f:204); SDIDEF default 0 (grinit:72).
Per-species SDImax source (SDIDEF fill / habitat table) still to locate — likely forest/habitat-driven like CI's R4SDI.

## dg_resid_sd (SIGMAR, blkdat.f:180-184) — MEASURED (active DATA, not the commented-out higher set):
0.3300, 0.2713, 0.3300, 0.3300, 0.3306, 0.3513, 0.3541, 0.3558, 0.3136, 0.1954, 0.3558, 0.6178  (sp12 RW=0.6178)
(NC DGSD=2.0 from grinit; dg_stddev_bound must also be set 2.0 — the BM/CI DGSD field-disconnect lesson.)

## Volume (grinit.f:85-94) — vollib-driven, NOT simple DATA: TOPD/DBHMIN/BFTOPD/SCFTOPD init 0, filled by the
volume library defaults (VOLKEY/vollib per species+forest, like CI's r4vol/VEQNNC). Chunk-8 = wire NC's VEQNNC/
vollib defaults (read from a live nct01.out VOLUME echo, as CI did). Not a species_coefficients.csv DATA read.

## SDImax — ECOCLASS/habitat-driven (RESOLVED source, chunk 2): SDIDEF(sp) is filled from ecocls.f's per-ecological-
class SDIMX table (ecocls.f:46 SITE/SDIMX(NENTRY)) OR from stand-input MAX_SDI (dbsstandin.f:798 RSTANDDATA(36)).
grinit SDIDEF=0 default; sdical weights XMAX=Σ SDIDEF(I)·BAXSP(I) (BA-weighted, sdical.f:170) then BAMAX=XMAX·
0.5454154·PMSDIU. So NC SDImax = the ecocls SDIMX-per-ecoclass table (like CI's R4SDI habitat SDImax). Chunk-2 =
port ecocls SITE/SDIMX + the ecological-class crosswalk (nct01 forest 371 → its ecoclass → SDIMX). Read ecocls.f.

## STILL TO MEASURE (chunk-1 completion) — remaining are STRUCTURAL (not flat DATA), deeper reads:
- dbh_max (find NC's diam-cap name); is_sprouting (which sp sprout — hardwoods BO/TO/OH + RW redwood); varmrt_varadj
  (varmrt.f); the forest-specific htdbh coeff arrays (htdbh.f SISKIY + per-forest); ht1/ht2/wykoff_ht2 (htgf.f); mort_bkgd (morts.f); sdi_max_default
  (sdical.f: BAMAX=XMAX·0.5454154·PMSDIU ⇐ per-sp XMAX SDI); volume stump/top_dib/dbh_min/scf_*/bf_* (grinit.f/
  vollib); htdbh_* (htdbh.f); varmrt_varadj; is_sprouting; dg_resid_sd (SIGMAR).
- species_translation.csv: full FIA-species→NC-12 crosswalk (dgf.f OSPMAP / the FIA map, NOT just the 12).

## SDImax detail (ecocls.f) — 90-entry PLANT-ASSOCIATION table (chunk-2 data extraction):
ecocls.f DATA (PA, SCIEN, SDIMX, SPC, SITE, NUMBR, IFLAG, FVSSEQ), NENTRY=90 plant associations (CDC411/CDC412/
CDC421/... = PSME-ABCO-PIJE etc.). Each PA → SDIMX (SDImax) + SITE + species. So NC SDImax = a 90-row PA→SDIMX
table + a PA-STRING crosswalk (the IE #143 habtyp pattern; FIA PV_CODE / plant-assoc string → PA row → SDIMX).
Chunk-2 = extract the 90-entry table (ecocls.f DATA blocks, I=1..90) + port the PA-string lookup + BA-weight
XMAX=Σ SDIDEF·BAXSP. Substantial (90 rows) — extract during chunk-2 implementation, validate SDI/BAMAX on nct01.

## PORT SCOPE (fully mapped 2026-08-11) — NC is a genuine ~8-chunk port, NOT a coefficient swap:
- ch1 species coeffs: flat DATA all MEASURED (this doc). ~16 of 41 species_coefficients.csv cols are flat; the rest
  (sdi_max/ht1/ht2/htdbh/volume) are STRUCTURAL (below) ⇒ the CSV can't be faithfully assembled standalone — build
  it alongside chunks 2/4/8.
- ch2 site/SDImax: 90-entry ecocls PA→SDIMX table + PA-string crosswalk (IE#143 pattern). REAL sub-port.
- ch3 DG: dgf.f ISCT diameter-sections + DGHAH habitat + 2 coeff sets + sp12 REDWOOD special CONSPP. REAL port.
- ch4 height: conifer site-index potential-ht (shared) + hardwood HD1-4 (sp5/7/8/11) + FOREST-specific htdbh (IFOR 1-7). REAL.
- ch5 crown / ch6 regent: per CI template + NC regent blend (SCALE=FNT/REGYR=5).
- ch7 mortality: REUSES easternmontana/utah mortality! (PMSC/PMD measured, EM/UT-form, Zeide self-thin) + sp12 RW special. MOSTLY SHARED.
- ch8 volume: vollib VEQNNC defaults (echo from nct01.out, like CI). REAL.
⇒ Reuse: mortality (EM/UT), conifer height (site path), linear bark (sp1-11). Variant-specific: DG, hardwood ht,
htdbh, ecocls SDImax, sp12 redwood (POWER bark + special mort). This is a multi-session implementation effort.

## NEXT
Finish reading the above arrays → assemble data/klamath/species_coefficients.csv + species_translation.csv →
wire src/variants/klamath/species.jl + site_index.jl → validate SITE INDEX on nct01 vs FVSnc_clean (chunk-2 gate).
Then NC-specific dg_coefficients.jl + diameter_growth.jl (port the ISCT-section/DGHAH/2-set dgf logic) → chunk-3.

## ═══ CHUNK-3 DG (nc/dgf.f) — DDS spec + coefficient arrays (MEASURED 2026-08-11) ═══
NC dgf: DGCONS entry loads site/forest-specific DGCON+DGDSQ per species; dgf loops species×trees computing
LN(DDS) into WK2. 5-yr rate: after DDS, TDDS=EXP(DDS); DDS=LOG(TDDS/2.0) (the /2 gives the 5-yr increment).
Clamp DDS≥-9.21 (default) / ≥-8.52 (sp2,6,9). Per-tree vars: ALD=ln(D); CRID=(ICR²/ln(D+1))/1000, =1.8 if D<2;
CR=ICR·0.01; BAL=(1-PCT/100)·BA; HOAVH=min(HT/AVH,1.5); PBA=PTBAA(pt); PBAL=PBA·(1-PCT/100) (→BAL if ≤0);
ALPBA=ln(PBA); PRD=ZRD(pt)/XMAXPT(pt) (point Zeide RD); PCCF(pt).

DGCONS setup (per species, IFOR=forest 1-7; nct01 forest 505=IFOR 1):
- CONSPP(default) = DGCON + COR + DGCCFA·ALRD + DGBA·ALBA ; CONSPP(sp12) = DGCON only.
- DGCON(default) = DGFOR(MAPLOC(IFOR,sp), sp) ; DGDSQ(default) = DGDS(MAPDSQ(IFOR,sp), sp).
- DGCON(12 RW) = -3.502444 + 0.415435·ln(SITEAR); DGCON(2,6,9) = DGLAT2(5,sp)+DGEL2·ELEV+DGSLP2·SLOPE+
  DGSLQ2·SLOPE²+DGSITE·SITEAR (ILAT=5).

DDS branches:
- DEFAULT (sp1,3,4,5,7,8,10,11): DDS = CONSPP + DGLD·ALD + CR·(DGCR+CR·DGCRSQ) + DGDSQ·D² + DGDBAL·BAL/ln(D+1)
  + DGPCCF·PCCF + DGHAH·HOAVH ; sp4(WF): DDS -= 0.15032.
- sp2,6,9 (SP/IC/RF): DDS = CONSPP + DGLD2·ALD + DGDSQ2·D²/1000 + DGCR2·CRID + DGDBA2·PBAL/ln(D+1)/100 + DGBA2·ALPBA.
- sp12 (RW): DGLT=EXP(CONSPP +0.185911·ln(D) -0.000073·D² -0.001796·PBAL -0.42078·PRD +0.589318·ln(CR·100)
  -0.000926·SLOPE·100 -0.002203·(SLOPE·100)·cos(ASPECT)); BRAT=nc_bratio(12,D); DDS=LN((D+DGLT)²·BRAT² - (D·BRAT)²)
  + COR + LN(COR2). (redwood uses the POWER bark BRATIO directly + a DIB-squared-diff DDS, then /2 for 5-yr.)

Forest-dependent DGCONS arrays (dgf.f):
- DGCCFA(12): -0.06784, 0,0,0,0,0,0,0,0,0,0,0  (only sp1).
- MAPLOC(7 forest,12 sp) col-major (sp-major rows), each row = sp's loc-class per forest 1-7:
  sp1 [1,1,1,2,3,3,2]; sp4 [1,1,1,2,3,3,2]; sp5 [1,1,1,1,1,1,1]... (sp2,3,6-12 all [1,1,1,1,1,1,1] except sp1/sp4).
  Actually rows: sp1=1,1,1,2,3,3,2 / sp2=1×7 / sp3=1,1,1,2,3,3,2 / sp4=1,1,1,2,3,3,2 / sp5-12=1×7.
- MAPDSQ(7,12): sp1=1,1,1,1,2,2,1 / sp4=1,1,1,1,2,2,1 / all others=1×7.
- DGFOR(6 locclass,12 sp): sp1[-2.00201,-2.19449,-1.84083,0,0,0] sp2[0×6] sp3[-2.54402,-2.41928,-2.75656,0,0,0]
  sp4[-1.88042,-2.06853,-1.69815,0,0,0] sp5[-1.69950,0,0,0,0,0] sp6[0×6] sp7[-2.68349,0×5] sp8[-0.94563,0×5]
  sp9[0×6] sp10[-4.6744,0×5] sp11[-0.94563,0×5] sp12[0×6].
- DGDS(4 idx,12 sp): sp1[-0.000328,-0.000248,0,0] sp2[0,0,0,0] sp3[-0.000313,0,0,0] sp4[-0.000356,-0.000268,0,0]
  sp5[-0.000875,0,0,0] sp6[0,0,0,0] sp7[-0.000338,0,0,0] sp8[-0.000373,0,0,0] sp9[0,0,0,0] sp10[-0.000728,0,0,0]
  sp11[-0.000373,0,0,0] sp12[0,0,0,0].
- STILL TO READ: DGLAT2(5,12)@182, DGSLP2@201, DGEL2@207, DGSITE@210, DGSLQ2, COR2 (sp2/6/9 + redwood), DGDSQ2/
  DGCR2/DGDBA2/DGBA2 (have from set-2), DGLD2 (have). + ELEV/SLOPE/ASPECT from stand; ICR/ITRE/PTBAA/ZRD/XMAXPT/
  PCCF stand-density arrays (shared engine — verify jl provides point-BA PTBAA + point-Zeide ZRD/XMAXPT).
NC dgf! port: nc_dgcons!(s) (DGCON/DGDSQ per sp) + dgf!(s) (3-branch DDS) + nc_bratio for the bark. Validate DDS
per-tree vs a live FVSnc_clean dgf DEBUG dump on nct01 (add WRITE at WK2(I)=DDS).

## Chunk-3 sp2/6/9 + redwood DGCON site arrays (dgf.f:182-217) — MEASURED (completes chunk-3 data):
DGCON(2,6,9) = DGLAT2(5,sp) + DGEL2·ELEV + DGSLP2·SLOPE + DGSLQ2·SLOPE² + DGSITE·SITEAR  (ILAT=5 ⇒ 5th col).
- DGLAT2(5,sp) [the ILAT=5 element per species]: sp1 0.1630, sp2 -0.4297, sp3 -0.1043, sp4 0.1434, sp5 -0.4297,
  sp6 0.0540, sp7 0.0540, sp8 0.1630, sp9 0.1434, sp10 -0.3995, sp11 0.0540, sp12 0.0000.
- DGSLP2(12): 0,0,0,0,0,0,0,0,0, 0.80370, 0,0  (only sp10).
- DGSLQ2(12): 0,0,0,-0.83400,0,0,0,0,-0.83400,0,0,0  (sp4, sp9).
- DGEL2(12): 0,0,0,0,0,0,0,0,-0.00700,0,0,0  (sp9).
- DGSITE(12): 0.47932, 0.01401, 0.56356, 0.47360, 0.20189, 0.01200, 0.32093, 0.00659, 0.00734, 1.10842, 0.00659, 0.
- COR2(12): the DGSCOR "COR2" calibration multiplier — DGCON += LN(COR2) when LDCOR2 & COR2>0 (dgf.f:508); redwood/
  sp2,6,9 also `DDS += COR + LN(COR2)`. Init COR2=1.0 (⇒ LN=0, no effect) until DGSCOR calibrates ⇒ port as 1.0
  baseline (calibration refinement later, the never-FFI-RNG class). COR (additive) likewise 0 baseline.
⇒ ALL chunk-3 DGCONS + DDS coefficients now MEASURED. Next = WRITE src/variants/klamath/diameter_growth.jl:
nc_dgcons!(s) [DGCON/DGDSQ per sp, 3 branches] + dgf!(s) [3-branch DDS → WK2] + nc_bratio(a,b,eqtype,d), integrate
with the shared DDS→DG + serial-corr engine (study CI diameter_growth.jl), validate DDS per-tree vs FVSnc_clean.

## ═══ CHUNK-4 HEIGHT (nc/htgf.f + findag.f + htcalc.f) — MEASURED 2026-08-11 ═══
Model = htgf driver → FINDAG (age solve) → HTCALC (site-height curve). SCALE=FINT/YR, XHT=XHMULT (MULTS kw, def 1).
**HTCALC(SINDX,ISPC,AG)→HGUESS** (site height at age AG); 6 species-group branches:
- CASE(1,3,12) OS/DF/RW: Z=2500/(SI−4.5); A=−0.954038+0.109757·Z; B=0.055818+0.0079224·Z; C=−0.0007338+0.0001977·Z;
  HGUESS = AG²/(A+B·AG+C·AG²) + 4.5.
- CASE(4,6,9) WF/IC/RF: X1=38.0202·AG^(−1.05213)·EXP(0.009557·AG); X2=101.842894·(1−EXP(−0.001442·AG^1.679259));
  HGUESS = (SI−69.91+X1·X2)/X1 + 4.5.
- CASE(5) MA: HGUESS = SI/(0.375 + 31.233/AG).
- CASE(7) BO: A=√AG−√50; HGUESS = (SI·(1+0.322·A) − 6.413·A)·0.80.
- CASE(8,11) TO/OH: HGUESS = SI/(0.204 + 39.787/AG)·0.85.
- CASE(2,10) SP/PP: HGUESS = (1.88·SI − 7.178)·(1−EXP(−0.025·AG))^(0.001·SI+1.64).
**FINDAG(H)→SITAGE,SITHT** (findag.f): AGMAX=200, HTMAX=300; if H≥300 SITAGE=200+(H−300)/0.10,SITHT=H. Else AG=2
step +2: HGUESS=HTCALC(AG); if HGUESS≥1 and (|HGUESS−H|≤TOLER=2 or H<HGUESS)→SITAGE=AG,SITHT=HGUESS; else if the
curve flattens (INCRNG: OLDHG≠0 & ΔHGUESS≥0.05 then <0.05)→lock SITAGE=AG; AG>AGMAX→SITAGE=AGMAX,SITHT=H.
**htgf DRIVER**: redwood(12) = special LTHTG(D,SINDX,DG10=DG/bark,H)·0.5·HGBND (see audit). DEFAULT = FINDAG(H)→
SITAGE; if H≥HTMAX→HTG=0.1; if SITAGE≥AGMAX→POTHTG=0.10; else AGP05=SITAGE+5, HGUESS2=HTCALC(SITAGE+FINT?),
POTHTG=HGUESS_next − SITHT (the site-curve height rise over the period); HTG=SCALE·XHT·POTHTG·EXP(HTCON). (Read
htgf.f:240-290 for the exact POTHTG age step — FINT vs 5 — before coding.) NOTE: HD1-4 appear UNUSED by HTCALC/
FINDAG (the height curves are the CASE formulas above, not HD1-4) — verify HD1-4 aren't used elsewhere in htgf.
⇒ chunk-4 port: nc_htcalc(si,sp,ag) + nc_findag(h,sp,si) + height_growth!(::Klamath) [redwood LTHTG + default
FINDAG/POTHTG]. All formulas measured; validate HTG per-tree / aggregate vs FVSnc_clean.

## HT1/HT2 (blkdat.f:167-178) — MEASURED (regent height-DBH model, chunk 6; were CSV placeholders):
- ht1: 4.78737,4.74961,4.78737,4.80268,4.73881,4.89619,4.80420,4.66181,4.83642,4.23251,4.66181,5.3401
- ht2: -7.31698,-7.19103,-7.31698,-8.40657,-9.44913,-12.55873,-9.92422,-8.33117,-7.04795,-8.31711,-8.33117,-15.9354
Used in nc/regent.f: BX=HT2, AX=HT1 (IABFLG=1) or AA (IABFLG=0) → Wykoff HT-DBH DK=BX/(ln(HK-4.5)-AX)-1 for small-tree DG.
CSV updated (cols 17/18). Still placeholder in CSV: dbh_max, st_htadj, st_break, wykoff_ht2, sdi_max_default (ecocls),
volume cols, htdbh cols, varmrt — filled as chunks 6/8 land.
