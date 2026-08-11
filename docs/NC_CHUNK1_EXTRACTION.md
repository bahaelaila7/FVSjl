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

## NEXT
Finish reading the above arrays → assemble data/klamath/species_coefficients.csv + species_translation.csv →
wire src/variants/klamath/species.jl + site_index.jl → validate SITE INDEX on nct01 vs FVSnc_clean (chunk-2 gate).
Then NC-specific dg_coefficients.jl + diameter_growth.jl (port the ISCT-section/DGHAH/2-set dgf logic) → chunk-3.
