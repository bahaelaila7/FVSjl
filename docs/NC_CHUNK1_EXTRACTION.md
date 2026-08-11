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

## STILL TO MEASURE (chunk-1 completion)
- bark: bratio.f BRKRAT(4,12) — the imap/DIB form; verify linear-encodable vs [0.80,0.99] clamp.
- dbh_max; small-tree st_* (regent.f DATA); ht1/ht2/wykoff_ht2 (htgf.f); mort_bkgd (morts.f); sdi_max_default
  (sdical.f: BAMAX=XMAX·0.5454154·PMSDIU ⇐ per-sp XMAX SDI); volume stump/top_dib/dbh_min/scf_*/bf_* (grinit.f/
  vollib); htdbh_* (htdbh.f); varmrt_varadj; is_sprouting; dg_resid_sd (SIGMAR).
- species_translation.csv: full FIA-species→NC-12 crosswalk (dgf.f OSPMAP / the FIA map, NOT just the 12).

## NEXT
Finish reading the above arrays → assemble data/klamath/species_coefficients.csv + species_translation.csv →
wire src/variants/klamath/species.jl + site_index.jl → validate SITE INDEX on nct01 vs FVSnc_clean (chunk-2 gate).
Then NC-specific dg_coefficients.jl + diameter_growth.jl (port the ISCT-section/DGHAH/2-set dgf logic) → chunk-3.
