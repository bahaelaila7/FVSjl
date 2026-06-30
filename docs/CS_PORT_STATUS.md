# CS (Central States) port — status & live baselines

Active campaign (see `docs/CS_GOAL.md` + `docs/CS_VARIANT_PORT_SCOPE.md`). Oracle = live
`FVScs` relinked via `test/harness/cs_oracle.sh`; canonical stand `tests/FVScs/cst01.key`.

## Chunk 1 — variant infra (IN PROGRESS)
- [x] `CentralStates <: AbstractVariant` registered (`variant_code`="CS", `nspecies`=96,
      `htg_period`=10, Zeide SDI, RNG 55329) — `src/variants/centralstates/centralstates.jl`;
      wired into `variant_from_code` + exports; suite holds 5391/2 (no SN/NE regression).
- [x] `test/harness/cs_oracle.sh` (relink live FVScs from `bin/FVScs_buildDir/*.o`).
- [ ] CS species roster + coefficient CSVs under `data/centralstates/` (extract from
      `cs/blkdat.f` + CS coefficient files — the mechanical bulk; mirrors `data/northeast/`).
- [ ] `centralstates/species.jl` (init_blockdata!), `site_index.jl`, crown — drive cst01
      **cycle-0** stand columns bit-exact.

## Live FVScs ground truth — cst01 stand-1 cycle-0 (1990), the cycle-0 target
```
yr   age  TPA  BA  SDI CCF TopHt QMD | TCuFt MCuFt SCuFt BdFt
1990  60  536  77  160 169  63    5.1 | 1517  1300   497  2903
```
NB the stand is S248112 (same inventory as net01/snt01) but with **CS species codes**
(tree 1 = `WN` not NE's `JP`) and CS coefficients. Stand metrics match net01's cycle-0
(TPA/BA/SDI/TopHt/QMD) EXCEPT **CCF 169** (vs NE 176 — CS crown-width) and the **volume**
(CS equations). So cycle-0 reduces to: read the CS species roster, CS crown width (CCF),
and CS volume equations on the shared density/stat code.

## Chunk 1b — species data extraction (IN PROGRESS)
- [x] **Roster extracted** (real): 96 CS species from `cs/blkdat.f` `DATA JSP/FIAJSP/PLNJSP`
      → `data/centralstates/species_roster.csv` (WN=8, SM=43 — matches cst01.tre). Shared
      `species_translation.csv` copied (already carries the `target_cs` column).
- [ ] Build `data/centralstates/species_coefficients.csv` (the loader's required roster+coef file).
      Per-species coefficient sources located:
      - `cs/blkdat.f`: `BKRAT` (bark_intercept=0/slope=BKRAT), `SIGMAR` (dg_resid_sd), `HHTMAX`
        (estab_hht_max), `XMIN`, `HT1`/`HT2` (HT-DBH curve → htdbh_coeffs).
      - `cs/grinit.f` / `sitset`: `SDIMAX`/`SDIDEF` (sdi_max_default), site_group / site index.
      - `cs/crown.f` / `cratet`: crown_bcr1..4 (crown-ratio coefficients).
      - `cs/varmrt.f` / mort tables: mort_bkgd_*, varmrt_varadj; `regent` min diam; dbh_min.
- [ ] `data/centralstates/{crown_width_*,site_species,htdbh_coeffs,...}.csv` (mirror data/northeast/).
- [ ] `centralstates/species.jl` (init_blockdata!: 96 sp, Zeide SDI, RNG 55329) + site/crown hooks.
- NOTE for cycle-0 stand columns (TPA/BA/SDI/QMD/TopHt): these are geometric (DBH/HT/TPA) and need
  only the roster to parse the .tre + the shared Zeide SDI; **CCF** needs CS crown-width; **volume**
  needs CS NVEL eq ids. So the first testable milestone is TPA/BA/SDI/QMD/TopHt, then CCF, then volume.

### Extraction progress (real values, from CS Fortran — no placeholders)
- [x] roster → `species_roster.csv` (96, cs/blkdat.f JSP/FIAJSP/PLNJSP)
- [x] HT-DBH curve → `htdbh_coeffs.csv` (96, cs/blkdat.f HT1/HT2; same curve form as NE)
- [x] bark_slope + dg_resid_sd → `_blkdat_extract.csv` (96, cs/blkdat.f BKRAT/SIGMAR; bark_intercept=0)
- [ ] estab_min_ht / estab_hht_max → cs/blkdat.f XMIN / HHTMAX (parser needs the full 96 — got 86,
      a continuation-line/comment split; fix the DATA-block collector)
- [ ] site_group + sdi_max_default → cs/grinit.f / cs/sitset.f (the species→site-group + SDImax tables)
- [ ] crown_bcr1..4 → cs/crown.f (the TWIGS crown-ratio BCR coefficients)
- [ ] mort_bkgd_intercept/dbh + varmrt_varadj → cs/varmrt.f (check if CS uses group constants like NE)
- [ ] regent_min_diam, dbh_min → likely constants (NE: 0.1, 5.0) — confirm vs cs/regent.f
- Then ASSEMBLE `species_coefficients.csv` (roster + all the above) + crown_width_* + site_species,
  add `centralstates/species.jl`, and drive cst01 cycle-0 (geometric cols TPA/BA/SDI/QMD/TopHt first —
  those don't read the coefficients, so they validate honestly; then CCF via crown-width, then volume).

### Remaining species-data extraction — finalized roadmap (per-species, from CS Fortran)
Extracted so far (real): roster, bark(BKRAT), dg_resid_sd(SIGMAR), HT-DBH(HT1/HT2), crown
BCR1-4(cs/crown.f), varmrt_varadj(VARADJ), sdi_max_default(SDICON). Remaining, with exact source:
- **estab_min_ht / estab_hht_max** — `cs/blkdat.f` `XMIN` / `HHTMAX` (parser must collect the full 96;
  it stopped at 86 on a comment/continuation split — fix the DATA-block collector to span comment lines).
- **site_group + ASITE/BSITE site coefficients** — `cs/sitset.f` (`SDICON` already pulled there). Trace
  the species→site-group mapping; CS site index may be per-species ASITE/BSITE (not the NE SICOEF-group
  pointer) — confirm the model before reusing NE `site_index.jl`.
- **mort_bkgd_intercept / mort_bkgd_dbh** — the background-mortality rate (NE had 4 groups). Trace
  `cs/varmrt.f` / `cs/morts.f` (cs/morts.f is the background-kill driver; the rate coefficients/groups
  are the per-species source).
- **regent_min_diam** — `cs/regent.f` (NE had 5 distinct values). **dbh_min** — constant 5.0 (NE).

ASSEMBLY: combine the per-species CSVs into `data/centralstates/species_coefficients.csv` (the loader's
required file, same 20 columns as `data/northeast/`), then add `centralstates/{species.jl, site_index.jl,
crown_ratio.jl}` reusing the NE structure with CS coefficients (htgf/crown/varmrt are 88-96% NE), and
drive `cst01` **cycle-0**: TPA/BA/SDI/QMD/TopHt are geometric (validate first, honest — they don't read
coefficients), then CCF via crown BCR, then volume via CS NVEL eq ids. CAUTION (re-trace discipline): a
"reuses NE" label is a hint — verify each CS routine's source for a CS-specific coefficient/branch.

### Re-trace verdict — CS site model is CS-SPECIFIC (do NOT reuse NE site_index)
`cs/sitset.f` converts site index via a per-species LINEAR model: `SITEAR(j) = -ASITE(j)/BSITE(j) +
(1/BSITE(j))·SITEAR(isisp)`, default site species ISISP=47, default SI 65. This is DIFFERENT from NE's
SICOEF 28×28 group-matrix (`site_group`) model — confirming the scope's `sitset` 23-48% flag. So
`centralstates/site_index.jl` must be CS-specific (ASITE/BSITE), NOT a copy of `ne_site_index_setup!`.
Extracted ASITE/BSITE (96, real) → `site_coef.csv`. The doctrine's re-trace discipline caught this:
"reuses NE" was a wrong assumption for site. (Crown/htgf/varmrt remain NE-shaped per BCR/curve form.)

### CORRECTION — extraction parser under-collects some DATA blocks (must fix before assembly)
The ad-hoc DATA-block collector returns short on several arrays: XMIN/HHTMAX → 86 (not 96),
ASITE/BSITE → 93/95 (not 96). So `site_coef.csv` was NOT written this round (the verdict that CS site is
ASITE/BSITE-based stands — that was read from the equation, not the count). ROOT: the collector's
break heuristics (stop on a line starting `C`/`DATA`) cut blocks that span comment-interrupted
continuations or have unusual spacing. FIX before any further extraction: replace it with a proper
Fortran fixed-form continuation parser (column-6 continuation, strip `Cxxx` full-line comments only,
accumulate until the matching closing `/`), then re-extract XMIN/HHTMAX/ASITE/BSITE to the full 96 and
RE-VERIFY the already-written files (BKRAT/SIGMAR/HT1/HT2/BCR/VARADJ/SDICON) are complete 96-counts too.

### RESOLVED — proper fixed-form parser; all CS coefficient data re-extracted + verified
Replaced the ad-hoc collector with `tools/fortran_data_extract.py` (proper Fortran fixed-form:
col-6 continuation, full-line C/*/! comments, N*v repeat expansion). RE-EXTRACTED every CS array and
OVERWROTE the data files — the old parser had CORRUPTED HHTMAX/XMIN/ASITE/BSITE (dropped the leading
value → mis-ordered; HHTMAX[1] read 27 instead of the raw 16). All now 96-count and spot-checked vs the
raw Fortran (HHTMAX[1]=16, ASITE[3]=-5.1489, BCR1[1]=4.0862, SDICON[1]=354, HT1[1]=4.4718 — all ✓).
Real CS species data complete for cycle-0: roster, bark, dg_resid_sd, estab(hhtmax/xmin), HT-DBH(HT1/HT2),
crown BCR1-4, site(ASITE/BSITE + sdi_max), varadj. Remaining (cycle-1+ only): mort_bkgd, regent_min_diam,
dbh_min(=5.0). NEXT: assemble species_coefficients.csv + centralstates/{species.jl, site_index.jl (CS
ASITE/BSITE), crown_ratio.jl (NE-shape + CS BCR)} → cst01 cycle-0.

### CHUNK 1 COMPLETE — cst01 cycle-0 ALL 6 stand columns BIT-EXACT vs live FVScs
TPA=536, BA=77, SDI=160, CCF=169, TopHt=63, QMD=5.1 (test/integration/test_cst01.jl;
suite 5399/2, no SN/NE regression). What landed:
- CS variant hooks: species.jl (blkdat init), site_index.jl (ASITE/BSITE SITSET + CS
  FORKOD lat/long/elev by IFOR). CS HT-DBH reuses the Southern Curtis-Arney+Wykoff dub
  (cs/htdbh.f ≡ sn/htdbh.f logic); CS crown reuses NE's TWIGS method (negated BCR4 in data).
- CS crown-width: cs/cwcalc.f is BYTE-IDENTICAL to ne/cwcalc.f (selects on the 2-char alpha,
  not species index) ⇒ CS reuses NE's crown_width_{equations,species}.csv verbatim. 92/96 CS
  alphas covered; the 4 gaps (3 blank pads + 'BW') are inert for cst01.
- SHARED BUG SURFACED (doctrine #3): the .sum ÷GROSPC scale-back (disply.f) — internal stats
  are per-stockable-acre, the .sum reports per-gross = stockable·GROSPC. FVSjl's summary_row
  ALREADY divides by gross_space, so it was correct; the apparent 10% gap was me comparing the
  internal per-stockable value (590) against the per-gross .sum target (536). GROSPC<1 is the
  first stand to exercise this (SN/NE test stands are GROSPC=1.0). Confirmed via the live .out
  'BASED ON STOCKABLE AREA' table = FVSjl's internal 590/85/63.4.
- Hopkins index needs lat/long ⇒ CS FORKOD (forest 905 → IFOR 1 → 37.95/91.77/10). Without it
  CCF read 143; with it, 169 exact.

NEXT (chunk 2): CS volume — wire CS NVEL equation ids into the R9 Clark + R9LOGS path so the
cyc0 .sum volume columns (Tcuft/Mcuft/Bdft 1517/1300/497...) and forest-type (503) come in.
Then chunk 3: cs/dgf.f (the one genuinely-new SN-family CS DG model) for cycle-1.
