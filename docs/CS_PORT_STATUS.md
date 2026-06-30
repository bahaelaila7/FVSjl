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
