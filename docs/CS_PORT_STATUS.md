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

### CHUNK 2 DONE — cst01 cycle-0 VOLUME columns ALL BIT-EXACT vs live FVScs
Tcuft=1517, Mcuft=1300, Scuft=497, Bdft=2903 (suite 5403/2). CS volume rides the eastern
NVEL Region-9 Clark cubic + R9LOGS board path (shared with NE) — cs/vols.f is byte-identical
to ne/vols.f; the R9 Clark cubic keys on FIA code (national). The CS-specific piece is the
merch standards: ported `_cs_merch` (cs/sitset.f:130-227 — softwoods=sp 1-7, eastern redcedar
sp 1 lower mins, IFOR rules; bf-equal). compute_volumes! + compute_volumes_ne! now route CS.
- SHARED-PATH FIX (caught by the doctrine's per-tree trace): the two BROKEN-TOP trees (d8.0,
  d10.4; live TRC HT 56/49) read ~4 cuft/acre low because CS never populated calib.bark_a/b
  ⇒ CFTOPK's bark_ratio fell to the clamp 0.80 instead of BKRAT. Added minimal cs_dgcons!
  (bark copy, mirrors ne_dgcons!) wired into setup_growth!; the full cs/dgf.f site-constant +
  serial-correlation port is chunk 3. With the bark copy, both broken-top trees match exactly
  (d8.0→8.6, d10.4→15.0) and all 4 volume columns close.
- The canonical cyc0 path is notre! → setup_growth! → compute_volumes! (setup_growth! does the
  CRATET dub + cs_dgcons!); the earlier manual dub+notre left bark=0.
- KNOWN cyc0 GAP (follow-up, not a stand/volume column): the trailing .sum forest-type field
  reads 999 vs live 503 — CS FORTYP classification (cs/fortyp* / forest_type.jl) is unported.

NEXT — chunk 3: cs/dgf.f (the one genuinely-new SN-family CS diameter-growth model: ln(DDS)
from DBH/site/crown/BA-percentile/QMD, ≥5") + the full cs_dgcons! site constant ⇒ cycle-1 DG.
Then the CS FORTYP forest-type field.

### cyc0 FORTYP done — full cst01 inventory row now 100% bit-exact
Forest type 503 (W.OAK-R.OAK-HICKORY) + size/stock class 2/2 match live. cs/fortyp.f and
cs/stkval.f are BYTE-IDENTICAL to ne's, and the stocking map is FIA-keyed (national) ⇒ CS
reuses NE's stocking_coeffs.csv + fia_stocking_map.csv verbatim. The only cyc0 .sum fields
not matching are the growth columns (period/accretion/mortality), which inherently need the
cycle 0→1 DG projection (cs/dgf.f, chunk 3). test_cst01.jl: 15/15. Cycle-0 is COMPLETE.

### CHUNK 3 IN PROGRESS — cs/dgf.f MODEL ported + validated bit-exact per-term; calibration COR open
The CS diameter-growth MODEL is done and PROVEN correct against a live debug stamp:
- Extracted the 12 DDS coefficient arrays (INTERC/VDBHC/DBHC/DBH2C/RDBHC/RDBHSQC/CRWNC/CRSQC/
  SBAC/BALC/SITEC + OBSERV) → data/centralstates/dg_coeffs.csv. Wrote `dgf!(s, ::CentralStates)`
  (src/variants/centralstates/diameter_growth.jl): the SN-family ln(DDS) regression + the QMD≥5/CR
  species caps + BAGE5/QMDGE5/BAL competition + the OB→IB bark conversion (WK2 = ln(IB DDS)).
- Updated cs_dgcons!: DGCON=0, ATTEN=OBSERV(ISPC), SMCON=0, bark copy. Wired diameter_growth!
  (the generic AbstractVariant application already dispatches dgf!) + calibration into setup_growth!.
- LIVE DEBUG STAMP (DEBUG/DGF keyword on cst01) confirms cs_dgf! is BIT-EXACT per-term: for the HI
  (ISPC=19) tree the live D1..D11 = (-0.647, 0, 1.084, -0.055, 0.250, -0.077, 2.005, -1.181, -0.098,
  0, 0.446) match jl's coefficients×inputs exactly; CR/PCT/SITE/BAGE5/QMDGE5/BAL all match.
- OPEN — the DG CALIBRATION (COR): cst01 trees carry MEASURED past DBH growth (read into diam_growth:
  WN d11.5→1.0, HI d6.5→2.3, …). FVS calibrates the per-species COR so the model reproduces the
  measured growth (hence live "CURR DIAM INCR" = the measured input, NOT the bare model). jl's
  calibrate_diameter_growth! currently yields COR only for sp47 (WO, goal 0.267) and 0 for sp8/19/43,
  so uncalibrated species predict ~0.7 where measured/live is ~2.3. NEXT: trace why the shared SN
  calibration doesn't accumulate COR for the other measured CS species (suspects: the measured-DG IDG
  conversion / GROWTH code for the CS .tre growth field; per-species obs threshold fnmin=5; or a
  cs/dgdriv.f specific the shared path misses). Validate per-species COR vs a live dgdriv debug stamp.
- Suite 5406/2 (no SN/NE regression, no cyc1 test added until the COR is bit-exact).

#### Refined COR diagnosis (fnmin / OLDRN residual)
Per-species valid measured-DG counts for cst01: sp47=5 (COR FIRES), sp8=3, sp19=4, sp43=4,
sp60=3 (all <5 ⇒ COR=0). The shared calibration's fnmin=5 species threshold exactly explains
which species calibrate. But live reproduces the measured growth for the <5 species too —
because a MEASURED tree grows by (≈) its measured increment in cycle 1 via its per-tree OLDRN
serial-correlation RESIDUAL (= ln(measured DDS) − ln(predicted DDS)), which is seeded even for
species below the COR threshold. The jl divergence is therefore in the OLDRN seeding for
measured trees in UNCALIBRATED species (snt01/net01 don't exercise this: their measured species
clear fnmin). NEXT: confirm vs a live dgdriv OLDRN stamp, then ensure measured trees in <fnmin
species carry their residual (not a BACHLO draw) so cycle-1 DG = measured input.

#### COR diagnosis CORRECTED (scale is fine; OLDRN residual is the locus)
Live DGDRIV calibration stamp (DEBUG) gives, for the HI ISPC=19 measured trees:
I=4 OBS.DG=0.6 TERM=8.4564, I=5 OBS.DG=0.7 TERM=9.926, I=10 OBS.DG=1.2 TERM=16.862.
TERM = DG·(2·BARK·WK3+DG)·SCALE. Back-solving I=4 with WK3≈7.25 (sensible backdated dbh of the
d7.9 tree) ⇒ SCALE=1.0. And jl's growth_fint=5 (the MEASUREMENT period, not the 10-yr cycle), so
jl's `5/dfint = 5/5 = 1.0` already equals live's YR/FINT = 10/10 = 1.0. ⇒ the SCALE is CORRECT;
the earlier "halved-scale" theory was WRONG (forcing 10/dfint=2.0 overshot the calibrated WO from
2.09→3.92). REVERTED — no change to simulate.jl.
True locus: OLDRN = RESLOG = ln(TERM) − WK2(backdated) is ~2× too small in jl (HI d6.5 OLDRN 0.587
vs the ~1.3 needed to reproduce the measured 2.3 via DGSCOR FRM). cs/dgdriv.f:412-422: a measured
tree gets OLDRN(I)=RESLOG when DGSD≥1.0 (met for HI's 4-5 trees); DGSCOR then folds OLDRN into FRM,
and DG=√(DSQ+DDS·FRM)−D. NEXT: instrument jl's calibrate_diameter_growth! to dump per-tree TERM +
backdated WK2 and diff vs live (I=4 HI TERM=8.4564) — that isolates whether the gap is the measured
DDS (TERM/WK3 backdating) or the backdated prediction WK2. The CS bark-in-dgf conversion (cs_dgf!
bark-converts WK2, unlike SN) is the prime suspect for a WK2 mismatch the shared calibration assumes.

#### DEFINITIVE diagnosis — calibration RESLOG goes the WRONG DIRECTION for CS (bark bookkeeping)
Per-tree OLDRN (=RESLOG=ln(TERM)−WK2_backdated), HI sp19, jl vs live (live = DGDRIV CUM.DEV stamp):
  measDG 0.6 (d7.9): jl +0.146  live +0.069
  measDG 0.7 (d8.0): jl −0.170  live +0.300
  measDG 1.2 (d8.2): jl −0.949  live +0.893
Live RESLOG INCREASES with measured DG (a faster-grown tree ⇒ larger positive residual = model
under-predicts). jl's DECREASES (even flips sign) ⇒ jl's TERM/WK3/WK2 bark bookkeeping is wrong for CS.
Analytic check (both-sides): cs/dgf.f WK2 = ln((D·BR+DIAGRI)²−(D·BR)²) = ln(inside-bark DDS). For the
matching residual, cs/dgdriv.f converts DG→DG_ib (line 329 DG=DG·BRATIO) THEN TERM=DG_ib·(2·BR·WK3+DG_ib)
= BR²·(2·WK3·DG+DG²) = inside-bark measured DDS — same bark space as WK2. So the shared FVSjl calibration
must, for CS: (1) convert measured DG to inside-bark via BKRAT, (2) TERM=DG_ib·(2·BR·WK3+DG_ib) with WK3 =
OUTSIDE backdated dbh, BR=BKRAT. The wrong-direction RESLOG means one of these is off (double bark, wrong
WK3 backdating basis, or the SN affine bark vs CS constant BKRAT). NEXT: instrument calibrate_diameter_
growth! to dump per-tree (WK3, TERM, WK2) and match live (I=4 HI: WK3≈7.255, TERM=8.4564, WK2=2.066,
RESLOG=0.069); fix the CS bark/backdate path. (Alternative path: proceed to chunk 4 height/mortality to
get the full cycle-1 .sum, then revisit — but the wrong-sign RESLOG will bias cycle-1 DG, so fix first.)

#### ★ CORRECTION — calibration + model are BIT-EXACT; the prior "2x wrong" was THREE misreads
Instrumented cs_dgf! (during calibration) and the calibration TERM directly. Findings:
1. cs_dgf! backdated WK2 matches live DGF stamp EXACTLY: i=3 d4.027 DDS 1.5808, i=4 d7.255 DDS
   2.0657, i=5 d7.247 DDS 2.0641 (== live 1.58079/2.066/2.064).
2. The calibration RESLOG matches live BIT-EXACT. jl per-tree RESLOG 0.0692/0.231/0.5931 ACCUMULATE
   to 0.0692/0.300/0.893 = live's CUM.DEV column. The earlier "wrong-direction RESLOG" was comparing
   jl's PER-TREE residual against live's CUMULATIVE-deviation column — a misread, not a bug.
3. The t.old_random values (0.146/−0.17/−0.95) are the BACHLO draws for the UNCALIBRATED sp19 (4<5
   trees), exactly as FVS does it (OLDRN=RESLOG only when DGSD≥1; else BACHLO) — not residuals.
4. The live "CURR DIAM INCR" = 2.3 I'd chased is the MEASURED INPUT echoed in the inventory treelist,
   NOT the cycle-1 projection. The live 2000 treelist shows the d6.5 HI tree TRIPLED: central DBH 7.2
   (projected DG 0.7), upper 8.2, lower 6.8. With proper tripling=true, jl's central DG = 0.654 vs
   live 0.7 (~7%); my earlier 1.26 used tripling=FALSE, which wrongly applies OLDRN to the lone record.
VERDICT: the CS DG MODEL + CALIBRATION are bit-exact. The only residual is a ~7% gap on the tripled
central DG (0.654 vs 0.7) = the DGSCOR serial-correlation / tripling-split detail (shared code; the
CS bark-converted WK2 may interact with the central-vs-triple OLDRN partition). NEXT: trace DGSCOR for
the central record vs a live growth-mode stamp, OR proceed to chunk 4 (height/mortality/crown) and
validate the full cycle-1 .sum (2000: TPA 518/BA 99/SDI 196/CCF 202/TopHt 68/QMD 5.9), where the small
DG residual surfaces aggregately. The earlier "halved-scale", "wrong-coefficient", "wrong-direction
RESLOG", and "2x-low DG" diagnoses are all SUPERSEDED — model+calibration verified correct.

#### Honest refinement — residual is the DGSCOR serial-correlation (not model/calibration)
Don't overclaim "DG bit-exact": the MODEL + CALIBRATION are bit-exact, but the projected cycle-1
DG still has a real residual in the DGSCOR serial-correlation/tripling step. Live's d6.5 HI tree
triples to 8.2/7.2/6.8 (TPA 3.82/9.17/2.29), TPA-weighted mean DG = 0.89 — ABOVE the base model
0.697 (the multiplicative OLDRN spread boosts the mean). jl's central/base DG = 0.654 (slightly
BELOW base 0.697 — a stronger FRM bias-correction than live's central 0.7≈base). So jl's mean will
land below live's 0.89. Locus = the DGSCOR FRM (SSIGMA/RHO + the central-vs-triple OLDRN partition)
for CS: the YR=10 AUTCOR variance multiplier (VMLT) and SIGMAR feed it, and the CS bark-converted
WK2 may interact. This is shared code (SN/NE bit-exact at YR=5/10), so suspect the CS YR=10 VMLT or
the uncalibrated-species BACHLO/RNG alignment. NEXT (two options): (a) a live growth-mode DGSCOR
stamp (DG/FRM per tree) to pin the central FRM + tripled spread, or (b) implement chunk 4 (height/
crown/mortality) and validate the full cycle-1 .sum (2000: TPA 518/BA 99/SDI 196/CCF 202/TopHt 68/
QMD 5.9) — the aggregate where the DG residual surfaces. Recommend (b): the .sum is the real target
and exercises DG+HTG+mort together; isolated per-tree DGSCOR tracing is high-effort/low-signal.

### ★ CHUNK 4 — full grow_cycle! RUNS for CS; cycle-1 .sum within ~4% of live
Ported the remaining growth-spine hooks; grow_cycle!(CentralStates) now runs end-to-end.
CS cyc1 (2000): TPA 512/BA 103/SDI 202/CCF 207/TopHt 68/QMD 6.1 vs live 518/99/196/202/68/5.9.
TopHt is BIT-EXACT (height growth correct); the rest is within ~3-4%.
- small_tree_growth!(::CentralStates): cs/regent.f — NE shape, XMIN=3.0, cs_balmod, MAPCS curve,
  DIAM floor = htdbh_db (budwidth/SNDBAL).
- mortality!: the generic AbstractVariant method serves CS. Supplied the background-mortality
  coeffs: RI=0.5/(1+exp(PMSC[g]+PMD[g]·D)) with g=IMAPCS[sp] (4 mortality groups, morts.f
  byte-identical CS/NE) → mort_bkgd_coeffs.csv. VARMRT efficiency: cs/varmrt.f ≡ ne/varmrt.f
  (only VARADJ DATA differs, CS supplies varmrt_varadj) ⇒ _varmrt_efftr! now dispatches
  Union{Northeast,CentralStates}.
- REMAINING cyc1 gaps (all small, to close for bit-exact): BA 103 vs 99 + QMD 6.1 vs 5.9 (~4%
  HIGH) = the DG over-growing — the DGSCOR serial-correlation/tripling residual (jl's tripled-mean
  DG runs above live's; the multiplicative OLDRN-spread for uncalibrated-species BACHLO draws is the
  prime suspect — validate the per-tree cyc1 DG/BACHLO RNG alignment vs a live growth-mode stamp).
  TPA 512 vs 518 (~1% over-kill) = mortality, likely downstream of the BA/DG diff (SDI-density kill).
  Suite 5406/2, no SN/NE regression.
NEXT: pin the DGSCOR/BACHLO cyc1-DG residual (now validatable aggregately via the .sum); then the
TopHt-exact spine + the bit-exact DG should bring BA/SDI/CCF/QMD/TPA to live. Then multi-cycle.

#### cyc1 residual refinement — per-tree DG matches; the pattern points to MORTALITY
Dumped jl's post-grow_cycle! tripled records for the d6.5 HI tree: jl 7.2/8.2/6.78 == live
7.2/8.2/6.8 — the DG/tripling DBH is RIGHT for this tree (so the DGSCOR worry may be smaller than
feared). The cyc1 .sum signature — QMD 6.1>5.9, BA 103>99, TPA 512<518 (FEWER trees but LARGER) —
is the classic "mortality over-kills the SMALL/suppressed trees" pattern: killing small stems
raises mean DBH/QMD and leaves BA high. So the prime cyc1 suspect is the MORTALITY distribution
(SDI-density kill + varmrt efficiency), not (only) DG. NEXT: compare jl's per-tree mortality (which
records lose TPA) vs the live 2000 treelist MORTAL-PER-ACRE column; check the CS SDImax (sp_sdi_def
from SDICON) and the density-mortality DR/D10 basis — a too-low CS SDImax over-thins. (DG per-tree
DBH validated against the live treelist; revisit DGSCOR only if mortality alone doesn't close it.)

#### cyc1 closure lead — jl over-kills ~6 TPA (24 vs live 18); per-tree mortality is the next probe
1990→2000 mortality: jl kills ~24 TPA (536→512), live ~18 (536→518). Live's per-tree MORTAL-PER-ACRE
(from the 2000 treelist) concentrates in a few records (3.4, 1.6, 1.4, 1.1) + many tiny; jl's total is
higher. The CS mortality path is: background RI=0.5/(1+exp(PMSC[g]+PMD[g]·D)) (g=IMAPCS, supplied) +
the SDI-density kill (shared SN/NE-validated code; CS SDImax from SDICON, sp_sdi_def). jl's stand
sp_sdi_def max = 726 (one species); the composition-weighted stand SDImax that the density kill uses
needs checking — the cyc0 SDI 160 vs SDImax drives how hard density mortality fires. NEXT (better
tooling than awk-on-treelist): (a) a live MORTS DEBUG stamp (per-tree RI/RN/RIP) to verify the CS
background RI + density rate per tree, and (b) dump jl's per-tree mort_pa vs the live treelist MORTAL
column (use FVSOut.db FVS_TreeList if a TREELIST DB run is set up, not the fixed-width .trl). The DG
per-tree DBH is validated (d6.5 HI triples match live), so isolate mortality first. This is the last
gap to cyc1 bit-exact; then multi-cycle + the .sum.save regression test.

#### cyc1 mortality — density PARAMETERS all match live (hypotheses ruled out)
Live MORTS DEBUG stamp + jl dumps confirm the CS density-mortality inputs are CORRECT:
- stand SDImax (SDICAL/CLMAXDEN, composition-weighted): jl 345.45 == live 345.449 (EXACT).
- per-species sp_sdi_def (SDICON): jl sp8=283/sp19=302/sp43=371/sp47=361 == live SDI MAX list.
- stand SDI cyc0: jl 175.57/gross = 159.6 ≈ live 160.
- ssigma (DGSCOR BACHLO spread) = SIGMAR for uncalibrated species (faithful; vardg back-derived so
  ssigma=SIGMAR exactly, vmlt=29.40 from autcor(10,10)).
- fint=10 cycle handling: NE-validated (NE also YR=10, bit-exact), so not the over-kill source.
RULED OUT: SDImax / sp_sdi_def / SDI / ssigma / fint-compounding. jl kills 23.8 TPA vs live ~18.
Remaining = a subtle coupled effect: either a small AGGREGATE DG over-prediction (the DGSCOR upper-
triple BACHLO mean for uncalibrated species pushing grown SDI up → more density kill — even though the
sampled d6.5 tree's DBH matched) OR the per-tree mortality DISTRIBUTION. NEXT (needs better tooling
than the stock DEBUG/MORTS, which dumps SDICAL but not per-tree XKILL): add a custom FVS morts.f stamp
dumping per-tree (D, RI, RN, RIP, XKILL) and diff vs jl's per-record mort_pa; and instrument jl's
grow_cycle! to expose the PRE-mortality grown BA/SDI (compare to the live MORTS 'OLDBA/RELDEN' = 84.35/
169.71) to settle DG-aggregate vs mortality-distribution. This is the final cyc1 gap — multi-session
per-tree work, as the SN/NE campaigns were.
