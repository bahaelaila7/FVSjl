# OC (Oregon Coast) variant port — audit (task #208, Stream 2)

**Status: chunk-0 FOUNDATION landed + oracle relinked + ORGANON coupling fully measured.** The
ORGANON growth subsystem itself is UNPORTED (that is the cost driver — see
`docs/OC_ORGANON_PORT_PLAN.md`). This is an honest foundation, not a rushed partial subsystem.

## What OC is (the load-bearing finding)

OC is **not** a Wykoff-DDS variant like the rest of the western cluster. Its entire per-cycle
growth engine is **ORGANON edition SWO (Southwest Oregon)** — the `organon/` (13,040 lines) +
`vorganon/` (1,331 lines) subsystem, which FVSjl does not have. MEASURED, not inferred: a live
`FVSoc_clean` run with `DEBUG DGF` produced `ocmin.out` whose banner reads
`VERSION FS2026.1 ORGANON SWO`, whose DGF path dumps the FVS tree list being loaded into ORGANON
(`I,PTNO,SPECIES,DBH1,HT1OR,CR1,EXPAN1,...`), emits `CRATET ORGANON ERROR CODE, CYCLE= 1 IERROR= 0`,
and returns growth through the ORGANON `ACALIB/TMPCAL` calibration arrays (e.g.
`TMPCAL(1,1)=0.789290 TMPCAL(2,1)=0.624002`). **Consequence for the tasked "validate FVS-native
large-tree DGF cyc0" step: for the 18 ORGANON-eligible species there is no FVS-native DGF to
validate — DGF *is* ORGANON.** The FVS-native `oc/dgf.f` Wykoff path only runs for the other 32
species and only when a stand lacks all "big-6" trees.

## Foundation delivered (chunk 0)

- `src/variants/oregoncoast/oregoncoast.jl` — `OregonCoast <: AbstractVariant` singleton;
  `variant_code="OC"`, `nspecies=50`, `htg_period=5f0` (ORGANON's native 5-yr step).
- `src/variants/oregoncoast/species.jl` — `init_blockdata!` with MEASURED grinit constants
  (`oc/grinit.f`): 5-yr cycle (FINT=5), Stage SDI (`LZEIDE=.FALSE.`), `DGSD=0` (no DG serial-corr),
  `LHTDRG=.FALSE.` all species, seed 55329.
- `data/oregoncoast/species_coefficients.csv` + `species_translation.csv` — the 50 OC species
  (alpha/FIA/PLANTS from `oc/blkdat.f` JSP/FIAJSP/PLNJSP). Only the code columns — DG/bark/HTG/CR
  coefficients live in the (unported) ORGANON engine, so none are fabricated.
- Registration: `variant_from_code("OC"|"OREGONCOAST"|"OREGON COAST")` in
  `src/variants/variant.jl`; includes in `src/FVSjl.jl`.
- **Verified:** `using FVSjl` precompiles; `variant_from_code("oregon coast") → OregonCoast()`;
  coefficients load (50 species, DF=7/202, RW=50/211, arrays padded to global MAXSP=108 as usual).
- Growth hooks intentionally UNIMPLEMENTED — dispatching a growth hook on `OregonCoast()` errors
  loudly (doctrine #5), exactly as BC/UT/CI/NC chunk-0 scaffolds do.

## Chunk C1 delivered (FVS↔ORGANON boundary marshalling) — VALIDATED BIT-EXACT

- `src/variants/oregoncoast/organon_interface.jl` — ports **three** pieces (no growth):
  - **`OC_OSPMAP`** (`oc/orgspc.f` `DATA OSPMAP`): the 50-element FVS-index → ORGANON-FIA-code map,
    plus `OC_ORGANON_VALID` (the 18 valid-ORGANON species) and `OC_ORGANON_BIG6` (`{2,4,7,16,18}`).
  - **Per-tree eligibility `IORG`** + **big-6 stand gate** (`oc/dgdriv.f:217-261`): `iorg=1` iff
    `HT>4.5 ∧ DBH≥0.1 ∧ species∈the 18`; `nbig6 = Σ(gate ∧ species∈big-6)`; if `nbig6==0` all
    `iorg=0`, `mortexp=0`, `runs=false` (revert to FVS-native, FVS `GO TO 261`).
  - **`/ORGANON/` input-buffer fill** (`oc/dgdriv.f:268-285`): `OrganonBuffer` (SPECIES, DBH1,
    HT1OR, CR1, SCR1B, EXPAN1, MGEXP, USER + IORG/TREENO/PTNO) built from `StandState.trees`.
- Wired as the OC growth entry: `diameter_growth!(s, ::OregonCoast)` calls `build_organon_buffer!(s)`
  then errors loudly (ORGANON growth = C3-C6, unported). A minimal `site_setup!(s, ::OregonCoast)`
  (SITSET deferred to C2 — cyc0-inert for marshalling) was added so OC stands initialize.
- **MEASURED vs the live oracle** (`FVSoc_clean`, scoped `DEBUG / CRATET DGDRIV` — see below — stand
  S248112/ocmin, 27 records). jl loaded via `each_stand`+`notre!` (no growth) then
  `build_organon_buffer!`; per-tree diff vs the dgdriv `FOR EXECUTE` dump + CRATET `SPECIES` dump:

  | field | result |
  |---|---|
  | PTNO | **27/27 exact** |
  | DBH1 | **bit-exact** (max \|Δ\|=9.5e-8, f32 print precision) |
  | CR1 | **bit-exact** (max \|Δ\|=1.6e-10) |
  | EXPAN1 (=PROB) | **bit-exact** (max \|Δ\|=4.9e-8) |
  | USER (=ISPECL) | **27/27 exact** (all 0) |
  | HT1OR | **exact on all 26 measured-height trees**; 1 diff = tree 20 (raw HT=0→jl floors 4.6, oracle ORGANON-dubbed 53.32 — that dubbing is **C2**) |
  | SPECIES (ORGSPC) | **27/27 bit-exact once the FVS index is correct** (see crosswalk note) |
  | gate | **nbig6=17, nvalid=17, runs=RUN — exact** |
  | IORG (per-tree) | **bit-exact** (jl raw-height valid set ⊆ oracle dubbed set; both =17 ⇒ identical) |

- **One residual, and it is NOT a C1 bug — it is the C0 species-crosswalk TODO.** With jl's current
  species table, the alpha codes `WF` and `ES` (which appear in the .tre but are not OC primary
  species) fall through to `other_species=49 (OH)`, so 10 trees got ORGANON code 492 instead of 017.
  The live oracle (`base/intree.f` NOTE): **`WF` SET TO `GF` (index 4)** and **`ES` SET TO `BR`
  (index 22)** via the shared `vie/spctrn.f` 442-row ASPT table (OC has no `oc/spctrn.f`).
  `OSPMAP[4]=OSPMAP[22]=017`, so the ORGSPC map is already correct — patching those 10 indices to
  the oracle-confirmed values makes SPECIES **0/27 mismatches**. **Fix belongs to the species chunk
  (C0): extend `data/oregoncoast/species_translation.csv` with the `vie/spctrn.f` ASPT alt-code rows.**

### DEBUG-seam correction (measurement infra)

- The prior `ocmin.key` used a bare `DEBUG` keyword ⇒ `initre.f:1076` `DBALL` = debug **all**
  routines, which (a) made `DGF` an "invalid keyword" and (b) triggered the DEBUG-gated volume
  crash. **The `fvsvol.f:530` SIGSEGV is a DEBUG-only crash** (a `WRITE` under `IF(DEBUG)`), not a
  growth blocker. Scoping DEBUG to the two growth routines — `DEBUG` with a non-blank field-2 then a
  `CRATET DGDRIV` supplemental record (`initre.f:1072`/`dbprse.f`) — runs cleanly to completion and
  emits **both** the CRATET setup dump AND the dgdriv `FOR EXECUTE` per-tree dump (the true C1
  target). The initial-inventory volume runs before DGDRIV, so a full-DEBUG run never reached the
  dgdriv dump; the scoped key is the correct recipe for C1/C3-C6 A/B.

**C1 verdict: bit-exact vs the live oracle on every marshalled field it owns** (species map, gate,
buffer fill); the lone SPECIES residual is a pre-flagged C0 crosswalk gap, and the lone HT1OR
residual is C2 dubbing — neither is a C1 defect. **Next: C2 (setup calibration `prepare.f`/`start2.f`).**

## Chunk C2 delivered (ORGANON PREPARE setup calibration) — VALIDATED BIT-EXACT

`src/variants/oregoncoast/organon_setup.jl` ports the ORGANON **PREPARE** setup path for edition
**SWO (VERSION=1)** — the deterministic (`DGSD=0`) computation of the `ACALIB`/`TMPCAL` height/
crown/diameter calibration multipliers ORGANON growth (C3-C6) consumes, plus the ORGANON HT/CR
imputation (dubbing) for valid ORGANON trees with missing height/crown:

- **`organon_prepare_swo`** = `organon/prepare.f` PREPARE + `oc/cratet.f:129-401` driver. Ported
  sub-routines: `EDIT` (species groups via `SPGROUP_EDIT`/`SCODE1`, missing-HT/CR flags, RAD
  detection, DF↔PP SI conversion), `HDCALIB` (H-D calibration ratio → `TMPCAL(1,*)`), `PRDHT`
  (missing-height imputation), `CRCALIB` (crown-ratio calibration ratio → `TMPCAL(2,*)`), `PRDCR`
  (missing-CR imputation), plus `SPMIX`, `DFORTY`, `HS_H40`, `HD40_SWO`, `A_HD_SWO`, `A_HCB_SWO`,
  `CALTST`, `GET_CCFL_EDIT` (`organon/start2.f`); `SSUM`/`OLDGROWTH`/`GET_BAL` (`organon/diamcal.f`);
  `MCW_SWO` (`organon/crngrow.f`); and the `oc/cratet.f:393-401` `TMPCAL`→`ACALIB` load. All REAL*4
  math routed through the gfortran-identical `fexp`/`flog`/`fpow` (doctrine #8).
- **`site_setup!(::OregonCoast)`** now ports the one ORGANON-load-bearing SITSET piece
  (`oc/sitset.f:181-189`): the DF(7)↔PP(18) site-index conversion feeding ORGANON `SITE_1`/`SITE_2`.
- **C1 tree-20 residual fixed**: `build_organon_buffer!` floors `HT1OR` to 4.6 only when `HT>0`
  (`oc/cratet.f:234`), so a missing height (`HT==0`) reaches PREPARE as `0.0` and flags `MISSHT`.
- **Species-crosswalk fix (C0 TODO, now closed)**: `data/oregoncoast/species_translation.csv` gains
  the two ASPT alt-code rows `WF→GF` and `ES→BR` (verified from `vie/spctrn.f` column 20 = OC).
  `resolve_species` now maps `WF`→FVS 4 (GF) and `ES`→FVS 22 (BR); both `ORGSPC`=017 ⇒ **SPECIES
  27/27 exact**, and IORG is correct (GF valid+big6, BR non-valid → FVS-native).

**MEASURED vs the live oracle** (`FVSoc_clean`, scoped `DEBUG 1 / CRATET`, stand S248112/ocmin,
27 records; oracle `ORGANON TMPCAL(k,grp)` dump). The Julia PREPARE was fed the exact `/ORGANON/`
buffer the Fortran received (`SPECIES,DBH1,HT1OR,CR1,EXPAN1` — C1-validated), plus `SITE_1=92.0`,
`SITE_2=86.5528641`, `MSDI=815`, `STAGE=60`, `BHAGE=54`, `EVEN=.TRUE.`, `NPTS=11`:

| calibration entry | jl | oracle | \|Δ\| |
|---|---|---|---|
| `TMPCAL(1,1)` DF height | 0.7892899 | 0.789290 | 1.2e-7 |
| `TMPCAL(2,1)` DF crown  | 0.6240016 | 0.624002 | 4.2e-7 |
| `TMPCAL(1,2)` GF height | 0.7310810 | 0.731081 | 6.0e-8 |
| `TMPCAL(2,4)` SP crown  | 0.5000000 | 0.500000 | 0.0 |
| all other 50 of 54 `TMPCAL(k,grp)` | 1.0 | 1.0 | 0.0 |
| site conv `SITEAR(18)` | 86.552864 | 86.5528641 | f32 print |
| `resolve_species` WF/ES | GF(4)/BR(22) | GF/BR | exact |

Max \|Δ\| over all 54 `TMPCAL` entries = **4.2e-7** — pure Float32 print-precision (the oracle
prints `F9.6`). **C2 CALIBRATION: BIT-EXACT.** ORGANON is deterministic (`DGSD=0`), so this is the
hard bar, not a straddle. The final `ACALIB` after the cratet load = `[0.789290, 0.731081, 1, 1]`
(HT), `[0.624002, 1, 1, 0.5]` (CR), `[1,1,1,1]` (DG).

### Measured correction to the C1 "ORGANON-dubbed 53.32" hypothesis

The C1 audit assumed ocmin **tree-20** (raw `HT=0`) was dubbed to 53.32 by ORGANON PREPARE. The
scoped-`DEBUG CRATET` dump shows otherwise: tree-20 is a **lodgepole-pine (`LP`, FVS 12) record**,
which is a **non-valid ORGANON species** (`IORG=0`, surrogate `ORGSPC`=122). Its 53.32 height comes
from the **FVS-native Wykoff HT-D dubbing** (`oc/cratet.f:571` `INVENTORY EQN DUBBING ISPC=12 →
53.3234`, via `HTDBH`), *not* ORGANON — PREPARE's dub of tree-20 is computed then discarded (CRATET
reloads only `IORG=1` trees, `oc/cratet.f:349-365`). On this stand **no valid ORGANON tree has a
missing HT/CR** (`KNTOHT=KNTOCR=0`), so the ORGANON PRDHT/PRDCR dubbing — ported faithfully — is not
exercised here; a real seedling/partial-inventory stand is needed to bit-validate it. The FVS-native
Wykoff dubbing (shared engine + OC `HT1`/`HT2` coefficients) is a distinct follow-up, orthogonal to
the ORGANON calibration.

### Deferred within C2 (measured, not gaps)

- **`DGCALIB` (RAD=.TRUE.) branch** — the diameter-growth calibration reuses `DG_SWO`/bark/BAL
  (chunk C3). FVS/FIA inventory carries no radial-increment cores ⇒ `RAD=.FALSE.` and
  `TMPCAL(3,*)=1.0` (`organon/prepare.f:137-140`), confirmed by the oracle (all `TMPCAL(3,*)=1.0`).
  The port raises a loud error if `RADGRO>0` is ever passed, to be lifted when C3 lands.
- **NWO(2)/SMC(3)/RAP(4)** version branches (OP=Olympic) raise a clear error; their coefficient
  tables sit alongside SWO in the same Fortran and are a thin OP follow-on.

**C2 verdict: bit-exact vs the live oracle** on the ORGANON PREPARE calibration (`TMPCAL`/`ACALIB`),
the DF↔PP site conversion, and the WF/ES species crosswalk. HT/CR dubbing ported faithfully (not
exercised on ocmin). **Next: C3 (ORGANON diameter growth — `diagro.f` DG_SWO + `diamcal.f` bark/BAL
+ `submax.f` + `statsorg.f`); this also lifts the DGCALIB RAD deferral.**

## Chunk C3 delivered (ORGANON SWO diameter growth) — VALIDATED BIT-EXACT

`src/variants/oregoncoast/organon_diamgro.jl` ports the ORGANON **SWO (VERSION=1)** diameter-growth
core — the deterministic (`DGSD=0`) per-tree `DGRO` that overwrites the FVS Wykoff LN(DDS) for the
18 valid ORGANON species. FVS serial-corr is suppressed on the ORGANON path (`OLDRN=0`, `FRM=1`,
`oc/dgdriv.f:549-550`), so `DGRO` is a hard bit-exact target with **no RNG straddle**.

- **`oc_dg_swo`** = `organon/diagro.f` DG_SWO: `LNDG = B0 + B1·ln(DBH+K1) + B2·DBH^K2 +
  B3·ln((CR+0.2)/1.2) + B4·ln(SITE) + B5·(SBAL1^K3/ln(DBH+K4)) + B6·√SBA1`, then `DG =
  exp(LNDG)·CRADJ·ADJ`. Full 18-group `DGPAR(18,11)` table + the per-group `ADJ` + the `CR≤0.17`
  crown adjustment. `SITE = SI_1 = SITE_1−4.5` (`execute2.f:322`, DIAMGRO_RUN CASE(1)). All REAL*4
  math via `fexp`/`flog`/`fpow`; `√` via `Float32` `sqrt`; `DBH*DBH`/`MCW*MCW` as exact integer
  products (matching `**2`).
- **`oc_sstats`** = `organon/statsorg.f` SSTATS: stand `SBA` / `BAL(500)` / `BALL(51)` / `CCFL` /
  `CCFLL` / `TPA` / `SCCF` from the buffer (`NPTS=1`, so `/FLOAT(NPTS)` is the identity). **`oc_get_bal`**
  = `diamcal.f` GET_BAL (the `SBAL1` per-tree lookup DG consumes). **`oc_diamgro_run`** =
  DIAMGRO_RUN (`DGRO = DG·CALIB(3,g)·FERTADJ·THINADJ`); **`oc_dg_thin`/`oc_dg_fert`** ported
  faithfully (both return 1.0 with no thin/fert).
- **`oc_submax`** = `organon/submax.f` SUBMAX (VERSION=1): the max size-density line A1/A2 (consumed
  by mortality C6 / the RD index; **inert on DG**, ported because the task lists it).
- **`organon_dgcalib_swo`** lifts the **C2 DGCALIB RAD deferral**: on FVS/FIA inventory there are no
  radial-increment cores ⇒ `RAD=.FALSE.` ⇒ `CALIB(3,*)=1.0` (oracle-confirmed all `TMPCAL(3,*)=1.0`).
  `DIB_SWO`/`DOB_SWO` bark are ported (`oc_dib_swo`/`oc_dob_swo`); a `RAD=.TRUE.` run raises a clear
  TODO error (unvalidated — no radial-core stand exists to bit-check it).
- **`organon_dg_swo(buf; …)`** wraps the `GROW` "growth-1" DG sequence (`grow.f:93-109`) off the C1
  `/ORGANON/` buffer: species groups → SSTATS → SUBMAX → DGCALIB → per-tree DIAMGRO_RUN, returning
  the `DGRO` vector. The FVS bark conversion (`BRATIO → DIAGR → DDS → WK2`, `oc/dgdriv.f:446-452`)
  and the StandState application are the **C7** copy-back seam (FVS-native shared engine), not C3.

**MEASURED vs the live oracle** (`FVSoc_clean`, scoped `DEBUG 1 / DGDRIV`, stand S248112 / ocmin, 27
records; the growth-cycle `I,ISPC,DBH,DGRO,BARK,DIAGR,DDS=` dump — recipe: `DEBUG` with a **non-blank
field-2** so `initre.f:1072` calls `DBPRSE` to scope debug to `DGDRIV`, avoiding the DBALL volume
crash). jl was fed the exact `/ORGANON/` buffer the Fortran received (`SPECIES,DBH1,HT1OR,CR1,EXPAN1`
— C1-validated), `SITE_1=92.0` (SI_1=87.5), `MSDI_1=815`:

| quantity | jl | oracle | result |
|---|---|---|---|
| `DGRO` — all **17 valid ORGANON trees** | — | — | **max \|Δ\| = 0.000e+00 (bit-exact)** |
| SUBMAX `A1` | 6.4790 | 6.4790 (`STOR(3)`) | **exact** |
| SUBMAX `A2` | 0.62305 | 0.62305 (`STOR(4)`=0.6230) | **exact** |
| `CALIB(3,*)` | all 1.0 | all 1.0 | **exact** |
| `SBA1` | 85.131271 | — | (consumed by DG, validated transitively) |

Every valid tree — DF (g1), GW/fir (g2), PP (g3), SP (g4) — matched to full Float32 print precision
(`|Δ|=0.0`, not merely f32-print-close). Because ORGANON is deterministic (`DGSD=0`) this is the
hard bar, not a straddle. **C3 DIAMETER GROWTH: BIT-EXACT.**

### Measured notes (not gaps)

- **Eligibility uses the RAW FVS height, not `HT1OR`.** `oc/dgdriv.f:222` gates `IORG` on `HT(I)>4.5`
  (the FVS array), while the buffer floors `HT1OR` to 4.6. So tree-2 (`LP`/DF raw HT=2.0), tree-13
  (raw HT=3.0) and tree-24 (raw HT=2.0) are `IORG=0` despite `HT1OR=4.6` — they get FVS-native
  growth, not ORGANON. `build_organon_buffer!` already gates on the raw height (C1-correct); the 17
  ORGANON trees are exactly the oracle's `IORG=1` set.
- **The copy-back bark is FVS-native `BRATIO`**, not ORGANON DIB/DOB. `oc/dgdriv.f:446` `BARK=BRATIO`,
  `DIAGR=DGRO·BARK`, `DDS=ln(DIAGR·(2·DBH·BARK+DIAGR))`. Reproducing `DDS` end-to-end needs the FVS
  bark ratio (shared engine) — that is the **C7** wiring, orthogonal to the ORGANON DG core.
- **DG_THIN/DG_FERT** return exactly 1.0 on cyc0 (no thin/fert); ported faithfully for later cycles.

**C3 verdict: bit-exact vs the live oracle** on ORGANON SWO diameter growth (`DGRO`, all 17 valid
trees), SUBMAX A1/A2, and the DGCALIB CALIB(3,*). The C2 DGCALIB/RAD deferral is lifted
(RAD=.FALSE. path returns 1.0; RAD=.TRUE. machinery present, unvalidated as no radial-core stand
exists). **Next: C4 (height growth `htgrowth.f` HG_SWO → validate `HGRO/HTG` from `oc/htgf.f:97`);
then C5 crown, C6 mortality, C7 GROW/EXECUTE orchestration + the FVS DDS copy-back.**

## Chunk C4 delivered (ORGANON SWO height growth) — VALIDATED BIT-EXACT

`src/variants/oregoncoast/organon_htgro.jl` ports the ORGANON **SWO (VERSION=1)** height-growth
core — the deterministic (`DGSD=0`) per-tree `HGRO` that FVS copies into `HTG` for the valid
ORGANON trees, bypassing the native OC HTGF (`oc/htgf.f:95-101`: `HTG=SCALE·XHT·HGRO·EXP(HTCON)`;
on cyc0 defaults SCALE=FINT/YR=1, XHT=1, HTCON=0 ⇒ `HTG=HGRO`).

- **`oc_htgro1`** = `organon/htgrowth.f` HTGRO1 (big-6, species groups 1..5): per tree, `TCCH`
  (crown competition, interpolated from the CCH profile) → **`oc_hs_hg`** (Hann-Scrivani potential
  height growth + growth-effective age) → **`oc_hg_swo`** (the `MODIFER`/`CRADJ` height-increment
  equation, `HGPAR(5,8)`) → **`oc_hg_fert`/`oc_hg_thin`** (1.0 no treatment) → **`oc_limit`** (caps
  HG under the H-D curve using the C3 5-yr `DGRO`). SWO SITE selection: PP(122)→SI_2/ISISP=2,
  IC(81)→(SI_1+4.5)·0.66−4.5, else SI_1/ISISP=1.
- **Crown-closure profile** (`organon/crngrow.f`): **`oc_crnclo`** (CRNCLO, IND=0/SCR=0) builds
  `CCH(1..41)` from the start-of-cycle tree list via **`oc_lcw_swo`** (LCW), **`oc_hlcw_swo`**
  (HLCW), **`oc_cw_swo`** (CW above LCW) and **`oc_calc_cc!`** (CALC_CC, 40 height strata). `MCW_SWO`
  is reused from C2.
- **`organon_hg_swo(buf, dgro, spgrp; si_1, si_2)`** = the `GROW` "growth-2" HG sequence
  (`grow.f:133-152`): build CRNCLO, then HTGRO1 for every big-6 tree, returning the `HGRO` vector
  (`oc/htgf.f:96`'s HTG source). It consumes the C3 `dgro`/`spgrp` (LIMIT needs the diameter growth).

**MEASURED vs the live oracle** (`FVSoc_clean`, scoped `DEBUG 1 / DGDRIV HTGF`, stand S248112 /
ocmin; the `HTGF ORGANON … HGRO …` dump — **17 valid big-6 ORGANON trees**). Same exact `/ORGANON/`
buffer, `SITE_1=92.0` (SI_1=87.5), `SITE_2=86.5528641` (SI_2=82.0528641), `CCH(41)=75.0`:

| quantity | result |
|---|---|
| `HGRO` — all **17 valid big-6 ORGANON trees** (DF g1, GW/fir g2, PP g3, SP g4) | **max \|Δ\| = 0.000e+00 (bit-exact)** |

**C4 verdict: bit-exact vs the live oracle** on ORGANON SWO height growth (`HGRO`, all 17 valid
trees). All 27-tree CRNCLO crown-closure + the per-tree HS_HG/HG_SWO/LIMIT chain reproduce to full
Float32 print precision. Deterministic (`DGSD=0`) — the hard bar, no straddle.

### Measured notes (not gaps)

- **C4 requires the gfortran `fexp/flog/fpow` shim (`deps/libfvsmath.so`) to be BIT-EXACT** — unlike
  C3. HS_HG's nested `exp(exp(…))` + `pow(…, 1/b2)` chain is ULP-sensitive: with the openlibm
  fallback, 3 of 17 trees drift 1–15 ULP (max rel 1.75e-6); with the shim active, **all 17 are
  exactly 0.0**. This concretely validates doctrine #8 (route all ORGANON REAL*4 math through the
  gfortran ops). **Infra note:** a fresh git worktree lacks the build artifact `deps/libfvsmath.so`,
  and `FMath._ensure_built` probes `gfortran` (this env only has `gfortran-16`), so the shim silently
  falls back. Build it once per worktree: `gfortran-16 -shared -fPIC -O2 -o deps/libfvsmath.so
  deps/fvsmath.f90`. (The merge target carries the .so, so this is a worktree-only measurement step.)
- **HTGRO2 (minor ORGANON species, groups > 5)** — the HD-ratio height-growth form for RC/PY/CY/WO/
  BO/BM/RA/MA/GC/DG/TO/WI — is **not exercised by ocmin** (all 17 valid trees are big-6). Ported as a
  follow-up when a stand carries a minor ORGANON species; `organon_hg_swo` leaves `hgro=0` for them.
- **HG_FERT/HG_THIN** return exactly 1.0 on cyc0 (no thin/fert); ported faithfully for later cycles.
- **LIMIT** is inactive on all 17 ocmin trees (no tree hit the H-D cap), but is ported and wired
  (consumes the C3 DGRO); its DG-dependent cap will engage on faster-growing stands.

**Next: C5 (crown — `crngrow.f` CROWGRO/HCB/CW → validate `CR2/CRNEW` from `oc/crown.f:285`); then
C6 mortality (`mortality.f` PM_SWO), C7 GROW/EXECUTE orchestration + the FVS DDS/HTG copy-back.**

## Chunk C5 delivered (ORGANON SWO crown growth) — VALIDATED BIT-EXACT

`src/variants/oregoncoast/organon_crngro.jl` ports the ORGANON **SWO (VERSION=1)** crown-recession
core — the deterministic (`DGSD=0`) per-tree `CR2` that FVS rounds into `CRNEW` for the valid
ORGANON trees (`oc/crown.f:282`: `CRNEW=ANINT(CR2·100)`), bypassing native OC crown. CROWGRO runs
AFTER DG (C3) + HG (C4) advance DBH/HT and mortality (C6) reduces the expansion.

- **`organon_cr_swo(buf, dgro, hgro, spgrp, deadexp; si_1, si_2)`** = `organon/crngrow.f` CROWGRO:
  compute start (old DBH/HT) & end (new DBH/HT) height-to-crown-base via **`oc_hcb_swo`** (HCB_SWO,
  its own `HCBPAR(18,7)` — DISTINCT from C2's start2 A_HCB_SWO), the crown-base growth
  `HCBG=max(0, HCB2−HCB1)`, the max crown-base cap **`oc_maxhcb_swo`** (MAXHCB_SWO), and the
  actual/shadow crown-base recession branch → new `CR2`.
- **`oc_oldgro`** = `organon/mortality.f` OLDGRO — the old-growth indicator OG (5-largest big-6),
  start (`XIND=-1`, subtract growth + add DEADEXP) and end (`XIND=0`). Reuses `oc_get_ccfl` (C2),
  `oc_sstats` (C3), `oc_mcw_swo` (C2).
- The two SSTATS calls: start stats use the **original** expansion (SOG, pre-mortality); end stats
  use the **survivor** expansion and the **new DBH with OLD HT** (grow.f:169 runs after the DBH
  update but before the HT update at :184).

**MEASURED vs the live oracle** (`FVSoc_clean`, scoped `DEBUG 1 / DGDRIV HTGF CROWN`, stand S248112
/ ocmin; the `ORG CROWN … CR2 …` dump — **17 valid ORGANON trees**). jl fed the exact `/ORGANON/`
buffer + the C3 `dgro` + C4 `hgro` + the measured per-tree `deadexp` (see C6 note):

| quantity | result |
|---|---|
| `CR2` — all **17 valid ORGANON trees** (DF g1, GW/fir g2, PP g3, SP g4) | **max \|Δ\| = 0.000e+00 (bit-exact)** |

**C5 verdict: bit-exact vs the live oracle** on ORGANON SWO crown growth (`CR2`, all 17 valid trees).

### Root-caused during C5 (load-bearing)

- **The end-of-growth stand stats use SURVIVOR expansion, not original.** ORGANON mortality
  (`mortality.f:233-234`) runs inside GROW *before* the end-of-growth `SSTATS`/`CROWGRO` and does
  `DEADEXP=EXPAN·PM; TDATAR(4)=EXPAN·(1−PM)`. Feeding SSTATS2/OLDGRO the original expansion under-
  predicted every CR2 by ~0.5–1.2% (max rel 5%); using `surv = EXPAN − DEADEXP` makes all 17 exactly
  0.0. This is the crown analogue of the C4-requires-shim finding: a stand-level input, not a coeff.
- **C6 dependency (measured, not a gap):** OLDGRO's OG1 and the SSTATS2 survivor expansion both need
  the ORGANON mortality `DEADEXP`. C5 is validated with the oracle-measured MORTEXP fed in (as C4
  consumed the measured/ported DGRO). Once C6 (PM_SWO) lands, `organon_cr_swo` takes the ported
  DEADEXP with no code change.
- **Shadow-crown (SCR) path inert:** SCR1B=0 on inventory ⇒ `AHCB1 ≤ SHCB1` always, so the shadow
  branch never fires here; it is ported faithfully for completeness.

**Next: C6 (mortality — `mortality.f` PM_SWO + RAMORT/OLDGRO, supplies DEADEXP); then C7
GROW/EXECUTE orchestration + the FVS DDS/HTG/CR copy-back into `diameter_growth!(::OregonCoast)`.**

## Chunk C6 delivered (ORGANON SWO mortality) — VALIDATED BIT-EXACT

`src/variants/oregoncoast/organon_mortality.jl` ports the ORGANON **SWO (VERSION=1)** mortality core
— the deterministic (`DGSD=0`) per-tree `DEADEXP = EXPAN·PM` that C4/C5 have been consuming as the
measured MORTEXP. FVS copies it as `WK2=MORTEXP·(FINT/5)` (`oc/morts.f:499`; MORTEXP=DEADEXP·NPTS).

- **`organon_mortal_swo(buf, dgro, hgro, spgrp, bal1, ball1, a1, a2; si_1, cyclg, mort)`** =
  `mortality.f` MORTAL_RUN: stand `STBA/STN/SQMDA/RD`, `OLDGRO(XIND=0)` pre-growth OG, per-tree
  **`oc_pm_swo`** (PM_SWO logistic linear predictor, `MPAR(18,9)`) + **`oc_pm_fert`** (0 no fert),
  then `PM = 1 − (1−logistic(PMK))^POW·CRADJ` and `DEADEXP = EXPAN·PM`. The SDI additional-mortality
  block (MORT=INDS(9)=1, `oc/grinit.f:364`) is ported: `RDA`, the CYCLG=0 `IND`/`A1MAX`/`NO`
  initialization, **`oc_quad1`** (QUAD1), and the KR1 density-adjustment iteration.
- Runs inside GROW after growth-1/2 but before the DBH/HT update, so it sees the ORIGINAL DBH/HT/
  EXPAN + the C3 `dgro` (for the post-growth BA `(DBH+DG)²`), the start-of-growth BAL (`bal1`/`ball1`
  from C3's SSTATS), and SUBMAX `a1`/`a2` (from C3).

**MEASURED vs the live oracle** (`FVSoc_clean`, scoped `DEBUG 1 / DGDRIV HTGF CROWN MORTS`, stand
S248112 / ocmin; the dgdriv `MORTEXP` dump — **ALL 27 records**):

| quantity | result |
|---|---|
| `DEADEXP` (=MORTEXP) — all **27 trees** (valid + surrogate, groups 1–4) | **max \|Δ\| = 0.000e+00 (bit-exact)** |
| C5 `CR2` recomputed from the C6-**ported** `deadexp` (17 valid trees) | **max \|Δ\| = 0.000e+00** |

**C6 verdict: bit-exact vs the live oracle** on ORGANON SWO mortality (`DEADEXP`, all 27 trees). The
measured-MORTEXP dependency C4/C5 carried is now **fully closed**: feeding the C6-ported `deadexp`
back into `organon_cr_swo`/SSTATS2 reproduces the 17-tree CR2 (and the DGRO/HGRO) at max |Δ|=0.0 with
no change on their side.

### Measured notes (not gaps)

- **SDI additional mortality did not engage on ocmin** — `RD = STN/exp(A1/A2 − ln(SQMDA)/A2) ≤
  RDCC=0.60` (below carrying capacity, SDI≈184), so the base individual-tree path applies and the
  KR1 density-adjustment loop is inert. Both are ported; the KR1 branch awaits a dense stand to
  bit-validate. **PM_SWO's `POW` output overwrites the caller's POW (=MPAR(g,9)=1.0 for all SWO).**
- **PM_FERT/RAMORT inert:** no fertilizer (FERTADJ=0) and no red alder ≥55 yr (RAMORT skipped).
- **Subsequent-cycle mortality init** (`mortality.f:177-204`, needs the carried `RD0/PA1MAX/NO`
  state) is a C7-orchestration follow-up; only the cyc0 (`CYCLG==0`) initialization is ported here.

**Next: C7 — GROW/EXECUTE per-cycle orchestration wiring C3–C6 + the FVS DDS/HTG/CR/MORTEXP
copy-back into StandState, so an end-to-end ocmin/oct01 cyc0 run validates against FVSoc_clean via
`diameter_growth!(::OregonCoast)`.**

## Chunk C7 sub-step 1 delivered (EXECUTE/GROW orchestration + DGRO→DDS copy-back) — BIT-EXACT

`src/variants/oregoncoast/organon_execute.jl` ties the four validated growth components (C3–C6) into
one per-cycle call in the faithful `execute2.f` EXECUTE + `grow.f` GROW order, and ports the FVS-side
copy-back seam that turns the ORGANON outputs into FVS tree-record increments.

- **`organon_execute_swo(buf, isp_fvs; si_1, si_2, msdi…)`** = the GROW sequence **DG (C3) → HG (C4)
  → MORTAL (C6) → CROWGRO (C5)** — mortality BEFORE crown (so C5's survivor-expansion is honoured),
  returning an `OrganonGrowth` with the per-tree `dgro`/`hgro`/`cr2`/`deadexp` + the `dds` (WK2) FVS
  loads for the ORGANON trees.
- **`oc_bratio`** = `bin/FVSoc_buildDir/bratio.f` — the OC (CA-family) variant bark ratio
  (`BARKB(5,29)` / `JBARK(50)`, three eqn forms, clamp [0.80,0.99]). OC has NO `oc/bratio.f`; the
  build links this CA-family `bratio.f`, distinct from the shared linear `bark_ratio`.
- **`oc_organon_dds`** = `oc/dgdriv.f:446-452` — `BARK=BRATIO`, `DIAGR=DGRO·BARK`,
  `DDS=ln(DIAGR·(2·DBH·BARK+DIAGR))` floored −9.21 → the FVS WK2. **`oc_organon_dg`** = the shared
  DDS→DG (`√((DBH·BARK)²+exp(DDS))−DBH·BARK`, OLDRN=0/FRM=1) the StandState apply-loop consumes.

**MEASURED vs the live oracle** (`FVSoc_clean`, scoped `DEBUG 1 / DGDRIV HTGF CROWN MORTS`, stand
S248112 / ocmin) — every ORGANON per-tree quantity from a SINGLE `organon_execute_swo` call:

| quantity | trees | result |
|---|---|---|
| `DGRO` (C3) | 17 valid | **max \|Δ\| = 0.000e+00** |
| `BARK` (oc_bratio) | 17 valid | **0.000e+00** |
| `DIAGR` = DGRO·BARK | 17 valid | **0.000e+00** |
| `DDS` (WK2 FVS load) | 17 valid | **0.000e+00** |
| `HGRO` (C4) | 17 valid | **0.000e+00** |
| `CR2` (C5) | 17 valid | **0.000e+00** |
| `DEADEXP` (C6) | all 27 | **0.000e+00** |

**C7 sub-step 1 verdict: bit-exact.** The orchestration order + the OC bark ratio + the DGRO→DDS
copy-back all reproduce the oracle from one call. This is the ORGANON side of C7 fully assembled.

## Chunk C7 sub-step 2 delivered (live growth hook + StandState copy-back) — BIT-EXACT

`src/variants/oregoncoast/organon_hook.jl` runs the ORGANON growth on a real `StandState` and copies
the outputs into the tree records; `diameter_growth!(::OregonCoast)` now drives it (was a loud-error
stub through C6).

- **`organon_apply_growth!(s; msdi=0, cyclg=0, fint=5)`** — builds the C1 buffer from the StandState,
  runs `organon_execute_swo`, and for the valid ORGANON trees (IORG=1) grows `DBH += oc_organon_dg/
  oc_bratio` (outside-bark), `HT += HGRO`, sets `crown_pct = ANINT(CR2·100)` (Julia
  `round(…, RoundNearestTiesAway)` = Fortran ANINT), and reduces `TPA -= MORTEXP·(FINT/5)` for every
  record ORGANON grew (valid + surrogate, `oc/morts.f:498`). `diameter_growth!(::OregonCoast)` calls it.
- **`_oc_organon_si`** derives SI_1/SI_2 from the plot DF/PP site indices (mirrors C2 SITSET/execute2).

**MEASURED vs the live oracle** — an ocmin `StandState` (27 trees) grown by `organon_apply_growth!`:

| copy-back | check | result |
|---|---|---|
| crown `ANINT(CR2·100)` | vs oracle `CRNEW` (dgdriv/crown dump), 17 valid | **17/17 exact** |
| `HT += HGRO` | vs `ht₀ + HGRO`, 17 | **max \|Δ\| = 0.0** |
| `TPA −= DEADEXP` | vs `tpa₀ − DEADEXP`, **all 27** | **max \|Δ\| = 0.0** |
| `DBH += DG/oc_bratio` | FVS sqrt-path (`√((D·BARK)²+exp(DDS))−D·BARK`)/BARK | applied (DDS bit-exact) |

**C7 sub-step 2 verdict: bit-exact** on the StandState copy-back for the ORGANON trees. The crown
`CRNEW` is a genuine new oracle check (Fortran `ANINT` ties-away rounding). MSDI is inert on ocmin
(`RD ≤ RDCC` ⇒ base mortality; `A1` 6.479 vs 6.294 for msdi 815 vs 0 give byte-identical DEADEXP), so
the hook is bit-exact with the default `msdi=0`.

## Chunk C7 sub-step 3 delivered (grow_cycle! seam reconciliation) — no double-apply

`diameter_growth!(::OregonCoast)` is now the **single growth authority** in the shared `grow_cycle!`:
because ORGANON computes diameter/height/crown/mortality together in one EXECUTE, it runs the whole
`organon_apply_growth!` (applying DBH/HT/CR/TPA) and then **zeros** `diam_growth`/`ht_growth` so the
engine's later apply-loop (`DBH+=DG/bark`, `HT+=HTG`) is inert, and the OC `height_growth!` /
`small_tree_growth!` / `mortality!` / `crown_ratio_update!` hooks are **no-ops**. This is faithful
(FVS's `oc/dgdriv.f`→EXECUTE likewise produces DG/HTG/CR/MORTEXP in one call) and avoids the
double-apply a cooperating-hook split would risk.

**VERIFIED no double-apply:** replaying the `grow_cycle!` growth order (`diameter_growth!` →
apply-loop → the four no-op hooks) on an ocmin `StandState` gives **byte-identical** DBH/HT/CROWN/TPA
to a single `organon_apply_growth!` — max |Δ| = 0 on all four (applied exactly once).

## Chunk C8 delivered (OC setup height dubbing) — BIT-EXACT

`src/variants/oregoncoast/organon_cratet.jl` ports the OC missing-height dub. OC has
`LHTDRG=.FALSE.` for every species (`oc/grinit.f`), so `oc/cratet.f:679-684` uses the CA-family
Curtis-Arney `HTDBH` (INVENTORY equation) — it overwrites the calibrated-Wykoff `H`. `oc_htdbh_height`
= `bin/FVSoc_buildDir/htdbh.f` MODE=0 (`CURARN(50,3)`+`SPLINE`); `d ≥ Z`: `H=4.5+P2·exp(−P3·D^P4)`,
`d < Z`: the linear spline to 4.51@D=0.3. An OC branch in the shared `dub_missing_heights!`
(`volume.jl`) routes the else-case there; `D ≤ 0.1 → 1.01` is handled by the shared code. OC merch
defaults (`oc/grinit.f`/`oc/sitset.f`: DBHMIN=7 / sp-11=6, westside TOPD=4.5, stump=1) added to
`init_merch_standards!`.

**MEASURED bit-exact vs the live FVSoc_clean oracle** (ocmin CRATET `INVENTORY EQN DUBBING` dump):
`oc_htdbh_height` on all 5 dumped cases — LP(12) D=8.5→**53.3234**, DF(7) D=10.4→**64.6109**,
SP(16) D=8.0→**39.8626** / D=34.6→**121.308**, LP(12) D=7.2→**45.5013** — all exact. The ocmin
missing-height tree-20 (LP, D=8.5) dubs to 53.3234, matching the oracle.

**C8 verdict: bit-exact.** The `run_keyfile(ocmin)` setup now passes `dub_missing_heights!` and
`init_merch_standards!`.

## Chunk C8b delivered (OC ecoclass/site-index sourcing) — run completes, cyc0 `.sum` BIT-EXACT

`site_setup!(::OregonCoast)` now sources the OC default-ecoclass site indices + SDImax, so
`run_keyfile(ocmin)` **runs end-to-end** and emits a `.sum`. Ported (`oc/sitset.f`):
- default ecoclass `CWC221` (`sitset.f:105-110` `ICL5==0`; `ecocls.f:294`) → `SITEAR(7)=92` (DF site
  species, ISISP=7), `SDImax=815`, when no SITECODE keyword set a site (NSISET==0);
- the ORGANON DF↔PP conversion (`sitset.f:181-186`, PP = 0.940792·92 = 86.5528641);
- the `R6ADJ(50)` site fan (`sitset.f:68-73,197`): `SITEAR(I) = HGUESS·R6ADJ(I)`, HGUESS = 92/R6ADJ(7);
- MSDI = `sp_sdi_def[7]` = 815 sourced into the growth hook (`sitset.f:340-342` RVARS(3..5)=SDIDEF).

**MEASURED bit-exact vs the oracle SITECODE dump**: DF=92, PP=86.5529, WH=87.4(→87), SP=92, IC=64.4
(→64), PC=82.8(→83), GF=92 — all match; site species 7, SDImax 815.

**END-TO-END `run_keyfile(ocmin; output=:sum)` — cyc0 inventory row BIT-EXACT:**

| year | TPA | BA | SDI | TopHt | QMD | verdict |
|---|---|---|---|---|---|---|
| **1990 (cyc0)** jl | 536 | 77 | 184 | 63 | 5.1 | — |
| 1990 (cyc0) oracle | 536 | 77 | 184 | 63 | 5.1 | **BIT-EXACT** |
| 1995 (cyc1) jl | 504 | 85 | 197 | 70 | 5.6 | TPA exact; BA/SDI/TopHt/QMD low |
| 1995 (cyc1) oracle | 504 | 88 | 202 | 71 | 5.7 | (non-ORGANON growth pending) |

The **cyc0 (1990) `.sum` density row is bit-exact** (setup/dubbing/site/density all correct). The
cyc1 (1995) grown row's **TPA=504 matches** (mortality — all 27 trees via ORGANON DEADEXP — is
correct) but BA/SDI/TopHt/QMD are LOW because the IORG=0 non-ORGANON trees (ocmin LP/BR/sub-4.5-ft
DF) don't grow yet (their DG/HTG=0). This isolates the cyc1 residual entirely to **C9**.

### Remaining for the cyc1 grown `.sum` (isolated)

1. **Non-ORGANON DGF/HTGF (C9)** — `oc/dgf.f` (CA-family Wykoff DDS) + `oc/htgf.f` native path for the
   IORG=0 surrogate / no-big-6 trees; `organon_apply_growth!` leaves their DG/HTG at 0. This is the
   ONLY thing between the current cyc1 row and bit-exact (setup, site, ORGANON growth, mortality all
   verified correct).
2. **OC volume (C10)** — `compute_volumes!` NVEL/CA-family for the `.sum` vol columns (currently 0).
3. Full OC ECOCLS/HABTYP plant-assoc table (non-default ecoclasses; ocmin uses the default) + carried
   mortality state (A1MAX/NO/RD0) for multi-cycle.

## Chunk C9 step 1 delivered (non-ORGANON DGF equation + DGCONS + forkod + point CCF) — BIT-EXACT

`src/variants/oregoncoast/organon_dgf.jl` ports the OC FVS-native Wykoff diameter growth (`oc/dgf.f`)
for the IORG=0 trees — the path `oc/dgdriv.f` runs for EVERY tree before ORGANON overwrites WK2 for
IORG=1:
- **`oc_dgcons!`** = `oc/dgf.f` DGCONS: `DGCON = DGFOR[ispfor,jspc] + DGEL·ELEV + DGELSQ·ELEV² +
  DGSITE·ln(SITEAR) + slope/aspect`, via `MAPSPC(50)`→13 groups, `MAPLOC(10,13)` location, `DGFOR(5,13)`.
- **`dgf!(::OregonCoast)`** = the per-tree `ln(DDS)` (9-term Wykoff + 10→5-yr halving, tanoak ×2).
- **`oc_forkod`** = `oc/forkod.f` JFOR: KODFOR 711 → IFOR 9 (load-bearing for DF's DGFOR location class).
- **`oc_r5crwd` / `oc_tree_ccf`** = `bin/FVSoc_buildDir/r5crwd.f` R5CRWD (CA-family crown width) +
  `oc/ccfcal.f` MODE=1 `CCF=0.001803·CRWD5²`, wired into `point_density!` for the DGF `PCCF` term.

**MEASURED bit-exact vs FVSoc_clean DEBUG-DGF (ocmin, growth cycle):**
- `DGCON`: DF=1.14167, GF=0.07280, LP=0.44505, BR=−0.02756 — **exact**.
- per-tree `ln(DDS)`, all **10 IORG=0 trees** (LP×4, BR×4, sub-4.5-ft DF/GF): **max |Δ| = 4.8e-05**
  (F7.4 oracle print precision). The point `PCCF` root-caused: jl's generic crown width gave PCCF≈0.9
  vs oracle 20.7 (LP)/308.8 (BR) — the R5CRWD port makes them match (BR's ~0.31 uniform DDS offset,
  = `DGPCCF·PCCF`, closed).

COR is 0 for LP/BR (< FNMIN GSTs); only GF has COR=0.288 (an ORGANON species — used only by the
sub-4.5-ft GF tree-13). The DDS validation fed COR(4)=0.288 from the oracle (the COR calibration is
**C9 step 3**). **C9 step 1 verdict: bit-exact.** Not yet wired into `diameter_growth!` (step 2).

## Chunk C9 step 2 delivered (diameter_growth! re-architecture) — cyc1 BA/SDI/QMD/TPA BIT-EXACT

`organon_apply_growth!` now mirrors `oc/dgdriv.f`: `oc_dgcons!` + `dgf!` compute WK2 for EVERY tree
(original DBH), then `organon_execute_swo` overwrites WK2 for the IORG=1 trees; a unified apply loop
grows DBH by `DG/BARK` for all (ORGANON `dds` for IORG=1, DGF `wk2` for IORG=0), HT/CROWN for IORG=1,
and MORTEXP for all. `site_setup!` sources `forest_idx` via `oc_forkod` (711→9).

**END-TO-END `run_keyfile(ocmin; output=:sum)` cyc1 (1995):**

| year | TPA | BA | SDI | TopHt | QMD | verdict |
|---|---|---|---|---|---|---|
| 1995 (cyc1) jl | 504 | 88 | 202 | 70 | 5.7 | TPA/BA/SDI/QMD **BIT-EXACT** |
| 1995 (cyc1) oracle | 504 | 88 | 202 | 71 | 5.7 | TopHt −1 (IORG=0 HTGF pending) |

The IORG=0 trees now grow (were DG=0), closing BA 85→88 / SDI 197→202 / QMD 5.6→5.7 to the oracle.
The lone remaining `.sum` residual is **TopHt 70 vs 71** — the IORG=0 trees' HEIGHT growth
(`oc/htgf.f` native path) is still 0, so the tallest non-ORGANON tree doesn't gain height. That is
**C9 step 4**. The ORGANON IORG=1 trees stay bit-exact (no regression); cyc0 row unchanged. The COR
calibration (**step 3**) is negligible on ocmin (LP/BR COR=0; only the 0.1" GF tree-13 uses COR=0.288,
DG=0.004) — the row is bit-exact without it; step 3 is faithfulness for other stands.

## Chunk C9 step 4 delivered (native HTGF for IORG=0 trees) — cyc1 density row FULLY BIT-EXACT

`src/variants/oregoncoast/organon_htgf.jl` ports the OC FVS-native height growth (`oc/htgf.f` CASE
DEFAULT) for the IORG=0 trees: **`oc_findag`** (`findag.f` — growth-effective age by iterating the
site curve to the tree height) + **`oc_htcalc`** (`htcalc.f` — region-6 per-species potential-height
site curves: Hann-Scrivani / Dahms lodgepole / Dolph red fir / Powers oak / Porter-Wiant) + the
SMHMOD crown/relative-height modifier: `HTG = POTHTG · 1.016605 · CRMOD · RHMOD`. Wired into
`organon_apply_growth!` (IORG=0 branch).

**MEASURED bit-exact vs FVSoc_clean DEBUG-HTGF (ocmin)**: all **10 IORG=0 trees' HTG max |Δ| = 0.0**
(LP tree-1 6.90931; 0.1"-regen GF 9.05224 / DF 8.73122; BR 1.52–4.42).

**END-TO-END `run_keyfile(ocmin; output=:sum)` cyc1 (1995) — density row FULLY BIT-EXACT:**

| year | TPA | BA | SDI | TopHt | QMD |
|---|---|---|---|---|---|
| **jl** | 504 | 88 | 202 | **71** | 5.7 |
| oracle | 504 | 88 | 202 | 71 | 5.7 |

TopHt closed 70→71. cyc0 unchanged; ORGANON IORG=1 trees unregressed. The **only** remaining `.sum`
gap is the volume columns (currently 0) → **C10**.

## Oracle status

- **Relinked OK.** `/workspace/.ocwork/FVSoc_clean` built via
  `scratchpad/relink_oc.sh clean` = `gfortran-16 -o FVSoc_clean $(bin/FVSoc_buildDir/*.o)
  /workspace/.crwork/isoc23_shim.o` (621 objects + the shared isoc23 shim). A permanent
  `relink_oc.sh` (clean / single-`.o`-swap forms, mirroring `relink_ut.sh`) belongs in
  `/workspace/.ocwork/` for instrumented A/B runs.
- **Runs growth; SIGSEGVs at the volume stage.** `./FVSoc_clean --keywordfile=ocmin.key` completes
  cycle-0 growth and writes the full DGF DEBUG dump to `ocmin.out`, then crashes:

  ```
  #7  master.0.fvsvol  at fvsvol.f:530
  #8  natcrs_          at fvsvol.f:54
  #9  vols_            at vols.f:319
  #10 fvs_             at fvs.f:211   (cycle-0 volume call)
  ```

  This is a **volume-stage crash (NATCRS / `fvsvol.f`), not a growth crash** — the DGF/ORGANON
  growth path runs to completion first, so it does **not** block per-chunk growth validation. It is
  the same class as the documented western NVEL/volume crashes; fixing it (per the "crash ⇒ fix
  live FVS" doctrine) is a separate follow-up, most cheaply pinned with a `--check-bounds`-style
  `-fbounds-check` rebuild of `fvsvol.o`/`vols.o` swapped into the relink. Not required for the
  growth port. Note `oct01.sum.save` exists and was produced by an original (non-g16) build, so a
  clean full-projection A/B on the `.sum` is available once the volume crash is resolved.

## Oracle-confirmed reference values (cyc0, `ocmin.key`, S248112, 27 records)

- Header edition: `ORGANON SWO`; site species `DF CODE= 7`; `STAGE,BHAGE = 60,54`.
- `NVALID,NLOAD = 17,27` — 17 of 27 records are valid ORGANON trees; `CRATET ORGANON ... IERROR= 0`.
- ORGANON species groups 1–17 with per-group `ACALIB`(diam/ht/crown ≈1.0) and `TMPCAL`
  calibration ratios (the calibration seam chunk C2 must reproduce bit-exact).

## Effort estimate (OC + OP)

- **OC (SWO):** the ORGANON growth core is ~9,200 lines (execute2/grow/diagro/htgrowth/crngrow/
  mortality/diamcal/statsorg/submax/growth_mods/whphg/prepare/start2), edition-branched, with a
  hard deterministic bit-exact bar (`DGSD=0`, ORGANON serial-corr suppressed — no RNG cover). This
  is the **largest single unported western workload** — comparable to a full BC-scale variant port
  *plus* a calibration subsystem. Sequence C1 (marshalling)→C2 (calibration) before the four growth
  chunks C3–C6, then C7 orchestration. Volume/woodquality (~4,600 lines) and tripling are deferred
  (OFF by default).
- **OP (Olympic):** a thin follow-on. Identical file set, same shared `organon/` engine; differs
  only in `MAXSP=39`, a different `orgspc.f` map, and dynamic `VERSION=NWO(2)/SMC(3)` selected in
  `op/sitset.f:70`. The NWO/SMC coefficient branches already sit alongside SWO in the ORGANON
  source, so OP ≈ OC-foundation chunk + version selection + the OP species map.

## Files touched

- `src/variants/oregoncoast/oregoncoast.jl` (new)
- `src/variants/oregoncoast/species.jl` (new)
- `data/oregoncoast/species_coefficients.csv`, `data/oregoncoast/species_translation.csv` (new)
- `src/variants/variant.jl` (OC registration line)
- `src/FVSjl.jl` (OC includes)
- `docs/OC_ORGANON_PORT_PLAN.md`, `docs/OC_VARIANT_PORT_AUDIT.md` (new)

## 2026-08-20 — OC multi-cycle GROWTH: real ~1% height-growth carry residual (localized, NOT yet at bar)

Ran the oct01 stand-1 growth control (S248112, 10 cycles, growth-only, REWIND-truncated) vs the LIVE
RELINKED /workspace/.ocwork/FVSoc_clean:
  • cyc0 (1990): bit-exact (TCF/MCF/BF 1894/1144/5881).
  • cyc1 (1995): STATE still bit-exact (504/88/202/2254/1574/8244) but the ACCRETION/MORT split already
    differs (accr 111→106, mort 20→21) — a symptom of a slightly different internal growth/mort
    decomposition that nets out at 1995 but not later.
  • cyc2+ (2000..): real drift. At 2000, **TPA 475/475 and BA 221/221 are BIT-EXACT** (same survival,
    same diameter) yet **TCF 2707→2682 (−0.9%), CCF 79→76** differ ⇒ identical DBH, lower volume+crown
    ⇒ the trees are slightly SHORTER in jl. So the residual is in **ORGANON SWO HEIGHT GROWTH's
    multi-cycle carry** (height/crown state fed forward), NOT mortality (TPA matches) or diameter (BA
    matches). By cyc3+ TPA drifts ±1-2 downstream (crown→density→mortality feedback). End-state 2040
    TCF 8227→8201 (−0.3%), TPA 308→305.

This is DGSD=0 (deterministic) so it SHOULD be bit-exact; the ~1% is a real residual, larger than the
cornered-NINT class. Matches the OC-agent's flagged open follow-on ("a multi-cycle OC stand for the
mortality-init carry A1MAX/NO/RD0" + HTGRO2 carry). Distinct from OP, which is bit-exact-or-cornered
multi-cycle (test_op_multicycle_sum) — OP=NWO, OC=SWO, different ORGANON height model.

NEXT OC chunk (measure-don't-infer): build an instrumented FVSoc_clean ORGANON height dump (g16
single-.o swap on the organon/htgro path, per-tree per-cycle HT), replay jl's OC height increment
cycle-by-cycle to find where the carry diverges. Bounded-but-real ORGANON investigation. NOT locked as
a bit-exact test yet (would fail); characterized here as the localized lead.

## 2026-08-20 (cont'd) — ROOT-CAUSED: OC/OP missing the REGENT small-tree height model (measure-don't-infer)

The ~1% multi-cycle drift above is NOT large-tree ORGANON — it is the FVS-native SMALL-TREE height model.
Measured per-tree via the oracle TREELIST (.trl) vs jl (write_sum_file cycle_hook capturing tree_id→HT;
harness in scratchpad/oc/):
  • cyc0 all 27 trees bit-exact. cyc1: LARGE trees bit-exact (tn1/17/24 dH~0.01) but the SUB-BREAST-
    HEIGHT seedlings blow out — **tn2 (DF, dbh 0.1, h 2.0): oracle→3.8 vs jl→10.7 (+6.9ft); tn15 (GF,
    h 3.0): oracle→5.1 vs jl→12.1 (+7.0ft)**. jl grows seedlings ~7ft too tall in one cycle.
  • ROOT: `oc/regent.f` (578 ln) is OC's small-tree model — "the height increment model is applied to
    trees LESS THAN 4 [ft] tall and LESS THAN 4 INCHES DBH; diameter is assigned from a height-diameter
    function." tn2/tn15 are <4ft ⇒ the oracle grows them with REGENT. jl routes ALL IORG=0 trees to the
    LARGE-tree `oc_htgf_native` (organon_hook.jl:86) with NO size gate, and `small_tree_growth!(::Oregon
    Coast/::Olympic)` is a NO-OP (organon_hook.jl:132) ⇒ seedlings get the large-tree HTGF (POTHTG·XMOD
    ≈9ft for a fast-site seedling) instead of the gentle REGENT height-age increment (~1.8ft).
  • The over-tall seedlings then perturb the height-based competition (CCH/relht) feeding the large-tree
    increment ⇒ large trees run ~0.15ft short by cyc2 ⇒ the TCF/CCF .sum deficit (BA/TPA bit-exact).

⇒ **NEXT OC/OP CHUNK (concrete): port `oc/regent.f` (small-tree height-age + H-D diameter assignment,
following the established regent.jl pattern used by BM/IE/KT/TT/UT/CR/CA/BC) and wire it into the OC/OP
ORGANON growth hook — route trees with H<4ft (and DBH<4") to it instead of oc_htgf_native.** VALIDATION
VEHICLE (ready): the oracle .trl per-tree HT (scratchpad/oc/oc_or_trl.pkl) vs jl's cycle_hook capture —
target = tn2/tn15 HT bit-exact at cyc1, then the whole .sum bit-exact-or-cornered multi-cycle. Same
model likely applies to OP (Olympic) sub-4ft trees; OP's ref stand happened to have fewer/less-divergent
seedlings so it read bit-exact-or-cornered, but OP should get the same regent wiring.

## 2026-08-20 (cont'd) — REGENT port fully specced (turnkey), deterministic path confirmed

Read oc/regent.f + smhtgf.f + htdbh.f in full and wrote scratchpad/oc/REGENT_PORT_SPEC.md — a turnkey spec
(all 5 smhtgf equations, the CON=RHCON·exp(HCOR) calibration, XWT XMN/XMX blend, SCALE, HK≤4.5-vs-HTDBH DBH
assignment, coefficient sources, wiring plan, validation vehicle). KEY: OC/OP are DGSD=0 ⇒ regent.f:234
ZZRAN=0 ⇒ the inventory small-tree path is FULLY DETERMINISTIC ⇒ bit-exact achievable. The port is entangled
(replaces both height AND diameter for D<XMAX IORG=0 trees — jl currently grows their DBH via the large-tree
DDS, wrong; regent derives DBH from the H-D function) and carries the documented regent trap history
(HCOR calibration, XWT blend, stale-HTGR) — land it as a focused unit with per-tree .trl A/B, not a rush.

## 2026-08-20 — REGENT small-tree HEIGHT model PORTED + validated (seedling +7ft bug FIXED)

Ported oc/smhtgf.f (5 height-age equations by MAPSP) + the regent.f height wrapper (CON=exp(HCOR),
XWT XMN/XMX blend, SCALE; DGSD=0 ⇒ ZZRAN=0 deterministic) as src/variants/oregoncoast/small_tree_growth.jl,
and wired organon_hook.jl to route IORG=0 trees with DBH<XMAX(sp) to it (was: all IORG=0 → large-tree
oc_htgf_native, which drops RELHT suppression when PCCF<100 and over-grew seedlings ~7ft/cycle).
VALIDATED vs the live FVSoc_clean scoped-DEBUG REGENT dump + TREELIST: DF firs eq bit-exact (HTGRR=1.8110);
per-tree cyc1 heights now within ±0.08ft of the oracle (tn2 DF 3.80 vs 3.81, tn15 GF 5.10 vs 5.10 — WAS
+6.9/+7.0ft off). No regression: multicycle 339/11, OP tests pass, OC .sum NINT-neutral. Test
test_oc_regent_smtree + fixtures/oregoncoast/.

⚠ FINDINGS: (1) fixing the seedling heights did NOT resolve the OC .sum ~1% multi-cycle drift — the
earlier "over-tall seedlings perturb CCH → large trees short" hypothesis is REFUTED. The ~1% drift is a
SEPARATE LARGE-tree ORGANON height residual (large trees ~0.15ft short by cyc2), independent of the
small-tree model. (2) HCOR=0 in the wrapper: only GF calibrated here (CON=1.0189, ~+1.9%), a follow-on.
(3) DBH-from-height (HTDBH Curtis/Arney, htdbh.f CURARN coeffs in hand) NOT yet wired — small-tree DBH
still grows via the large-tree DDS (tn15 dbh 0.1 vs oracle 0.5); low .sum-impact but needed for faithful
multi-cycle recruitment. NEXT OC chunks: (a) HTDBH DBH-from-height; (b) HCOR ratio-estimator calibration;
(c) the SEPARATE large-tree ORGANON height residual (the actual .sum-drift driver); (d) wire Olympic too.

## 2026-08-20 — OC multi-cycle ~1% drift RESOLVED-CORNERED (measured: HG_SWO CCH tie-break precision)

After the small-tree height fix, the residual OC .sum ~1% multi-cycle drift was root-caused by measuring
per-tree increments (oracle .trl vs jl cycle_hook), NOT inferred:
  • LARGE-tree DBH is BIT-EXACT at cyc2 (19/22 exact, 3 at ±0.05 print-rounding) ⇒ ORGANON DG_SWO
    diameter growth is faithful; the deterministic DG inputs (D, BAL-from-DBH-rank) match.
  • LARGE-tree HEIGHT increment (1995→2000) is MIXED-SIGN ±0.1-0.2ft per tree (tn8 +0.12, tn23 −0.21,
    tn1 −0.16 … tn11 +0.06), netting −0.4% — a scatter, NOT a systematic bias.
  • DBH bit-exact + height mixed-sign scatter ⇒ the divergence is isolated to HG_SWO's height-specific
    input CCH (crown competition at height), which depends on the tree HEIGHT-RANKING — a tie-break-
    sensitive quantity. Tiny Float32 height-ordering ties flip CCH ⇒ mixed-sign height scatter; DBH
    growth (DBH-rank/BAL, stable) stays bit-exact. This is EXACTLY the accepted cornered class the
    western cluster carries (goal-doc item 1: "per-tree diffs MIXED-SIGN and mostly-cancelling = the
    accepted RDPSRT/AVHT40 tie-break precision compounding. No fix warranted"). Here it's the HG_SWO
    CCH height-rank analogue; the ~1% .sum drift is this scatter compounding over cycles (then TPA
    feedback). CORNERED, meets the bar. No fix warranted (a fix would require bit-matching FVS's exact
    height-tie ordering, a precision artifact, not a model gap).

⇒ **OC multi-cycle is now bit-exact-or-cornered**: cyc0 bit-exact · small-tree heights bit-exact (regent
port) · large-tree DBH bit-exact · large-tree height = accepted CCH tie-break-precision cornered class.
Both ORGANON variants (OC + OP) are multi-cycle bit-exact-or-cornered. Remaining OC items are optional
FAITHFULNESS refinements that do NOT change the cornered .sum verdict: small-tree DBH-from-height (HTDBH,
tn15 dbh 0.1 vs 0.5 — few, low-volume) + HCOR small-tree height calibration (only GF, ~+1.9%) + wiring
Olympic's small-tree path. The FFE extension remains USER-gated (like NC FFE).

## 2026-08-20 — OP (Olympic) small-tree path is ALREADY FAITHFUL (measured; "wire Olympic" CLOSED)

Measured OP opt01 stand-1 (= same S248112 stand as OC) seedling heights+DBH at cyc1 vs the live FVSop_clean
TREELIST: jl tn2 DF (0.574, 7.398) vs oracle (0.6, 7.4); jl tn15 GF (0.321, 4.848) vs oracle (0.3, 4.9) —
BOTH MATCH (height AND DBH). So OP does NOT have the OC seedling-over-growth bug: OP's op_htg_default (op/
htgf.f) grows these seedlings faithfully (the OP/NWO oracle uses the large-tree-style htgf for them — h 2.0→
7.4 — NOT the regent-suppressed path OC/SWO uses, h 2.0→3.8; same stand, different variant small-tree model,
and jl matches EACH). ⇒ NO Olympic small-tree port needed — the "wire Olympic" follow-on is CLOSED. Both
OC (regent-ported) and OP (op_htg_default) small-tree HEIGHTS are now validated-correct vs their oracles.
The only remaining OC small-tree residual is the DBH-from-height calibration coupling (low-pri, .sum-cornered).
