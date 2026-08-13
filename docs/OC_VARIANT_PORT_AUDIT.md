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
