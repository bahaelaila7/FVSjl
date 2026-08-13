# OC/OP ORGANON coupling — port plan (task #208, Stream 2)

**Bottom line up front.** OC (Oregon Coast / SW Oregon) and OP (Olympic) are NOT Wykoff-DDS
variants. Their per-cycle growth engine is **ORGANON** — the `organon/` (13,040 lines) +
`vorganon/` (1,331 lines) subsystem that the FVSjl port does **not** yet have. MEASURED from a
live `FVSoc_clean` DEBUG-DGF run: the `.out` header reads `VERSION FS2026.1 ORGANON SWO`, the DGF
path loads the FVS tree list into ORGANON, `CRATET ORGANON ERROR CODE, CYCLE= 1 IERROR= 0` is
emitted, and diameter/height/crown come back through the ORGANON `ACALIB/TMPCAL` calibration
arrays. There is effectively **no FVS-native large-tree DGF to validate independently** for the
18 ORGANON-eligible species — DGF *is* ORGANON. This is the real cost driver of the OC/OP port.

All citations are `file:line` from `/workspace/ForestVegetationSimulator`.

## 1. Model overview — the coupling architecture

A single logical `LORGANON` (`common/ORGANON.F77:6`) gates all coupling. It is **hardcoded
`.TRUE.`** in `oc/grinit.f:342` — no keyword disables it in normal OC operation. A per-tree flag
`IORG(I)` (0/1, `ORGANON.F77:78`) marks each record as a "valid ORGANON tree". The single ORGANON
growth entry per cycle is `CALL EXECUTE(...)` at `oc/dgdriv.f:370-377` (`organon/execute2.f:24`);
the single one-time setup entry is `CALL PREPARE(...)` at `oc/cratet.f:269` (`organon/prepare.f:31`).

**Two gates decide who ORGANON grows:**
- **Per-tree eligibility** (`oc/dgdriv.f:218-248`): `IORG(I)=1` iff `HT>4.5 AND DBH>=0.1` AND the
  FVS species index ∈ `{2,3,4,7,8,16,18,24,27,30,31,34,36,37,38,39,42,46}` (18 of 50 species).
- **Big-6 stand gate** (`oc/dgdriv.f:217,254-261`): the ORGANON engine only runs if the stand has
  ≥1 "big-6" tree (`{2,4,7,16,18}` = IC, GF/WF, DF, SP, PP). If `NBIG6==0`, `IORG` is forced 0 for
  all trees, `MORTEXP=0`, `EXECUTE` is skipped (`GO TO 261`), and the stand reverts to FVS-native.

**Which component comes from where (per normal cycle, big-6 present):**

| Component | Valid ORGANON tree (`IORG=1`) | Non-valid / no big-6 |
|---|---|---|
| Diameter | ORGANON `DGRO(I)` → bark-adj DDS overwrites FVS `WK2(I)`; FVS serial-corr suppressed (`OLDRN=0,FRM=1`), so ORGANON DG is deterministic (`oc/dgdriv.f:444-553`) | `oc/dgf.f` FVS Wykoff LN(DDS) + `DGSCOR` error |
| Height | ORGANON `HGRO(I)` → `HTG(I)`, native HTGF bypassed (`oc/htgf.f:95-101`) | FVS-native HTGF (`oc/htgf.f:105+`) |
| Crown | ORGANON `CR2(I)` → `CRNEW(I)`, native crown bypassed (`oc/crown.f:281-288`) | FVS-native crown |
| Mortality | ORGANON `MORTEXP(I)*(FINT/5)` for the **whole stand** whenever `SMORMT>0`; FVS SDI/background skipped (`oc/morts.f:498-504,573,750`) | FVS-native SDI/background mortality |
| CR/HT dubbing (setup) | ORGANON `PREPARE` (`oc/cratet.f:269`) | FVS dubbing |
| Regen/ingrowth | new seedlings forced `IORG=0` → FVS-native (`oc/regent.f:186`) | FVS-native |

**Control defaults (MEASURED, `oc/grinit.f`):** `VARACD='OC'`, `MAXSP=50`, seed 55329,
`FINT=FINTH=FINTM=5` (5-yr cycle — ORGANON's native step), `LZEIDE=.FALSE.` (Stage SDI), `DGSD=0`
(no DG serial-corr), `LHTDRG=.FALSE.` all species. ORGANON `VERSION=1` (SWO) hardcoded
(`oc/grinit.f:345`); software edition = 9.1 (`execute2.f:1093`). Tripling is OFF by default
(`INDS(5)=0`, `oc/grinit.f:360`). ORGANON volumes OFF by default (`LORGVOLS=.FALSE.`,
`oc/grinit.f:343`) — FVS volume used unless keyword `ORGVOLS`.

## 2. Data flow FVS ↔ ORGANON

**FVS-facing boundary = the `/ORGANON/` common block** (`common/ORGANON.F77:28-35`), a fixed
2000-record structure. FVS fills inputs before `EXECUTE`: `SPECIES,DBH1,HT1OR,CR1,SCR1B,EXPAN1,
MGEXP,USER` (`oc/dgdriv.f:268-285`). ORGANON returns `DGRO,HGRO,CRCHNG,SCRCHNG,MORTEXP,NTREES2,
DBH2,HT2OR,CR2,SCR2B,EXPAN2,OCC,OAHT`.

**Internal ORGANON tree structure = argument-passed arrays, NOT common blocks.** Inside
`EXECUTE`/`GROW`: `TDATAI(2000,3)` (int: FIA species, species-group, user code), `TDATAR(2000,8)`
(real: DBH,HT,CR,EXPAN,…), results in `GROWTH(2000,4)` (`execute2.f:221-228,429-430`). This is the
seam a Julia port targets: fill `TDATAI/TDATAR` from `StandState.trees`, run the ORGANON growth
sequence, copy `GROWTH` back.

**Copy-back into FVS:** DBH `DGRO→WK2→DG` (`oc/dgdriv.f:446-557`); HT `HGRO→HTG` (`oc/htgf.f:96`);
CR `CR2→CRNEW` (`oc/crown.f:282`); mort `MORTEXP→WK2` (`oc/morts.f:499`).

## 3. `orgspc.f` species map (the 50→18 collapse)

`ORGSPC(INSPEC,OUTSPC)` (`oc/orgspc.f:1`) maps FVS index → ORGANON FIA code via the 50-element
`DATA OSPMAP` (`:55-65`). Every FVS species maps to *some* FIA code (so ORGANON gets correct stand
density even for surrogate trees), but only **18 map to a valid ORGANON code and get ORGANON
growth**: 2=IC(081) 3=RC(242) 4=GF(017) 7=DF(202) 8=WH(263) 16=SP(117) 18=PP(122) 24=PY(231)
27=CY(805) 30=WO(815) 31=BO(818) 34=BM(312) 36=RA(351) 37=MA(361) 38=GC(431) 39=DG(492) 42=TO(631)
46=WI(920). The other 32 pass to ORGANON under a surrogate FIA code but grow FVS-native.

## 4. ORGANON routine inventory (grouped)

**Growth-critical core to port (~9,200 lines):**

| File | Lines | Role |
|---|---|---|
| `organon/execute2.f` | 1275 | `EXECUTE` per-cycle driver (marshals FVS↔ORGANON; EDIT_RUN/SSTATS/GROW); `GET_ORGRUN_EDITION`=9.1 |
| `organon/grow.f` | 220 | `GROW` — sequences diameter→height→mortality→crown |
| `organon/diagro.f` | 688 | `DIAMGRO_RUN`+`DG_SWO/NWO/SMC/RAP`, `DG_THIN/FERT`, `GET_BAL_RUN` — diameter growth |
| `organon/htgrowth.f` | 1028 | `HTGRO1/2`+`HG_SWO/NWO/SMC/RAP`, `HD_*`, `LIMIT`, `RAGEA/RAH40` — height growth |
| `organon/crngrow.f` | 1517 | `CROWGRO`+`HCB_*`, `MCW/LCW/CW_*`, `CALC_CC`, `CRNCLO` — crown recession/width/closure |
| `organon/mortality.f` | 800 | `MORTAL_RUN`+`PM_SWO/NWO/SMC/RAP`, `RAMORT`, `OLDGRO`, `QUAD1/2` — mortality |
| `organon/diamcal.f` | 632 | `DGCALIB`, `DIAMGRO`, `MORTAL`, `PHTS`, `DIB/DOB_*` bark, `SSUM`, `GET_BAL` |
| `organon/statsorg.f` | 135 | `SSTATS`, `HTFORTY`, `RASITE` — per-cycle stand stats |
| `organon/submax.f` | 140 | `SUBMAX` — max-SDI density limit for mortality |
| `organon/growth_mods.f` | 75 | `GG_MODS` (genetic worth), `SNC_MODS` (Swiss needle cast) modifiers |
| `organon/whphg.f` | 118 | `SITECV_F/SITEF_C/SITEF_SI` — W-hemlock site/potential-height |
| `organon/prepare.f` | 795 | **setup** `PREPARE`+`EDIT` — one-time HT/CR imputation (from `oc/cratet.f:269`) |
| `organon/start2.f` | 1087 | **setup** `HDCALIB/CRCALIB/PRDHT/PRDCR`, `HD40_*` — start-of-run H-D & crown calibration |

**Conditional (tripling OFF by default in OC):** `tripleorg.f` (460), `orgtrip.f` (23).

**Interface (small):** `varget.f` (199), `varput.f` (199) control-buffer pack/unpack; `getorgv.f`
(29) Event-Monitor `OCC/OAHT`; `vorganon/orgtab.f` (291) settings-table print; `vorganon/orin.f`
(451) + `vorganon/org_intree.f` (589) native `.INP` treelist reader (**alternate input path — NOT
the FVS/FIA path; skip for the FVS-coupled port**).

**Deferrable output-side (~4,600 lines):** volume `orgvol.f` (451) + `vols.f` (612) + `voleqns.f`
(1361); woodquality `woodqual.f` (228) + `woodq2.f` (878). OFF by default (`LORGVOLS=.FALSE.`); FVS
volume is used. `orgfert.f` (90) is a **doc-only** note file (no compiled subroutine) — ignore.

## 5. Chunk decomposition (dependency-ordered)

Each chunk validated **bit-exact vs the live oracle** via the FVS DEBUG seams in §7. Reuse the
shared jl engine; the ORGANON step slots in behind the existing `diameter_growth!`/`height_growth!`
/`crown`/`mortality!` hooks on `OregonCoast()`, gated by the `IORG`/big-6 logic ported into
`oc/dgdriv.f`'s Julia analogue.

- **C0 — Foundation (DONE).** `OregonCoast <: AbstractVariant`, MAXSP=50 species table, grinit
  constants, registration, include, oracle relinked. See `docs/OC_VARIANT_PORT_AUDIT.md`.
- **C1 — orgspc species map + IORG/big-6 gate + `/ORGANON/` buffer marshalling.** Port
  `oc/orgspc.f` `OSPMAP`, the `oc/dgdriv.f:210-285` fill loop, and the two gates. No growth yet —
  validate the *tree list handed to ORGANON* against the `I,PTNO,SPECIES,DBH1,HT1OR,CR1,EXPAN1`
  DEBUG dump (already captured). Small, high-leverage.
- **C2 — Setup calibration (`prepare.f` + `start2.f`).** HT/CR imputation + H-D/crown start
  calibration (`ACALIB/TMPCAL`). Validate against the `ORGANON ACALIB(k,sp)/TMPCAL(k,sp)` dump and
  `CRATET ORGANON ERROR CODE`. This is where the calibration ratios (e.g. `TMPCAL(1,1)=0.789290`)
  are set; getting these bit-exact is prerequisite to any growth match.
- **C3 — Diameter growth (`diagro.f` DG_SWO + `diamcal.f` bark/BAL + `submax.f` + `statsorg.f`).**
  The largest single payoff. Validate per-tree `DGRO/BARK/DIAGR/DDS` from `oc/dgdriv.f:449` DEBUG.
- **C4 — Height growth (`htgrowth.f` HG_SWO).** Validate `HGRO/HTG` from `oc/htgf.f:97`.
- **C5 — Crown (`crngrow.f` CROWGRO/HCB/CW).** Validate `CR2/CRNEW` from `oc/crown.f:285`.
- **C6 — Mortality (`mortality.f` PM_SWO + RAMORT/OLDGRO).** Whole-stand `MORTEXP`; validate from
  `oc/morts.f:502`.
- **C7 — `GROW`/`EXECUTE` orchestration + `growth_mods.f`/`whphg.f`.** Wire C3–C6 into the
  per-cycle sequence; run full-cycle A/B vs `oct01.sum.save`.
- **C8 (deferred) — ORGANON volumes/woodquality** (`orgvol/vols/voleqns/woodqual/woodq2`) — only
  needed for the `ORGVOLS` keyword; FVS volume covers the default path.
- **C9 (deferred) — tripling** (`tripleorg.f`) — OFF by default; port only if a keyword enables it.

**SWO-first, NWO/SMC nearly free:** every version branch (`DG_NWO/SMC/RAP`, `HG_*`, `PM_*`,
`HD40_*`) sits alongside SWO in the same files. Porting OC=SWO lays all the scaffolding; OP just
selects `VERSION=2/3` and swaps the coefficient branch. See §6.

## 6. OP (Olympic) — commonality with OC

OP has an **identical file set** to OC (`blkdat,cratet,crown,dgdriv,dgf,grinit,grohed,htgf,morts,
orgspc,regent,sitset`) and **shares the exact same `organon/`+`vorganon/` engine** and `/ORGANON/`
common block. Differences: (a) `MAXSP=39` (`op/common/PRGPRM.F77:12`) vs 50, and a different
`orgspc.f` map; (b) ORGANON `VERSION` is **dynamic** — `op/grinit.f:343` default NWO, then
`op/sitset.f:70` sets `VERSION=IMODTY` (2=NWO, 3=SMC) from the model-type keyword (OC is fixed
SWO=1); banner "ORGANON NWO&SMC" (`op/grohed.f:34`). Same control style (`FINT=5`, `LZEIDE=.FALSE.`,
`DGSD=0`, `LHTDRG=.FALSE.`). **OP after OC ≈ the OC foundation chunk + the NWO/SMC coefficient
branches (already present in the ported ORGANON source) + version-selection in `sitset.f`.**

## 7. Validation strategy — the DEBUG seams

**No DEBUG hooks exist inside the ORGANON core** (`organon/`, `vorganon/` — no `DBCHK`). All
instrumentation is FVS-wrapper-side, gated by the standard `DEBUG` keyword + `CALL DBCHK`:
- `oc/dgdriv.f:116` `DBCHK(DEBUG,'DGDRIV',...)`; ORGANON inputs (`:282`), error codes (`:381-429`),
  per-tree `DGRO/BARK/DIAGR/DDS` (`:449`), `HGRO/CR2/MORTEXP/PROB` (`:460`).
- `oc/htgf.f:97` `HGRO/HTG`; `oc/crown.f:285` `CR2/CRNEW`; `oc/morts.f:502` `MORTEXP`;
  `oc/cratet.f:280+` `PREPARE` errors + dubbed CR/HT counts, plus the `ORGANON ACALIB/TMPCAL` dump.

**Per-chunk A/B protocol:** run `FVSoc_clean --keywordfile=<min>.key` with `DEBUG\nDGF\nEND`
(scope DEBUG to the relevant routine), diff the Julia per-tree output against the `.out` dump.
Prove real-vs-cornered by per-tree `DGRO`/`DDS`, never a `.sum` divergence (doctrine, cf #206).
Because `DGSD=0` and ORGANON serial-corr is suppressed, ORGANON DG is **deterministic** — no RNG
straddle on the ORGANON path, so bit-exact is a hard target (unlike the Wykoff variants' OLDRN
straddle). The FVS-native fallback path (32 species) reuses the shared Wykoff engine and *does*
carry the usual `DGSCOR` behaviour.

**Oracle:** `/workspace/.ocwork/FVSoc_clean` (relink via `relink_oc.sh`, gfortran-16 + isoc23
shim). Growth runs and the DGF dump is fully available; the run SIGSEGVs later in the FVS-native
volume path (`fvsvol.f:530` via `natcrs_ ← vols.f:319`) — a **volume-stage crash that does NOT
block growth validation**. See `docs/OC_VARIANT_PORT_AUDIT.md` for the crash detail and the
`build_g16`-style single-`.o`-swap instrumentation path.

## 8. Realistic effort estimate

The ORGANON growth core is ~9,200 lines of dense, edition-branched Fortran with heavy calibration
coupling (C2 must be bit-exact before C3–C6 can match). This is comparable in size to porting a
whole new mid-size variant *plus* a calibration subsystem — realistically **the largest single
chunk of unported western work**, on the order of the BC full port but with a stiffer bit-exact bar
(deterministic, no RNG cover). Recommended sequencing: C1→C2 first (marshalling + calibration
proven bit-exact) before committing to the four growth chunks. OC (SWO) is the pilot; OP (NWO/SMC)
is a thin follow-on once the engine + branch scaffolding exist.
