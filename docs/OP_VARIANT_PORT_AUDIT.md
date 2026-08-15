# OP (Olympic) ORGANON variant port — audit

Stream 2 [#208], the thin ORGANON follow-on to OC (Oregon Coast). Branch: kt-variant-port worktree.
Oracle: `/workspace/.opwork/FVSop_clean` (+ `FVSop_dbg` for instrumented dumps), built with
`/workspace/.opwork/build_g16.sh` (mirrors `.ocwork/relink_oc.sh`: gfortran-16 relink of the
prebuilt `FVSop_buildDir/*.o` + the shared `.crwork/isoc23_shim.o`).

## Oracle build

- `FVSop_buildDir` already carries 616 prebuilt `*.o`. `build_g16.sh clean` relinks `FVSop_clean`
  (production-equivalent). `build_g16.sh dbg <swap.o>` relinks `FVSop_dbg` swapping one instrumented
  object (the swap object must be named identically to the buildDir object it replaces, e.g.
  `diagro.o`).
- Reference `.sum` generated from the live relink via the inline keyfile `opdbg.key` (stand S248112,
  the same stand as OC's ocmin, NUMCYCLE=1). Growth columns match `tests/FVSop/opt01.sum.save` case
  1 (UNTHINNED CONTROL) exactly: 1990 `536 77 184 63 5.1`, 1995 `499 98 219 71 6.0`. Only MERCH/BdFt
  volume differs (opt01.sum.save is a STALE 20260413 build; the relink is 20260401) — the growth path
  is bit-identical across builds, so the relinked oracle is the authoritative target (per the EC
  `.sum.save`-is-stale lesson).
- Like OC, a full multi-cycle run later SIGSEGVs in the FVS-native volume path (fvsvol.f) — a
  VOLUME-stage crash that does NOT block growth validation.

## MEASURED OP DGF structure (not assumed)

Instrumented the live `FVSop_dbg` (scoped DEBUG + a `WRITE` inside `DG_NWO_RUN`, buildDir restored
pristine afterward — marker count 0):

- **ORGANON version = 2 (NWO / Northwest Oregon).** `grohed.f` header = "ORGANON NWO&SMC";
  `sitset.f` defaults `IMODTY=2` and sets `VERSION=IMODTY`; the run logs "MODEL TYPE NOT RECOGNIZED.
  BEING SET TO 2= NORTHWEST OREGON". (SMC=3 keyword-selectable, RAP=4 exists — neither is the default.)
  This is the key OC↔OP difference: OC is SWO (VERSION=1), OP is NWO (VERSION=2).
- **MAXSP = 39** (op/orgspc.f species order; PRGPRM.F77), vs OC's 50.
- **Species map:** op/orgspc.f `OSPMAP` (FVS 1..39 → ORGANON FIA code) composed with
  organon/execute2.f `SPGROUP_RUN` `SCODE2` (VERSION 2/3: FIA code → NWO group 1..11 =
  DF,GF,WH,RC,PY,MD,BL,WO,RA,PD,WI). Ported verbatim in `species.jl`.
- **Reuse verdict:** the shared ORGANON *engine structure* (organon/diagro.f `DIAMGRO_RUN` driver,
  `GET_BAL`/`SSTATS` stand stats, `SUBMAX`, the ln(DDS) form, `CRADJ`, `DGCALIB`, `PREPARE`) is
  version-agnostic and REUSABLE. What NWO needs on top is DATA-only: the `DG_NWO_RUN` `DGPAR(11,11)`
  coefficient table + the NWO species-group `ADJ` multiplier + the OP species map. The OC Julia port
  hard-coded VERSION=1 and raises on VERSION≠1, so OP carries its NWO data alongside. NOT a pure
  clone of OC-SWO — the coefficient tables genuinely differ (Zumrawi & Hann 1993 for DF/GF vs the SWO
  Hann & Hanus 2002 set) — but the port cost is data, not new machinery.

## Validated chunks (this session)

### DG_NWO diameter-growth core — BIT-EXACT
`src/variants/olympic/organon_diamgro_nwo.jl` (`op_dg_nwo` / `op_diamgro_run_nwo`) +
`src/variants/olympic/species.jl` (`op_spgroup_nwo`).

Validation: instrumented `DG_NWO_RUN` to dump `ISPGRP,DBH,CR,SITE,SBAL1,SBA1,DG` as `Z8.8` hex
(lossless Float32 round-trip), stand S248112, 73 tree-calls. Re-computing DG in Julia from the
hex-exact inputs (via the registered `FVSjl.op_dg_nwo`):

    73 records   BIT-EXACT = 72   1-ULP = 1   worst-ULP = 1

The single 1-ULP tree is the documented irreducible libm `exp`/`log`/`pow` intrinsic difference
(gfortran-16 vs Julia), consistent with the whole-codebase ULP floor. Growth is deterministic
(DGSD=0 on the ORGANON path), so this is a genuine bit-exact target, not an RNG straddle. Earlier
decimal-input dump gave 49/73 exact + 24 within 1 ULP; switching to hex inputs collapsed those 24 to
exact, proving the residual was debug-print truncation, not a port error.

### Foundation
`src/variants/olympic/olympic.jl` — `Olympic <: AbstractVariant` (VARACD "OP", MAXSP=39,
`organon_version`=2, htg_period=5), registered in `variant.jl` (`variant_from_code("OP")`) and
`FVSjl.jl` includes. `site_setup!` carries the op/sitset.f Nigh(1995) DF(16)↔WH(19) site-index
conversion (analogous to OC's DF↔PP). `using FVSjl` loads clean and precompiles.

## buildDir-pristine confirmation

`FVSop_buildDir/diagro.f`: marker count 0, `diff` vs `organon/diagro.f` = identical. `diagro.o`
untouched (Jun 4 timestamp; swaps used `.opwork/diagro.o`). Clean relink from the restored tree
succeeds.

## ORGANON NWO growth ENGINE — VALIDATED bit-exact-or-1-ULP (this session)

`src/variants/olympic/organon_nwo.jl` ports the full NWO growth engine (C1 marshalling, C2 PREPARE,
C3 DG, C4 HG, C5 crown, C6 mortality, C7 execute), reusing OC's version-agnostic helpers
(`OrganonBuffer`, `oc_get_bal`, `oc_get_ccfl`, `oc_caltst`, `oc_quad1`, `oc_spmix`, `oc_dforty`,
`oc_oldgrowth`, `oc_oldgro`, `oc_pm_fert`, `oc_hg_fert`, `oc_hg_thin`) + the NWO coefficient tables
and NWO structural differences (all cited `organon/*.f`).

**Validation.** The live `FVSop_clean` was run on `opdbg.key` (stand S248112, DEBUG 1 / DGDRIV HTGF
CROWN MORTS) → `/workspace/.opwork/opdbg.out`, the per-tree ORGANON dump. Feeding the engine the
exact `/ORGANON/` FOR EXECUTE buffer (27 trees) + SI_1=93.5, SI_2=83.17, ACALIB(1,1)=0.795231164,
ACALIB(2,1)=0.667298734, MSDI=950 reproduces (test `test/unit/test_op_organon_nwo.jl`, 38 assertions):

| Chunk | quantity | result |
|---|---|---|
| **C2 PREPARE** (VERSION=2) | ACALIB(1,1) height calib | **0.795231164 — Δ=0.0 (bit-exact)** |
| | ACALIB(2,1) crown calib | **0.667298734 — Δ=0.0 (bit-exact)** |
| **C3 DG_NWO** | DGRO, 7 IORG=1 (DF) trees | **7/7 bit-exact (maxΔ=0)** |
| **C4 HG_NWO** (B_HG+HG_NWO) | HGRO, all 27 buffer trees | **25/27 bit-exact, 2×1-ULP (maxΔ=9.5e-7)** |
| **C5 CROWGRO** (NWO CALIB(2)) | CR2, all 27 buffer trees | **18/27 bit-exact, rest ≤2-ULP (maxΔ=1.2e-7)** |
| **C6 PM_NWO** | MORTEXP/DEADEXP, all 27 | **27/27 bit-exact (maxΔ=0)** |

ORGANON is DGSD=0 (deterministic) ⇒ genuine bit-exact targets; the HGRO/CR2 residuals are the
documented irreducible gfortran-16↔Julia libm exp/log/pow ULP (crown chains exp/log through
HCB→PCR→CALIB2, accumulating the ULP — same as OC). Key measured NWO facts: SI_1=98 (DF), SI_2=87.67
(WH, via the DF→WH −0.432+0.899·SI conversion); MSDI=950; **big-6 stand gate counts only GF(3)/DF(16),
valid ORGANON species = {3,16,18,19,21,22,23,28,33,34,37}** (op/dgdriv.f) — so in S248112 only the 7
DF trees are ORGANON-grown, the 20 LP/PP/SP/ES/WF trees grow FVS-native (see below). PM_NWO has a
DISTINCT per-group form (DF: √DBH & CR^0.25; GF: BAL/DBH). CROWGRO applies CALIB(2) (SWO does not).

## What remains for the full end-to-end cyc0 `.sum` (NOT ORGANON-specific)

The ORGANON NWO engine is complete + validated, but the reference `.sum` needs three more subsystems,
none ORGANON and all large (measured, not assumed):

1. **FVS-native Wykoff DGF/HTGF/crown for the 20 non-ORGANON trees.** In S248112 only 7 of 27 trees
   are valid-ORGANON (DF); LP/PP/SP/ES/WF grow via op/dgf.f (Wykoff LN(DDS)) + op/htgf.f — a
   separate variant-scale port (OC's chunk C9 analogue). The `.sum` growth columns cannot be
   bit-exact without it. `op_build_organon_buffer!` + `op_execute_nwo` handle the ORGANON side; the
   StandState growth hook (`diameter_growth!(::Olympic)`) and species/blockdata loader are NOT yet
   wired (`coefficients(::Olympic)` still errors) pending this.
2. **Ecoclass site-index fan** (op/sitset.f ECOCLS/SICHG/HTCALC) — plant association `CHS133` →
   DF SITEAR(16)=98 → the per-species fan. ORGANON only consumes SITE_1(DF)/SITE_2(WH), but the
   FVS-native trees + volume need the full fan.
3. **R6 volume** (`compute_volumes_op!`) — OP has no dedicated volume routine; it uses the shared
   Region-6 Behre/Flewelling path (like SO/WC/PN chunk 8). The `.sum` TCuFt/MCuFt/BdFt columns
   (1990: 1472/972/5003; 1995: 2256/1745/9249) need it.

`op_prepare_nwo`'s internal HT dub for a NON-valid tree (e.g. tree20 LP → 49.84) is intentionally NOT
written back (IORG=0 trees get FVS HTDBH → 62.39, item 1); only the ACALIB it produces is used, and
that is bit-exact.

## buildDir-pristine confirmation (this session)

`FVSop_buildDir/diagro.f`: marker count 0, `diff` vs `organon/diagro.f` = identical. Only the ORACLE
was run (writes to `/workspace/.opwork/opdbg.out`); no Fortran/buildDir object was modified.

Do NOT merge until the orchestrator re-runs and verifies.

## CHUNK 1 — FVS-native Wykoff DGF/HTGF for the non-ORGANON species — VALIDATED per-tree

Ports the FVS-native large-tree growth (op/dgf.f + op/htgf.f + op/findag.f + op/htcalc.f) for the
20 non-ORGANON trees (op is ORGANON only for {3,16,18,19,21,22,23,28,33,34,37}; everything else —
LP/PP/SP/ES/WF — grows FVS-native). Files:
- `src/variants/olympic/diameter_growth.jl` — `op_dgcon` (ENTRY DGCONS), `op_dg_dgdsq`, `op_dgf_dds`
  (DEFAULT ln(DDS)), `op_dgcons!`/`dgf!(::Olympic)`. op/dgf.f MAPSPC (39→20 groups) + all DATA tables.
- `src/variants/olympic/height_growth.jl` — `op_htcalc` (PN edition: DF/WO King, SS/RC Farr, Curtis
  "misc" EXCLUDES 16:18), `op_findag` (findag.f byte-identical to WC), `op_htg_default` (5-YEAR step,
  AGP10=SITAGE+5), `height_growth!(::Olympic)`. Reuses `op_bratio` (organon_nwo.jl).
- `data/olympic/{species_coefficients,species_translation}.csv` — op/blkdat.f 39-species table +
  op/bratio.f bark ⇒ `coefficients(::Olympic)` now returns a valid struct (was: errored).

**Validation — RUN vs live FVSop_clean, stand S248112, cyc0** (`/workspace/.opwork/opdbg.out`, DEBUG
DGF/HTGF dump). Feeding the exact oracle per-tree inputs (ELEV=7, SLOPE=0.30, ASPECT=5.49779, IFOR=6,
BA=62.5335, AVH=63.4388; SITEAR WF=97.98189/ES=139.15092/LP=98.22826/SP=PP=139.15092):

| Chunk | quantity | result |
|---|---|---|
| **DGCONS** | DGCON(ISPC) for WF/ES/LP/SP/PP | **5/5 bit-exact to the F9.5 9030 print** |
| **DGF DEFAULT** | LN(DDS), 19 non-ORGANON trees | **17/19 bit-exact to the F7.4 9001 print; 2 sub-0.1" trees input-print-limited** |
| **HTGF DEFAULT** | HTG(I), 16 trees with a format-901 dump | **16/16 bit-exact to the F9.2 print** |

Test `test/unit/test_op_native_growth.jl` (41 assertions, all pass; wired into `test/runtests.jl`).
The two sub-0.1" DGF trees are DBH=0.1 input records that DGDRIV back-dates to an internally-computed
diameter the DGF DEBUG dump prints only to F11.4 (=0.0828); for a sub-0.1" tree ln(D) is hyper-
sensitive to that 5th digit, so they are INPUT-print-limited (the ORGANON decimal-vs-hex lesson), not
a formula error — the other 17 are exact at the oracle's full print precision. op is DGSD>0 (OLDRN
serial-corr on this path) so multi-cycle is straddle-class; the cyc0 DDS is a deterministic bit-exact
target.

**Growth hook.** `dgf!(::Olympic)` and `height_growth!(::Olympic)` are defined and dispatch; both
coexist with the ORGANON-NWO engine via the inline op/dgdriv.f IORG gate (sp∈valid ∧ HT>4.5 ∧
DBH≥0.1 ⇒ ORGANON, skipped in the FVS-native path). `StandState(Olympic())` now constructs.

**End-to-end opt01 `.sum` — STILL BLOCKED** (documented, not reached). Reaching the cyc0 `.sum`
needs the OTHER chunks, none ported yet:
1. **Site-index fan** (op/sitset.f ECOCLS/SICHG/HTCALC) — the DGF/HTGF above CONSUME `SITEAR(ISPC)`
   (WF=97.98189, ES/SP/PP=139.15092, LP=98.22826, DF=98.00) which the fan produces from the ecoclass
   `CHS133`. The port fed these as measured inputs; a live run needs the fan + a forest_idx/`IFOR=6`
   loader.
2. **DGDRIV calibration + ORGANON coexistence driver** — no OP method wires the FVS-native `dgf!`
   into the shared DGDRIV calibration loop AND the ORGANON `op_build_organon_buffer!`/`op_execute_nwo`
   into one per-cycle growth step (the generic `diameter_growth!(::AbstractVariant)` calls only `dgf!`).
3. **R6 volume** (`compute_volumes_op!`) — the `.sum` TCuFt/MCuFt/BdFt columns.
4. **Crown + mortality** wiring for the FVS-native trees.

So the delivered bar is the stated fallback: **cyc0 large-tree DGF/HTGF bit-exact-or-input-print-limited
per-tree vs the live oracle**, plus `coefficients(::Olympic)` no longer erroring and both growth hooks
dispatching.
