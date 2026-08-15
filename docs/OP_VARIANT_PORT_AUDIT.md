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

## What remains for the next OP chunk (mirrors OC C1-C10)

The DG_NWO *equation* is validated in isolation; a full end-to-end cyc0 `.sum` beachhead still needs
the shared ORGANON engine wired for NWO (mostly reuse of OC's version-agnostic plumbing):

1. **ORGANON marshalling (C1):** `org_intree` FVS→ORGANON buffer for MAXSP=39 (species → FIA →
   ISPGRP into TDATAI(,2)); IORG/big-6 classification.
2. **PREPARE setup (C2):** the VERSION=2 branch of organon/prepare.f (ACALIB/TMPCAL, HT/CR dubbing,
   the SI_1/SI_2 marshalling — SITE feeds DG_NWO's `ln(SITE)`; WH group 3 uses SITE_2).
3. **Stand-stat plumbing:** reuse OC's version-agnostic `SSTATS`/`GET_BAL`/`SUBMAX` (DBH²·EXPAN) to
   produce the `sbal1`/`sba1` that DG_NWO consumes — validated here via the oracle's values.
4. **DGDRIV COR calibration:** the FVS-native DGDRIV self-calibration from observed treelist DG
   (fort.16 showed nonzero COR, e.g. ISPC=2 COR=0.0338) — shared with OC's `dg_cor`.
5. **Height / crown / mortality / volume:** NWO branches of htgrowth.f (HG_NWO), crngrow.f,
   mortality.f (op/morts.f has explicit NWO tables), orgvol.f — each a bit-exact chunk.

Do NOT merge until the orchestrator re-runs and verifies.
