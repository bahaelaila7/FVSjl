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
