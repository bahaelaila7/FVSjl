# Three-way scenario harness (Fortran ↔ FVSjulia ↔ FVSjl)

Broader-coverage regression testing for the Southern variant. The guiding rule
(per the project owner): **never validate FVSjl against a possibly-faulty oracle.**
So every scenario is first confirmed `FVSjulia (Oracle A) == live Fortran FVSsn`
before any FVSjl comparison is trusted.

## Pieces

| file | role |
|------|------|
| `gen_scenarios.sh`  | derive a matrix of valid single-stand `.key` files from `snt01.key`, varying site index, site species, ecological unit / forest, cycle count, fire, and thinning. Copies the shared `snt01.tre`. |
| `validate_oracle.sh`| Leg 1+2 — run each scenario through the rebuilt Fortran `FVSsn` and FVSjulia, diff the full `.sum` (whitespace- and timestamp-insensitive). Proves the oracle is trustworthy. |
| `fvsjl_cycle0.jl`   | Leg 3 helper — initialize FVSjl on a key and print the cycle-0 stand summary (`TPA BA SDI CCF TopHt QMD`, per-acre quantities ÷ gross_space). |
| `three_way.sh`      | run all three legs and report. |

The Fortran ground truth is built on demand by FVSjulia's `tests/fortran_baseline.sh`
(rebuilds `/tmp/FVSsn_new` + a glibc shim from `bin/FVSsn_buildDir/*.o`).

## Run

```bash
bash test/harness/gen_scenarios.sh        # (re)generate scenarios/*.key + *.tre
bash test/harness/three_way.sh            # Fortran vs FVSjulia vs FVSjl
```

## Scenario coverage

- **Baseline (12):** site lo/hi, site-species LP/YP, ecounit M221/232, forest 808,
  3/20 cycles, fire, THINBTA — `gen_scenarios.sh`.
- **Species (36):** homogeneous stands of a broad SN species set (pines, oaks,
  hickories, maples, gums, ashes, cypress, cedar, elm, basswood…) — exercise
  **~29 distinct FIA forest types** (103, 141, 142, 161-168, 181, 501-520, 601-607,
  702-706, 801-809, 997) — `gen_species_scenarios.sh`.
- **Fire (3):** scheduled SIMFIRE at early/mid/late years with varied
  intensity/%-area.

- **Mixes (16):** oak-pine / conifer / hardwood compositions that reach the forest
  types homogeneous stands can't, so **all 8 DGF forest-type groups** (the
  `_dgf_forest_group` codepaths: lohd/nohd/okpn/sfhp/uphd/upok/ylpn/none) are hit.
- **SITSET (8):** vary SITECODE across the 9 site-index master-group representatives
  → each SITSET A/B/C/D site-index codepath.
- **SPCTRN (4):** foreign (non-SN) species codes (BF/JP/NS/QA) → the species-
  translation crosswalk codepath.

Goal is **codepath coverage**, not enumeration: every distinct forest-type-group
and species-resolution branch the engine takes is exercised and Fortran-validated.

## Status

- **Leg 1+2 (oracle):** baseline 12/12 + fire 3/3 — FVSjulia matches live Fortran
  byte-for-byte. (FVSjulia + FFE is slow; allow ~5 min/scenario.)
- **FORTYP port:** FVSjl `compute_forest_type!` matches the Fortran `.sum` FORTYP on
  **34/34** runnable species scenarios across the ~29 types.
- **Leg 3:** FVSjl cycle-0 stand summary matches the Fortran `.sum` row 1 on the
  baseline matrix (init / NOTRE / stand statistics).

Leg 3 is limited to cycle 0 until FVSjl gains the volume + `.sum` writer (C5) and
fire (C7); once those land this harness extends to the full multi-cycle `.sum` and
the fire/DB tables with no change to the oracle-validation flow.
