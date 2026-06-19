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

## Status / coverage

- **Leg 1+2:** 12/12 scenarios — FVSjulia matches live Fortran byte-for-byte
  (baseline, site lo/hi, site-species LP/YP, ecounit M221/232, forest 808,
  3/20 cycles, simulated fire, THINBTA).
- **Leg 3:** 12/12 — FVSjl cycle-0 stand summary matches the Fortran `.sum` row 1
  (init / NOTRE / stand statistics validated across the matrix).

Leg 3 is limited to cycle 0 until FVSjl gains the volume + `.sum` writer (C5) and
fire (C7); once those land this harness extends to the full multi-cycle `.sum` and
the fire/DB tables with no change to the oracle-validation flow.
