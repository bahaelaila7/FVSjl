# FVSjl

An idiomatic, maintainable, thread-safe Julia reimplementation of the USFS
**Forest Vegetation Simulator** — Southern (SN) variant first, structured so other
variants (CS/NE/LS) plug in cleanly.

It is a **drop-in replacement** for the Fortran `FVSsn`: same `.key`/`.tre` inputs,
same SQLite / `.sum` outputs, and — in the default `faithful=true` mode —
**bit-exact** results, validated against the faithful port at `/workspace/FVSjulia`
and the live Fortran build.

This is the *modernized* successor to that faithful 1:1 port. Where the port was a
direct transliteration of Fortran (globals + gotos), FVSjl is built around:

- **explicit state, no globals** — one `StandState` passed in, so stands run on
  separate threads with no contention;
- **pure numeric kernels + stateful orchestration** — easy to test;
- **Structure-of-Arrays tree data + preallocated scratch** — no hotpath allocation,
  autovectorizable;
- **variants via multiple dispatch** — add a variant without touching the engine.

See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md), [`docs/PORTING.md`](docs/PORTING.md),
and [`docs/DIVERGENCES.md`](docs/DIVERGENCES.md).

## Status

The SN natural-process core (growth, mortality, density, volume, regeneration/establishment, stump
sprouting), the management/disturbance keywords (all thinning methods + modifiers, YARDLOSS, FIX*/MULT
multipliers, event monitor), structural-stage classification, ECON, modern IO (YAML/CSV ↔ .key/.tre),
and the **full FFE fire/fuels/carbon extension** (both CARBCALC methods, crown-lift, the `.out` CARBREPT
report, and all 9 FFE DBS tables) are ported and validated against the live Fortran oracle. Suite:
**4392 pass + 10 `@test_broken`** (the `@test_broken` track the one known non-bit-exact residual — the
FFE dead-pool one-cycle crown-lift timing lag).

See `docs/DECISION_FLOW.md` (+ the interactive `docs/decision_flow*.html` call-graphs),
`docs/REMAINING_WORK.md`, and `docs/TOLERANCE_AND_COVERAGE_AUDIT.md` for the per-subsystem status,
what remains, and the test-tolerance taxonomy. Each chunk is validated against the oracle before the next.

## Run

```julia
using FVSjl
FVSjl.main(["--keywordfile=/path/to/stand.key"])   # (CLI wired up in a later chunk)
```

## Test

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
