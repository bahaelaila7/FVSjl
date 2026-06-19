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

Under active construction, chunk by chunk; each chunk is validated against the
oracle before the next begins. See the project plan for the chunk breakdown.

## Run

```julia
using FVSjl
FVSjl.main(["--keywordfile=/path/to/stand.key"])   # (CLI wired up in a later chunk)
```

## Test

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
