# FVSjl Architecture

A from-scratch, idiomatic Julia reimplementation of the Forest Vegetation
Simulator. It is a drop-in replacement for the live Fortran FVS — the **Southern
(`FVSsn`)**, **Northeast (`FVSne`)**, **Central States (`FVScs`)**, and **Lake States
(`FVSls`)** variants — and, in the default `faithful=true` mode, bit-exact to it.

This document explains the *shape* of the code. For how to translate a new
subroutine or variant see [PORTING.md](PORTING.md); for places where we knowingly
differ from Fortran see [DIVERGENCES.md](DIVERGENCES.md).

## The core idea: explicit state, no globals

Fortran FVS keeps ~52 COMMON blocks (~1,300 globals) mutated by every subroutine.
That is the root cause of the faithful port being unreadable and single-threaded.

FVSjl replaces every COMMON block with a plain `mutable struct`, all composed onto
a single `StandState{V}` (see `src/core/state.jl`) that is passed as the first
argument to nearly every function:

| COMMON block(s)              | Struct (field on StandState) |
|------------------------------|------------------------------|
| `CONTRL` / `CONCHR`          | `Control`   (`.control`)     |
| `ARRAYS`                     | `TreeList`  (`.trees`)       |
| `PLOT` / `PLTCHR`            | `PlotData`  (`.plot`)        |
| `PLTCHR` species tables + coeffs | `SpeciesData` (`.species`) |
| `CALCOM` / `HTCAL`           | `Calibration` (`.calib`)     |
| `PDEN`                       | `Density`   (`.density`)     |
| `OUTCOM` / `SUMTAB`          | `OutputState` (`.out`)       |
| `WORKCM` + per-tree work cols| `Scratch`   (`.scratch`)     |
| `RANCOM` / `ESRNCM`          | `FVSRng`    (`.rng`)         |
| `ESCOMN` / `ESCOM2` / ...    | `Establishment` (`.estab`)   |
| `DBSCOM`                     | `DbsState`  (`.dbs`)         |
| `FMCOM` / `FMFCOM` / ...     | `FireState` (`.fire`, lazy)  |
| `ECNCOM` / ...               | `EconState` (`.econ`, lazy)  |

Consequences:
- **No globals** (requirement #1).
- **One state per stand/thread** → no shared mutable state → stands run on
  separate threads with zero contention (requirement #5).
- **Testable** — construct a state, call a function, assert (requirement #6).

## Two-layer functions: pure kernels + stateful orchestration

Numeric formulas are written as **pure functions** of their inputs (coefficients +
a few scalars), returning values — independently testable, allocation-free, and
trivially autovectorizable. The **orchestration** layer loops over trees, calls
kernels, and writes results back into `StandState`.

```julia
# pure kernel (testable in isolation)
@inline ln_dds(c, dbh, ba, ...) = c.intercept + c.ldbh*log(dbh) + ...

# stateful orchestration
function diameter_growth!(s::StandState, ::Southern)
    t = s.trees
    @inbounds for i in 1:t.n
        s.scratch.wk[2, i] = ln_dds(coef(s, t.species[i]), t.dbh[i], s.plot.basal_area, ...)
    end
end
```

## Data layout (requirements #3, #4)

`TreeList` is a Structure-of-Arrays: each attribute is a `Vector` preallocated to
`MAXTRE`, with `n` active records. This is cache-friendly and lets LLVM
auto-vectorize the tight per-tree loops. **We do not** hand-vectorize branchy
per-tree logic — that fights both readability and bit-exactness. Scratch space
lives in `Scratch` so the cycle loop never allocates.

## Variants (requirement #10)

`AbstractVariant` singletons (`Southern`, `Northeast`, and later others such as
`CentralStates` / `LakeStates`) are the type parameter of `StandState`, so variant
hooks (`diameter_growth!`, `mortality!`, ...) dispatch at zero cost and devirtualize
(trim-friendly). The base `engine/` is variant-agnostic and calls the generic
hooks declared in `src/variants/variant.jl`; each variant supplies methods under
`src/variants/<variant>/` (e.g. `southern/`, `northeast/`). Adding a variant never
edits the engine. **Four are complete and validated: SN (`FVSsn`), NE (`FVSne`), CS
(`FVScs`), and LS (`FVSls`).** NE carries its own structurally-new pieces (R9 Clark
volume + BAL-potential competition); CS and LS share an SN-family ln(DDS) diameter-growth
model (`centralstates/` + `lakestates/`) plus LS-specific site/volume/FFE/sprout coefs.
Shared `engine/io/core` code gates variant behavior on the variant/coefficients rather
than hardening to any one, so each added variant kept the earlier ones bit-exact.

## Faithful mode & divergences (requirements #8, #9)

`control.faithful` (default `true`) selects bit-exact Fortran behaviour. Suspected
FVS bugs are fixed only behind `!faithful`, marked in-code with `# DIVERGENCE:` and
logged in [DIVERGENCES.md](DIVERGENCES.md). Float arithmetic is `Float32` (FVS is
`REAL*4`); `nint` rounds half-away-from-zero to match Fortran `NINT`. The RNG
(`FVSRng`) reproduces the Park–Miller streams exactly and is swappable later.

## Trim-friendliness (size, deferred)

Modernization does not shrink the Julia runtime image; a small standalone binary
needs `juliac --trim` (Julia ≥1.12). That requires exactly this style: type-stable,
no globals, no reflection, devirtualized dispatch. We write to be trim-ready and
verify the small-binary pass as later follow-up work, not a current gate.

## Testing

For **Southern**, `test/oracle/oracle.jl` runs any `.key` through **Oracle A** (the
1:1 faithful port at `/workspace/FVSjulia`, always available) and through FVSjl, then
diffs `.sum`/DB output with Float32 tolerance; **Oracle B** (live Fortran `FVSsn`) is
the final confirmation. For **Northeast** there is *no* faithful-port oracle (no
FVSjulia NE) — it is validated directly against **live `FVSne`** (relinked from the
build objects; see `test/harness/ne_oracle.sh`). Per-kernel unit tests live in
`test/unit/`, full-run parity in `test/integration/`.
