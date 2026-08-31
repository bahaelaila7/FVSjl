# FVSjl

An idiomatic, thread-safe Julia reimplementation of the USFS **Forest Vegetation
Simulator (FVS)** — a **drop-in replacement** for the live Fortran FVS. It reads the
same `.key` / `.tre` inputs and writes the same `.sum` / SQLite outputs, and in the
default `faithful = true` mode it is **bit-exact** to the Fortran (barring only
single-precision ULP and a set of documented, named numerical divergences).

The validation doctrine is **bit-exact-or-cornered vs the live relinked Fortran oracle**
(`FVS{v}_g16`, gfortran-16), measured per subsystem/chunk, never inferred. "Cornered"
means a divergence that was *measured* and reduced to a named floating-point primitive
(e.g. the `#206` OLDRN serial-correlation growth straddle, an RDPSRT self-thin tie-break,
a DGSCOR/volume ULP) — not an unexplained difference. The full method is
**[docs/DOCTRINE.md](docs/DOCTRINE.md)**.

## Validated variants (24)

All 24 FVS geographic variants are ported and validated as bit-exact-or-cornered
drop-ins for their live Fortran counterparts:

| Cluster | Variants |
|---|---|
| **Eastern** | Southern (SN), Northeast (NE), Central States (CS), Lake States (LS) |
| **Western Rockies** | Central Rockies (CR), Kootenai (KT), Inland Empire (IE), Eastern Montana (EM), Blue Mountains (BM), Teton (TT), Utah (UT), Central Idaho (CI) |
| **Pacific / coastal** | Central California (CA), East Cascades (EC), West Cascades (WC), West Sierra (WS), South Central Oregon (SO), Klamath/NC, Pacific Northwest (PN), Oregon Coast (OC, ORGANON), Olympic (OP, ORGANON) |
| **Other** | British Columbia (BC), Ontario (ON, Penner), Southeast Alaska (AK) |

Each variant has diameter/height growth, mortality, density, crown, and volume, and
most also carry FFE / ECON / mistletoe / Climate / establishment, validated per
subsystem against the freshly-relinked live oracle. Every species of a variant is
exercised by an all-species coverage test, and the reference multi-stand scenarios
(control · thinning · shelterwood · ECON · FFE fire · bare-ground planting) are
validated end-to-end against the live binary.

## Extensions

The natural-process core (diameter/height growth, mortality, density, crown, volume,
regeneration/establishment, stump sprouting) and the management & disturbance keywords
(all thinning methods + modifiers, harvest, fertilization, the Event Monitor) are joined
by the full extension suite — all ported and validated bit-exact-or-cornered vs the live
Fortran oracle:

- **FFE** fire / fuels / snags / **carbon**
- **ECON** economic analysis
- **FVS-Climate**
- **Dwarf mistletoe** (MISTOE)
- **Western Root Disease** (WRD)
- **Insect & beetle models** — Douglas-fir beetle (DFB), tussock moth (DFTM), mountain
  pine beetle (MPB), western spruce budworm (WSBWE), and the **WWPB** western pine beetle
  with the **PPE** parallel-processing landscape harness
- **COVER** canopy/shrub cover
- **DBS** database output tables (FFE, carbon, structure, mistletoe, root disease, …)
- **Establishment / regeneration** (AUTOES, PLANT, NATURAL, SPROUT, site prep, …)
- **Modern readable I/O** — YAML ⇄ `.key`, CSV ⇄ `.tre` (lossless, either direction)

## FIA behaviour-compat validation

Beyond the reference scenarios, FVSjl is validated against **FVS-ready FIA plots** at
full population scale — every stand projected the full horizon, all 10 `.sum` columns
per cycle vs freshly-relinked live FVS.

- **Eastern four — EXHAUSTIVE**: 1,468,376 stands, **99.99% bit-exact-or-cornered**; the
  residuals all classify to named cornered primitives, and the handful of `live_crash`
  cases are where **live FVS itself** SIGFPEs on extreme FIA geometry while FVSjl runs
  clean.
- **Western / non-eastern — in progress**: a full-population sweep of the 15 FVS-ready
  non-eastern variants (688,903 stands), run under the same per-cycle live-oracle
  comparison. It is being run **clean on final code** under the cap-and-fix discipline
  of [docs/DOCTRINE.md](docs/DOCTRINE.md) (pause at ~100 unexplained divergences, dig
  each to a named primitive or a real bug, fix upstream-first, resume). See
  [docs/PORT_STATUS.md](docs/PORT_STATUS.md) for the current per-variant state.

## Install & run

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Run a stand from the command line (the CLI picks the output format and variant):

```bash
# A .key defaults to Southern; pass --variant <CODE> for any other variant.
julia --project bin/fvsjl-run.jl  stand.key                 # → stand.sum (SN)
julia --project bin/fvsjl-run.jl  stand.key  --variant IE   # Inland Empire (any of the 24 codes)
julia --project bin/fvsjl-run.jl  stand.yaml --output csv -o out.csv
```

A YAML stand carries its own `variant:`; a `.key` defaults to SN unless you pass
`--variant`. Or call the library directly:

```julia
using FVSjl
txt = run_keyfile("stand.key"; variant = FVSjl.InlandEmpire())   # returns the .sum text
```

Convert between the legacy fixed-column forms and the modern readable forms (lossless,
either direction — inferred from the extensions):

```bash
julia --project bin/fvsjl-translate.jl  stand.key  stand.yaml   # .key ⇄ .yaml
julia --project bin/fvsjl-translate.jl  stand.tre  stand.csv    # .tre ⇄ .csv
```

Export **FVS-ready FIA plots** (by `STAND_CN`) from an FIA SQLite database to standalone
stand files that run with no database:

```bash
# one CN → out/<CN>.key + out/<CN>.tre  (or --format yaml for .yaml + .csv)
julia --project bin/fvsjl-fia-export.jl  fia.db  163384065010854  out/
julia --project bin/fvsjl-run.jl         out/163384065010854.key
```

Worked examples (thinning, multi-stand, multi-scenario, semantic YAML, FIA export) are in
[`examples/`](examples/). The three CLI tools and every workflow are in
[docs/TOOLS.md](docs/TOOLS.md).

## Documentation

**Current status & method**

- **[docs/PORT_STATUS.md](docs/PORT_STATUS.md)** — the canonical current port &
  validation status: every variant, extension, and the FIA full-population sweeps.
- **[docs/DOCTRINE.md](docs/DOCTRINE.md)** — the validation doctrine and the fix-execution
  method: bit-exact-or-cornered vs the live oracle, how to trace and fix a divergence,
  and the FIA full-population sweep discipline.

**Using FVSjl**

- **[docs/TOOLS.md](docs/TOOLS.md)** — the three command-line tools and every workflow: run a
  stand, convert `.key ↔ .yaml` / `.tre ↔ .csv`, and export FVS-ready FIA plots (by CN) to
  standalone files. Start here to *do* something.
- **[docs/KEYWORDS.md](docs/KEYWORDS.md)** — every keyword FVSjl recognizes, each with a full
  plain-English explanation of what it does, its parameters (position/units/defaults/codes),
  and a worked `.key` + `.yaml` example.
- **[docs/FORMATS.md](docs/FORMATS.md)** — the input formats: `.key`/`.yaml` (two YAML
  flavors — an order-preserving keyword stream *and* a declarative `fvs-stand/v1`) and
  `.tre`/`.csv`, plus how the keyword and tree files pair up, species groups, and variants.

**Understanding the model**

- **[docs/DECISION_FLOW.md](docs/DECISION_FLOW.md)** — a bird's-eye map of how a run
  flows from input to output (the major routines and the conditions that gate them).
- **[docs/decision_flow_fvsjl.html](docs/decision_flow_fvsjl.html)** — interactive
  call-graph of **FVSjl** (click a node to expand its callees, hover for the code excerpt;
  shape = scope, fill = port status). Open in a browser.
- **[docs/decision_flow.html](docs/decision_flow.html)** — the same view of the **FVS
  Fortran** semantics (the oracle FVSjl is validated against), for side-by-side reading.

**Design**

- **[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)** — how the code is shaped: explicit
  state (no globals), pure kernels + stateful orchestration, Structure-of-Arrays tree
  data, and variants via multiple dispatch.

## Why a rewrite

Where the Fortran keeps ~52 COMMON blocks (~1,300 globals) mutated by every subroutine,
FVSjl is built around:

- **explicit state, no globals** — one `StandState` passed in, so stands run on separate
  threads with zero contention;
- **pure numeric kernels + stateful orchestration** — easy to test in isolation;
- **Structure-of-Arrays tree data + preallocated scratch** — no hot-path allocation,
  autovectorizable;
- **variants via multiple dispatch** — add a variant without touching the engine.

## Test

```bash
julia --project=. test/runtests.jl
```

The hard cross-variant gate is `test/integration/test_multicycle.jl` — **339 pass / 11
broken, byte-identical** to the Fortran oracle (the 11 broken are the documented cornered
set: ULP-class / FVS-UB / eigensolver-tie / accepted named primitives). Every variant is
validated against its freshly-relinked live Fortran build; the per-variant and per-subsystem
detail is in [docs/PORT_STATUS.md](docs/PORT_STATUS.md).
