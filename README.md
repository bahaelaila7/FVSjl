# FVSjl

An idiomatic, thread-safe Julia reimplementation of the USFS **Forest Vegetation
Simulator (FVS)** — a **drop-in replacement** for the live Fortran FVS. It reads the
same `.key` / `.tre` inputs and writes the same `.sum` / SQLite outputs, and in the
default `faithful = true` mode it is **bit-exact** to the Fortran (barring only
single-precision ULP and a few documented numerical divergences).

Two variants are validated, bit-exact drop-ins for their live Fortran counterparts:

| Variant | Replaces | Tag |
|---------|----------|-----|
| **Southern (SN)** | `FVSsn` | `FVSsn-complete` |
| **Northeast (NE)** | `FVSne` | `FVSne-complete` |

The natural-process core (diameter/height growth, mortality, density, crown, volume,
regeneration/establishment, stump sprouting), the management & disturbance keywords (all
thinning methods + modifiers, harvest, fertilization, the event monitor), the full **FFE
fire / fuels / carbon** extension, the **ECON** extension, and **modern readable I/O**
(YAML ⇄ `.key`, CSV ⇄ `.tre`) are all ported and validated against the live Fortran oracle.

## Install & run

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Run a stand from the command line (the CLI picks the output format and variant):

```bash
# Southern is the default for a .key; pass --variant NE to run as Northeast.
julia --project bin/fvsjl-run.jl  stand.key                 # → stand.sum (SN)
julia --project bin/fvsjl-run.jl  stand.key  --variant NE   # run as Northeast
julia --project bin/fvsjl-run.jl  stand.yaml --output csv -o out.csv
```

A YAML stand carries its own `variant:`; a `.key` defaults to SN unless you pass
`--variant`. Or call the library directly:

```julia
using FVSjl
txt = run_keyfile("stand.key"; variant = FVSjl.Northeast())   # returns the .sum text
```

Convert between the legacy fixed-column forms and the modern readable forms (lossless,
either direction — inferred from the extensions):

```bash
julia --project bin/fvsjl-translate.jl  stand.key  stand.yaml   # .key ⇄ .yaml
julia --project bin/fvsjl-translate.jl  stand.tre  stand.csv    # .tre ⇄ .csv
```

Worked examples (thinning, multi-stand, multi-scenario, semantic YAML) are in
[`examples/`](examples/).

## Documentation

**Using FVSjl**

- **[docs/KEYWORDS.md](docs/KEYWORDS.md)** — every keyword FVSjl recognizes, grouped by
  purpose, with its named parameters in both the legacy `.key` and modern `.yaml` forms.
- **[docs/FORMATS.md](docs/FORMATS.md)** — the input formats: `.key`/`.yaml` (two YAML
  flavors — an order-preserving keyword stream *and* a declarative `fvs-stand/v1`) and
  `.tre`/`.csv`, plus how the keyword and tree files pair up and the translate CLI.

**Understanding the model**

- **[docs/DECISION_FLOW.md](docs/DECISION_FLOW.md)** — a bird's-eye map of how a run
  flows from input to output (the major routines and the conditions that gate them).
- **[docs/decision_flow_fvsjl.html](docs/decision_flow_fvsjl.html)** — interactive
  call-graph of **FVSjl** (both variants + FFE/econ; click a node to expand its callees,
  hover for the code excerpt; shape = scope, fill = port status). Open in a browser.
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

Current suite: **5391 pass / 2 broken** (the two broken are an accepted SN COMPRESS
eigensolver divergence and a NOHTDREG ULP residual — both documented). Both variants are
validated against the live Fortran builds; SN is additionally checked against the 1:1
faithful port at `/workspace/FVSjulia` (there is no faithful-port oracle for NE — it is
validated directly against live `FVSne`).
