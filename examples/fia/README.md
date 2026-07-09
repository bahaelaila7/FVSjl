# Example — FVS-ready FIA plots → standalone stand files

This folder shows how to turn **FVS-ready FIA plots (identified by their `STAND_CN`)** into
self-contained stand files you can run, convert, and share — with **no database** needed to
run them. The tool is [`bin/fvsjl-fia-export.jl`](../../bin/fvsjl-fia-export.jl); the full
reference is [docs/TOOLS.md](../../docs/TOOLS.md).

## What's shipped here (no database required)

One real FVS-ready plot, `STAND_CN = 163384065010854` (an oak/hickory stand, ecoregion 221Ja),
exported to **both** input forms so you can inspect the format:

| file | what it is |
|------|------------|
| `163384065010854.key`  | legacy keyword file — `STDINFO`/`DESIGN`/`SITECODE`/`INVYEAR`/`GROWTH` materialized from `FVS_STANDINIT_COND` |
| `163384065010854.tre`  | legacy fixed-column tree records (from `FVS_TREEINIT_COND`) |
| `163384065010854.yaml` | the same stand as modern keyword YAML (self-describing — carries `variant: SN`) |
| `163384065010854.csv`  | the same trees as a named-column CSV |
| `163384065010854.fvsjl.sum` | the `.sum` from running the exported stand |

Run the shipped stand directly (it reads the `.tre`/`.csv` beside it):

```bash
julia --project=. bin/fvsjl-run.jl examples/fia/163384065010854.key    # legacy
julia --project=. bin/fvsjl-run.jl examples/fia/163384065010854.yaml   # modern, no --variant needed
```

## Generate your own from a database

Point the export tool at an FVS-ready FIA SQLite database (the `FVS_STANDINIT_COND` /
`FVS_TREEINIT_COND` tables — e.g. `SQLITE_FIADB_ENTIRE.db`) and give it CNs:

```bash
# one CN → out/<CN>.key + out/<CN>.tre
julia --project=. bin/fvsjl-fia-export.jl SQLITE_FIADB_ENTIRE.db 163384065010854 out/

# a file of CNs (one per line) → modern yaml/csv, with a fidelity check vs the DB reader
julia --project=. bin/fvsjl-fia-export.jl SQLITE_FIADB_ENTIRE.db @cns.txt out/ --format yaml --validate
```

Or run the guided demo, which does export → run → convert → validate:

```bash
bash examples/fia/export_and_run.sh                                  # shows the shipped stand only
bash examples/fia/export_and_run.sh SQLITE_FIADB_ENTIRE.db 163384065010854,…   # full round trip
```

## Fidelity note

The tree list and direct-measurement stand fields are reproduced **exactly**. The one thing
an `STDINFO` card cannot carry is the plot's `LATITUDE`/`LONGITUDE` (there are no such card
fields), so a materialized stand falls back to the national-forest-table default — a
Hopkins-index difference that can move the crown competition factor (`CCF`) by ~1 at cycle 0.
`--validate` reports this precisely. For **byte-exact** FIA reproduction, run the CN through
the `DATABASE` keyword block instead (see [docs/TOOLS.md §3](../../docs/TOOLS.md#3-run-directly-from-fvs-ready-fia-plots-cns));
for **portable** files, use the export.
