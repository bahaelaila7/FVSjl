# FVSjl command-line tools & workflows

FVSjl ships three CLI drivers under `bin/`. They cover the whole lifecycle of a stand:
get it in (from files or from FVS-ready FIA plots), run it, and move it between the legacy
and modern input forms. The engine reads **either** form directly, so conversion is only
for interoperability (feeding stock FVS, or modernizing a legacy stand) — never a required
step.

| tool | what it does |
|------|--------------|
| [`bin/fvsjl-run.jl`](#1-run-a-stand--fvsjl-run) | run a stand (`.key`/`.yaml` + `.tre`/`.csv`) and write its `.sum`/CSV summary |
| [`bin/fvsjl-translate.jl`](#2-convert-between-forms--fvsjl-translate) | convert `.key ↔ .yaml` and `.tre ↔ .csv`, both ways, losslessly |
| [`bin/fvsjl-fia-export.jl`](#4-export-fvs-ready-fia-plots-to-files--fvsjl-fia-export) | export FVS-ready **FIA plots (by CN)** to standalone stand files |

Two supporting references: **[FORMATS.md](FORMATS.md)** (the file formats — both YAML
flavors, the tree `.tre`/`.csv`, the `.sum`) and **[KEYWORDS.md](KEYWORDS.md)** (every
keyword). Runnable examples are under [`../examples/`](../examples/) (see
[`examples/convert_and_run.sh`](../examples/convert_and_run.sh) and
[`examples/fia/`](../examples/fia/)).

Every command below assumes you are in the project root and prefixes Julia with
`julia --project=.` (abbreviated `$JL`).

---

## The four input files

A stand is **keywords** (what to run) + **trees** (the inventory). Each comes in a legacy
fixed-column form and a modern readable form; the four are interchangeable:

| role | legacy | modern | who reads it |
|------|--------|--------|--------------|
| keywords | `stand.key` | `stand.yaml` | `run_keyfile`, `each_stand` (auto-detect by extension) |
| trees | `stand.tre` | `stand.csv` | found beside the keyword file by base name; `read_tree_records` auto-detects |

The tree file is located by **base name**: `run_keyfile("x.yaml")` reads `x.csv` (or `x.tre`)
sitting next to it. So a modern run is `x.yaml` + `x.csv`; a legacy run is `x.key` + `x.tre`;
mixing (`x.yaml` + `x.tre`) also works.

> **Does the driver understand YAML/CSV directly?** **Yes.** `_keyword_reader` dispatches on
> the extension — `.yaml`/`.yml` → the YAML reader, else the `.key` lexer — and
> `read_tree_records` dispatches `.csv` → the CSV reader, else `.tre`. So
> `run_keyfile("stand.yaml")` and `bin/fvsjl-run.jl stand.yaml` run a modern stand with **no
> conversion step**. The YAML also carries its own `variant:` / `output_format:`, so a
> `.yaml` is self-describing (a `.key` is not — it needs `--variant` at run time).

---

## 1. Run a stand — `fvsjl-run`

```bash
$JL bin/fvsjl-run.jl <stand.{key,yaml}> [--variant SN|NE|CS|LS] [--output sum|csv] [-o outfile]
```

Runs the stand and writes its summary to stdout (or `-o` file). A `.key`/`.tre` carries
neither a variant nor an output preference, so those are flags; a `.yaml`'s `variant:`/
`output_format:` supply them when the flag is omitted (an explicit flag always wins).

```bash
# legacy .key (reads thinba.tre beside it), .sum to stdout — the defaults:
$JL bin/fvsjl-run.jl examples/thinba/thinba.key

# modern: YAML keywords + CSV trees, no conversion, modern CSV summary to a file:
$JL bin/fvsjl-run.jl examples/semantic/thinsdi.yaml --output csv -o thinsdi.sum.csv

# a .key with no variant marker, run as Northeast (or CS / LS):
$JL bin/fvsjl-run.jl net01.key --variant NE
```

The library entry point is `run_keyfile(path; variant=…, output=…)`; a multi-stand file
yields one `.sum` block per stand.

---

## 2. Convert between forms — `fvsjl-translate`

```bash
$JL bin/fvsjl-translate.jl <src> <dst> [tree-format] [--flat]
```

Direction is inferred from the extensions. The round-trip is **lossless** (a re-converted
file simulates byte-identically):

```bash
# keywords: legacy .key  ->  order-aware hierarchical YAML   (and back)
$JL bin/fvsjl-translate.jl examples/thinba/thinba.key   thinba.yaml
$JL bin/fvsjl-translate.jl thinba.yaml                  thinba.key      # --flat for a flat list

# trees:    legacy .tre   ->  named-column CSV               (and back)
$JL bin/fvsjl-translate.jl examples/thinba/thinba.tre   thinba.csv
$JL bin/fvsjl-translate.jl thinba.csv                   thinba.tre

# the SEMANTIC YAML (format: fvs-stand/v1) unravels to a .key for stock FVS:
$JL bin/fvsjl-translate.jl examples/semantic/thinsdi.yaml  thinsdi.key
```

`tree-format` (3rd arg) overrides the `.tre` FORTRAN FORMAT for `.tre ↔ .csv` on a stand
whose `TREEFMT` is non-default (default is the SN layout). The two YAML *flavors* — the
order-preserving *keyword-stream* form and the declarative *semantic* (`fvs-stand/v1`) form
— are both read directly and both convert to `.key`; see [FORMATS.md](FORMATS.md).

> **Does the conversion tool need to be run before a stand can be used?** **No.** The engine
> reads `.yaml`/`.csv` natively (§1). Use `fvsjl-translate` only to hand a modern stand to
> **stock FVS** (which needs `.key`/`.tre`), or to modernize a legacy stand for readability.

---

## 3. Run directly from FVS-ready FIA plots (CNs)

FVSjl's `DATABASE`/`DSNin` keyword block reads a stand + tree list straight out of an
FVS-ready FIA SQLite database at run time — no intermediate files. This is the authoritative
path (it consumes FVS's own FVS-ready tables). Point `DSNin` at the database and give the two
`StandSQL`/`TreeSQL` queries; `%StandID%` is filled from `STDIDENT`:

```text
STDIDENT
232388261010854
DATABASE
DSNin
/path/to/SQLITE_FIADB_ENTIRE.db
StandSQL
SELECT * FROM FVS_STANDINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
TreeSQL
SELECT * FROM FVS_TREEINIT_COND WHERE STAND_CN = '%StandID%'
EndSQL
END
NUMCYCLE  10
ECHOSUM
PROCESS
STOP
```

Run it like any other keyfile: `$JL bin/fvsjl-run.jl fromcn.key`. To run *many* CNs this way
without hand-writing a query per CN — or to get **portable files** that no longer need the
database — use the export tool below.

---

## 4. Export FVS-ready FIA plots to files — `fvsjl-fia-export`

```bash
$JL bin/fvsjl-fia-export.jl <fia.db> <CN | @cnfile | CN1,CN2,…> [outdir] \
      [--variant SN|NE|CS|LS] [--format key|yaml] [--numcycle N] [--validate]
```

Given an FVS-ready FIA SQLite database and one or more stand **CNs**, it writes for **each
CN** a self-contained stand that runs with **no database**:

- `<outdir>/<CN>.key` + `<CN>.tre`  — legacy fixed-column (**default**)
- `<outdir>/<CN>.yaml` + `<CN>.csv` — modern readable (`--format yaml`)

The setup cards (`STDIDENT`/`STDINFO`/`DESIGN`/`SITECODE`/`INVYEAR`/`GROWTH`/`NUMCYCLE`) are
materialized from the `FVS_STANDINIT_COND` columns — the same fields the native reader
consumes — and the tree records come straight from `FVS_TREEINIT_COND`. The output is
portable: run it with `fvsjl-run`, convert it with `fvsjl-translate`, share it, archive it.

```bash
# one CN → a standalone legacy stand you can run anywhere:
$JL bin/fvsjl-fia-export.jl SQLITE_FIADB_ENTIRE.db 163384065010854 out/
$JL bin/fvsjl-run.jl out/163384065010854.key            # runs with no DB

# a list of CNs → modern yaml/csv, self-describing variant, with a fidelity check:
$JL bin/fvsjl-fia-export.jl SQLITE_FIADB_ENTIRE.db @cns.txt out/ --format yaml --validate

# a comma-separated list, run as Central States, 20 cycles:
$JL bin/fvsjl-fia-export.jl fia.db 100,101,102 out/ --variant CS --numcycle 20
```

**Performance.** The FVS-ready master (e.g. `SQLITE_FIADB_ENTIRE.db`, ~1.5 M rows) is **not
indexed on `STAND_CN`**, so the tool first builds a small **in-memory indexed working set**
of just the requested CNs (one scan of the master, then O(log n) lookups). The master is
opened **read-only** (`mode=ro&immutable=1`) and never modified.

**Fidelity (`--validate`).** For each export the tool also runs the same CN via the
`DATABASE` reader (§3) and compares the `.sum`. The tree list and the direct-measurement
fields are **exact**; the one thing an `STDINFO` card cannot carry is the stand's
`LATITUDE`/`LONGITUDE` (there are no such card fields), so the reconstructed stand falls back
to the national-forest-table default — a Hopkins-index difference that can nudge the crown
competition factor (`CCF`) by ~1 at cycle 0. `--validate` reports exactly that:

```
validate 163384065010854: cycle-0 faithful, differs only in [CCF 114.0≠113.0, …]
   — crown-model ULP (STDINFO can't carry the DB lat/long); later cycles = model ULP
```

i.e. the exported stand reproduces the initial inventory faithfully (bit-exact trees, ±1-unit
crown ULP); multi-cycle divergence past that is the model's own dense-phase single-precision
tail, identical to running the `DATABASE` reader against live FVS. For **byte-exact** FIA
reproduction, use the `DATABASE` reader form (§3) directly; for **portable** files, use the
export.

A worked, runnable demo is in [`examples/fia/export_and_run.sh`](../examples/fia/export_and_run.sh).

---

## Cheat-sheet

```bash
# RUN
$JL bin/fvsjl-run.jl stand.key                          # legacy
$JL bin/fvsjl-run.jl stand.yaml                          # modern (self-describing)
$JL bin/fvsjl-run.jl stand.key --variant NE --output csv -o s.csv

# CONVERT (lossless, both ways)
$JL bin/fvsjl-translate.jl stand.key  stand.yaml         # keywords .key <-> .yaml
$JL bin/fvsjl-translate.jl stand.tre  stand.csv          # trees    .tre <-> .csv

# FROM FIA CNs
$JL bin/fvsjl-fia-export.jl fia.db 163384065010854 out/                 # -> out/<CN>.key + .tre
$JL bin/fvsjl-fia-export.jl fia.db @cns.txt out/ --format yaml --validate  # -> .yaml + .csv
$JL bin/fvsjl-run.jl out/163384065010854.key            # then run the exported stand
```
