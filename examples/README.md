# FVSjl examples (Southern variant)

Each example lives in its **own folder** and is provided in several equivalent forms:

```
thinba/      thinba.key   thinba.tre   thinba.yaml   thinba.csv
thinsdi/     thinsdi.key  thinsdi.tre  thinsdi.yaml  thinsdi.csv
multistand/  multistand.key  …  .tre  …  .yaml  …  .csv
semantic/      thinsdi.yaml  multistand.yaml   (+ .tre/.csv)   ← the SEMANTIC form
multiscenario/ stand.yaml    (+ stand.tre/.csv)              ← N scenarios, one stand
fia/           163384065010854.{key,tre,yaml,csv}           ← a real FVS-ready FIA plot, exported
convert_and_run.sh                            ← run + convert every form
fia/export_and_run.sh                         ← FIA CNs → standalone files, then run/convert
```

> **From FVS-ready FIA plots.** The `fia/` folder ships one real FIA plot (`STAND_CN
> 163384065010854`) exported to standalone `.key`/`.tre` + `.yaml`/`.csv` by
> [`bin/fvsjl-fia-export.jl`](../bin/fvsjl-fia-export.jl) — runnable with no database. To turn
> your own list of CNs into files, see [`fia/README.md`](fia/README.md) and the CLI guide
> [`../docs/TOOLS.md`](../docs/TOOLS.md).

> **Multi-scenario runs.** `multiscenario/stand.yaml` projects one stand under four
> treatments (control / light thin / heavy thin / thin-twice) — every scenario reads the
> same `stand.csv`/`stand.tre`. The semantic converter emits a `REWIND 2` before each
> scenario after the first so stock FVS re-reads that inventory; FVSjl re-reads implicitly.
> Verified: the converted `.key` reproduces **live FVSsn** within the ±1-cuft single-precision
> tail, and all input forms give byte-identical FVSjl output.

| file | what it is |
|------|------------|
| `*.key` | legacy fixed-column keyword file (**the runnable form**) |
| `*.tre` | legacy fixed-column tree records |
| `*.yaml` | keyword file — **two flavors**: order-aware *keyword-stream* (here) and *semantic* (`semantic/`) |
| `*.csv` | tree records as a named-column table |
| `*.fvsjl.sum` | the `.sum` summary from **FVSjl**, run from the `.yaml` + `.csv` |
| `*.fvssn.sum` | the `.sum` summary from **live FVSsn**, run from the equivalent `.key` + `.tre` |

> The two `.sum` files in each folder are committed side by side so you can see FVSjl and
> the original Fortran produce the **same output** — they agree to within the ±1-cuft
> single-precision tail (compare the data rows; the `-999` header timestamp differs per
> run). The column-by-column `.sum` format is documented in
> [`../docs/FORMATS.md`](../docs/FORMATS.md#4-the-sum-summary-output).

> **Two YAML flavors.** The `thinba/thinsdi/multistand` folders hold the **keyword-stream**
> YAML (an order-preserving image of the `.key`, documented below). The `semantic/` folder
> holds the **semantic** YAML (`format: fvs-stand/v1`) — a stand described by *intent*
> (`invyr`, `numcycle`, `treatments`, `treelist`) that the converter unravels into keyword
> order. Full reference for both, the tree CSV/TREEFMT, and the species codes:
> [`../docs/FORMATS.md`](../docs/FORMATS.md); every keyword + parameter:
> [`../docs/KEYWORDS.md`](../docs/KEYWORDS.md).

## Run & convert in one go

```bash
bash examples/convert_and_run.sh      # runs the examples + converts .key↔.yaml and .tre↔.csv
```

## Running an example

From the `FVSjl` directory (the `.tre` next to the `.key` is read automatically):

```bash
julia --project=. -e 'using FVSjl; print(run_keyfile("examples/thinba/thinba.key"))'
```

That prints the `.sum` summary to stdout. Save it instead with:

```bash
julia --project=. -e 'using FVSjl; write("examples/thinba/thinba.sum", run_keyfile("examples/thinba/thinba.key"))'
```

`multistand/` holds three stands (control / THINBBA / THINSDI) in one keyword file,
separated by `PROCESS` — running it produces three stand blocks in the `.sum`.

## The structured YAML form — order-aware hierarchy

Each keyword is a block whose parameters are **named**, and values keep their natural
type (**numbers unquoted; code strings quoted**). The whole keyword stream stays an
**ordered sequence** (order is significant in FVS); the YAML only *groups* it into
named sections for readability. The shape is a `stand:` map → an ordered list of
section blocks → ordered keyword entries:

```yaml
stand:
  - setup:                       # ── grouped, but still in input order ──
      - STDINFO:
          forest_code: 80106
          habitat: "231Dd"       # a code string → quoted
          stand_age: 60          # numbers → unquoted
      - NUMCYCLE: { cycles: 6 }
  - treatments:                  # ── ORDER MATTERS within a section ──
      - THINBBA:
          year: 2010
          residual_basal_area: 80
```

Flattening the blocks top-to-bottom reproduces the **exact** keyword order — the
grouping never reorders. A keyword with no schema, or a free-form continuation line,
falls back to `keyword:`/`params:` or a verbatim `raw:` entry (still lossless). See
[`../docs/KEYWORDS.md`](../docs/KEYWORDS.md) for the section list, the order-significant
relationships (SPGROUP→THIN, same-cycle order, COMPUTE def-before-use), and worked
side-by-side examples.

Regenerate the YAML from any `.key` with either tool (both emit the hierarchical form;
add `--flat` for the legacy single `keywords:` list):

```bash
julia --project=. bin/fvsjl-translate.jl examples/thinba/thinba.key examples/thinba/thinba.yaml
python3 examples/key_to_structured_yaml.py examples/thinba/thinba.key > examples/thinba/thinba.yaml
```

### Parameter schema (the keywords used here)

A full reference for **every** supported keyword (purpose + named parameters) is in
[`../docs/KEYWORDS.md`](../docs/KEYWORDS.md).


| keyword | named parameters (in order) |
|---------|------------------------------|
| `STDIDENT` | `id` (the stand label on the following line) |
| `DESIGN` | `basal_area_factor`, `fixed_plot_area_inverse`, `break_dbh`, `number_of_plots`, `nonstockable_code`, `sample_weight`, `stockable_proportion` |
| `STDINFO` | `forest_code`, `habitat`, `stand_age`, `aspect`, `slope`, `elevation`, `stand_origin` |
| `SITECODE` | `site_species`, `site_index` |
| `INVYEAR` | `year` |
| `NUMCYCLE` | `cycles` |
| `THINBBA` | `year`, `residual_basal_area`, `cut_efficiency`, `dbh_min`, `dbh_max`, `species`, `plot` |
| `THINSDI` | `year`, `residual_sdi`, `cut_efficiency`, `dbh_min`, `dbh_max`, `species`, `plot` |
| `TREEFMT` | `format` (the Fortran FORMAT string) |
| `TREEDATA` / `PROCESS` / `STOP` | *(no parameters)* |

(Only present fields are emitted — e.g. `DESIGN 11.0 1.0` lists just `number_of_plots`
and `nonstockable_code`.)

## The two thinning keywords

- **`THINBBA`** — thin from below to a residual **basal area** (`residual_basal_area`).
  `THINBBA  2010  80` cuts to 80 sq ft/ac in 2010 (≈ −296 TPA in this stand).
  (FVS has no keyword named `THINBA`; basal-area thinning is `THINBBA` "from below" /
  `THINABA` "from above".)
- **`THINSDI`** — thin to a residual **Stand Density Index** (`residual_sdi`).
  `THINSDI  2010  200` cuts to SDI 200 (≈ −104 TPA) — a different removal than the BA
  target, since SDI weights large stems differently than basal area.

Both fire on the **2010** row of the `.sum`, shown by the non-zero *removed* columns
(removed TPA / cuft / bdft).
