# FVSjl examples (Southern variant)

Each example lives in its **own folder** and is provided in several equivalent forms:

```
thinba/      thinba.key   thinba.tre   thinba.yaml   thinba.csv
thinsdi/     thinsdi.key  thinsdi.tre  thinsdi.yaml  thinsdi.csv
multistand/  multistand.key  …  .tre  …  .yaml  …  .csv
```

| file | what it is |
|------|------------|
| `*.key` | legacy fixed-column keyword file (**the runnable form**) |
| `*.tre` | legacy fixed-column tree records |
| `*.yaml` | **structured** keyword file — one block per keyword, *named* parameters |
| `*.csv` | tree records as a named-column table |

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

## The structured YAML form

Unlike a flat positional dump, each keyword is a block whose parameters are **named**,
and values keep their natural type — **numbers are numbers (unquoted); only genuine
strings are quoted**. Keyword order is preserved (it is significant in FVS), so the
file is a YAML *list*:

```yaml
keywords:
  - STDINFO:
      forest_code: 80106
      habitat: "231Dd"      # a code string → quoted
      stand_age: 60         # numbers → unquoted
      aspect: 315
      slope: 30
      elevation: 7
  - THINBBA:
      year: 2010
      residual_basal_area: 80
```

Regenerate the structured YAML from any `.key` with the bundled tool:

```bash
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
