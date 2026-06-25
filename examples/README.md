# FVSjl examples — Southern variant thinning, in legacy *and* modern forms

Each example is provided in **four equivalent files**:

| form | keywords | trees |
|------|----------|-------|
| legacy (stock FVS) | `*.key` (fixed-column) | `*.tre` (fixed-column) |
| modern (FVSjl)     | `*.yaml` (readable)    | `*.csv` (named columns) |

The `.yaml`/`.csv` are produced from the `.key`/`.tre` with the translator and are
**losslessly equivalent** — they parse to the identical keyword stream and tree list
(verified: `read_keywords_yaml(x.yaml)` ≡ `read_keyfile_records(x.key)`, and
`read_trees_csv(x.csv)` ≡ `read_tree_records(x.tre)`). Convert either way:

```bash
julia --project bin/fvsjl-translate.jl sn_thinba.key  sn_thinba.yaml   # key  → yaml
julia --project bin/fvsjl-translate.jl sn_thinba.yaml sn_thinba.key    # yaml → key
julia --project bin/fvsjl-translate.jl sn_thinba.tre  sn_thinba.csv    # tre  → csv
```

Run the legacy form through the engine (reads `<stem>.tre` automatically):

```julia
using FVSjl
print(run_keyfile("examples/sn_thinba.key"))   # .sum text to stdout (+ .sum / .db files)
```

## The two examples

Both start from the same SN stand (FIA-coded inventory in the `.tre`/`.csv`), a 5-cycle
projection from 1990, and schedule one thinning in **2010**.

### `sn_thinba` — THINBBA (thin from below to a residual basal area)
`THINBBA  2010  80` → thin from below until basal area is reduced to **80 sq ft/ac**.
Fields: `1`=year/cycle, `2`=residual BA, `3`=cut efficiency, `4`/`5`=min/max DBH,
`6`=species (0=all), `7`=plot. In the run this removes ~296 TPA at 2010.

(FVS has no keyword literally named `THINBA`; basal-area thinning is `THINBBA`
"from below" and `THINABA` "from above" — swap the keyword to thin from above.)

### `sn_thinsdi` — THINSDI (thin to a residual Stand Density Index)
`THINSDI  2010  200` → thin until the Stand Density Index is reduced to **200**.
Same field layout, with field `2` = residual SDI. In the run this removes ~104 TPA at
2010 — a *different* removal than the BA target, since SDI weights large stems
differently than basal area.

In the `.sum` output the cut shows up as non-zero **removed** columns
(removed TPA / cuft / bdft) on the 2010 row.
