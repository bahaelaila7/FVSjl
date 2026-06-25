# FVSjl keyword reference

Every keyword FVSjl recognizes, grouped by purpose, with its **named parameters**
(the names used by the structured-YAML form — see `examples/`). In a legacy `.key`
file parameters are positional fixed-column fields; the *field number* in parentheses
is the 1-based parameter slot (cols `10·n+1 … 10·n+10`). A blank date field defaults
to **cycle 1**; a date ≥ 1000 is a calendar year, < 1000 a cycle number.

Structured-YAML form:

```yaml
- THINBBA:
    year: 2010
    residual_basal_area: 80
```

---

## Stand identification & run setup

| keyword | purpose | named parameters |
|---|---|---|
| **STDIDENT** | Stand id/label. | `id` — the label, taken from the line that follows. |
| **STDINFO** | Stand attributes. | `forest_code`(1), `habitat`(2), `stand_age`(3), `aspect`(4, deg), `slope`(5, %), `elevation`(6, ×100 ft), `stand_origin`(9). |
| **DESIGN** | Plot/inventory design. | `basal_area_factor`(1), `fixed_plot_area_inverse`(2), `break_dbh`(3), `number_of_plots`(4), `nonstockable_code`(5), `sample_weight`(6), `stockable_proportion`(7). |
| **SITECODE** | Site index. | `site_species`(1, 0=all), `site_index`(2). |
| **INVYEAR** | Inventory (start) year. | `year`(1). |
| **NUMCYCLE** | Number of projection cycles. | `cycles`(1). |
| **TIMEINT** | Cycle length. | `cycle`(1, 0=all), `length`(2, yrs). |
| **CYCLEAT** | Insert an extra cycle boundary. | `year`(1). |
| **TREEFMT** | Custom tree-record FORTRAN FORMAT. | `format` — the FORMAT string on the following (up to two 80-col) lines. |
| **TREEDATA** | Read the tree records (`<stem>.tre`, or the self-describing `<stem>.csv`). | *(none)* |
| **NOTREES** | Start with no trees (bare stand). | *(none)* |
| **TFIXAREA** | Total fixed-plot area (small-tree expansion). | `area`(1). |
| **MANAGED** | Mark the stand as managed (affects some defaults). | `date`(1, optional). |
| **RESETAGE** | Reset stand age. | `date`(1), `age`(2). |
| **PROCESS** | End this stand's keywords; project it. | *(none)* |
| **STOP** | End the run. | *(none)* |

## Density limits

| keyword | purpose | named parameters |
|---|---|---|
| **SDIMAX** | Per-species maximum SDI (self-thinning limit) + the mortality SDI band. | `species`(1), `sdimax`(2), `pct_lo`(5), `pct_hi`(6). |
| **SDICALC** | Select the SDI computation method (Reineke vs Zeide summation). | `method`(1). |
| **BAMAX** | Pin the stand maximum basal area (derives per-species SDImax). | `bamax`(1). |

## Growth calibration & modifiers

| keyword | purpose | named parameters |
|---|---|---|
| **GROWTH** | Declare the type/period of the input measured growth used to calibrate. | `idg`(1, diameter data type), `ihtg`(2, height type), `dg_period`(3), `htg_period`(4). |
| **NOCALIB** | Disable the self-calibration (use uncalibrated coefficients). | `species`(1, 0=all). |
| **READCORD / REUSCORD** | Read / reuse a large-tree **diameter**-growth correction (COR). | `species`(1), `value`(2). |
| **READCORH / REUSCORH** | Read / reuse a large-tree **height**-growth correction. | `species`(1), `value`(2). |
| **READCORR / REUSCORR** | Read / reuse a small-tree **height**-growth correction. | `species`(1), `value`(2). |
| **BAIMULT** | Diameter (basal-area-increment) growth multiplier. | `date`(1), `species`(2), `multiplier`(3), `dbh_min`(4), `dbh_max`(5). |
| **HTGMULT** | Height-growth multiplier. | `date`(1), `species`(2), `multiplier`(3), `dbh_min`(4), `dbh_max`(5). |
| **CRNMULT** | Crown-ratio/-width multiplier. | `date`(1), `species`(2), `multiplier`(3). |
| **REGDMULT / REGHMULT** | Regeneration diameter / height multipliers. | `date`(1), `species`(2), `multiplier`(3). |
| **DGSTDEV** | Scale the stochastic diameter-growth standard deviation. | `value`(1). |
| **SERLCORR** | Serial-correlation (AR/MA) parameters for stochastic DG. | `ar`(1, BJPHI), `ma`(2, BJTHET). |
| **RANNSEED** | Set the random-number seed. | `seed`(1). |
| **FIXDG** | One-shot fixed diameter increment. | `date`(1), `species`(2), `value`(3), `dbh_min`(4), `dbh_max`(5). |
| **FIXHTG** | One-shot fixed height increment. | `date`(1), `species`(2), `value`(3), `dbh_min`(4), `dbh_max`(5). |
| **FIXMORT** | One-shot forced mortality rate. | `date`(1), `species`(2), `rate`(3), `dbh_min`(4), `dbh_max`(5), `option`(5/6). |
| **MORTMULT** | Mortality-rate multiplier (with DBH window). | `date`(1), `species`(2), `multiplier`(3), `dbh_min`(4), `dbh_max`(5). |
| **TREESZCP** / **SIZCAP** | Per-species size cap (max DBH/height). | `species`(1), `dbh_cap`(2). |

## Thinning & harvest

All `THIN*` keywords share the layout: `year`(1), then a residual/target(2),
`cut_efficiency`(3), `dbh_min`(4), `dbh_max`(5), `species`(6), `plot`(7).

| keyword | purpose | target parameter (2) |
|---|---|---|
| **THINBBA** | Thin **from below** to a residual basal area. | `residual_basal_area` |
| **THINABA** | Thin **from above** to a residual basal area. | `residual_basal_area` |
| **THINBTA** | Thin **from below** to a residual trees/acre. | `residual_tpa` |
| **THINATA** | Thin **from above** to a residual trees/acre. | `residual_tpa` |
| **THINSDI** | Thin to a residual Stand Density Index. | `residual_sdi` |
| **THINCC** | Thin to a residual crown-competition factor. | `residual_ccf` |
| **THINHT** | Thin by height. | `height` |
| **THINQFA** | Thin to a target TPA/BA/SDI by Q-factor. | per method |
| **THINRDEN** | Thin to a relative density (SDI-line). | `relative_density` |
| **THINDBH** | Remove a proportion within a DBH window. | `dbh_min`(2), `dbh_max`(3), `cut_efficiency`(4), `residual_tpa`(6), `species`(7). |
| **THINPT / SETPTHIN / THINPRSC / THINAUTO** | Point/prescription/automatic thinning variants. | per method |
| **SPECPREF** | Species preference order for removal. | `species`(1), `preference`(2). |
| **LEAVESP / SPLEAVE** | Protect a species from cutting (leave). | `species`(1). |
| **CUTEFF** | Default proportion of selected trees removed. | `proportion`(1). |
| **MINHARV** | Minimum harvest threshold (else no cut). | `min_volume`(1). |
| **SALVAGE** | Salvage dead/damaged trees. | `date`(1), `species`(2), … |
| **YARDLOSS** | Logging-damage / yarding loss fractions. | per field |

## Establishment & regeneration

| keyword | purpose | named parameters |
|---|---|---|
| **ESTAB** | Open an establishment packet (then `PLANT`/`NATURAL`/`SPROUT` … `END`). | `disturbance_date`(1). |
| **PLANT** | Schedule planting. | `year`(1), `species`(2), `tpa`(3), `survival_pct`(4), `age`(5), `height`(6), `shade`(7). |
| **NATURAL** | Schedule natural regeneration. | `year`(1), `species`(2), `tpa`(3), `survival_pct`(4), `age`(5), `height`(6), `shade`(7). |
| **SPROUT / NOSPROUT** | Enable / disable stump sprouting. | *(flags)* |

> Note: FVSjl's establishment fires only for **scheduled** ESTAB/PLANT/NATURAL
> activities (+ sprouting). FVS's default-on **auto-establishment** (natural ingrowth
> by stocking each cycle, switched off by `NOAUTOES`) is not yet wired — so to match
> a stock-FVS run that auto-regenerates, add `NOAUTOES`, or schedule the regen
> explicitly.

## Volume & merchandising

| keyword | purpose | named parameters |
|---|---|---|
| **VOLUME** | Cubic merch standards (stump/top dia/min DBH). | `date`(1), `species`(2), `cf_stump`, `cf_top`, `cf_min_dbh`, … |
| **BFVOLUME** | Board-foot merch standards + equation. | `date`(1), `species`(2), `bf_stump`, `bf_top`, `bf_min_dbh`. |
| **VOLEQNUM** | Override the per-species volume equation number. | `species`(1), `equation`(2). |
| **MCDEFECT / BFDEFECT** | Per-species cubic / board defect curves by DBH. | `date`(1), `species`(2), then the DBH-class % values. |
| **MCFDLN / BFFDLN** | Log-linear form-model coefficients (B0/B1). | `species`(1), `b0`(2), `b1`(3). |

## Species groups

| keyword | purpose | named parameters |
|---|---|---|
| **SPGROUP** | Define a named species group (referenced as `−N` elsewhere). | `group`(1), then member species. |

## Site & treatments

| keyword | purpose | named parameters |
|---|---|---|
| **SETSITE** | Scheduled mid-run site/habitat/SDImax change. | `date`(1), `habitat`(2), `bamax`(3), `species`(4), `site_index`(5), `si_flag`(6), `sdimax`(7). |
| **FERTILIZ** | Schedule a fertilizer application. | `date`(1), `nitrogen`(2). |

## Output & database

| keyword | purpose | named parameters |
|---|---|---|
| **ECHOSUM** | Write the per-cycle summary to the `.sum` file. | *(none)* |
| **DATABASE … END** | Open a database-output block (SQLite). | *(block; see below)* |
| **DSNOUT** | Output database filename (inside `DATABASE`). | `filename` — on the following line. |
| **SUMMARY** | Emit the `FVS_Summary` table. | `level`(1, optional). |
| **TREELIDB** | Emit the per-tree `FVS_TreeList` table. | `level`(1, optional). |
| **CUTLIST** | Emit the removed-tree (cut) list. | `level`(1). |
| **COMPUTE … END / COMPUTDB** | Event-monitor user variables (and DB output). | *(expression block)* |
| **STRCLASS** | Stand structural-stage classification output. | `date`(1). |

## Event monitor

| keyword | purpose | named parameters |
|---|---|---|
| **IF … THEN … {activities} … ENDIF** | Conditionally schedule activities. | the condition expression follows `IF`; the bracketed activities run when it is true. |
| **COMPUTE … END** | Define user variables / expressions evaluated each cycle. | *(expression block)* |

## Fire & fuels (FFE)

| keyword | purpose | named parameters |
|---|---|---|
| **FMIN … END** | Open the Fire-and-Fuels (FFE) keyword block. | *(block)* |
| **SIMFIRE** | Schedule a simulated fire. | `date`(1), … |
| **POTFIRE / POTFLAME / FLAMEADJ** | Potential-fire / flame outputs & adjustments. | per field |
| **CARBREPT / CARBCALC** | Stand carbon report / carbon-calc method. | `date`(1) / `method`(1). |
| **TOPKILL** | Top-kill (broken-top) damage event. | `date`(1), `species`(2), … |
| **HTGSTOP** | Stop height growth (top damage). | `date`(1), `species`(2), … |

## Economics (ECON)

| keyword | purpose | named parameters |
|---|---|---|
| **ECON … END** | Open the economic-analysis block. | *(block)* |
| **ANNUCST** | Annual cost. | `date`(1), `amount`(2), … |
| **HRVRVN** | Harvest revenue (per unit volume). | `species`(1), `value`(2), … |
| **HRVVRCST** | Variable harvest cost. | per field |
| **TCONDMLT** | Tree-condition multiplier. | per field |

## Compression & tripling

| keyword | purpose | named parameters |
|---|---|---|
| **COMPRESS** | Cluster tree records (PC-score) to speed projection. | `date`(1), … |
| **NOTRIPLE** | Disable record tripling. | *(none)* |
| **NUMTRIP** | Number of tripled records. | `count`(1). |

---

### Notes

- Parameter **field positions** come from the FVSjl keyword handlers
  (`src/engine/keyword_dispatch.jl`); the `# OPTION n — NAME (initre.f:…)` comment
  on each handler is the authoritative field-by-field source.
- The structured-YAML schema (`src/io/yaml_keywords.jl` `_KW_SCHEMA`) gives the
  machine-readable named parameters for the keywords most used in stand setup and
  thinning; extend it to add named params for more keywords.
- Keywords accepted but treated as `.sum`-inert no-ops in some variants (e.g.
  certain calibration cards) still parse without error.
