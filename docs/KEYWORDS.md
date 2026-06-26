# FVSjl keyword reference

Every keyword FVSjl recognizes, grouped by purpose, with its **named parameters**.
A keyword can be written in **two equivalent forms** — FVSjl reads either and runs
them identically:

1. **Legacy `.key`** — fixed-column Fortran cards. The keyword name is in cols 1–8;
   each parameter is a 10-column field (field `n` at cols `10·n+1 … 10·n+10`). The
   *field number* in parentheses below is that 1-based slot.
2. **Modern `.yaml`** — an **order-aware, hierarchical** document (see next section).
   Parameters are **named** and keep their natural type (numbers unquoted, codes quoted).

A blank date field defaults to **cycle 1**; a date ≥ 1000 is a calendar year, < 1000 a
cycle number. Convert between the forms with `bin/fvsjl-translate.jl` (or
`examples/key_to_structured_yaml.py`); the round-trip is lossless.

---

## The YAML form — order-aware hierarchy

> **Order matters in FVS.** The keyword stream is an *ordered sequence*: activities
> schedule in input order, later keywords override earlier ones, and several keywords
> depend on the state a *prior* keyword established. The YAML form therefore stays an
> ordered list end-to-end — the hierarchy only **groups** that ordered stream into
> named, readable sections; it never reorders it.

The default (hierarchical) shape is a `stand:` map whose value is an **ordered list of
section blocks**, each block an **ordered list of keyword entries**:

```yaml
stand:
  - setup:               # ── stand setup & inventory ──
      - STDIDENT: { id: "my stand" }
      - STDINFO:  { forest_code: 80106, habitat: "231Dd", stand_age: 60 }
      - NUMCYCLE: { cycles: 5 }
  - species_groups:      # ── define groups BEFORE the THIN that names them ──
      - SPGROUP: { group: "MAPHIC" }
      - raw: "SM HI"     # the member-species list (carried verbatim)
  - treatments:          # ── ORDER MATTERS: same-cycle activities run top-to-bottom ──
      - THINBBA: { year: 2005, residual_basal_area: 90 }
      - THINDBH: { year: 2015, dbh_min: 0, dbh_max: 12, cut_efficiency: 0.8 }
  - setup:
      - TREEDATA: {}
      - PROCESS:  {}
      - STOP:     {}
```

**Flattening the blocks top-to-bottom reproduces the exact keyword order** — grouping
is strictly order-preserving (a section may appear more than once, e.g. `setup` above,
precisely so the original interleaving is kept). Section names are *labels only* on
read; the list order is the source of truth.

Each entry is one of:

| entry form | when used | example |
|---|---|---|
| `KEYWORD: { named: …, params: … }` | keyword is in the schema and every present field maps to a name | `THINBBA: { year: 2005, residual_basal_area: 90 }` |
| `keyword: "NAME"` + `params: [ … ]` | no schema, or a field has no name (positional, still lossless) | `keyword: "VOLEQNUM"`, `params: ["0", "203.0"]` |
| `raw: "…"` | a free-form continuation line (an `IF` condition, a `COMPUTE` body line, SPGROUP members, the `TREEFMT` format) | `raw: "(FRAC(CYCLE/2.0) EQ 0.0)"` |

A legacy **flat** form is also accepted (and emitted with `--flat`): a single ordered
`keywords:` list of those same entries — no `stand:`/sections.

### Order-significant relationships you must preserve

These are the places where *relative order changes the result*. The hierarchical form
keeps each in a single contiguous, ordered block:

| relationship | rule | section |
|---|---|---|
| **SPGROUP → THIN** | a `SPGROUP` must come **before** any thinning that references the group (via a negative species code). | `species_groups` then `treatments` |
| **same-cycle activities** | two activities scheduled in the **same cycle/year** are applied in **input order** (e.g. a `THINBBA` then a `THINDBH` in 2005). | `treatments` |
| **COMPUTE def-before-use** | a `COMPUTE` variable must be **defined before** any later keyword or condition that uses it; variables defined in an earlier cycle persist. | `event_monitor` |
| **IF … THEN … ENDIF** | the condition and its guarded activities form one ordered block; the activities fire in order when the condition holds. | `event_monitor` |
| **ESTAB … END** | `PLANT` / `NATURAL` / `SPROUT` cards live **inside** the `ESTAB … END` packet, in order. | `regeneration` |
| **override keywords** | a later `SETSITE`, multiplier, or merch-standard card **overrides** an earlier one for the overlapping date/species window. | various |

Everything else is grouped purely for readability; reordering *across* unrelated
sections has no effect, but the converter never does so — it preserves the input order.

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

## Worked examples (`.key` ↔ `.yaml`, side by side)

Each example shows the legacy card(s) and the equivalent hierarchical-YAML entry.
All are drawn from `test/keyword_coverage/scenarios/`.

### Thinning — `THINABA` (from above to a residual basal area)

```text                                  # .key
THINABA       2005.0     100.0
```
```yaml                                  # .yaml  (in the `treatments` section)
- THINABA:
    year: 2005
    residual_basal_area: 100
```
Thin **from above** in 2005 down to a residual 100 sq ft/ac of basal area. (`THINBBA`
is the same card "from below"; `THINBTA`/`THINATA` target trees-per-acre, `THINSDI`
an SDI, `THINCC` a CCF.)

### `THINDBH` — remove within a DBH window, by efficiency

```text
THINDBH       2015.0       0.0      12.0      0.80
```
```yaml
- THINDBH:
    year: 2015
    dbh_min: 0
    dbh_max: 12
    cut_efficiency: 0.8
```
Remove 80 % (`cut_efficiency`) of trees with DBH in [0, 12) in 2015.

### `SPGROUP` + group thin — **order-significant**

```text
SPGROUP   MAPHIC
SM HI
THINDBH       2005.0       0.0      99.0      1.00        -1     100.0
```
```yaml
- species_groups:
    - SPGROUP:
        group: "MAPHIC"
    - raw: "SM HI"            # member species (SM + HI), carried verbatim
- treatments:
    - keyword: "THINDBH"      # species = -1 references the group above
      params: ["2005.0", "0.0", "99.0", "1.00", "-1", "100.0"]
```
The `SPGROUP` block **must precede** the `THINDBH` whose species code `-1` references
the group. (Here `THINDBH` keeps the positional form because its field-5 species code
has no schema name — still lossless.)

### `SPECPREF` / `LEAVESP` — removal preference & protected species

```text
LEAVESP           63
```
```yaml
- LEAVESP:
    species: 63
```
Protect species 63 (white oak) from cutting during thins. `SPECPREF` similarly sets a
removal-preference order. These come **before** the thin they modify.

### Establishment — `ESTAB … PLANT … END` packet (order-significant)

```text
ESTAB         2000.0
PLANT         2000.0      LP       300.0      90.0
END
```
```yaml
- regeneration:
    - ESTAB:
        disturbance_date: 2000
    - PLANT:
        year: 2000
        species: "LP"          # alpha species code → quoted
        tpa: 300
        survival_pct: 90
    - END: {}
```
Open an establishment packet for the 2000 disturbance, plant 300 loblolly-pine TPA at
90 % survival, then close with `END`. `PLANT`/`NATURAL`/`SPROUT` cards live **inside**
the packet, in order.

### `VOLEQNUM` — override the cubic volume equation

```text
VOLEQNUM           0     203.0
```
```yaml
- VOLEQNUM:
    species: 0               # 0 = all species
    equation: 203
```

### `VOLUME` — cubic merchandising standards

```text
VOLUME        2000.0       0.0       1.0       4.0       8.0
```
```yaml
- VOLUME:
    date: 2000
    species: 0
    cf_stump: 1
    cf_top: 4
    cf_min_dbh: 8
```

### `MCDEFECT` — cubic defect curve by DBH class

```text
MCDEFECT      2000.0      22.0       5.0      10.0      15.0      20.0
```
```yaml
- keyword: "MCDEFECT"
  params: ["2000.0", "22.0", "5.0", "10.0", "15.0", "20.0"]
```
`date`(1), `species`(2), then the per-DBH-class defect percentages.

### `COMPUTE … END` — Event-Monitor user variable (**def-before-use**)

```text
COMPUTE       0
MYBA = BBA
END
```
```yaml
- event_monitor:
    - keyword: "COMPUTE"
      params: ["0"]
    - raw: "MYBA = BBA"      # the assignment body, carried verbatim
    - END: {}
```
Define `MYBA` = the before-thin basal area each cycle. The definition **must appear
before** any later keyword/condition that reads `MYBA`.

### `IF … THEN … ENDIF` — conditional scheduling (**ordered block**)

```text
IF
(FRAC(CYCLE/2.0) EQ 0.0)
THEN
THINBBA          0.0     110.0
ENDIF
```
```yaml
- event_monitor:
    - IF: {}
    - raw: "(FRAC(CYCLE/2.0) EQ 0.0)"   # the condition expression
    - THEN: {}
    - THINBBA:
        year: 0                          # 0 = current cycle
        residual_basal_area: 110
    - ENDIF: {}
```
Every even cycle, thin to 110 sq ft/ac. The condition and its guarded activities are
one contiguous ordered block.

### `SERLCORR` — serial-correlation parameters

```text
SERLCORR
```
```yaml
- SERLCORR: {}
```

### `TIMEINT` — cycle length

```text
TIMEINT          0      10.
```
```yaml
- TIMEINT:
    cycle: 0           # 0 = all cycles
    length: 10
```
Set every cycle to 10 years.

### `COMPRESS` — cluster records to speed projection

```text
COMPRESS      2000.0
```
```yaml
- COMPRESS:
    date: 2000
```

### `THINAUTO` — automatic thinning

```text
THINAUTO      2000.0
```
```yaml
- THINAUTO:
    year: 2000
```

### Database output — `DATABASE … SUMMARY/CUTLIST … END`

```text
DATABASE
SUMMARY            2
CUTLIST            2
END
```
```yaml
- output:
    - DATABASE: {}
    - SUMMARY:
        level: 2
    - CUTLIST:
        level: 2
    - END: {}
```
Open the SQLite-output block and request the `FVS_Summary` and `FVS_CutList` tables,
then close with `END`.

---

### Notes

- Parameter **field positions** come from the FVSjl keyword handlers
  (`src/engine/keyword_dispatch.jl`); the `# OPTION n — NAME (initre.f:…)` comment
  on each handler is the authoritative field-by-field source.
- The named-parameter schema (`src/io/yaml_keywords.jl` `_KW_SCHEMA`) and the
  section grouping (`_KW_SECTION` in the same file) are the machine-readable source
  for the YAML form; extend them to add named params / a section for more keywords.
  A keyword without a schema entry still round-trips via the positional
  `keyword:`/`params:` form, and an unrecognized section falls into `other`.
- The hierarchical grouping is **order-preserving by construction**: the writer emits
  consecutive same-section records as one block and never reorders, so flattening the
  blocks reproduces the exact `.key` keyword sequence. The reader treats section names
  as labels only — list order is authoritative. Both directions are lossless w.r.t.
  the dispatch-relevant record (name, values, presence, field text); the original
  card's exact column padding is intentionally not preserved (it carries no meaning
  the handlers use).
- Keywords accepted but treated as `.sum`-inert no-ops in some variants (e.g.
  certain calibration cards) still parse without error.
