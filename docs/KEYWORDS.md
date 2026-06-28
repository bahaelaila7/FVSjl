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

> There is also a **second, *semantic* YAML flavor** (`format: fvs-stand/v1`) that
> describes a stand by intent (`invyr`, `numcycle`, `treatments`, `treelist`…) rather
> than mirroring the keyword stream, and is unraveled into the canonical keyword order
> for you. That flavor — plus the tree **`.csv`** format and how the keyword and tree
> files pair up — is documented in **[FORMATS.md](FORMATS.md)**. The form described
> below is the order-preserving image of a `.key`.

---

## Contents

- [The YAML form — order-aware hierarchy](#the-yaml-form--order-aware-hierarchy)
- [Order-significant relationships](#order-significant-relationships-you-must-preserve)
- **Keyword groups** (every keyword + its parameters):
  - [Stand identification & run setup](#stand-identification--run-setup) — STDIDENT · STDINFO · DESIGN · SITECODE · INVYEAR · NUMCYCLE · TIMEINT · CYCLEAT · TREEFMT · TREEDATA · NOTREES · TFIXAREA · MANAGED · RESETAGE · PROCESS · STOP
  - [Density limits](#density-limits) — SDIMAX · SDICALC · BAMAX
  - [Growth calibration & modifiers](#growth-calibration--modifiers) — GROWTH · NOCALIB · READCORD/H/R · REUSCORD/H/R · BAIMULT · HTGMULT · CRNMULT · REGDMULT · REGHMULT · DGSTDEV · SERLCORR · RANNSEED · FIXDG · FIXHTG · FIXMORT · MORTMULT · MORTMSB · TREESZCP · NOHTDREG · HTGSTOP · TOPKILL
  - [Thinning & harvest](#thinning--harvest) — THINBBA · THINABA · THINBTA · THINATA · THINSDI · THINCC · THINHT · THINQFA · THINRDEN · THINDBH · THINPT · SETPTHIN · THINAUTO · SPECPREF · LEAVESP · SPLEAVE · CUTEFF · MINHARV · SALVAGE · YARDLOSS
  - [Establishment & regeneration](#establishment--regeneration) — ESTAB · PLANT · NATURAL · SPROUT · NOSPROUT
  - [Volume & merchandising](#volume--merchandising) — VOLUME · BFVOLUME · VOLEQNUM · MCDEFECT · BFDEFECT · MCFDLN · BFFDLN
  - [Species groups](#species-groups) — SPGROUP
  - [Site & treatments](#site--treatments) — SETSITE · FERTILIZ
  - [Output & database](#output--database) — ECHOSUM · DATABASE · DSNOUT · SUMMARY · TREELIDB · CUTLIST · COMPUTDB · STRCLASS
  - [Event monitor](#event-monitor) — IF · COMPUTE
  - [Fire & fuels (FFE)](#fire--fuels-ffe) — FMIN · SIMFIRE · FLAMEADJ · MOISTURE · DROUGHT · the FUEL\* · the SNAG\* · the POTF\* · SALVAGE/SALVSP · PILEBURN · FUELMOVE · FUELTRET · FUELMODL · DEFULMOD · FIRECALC · CARBREPT · CARBCALC · CANCALC · SOILHEAT
  - [Economics (ECON)](#economics-econ) — ECON · STRTECON · ANNUCST · HRVVRCST · HRVRVN · TCONDMLT
  - [Compression & tripling](#compression--tripling) — COMPRESS · NOTRIPLE · NUMTRIP
- [Worked examples (`.key` ↔ `.yaml`)](#worked-examples-key--yaml-side-by-side)

> **Conventions.** In each table the number in `(n)` after a parameter is its 1-based
> **field position** on the `.key` card (field `n` at cols `10·n+1 … 10·n+10`); the same
> name is the YAML key. `date`(1) blank ⇒ cycle 1; ≥ 1000 = calendar year, < 1000 = cycle.
> `species` accepts an alpha/FIA/PLANTS code, `0`/`ALL`, or `−N` for [SPGROUP](#species-groups) *N*.
> A keyword absent from these tables but recognized by FVSjl is listed in its group's notes.
> The **structures** (the two YAML flavors, the tree `.tre`/`.csv`) are in [FORMATS.md](FORMATS.md).

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
| **MORTMSB** | Alternate "mature-stand breakup" mortality (a DBH-range kill, `msbmrt.f`). | `date`(1), `species`(2), `mort_rate`(3), `dbh_min`(4), `dbh_max`(5). |
| **TREESZCP** / **SIZCAP** | Per-species size cap (max DBH/height). | `species`(1), `dbh_cap`(2). |
| **NOHTDREG** | Control the HT-DBH (`LHTDRG`) calibration: field 1 = `0` suppress / `1` invoke. | `mode`(1). |
| **HTGSTOP** | Top-damage that **scales height growth** (`htgstp.f`). | `date`(1), `species`(2), `dbh_min`(4), `dbh_max`(5), `severity`. |
| **TOPKILL** | Top-damage that **top-kills** a fraction of height (broken top). | `date`(1), `species`(2), `dbh_min`(4), `dbh_max`(5), `severity`. |

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

The FFE keywords live **inside an `FMIN … END` block** (open FFE with `FMIN`, close with
`END`). They configure the fire-and-fuels extension: a scheduled fire, fuel loadings and
dynamics, the snag (standing-dead) sub-model, and the potential-fire report. (`date`(1)
follows the usual rule; `species` accepts a code / `0`=all / `−N`=group.)

**Fire event & behavior**

| keyword | purpose | named parameters |
|---|---|---|
| **FMIN … END** | Open / close the FFE keyword block. | *(block)* |
| **SIMFIRE** | Schedule a simulated fire + its weather. | `date`(1), `wind`(2, mph), `moisture_code`(3, 1–4), `temp`(4, °F), `mort_code`(5, 0/1), `pct_burn`(6, %), `season`(7, 1–4). |
| **FLAMEADJ** | Adjust flame length & crown burn. | `flame_mult`(2), `crown_burn`(4, %). |
| **FIRECALC** | Fire-behavior method. | `method`(1) — `0` old-FM-logic (the SN default & only ported path; 1/2 warn). |
| **MOISTURE** | Override the 7 fuel-moisture values for a date. | `date`(1), then `1hr`(2), `10hr`(3), `100hr`(4), `3in+`(5), `duff`(6), `live_woody`(7), `live_herb`(8, blank ⇒ live-woody) — all %. |
| **DROUGHT** | Drought years (UT/CR/LS only). | *(recognized; **no-op in SN**)* |
| **CANCALC** | Canopy base-ht / bulk-density for the crown-fire model. | *(recognized; **no-op in SN** — no crown fire)* |

**Fuel loadings & dynamics**

| keyword | purpose | named parameters |
|---|---|---|
| **FUELINIT** | Set initial **hard** fuel loadings (tons/ac by size class). | size-class loadings (`<.25`, `.25–1`, `1–3`, `3–6`, `6–12`, `>12`, litter, duff). |
| **FUELSOFT** | Set initial **soft** (decayed) fuel loadings. | same size classes as FUELINIT. |
| **FUELMODL** | Force standard fire-behavior fuel model(s). | `date`(1), then up to **3** `(model#, weight)` pairs (fields 2–7; weights normalized). |
| **DEFULMOD** | Define / alter a custom fuel model. | model number + the model's fuel-load / depth parameters (on following lines). |
| **FUELDCAY** | Set per-pool fuel **decay rates**. | per size-class decay rate. |
| **FUELMULT** | Multiply the fuel decay rates. | per size-class multiplier. |
| **FUELPOOL** | Assign a species to a fuel **decay class**. | `species`(1), `decay_class`(2, 1–4). |
| **FUELMOVE** | Move fuel between size pools. | `date`(1), `from`(2, 0–11), `to`(3), `amount`(4), `proportion`(5), `leave`(6), `target_final`(7). |
| **FUELTRET** | Fuel-treatment depth adjustment. | `date`(1), `treatment`(2, 0–2), `harvest`(3, 1–3), `depth_mult`(4, −1 ⇒ table). |
| **DUFFPROD** | Proportion of decayed wood that becomes **duff**. | `proportion`(per pool). |
| **PILEBURN** | Jackpot / pile burn. | `date`(1), `type`(2), `affect`(3, %), `treat`(4, %), `consumption`(5, %), `mortality`(6, %). |

**Snags (standing dead) & salvage**

| keyword | purpose | named parameters |
|---|---|---|
| **SNAGINIT** | Add user snags. | `species`(1), `dbh`(2), `ht_at_death`(3), `cur_ht`(4), `age`(5), `density`(6, stems/ac). |
| **SNAGFALL** | Per-species snag **fall** rates. | `species`(1), then the fall-rate parameters. |
| **SNAGDCAY** | Per-species snag **decay rate** `DECAYX` (hard→soft). | `species`(1), `decayx`(2). |
| **SNAGBRK** | Per-species snag **height-loss** (`YRS50`/`YRS30` → HTX). | `species`(1), `yrs_to_50pct`(2), `yrs_to_30pct`(3). |
| **SNAGPSFT** | Per-species fraction **soft at creation**. | `species`(1), `prop_soft`(2, 0–1). |
| **SNAGPBN** | Post-burn snag-fall parameters. | `pb_soft`(1), `pb_small`(2), `pb_time`(3, yrs), `pb_size`(4), `pb_scorch`(5). |
| **SALVAGE** | Remove (salvage) snags. | `date`(1), `dbh_min`(2), `dbh_max`(3), `max_age`(4), `ok_soft`(5, 0=all/1=hard/2=soft), `prop_removed`(6), `prop_left`(7). |
| **SALVSP** | Salvage species cut/leave list. | `date`(1), `species`(2), `flag`(3, <1 cut-list / ≥1 leave-list). |
| **FMORTMLT** | Fire-caused mortality multiplier. | `date`(1), `species`(2), `multiplier`(3), `dbh_min`(4), `dbh_max`(5). |

**Reports & potential fire**

| keyword | purpose | named parameters |
|---|---|---|
| **CARBREPT** | Request the FFE Stand **Carbon** report. | *(flag)* |
| **CARBCALC** | Carbon method & units. | `method`(1, 0=FFE/1=Jenkins), `units`(2, 0=US-t/ac / 1=t/ha / 2=t/ac). |
| **POTFIRE** / **POTFLAME** | Request the **Potential Fire** report. | *(flag)* |
| **POTFMOIS** | Potential-fire moisture by scenario. | `scenario`(1, 1=severe/2=moderate), then 7 moisture % (2–8). |
| **POTFWIND / POTFTEMP / POTFSEAS / POTFPAB** | Potential-fire wind / temp / season / % area burned. | `severe`(1), `moderate`(2). |
| **SOILHEAT** | Request the soil-heating report. | *(recognized; report not emitted)* |

> Report-only FFE keywords (`BURNREPT`, `FUELOUT`, `SNAGSUM`, `MORTREPT`, …) are
> **recognized no-ops** — the equivalent data is available via the [DBS](#output--database)
> tables. Any FFE keyword *not* listed here warns at parse time (it may change results).

## Economics (ECON)

The ECON keywords live **inside an `ECON … END` block**. They drive the discounted
cost/revenue analysis (`FVS_EconSummary`) and the log-graded harvest value
(`FVS_EconHarvestValue`).

| keyword | purpose | named parameters |
|---|---|---|
| **ECON … END** | Open / close the economic-analysis block. | *(block)* |
| **STRTECON** | Analysis start + discount rate. | `start`(1, year/delay), `discount_rate`(2, **%**), `known_sev`(3), `compute_sev`(4, flag). |
| **ANNUCST** | Annual management cost. | `amount`(1, $/ac/yr). |
| **HRVVRCST** | Variable harvest **cost** by DBH class. | `amount`(1), `dbh_class`(2), `dbh_min`(3), `dbh_max`(4). |
| **HRVRVN** | Harvest **revenue** by species + DBH. | `amount`(1, $/unit), `unit`(2, 1=TPA / 2=BF·1000 / 3=Ft³·100 / 4=BF·1000 per **log** / 5=Ft³·100 per **log**), `dbh_min`(3), `species`(4). |
| **TCONDMLT** | Tree-condition (point-weight) multiplier. | per field *(faithful single-point)*. |

> The log-graded units (`HRVRVN … 4`/`5`) feed the per-DIB-class `FVS_EconHarvestValue`
> detail table; the per-tree units (1/2/3) feed the discounted summary. See
> [docs/FORMATS.md](FORMATS.md) for how a harvest is valued.

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
