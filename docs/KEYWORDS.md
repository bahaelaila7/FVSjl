# FVSjl keyword reference

Every keyword FVSjl recognizes, grouped by purpose, with — for **each keyword** — a full
explanation of *what it does* (the underlying forestry/model concept, not just a restatement
of the name), its **named parameters** (field position, units, defaults, codes), and a
**worked example** in both input forms.

A keyword can be written in **two equivalent forms** — FVSjl reads either and runs them
identically:

1. **Legacy `.key`** — fixed-column Fortran cards. The keyword name is in cols 1–8; each
   parameter is a 10-column field (field `n` at cols `10·n+1 … 10·n+10`). The *field number*
   in parentheses below is that 1-based slot.
2. **Modern `.yaml`** — an **order-aware, hierarchical** document (see next section).
   Parameters are **named** and keep their natural type (numbers unquoted, codes quoted).

A blank date field defaults to **cycle 1**; a date ≥ 1000 is a calendar year, < 1000 a cycle
number. Convert between the forms with `bin/fvsjl-translate.jl` (or
`examples/key_to_structured_yaml.py`); the round-trip is lossless.

> **Variants.** This reference applies to **all four ported variants** — **SN** (Southern),
> **NE** (Northeast), **CS** (Central States), and **LS** (Lake States). The keyword grammar
> is *shared*: the variants differ in the underlying model coefficients (growth, volume,
> site, forest type) and species list, not in which keywords exist. Stock FVS has no variant
> keyword (you pick the `FVSsn`/`FVSne`/`FVScs`/`FVSls` binary); in FVSjl the variant is a run
> option (`--variant SN|NE|CS|LS`, the YAML `variant:` key, or `run_keyfile(...; variant=…)`) —
> see [FORMATS.md](FORMATS.md). A handful of keywords are recognized but `.sum`-inert in a
> given variant (noted inline below). Species codes and defaults in the examples are the
> Southern (SN) set unless stated; the *layout* is identical across variants.

> There is also a **second, *semantic* YAML flavor** (`format: fvs-stand/v1`) that describes a
> stand by intent (`invyr`, `numcycle`, `treatments`, `treelist`…) rather than mirroring the
> keyword stream, and is unraveled into the canonical keyword order for you. That flavor —
> plus the tree **`.csv`** format and how the keyword and tree files pair up — is documented
> in **[FORMATS.md](FORMATS.md)**. The form shown in each YAML example below is the
> order-preserving image of a `.key`.

---

## Contents

- [The YAML form — order-aware hierarchy](#the-yaml-form--order-aware-hierarchy)
- [Order-significant relationships](#order-significant-relationships-you-must-preserve)
- **Keyword groups** (every keyword documented with explanation + parameters + example):
  - [Stand identification & run setup](#stand-identification--run-setup) — STDIDENT · STDINFO · DESIGN · SITECODE · INVYEAR · NUMCYCLE · TIMEINT · CYCLEAT · TREEFMT · TREEDATA · NOTREES · TFIXAREA · MANAGED · RESETAGE · PROCESS · REWIND · STOP
  - [Density limits](#density-limits) — SDIMAX · SDICALC · BAMAX
  - [Growth calibration & modifiers](#growth-calibration--modifiers) — GROWTH · NOCALIB · READCORD/H/R · REUSCORD/H/R · BAIMULT · HTGMULT · CRNMULT · REGDMULT · REGHMULT · DGSTDEV · SERLCORR · RANNSEED · FIXDG · FIXHTG · FIXMORT · MORTMULT · MORTMSB · TREESZCP · NOHTDREG · HTGSTOP · TOPKILL
  - [Thinning & harvest](#thinning--harvest) — THINBBA · THINABA · THINBTA · THINATA · THINSDI · THINCC · THINHT · THINQFA · THINRDEN · THINDBH · THINPT · SETPTHIN · THINAUTO · THINPRSC · SPECPREF · LEAVESP · SPLEAVE · CUTEFF · MINHARV · SALVAGE · YARDLOSS
  - [Establishment & regeneration](#establishment--regeneration) — ESTAB · PLANT · NATURAL · SPROUT · NOSPROUT · NOAUTOES
  - [Species groups](#species-groups) — SPGROUP
  - [Site & treatments](#site--treatments) — SETSITE · FERTILIZ
  - [Volume & merchandising](#volume--merchandising) — VOLUME · BFVOLUME · VOLEQNUM · MCDEFECT · BFDEFECT · MCFDLN · BFFDLN
  - [Output & database](#output--database) — ECHOSUM · DATABASE · DSNOUT · SUMMARY · TREELIDB · CUTLIST · COMPUTDB · STRCLASS
  - [Event monitor](#event-monitor) — IF · COMPUTE
  - [Fire & fuels (FFE)](#fire--fuels-ffe) — FMIN · SIMFIRE · FLAMEADJ · FIRECALC · MOISTURE · DROUGHT · CANCALC · FUELINIT · FUELSOFT · FUELMODL · DEFULMOD · FUELDCAY · FUELMULT · FUELPOOL · FUELMOVE · FUELTRET · DUFFPROD · PILEBURN · SNAGINIT · SNAGFALL · SNAGDCAY · SNAGBRK · SNAGPSFT · SNAGPBN · SALVAGE · SALVSP · FMORTMLT · CARBREPT · CARBCALC · POTFIRE · POTFMOIS · POTFWIND/TEMP/SEAS/PAB · SOILHEAT
  - [Economics (ECON)](#economics-econ) — ECON · STRTECON · ANNUCST · HRVVRCST · HRVRVN · TCONDMLT
  - [Compression & tripling](#compression--tripling) — COMPRESS · NOTRIPLE · NUMTRIP

> **Conventions.** In each entry the number in `(n)` after a parameter is its 1-based **field
> position** on the `.key` card (field `n` at cols `10·n+1 … 10·n+10`); the same name is the
> YAML key. `date`(1) blank ⇒ cycle 1; ≥ 1000 = calendar year, < 1000 = cycle. `species`
> accepts an alpha/FIA/PLANTS code, `0`/`ALL`, or `−N` for [SPGROUP](#species-groups) *N*.
> The **structures** (the two YAML flavors, the tree `.tre`/`.csv`) are in [FORMATS.md](FORMATS.md).

---

## The YAML form — order-aware hierarchy

> **Order matters in FVS.** The keyword stream is an *ordered sequence*: activities schedule
> in input order, later keywords override earlier ones, and several keywords depend on the
> state a *prior* keyword established. The YAML form therefore stays an ordered list end-to-end
> — the hierarchy only **groups** that ordered stream into named, readable sections; it never
> reorders it.

The default (hierarchical) shape is a `stand:` map whose value is an **ordered list of section
blocks**, each block an **ordered list of keyword entries**:

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

**Flattening the blocks top-to-bottom reproduces the exact keyword order** — grouping is
strictly order-preserving (a section may appear more than once, e.g. `setup` above, precisely
so the original interleaving is kept). Section names are *labels only* on read; the list order
is the source of truth.

Each entry is one of:

| entry form | when used | example |
|---|---|---|
| `KEYWORD: { named: …, params: … }` | keyword is in the schema and every present field maps to a name | `THINBBA: { year: 2005, residual_basal_area: 90 }` |
| `keyword: "NAME"` + `params: [ … ]` | no schema, or a field has no name (positional, still lossless) | `keyword: "VOLEQNUM"`, `params: ["0", "203.0"]` |
| `raw: "…"` | a free-form continuation line (an `IF` condition, a `COMPUTE` body line, SPGROUP members, the `TREEFMT` format, a `READCOR*` value block) | `raw: "(FRAC(CYCLE/2.0) EQ 0.0)"` |

A legacy **flat** form is also accepted (and emitted with `--flat`): a single ordered
`keywords:` list of those same entries — no `stand:`/sections.

### Order-significant relationships you must preserve

These are the places where *relative order changes the result*. The hierarchical form keeps
each in a single contiguous, ordered block:

| relationship | rule | section |
|---|---|---|
| **SPGROUP → THIN** | a `SPGROUP` must come **before** any thinning that references the group (via a negative species code). | `species_groups` then `treatments` |
| **same-cycle activities** | two activities scheduled in the **same cycle/year** are applied in **input order** (e.g. a `THINBBA` then a `THINDBH` in 2005). | `treatments` |
| **COMPUTE def-before-use** | a `COMPUTE` variable must be **defined before** any later keyword or condition that uses it; variables defined in an earlier cycle persist. | `event_monitor` |
| **IF … THEN … ENDIF** | the condition and its guarded activities form one ordered block; the activities fire in order when the condition holds. | `event_monitor` |
| **ESTAB … END** | `PLANT` / `NATURAL` / `SPROUT` cards live **inside** the `ESTAB … END` packet, in order. | `regeneration` |
| **FMIN … END / ECON … END** | the FFE fire/fuels cards and the economics cards live **inside** their block, closed by `END`. | `fire_and_fuels` / `economics` |
| **override keywords** | a later `SETSITE`, multiplier, or merch-standard card **overrides** an earlier one for the overlapping date/species window. | various |

Everything else is grouped purely for readability; reordering *across* unrelated sections has
no effect, but the converter never does so — it preserves the input order.

---
## Stand identification & run setup

#### STDIDENT
**What it does:** Assigns the stand its identifying label and descriptive title. Unlike ordinary keywords, STDIDENT carries no fields on its own card — the identifier is read from the *next* line of the keyword file. The first whitespace-delimited token (up to 26 characters) becomes the stand ID (`NPLT`, the key that ties the run to its inventory and labels every `.sum`/database row); the remainder of the line (up to 72 characters) becomes the free-text title (`TITLE`) that is carried into the CSV/`.sum` output. It has no effect on the projection itself — it is pure bookkeeping.
**Parameters:**
- `id` — the label, taken from the line that follows the keyword. First token (≤26 chars) = stand ID; the rest of the line (≤72 chars) = descriptive title.
**Example:**
```text
STDIDENT
LOBLOLLY01                Loblolly pine demo, Talladega NF
```
```yaml
- STDIDENT: { id: "LOBLOLLY01" }
- raw: "LOBLOLLY01                Loblolly pine demo, Talladega NF"
```
*Names the stand `LOBLOLLY01` and titles it "Loblolly pine demo, Talladega NF".*

#### STDINFO
**What it does:** Supplies the stand's location, site, and topography attributes. Field 1 is the FVS location/forest code (`KODFOR`), which selects the region-8/9 volume equations, default latitude/longitude/elevation, and (in SN) certain growth overrides; a code that resolves to no recognized Southern forest is remapped to Talladega NF (80106), and Savannah River / Fort Bragg / reservation "pseudo-codes" are collapsed to their canonical 5-digit forest. Field 2 is the ecological-unit / habitat code (decoded by HABTYP into the PCOM ecological classification). The remaining fields set physical stand geometry — age, aspect, slope, elevation — that feed the topographic modifiers in the growth and site models.
**Parameters:**
- `forest_code`(1) — FVS location/national-forest code (`KODFOR`); drives volume equations and default geography. Blank/foreign SN codes default to 80106 (Talladega NF).
- `habitat`(2) — ecological-unit / habitat-type code (numeric → SNECU index, or an alpha unit like `231Dd`).
- `stand_age`(3) — stand age in years at inventory.
- `aspect`(4) — aspect in degrees (0–360); stored internally as radians.
- `slope`(5) — slope in percent; stored internally as a fraction (÷100).
- `elevation`(6) — elevation in hundreds of feet (e.g. 12 = 1200 ft); only applied when > 0.
- `stand_origin`(9) — 0 = natural, 1 = plantation; any out-of-range value is treated as 0.
**Example:**
```text
STDINFO      80106.      231.       60.      180.       35.       12.
```
```yaml
- STDINFO: { forest_code: 80106, habitat: "231", stand_age: 60, aspect: 180, slope: 35, elevation: 12 }
```
*Places a 60-year-old stand on Talladega NF (ecological unit 231), on a 35 % west-facing (180°) slope at 1200 ft.*

#### DESIGN
**What it does:** Describes the inventory sampling design so FVS can expand the tree records to per-acre values. The **basal-area factor** (BAF) is the prism/angle-gauge constant: every tree tallied "in" on a variable-radius point represents BAF square feet of basal area per acre, so a larger tree stands for fewer trees/acre. The **break DBH** (BRK) is the threshold below which trees are tallied on a fixed small-tree subplot instead of the variable-radius point; `fixed_plot_area_inverse` (FPA) is the inverse of that subplot's area (e.g. 300 ⇒ 1/300-acre plots). The remaining fields count the sample points, flag nonstockable points, weight the sample, and give the stockable (gross) proportion of the stand.
**Parameters:**
- `basal_area_factor`(1) — BAF for the variable-radius (prism) plot, ft²/ac per tallied tree. Default 40.
- `fixed_plot_area_inverse`(2) — inverse fixed-subplot area (FPA), e.g. 300 = 1/300 acre. Default 300.
- `break_dbh`(3) — DBH (in) below which trees are on the fixed subplot, above which on the variable plot (BRK). Default 5.0.
- `number_of_plots`(4) — number of sample points/plots (IPTINV), rounded to integer.
- `nonstockable_code`(5) — count of nonstockable points (NONSTK), rounded to integer.
- `sample_weight`(6) — per-stand sample weight (SAMWT); if ≤ 0 and points > 0, it defaults to the plot count.
- `stockable_proportion`(7) — gross stockable proportion of the stand (GROSPC, 0–1). A value entered as 1 < g ≤ 100 is read as a percent and divided by 100. Default 1.0.
**Example:**
```text
DESIGN        40.      300.        5.        1.
```
```yaml
- DESIGN: { basal_area_factor: 40, fixed_plot_area_inverse: 300, break_dbh: 5, number_of_plots: 1 }
```
*A single-point cruise: BAF-40 prism for trees ≥ 5" DBH, a 1/300-acre fixed subplot for smaller trees.*

#### SITECODE
**What it does:** Sets the **site index** — the height (in feet) that dominant/codominant trees of the site species reach at a reference base age — which is the model's core productivity measure: higher site index means faster diameter and height growth. Field 1 selects the species (0/blank = all species, −N = species group N, or a species sequence index), field 2 is the site-index value. The first SITECODE card whose field 1 is a real species (1…MAXSP) also fixes the *site species* (`ISISP`), the species whose site index represents the stand. Note FVS treats a field-2 value ≤ 5 as "not supplied" (SITEAR is only overwritten when the value exceeds 5).
**Parameters:**
- `site_species`(1) — species selector: 0/blank = all species, −N = SPGROUP N, else a species index; the first ≥1 value also sets the site species.
- `site_index`(2) — site index (ft at base age) applied to the selected species; in FVS a value ≤ 5 is ignored.
**Example:**
```text
SITECODE       0.       70.
```
```yaml
- SITECODE: { site_species: 0, site_index: 70 }
```
*Sets a site index of 70 ft for all species in the stand.*

#### INVYEAR
**What it does:** Records the calendar year of the inventory, i.e. the year the tree list was measured. This becomes the run's start year (`cycle_year[1]`) — the anchor from which every projection cycle boundary, scheduled activity date, and reported summary year is counted. All later "cycle number" dates (dates < 1000) are resolved relative to this year.
**Parameters:**
- `year`(1) — the inventory / start calendar year, rounded to an integer.
**Example:**
```text
INVYEAR     2000.
```
```yaml
- INVYEAR: { year: 2000 }
```
*Anchors the projection to a year-2000 inventory.*

#### NUMCYCLE
**What it does:** Sets how many projection cycles FVS runs — the stand is grown forward one cycle at a time, and this is the count of cycles (hence the length of the projection = cycles × cycle length). The value is accepted only if it is between 1 and MAXCYC (40); anything outside that range leaves the previous cycle count unchanged.
**Parameters:**
- `cycles`(1) — number of projection cycles, integer 1…40 (MAXCYC).
**Example:**
```text
NUMCYCLE       5.
```
```yaml
- NUMCYCLE: { cycles: 5 }
```
*Projects the stand for 5 cycles.*

#### TIMEINT
**What it does:** Sets the length (in years) of a projection cycle. Field 1 is a cycle index: 0 or blank means "apply this length uniformly to all cycles" (setting the global `FINT`/year), while a positive N sets only cycle N's length (stored as an individual override that the schedule builder cumulates into the boundary years). Cycle length matters because the growth models scale by it — diameter and height increments scale by FINT/5, and the autocorrelation / age bookkeeping key off it — so a 10-year cycle grows roughly twice the increment of a 5-year cycle. If field 2 is left blank the length defaults to 10 years. (SN is calibrated to 5-year cycles, NE to 10-year.)
**Parameters:**
- `cycle`(1) — cycle index to set (0/blank = all cycles uniformly, else the individual cycle N), absolute value taken.
- `length`(2) — cycle length in years. Default 10.
**Example:**
```text
TIMEINT        0.       10.
```
```yaml
- TIMEINT: { cycle: 0, length: 10 }
```
*Makes every cycle 10 years long.*

#### CYCLEAT
**What it does:** Requests an extra cycle boundary at a specific calendar year. FVS inserts the year as a new boundary strictly *inside* the run (it does not extend the end or move the start), which increases the effective cycle count and forces a reporting/growth-recomputation point at that year — useful to line the output up with a management event or a target reporting year. Duplicate and non-positive years are ignored.
**Parameters:**
- `year`(1) — the calendar year at which to insert a cycle boundary (must be positive and inside the run), rounded to integer.
**Example:**
```text
CYCLEAT     2013.
```
```yaml
- CYCLEAT: { year: 2013 }
```
*Adds a cycle break at 2013 so the projection reports (and re-grows) exactly at that year.*

#### TREEFMT
**What it does:** Supplies a custom FORTRAN FORMAT statement describing the column layout of the fixed-column `.tre` tree-record file, overriding the built-in default format. Like STDIDENT it carries no fields on its own card; the FORMAT string is read from the next (up to two) 80-column lines and concatenated. This lets non-standard tree files be read without reformatting. It is ignored when the companion tree file is a self-describing `.csv` (named columns).
**Parameters:**
- `format` — the FORTRAN FORMAT string, read from the following one or two ≤ 80-column lines (concatenated to 160 chars).
**Example:**
```text
TREEFMT
(I4,I7,F6.0,I1,A3,F4.1,F3.1,2F3.0,2(1X,2F3.1),I1,3(F3.0,I1),6F3.0,
 F2.0,2I1,I2,2I3,2I1,F3.0)
```
```yaml
- TREEFMT: { format: "(I4,I7,F6.0,I1,A3,F4.1,F3.1,2F3.0,2(1X,2F3.1),I1,3(F3.0,I1),6F3.0,F2.0,2I1,I2,2I3,2I1,F3.0)" }
```
*Reads the `.tre` file with an explicit column layout instead of the default tree-record format.*

#### TREEDATA
**What it does:** Triggers reading of the stand's tree records. FVSjl loads `<stem>.tre` (fixed-column, parsed with the current TREEFMT) or the self-describing `<stem>.csv` (named columns, no TREEFMT needed), where `<stem>` is the keyword file's path with the extension stripped. If no TREEDATA keyword appears at all, FVS still reads the tree file once at the end of keyword processing — so TREEDATA is mainly explicit control of *when* the inventory is read; NOTREES suppresses the default read.
**Parameters:** None.
**Example:**
```text
TREEDATA
```
```yaml
- TREEDATA: {}
```
*Reads the inventory tree list (`<stem>.tre` or `<stem>.csv`).*

#### NOTREES
**What it does:** Declares that the stand starts with no trees — a bare stand. It suppresses the automatic tree-file read that would otherwise happen at the end of keyword processing, so the projection begins from an empty tree list (typically to be populated by scheduled establishment: PLANT / NATURAL / ESTAB). Used to model afforestation or planting on open ground.
**Parameters:** None.
**Example:**
```text
NOTREES
```
```yaml
- NOTREES: {}
```
*Starts the stand bare (no inventory), to be populated by planting/natural regeneration.*

#### TFIXAREA
**What it does:** Sets the total fixed-plot area (`TFPA`) used when expanding the small-tree (DBH < break) sample to a per-acre basis. When present, the small-tree expansion factor becomes 1/TFPA rather than the default derived from the fixed-plot inverse; this corrects the trees-per-acre representation of the small-tree tally when the fixed subplots cover a known total area.
**Parameters:**
- `area`(1) — total fixed-plot area (TFPA), used as the small-tree expansion 1/TFPA.
**Example:**
```text
TFIXAREA     300.
```
```yaml
- TFIXAREA: { area: 300 }
```
*Expands the small-tree tally by 1/300 for the per-acre estimate.*

#### MANAGED
**What it does:** Flags the stand as managed/planted, which switches on the "planted/managed" diameter-growth term in the DGF large-tree growth model — a per-species increment `dg_planted[sp]` added to ln(DDS), so managed stands grow faster in diameter than the unmanaged default. A bare `MANAGED` card (or any field-2 value other than 0) sets the managed flag; an explicit field-2 value of 0 marks the stand unmanaged. A *dated* MANAGED (field 1 > 0), which schedules the change for a future cycle, is not yet ported (MANAGED is normally a one-time setup flag).
**Parameters:**
- `date`(1) — optional date; a dated (> 0) form is deferred (not applied). Blank/0 = immediate.
- `flag`(2) — 0 sets unmanaged; any other value (or a bare card) sets managed.
**Example:**
```text
MANAGED
```
```yaml
- MANAGED: {}
```
*Marks the stand as managed, enabling the planted-tree diameter-growth term.*

#### RESETAGE
**What it does:** Rebases the stand age so that at the activity's date the stand age equals field 2 — i.e. `age = age − date + start`. This is used when the inventoried age is wrong or when a regeneration event effectively resets the stand clock. Because the reset is applied after the cycle's summary row is written, the reset year's own `.sum` row keeps the old age and the new age takes effect on the following row. Field 1 may be a calendar year or a 1-based cycle number resolved against INVYEAR.
**Parameters:**
- `date`(1) — the date to rebase at: a calendar year (≥ 1000) or a cycle number (< 1000, resolved from INVYEAR + cycle length). Default cycle 1.
- `age`(2) — the stand age (years) the stand should have at that date. Default 0.
**Example:**
```text
RESETAGE    2010.       40.
```
```yaml
- RESETAGE: { date: 2010, age: 40 }
```
*Sets the stand age to 40 years as of 2010, rebasing the age track from then on.*

#### PROCESS
**What it does:** Terminates the current stand's keyword block and runs its projection. Every keyword before PROCESS configures the stand; PROCESS closes that configuration and hands it to the growth engine. If the stand had no explicit TREEDATA and is not NOTREES, the tree file is read at this point before projecting. Multiple stands in one keyword file are separated by PROCESS cards.
**Parameters:** None.
**Example:**
```text
PROCESS
```
```yaml
- PROCESS: {}
```
*Closes this stand's keywords and projects it.*

#### REWIND
**What it does:** In stock FVS, rewinds an input unit so the next stand re-reads it from the beginning — `REWIND 2` re-reads the tree-data unit, which is the mechanism for giving several stands (scenarios) the same inventory. In FVSjl it is a recognized no-op: the engine re-reads the shared tree data implicitly for each stand, so REWIND has no effect on results. The semantic YAML converter auto-emits it between stands for fidelity.
**Parameters:**
- `unit`(1) — the Fortran input unit to rewind (e.g. 2 = tree data). Ignored by FVSjl (no-op).
**Example:**
```text
REWIND         2.
```
```yaml
- REWIND: { unit: 2 }
```
*Would re-read the tree-data unit for the next stand (a no-op in FVSjl, which re-reads implicitly).*

#### STOP
**What it does:** Ends the entire run. STOP marks the end of the keyword stream — no further stands are read. It is the run-level terminator, as opposed to PROCESS, which only ends one stand and continues to the next.
**Parameters:** None.
**Example:**
```text
STOP
```
```yaml
- STOP: {}
```
*Ends the keyword run.*

## Density limits

#### SDIMAX
**What it does:** Overrides the maximum **Stand Density Index** (SDImax) and the self-thinning mortality band. SDI (Reineke) is `TPA·(QMD/10)^1.605` — the number of trees per acre a stand would have if its quadratic mean diameter were 10 inches, a scale-independent density measure; SDImax is the biological upper limit for the species/site at which density-dependent (self-thinning) mortality drives the stand. Field 1 selects the species (blank/0 = all, −N = group, else a species *sequence index* — a numeric field is the species index, not an FIA code), field 2 (> 0) pins that species' SDImax as a flagged user value the site setup then leaves untouched. Fields 5/6 set the lower/upper self-thinning percents (PMSDIL/PMSDIU): below the lower percent of SDImax there is no density mortality, and it ramps to full at the upper percent.
**Parameters:**
- `species`(1) — species selector: blank/0 = all, −N = SPGROUP N, else a species sequence index (SPDECD; the numeric field is the index, e.g. `11`, not an FIA code).
- `sdimax`(2) — the maximum SDI (SDIDEF) for the selected species; applied only when > 0.
- `pct_lo`(5) — lower self-thinning percent PMSDIL (floored at 10); stored as the fraction (÷100). FVS default 55.
- `pct_hi`(6) — upper self-thinning percent PMSDIU (capped at 95); stored as the fraction (÷100). FVS default 85.
**Example:**
```text
SDIMAX        11.      480.                            55.       85.
```
```yaml
- SDIMAX: { species: 11, sdimax: 480, pct_lo: 55, pct_hi: 85 }
```
*Caps species 11's SDImax at 480, with self-thinning mortality ramping between 55 % and 85 % of that.*

#### SDICALC
**What it does:** Selects how Stand Density Index is computed and the minimum DBH each method counts. The **Reineke** method uses the single stand quadratic-mean-diameter formula `TPA·(QMD/10)^1.605`; the **Zeide** (summation) method instead sums each tree's `(d/10)^1.605` over the stand, which handles irregular diameter distributions better. Field 3 chooses the method (≥ 1 ⇒ Zeide, else Reineke), and fields 1/2 give the minimum DBH included by each method. The chosen method drives *both* the reported `.sum` SDI column and the SDImax self-thinning mortality, matching FVS's shared `LZEIDE` flag.
**Parameters:**
- `dbh_stage`(1) — minimum DBH (in) for the Reineke/Curtis SDI (DBHSTAGE). Default 0.
- `dbh_zeide`(2) — minimum DBH (in) for the Zeide summation SDI (DBHZEIDE). Default 0.
- `method`(3) — SDI method: ≥ 1 ⇒ Zeide summation, else Reineke (LZEIDE).
**Example:**
```text
SDICALC        1.        3.        1.
```
```yaml
- SDICALC: { dbh_stage: 1, dbh_zeide: 3, method: 1 }
```
*Uses the Zeide summation SDI (counting trees ≥ 3" DBH) for both the report and self-thinning.*

#### BAMAX
**What it does:** Pins the stand's maximum **basal area** (BAMAX, ft²/acre — basal area is the cross-sectional area of tree stems at breast height summed per acre, a standard stocking measure). When set (> 0), the site setup derives every species' default SDImax from it as `BAMAX / (0.5454154·PMSDIU)` instead of from the per-species SDI constants, so the whole stand's density-driven self-thinning mortality is keyed to the user's maximum basal area. Without the keyword the SDImax stays dynamic (per-species defaults).
**Parameters:**
- `bamax`(1) — the stand maximum basal area (ft²/ac); applied only when > 0.
**Example:**
```text
BAMAX        300.
```
```yaml
- BAMAX: { bamax: 300 }
```
*Fixes the stand's maximum basal area at 300 ft²/ac, from which every species' SDImax is derived.*
## Growth calibration & modifiers

These keywords adjust FVS's growth, mortality, and calibration machinery. Two broad
families live here: **calibration controls** (how FVS tunes its diameter/height-growth
equations to your inventory, and the large-/small-tree correction terms `COR`) and
**modifiers** (multipliers, one-shot scalers, forced mortality, size caps, and top-damage
events that bend the projection).

> **Shared DBH-window convention.** Several cards carry a `dbh_min`(4)/`dbh_max`(5) pair
> that restricts the effect to trees inside it. FVS uses a **half-open** window `[dbh_min,
> dbh_max)` for the mortality/one-shot cards (`MORTMULT`, `FIXDG`, `FIXHTG`, `FIXMORT`)
> and a **closed** window `[dbh_min, dbh_max]` for `CRNMULT`. Defaults are `dbh_min = 0`,
> `dbh_max = 99999` (all trees). The plain growth multipliers `BAIMULT`/`HTGMULT`/
> `REGDMULT`/`REGHMULT` are **not** DBH-windowed. `HTGSTOP`/`TOPKILL` use a **height**
> window, not a DBH window.

> **Species field.** As everywhere, `species` accepts an alpha/FIA/PLANTS code, a numeric
> species **sequence index**, `0`/blank = all species, or `−N` for [SPGROUP](#species-groups) *N*.

> **Multiplier precedence.** For the persistent multipliers (`BAIMULT`, `HTGMULT`,
> `MORTMULT`, `CRNMULT`, `REGDMULT`, `REGHMULT`) the most recent card dated on or before
> the cycle wins, and a species-specific card beats an all-species (`species = 0`) card of
> the same date. The multiplier persists from its date onward; there is no useful default
> for the multiplier field itself — a blank reads `0.0` (which zeroes growth), so always
> supply it.

#### GROWTH
**What it does:** Declares the *type and measurement period* of the input growth data that
FVS uses when it self-calibrates its diameter- and height-growth equations to your stand
(the `LSTART`/calibration pass). It does not itself add growth — it tells the calibrator
how to read the DG/HTG fields on the tree records (an increment vs. a past diameter) and
over how many years they were measured, so the fitted correction (`COR`) is on the right
basis. FVSjl currently honors the default reading (the DG field is a 5-year increment,
which is bit-exact); the past-DBH interpretations (`idg` = 1 or 3) and non-5-year periods
are captured but deferred (the WK3 past-DBH calibration chunk).
**Parameters:**
- `idg`(1) — diameter-data type: `0` none/increment, `1` or `3` the DG field is a *past DBH* (→ PDBH), `2` increment. Default `0`.
- `dg_period`(2) — FINT, the diameter-growth measurement period in years. Default `5`.
- `ihtg`(3) — height-data type (analogous to `idg`). Default `0`.
- `htg_period`(4) — FINTH, the height-growth measurement period in years. Default `5`.
- `mort_period`(5) — FINTM, the mortality-observation period in years. Default `5`.
**Example:**
```text
GROWTH           2.0       5.0       0.0
```
```yaml
- GROWTH: { idg: 2, dg_period: 5, ihtg: 0 }
```
*Declares the diameter-growth field as a 5-year increment and supplies no measured height growth.*

#### NOCALIB
**What it does:** Turns **off** FVS's diameter-growth *self-calibration* (`LDGCAL`) for a
species. Normally FVS compares the growth implied by your inventory's DG measurements to
its model prediction and fits a per-species scale correction (the `COR` factor) so the
projection starts consistent with your data; `NOCALIB` skips that fit and uses the raw,
uncalibrated coefficients instead. Use it when your DG measurements are unreliable or you
want the pure regional model. (SN also clears the height self-calibration flag, but FVSjl
does no large-tree height self-calibration, so the height side is already inert.)
**Parameters:**
- `species`(1) — species to un-calibrate: `0`/blank = all, `−N` = SPGROUP *N*, else a code. Default all.
**Example:**
```text
NOCALIB          0.0
```
```yaml
- NOCALIB: { species: 0 }
```
*Disables diameter-growth self-calibration for every species (project on the raw coefficients).*

#### READCORD
**What it does:** Reads a table of per-species **large-tree diameter-growth corrections**
(`COR2`). Each value is added as `ln(COR2)` to the diameter-growth constant `DGCON`
*before* calibration (dgf.f), so it directly scales the large-tree DDS (diameter-squared
growth) of that species — a way to hand-tune large-tree diameter growth up or down when
you have local knowledge the regional model lacks. The correction terms are read from the
**continuation line(s)** that follow the keyword, `8F10.0` (eight 10-column fields per
line), one value per species in species-index order; a blank field reads `0.0` and (being
`≤ 0`) applies no correction.
**Parameters:**
- *(continuation block)* — up to `MAXSP` values, 8 per line, `value[i]` = the COR2 multiplier for species `i`. Blank ⇒ no correction for that species.
**Example:**
```text
READCORD
       1.0       1.0      1.15       1.0       1.0       1.0       1.0       1.0
```
```yaml
- READCORD: {}
- raw: "       1.0       1.0      1.15       1.0       1.0       1.0       1.0       1.0"
```
*Boosts large-tree diameter growth of species 3 by 15 % (COR2 = 1.15) and leaves the other species unchanged.*

#### REUSCORD
**What it does:** Re-enables the previously-read `READCORD` diameter corrections **without
reading a new block** — the multi-stand carry-over form. In a run that projects several
stands, one stand's `READCORD` table can be reused by a later stand with a bare `REUSCORD`
card, avoiding re-typing the values.
**Parameters:** None.
**Example:**
```text
REUSCORD
```
```yaml
- REUSCORD: {}
```
*Re-applies the last-read large-tree diameter-growth corrections to this stand.*

#### READCORH
**What it does:** Reads a table of per-species **large-tree height-growth corrections**
(`HCOR2`). Each value is added as `ln(HCOR2)` to the height-growth constant `HTCON`
before calibration (htgf.f), scaling that species' large-tree height increment. Same
continuation-block layout as `READCORD` (`8F10.0`, one value per species, blank ⇒ no
correction).
**Parameters:**
- *(continuation block)* — up to `MAXSP` values, 8 per line, `value[i]` = the HCOR2 multiplier for species `i`. Blank ⇒ no correction.
**Example:**
```text
READCORH
       1.0      0.90       1.0       1.0       1.0       1.0       1.0       1.0
```
```yaml
- READCORH: {}
- raw: "       1.0      0.90       1.0       1.0       1.0       1.0       1.0       1.0"
```
*Reduces large-tree height growth of species 2 by 10 % (HCOR2 = 0.90).*

#### REUSCORH
**What it does:** Re-enables the previously-read `READCORH` height corrections without
reading a new block (the multi-stand carry-over form for large-tree height corrections).
**Parameters:** None.
**Example:**
```text
REUSCORH
```
```yaml
- REUSCORH: {}
```
*Re-applies the last-read large-tree height-growth corrections to this stand.*

#### READCORR
**What it does:** Reads a table of per-species **small-tree height-growth corrections**
(`RCOR2`). The value is the small-tree height constant `RHCON` — a multiplier on the
seedling/sapling (REGENT) height-growth constant (regent.f), so it tunes the height growth
of the smallest trees, below the large-tree model's DBH threshold. Same continuation-block
layout as the other `READCOR*` cards.
**Parameters:**
- *(continuation block)* — up to `MAXSP` values, 8 per line, `value[i]` = the RCOR2 multiplier for species `i`. Blank ⇒ no correction.
**Example:**
```text
READCORR
       1.1       1.1       1.0       1.0       1.0       1.0       1.0       1.0
```
```yaml
- READCORR: {}
- raw: "       1.1       1.1       1.0       1.0       1.0       1.0       1.0       1.0"
```
*Raises small-tree height growth of species 1 and 2 by 10 % (RCOR2 = 1.1).*

#### REUSCORR
**What it does:** Re-enables the previously-read `READCORR` small-tree height corrections
without reading a new block (the multi-stand carry-over form).
**Parameters:** None.
**Example:**
```text
REUSCORR
```
```yaml
- REUSCORR: {}
```
*Re-applies the last-read small-tree height-growth corrections to this stand.*

#### BAIMULT
**What it does:** Multiplies **large-tree diameter growth** — specifically the diameter-
squared / basal-area increment (`DDS`) — of the selected species from the given date
onward. FVS applies it as `ln(multiplier)` added to `ln(DDS)` (dgdriv.f), so a multiplier
of `1.20` grows diameters 20 % faster in the DDS sense. Use it to nudge diameter growth for
a species or the whole stand (e.g. to reflect a fertilization or genetic-gain effect the
base model doesn't capture). Not DBH-windowed.
**Parameters:**
- `date`(1) — cycle/year the multiplier takes effect. Default cycle 1.
- `species`(2) — target species (`0` = all). Default all.
- `multiplier`(3) — the DDS multiplier (no default — supply it).
**Example:**
```text
BAIMULT       2010.0       11.0       1.2
```
```yaml
- BAIMULT: { date: 2010, species: 11, multiplier: 1.2 }
```
*From 2010 on, increases the basal-area increment (DDS) of species 11 by 20 %.*

#### HTGMULT
**What it does:** Multiplies **large-tree height growth** of the selected species from the
given date onward. The factor scales the computed height increment directly (htgf.f), so
`0.85` slows height growth by 15 %. Companion to `BAIMULT` for the height dimension. Not
DBH-windowed.
**Parameters:**
- `date`(1) — cycle/year the multiplier takes effect. Default cycle 1.
- `species`(2) — target species (`0` = all). Default all.
- `multiplier`(3) — the height-increment multiplier (no default).
**Example:**
```text
HTGMULT       2010.0        0.0      0.85
```
```yaml
- HTGMULT: { date: 2010, species: 0, multiplier: 0.85 }
```
*From 2010 on, reduces every species' height growth by 15 %.*

#### CRNMULT
**What it does:** Multiplies the **crown-ratio change** of the selected species (sn/crown.f),
adjusting how quickly crowns recede or expand as the stand develops. Because crown ratio
feeds diameter and height growth, this is an indirect growth modifier as well as a crown
report adjustment. Unlike `BAIMULT`/`HTGMULT`, `CRNMULT` **is** DBH-windowed, using a
**closed** `[dbh_min, dbh_max]` window.
**Parameters:**
- `date`(1) — cycle/year it takes effect. Default cycle 1.
- `species`(2) — target species (`0` = all). Default all.
- `multiplier`(3) — the crown-ratio-change multiplier (no default).
- `dbh_min`(4) — window lower bound (closed). Default `0`.
- `dbh_max`(5) — window upper bound (closed). Default `99999`.
**Example:**
```text
CRNMULT       2005.0        0.0       1.1       5.0      20.0
```
```yaml
- CRNMULT: { date: 2005, species: 0, multiplier: 1.1, dbh_min: 5, dbh_max: 20 }
```
*From 2005, increases the crown-ratio change by 10 % for all trees with DBH in [5, 20].*

#### REGDMULT
**What it does:** Multiplies **small-tree (regeneration) diameter growth** for the selected
species from the given date — the seedling/sapling analogue of `BAIMULT`, applied inside
the REGENT small-tree model rather than the large-tree DDS model. Useful for tuning early
diameter development of regenerating cohorts. Not DBH-windowed.
**Parameters:**
- `date`(1) — cycle/year it takes effect. Default cycle 1.
- `species`(2) — target species (`0` = all). Default all.
- `multiplier`(3) — the small-tree diameter-growth multiplier (no default).
**Example:**
```text
REGDMULT      2005.0        0.0       1.15
```
```yaml
- REGDMULT: { date: 2005, species: 0, multiplier: 1.15 }
```
*From 2005, increases small-tree diameter growth of all species by 15 %.*

#### REGHMULT
**What it does:** Multiplies **small-tree (regeneration) height growth** for the selected
species from the given date — the seedling/sapling analogue of `HTGMULT`, applied in the
REGENT small-tree height model. Not DBH-windowed.
**Parameters:**
- `date`(1) — cycle/year it takes effect. Default cycle 1.
- `species`(2) — target species (`0` = all). Default all.
- `multiplier`(3) — the small-tree height-growth multiplier (no default).
**Example:**
```text
REGHMULT      2005.0       11.0       0.9
```
```yaml
- REGHMULT: { date: 2005, species: 11, multiplier: 0.9 }
```
*From 2005, slows small-tree height growth of species 11 by 10 %.*

#### DGSTDEV
**What it does:** Sets `DGSD`, the number of standard deviations to which the **stochastic
serial-correlation diameter-growth variation** is bounded. FVS adds a random, serially-
correlated perturbation to each tree's diameter growth (the `DGSCOR` machinery); `DGSD`
caps how far that perturbation can stray from the mean. The default is `2.0`; setting
`DGSD < 1` turns the random variation **off** entirely, giving deterministic diameter
growth — the usual choice for reproducible, non-Monte-Carlo runs.
**Parameters:**
- `value`(1) — DGSD, the standard-deviation bound. Default `2.0`; `< 1` ⇒ deterministic DG.
**Example:**
```text
DGSTDEV          0.0
```
```yaml
- DGSTDEV: { value: 0 }
```
*Turns off the stochastic diameter-growth variation (deterministic growth).*

#### SERLCORR
**What it does:** Sets the two ARMA(1,1) **serial-correlation** parameters of the stochastic
diameter growth. "Serial correlation" means a tree that grew faster than predicted one
cycle tends to keep doing so the next — the perturbations are correlated through time
rather than independent. Field 1 is the autoregressive term `BJPHI` and field 2 the
moving-average term `BJTHET`; together they define the autocorrelation series (`autcor`),
so changing them alters the per-cycle variance/covariance of the DG perturbations.
**Parameters:**
- `ar`(1) — BJPHI, the AR(1) coefficient. Default `0.74`.
- `ma`(2) — BJTHET, the MA(1) coefficient. Default `0.42`.
**Example:**
```text
SERLCORR         0.6       0.3
```
```yaml
- SERLCORR: { ar: 0.6, ma: 0.3 }
```
*Sets a weaker serial correlation (AR 0.6, MA 0.3) in the stochastic diameter growth.*

#### RANNSEED
**What it does:** Reseeds the run's main random-number stream (`RANSED`) — the source for
all stochastic draws (DG perturbations, mortality, top-damage, fire). A non-zero seed makes
the whole random projection **reproducible**; the same seed always yields the same run.
Set at keyword-read time, before any draw. A field-1 value of `0` requests a clock-based
seed (`GETSED`) which is intentionally *not* reproduced here (it would be non-deterministic);
a blank field restarts the stream from its saved seed.
**Parameters:**
- `seed`(1) — RNG seed. Non-zero ⇒ install that seed (forced odd); `0` ⇒ clock seed (not reproduced, treated as no-op); blank ⇒ restart from the saved seed.
**Example:**
```text
RANNSEED     55329.0
```
```yaml
- RANNSEED: { seed: 55329 }
```
*Seeds the random-number stream with 55329 for a reproducible stochastic run.*

#### FIXDG
**What it does:** A **one-shot** diameter-growth scaler. In the single cycle whose year
range contains the keyword date, it multiplies the diameter increment (`DG`) of every
matching tree by `value`, then reverts — unlike `BAIMULT`, which persists. It fires after
all growth is computed but before mortality, so the mortality model sees the scaled growth,
and the tripled upper/lower records get the same factor. Use it to model a one-time growth
shock (e.g. a drought year, a release treatment) confined to one cycle. DBH-windowed
(`[dbh_min, dbh_max)`); multiple `FIXDG` cards in the same cycle compound in keyword order.
**Parameters:**
- `date`(1) — the cycle whose window contains this year gets scaled. Default cycle 1.
- `species`(2) — target species (`0` = all). Default all.
- `value`(3) — the one-cycle DG multiplier (no default).
- `dbh_min`(4) — window lower bound (half-open). Default `0`.
- `dbh_max`(5) — window upper bound (half-open). Default `99999`.
**Example:**
```text
FIXDG         2015.0        0.0       0.7       0.0      12.0
```
```yaml
- FIXDG: { date: 2015, species: 0, value: 0.7, dbh_min: 0, dbh_max: 12 }
```
*In the 2015 cycle only, cuts diameter growth to 70 % for all trees under 12" DBH (a one-year drought shock).*

#### FIXHTG
**What it does:** The one-shot height-growth analogue of `FIXDG`: in the cycle containing
the date, multiply the height increment (`HTG`) of every matching tree by `value`, then
revert. Same firing point (after growth, before mortality) and same tripled-record handling.
DBH-windowed (`[dbh_min, dbh_max)`); same-cycle cards compound in order.
**Parameters:**
- `date`(1) — the cycle whose window contains this year gets scaled. Default cycle 1.
- `species`(2) — target species (`0` = all). Default all.
- `value`(3) — the one-cycle HTG multiplier (no default).
- `dbh_min`(4) — window lower bound (half-open). Default `0`.
- `dbh_max`(5) — window upper bound (half-open). Default `99999`.
**Example:**
```text
FIXHTG        2015.0       11.0       0.5
```
```yaml
- FIXHTG: { date: 2015, species: 11, value: 0.5 }
```
*In the 2015 cycle only, halves the height growth of species 11.*

#### FIXMORT
**What it does:** A **one-shot forced-mortality override** applied in the cycle containing
the date, after the normal BA/self-thinning mortality (morts.f). It sets or adjusts the
periodic mortality rate for matching trees to a user value — a way to inject a specific
kill (a known insect/disease event, a sanitation cut modeled as mortality) instead of
relying on the background mortality model. How the given rate combines with the model's own
rate is controlled by the `combine` option. The rate is clamped to `[0,1]` unless
`combine = 3` (multiply). DBH-windowed (`[dbh_min, dbh_max)`).
**Parameters:**
- `date`(1) — cycle/year the override fires. Default cycle 1.
- `species`(2) — target species (`0` = all). Default all.
- `mort_rate`(3) — the mortality rate/period, meaning set by `combine`.
- `dbh_min`(4) — window lower bound (half-open). Default `0`.
- `dbh_max`(5) — window upper bound (half-open). Default `999`.
- `combine`(6) — how to apply the rate: `0` replace the model rate (default), `1` add, `2` take the max, `3` multiply. (When absent, rate is capped at 1.)
- `point_flag`(7) — point/size-concentration flag; the concentration path (`≥ 10`) is deferred. Default `0`.
**Example:**
```text
FIXMORT       2020.0        0.0      0.25       6.0      99.0       2.0
```
```yaml
- FIXMORT: { date: 2020, species: 0, mort_rate: 0.25, dbh_min: 6, dbh_max: 99, combine: 2 }
```
*In 2020, forces at least a 25 % mortality (`combine = 2` = max with the model rate) on all trees ≥ 6" DBH.*

#### MORTMULT
**What it does:** Multiplies the **background mortality rate** of matching trees from the
given date onward (morts.f) — a persistent scaler, unlike the one-shot `FIXMORT`. A
multiplier above 1 kills faster, below 1 slower; combined with the DBH window it can, e.g.,
elevate mortality only in a small-diameter class. Use it to represent a chronic stressor
(ongoing disease pressure) or to reduce mortality where the base model over-predicts.
DBH-windowed (`[dbh_min, dbh_max)`).
**Parameters:**
- `date`(1) — cycle/year it takes effect. Default cycle 1.
- `species`(2) — target species (`0` = all). Default all.
- `multiplier`(3) — the mortality-rate multiplier (no default).
- `dbh_min`(4) — window lower bound (half-open). Default `0`.
- `dbh_max`(5) — window upper bound (half-open). Default `99999`.
**Example:**
```text
MORTMULT      2005.0       11.0       2.0       0.0       5.0
```
```yaml
- MORTMULT: { date: 2005, species: 11, multiplier: 2.0, dbh_min: 0, dbh_max: 5 }
```
*From 2005 on, doubles the mortality rate of species 11 for trees under 5" DBH.*

#### MORTMSB
**What it does:** Enables the alternate **"mature-stand breakup"** mortality model — a
DBH-range kill routine (`msbmrt.f`) that replaces the normal density-driven mortality once
the stand's quadratic mean diameter reaches a threshold, modeling the accelerating breakup
of an old, large-tree stand. The six fields set the QMD trigger, the slope of the breakup
mortality relationship, an efficiency, the DBH window it acts on, and a method flag. Field
validation mirrors FVS exactly: any out-of-range field prints an error and **resets all six
parameters to their MSB-off defaults** (so a bad card disables the alternate model rather
than partially applying it).
**Parameters:**
- `qmd`(1) — QMDMSB, the quadratic-mean-diameter threshold that triggers breakup mortality (`> 0`). Default `999` (off).
- `slope`(2) — SLPMSB, the breakup-mortality slope, must be in `[−10, −1.605]`.
- `efficiency`(3) — EFFMSB, the mortality efficiency in `(0, 1]`. Default `0.90`.
- `dbh_min`(4) — DLOMSB, lower DBH bound (`≥ 0`). Default `0`.
- `dbh_max`(5) — DHIMSB, upper DBH bound (`≥ 0`, must exceed `dbh_min`). Default `999`.
- `method`(6) — MFLMSB, the method flag ∈ `{1, 2, 3}`. Default `1`.
**Example:**
```text
MORTMSB         20.0      -3.0      0.90       5.0      40.0       1.0
```
```yaml
- MORTMSB: { qmd: 20, slope: -3.0, efficiency: 0.90, dbh_min: 5, dbh_max: 40, method: 1 }
```
*Switches on mature-stand-breakup mortality once QMD reaches 20", applied to trees of 5–40" DBH with slope −3.0.*

#### TREESZCP  *(SIZCAP)*
**What it does:** Imposes a per-species **maximum tree size** (a diameter cap, and
optionally a height cap). Once a tree reaches the cap, FVS stops it from growing past the
limit — the cap governs diameter growth (dgbnd), height growth (htgf), and, unless the
no-mortality flag is set, applies an annual mortality rate to the capped trees (a size-cap
mortality floor in morts.f). Use it to keep a species from growing to unrealistic sizes on
a good site, or to model a species that "tops out". There is **no date field** — the cap
holds for the whole projection. (`SIZCAP` is the internal array name; the recognized
keyword is `TREESZCP`.)
**Parameters:**
- `species`(1) — target species: `0` = all, `−N` = SPGROUP *N*, else a code. Default all.
- `dbh_cap`(2) — maximum DBH (in). Default `999`.
- `mort_rate`(3) — annual mortality rate applied to trees at the cap. Default `1`.
- `no_mort_flag`(4) — IDMFLG; if set, suppress the cap mortality (grow-cap only). Default `0`.
- `ht_cap`(5) — maximum height (ft); `> 0` to activate. Default `999`.
**Example:**
```text
TREESZCP        11.0      30.0      0.05       0.0     120.0
```
```yaml
- TREESZCP: { species: 11, dbh_cap: 30, mort_rate: 0.05, no_mort_flag: 0, ht_cap: 120 }
```
*Caps species 11 at 30" DBH and 120 ft tall, applying a 5 % annual mortality to trees that reach the cap.*

#### NOHTDREG
**What it does:** Controls the per-species **height–diameter (HT-DBH) regression
calibration** (`LHTDRG`) — despite the "NO" in the name, it is a toggle, not just a
suppressor. Field 2 > 0 **invokes** the calibration (fit each species' HT-DBH relationship
to your inventory and use it for the small-tree DBH-from-height inversion); blank/zero
**suppresses** it, which is the SN default. FVSjl faithfully models the default/suppress
path (regent.f's inventory-inverse small-tree DBH branch); the invoke form — which would
switch non-Wykoff species to the Wykoff HT-DBH branch and run the ≥3-observation regression
fit — is unported and flagged rather than silently mis-run.
**Parameters:**
- `species`(1) — target species: `0`/blank = all, `−N` = SPGROUP *N*, else a code. Default all.
- `invoke`(2) — `> 0` invoke the HT-DBH calibration (`LHTDRG = TRUE`); blank/`0` suppress (the default). Default suppress.
**Example:**
```text
NOHTDREG         0.0       0.0
```
```yaml
- NOHTDREG: { species: 0, invoke: 0 }
```
*Suppresses the HT-DBH regression calibration for all species (the SN default; the modeled path).*

#### HTGSTOP
**What it does:** A scheduled **top-damage event** that *scales height growth* (htgstp.f,
activity 110) — a way to model a top-damaging disturbance (ice, wind, browse) that stunts
height development without necessarily killing the top. In the cycle containing the date,
for each matching-species tree whose **height** falls in the window `(HT1, HT2]`, the height
increment is multiplied by a kill proportion `PKIL` drawn from `BACHLO(AVEPRB, STDPBR)`
(deterministic and equal to `AVEPRB` when `STDPBR ≤ 0`); a per-tree probability `PRB < 1`
can let some trees escape via an RNG draw. Runs after tripling/mortality, before the height
increment is applied.
**Parameters:**
- `date`(1) — cycle/year the event fires. Default cycle 1.
- `species`(2) — target species: `0`/negative = all/group. Default all.
- `ht_min`(3) — HT1, lower height bound (ft, exclusive). Default `0`.
- `ht_max`(4) — HT2, upper height bound (ft, inclusive). Default `9999`.
- `prob`(5) — PRB, per-tree probability the tree is affected. Default `1` (all).
- `mean_severity`(6) — AVEPRB, mean of the kill proportion. Default `0`.
- `sd_severity`(7) — STDPBR, standard deviation of the kill proportion (`≤ 0` ⇒ deterministic). Default `0`.
**Example:**
```text
HTGSTOP       2010.0        0.0       0.0      30.0       1.0       0.5       0.0
```
```yaml
- HTGSTOP: { date: 2010, species: 0, ht_min: 0, ht_max: 30, prob: 1.0, mean_severity: 0.5, sd_severity: 0 }
```
*In 2010, halves (severity 0.5) the height growth of every tree up to 30 ft tall.*

#### TOPKILL
**What it does:** A scheduled **top-damage event** that *top-kills* part of a tree's height
(htgstp.f, activity 111) — a stronger disturbance than `HTGSTOP`: it removes a fraction of
the standing height rather than just slowing growth. In the cycle containing the date, for
each matching tree whose height is in `(HT1, HT2]`, the height is reduced to `H·(1 − PKIL)`
(with `PKIL` capped at 0.8). For a tall, large tree (`H ≥ 25`, `D ≥ 6`) whose implied broken
top diameter is `≥ 4"`, FVS records a **permanent broken top** (sets NORMHT/truncation) and
cuts the crown ratio accordingly. `PKIL` and the per-tree escape probability work exactly as
in `HTGSTOP` (`BACHLO(AVEPRB, STDPBR)`, `PRB`).
**Parameters:**
- `date`(1) — cycle/year the event fires. Default cycle 1.
- `species`(2) — target species: `0`/negative = all/group. Default all.
- `ht_min`(3) — HT1, lower height bound (ft, exclusive). Default `0`.
- `ht_max`(4) — HT2, upper height bound (ft, inclusive). Default `9999`.
- `prob`(5) — PRB, per-tree probability the tree is affected. Default `1`.
- `mean_severity`(6) — AVEPRB, mean fraction of height killed (capped at 0.8). Default `0`.
- `sd_severity`(7) — STDPBR, standard deviation of the kill fraction (`≤ 0` ⇒ deterministic). Default `0`.
**Example:**
```text
TOPKILL       2010.0        0.0      25.0    9999.0       1.0       0.3       0.0
```
```yaml
- TOPKILL: { date: 2010, species: 0, ht_min: 25, ht_max: 9999, prob: 1.0, mean_severity: 0.3, sd_severity: 0 }
```
*In 2010, tops-kills 30 % of the height of every tree taller than 25 ft, recording a permanent broken top on large trees.*
## Thinning & harvest

Thinning keywords schedule a **cut** on a given cycle. FVSjl runs them in two passes each cutting cycle (mirroring FVS `cuts.f`): first the **modifiers** (`SPECPREF`, `SPLEAVE`/`LEAVESP`, `MINHARV`, `SETPTHIN`, `TCONDMLT`) set the state the methods read, then the **methods** (`THIN*`) actually remove trees. Same-cycle activities apply in input order.

### What "thin from below" vs "thin from above" means

Every removal method works against a **residual target** — a density (basal area, trees/acre, SDI, relative density, or crown cover) you want the stand *left at* after the cut. The amount to remove is `current stocking − target`. The **direction** decides which trees are taken first to reach that target:

- **From below (low thinning)** — trees are ranked by size and the **smallest-diameter (or shortest) trees are removed first**, working *up* the size distribution until the residual target is met. Silviculturally this imitates natural suppression mortality: it clears out the overtopped, small, low-vigor stems and concentrates future growth on the largest dominant/codominant crop trees. It is the usual "improvement/tending" thin. Keywords: **THINBBA**, **THINBTA** (and `THINSDI`/`THINCC`/`THINRDEN`/`THINPT` with direction = 1).
- **From above (high/crown thinning)** — trees are ranked by size and the **largest-diameter trees are removed first**, working *down* until the residual is met. This harvests the most valuable large stems now (financial maturity, an overstory-removal/shelterwood cut, or a deliberate shift of size/species composition) and releases the smaller trees left behind. Keywords: **THINABA**, **THINATA** (and the index thins with direction = 2).
- **Throughout (proportional)** — direction = 0 on the index thins removes the same *fraction* from every eligible tree, preserving the diameter distribution.

The direction and the target *metric* (BA vs TPA) are baked into the keyword **name** for the four density thins (`B`/`A` = below/above; `BA`/`TA` = basal area / trees-per-acre target). The index thins carry direction in a field.

### The shared `THIN*` card layout

A THIN card is `KEYWORD` in cols 1–8, then up to seven 10-column fields. `year`(1) is always the date (blank ⇒ cycle 1; ≥ 1000 = calendar year, < 1000 = cycle). Fields 2–7 are method parameters, and **which slot means what depends on the keyword family** (FVSjl's YAML schema names the canonical slots `year, target, cut_efficiency, dbh_min, dbh_max, species, plot`, but the engine reads three distinct layouts):

| field | **THINBBA/ABA/BTA/ATA** (density thins) | **THINDBH / THINHT** (class-residual) | **THINSDI/CC/RDEN/PT** (index thins) |
|---|---|---|---|
| 1 | date | date | date |
| 2 | **residual target** (BA or TPA) | class lower bound (DBH / HT) | **residual target** (SDI / cover % / RD / metric) |
| 3 | **cut_efficiency** | class upper bound | **cut_efficiency** |
| 4 | dbh_min (`≥`) | cut_efficiency | species (SPDECD) |
| 5 | dbh_max (`<`) | species (SPDECD) | dbh_min (`≥`) |
| 6 | ht_min *(YAML calls this `species`)* | residual TPA in class | dbh_max (`<`) |
| 7 | ht_max *(YAML calls this `plot`)* | residual BA in class | direction (0 throughout / 1 below / 2 above) |

- **target**(2) — the density the stand (or class) is thinned *down to*. Its units are the metric named by the keyword. A target ≥ current stocking cancels the cut.
- **cut_efficiency**(3 or 4) — the fraction (0–1) of each **selected** tree's trees-per-acre actually removed (see CUTEFF). Blank ⇒ the global `CUTEFF` default (1.0). Where a residual target binds, the method computes its own efficiency (`remove/stock`) so the residual is hit exactly; efficiency only *caps* the cut when there is no target.
- **dbh_min**(4)/**dbh_max**(5) — restrict the cut to a diameter window, inches; `≥ dbh_min` and `< dbh_max`. Default 0 / no-upper-limit (a blank or `≤ dbh_min` max means "no upper bound"). Trees outside the window are untouched.
- **species** (SPDECD-decoded: `0`/blank = all, a positive integer = species **sequence index**, an alpha/FIA code, or `−N` for [SPGROUP](#species-groups) *N*) — only the class-residual and index thins read a species field; the four density thins ignore species selection (bias them with `SPECPREF`, protect species with `SPLEAVE`).
- **ht_min/ht_max** (density thins only) — an optional total-height window, feet, applied like the DBH window. *Caveat: the current YAML schema labels these slots `species`/`plot`; on a density thin they are the height window, so avoid writing `species:`/`plot:` there.*

> **Note on `SALVAGE`.** The base-model `SALVAGE` thinning keyword is *abandoned* in FVS (`cuts.f`) and is a recognized no-op in FVSjl. The active salvage operation is the **Fire & Fuels (FFE)** `SALVAGE` keyword (removes standing dead / snags) — documented under [Fire & fuels](#fire--fuels-ffe).

#### THINBBA
**What it does:** Thins **from below to a residual basal area**. Trees in the (optional) DBH/height window are ranked smallest-diameter-first and whole records are removed (each at `cut_efficiency`) until the stand's live basal area drops to the field-2 target; the last record is a partial cut that lands exactly on the target. This is the classic low-thinning-to-a-BA improvement cut — it strips the suppressed understory and leaves the largest crop trees.
**Parameters:**
- `year`(1) — date; blank ⇒ cycle 1, ≥ 1000 = calendar year, < 1000 = cycle number.
- `residual_basal_area`(2) — target residual BA, ft²/acre. A value ≥ current BA cancels the cut.
- `cut_efficiency`(3) — fraction 0–1 of each selected tree removed; blank ⇒ global CUTEFF (1.0).
- `dbh_min`(4) — smallest DBH eligible to cut, inches (`≥`). Default 0.
- `dbh_max`(5) — largest DBH eligible, inches (`<`); blank/`≤ dbh_min` ⇒ no upper limit.
- `ht_min`(6) — smallest total height eligible, ft (`≥`). Default 0. *(YAML slot named `species`.)*
- `ht_max`(7) — largest total height eligible, ft (`<`); blank ⇒ no upper limit. *(YAML slot named `plot`.)*
**Example:**
```text
THINBBA           2010        70.
```
```yaml
- THINBBA: { year: 2010, residual_basal_area: 70 }
```
*Thins the stand from below in 2010 down to 70 ft²/acre residual basal area, all species, all sizes.*

#### THINABA
**What it does:** Thins **from above to a residual basal area**. Identical to THINBBA but trees are ranked largest-diameter-first, so the biggest stems are harvested until BA falls to the target. Used for overstory removal / harvesting the most valuable large trees while releasing the smaller understory.
**Parameters:**
- `year`(1) — date.
- `residual_basal_area`(2) — target residual BA, ft²/acre.
- `cut_efficiency`(3) — fraction 0–1 removed per selected tree; blank ⇒ CUTEFF (1.0).
- `dbh_min`(4) / `dbh_max`(5) — DBH window, inches (`≥`/`<`); defaults 0 / no limit.
- `ht_min`(6) / `ht_max`(7) — optional height window, ft *(YAML slots `species`/`plot`)*.
**Example:**
```text
THINABA           2020        40.        1.0       12.
```
```yaml
- THINABA: { year: 2020, residual_basal_area: 40, cut_efficiency: 1.0, dbh_min: 12 }
```
*In 2020 removes the largest trees ≥ 12 in DBH down to a 40 ft²/acre residual — a from-above harvest of the big stems.*

#### THINBTA
**What it does:** Thins **from below to a residual trees-per-acre**. Same as THINBBA but the target metric is stem count: smallest trees are removed first until live TPA reaches the field-2 target. A density-control thin expressed in stems rather than basal area.
**Parameters:**
- `year`(1) — date.
- `residual_tpa`(2) — target residual trees/acre.
- `cut_efficiency`(3) — fraction 0–1 removed per selected tree; blank ⇒ CUTEFF (1.0).
- `dbh_min`(4) / `dbh_max`(5) — DBH window, inches.
- `ht_min`(6) / `ht_max`(7) — optional height window, ft *(YAML slots `species`/`plot`)*.
**Example:**
```text
THINBTA           2015       300.
```
```yaml
- THINBTA: { year: 2015, residual_tpa: 300 }
```
*Reduces the stand from below in 2015 to 300 trees/acre residual.*

#### THINATA
**What it does:** Thins **from above to a residual trees-per-acre** — largest stems removed first until TPA reaches the target. Used when you want to cut the biggest trees but express the residual as stem count.
**Parameters:**
- `year`(1) — date.
- `residual_tpa`(2) — target residual trees/acre.
- `cut_efficiency`(3) — fraction 0–1 removed per selected tree; blank ⇒ CUTEFF (1.0).
- `dbh_min`(4) / `dbh_max`(5) — DBH window, inches.
- `ht_min`(6) / `ht_max`(7) — optional height window, ft *(YAML slots `species`/`plot`)*.
**Example:**
```text
THINATA           2025       150.
```
```yaml
- THINATA: { year: 2025, residual_tpa: 150 }
```
*In 2025 harvests the largest trees down to 150 trees/acre.*

#### THINSDI
**What it does:** Thins to a residual **Stand Density Index** (Zeide summation SDI in SN: `Σ tpa·(DBH/10)^1.605`). Field 7 sets the direction — 0 removes the same fraction from every eligible tree (proportional), 1 thins from below, 2 from above — until the class SDI reaches the field-2 target. SDI is a size-independent density measure, so this targets stocking directly.
**Parameters:**
- `year`(1) — date.
- `residual_sdi`(2) — target residual SDI.
- `cut_efficiency`(3) — fraction 0–1; blank ⇒ CUTEFF (1.0). Bounds the from-below/above cut; the proportional path derives its own efficiency to hit the target.
- `species`(4) — SPDECD species/group filter (0 = all).
- `dbh_min`(5) / `dbh_max`(6) — DBH window, inches.
- `direction`(7) — 0 = throughout (proportional), 1 = from below, 2 = from above. Default 0.
**Example:**
```text
THINSDI           2030       200.         0.        0.        0.        1.
```
```yaml
- THINSDI: { year: 2030, residual_sdi: 200, species: 0, direction: 1 }
```
*Thins from below in 2030 to a residual SDI of 200, all species.*

#### THINCC
**What it does:** Thins to a residual **crown/canopy cover percent**. The class metric is crown-cover area (`Σ tpa·CW²`, forest-grown crown widths); the target cover % is converted to an equivalent crown area and trees are removed (proportional/below/above per field 7) until cover drops to the target. A target ≥ 100 % cancels the cut. Used to manage canopy closure / light directly.
**Parameters:**
- `year`(1) — date.
- `residual_ccf`(2) — target residual canopy cover, percent (0–100). *(Schema name is `residual_ccf`; the value is a cover percent, not a CCF.)*
- `cut_efficiency`(3) — fraction 0–1; blank ⇒ CUTEFF (1.0).
- `species`(4) — SPDECD filter (0 = all).
- `dbh_min`(5) / `dbh_max`(6) — DBH window, inches.
- `direction`(7) — 0 throughout / 1 below / 2 above. Default 0.
**Example:**
```text
THINCC            2035        60.         0.        0.        0.        1.
```
```yaml
- THINCC: { year: 2035, residual_ccf: 60, direction: 1 }
```
*Opens the canopy from below in 2035 to 60 % residual crown cover.*

#### THINHT
**What it does:** Removes trees within a **total-height class** down to a residual density. Like THINDBH but the class is bounded by height (fields 2–3) instead of diameter; within that height band it cuts `cut_efficiency` of each eligible record in record order until the residual TPA (field 6) and/or BA (field 7) target is met. Useful for taking out a specific height stratum.
**Parameters:**
- `year`(1) — date.
- `ht_min`(2) — lower height bound, ft (`≥`).
- `ht_max`(3) — upper height bound, ft (`<`); blank/`≤ ht_min` ⇒ no upper limit.
- `cut_efficiency`(4) — fraction 0–1; blank ⇒ CUTEFF (1.0).
- `species`(5) — SPDECD filter (0 = all).
- `residual_tpa`(6) — residual trees/acre left in the class.
- `residual_ba`(7) — residual basal area left in the class (ft²/acre); if > 0 the target is BA, else TPA.
**Example:**
```text
THINHT            2012         0.       25.         .8        0.       50.
```
```yaml
- THINHT: { year: 2012, ht_min: 0, ht_max: 25, cut_efficiency: 0.8, residual_tpa: 50 }
```
*In 2012 cuts 80 % of the stems under 25 ft tall, leaving 50 trees/acre in that height class.*

#### THINQFA
**What it does:** Thins a diameter distribution to a target **Q-factor** (a De Liocourt reverse-J: each smaller 1-inch-class holds `Q` times the trees of the next-larger class). It back-computes a per-DBH-class residual TPA/BA/SDI from the overall target and the Q ratio, then thins each class to it — shaping an uneven-aged/selection structure rather than a single flat residual. This is a **two-record** keyword: the second line is the target-units switch.
**Parameters (record 1):**
- `year`(1) — date.
- `dbh_min`(2) — smallest DBH class, inches. Default 0.
- `dbh_max`(3) — largest DBH class, inches. Default 24.
- `species`(4) — SPDECD filter (0 = all).
- `q_factor`(5) — the between-class ratio Q. Default 1.4.
- `class_width`(6) — diameter class width, inches. Default 2.
- `target`(7) — overall residual target (interpreted per record 2).
**Parameters (record 2):**
- units switch — one integer: `≤ 0` ⇒ target is basal area, `= 1` ⇒ trees/acre, `> 1` ⇒ SDI.
**Example:**
```text
THINQFA           2040        2.       20.         0.       1.3        2.       60.
1
```
```yaml
- keyword: "THINQFA"
  params: ["2040", "2", "20", "0", "1.3", "2", "60"]
- raw: "1"
```
*In 2040 shapes a reverse-J distribution across the 2–20 in classes (Q = 1.3, 2-in classes) to a 60 trees/acre target (the trailing `1` = TPA units).*

#### THINRDEN
**What it does:** Thins to a residual **Curtis relative density** (RD, a density index linearised about QMD). Direction (field 7) selects proportional/below/above; trees are removed until the class RD reaches the field-2 target. An alternative density metric to SDI, common in eastern management.
**Parameters:**
- `year`(1) — date.
- `residual_relsdi`(2) — target residual relative density.
- `cut_efficiency`(3) — fraction 0–1; blank ⇒ CUTEFF (1.0).
- `species`(4) — SPDECD filter (0 = all).
- `dbh_min`(5) / `dbh_max`(6) — DBH window, inches.
- `direction`(7) — 0 throughout / 1 below / 2 above. Default 0.
**Example:**
```text
THINRDEN          2018        50.         0.        0.        0.        1.
```
```yaml
- THINRDEN: { year: 2018, residual_relsdi: 50, direction: 1 }
```
*Thins from below in 2018 to a residual Curtis relative density of 50.*

#### THINDBH
**What it does:** Removes a proportion of trees within a **DBH class** down to a residual density. Within the diameter window (fields 2–3) and species filter it cuts `cut_efficiency` of each eligible record in record order until the residual TPA (field 6) and/or BA (field 7) is reached; with no residual target it simply removes `cut_efficiency` of every tree in the window. The workhorse "cut this size range of this species by X%" operation.
**Parameters:**
- `year`(1) — date.
- `dbh_min`(2) — lower DBH bound, inches (`≥`).
- `dbh_max`(3) — upper DBH bound, inches (`<`); blank/`≤ dbh_min` ⇒ no upper limit.
- `cut_efficiency`(4) — fraction 0–1 of each selected tree removed; blank ⇒ CUTEFF (1.0).
- `species`(5) — SPDECD filter (0 = all, `−N` = SPGROUP N).
- `residual_tpa`(6) — residual trees/acre left in the class (0 = cut to efficiency, no TPA floor).
- `residual_ba`(7) — residual basal area left in the class, ft²/acre (if > 0 the target is BA).
**Example:**
```text
THINDBH           2015         0.       12.         .8
```
```yaml
- THINDBH: { year: 2015, dbh_min: 0, dbh_max: 12, cut_efficiency: 0.8 }
```
*In 2015 removes 80 % of the trees under 12 in DBH (all species, no residual floor).*

#### THINPT
**What it does:** A **point-level** thin: each sample point is thinned independently to a residual, in the metric and on the point(s) named by a preceding `SETPTHIN` card this cycle (TPA, BA, SDI, crown cover, or Curtis RD). THINPT itself supplies the residual value, DBH window, and direction; the per-point stocking is scaled by the point count. Used to model treatments applied point-by-point rather than to the whole-stand average. (Requires a same-cycle `SETPTHIN`; without one it is a no-op.)
**Parameters:**
- `year`(1) — date.
- `residual`(2) — target residual, in the SETPTHIN metric.
- `cut_efficiency`(3) — fraction 0–1; blank ⇒ CUTEFF (1.0).
- `species`(4) — SPDECD filter (0 = all).
- `dbh_min`(5) / `dbh_max`(6) — DBH window, inches.
- `direction`(7) — 0 throughout / 1 below / 2 above.
**Example:**
```text
SETPTHIN          2020         0.        2.
THINPT            2020        40.         0.        0.        0.        0.        1.
```
```yaml
- SETPTHIN: { year: 2020, point: 0, metric: 2 }
- THINPT:   { year: 2020, residual: 40, direction: 1 }
```
*SETPTHIN selects all points with a basal-area metric; THINPT then thins each point from below to 40 ft²/acre residual BA in 2020.*

#### SETPTHIN
**What it does:** Sets the point-thinning **prescription** that a same-cycle `THINPT` reads: which sample point to treat (0 = every point) and which density **metric** the THINPT residual is expressed in. It removes no trees itself — it only arms THINPT.
**Parameters:**
- `year`(1) — date.
- `point`(2) — sample-point number to thin; `0` = all points. Default 0.
- `metric`(3) — residual metric code: `1` = trees/acre, `2` = basal area, `3` = SDI (Zeide), `4` = crown cover, `5` = Curtis relative density.
**Example:**
```text
SETPTHIN          2020         3.        1.
```
```yaml
- SETPTHIN: { year: 2020, point: 3, metric: 1 }
```
*Arms a trees-per-acre-metric point thin on point 3 in 2020 (paired with a THINPT that supplies the residual TPA).*

#### THINAUTO
**What it does:** **Automatic density-triggered** thinning. When the stand's stocking reaches `autmax` % of full/normal stocking (Reineke SDImax basis in SN/NE) the stand is thinned **from below** down to `autmin` % of full stocking; below the trigger it does nothing. Lets a single card impose a repeating "thin whenever it gets too dense" regime keyed to natural stocking rather than a fixed year.
**Parameters:**
- `year`(1) — date the automatic rule becomes active.
- `autmin`(2) — lower stocking bound, percent of full stocking (the residual). Default 45.
- `autmax`(3) — upper stocking bound, percent — the trigger. Default 60.
- `cut_efficiency`(4) — fraction 0–1 removed per selected tree; blank ⇒ CUTEFF (1.0).
**Example:**
```text
THINAUTO          2000        50.       70.
```
```yaml
- THINAUTO: { year: 2000, cut_efficiency: 50 }
```
*From 2000 on, whenever stocking reaches 70 % of full it is thinned from below back to 50 %.* (In the YAML schema fields 2–3 map to `cut_efficiency`; write the percents positionally if you need both `autmin` and `autmax`.)

#### THINPRSC
**What it does:** A **prescription** thin: it removes exactly the tree records the user **pre-marked** with a cut code in the input tree data (cut code ≥ 2), at `cut_efficiency`, across all species and diameters. Rather than a density rule, this replays a marked-tree cruise/marking — the model cuts the stems the field crew flagged.
**Parameters:**
- `year`(1) — date.
- `cut_efficiency`(2) — fraction 0–1 of each marked record removed; blank ⇒ CUTEFF (1.0). (A record left at ≤ 0.0005 residual is removed whole.)
**Example:**
```text
THINPRSC          2010         1.0
```
```yaml
- keyword: "THINPRSC"
  params: ["2010", "1.0"]
```
*In 2010 removes every tree record marked with a cut code ≥ 2 in the tree list, at full efficiency.*

#### SPECPREF
**What it does:** Sets a per-species **removal preference** that biases the ordering of any thin: the preference value is added to each tree's cut-priority weight, so a species with a higher value is removed *before* others of the same size. It does not by itself remove trees — it steers *which* species the density thins take first (e.g. "when thinning, take the pines before the oaks"). Persists until changed.
**Parameters:**
- `year`(1) — date the preference takes effect.
- `species`(2) — SPDECD species/group: a code, a positive sequence index, `0` = all, `−N` = SPGROUP N.
- `preference`(3) — integer priority added to the removal weight; higher ⇒ removed earlier.
**Example:**
```text
SPECPREF          2010        LP         5.
```
```yaml
- SPECPREF: { year: 2010, species: "LP", preference: 5 }
```
*From 2010, loblolly pine gets a +5 removal preference, so thins take loblolly ahead of equally-sized trees of other species.*

#### LEAVESP
**What it does:** Marks a species as **protected (leave)** — a "leave" species is excluded from *every* thin: its trees count toward neither the class stocking nor the removal, so no cut touches them. `LEAVESP` is an alias of `SPLEAVE`. Use it to hold a species entirely out of harvest (e.g. protect a hardwood or a rare species while thinning conifers).
**Parameters:**
- `year`(1) — date the flag takes effect.
- `species`(2) — SPDECD species/group: a code/index, `0` = reset all species to *not* protected, `−N` = SPGROUP N.
- `leave_flag`(3) — `> 0` ⇒ protect (leave); `≤ 0` ⇒ clear the flag. Default protects when the species field is given.
**Example:**
```text
LEAVESP           2010        WO         1.
```
```yaml
- LEAVESP: { year: 2010, species: "WO" }
```
*Protects white oak from all thinning starting in 2010.*

#### SPLEAVE
**What it does:** Identical to `LEAVESP` (same handler, cut code 206): flags a species as a **leave species** that is excluded from every thinning operation. Provided as the alternate spelling FVS accepts.
**Parameters:**
- `year`(1) — date.
- `species`(2) — SPDECD species/group; `0` = reset all, `−N` = SPGROUP N.
- `leave_flag`(3) — `> 0` protect, `≤ 0` clear.
**Example:**
```text
SPLEAVE           2010        -2         1.
```
```yaml
- SPLEAVE: { year: 2010, species: -2 }
```
*Protects every species in [SPGROUP](#species-groups) 2 from thinning starting in 2010.*

#### CUTEFF
**What it does:** Sets the **default cut efficiency** — the proportion of the **selected** trees that a thin actually removes when its own `cut_efficiency` field is left blank. "Selected trees" are the stems a thinning rule picked for removal: those inside the DBH/height window and passing the species filter (and not protected by `SPLEAVE`), ranked by the method. An efficiency of 1.0 removes the whole selection to the target; 0.8 models an imperfect operation that removes 80 % of each marked stem's trees-per-acre and misses 20 %. It also supplies the default proportion for the blank-field forms of THINPRSC/THINAUTO (and TOPKILL/HTGSTOP damage proportions).
**Parameters:**
- `proportion`(1) — default cut/affect fraction, 0–1. FVS default is 1.0.
**Example:**
```text
CUTEFF             0.9
```
```yaml
- CUTEFF: { proportion: 0.9 }
```
*Makes every later thin that omits its efficiency field remove 90 % of the trees it selects instead of 100 %.*

#### MINHARV
**What it does:** Sets **minimum-harvest thresholds** that gate a cut: after a cutting cycle's thins are computed, if the total removal falls **below any** threshold that is set (basal area, total cubic, merch cubic, sawlog cubic, or board feet), the *entire* cut is **canceled** and the stand restored to pre-thin. This models "don't bother logging unless the sale is big enough to be worthwhile." Thresholds persist across cycles once set; defaults of 0 leave the gate off.
**Parameters:**
- `year`(1) — date the thresholds take effect.
- `ba_min`(2) — minimum removed basal area, ft²/acre.
- `tcf_min`(3) — minimum removed total cubic volume, ft³/acre.
- `cf_min`(4) — minimum removed merchantable cubic volume, ft³/acre.
- `scf_min`(5) — minimum removed sawlog cubic volume, ft³/acre.
- `bf_min`(6) — minimum removed board-foot volume, bd ft/acre.
**Example:**
```text
MINHARV           2010        20.
```
```yaml
- MINHARV: { year: 2010, min_volume: 20 }
```
*From 2010, any scheduled cut that would remove less than 20 ft²/acre of basal area is canceled outright.* (The YAML `min_volume` maps to field 2 / `ba_min`; set the other minimums positionally.)

#### SALVAGE  *(base model — inert)*
**What it does:** The base-model `SALVAGE` thinning keyword is **abandoned** in FVS (`cuts.f`) and is a recognized **no-op** in FVSjl — accepted without error so it doesn't fall through as an unknown keyword, but it removes nothing. The functioning salvage operation (removing standing dead trees / snags after mortality or fire) lives in the Fire & Fuels Extension; see the FFE `SALVAGE`/`SALVSP` keywords under [Fire & fuels](#fire--fuels-ffe).
**Parameters:** *(none acted on — inert in the base model.)*
**Example:**
```text
SALVAGE           2010
```
```yaml
- keyword: "SALVAGE"
  params: ["2010"]
```
*Recognized but has no effect in the base model; use the FFE SALVAGE keyword (inside an FMIN block) to actually salvage dead trees.*

#### YARDLOSS
**What it does:** Models **yarding / logging loss** — a fraction of the harvested volume that is left on site during extraction instead of being hauled out. It scales the **reported** removed merchantable, sawlog, and board-foot volumes by `(1 − PRLOST)` (total cubic and basal area still reflect the full physical removal), and routes the lost material into the fuel/down-wood pools, split into downed vs standing pieces. Used to make harvest-yield reporting reflect real extraction losses.
**Parameters:**
- `year`(1) — date the loss fractions take effect.
- `prlost`(2) — proportion of the harvested merch/saw/board volume lost in yarding (left on site), 0–1. Default 0 (inactive).
- `prdsng`(3) — of that yarding loss, the proportion left as **downed** material; the remainder `(1 − prdsng)` becomes standing snags. 0–1.
- `prcrwn`(4) — the crown proportion of the loss, 0–1.
**Example:**
```text
YARDLOSS          2015        0.1        0.7        0.5
```
```yaml
- keyword: "YARDLOSS"
  params: ["2015", "0.1", "0.7", "0.5"]
```
*From 2015, 10 % of harvested merchantable/sawlog/board volume is lost in yarding (70 % of it downed, 30 % left as snags), reducing reported yields and adding to the down-wood pools.*
## Establishment & regeneration

The Southern (SN) variant — like NE/CS/LS — implements a **partial, keyword-driven** establishment model: FVSjl does not spontaneously add ingrowth on its own. New trees enter the stand only where you explicitly schedule them inside an **establishment packet**. A packet is the block of cards bracketed by `ESTAB … END`: `ESTAB` opens it and carries the *date of the triggering disturbance* (a harvest, fire, site-prep, etc.); the `PLANT`, `NATURAL`, and `SPROUT` cards inside it each schedule one regeneration activity relative to that disturbance; `END` closes the packet. The bracket matters because the sub-keywords (`PLANT`/`NATURAL`/`SPROUT`/`NOSPROUT`/`TALLY`) are *establishment-extension* keywords read by the ESIN parser only **inside** a packet — a bare `SPROUT` at top level is an "INVALID KEYWORD" in stock FVS (FVSjl still tolerates the standalone `NOSPROUT`/`NOAUTOES` disable forms).

#### ESTAB
**What it does:** Opens an establishment packet and reads cards until `END`. Field 1 is the *date of the disturbance* that this regeneration responds to (`IDSDAT`), which anchors the timing of every `PLANT`/`NATURAL` inside (age/height are grown from it via the ESNUTR/ESSUBH delay-and-generation-time logic). At `END`, FVSjl schedules a `TALLY` establishment trigger (activity 427) at the disturbance date and marks the packet's date sentinel so downstream defaulting applies. The packet is what makes regeneration *conditional on a disturbance year* rather than a free-standing planting.
**Parameters:**
- `date`(1) — calendar year (≥ 1000) or cycle number (< 1000) of the disturbance the packet responds to; blank ⇒ treated as unset (defaulted downstream). No other fields.
**Example:**
```text
ESTAB         2000.0
PLANT         2000.0      LP       300.0      90.0
NATURAL       2000.0      WO        50.0      80.0       3.0       1.0
END
```
```yaml
- regeneration:
    - ESTAB: { date: 2000 }
    - PLANT: { year: 2000, species: "LP", tpa: 300, survival_pct: 90 }
    - NATURAL: { year: 2000, species: "WO", tpa: 50, survival_pct: 80, age: 3, height: 1 }
    - END: {}
```
*Opens a packet keyed to a year-2000 disturbance, plants 300 loblolly pine/ac at 90 % survival, and adds 50 natural white-oak seedlings/ac, then closes the packet.*

#### PLANT
**What it does:** Schedules a **planting** — an operator-established cohort (`MANAGD`/artificial regeneration) of a chosen species at a given date. FVSjl creates the seedling records at that cycle: it spreads `tpa · survival_pct/100` across the inventory point replicates, drawing an established height per replicate from the species' establishment height-at-age curve (ESSUBH) plus a BACHLO random draw, then back-solves DBH. A planted cohort also flips the managed-stand flag that feeds the `dg_planted` diameter-growth term. Species must resolve (alpha code like `LP`, an FIA number, or a numeric sequence index) or the card silently plants nothing.
**Parameters:**
- `year`(1) — planting year (≥ 1000) or cycle number (< 1000); blank ⇒ cycle 1.
- `species`(2) — species to plant: 2-letter alpha code (e.g. `LP`), FIA number, or numeric sequence index (SPDECD-decoded). Required; unresolved ⇒ no trees.
- `tpa`(3) — trees per acre established (before survival). **Must be > 0** or the card is skipped (esin.f:143).
- `survival_pct`(4) — percent surviving to the tally, applied as `tpa·survival/100`. Out-of-range (< 0.001 or > 100) ⇒ defaults to **100 %** (esin.f:149).
- `age`(5) — age of the planting stock at establishment, years; default 0 (seedling). Feeds the ESSUBH base-height computation.
- `height`(6) — override established height, feet; default 0 ⇒ model-derived height.
- `shade`(7) — shade/tolerance code carried for the seedling; default 0.
**Example:**
```text
ESTAB         2005.0
PLANT         2005.0      LP       435.0      85.0       1.0
END
```
```yaml
- regeneration:
    - ESTAB: { date: 2005 }
    - PLANT: { year: 2005, species: "LP", tpa: 435, survival_pct: 85, age: 1 }
    - END: {}
```
*Plants 435 loblolly pine/ac of age-1 stock in 2005 at 85 % survival (≈ 370 surviving seedlings/ac).*

#### NATURAL
**What it does:** Schedules **natural regeneration** — seedlings arriving on their own (seed-fall/advance regeneration) after the disturbance, as opposed to operator planting. It uses the identical field layout and creation path as `PLANT`, but does **not** set the managed-stand flag, and its established heights are drawn on the natural-height acceptance window (estab.f) rather than the planted path. Use it to represent expected volunteer ingrowth of a species you did not plant.
**Parameters:**
- `year`(1) — regeneration year (≥ 1000) or cycle number (< 1000); blank ⇒ cycle 1.
- `species`(2) — species (alpha/FIA/numeric index, SPDECD-decoded). Required.
- `tpa`(3) — trees per acre; **must be > 0** or skipped.
- `survival_pct`(4) — percent surviving; out-of-range ⇒ 100 %.
- `age`(5) — seedling age at establishment, years; default 0.
- `height`(6) — override height, feet; default 0 ⇒ model-derived.
- `shade`(7) — shade/tolerance code; default 0.
**Example:**
```text
ESTAB         2000.0
NATURAL       2000.0      WO        80.0      75.0       2.0
NATURAL       2000.0      RM       120.0     100.0
END
```
```yaml
- regeneration:
    - ESTAB: { date: 2000 }
    - NATURAL: { year: 2000, species: "WO", tpa: 80, survival_pct: 75, age: 2 }
    - NATURAL: { year: 2000, species: "RM", tpa: 120, survival_pct: 100 }
    - END: {}
```
*Adds 80 white-oak/ac (age 2, 75 % survival) and 120 red-maple/ac (100 % survival) as natural ingrowth after a year-2000 disturbance.*

#### SPROUT
**What it does:** Enables (and parameterizes) **stump sprouting** — vegetative regeneration where cut or killed hardwood stumps re-sprout. It is an establishment sub-keyword read inside an `ESTAB` packet (esin.f option 26). By itself it sets the sprouting flag `LSPRUT = true`; its fields optionally scale the sprout **count** and **height** and restrict which parent-stump diameters sprout. Species selection is variant-aware: FVSjl validates the chosen species (or every member of a species group) against the variant's known sprouters (`is_sprouting`/ISPSPE) — a species that cannot sprout, a group where not all members sprout, or the −999 "no-species" sentinel instead **disables** sprouting. A bare `SPROUT <date>` (blank species ⇒ 0 ⇒ ALL) enables all-species sprouting.
**Parameters:**
- `date`(1) — activation date (year or cycle).
- `species`(2) — SPDECD selector: `0`/blank = all sprouting species, `−N` = SPGROUP group *N*, `> 0` = a single species; must be a valid sprouter or sprouting is turned off.
- `smult`(3) — sprout-count multiplier applied to the parent; default **1**.
- `hmult`(4) — sprout-height multiplier; default **1**.
- `dbh_lo`(5) — lower parent-stump DBH bound for the override window, inches; default **0**.
- `dbh_hi`(6) — upper parent-stump DBH bound, inches; default **999**.
**Example:**
```text
ESTAB         2010.0
SPROUT        2010.0      WO         3.0       1.0       0.0      20.0
END
```
```yaml
- regeneration:
    - ESTAB: { date: 2010 }
    - SPROUT: { year: 2010, species: "WO", smult: 3, hmult: 1, dbh_lo: 0, dbh_hi: 20 }
    - END: {}
```
*Enables stump sprouting in 2010, tripling the sprout count of white-oak stumps 0–20 in DBH (with default sprout height).*

#### NOSPROUT
**What it does:** Disables stump sprouting (`LSPRUT = false`, esin.f option 27). It is the explicit "off" switch for the sprouting process — the counterpart to `SPROUT`. FVSjl accepts it both **inside** an `ESTAB` packet (its native ESIN position) and as a standalone top-level card.
**Parameters:** *(none)* — a date field, if present, is ignored.
**Example:**
```text
NOSPROUT
```
```yaml
- regeneration:
    - NOSPROUT: {}
```
*Turns off all stump sprouting for the run.*

#### NOAUTOES
**What it does:** Disables **automatic establishment** — the model's would-be spontaneous ingrowth/auto stump-sprouting (estab.f `LFLAG`). In FVSjl this clears the same sprouting/auto-establishment flag (`LSPRUT = false`), so no trees are added except those you schedule with `PLANT`/`NATURAL`/`SPROUT`. Note the SN variant is already a *partial, keyword-driven* model with no full automatic ingrowth, so `NOAUTOES` chiefly documents intent and suppresses any auto stump-sprouting; the deliberately scheduled regeneration in a packet is unaffected.
**Parameters:** *(none)*
**Example:**
```text
NOAUTOES
```
```yaml
- regeneration:
    - NOAUTOES: {}
```
*Suppresses automatic establishment so only explicitly scheduled regeneration enters the stand.*

> **Order-significance.** Within a packet the cards apply in the order written, and the bracket is mandatory: `ESTAB` first, then any mix of `PLANT`/`NATURAL`/`SPROUT`/`NOSPROUT`, then `END`. A `SPROUT`/`NOSPROUT` that both appear resolve to whichever came **last**. `PLANT`/`NATURAL` cards with the same year accumulate (multiple species/cohorts). A `SPROUT` naming a `−N` group requires that `SPGROUP` group to have been defined earlier in the file.

## Species groups

A **named species group** lets you refer to several species at once with a single negative code. You define the group with `SPGROUP`, and thereafter any keyword's species field written as `−N` means "every member of the *N*-th group defined" (groups are numbered in definition order). This is the mechanism behind e.g. thinning "all oaks and hickories" with one card.

#### SPGROUP
**What it does:** Defines a species group (vbase/initre.f:4726). The keyword card's field 1 is the (optional) group **name**; the group's members are listed on the **next record** as space-separated species (alpha codes like `SM HI`, FIA numbers, or numeric indices). `ALL`/`0`/duplicates are skipped; up to 90 members, up to 30 groups total; a blank name auto-generates `GROUPnn`. Once defined, the group is referenced from any species field by its negative index `−N` (first group = −1, second = −2, …). Because the reference is by definition order, **every `SPGROUP` must be declared before the keyword that names it.**
**Parameters:**
- `group`(1) — group name (optional); blank ⇒ auto-named `GROUP01`, `GROUP02`, …
- *(next record)* — the member list: space-separated species codes (alpha/FIA/numeric); `ALL`, `0`, and duplicates are dropped; max 90 members.
**Example:**
```text
SPGROUP   OAKHIC
WO RO BO SO HI
THINDBH       2010.0       0.0      99.0      1.00        -1     100.0
```
```yaml
- species_groups:
    - SPGROUP: { group: "OAKHIC" }
    - raw: "WO RO BO SO HI"        # member species, carried verbatim
- treatments:
    - keyword: "THINDBH"           # species -1 references OAKHIC above
      params: ["2010.0", "0.0", "99.0", "1.00", "-1", "100.0"]
```
*Defines the group `OAKHIC` (white/red/black/scarlet oak + hickory) and thins that group across all diameters in 2010; the `SPGROUP` block must come before the `THINDBH −1` that references it.*

## Site & treatments

These keywords change the **site conditions or growth environment** partway through a run, rather than the tree list directly. Schedule them at a date and they fire in the cycle that contains that date.

#### SETSITE
**What it does:** Schedules a **mid-run change to site parameters** (initre.f:13800, activity 120). You would schedule it to represent a site becoming more/less productive over time, or to correct a site index for a species after an event — anything that should change the growth potential without editing the tree records. When the scheduled cycle arrives, FVSjl updates the per-species **site index** (directly or as a percent change), and optionally the maximum basal area (`BAMAX`) and maximum SDI (`SDImax`), then re-seeds the site-dependent diameter-growth constants (the FVS `RCON`) so subsequent growth uses the new site. The habitat-type field (2) is parsed but not yet wired in the SN family (SN growth keys off forest type, not habitat) — a documented gap.
**Parameters:**
- `year`(1) — date of the change (calendar year ≥ 1000 or cycle number < 1000); blank ⇒ cycle 1.
- `habitat`(2) — habitat-type code; parsed but currently ignored in SN (default 0).
- `bamax`(3) — new maximum basal area (sq ft/ac); 0 ⇒ leave unchanged. Also updates the BAMAX-derived SDImax default.
- `species`(4) — species the change applies to (SPDECD: `0`/blank = all, `−N` = group, else a code).
- `site_index`(5) — new site index (feet); interpretation set by `si_flag`. Clamped to ≥ 1.
- `si_flag`(6) — 0 ⇒ field 5 is a **direct** site-index value; non-zero ⇒ field 5 is a **percent change** applied to the current site index.
- `sdimax`(7) — new maximum SDI; 0 ⇒ leave unchanged (overrides the BAMAX-derived value if both given).
**Example:**
```text
SETSITE       2020.0                             0        75.0
```
```yaml
- site:
    - SETSITE: { year: 2020, site_index: 75, si_flag: 0 }
```
*In 2020, sets the site index of all species to 75 (direct value) and rebuilds the diameter-growth constants for the new site.*

#### FERTILIZ
**What it does:** Schedules a **fertilizer application** (ffin.f/ffert.f), modeling a 200-lb-nitrogen dose. When it fires, FVSjl applies a multiplicative boost to each tree's squared-diameter change (DDS) and height growth for up to **10 years** after application, scaled by the application efficacy and prorated to the cycle years falling inside that window; the effect carries across cycles until it expires. Only nitrogen at the single representable amount (200 lb) is modeled — the phosphorus/potassium fields are read but ignored — and only the efficacy multiplier actually feeds the response (activity 260). (SN lies outside the model's calibrated Douglas-fir/grand-fir range; stock FVS warns but still applies the species-agnostic factor, and FVSjl matches that.)
**Parameters:**
- `year`(1) — application date (year or cycle); blank ⇒ cycle 1.
- `nitrogen`(2) — N applied; forced to **200 lb** (the only representable amount).
- *(field 3 = P, field 4 = K)* — parsed but **ignored**.
- `efficacy`(5) — application efficacy multiplier scaling the response; default **1**.
**Example:**
```text
FERTILIZ      2015.0     200.0       0.0       0.0       1.0
```
```yaml
- site:
    - FERTILIZ: { year: 2015, nitrogen: 200 }
```
*Applies 200 lb N/ac in 2015 at full efficacy, boosting diameter and height growth for the following 10 years.*

*Both `SETSITE` and `FERTILIZ` are scheduled activities: a date that lands mid-cycle fires in the cycle **containing** it (OPCYCL bucketing), and the change persists from the cycle start — so place them at the year you want the effect to begin.*
## Volume & merchandising

FVSjl reports four stem volumes per tree, and these keywords control how each is
computed. **Total cubic** is the whole-stem cubic-foot wood volume. **Merch(antable)
cubic** and **sawtimber cubic** are the cubic-foot volumes of the *merchantable* portion
of the stem — the piece left after you cut a **stump** off the bottom and stop at a
**top diameter inside bark** (DIB), and only if the tree is at least a **minimum DBH**.
Those three numbers (stump height, top DIB, min DBH) are the *cubic merchandising
standards*. **Board feet** is a lumber-yield measure (Scribner or International ¼″ log
rule) computed on the sawtimber-sized part of the stem — always larger minimums than
cubic. A **volume equation number** picks which taper/volume model (an NVEL string such
as a Clark equation) turns DBH + height into that cubic profile. A **defect curve**
deducts a percentage of the *gross* volume as cull/rot, given as a fraction at each
5″ DBH class (interpolated between classes); cubic defect (`MCDEFECT`) trims the
pulpwood/topwood part while board defect (`BFDEFECT`) trims board feet *and* sawtimber
cubic. The **log-linear form model** (`MCFDLN`/`BFFDLN`) is an alternative defect source:
a two-parameter correction `VOLCOR = exp(B0 + B1·ln V)` whose implied shrinkage competes
(as a max) with the DBH curve and the per-tree input defect.

#### VOLUME
**What it does:** Overrides the **cubic** merchandising standards — for the merch-cubic
product and, separately, the sawtimber-cubic product — for one species, a species group,
or all species (`volkey.f`, activity 217). It resets, per species, the minimum DBH a tree
must reach to yield merch cubic, the top diameter inside bark the merchantable bole is
cut back to, and the stump height, plus the three sawtimber-cubic counterparts. These
standards feed the Clark taper integration (and the broken-top `CFTOPK` re-fit), so
lowering the top DIB or min DBH makes more of each stem merchantable and raises reported
`Mcuft`/`Scuft`. Field 7 (`method`) selects the cubic model — 6 = Clark (the default and
the only model SN/NE use), 5 = the Central-States Gevorkiantz/DVEE model (CS variant
only). Undated cards take effect immediately; a dated card is scheduled to its cycle.
**Parameters:**
- `year`(1) — schedule year (≥1000 = calendar year, <1000 = cycle number); blank ⇒ cycle 1.
- `species`(2) — `0`/`ALL` = every species, a positive alpha/FIA/PLANTS code = one species, `−N` = [SPGROUP](#species-groups) *N*. Default `0`.
- `dbh_min`(3) — DBHMIN, minimum DBH (in) for a tree to yield **merch cubic**. SN default ≈ 4″ (softwoods) / 4–6″ (some hardwoods); higher on North-Carolina forests.
- `top_diam`(4) — TOPD, **merch-cubic** top diameter inside bark (in). Default 4.0″.
- `stump`(5) — STMP, **merch-cubic** stump height (ft). Default 0.5.
- `form_class`(6) — FRMCLS, Girard form class. Not read by the Clark taper path (inert in SN/NE/CS/LS).
- `method`(7) — METHC, cubic-volume method: 6 = Clark (default), 5 = CS Gevorkiantz/DVEE (Central-States only; ignored elsewhere).
- `scf_min_dbh`(8) — SCFMIND, minimum DBH (in) to be **sawtimber**. SN default ≈ 10″ (softwood) / 12″ (hardwood).
- `scf_top_dib`(9) — SCFTOPD, **sawtimber-cubic** top diameter inside bark (in). Default ≈ 7″ (softwood) / 9″ (hardwood).
- `scf_stump`(10) — SCFSTMP, **sawtimber-cubic** stump height (ft). Default 1.0.
**Example:**
```text
VOLUME      2010.0       0.0       5.0       4.0       0.5                          11.0       8.0       1.0
```
```yaml
- VOLUME:
    year: 2010
    species: 0            # 0 = all species
    dbh_min: 5.0          # merch-cubic min DBH (in)
    top_diam: 4.0         # merch-cubic top DIB (in)
    stump: 0.5            # merch-cubic stump height (ft)
    scf_min_dbh: 11.0     # sawtimber-cubic min DBH (in)
    scf_top_dib: 8.0      # sawtimber-cubic top DIB (in)
    scf_stump: 1.0        # sawtimber-cubic stump height (ft)
```
*From 2010 on, every species must reach 5″ DBH for merch cubic and 11″ DBH for sawtimber, with sawtimber cut back to an 8″ top — a larger-log utilization standard (`form_class`/`method` fields 6–7 left blank keep their defaults).*

#### BFVOLUME
**What it does:** Overrides the **board-foot** merchantability standards — the minimum
DBH, top diameter inside bark, and stump height used for the Scribner/International board
computation — for a species, group, or all species (`volkey.f`, activity 218). Board
minimums are normally set equal to the sawtimber-cubic minimums ("bf-equal") in the
eastern variants, so `BFVOLUME` is what you use to make the board rule *differ* from the
sawtimber-cubic rule. When the board standards (or the board equation via `VOLEQNUM`)
diverge from the sawtimber ones, FVS recomputes board feet from a separate board taper
call rather than riding the sawtimber call. The FVSjl SN Clark path carries only the
three standards; the `form_class`/`method` fields are read but unused by that taper.
**Parameters:**
- `year`(1) — schedule year; blank ⇒ cycle 1.
- `species`(2) — `0`/`ALL`, alpha/FIA code, or `−N` for SPGROUP *N*. Default `0`.
- `bf_min_dbh`(3) — BFMIND, minimum DBH (in) for board feet. Default = the sawtimber-cubic min (≈ 10″ softwood / 12″ hardwood).
- `bf_top_dib`(4) — BFTOPD, board-foot top diameter inside bark (in). Default ≈ 7″.
- `bf_stump`(5) — BFSTMP, board-foot stump height (ft). Default 1.0.
- `form_class`(6) — FRMCLS, form class; read but not used by the Clark taper path.
- `method`(7) — METHB, board-volume method; read but not used by the Clark taper path.
**Example:**
```text
BFVOLUME    2010.0       0.0      12.0       8.0       1.0
```
```yaml
- BFVOLUME:
    year: 2010
    species: 0
    bf_min_dbh: 12.0      # board min DBH (in)
    bf_top_dib: 8.0       # board top DIB (in)
    bf_stump: 1.0         # board stump height (ft)
```
*Requires a 12″ DBH and cuts boards back to an 8″ top for every species from 2010 — a stricter board rule than the sawtimber-cubic standard, so board feet are recomputed from their own taper call.*

#### VOLEQNUM
**What it does:** Overrides the NVEL cubic volume-equation identifier for a species, so
that species' cubic profile is built from a different taper/volume model than the one
`VOLEQDEF`/the location default assigned (`initre.f:5061`). The equation id is the 10-character
NVEL string (a Clark equation such as `841CLKE318`, where the digits encode
species/region/model). The override is stored and applied **after** the defaults are set,
so it wins. Because board feet ride the sawtimber cubic call by default, changing the
cubic equation also changes board feet unless a separate board equation/standard splits
them (see `BFVOLUME`). This keyword has **no date field** — field 1 is the species.
**Parameters:**
- `species`(1) — the species to override: alpha/FIA/PLANTS code, or `−N` for SPGROUP *N* (applies the id to every member).
- `equation`(2) — the NVEL volume-equation id, a text token up to 10 characters (e.g. `841CLKE318`); a bare number is accepted too. Blank ⇒ the card is ignored.
**Example:**
```text
VOLEQNUM         131  841CLKE318
```
```yaml
- VOLEQNUM: { species: 131, equation: "841CLKE318" }   # 131 = loblolly pine (FIA); NVEL Clark cubic eqn
```
*Forces loblolly pine to use NVEL Clark equation `841CLKE318` instead of its location default, changing that species' cubic (and, by default, board) volumes.*

#### MCDEFECT
**What it does:** Sets a per-species **cubic** defect curve (`CFDEFT`, `sdefet.f`,
activity 215): the fraction of *gross* cubic volume treated as cull/rot at DBH classes
5, 10, 15, 20 and 25″, interpolated linearly between classes and held flat above 25″
(out to 30/35/40″) and anchored at 0% at DBH 0. In the volume path (`vols.f`) the curve
value is turned into an integer percent `ICDF` and cuts the **pulpwood/topwood** cubic
(the merch-cubic-minus-sawtimber part, `MCFV − SCFV`); the sawtimber part is left for the
board step. `ICDF` is the **largest** of three sources: this DBH curve, the per-tree
input defect, and the `MCFDLN` log-linear form model — so `MCDEFECT` only raises the cull,
never lowers it. Blank fields are ALGSLP-filled (clamped/interpolated from the supplied
points), **not** taken as 0. An undated card applies now; a dated card is scheduled to
its cycle.
**Parameters:**
- `year`(1) — schedule year; blank/undated ⇒ apply immediately.
- `species`(2) — `0`/`ALL`, alpha/FIA code, or `−N` for SPGROUP *N*. Default `0`.
- `defect_5`(3) — cull **fraction** (0.0–1.0) at DBH 5″.
- `defect_10`(4) — cull fraction at DBH 10″.
- `defect_15`(5) — cull fraction at DBH 15″.
- `defect_20`(6) — cull fraction at DBH 20″.
- `defect_25`(7) — cull fraction at DBH 25″ (extended flat to 30/35/40″). Default curve is all 0 (no defect).
**Example:**
```text
MCDEFECT    2010.0       0.0      0.02      0.05      0.10      0.15      0.20
```
```yaml
- MCDEFECT: { year: 2010, species: 0, defect_5: 0.02, defect_10: 0.05, defect_15: 0.10, defect_20: 0.15, defect_25: 0.20 }
```
*Deducts 2%–20% of gross topwood/pulpwood cubic across the 5″→25″ DBH classes (rising with tree size) for all species from 2010.*

#### BFDEFECT
**What it does:** Sets a per-species **board-foot** defect curve (`BFDEFT`, `sdefet.f`,
activity 216) with the same 5/10/15/20/25″ layout, interpolation, and flat extension as
`MCDEFECT`. The difference is what it reduces: the resulting percent `IBDF` cuts **both
board feet and sawtimber cubic** (`vols.f`), and because reported merch cubic includes
the surviving sawtimber, a board defect also lowers merch cubic. An input board-defect is
applied to sawtimber cubic even for trees too small to yield board feet, while the
curve/form contributions only update where board feet exist. `IBDF` is again the max of
the DBH curve, the per-tree input defect, and the `BFFDLN` form model. Board defect
fractions are typically higher than cubic ones because sawtimber grading is stricter.
**Parameters:**
- `year`(1) — schedule year; blank/undated ⇒ apply immediately.
- `species`(2) — `0`/`ALL`, alpha/FIA code, or `−N` for SPGROUP *N*. Default `0`.
- `defect_5`(3) — cull **fraction** (0.0–1.0) at DBH 5″.
- `defect_10`(4) — cull fraction at DBH 10″.
- `defect_15`(5) — cull fraction at DBH 15″.
- `defect_20`(6) — cull fraction at DBH 20″.
- `defect_25`(7) — cull fraction at DBH 25″ (extended flat to 30/35/40″). Default curve is all 0.
**Example:**
```text
BFDEFECT    2010.0      -2.0      0.05      0.10      0.15      0.20      0.25
```
```yaml
- BFDEFECT: { year: 2010, species: -2, defect_5: 0.05, defect_10: 0.10, defect_15: 0.15, defect_20: 0.20, defect_25: 0.25 }
```
*Culls 5%–25% of board feet (and sawtimber cubic) with DBH for the species in group 2 from 2010 — a steeper defect than the cubic curve, reflecting stricter sawlog grading.*

#### MCFDLN
**What it does:** Sets the two **cubic** log-linear form-model coefficients B0/B1
(`CFLA0`/`CFLA1`, `sdefln.f`, option 39), activating an alternate cubic defect source.
For each tree the model computes a corrected volume `VOLCOR = exp(B0 + B1·ln(V))` from the
topwood cubic `V = MCFV − SCFV`, and the implied shrinkage `(V − VOLCOR)/V` (as a percent)
competes — as a maximum — with the `MCDEFECT` DBH curve and the per-tree input defect for
`ICDF`. The defaults B0 = 0, B1 = 1 give `VOLCOR = V`, i.e. no correction, so the model is
inert until you set it. This keyword has **no date field** — field 1 is the species — and
only the fields you supply are written (a blank B0 or B1 keeps its default).
**Parameters:**
- `species`(1) — `0`/`ALL`, alpha/FIA code, or `−N` for SPGROUP *N*.
- `b0`(2) — B0, intercept of the log-linear form model. Default 0.
- `b1`(3) — B1, slope on `ln(V)`. Default 1.
**Example:**
```text
MCFDLN           0.0  -0.10536       1.0
```
```yaml
- keyword: "MCFDLN"
  params: ["0.0", "-0.10536", "1.0"]   # species=0, B0=-0.10536, B1=1.0
```
*With B1 = 1 and B0 = −0.10536, `VOLCOR = V·e^-0.10536 ≈ 0.90 V`, i.e. a ~10% topwood-cubic reduction for all species — applied only where it exceeds the tree's other defect sources.*

#### BFFDLN
**What it does:** Sets the two **board-foot** log-linear form-model coefficients B0/B1
(`BFLA0`/`BFLA1`, `sdefln.f`, option 40) — the board-side analogue of `MCFDLN`. It
computes `VOLCOR = exp(B0 + B1·ln(BFV))` from the tree's board feet, and the implied
percent reduction competes (as a max) with the `BFDEFECT` DBH curve and the per-tree input
defect for `IBDF`, which then trims board feet and sawtimber cubic. Defaults B0 = 0,
B1 = 1 make it a no-op. No date field — field 1 is the species; only supplied fields are
written.
**Parameters:**
- `species`(1) — `0`/`ALL`, alpha/FIA code, or `−N` for SPGROUP *N*.
- `b0`(2) — B0, intercept of the board log-linear form model. Default 0.
- `b1`(3) — B1, slope on `ln(BFV)`. Default 1.
**Example:**
```text
BFFDLN          -1.0  -0.22314       1.0
```
```yaml
- keyword: "BFFDLN"
  params: ["-1.0", "-0.22314", "1.0"]   # SPGROUP 1, B0=-0.22314, B1=1.0
```
*Applies `VOLCOR = BFV·e^-0.22314 ≈ 0.80 BFV`, a ~20% board-foot (and sawtimber-cubic) reduction for every species in group 1, taken whenever it beats that species' other board-defect sources.*
## Output & database

#### ECHOSUM
**What it does:** `ECHOSUM` requests that FVS write the stand-summary table (the `.sum` file) to output. The `.sum` file is FVS's canonical fixed-column projection report: one row per cycle giving the start-of-period and after-treatment stand statistics — year, age, trees/acre, basal area, SDI, CCF, top height, QMD, total and merchantable cubic and board-foot volumes, periodic and mean annual increment, mortality, and the sampling weight. FVSjl always computes and can emit this table (see `src/io/summary.jl`), so `ECHOSUM` is recognized as a report-control no-op — the summary is produced regardless — and it is preserved losslessly through the `.key`↔`.yaml` round-trip. The modern CSV form of the same table carries the extra `Title` column that the fixed-column `.sum` drops.
**Parameters:** *(none)* — a bare keyword that only toggles emission of the summary table.
**Example:**
```text
ECHOSUM
```
```yaml
- ECHOSUM: {}
```
*Requests the per-cycle `.sum` stand-summary table for the stand.*

#### DATABASE
**What it does:** `DATABASE` opens the DBS (DataBase Structured output) block, read line-by-line until a closing `END` (dbsin.f). It directs FVS to emit its outputs into a SQLite database file instead of / in addition to the text reports, and it can also pull the stand and tree list *in* from a FIA "FVS-ready" SQLite database. Inside the block, sub-keywords select which tables are written (`SUMMARY`→FVS_Summary, `TREELIDB`→FVS_TreeList, `COMPUTDB`→FVS_Compute, plus `CUTLIST` for FVS_CutList) and name the output/input files (`DSNOUT`/`DSNIN`); `StandSQL`/`TreeSQL` … `EndSQL` provide the SQL queries used to read an input database. Each table shares a `CaseID`/`StandID` key so a run's stand-, tree-, cut-, carbon-, fuels-, and fire-event tables all join together.
**Parameters:** block contents (read until `END`):
  - `DSNOUT` — the *next raw line* is the output SQLite filename.
  - `DSNIN` — the next raw line is the input SQLite filename to load the stand/tree data from.
  - `STANDSQL` / `TREESQL` — SQL query text over the following raw lines up to a line reading `EndSQL`; used with `DSNIN` to load the FIA stand + tree list.
  - `SUMMARY` — enable the FVS_Summary table (the `.sum` columns, per cycle).
  - `TREELIDB` — enable the FVS_TreeList per-tree snapshot table.
  - `COMPUTDB` — enable the FVS_Compute event-monitor-variable table.
**Example:**
```text
DATABASE
DSNOUT
FVSOut.db
SUMMARY
TREELIDB
COMPUTDB
END
```
```yaml
- output:
    - DATABASE: {}
    - DSNOUT: { filename: "FVSOut.db" }
    - SUMMARY: {}
    - TREELIDB: {}
    - COMPUTDB: {}
    - END: {}
```
*Writes the run's per-cycle summary, per-tree list, and COMPUTE variables to the SQLite file `FVSOut.db`.*

#### DSNOUT
**What it does:** `DSNOUT` names the output SQLite database file for the `DATABASE` block. It is a sub-keyword of `DATABASE` (not a top-level keyword): FVS reads the filename from the raw line that *follows* the `DSNOUT` card. Every DBS table selected in the block (FVS_Summary, FVS_TreeList, FVS_CutList, FVS_Compute, and the FFE fire/fuels/carbon tables) is written into this one file. Without `DSNOUT` the DBS block selects tables but has no destination file.
**Parameters:** filename (next raw line) — path to the SQLite database to create/append.
**Example:**
```text
DATABASE
DSNOUT
/runs/stand42_out.db
SUMMARY
END
```
```yaml
- output:
    - DATABASE: {}
    - DSNOUT: { filename: "/runs/stand42_out.db" }
    - SUMMARY: {}
    - END: {}
```
*Sends the FVS_Summary table to the SQLite file `/runs/stand42_out.db`.*

#### SUMMARY
**What it does:** `SUMMARY` (inside a `DATABASE` block) enables the **FVS_Summary** table, the SQLite mirror of the text `.sum` file. Its columns are exactly the `.sum` columns FVSjl already computes for each cycle — year, age, TPA, BA, SDI, CCF, top height, QMD, total/merch cubic and board volumes, PAI/MAI, mortality, and sample weight — one row per period, keyed by `CaseID`/`StandID`. It is the DBS analog of `ECHOSUM`: `ECHOSUM` controls the text report, `SUMMARY` controls the database table.
**Parameters:** *(none)* — bare sub-keyword; sets the internal `ISUMARY=1` flag.
**Example:**
```text
DATABASE
DSNOUT
out.db
SUMMARY
END
```
```yaml
- output:
    - DATABASE: {}
    - DSNOUT: { filename: "out.db" }
    - SUMMARY: {}
    - END: {}
```
*Enables the per-cycle FVS_Summary stand-statistics table in the SQLite output.*

#### TREELIDB
**What it does:** `TREELIDB` (inside a `DATABASE` block) enables the **FVS_TreeList** table: a per-tree snapshot for each cycle. Each row is one tree record with its identity and full state — tree id and index, species in three code systems (FVS/PLANTS/FIA), trees-per-acre and mortality-per-acre, DBH, diameter growth, height, height growth, crown ratio and crown width, BA percentile, point BAL, the four volumes (total cubic, merch cubic, sawlog cubic, board feet), truncated/broken-top height, the two merch-height-to-top-diameter values, and tree age. This is the database form of the FVS `TREELIST` text report.
**Parameters:** *(none)* — bare sub-keyword; sets `dbs_treelist = true`.
**Example:**
```text
DATABASE
DSNOUT
out.db
SUMMARY
TREELIDB
END
```
```yaml
- output:
    - DATABASE: {}
    - DSNOUT: { filename: "out.db" }
    - SUMMARY: {}
    - TREELIDB: {}
    - END: {}
```
*Adds a per-tree, per-cycle FVS_TreeList snapshot table alongside the summary in `out.db`.*

#### CUTLIST
**What it does:** `CUTLIST` requests the list of trees **removed** by thinning/harvest — the FVS_CutList table (dbscuts.f). Whenever a cut fires in a cycle, each removed record is captured and written as a row with the same fields as the tree list — tree id/index, species, removed TPA, DBH and its growth, height and height growth, crown ratio/width, BA percentile, point BAL, the four volumes, truncated height, the two merch heights, and tree age — keyed by cut year and period length. As a top-level keyword it enables the DBS FVS_CutList table (`dbs_cutlist = true`); the plain text-report form is a recognized report no-op.
**Parameters:** *(none)* — bare keyword; toggles emission of the FVS_CutList table.
**Example:**
```text
CUTLIST
THINBBA           2010      80.
```
```yaml
- output:
    - CUTLIST: {}
- treatments:
    - THINBBA: { year: 2010, residual_basal_area: 80 }
```
*Enables the FVS_CutList table so every tree removed by the 2010 thin-to-basal-area is recorded.*

#### COMPUTDB
**What it does:** `COMPUTDB` (inside a `DATABASE` block) enables the **FVS_Compute** table (dbscmpu.f): the per-cycle values of the event-monitor `COMPUTE` variables. Its schema is *dynamic* — one REAL column per named COMPUTE variable, in declaration order, created on first write — so the columns depend on the run's `COMPUTE` blocks. Each row is one growing cycle (`Year`), with a variable written NULL in cycles before it becomes active. This lets the user's derived quantities (e.g. a computed density index or ratio) be queried alongside the standard summary.
**Parameters:** *(none)* — bare sub-keyword; sets `dbs_compute = true`.
**Example:**
```text
DATABASE
DSNOUT
out.db
COMPUTDB
END
COMPUTE            0
STOCK = BBA / BSDI
END
```
```yaml
- output:
    - DATABASE: {}
    - DSNOUT: { filename: "out.db" }
    - COMPUTDB: {}
    - END: {}
- event_monitor:
    - COMPUTE: { date: 0 }
    - raw: "STOCK = BBA / BSDI"
    - END: {}
```
*Records the per-cycle value of the user variable `STOCK` into the FVS_Compute table.*

#### STRCLASS
**What it does:** `STRCLASS` activates the stand **structural-stage classification** (SSTAGE, ksstag.f) and can override its classification thresholds. The classifier assigns the stand to a structural stage (e.g. stand-initiation, stem-exclusion, understory-reinitiation, old-forest single/multi-story) from canopy cover, the size (DBH) of the uppermost stratum, and gap/sawtimber percentages. Turning it on makes six event-monitor variables available for conditional logic — `BSCLASS`/`ASCLASS` (structural class before/after thin), `BSTRDBH`/`ASTRDBH` (uppermost-stratum DBH), and `BCANCOV`/`ACANCOV` (canopy cover). FVSjl always computes the class when `STRCLASS` is present; the separate per-cycle text report is deferred.
**Parameters:**
- `print`(1) — print code; 0 = compute but don't print, else compute and print (FVSjl always computes). Default: report deferred.
- `GAPPCT`(2) — gap-percent threshold. Default: internal.
- `SSDBH`(3) — stem-exclusion DBH threshold. Default: internal.
- `SAWDBH`(4) — sawtimber DBH threshold. Default: internal.
- `CCMIN`(5) — minimum canopy cover. Default: internal.
- `TPAMIN`(6) — minimum trees/acre. Default: internal.
- `PCTSMX`(7) — max small-tree percent. Default: internal.
**Example:**
```text
STRCLASS           1
```
```yaml
- output:
    - STRCLASS: { print: 1 }
```
*Turns on structural-stage classification (default thresholds), enabling the BSCLASS/BCANCOV event-monitor variables.*

## Event monitor

The **event monitor** (evmon.f) is FVS's per-cycle expression evaluator. Each cycle it can read *stand variables* — quantities describing the current stand — and evaluate arithmetic/logical expressions over them, so management can be made *conditional* on the stand's state rather than fixed to calendar years. Example stand variables include `BBA` (before-thin basal area per acre), `BSDI` (before-cut Reineke stand density index), `TPA` (trees per acre), `AGE` (stand age = inventory age + elapsed years), `CYCLE`, `YEAR`, the structural-stage variables from `STRCLASS` (`BSCLASS`/`BSTRDBH`/`BCANCOV` and their `A…` after-thin counterparts), and the constants `YES`/`NO`/`ALL`. Expressions support `+ - * /`, `**` exponentiation, the comparators `GT GE LT LE EQ NE`, and functions `FRAC INT MOD EXP SQRT ALOG ABS RANN`. Two constructs drive it: `IF…THEN…{activities}…ENDIF` schedules activities only in cycles where a condition is true, and `COMPUTE…END` defines named user variables (evaluated def-before-use) that later conditions can read.

#### IF
**What it does:** `IF` begins an event-monitor conditional. The lines after `IF`, up to a line reading `THEN`, form the condition expression (over the stand variables above); the keywords after `THEN`, up to `ENDIF`, are the activities to schedule if the condition is true. Each cycle the monitor re-evaluates the condition and, when it holds, fires the enclosed activities (e.g. a thinning) with their year filled in at that time. This is how a prescription like "thin to a target basal area whenever the stand exceeds a density threshold" is expressed. The enclosed activities carry no date of their own (field 1 blank) — the trigger cycle supplies the year — and unrecognized keywords inside the block (such as a nested `COMPUTE`) are currently skipped.
**Parameters:**
- condition (raw lines before `THEN`) — an event-monitor boolean expression, e.g. `BBA GT 150`.
- activities (keyword cards between `THEN` and `ENDIF`) — one or more schedulable activities; supported inside `IF` are the thinning keywords (`THINBTA` `THINATA` `THINBBA` `THINABA` `THINPRSC` `THINDBH`, …) and `SPECPREF`. Their field 1 (date) is blank/filled at fire time; fields 2-7 are the activity's own parameters.
- `ENDIF` — closes the block.
**Example:**
```text
IF
BBA GT 150
THEN
THINBBA                      100.       0.       0.       0.       999.
ENDIF
```
```yaml
- event_monitor:
    - IF: {}
    - raw: "BBA GT 150"
    - THEN: {}
    - THINBBA: { year: 0, residual_basal_area: 100, cut_efficiency: 0, dbh_min: 0, dbh_max: 999 }
    - ENDIF: {}
```
*Each cycle, if before-thin basal area exceeds 150 ft²/ac, thins from below to a residual 100 ft²/ac.*

#### COMPUTE
**What it does:** `COMPUTE` opens an event-monitor **user-variable** block (initre.f option 33 → EVUSRV), read as `NAME = expression` assignments until a closing `END`. Field 1 is the start year/date (`IDT`, default 1): the block fires at that date only (`IDT=0` = every cycle; a cycle number `0<IDT<1000` = that 1-based cycle; else the calendar-year cycle that contains it), *not* automatically every cycle — so `COMPUTE MYCYC = CYCLE` with the default date freezes `MYCYC` at 1. Assignments are evaluated in order so a later definition may reference an earlier one (def-before-use), and each defined name persists in the event-monitor variable table for use in subsequent `IF` conditions (and, with `COMPUTDB`, is written to FVS_Compute). This lets you build derived quantities — ratios, indices, flags — from the built-in stand variables.
**Parameters:**
- `date`(1) — start year/cycle `IDT` for the block. Default: 1 (cycle 1). Codes: `0` = all cycles; `1…999` = that cycle number; `≥1000` = a calendar year.
- block body — `NAME = expression` lines (evaluated in order, def-before-use), terminated by `END`.
**Example:**
```text
COMPUTE            0
DENSE = BBA / BSDI
HEAVY = DENSE GT 1.2
END
IF
HEAVY GT 0
THEN
THINBBA                      120.       0.       0.       0.       999.
ENDIF
```
```yaml
- event_monitor:
    - COMPUTE: { date: 0 }
    - raw: "DENSE = BBA / BSDI"
    - raw: "HEAVY = DENSE GT 1.2"
    - END: {}
    - IF: {}
    - raw: "HEAVY GT 0"
    - THEN: {}
    - THINBBA: { year: 0, residual_basal_area: 120, cut_efficiency: 0, dbh_min: 0, dbh_max: 999 }
    - ENDIF: {}
```
*Defines `DENSE` (a BA-to-SDI ratio) and the flag `HEAVY` every cycle, then thins to 120 ft²/ac basal area whenever `HEAVY` is set.*
## Fire & fuels (FFE)

The **Fire and Fuels Extension (FFE)** bolts a fire/fuels/carbon model onto the growth simulation. All of its keywords live **inside an `FMIN … END` block** — `FMIN` opens the extension (and activates it for the stand), the keywords in between configure it, and `END` closes the block (exactly like the `ESTAB … END` regeneration packet). Once active, the FFE tracks, every cycle:

- **Surface fuel loadings** — dead woody debris on the ground, split into **size classes** by piece diameter: **1-hour** (`< ¼″`, the fast-drying fine fuels), **10-hour** (`¼–1″`), **100-hour** (`1–3″`), and **1000-hour** (`> 3″`, further binned 3–6″, 6–12″, 12–20″, 20–35″, 35–50″, `> 50″`), plus the **litter** (freshly-fallen needles/leaves) and **duff** (the decomposed organic layer beneath the litter) pools. Each pool exists in two states: **hard** (sound wood) and **soft** (rotten/decayed wood), which burn and decay differently.
- **Snags** — standing dead trees. The model ages each snag cohort, letting it lose height (**breakage**), decay from hard to soft, and eventually **fall** over into the down-wood pools, where it feeds the surface fuels.
- **Fuel dynamics** — fuels are added (tree mortality, crown-lift, litterfall), decayed (a fraction of each decayed pool passing to duff, the rest respired to the air), moved between pools, and treated (piling, crushing) year over year.
- **Fire behavior & effects** — a **scheduled fire** (`SIMFIRE`) or the hypothetical **potential fire** (`POTFIRE`) is run through a **fuel model** (a standard Rothermel fire-spread model that maps the fuel bed to spread rate and intensity) at a given **fuel moisture**, wind, and temperature, yielding a **flame length**, a **scorch height** (how high the crown is heated), fuel **consumption**, and tree **mortality** (larger/thicker-barked trees survive; small and thin-barked trees are killed and become snags).
- **Carbon** — the FFE rolls the live biomass, standing/down dead wood, forest floor, and (optionally) removed/released carbon into a per-cycle **Stand Carbon Report**.

**Variant note (SN).** Several FFE keywords are **recognized but inert in the Southern (SN) variant**: `DROUGHT` (drought-year fuel effects apply only to UT/CR/LS), `CANCALC` (canopy base-height / bulk-density for the crown-fire model — SN runs no crown fire), and `SOILHEAT` plus the family of report-only keywords (`BURNREPT`, `FUELOUT`, `SNAGSUM`, `MORTREPT`, `FUELREPT`, `DWDVLOUT`, …) whose text tables FVSjl does not print (the equivalent data is emitted via the DBS tables). `FIRECALC` method 0 (the SN default) is the only ported fire-behavior path. Any FFE keyword *not* handled warns at parse time, because ignoring a model keyword would silently change results.

### Fire event & behavior

#### FMIN
**What it does:** Opens the Fire and Fuels Extension keyword block and activates the FFE for the stand. Everything up to the matching `END` is read as FFE configuration (the fire event, fuel loadings, snag parameters, and report requests). With no interior keywords the FFE still runs its fuel/snag dynamics from the default loadings; it just has no scheduled fire. This is the FFE analogue of the `ESTAB` and `ECON` block openers.
**Parameters:** *(none — a bare block opener; `END` closes the block)*
**Example:**
```text
FMIN
SIMFIRE         2000     10.00         1      50.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SIMFIRE: { date: 2000, wind: 10, moisture_code: 1, temp: 50 }
    - END: {}
```
*Opens the FFE, schedules one fire in 2000, then closes the block.*

#### SIMFIRE
**What it does:** Schedules a simulated fire and specifies the weather it burns under. In the fire's cycle the FFE runs the fuel bed through the fuel model at the given wind/moisture/temperature to compute flame length, scorch height, fuel consumption, and fire-caused tree mortality (killed trees become snags). Each `SIMFIRE` card is its own event, so repeating the keyword schedules several fires. A bare `SIMFIRE` (no fields) fires in cycle 1 (the inventory year), not "never".
**Parameters:**
- `date`(1) — calendar year, or a cycle number if `≤ MAXCYC`; blank ⇒ cycle 1.
- `wind`(2) — 20-ft wind speed, mi/hr. Default `20`.
- `moisture_code`(3) — fuel-dryness scenario, `1`=very dry … `4`=very moist. Default `1`. A code outside 1–4 (e.g. `9`) is a no-op that leaves moisture at the last (moderate, model-3) value.
- `temp`(4) — air temperature, °F (drives scorch). Default `70`.
- `mort_code`(5) — `0` = no fire mortality, `1` = compute mortality. Default `1`.
- `pct_burn`(6) — percent of the stand area that burns, %. Default `100`.
- `season`(7) — burn season, `1`=spring/dormant … `4`. Default `1`.
**Example:**
```text
FMIN
SIMFIRE         2015     15.00         2      75.0         1      60.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SIMFIRE: { date: 2015, wind: 15, moisture_code: 2, temp: 75, mort_code: 1, pct_burn: 60 }
    - END: {}
```
*Burns 60 % of the stand in 2015 under a 15 mph wind, moderately-dry fuels, 75 °F, with mortality computed.*

#### FLAMEADJ
**What it does:** Adjusts the fire's computed flame length and the fraction of crown material consumed, letting the user tune fire intensity without changing the fuel bed. The flame multiplier scales the modeled flame length (and hence scorch and mortality); the crown-burn fraction sets what proportion of the tree crowns is consumed by the fire.
**Parameters:**
- `date`(1) — date/cycle the adjustment applies to.
- `flame_mult`(2) — multiplier on the computed flame length (`FLMULT`). Default `1.0` (no change).
- `crown_burn`(4) — percent of crown burned; stored as a fraction (`×0.01`). A value `≤ −1` means 0. Default: model-computed.
**Example:**
```text
FMIN
FLAMEADJ        2015                   1.50                 25.0
SIMFIRE         2015     15.00         2      75.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FLAMEADJ: { date: 2015, flame_mult: 1.5, crown_burn: 25 }
    - SIMFIRE: { date: 2015, wind: 15, moisture_code: 2, temp: 75 }
    - END: {}
```
*Amplifies the 2015 fire's flame length by 1.5× and burns 25 % of the crowns.*

#### FIRECALC
**What it does:** Selects the fire-behavior calculation logic. The SN default (`0`, "old FM logic") maps the fuel bed to weighted standard fuel models — this is the only path FVSjl implements, and it faithfully matches live FVS for SN. Methods `1` (new FM logic) and `2` (modelled loads → custom model 89), along with their surface-area-to-volume / bulk-density / heat-content overrides, are alternative behavior models that are **not ported**; requesting them warns and falls back to method 0.
**Parameters:**
- `date`(1) — date/cycle.
- `method`(2) — `0` = old-FM-logic (default, only ported path); `1`/`2` warn and fall back.
- *(fields 4–9)* — `USAV`/`UBD`/heat-content overrides for methods 1/2 — not ported.
**Example:**
```text
FMIN
FIRECALC                     0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FIRECALC: { method: 0 }
    - END: {}
```
*Explicitly selects the (default) old-FM-logic fire-behavior path — a no-op in SN.*

#### MOISTURE
**What it does:** Overrides the seven fuel-moisture values used by a scheduled fire on a given date, replacing the moisture-code table lookup with explicit percentages. Moisture is the water content (as % of dry weight) of each fuel class; drier fuels spread fire faster and burn more completely. The seven values cover the dead 1-/10-/100-hr and `3″+` woody classes, the duff, and the live woody and live herbaceous fuels.
**Parameters:**
- `date`(1) — date/cycle the override applies to; blank ⇒ cycle 1.
- `1hr`(2), `10hr`(3), `100hr`(4), `3in+`(5), `duff`(6), `live_woody`(7), `live_herb`(8) — fuel moisture, %. Blank ⇒ 0, **except** a blank `live_herb`(8) defaults to the `live_woody`(7) value.
**Example:**
```text
FMIN
MOISTURE        2015       4.0       5.0       8.0      12.0      40.0      90.0
SIMFIRE         2015     15.00
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - MOISTURE: { date: 2015, "1hr": 4, "10hr": 5, "100hr": 8, "3in_plus": 12, duff: 40, live_woody: 90 }
    - SIMFIRE: { date: 2015, wind: 15 }
    - END: {}
```
*Sets explicit dry-fuel moistures for the 2015 fire (live herb defaults to the 90 % live-woody value).*

#### DROUGHT  *(no-op in SN)*
**What it does:** Marks a span of drought years that (in the UT/CR/LS variants) dry the fuel model and raise fire intensity. **Recognized but a no-op in SN** — the Southern FFE ("OZ-FFE") does not apply the drought adjustment (`fmvinit.f`: "not used in OZ-FFE"). Included so an SN key using it parses cleanly.
**Parameters:**
- `year_begin`(1), `year_end`(2) — first/last drought year (`IDRYB`/`IDRYE`). *Inert in SN.*
**Example:**
```text
FMIN
DROUGHT         2005      2007
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - DROUGHT: { year_begin: 2005, year_end: 2007 }
    - END: {}
```
*Declares a 2005–2007 drought — recognized but has no effect in the SN variant.*

#### CANCALC  *(no-op in SN)*
**What it does:** Configures the canopy base-height and canopy bulk-density options that feed the **crown-fire** model (whether fire climbs from the surface into the canopy). **Recognized but a no-op in SN** — the Southern variant runs no crown-fire model, so these settings have nothing to act on.
**Parameters:** *(canopy base-height / bulk-density option codes)* — *inert in SN.*
**Example:**
```text
FMIN
CANCALC            1
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - CANCALC: { method: 1 }
    - END: {}
```
*Selects a canopy-fire option — recognized but inert in SN (no crown fire).*

### Fuel loadings & dynamics

#### FUELINIT
**What it does:** Sets the stand's initial **hard** (sound) surface-fuel loadings, in tons/acre, by piece-size class — overriding the FFE's modeled starting fuel bed. Use it to initialize a stand to a measured or known fuel load. Any field left blank (or `−1`) keeps the model default for that size class; the lumped `< 1″` field is split between the two finest classes unless they are given explicitly.
**Parameters:**
- `lt1in`(1) — `< 1″` fuels lumped (split into `<¼″` + `¼–1″`), tons/ac.
- `1to3in`(2), `3to6in`(3), `6to12in`(4), `12to20in`(5) — the coarser woody classes.
- `litter`(6), `duff`(7) — forest-floor litter and duff.
- `lt025in`(8), `025to1in`(9) — the `<¼″` (1-hr) and `¼–1″` (10-hr) fine classes explicitly.
- `20to35in`(10), `35to50in`(11), `gt50in`(12) — the largest woody classes.
- All: tons/ac; blank/`−1` ⇒ keep default.
**Example:**
```text
FMIN
FUELINIT                 2.0       1.5       1.0                 3.5       8.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FUELINIT: { "1to3in": 2.0, "3to6in": 1.5, "6to12in": 1.0, litter: 3.5, duff: 8.0 }
    - END: {}
```
*Initializes the hard woody, litter (3.5 t/ac) and duff (8 t/ac) fuel pools to measured values.*

#### FUELSOFT
**What it does:** Sets the stand's initial **soft** (rotten/decayed) surface-fuel loadings, in tons/acre, for the nine woody size classes directly (classes 1–9). Soft fuels are the decayed counterpart of the hard pools — they hold more moisture, decay faster, and consume differently in a fire. Blank/`−1` keeps the model default.
**Parameters:**
- `s1`(1)…`s9`(9) — soft-fuel loading for woody size classes 1 through 9 (finest → coarsest), tons/ac; blank/`−1` ⇒ keep default.
**Example:**
```text
FMIN
FUELSOFT       0.30      0.50      0.80      0.60
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FUELSOFT: { s1: 0.30, s2: 0.50, s3: 0.80, s4: 0.60 }
    - END: {}
```
*Seeds the four finest soft (rotten) woody-fuel classes with measured loadings.*

#### FUELMODL
**What it does:** Forces the fire to use one or more specified **standard fire-behavior fuel models** (instead of the FFE auto-selecting them from cover type and fuel load). A fuel model is a parameterized fuel bed (loadings, depth, surface-area ratios, moisture of extinction) that the Rothermel equations turn into spread rate and intensity. Up to three models can be blended with weights (normalized to sum 1) to represent a mixed fuel bed.
**Parameters:**
- `date`(1) — date/cycle; blank ⇒ cycle 1.
- `(model#, weight)` pairs — fields 2–7 hold up to **3** pairs. `model#` is a standard model 1–53; a blank/`≤0` weight defaults to `1`. Weights are normalized to sum 1. No valid model ⇒ auto-selection.
**Example:**
```text
FMIN
FUELMODL        2010         9      0.70        10      0.30
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FUELMODL: { date: 2010, models: [ {model: 9, weight: 0.70}, {model: 10, weight: 0.30} ] }
    - END: {}
```
*From 2010, forces a 70/30 blend of standard fuel models 9 and 10.*

#### DEFULMOD
**What it does:** Defines a new custom fuel model or alters an existing standard one, so the user can supply a site-specific fuel bed. It starts from a chosen standard model and overrides the fuel loadings, surface-area-to-volume ratios, bed depth, and moisture of extinction. Because a model has many parameters, this keyword reads a **supplemental record** (the next line, seven 10-column fields) for the remaining values. Any parameter left `−1`/blank keeps the base model's value.
**Parameters:**
- `date`(1) — date/cycle.
- `model#`(2) — the fuel model number to define/alter.
- fields 3–7 — dead 1-/10-/100-hr SAV, live SAV, dead 1-hr load (`PRMS(2–6)`).
- **supplemental line** (7×`F10`) — dead 10-/100-hr load, live-woody load, bed depth, moisture of extinction, live-herb SAV, live-herb load (`PRMS(7–13)`).
- All overrides: `−1`/blank ⇒ keep the standard model's value.
**Example:**
```text
FMIN
DEFULMOD        2010        14      1800      1600      1500      1500      2.5
        3.0       0.5       1.0       0.20     30.0      1600      0.30
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - DEFULMOD: { date: 2010, model: 14, dead_1hr_sav: 1800, dead_10hr_sav: 1600, dead_100hr_sav: 1500, live_sav: 1500, dead_1hr_load: 2.5, supplemental: [3.0, 0.5, 1.0, 0.20, 30.0, 1600, 0.30] }
    - END: {}
```
*Redefines fuel model 14 with custom loadings, a 1.0-ft depth and 0.20 moisture of extinction.*

#### FUELDCAY
**What it does:** Sets the annual **decay rate** of specific fuel size classes for one fuel decay class (the species-group decay speed, 1 = slow-rotting … 4 = fast). Decay is the fraction of a pool that breaks down each year; the FFE moves some of it to duff (see `DUFFPROD`) and respires the rest. Setting the decay-class ID to `≥ 5` copies decay-class 4's rates onto **all** classes. Woody/litter rates are capped at 1.0.
**Parameters:**
- `decay_class`(1) — decay-class ID 1–4 (`≥5` ⇒ apply class-4 rates to all).
- `litter`(2) — decay rate for the litter pool (size 10).
- `duff`(3) — decay rate for the duff pool (size 11).
- `s1`(4), `s2`(5), `s3`(6) — decay rates for woody size classes 1–3.
- `s4to9`(7) — decay rate applied to woody size classes 4–9.
**Example:**
```text
FMIN
FUELDCAY           2      0.50      0.02      0.30      0.15      0.08      0.04
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FUELDCAY: { decay_class: 2, litter: 0.50, duff: 0.02, s1: 0.30, s2: 0.15, s3: 0.08, s4to9: 0.04 }
    - END: {}
```
*Sets per-size-class annual decay rates for decay class 2.*

#### FUELMULT
**What it does:** Scales the total fuel decay rate of every size class by a per-decay-class multiplier — a quick way to speed up or slow down overall fuel breakdown without re-specifying each rate. One multiplier per decay class (1–4); a blank field leaves that class unchanged. The resulting rates are capped at 1.0.
**Parameters:**
- `class1`(1), `class2`(2), `class3`(3), `class4`(4) — decay-rate multiplier for decay classes 1–4; blank ⇒ leave as-is.
**Example:**
```text
FMIN
FUELMULT        1.20      1.20      1.50      1.50
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FUELMULT: { class1: 1.20, class2: 1.20, class3: 1.50, class4: 1.50 }
    - END: {}
```
*Speeds up decay by 20 % for the slow classes and 50 % for the fast classes.*

#### FUELPOOL
**What it does:** Assigns a species to a fuel **decay class**, controlling how fast that species' dead wood rots once it becomes surface fuel or a fallen snag. The four decay classes range from slow-decaying (durable species) to fast-decaying. This routes each species' contributed debris to the right decay-rate set (`FUELDCAY`/`FUELMULT`).
**Parameters:**
- `species`(1) — species (alpha/FIA code, `0`/blank = all, `−N` = SPGROUP *N*).
- `decay_class`(2) — decay class `1`–`4`.
**Example:**
```text
FMIN
FUELPOOL          LP         3
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FUELPOOL: { species: "LP", decay_class: 3 }
    - END: {}
```
*Puts loblolly-pine dead wood into decay class 3.*

#### FUELMOVE
**What it does:** Transfers fuel from one size pool to another on a given date — used to model mechanical treatments (chipping, crushing, redistribution) that reshuffle the fuel bed among size classes. The transfer can be an absolute amount, a proportion of the source pool, with a floor to leave behind and a ceiling target for the destination.
**Parameters:**
- `date`(1) — date/cycle.
- `from`(2) — source size pool, 0–11. Default `6`.
- `to`(3) — destination size pool, 0–11. Default `11`.
- `amount`(4) — absolute tons/ac to move. Default `0`.
- `proportion`(5) — fraction of the source pool to move. Default `0`.
- `leave`(6) — amount to leave in the source (`Z`). Default `9999`.
- `target_final`(7) — cap on the destination final amount (`Q`). Default `0`.
**Example:**
```text
FMIN
FUELMOVE        2012         3        11                0.50
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FUELMOVE: { date: 2012, from: 3, to: 11, proportion: 0.50 }
    - END: {}
```
*In 2012 moves half of size-pool 3 into the duff pool (11).*

#### FUELTRET
**What it does:** Adjusts the surface-fuel **bed depth** after a fuel treatment (piling, lop-and-scatter, crushing), which changes fire spread and intensity. The depth multiplier can be given directly or looked up from the built-in `DPMULT` table keyed by treatment type and harvest type. A shallower bed (multiplier `< 1`) slows fire; a fluffed bed (`> 1`) speeds it.
**Parameters:**
- `date`(1) — date/cycle.
- `treatment`(2) — treatment type `0`–`2`. Default `0`.
- `harvest`(3) — harvest type `1`–`3`. Default `1`.
- `depth_mult`(4) — depth multiplier; `−1` (or blank) ⇒ use the `DPMULT` table: treatment 0 → 1.0/1.3/1.6 by harvest type; treatment 1 → 0.83; treatment 2 → 0.75.
**Example:**
```text
FMIN
FUELTRET        2012         2         1
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FUELTRET: { date: 2012, treatment: 2, harvest: 1 }
    - END: {}
```
*In 2012 applies treatment type 2 (a 0.75× depth reduction from the table).*

#### DUFFPROD
**What it does:** Sets the proportion of each decayed fuel pool that becomes **duff** (versus being respired to the air) as it breaks down — the split between building up the forest-floor organic layer and carbon loss. Specified per size class of one decay class (ID `≥ 5` copies class-4 values to all). Values are clamped to [0, 1]; the default is 0.02 (2 % to duff).
**Parameters:**
- `decay_class`(1) — decay-class ID 1–4 (`≥5` ⇒ apply class-4 values to all).
- `litter`(2) — proportion for the litter pool (size 10).
- `s1`(3), `s2`(4), `s3`(5) — proportions for woody size classes 1–3.
- `s4to9`(6) — proportion for woody size classes 4–9.
- `all`(7) — set all sizes 1–10 at once.
- All: fraction in [0, 1]; default `0.02`.
**Example:**
```text
FMIN
DUFFPROD           1                                            0.10
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - DUFFPROD: { decay_class: 1, all: 0.10 }
    - END: {}
```
*Routes 10 % of all decayed woody material in decay class 1 into the duff pool.*

#### PILEBURN
**What it does:** Simulates a **jackpot or pile burn** — igniting concentrations of activity fuels (slash piles) rather than a broadcast surface fire. It consumes a share of the affected fuel and can kill a fraction of the trees near the piles. The parameters control how much of the stand is affected, how much fuel is treated/consumed, and the resulting tree mortality.
**Parameters:**
- `date`(1) — date/cycle.
- `type`(2) — pile/burn type. Default `1`.
- `affect`(3) — percent of stand area affected, %. Default `70`.
- `treat`(4) — percent of fuel treated/piled, %. Default `10`.
- `consumption`(5) — percent of piled fuel consumed (`FULCON`), %. Default `80`.
- `mortality`(6) — percent tree mortality near piles (`TRMORT`), %. Default `0`.
**Example:**
```text
FMIN
PILEBURN        2013         1      70.0      15.0      90.0       5.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - PILEBURN: { date: 2013, type: 1, affect: 70, treat: 15, consumption: 90, mortality: 5 }
    - END: {}
```
*In 2013 burns piles over 70 % of the stand, consuming 90 % of the piled fuel with 5 % tree mortality.*

### Snags & salvage

#### SNAGINIT
**What it does:** Adds user-specified **snags** (standing dead trees) to the stand at the start of the run, letting you initialize known standing dead wood that the tree list does not carry. Each card describes one snag cohort by species, size, height at death and current top height (a snag whose current height is below its death height has already broken), age since death, and density. These snags then age, decay, break, and fall through the normal FFE snag dynamics.
**Parameters:**
- `species`(1) — species (SPDECD code); `0`/`−999` ⇒ ignored.
- `dbh`(2) — DBH at death, inches.
- `ht_at_death`(3) — total height at death (`HTDEAD`; sets taper & dead volume), ft.
- `cur_ht`(4) — current top height (`HTIH`; the present, possibly-broken top), ft.
- `age`(5) — years since death.
- `density`(6) — snag density, stems/ac.
**Example:**
```text
FMIN
SNAGINIT          LP      14.0      60.0      50.0       8.0      12.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SNAGINIT: { species: "LP", dbh: 14.0, ht_at_death: 60.0, cur_ht: 50.0, age: 8.0, density: 12.0 }
    - END: {}
```
*Adds 12 stems/ac of 14″ loblolly-pine snags, dead 8 years and already broken from 60 ft to a 50-ft top.*

#### SNAGFALL
**What it does:** Sets per-species snag **fall-rate** parameters — how quickly standing dead trees topple to the ground and become down wood. `FALLX` scales the rate of falling; `ALLDWN` is the snag age by which the last 5 % have fallen (i.e. essentially all are down). Faster fall shortens the standing-dead phase and feeds the surface fuels sooner.
**Parameters:**
- `species`(1) — species (SPDECD; `0`/blank = all, `−N` = SPGROUP *N*).
- `fallx`(2) — rate-of-fall correction (`FALLX`), clamped `≥ 0.001`.
- `alldwn`(3) — snag age (yr) by which the last 5 % fall (`ALLDWN`), clamped `≥ 0`.
**Example:**
```text
FMIN
SNAGFALL          LP      1.50      25.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SNAGFALL: { species: "LP", fallx: 1.50, alldwn: 25.0 }
    - END: {}
```
*Makes loblolly-pine snags fall 1.5× faster, with all down by age 25.*

#### SNAGDCAY
**What it does:** Sets a per-species snag **decay-rate multiplier** (`DECAYX`), controlling how fast a hard snag decays into a soft (rotten) snag. Soft snags fall and break differently and count differently for fuel/carbon. Larger `DECAYX` decays snags faster. Defaults are roughly 0.07/0.21/0.35 across the decay classes. *(In SN snags do not lose height — `HTX = 0` — so the height-loss coupling is inert; `DECAYX` still drives the hard→soft transition, which is fully honored.)*
**Parameters:**
- `species`(1) — species (SPDECD; `0`/blank = all, `−N` = SPGROUP *N*).
- `decayx`(2) — snag decay-rate multiplier (`DECAYX`), clamped `≥ 0`.
**Example:**
```text
FMIN
SNAGDCAY          LP      0.25
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SNAGDCAY: { species: "LP", decayx: 0.25 }
    - END: {}
```
*Speeds loblolly-pine snag hard→soft decay to a 0.25 rate multiplier.*

#### SNAGBRK
**What it does:** Sets per-species snag **breakage** (height-loss) rates — how a standing dead tree loses its top over time, expressed as the years to lose 50 % of its height and the years to reach 30 % of its original height, separately for hard and soft snags. These are converted to the FFE's `HTX` height-loss coefficients. **By default (empty) `HTX = 0` in SN and NE — snags keep their full height** (no breakage); this keyword turns breakage on.
**Parameters:**
- `species`(1) — species (SPDECD).
- `yrs50_hard`(2), `yrs50_soft`(3) — years to lose 50 % of height, for hard / soft snags (`YRS50`).
- `yrs30_hard`(4), `yrs30_soft`(5) — years to reach 30 % of original height, for hard / soft snags (`YRS30`).
**Example:**
```text
FMIN
SNAGBRK           LP      20.0      10.0      40.0      25.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SNAGBRK: { species: "LP", yrs50_hard: 20, yrs50_soft: 10, yrs30_hard: 40, yrs30_soft: 25 }
    - END: {}
```
*Turns on loblolly-pine snag breakage: hard snags lose 50 % of height by 20 yr, soft ones by 10 yr.*

#### SNAGPSFT
**What it does:** Sets the per-species proportion of newly-created snags that start out **soft** (already partly decayed) at the moment of death, rather than hard. Fire-killed or disease-killed trees can enter the snag pool already rotten; this fraction seeds them directly into the soft state, changing subsequent fall/decay/fuel behavior.
**Parameters:**
- `species`(1) — species (SPDECD; `0`/blank = all, `−N` = SPGROUP *N*).
- `prop_soft`(2) — proportion soft at creation (`PSOFT`), clamped [0, 1].
**Example:**
```text
FMIN
SNAGPSFT          LP      0.20
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SNAGPSFT: { species: "LP", prop_soft: 0.20 }
    - END: {}
```
*Makes 20 % of new loblolly-pine snags start in the soft (decayed) state.*

#### SNAGPBN
**What it does:** Sets the **post-burn** snag-fall parameters — how fire-killed trees fall in the years after a fire, which differs from ordinary snag fall (fire-weakened boles fall faster, and small/scorched stems fall soonest). It controls the soft fraction, the small-stem fraction, the time window over which post-burn fall is accelerated, and size/scorch thresholds.
**Parameters:**
- `pb_soft`(1) — post-burn soft fraction (`PBSOFT`), clamped [0, 1].
- `pb_small`(2) — post-burn small-stem fraction (`PBSMAL`), clamped [0, 1].
- `pb_time`(3) — accelerated-fall window (`PBTIME`), years, min `1`.
- `pb_size`(4) — size threshold (`PBSIZE`), min `0`.
- `pb_scorch`(5) — scorch threshold (`PBSCOR`), min `0`.
**Example:**
```text
FMIN
SNAGPBN        0.50      0.80       5.0       6.0      50.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SNAGPBN: { pb_soft: 0.50, pb_small: 0.80, pb_time: 5, pb_size: 6, pb_scorch: 50 }
    - END: {}
```
*Accelerates post-fire snag fall over a 5-year window, dropping 80 % of small (< 6″) stems.*

#### SALVAGE (FFE)
**What it does:** Schedules a **salvage** operation that removes standing snags (dead-tree harvest) on a given date, within a DBH window and up to a maximum snag age, optionally restricted to hard or soft snags. A fraction is removed from the stand; a further proportion of the removed material can be **left on site as down wood** (adding to surface fuels) rather than hauled off. *(This is the functioning salvage keyword — the base-model `SALVAGE` under [Thinning](#thinning--harvest) is inert.)*
**Parameters:**
- `date`(1) — date/cycle; blank ⇒ cycle 1.
- `dbh_min`(2) — minimum snag DBH, in. Default `0`.
- `dbh_max`(3) — maximum snag DBH, in. Default `999`.
- `max_age`(4) — maximum snag age to salvage, yr. Default `5`.
- `ok_soft`(5) — `0` = all, `1` = hard only, `2` = soft only. Default `1`.
- `prop_removed`(6) — fraction of eligible snags removed (`PROP`). Default `0.9`.
- `prop_left`(7) — proportion of removed material left as down wood (`PROPLV`). Default `0`.
**Example:**
```text
FMIN
SALVAGE         2017      10.0     999.0       3.0         1      0.80      0.10
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SALVAGE: { date: 2017, dbh_min: 10, dbh_max: 999, max_age: 3, ok_soft: 1, prop_removed: 0.80, prop_left: 0.10 }
    - END: {}
```
*In 2017 salvages 80 % of hard snags ≥ 10″ and ≤ 3 yr old, leaving 10 % of the removed wood on site.*

#### SALVSP
**What it does:** Defines a species cut/leave list that scopes which species a co-scheduled FFE `SALVAGE` acts on — a cut-list (salvage only these species) or a leave-list (salvage all *except* these). This mirrors the `LEAVESP`/`SPECPREF` species filtering for live-tree thinning, but for the snag salvage.
**Parameters:**
- `date`(1) — date/cycle; blank ⇒ cycle 1.
- `species`(2) — species (SPDECD; `0`/`ALL`, index, or `−N` group).
- `flag`(3) — `< 1` ⇒ cut-list (salvage these), `≥ 1` ⇒ leave-list (spare these).
**Example:**
```text
FMIN
SALVSP          2017        LP         0
SALVAGE         2017
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SALVSP: { date: 2017, species: "LP", flag: 0 }
    - SALVAGE: { date: 2017 }
    - END: {}
```
*Restricts the 2017 salvage to loblolly-pine snags only (a cut-list).*

#### FMORTMLT
**What it does:** Multiplies the **fire-caused mortality** probability for a species within a DBH window — a knob to calibrate how many trees a fire kills (independently of ordinary background mortality). It is applied as `PMORT = PMORT · FMORTMLT` when the fire computes per-tree kill. **Note its field order differs from the growth `*MULT` keywords: the multiplier comes before the species.**
**Parameters:**
- `date`(1) — date/cycle.
- `multiplier`(2) — fire-mortality multiplier (applied to `PMORT`).
- `species`(3) — species (SPDECD; `0` = all, `−N` = group).
- `dbh_min`(4), `dbh_max`(5) — DBH window `[d1, d2)` the multiplier applies within. Defaults `0` / `99999`.
**Example:**
```text
FMIN
FMORTMLT        2015      1.50         0       0.0      12.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - FMORTMLT: { date: 2015, multiplier: 1.50, species: 0, dbh_min: 0, dbh_max: 12 }
    - END: {}
```
*From 2015, raises fire-caused mortality 1.5× for all species with DBH under 12″.*

### Reports & potential fire

#### CARBREPT
**What it does:** Requests the per-cycle **Stand Carbon Report**, which tallies the stand's carbon pools — live above/below-ground biomass, standing dead (snags), down dead wood, forest floor (litter/duff), and (with a fire) released and removed carbon. It is a flag: it only turns on report-row collection; the FFE fuel/snag dynamics run regardless.
**Parameters:** *(none — a flag)*
**Example:**
```text
FMIN
CARBREPT
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - CARBREPT: {}
    - END: {}
```
*Turns on the per-cycle Stand Carbon Report.*

#### CARBCALC
**What it does:** Selects the carbon-accounting **method** and reporting **units** for the Stand Carbon Report. The method chooses between the FFE fuel-based biomass and the Jenkins national biomass equations (the default, and the set FVSjl's live pools implement). The units field converts the reported values between US and metric conventions.
**Parameters:**
- `method`(1) — `0` = FFE fuel-based, `1` = Jenkins national biomass (default). Clamped to 0–1.
- `units`(2) — `0` = US tons/ac (default), `1` = metric t/ha, `2` = metric t/ac. Clamped to 0–2.
**Example:**
```text
FMIN
CARBREPT
CARBCALC           1         1
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - CARBREPT: {}
    - CARBCALC: { method: 1, units: 1 }
    - END: {}
```
*Reports carbon by the Jenkins method in metric tonnes/hectare.*

#### POTFIRE (POTFLAME)
**What it does:** Requests the **Potential Fire report**, which characterizes how a fire *would* behave in each cycle under specified severe and moderate weather **without actually burning the stand**. It reports the hypothetical flame length, scorch, mortality, and fuel consumption, so the user can gauge fire hazard as the stand and fuels evolve. `POTFLAME` is an accepted alias.
**Parameters:** *(none — a flag; the scenario weather is set by the `POTF*` keywords below)*
**Example:**
```text
FMIN
POTFIRE
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - POTFIRE: {}
    - END: {}
```
*Turns on the per-cycle Potential Fire hazard report.*

#### POTFMOIS
**What it does:** Sets the fuel **moisture** values for one Potential Fire weather scenario (severe or moderate), overriding the built-in scenario defaults. It supplies the seven moisture percentages (dead 1-/10-/100-hr, `3″+` woody, duff, live-woody, live-herb) used when computing that scenario's potential fire behavior.
**Parameters:**
- `scenario`(1) — `1` = severe, `2` = moderate. Default `1`.
- fields 2–8 — the 7 moisture % (1hr, 10hr, 100hr, `3in+`, duff, live-woody, live-herb) for that scenario. A blank field uses the scenario's `FMMOIS` default; a blank live-herb(8) ⇒ the (resolved) live-woody value.
**Example:**
```text
FMIN
POTFIRE
POTFMOIS           1       3.0       4.0       6.0      10.0      30.0      80.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - POTFIRE: {}
    - POTFMOIS: { scenario: 1, "1hr": 3, "10hr": 4, "100hr": 6, "3in_plus": 10, duff: 30, live_woody: 80 }
    - END: {}
```
*Sets dry fuel moistures for the severe potential-fire scenario.*

#### POTFWIND / POTFTEMP / POTFSEAS / POTFPAB
**What it does:** Each of these sets one weather/effect parameter for the two Potential Fire scenarios (a **severe** value in field 1 and a **moderate** value in field 2): `POTFWIND` the 20-ft wind speed (mi/hr, higher ⇒ faster spread and longer flames), `POTFTEMP` the air temperature (°F, feeds scorch height), `POTFSEAS` the burn-season code (shapes live-fuel moisture/mortality assumptions), and `POTFPAB` the percent of stand area burned (the fraction assumed to burn when the potential fire's effects are tallied). They only take effect when `POTFIRE` is on.
**Parameters (each keyword):**
- `severe`(1) — the value for the severe scenario.
- `moderate`(2) — the value for the moderate scenario.
**Example:**
```text
FMIN
POTFIRE
POTFWIND        20.0       6.0
POTFTEMP        90.0      70.0
POTFPAB        100.0      50.0
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - POTFIRE: {}
    - POTFWIND: { severe: 20, moderate: 6 }
    - POTFTEMP: { severe: 90, moderate: 70 }
    - POTFPAB: { severe: 100, moderate: 50 }
    - END: {}
```
*Sets 20/6 mph winds, 90/70 °F temperatures, and 100/50 % area-burned for the severe/moderate potential-fire scenarios.*

#### SOILHEAT  *(no-op)*
**What it does:** Requests the **soil-heating report** — how deeply a fire heats the soil (relevant to root and seed-bank survival). **Recognized but a no-op:** FVSjl does not print the soil-heating text report, so the keyword is accepted (so keys parse) but produces no output, like the other report-only FFE keywords.
**Parameters:** *(none — a report flag; not emitted)*
**Example:**
```text
FMIN
SOILHEAT
END
```
```yaml
- fire_and_fuels:
    - FMIN: {}
    - SOILHEAT: {}
    - END: {}
```
*Requests the soil-heating report — recognized but not emitted by FVSjl.*

> **Report-only FFE keywords.** `BURNREPT`, `FUELOUT`, `SNAGOUT`, `SNAGSUM`, `MORTREPT`, `FUELREPT`, `MOREOUT`, `LANDOUT`, `STATFUEL`, `MORTCLAS`, `SNAGCLAS`, `DWDVLOUT`, `DWDCVOUT`, `FMODLIST`, `FUELFOTO`, `CANFPROF`, `SVIMAGES`, `CARBCUT` (and `SOILHEAT`, above) are **recognized no-ops** — the text tables are not printed, but the equivalent data is emitted via the [DBS](#output--database) tables (`FVS_Carbon`, `FVS_Fuels`, `FVS_SnagSum`, `FVS_Down_Wood_Vol`/`_Cov`, `FVS_BurnReport`, `FVS_Consumption`, `FVS_Mortality`, `FVS_PotFire`). Any FFE keyword *not* documented above **warns at parse time**, because it is a model/override keyword whose omission would change results.
## Economics (ECON)

The ECON keywords all live **inside an `ECON … END` block** and drive FVS's discounted cash-flow analysis of a stand's management. The block produces two outputs: `FVS_EconSummary` (the discounted **present net value** and its components) and, for the log-graded revenue units, `FVS_EconHarvestValue` (a per-species, per-DIB-class detail table). A **discount rate** is the annual interest rate used to bring future dollars back to today: a cost or revenue of `amt` occurring `t` years out is worth `amt / (1+rate)^t` now. Present net value (PNV) is discounted revenue minus discounted cost over the horizon (costs booked at the start of each year, revenues at the end). **SEV (soil expectation value / Faustmann land value)** is the present value of an infinite series of identical rotations — `net·(1+rate)^R / ((1+rate)^R − 1)` for rotation length `R` — i.e. what the bare land is worth if managed this way forever.

#### ECON
**What it does:** Opens (and, with its matching `END`, closes) the economic-analysis block. Everything between `ECON` and `END` — `STRTECON`, `ANNUCST`, `HRVVRCST`, `HRVRVN` — is parsed into the stand's `EconState` and marks the analysis active. The block itself takes no fields; it is purely a delimiter that turns on discounting and collects the cost/revenue tables the discounting core (`econ_stand_pnv`) later values. Without it, no economic output is produced.
**Parameters:** *(none — block delimiter; must be terminated by an `END` card)*
**Example:**
```text
ECON
STRTECON       2000       4.0
ANNUCST        2.50
HRVVRCST      12.00         1       0.0      12.0
HRVVRCST      20.00         1      12.0     999.0
HRVRVN       180.00         4       6.0     ALL
HRVRVN       350.00         4      12.0     ALL
HRVRVN       525.00         4      18.0     ALL
END
```
```yaml
- economics:
    - ECON: {}
    - raw: "STRTECON       2000       4.0"
    - raw: "ANNUCST        2.50"
    - raw: "HRVVRCST      12.00         1       0.0      12.0"
    - raw: "HRVVRCST      20.00         1      12.0     999.0"
    - raw: "HRVRVN       180.00         4       6.0     ALL"
    - raw: "HRVRVN       350.00         4      12.0     ALL"
    - raw: "HRVRVN       525.00         4      18.0     ALL"
    - END: {}
```
*Opens an economic analysis discounting at 4%, with a $2.50/ac/yr fixed cost, per-tree cutting costs, and log-graded board-foot revenue in three DIB classes.*

#### STRTECON
**What it does:** Sets the analysis start and the discount rate for the ECON block. FVSjl honors the **discount rate** (field 2, entered as a **percent** — e.g. `4.0` = 4%/yr, stored internally as `0.04`); this is the rate applied to every cost and revenue in `econ_stand_pnv`/`econ_present_value`. The start-year/delay (field 1) and the two SEV fields (a known end-of-horizon soil expectation value, and a flag to have FVS compute SEV) are recognized but not modeled by the discounting core in FVSjl (the SEV/forest-value kernels exist — `econ_sev`, `econ_forest_value` — but are not wired to these fields). The default discount rate is 0 (undiscounted).
**Parameters:**
- `start`(1) — analysis start year or delay (years). Recognized, not modeled. Default: block start.
- `discount_rate`(2) — annual discount/interest rate, in **percent** (divided by 100 internally). Units: %/yr. Default: 0.
- `known_sev`(3) — a known end-of-horizon soil expectation value ($/ac). Recognized, not modeled.
- `compute_sev`(4) — flag requesting FVS compute SEV. Recognized, not modeled.
**Example:**
```text
STRTECON       2000       4.0
```
```yaml
- raw: "STRTECON       2000       4.0"
```
*Starts the economic analysis in 2000 and discounts all future cash flows at 4% per year.*

#### ANNUCST
**What it does:** Adds a fixed **annual management cost** ($/ac/yr) that accrues every year of the analysis horizon, discounted at the start of each year. Multiple `ANNUCST` cards accumulate (they are summed into a single per-year cost). This is the recurring overhead — administration, property tax, road maintenance — as opposed to the per-harvest variable costs from `HRVVRCST`. Default is 0 (no annual cost).
**Parameters:**
- `amount`(1) — annual cost per acre per year ($/ac/yr). Accumulates across cards. Default: 0.
**Example:**
```text
ANNUCST        2.50
```
```yaml
- raw: "ANNUCST        2.50"
```
*Charges $2.50 per acre every year of the analysis, discounted as an annual cost.*

#### HRVVRCST
**What it does:** Defines a **variable harvest cost by DBH class** — the per-unit cost of cutting trees whose DBH falls in `[dbh_lo, dbh_hi)`. The cost is charged as `amount × volume` for each removed tree, where the volume interpretation depends on the unit code (per tree, per thousand board feet, or per hundred cubic feet). These costs flow into each harvest's undiscounted cost stream (`econ_value_harvest`) and then into PNV. Multiple cards define successive DBH classes (e.g. cheaper per-tree cost for small stems, higher for large); every card applies to all species.
**Parameters:**
- `amount`(1) — cost per unit ($/unit). Required. Default: 0.
- `unit`(2) — unit code: **1** = per tree (TPA), **2** = per MBF (thousand board feet, `bdft·tpa/1000`), **3** = per CCF (hundred cubic feet, `cuft·tpa/100`). Default: 0.
- `dbh_lo`(3) — lower DBH bound of the class, inclusive (inches). Default: 0.
- `dbh_hi`(4) — upper DBH bound, exclusive (inches). Default: 999 (blank or ≤0 ⇒ 999, i.e. no upper limit).
**Example:**
```text
HRVVRCST      12.00         1       0.0      12.0
HRVVRCST      20.00         1      12.0     999.0
```
```yaml
- raw: "HRVVRCST      12.00         1       0.0      12.0"
- raw: "HRVVRCST      20.00         1      12.0     999.0"
```
*Charges $12 to cut each tree under 12" DBH and $20 to cut each tree 12" and larger.*

#### HRVRVN
**What it does:** Defines **harvest revenue by species and DBH/DIB class** — the per-unit price received for removed trees of a given species (or all species) at or above a diameter. For the per-tree/whole-tree units (1/2/3) the revenue is `amount × volume` and flows into PNV via `FVS_EconSummary`. For the **log-graded units (4/5)** each stem is bucked into logs, and each log is bucketed into the DIB class of its small-end inside-bark diameter and valued at that class's price; these populate the `FVS_EconHarvestValue` detail table and are **report-only** (they do not feed `FVS_EconSummary.Revenue`/PNV). Successive `HRVRVN` cards with the same unit and species define a ladder of DIB classes by their lower bounds (a class's upper bound is the next-larger class's lower bound, or 999.9 for the top class). A species-specific price overrides an `ALL` price for the same class.
**Parameters:**
- `amount`(1) — revenue (price) per unit ($/unit). Required. Default: 0.
- `unit`(2) — unit code: **1** = per tree (TPA); **2** = BF_1000, per thousand board feet of the **whole tree**; **3** = FT3_100, per hundred cubic feet of the **whole tree**; **4** = BF_1000_LOG, per thousand board feet **log-graded** (each stem bucked into logs, valued by log small-end DIB class — Scribner board feet, `FVS_EconHarvestValue`); **5** = FT3_100_LOG, per hundred cubic feet **log-graded** (per-log cubic bucking by DIB class). Units 4/5 are report-only. Default: 0.
- `dbh_lo`(3) — the class lower bound, inclusive: DBH for whole-tree units, or log-end **DIB** class bound for the log-graded units (inches). The upper bound is always taken from the next class (internally 999). Default: 0.
- `species`(4) — species the price applies to: an alpha/FIA/PLANTS code, or `ALL`/blank = every species. Default: ALL.
**Example:**
```text
HRVRVN       180.00         4       6.0     ALL
HRVRVN       350.00         4      12.0     ALL
HRVRVN       525.00         4      18.0     ALL
HRVRVN       600.00         4      18.0     SM
```
```yaml
- raw: "HRVRVN       180.00         4       6.0     ALL"
- raw: "HRVRVN       350.00         4      12.0     ALL"
- raw: "HRVRVN       525.00         4      18.0     ALL"
- raw: "HRVRVN       600.00         4      18.0     SM"
```
*Prices log-graded board feet (BF_1000_LOG) at $180/$350/$525 per MBF for the 6"/12"/18" DIB classes, with sugar maple in the top class overriding to $600/MBF.*

#### TCONDMLT
**What it does:** A **cut modifier** (not a thinning method itself) that adjusts each tree's cut-priority weight during a thinning by its **tree-condition class and special status**. FVS ranks trees for removal by a weight `WK2 = size + species-preference + TCWT·IMC + SPCLWT·ISPECL + point-density terms`; `TCONDMLT` supplies the `TCWT` (management-condition-code multiplier), `SPCLWT` (special-status multiplier), and the optional per-point density weights (basal area, CCF, TPA). A larger weight means the tree is removed earlier, so positive/negative multipliers bias cutting toward or away from trees in certain condition or status classes. All weights default to 0, which leaves only the base size + species-preference ranking (i.e. `TCONDMLT` is inert until set). The point-density terms are faithful for a single point and only matter for multi-point, point-density thinning.
**Parameters:**
- `date`(1) — schedule year, or cycle number if < 1000 (blank ⇒ cycle 1).
- `TCWT`(2) — weight multiplying the tree's management-condition code `IMC`. Default: 0.
- `SPCLWT`(3) — weight multiplying the tree's special-status flag `ISPECL`. Default: 0.
- `PBAWT`(4) — weight on the point's basal area (`PTBAA`). Default: 0.
- `PCCFWT`(5) — weight on the point's crown competition factor (`PCCF`). Default: 0.
- `PTPAWT`(6) — weight on the point's trees per acre (`PTPA`). Default: 0.
**Example:**
```text
TCONDMLT       2000       5.0      10.0
```
```yaml
- raw: "TCONDMLT       2000       5.0      10.0"
```
*Biases the year-2000 thinning to remove poorer-condition and special-status trees first, by adding 5×condition-code + 10×special-status to each tree's cut-priority weight.*

## Compression & tripling

FVS represents a stand as a list of tree **records**, each carrying a species, diameter, height, and a trees-per-acre (TPA) expansion factor. Two mechanisms manage that list. **Record tripling** splits each record into three — a central record plus a slightly larger and a slightly smaller one, with the TPA divided roughly 0.60/0.25/0.15 — so a single inventory record represents a small spread of diameters. This better captures within-record variance and, crucially, drives the stochastic mortality and diameter-growth serial-correlation machinery in the early cycles. **Compression** goes the other way: when the record list grows too large (from tripling and regeneration), similar records are clustered — via a principal-component score over their attributes — and merged into representative classes to keep the projection fast.

#### COMPRESS
**What it does:** Schedules **tree-record compression** at a given date: the record list is clustered down to a target number of classes to speed the projection, with a fraction of the classes found by simple attribute breaks (Method 1) and the remainder by principal-component splitting (Method 2, `comprs.f`). Similar records are merged into representative classes (species, DBH, height, and TPA-weighted), reducing record count with minimal loss of stand structure. In FVSjl the full IBM-EIGEN PCA clustering algorithm is ported and validated bit-exact against Fortran. Compression also suppresses tripling in the following cycle (the merged records are already the working set).
**Parameters:**
- `date`(1) — schedule year, or cycle number if < 1000. Default: 1.
- `target`(2) — target number of records/classes to compress to. Default: MAXTRE/2 = 1500.
- `pn1`(3) — percent of the target classes to form by Method 1 (attribute breaks); the rest by Method-2 PC splitting. Units: %. Default: 50.
**Example:**
```text
COMPRESS       2.0     200.0      60.0
```
```yaml
- keyword: "COMPRESS"
  params: ["2.0", "200.0", "60.0"]
```
*At cycle 2, compresses the tree-record list to 200 classes, 60% of them by Method-1 attribute breaks and the rest by principal-component splitting.*

#### NOTRIPLE
**What it does:** **Disables record tripling** entirely for the stand by setting the tripling-cycle count (`ICL4`) to 0. With tripling off, each inventory record is projected as a single record with no diameter spread, so the deterministic early-cycle variance expansion does not happen. This is used when the input already resolves diameter distribution finely, or to reduce record count and runtime, at the cost of the extra within-record variance tripling would provide. Takes no fields.
**Parameters:** *(none — sets the tripling count `ICL4` to 0)*
**Example:**
```text
NOTRIPLE
```
```yaml
- NOTRIPLE: {}
```
*Turns off record tripling so every tree record is projected as a single record with no diameter spread.*

#### NUMTRIP
**What it does:** Sets **how many early cycles use record tripling** (the `ICL4` count). By default the first 2 cycles triple each record (spreading diameters and driving stochastic mortality); after that, growth switches to the stochastic serial-correlation path with the un-tripled list. `NUMTRIP` overrides that count — e.g. `NUMTRIP 3` extends tripling to the first 3 cycles, `NUMTRIP 0` is equivalent to `NOTRIPLE`. Use it to tune how long the variance-expansion tripling stays active.
**Parameters:**
- `count`(1) — number of initial cycles to apply record tripling (`ICL4`). Default: 2.
**Example:**
```text
NUMTRIP           3
```
```yaml
- NUMTRIP: { count: 3 }
```
*Applies record tripling for the first 3 cycles instead of the default 2 before switching to the stochastic serial-correlation growth path.*
---

### Notes

- Every keyword's fields, units, defaults, and codes above are taken from the FVSjl keyword
  handlers in `src/engine/keyword_dispatch.jl` — each handler's `# OPTION n — NAME
  (initre.f:…)` comment is the authoritative field-by-field source — cross-checked against the
  FVS Fortran (`initre.f`, `cuts.f`, the `fm*.f` FFE, `volkey.f`/`sdefet.f`/`sdefln.f` volume,
  `esin.f`/`estab.f` establishment, `econ*` economics, `comprs.f`).
- The named-parameter schema (`src/io/yaml_keywords.jl` `_KW_SCHEMA`) and the section grouping
  (`_KW_SECTION` in the same file) are the machine-readable source for the YAML form; extend
  them to add named params / a section for more keywords. A keyword without a schema entry
  still round-trips via the positional `keyword:`/`params:` form, and an unrecognized section
  falls into `other`. Where an example uses the positional `keyword:`/`params:` form, that
  keyword has no (or partial) named schema yet — the positional form is exact and lossless.
- The hierarchical grouping is **order-preserving by construction**: the writer emits
  consecutive same-section records as one block and never reorders, so flattening the blocks
  reproduces the exact `.key` keyword sequence. The reader treats section names as labels only
  — list order is authoritative. Both directions are lossless w.r.t. the dispatch-relevant
  record (name, values, presence, field text); the original card's exact column padding is
  intentionally not preserved (it carries no meaning the handlers use).
- Keywords accepted but treated as `.sum`-inert no-ops in some variants (certain FFE cards in
  SN — `DROUGHT`/`CANCALC`/`SOILHEAT` and the report-only family; the base-model `SALVAGE`;
  certain calibration cards) still parse without error, and are noted inline where relevant.
- Runnable end-to-end examples live in [`../examples/`](../examples/); the input **structures**
  (both YAML flavors, the tree `.tre`/`.csv`) are in [FORMATS.md](FORMATS.md).
