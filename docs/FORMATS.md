# FVSjl input formats

FVSjl reads a stand from **keyword** input (what to run) plus **tree** input (the
inventory). Each comes in a legacy fixed-column Fortran form and a modern readable
form, and the engine consumes the parsed records, never the file — so the forms are
fully interchangeable.

| Role     | Legacy (stock FVS)        | Modern readable                                  |
|----------|---------------------------|--------------------------------------------------|
| Keywords | `.key` fixed-column cards | `.yaml` — **two** flavors (see below)            |
| Trees    | `.tre` fixed-column       | `.csv` named-header columns                      |

Convert any of these with the CLI (direction inferred from the extensions):

```
julia --project bin/fvsjl-translate.jl  <src>  <dst>  [tree-format]  [--flat]
#   .key  ↔ .yaml          .tre ↔ .csv
```

There are **two YAML keyword flavors**, distinguished by the top-level `format:` key:

1. **Keyword-stream YAML** (no `format:`, or `stand:`/`keywords:` as an ordered list) —
   an order-preserving image of the `.key` keyword stream, grouped into sections. It
   round-trips a `.key` **losslessly**. Documented in [KEYWORDS.md](KEYWORDS.md).
2. **Semantic stand YAML** (`format: fvs-stand/v1`) — a declarative description of a
   stand by *intent* (`invyr`, `numcycle`, `treatments`, `treelist`…). The converter
   **unravels** it into keyword records in the canonical FVS order. Documented here.

Both are read by the same entry points (`read_keyword_records`, `each_stand`,
`run_keyfile`); the discriminator picks the reader automatically.

## Contents

- [1. Semantic stand YAML (`fvs-stand/v1`)](#1-semantic-stand-yaml--format-fvs-standv1)
  - [Shape](#shape) · [Worked example](#worked-example) · [Emission order](#emission-order-how-the-stand-is-unraveled)
  - [Field → FVS-keyword reference](#field--fvs-keyword-reference) · [`raw_keywords:` escape hatch](#the-raw_keywords-escape-hatch) · [`treelist:` block](#the-treelist-block)
- [2. Tree input — `.tre` and `.csv`](#2-tree-input--tre-and-csv)
  - [Tree CSV schema](#tree-csv-schema) · [The TREEFMT FORMAT string](#the-treefmt-format-string)
- [3. Species codes (Southern variant)](#3-species-codes-southern-variant)
- [4. The `.sum` summary output](#4-the-sum-summary-output)
- [Running & converting — commands](#running--converting--commands)

For the **keywords** themselves (every keyword + every parameter), see
[KEYWORDS.md](KEYWORDS.md). Runnable examples are in [`../examples/`](../examples/).

---

## 1. Semantic stand YAML — `format: fvs-stand/v1`

The keyword-stream form still mirrors `.key` ordering. The **semantic** form does not:
a stand is a self-contained map of named settings, and the tool emits the keywords in
the right order for you. Keys mirror **FVS keyword/parameter names** (so an FVS user
maps them at a glance).

### Shape

```yaml
format: fvs-stand/v1     # REQUIRED discriminator (anything starting `fvs-stand`)
variant: SN              # which FVS model — SN = Southern (FVSsn), NE = Northeast. Default SN.

stand:                   # ONE stand …
  invyr: 1990
  ...

# — or —

stands:                  # … or MANY (each map is one self-contained stand)
  - { invyr: 1990, ... }
  - { invyr: 1995, ... }
```

A single `STOP` is appended after the last stand; each stand gets its own `PROCESS`.

**`variant:`** names the FVS geographic model to run the file as — `SN` (Southern, the
`FVSsn` binary) or `NE` (Northeast). In stock FVS the variant is *which binary you run*,
not a keyword; making the config carry it keeps a YAML file self-describing (the engine
reads it and dispatches to that variant). It is **not** written into a converted `.key`
(stock FVS has no variant keyword — you pick the `FVSsn`/`FVSne` binary). Resolution
order: an explicit `run_keyfile(...; variant=…)` argument wins; else the file's `variant:`;
else `SN`. The keyword-stream YAML accepts the same top-level `variant:` key.

### Worked example

```yaml
format: fvs-stand/v1
variant: SN
stand:
  stdident: "SN THINSDI — thin to residual SDI 200 in 2010"
  noautoes: true
  design:   { plots: 11, nsc: 1 }
  stdinfo:  { forest: 80106, habitat: "231Dd", age: 60, aspect: 315, slope: 30, elev: 7 }
  sitecode: { species: 63, index: 60.0 }
  invyr:    1990
  numcycle: 6
  echosum:  true
  treatments:
    - thinsdi: { year: 2010, sdi: 200 }
  treelist:
    format: "(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,...)"
```

unravels (`fvsjl-translate stand.yaml stand.key`) to exactly:

```
STDIDENT
SN THINSDI — thin to residual SDI 200 in 2010
NOAUTOES
DESIGN                                  11        1
STDINFO   80106     231Dd     60        315       30        7
SITECODE  63        60.0
INVYEAR   1990
NUMCYCLE  6
ECHOSUM
THINSDI   2010      200
TREEFMT
(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,...)
TREEDATA
PROCESS
STOP
```

The numeric field **text** is the natural value (`11`, not `11.0`); FVS parses both to
the same number, so the simulation is byte-identical. The semantic form is *semantic*,
not a byte-image of a particular `.key`.

### Emission order (how the stand is "unraveled")

Keys may appear in any order in the YAML (it is a map); the converter emits the
keywords in this fixed canonical order, which satisfies every hard FVS ordering
constraint (STDIDENT first; species groups before the treatments that reference them;
COMPUTE before IF; tree list and PROCESS last):

```
stdident → control flags → design → stdinfo → sitecode/setsite → invyr →
numcycle/timeint/cycleat → density (sdimax/bamax/sdicalc) → spgroup →
growth modifiers → output (echosum/summary) → raw_keywords → treatments →
regeneration → event_monitor → treelist (TREEFMT/TREEDATA) → PROCESS
```

`treatments` is emitted in the **author order of the list** — same-year cut activities
run in input order in FVS, so the list order is preserved verbatim.

### Field → FVS-keyword reference

**Scalars** — `key: value` → the keyword with `value` in the noted field:

| YAML key   | FVS keyword | field | YAML key    | FVS keyword | field |
|------------|-------------|-------|-------------|-------------|-------|
| `invyr`    | `INVYEAR`   | 1     | `managed`   | `MANAGED`   | 1     |
| `numcycle` | `NUMCYCLE`  | 1     | `numtrip`   | `NUMTRIP`   | 1     |
| `bamax`    | `BAMAX`     | 1     | `rannseed`  | `RANNSEED`  | 1     |
| `sdicalc`  | `SDICALC`   | 1     | `tfixarea`  | `TFIXAREA`  | 1     |

**Flags** — `key: true` emits the bare card; `false`/absent emits nothing:

| YAML key   | FVS keyword | YAML key   | FVS keyword |
|------------|-------------|------------|-------------|
| `noautoes` | `NOAUTOES`  | `notriple` | `NOTRIPLE`  |
| `screen`   | `SCREEN`    | `echosum`  | `ECHOSUM`   |
| `summary`  | `SUMMARY`   | `sprout`   | `SPROUT`    |

**Maps** — `key: { param: value, … }` → one card; params map to FVS field positions:

| YAML key   | FVS keyword | params (→ field #)                                                        |
|------------|-------------|---------------------------------------------------------------------------|
| `design`   | `DESIGN`    | `baf`1 `fixed_plot`2 `break_dbh`3 `plots`4 `nsc`5 `sample_weight`6 `stockable`7 |
| `stdinfo`  | `STDINFO`   | `forest`1 `habitat`2 `age`3 `aspect`4 `slope`5 `elev`6 `origin`9          |
| `sitecode` | `SITECODE`  | `species`1 `index`2                                                       |
| `estab`    | `ESTAB`     | scalar or `{ date }` → field 1                                            |

**`stdident:`** `"<text>"` → `STDIDENT` followed by the id text on its own line.

**List-valued keys** — each list entry is a card of the named keyword:

| YAML key       | entries are `{<tag>: {params}}` for tags …                                        |
|----------------|-----------------------------------------------------------------------------------|
| `treatments`   | `thinbba thinaba thinbta thinata thinsdi thincc thinrden thindbh thinht thinauto salvage` |
| `regeneration` | `plant natural`                                                                   |
| `growth`       | `baimult htgmult mortmult fixmort fixdg`                                           |
| `sdimax`       | `{ species, sdimax, pct_lo, pct_hi }` → `SDIMAX` (fields 1,2,5,6)                  |
| `setsite`      | `{ year, habitat, bamax, species, index, flag, sdimax }` → `SETSITE` (1–7)        |
| `timeint`      | `{ cycle, length }` → `TIMEINT` (1,2)                                              |
| `cycleat`      | `[year, …]` → one `CYCLEAT` per year                                               |
| `spgroup`      | `{ name, species: [codes…] }` → `SPGROUP <name>` + a member line                  |

**Treatment params** (FVS-term keys; `year` is field 1, blank ⇒ cycle 1):

| tag(s)                          | params (→ field #)                                              |
|---------------------------------|----------------------------------------------------------------|
| `thinbba` `thinaba`             | `ba`2 `eff`3 `dmin`4 `dmax`5 `species`6 `plot`7                 |
| `thinbta` `thinata`             | `tpa`2 `eff`3 `dmin`4 `dmax`5 `species`6 `plot`7               |
| `thinsdi`                       | `sdi`2 `eff`3 `dmin`4 `dmax`5 `species`6 `plot`7               |
| `thincc`                        | `ccf`2 `eff`3 `dmin`4 `dmax`5 `species`6 `plot`7               |
| `thinrden`                      | `rsdi`2 `eff`3 `dmin`4 `dmax`5 `species`6 `plot`7             |
| `thindbh`                       | `dmin`2 `dmax`3 `eff`4 `tpa`6 `species`7                        |
| `thinht`                        | `hmin`2 `hmax`3 `eff`4 `tpa`6 `species`7                        |
| `thinauto`                      | `eff`2                                                          |
| `salvage`                       | `dmin`2 `dmax`3 `eff`4 `species`5                               |
| `plant` `natural`               | `species`2 `tpa`3 `survival`4 `age`5 `height`6 `shade`7         |
| `baimult` `htgmult` `mortmult`  | `species`2 `mult`3 `dmin`4 `dmax`5                              |
| `fixmort`                       | `species`2 `rate`3 `dmin`4 `dmax`5 `option`6                    |
| `fixdg`                         | `species`2 `value`3 `dmin`4 `dmax`5                             |

### The `raw_keywords:` escape hatch

The schema above models the common keywords. **Anything else** — the full keyword
space (FFE/`FMIN`, `ECON`, the event monitor, any control card) — rides along in a
`raw_keywords:` list, written in the **keyword-stream entry form** (named `{KW: {…}}`,
positional `{keyword:, params:}`, or verbatim `{raw:}`):

```yaml
stand:
  invyr: 1990
  treatments:
    - thin_sdi: { year: 2010, residual_sdi: 200 }   # modeled
  raw_keywords:                                      # not (yet) modeled → verbatim
    - FMIN: {}
    - SIMFIRE: { year: 2005 }
    - { raw: "Econ" }
    - STRTECON: { ... }
    - { raw: "End" }
```

`raw_keywords` are emitted **after** the modeled non-treatment cards and **before** the
tree list / PROCESS — the right place for extension blocks (`ECON`, FFE/`FMIN`) and the
event monitor, which need only to precede `PROCESS`. (A keyword with a same-stand
ordering dependency on the modeled treatments, e.g. an `SPGROUP` a later `THIN` names,
should use the modeled `spgroup:` key or the keyword-stream YAML form instead.)

### The `treelist:` block

```yaml
treelist:
  format: "(I4,T1,I7,F6.0,...)"   # optional: → TREEFMT + the FORMAT line(s), for a .tre companion
# treelist: {}                    # just emit TREEDATA (inherit the format / use a .csv companion)
```

`format:` emits `TREEFMT` and the FORMAT string (auto-split if > 72 cols, since
`kw_treefmt!` reads two ≤80-col lines); `treelist:` then emits `TREEDATA`. The inventory
itself is **not** inline — it always comes from the companion file resolved by the
keyfile's base name (covered next). So `treelist.format` selects only the `.tre` *layout*,
never the *file*. (Tree records are never embedded in the YAML; there is no `data:`/`file:`
key — every stand reads the one companion, which is what makes a multi-scenario run work.)

### Multi-scenario runs (several stands, one inventory)

A `stands:` list describes **N scenarios of the same stand** — each entry is the same
inventory under a different treatment. They all read the **one** companion tree file
(`<keyfile>.csv`/`.tre`, by base name). Stock FVS reads that file sequentially, so on its
own only the first scenario would get trees; the converter therefore emits a **`REWIND 2`**
before every scenario after the first (this rewinds the tree-data unit so the next stand
re-reads it — exactly the pattern in FVS's own `snt01.key`). FVSjl re-reads the file for
every stand implicitly, so `REWIND` is a no-op there and the stream is correct for both —
**the converted `.key` reproduces live FVSsn** (verified: `examples/multiscenario/`, all
scenarios start from the identical 536-TPA inventory, FVSjl vs FVSsn within the ±1-cuft
single-precision tail). To give scenarios *different* inventories you'd use separate
keyfiles, each with its own companion.

---

## 2. Tree input — `.tre` and `.csv`

`TREEDATA` tells FVS to load the inventory. The records come from a file **resolved by
the keyword file's base name** (input path minus extension):

* `<base>.csv` — used if present (a **self-describing** named-header CSV; no `TREEFMT`
  needed), else
* `<base>.tre` — the fixed-column Fortran file, sliced by the current `TREEFMT`
  (defaults to the SN layout, `sn/blkdat.f`).

So `thinsdi.yaml` (or `.key`) loads `thinsdi.csv` if it exists, otherwise `thinsdi.tre`.
This is the only coupling between the keyword file and the tree file — **the tree form
is chosen independently of the keyword form.** Any combination works:

| keyword file   | tree companion   | notes                                           |
|----------------|------------------|-------------------------------------------------|
| `.key`         | `.tre`           | stock-FVS layout                                |
| `.yaml` (either)| `.csv`          | fully modern; `TREEFMT` is ignored for CSV      |
| `.yaml`        | `.tre`           | needs `treelist.format` / a `TREEFMT` keyword   |

### Tree CSV schema

A header row (required, must match exactly) then one row per tree. 25 columns, in this
order — each maps to a fixed-column field of the `.tre` record (`intree.f`):

| CSV column    | FVS var  | type  | meaning                          |
|---------------|----------|-------|----------------------------------|
| `plot`        | ITREI    | int   | point/plot number                |
| `id`          | IDTREE   | int   | tree id                          |
| `tpa`         | PROB     | float | trees per acre represented       |
| `history`     | ITH      | int   | history/status code              |
| `species`     | CSPI     | code  | species (FIA/PLANTS/alpha code)  |
| `dbh`         | DBH      | float | diameter at breast height        |
| `diam_growth` | DG       | float | measured diameter growth         |
| `height`      | HT       | float | total height                     |
| `top_height`  | THT      | float | height to a broken/dead top      |
| `ht_growth`   | HTG      | float | measured height growth           |
| `crown_pct`   | ICR      | int   | crown ratio (percent)            |
| `damage1..6`  | IDAMCD   | int×6 | 3 (agent, severity) damage pairs |
| `mort_code`   | IMC1     | int   | mortality code                   |
| `cut_code`    | KUTKOD   | int   | cut/removal code                 |
| `pest1..5`    | IPVARS   | int×5 | pest-extension variables         |
| `birth_age`   | ABIRTH   | float | age at birth, if supplied        |

The CSV is a **lossless** image of the parsed `TreeRecord` (numbers compact but exact;
the species code stored stripped). `.tre ↔ .csv` round-trips through
`read_tree_records` / `write_trees_csv` reproduce the same records, so the engine yields
identical results from either.

```
# trees.csv
plot,id,tpa,history,species,dbh,diam_growth,height,top_height,ht_growth,crown_pct,damage1,damage2,damage3,damage4,damage5,damage6,mort_code,cut_code,pest1,pest2,pest3,pest4,pest5,birth_age
1,1,5,1,131,9.1,0,46,0,0,40,0,0,0,0,0,0,0,0,0,0,0,0,0,0
```

The `.tre` FORMAT (the column layout) is **only** needed for fixed-column `.tre` files;
a `.csv` carries its own column meaning in the header, so `TREEFMT` / `treelist.format`
is irrelevant when a CSV companion is used. Override the `.tre` layout for translation
with the CLI's `[tree-format]` argument; the SN default is in `DEFAULT_TREE_FORMAT`
(`src/io/treedata.jl`).

### The TREEFMT FORMAT string

A `.tre` file is fixed-column: each tree is one line, and a tree's fields sit in
specific **column ranges**. A Fortran **FORMAT string** says where. It is the same
string FVS uses on a `TREEFMT` card (the `treelist.format` value in the semantic YAML),
and it is parsed by `src/io/treedata.jl`. The default SN layout is:

```
(I4,T1,I7,F6.0,I1,A3,F4.1,F3.1,2F3.0,F4.1,I1,3(I2,I2),2I1,I2,2I3,2I1,F3.0)
```

The descriptors used:

| descriptor | meaning |
|------------|---------|
| `Iw`       | integer in the next `w` columns |
| `Fw.d`     | real in the next `w` columns; if the text has no `.`, `d` implied decimals (e.g. `F4.1` reads `"  91"` as `9.1`) |
| `Aw`       | character string in the next `w` columns (the species code) |
| `Tn`       | **tab**: jump to absolute column `n` (no field read) |
| `n(...)`   | repeat the group `n` times — e.g. `3(I2,I2)` = the 6 damage fields |
| `nFw.d`    | repeat a descriptor `n` times — e.g. `2F3.0` = two `F3.0` floats |

The descriptors are consumed **left to right**, advancing a column cursor, and produce
the 25 fields in the [Tree CSV schema](#tree-csv-schema) order (`plot, id, tpa, …`). So
the default reads `plot` from cols 1–4 (`I4`), then `T1` jumps **back to column 1** and
`id` is read from cols 1–7 (`I7`) — the SN layout deliberately **overlaps** the 4-col
plot and the 7-col id (a small id sits inside the plot's columns; a wide id fills the
full 7). `tpa` is then `F6.0` from cols 8–13, and so on. To author a `.tre` by hand,
match each value to its descriptor's columns; or just use the `.csv` form, which needs
no FORMAT at all.

---

## 3. Species codes (Southern variant)

The `species` field (`.tre`/`.csv` column 5) accepts **any** of three codes — FVSjl
resolves them all to the same internal species. Keywords that take a `species` argument
accept the **FVS alpha** code (or `0`/`ALL`, or `−N` for an [SPGROUP]). The 90 Southern
species (from `data/southern/species_coefficients.csv`, the source of truth):

| # | FVS | common name | FIA | PLANTS |
|---|-----|-------------|-----|--------|
| 1 | FR | fir | 010 | ABIES |
| 2 | JU | eastern redcedar | 057 | JUNIP |
| 3 | PI | spruce | 090 | PICEA |
| 4 | PU | sand pine | 107 | PICL |
| 5 | SP | shortleaf pine | 110 | PIEC2 |
| 6 | SA | slash pine | 111 | PIEL |
| 7 | SR | spruce pine | 115 | PIGL2 |
| 8 | LL | longleaf pine | 121 | PIPA2 |
| 9 | TM | Table Mountain pine | 123 | PIPU5 |
| 10 | PP | pitch pine | 126 | PIRI |
| 11 | PD | pond pine | 128 | PISE |
| 12 | WP | eastern white pine | 129 | PIST |
| 13 | LP | loblolly pine | 131 | PITA |
| 14 | VP | Virginia pine | 132 | PIVI2 |
| 15 | BY | baldcypress | 221 | TADI2 |
| 16 | PC | pondcypress | 222 | TAAS |
| 17 | HM | eastern hemlock | 260 | TSUGA |
| 18 | FM | Florida maple | 311 | ACBA3 |
| 19 | BE | boxelder | 313 | ACNE2 |
| 20 | RM | red maple | 316 | ACRU |
| 21 | SV | silver maple | 317 | ACSA2 |
| 22 | SM | sugar maple | 318 | ACSA3 |
| 23 | BU | buckeye | 330 | AESCU |
| 24 | BB | birch | 370 | BETUL |
| 25 | SB | sweet birch | 372 | BELE |
| 26 | AH | American hornbeam | 391 | CACA18 |
| 27 | HI | hickory | 400 | CARYA |
| 28 | CA | catalpa | 450 | CATAL |
| 29 | HB | hackberry | 460 | CELTI |
| 30 | RD | eastern redbud | 471 | CECA4 |
| 31 | DW | flowering dogwood | 491 | COFL2 |
| 32 | PS | common persimmon | 521 | DIVI5 |
| 33 | AB | American beech | 531 | FAGR |
| 34 | AS | ash | 540 | FRAXI |
| 35 | WA | white ash | 541 | FRAM2 |
| 36 | BA | black ash | 543 | FRNI |
| 37 | GA | green ash | 544 | FRPE |
| 38 | HL | honeylocust | 552 | GLTR |
| 39 | LB | loblolly-bay | 555 | GOLA |
| 40 | HA | silverbell | 580 | HALES |
| 41 | HY | American holly | 591 | ILOP |
| 42 | BN | butternut | 601 | JUCI |
| 43 | WN | black walnut | 602 | JUNI |
| 44 | SU | sweetgum | 611 | LIST2 |
| 45 | YP | yellow-poplar | 621 | LITU |
| 46 | MG | magnolia | 650 | MAGNO |
| 47 | CT | cucumbertree | 651 | MAAC |
| 48 | MS | southern magnolia | 652 | MAGR4 |
| 49 | MV | sweetbay | 653 | MAVI2 |
| 50 | ML | bigleaf magnolia | 654 | MAMA2 |
| 51 | AP | apple | 660 | MALUS |
| 52 | MB | mulberry | 680 | MORUS |
| 53 | WT | water tupelo | 691 | NYAQ2 |
| 54 | BG | blackgum | 693 | NYSY |
| 55 | TS | swamp tupelo | 694 | NYBI |
| 56 | HH | eastern hophornbeam | 701 | OSVI |
| 57 | SD | sourwood | 711 | OXAR |
| 58 | RA | redbay | 721 | PEBO |
| 59 | SY | sycamore | 731 | PLOC |
| 60 | CW | cottonwood | 740 | POPUL |
| 61 | BT | bigtooth aspen | 743 | POGR4 |
| 62 | BC | black cherry | 762 | PRSE2 |
| 63 | WO | white oak | 802 | QUAL |
| 64 | SO | scarlet oak | 806 | QUCO2 |
| 65 | SK | southern red oak | 812 | QUFA |
| 66 | CB | cherrybark oak | 813 | QUPA5 |
| 67 | TO | turkey oak | 819 | QULA2 |
| 68 | LK | laurel oak | 820 | QULA3 |
| 69 | OV | overcup oak | 822 | QULY |
| 70 | BJ | blackjack oak | 824 | QUMA3 |
| 71 | SN | swamp chestnut oak | 825 | QUMI |
| 72 | CK | chinkapin oak | 826 | QUMU |
| 73 | WK | water oak | 827 | QUNI |
| 74 | CO | chestnut oak | 832 | QUPR2 |
| 75 | RO | northern red oak | 833 | QURU |
| 76 | QS | Shumard oak | 834 | QUSH |
| 77 | PO | post oak | 835 | QUST |
| 78 | BO | black oak | 837 | QUVE |
| 79 | LO | live oak | 838 | QUVI |
| 80 | BK | black locust | 901 | ROPS |
| 81 | WI | willow | 920 | SALIX |
| 82 | SS | sassafras | 931 | SAAL5 |
| 83 | BD | basswood | 950 | TILIA |
| 84 | EL | elm | 970 | ULMUS |
| 85 | WE | winged elm | 971 | ULAL |
| 86 | AE | American elm | 972 | ULAM |
| 87 | RL | slippery elm | 975 | ULRU |
| 88 | OS | other softwood | 299 | 2TN |
| 89 | OH | other hardwood | 998 | 2TB |
| 90 | OT | other tree | 999 | 2TREE |

The **FIA** column is the standard USDA-FS Forest Inventory & Analysis species code; the
**PLANTS** column is the USDA PLANTS symbol (the scientific-name identifier). The common
names are a convenience gloss — the codes are authoritative. (Common-name spelling is not
parsed; only the FVS/FIA/PLANTS codes are accepted in the `species` field.)

---

## 4. The `.sum` summary output

Running a stand (`run_keyfile`, or stock FVS) produces the **`.sum`** summary — one block
per stand (scenario), each a `-999` header line followed by **one row per cycle**. FVSjl
emits the identical SUMOUT/IOSUM layout as FVSsn, so the two are directly comparable; each
example folder ships both for reference:

| file | produced by |
|------|-------------|
| `<name>.fvsjl.sum` | **FVSjl**, run from the `.yaml` + `.csv` (`run_keyfile`) |
| `<name>.fvssn.sum` | **live FVSsn** (the Fortran binary), run from the equivalent `.key` + `.tre` |

They match to within the **±1-cuft single-precision tail** (volume columns may differ by 1
from Float32 rounding; the data otherwise agrees). The `-999` header timestamp differs per
run — compare the data rows.

**Header line** — `-999  <#periods>  <StandID>  <MgmtID>  <SampleWt>  <Variant>  <Date>  <Time>  …  <#plots>`
(`SampleWt` in Fortran `E15.7`, e.g. `0.1100000E+02` = 11.0).

**Data row** — one per cycle, 29 fields in four groups (a field is `0` when not applicable,
e.g. the removal columns outside a thinning cycle):

| # | field | meaning |
|---|-------|---------|
| | | **Identity + start-of-period stand** (live, *before* any thinning this cycle) |
| 1 | Year | calendar year of the cycle boundary |
| 2 | Age | stand age (yrs) |
| 3 | Tpa | live trees per acre |
| 4 | BA | basal area (ft²/ac) |
| 5 | SDI | stand density index |
| 6 | CCF | crown competition factor |
| 7 | TopHt | top height (ft) |
| 8 | QMD | quadratic mean diameter (in) |
| 9 | TCuFt | total cubic-foot volume / ac |
| 10 | MCuFt | merchantable cubic ft / ac |
| 11 | SCuFt | sawtimber cubic ft / ac |
| 12 | BdFt | board feet / ac (Scribner) |
| | | **Removed by treatment this period** (0 if no cut) |
| 13 | RTpa | removed trees / acre |
| 14 | RTCuFt | removed total cubic ft |
| 15 | RMCuFt | removed merch cubic ft |
| 16 | RSCuFt | removed sawtimber cubic ft |
| 17 | RBdFt | removed board feet |
| | | **After-treatment stand** (= the start columns when no cut occurred) |
| 18 | ATBA | basal area after treatment |
| 19 | ATSDI | SDI after treatment |
| 20 | ATCCF | CCF after treatment |
| 21 | ATTopHt | top height after treatment |
| 22 | ATQMD | QMD after treatment |
| | | **Growth (this period) + classification** |
| 23 | PrdLen | period length (yrs) |
| 24 | Accret | periodic annual accretion (cuft/ac/yr) |
| 25 | Mort | periodic annual mortality (cuft/ac/yr) |
| 26 | MAI | mean annual increment (cuft/ac/yr) |
| 27 | ForType | forest-type code |
| 28 | SizeCls | size class |
| 29 | StkCls | stocking class |

For a multi-scenario file the blocks appear in scenario order; e.g.
`examples/multiscenario/stand.fvsjl.sum` has four blocks, all starting from the identical
inventory (row 1 of each), then diverging by treatment.

---

## Running & converting — commands

The engine reads a `.key` **or** a `.yaml` (either flavor) directly; the tree companion
(`.tre`/`.csv`) is found by base name. A runnable script that exercises all of the below
is [`../examples/convert_and_run.sh`](../examples/convert_and_run.sh).

**Run a stand** (prints the `.sum`):

```bash
# legacy .key (+ thinba.tre)             |  semantic YAML (+ thinsdi.csv)
julia --project=. -e 'using FVSjl; print(run_keyfile("examples/thinba/thinba.key"))'
julia --project=. -e 'using FVSjl; print(run_keyfile("examples/semantic/thinsdi.yaml"))'
# save it instead of printing:
julia --project=. -e 'using FVSjl; write("out.sum", run_keyfile("examples/semantic/multistand.yaml"))'
```

**Convert** (direction inferred from the extensions; `bin/fvsjl-translate.jl <src> <dst>`):

```bash
JL='julia --project=.'

# keywords: .key  ->  .yaml   (and back).  Works for BOTH yaml flavors on the way back.
$JL bin/fvsjl-translate.jl examples/thinba/thinba.key      thinba.yaml      # keyword-stream YAML
$JL bin/fvsjl-translate.jl examples/thinba/thinba.key      thinba.yaml --flat   # legacy flat list
$JL bin/fvsjl-translate.jl thinba.yaml                     thinba.key       # YAML -> .key (stock FVS)
$JL bin/fvsjl-translate.jl examples/semantic/thinsdi.yaml  thinsdi.key      # SEMANTIC YAML -> .key

# trees: .tre  ->  .csv   (and back).  Pass the FORMAT for a non-default .tre layout:
FMT='(T24,I4,T1,I4,T31,F2.0,I1,A3,F3.1,F2.1,T45,F3.0,T63,F3.0,T60,F3.1,T48,I1,T52,I2,T66,5I1,T54,7I1,T75,F3.0)'
$JL bin/fvsjl-translate.jl examples/thinba/thinba.tre      thinba.csv "$FMT"
$JL bin/fvsjl-translate.jl thinba.csv                      thinba.tre "$FMT"
```

> A semantic-YAML → `.key` conversion **unravels** the declarative stand into keyword
> cards in the canonical order; `.key` → YAML emits the order-preserving keyword-stream
> form (there is no `.key` → *semantic* direction — the semantic form is for authoring).
