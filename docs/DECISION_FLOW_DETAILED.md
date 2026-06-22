# FVS-Southern decision flow — branch-level detail

Companion to [`DECISION_FLOW.md`](DECISION_FLOW.md) (which is the bird's-eye
input→output map). This document keeps the same spine but expands each hot-path
routine into its **decision branches**, with the gating condition and the FVSjl
port status of each branch. The purpose is branch-completeness: a routine can be
"ported" yet silently miss a branch the test data never exercises (this is how
the HABTYP-index, BAMAX, and size-cap gaps slipped in). Each branch listed here
is something to confirm, and each ⛔ is a scenario to add.

Legend: ✅ ported & validated · 🟡 partial / simplified · ⛔ not ported (in-scope)
· 🧊 out of scope by the plan / no scenario · ⚪ no-op in the oracle too (gated
off) — FVSjl is correct to omit it.

Fortran refs are `file.f`; FVSjl refs are `src/...`.

---

## INITRE — input (`initre.f`) → `initialize()`

| branch / condition | effect | status |
|---|---|---|
| keyword loop: dispatch each card by name | sets all run options | ✅ `engine/keyword_dispatch.jl` |
| STDINFO field 1 (forest code) | lat/long/elev defaults (FORKOD) | ✅ |
| STDINFO field 2 = habitat (`HABTYP`/`HBDECD`) | numeric ⇒ index into 320-entry SNECU table; alpha ⇒ uppercased exact match; else default #122 "231DD" | ✅ `variants/southern/habitat.jl` |
| SITECODE / site-species (`SITSET`) | per-species site index array | ✅ `site_index.jl` |
| DESIGN / sample design | TPA expansion factors | ✅ |
| INVYEAR / NUMCYCLE / TIMEINT | cycle calendar | ✅ |
| thinning/harvest keywords (THIN*/SALVAGE/SPECPREF/…) | schedule CUTS + set cut modifiers | 🟡 5 of ~17 methods + 0 of 6 modifiers ported — see the destructured **CUTS** section below for the per-keyword audit |
| MSB / SIZECAP / MORTMULT / FIXMORT / FIXDG / FIXHTG / HTGSTOP / TOPKILL / FFERT | option activities | ⛔ keyword paths not wired (defaults = no-op) |
| BAMAX (SETSITE basal-area max) keyword | sets LBAMAX + BAMAX | 🟡 BAMAX honored in MORTS; keyword path partial |

## NOTRE / SETUP — build records (`notre.f`, `setup.f`)

| branch | effect | status |
|---|---|---|
| expand tree records by sample design | live `TreeList` (SoA) | ✅ `notre!` |
| RNG seed (`GETSED`) | faithful LCG seed | ✅ `FVSRng` |
| DGCONS / site-dep DG constants | `calib.dg_const`, `atten` | ✅ `dgcons!` |
| SDICHK over-dense reset (`sdichk.f`) | if TPROB>(PMSDIU+0.05)·max ⇒ reset all SDIDEF | ✅ `sdi_max_check!` |
| LSTART DG calibration (`DGDRIV` cal pass) | per-species COR from input DG | ✅ `calibrate_diameter_growth!` |

---

## Per-cycle REPORT step (before growing)

| routine | branches | status |
|---|---|---|
| `CRATET` | dub HT=0 trees; resolve broken-top NORMHT; LSTART dub-crowns; ITRN==0 / TPROB≤0 early returns | ✅ `dub_missing_heights!` (init); per-cycle re-dub only matters with regen |
| `VOLS` | per-tree cuft (R8 Clark) / merch MCF (D≥DBHMIN) / sawtimber SCF / board-ft (Scribner); topkill CFTOPK/BFTOPK | ✅ `compute_volumes!` |
| `CWIDTH` | crown-width by species eq | ✅ |
| `STATS`/`DISPLY`/`SUMOUT` | stand stats + `.sum` row | ✅ `io/summary.jl` |
| `EXTREE`/`CVGO`/`MISPRT`/`RDPR`/`BRPR` | tree list, cover, mistletoe/down-wood/snag reports | 🧊 |

---

## GROW step — `GRINCR` then `GRADD`

### `CUTS` — scheduled thinning / harvest (`cuts.f` + `cutstk.f`) → `cuts!`

Runs at the top of `GRINCR` (before `DGDRIV`). The keyword dispatch (`cuts.f` computed-GOTO
on `ICFLAG`/`IACTK`) is **destructured below into atomic methods + modifiers**, each with port
status. `ICFLAG = IACTK − 220` for cut activities; modifiers use `IACTK` 201–206. ⛔ = not
ported (silently ignored today — a real gap, not a no-op). ⚠ = parsed but wrong/partial.

**Cut methods** (each removes TPA/BA over an eligibility class, ranked by `RDPSRT`):

| keyword | ICFLAG | cuts.f label | semantics | status |
|---|---|---|---|---|
| THINBTA | 3 | label_200 | from below to residual **TPA** | ✅ `_thin_sorted!` |
| THINATA | 4 | label_225 | from above to residual **TPA** | ✅ `_thin_sorted!` |
| THINBBA | 5 | label_250 | from below to residual **BA** | ✅ `_thin_sorted!` |
| THINABA | 6 | label_275 | from above to residual **BA** | ✅ `_thin_sorted!` |
| THINDBH | 8 | label_325 | proportional cut of a **DBH class** to residual TPA/BA | ✅ `_thindbh!` |
| THINAUTO | 1 | label_150 | auto-thin to FULSTK when stocking exceeds a trigger | ⛔ **used by no test** |
| THINPRSC | 7 | label_300 | **prescription** thin (per-DBH-class residual table) | ⛔ **USED BY snt01/sn (stand 3)** |
| xSALVAGE | 9 | label_300 | **salvage** dead/damaged trees | ⛔ **USED BY snt01/sn (stand 4)** |
| THINSDI | 10/14/16 | label_400 | thin to a target **SDI** | ⛔ |
| THINHT | 12 | label_325 | thin a **height** class | ⛔ |
| THINMIST | 13 | label_450 | mistletoe (DMR-based) thin | ⛔ |
| THINRDEN | 14 | label_400 | **relative-density** thin | ⛔ |
| SETPTHIN/THINPT | 15 | label_475 | **point** (plot-specific) thin | ⛔ |
| THINQFA | 17 | label_350 | **Q-factor** (per-class diameter-dist) thin | ⛔ |

**Cut modifiers** (set state the cut loop reads — do NOT remove trees themselves):

| keyword | IACTK | cuts.f label | effect on the cut | status |
|---|---|---|---|---|
| SPECPREF | 201 | label_1200 | per-species **cut preference** → reorders RDPSRT (which species go first) | ✅ `_apply_specpref!` + IORDER term in the `_thin_sorted!` key — cut_specpref year-2000 removal bit-exact vs Oracle A (single-species path; the species-GROUP path needs SPGROUP) |
| TCONDMLT | 202 | — | thin-**condition multiplier** (scales the trigger) | ⛔ |
| YARDLOSS | 203 | label_1325 | **yarding loss** → removed-volume accounting | ⛔ **USED BY snt01/sn (stand 4)** |
| SPLEAVE/LEAVESP | 206 | label_1340 | **leave** named species (exclude from cut) | ⛔ |
| SPGROUP | (125) | — | define **species groups** (referenced by SPECPREF/LEAVESP) | ⛔ |
| CUTEFF | (52) | — | default **cutting efficiency** | ⛔ |

| shared cut-loop branch | effect | status |
|---|---|---|
| `CLSSTK` class stocking (TPA jtyp=1 / BA jtyp=2) over eligibility window | budget | ✅ `_clsstk` |
| `RDPSRT` size rank (−DBH below / +DBH above) **+ IORDER[sp]** (SPECPREF) then whole-record removal ×cuteff | selection | ✅ (⚠ tie-break/stable-sort vs oracle not yet reconciled; TCONDMLT/point/density weights default 0) |
| `TREDEL` compact removed (PROB≤0) records | RNG alignment | ✅ `compact_live!` |
| **post-thin DGSCOR traversal order** on the compacted set | stochastic draw alignment | ⚠ diverges after a thin (s29 RNG splits at ICYC=4) — see [[fvsjl-c5-sum-state]] |
| removed-volume columns (rem_tpa/cuft/mcuft/scuft/bdft) | `.sum` reporting | ✅ |

### `DGDRIV` — diameter growth (`dgdriv.f`) → `diameter_growth!`

| branch / condition | effect | status |
|---|---|---|
| `LSTART` (label_100) | **calibration pass** (once): predict DDS for input trees, derive per-species COR, set OLDRN | ✅ `calibrate_diameter_growth!` |
| ↳ IDG=1 or 3 | convert input DG to inside-bark | ✅ |
| ↳ poor-sample abort (csnxx<0) | skip species calibration | 🟡 unverified — needs a sparse-species key to exercise |
| ↳ DGSD≥1 ⇒ set OLDRN residuals (BACHLO draws) | reproducible per-tree noise | ✅ |
| normal cycle: `MULTS(7)` cov, `AUTCOR` | serial-correlation params (σ,ρ) | ✅ `serial_correlation.jl` |
| `DGF(DBH)` | predicted ln(DDS) per tree | ✅ `dgf!` |
| `MULTS(1)` XDMULT | DG multiplier keyword | 🟡 (no MULT keyword in tests) |
| ICYC==1 special | first-cycle DG handling | ✅ |
| `LDGCAL[sp]` | apply species COR or not | ✅ |
| **`LTRIP` true** ⇒ deterministic tripling DG (central/upper/lower × `MISDGF`) | 3 weighted DGs | ✅ `triple_records!` stash |
| **`LTRIP` false** ⇒ `DGSCOR` serial-correlated DG (ssigma, frm, rho, OLDRN) × `MISDGF` | stochastic single DG | ✅ (this is the cyc3+ OLDRN tail source) |
| `MISDGF` mistletoe DG reduction | ×DG | ⚪ no-op without mistletoe |
| **Fort Bragg (IFOR=20)**: dg5 special DG (sp 8/13), ATTEN override (2056/689), special bark (sp 5,6,8,11,13) | longleaf/loblolly growth | ✅ `dgf!`/`dgcons!` (dgf.f:515) — s30 bit-exact |
| **bark source** (BRATIO) for DG / volume / mortality | inside-bark conversion | ✅ unified to one per-stand `calib.bark_a/bark_b` (was two duplicate tables) |

### `HTGF` — height growth (`htgf.f`) → `height_growth!`

| branch / condition | effect | status |
|---|---|---|
| `MULTS(2)` XHMULT | HTG multiplier keyword | 🟡 |
| PROB≤0 ⇒ skip | dead record | ✅ |
| `HTCALC` mode 0 | back out tree AGE from current HT on Chapman-Richards curve | ✅ |
| `HTCALC` mode 9 | 5-yr HT increment from that age | ✅ |
| htmax−hti ≤ 1 ⇒ tiny floor, goto apply | near-max-height cap | ✅ |
| crown-ratio modifier `hgmdcr` (≤1) | scale HTG | ✅ |
| relative-height modifier (AVH>0, relht≤1.5) | scale HTG | ✅ |
| `htgmod` clamp [0.1, 2.0] | bound modifier | ✅ |
| HTG floor 0.1 | min growth | ✅ |
| **HT+HTG > SIZCAP[sp,4] ⇒ cap HTG** | species height cap | ⛔ not in FVSjl; ⚪ no-op now (SIZCAP=999), needed for the SIZECAP keyword |
| `LTRIP` ⇒ repeat caps for upper/lower records | tripled HTG | ✅ |
| HTCONS entry: HTCON from HCOR2 calibration | per-species HT calib | ✅ (`htg_cor`, =0 for snt01) |

### `REGENT` — small-tree growth (`regent.f`) → `small_tree_growth!`

| branch / condition | effect | status |
|---|---|---|
| `LSTART` (label_40) | calibration of small-tree HT | ✅ |
| `lestb` (establishment mode) | grow newly-established trees only (i≥itrnin), random CR draw, FINT≤5 split | ⛔ **C4 regen** |
| d ≥ xmx (=3") ⇒ skip (large tree) | hand off to DGF | ✅ |
| MANAGD==1 ⇒ ddum=1 | managed-stand modifier | 🟡 |
| `lskiph` height-skip vs HTCALC blend | height path select | ✅ |
| htmax−h ≤ 1 ⇒ floor | near-max | ✅ |
| [xmn,xmx] weight blend with large-tree HTG (xwt; d≤xmn or lestb ⇒ xwt=0) | small/large blend | ✅ |
| htgr floor 0.1 | min | ✅ |
| DGSD≥1 ⇒ BACHLO ±noise (uses ESRANN in estab) | reproducible noise | ✅ (main RANN; ESRANN path is estab ⛔) |
| HTDBH inverse (height→dbh) | derive DG from HT growth | ✅ `_htdbh_dbh` (incl. **Fort Bragg IFOR=20 p2/p3/p4 override** for sp 6/8/11/13 — htdbh.f:145) |

### `MORTS` — mortality (`morts.f`) → `mortality!`  *(see DECISION_FLOW.md §3)*

| branch / condition | effect | status |
|---|---|---|
| RMSQD==0 ⇒ reset line; PMSDIL/U rescale | init | ✅ |
| ICYC>1 & |t−TPAMRT|>1 ⇒ reset self-thinning line | per-cycle recompute | ✅ ported (mortality.jl tpa_mort) — resets after thinning/ingrowth; closed stands persist |
| SDIMAX<5 ⇒ background (Hamilton) only | sparse stand | ✅ |
| t>t85d0 ⇒ tn10=t85d10 (over-dense) | strong self-thin | ✅ |
| t55d0<t≤t85d0 ⇒ solve self-thinning line (iterate treeit) | intermediate | ✅ |
| t≤t55d10 ⇒ tn10=t (none) | low density | ✅ |
| per-tree rip = Hamilton ri or rn; XMMULT window (MORTMULT) | rate | ✅ (XMMULT default 1) |
| `VARMRT` distribute (t−tn10) toward suppressed | kill assignment | ✅ |
| QMD-convergence: recompute d10n, re-iterate ≤10 | end-QMD fixed point | ✅ |
| **MSB alternate mortality** (d10>QMDMSB ⇒ MSBMRT) | extra mortality | ⛔ keyword-only (QMDMSB=999 default) |
| **SIZE-CAP mortality** (d+g≥SIZCAP[,1]) | cap big trees | ⛔ keyword-only (SIZCAP=999 default) |
| **BAMAX enforcement** (scale kills until BA≤BAMAX) | density BA cap | ✅ (commit aedecd1 — was the multi-cycle gap) |
| FIXMORT keyword | forced mortality | ⛔ keyword option |
| TPAMRT = surviving TPA | next-cycle reset basis | 🟡 |

### `TRIPLE` + `REASS` — record tripling (`triple.f`) → `triple_records!`

| branch | effect | status |
|---|---|---|
| only when `LTRIP` & ITRN>0 (cyc1-2) | split each live tree into 3 | ✅ TRIPLE_CYCLE_LIMIT=2 |
| central 0.60 / upper 0.25 / lower 0.15 weights | record probs | ✅ |
| weight<0.2 ⇒ break (skip degenerate) | guard | ✅ |
| ITRN*=3, IREC1*=3, REASS reindex | record bookkeeping | ✅ |

### `UPDATE` — apply growth (`update.f`) → inline in `grow_cycle!`

| branch | effect | status |
|---|---|---|
| DBH += DG/bark (bark at pre-growth DBH) | outside-bark DBH | ✅ |
| HT += HTG | height | ✅ |
| NORMHT>0 ⇒ NORMHT grows by HTG·100 (broken-top) | topkill cubic keeps growing | ✅ |
| wki cap (≤PROB) | mortality apply | ✅ |

### `SDICAL` / `SDICLS` — stand max SDI (`sdical.f`)

| branch | effect | status |
|---|---|---|
| BA-weighted SDIDEF over live trees ⇒ xmax | stand SDImax | ✅ `stand_sdimax` |
| LSTART ⇒ include dead (IMC==7) records | initial xmax | 🟡 |
| `CLMAXDEN(SDIDEF,xmax)` climate SDImax reduction | shrink for climate | ⚪ `IF(.NOT.LCLIMATE)RETURN` — no-op without CLIMATE ext; oracle stubs it; FVSjl correct |
| `!LBAMAX` ⇒ BAMAX = xmax·0.5454154·PMSDIU; else back-solve xmax from BAMAX | BA cap | ✅ (computed in `mortality!`) |
| SDICLS ⇒ SDI class / stand stage bounds | class | 🟡 (class only; not in .sum math) |

### `CROWN` — crown-ratio update (`crown` / crownw) → `crown_ratio_update!`

| branch | effect | status |
|---|---|---|
| LSTART ⇒ dub initial crowns | init CR | ✅ |
| ITRN==0 / TPROB≤0 ⇒ early return | empty | ✅ |
| CRNMLT/DLOW/DHI/ICFLG keyword (per sp / group / all) | CR multipliers | ⛔ keyword path (defaults inert) |
| relsdi = SDIAC/SDIDEF·10 clamp[1,12] | density driver | ✅ |
| acrnew via MCREQN form (5 eqn types `imceqn`) | mean CR | ✅ |
| Weibull draw at diameter percentile | per-tree CR | ✅ |
| ±1%/yr change limit; crown-length cap; clamp[10,95] | bound CR | ✅ |

---

## GRADD-only branches not in the growth core

| routine | role | status |
|---|---|---|
| `MPBCUP`/`DFBWIN`/`MISTOE`/`TMCOUP`/`BWECUP` | insect/disease record edits | 🧊 |
| `FMMAIN` | **FFE fire effects + fire mortality** | 🧊 C7 (the `s10_fire`/`fire_*` divergence) |
| `BRTREG`/`RDTREG`/`CLAUESTB`/`ESNUTR` | sprout/planted/natural regen + nutrient hook | ⛔ C4 regen |
| `HTGSTP` | HTGSTOP/TOPKILL keyword height edits | ⛔ keyword (no-op otherwise) |
| `CVGO`/`CVBROW`/`CVCNOP` | canopy cover | 🧊 |

---

## Reading the status at a glance

- **All ✅ in the growth core** ⇒ the C3/C4 hot path is branch-complete for the
  current (uncalibrated-keyword, no-extension) scenario family. The remaining
  numeric residual is the 🟡 serial-correlation OLDRN tail in `DGSCOR`, not a
  missing branch.
- **Every ⛔** is keyword- or chunk-gated. To turn one live and validate it, add
  a scenario that sets its keyword (e.g. `SIZECAP`, `MORTMULT`, `MATUREW` for
  MSB) or its chunk's inputs (regen: `NATURAL`/`PLANT`; fire: FFE keywords), then
  diff against Fortran.
- **⚪** rows are *correctly* absent (gated off in the oracle too); do not "fix"
  them without first enabling their gate (e.g. the CLIMATE extension).
