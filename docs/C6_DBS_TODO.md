# C6 (DBS / SQLite output) tracker

Items deferred to C6 because they are consumed by the DBS database I/O + output
tables, not the natural-process `.sum`. Carried here so they are not forgotten when
C6 is built. Each is fully scoped in its source tracker.

## From the C3/C4/C5 semantic audit (`SEMANTIC_AUDIT_C3C4C5.md`)

| item | what | validation target (C6) |
|---|---|---|
| **G1 — CULL / DEFECT** volume reduction | cull%/defect% reduce merch/board volume; come only from DBS *database* input (text `.tre` reads 25 fields, both = 0) | a DBS-input scenario with cull/defect trees → cuft/bdft vs live Fortran |
| **G1 — DECAYCD / WDLDSTEM defaults** | `intree` DECAYCD=3-for-dead / WDLDSTEM=1-woodland; affect dead-tree carbon + woodland multi-stem | DBS Carbon + per-tree state vs oracle |
| **G2 — biomass / carbon** | JENKINS + WOODDEN (Jenkins 2003 AGB + wood density, WDBKWT 2677-species table) fill ABVGRD/MERCH/CUBSAW/FOLI bio + carbon | DBS Carbon tables + per-tree biomass dump vs Oracle A |

## From the natural-process tracker (`NATURAL_PROCESS_TODO.md` §C)

| item | what | C6 table |
|---|---|---|
| PCTILE / DIST / COMP | per-attribute percentile + DBH-class distribution + species composition | DBS Compute / `.out` detail |
| SSTAGE | stand-structure-stage code | DBS / `.out` |
| SDICLS detail (SDIBC/SDIAC) | SDI class column (the `.sum` class code already matches) | `.out` detail |
| SILFTY | silvicultural forest type code | `.out` detail |

## Native DBS tables (the C6 chunk proper)
The 18+ base DBS output tables (Summary, TreeList, Compute, Carbon, …) — see the
existing FVSjulia sndb table set as the parity target. (Fire/econ DBS tables are C7/C8.)

**C6 STARTED (2026-06-24):**
- ✅ **FVS_Summary** (`src/io/dbs_output.jl`, `write_dbs_summary!`) + **FVS_Cases** — the
  per-cycle summary written to the DSNOUT SQLite db. The **DATABASE** block keyword is ported
  (`kw_database!`: `DSNOUT` filename on the next line, `SUMMARY` ⇒ ISUMARY); `run_keyfile`
  collects the `SummaryRow`s (via `write_sum_file`'s `collect_rows`) and appends them. **Bit-
  exact vs the live Fortran `FVSOut.db`** (`test_dbs_summary.jl`, `dbs_summary.key`: every
  column — Tpa/BA/SDI/CCF/TopHt/QMD/MCuFt/BdFt — matches, only ±1 cuft Float32 noise). This is
  the same data the text `.sum` carries, into a database (the "same SQLite outputs" goal).
- ✅ **FVS_TreeList** (`write_dbs_treelist!` + `treelist_snapshot`) — the per-cycle, per-tree
  detail table (dbstrls.f). `write_sum_file` gained a `cycle_hook` that snapshots the start-of-
  cycle tree list (DATABASE `TREELIDB` ⇒ `dbs_treelist`); columns map directly from the
  `TreeList` struct (species FVS/PLANTS/FIA, TPA, DBH, DG, Ht, HtG, PctCr, CrWidth, BAPctile,
  PtBAL, TCuFt/MCuFt/SCuFt/BdFt, TruncHt, the merch-top heights Ht2TDCF/BF, TreeAge). **Per-tree
  TPA = `t.tpa/gross_space`** (= Fortran `PROB/GROSPC`; the validation caught a missing /g).
  Validated: each cycle's Σ(TPA) and Σ(TCuFt·TPA) reconstruct the Fortran-bit-exact `.sum` stand
  TPA / cubic volume (`test_dbs_treelist.jl`). The exact per-tree row set differs from Fortran
  only by the tripling/COMCUP record PARTITION (same totals). Not-yet-filled cols (nullable):
  MortPA, TreeVal/SSCD/PtIndex, MistCD, MDefect/BDefect split, EstHt, ActPt.
- ✅ **FVS_Compute** (`write_dbs_compute!`, dbscmpu.f) — DONE. Dynamic schema (one REAL col per
  COMPUTE var, declaration order); `snapshot_compute!` evaluates the active COMPUTE defs at each
  GROWING cycle's start (`write_sum_file`'s `compute_collect`; the final inventory cycle gets no
  row, matching Fortran — the event monitor runs only during growth). DATABASE `COMPUTDB` ⇒
  `dbs_compute`. **Bit-exact (Float32) vs live Fortran FVSOut.db** (`test_dbs_compute.jl`):
  MYBA=BBA and MYSDI=BSDI match every cycle. The earlier "Fortran writes only 1 row" finding was a
  **misdiagnosis** — the cadence is date-gated by the COMPUTE date; `COMPUTE 0` (active every
  cycle) writes every growing cycle. This also validated the **BSDI = raw Reineke SDIBC fix**
  end-to-end through the DBS path (MYSDI 202.94, not the BA-77 of the old copy-paste bug). NOTE: a
  bare `TPA` event var exists in FVSjl but not in stock Fortran (Fortran leaves `MYTPA=TPA` null),
  so it is not a valid parity column.
- ✅ **FVS_InvReference** (`write_dbs_invref!`, dbsinvref.f) — DONE. A once-per-stand dump of the
  variant's 90-species master list: FVS/PLANTS/FIA codes, SDI method (`zeide_sdi` → "ZEIDE"),
  per-species SDImax + site index, cubic/board volume-equation ids and merch specs (min DBH / top
  dia / stump for total / sawtimber / board). All from engine state after `compute_volumes!`;
  emitted whenever the DSNOUT db is active. **Bit-exact vs live Fortran FVSOut.db — 0 mismatches
  across 90 species × 19 columns** (`test_dbs_invref.jl`); incl. the 9 woodland species' CFMinDBH=6.
  One fix: SiteIndex/SDIMax use FVS `trunc(x+0.5)` (round-half-up), not Julia `round` (half-even).
- ⛔ The remaining ~14 tables + the FVS_Cases full schema. Findings:
  - **FVS_Carbon** needs belowground / forest-floor / shrub-herb carbon pools FVSjl doesn't
    compute (only aboveground-live via Jenkins + the FFE dead pools) — the G2 biomass chunk.

> Next: TreeList (per-tree records) + Compute (event-monitor vars), then G1/G2-dependent
> Carbon. The Summary writer establishes the pattern (CREATE TABLE IF NOT EXISTS + prepared
> INSERT, gated on the DATABASE keyword, validated vs the Fortran .db).
