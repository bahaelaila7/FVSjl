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
- ⛔ The remaining ~16 tables (Compute, Carbon, …) + the FVS_Cases full schema.

> Next: TreeList (per-tree records) + Compute (event-monitor vars), then G1/G2-dependent
> Carbon. The Summary writer establishes the pattern (CREATE TABLE IF NOT EXISTS + prepared
> INSERT, gated on the DATABASE keyword, validated vs the Fortran .db).
