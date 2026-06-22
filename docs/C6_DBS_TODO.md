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

> When C6 starts: do G1/G2 first (they make VOLS branch-complete), then the Compute/
> Carbon tables that consume them, then the rest of the DBS schema.
