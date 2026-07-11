# Pillar-1 stratified FIA sample manifests (per variant)

Deterministic, reproducible stratified samples of real FIA plots from the read-only master
`SQLITE_FIADB_ENTIRE.db` (`FVS_STANDINIT_COND`, `VARIANT` column). Each `<v>_manifest.txt` is
`STAND_CN<TAB>VARIANT`, 500 stands/variant (2000 total — materially larger than the 162-stand
modernization baseline). These are the Pillar-1 sample manifests for the FIA/FVS behaviour-compatibility
campaign (docs/FIA_FVS_COMPAT_AUDIT.md).

## Strata & method
Stands are ordered by `(ECOREGION, LOCATION, STAND_CN)` and every K-th is taken (even stride), so the
sample spreads across ecoregions (the ecological unit driving DG EUT terms + species/geography) and
national forests. No RNG — fully reproducible.

| variant | population | sampled | distinct ECOREGION | distinct LOCATION |
|---------|-----------:|--------:|-------------------:|------------------:|
| SN      |    637,641 |     500 |                170 |                76 |
| NE      |    178,149 |     500 |                111 |                 6 |
| CS      |    255,952 |     500 |                 93 |                 3 |
| LS      |    400,649 |     500 |                 96 |                 8 |

## Regenerate
```
julia --project=. test/harness/fia/extract_sample.jl <SN|NE|CS|LS> 500 test/harness/fia/manifests/<v>_manifest.txt
```
`extract_sample.jl` is the extraction script (read-only on the DB). Run the projection differential over a
manifest with `ledger_fia.jl <manifest> <V> none` (Pillar-2) or `manage_fia.jl <manifest> <V> <regime>` (Pillar-3).
