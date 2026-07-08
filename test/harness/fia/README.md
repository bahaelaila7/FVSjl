# FIA-DB validation harness (dev-only; not part of the suite)

`validate_fia.jl` runs FVSjl's native FIA reader vs live FVS on real stands from an
FIA "FVS-ready" SQLite database, per stand, and diffs the `.sum` (cycle-0 bit-exact +
per-cycle mean |rel diff|).

Requires (NOT in CI): the FIA DB (e.g. /workspace/SQLite_FIADB_ENTIRE.db, opened
READ-ONLY) and the relinked live binary /tmp/FVS{ls,sn,ne,cs}_new. Stand list =
tab-separated `STAND_CN<TAB>tag` lines.

The self-contained regression lock lives in test/integration/test_fia_reader.jl
(fixture test/fixtures/fia/ls_sample.db + golden live .sum) and DOES run in the suite.
