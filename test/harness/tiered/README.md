# Tiered integration suite

The gate that decides whether FVSjl still matches the live Fortran FVS — across **variants × regimes × outputs** —
and the **safety net for refactoring** once variants are certified. It replaces reliance on the SN-only
`test_multicycle` smoke gate (10 SN scenarios, 5 `.sum` columns).

## Tiers

| Tier | What it proves | How to run | When |
|---|---|---|---|
| **Fast** (oracle) | jl output == committed **live-oracle goldens** at print precision: every `.sum` field of every row, every DBS table present/absent, every column of every row (numeric-exact at FVS's stored precision, text exact); per variant×regime **signed tally** (one-directional material bias fails) | `julia --project -e 'using Test, FVSjl; include("test/integration/test_tiered.jl")'` (also in `Pkg.test`) | every merge |
| **Quick** | the fast tier on 3 stands × 5 key regimes per variant (`none, thinbba, plant_cyc, simfire, mistletoe`) | `TIERED=quick` before either command | while iterating |
| **Snapshot** (jl vs jl) | **bit-identity** with jl's own blessed outputs: `.sum` rows, every DBS table at full precision, and the full per-record state (every TreeList + Scratch vector, raw Float32 bits) at every cycle | `julia --project=. test/harness/tiered/snapshot.jl check` | every refactor step (quick) / before merge (full) |
| **Coverage gate** | every ported function the fixtures exercised at baseline is still exercised; zero-coverage functions are *listed*, never silent | `test/harness/tiered/coverage_check.sh` | before merge |
| **Thread determinism** | snapshot results identical at 1 and N threads | `TIERED_THREADS=1 … snapshot.jl check` and `TIERED_THREADS=8 … snapshot.jl check` | before merge |
| **Mutation** | the snapshot net catches realistic subtle changes (1-ULP coefficient, a*b*c association, Float64 intermediate, tie-break flip, libm swap) — reports the kill rate | `julia --project=. test/harness/tiered/mutation.jl` | when changing the suite |
| **Slow** | (1) regenerated goldens from the *current* oracles are identical to the committed ones (golden drift ⇒ oracle changed); (2) population regime sweeps (400/120 stratified stands) + signed tallies | `test/harness/tiered/slow_tier.sh [V …]` | before master merge, nightly |

Why both oracle and snapshot tiers: the oracle tier allowlists known residuals (`KNOWN_RESIDUALS.toml`), so a change
that shifts an already-allowlisted OPEN/CORNER cell, or drifts below print precision, is invisible to it. The
snapshot tier sees every bit.

## Fixtures (`test/fixtures/tiered/<v>/`, git-tracked)

Built by `test/harness/tiered/make_fixtures.jl <VARIANT> <K>` from the **live oracle** (`BIN` in
`test/harness/fia/ledger_fia.jl`; the FIA-DB path reads trees via DBSTREESIN, so it is unaffected by the g16
`-std=legacy` TREEFMT back-tab read bug):

* `stands.db` — K stratified FIA stands (extract_sample.jl stride over ECOREGION/LOCATION), self-contained;
* `<cn>_<regime>.key` — one keyfile per stand × regime (`none, thinbba, salvage, plant_cal, plant_cyc, simfire,
  mistletoe, rootdis, climate, econ, cover`), regime cards = `ledger_fia.jl regime_block` verbatim, plus per-regime
  report requests (TREELIST/CUTLIST/STRCLASS/FFE reports + DATABASE toggles, the two-DATABASE-block emission recipe);
  `plant_cal` uses a calendar PLANT date, `plant_cyc` the cycle-number form (both real user forms);
* `<cn>_<regime>.live.sum`, `<cn>_<regime>.<Table>.csv`, `<cn>_<regime>.tables` — the goldens;
* `PROVENANCE.toml` — oracle path + sha256, build script and whether `main.f` was built `-std=legacy`, FVS source
  commit + dirty flag, generator commit, date, stand list;
* `SNAPSHOT.tsv` — the blessed bit-identity manifest (one SHA-1 per stand×regime×{sum, table, cycle}).

WPBR is not in the matrix: the standard `FVS{v}_g16` oracles are not linked with BRUST; `test_wpbr` carries its own
live goldens.

## Allowlist rules (`test/fixtures/tiered/KNOWN_RESIDUALS.toml`)

Each entry: `variant/stand/regime/file/col/year` (`"*"` wildcards) + `status = "CORNER"` with `proof` (named FP
primitive + per-record evidence) **or** `status = "OPEN"` with `issue` (the tracked bug/owner), + `max_cells` (the
measured count). The fast tier FAILS on: any unlisted mismatch; an entry matching **more** cells than `max_cells`
(a grouped entry cannot hide spread); an entry matching **nothing** in a full run (unexpected pass — remove it); an
entry without status/proof/issue; a missing fixture (never a skip); a jl crash (cases run in crash-isolated worker
subprocesses — a SIGSEGV is recorded as a `CRASH` mismatch, not a lost run). Entries are drafted from measurement
(`report.jl` → `classify.jl` → `allowlist_draft.jl`), never by hand-waving; every `NEW:` tag is a finding.

## Refactor workflow (after a variant is certified)

1. **Bless at certified master:** `git checkout master && julia --project=. test/harness/tiered/snapshot.jl bless
   "certified <what>"` → commit the `SNAPSHOT.tsv` files. Record the coverage baseline:
   `test/harness/tiered/coverage_check.sh baseline` → commit.
2. **Refactor** on a branch, in small steps. After each step: `TIERED=quick julia --project=.
   test/harness/tiered/snapshot.jl check` (minutes). Any `DIFF` ⇒ `snapshot.jl detail <V> <stem>` names the first
   differing cycle / variable / record with blessed-vs-now Float32 bits.
3. **Before merge:** full `snapshot.jl check` at `TIERED_THREADS=1` **and** `=N`; the oracle fast tier
   (`test_tiered.jl`, full); `coverage_check.sh`; full `Pkg.test`.
4. **Re-bless only intentionally** — a real model fix that changes outputs: fix → oracle fast tier improves (tighten or
   remove the KNOWN_RESIDUALS entries) → `snapshot.jl bless "<why>"` (the manifest records commit + reason).

## Files

| File | Role |
|---|---|
| `tiered_common.jl` | regimes, keyfile builder, fixed-width `.sum` layout (FORMAT 9014, 29 fields), DBS table I/O |
| `tiered_runner.jl` | run a case, compare vs goldens, signed tally, allowlist, crash-isolated `run_isolated` |
| `tiered_worker.jl` | worker subprocess (compare / snapshot modes) |
| `make_fixtures.jl` | build a variant's fixture from the live oracle |
| `report.jl`, `classify.jl`, `classify_fns.jl`, `allowlist_draft.jl` | measure, classify, and draft allowlist entries |
| `snapshot.jl` | bless / check / detail / dump (bit-identity tier) |
| `coverage_gate.jl`, `coverage_check.sh` | coverage baseline + gate |
| `mutation.jl` | mutation kill-rate check |
| `slow_tier.sh` | golden-drift + population sweeps |
