# FVSjl — certification status

_Last updated 2026-10-05 (master `f0272193`). Roadmap and remaining work: **[ROADMAP.md](ROADMAP.md)**.
Method: **[DOCTRINE.md](DOCTRINE.md)**. Older narrative status (FIA sweeps, extensions, PPE): [PORT_STATUS.md](PORT_STATUS.md)._

## What "certified" means here

A variant is certified when the **tiered suite** (`test/integration/test_tiered.jl`,
[test/harness/tiered/README.md](../test/harness/tiered/README.md)) reports **zero unexplained cells** for it: every
cell of the `.sum` and of every DBS table that live FVS writes, for every stand × regime in its fixture, equals the
live Fortran oracle at print precision — or is an allowlist entry (`test/fixtures/tiered/KNOWN_RESIDUALS.toml`) whose
cause is measured down to a named primitive (`CORNER`). `OPEN` entries are tracked bugs; a variant with OPEN entries
is not certified. The full `Pkg.test()` must be green on the exact master tip.

Fixtures: 10–12 stratified stands per variant (FIA where the variant has FIA plots; FVS's own test databases for BC;
converted neighbour-variant FIA stands for OC/OP/ON), each run under the regimes its oracle build supports — none,
thinbba, salvage, plant (2 date forms), simfire, mistletoe, climate, rootdis, econ, cover — with goldens produced by
the live Fortran (`FVS{v}_g16`, or the private crash-fixed builds listed below).

## Current state (24 variants)

"Off" = tiered cells that differ from live; "compared" = every golden cell checked (`.sum` rows × fields + Σ table
rows × columns; script `/workspace/.postswap/denom.jl`). A crash counts as one off cell, so crash-heavy rows overstate
exactness.

| Group | Variant | Cases | Cells compared | Off (2026-10-04 AM) | Off (master now) | % exact | Where the rest is |
|---|---|---:|---:|---:|---:|---:|---|
| Core | SN | 132 | 466,237 | 62 | **0** | 100% | — |
| | IE | 132 | 971,091 | 926 | **0** | 100% | — |
| | EM | 132 | 690,049 | 698 | **0** | 100% | — |
| | BM | 132 | 280,949 | 9,350 | **0** | 100% | — |
| Western | TT | 110 | 511,958 | 53,551 | **0** | 100% | — |
| | NC | 110 | 245,530 | 23,002 | 1 | 99.9996% | rootdis |
| | PN | 110 | 509,215 | 7,888 | 1 | 99.9998% | climate |
| | EC | 110 | 314,151 | 24,789 | 7 | 99.998% | rootdis / FFE |
| | WC | 110 | 341,295 | 1,757 | 10 | 99.997% | simfire |
| | SO | 110 | 819,712 | 56,238 | 41 | 99.995% | FFE |
| | UT | 110 | 517,866 | 43,333 | 77 | 99.985% | FFE, rootdis |
| | KT | 110 | 665,474 | 82,984 | 97 | 99.985% | FFE |
| | CR | 110 | 431,194 | 62,672 | 107 | 99.975% | mistletoe |
| | CI | 110 | 214,008 | 3,020 | 179 | 99.916% | salvage/simfire down wood |
| | AK | 90 | 431,138 | 2,283 | 272 | 99.937% | FFE |
| | CA | 100 | 434,154 | 14,475 | 670 | 99.846% | FFE |
| | WS | 110 | 560,098 | 29,758 | 1,217 | 99.783% | FFE |
| Eastern | NE | 70 | 461,313 | not measured | 2,413 | 99.48% | none-regime per-record drift, FFE |
| | LS | 70 | 392,704 | not measured | 4,276 | 98.91% | none, FFE, plant |
| | CS | 70 | 427,228 | not measured | 12,952 | 96.97% | none-regime per-record drift |
| | ON | 60 | 357,738 | not measured | 19,106 + 40 crashes | 94.66% | growth drift; FFE + establishment unported |
| ORGANON | OC | 100 | 192,907 | not measured | 15,700 | 91.86% | startup DG calibration unported |
| | OP | 100 | 173,264 | not measured | 19,545 | 88.72% | Float32 math in the FVS-native paths |
| Canada | BC | 90 | 1,968,465 | not measured | 467,601 (master) · **57** (pending `integ-1005e`) | 76.2% · 99.997% | post-fire down-wood carbon |

Totals on master: core 4 — 0 / 2.41M; 13 original western — 2,679 / 6.0M (99.955%).

CS/LS/NE/OC/OP/ON entered the tiered suite on 2026-10-05; before that they were validated only by shipped-test goldens,
keyword-coverage tests and (CS/LS/NE) the FIA full-population `.sum` sweep, which a per-cell tiered comparison is much
stricter than.

## Pending integration

| Branch | Content | State |
|---|---|---|
| `integ-1005e` | `bc-port`: BC FFE port, Kozak volume, metric DBS/keyword/ESTAB conversions, BC growth/crown/REGENT/sprouting, MISTOFF, shared COMPRESS/RDTRP fixes — BC 467,601 → 57 | measuring, then suite |
| `east-resid` | CS/LS/NE residual round 1 | final measurement |
| `oc-op-on` | ON establishment, OC/OP site setup and math | in progress (ON first) |
| `west-tail` | western tail + FLAMEADJ crown fire, removed carbon, activity fuels | wrapping up |

## How work lands (the gate)

1. A fix is measured against the live oracle and committed with the Fortran file:line it follows, plus a live-fixture
   regression test.
2. Branches are combined on an `integ-*` branch with `git merge --no-ff` (never cherry-pick or rebase); conflicts that
   are two ports of one Fortran mechanism are resolved to one implementation.
3. Tiered ALL is measured on the combined tip, the allowlist is reconciled from that measurement (stale entries
   dropped; caps raised only with the measured reason), and the TOML is parsed with `tomllib`.
4. The full `Pkg.test()` runs on that exact commit; failures are compared by name with the previous green run. Any new
   failure is traced to a cause before merging — a test is changed only when a live measurement shows the old
   expectation was wrong.
5. Master only moves by `--no-ff` merge of a green integration branch.

Recent examples of what step 4 caught before master: an FMMAIN pass running ahead of the mistletoe spread (snag
random draws 275 early), ON emptied-record merch volume overwritten by a generic PROB=0 rule, a closure that boxed
captured variables and made SN/NE/CS/LS mortality allocate 471 KB/cycle, and an out-of-bounds species-table read in
KT volume setup.

## Integration history since the 2026-09-28 → 10-04 outage

| Master | Date | What landed | Suite (pass / fail / err / broken) |
|---|---|---|---|
| `279c33b0` | 10-04 | BM regime residuals (`bm-resid`) | 123,953 / 0 / 0 / 23 |
| `ca2483b9` | 10-04 | western shared round 2 | 123,957 / 0 / 0 / 23 |
| `fbf0d819` | 10-04 | CORE round 1 | 124,365 / 0 / 0 / 23 |
| `cb83ec2b` | 10-04 | western crash sweep (180 → 0 crashes) | 138,198 / 0 / 0 / 22 |
| `6c98a149` | 10-05 | KT/WS/CA/AK + TIMEINT 5 | 138,331 / 0 / 0 / 22 |
| `e869f8f5` | 10-05 | western round 3, KT/WS/CA/AK/BC round 2, CORE round 2 | 139,448 / 0 / 0 / 21 |
| `bb2f1f8a` | 10-05 | western round 4 + FMMAIN after MISTOE | 139,495 / 0 / 0 / 20 |
| `f0272193` | 10-05 | all 24 variants in the tiered suite + CORE round 3 (core at 0) | 139,553 / 0 / 0 / 19 |

## Oracles

Shared live oracles: `/workspace/.{v}work/FVS{v}_g16` (gfortran-16 rebuilds of each variant's buildDir; CR
`FVScr_clean`, KT `FVSkt_clean`). Private relinks used for goldens where the shared binary crashes — each is the shared
build's objects plus a call-site patch, documented in `/workspace/ORACLE_SOURCE_AUDIT_2026-09-19.md`:

| Oracle | Patch | Audit § |
|---|---|---|
| `/workspace/.bcwork/dbfix/FVSbc_dbfix` | DBSTREESIN 34-argument call (every DB tree read segfaulted); CFTOPK SCF argument (FFE/top-kill SIGFPE) | 7 |
| `/workspace/.ocwork/tiered/FVSoc_tiered`, `/workspace/.opwork/tiered/FVSop_tiered` | ORGANON mortality underflow trap after heavy thins | 8 |
| `/workspace/.onwork/tiered/FVSon_tiered` | FMSVL2 / CFTOPK argument lists (SIGSEGV, garbage top-kill volume) | 8 |

Shared oracle binaries and buildDir sources are never modified without explicit approval.
