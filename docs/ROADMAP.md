# FVSjl — roadmap

_Last updated 2026-10-05. Current numbers: **[CERTIFICATION_STATUS.md](CERTIFICATION_STATUS.md)**. Method:
**[DOCTRINE.md](DOCTRINE.md)**. (Supersedes the per-cluster scoping docs `WESTERN_UNPORTED_VARIANTS_ROADMAP.md`,
`WESTERN_ROLLOUT_GOAL.md`, `CS_GOAL.md`/`LS_GOAL.md`/`NE_GOAL.md`, whose port phases are complete.)_

## Goal

Every one of the 24 FVS variants certified as a drop-in replacement: the tiered suite shows zero unexplained cells
(every residual fixed or measured to a named primitive) and the full test suite is green on master. After that, the
code is refactored to a cleaner design with the tiered suite as a lock-step safety net.

## Phases

| # | Phase | State |
|---|---|---|
| 1 | Port every variant's base model + extensions | **done** |
| 2 | Tiered suite with live goldens for every variant | **done** 2026-10-05 (CS/LS/NE/ON/OC/OP were the last six; BC on 10-04) |
| 3 | Drive every variant's tiered cells to zero or to named corners | **in progress** — see below |
| 4 | Certification pass: re-audit every CORNER entry against the live oracle, close or justify each OPEN entry, full suite on the final tip | not started |
| 5 | Refactor to a cleaner design (jl-vs-jl bit-identity snapshots + oracle fast tier + coverage + mutation + thread-determinism as the net) | after phase 4 |

## Phase 3 — remaining work, in priority order

Priority set by the user on 2026-10-05: **eastern variants first** (CS, LS, NE, ON).

| Priority | Variants | Off / compared | What is known | Estimate |
|---|---|---|---|---|
| 1 | CS | 12,952 / 427k | Mostly `none`-regime per-record drift from cycle 1, starting at the calibration-initialised serial-correlation residual (dgdriv.f OLDRN) | 2–3 days |
| 1 | LS, NE | 4,276 / 393k, 2,413 / 461k | Same OLDRN drift; eastern FFE reports; LS plant summary and FVS54 | 2–3 days |
| 1 | ON | 19,106 / 358k + 40 crashes | ON FFE and ON establishment were unported (establishment now ported on its branch); then growth drift on dense stands | 3–5 days |
| 2 | OC, OP | 15,700 / 193k, 19,545 / 173k | OC startup DG calibration (oc/dgdriv.f) never runs; OP FVS-native paths use Julia math where op/*.f call glibc expf/logf/powf | 3–5 days, in parallel |
| 3 | Western tail (WS, CA, AK, CI, CR, KT, UT, SO, WC, EC, NC, PN) | 2,679 / 6.0M | Mostly FFE salvage/simfire down-wood/flame cells; CR mistletoe DMR | 1–2 days |
| 3 | BC | 57 / 1.97M (pending merge) | post-fire down-wood carbon (+40% ten years after a burn) | <1 day |

Shared items found during phase 3, owned by the western-tail round: FVS_Carbon `Total_Removed_Carbon` never written
(BIOREM(2) harvest terms); FLAMEADJ with FLMULT≠1 skipped the crown-fire path; SO activity-fuel consumption; EC/NC/WC
root-disease cycle-1 WK1; EM NOAUTOES+PLANT planted birth heights; remaining `NORMHT*0.01` sites.

## Timeline

| | Estimate (2026-10-05) |
|---|---|
| Phase 3 complete | 4–7 days |
| Phase 4 (certification pass) | +1–2 days |
| **Certification** | **~5–9 days from 2026-10-05** |
| Phase 5 refactor | +1–2 weeks |

History of the estimate: ~2–3 weeks (2026-09-28, before the API outage), 1.5–2.5 weeks (2026-10-04 restart),
5–9 days (2026-10-05). It shortened because BC (467,601 → 57) and the core four (→ 0) closed far faster than planned.

## Risks

1. **ON and the ORGANON variants (OC/OP)** are the least predictable: ON's FFE is a full subsystem port, and ORGANON
   has been hard every time it was touched.
2. **Integration cost.** Several work streams touch the shared engine (FFE, mortality, volumes), so each merge brings
   out cross-branch interactions. Each is caught by the tiered + full-suite gate before master, but each costs an
   extra ~1 h suite run (more when the machine is saturated by concurrent sweeps).
3. **Diminishing returns.** What remains is increasingly last-digit Float32 drift and rarely exercised paths; each fix
   moves fewer cells, and some residuals may end as measured CORNERs rather than fixes.
4. **Coverage beyond the fixtures.** The tiered stands are a stratified sample. Phase 4 includes re-running the FIA
   full-population sweeps on the final engine for the variants that have FIA plots.
