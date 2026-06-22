# Natural-process completeness tracker

Items in **C3 (growth) / C4 (mortality+density+regen) / C5 (volume+biomass+stats)**
that are natural-process (NOT management/disturbance/keyword-option) and are still
unported, partial, or unvalidated. Driven from the FVS **decision flow** as the
completeness oracle — NOT from passing tests. A row is "done" only when the FVS
logic is **genuinely ported from the Fortran** (Oracle A `/workspace/FVSjulia/src`)
AND a Fortran-validated coverage scenario exercises its branch.

> **Discipline (user rule):** port the real FVS code for each item; use tests only to
> *catch and fix bugs*, never to back-fill logic that just makes a test pass.

Legend: ⛔ unported · 🟡 partial/simplified · 🔬 ported-but-accuracy-tail · ✅ done

## A. Unported natural DYNAMICS (real mechanisms, no current coverage)

| item | FVS source | what it does | why untested today | status |
|---|---|---|---|---|
| **COMCUP / COMPRS** | `base/comprs.f` (+EIGEN!), `base/comcup.f` | record compression when live records overflow `MAXTRE÷3` (auto) or on COMPRESS keyword | snt01 caps at 243 records → never triggers | ⛔ |
| **ESFLTR** | `base/esfltr.f` | establishment filter — flag "best tree" records for the estab model | regen scenarios don't hit it | ⛔ |
| **NPTIDS>1 replication** | `base/estab.f` (nn-loop over IPTIDS) | establishment replicated per inventory point (multi-point) | bare scenarios are single-point | 🟡 (single-point done) |

## B. Natural ACCURACY tails (ported, not bit-exact)

| item | where | symptom | status |
|---|---|---|---|
| per-species DG calibration | `sn/dgf.f` calib | sp33/65, single-species `all_*`, s31 loblolly off by ~1 | 🔬 |
| Float32 / sawtimber-threshold | volume | ±1 ulp in cuft/bdft; one boundary tree flips sawtimber class | 🔬 (irreducible) |

## C. Classification / reporting (in C4/C5 but OUTPUT, not dynamics)

| item | FVS source | feeds | status |
|---|---|---|---|
| SSTAGE | `base/sstage.f` | stand-structure-stage code (.out/DBS) | ⛔ |
| SDICLS → SDIBC/SDIAC | `base/sdicls.f` | SDI class column (.out); `.sum` class code already matches | 🟡 |
| PCTILE / DIST / COMP | `base/pctile.f`/`dist.f` | `.out` detail + DBS Compute tables | 🟡 |
| SILFTY | `base/silfty.f` | silvicultural forest type (.out) | 🟡 |

## D. Done this stretch (kept here so docs stay honest)

- ✅ MORTS self-thinning **line-reset** (`|tt−TPAMRT|>1`) — mortality.jl:229-233 (was 🟡)
- ✅ post-thin **DGSCOR record order** — TRIPLE interleaved append (was ⚠-diverges)
- ✅ regen **birth-cycle volume** (cyc1 cuft=0), regen **point-BA** (single inventory point)
- ✅ Fort Bragg **KODFOR→81110** volume remap
- ✅ **multi-stand driver** (`each_stand`/`run_keyfile`): TREFMT persistence + default INTREE

## Order of work (least-dependent first, each: port → scenario → validate)

1. **COMPRS / COMCUP** — dense/long-stand scenario that overflows records; port `comprs.f`.
2. **NPTIDS>1** — multi-point establishment scenario; port the per-point estab loop.
3. **ESFLTR** — regen scenario; port the establishment filter.
4. Classification (SSTAGE/SDICLS/SILFTY) + PCTILE/DIST/COMP — belongs with C5-detail/C6 output.
5. DG calibration tail — pick one calibrated outlier species, reconcile per-tree.
