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
| **COMCUP zero-PROB delete** | `base/comcup.f` (top, before COMPRS) | every cycle, TREDEL records with `PROB ≤ 1e-5` (suppressed trees whose expansion → ~0) so the next cycle's DGSCOR doesn't draw an extra deviate for a dead record | snt01 live records never reach 1e-5 → bit-exact without it; dense/long stands trigger it | ✅ `comcup!` (trees.jl) in grow_cycle!; dense_long 30-cyc scenario validates ±1 vs oracle (test_longrun.jl) |
| **ESFLTR** | `estb/esfltr.f` | sets BAAINV/TPAAINV (understory/overstory densities) + IESTAT best-tree flags that feed the AUTO-establishment probability | SN establishment is keyword-driven (PLANT/NATURAL give TPA explicitly) → the auto-estab probability path is bypassed; Oracle A stubs ESFLTR and stays bit-exact | ✅ confirmed no-op for SN (matches oracle; would only matter for auto-ingrowth, which SN doesn't use) |
| **NPTIDS>1 replication** | `base/estab.f` (nn-loop over IPTIDS) | establishment replicated per inventory point (multi-point) | bare scenarios are single-point | ✅ per-point loop + plot_id=nn (point_ba×PI→ba_v); bare_multipoint/mp3 match oracle |

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

1. ✅ **COMCUP zero-PROB delete** — DONE (dense_long 30-cyc scenario, ±1 vs oracle).
2. ✅ **NPTIDS>1** — DONE (bare_multipoint 5-pt / bare_mp3 3-pt match oracle).
3. ✅ **ESFLTR** — confirmed no-op for SN (auto-estab probability path, bypassed by keyword PLANT/NATURAL).
4. Classification (SSTAGE/SDICLS/SILFTY) + PCTILE/DIST/COMP — belongs with C5-detail/C6 output.
5. DG calibration tail — pick one calibrated outlier species, reconcile per-tree.

> **Natural DYNAMICS are now complete** (COMCUP + NPTIDS>1 ported; ESFLTR confirmed no-op for SN).
> Remaining natural items are the DG-accuracy tail (#5) and classification/reporting (#4 → C5-detail/C6),
> which are output, not simulation dynamics.
