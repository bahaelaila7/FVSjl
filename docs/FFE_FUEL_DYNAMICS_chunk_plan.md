# FFE surface-fuel dynamics — chunk plan (the grown-cycle carbon remainder)

The **Stand Carbon Report** is now bit-exact at the **inventory cycle** (every column: Jenkins live
above/merch/below, snag, DDW ×0.5, forest floor ×0.37, shrub/herb via FULIV2, total — see
`test_carbon.jl` + `carbon.jl`). The remaining gap is the **grown-cycle** FFE pools: the Fortran
report decays the down-wood/litter/duff each cycle (e.g. DDW 3.8→2.5, Floor 9.1→6.6 from 1990→1995)
and accumulates new dead wood from tree mortality. That is a distinct FFE **surface-fuel-dynamics
subsystem**, not part of the carbon report itself — this plan scopes it from the Fortran flow.

## What FVSjl already has
- `fmcba!` (fire/fmcba.jl) — the per-cycle FFE cover-type + **initial** dead-fuel loading (FUINI) and
  live herb/shrub loading (FULIV / FULIV2 override). It loads `fire.cwd` once (`fuels_init`) and
  re-sets `fire.flive` each call. ✅ Correct for the inventory cycle.
- `fuel_model.jl` (FMCFMD3 dynamic fuel model), `fmburn.jl` (the fire-event behavior), `snag.jl`
  (FMSFALL snag falldown + hard→soft decay). These run **only on a fire event**.

## What's missing (the subsystem to port), most-upstream first
1. **Per-cycle down-wood decay** — the `fire.cwd` woody/litter/duff pools decay by size-class- and
   decay-class-specific annual rates each cycle (fire/base: trace where CWD is multiplied down between
   cycles — start at `fmcfir.f` (373 ln) and `fmsadd.f` (447 ln); confirm the decay-rate source/table
   and whether it is in `fmcblk`/`fmvinit`). Validation: the report's DDW + Floor columns at 1995/2000.
2. **Mortality → dead-wood accumulation** — trees that die each cycle add to the down-wood pools
   (FMSADD). Couples to the existing periodic-mortality + snag path. Validation: DDW growth offsets the
   decay in mixed cycles; snag→DDW transfer on falldown (`snag.jl` already has falldown timing).
3. **Per-cycle FFE driver** — call the fuel update (`fmcba!` for live/cover + the new decay/accum) every
   cycle for FFE-active stands, from `grow_cycle!` (currently only `fmburn!` runs, and only on a fire
   year). Order vs growth/mortality must match the Fortran FFE main (`fmmain.f`) sequence.
4. **`stand_carbon_report` per-cycle emission + the `.out` report WRITER** — byte-exact like
   `write_structure_report` (SSTAGE): FORMAT headers + the per-cycle rows. Only do this AFTER 1-3 so
   the grown-cycle rows are correct; the inventory row is already bit-exact.

## ⚠ Mandatory regression gate (the lesson)
Wiring the fuel update into `grow_cycle!` changes `fire.cwd`/`fire.flive` at **fire time**, which
`fmburn!` reads — so it can move the validated fire results. **snt01 is ecounit 231Dd**, so the
already-committed FULIV2 override ALSO changes its `flive`; neither FULIV2 nor this chunk is exercised
by the current cycle-0/1 `test_snt01.jl` (the fire fires late). Before/after this chunk, run the full
snt01 stand-4 (and `fire_early`) to a fire cycle and diff the post-fire `.sum` (TPA/BA/mortality) vs
the Fortran baseline — the fire path must stay within its existing residual. Add that as a committed
fire-stand test (it is currently a coverage gap independent of this chunk).

## Validation target
The committed `carbon_jenkins.report.save` already has the 1995/2000 rows (DDW 2.5/…, Floor 6.6/…,
Shb 1.0/…, Total 129.9/…). Each chunk above is validated by another report column/cycle reconciling
bit-exact — same method as the inventory cycle. See [[fvsjl-ground-truth-binary-limits]] (the report
is a `.out` text report the stripped binary prints, so this is bit-exact-validatable).
