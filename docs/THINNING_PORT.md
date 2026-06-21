# Thinning / CUTS port plan (FVSjl)

The audit found thinning is **unported** in FVSjl (no `THIN*` handler, no removal in
`grow_cycle!`). It is the top upstream gap. This is the goal-directed plan so the
implementation does not run in circles — target, algorithm, and integration points
are all captured below from the oracle (`base/cuts.f` / `FVSjulia/src/base/cuts.jl`).

## Validation target (committed: `test/harness/scenarios/s28_thindbh_cut.key`)
`THINDBH` at 2000 (date), DBHlo 0, DBHhi 999, cuteff 1, species 0(all), residTPA 0,
residBA 40 — i.e. thin all trees to residual basal area 40. Oracle result:

| cyc (yr) | TPA | BA | SDI | QMD | Tcuft | removed TPA |
|---|---|---|---|---|---|---|
| 1 (1995) | 507 | 103 | 202 | 6.1 | 1881 | 0 |
| 2 (2000) | 470 | 126 | 237 | 7.0 | 2481 | **334** (the thin) |
| 3 (2005) | 132 | 47  | 85  | 8.1 | 964  | 0 |
| 4 (2010) | 128 | 59  | 100 | 9.1 | 1243 | 0 |

FVSjl currently reproduces the *unthinned* baseline (439/147 at cyc3) → the gap.

## Algorithm (THINDBH = icflag 8, the lspecl + lbarea path)
`cuts.f` is a big multi-method router (THINDBH/BTA/ATA/BBA/ABA/CC/SDI/QFA/…). Start
with **THINDBH** (Fortran-baseline ground truth via snt01 block 2 + s28).

1. **Param setup (label_325→355):** valmin=DBHlo, valmax=DBHhi, cuteff, species,
   ctpa=residTPA, cba=residBA. `cstock = CLSSTK(jtyp, species, DBHlo, DBHhi, 0,999)`
   (jtyp=2 if cba>0 ⇒ basal-area class stock, else 1 ⇒ TPA). `rstock = ctpa+cba`,
   `remove = cstock − rstock`, `cuteff = min(1, remove/cstock)`. (`CLSSTK` = sum of
   the chosen stock over trees in the DBH/HT/species class — needs porting.)
2. **Priority (label_700):** for THINDBH `lspecl=true` ⇒ `WK2[i]=1` for trees with
   DBH∈[valmin,valmax) and species included and not `LEAVESP`, else 0. (No sort for
   lspecl — trees processed in record order; `RDPSRT` sort is only for non-lspecl.)
3. **Trial removal loop (label_900→1100):** walk records; for each eligible tree
   `prem = WK4[it]*cuteff` (WK4 = pre-thin PROB copy), `cut_v = prem·d²·0.005454154`
   (lbarea) or `prem` (TPA). Accumulate `totcut`; when `remove−totcut−cut_v < 0` the
   tree is the **last** one and `prem` is scaled so the cut hits `remove` exactly.
   Reduce the tree's PROB by `prem` (→ `TREMOV`), accrue removed volume.
4. **TREMOV/TREDEL:** subtract `prem` from PROB; drop records that reach 0 (`TREDEL`)
   and compact. Accumulate `ONTREM`/removed-volume for the `.sum` removed columns.
5. **Post-thin:** recompute density (`DENSE`) and after-treatment SDI/stats
   (`SDICAL`/`SDICLS` → the ATBA/ATSDIX/ATTPA `.sum` "after treatment" fields).

## Cycle-order integration (the tricky part)
In the oracle, `CUTS` runs in **GRINCR at the start of the cycle, before growth**,
and the cycle's `.sum` row reports the removed columns for that cycle. In FVSjl the
host loop is `summary_row` then `grow_cycle!`. So:
- Add a `cuts!(s)` step at the **top of `grow_cycle!`** (after `compute_density!`,
  before `diameter_growth!`), running any thin scheduled for the current year.
- The removed-volume columns belong to the **current** cycle's row, so `summary_row`
  must read them after `cuts!` — i.e. either run `cuts!` before `summary_row` for the
  cycle, or stash removed totals on the state for the next `summary_row`. Match the
  oracle: the 2000 row shows removed=334 with the *pre-thin* live stats (470/126),
  and the 2005 row shows the *post-thin* grown stand. So removed is reported on the
  thin cycle's row alongside pre-thin live stats.
- **Tripling interaction:** by cyc2 the records are already tripled (ITRN≈243); CUTS
  operates on the tripled set. Port must cut tripled records (per-record PROB), which
  is why uniform-proportional ≠ exact (the 334 vs naive-321 gap is the record-order
  trial-loop + tripling). Validate against the oracle's per-record cut, not a formula.

## State changes
- A `ScheduledActivity` list (date + icflag + 6 params) on `Control` (or a new
  `schedule` field), populated by the `THIN*` keyword handlers in
  `keyword_dispatch.jl`. Update the `Control`/`StandState` constructors.
- `WK4`-equivalent pre-thin PROB scratch + removed-volume accumulators on the state.

## Milestones (each validated before the next)
1. **Scheduling**: parse `THINDBH` (+ date) into a `ScheduledActivity`; unit-test the
   parse. (No behaviour change; suite stays green.)
2. **CLSSTK + cuts! for THINDBH/lbarea**: implement class-stock + the trial loop +
   TREMOV; wire `cuts!` into `grow_cycle!`. Validate `s28` cyc3 (132/47/8.1/964) and
   removed=334 vs the oracle, then wire `s28` into the multi-cycle test.
3. **Removed-volume `.sum` columns** + post-thin after-treatment stats.
4. **Other methods**: THINBTA/ATA (CLSSTK from-below/TPA), THINBBA/ABA, THINDBH
   multi-line DBH-class schedule (snt01 block 2 → match the Fortran baseline `.sum`),
   then THINCC/THINSDI/THINQFA.
5. Re-validate the whole C10 matrix (the thinning scenarios become live).

## Helpers to port (dependencies)
`CLSSTK` (class stocking), `RDPSRT` (priority sort, already used elsewhere?),
`TREMOV`/`TREDEL` (record removal+compaction), and for later methods `CUTQFA`,
`SDICLS`, `CCCLS`, `RDCLS`. CLSSTK + TREMOV/TREDEL are the immediate needs.
