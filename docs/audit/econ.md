# Bandaid audit — Economic analysis (ECON extension)

**Module file:** `/workspace/FVSjl/src/engine/econ.jl`
**FVS source checked:** `econ/eccalc.f`, `econ/echarv.f`, `econ/ecin.f`, `econ/ecstatus.f`,
`econ/ecinit.f`, `common/ECNCOM.F77` (build copy `bin/FVSsn_buildDir/ECNCOM.F77`).
**Adjacent wiring read:** `src/core/state.jl`, `src/engine/keyword_dispatch.jl`,
`src/engine/simulate.jl`, `test/unit/test_econ.jl`.

## Summary

The self-contained financial kernels are faithful to `eccalc.f`:

- `econ_present_value` ↔ `computePV` (eccalc.f:991-1002): `amt<=0 → 0`, else `amt/(1+rate)^time`. Exact.
- `econ_pnv` ↔ `computePNV` (eccalc.f:1006-1021): cost@`i-1`, rev@`i`, loop `i=1:n`. Exact.
- `econ_bc_ratio` ↔ eccalc.f:683 (`revDisc/costDisc`). Core exact (guard threshold flagged below).
- `econ_rate_of_return` ↔ eccalc.f:686-687. Exact.
- `econ_sev` ↔ the perpetual-rotation factor in `calcAppreSev` (eccalc.f:1102-1104): with no
  appreciation, `calcAppreSev = undiscAmt·factor/(factor−1)/(1+rate)^discTime`,
  i.e. `PV(event)·(1+rate)^endTime/((1+rate)^endTime−1)` — matches `net_rotation·f/(f−1)`.
- `econ_forest_value` ↔ eccalc.f:649-655: `discSev=PV(sevInput,endTime)`, `fv=pnv+discSev`,
  `reprod=fv−sevInput`. Exact.
- ECON unit constants (TPA=1/BF_1000=2/FT3_100=3) match ECNCOM.F77 (param block ~line 21; the
  `:19` cite is a couple lines off but the values are correct).
- `harvest_value` DBH band `dbh_lo <= dbh < dbh_hi` and `sp==0 ⇒ all` match echarv.f:55,63 (`>= lo
  .and. < hi`); the `/1000` and `/100` unit conversions match valueHarvest (eccalc.f:284/290) and
  echarv.f:97/129 (raw volume × price, conversion applied to one factor).

Four concerns flagged below.

---

## FLAG 1 — BANDAID: `econ_stand_pnv` under-discounts every harvest by one year

- **jl symbol/line:** `econ_stand_pnv`, `src/engine/econ.jl:168-172`
  ```julia
  for (yr, cost, rev) in ec.harvests
      t = Int(yr) - base
      disc_cost += econ_present_value(cost, max(0, t - 1), rate)
      disc_rev  += econ_present_value(rev,  max(0, t),     rate)
  end
  ```
  with `base = ec.base_year`, and `base_year` set to the **first econ-active cycle year**
  (simulate.jl:212), harvests recorded with `yr = current_cycle_year` (the cycle-start year).
- **Claim:** docstring (econ.jl:152) says "Present net value … (eccalc.f)", i.e. claims faithfulness
  to the eccalc.f discounting convention.
- **FVS source checked:** eccalc.f:101 `beginTime = beginAnalYear - startYear + 1`; a cycle's
  harvest is booked into `undiscCost(beginTime)/undiscRev(beginTime)` (eccalc.f:328-329); calcEcon
  discounts `costDisc += computePV(undiscCost(i), i-1)` / `revDisc += computePV(undiscRev(i), i)`
  (eccalc.f:628-629). So a harvest in calendar year `Y` has analysis index `i = Y - startYear + 1`,
  giving **cost discounted `Y-startYear` years and revenue `Y-startYear+1` years**. For the very
  first econ year (`Y = startYear = beginAnalYear`), `beginTime=1`: cost@0, **revenue@1 year**.
- **What jl does instead:** with `base = startYear`, jl uses cost@`max(0, (Y-base)-1)` and
  rev@`Y-base`. That is **one year short of FVS for both streams** (cost `t-1` vs FVS `t`; rev `t`
  vs FVS `t+1`), and at the base year it collapses the cost-at-start / revenue-at-end intra-year
  split (FVS = 0 / 1; jl = 0 / 0).
- **Why it reads as a bandaid:** the only check is `test/unit/test_econ.jl:155`
  `@test pnv.disc_rev ≈ econ_present_value(rev, 2000 - 1990, 0.04f0)` — it asserts jl's own
  `PV(rev, 10)` convention. The faithful eccalc.f value for a year-2000 harvest with start 1990 is
  `PV(rev, 11)`. The test validates the formula against itself, not against FVS.
- **Faithfulness impact:** every harvest revenue (and cost beyond the base year) is discounted one
  year too little, biasing PNV/B-C/RRR upward. Blast radius is currently limited: `econ_stand_pnv`
  is only exercised by the unit test, not wired into any `.out`/DBS writer yet — but it is the
  intended stand-level aggregator and contradicts the cited source.

---

## FLAG 2 — GAP: log-graded revenue units silently valued at zero

- **jl symbol/line:** `harvest_value`, `src/engine/econ.jl:104-106`
  ```julia
  vol = r.unit == ECON_TPA     ? tpa :
        r.unit == ECON_BF_1000 ? bdft * tpa / 1000f0 :
        r.unit == ECON_FT3_100 ? cuft * tpa / 100f0  : 0f0
  ```
  Comment (econ.jl:96-97): "Log-graded units are handled by the log-bucking layer (TODO)."
- **FVS source checked:** ECNCOM.F77 defines `BF_1000_LOG = 4, FT3_100_LOG = 5`. echarv.f:72-91
  (BF) and 105-123 (FT3) bucket per-log volume into `revVolume(spId,BF_1000_LOG/FT3_100_LOG,·)`,
  and valueHarvest (eccalc.f:305-326, with `MAX_REV_UNITS=5`) values them like the whole-tree
  units. So units 4/5 are real revenue, not a no-op.
- **Faithfulness impact:** an HRVRVN keyword that specifies log-graded units would contribute **zero**
  revenue in FVSjl. Documented as TODO and outside the sn.key path, so a true GAP rather than a
  silent bug, but it must be closed before any log-bucked keyfile is trusted.

---

## FLAG 3 — GAP: discount rate is a hardcoded 4%, never read from input

- **jl symbol/line:** `EconState` default `discount_rate = 0.04f0` (`src/core/state.jl:639`),
  consumed by every econ.jl kernel via `ec.discount_rate` (e.g. econ_stand_pnv:157). `kw_econ!`
  (keyword_dispatch.jl:1261-1290) parses ANNUCST/HRVVRCST/HRVRVN but **never sets a rate**.
- **FVS source checked:** ecinit.f:15 `discountRate = 0.0` (default), and the live value is set
  from the **STRTECON** event-monitor activity: ecstatus.f:49 `discountRate = strtParms(1)` /
  eccalc.f:151, then eccalc.f:91/153 `rate = discountRate/100.0`. The rate is *not* a field on the
  ECON/ANNUCST keyword lines (ecin.f has no `discountRate=` assignment).
- **Faithfulness impact:** because STRTECON is unported, `ec.discount_rate` is permanently the
  magic `0.04`. Every PNV/SEV/forest-value the (otherwise faithful) kernels produce uses 4%
  regardless of the keyfile's specified rate. The discounting math is correct; its single most
  important input is not wired to FVS. GAP (the `0.04` constant has no FVS basis — FVS defaults 0.0).

---

## FLAG 4 — GAP (low): B-C ratio / RRR guards drop FVS's `NEAR_ZERO` threshold

- **jl symbol/line:** `econ_bc_ratio` (econ.jl:43-44, guard `disc_cost > 0f0`) and
  `econ_rate_of_return` (econ.jl:52-55, guard `disc_cost > 0 && disc_rev > 0`).
- **FVS source checked:** eccalc.f:58 `real, parameter :: NEAR_ZERO = 0.01`; eccalc.f:681
  `if (costDisc > NEAR_ZERO)` gates bcRatio, eccalc.f:685 `if (revDisc > NEAR_ZERO)` gates rrr.
- **Faithfulness impact:** for a discounted cost or revenue in `(0, 0.01]`, FVS leaves the metric
  blank/uncalculated, whereas FVSjl returns a (potentially very large) ratio / rate. Negligible for
  realistic stands, but it is a defined source threshold the port silently omits.

---

## Decisions reviewed but found faithful (not individually flagged)

computePV, computePNV, bc-ratio core, rrr core, SEV factor, forest/reprod value, unit constants,
DBH half-open band + species-0 wildcard, the `/1000`–`/100` unit conversions,
`econ_value_harvest`/`record_harvest!` accumulation, and the `econ_stand_pnv` annual-cost loop
(`t=0…N-1 ⇒ PV(annCost, 0…N-1)` matches eccalc.f `i=1…endTime ⇒ i-1`).
