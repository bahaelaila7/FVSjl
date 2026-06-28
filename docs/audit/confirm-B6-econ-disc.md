# Confirm B6 — ECON harvest discounting off by one year vs eccalc.f convention

**Verdict: CONFIRMED** (semantic mismatch; both harvest cost and harvest revenue are
discounted exactly one year too little relative to the FVS `eccalc.f` convention).

- **FVSjl path:** `src/engine/econ.jl:156-174` (`econ_stand_pnv`), the harvest loop at
  `econ.jl:168-172`; flow plumbing at `src/engine/simulate.jl:210-215`.
- **FVS source:** `econ/eccalc.f` (build copy `bin/FVSsn_buildDir/eccalc.f`), lines
  101-102, 328-329, 380-387, 624-634, 991-1002.

---

## 1. What the Julia computes (econ.jl:156-172)

```julia
function econ_stand_pnv(ec::EconState, end_year::Integer)
    rate = ec.discount_rate
    base = Int(ec.base_year)
    ...
    # annual management cost: time t = 0 … end-base-1
    if ec.ann_cost > 0f0
        for t in 0:(end_year - base - 1)
            disc_cost += econ_present_value(ec.ann_cost, t, rate)        # ann cost @ t
        end
    end
    # harvest cash flows
    for (yr, cost, rev) in ec.harvests
        t = Int(yr) - base
        disc_cost += econ_present_value(cost, max(0, t - 1), rate)       # harvest cost @ t-1
        disc_rev  += econ_present_value(rev,  max(0, t),     rate)       # harvest rev  @ t
    end
    ...
```

`base = ec.base_year`, set once to the **first econ-active cycle year** (the cycle-start
year, `current_cycle_year`) at `simulate.jl:212`; each harvest is pushed with
`yr = current_cycle_year` (the cycle-start year) at `simulate.jl:214`.
`econ_present_value(amt, time, rate) = amt/(1+rate)^time` (econ.jl:20-21) — exact match to
`computePV` (eccalc.f:991-1002).

So for a harvest in calendar year `Y` with `t = Y - base`:
- **harvest cost** discounted `max(0, t-1)` years,
- **harvest revenue** discounted `max(0, t)` years.

## 2. What the FVS Fortran computes (eccalc.f)

- `beginTime = beginAnalYear - startYear + 1` (eccalc.f:101). `startYear = econStartYear`
  (eccalc.f:88). Time index is **1-based**: the start year maps to index 1.
- A cycle's harvest is booked into the year-indexed accumulators at `beginTime`
  (the cycle's first year, where harvests always occur):
  `undiscCost(beginTime) += harvCst` (eccalc.f:328),
  `undiscRev(beginTime)  += harvRvn` (eccalc.f:329).
- `calcEcon` discounts each year index `i` (eccalc.f:625-629):
  ```
  costDisc += computePV(undiscCost(i), i-1, rate)   ! costs accrue at BEGINNING of year
  revDisc  += computePV(undiscRev(i),  i,   rate)   ! revenues accrue at END   of year
  ```

For a harvest in calendar year `Y` the analysis index is `i = Y - startYear + 1`, so:
- **harvest cost** discounted `i-1 = Y - startYear` years,
- **harvest revenue** discounted `i = Y - startYear + 1` years.

In 0-based offset `t = Y - startYear` (≡ Julia's `t` when `base = startYear`):
FVS gives **cost @ t**, **rev @ t+1**.

## 3. The precise semantic divergence

| stream        | FVS exponent          | Julia exponent        | error      |
|---------------|-----------------------|-----------------------|------------|
| harvest cost  | `t`   (`= i-1`)       | `max(0, t-1)`         | one short  |
| harvest rev   | `t+1` (`= i`)         | `max(0, t)`           | one short  |

Both harvest streams are under-discounted by exactly one year → present values too high
(less discounting), inflating PNV. At the base year (`t=0`) FVS splits cost@0 / rev@1
(cost at start, revenue at end of the year); Julia collapses this to cost@0 / rev@0,
losing the intra-year start/end distinction that the eccalc.f comments call out
(eccalc.f:628-629, "Costs accrue at beginning of year" / "Revenues accrue at end").

**Internal-inconsistency corroboration (static, no test):** the *annual* cost in the very
same function discounts year `Y` at exponent `t = Y-base` (econ.jl:163-164), faithful to
FVS `computePV(undiscCost(j), j-1)` since FVS index `j = (Y-startYear)+1` ⇒ exponent
`j-1 = Y-startYear = t`. A *harvest* cost in the same year `Y` is, per FVS, equally a
"cost at beginning of year" and must use the identical exponent `t`. Julia instead uses
`t-1` for the harvest cost — so two cash flows that FVS treats with one convention are
discounted differently. This contradiction is internal to the Julia and does not depend
on any FVS run.

## 4. Faithful fix (do NOT apply)

```julia
for (yr, cost, rev) in ec.harvests
    t = Int(yr) - base                                   # = i-1, the FVS 0-based offset
    disc_cost += econ_present_value(cost, t,     rate)   # eccalc.f:628  cost @ i-1
    disc_rev  += econ_present_value(rev,  t + 1, rate)   # eccalc.f:629  rev  @ i
end
```

`yr >= base` always (harvests occur on/after the first econ cycle), so `t >= 0` and the
`max(0, …)` guards become unnecessary; `t+1 >= 1`.

**Base-alignment caveat (separate concern):** this fix assumes `ec.base_year == startYear`
(FVS `econStartYear`). Julia sets `base_year` to the first *econ-active cycle year*
(simulate.jl:212), which equals `econStartYear` only when ECON activates at the simulation
start and is not re-initialized by `StrtEcon`/`ECON_START_YEAR` (eccalc.f:127-139). If ECON
can start later, the `base` definition is a second, independent divergence and must be
pinned to `econStartYear` for the discounting to match; it is not part of this off-by-one
but would compound it. Do not let a `base_year = startYear - 1` style tweak be used to
"absorb" the off-by-one — that would be a second bandaid masking this one.

## 5. Upstream rank: **LEAF**

`econ_stand_pnv` outputs `pnv / disc_cost / disc_rev` consumed only by the ECON report
tables and the Event-Monitor registers `EV_PNV / EV_DISCCOST / EV_DISCREVN`
(eccalc.f:630-647). They do **not** feed growth, mortality, density, regen, or any tree
state, so a fix changes only the ECON report/event values, not the simulated stand.
(Could rise to MID only if a keyword file branches activity scheduling on a PNV/DiscCost
event variable; not the case in the SN econ scenarios.)

## 6. Reachability

- **Unit:** `test/unit/test_econ.jl:138-156` exercises this exact path but asserts the
  Julia's *own* convention — `@test pnv.disc_rev ≈ econ_present_value(rev, 2000-1990, 0.04)`
  (i.e. `PV(rev, 10)`), whereas the faithful eccalc.f value for a 2000 harvest with start
  1990 is `PV(rev, 11)`. The test validates the formula against itself; it will **not**
  move when the bug is present and **will need updating** (to `PV(rev, 11)` / `PV(cost,10)`)
  when the fix lands — its current passing is not evidence of faithfulness.
- **Integration:** the sn.key ECON scenario (per project memory, "both econ tables match
  Fortran") touches the path. Any apparent table match must be treated as non-evidence per
  method; if it genuinely matches it is most plausibly because the scenario's harvest sits
  at/near the base year and the integer-rounded report column absorbs the one-year
  discount, or because a compensating `base_year` offset hides it — both worth a live
  differential to settle (see below).

**Net:** fixing this is a silent correctness change to the ECON tables. It will require the
self-referential unit assertion to be re-pinned to the eccalc.f exponents, and should be
checked against an instrumented `eccalc.f` run for a harvest year strictly greater than
`econStartYear` (where `t > 0` makes the divergence unmaskable).

## 7. Masked-bug watch

If, after this fix, the sn.key ECON tables *regress* away from live FVS, that is a signal
the **base_year definition** (Section 4 caveat) is a second, compensating divergence —
not a reason to revert. The correct end state is: `base = econStartYear` AND harvest
`cost @ t`, `rev @ t+1`. Both must be faithful simultaneously.
