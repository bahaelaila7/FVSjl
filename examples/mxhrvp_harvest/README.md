# MXHRVP — multistand harvest scheduling

The PPE **multistand treatment scheduling** feature (MXHRVP): instead of treating each
stand on its own schedule, you define a landscape harvest **policy** and the scheduler
picks *which* stands to cut each master cycle to meet a landscape **target**, in
**priority** order. Every `hv*.f` kernel is ported and validated bit-exact-or-cornered vs
the historical `FVSppe` oracle (see `docs/PORT_STATUS.md` → PPE MXHRVP).

```
mxhrvp_harvest/
  member_stand.key   member_stand.tre   ← one landscape member stand (EastCascades)
  run_harvest.jl                         ← the worked demo (run it)
```

## Run it

```bash
JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project=. examples/mxhrvp_harvest/run_harvest.jl
```

## The entry point

```julia
using FVSjl: ppe_run_landscape_harvest!, PPEStand, EastCascades
stands = [PPEStand("member_stand.key"; area = 11.0) for _ in 1:3]   # a 3-stand, 11-ac-each landscape
res = ppe_run_landscape_harvest!(stands; variant = EastCascades(),
          labels = ["ALL","ALL","ALL"],           # each stand's policy label set (SPLABEL)
          mslabel = "ALL",                          # MSPLABEL — the policy's candidate filter
          target_expr = "1000",                     # TARGET — the landscape resource flow
          priority_expr = "BBA",                    # PRIORITY — the per-stand ranking key
          credit_expr = "BBA",                      # CREDIT/UNITS — a stand's contribution if cut
          master_years = [1990, 2000, 2010])        # the master-cycle boundaries
```

`res` is one row per master cycle — the **HVSEL selection table**: `.selected` (which
stands were cut), `.priority` / `.credit_sel` per stand, `.target`, `.selected_resource`
(the total harvested toward the target), `.pct_of_target`, and `.hvpart` (the partial-cut
fraction of the one status-4 stand, if any). It mirrors what the Fortran `FVSppe` prints as
the per-policy selection table + the post-harvest **COMPOSITE YIELD STATISTICS**.

## What the demo shows

| # | policy | what it demonstrates |
|---|---|---|
| **(a)** | `TARGET=1000`, `PRIORITY=BBA`, `CREDIT=BBA` | a basic landscape target-flow schedule — stands cut in priority order toward a 1000-unit basal-area flow (selected resource 232 → 327 → 433 over the three cycles). |
| **(b)** | `CREDIT=SELECTED*BBA`, `TARGET=150`, `lprtct=true` | the **EXACT** partial cut — when only part of the marginal stand is needed, it is *partially* selected (status 4) and `hvpart` is the fraction (0.938 → 0.378 → 0.04 here). |
| **(c)** | `hvsel!(…; lhvmxc=true, hvmxcc=60)` | the **MXCLRCUT** max-contiguous-clearcut constraint — a would-be-clearcut stand is not cut if it would create more than `HVMXCC` contiguous clear-cut acres with the already-selected clearcut stands. Shown at the `hvsel!` kernel level (three would-be-clearcut stands in a star, 50 ac each, 60-ac cap → only the first is cut). |

> **Notes.** The policy expressions (`target_expr` / `priority_expr` / `credit_expr`) are
> evaluated by the ported Event Monitor over the landscape statistics (BBA, before-thin
> TPA/CuFt/BdFt, accretion, mortality, …) and the per-stand `SELECTED` variable — the same
> algebra as a `.key` `MSPOLICY … MSPLABEL/TARGET/PRIORITY/CREDIT … END` block. The
> selection outcomes and cyc-0 numerics are bit-exact vs the oracle; later-cycle credit
> values are cornered to the EastCascades growth residual (`#207`). (c) is shown at the
> kernel level because the coordinator-level `MXCLRCUT`/`DEFCCUT` clearcut-sign + neighbor
> wiring is a documented follow-up.
