# Divergences from Fortran FVS

This file records **every** place where FVSjl knowingly differs from the reference
Fortran FVS, in two categories:

1. **Suspected FVS bugs** — behaviour we believe is wrong given the code's intent or
   comments. These are fixed only behind the non-default `faithful=false` path and
   listed here so they can be raised with the FVS maintainers for confirmation.
2. **Unavoidable numerical drift** — e.g. a transcendental function whose last ulp
   differs between the Julia and Fortran math libraries. Documented, not chased.

In `faithful=true` mode (the default, and what the test suite runs) FVSjl is
intended to be bit-exact to Fortran, so this list has no behavioural effect there.

Format per entry:

> ### <short title>
> - **Where (FVSjl):** file:func
> - **Where (Fortran):** file.f:line
> - **Category:** suspected-bug | numeric-drift
> - **Gate:** `faithful=false` (bug fixes only) | always (drift)
> - **Description / evidence:**
> - **Status:** open | reported | confirmed

---

_No behavioural (faithful-mode) divergences recorded yet — snt01.sum is bit-exact._

---

# Semantic coverage vs the FVS call graph

A different question from the behavioural drift above: does FVSjl **implement every
semantic** in the FVS call graph (`docs/decision_flow.html`, extracted from the
faithful oracle)? Audit method:

1. **Exercised path = empirically verified.** snt01.sum is bit-exact and the C10
   matrix is bit-exact at cyc0/1 and BA-exact every cycle. Any divergence on a
   routine the scenarios actually run would surface as a numeric diff, so the
   exercised semantics are covered by construction.
2. **Risky in-scope routines spot-checked** (routines reachable & not extension-
   gated, where a silent gap could still hide): `DGBND` large-tree growth bound —
   **covered** (`bark_and_bounds.jl:dg_bound`, applied in both the REGENT and DGDRIV
   paths, diameter_growth.jl:465-482); `DAMPRO`/`BASDAM` damage — **covered**
   (input-time in INITRE; snt01's top-killed tree is bit-exact); `FORTYP`/`SILFTY`,
   DG calibration, `MAICAL` — all present.

**Result: no divergence in the semantics any current scenario exercises**, with one
exception (item 1 below). Everything else that is missing is a *conditionally-gated*
branch no test scenario triggers — structurally divergent, not yet validated.

### 1. DGSCOR per-record serial-correlation (the one *live* divergence)
- **Where:** `variants/southern/diameter_growth.jl:dgscor!` vs `base/dgscor.f`.
- **Gate:** always, but only visible in untripled cycles (cyc ≥ 3).
- **Evidence:** TPA/BA/SDI/QMD/TopHt stay bit-exact; only the nonlinear cubic-volume
  sum drifts ~1–2% because the per-record growth draw distributes differently while
  the aggregate matches. Needs exact per-record RNG-order matching to close.

### Conditionally-gated branches absent from FVSjl (diverge only if triggered)
Each is present in the oracle graph but not (yet) in FVSjl; the **gate** is what a
scenario must set to make it execute and become Fortran-validatable.

| missing semantics (FVS routines) | gate that triggers it | chunk |
|---|---|---|
| record compression `COMPRS/COMCUP/CMRANG/MEANSD` | `COMPRESS` keyword or >~3000 records (tests sit ~243) | — |
| establishment/regen `ESTAB/ESIN/ESUCKR/ESPLT*/…` (~20) | `NATURAL`/`PLANT`/`SPROUT`/AUTOES (tests are NOAUTOES) | C4 |
| MSB mortality `MSBMRT`, size-cap mortality | `MATUREW` / `SIZECAP` (defaults QMDMSB=999, SIZCAP=999) | C4 |
| `FFERT`, `HTGSTP` (HTGSTOP/TOPKILL), `MULTS`/MORTMULT/FIXMORT/FIXDG/FIXHTG | the respective keyword (defaults inert) | — |
| event-monitor expression eval `EV*`/`ALGEVL` | `COMPUTE`/`IF` keywords (basic IF/THEN works; full expr partial) | — |
| stand-structure stage `SSTAGE/ISSTAG/KSSTAG/COVOLP` | structure-class output (not in .sum growth math) | — |
| fire `FM*` (~100), econ `EC*`, insects `MPB/DFB/TM/BWE`, mistletoe `MIS*` | FFE / econ / pest keywords | C7/C8 |
| alt-variant volume eqs `FVSBRUCEDEMARS/FVSSIERRALOG/FVSOLDGRO/HANNBARE…` | a non-SN volume-equation number | — |

**Limit of this audit:** the exercised path is verified empirically and the risky
in-scope routines individually; the gated rows are identified *structurally* (in the
oracle graph, absent in FVSjl), not yet branch-audited. Turning on each gate with a
scenario is how to convert "structurally diverges" into "validated" — that is the
test backlog implied by `docs/DECISION_FLOW_DETAILED.md`.
