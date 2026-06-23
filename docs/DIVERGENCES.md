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
| MSB mortality `MSBMRT` | `MATUREW` (default QMDMSB=999) | C4 |
| size-cap (SIZCAP): DG bound + size-cap mortality + HT cap | `TREESZCP` — ✅ ported (kw_treeszcp!; sn/morts.f:692 floor + htgf.f:286 HT cap, test_treeszcp.jl). Residuals: cap mid-cycle TPA/BA carry the regen response to cap-driven mortality (QMD bit-exact, endpoint matches); htcap TopHt ≤4' declining-stand drift (TPA/BA/QMD bit-exact) | — |
| BAIMULT/HTGMULT/MORTMULT/REGHMULT/REGDMULT (`MULTS` 1/2/4/3/6) | the keyword — ✅ ported (active_multiplier; bit-exact vs Fortran, test_multipliers.jl) | — |
| FIXDG/FIXHTG | the keyword — ✅ ported (one-shot DG/HTG scaler + DBH window + tripled-record scaling; bit-exact, test_fix_scalers.jl) | — |
| HTGSTOP/TOPKILL | the keyword — ✅ ported (htgstp! + crown negative-ICR bypass; deterministic + stochastic bit-exact through firing cycle, test_htgstp.jl) | — |
| FIXMORT | the keyword — ✅ normal path ported (replace/add/max/mult + TPAMRT ordering; bit-exact, test_fixmort.jl); point/size concentration deferred | — |
| `FFERT` | the keyword (defaults inert) | — |
| event-monitor expression eval `EV*`/`ALGEVL` | `COMPUTE`/`IF` keywords (basic IF/THEN works; full expr partial) | — |
| stand-structure stage `SSTAGE/ISSTAG/KSSTAG/COVOLP` | structure-class output (not in .sum growth math) | — |
| fire `FM*` (~100), econ `EC*`, insects `MPB/DFB/TM/BWE`, mistletoe `MIS*` | FFE / econ / pest keywords | C7/C8 |
| alt-variant volume eqs `FVSBRUCEDEMARS/FVSSIERRALOG/FVSOLDGRO/HANNBARE…` | a non-SN volume-equation number | — |

**Limit of this audit:** the exercised path is verified empirically and the risky
in-scope routines individually; the gated rows are identified *structurally* (in the
oracle graph, absent in FVSjl), not yet branch-audited. Turning on each gate with a
scenario is how to convert "structurally diverges" into "validated" — that is the
test backlog implied by `docs/DECISION_FLOW_DETAILED.md`.

## Branch-level differential (`tools/branch_diff.js`)

The call-graph coverage above sees whole *routines*; it cannot see a branch missing
*inside* a routine that is otherwise ported (the class that hid the BAMAX gap). To
catch those, `branch_diff.js` fingerprints every distinctive numeric constant
(4-significant-figure value) in each C3/C4/C5 oracle routine and flags any whose
value appears nowhere in FVSjl (`src/*.jl` + the coefficient CSVs). A flagged
constant ⇒ a candidate un-ported formula/branch. After clearing false positives
(printf format specs, `f0`-suffixed integers, CSV precision), the real findings:

| finding | where | gate | status |
|---|---|---|---|
| **Fort Bragg special equations** — `dgf` dg5 growth (sp 8,13), `bratio` bark (sp 5,6,8,11,13), DG ATTEN override (sp 8=2056,13=689), AND KODFOR→81110 remap for volume | `dgf.f:515-537`, `bratio.f:106-118`, `dgf.f:636`, `forkod.f:137-140` | `IFOR==20` ⇔ forest code 701 (Fort Bragg → mapped to Uwharrie 81110) | ✅ **PORTED** — dgf/bratio/atten (s30 growth bit-exact) + KODFOR=81110 remap so VOLEQDEF gets region 8 (was zero-volume; s30 .sum now ±1 ulp). `test_fortbragg_coverage.jl` |
| REGENT establishment-mode crown draw (`0.89722−0.0000461·PCCF + 0.07985·ran`) | `regent.f:107-113` | establishment (`lestb`) — regen keywords | un-ported (C4 regen) |
| SDICLS relative-density quadratics (3 species-group eqns) | `sdical.f:320-345` | stand-stage classification (not the .sum SDI number) | partial (SDICLS) |
| `cfvol` merch-height conversion (8.3333) | `cfvol.f:37` | alternate cubic-volume equation (snt01 uses R8 Clark) | un-ported, vol-eq-gated |
| carbon factors (0.47/0.501/…) | `vols.f:157` | carbon output | C6/DBS |
| CRNMLT DBH bound default 99.0 | `crown.f` | `CRNMLT` keyword | gated, inert |

Everything else in the C3/C4/C5 core fingerprints clean. None of the above is on the
currently-validated path (consistent with snt01 being bit-exact).

## Test-coverage audit vs the decision flow (blindspots)

Auditing the test suite against the decision flow found, and closed, several
coverage blindspots — and in doing so corrected a wrong assumption:

1. **Physiography** — 133/135 scenarios used ecounit `231Dd`; only 3 of 11 DGF
   physiography branches were exercised. **Closed**: 8 numeric-index scenarios
   (s12–s19), all bit-exact at cyc0+cyc1 → every `dg_phys_*` column validated. No bug.
2. **Multi-cycle depth** — the suite only checked cyc0/cyc1. **Closed**:
   `test_multicycle.jl` regresses all 11 cycles vs a Fortran-oracle golden (BA kept
   tight as the BAMAX sentinel). Suite 187→627 assertions.
3. **Forest codes** — only 2 exercised. **Closed (partial)**: s20–s22 (forests
   802/806/809) bit-exact. Special mappings (905/908/701) still untested.
4. **★ Thinning/CUTS — WRONG ASSUMPTION CORRECTED.** It was marked "ported" from
   old FVSjulia-era notes, but **FVSjl has no thinning at all**: no `THIN*` keyword
   handler in `keyword_dispatch.jl`, no removal step in `grow_cycle!`. The snt01
   suite only validates the *unthinned* first stand; the `s11_thinbta` scenario
   silently no-ops in both engines (oracle had non-cutting params, FVSjl ignores the
   keyword), so it never tested anything. **CUTS / TREMOV / TREDEL / SDICLS-after-
   treatment / the removed-volume `.sum` columns are all un-ported.** This is the
   single biggest gap and is **upstream** in the cycle (GRINCR, before growth).

**Revised porting order (upstream→downstream):** **CUTS/thinning** (biggest, most
upstream, real scenarios exercise it) → Fort Bragg `IFOR==20` → REGENT establishment
→ SDICLS density → `cfvol`/carbon.
