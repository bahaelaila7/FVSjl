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

### Board feet derived from the shared sawtimber R8 Clark call (not a separate BF taper)
- **Where (FVSjl):** `engine/volume.jl:compute_volumes!` (`bf = v[10]` from the single per-tree `_R8CLARK_VOL` call)
- **Where (Fortran):** `base/vols.f` / `base/fvsvol.f` — board feet uses its own BFMIND/BFTOPD/BFSTMP merch limits
- **Category:** modelling-simplification (architecture)
- **Gate:** always (no keyword default exercises it; snt01 + all default scenarios are bit-exact)
- **Description / evidence:** FVSjl computes all four per-tree volumes (total/merch/sawtimber
  cubic + board feet) from ONE R8 Clark taper call parameterised by the *cubic/sawtimber*
  merch standards, reading board feet as `v[10]`. By default this is exact because the SN
  merch table has `scf_top_dib == bf_top_dib`, `scf_min_dbh == bf_min_dbh`, and
  `scf_stump == bf_stump` for every species (sawtimber and board-foot standards coincide).
  The **VOLUME** keyword's merch-top/stump params (`TOPD/SCFTOPD/SCFSTMP`) and the entire
  **BFVOLUME** keyword can break that coincidence; when they do, FVSjl's board feet ride the
  changed sawtimber call while Fortran's board feet stay on the (separate) BF standards, so
  the board-foot `.sum` column diverges (observed up to ~7% on a `SCFTOPD 7→9` override).
  The **VOLUME `DBHMIN` gate is exact** (it only gates merch cubic, never the taper call) and
  is validated bit-exact vs Fortran (`test_volume_override.jl`). BFVOLUME currently only
  affects topkilled-tree board feet (via `bftopk`, which does read the per-stand BF standards).
- **Full-fix roadmap (investigated, deferred):** Fortran's structure is `fvsvol.f:257` — it sets
  `BFPFLG=1` when the BF standards equal the sawtimber standards (the default → board feet from
  the one sawtimber call), else `BFPFLG=0` and a **second** `VOLINITNVB` call is made with
  `MTOPP=BFTOPD, STUMP=BFSTMP, PROD='01'` to get board feet (`fvsvol.f:362-503`). A prototype
  second `_R8CLARK_VOL` call made board feet match — **but** it surfaced a further coupling:
  the BF call's Region-8 "≥10 ft of product" rule (`fvsvol.f:499`, `HT1PRD<10 ⇒ TVOL(4)=0`)
  zeros the **sawtimber cubic** too, so a BFTOPD override also drops the merch-cubic / sawtimber
  `.sum` columns. Replicating that needs the un-zeroed BF-top log height (FVSjl's `_R8CLARK_VOL`
  already zeros `sawHt` at a slightly different `merchL+stump+trim`≈9.5 threshold), so the full
  port must thread `HT1PRD` out of the taper model and apply the 10-ft rule across both products.
- **Status:** open — needs the BFPFLG=0 second call **plus** the Region-8 10-ft product coupling; deferred.

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
- **Gate:** always, but only visible from the first UNTRIPLED grow (snt01 cycle 2,
  the 2005 .sum row; tripling covers cycles 0–1).
- **Evidence (measured 2026-06-23, base.trl per-record vs FVSjl):** every linear
  aggregate is bit-exact EVERY cycle — TPA/BA/SDI/QMD/TopHt match to the .sum print.
  The only residual is the nonlinear cubic-volume sum: **±1 cuft on ~3026 (≈0.03%)**,
  first at 2005, slowly compounding (cuft ±1, bdft ±4–18 by 2040). NOT ~1–2% (that was
  a stale pessimistic estimate). Reconstructed the Fortran stand cuft from base.trl
  (243 records × TPA·TotCuFt = 3026.3 ✓) and compared per-record to FVSjl: per-record
  HEIGHTS agree exactly and per-record DBH agrees **to the .trl's 0.1" print precision**
  (max|Δdbh| = 0.050" = exactly half the 0.1" quantum, uniform across all 5 species;
  mean 0.027" = the rounding mean). So the per-record growth matches at print precision;
  the ±1 cuft is a *sub-0.1"-per-record* effect in the nonlinear volume, below .trl
  resolution. **To diagnose further requires SOURCE-LEVEL high-precision Fortran dumps**
  (a custom WRITE of per-record DBH/DG/OLDRN in dgscor.f / the cycle-2 tree loop) to see
  whether it is a closable within-species draw-assignment/serial-correlation ordering
  difference for near-identical records, or irreducible Float32 propagation. At ≈0.03%
  with per-record agreement already at print precision, this is at/near the Float32 floor.
- **Resolution (2026-06-23): IRREDUCIBLE NUMERIC DRIFT (category 2), not a bug.** Verified
  `dgscor!` line-by-line against base/dgscor.f — it is a bit-exact transliteration (BACHLO·
  RHOCP + RHO·OLDRN, the |FRM|>DGSD·SSIG bounded REDRAW, the DDS>4/>5 taper, OLDRN=FRM,
  EXP(FRM)). So the DGSCOR kernel arithmetic is NOT the source. The residual is upstream
  transcendental ulp (exp/log/sqrt in DGF's WK2, autcor's SSIG/RHO) propagating through the
  **bounded redraw**: a sub-ulp OLDRN difference can flip whether a given tree re-draws, which
  reshuffles the per-record assignment of an otherwise-identical draw multiset (aggregate BA
  averages out, nonlinear cuft drifts ≈0.03%) and the AR(1) carry compounds it slowly. Closing
  it would require bit-identical transcendental libs Julia↔Fortran — out of scope. ⇒ snt01 is
  effectively bit-exact: every linear aggregate exact every cycle, only provably-irreducible
  Float32 noise in the nonlinear volume. Same class as the "single off-by-1 BdFt" drift.

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
| CRNMLT DBH bound default 99.0 | `crown.f` | `CRNMLT` keyword | ✅ ported (active_crn_mult; the multiplier IS applied — the "inert" note referred only to the bound default, not the multiplier) |

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
