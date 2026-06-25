# Faithful artifacts — code kept verbatim for bit-exactness, not because it is the best idiom

FVSjl is an idiomatic Julia rewrite, but a handful of places deliberately reproduce the **exact
implementation** of the Fortran FVS rather than the canonical CS/math idiom. They are kept ONLY because
FVS's downstream numerics are sensitive to their exact behavior, so modernizing them would break
byte-for-byte agreement with the oracle — even though the modern version would be cleaner and just as
*correct* in the forest-science sense.

This file catalogs every such artifact: what it is, the canonical alternative, why it's kept, and
**whether it is bit-validated** (we KNOW it matters) or merely conservative.

Three categories of "non-canonical":
- **FAITHFUL ARTIFACT** — verbatim port; changing it demonstrably breaks bit-exactness vs the oracle.
  Keep until/unless the bit-exact-vs-Fortran goal is dropped.
- **IRREDUCIBLE** — cannot be made bit-exact by any choice (transcendental last-ULP); not "kept", just
  inherent. Documented in `DIVERGENCES.md`.
- **BEST-EFFORT STAND-IN** — a modern approximation used where the oracle is unavailable (binary-blocked);
  NOT validated against FVS and may not even match it. A candidate to revisit, not a faithful artifact.

---

## 1. RNG — Park–Miller multiplicative LCG  ·  FAITHFUL ARTIFACT
`src/core/rng.jl` (`rann!`, `esrann!`, `ranseed!`, `rannput!`)

- **What:** a Park–Miller minimal-standard LCG, multiplier `16807`, modulus `2³¹−1`, two independent
  streams (main + establishment, seeded 55329). Each `rann!` is one LCG step returning a Float32 in (0,1).
- **Canonical:** Julia's `Random` (Xoshiro256++ / MersenneTwister) with `rand()`.
- **Why kept:** every stochastic result (diameter-growth serial correlation, mortality kill distribution,
  tripling perturbation, establishment heights) consumes this stream in a specific order; matching FVS
  bit-for-bit requires the identical LCG sequence. A modern PRNG would give statistically-equivalent but
  byte-different results.
- **Validated:** yes — the whole snt01 growth/mortality trajectory is bit-exact only with this stream.
- **Cost of keeping:** the LCG is a known weak generator (short period 2³¹, serial correlation); fine for
  this purpose because we are *reproducing* FVS, not generating high-quality randomness.

## 2. Normal draw — BACHLO composite rejection (Tocher 1963)  ·  FAITHFUL ARTIFACT
`src/core/rng.jl` (`bachlo`)

- **What:** draws `N(xbar, stdev)` by Batchelor/Tocher's composite-rejection method — a specific
  accept/reject loop over uniforms from the LCG.
- **Canonical:** `randn()` (Ziggurat) scaled and shifted.
- **Why kept:** it consumes a *specific number* of LCG uniforms per draw and has specific rejection
  bounds, so it both (a) produces the exact deviate and (b) advances the shared stream by the exact
  amount — every later draw depends on that alignment. `randn` would desynchronize the entire stream.
- **Validated:** yes (coupled to #1).

## 3. Sort — `_rdpsrt!` Singleton/Sedgewick indirect quicksort (with `@goto`)  ·  FAITHFUL ARTIFACT
`src/engine/cuts.jl` (`_rdpsrt!`) — the ONLY `@goto`/`@label` left in the codebase.

- **What:** a verbatim port of RDPSRT (rdpsrt.f): middle-pivot Hoare-partition quicksort with an explicit
  64-deep stack, sorting an index permutation, descending. Uses `@goto` to mirror the Fortran control flow.
- **Canonical:** `sortperm(key; rev=true)` or `sort!`.
- **Why kept:** quicksort is **unstable** and this one's tie-order among equal keys is an *artifact* of the
  partition mechanics (no semantic meaning). But FVS's downstream IS sensitive to it: which of two
  equal-DBH lineages a from-below cut removes, and the post-thin physical record order that seeds the
  per-record DGSCOR/serial-correlation RNG traversal and tripling. A stable sort or a different quicksort
  would be equally *correct* but byte-different. (See the long comment at the function head.)
- **Validated:** yes — snt01 stands 2/3 post-thin trajectory depends on the exact tie-order.
- **Note:** this is the one place we reproduce an implementation *accident*, not a model semantic.

## 4. Serial correlation — AUTCOR ARMA(1,1) + DGSCOR bounded-normal recurrence  ·  FAITHFUL ARTIFACT
`src/variants/southern/serial_correlation.jl`

- **What:** `autcor` computes ARMA(1,1) variance/covariance multipliers; `dgscor` perturbs each tree's
  `ln(DDS)` by `frm = exp(BACHLO(0,ssigma)·rhocp + rho·OLDRN_prev)`, a per-tree auto-correlated error with
  a **bounded redraw** (reject-and-redraw until the deviate is in range).
- **Canonical:** a standard AR(1)/ARMA simulation from a library, or vectorized `randn`-based noise.
- **Why kept:** the exact draw order, the bounded-redraw loop (which consumes a variable number of LCG
  uniforms), and the `OLDRN` carry-over across cycles must match. This is the source of the documented,
  irreducible DGSCOR ±0.03% cubic-volume drift (the bounded redraw amplifies transcendental-ulp).
- **Validated:** yes (the per-record DG is bit-exact for uncalibrated species; the residual tail is the
  transcendental component, see #9).

## 5. Rounding — Fortran NINT / add-0.5-and-truncate  ·  FAITHFUL ARTIFACT
`src/core/parameters.jl` (`nint`), `src/io/summary.jl` (`di`/`dt`), `src/engine/volume.jl`,
`src/variants/southern/crown_ratio.jl`, `src/engine/establishment.jl`, `src/engine/compress.jl`

- **What:** two distinct Fortran rounding idioms, both reproduced exactly:
  - `nint(x) = round(x, RoundNearestTiesAway)` — Fortran `NINT` (round half AWAY from zero).
  - `trunc(Int, x + 0.5)` — the `.sum`/report integer conversion (add 0.5, truncate toward zero) and the
    crown-ratio/normht `INT(x+0.5)` casts.
- **Canonical:** Julia's `round(x)` (round half to **even**, banker's rounding).
- **Why kept:** Julia's default ties-to-even differs from Fortran's ties-away exactly on the `.5`
  boundary; the `.sum` columns, crown-ratio integers (ICR), NORMHT, and Scribner board-foot all depend on
  the Fortran rule. One wrong tie-break flips a printed integer or a class boundary.
- **Validated:** yes (the `.sum` integer columns + crown-ratio are bit-exact).

## 6. Float32 / REAL*4 arithmetic throughout  ·  FAITHFUL ARTIFACT (pervasive)
the whole `src/` numeric engine

- **What:** every numeric kernel uses `Float32`, including *forcing* Float32 intermediates where Julia
  would otherwise promote to `Float64` (literals written `…f0`, explicit `Float32(...)` casts, e.g.
  `src/io/summary.jl:203` "Computed in Float32 to match FVS REAL*4").
- **Canonical:** `Float64` everywhere (Julia's default; faster on most hardware, more accurate).
- **Why kept:** Fortran FVS is `REAL*4`; the accumulated single-precision rounding is part of the answer.
  Promoting to Float64 changes the last bits of every transcendental and accumulation → different `.sum`
  values and different DGSCOR redraw outcomes.
- **Validated:** yes (foundational — the entire bit-exact result rests on it).
- **Cost:** Float32 transcendentals are not faster on x86, and Float32 is less accurate; we accept both to
  match the oracle.

## 7. Lookup search order — sequential / binary search matching Fortran  ·  FAITHFUL ARTIFACT (minor)
`src/engine/r8clark_vol.jl` (DIBMEN sequential, `_R8CF` binary search), `src/engine/volume_equations.jl`
(`_VOL_SNFIA` binary search)

- **What:** the volume-coefficient tables are searched with the Fortran's exact search (sequential for
  DIBMEN, binary for the keyed Clark/FIA tables).
- **Canonical:** a `Dict` keyed by (geoa, species).
- **Why kept:** when a key appears more than once or a boundary is ambiguous, the *which-row-wins* must
  match the Fortran's search; a `Dict` could resolve duplicates differently. (Low risk, but kept for
  safety since it's volume-critical.)
- **Validated:** yes (volume is bit-exact).

## 8. Fixed-iteration numerical solvers  ·  FAITHFUL ARTIFACT
the Pretzsch self-thinning line solver (`src/variants/southern/mortality.jl` `_pretzsch_tn10`, fixed
100-iteration loop + ±5 convergence band), the height↔DBH inversion (`_htdbh_dbh`), the establishment /
small-tree bounded-height redraw loops (`src/engine/establishment.jl`, `src/variants/southern/small_tree_growth.jl`)

- **What:** Newton/secant/bisection-style loops with FVS's exact iteration cap and exact convergence
  criterion (e.g. "stop when |diff| ≤ 5").
- **Canonical:** `Roots.jl` / a library solver to a tight tolerance.
- **Why kept:** two reasons — (a) the *converged value* is the FVS value only at FVS's tolerance (a tighter
  solve gives a slightly different number); (b) the redraw loops consume LCG uniforms, so their exact trip
  count keeps the RNG stream aligned (#1/#2).
- **Validated:** yes.

## 9. Transcendental functions (exp / log / sqrt / `^`)  ·  IRREDUCIBLE (not "kept")
everywhere

- **What:** Julia's libm and Fortran's may differ in the last ULP of `exp`/`log`/`sqrt`/`pow`; and the
  *evaluation order* of an expression (we mirror the Fortran term order rather than Horner-factoring) is
  preserved to keep Float32 rounding identical.
- **Why it's not an "artifact we keep":** there is no choice that makes these bit-identical across
  toolchains; the residual is the documented "transcendental-ulp" tail (`DIVERGENCES.md`: the DGSCOR
  ±0.03% drift, the occasional ±1 board-foot at a `.sum` boundary). We DO preserve the Fortran expression
  order (don't Horner-factor or rearrange) to minimize it — that part IS a faithful choice.
- **Validated:** the residual is bounded and characterized, not eliminable.

## 10. Normal CDF `_normal_cdf` (FMPOFL torching probability)  ·  BEST-EFFORT STAND-IN ⚠
`src/engine/fire/fmburn.jl` (`_normal_cdf`)

- **What:** a standard normal CDF via **Abramowitz & Stegun 26.2.17**, used by `torching_probability`
  (FVS_PotFire).
- **Honest status:** FVS's `FMPOFL_NPROB` uses a *different* rational approximation, so this is **NOT a
  faithful artifact** — it would not be bit-exact even if it could be checked. It isn't checked: FVS_PotFire
  is **binary-blocked** (the stripped DBS build can't emit the table to diff). So this is a best-effort
  stand-in, validated only by semantics (monotonic in flame, ∈[0,1]).
- **To revisit:** if a fuller FVS binary becomes available, port the exact NPROB rational form and validate.

---

## Summary table

| # | Artifact | File | Canonical alternative | Status |
|---|---|---|---|---|
| 1 | Park–Miller LCG RNG | core/rng.jl | `Random` (Xoshiro/MT) | FAITHFUL (validated) |
| 2 | BACHLO normal rejection | core/rng.jl | `randn` (Ziggurat) | FAITHFUL (validated) |
| 3 | `_rdpsrt!` quicksort + `@goto` | engine/cuts.jl | `sortperm`/`sort!` | FAITHFUL (validated) |
| 4 | AUTCOR/DGSCOR ARMA recurrence | southern/serial_correlation.jl | ARMA lib / `randn` noise | FAITHFUL (validated) |
| 5 | NINT / add-0.5-truncate rounding | parameters.jl, io/summary.jl, … | `round` (ties-to-even) | FAITHFUL (validated) |
| 6 | Float32 / REAL*4 everywhere | all of src/ | Float64 | FAITHFUL (validated, pervasive) |
| 7 | Fortran table search order | r8clark_vol.jl, volume_equations.jl | `Dict` | FAITHFUL (validated, minor) |
| 8 | Fixed-iteration solvers | mortality.jl, establishment.jl, … | `Roots.jl` to tight tol | FAITHFUL (validated) |
| 9 | exp/log/sqrt/`^` + expr order | everywhere | (none) | IRREDUCIBLE (bounded) |
| 10 | `_normal_cdf` (A&S 26.2.17) | fire/fmburn.jl | exact FVS NPROB rational | BEST-EFFORT STAND-IN ⚠ |

**If the bit-exact-vs-Fortran goal is ever dropped**, items 1–8 can be replaced by their canonical
counterparts with no loss of forest-science correctness (only loss of byte-for-byte agreement). Item 9 is
inherent. Item 10 should be replaced by the exact FVS form once a validating binary exists. Until then, the
discipline is: **don't "clean up" any of items 1–8 — the messiness is load-bearing for the oracle match.**

---

## Counter-example — where we DID modernize and accept divergence: COMPRESS
`src/engine/compress.jl`

For contrast: the `COMPRESS` keyword's record-clustering uses a PC-score (eigenvector) projection. The
Fortran uses a 1966 IBM-SSP power-iteration eigensolver; FVSjl substitutes
`LinearAlgebra.eigen`. This is the deliberate INVERSE choice — we picked the modern, correct eigensolver
and **accepted** that the exact class *partition* (and therefore the later-cycle trajectory) is not
bit-identical to FVS, because (a) reproducing a 1966 power-iteration's convergence path is high-cost,
low-value, and (b) the compression-cycle aggregate (`.sum` row) is still conserved and matches. It is
flagged 🟡 in `MANAGEMENT_DISTURBANCE_TODO.md`. This is the one place the project knowingly traded
bit-exactness for a canonical algorithm — every item 1–8 above is where it did NOT, and why.
