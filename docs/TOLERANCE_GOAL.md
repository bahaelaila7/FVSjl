# ACTIVE GOAL — FVSjl test-tolerance closure campaign

## Mission
Drive EVERY numerical tolerance in the FVSjl test suite to one of exactly two states:
1. **BIT-EXACT** — the assertion is `==` (or a comparison of FVSjl's *formatted* `.sum`/report
   output against the live-Fortran golden, which is by construction exact), tolerance 0; OR
2. **PROVEN-IRREDUCIBLE ULP** — the residual is CORNERED to the exact low-level Float32
   operation that produces it (e.g. a named non-associative sum order, a `^(7/6)` transcendental,
   a print half-width against a rounded oracle field), with the traced root documented IN the test,
   and the bound set to exactly that irreducible width — not a loosened multiple of it.

No empirical bounds survive: no "measured floor × 1.5", no percentages (2.5%, 3%, `1+0.002·x`),
no multi-unit integer slack (`≤ 4`, `CCF ≤ 10`, `atol=8`) that merely *covers a residual class*
without cornering it. "Documented accepted class" is NOT "proven ULP for this exact comparison."

The **only** permitted non-passing assertions are `@test_broken` whose root is a documented
eigensolver-class / genuinely-irreducible divergence (COMPRESS s22, the R8-VOLUME s32 leak, the
WK3/DGSCOR sp33/65 tail) — and even those must carry a precise both-sides traced verdict.

## THE DOCTRINE (imprinted — apply to every fix)
1. **TRACE LOGIC, NOT RUNTIME.** Trace the residual in BOTH FVSjl and the live FVS Fortran; match
   SEMANTICS. Find WHY the bits differ before touching a bound.
2. **UPSTREAM-FIRST.** Fix the most-upstream cause. A grown-cycle BA tolerance usually traces to a
   DGF / mortality / volume op — fix THAT (bit-exact), which collapses many downstream bounds at once.
3. **REGRESSION = MASKED-BUG SIGNAL.** A semantically-certain fix that regresses another test likely
   unmasked a hidden bug — NOTE it, keep going, don't revert a faithful fix to keep a stale bound green.
4. **FIX/CORNER FAITHFULLY FIRST + verify vs LIVE FVS, THEN tighten the test.** Never loosen a bound
   to pass; never declare ULP without cornering the op. A green test on a loosened bound is a lie.
5. **TOLERANCES = 0 EXCEPT PROVEN ULP.** This IS the campaign. For each bound: either drive it to `==`
   by fixing the FVSjl op to bit-exactness, OR corner the exact Float32 op and prove the bound is the
   irreducible width (prefer comparing FORMATTED output with `==` — the `.sum` prints fixed decimals,
   so a print-rounded field is exactly matched by comparing the rendered string / rounded value).
6. **DOCUMENT EVERY verdict + both-sides logic** in the test comment and this file — so we never
   re-litigate. Each closed item: BIT-EXACT (op fixed) or ULP (op named + width proven).
7. **VARIANT-AWARE — DO NOT HARDEN.** Shared `src/engine|io|core` gates on variant/coefficients;
   variant-specific code under `src/variants/<v>/`. Every fix must keep **SN, NE, CS, AND LS bit-exact**
   (all four are validated drop-ins: tags FVSsn-complete … FVSsn+ne+cs+ls-done).
8. **DECONFOUND THE ELEMENTARY OPS VIA A FORTRAN COMPANION (FFI).** Isolate the SMALLEST set of
   known compounding low-level operations — the not-correctly-rounded transcendentals `exp`/`log`/
   `pow(**)` (NOT `sqrt`, which is IEEE-correctly-rounded and already bit-identical) and the
   eigensolver (COMPRESS EIGEN/Jacobi) — into a gfortran companion library (compiled from / linking
   the SAME libm as `bin/FVS*_buildDir`) and call the REAL ones from Julia through `ccall`/FFI. Keep
   the pure-Julia implementations in place but RENAMED `xxx_julia`; the default `xxx` dispatches to the
   Fortran op. Purpose: REMOVE op-level ULP as a confound so any residual that survives is a genuine
   SEMANTIC (logic) mismatch, not a libm rounding artifact — this is what lets doctrine #1–#5 actually
   bite (no more "irreducible-ULP" fog masking real bugs like the QMDGE5 cap). Prove the premise first
   (bit-compare Julia `exp/log/pow` vs gfortran `expf/logf/powf` over the real input ranges); only wire
   the FFI for the ops that actually differ. Keep the companion minimal, documented, and variant-safe.
9. **EXPOSE, DON'T HIDE — every non-bit-exact residual is `@test_broken`, never a green `tol>0`.** The
   END STATE: the ONLY assertions that may remain `@test_broken` (or carry a residual) are ULPs cornered to
   ONE fundamental, ISOLATED, PORTABLE, FVS-SEMANTICS-FREE numeric primitive — the eigensolver (EIGEN/Jacobi),
   the COR / ARMA serial-correlation recurrence, or a transcendental (`exp`/`log`/`pow`) — i.e. irreducible
   only because a low-level primitive rounds differently and we have not or cannot FFI it. EVERYTHING ELSE —
   every semantic/logic/ordering difference — must be driven to BIT-EXACT (`==` or rendered-`==`). Crucially,
   a PASSING `tol>0` assertion HIDES a non-bit-exact residual inside the green suite (a lie by omission). So
   convert every surviving `tol>0` into a `@test_broken` with a bound TIGHTER than the current one — ideally
   `@test_broken ==` / rendered-`==`, or `@test_broken all(rows bit-exact)` — so the residual is VISIBLE as
   broken until it is actually closed. GREEN ⇔ bit-exact; BROKEN ⇔ a documented, cornered, still-open residual
   (with its primitive named). No residual passes silently. Sum-order accumulations that are NOT one portable
   primitive do NOT get a free pass — they become `@test_broken` too, until matched to FVS's loop order or
   proven to reduce to a named primitive.

## Oracle & runner
- Oracle = live Fortran per variant: `/tmp/FVS{sn,ne,cs,ls}_new` (relink from `bin/FVS*_buildDir/*.o`
  + the glibc shim, via `test/harness/*_oracle.sh`). Validate every fix against the freshly-relinked
  live binary, not a possibly-stale golden.
- Runner: `julia --project=. test/runtests.jl` (NOT Pkg.test). Baseline before starting: **7658 pass /
  2 broken**. No fix may regress SN/NE/CS/LS or drop a pass.

## The work list
The full per-line inventory is **docs/TOLERANCE_AUDIT.md** (the checklist — tick each as BIT-EXACT or
ULP-PROVEN with its traced root). Ordered by leverage (upstream-first):

1. **Re-measure post-QMDGE5-fix** — the CS-family DGF fix this session may already have made the CS
   all-species (`test_allspecies.jl` CS rows) and `test_cst01.jl` grown/later bounds bit-exact or much
   tighter. Re-run, tighten to the true floor, corner the remainder. (Cheapest wins first.)
2. **Volume / board-foot "Scribner Float32 noise"** (`1+0.002·x` / `1+0.005·x`, `BdFt atol`) across
   ~12 keyword tests + net01/cst01/dbs — corner to the exact Scribner/Clark round or fix the op; ideally
   compare the `.sum`-rendered integer with `==`.
3. **The ~69 `±1/±2` per-column bounds** — for each, either prove it's a specific print/sum-order ULP
   (then compare rendered output `==`) or fix the op. This is the biggest population.
4. **Grown-cycle percentage bounds** — `test_allspecies` `_ALLSP_TOL`, `test_timeint` 3%,
   `test_multicycle` rtol, `test_carbon` `0.005·v+0.1` — trace each residual to its op.
5. **FFE fire/carbon loosest atols** — `test_carbon` FFE rows, `test_lst01_ffe` scorch/flame,
   `test_fire` BA/TPA — the fire kill-distribution + transcendental scorch: corner or fix.
6. **@test_broken** — re-verify each root is genuinely irreducible; document the both-sides verdict.

## Re-trace discipline
An "accepted / ULP / documented" label on a bound may be a misread — re-verify against the live source
before trusting it. Prior campaigns repeatedly found REAL bugs behind "irreducible-RNG"/"ULP" floors
(the LS QMDGE5 cap this very session was hiding behind a "terminal tripling-spread" label). Trace every
floor to ground.

## Off-switch
This reminder fires from the Stop hook until every tolerance is closed. To silence it (only when the
audit checklist is 100% BIT-EXACT or PROVEN-ULP-with-traced-root):
`touch /workspace/FVSjl/docs/TOLERANCE_COMPLETE`
