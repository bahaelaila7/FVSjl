# ACTIVE GOAL — Fix live-FVS crashes on real FIA inventory (minimal upstream source patches)

## Mission
Live FVS (the Fortran reference, `/workspace/FVSjl/tmp/oracles/FVS{sn,ne,cs,ls}_new`) **crashes** (SIGFPE:
divide-by-zero / invalid arithmetic) on a small set of real FIA stands. These are genuine FVS SOURCE bugs —
unguarded arithmetic on degenerate-but-legal inputs (zero density, off-curve trees, exhausted mortality pools,
etc.). Per the standing directive "on any live-FVS crash, trace + patch + fix live FVS for maintainer submission;
never tolerate as live_crash or replicate in FVSjl" ([[feedback-crash-means-fix-live-fvs]]), root-cause every
crash class and propose the **most minimal, semantically-plausible** Fortran fix — a guard that preserves the
routine's intended behaviour on all non-degenerate inputs and does the mathematically-correct thing on the
degenerate one — document the reasoning, apply the patch to the FVS source (+ a maintainer-submittable
`docs/patches/*.patch`), and verify the crashing stand now runs clean and produces a sensible `.sum`.

## Measured scope (2026-07-18, current tmp/oracles binaries)
12 currently-crashing FIA stands across 6 distinct SIGFPE sites (list: `docs/fvs_livecrash_stands.txt`):
- `cs/varmrt.f:162`  ADJUST=TEMKIL/TEMSUM — div-by-zero when the mortality search has no killable TPA left.
- `cs/grincr.f:437` and `:449` (3+1 stands) — root TBD.
- `cs/fvs.f:197` (1) — root TBD (main driver).
- `ls/dgdriv.f:134` (3) and `:353` (2) — diameter-growth calibration, root TBD.
- `ls/htdbh.f:336` (1) — height↔DBH inverse, root TBD.
(The 60 rows the sweep recorded as `live_crash` are partly STALE — only these 12 still crash on the current
build; the rest ran clean on re-test. FVSjl runs all of them clean — it is the correct side.)

## DOCTRINE
1. **BOTH-SIDES-TRACE the MATH.** Read the crashing routine; identify the exact unguarded operation (division,
   log/sqrt of a non-positive, power of a negative, out-of-domain root). Confirm the degenerate input value by
   reproduction/instrumentation — never guess the trigger.
2. **MINIMAL + SEMANTICALLY PLAUSIBLE.** The fix guards ONLY the degenerate case and returns the mathematically/
   silviculturally correct result there (e.g. no killable TPA ⇒ apply zero additional mortality; off-curve tree
   ⇒ the documented HTMAX/EDH=0.1 fallback). It must be a NO-OP on every input that already worked. Prefer the
   smallest diff (one guard) over refactors. Match the FVS coding idiom + fixed-form columns.
3. **DON'T BREAK A WORKING STAND.** After patching, re-run BOTH the crashing stand (must run clean + sensible
   `.sum`) AND a sample of normal stands (byte-identical `.sum` — the guard must not perturb the non-degenerate
   path). Instrument→build→run→**RESTORE source + rebuild clean .o**→verify oracle pristine (debug-FVS discipline).
4. **DOCUMENT EVERY VERDICT.** One slice per crash class in `docs/FVS_LIVECRASH_AUDIT.md` (site, trigger value,
   root cause, the guard + why it's correct, validation) + a `docs/patches/*.patch`. Add to `FVS_SOURCE_BUGS.md`.
5. **BUILD CAVEAT (honest).** In-container gfortran (12.2.0) ≠ the official oracle build (15.2.1); a relinked
   oracle is a VALIDATION vehicle, not the shipped artifact. The deliverable is the SOURCE PATCH + reasoning +
   in-container validation. If a specific `.o` cannot be recompiled (ABI/`.mod`), say so and validate via the
   build-flag or a source-only path; never claim an unvalidated fix.
6. **VARIANT-SCOPE.** A shared routine (varmrt/grincr/dgdriv/htdbh live in per-variant buildDirs but are often
   textually identical across variants) ⇒ apply the same guard to every variant that carries the routine, and
   note which variants were latent-vulnerable.

## Done-state
Every one of the 6 crash classes: root-caused (trigger value confirmed) + minimally patched + the crashing
stands run clean with a sensible `.sum` + normal stands byte-identical + a documented `docs/patches/*.patch` +
an audit slice. `docs/FVS_SOURCE_BUGS.md` lists each. Then all 12 stands (and a re-sweep of the 60 recorded
live_crash) are crash-free.

## Off-switch
`touch /workspace/FVSjl/docs/FVS_LIVECRASH_COMPLETE` (USER's call). Working log: `docs/FVS_LIVECRASH_AUDIT.md`.
