# ACTIVE GOAL — FVSjl Lake States (LS) variant port

## Mission
Port FVS **Lake States (LS)** into FVSjl as a *semantically faithful* drop-in for the live
Fortran **FVSls** — the same standard that made Southern, Northeast, and Central States validated
drop-ins (tags `FVSsn-complete`, `FVSne-complete`, `FVSsn+ne`, `FVScs` / `FVSsn+ne+cs-done`).
Bit-exact end-to-end barring only ULP single-precision and documented eigensolver-class divergences.
Full scope + per-routine NE-reuse table + chunk order: **`docs/LS_VARIANT_PORT_SCOPE.md`**.

## Oracle (LS-specific — read carefully)
- Sole oracle = **live FVSls**: relink from `bin/FVSls_buildDir/*.o` + the glibc shim, exactly like
  NE/CS. Harness = **`test/harness/ls_oracle.sh`** (BUILDDIR=`FVSls_buildDir`, BIN=`/tmp/FVSls_new`) —
  already created and producing `lst01.sum`. There is **NO Oracle-A / FVSjulia for LS** — no 1:1
  transliteration exists; validate against the live binary only.
- Canonical test stand: **`/workspace/ForestVegetationSimulator/tests/FVSls/lst01.key`** (+ `lst01.tre`).
  Treat any committed `.sum.save` as potentially STALE — validate cyc0+ against the freshly-relinked
  live binary. Live cycle-0 (1990) ground truth for lst01 (forest 401, 10-yr cycles):
  `TPA 536 · BA 77 · SDI 160 · CCF 171 · TopHt 63 · QMD 5.1 · TCuFt 1551 · MCuFt 1338 · SCuFt 480 · BdFt 1887`.
- Canonical runner: `julia --project=. test/runtests.jl` (NOT `Pkg.test`). Record the current suite
  baseline in the scope doc before starting; adding LS must NOT regress SN/NE/CS.

## The leverage (why LS is another cheap variant)
- **~90%+ of the FVSls source is the shared base** already ported variant-agnostically (engine, keyword
  dispatch, density, `grow_cycle!`/`mortality_and_fire!`, FFE base + DBS tables, ECON, IO).
- The variant surface is **~30 routines**; LS is an **eastern** variant so most are near-NE/CS. The
  scope doc holds the measured reuse table (a background scan is producing it). Known so far:
  **`htcalc.f` is IDENTICAL to NE**; the genuinely-LS-specific models are **`dgf.f`** (an SN-family
  ln(DDS) per-species regression, ~508 lines — extend the Southern `dgf!` framework, NOT NE's
  BAL-potential iteration), **`crown.f`**, and the **volume/merch** path (confirm R9 Clark vs LS-specific).
- **68 species** (LS `MAXSP = 68`; SN 90, NE 108, CS 96) — the LS coefficient CSVs are the main data lift.

## Infra constants (confirmed from FVSls_buildDir)
- `MAXSP = 68` (PRGPRM.F77) · `YR = 10.0` (blkdat.f:71 — 10-yr cycle, like NE/CS) ·
  RNG seed `S0/SS = 55329` (blkdat.f:302 — same as SN/NE/CS) · `LZEIDE = .TRUE.` (grinit.f:126 —
  Zeide SDI, like NE/CS). Variant designator `"LS"`; lst01 forest-code default 401.

## Ordered work (most-upstream / least-dependent first)
1. **Variant infra:** `LakeStates <: AbstractVariant`, `variant_code(::LakeStates)="LS"`,
   `nspecies→68`, `htg_period→10`, Zeide SDI, LS FORKOD/site defaults, LS species CSVs, registry entry.
   ⇒ drive `lst01` **cycle-0** stand columns bit-exact (TPA/BA/SDI/CCF/QMD/TopHt) + `test/integration/test_lst01.jl`.
2. **Volume:** wire LS volume equation ids into the existing R9 Clark + R9LOGS path (or LS-specific if it
   differs — confirm from vollib). ⇒ cycle-0 `.sum` volume columns (1551/1338/480/1887).
3. **Diameter growth:** `ls/dgf.f` (the SN-family LS DDS model) + LS competition/`balmod`; `dgdriv`
   calibration shared/near-NE. ⇒ cycle-1 DG.
4. **Height growth** (reuse NE `htgf`/`htcalc` + LS coefs) · **crown** (`ls/crown.f`) · **mortality**
   (`varmrt`) · `htdbh`/`sitset`/`forkod`. ⇒ drive `lst01` to cycle-1+ vs live.
5. LS-active shared branches (un-gate, don't harden) → FFE (`fire/ls/` if present), establishment,
   sprouting, thinning — validate each at its key cycle vs live.

## THE DOCTRINE (imprinted — apply to every verification and fix)
1. **TRACE LOGIC, NOT RUNTIME.** Trace the logic in BOTH FVSjl and the FVS LS Fortran; match SEMANTICS.
2. **UPSTREAM-FIRST.** Work the most-upstream, least-dependent issue first.
3. **REGRESSION = MASKED-BUG SIGNAL.** A semantically-certain fix that regresses other tests likely
   unmasked a hidden bug — NOTE it, keep going, circle back. Do NOT revert a faithful fix to keep a
   stale test green (a `.sum.save` can be STALE — the oracle, not your fix, may be wrong).
4. **PORT SEMANTICS FAITHFULLY FIRST + double-check vs LIVE FVSls, THEN write the test.** Never test-first.
   A green test on a degenerate scenario (calibration never fired, keyword path not triggered) is vacuous.
5. **ALL TEST TOLERANCES = 0 EXCEPT PROVEN ULP.** Do not use a global cutoff to declare a residual "ULP".
   Corner the exact low-level Float32 operation that produces the differing ULP (e.g. non-associative
   sum order) before accepting it. Document each accepted ULP with its traced root.
6. **DOCUMENT EVERY verdict + both-sides logic** (docs/LS_*), so we never run in circles.
7. **VARIANT-AWARE — DO NOT HARDEN.** In shared `src/engine|io|core/*` gate variant behavior on the
   variant/coefficients; LS-specific code under `src/variants/lakestates/`. Adding LS must keep **SN,
   NE, AND CS bit-exact**. Reuse the NE/CS routine where LS ≈ it; write LS-specific only where it differs.

## Re-trace discipline
A "done" / "ported" / "reuses NE" label can be a misread — re-verify against the LS Fortran source
before assuming scope (% similarity is a hint, not a guarantee — a 96%-similar routine can hide an
LS-specific coefficient or branch). Prior campaigns repeatedly caught real bugs behind stale
already-done flags and behind "irreducible-RNG" floors (trace every floor to ground).

## Off-switch
This reminder fires from a Stop hook until the LS port is complete. To silence it:
`touch /workspace/FVSjl/docs/LS_COMPLETE`
