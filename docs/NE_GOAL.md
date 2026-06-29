# ACTIVE GOAL — FVSjl Northeast (NE) variant port

## Mission
Port FVS **Northeast (NE)** into FVSjl as a *semantically faithful* drop-in for the live
Fortran **FVSne** — the same standard that made Southern a validated drop-in (tag
`FVSsn-complete`). Bit-exact end-to-end barring only ULP single-precision and documented
eigensolver-class divergences.

## Oracle (NE-specific — read carefully)
- Sole oracle = **live FVSne**: `bash test/harness/ne_oracle.sh <key> <outdir>` relinks
  `/tmp/FVSne_new` from `bin/FVSne_buildDir/*.o` + the glibc shim (shipped `bin/FVSne` is
  GLIBC-broken). There is **NO Oracle-A / FVSjulia for NE** — no 1:1 transliteration exists.
- `tests/FVSne/net01.sum.save` is **STALE** (source 20250930 vs current 20260401): cyc-0
  matches, cyc1+ DIVERGE. Validate cyc1+ against the LIVE `FVSne_new`; regenerate the `.save`
  before promoting cyc1+ asserts.
- Canonical runner: `julia --project=. test/runtests.jl` (NOT `Pkg.test` — the NTuple{9}
  BoundsError is a stale-sandbox-manifest artifact). Suite baseline 5172 pass / 2 broken
  (the 2 broken = accepted SN COMPRESS eigensolver + NOHTDREG ULP — not NE).

## Current state (continue from here)
- Variant infra DONE: `Northeast<:AbstractVariant`, MAXSP 90→108 + per-variant `nspecies`,
  tolerant CSV loader, variant-dispatched `site_setup!`, NE FORKOD defaults (IFOR=2).
- net01 **cycle-0** bit-exact on 6 stand columns (TPA 536 / BA 77 / SDI 160 / CCF 176 /
  QMD 5.1 / TopHt 63) — `test/integration/test_net01.jl`.
- R9 Clark CUBIC volume PORTED as WIP (`src/engine/r9clark_vol.jl`, NVEL R9 — NOT R8 Clark)
  and validated PER-TREE <1% vs the live `.trl` — but NOT yet wired into `compute_volumes!`.

## Ordered work (most-upstream / least-dependent first)
1. Wire R9 volume into `compute_volumes!` via variant dispatch (+ the d<1→0 guard); close the
   <1% residual; add board feet (R9LOGS / r9bdft Scribner, vol[2]). ⇒ cycle-0 .sum volume cols.
2. Diameter growth (`ne/dgf.f`) — incl. the structurally-NEW **BAL competition** (badist/balmod).
3. Height growth (`htgf`) · crown ratio (cwcalc species selection) · mortality (`morts`).
4. ⇒ drive net01 to cycle-1+, validated vs live `FVSne_new`; then NE-active shared branches + FFE.

## THE DOCTRINE (imprinted — apply to every verification and fix)
1. **TRACE LOGIC, NOT RUNTIME.** Trace the logic in BOTH FVSjl and the FVS NE Fortran; match SEMANTICS.
2. **UPSTREAM-FIRST.** Work the most-upstream, least-dependent issue first.
3. **REGRESSION = MASKED-BUG SIGNAL.** A semantically-certain fix that regresses other tests likely
   unmasked a hidden bug — NOTE it, keep going, circle back. Do NOT revert a faithful fix to keep a
   stale test green (and remember the NE `.save` can be STALE — the oracle, not your fix, may be wrong).
4. **PORT SEMANTICS FAITHFULLY FIRST + double-check vs LIVE FVSne, THEN write the test.** Never test-first.
5. **DOCUMENT EVERY verdict + both-sides logic** (docs/NE_*), so we never run in circles.
6. **VARIANT-AWARE — DO NOT HARDEN.** In shared `src/engine|io|core/*` gate variant behavior on the
   variant/coefficients; SN code lives under `src/variants/southern/`, NE under `src/variants/northeast/`.
   Match the FVS **base** routine's general semantics — adding NE must keep SN bit-exact (don't harden
   base code to either variant).

## Re-trace discipline
A "done" / "ported" / "validated" label can be a misread — re-verify against the NE Fortran source before
assuming scope. The SN campaign repeatedly caught real bugs hiding behind stale already-done flags.

## Off-switch
This reminder fires from a Stop hook until the NE port is complete. To silence it:
`touch /workspace/FVSjl/docs/NE_COMPLETE`
