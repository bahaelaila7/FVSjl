# ACTIVE GOAL — FVSjl structural-backlog campaign (Campaign 2)

## Mission
Make **FVSjl (Southern)** a faithful drop-in for **FVSsn** — *semantically faithful*, not merely bit-exact.
The flag audit→fix campaign is COMPLETE (all 67 flags resolved + 3 bandaid follow-ups closed). This campaign
works the **deferred structural backlog** — the infrastructure/feature items that flag-fixes can't reach.
The only accepted divergences are ULP floating-point and the COMPRESS eigensolver.

## Current task
The ordered backlog in `docs/audit/BACKLOG.md`, **most-upstream / least-dependent / easiest first**:
1. Single-canopy structure-stage (NTREES≤1)            ← START HERE (bounded, independent)
2. Per-point density layer → 2a multi-point pccf, 2b TCONDMLT point weights
3. Snag soft-decay (DECAYX) → 3a SNAGDCAY/SNAGBRK, 3b fire-snag hard→soft
4. density notre! FINT/FINTM dead inflation
5. COMPRESS #29 post-compression DGSCOR residual (greens the suite)
6. NOHTDREG/LHTDRG HT-DBH calibration subsystem
7. log-graded HRVRVN revenue (ecvol log-bucking)
8. FFE phasing #28 co-refactor (+ fire-carbon released-from-fire value)

Method (principle #4): port the semantics FAITHFULLY + double-check vs LIVE FVS (binary at
`/workspace/FVSjl/tmp/FVSsn_full`; `LiveFVS.run_key`), THEN write a scenario that exercises it. Document each
item's verdict + both-sides logic in `docs/audit/` (INDEX.md / BACKLOG.md). Keep the suite green (or move it in
a verified-faithful way). Within an infra item, build the layer first, then its small consumers.

## The 5 principles (IMPRINTED — apply to every verification and fix)
1. **TRACE LOGIC, NOT RUNTIME.** Trace the logic in BOTH FVSjl and the FVS Fortran; match SEMANTICS.
2. **UPSTREAM-FIRST.** Work the most-upstream, least-dependent issue first.
3. **REGRESSION = MASKED-BUG SIGNAL.** A semantically-certain fix that regresses other tests likely unmasked a
   hidden bug — NOTE it, keep going, circle back. Do NOT revert a faithful fix to keep a stale test green.
4. **PORT THE SEMANTICS FAITHFULLY FIRST + double-check vs live, THEN write the test.** Never test-first.
5. **DOCUMENT EVERYTHING** in `docs/audit/` so we never run in circles.
6. **VARIANT-AWARE — DO NOT HARDEN TO SOUTHERN.** In shared/base modules (`src/engine/*`, `src/io/*`,
   `src/core/*`) gate variant-specific behavior on the variant/coefficients; SN-only code lives under
   `src/variants/southern/`. Match the FVS **base** routine's general semantics.

## Re-trace discipline (this campaign's hard-won lesson)
A "deferred" or "cleared" label can be a misread — the flag campaign caught real bugs (event-monitor `**`,
redcedar TFALL, compress decay-averaging) AND stale already-done flags (FMDYN) by re-tracing against SOURCE.
Re-verify each backlog item against the FVS Fortran before assuming its scope.

## Bar for "done" with any item
Traced in both implementations · semantics matched to FVS source (file:line cited) · live-validated where an
oracle exists · suite stays bit-exact OR moves verified-faithful · verdict + fix logged in `docs/audit/`.

## Off-switch
This reminder fires from a Stop hook until the campaign is complete. To silence it:
`touch /workspace/FVSjl/docs/audit/CAMPAIGN_DONE`
