# ACTIVE GOAL — FVSjl FIA/FVS behaviour-compatibility campaign (4 variants)

## Mission
FVSjl is a validated bit-exact-or-cornered drop-in for FVS SN+NE+CS+LS on the *curated* test scenarios
(the modernization + tolerance campaigns are CLOSED — `docs/MODERNIZATION_COMPLETE`,
`docs/TOLERANCE_COMPLETE`; floor 38527/143/0). This campaign proves that fidelity holds on **real FIA
inventory at scale, over full multi-cycle projections and under management**, for every variant — i.e.
FVSjl reproduces live FVS's *behaviour*, not just cycle-0 inventory, on the plots foresters actually run.

Oracle = **freshly-relinked live FVS** per variant (`/tmp/FVS{sn,ne,cs,ls}_new`, via
`test/harness/*_oracle.sh`). Inputs = real plots from `SQLITE_FIADB_ENTIRE.db` (read-only, NEVER modified)
through the native `DATABASE`/`DSNIN` reader. Every claim is measured vs the live binary, not a stale golden.

Baseline already achieved (modernization #85/#94): 162 real cross-variant stands, **cycle-0** inventory
bit-exact on 10 `.sum` cols (159/162; 7/10 cols perfect on all 162). This campaign extends that from
*one inventory row* to *the whole projected trajectory + management*, and from 162 to a stratified sample.

## The FOUR PILLARS (each with a measurable done-state)

### 1. SCALE & STRATIFICATION
- A documented, reproducible **stratified sample of real FIA plots per variant** — spanning forest types,
  stand structures (even/uneven-aged, seedling→sawtimber), site classes, and the geographies each variant
  covers — extracted deterministically from `SQLITE_FIADB_ENTIRE.db`.
- Done-state: a per-variant plot manifest (plot IDs + strata) + an extraction script that regenerates it;
  materially larger than the 162-stand baseline.

### 2. MULTI-CYCLE PROJECTION COMPATIBILITY
- Each sampled plot **projected the full default horizon** (not cycle-0 only) through both live FVS and
  FVSjl; **every cycle's `.sum` row bit-exact-or-cornered** vs live on all 10 columns (TPA/BA/SDI/CCF/
  TopHt/QMD/TCuFt/MCuFt/SCuFt/BdFt), for all 4 variants.
- Done-state: a projection differential over the sample; per-variant pass rate documented; every residual
  bit-exact or cornered to a named primitive.

### 3. MANAGEMENT-SCENARIO COMPATIBILITY
- Real plots run under **standard silvicultural regimes** (thinning by BA/TPA/DBH, salvage, planting/
  regeneration, prescribed fire/SIMFIRE) match live FVS across the projection — the keyword behaviour,
  validated on real inventory rather than synthetic stands.
- Done-state: a management-scenario differential over a plot subset per variant; bit-exact-or-cornered.

### 4. DIVERGENCE TAXONOMY & CORNERING
- Every non-bit-exact plot/cycle **root-caused BOTH-SIDES-TRACED** (FVS source + FVSjl) and either FIXED
  (without regressing the 38527/143/0 floor) or cornered to a named primitive (ULP-class direct/compounded,
  FVS bug, FVS-UB) in `docs/FIA_FVS_COMPAT_AUDIT.md` + `docs/FVS_SOURCE_BUGS.md`.
- Done-state: a documented divergence taxonomy; no unexplained divergence remains.

## DOCTRINE (carried from the closed campaigns — apply to EVERY change)
1. **NEVER REGRESS THE FLOOR.** `julia --project=. test/runtests.jl` = 38527 pass / 143 broken / 0 fail.
   Every fix keeps it. The 143 broken stay the exact documented cornered set (see MODERNIZATION_COMPLETE).
2. **VALIDATE VS FRESHLY-RELINKED LIVE FVS**, never a stale golden. Regenerate the oracle per run.
3. **BOTH-SIDES-TRACE.** A divergence verdict requires reading the FVS source AND the FVSjl path — never
   infer a fix from test pass/fail (the s32 prod=="01" lesson). Bit-exact-on-a-test ≠ faithful.
4. **BIT-EXACT OR CORNERED-TO-A-NAMED-PRIMITIVE.** Float32 op-order is semantics; a residual is acceptable
   ONLY when named (ULP direct/compounded, FVS bug, FVS-UB) — not a padded tolerance that hides logic gaps.
5. **VARIANT-SAFE.** Shared `src/engine|io|core` gates on variant/coefficients; specifics under
   `src/variants/<v>/`. Every change keeps SN+NE+CS+LS correct.
6. **MEASURE, DON'T GUESS.** Per-tree/per-cycle differentials, debug-FVS stamps (instrument .f → build →
   run → RESTORE source + rebuild clean .o → verify oracle pristine). SQLITE_FIADB_ENTIRE.db is read-only.
7. **INCREMENTAL + REVERSIBLE + DOCUMENTED.** One slice at a time; validate before the next; record each
   slice in `docs/FIA_FVS_COMPAT_AUDIT.md` (plots covered, pass rate, divergences found/cornered/fixed).

## Infrastructure
- FIA reader + translator: `read_fia_database` / raw-FIA translator (modernization #86/#87), validated.
- Harness: `test/harness/fia/` (validate_fia.jl, validate_fia_cols.jl, sweep/). Extend for multi-cycle +
  management differentials.
- Oracle relink: `test/harness/{sn,ne,cs,ls}_oracle.sh`. Suite runner: `julia --project=. test/runtests.jl`.
- Working checklist: `docs/FIA_FVS_COMPAT_AUDIT.md` (create as the first slice lands).

## Off-switch
Fires from the Stop hook until all four pillars are met AND the suite is at the floor. To silence:
`touch /workspace/FVSjl/docs/FIA_FVS_COMPAT_COMPLETE`
