# FVSjl porting & validation doctrine

The non-negotiable rules for porting FVS (Fortran) to FVSjl (Julia) and proving it correct.
Consolidated from the `feedback-*.md` memory notes and the practices used throughout the port.
When any rule below conflicts with convenience, the rule wins.

## 1. Validate BIT-EXACT-OR-CORNERED vs the LIVE relinked Fortran oracle
Every chunk is validated against a freshly-relinked live Fortran build — `FVS{v}_g16`
(gfortran-16 `-std=legacy -w -fno-automatic` + isoc23 shim) or `FVS{v}_clean`, under
`/workspace/.{v}work`. **Never** validate against a prior Julia reimplementation
(FVSjulia), against a port comment, or against a test that doesn't exercise the real
semantic. After a container restart, RELINK the oracles (do not trust a stale binary).

## 2. MEASURE — never infer
- Instrument the oracle (single-`.o` dump swap on the g16 build) and A/B the actual
  numbers at the divergence point, byte-first.
- Read the **FVS Fortran as ground truth** — not the Julia code's own comments (which
  have been wrong: the NC/UT/OC episodes), not a test's expectations.
- "measure-don't-infer" applies to your OWN fix too: implement it, then A/B it.
- The full end-to-end `run_keyfile` `.sum` diff is REQUIRED. Controlled-input
  per-subsystem tests pass while the full engine hides crash / latitude / straddle bugs.
- **Per-record / per-tree treelist diffs are INVALID after record TRIPLING** — the
  tripled records don't correspond 1:1 across the two engines. Validate on the `.sum`
  aggregates, or on the pre-split (pre-tripling) window.

## 3. "Cornered" is a MEASURED verdict, not a label of convenience
A divergence may be called *cornered* only when it has been measured and reduced to a
**named floating-point primitive**: the #206 OLDRN serial-correlation growth straddle
(active when DGSD≥1), the RDPSRT unstable-quicksort self-thin tie-break, a DGSCOR /
volume ULP, or the FVS-flagged >1000-TPA-per-record mega-record instability. An
unexplained difference is **not** cornered — it is a needs-dig.
⚠ A "cornered" label can MASK a real bug (proven repeatedly: NSTORE = 3 real bugs, EC
TCuFt, the NC "extreme-dense" residual = a small-tree-growth bug). Re-audit corners
against the live oracle; do not inherit a stale cornered verdict.

## 4. Float & RNG rules
- Use a glibc libm `ccall` for Float32 transcendentals (exp/log/pow) so they match
  gfortran bit-for-bit. `src/core/fmath.jl` provides `fexp/flog/fpow`.
- **NEVER FFI an RNG.** Port the RNG faithfully (BMRANN/ESRANN/etc.) so the draw
  sequence matches; an FFI'd RNG is not reproducible or portable.

## 5. The hard gate
`test/integration/test_multicycle.jl` MUST read **339 pass / 11 broken, BYTE-IDENTICAL**
after every commit. It is the guard that a change didn't perturb the validated engine.
Run it under DEFAULT bounds. `--check-bounds=yes` is the crash-repro tool, not the gate.

A new variant / extension / subsystem is landed **ADDITIVE + INERT**: it adds code
but no `simulate.jl` seam (or a seam gated off) until its effect is validated, so the
gate stays byte-identical the whole time. Wire the live seam only once the effect is
proven bit-exact-or-cornered.

## 6. Commit, reuse & merge discipline
- REUSE the shared engine — dispatch a variant into the shared kernels/drivers rather
  than duplicating them (a bug fixed once in a shared driver fixes every variant; a
  divergence found in one often lurks in its siblings — audit the family).
- Commit ONLY validated chunks. A branch with an unverified fix is marked WIP and NOT merged.
- Each fix must be variant-DISPATCH-GUARDED so it cannot regress other variants.
- Verify "merged" by RUNNING the ref `.sum`, not by `git cherry` (patch-id gives
  false negatives under different SHAs).

## 7. Trust nothing on an agent's say-so — re-verify through the committed path
An agent's "bit-exact COMPLETE" self-report is NOT verification. Before merging any
agent's work: run its OWN test through the COMMITTED engine path (not a side driver),
and run the gate, in the main session. (The OC crown episode: the agent validated a
side per-tree driver, its committed test failed 3/6.)

## 8. Bias vs straddle needs a SIGNED ≥10-stand sign-tally
A single stand cannot distinguish a systematic bias from a floating-point straddle.
Run ≥10 stands, tally the SIGN of (jl − live) at the divergence cycle: one-directional
⇒ real bug; mixed ⇒ cornered. Magnitude-only ledgers cannot make this call.

## 9. Sweep bugs are REGIME-specific, not variant-specific
Mature reference stands are bit-exact-or-cornered, but LARGE systematic bugs live in
the SEEDLING / WOODLAND / BARE-ESTAB regimes that ref stands don't represent. A "full
sweep" must STRATIFY by regime. This is why full-population FIA sweeping finds bugs
that sampled/ref-stand validation structurally cannot.

## 10. Source archaeology before an "un-portable" verdict
"deleted ≠ absent" — check `git log --all --diff-filter=D` before declaring source
absent. Build the historical **self-consistent** tree (all at one revision) before an
"un-oracle-able" verdict — don't link old code against a modern base. (The PPE/FVSppe
recovery.)

## 11. Persistent infrastructure
Durable state lives on `/workspace` (btrfs, survives container restart), NEVER `/tmp`.
`JULIA_DEPOT_PATH=/workspace/.julia_depot`. Oracles under `/workspace/.{v}work`.
Long-running jobs must be cursor-checkpointed + resumable (see the FIA sweep harness).
Do NOT read the live sweep DB while a sweep writes it.

## 12. Scope
Port EVERY column / subsystem bit-exact, including report-only and cosmetic output.
When the roadmap order is clear, don't stop-and-ask — proceed and validate. The
off-switch is the user's call.

## Executing a fix — the method (how, not just what)
Source: `feedback-audit-fix-doctrine`, `feedback-verify-semantics-from-fvs-code`,
`feedback-decision-flow-is-completeness-oracle`, `feedback-test-must-exercise-the-semantic`.

1. **UPSTREAM-FIRST / least-dependent first.** Work the most-upstream, least-dependent
   divergence first — dependents inherit the corrected upstream behavior, so a downstream
   symptom often vanishes once its root is fixed. (The NC "under-mortality" was really an
   upstream small-tree-GROWTH starvation; don't patch the symptom.)
2. **MAP THE FVS SEMANTICS HONESTLY FIRST, before attempting any fix.** Read the exact FVS
   subroutine, trace the data flow (which call sets which value, what gates each block, the
   restore/override order), and implement ONLY what the FVS code does.
3. **TRACE LOGIC, NOT RUNTIME — read the code, don't fit the test.** A fix must match FVS
   SEMANTICS in both the Fortran and the jl path. NEVER infer a fix from test pass/fail.
   "It made the suite pass" is how bandaids + masked bugs are born (the s32 `prod=="01"`
   proxy passed the tested cases but diverged on untested ones). **Bit-exact-on-a-test ≠
   faithful** — bit-exactness on tested scenarios is necessary but NOT sufficient.
4. **PORT THE SEMANTICS FAITHFULLY FIRST, THEN write the test.** Never test-first on a gap
   (it invites fitting the code to the test). And the test scenario MUST actually fire the
   ported branch — if the ported code can be deleted and the test still passes, the test is
   vacuous. Prefer an empirical differential against the live oracle to settle "does X matter".
5. **REGRESSION FROM A FAITHFUL FIX = MASKED-BUG SIGNAL, not a reason to revert.** If a
   semantically-certain fix regresses other tests/stands, that is a GOOD sign it unmasked a
   hidden bug. NOTE it, keep going upstream-first, circle back. Do NOT revert a faithful fix
   to keep a stale test green. (This is why a fix invalidates the whole affected population,
   not just the flagged stands — re-sweep the passers too; see the sweep discipline below.)
6. **SERIALIZE THE GATE.** One fix → run the suite/gate → next fix, so any regression stays
   attributable to the change that caused it. Each fix is variant-dispatch-guarded.
7. **CHUNK large Fortran files.** Skeleton + a saved chunk plan; execute one chunk per
   session (read only that chunk's source), mark verified; don't re-read covered sections.
8. **DOCUMENT every verdict** (the both-sides logic, the faithful fix, any noted masked-bug
   regression) so the work never runs in circles.
9. **Drive completeness from the FVS decision flow, NOT from tests.** Destructure any coarse
   node into its atomic branches, each with a port status, so untested gaps can't hide.

## FIA full-population sweep discipline (recovered from the eastern 1.47M sweep)
The eastern SN/NE/CS/LS sweep reached a *trustworthy* 99.99% because it ran this way.
Deviating (as happened on the western sweep — DIGCAP raised, backlog ballooned to ~110k)
produces provisional, untrustworthy numbers. Source: `docs/FIA_FVS_COMPAT_GOAL.md`,
`docs/fia_ledger.README.md`, `docs/fia_divergence_taxonomy.md`, `run_expand_loop.sh`.

1. **CAP-AND-FIX: DIGCAP ≈ 100 (never raise it for fresh clusters).** The loop PAUSES when
   the dig-queue reaches ~100–200 material discrepancies; you then root-cause and fix them
   BEFORE the sweep continues. One slice at a time; validate before the next. (The one time
   the east raised the cap 100→500, it was for a batch that was a KNOWN already-classified
   primitive — never for fresh unexplained divergence.)
2. **STATUS-FLIP LEDGER: re-run a FIXED stand set after EVERY fix, diff BOTH directions.**
   Keep a committed per-stand ledger (e.g. 1000 stratified stands/variant). After any fix,
   re-run it and diff: a `needs_dig` that should now clear is only half — a previously
   `bit_exact` stand that REGRESSED is the dangerous half. The 339/11 gate guards only the
   ref stands, NOT the swept population, so a fix is not proven safe until the passers are
   re-swept and shown unchanged.
3. **A fix invalidates the population, not just the flagged stands.** A DB that mixes stands
   swept across several code versions is PROVISIONAL — trustworthy numbers need a coherent
   sweep on final code. Use a stale/mixed DB only as bug-finding signal, never as a verdict.
4. **MATERIAL gate for "dig-worthy":** `>1 unit AND ≥1% rel` — a ±1 straddle is never counted
   as a real divergence.
5. **Per-stand both-sides-trace to a named primitive; no unexplained divergence remains**
   (Pillar 4). Cluster-sampling finds bugs but does NOT discharge the endpoint — every
   non-bit-exact stand must end classified (fixed, cornered-to-a-named-primitive, or an FVS
   bug where jl is the correct side).
6. **Don't quote a pass-rate as final while unexplained divergence stands** — the eastern
   README calls that "PREMATURE." In-progress numbers are provisional; say so.

---
_The per-topic detail behind each rule lives in the `feedback-*.md` auto-memory notes._
