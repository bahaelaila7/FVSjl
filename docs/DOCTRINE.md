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

## 6. Commit & merge discipline
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

---
_The per-topic detail behind each rule lives in the `feedback-*.md` auto-memory notes._
