# ACTIVE GOAL — FVSjl modernization campaign

## Mission
Make **FVSjl** a *modern, maintainable, idiomatic-Julia, allocation-free-in-the-memory-path,
massively-parallelizable (Struct-of-Arrays)* implementation that is a **100% drop-in replacement for
FVS SN + NE + CS + LS** — *barring only* irreducible ULP-class divergences, FVS undefined behavior,
and known FVS output bugs (all already documented + cornered by the completed tolerance campaign).

The tolerance campaign is DONE (`docs/TOLERANCE_COMPLETE`): the suite is **bit-exact or
cornered-to-a-named-primitive**, all four variants bit-exact drop-ins on the exercised scenarios.
**That end state is now the CORRECTNESS FLOOR this campaign must never regress** — every refactor,
allocation-removal, SoA change, and parallelization must keep it.

## The FOUR PILLARS (each with a measurable done-state)

### 1. CORRECTNESS & COVERAGE — the "100% drop-in" claim, made defensible for ALL FOUR variants
- Every test stays BIT-EXACT (`==` / rendered-`==` vs live Fortran) or a documented `@test_broken`
  cornered to a named primitive. Suite must stay green (6885 pass / 128 broken / 0 fail baseline; a
  dropped pass or a new/undocumented broken = regression).
- **Close the coverage gap**: SN is deeply exercised; NE/CS/LS cover the core spine but are narrower.
  Broaden the NE/CS/LS scenario/keyword harness toward SN's breadth so "100% drop-in" is defensible
  for each variant — every NEW scenario validated BIT-EXACT-or-cornered vs the freshly-relinked live
  binary (`/tmp/FVS{sn,ne,cs,ls}_new`), not a stale golden.
- Done-state: a documented per-variant coverage matrix; NE/CS/LS materially broadened; all bit-exact-or-cornered.

### 2. ALLOCATION-FREE MEMORY PATH
- The per-cycle simulation hot path — `grow_cycle!` and everything it transitively calls (diameter/
  height growth, mortality, crown, volume, FFE fire/fuel/snag/carbon, summary accumulation) —
  allocates **ZERO** heap after warmup, or a documented+justified floor.
- Kill ad-hoc allocations: per-cycle `Dict`/temporary `Vector`/comprehensions in the hot path (e.g. the
  snag-binning `Dict`s in `book_mortality_snags!`, per-cycle `Float32[...]` temporaries). Replace with
  preallocated, reused `Scratch`/state buffers sized to `MAXTRE`.
- Done-state: `@allocated`/`BenchmarkTools` on a representative multi-cycle run per variant ≈ 0 bytes/cycle
  (documented floor), measured + recorded in the audit.

### 3. SoA + MASSIVE PARALLELISM
- Tree/snag/plot state is **Struct-of-Arrays** (columnar `Vector`s, preallocated to `MAXTRE`), not
  array-of-structs. Verify/complete the SoA layout; no per-record heap objects in the hot path.
- **No shared mutable global state**; per-stand RNG; per-stand scratch. The engine must be safe to run
  many stands concurrently.
- Demonstrate massive parallelism: a `Threads.@threads` (and/or `Distributed`) run over many stands is
  **BIT-IDENTICAL** to the serial run and scales. That bit-identity is the correctness test for parallelism.
- Done-state: a parallel multi-stand driver that is bit-identical to serial across all 4 variants; SoA
  verified; a documented "no shared mutable state" audit.

### 4. IDIOMATIC & MAINTAINABLE JULIA
- Hot-path functions are **type-stable** (`@code_warntype` / `JET` clean — no `Any`, no red flags).
- Idiomatic, consistent style; no dead code; clear module boundaries (`core`/`engine`/`io`/`variants/<v>`).
- Done-state: type-stability audit on the hot path clean; a style/dead-code pass recorded.

## THE DOCTRINE (imprinted — apply to EVERY change)
1. **BIT-EXACT IS THE FLOOR — NEVER REGRESS.** Re-run `julia --project=. test/runtests.jl` after every
   change. Baseline = **6885 pass / 128 broken / 0 fail / 0 error**. A dropped pass, a new fail, or a new/
   silently-flipped `@test_broken` = REVERT (or fix before proceeding). The 128 broken must stay the exact
   documented cornered set. All four variants stay bit-exact.
2. **FLOAT32 OP-ORDER IS SEMANTICS.** A refactor that reorders a Float32 sum/accumulation, changes an op
   sequence, or alters a loop/iteration order BREAKS bit-exactness (the tolerance campaign proved this
   over and over — non-associative sums, DGSCOR order, the cone-split, the species-sorted SDI sum).
   "Cleaner" / "more idiomatic" is NEVER worth a single bit. If an idiomatic rewrite changes the op order,
   either prove it bit-exact or don't do it.
3. **MEASURE, DON'T GUESS.** Find real allocations/hotspots with `@allocated`, `BenchmarkTools.@ballocated`,
   `Profile`, `--track-allocation=user`, `@code_warntype`/`JET`. No speculative optimization. Every
   allocation-removal must be PROVEN (bytes→0) AND bit-exact-preserving.
4. **VARIANT-SAFE — DO NOT HARDEN.** Shared `src/engine|io|core` gates on variant/coefficients; variant
   specifics under `src/variants/<v>/`. Every change keeps SN + NE + CS + LS bit-exact.
5. **PARALLELISM = NO SHARED MUTABLE STATE + PER-STAND RNG.** The bit-identity of serial-vs-parallel is the
   test. Global mutable state, shared scratch across stands, or a shared RNG = a bug to fix, not a feature.
6. **SoA FIRST.** Columnar `Vector`s preallocated to `MAXTRE`, reused, no per-record heap objects; enables
   both allocation-free and SIMD/parallel.
7. **INCREMENTAL + REVERSIBLE.** One slice at a time — validate (suite green + the pillar's metric) before
   the next. Never leave the tree broken across a stop.
8. **DOCUMENT EACH SLICE** in `docs/MODERNIZATION_AUDIT.md`: what changed, the before/after metric
   (allocated bytes / type-stability / parallel bit-identity / coverage added), and the bit-exact re-verify.
9. **RE-TRACE DISCIPLINE.** An allocation/type-instability/parallelism claim may be a misread — verify with
   the actual tool output, not intuition. (Prior campaign lesson: measurement beats labels every time.)

## Oracle & runner
- Live Fortran per variant: `/tmp/FVS{sn,ne,cs,ls}_new` (relink from `bin/FVS*_buildDir/*.o` + glibc shim,
  via `test/harness/*_oracle.sh`). Validate coverage additions vs the freshly-relinked binary.
- Suite runner: `julia --project=. test/runtests.jl` (NOT `Pkg.test`). Perf: `BenchmarkTools`, `@allocated`,
  `--track-allocation=user`, `Profile`; type: `@code_warntype` / `JET.@report_opt`.

## Working checklist
`docs/MODERNIZATION_AUDIT.md` — per-pillar task list + baseline metrics; tick each slice as it lands
(bit-exact re-verified + metric recorded).

## Off-switch
This reminder fires from the Stop hook until ALL FOUR pillars are met AND the suite is at the bit-exact-or-
cornered floor. To silence it (only when the audit checklist is 100% done, every slice bit-exact-verified):
`touch /workspace/FVSjl/docs/MODERNIZATION_COMPLETE`
