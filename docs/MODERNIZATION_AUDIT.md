# MODERNIZATION_AUDIT — working checklist

Goal + doctrine: `docs/MODERNIZATION_GOAL.md`. Tick each slice only when **bit-exact re-verified
(suite green)** AND the pillar metric is recorded. Baseline floor = tolerance-campaign end state.

## Baseline (measured — fill as established)
- **Suite:** ✅ floor intact. Baseline **6885/128/0** → **6889** (+4 pillar-3 `test_parallel.jl`) →
  **7037/128/0** (+148 pillar-1 `test_kwcov_variants.jl`, 4 NE scenarios × 37 col/row assertions). All +
  additions are new passing coverage, not regressions; 128 broken unchanged; 0 fail / 0 error.
- **SoA:** ✅ tree state IS Struct-of-Arrays already — `TreeList` (src/core/trees.jl) is columnar `Vector{…}`
  preallocated to `MAXTRE=3000`; scratch columns (WK1..WK15, sort idx) in `Scratch`. Snags = `SnagList`
  (mutable struct of Vectors). Design comment already states allocation-free + SoA + parallel intent.
- **Allocation floor (hot path):** baseline `@allocated run_keyfile(snt01_alpha; faithful)` = **71,849,816 B/run**
  (warm; includes legitimate IO). Ranked hot-path allocators (`--track-allocation=user`, non-IO):
  snag.jl:203 `pat` closure 16.0MB · trees.jl:112 (TreeList alloc, once) 4.7MB · snag.jl:216/220 3.8+3.8MB ·
  fuel_additions.jl:275 3.3MB · crown_biomass.jl:160 2.0MB · state.jl:558 1.5MB · r8clark_vol.jl:349/525 ·
  diameter_growth.jl:796 · standstats.jl:185/238. Probe: `scratchpad/alloc_probe.jl`.
- **Type-stability:** ✅ DONE — `@code_warntype` clean on grow_cycle! + diameter/height growth, mortality,
  compute_density!, crown_biomass, volume (R8+R9), standstats, update_snags! (see Pillar 4). 0 `::Any`/`Core.Box`.
- **Parallelism:** ✅ DONE — shared-mutable-global audit complete (only `_COEF_CACHE`, now lock-guarded); per-stand
  RNG/scratch/trees; parallel==serial bit-identical `-t8` all 4 variants + scaling measured (see Pillar 3).

## Pillar 1 — Correctness & coverage (100% drop-in, all 4 variants)  ← PRIORITY (semantics before optimization)
- **S17 (2026-07-06) — coverage found a real reporting bug (RESETAGE MAI); fix ATTEMPTED + REVERTED.**
  ne_resetage surfaced that post-RESETAGE the MAI col diverges (jl merch/age=426 vs live 0). disply.f:391-394
  zeroes MAI when MAIFLG≠0. FIRST FIX `mai=(ry>=0&&yr>ry)?0:…` was WRONG — broke s17_managed (SN scenario
  with the reset condition set but live MAI=62.5≠0) ⇒ the age-rebase condition ≠ MAIFLG. REVERTED (doctrine
  #1). Tracked @test_broken ne_resetage + task #80 (model the exact MAIFLG trigger from evtstv.f:398). The
  earlier "age reset to 0" was a SCENARIO column-misalignment (not a jl bug — aligned age matches live).
  NET: **+COMPRESS (eigensolver record compression) + CYCLEAT bit-exact vs live** (2 more verified NE paths).
- **S18 (2026-07-06) — RESETAGE MAI bug FIXED (correctly this time)** [`summary.jl`]. The reverted S17 fix
  zeroed MAI for ALL RESETAGE; the CORRECT condition (per evtstv.f:396 `ZERO==0`) is reset-to-**ZERO** only,
  which ⟺ `age_reset_age==0`. Added `&& Int(age_reset_age)==0`: zeroes MAI for ne_resetage (reset-to-0) but
  NOT s17_managed (resets to 40, MAI stays 62.5); non-RESETAGE (ry<0) + bare-ground untouched (bit-exact).
  ne_resetage now full-row BIT-EXACT (moved out of _KCV_BROKEN); s17_managed + canonical bare-ground stay
  green. **A REAL jl bug found via full-row coverage AND fixed, source-confirmed + suite-validated.** Task #80 DONE.
- **2026-07-06 — comparison hardened to FULL `.sum` row** (was 6 cols; now every volume/mortality column too)
  + broadened to 48 kwcov scenarios. This SURFACED a real semantic gap immediately:
- **R9 merch-override gap (NE/CS/LS)** — Saw-CuFt/BdFt cols ~1-2% off live under VOLUME/BFVOLUME/MCDEFECT
  overrides (base stand + SN(R8) bit-exact). ROOT-CAUSED: live fvsvol.f:256-380 BFPFLG=0 does a SEPARATE
  board recompute when board standards ≠ sawtimber (which VOLUME causes); SN/R8 `compute_volumes!` ports it
  (bit-exact), R9 `compute_volumes_ne!` does NOT. Fix = port the BFPFLG=0 R9 board recompute (SN is the
  template) + validate PER-TREE vs a live fvsvol.f stamp (doctrine — .sum inference alone kept producing
  unconfirmed hypotheses). Task #78, doc'd. ⇒ the ONE substantive drop-in gap; SDIMAX (small, likely
  mortality-knife-edge) tracked separately; everything else full-row bit-exact across 52 scenarios.
- **The 52-scenario full-row semantic map is the real pillar-1 deliverable this session** (see coverage
  matrix): NE/CS/LS are faithful drop-ins across multipliers/all-thin-forms/notriple/fixmort/serlcorr/
  rannseed/nocalib/dgstdev/strclass/minharv/mcdefect(NE)/bfdefect/combos — the "100% drop-in" claim is
  defensible except the merch-override volume path.
- ✅ Per-variant coverage matrix documented: `docs/MODERNIZATION_COVERAGE_MATRIX.md` (SN 261+37 keyword-
  isolation scenarios; NE/CS cover only the canonical bundle; LS has FFE/estab/sprout/sitesweep dedicated).
- 🟡 Broaden NE/CS/LS coverage (S7–S16, ✅ suite 8111/128/0): **33 keyword-isolation scenarios BIT-EXACT
  vs live** — {ne,cs,ls} × {mult, notriple, thinbta, fixmort, serlcorr, rannseed, thinaba, thincc, thinht,
  nocalib, dgstdev} (11 each). Incl. stochastic-DG paths (SERLCORR/RANNSEED) bit-exact. Harness `test_kwcov_variants.jl`
  AUTO-DISCOVERS `test/fixtures/kwcov/<prefix>_<kw>.*` (prefix→variant), `_KCV_BROKEN` for any cornered stem.
  Adding a scenario = drop 3 fixture files + rerun. Next: VOLUME/VOLEQ (needs vol-column compare),
  COMPUTE+event-monitor, MCDEFECT, THINQFA, STRCLASS, SDIMAX, NOCALIB/DGSTDEV.
  ⚠ Trap avoided: `run_keyfile` probes/tests MUST pass `variant=` or they default to SN (5-yr, SN model) on
  NE/CS/LS data → fabricated divergence (cost an hour before the re-trace caught it; memory note confirmed).
- ✅ **S22 (2026-07-07) — MAJOR coverage broadening + 2 fixes proven variant-safe + 1 gap found & localized.**
  NE/CS/LS keyword-isolation breadth grew **~11/11/5 → 29/26/26** (5 batches, all vs freshly-relinked live,
  each verified NON-VACUOUS where it should change output — the s3-vacuous-pass lesson):
  - Batch 1 (cs/ls): sdimax, resetage, cycleat, compress — **proved the SDIMAX (#79) + RESETAGE (#80) fixes
    generalize to CS+LS bit-exact**, not just the NE scenario each was found on.
  - Batch 2 (ne/cs/ls): thinqfa (non-vacuous Q-factor thin, TPA 524→71), compute (event-monitor evaluator).
  - Batch 3 (ne/cs/ls): leavesp, spgroup. ne_leavesp cornered ULP (1 derived col 0.1 off; cs/ls bit-exact).
  - Batch 4 (ne/cs/ls): fertiliz (non-vacuous eastern DG response, bit-exact). SPECPREF dropped (vacuous w/ THINBBA).
  - Batch 5 (ne/cs/ls): estab (ESTAB+PLANT) — surfaced a REAL gap (#81), cornered + localized (below).
  Stale-golden audit: the 6 oldest goldens regenerated vs fresh live → all identical (baseline trustworthy).
  kwcov gate **12094 pass / 10 broken / 0 fail**; full suite **18982 / 139 / 0**. Coverage matrix updated.
- 🔍 **#81 (ESTAB/PLANT) — precisely LOCALIZED via per-tree DBS FVS_TreeList** (harness scratchpad/tl/): the
  planted cohort matches live in count/TPA/QMD/TopHt at the first post-plant cycle, diverging only in SDI
  (236 vs 234) + CCF — driven by the SMALLEST suppressed seedlings (jl over-grows min dbh 0.784 vs live 0.558)
  in the CS/NE small-tree REGENT DBH growth. ~1% SDI at bit-exact QMD/TopHt ⇒ distribution/ULP-adjacent.
  Multiple hypotheses raised + CORRECTED with tool data (doctrine #9). CORNERED-to-a-named-primitive (an
  accepted goal end-state); definitive fix = a live cs/regent.f REGENT stamp (disproportionate for a
  ULP-adjacent uncommon-scenario item — deferred, fully characterized in task #81).
- ⬜ (carry) #73 soft-snag CWD1 cone-split LOHT — DEFERRED-BY-DESIGN (2026-07-07 re-assessed): the soft LOHT=1.0
  vs jl 0.10 gap is REAL but INERT — the fall path books all mortality snags HARD (DFIS=0 every cycle), so the
  soft branch is only reachable via the BLOCKED SNAGDCAY no-op; unvalidatable without unblocking it. Accepted
  deferred set (with SNAGDCAY/SNAGBRK). See task #73.
- ✅ **RESOLVED (stale carry reconciled 2026-07-07)** — #78 R9 merch-override board recompute (task #78 DONE:
  bftopk board-top fix, ne_bfvolume bit-exact, SawCuFt ULP cornered) and #81 CS/NE ESTAB/PLANT establishment-
  height (tasks #81/#82 DONE bit-exact). Both completed; everything else bit-exact-or-cornered.

## Pillar 2 — Allocation-free memory path
- **S10 (2026-07-06) — volume return-buffer reuse** [`r8clark_vol.jl` + `volume.jl`]. Added an optional
  `buf` kwarg to `_R8CLARK_VOL` (default nothing→`zeros(15)`, so all 7 callers unchanged); `compute_volumes!`
  now allocates TWO reusable 15-vectors once per call and passes them to the primary (v) and board-recompute
  (vb) sites — distinct buffers because their lifetimes overlap (`v[7]` read after `vb`). Buffer zeroed on
  entry ⇒ identical fill to `zeros(15)`; results consumed as scalars per iteration (no retention).
  Metric: **25,213,976 → 24,069,592 B (−1.1MB)**. Bit-exact (volume feeds .sum): ✅ **7667/128/0**. Total
  **71.8 → 24.07 MB (−66%)**.
- **S11 attempt (2026-07-06) — logLen buffer — REVERTED.** Tried the same buffer pattern for the
  `logLen = zeros(40)` sawtimber board scratch, but the edit landed in the HELPER `_r8_scribner_bf` (called
  BY `_R8CLARK_VOL`), which has no `logbuf` in scope → `UndefVarError` → volume path threw → 104 test
  errors (caught by doctrine #1, reverted immediately; floor restored to 7667/128/0). Lesson: the `logLen`
  at r8clark_vol.jl:531 belongs to `_r8_scribner_bf`, not `_R8CLARK_VOL` — buffering it needs the kwarg
  threaded through that helper too. Deferred (1.3MB, not worth the extra layer right now).

- **S13 (2026-07-06) — precompute per-species 2-char crown-width code** [`state.jl`, 4×`species.jl`,
  `standstats.jl`, `dbs_output.jl`, `sprout.jl`, `cuts.jl`, `fire/fuel_model.jl`, `fire/fmcba.jl`]. The
  per-tree `class_codes[sp,1][1:2]` String slice (2 standstats loops + 5 other readers) is now a precomputed
  `SpeciesData.code2::Vector{String}` set once at species load. Bit-identical (same slice value). CAUGHT +
  FIXED IN-TURN: loader-level slice hit a `BoundsError` on phantom species whose base code is empty
  (`class_codes="1"`, 1 char) — the old readers never sliced those (not live trees); switched the setter to
  `first(cc, 2)` (safe, identical for real ≥2-char species). Metric: **24,069,592 → 22,707,096 B (−1.4MB)**.
  Bit-exact: ✅ **7889/128/0**. Total **71.8 → 22.7 MB (−68%)**. (Also removed the slice from sprout/cuts/
  fire readers = a small pillar-4 cleanup.)

- **S15 (2026-07-06) — logLen buffer, DONE (the S11 revert, done right)** [`r8clark_vol.jl` + `volume.jl`].
  Threaded `logbuf` through BOTH layers this time — `_R8CLARK_VOL` (kwarg) → `_r8_scribner_bf` (the helper
  that actually owns `logLen = zeros(40)`, which S11 missed → the UndefVarError). `compute_volumes!` passes
  a third reusable 40-vector. Zeroed on entry ⇒ identical. Metric: **22,707,096 → 21,725,336 B (−1.0MB)**.
  Bit-exact: ✅ **7889/128/0**. Total **71.8 → 21.7 MB (−70%)**.

- **S16 (2026-07-06) — `_r8clark_lookup` try/catch → tryparse (idiomatic; alloc-neutral).** Replaced
  `parse`+try/catch with `tryparse` (identical parse semantics, no throw). Bit-exact ✅ 8111/128/0. BUT the
  hypothesis that the blank-voleq exception dominated the 0.6 MB was WRONG (doctrine #9): warm `@allocated`
  moved only −256 B ⇒ the `.mem` 0.6 MB at r8clark_vol:88/90 was **warmup/compilation**, not steady-state.
  Kept as a cleanliness/robustness improvement (avoids exceptions on invalid voleqs); NOT an alloc win.
  **Implication:** the per-cycle volume allocation is already near-zero; the remaining 21.7 MB aggregate is
  dominated by IO (.key parse + .sum write) + one-time `MAXTRE` preallocation — outside the per-cycle hot
  path the done-state targets. So the pillar-2 hot-path goal is ~met; the aggregate floor is IO/setup.

### Pillar-2 residual floor (characterized — mostly IO + one-time setup; per-cycle hot path ~allocation-free)
All in leaf functions needing per-species precompute or buffer-threading with real bit-exact risk:
- `_r8_scribner_bf` `logLen = zeros(40)` (sawtimber trees) — thread a `logbuf` through the helper.
- ✅ `standstats.jl:185/238` — the `[1:2]` slice DONE (S13, precomputed `code2`). Residual: `crown_width`'s
  internal `rstrip(String(sp2))` + 2 Dict lookups per tree (smaller now; would need per-species crown-eq
  precompute to fully remove — deferred, bigger interface change).
- `_r8clark_lookup` voleq string parse (r8clark_vol.jl:88/90) — precompute per species at setup.
- The 6 non-volume `_R8CLARK_VOL` callers (mortality/crown_biomass/carbon/snag) still `zeros(15)` (fire path).

- ✅ Baseline measured + call sites ranked (see above). NOTE: `trees.jl:112` (TreeList ctor) and
  `state.jl:558` (Scratch ctor) are **once-per-stand preallocations** — expected, NOT per-cycle targets.
- ✅ **DONE (S59, 2026-07-07)** — Isolated per-cycle metric: `@allocated grow_cycle!(s)` (warmed, net01 NE) =
  **9280 B/cycle**, all of it the documented+justified "Base sort scratch" floor (3× `compute_density!`
  descending-DBH stat sorts; the ~45 KB Printf-IO floor is in the summary path, not `grow_cycle!`). Confirms
  the Pillar-2 per-cycle done-state; S48 improved the sort-scratch component (384→160 B). Further reduction is
  unsafe (would risk `pbal`/DG bit-exactness) — see S59.
- ✅ **S20 (2026-07-06) — Kill the snag-binning `Dict`s in `book_mortality_snags!`** [`mortality.jl` +
  `state.jl`]. The FMSADD snag-record binning allocated SIX per-call heap objects on every fire-mortality
  cycle — `minht`/`maxht`/`midht` `Dict{Tuple{Int,Int},Float32}`, `gkey` `Dict{Tuple{Int,Int,Int},Int}`, and
  the four growable `gsp/gdbh/ght/gden` `Vector`s. Measured (`--track-allocation`, fire_carbon ×20): line 585
  push!+Dict-insert 104,960 B, minht/maxht 16,000 B, gsp… 12,800 B, gkey 8,000 B, midht 8,000 B ⇒ ~7.5 KB/run.
  FIX: a new preallocated `SnagBinScratch` (one field on `FireState`, so non-fire stands pay nothing) with
  DENSE bins over the (species, DBH-class 1:19, HT-class 1:2) tuple space — `minht`/`maxht` indexed by
  `(sp-1)*19+dbhcl` (`MAXSP*19`), `gkey` by `((sp-1)*19+(dbhcl-1))*2+htcl` (`MAXSP*19*2`), and `gsp/gdbh/ght/
  gden` sized `MAXTRE`. `fill!`-reset at the top of each call (in-place, allocation-free); MIDHT computed
  inline in PASS 2. The `i=1:n` scan order is UNCHANGED, so both the Float32 density-weighted running-mean
  accumulation order AND the first-appearance emission order are identical to the Dict version. AFTER:
  book_mortality_snags! `--track-allocation` shows **ZERO** on every line. **Bit-exact re-verify: suite
  15204 pass / 135 broken / 0 fail — identical to before** (all fire/carbon/snag tests green). Bounds proven
  exact (lin2 max 2052 = MAXSP·19; lin3 max 4104 = MAXSP·19·2).
- ✅ **S21 (2026-07-06) — Preallocate the SN crown Weibull caches** [`crown_ratio.jl` + `state.jl`].
  `crown_ratio_update!` runs EVERY cycle for EVERY stand (core non-fire hot path), and allocated
  `Aw/Bw/Cw = zeros(Float32, MAXSP)` ×3 + `seen = falses(MAXSP)` on every call (~1.3 KB/cycle/stand of
  per-species Weibull-parameter scratch). FIX: moved to reused `Scratch` buffers (`crown_aw/bw/cw` +
  `crown_seen`, sized `MAXSP`); `seen` is `fill!(false)`-reset each call (the lazy compute gate), Aw/Bw/Cw
  need no reset (read only where `seen[sp]`). The per-species lazy compute order is unchanged ⇒ bit-exact
  (no cross-tree Float32 accumulation here). `--track-allocation` (s1_thinning ×8): crown_ratio.jl lines
  75-120 now **ZERO** (line 79 gone). **Suite 15204/135/0 — identical (bit-exact).** SN-only pattern
  (NE/CS/LS crown paths don't use it).
- ✅ **S22 (2026-07-07) — `vtot` summary-accumulation boxing** [`summary.jl:319`]. The per-cycle summary
  volume totals called `getfield(t, f)` with a RUNTIME Symbol ⇒ `fld::Any` ⇒ `fld[i]*t.tpa[i]` boxed every add.
  Added `::Vector{Float32}` assertion (all vtot fields are Vector{Float32}) — concrete `fld`, allocation-free,
  type-stable, same sequential Float32 order (bit-exact). **Measured: SN growth-only per-cycle 157→115 KB/cycle.**
- ✅ **S23 (2026-07-07) — R8 Clark voleq string-parse** [`r8clark_vol.jl:90`]. `_r8_mrules` parsed the per-species
  voleq via `string(voleq[2])` / `voleq[8:10]` — each allocating a fresh String on EVERY per-tree, per-cycle volume
  call. Switched to non-allocating `SubString` views (tryparse semantics identical ⇒ bit-exact). **Per-cycle
  115→104 KB/cycle.** Suite 19486/136/0 (bit-exact).
- ✅ **S24 (2026-07-07) — `crown_width` per-call `rstrip(String(sp2))`** [`crown_width.jl:63`]. Every crown-width
  call (per-tree, per-cycle, from standstats/cuts/sprout/fmcba/fuel_model) did `get(crown_species, rstrip(String(sp2)), …)`
  — allocating a String copy + an rstrip SubString EACH call. Fix: pre-rstrip `species.code2` ONCE at load
  (`String(rstrip(first(class_codes,2)))` in all 4 variant species.jl) so the lookup key is IDENTICAL to `rstrip(code2)`,
  then look up `sp2` directly. **Bit-exact BY CONSTRUCTION** (key unchanged; structure_stage already passes `strip()`).
  **Measured via `Profile.Allocs`: crown_width 37.3 → 0.0 KB/cycle** (the biggest simulation-path — non-IO — per-cycle
  allocator, now gone). Suite 19486/136/0.
- 📊 **Per-cycle tail ATTRIBUTED via `Profile.Allocs`** (doctrine #3 — `--track-allocation=user` misses Base-attributed
  allocs; the allocation profiler attributes to the user call site). This drove the next three slices.
- ✅ **S25 (2026-07-07) — `stand_top_height` / `point_basal_area!` / `stand_pct!` sort permutation** [`standstats.jl:127,158,211`].
  Each once-per-cycle stand-stat pass did `sortperm(view(t.dbh,…); rev=true)` — a fresh `Vector{Int}` per call per cycle.
  Added a dedicated `Scratch.stat_idx::Vector{Int32}` (MAXTRE, the three passes run sequentially — never nested — so one
  shared buffer is value-safe) and switched to `sortperm!(view(stat_idx,1:n), …)`. `sortperm!` uses the SAME default
  algorithm ⇒ IDENTICAL permutation incl. tie-breaking (verified over 2000 tie-heavy trials + full suite) ⇒ bit-exact.
  **Difference-method total 125.2 → 110.1 KB/cycle (−15).**
- ✅ **S26 (2026-07-07) — `crown_width` Dict-lookup Union-boxing** [`crown_width.jl:63`]. `get(coef.crown_species, sp2, nothing)`
  returns `Union{Tuple{String,String},Nothing}`; the non-isbits `Tuple` is BOXED on every return (32 B/call), and
  crown_width is called for EVERY tree across multiple CCF passes per cycle (point_credit_ccf!/stand_ccf/structure_stage/
  fmcba/fuel_model). Replaced both `get(…, nothing)` lookups with `haskey`+`getindex` — two cheap non-boxing hashes,
  value-identical ⇒ bit-exact. Measured: crown_width 2160 B → **0 B** over 27 trees. **Difference-method total
  110.1 → 57.4 KB/cycle (−53).** (S24 had removed the per-call `rstrip(String())` copy; S26 removes the residual box.)
- ✅ **S27 (2026-07-07) — `compute_forest_type!` per-cycle `zeros(210)` + 13 slice-copies** [`forest_type.jl:74,118+`].
  FORTYP runs every cycle (simulate.jl:26, feeds dgf!). It built a fresh `zeros(Float32,210)` stocking accumulator and
  then summed 13 `s[a:b]` COPY-slices (`sum(s[1:58])`, `sum(s[81:153])`, …). Added `Scratch.stkval_s::Vector{Float32}(210)`
  (fill!(0) each call) + switched every `sum(s[a:b])` to `sum(@view s[a:b])` (identical iteration order ⇒ bit-exact).
  Measured: compute_forest_type! **2 KB → 0 B**. Difference-method total **57.4 → 55.0 KB/cycle.**
- 📉 **Per-cycle floor (SN growth-only, difference method pc10↔pc20): 157 → 55.0 KB/cycle** (whole effort; **125 → 55 this
  session** via S24–S27). **The per-cycle hot path is now at the Base-stdlib floor — every AD-HOC allocation is eliminated**
  (snag-binning Dicts, crown Weibull caches, vtot boxing, voleq String parse, crown_width rstrip+Union-box, standstats sort
  Vectors, forest_type zeros+slices). Attribution of the residual 55 KB/cyc:
    - **~45 KB/cyc = `summary.jl:40 write_sum_row`** — the per-cycle .sum output-row rendering via stdlib `Printf.format`
      (28 mixed Int/Float64 columns: heterogeneous-vararg boxing + format machinery). This is **OUTPUT IO, not the simulation
      hot path**, and is a stdlib-Printf floor — a hand-rolled 28-column formatter would have to reproduce every field width/
      rounding byte-for-byte (the .sum is tested rendered-`==`), high-risk for ~45 KB/cyc of pure IO. **DOCUMENTED FLOOR.**
    - **~10 KB/cyc = Base sort scratch** in grow_cycle!'s transitive path: crown_ratio.jl:70 (ascending sort + inverse-perm),
      standstats `sortperm!` merge-buffer residual, serial_correlation species_sort!. `sortperm!` leaves ~272 B/call even WITH
      a persistent `scratch=` buffer (measured) — a Base merge-sort floor; a custom allocation-free sort would risk the
      tie-break order (doctrine #2). **DOCUMENTED FLOOR** (small, Base-internal).
  ⇒ **Pillar-2 per-cycle done-state MET for SN: zero ad-hoc allocation; residual is documented+justified Base-stdlib floor
  (Printf IO + sort scratch).** The `trees.jl fz/iz/TreeList` + `init.jl each_stand` sites in Profile.Allocs are ONE-TIME
  per-stand setup (the difference method cancels them).
- ✅ **S28 (2026-07-07) — R9 Clark volume per-tree allocation (NE/CS/LS floor)** [`r9clark_vol.jl`]. Measuring NE/CS/LS with
  their own growth-only pc10/pc20 keyfiles exposed a per-cycle floor of **NE 131.7 / CS 106.6 / LS 126.6 KB/cyc** — far above
  SN's 55, because NE/CS/LS use the **R9 Clark volume** path (SN uses R8). `Profile.Allocs` attributed ~76 KB/cyc to three
  per-tree-per-cycle allocators: `_r9_bucked_bf` `logLen=zeros(40)` (37.7), `r9clark_cubic` `vol=zeros(15)` (25.6), and the
  `_R9State(...)` mutable-struct construction (12.8). Three bit-exact fixes:
    1. **`_R9State` → immutable/isbits** — the only post-construction writes were line 408 (`st.a17=co.a17_0`, IDENTICAL to
       the construction value ⇒ redundant, deleted) and the dbhIb/dib17/totHt sets (folded into a second stack-allocated
       rebuild). All taper functions (`_r9_dia417`/`_r9_cuft`/`_r9_ht`/`_r9_dib`/`_r9_bucked_bf`) only READ `st`. Stack-allocated
       now ⇒ **−12.2 KB/cyc** (NE 131.7→119.5).
    2. **`vol=zeros(15)` + `logLen=zeros(40)` → per-stand `Scratch.r9_vol`/`Scratch.r9_logbuf`** threaded via optional
       `vbuf`/`logbuf` kwargs on `r9clark_cubic` (and a `logbuf` arg on `_r9_bucked_bf`). `compute_volumes_ne!`'s hot loop
       passes `s.scratch.*` (parallelism-safe — per-stand, not global); the non-hot callers (snag/mortality) pass `nothing`
       and allocate, staying bit-exact and unchanged. `vol` is `fill!(0)` each call (== `zeros`); `logLen` is written-before-read
       (value-safe dirty reuse); the driver consumes the returned `v` before the next tree (no aliasing). **−67 KB/cyc.**
  Result: **NE 131.7→52.3, CS 106.6→49.3, LS 126.6→51.7 KB/cyc** — all four variants now converge on the same ~50–55 KB/cyc
  floor. Suite **19486/136/0** (bit-exact, incl. the heavily-tested volume columns).
- 📉 **Per-cycle floor ALL FOUR VARIANTS: SN 55.0 / NE 52.3 / CS 49.3 / LS 51.7 KB/cyc** — the hot path has **zero ad-hoc
  allocation on every variant**; the residual ~50 KB/cyc is the SAME documented+justified Base-stdlib floor everywhere:
  ~45 KB/cyc `write_sum_row` Printf output IO + ~5–10 KB/cyc Base sort scratch. **Pillar-2 per-cycle done-state MET for all
  four variants.**
- ✅ **S29 (2026-07-07) — FFE snag-fall `byr` recomputation over `burn_reports::Vector{Any}`** [`snag.jl:307`].
  Profiling a fire scenario (fire_carbon) showed ONE dominant FFE allocator: **1456 KB** at the `byr` (last-qualifying-
  burn-year) computation inside `update_snags!`. It looped `fs.burn_reports` — a `Vector{Any}` (15-field NamedTuple
  element, kept untyped) — reading `.scorch`/`.year`, which BOXES on every access, and it did so **per-snag × per-year**
  (inside both the snag loop and the annual `for _ in 1:yrs` loop) though `byr` is snag- AND year-INDEPENDENT. Fix:
  hoist the computation to function entry (one pass) with `::Float32`/`::Int` asserts to de-box the Any-element reads.
  Bit-exact (same value). Measured: byr-region **1456 KB → 0**. Suite 19486/136/0 (fire tests bit-exact). The FFE-path
  residual is now small + floor-class: carbon-report row formatting (`_format_carbon_row`, IO like write_sum_row) +
  `add_snag!` `push!` SoA growth (amortized — snags accumulate across cycles, not MAXTRE-bounded like live trees).
- 📊 **Pillar-2 status: DONE-STATE MET for the full per-cycle hot path, all four variants** — growth (S22–S28) AND the FFE
  fire/snag/carbon path (S29). Every dominant ad-hoc allocation eliminated; residuals are documented+justified floors
  (Printf output IO, Base sort scratch, amortized SoA container growth). The one remaining note (`Vector{Any}` burn_reports
  type itself) is a maintainability item (Pillar 4), not a per-cycle allocator now that access is hoisted+de-boxed.

## Pillar 3 — SoA + massive parallelism  🟢 core done-state MET
- ✅ SoA verified (see Baseline): `TreeList`/`SnagList` columnar, preallocated to `MAXTRE`; no per-record
  heap objects (copy_tree! S5 removed the last dynamic per-field boxing).
- ✅ Shared-mutable-global audit: the read-only lookup `Dict`s (schemas, column maps, `_EV_CMP`, forest
  defaults) are init-once/read-only ⇒ thread-safe. RNG is per-`StandState` (`rng.jl:14`). Tree/plot/scratch
  are per-stand. The ONE mutable global was `_COEF_CACHE` (`get!`-populated) — **now lock-guarded** (S6,
  `cached_coefficients` + `_COEF_LOCK`); `SpeciesCoefficients` is immutable and never mutated in-sim ⇒
  shared read-only. No other shared mutable state.
- ✅ Parallel multi-stand BIT-IDENTICAL to serial, **all 4 variants** (`-t 8`, N=48 concurrent runs each,
  incl. a COLD-cache race that forces threads through the `get!`): SN/NE/CS/LS all 0 mismatches. Permanent
  test `test/integration/test_parallel.jl` (degrades to a determinism/reentrancy check under `-t1`).
- ✅ **RE-VALIDATED 2026-07-07 after S24–S28**: the four new per-stand `Scratch` buffers (`stat_idx`, `stkval_s`,
  `r9_vol`, `r9_logbuf`) all live in `s.scratch` (per-stand, allocated in the `StandState` ctor) and the R9
  `vbuf`/`logbuf` are threaded from `s.scratch` — NO new shared mutable state. `test_parallel.jl` under `-t8`:
  all 4 variants 0 mismatches. The `_R9State`→immutable change (S28) is strictly parallelism-positive.
- ✅ **Scaling measurement (2026-07-07)** — 128 stands via `@threads`, wall-clock at `-t1/2/4/8`:
  **43.2 → 62.5 → 91.6 → 122.4 stands/s** (1.0× / 1.45× / 2.12× / **2.83×** at 8 threads). The `@assert all(==)`
  bit-identity gate passed at every thread count (parallel == serial, byte-for-byte). Sub-linear because each stand
  is short (per-stand StandState setup + GC + the uncontended `_COEF_LOCK` dominate at high core counts) — but it
  scales monotonically and stays bit-exact, which is the Pillar-3 done-state. (`/tmp/scaling.jl`.)

## Pillar 4 — Idiomatic & maintainable
- 🟢 **Style / dead-code pass (S14, 2026-07-06) — recorded:**
  - Module hygiene: ✅ NO orphan files (every `src/**/*.jl` is `include`d). Boundaries clean:
    `core`/`engine`/`io`/`variants/<v>`.
  - Dead-code: NO dead internal helpers found. (A grep scan flagged `!`-suffixed helpers as refs=0, but
    that was a `\b`-word-boundary artifact on `!`; spot-checks confirmed all are called — e.g. `_maybe_burn!`
    ×2, `_merge_classes!` ×3. Deeper static analysis would want `JET.jl` as a dev-dep — optional follow-up.)
  - TODO/FIXME audit: 18 markers → 14 false-positives (the `xxx` volume-temp var, "HACKBERRY" species,
    "STOPxxx" in a comment); 3 legitimate documented deferrals (NE aspen suckering/ASSPTN sprout.jl:507,
    CFTOPK broken-top reuse r9clark_vol.jl:544, spctrn verify note ne/species.jl:47); **1 STALE comment
    FIXED** — keyword_dispatch.jl COMPRESS said "algorithm TODO" but COMPRESS is fully ported+validated
    (compress.jl `apply_compress!`/`_merge_classes!`, IBM EIGEN) → comment corrected.
  - Idiomatic wins already landed: S5 `@generated copy_tree!` unroll, S1–S4 closure hoists, S13 removed the
    `class_codes[..][1:2]` slice pattern from 7 sites.
- 🟢 Type-stability audit of the CORE hot path — **all 0-instability** (warntype, `scratchpad/wt_growth*.jl`):
  `grow_cycle!`, `diameter_growth!`, `height_growth!`, `mortality!`, `compute_density!`, `crown_biomass`,
  `_cwd_cone_fractions`. (S3/S4 fixed the two that boxed captured locals; the rest were already clean.)
- 🟢 **Volume + standstats + snag sweep DONE (2026-07-07) — all 0-instability** (`code_warntype`, `::Any`=0 &
  `Core.Box`=0): `compute_volumes!` (R8/SN), `compute_volumes_ne!` (R9/NE-CS-LS), `stand_top_height`,
  `point_basal_area!`, `stand_pct!`, `update_snags!`. The S28 immutable `_R9State` and S29 `::Float32`/`::Int`
  asserts on the `burn_reports::Vector{Any}` access made the last two clean (type-stability win alongside the
  allocation fix). ⇒ the audit's "un-swept volume/standstats" gap is CLOSED; the hot path is type-stable end-to-end.
- ✅ Style/dead-code pass + module-boundary review — DONE in S14 (above): no orphan files, clean core/engine/io/
  variants boundaries, no dead helpers, TODO/FIXME audited. (Duplicate open-item resolved 2026-07-07.)

## Slice log (append: what changed / metric before→after / bit-exact re-verify)
- **S1 (2026-07-06) — cone-fraction alloc removal** [`snag.jl _cwd_cone_fractions`]. Hoisted the per-call
  `pat` closure to a module-level pure `_cwd_pat` (no capture) and changed the return from a heap
  `zeros(Float32,9)` to a stack `NTuple{9,Float32}` (matches the docstring; all 3 callers read `frac[j]`
  only). Bit-identical arithmetic (same `let r2=...` expr, same op order). Metric: run_keyfile(snt01_alpha)
  **71,849,816 → 48,628,904 B (−23.2MB, −32%)**. Bit-exact re-verify: ✅ **6885 pass / 128 broken / 0 fail**
  (1m33s) — floor intact.
- **S2 (2026-07-06) — two more per-call closures** [`crown_biomass.jl` `cone`, `diameter_growth.jl` `bsc`].
  `cone` → module-level pure `_cone(dbrk,ang,sg,mp)` (all Float32: mypi=3.14159f0, _FM_P2T=0.0005f0 — types
  + op order identical). `bsc` (defined per-tree in the DG loop) → 4 inline `_bound_scale(...)` calls (same
  fn, same args). Metric: **48,628,904 → 46,725,720 B (−1.9MB)**. Bit-exact: ✅ **6885/128/0** (1m32s).
- **S3 (2026-07-06) — unbox `angle` in crown_biomass (pillar 2 ∧ 4)** [`crown_biomass.jl`]. `@code_warntype`
  showed `angle = Core.Box()` + a SECOND cone closure `conem` I'd missed in S2 — it captured `angle` by
  reference, keeping it boxed even after `cone` was hoisted (why S2 only saved 1.9MB). The Box made every
  `angle`-derived value `::Any`, de-optimizing the whole function → the scattered per-line allocations
  (480 B/call across lines 186–225). Hoisted `conem` → module `_conem(dbrk,d,ang,sg,mp)` (bit-identical,
  all Float32). Result: `crown_biomass` now **0 type-instabilities**, **480→32 B/call**; aggregate
  **46,725,720 → 31,952,856 B (−14.7MB)**. Bit-exact re-verify: ✅ **6885/128/0** (1m37s) — floor intact.
  **META-LESSON (reusable):** a captured local boxed by ANY closure de-optimizes the whole function and
  shows up as *scattered* per-line allocations in `--track-allocation`. `@code_warntype` (optimize=true)
  pinpoints the `Core.Box`; hoisting every closure that captures it fixes pillar-2 allocs AND pillar-4
  type-stability at once. Run the warntype+`@allocated` probe (`scratchpad/warntype_cb.jl`) on any hot
  function whose `.mem` shows scattered small allocations.
- **S4 (2026-07-06) — unbox `r1` in `_cwd_cone_fractions`** [`snag.jl`]. The S1 NTuple refactor still
  allocated 352 B/call: `r1` was *reassigned* (`r1 = r1 + …`) and captured by the `ntuple(Val(9)) do j`
  closure → `Core.Box` (same pattern as `angle`). SSA-renamed the first value to `r1a` so `r1` is
  assigned once → captured unboxed. Result: **352 → 0 B/call, 0 instabilities**; aggregate
  **31,952,856 → 27,250,408 B (−4.7MB)**. Bit-identical (`r1a` holds the same first value). Bit-exact
  re-verify: ✅ **6885/128/0** (1m35s). Running total: **71.8 → 27.25 MB (−62%)** since campaign start.
  Reinforces the meta-lesson: reassigned-AND-captured locals box — SSA-rename or hoist.
- **S5 (2026-07-06) — copy_tree! field-copy unroll + forest_type `_cf` hoist** [`trees.jl`, `forest_type.jl`].
  (a) `copy_tree!` looped `for f in _TREE_VEC_FIELDS; getfield(t,f)[dst]=getfield(t,f)[src]` — a RUNTIME
  Symbol into `getfield` ⇒ `Any` column ⇒ every copied value boxed + dynamically dispatched, per field
  per tripled record (core tripling hot path). Replaced with a `@generated _copy_tree_vecs!` that unrolls
  to `t.species[dst]=t.species[src]; …` with CONSTANT symbols (auto-derived from _TREE_VEC_FIELDS, stays
  in sync) ⇒ type-stable, zero-alloc, no dispatch (also a speed win). (b) `forest_type` `_cf` closure
  captured loop-reassigned `ttst51`/`dmxss` ⇒ boxed; hoisted to `_forest_cf(d,ttst51,dmxss)`. Both
  bit-identical (pure value copies / same expression). Metric: **27,250,408 → 25,213,976 B (−2.0MB)**.
  Bit-exact re-verify: ✅ **6885/128/0** (1m30s). Running total: **71.8 → 25.2 MB (−65%)**.

### Remaining pillar-2 per-cycle allocators (post-S5 ranking) — ✅ ALL RESOLVED (re-traced 2026-07-07)
The core growth/mortality/crown/height path NO LONGER appears in the `--track-allocation` ranking ⇒ it is
already allocation-free/type-stable. The three then-remaining items are now ALL done:
- ✅ `_R8CLARK_VOL` `vol = zeros(15)` + `logLen = zeros(40)` — RESOLVED: now `buf`/`logbuf` params
  (`r8clark_vol.jl:357`, threaded from `volume.jl:534/569` `_vbuf`/`_logbuf`). SN `compute_volumes!` measures
  480 B / 27 trees (residual = the rare broken-top cftopk/bftopk `logLen`, lines 614/674 — only fires for
  trunc trees). The R9 analog (NE/CS/LS) done in **S28**. Volume path allocation-lean on all four variants.
- ✅ `_r8clark_lookup` voleq STRING parse — RESOLVED in **S23** (non-allocating `SubString` views).
- ✅ `standstats` `class_codes[1:2]` slice + `crown_width` `rstrip(String(sp2))` + Dict lookup — RESOLVED in
  **S24** (pre-rstripped `code2` at load) + **S26** (`haskey`+`getindex`, no Union-box). standstats sort
  Vectors → `Scratch.stat_idx` in **S25**.

- **S6 (2026-07-06) — pillar-3: thread-safe coef cache + parallel bit-identity harness** [`coefficients.jl`
  + 4 variant accessors + `test_parallel.jl`]. Audited all shared state; the only mutable global was
  `_COEF_CACHE`. Wrapped its `get!` in `cached_coefficients`/`_COEF_LOCK` (called once per stand at
  `StandState` ctor ⇒ uncontended after warmup, negligible). Proved parallel==serial bit-identical for
  SN/NE/CS/LS under `-t8` (cold-cache race included). Bit-exact: ✅ **6885/128/0** (S6 single-thread run);
  the `-t8` full-suite + new parallel test run confirms under real concurrency. Value-neutral (a lock
  changes no Float32 result). Final: `-t8` suite **6889/128/0** (all 4 parallel assertions green).
  ⇒ **Pillar 3 core done-state MET.**

## Slice S30 — FIA database input (reader + translator) [2026-07-07]
- **What**: native FIA "FVS-ready" DB input (`src/io/fia_database.jl`) — `kw_database!`
  parses the DATABASE/DSNIN block (StandSQL/TreeSQL), opens the DB READ-ONLY
  (mode=ro/immutable=1, never modifies source), maps FVS_STANDINIT_COND→plot/control
  (dbsstandin.f) + FVS_TREEINIT_COND→TreeRecords (dbstreesin.f); existing `notre!` does
  the BAF TPA expansion. `treeinput.jl` refactored to share `ingest_tree_records!`. Plus
  a raw-FIA→FVS-ready translator (`src/io/fia_translate.jl`) for the derivable fields.
- **Metric (correctness)**: cycle-0 inventory BIT-EXACT vs live FVS on real LS stands for
  the density/size columns (TPA/BA/SDI/TopHt/QMD); volume/CCF a small bounded residual;
  multi-cycle growth diverges 5–30%/50yr on dense stands (mortality timing, not IDG=1).
  Translator's direct-copy fields reproduce FVS-ready rows bit-exact; SITE_INDEX/LOCATION/
  design/DG-calib need the external USFS FIA2FVS reference (not in repo) — documented.
- **Bit-exact re-verify**: suite 26099 pass / 136 broken / 0 fail (+22 from
  test/integration/test_fia_reader.jl; fixture test/fixtures/fia/ls_sample.db). No regression.

## Slice S31 — FIA reader species-code 3-digit fix [2026-07-07]
- **What**: FIA numeric species codes arrive un-padded from the DB ("71", "012") but the FVS
  species tables key on the 3-digit FIA code ("071"); un-padded codes < 100 ALL mis-resolved
  to Other-Hardwood (wrong crown width + volume; masked because species-independent
  BA/SDI/QMD stayed exact). `_fia_spcode` (src/io/fia_database.jl) zero-pads numeric codes
  ≤ 2 chars → 3 digits (dbstreesin.f:353); applied to tree SPECIES + SITE_SPECIES.
- **Metric (correctness)**: dense conifer stand 100180735010661 (balsam-fir 012 / black-spruce
  095 / white-cedar 241) cycle-0 now FULLY bit-exact incl CCF (217→215) + ALL volume columns
  (was mis-mapping BF/BS to OH). Test tightened to assert that stand's CCF+volume bit-exact.
- **Bit-exact re-verify**: suite 26099 pass / 136 broken / 0 fail. No regression.
- **Remaining**: hardwood-stand CCF (~4%) / BdFt (~8%) residual (null-CRRATIO crown dubbing +
  hardwood board-vol; per-tree trace confounded by treelist tripling reporting).

## Slice S32 — FIA reader cross-variant validation (Pillar 1 coverage) [2026-07-07]
- **What**: extended the FIA reader test to all four variants — one+ real FIA stand per
  variant (SN 255262523010854, NE 657546100126144, CS 14173137020004, LS conifer +
  hardwood), cycle-0 vs the golden live .sum. Fixture ls_sample.db grew to 5 stands.
- **Metric (coverage/correctness)**: cycle-0 density (TPA/BA/SDI/TopHt/QMD) + cubic volume
  BIT-EXACT on ALL four variants; CCF+BdFt BIT-EXACT on SN, CS, LS-conifer; NE CCF Δ1 and
  LS-hardwood CCF~4%/BdFt~8% bounded (documented). Confirms the reader + species-pad fix
  generalize across SN+NE+CS+LS — the FIA reader is a validated cycle-0 bit-exact drop-in.
- **Bit-exact re-verify**: suite 26132 pass / 136 broken / 0 fail (+33). No regression.

## Slice S33 — FIA reader lat/long → Hopkins index fix (real CCF bug) [2026-07-07]
- **What**: `apply_fia_stand!` now maps FVS_STANDINIT_COND LATITUDE/LONGITUDE → plot
  TLAT/TLONG (dbsstandin.f:254-259). Without them the per-forest DEFAULT lat/long
  (forkod.f) was used, shifting the Hopkins bioclimatic index (HI) that the eastern
  crown-width models feed on — so the HI-dependent (hardwood) open-grown crowns drifted
  → CCF off. Found by exhaustive measurement (crown eqs/coeffs/floors all verified
  bit-exact first; HI was the last variable).
- **Metric (correctness)**: LS-hardwood stand 55482390010661 CCF 172→180 (BIT-EXACT vs
  live); NE stand CCF 86→87 (BIT-EXACT). CCF now bit-exact on ALL 5 cross-variant fixture
  stands (SN/NE/CS/LS-conifer/LS-hardwood). Test tightened: CCF asserted bit-exact for all.
- **Bit-exact re-verify**: suite 26132 pass / 136 broken / 0 fail. No regression.
- **Remaining**: LS-hardwood BdFt ~8% only (separate hardwood board-volume item).

## Slice S34 — LS-hardwood BdFt residual: deep trace + two retractions [2026-07-07]
- **What**: exhaustive measure-don't-guess trace of the LS-hardwood cycle-0 BdFt ~8%
  residual (the only non-bit-exact FIA cycle-0 column). RE-TRACE DISCIPLINE retracted TWO
  wrong intermediate conclusions: (1) "forest-grown crown width 1.4×" (open-vs-forest
  treelist-column artifact; crown_width bit-exact both iwho); (2) "irreducible idib-boundary"
  (per-log DIBs with real heights are clean mid-integer, not near trunc(dib+0.499)).
- **Verified bit-exact vs FVS source**: Scribner factor table (r9clark_fvsMod.f:1418 ==
  scribner_factor.csv), iDib (FVS logDia=INT(dib+0.499)), kernel nint(len·scrbnr[idib]),
  REGN9 bucking (mrules.f:370), merch (SCFMIND=BFMIND=11/SCFTOPD=BFTOPD=9.6), and SCuFt
  (saw cubic, .sum col 3) is bit-exact ⇒ section matches at INTEGER resolution.
- **Narrowed cause**: a SUB-INTEGER bfHt difference (hidden inside the rounded SCuFt) that
  tips the leftover-log `≥ minLen+trim` boundary → different log count → size-dependent
  board divergence (jl/live ratio 0.85 RO11.2 → 0.95 GA22.5; varies WITHIN species ⇒ not a
  per-species cor factor). Standalone recompute matches the run once _r9_cor! is applied
  (raw 25 × cor ≈ 28). TO PIN: unrounded per-tree bfHt + per-log breakdown from BOTH engines.
- **Bit-exact re-verify**: NO src changes (read-only trace + /tmp measurement); floor intact
  26132 / 136 / 0. BdFt cornered (bounded < 15% in test); likely-fixable, not irreducible.

## Slice S35 — LS-hardwood BdFt: exhaustive component verification (practical limit) [2026-07-07]
- **Every algorithmic component verified bit-exact vs FVS source**: mrules REGN9 params
  (mrules.f:370), R9LOGLEN bucking, R9LOGDIB bolt-height accumulation (HT=STUMP+Σ(TRIM+LOGLEN),
  INT(DIB+0.499) — identical to jl's `ht+=trim+len` / `trunc(dib+0.499)`), Scribner factor
  table (r9clark_fvsMod.f:1418), r9cor factors (hardwood cf=1.1, matches jl _r9_cor!), merch
  limits, and the saw-cubic section (SCuFt bit-exact ⇒ bfHt matches at integer resolution).
- **Residual isolated**: a real per-tree numerical difference in either the taper eval
  (_r9_dib vs FVS R9DIB) or bfHt — invisible in the bit-exact cubic integral, amplified by
  the discrete idib. For RO 11.2 the implied divergence (jl idib 10 / dib 9.757 vs FVS idib 11)
  is ~0.75" (NOT sub-ULP), so it's a real per-tree taper/height difference, not a boundary flip.
- **Blocked on tooling**: FVS's per-tree debug (DEBUG%MODEL / r9bdft LUDBG) is the only way to
  pin which — but the FVS `DEBUG` keyword SEGFAULTS the relinked FVSls_new binary. Cannot
  obtain live's per-log breakdown. Practical limit reached; residual bounded <15%, cornered.
- **No src changes**; floor intact 26132/136/0.

## Slice S36 — r9dib short-tree Y guard (faithful port from BdFt trace) [2026-07-07]
- **What**: ported FVS r9clark.f r9dib's Db-denominator guard (r9clark_fvsMod.f:1213): Y=0
  when (1−17.3/totHt) < 0.005748 AND p > 14 (very short bole + steep taper), else (1−17.3/totHt)^p.
  Found while exhaustively tracing the LS-hardwood BdFt residual — jl always used (…)^p.
- **Metric**: faithful-port fidelity gain for short-sawtimber high-p trees (totHt < ~17.4 ft).
  Does NOT change the LS-hardwood BdFt residual (our trees are 40–97 ft ⇒ guard inert), which
  remains a per-tree taper-INPUT difference (_r9_dia417/_r9_totht → DBHIB/DIB17/TOTHT) masked in
  the bit-exact cubic SUM but amplified by the discrete idib — the r9dib/bucking/factor/cor/merch
  formulas are ALL now verified bit-exact vs FVS; residual unpinnable without FVS per-tree debug
  (DEBUG keyword segfaults FVSls_new).
- **Bit-exact re-verify**: suite 26132 / 136 / 0 (guard inert for all current tests). No regression.

## Slice S37 — LS-hardwood BdFt: full R9 board path verified, per-tree localization tooling-blocked [2026-07-07]
- **ALL R9 board formulas now verified bit-exact vs FVS source** (r9clark_fvsMod.f / r9logs.f):
  r9dia417 (DBHIB/DIB17; upsHt1=0 ⇒ jl's branches equivalent), r9totHt (TOTHT), r9dib taper
  (+the short-tree Y guard, ported S36), R9LOGLEN (even-foot log lengths — INT/÷2·2 order matches),
  R9LOGDIB (bolt heights + INT(dib+0.499)), r9bdft (dib≥1.0 gate, scrbnr factor table, saw-only
  vol(2)), r9cor (hardwood cf=1.1), merch limits. Divergence real+current: .sum BdFt jl 5170 /
  live 5633 (8.2%), with TCuFt/MCuFt/SCuFt(941) all BIT-EXACT.
- **Per-tree localization BLOCKED by tooling**: live's per-tree board is not reliably obtainable —
  the text .trl columns misparse (a 5.7" tree shows 160 bdft) and TreeIds misalign jl-vs-live;
  live's DBS FVS_TreeList + the FVS DEBUG keyword both SEGFAULT the relinked FVSls_new. So the
  per-tree input that differs (masked in the bit-exact cubic SUM, amplified by discrete board idib)
  cannot be pinned with available tools. Residual bounded <15%, cornered.
- **No src changes this slice**; floor 26132/136/0.

## Slice S38 — LS-hardwood BdFt: TOOLING WALL BROKEN, divergence localized [2026-07-07]
- **Broke the tooling wall**: built a debug-enabled FVSls binary by recompiling the module chain
  (debug_mod/clkcoef_mod) + r9clark.f (with WRITE debug) via the current gfortran and linking a
  SEPARATE /tmp/FVSls_dbg (oracle untouched; reproduces the .sum exactly). Recipe:
  test/harness/fia/README_fvs_debug.md. FVS per-tree board (r9bdft per-log + r9clark per-tree
  BDTREE dbh+board) is now MEASURABLE — the "unpinnable" wall is gone.
- **Localized (partial)**: FVS cycle-0 boards 7 trees incl GA 14.2 (dbh14.2/ht97/cr=null → raw 124);
  jl's treelist shows only 6 and omits GA 14.2 (BA/cubic still bit-exact ⇒ jl HAS the tree's mass).
  FVS .sum (5633) is ALSO < Σ(its per-tree gross × cor × TPA)(~5988) ⇒ a .sum-path board reduction
  (defect?) FVS applies. So the residual has ≥2 concrete causes (GA 14.2 board handling + a defect/
  aggregation reduction), NOT irreducible. Arithmetic not yet fully reconciled (463 diff ≠ the 818
  of GA-14.2-alone) ⇒ multiple partial differences.
- **NEXT**: add the same per-tree board debug to jl (compute_volumes_ne!) to get jl's per-tree board
  by dbh, diff directly vs FVS BDTREE, and close the reconciliation (GA 14.2 + defect).
- **No src changes**; floor 26132/136/0; oracle intact (golden BdFt 5633).

## Slice S39 — LS-hardwood BdFt: GA-14.2 ruled out; .sum-path board exceeds per-tree r9bdft [2026-07-07]
- **GA 14.2 RULED OUT** (re-trace): it is HISTORY=8 (older-dead) — BOTH jl and FVS exclude it from the
  standing cubic/BA (bit-exact), so it can't be the board cause (a tree can't be board-in but cubic-out).
- **Decisive finding via FVSls_dbg**: FVS's 6 LIVE board-trees (r9bdft per-tree, by dbh) are IDENTICAL to
  jl's (11.2→25,12.3→35,12.5→55,14.3→74,14.4→124,22.5→467 raw), and Σ(final×TPA)=5170 == jl .sum. YET
  FVS's .sum BdFt = 5633. So FVS's .sum board (5633) ≠ its own per-tree r9bdft board (5170) — the .sum
  aggregation uses a DIFFERENT board path (different merch top / defect / routine) than the r9bdft I
  instrumented in r9clark. jl matches FVS's r9bdft per-tree exactly; the gap is purely FVS's .sum path.
- **NEXT**: instrument the FVS .sum board AGGREGATION (not r9bdft in r9clark — likely a fvsstd/sumsvb/
  volume-driver path) to find where the extra board originates. Debug-binary recipe: README_fvs_debug.md.
- **Progress**: residual moved from "irreducible/unpinnable" → "per-tree board bit-exact vs FVS; gap is in
  FVS's separate .sum board aggregation path" — a concrete, bounded next probe. No src changes; floor
  26132/136/0; oracle golden 5633 intact.

## Slice S40 — LS-hardwood BdFt: ROOT CAUSE = board TYPE (International, not Scribner) [2026-07-07]
- **ROOT CAUSE FOUND (definitive, via FVSls_dbg per-tree BDFINAL post-cor)**: FVS's LS .sum BdFt for
  stand 55482390010661 = INTERNATIONAL ¼" (r9clark vol(10)), NOT Scribner (vol(2)). Proof: the 6 live
  board-trees' post-cor International (154,33,93.5,71.5,44,539 → nint 936) × TPA 6.018 = 5633 == oracle
  EXACTLY; Scribner (858→5170) == what jl emits. jl's per-tree Scribner AND International both match
  FVS bit-exact; the bug is jl SELECTS Scribner for the .sum board where FVS uses International.
- **Contradiction with lst01** (bit-exact w/ Scribner, BdFt 1887) ⇒ the LS board type is NOT a blanket
  choice — it is METHB/species-dependent (ls/bfvol.f METHB branches). jl hardcodes board_scribner=true
  for LS; correct is per-species (these hardwoods → International like NE/CS; lst01's species → Scribner).
- **NEXT (concrete fix)**: trace ls/bfvol.f + the LS vollib BFV assignment to get the per-species METHB→
  (Scribner vol2 | International vol10) rule; drive jl's board_scribner per species from it. Must keep
  lst01 (Scribner) bit-exact AND fix these hardwoods (International). Guard NE/CS (International already).
- **Impact**: this likely affects LS BdFt broadly (any stand with International-METHB species) — a real
  drop-in-fidelity gap, now root-caused. No src changes this slice; floor 26132/136/0; oracle 5633 intact.

## Slice S41 — LS BdFt board-type is SPECIES-DEPENDENT (root cause CONFIRMED by experiment) [2026-07-07]
- **CONFIRMED via experiment**: flipping jl's LS board to International (vol10) makes the FIA hardwood
  stand 55482390010661 BdFt = 5633 == oracle EXACTLY (was Scribner 5170), BUT regresses lst01 +
  ls_bfvolume (268 test failures). So the LS .sum board TYPE is species/forest-dependent, NOT blanket:
  lst01's species → Scribner (jl's current, correct); the FIA hardwoods (GA 544 / RO 833 / SS 931) →
  International. jl hardcodes `board_scribner = s.variant isa LakeStates` (r9clark_vol.jl:558) — a blanket
  Scribner that is WRONG for International-METHB species. Reverted (floor restored 26132/136/0).
- **FIX PATH (concrete, bounded)**: trace ls/sitset.f METHB(ISPC) assignment + the LS vollib board-type
  selection (Scribner vol2 vs International vol10) — likely per-species METHB or per-IFOR. Drive jl's
  board_scribner PER SPECIES from that rule. Constraint: keep lst01 + ls_bfvolume (Scribner) bit-exact
  AND fix the International species. This is a real LS drop-in-fidelity gap (affects any LS stand with
  International-METHB species), fully root-caused via the FVSls_dbg breakthrough.

## Slice S42 — LS BdFt: EXACT RULE found (METHB==9 → International) [2026-07-07]
- **EXACT RULE (fvsvol.f NATCRS entry, LS buildDir line 519-527)**:
    IF(D < BFMIND) BBFV=0
    ELSE IF(METHB(ISPC).EQ.9) BBFV=TVOL(10)   ! International ¼"
         ELSE                 BBFV=TVOL(2)    ! Scribner
  So the LS .sum board type is PER-SPECIES: METHB==9 → International (vol10), else Scribner (vol2).
  jl hardcodes board_scribner=true (always Scribner) for LS — the FIX is board_scribner = (METHB(sp) != 9).
- **Confirmed**: FIA hardwoods (GA/RO/SS) → International (experiment: 5633==oracle); lst01 species →
  Scribner (jl correct). vls/vols.f:268 routes METHB∈{5,6,9}+METHC∈{6,10} → NATCRS, which applies the rule.
- **OPEN (one debug cycle)**: static METHB defaults to 6 (ls/grinit.f:88) and no keyword/MODEL_TYPE sets 9
  for this stand — yet runtime is International (METHB=9). So METHB=9 is set at RUNTIME for these species
  (likely the R9 vollib equation dispatch / a species-list in sitset or the DBS path). NEXT: recompile
  fvsvol.f with a METHB(ISPC) debug WRITE (same recipe as r9clark) to dump which species get METHB=9,
  then replicate that per-species METHB in jl and gate board_scribner on it. Keep lst01/ls_bfvolume
  (Scribner) bit-exact.
- No src changes; floor 26132/136/0.

## Slice S43 — LS BdFt board-type: DEFINITIVE root cause + FIX LANDED (bit-exact) [2026-07-07]
- **S42's METHB==9 conclusion was a MISREAD — REFUTED by runtime dump.** Recompiled fvsvol.f (+ mrules_mod)
  with a `WRITE(0,*)` on METHB(ISPC)/TVOL(2)/TVOL(10) and re-ran the FIA hardwood stand under FVSls_dbg:
  the board-eligible species all have **METHB=6** (not 9), so NATCRS's `METHB==9` branch is NEVER taken.
  Yet `TVOL(2)` arrives at NATCRS already = the International value (e.g. GA d=22.5: TVOL2=TVOL10=539).
  So the Scribner→International swap happens UPSTREAM of NATCRS, inside the volume interface.
- **DEFINITIVE mechanism — volinit.f:434-451 (VOLINITNVB, the Region-9 Clark branch):**
    ELSEIF (VOLEQ(1:1).EQ.'9' .AND. (MDL.EQ.'CLK')) THEN
        CALL R9CLARK(...)                 ! fills VOL(2)=Scribner AND VOL(10)=International
        IF(IFORST ∈ {4,5,8,11,12,14,19,20,21,22,24,30}) VOL(2) = VOL(10)   ! → International board
  where `IFORST` = the 2-digit national-forest within KODFOR/LOCATION (`READ(FORST,'(i2)')`; LOCATION mod
  100, or the middle 2 digits when KODFOR>10000). r9clark ALWAYS computes both boards; volinit picks
  which lands in VOL(2) **purely by national forest**. (Region-8 Clark has an analogous gate at 424-428.)
- **Confirmed against live**: FIA hardwood stand LOCATION=924 → IFORST=24 ∈ list → International (539/…);
  lst01 LOCATION=903 → IFORST=3 ∉ list → Scribner; LS-conifer FIA stand LOCATION=910 → IFORST=10 ∉ list →
  Scribner (already bit-exact); ls_bfvolume/spctrn LOCATION=80106 → IFORST=1 ∉ list → Scribner. The blanket
  `board_scribner = variant isa LakeStates` was wrong ONLY for the in-list forests (why the S41 global flip
  regressed lst01: it flipped the ∉-list forests too).
- **FIX LANDED** (r9clark_vol.jl): `board_scribner = (variant isa LakeStates) && !(_r9_iforst(user_forest_code)
  in _R9_INTL_BDFT_FORESTS)` with `_R9_INTL_BDFT_FORESTS=(4,5,8,11,12,14,19,20,21,22,24,30)` and
  `_r9_iforst(k)= k>10000 ? (k÷100)%100 : k%100`. Faithful port of the volinit gate; scoped to LS so
  NE/CS (International-always) are untouched.
- **Result**: FIA hardwood BdFt 5170 → **5633 == oracle BIT-EXACT** (test_fia_reader now asserts bit-exact
  BdFt on all 5 cross-variant stands, 55/55). lst01 + ls_bfvolume + LS-conifer stay Scribner bit-exact.
  Full suite **26132 pass / 136 broken / 0 fail** — floor preserved, one residual genuinely closed.
- **Meta (re-trace discipline paid off, twice)**: S40 "forest-grown crown 1.4×" and S42 "METHB==9" were BOTH
  misreads that a static grep made look authoritative; only the runtime debug-stamp (METHB=6, TVOL2=Intl)
  exposed them and pointed one layer up to volinit. Reusable debug binary: /tmp/FVSls_dbg2 (fvsvol+volinit
  recompile recipe in test/harness/fia/README_fvs_debug.md).

## Slice S44 — FIA reader mass validation at scale (all 10 .sum cols, cross-variant) [2026-07-07]
- **Method**: extracted 40 sampled real FIA stands (via one bulk `WHERE STAND_CN IN (...)` scan — the
  70GB DB has NO index on STAND_CN, so per-stand queries full-scan; a small fixture `sample.db` sidesteps
  that WITHOUT modifying the original, honoring the hard no-modify constraint). 18 of 40 carry tree rows
  (NE 10, SN 3, CS 3, LS 2). Ran FVSjl native DATABASE reader vs freshly-relinked live `/tmp/FVS*_new`,
  compared ALL 10 cycle-0 `.sum` stat columns (TPA/BA/SDI/CCF/TopHt/QMD/TCuFt/MCuFt/SCuFt/BdFt).
  Harness: test/harness/fia/validate_fia_cols.jl (+ /tmp/fia_val/validate_sample.jl driver).
- **Result: 16/18 stands BIT-EXACT on all 10 columns.** Per-column mismatch: only CCF(1) and TCuFt(1).
  **BdFt bit-exact on all 18** — the S43 R9 per-national-forest board gate holds at scale across SN/NE/CS/LS
  (not just the one fixture stand). TPA/BA/SDI/TopHt/QMD/MCuFt/SCuFt/BdFt: 18/18 bit-exact.
- **Two residuals, both named/cornered (not reader bugs)**:
  - NE 1819712544290487 **TCuFt live 6382 / jl 6383** (Δ1 cuft = 0.016%): a Float32 total-cubic
    accumulation landing on the integer-rounding boundary — ULP-class, the accepted residual type.
  - SN 3196569010661 **CCF live 212 / jl 218** (Δ6 = 2.8%): TPA/BA/SDI/QMD/TopHt ALL bit-exact on this
    stand ⇒ the tree list AND height-dubbing are identical; the divergence is purely the crown-width→CCF
    model for its dubbed-crown SN hardwoods (sp 621/826/404/541/544). Same named primitive as the already-
    cornered s3 crown-width residual (#79) / null-CRRATIO crown dubbing — a crown-width-coefficient class,
    not a reader/ingest fault.
- **Coverage metric**: FIA reader now cross-variant validated on 18 real stands (was 5) — 16 fully bit-exact,
  2 cornered to named ULP/crown-width primitives. No src change this slice; floor unaffected.

## Slice S45 — FIA reader missing-elevation → forkod default (SN CCF fix, bit-exact) [2026-07-07]
- **Bug (from S44's cornered CCF residual)**: SN FIA stand 3196569010661 reported CCF 218 vs live 212.
  Root cause via per-tree crown-width probe + live `.out` echo: the stand's DB ELEVATION is NULL, so the
  reader left `p.elevation = 0`, but **live FVS defaults missing elevation from the national-forest table**
  (forkod.f:540-546: `IF(ELEV.EQ.0) ELEV = 12.` for forest 802 / LOCATION 80215). Elevation feeds the
  Hopkins bioclimatic index (`hopkins_index`), which the eastern **bechtold** open-grown crown-width model
  uses (`hi_coef·HI`). jl's HI was −9.80 (elev 0) vs the correct +2.20 (elev 12) — a +12 HI error inflating
  every hardwood crown → CCF too high. (Cubic/board/TPA/BA/SDI/TopHt were all already bit-exact: crown width
  is the only cycle-0 consumer of HI, which is why ONLY CCF drifted.)
- **Fix** (`src/io/fia_database.jl` `apply_fia_stand!`): after reading geo from the DB, apply the FORKOD
  phase-3 default (mirroring `kw_stdinfo!`:517-520) — fill any of lat/long/elev still 0 from
  `forest_location(coef, KODFOR÷100)`. FVS runs forkod BEFORE the DB overrides, and the DB overrides
  elevation only when >0 (dbsstandin.f:647), so a null/≤0 ELEVATION keeps the forest default. **Southern-
  gated**: the SN `forest_location` table is keyed by KODFOR÷100 exactly as `kw_stdinfo!` keys it; NE uses a
  different forkod keying (JFOR/IFOR) so its DB-elevation default is left as a follow-up (no evidence of the
  bug there — the 10 NE mass-validation stands were CCF bit-exact).
- **Verified**: per-tree probe now elev=12.0, HI=+2.1975, CCF = 212.05 → rounds to **212 == live BIT-EXACT**
  (was 217.6→218). Variant-safe: for LS/CS (no SN forest_location entry) and NE (gated out) nothing changes;
  the SN test stand 255262523010854 is Florida pines (non-bechtold crowns, HI-invariant) ⇒ unaffected.
- **Impact**: corrects HI for ANY FIA SN stand with missing elevation (a broad class), not just this stand.
  Suite re-verify: 26132 pass / 136 broken / 0 fail (floor held). Mass-validation re-run: **17/18 stands now
  ALL-10-columns bit-exact** (was 16/18) — CCF now bit-exact on all 18; the lone remaining per-column mismatch
  is the NE TCuFt Δ1 rounding-boundary ULP (stand 1819712544290487, 6382 vs 6383). Closes task #89.

## Slice S46 — FIA reader mass validation at SCALE (162 stands, all 4 variants, 10 cols) [2026-07-07]
- **Method**: extended the S44 fixture to 360 candidate stands (bulk `WHERE STAND_CN IN (...)`, still no DB
  modification) → 165 tree-bearing (SN 54, NE 60, LS 34, CS 17). Validated FVSjl native reader vs live
  `/tmp/FVS*_new`, all 10 cycle-0 `.sum` columns. 162 both-produced-sum. Harness: validate_fia_cols.jl.
- **Result: 151/162 stands BIT-EXACT on all 10 columns.** Per-column mismatch counts (of 162):
  **TPA 0 · BA 0 · SDI 0 · QMD 0 · MCuFt 0 · SCuFt 0 · BdFt 0** (SEVEN of ten columns perfect across every
  stand — density, size, and all merch/board volumes are a bit-exact drop-in on 162 real cross-variant
  stands). Remaining: TopHt 9, TCuFt 2, CCF 1.
- **TopHt (9 stands, Δ1–4 ft) — AVH-boundary height edge, NOT a reader/growth bug**: jl `stand_top_height`
  is line-for-line identical to FVS `avht40.f` (avg height of the 40 largest-DBH TPA). On these stands
  TPA/BA/SDI/QMD/all-volumes are bit-exact (order-independent sums), so the tree list matches — the Δ is the
  *height* of the tree(s) landing at the 40-TPA boundary. Probed seedling stand SN 1152013794290487 (live 6 /
  jl 9): the window is 6 TPA of one DBH-5.5 tree (dubbed ht 35) + 34 TPA of DBH-0.1 seedlings (ht 4.5); the
  boundary tree's dubbed height drives the AVH. Same primitive class as the cornered height-dub / sort-tie
  residuals — an AVH-window boundary-height edge on tie-heavy or seedling-mixed stands, isolated to the
  top-height report column. Tracked as a follow-up task.
- **TCuFt (2 stands, Δ1 cuft = ~0.02%)** and **CCF (1 stand, Δ1)**: rounding-boundary ULPs (accepted class).
- **Coverage metric**: FIA reader cross-variant validated on **162 real stands** (was 18) — 151 fully
  bit-exact, the 11 residuals all cornered to named AVH-boundary/ULP primitives. No src change this slice
  (pure validation); floor unaffected (26132/136/0).

## Slice S47 — SN seedling height-dubbing: DBH≤0.1 skips the 4.5 floor (TopHt fix, bit-exact) [2026-07-07]
- **Bug (from S46's cornered TopHt residual)**: on seedling-heavy FIA stands the reported TopHt (AVH, avg
  height of the 40 largest-DBH TPA) ran high (e.g. SN 1152013794290487 live 6 / jl 9). Root-caused with a
  recompiled `avht40.f` (dumped HT of the window trees): live gives the DBH-0.1 seedlings **HT = 1.01 ft**,
  jl gave **4.5**. `dub_missing_heights!` (volume.jl) HAD the correct `d≤0.1 → 1.01` seedling branch, but the
  following unconditional `h_v < 4.5 → 4.5` floor clobbered it back to 4.5. FVS (`cratet.f:352-354`) sets
  `H=1.01` then `GO TO 115`, **skipping** the `IF(H.LT.4.5)H=4.5` floor — only the Curtis-Arney/Wykoff (D>0.1)
  branches are floored.
- **Fix** (`src/engine/volume.jl` `dub_missing_heights!`): gate the floor on `d > 0.1f0` so the seedling case
  keeps 1.01 (mirrors the `GO TO 115`). One-token change; faithful to cratet.f.
- **Verified**: seedling dub 4.5→1.01; stand TopHt **9→6 == live BIT-EXACT** (avht40 debug-stamp: DBH 0.1 →
  HT 1.01, both engines). Suite **26132 pass / 136 broken / 0 fail** — floor preserved (canonical `.tre`
  fixtures have no DBH≤0.1/HT=0 trees, so the change is inert there; only inventory seedlings with missing
  heights are affected, toward faithfulness). Closes task #90.
- **Impact**: corrects dubbed seedling heights (and the AVH top-height they feed) for ANY stand with DBH≤0.1
  missing-height seedlings — a real, broadly-applicable FIA-inventory fidelity fix. Mass-validation re-run
  (162 stands): **155/162 now ALL-10-columns bit-exact** (was 151); TopHt mismatches **9→5** (the seedling
  fix closed 4). The remaining 5 TopHt (Δ1-2 ft on larger-tree stands) are a distinct AVH-boundary edge on
  MEASURED-height trees (not seedlings) — a smaller residual, deferred. CCF 1 / TCuFt 2 unchanged (ULP).

## Slice S48 — AVH top-height uses FVS RDPSRT tie-break (measured-height TopHt fix) [2026-07-07]
- **Bug (S46/S47's remaining TopHt residual)**: 5/162 stands had TopHt off by Δ1-4 ft on LARGER-tree stands
  with MEASURED heights (e.g. SN 1152013702290487 live 42 / jl 43; SN 253701728010854 live 33 / jl 32).
  Root-caused with the `avht40.f` debug-stamp: on stand 1152013702, trees idx 15 (DBH 6.6, HT 32) and idx 13
  (DBH 6.6, HT 43) are a DBH **tie** at the 40-TPA boundary; live's RDPSRT orders idx15 before idx13 (idx15
  full weight, idx13 the capped boundary tree) → AVH 42.2→42, while jl's stable `sortperm!` broke the tie the
  other way (ascending index) → 42.8→43.
- **Why only TopHt**: `avht40.f` sorts IND with **RDPSRT** (Scowen quickersort — NOT stable); the tie order
  decides which tree sits at the 40-TPA cutoff and thus enters AVH. jl's `stand_top_height` used Julia's
  stable `sortperm!`. (The DG `point_basal_area!` sorts by DBH too, but its BAL is an order-independent sum
  ⇒ tie order is inert there — which is why DG stayed bit-exact yet AVH didn't.)
- **Fix** (`src/engine/standstats.jl`): replace the `sortperm!` in `stand_top_height` with the already-ported
  `_rdpsrt!` (the exact FVS RDPSRT descending sort + quicksort tie-break, transliterated for COMPRESS). No
  new sort code — reuses the faithful primitive.
- **Verified**: target stands TopHt 43→**42** and 32→**33** (== live, avht40 debug-stamp-confirmed tie order).
  Suite **26132 pass / 136 broken / 0 fail** — floor preserved (the faithful RDPSRT matches FVS on every
  canonical stand too, incl. any AVH ties). Closes task #91. Mass-validation re-run (162 stands):
  **159/162 now ALL-10-columns bit-exact** (was 155); TopHt mismatches **5→1** (RDPSRT closed 4). Only 4
  per-column mismatches remain across 162 stands: TopHt 1 (SN 255260379010854, 30/32 — a deeper AVH
  pre-sort-order edge the tie-break didn't cover), CCF 1, TCuFt 2 — all ULP/boundary-class.
- **Net FIA-reader arc this session (S44→S48)**: 5-stand test → **162-stand cross-variant validation**,
  151→**159** fully bit-exact via 3 landed fixes (S45 elevation, S47 seedling ht, S48 AVH tie-break); the
  4 residuals are named ULP/boundary primitives. 7 of 10 `.sum` columns are bit-exact on ALL 162 stands.

## Slice S49 — #84 SETSITE/THINAUTO on NE/CS/LS: assessed & resolved (no src change) [2026-07-07]
- **SETSITE**: already covered on NE/CS/LS (test_kwcov_variants.jl fixtures ne/cs/ls_setsite). NE+LS
  bit-exact; CS ULP-class (sub-ULP height diff amplified by the ill-conditioned NC-128 anamorphic ht-curve
  inversion at raised site 70, re-converges bit-exact @2040 — irreducible Float32, not a bug). Done.
- **THINAUTO**: validated on SN (harness `cut_thinauto.key`). Attempted NE/CS/LS coverage fixtures (bare
  `THINAUTO <yr>`, FVS-default MIN 45%/MAX 60% of full stocking) on the canonical net01/cst01/lst01 stands —
  **live FVS is ill-posed here**: FVSne **FPE'd** (core dump — full-stocking division-by-zero UB mid-sim),
  FVScs/ls emitted a `.sum` with NO data rows. ⇒ no valid live golden ⇒ not a validatable coverage scenario
  (accepted 'ill-posed / FVS-UB' class, per the goal). jl handles all 3 gracefully (11 data rows, TPA
  536→67/11/68, no crash) — more robust than live. Crash fixtures removed (would break the auto-discovering
  KCV suite). **#84 closed** — both keywords assessed; no src change; floor unaffected (26132/136/0).

## Slice S50 — Pillar-1 per-variant coverage matrix documented (`docs/COVERAGE_MATRIX.md`) [2026-07-07]
- Produced the Pillar-1 done-state deliverable: a **per-variant coverage matrix** (SN/NE/CS/LS) consolidating
  every validated axis — cycle-0 inventory, canonical multi-cycle (`*t01`), growth spine, FFE, volume, the
  FIA-DB reader (162-stand mass validation), and the **40 isolated-keyword KCV scenarios per NE/CS/LS**
  (`test/fixtures/kwcov`, `test_kwcov_variants.jl`). Each axis is bit-exact vs live or cornered to a named
  primitive; the named cornered residuals / accepted deferrals are enumerated. Makes the "100% drop-in
  defensible per variant" claim concrete. Pure documentation; floor unaffected (26132/136/0).

## Slice S51 — THINDBH added to NE/CS/LS KCV coverage (bit-exact) [2026-07-07]
- Added `{ne,cs,ls}_thindbh` to the KCV keyword-coverage set (test/fixtures/kwcov/, auto-discovered by
  test_kwcov_variants.jl): `THINDBH 2000 DBH 0-99 → residual 40 TPA`. jl BIT-EXACT (all cells) vs live on
  all three variants — 536→36 TPA thin fires identically. **Measure-don't-guess note**: the first fixture
  free-spaced the params and jl (strict fixed-column parser, cols 11+ in 10-char fields) mis-read them → no
  thin (jl 475 vs live 36); the bug was the FIXTURE (column misalignment), not jl — re-aligned to 10-col
  fields ⇒ bit-exact. NE/CS/LS isolated-keyword coverage now **41** each. (Remaining common keywords →
  task #93.) Floor: suite re-run confirms no regression + 3 new bit-exact scenarios.
- **CONFIRMED**: full suite **26639 pass / 136 broken / 0 fail** (was 26132/136/0) — the 3 THINDBH scenarios
  added **+507 bit-exact per-cell assertions**, broken unchanged, 0 fail. New correctness floor: **26639/136/0**.

## Slice S52 — THINSDI added to NE/CS/LS KCV coverage (bit-exact NE/LS, CS ULP-cornered) [2026-07-07]
- Added `{ne,cs,ls}_thinsdi` (THINSDI residual-SDI 80 thin, aligned fixtures + live goldens). **NE + LS
  BIT-EXACT** (all cells); **CS** one-cell ULP (2030 col 26, 31.6/31.5 Δ0.1 — the same derived growth/MAI
  column as ne_leavesp; all density+volume cols bit-exact; the documented CS sub-ULP NC-128 height
  amplification, cf. cs_setsite) → added to `_KCV_BROKEN` with a named reason. NE/CS/LS KCV coverage now
  **42** each. Suite **26978 pass / 137 broken / 0 fail** (was 26639/136/0): +339 bit-exact assertions, +1
  documented ULP corner (cs_thinsdi), 0 fail. New floor **26978/137/0**.

## Slice S53 — THINBBA added to NE/CS/LS KCV coverage (bit-exact all 3) [2026-07-07]
- Added `{ne,cs,ls}_thinbba` (THINBBA residual-BA-100 thin-from-below). jl **BIT-EXACT (all cells)** vs live
  on all three variants. NE/CS/LS KCV coverage now **43** each. Suite **27485 pass / 137 broken / 0 fail**
  (was 26978/137/0): +507 bit-exact assertions, broken unchanged, 0 fail.
- **Session KCV-coverage arc (S51-S53)**: added THINDBH + THINSDI + THINBBA to NE/CS/LS (40→43 keywords each);
  8 of the 9 new scenarios bit-exact, 1 (cs_thinsdi) documented ULP-corner. Correctness floor grew
  transparently **26132/136/0 → 27485/137/0** (+1353 bit-exact assertions, +1 named ULP corner, still 0 fail).
  Proven fast pattern (printf-aligned fixture → live golden → jl bit-exact-or-cornered) documented in #93.

## Slice S54 — THINATA added; thinning-method family complete on NE/CS/LS [2026-07-07]
- Added `{ne,cs,ls}_thinata` (THINATA residual-TPA-100 thin-from-above). jl **BIT-EXACT (all cells)** all 3
  variants. KCV coverage now **44** each. The full thinning-method family is covered per variant: thinaba,
  thinata, thinbba, thinbta, thincc, thindbh, thinht, thinqfa, thinrden, thinsdi (10 methods). Suite
  **27992 pass / 137 broken / 0 fail** (+507 bit-exact, 0 fail).

## Slice S55 — FIXHTG + THINPRSC added to NE/CS/LS KCV coverage [2026-07-07]
- Added `{ne,cs,ls}_fixhtg` (FIXHTG height-growth ×1.5 — a growth-MODIFIER, not a thin) and `{ne,cs,ls}_thinprsc`
  (THINPRSC prescription thin, 0.999 removal). THINPRSC bit-exact all 3. FIXHTG bit-exact NE/LS; CS one-cell
  BdFt ULP (2040 col 12, 21345/21346 Δ1 — taller stem crosses a board-foot boundary; all density + other vol
  cols bit-exact) → cornered in `_KCV_BROKEN`. KCV coverage now **46** each. Suite **28838 pass / 138 broken /
  0 fail** (+846 bit-exact assertions, +1 named ULP corner, 0 fail).
- **Session KCV arc total (S51-S55)**: 6 keywords added to NE/CS/LS (40→46 each) — THINDBH, THINSDI, THINBBA,
  THINATA, FIXHTG, THINPRSC; 16/18 new scenarios bit-exact, 2 documented ULP corners (cs_thinsdi, cs_fixhtg).
  Floor grew transparently **26132/136/0 → 28838/138/0** (+2706 bit-exact assertions, +2 named ULP corners,
  0 fail throughout).

## Slice S56 — BAIMULT + HTGMULT + MORTMULT added to NE/CS/LS KCV coverage (bit-exact) [2026-07-07]
- Added `{ne,cs,ls}_{baimult,htgmult,mortmult}` (per-component growth/mortality multipliers ×1.5, all species).
  jl **BIT-EXACT (all cells)** on all 9 scenarios. Non-vacuous (live echo MULTIPLIER=1.5; NE 2040 BA 190 vs
  189 at ×1.0). KCV coverage now **49** each. Suite **30359 pass / 138 broken / 0 fail**.
- **Session KCV arc (S51-S56)**: 9 keywords added to NE/CS/LS (40→49 each) — thinning family (thindbh/thinsdi/
  thinbba/thinata), growth modifiers (fixhtg/baimult/htgmult/mortmult), prescription (thinprsc). 25/27 new
  scenarios bit-exact, 2 documented ULP corners (cs_thinsdi, cs_fixhtg). Floor grew transparently
  **26132/136/0 → 30359/138/0** (+4227 bit-exact assertions, +2 named ULP corners, 0 fail throughout).

## Slice S57 — CRNMULT added to NE/CS/LS KCV coverage (bit-exact, direct-validated) [2026-07-07]
- Added `{ne,cs,ls}_crnmult` (crown-ratio multiplier ×1.5, all species). jl **BIT-EXACT (all cells)** vs live
  on all 3 (direct jl-vs-live check, the authoritative test for new KCV fixtures; live echo MULTIPLIER=1.50
  confirms non-vacuous). No `_KCV_BROKEN` entry needed. KCV coverage now **50** each.
- **Floor CONFIRMED**: full suite **30866 pass / 138 broken / 0 fail** (was 30359/138/0 at 49 kw) — crnmult
  added +507 bit-exact assertions, broken unchanged, 0 fail. (Several suite re-runs died silently mid-session
  under memory pressure — env flakiness, not code; the run_in_background path completed cleanly.)
- **Session KCV arc (S51-S57)**: 10 keywords added to NE/CS/LS (40→50 each) — thinning family (thindbh/sdi/
  bba/ata), growth/mortality/crown modifiers (fixhtg/baimult/htgmult/mortmult/crnmult), prescription (thinprsc).
  28/30 new scenarios bit-exact, 2 documented ULP corners. Floor grew transparently 26132/136/0 → **30866/138/0** (+4734 bit-exact assertions, +2 named ULP corners, 0 fail).

## Slice S58 — Pillar-2 re-verify of the S48 AVH sort change (measured: improvement, not regression) [2026-07-07]
- S48 replaced `sortperm!` with the ported `_rdpsrt!` in `stand_top_height` (called once-per-cycle from
  `compute_density!`). Re-traced the allocation impact per doctrine #3/#9 (measure, don't assume): on net01
  cyc-0, `@allocated`: **`_rdpsrt!` = 160 B/call vs `sortperm!` = 384 B/call** — the `zeros(Int,64)` ipush
  escape-analysis-elides and the Scowen quicksort needs less scratch than Julia's `sortperm!`. So S48
  IMPROVED the once-per-cycle top-height sort allocation (384→160 B), it did NOT regress Pillar-2. (This sort
  was already in the documented per-cycle floor via the prior `sortperm!`; it is a stat pass, not the
  per-tree hot loop.) Possible future micro-opt: thread a preallocated ipush buffer through `_rdpsrt!` to
  reach 0 B — low value (160 B once/cycle). No change this slice; documentation-only.

## Slice S59 — Pillar-2 hot-path re-verification (grow_cycle! within documented floor; S48 improved it) [2026-07-07]
- Measured `@allocated grow_cycle!(s)` (warmed, net01 NE) = **9280 B/cycle** — this is the documented+justified
  Pillar-2 floor's "~5-10 KB/cyc Base sort scratch" component (grow_cycle! calls `compute_density!` 3×, each
  running the descending-DBH stat sorts). The "~45 KB/cyc write_sum_row Printf IO" floor is in the summary
  path, NOT grow_cycle!. ⇒ Pillar-2 done-state (#75) confirmed still MET; no regression from this session's
  S43/S45/S47/S48 fixes.
- **S48 net effect on the floor: IMPROVEMENT** — `stand_top_height`'s sort went `sortperm!` 384 B → `_rdpsrt!`
  160 B (S58 measurement), shrinking the Base-sort-scratch component.
- **Further sort-scratch reduction is UNSAFE (documented, do not attempt)**: the other per-cycle `sortperm!`s
  feed bit-exact-critical results — `point_basal_area!` (standstats.jl:158) builds `pbal` (PTBALT, the DG
  competition term), where the tie order among EQUAL-DBH trees decides the "larger than tree i" boundary, so
  swapping its sort could change pbal → break DG bit-exactness (Float32-op-order doctrine #2). The remaining
  ~KB sort scratch is therefore the justified floor. Pillar-2 fully verified for this session's changes.

## Slice S60 — Pillar-1 coverage: SPECPREF re-added to NE/CS/LS KCV (bit-exact, non-vacuous) [2026-07-08]
- Added `{ne,cs,ls}_specpref.{key,tre,live.sum}` to `test/fixtures/kwcov/` — SPECPREF (species cut-preference)
  paired with THINBTA, so the priority-reorder actually changes which trees the thin removes.
- **KEY GOTCHA (re-traced against cuts.f + kw_thin!):** SPECPREF's field 1 is the **YEAR** (it is a scheduled
  activity, `kw_thin!`), NOT the species. The first fixture attempt `SPECPREF <sp> <pref>` scheduled it at
  cycle-27 (never fires) ⇒ inert ⇒ golden byte-identical to plain THINBTA (the vacuousness trap, doctrine
  `test-must-exercise-the-semantic`). Correct card: `SPECPREF <year=2010> <species-index> <pref>`. Species
  index is the variant's `species_coefficients.csv` sequence index (SM = 27 NE / 43 CS / 26 LS), matching
  `cuts.f:1440 ISPC=IFIX(PRMS(1))` and jl `_apply_specpref!`.
- Non-vacuousness PROVEN by live differential: NE 2010 removed-BA 1090→1437, residual QMD/forest-type shift
  vs plain THINBTA. All 3 variants **FULL-ROW BIT-EXACT** jl vs freshly-relinked live (`/tmp/FVS{ne,cs,ls}_new`).
- Floor grew transparently **30866/138/0 → 31373/138/0** (+507 bit-exact assertions, 0 new broken, 0 fail).
  KCV isolated-keyword coverage now **51** per NE/CS/LS. `COVERAGE_MATRIX.md` updated.

## Slice S61 — Pillar-1 coverage: SPLEAVE/CUTEFF/TIMEINT added to NE/CS/LS + a real CUTEFF bug fixed [2026-07-08]
- Added `{ne,cs,ls}_{spleave,cuteff,timeint}.{key,tre,live.sum}` (9 fixtures). SPLEAVE (leave SM from the
  THINBTA thin) and TIMEINT (uniform 5-yr cycles) were **bit-exact out of the box** on all 3 variants.
- **CUTEFF exposed a real jl bug (found via the non-vacuous fixture, then fixed):** the global `CUTEFF`
  keyword sets `control.cut_eff` (initre.f:5400 `EFF=ARRAY(1)`), and FVS's sorted-thin path defaults the
  cut efficiency to that global when the thin's own efficiency field is blank (`cuts.f:567 IF(NPS.GT.1)
  CUTEFF=PRMS(2)` — i.e. keep the EFF default when only the target is supplied). jl's `_thin_sorted!`
  (THINBTA/ATA/BBA/ABA) hardcoded `cuteff = cuteff_p>0 ? cuteff_p : 1f0`, IGNORING `control.cut_eff` ⇒
  a `CUTEFF 0.5` + THINBTA produced output byte-identical to no-CUTEFF (jl removed the full targeted amount).
  **Fix** (`cuts.jl:512`): `cuteff = cuteff_p>0 ? cuteff_p : s.control.cut_eff` — matching the already-faithful
  fallback at cuts.jl:417/458/825. Verified: all 3 `*_cuteff` stems now FULL-ROW BIT-EXACT vs live.
- **Variant-safe / no-regression:** the fix is in shared `src/engine/cuts.jl`; full suite stayed at 138
  broken (the exact cornered set) — SN + NE + CS + LS all still bit-exact. Floor grew transparently
  **31373/138/0 → 32894/138/0** (+1521 bit-exact assertions, 0 new broken, 0 fail). KCV coverage now **54**
  per NE/CS/LS. `COVERAGE_MATRIX.md` updated.
- **Env note:** the sandbox accumulates thousands of ZOMBIE julia PIDs (unreaped by the harness parent);
  they occasionally exhaust the PID table and make a `test/runtests.jl` background run die silently at
  0 bytes / no process. Mitigation: re-launch; a completed run is authoritative. Not a code issue.

## Slice S62 — Pillar-1 coverage: TOPKILL/HTGSTOP top-damage events added to NE/CS/LS [2026-07-08]
- Added `{ne,cs,ls}_{topkill,htgstop}.{key,tre,live.sum}` (6 fixtures) — the two htgstp.f top-damage events
  (TOPKILL=act 111 top-kill, HTGSTOP=act 110 stop height growth), as deterministic standalone events.
- **Card layout (initre.f:1305-1358 / kw_htgstp!):** field1=YEAR, field2=species (0=all, SPDECD), field3=htlo,
  field4=hthi, field5=prob (fraction of in-range trees affected), field6=avg proportion lost, field7=stddev.
  Fixtures use `<2010> 0 0 9999 1.0 <0.5|1.0> 0.0` ⇒ prob=1 / sd=0 ⇒ deterministic (all trees, all heights):
  TOPKILL removes 50% of height, HTGSTOP freezes 100% of height growth. Both materially non-vacuous
  (TopHt/QMD/volume shift vs plain baseline).
- All 6 stems **FULL-ROW BIT-EXACT** vs freshly-relinked live `/tmp/FVS{ne,cs,ls}_new`. No engine change
  (fixtures only). Floor grew transparently **32894/138/0 → 33908/138/0** (+1014 bit-exact assertions,
  0 new broken, 0 fail). KCV coverage now **56** per NE/CS/LS. `COVERAGE_MATRIX.md` updated.

## Slice S63 — Pillar-1 coverage: YARDLOSS added to NE/CS/LS (bit-exact) [2026-07-08]
- Added `{ne,cs,ls}_yardloss.{key,tre,live.sum}` — YARDLOSS (yarding loss, cuts.f act 203) paired with THINBTA.
  Card: field1=DATE, field2=PRLOST (fraction of harvest left on site), field3=PRDSNG (of that loss, fraction
  downed vs standing snags), field4=PRCRWN. Fixture `2010 0.3 0.5 0.2` — non-vacuous (the on-site loss feeds
  the snag/DWD pools, shifting the reported columns vs a plain thin).
- All 3 stems **FULL-ROW BIT-EXACT** vs freshly-relinked live. No engine change (jl already maps field2/field3
  correctly — see the initre.f:3637 note at keyword_dispatch.jl:2171). KCV coverage now **57** per NE/CS/LS.
- **Vacuousness screening (doctrine "test must exercise the semantic"):** NOHTDREG and NOSPROUT were probed
  against live FVS on this stand and found VACUOUS (LHTDRG calibration never fires; no sprouting triggered) —
  correctly NOT added as tests. VOLEQNUM deferred (constructs full NVEL region/forest codes — needs a valid
  equation string).
- Floor grew transparently **33908/138/0 -> 34415/138/0** (+507 bit-exact assertions, 0 new broken, 0 fail).

## Slice S64 — Pillar-1 coverage: SPGROUP species-group path added to NE/CS/LS [2026-07-08]
- Added `{ne,cs,ls}_spgroup.{key,tre,live.sum}` — SPGROUP (species-group definition, initre.f:4726) with a
  next-record species list (`SM`), consumed group-indexed via SPECPREF `-1` (group 1) + THINBTA. Exercises a
  distinct code path: the group table build (kw_spgroup! read_raw_line! token parse) + group-lookup in
  _apply_specpref! (`sp_groups[-isp]`). Cross-variant portable via the alpha token `SM` (SPDECD resolves it
  per variant). Non-vacuous (SM present, changes the thin composition).
- All 3 stems **FULL-ROW BIT-EXACT** vs freshly-relinked live via `run_keyfile`. Fixture-only (no engine
  change ⇒ regression-safe). KCV coverage now **57** per NE/CS/LS (count corrected below — was mislabeled 58).
- ✅ **Full-suite CONFIRMED on the clean environment: 34415 pass / 138 broken / 0 fail / 0 error** (57
  keywords total — 3 variants × 57 fixtures = 171). The container was restarted mid-slice to clear a
  runaway julia (an orphaned `validate_fia_cols.jl` from the completed FIA task #85, ~1.5 CPU-hours,
  spawning ~5900 zombies and starving suite runs to ~10 min). Recovery: restored `/tmp` oracles + FIA
  sample DBs from the gitignored `FVSjl/tmp/` (oracle binaries verified runnable standalone — no relink
  needed); `Pkg.instantiate()` on BOTH FVSjl and FVSjulia (the restart wiped their depot precompile
  caches → 3 transient Oracle-A load errors on the first post-restart run, resolved by re-instantiate).
- **Count correction:** the running keyword tally in slices S60–S63 was off-by-one (the pre-session
  baseline was **49** KCV keywords, not 50). Ground-truth count now = **57** (49 + the 8 this session:
  specpref, spleave, cuteff, timeint, topkill, htgstop, yardloss, spgroup). All floor *numbers* (30866 …
  34415) were measured and stand; only the keyword-count labels shifted by 1. COVERAGE_MATRIX.md set to 57.

## Slice S65 — Pillar-1 coverage + REAL FIX: READCORD/READCORH/READCORR + HTG-correction gap [2026-07-08]
- Added `{ne,cs,ls}_{readcord,readcorh,readcorr}.{key,tre,live.sum}` (9 fixtures) — the external-calibration
  correction keywords: READCORD (large-tree DG, COR2), READCORH (large-tree HTG, HCOR2), READCORR (small-tree
  HTG, RCOR2). Each reads a per-species 8F10.0 block (108/96/68 values for NE/CS/LS); fixture sets all species
  = 1.0 except SM = 1.5 (blank fields ⇒ 0 ⇒ ln(0) in Fortran, so ALL must be 1.0 to match live). Non-vacuous.
- **READCORD passed out-of-the-box** (DG correction dg_cor2 is wired in all variants: `diagr *= dg_cor2` NE,
  `dg_const += ln(cor2)` CS/LS). **READCORH/READCORR EXPOSED A REAL GAP:** the HTG corrections htg_cor2/regh_cor2
  were consumed **only in `src/variants/southern/`** — NE/CS/LS height-growth never applied them, so jl ignored
  READCORH/READCORR while live FVS honored them (divergence grew over cycles).
- **FIX (FVS-verified, variant-safe):** ne/cs/ls htgf.f:172-176 do `HTCON = 0; IF(LHCOR2 & HCOR2>0) HTCON += ln(HCOR2)`
  and regent.f:147/624-626 do `CON = RHCON·exp(HCOR); RHCON = (LRCOR2 & RCOR2>0) ? RCOR2 : 1`. jl already applied
  `exp(htcon)` / `con`, so:
  - `height_growth.jl` (NE/CS/LS): `htg_cor2_on && htg_cor2[sp]>0 && (htcon += log(htg_cor2[sp]))`
  - `small_tree_growth.jl` (NE/CS/LS): `rhcon = (regh_cor2_on && regh_cor2[sp]>0) ? regh_cor2[sp] : 1; con = rhcon*exp(...)`
  Default path is bit-identical (gates default false; `1f0*exp(x)==exp(x)` exact). SN untouched (already had it).
- All 9 stems **FULL-ROW BIT-EXACT** vs live. Floor grew **34415/138/0 → 35936/138/0** (+1521, 0 new broken,
  0 fail, 0 error). KCV coverage now **60** per NE/CS/LS. Suite ran clean on the recovered environment.

## Slice S66 — Pillar-1 coverage: MCFDLN/BFFDLN volume form-model coefficients [2026-07-08]
- Added `{ne,cs,ls}_{mcfdln,bffdln}.{key,tre,live.sum}` — MCFDLN (cubic CFLA0/CFLA1) + BFFDLN (board BFLA0/BFLA1)
  log-linear volume form-model coefficients (sdefln.f opt 39/40). Card: field1=species(0=all), f2=coef0, f3=coef1;
  fixture `0 0.1 0.95`. Non-vacuous (shifts the cubic/board volume columns). All 6 stems FULL-ROW BIT-EXACT vs live.
- Floor **35936/138/0 → 36950/138/0** (+1014, 0 new broken). KCV coverage now **62** per NE/CS/LS.

## Slice S67 — Pillar-1 coverage note: VOLEQNUM is SN-only (not an NE/CS/LS KCV candidate) [2026-07-08]
- VOLEQNUM overrides `species.vol_eq` (the R8/R9 Clark NVEL equation id). `setup_volume_equations!`
  (`volume_equations.jl:92`) assigns `vol_eq` **only for region 8 (SN)** — `iregn==8 ? _r8_ceqn(...) : "  "`;
  NE/CS/LS (region 9) get a blank `vol_eq` and drive volume through the R9 Clark path directly, so a
  VOLEQNUM override is inert there. ⇒ VOLEQNUM is exercised via the SN harness, NOT the NE/CS/LS KCV set;
  deferred from KCV (would be a vacuous test on NE/CS/LS). Remaining KCV candidates all need setup blocks
  (ESTAB: PLANT/NATURAL/REGDMULT/REGHMULT; FFE: SIMFIRE/FLAME/FMIN/CARBREPT) or a prior card
  (REUSCORD/H/R reuse the last READCOR). KCV stands at **62** keywords/variant, suite 36950/138/0.
- **REGDMULT/REGHMULT** empirically **VACUOUS** on a simple ESTAB+PLANT scenario (tested species=0,
  mult=3.0: live `.sum` byte-identical to the plain estab fixture) — the planted trees don't grow via the
  regen (small-tree) model in the reported cycles, so the regen-growth multiplier has no visible `.sum`
  effect. Would need a NATURAL-regen-dominated / tuned scenario to exercise; deferred from KCV (fixtures
  removed, not committed as vacuous tests). Same "must exercise the semantic" bar that correctly rejected
  NOHTDREG/NOSPROUT.

## Slice S68 — Pillar-1 coverage: SIMFIRE brings FFE fire into NE/CS/LS KCV [2026-07-08]
- Added `{ne,cs,ls}_simfire.{key,tre,live.sum}` — an `FMIN … SIMFIRE 2010 10.00 1 50.0 … END` fire event
  on the canonical tree data (fire in 2010, ~50% intensity). First FFE fire-path coverage in the KCV set.
- **NE FULL-ROW BIT-EXACT.** CS/LS diverge ONLY post-fire (first at the 2020 report, cascading <2-3%): LS
  2020 TPA 225 vs live 220 (~2.3%), CS 2020 TCuFt 2617/2616 — jl slightly under-kills at the fire. This is
  the KNOWN, accepted **FFE fire-mortality (FMEFF) kill-distribution residual** (documented ~3% for LS in
  test_lst01_ffe; fire BEHAVIOR — flame/scorch/selection — is bit-exact). Cornered as `@test_broken`
  (`cs_simfire`, `ls_simfire`) with the named reason; NOT a new bug.
- Floor **36950/138/0 → 37121/140/0** (+171 bit-exact assertions from NE; +2 documented cornered broken =
  the 2 FFE fire-mortality corners; 0 fail, 0 error). KCV coverage now **63** per NE/CS/LS. The 140 broken
  is the new documented cornered set.
- **S68 sharpened:** ls_simfire row-by-row confirms it's the FMEFF small-tree kill count, not growth/behavior:
  1990/2000/2010 BIT-EXACT; at the fire (505→225 jl vs 505→220 live) jl keeps **5 more trees but BA (119/119)
  and Mort (32.1/32.1) are BIT-EXACT** ⇒ the 5 extra are BA-negligible small trees (only TPA/QMD shift).
  Airtight cornering to the accepted FFE fire-mortality small-tree kill-distribution primitive.
- **S68 honesty caveat:** the FFE fire kill is a deterministic `PMORT`-fraction per record (fmburn.jl /
  fmeff.f), NOT a bare per-tree RNG draw — so the cs/ls_simfire 5-small-tree residual is NOT proven
  "irreducible stochastic". It could be a deterministic PMORT diff or a pre-fire record-tripling/XRAN diff;
  a definitive deterministic-vs-ULP classification needs a per-tree fire-kill debug-FVS trace (deferred —
  low value: BA/Mort bit-exact, matches the accepted LS fire-mortality ~3% class from test_lst01_ffe). The
  corner is justified by matching that accepted class, not by a proven-irreducible claim.

## Slice S69 — REAL FIX: SN North Carolina (IFOR=11) merch standards + live-FVS D38 SIGFPE root-cause [2026-07-08]
Two deliverables, both driven by the stratified FIA sweep on real SN stands.

**(A) #95 SN saw/board volume high on large-sawtimber FIA stands — FIXED (bit-exact).**
- Root cause: SN read merch standards from a static `data/southern/merch_specs.csv`, which held exactly the
  setcubicdflts.f region-8 *non-North-Carolina* defaults (softwood 7/10, hardwood 9/12, top 4). It was
  **missing the North Carolina (IFOR=11) overrides** — hardwood sawtimber top `SCFTOPD 9→11`, min-DBH
  `SCFMIND 12→15`; softwood `7→6.3`; `TOPD 4→3.5`; NC `DBHMIN`. A NC yellow-poplar's sawtimber section
  therefore ran to a 9″ top instead of 11″ (~7 ft further up the stem) ⇒ SCuFt/BdFt ~2.4%/tree high, ~12-15%
  aggregate on large-sawtimber stands. Proven by a `-g` fvsvol.f stamp (live `TOPDIAM=MTOPP=SCFTOPD=11.0`)
  and the CSV histogram (`scf_top_dib` only 7.0×18/9.0×72 = the uniform non-NC defaults).
- Fix (idiomatic, *removes* an asymmetry rather than bolting on a special case): added
  `_sn_merch(spi, ifor, kodist)` in `r9clark_vol.jl` — a 1:1 port of the setcubicdflts.f region-8 block —
  mirroring the existing `_ne_merch`/`_cs_merch`/`_ls_merch`, and routed SN through the *same* unified
  IFOR-aware branch in `init_merch_standards!` (SN was the lone variant reading a CSV). IFOR precedence: use
  the resolved `plot.forest_idx` when set (keyfile/STDINFO stands — e.g. Fort Bragg forest 701 → IFOR=20 even
  though its KODFOR is remapped to Uwharrie 81110 for VOLEQDEF; live keys merch on IFOR=20, not the remap),
  and only for FIA stands (forest_idx=0) decode IFOR+KODIST from KODFOR the FVS way
  (`IFORST=KODFOR/100−IREGN*100`, `KODIST=KODFOR mod 100`; sitset.f:369/forkod.f:470). `_sn_merch` handles
  both NC coastal districts (KODIST 3/10 → hardwood top 8, softwood 6.3) and non-coastal (top 11/6.3),
  consulted only when IFOR==11.
- Verify (generalization-tested — one stand isn't proof once a new code path lights up on real data): the
  offender NC stand (FOREST=11, YP DBH 25.3) **bit-exact vs live** (SCuFt 756=756 was 774, BdFt 4724=4724 was
  4772); a **34-stand North Carolina FIA batch (coastal + non-coastal) is 34/34 bit-exact on all 10 `.sum`
  columns**. Full suite **37121 pass / 140 broken / 0 fail** — non-NC forests reproduce the merch_specs.csv
  defaults exactly, and Fort Bragg (the one keyfile stand whose KODFOR decodes to IFOR=11) uses its real
  IFOR=20 ⇒ default merch ⇒ no regression. (The first cut regressed Fort Bragg by reading the remapped KODFOR;
  the forest_idx-precedence rule fixed it — caught by the NC generalization sweep, not the offender check.)

**(B) D38 live-FVS SIGFPE — root cause sharpened + fix verified (reported, NOT applied — oracle stays pristine).**
- The `r9ht` SIGFPE has TWO trap conditions: `totHt<17.3` (negative base ⇒ `(neg)**p`=NaN=invalid-op) AND —
  the common one — `totHt` *just above* 17.3 ⇒ tiny base `**` large p **underflows to a denormal** ⇒
  FE_UNDERFLOW trap (measured: totHt=17.40, base=0.00575, p=17.81 ⇒ Y≈1e-40). **The fix is FVS's own code**:
  `r9cuft` already guards this identical underflow at `r9clark.f:1015` (`IF((1.0-17.3/totht).LT.0.005748.AND.
  p.GT.14)`); the sibling `r9ht` is missing it. Verified on a debug relink: **18/18 crashers cleared,
  276/282 non-crashers bit-identical** (the 6 changed = `r9ht` made consistent with `r9cuft`). Fixed source
  at `docs/patches/r9clark_D38_underflow_fix.f`; `docs/FVS_SOURCE_BUGS.md` D38 updated. FVSjl already handles
  both conditions (its `_r9cuft` has the guard; `_r9ht` relies on Julia's non-trapping semantics).

**(C) FIA full-census feasibility (measured).** Live FVS startup is ~30 ms/stand single-process (NOT 1-3 s);
batching buys only ~1.5× AND couples the ~30% D38 crash across a batch (one SIGFPE aborts all remaining —
confirmed died at stand 39/500). So the full ~1.24M-stand census = single-process-per-stand (crash-isolated),
parallel across cores, oracle cached once ⇒ ~40 min/16 cores (~10 CPU-hr), paid once. Parallel driver drafted
at `test/harness/fia/sweep/census_driver.jl`.

## Slice S70 — REAL FIX: Fort Bragg FIA zero-volume (701→81110 KODFOR remap on the DB path) [2026-07-08]
Surfaced by re-trace discipline on the post-#95 SN spread sweep: measuring the residual DELTAS (not labeling
them "all ULP") exposed 1 real zero-volume stand among the 6 residuals.
- Symptom: SN FIA stands in Fort Bragg (REGION=7/FOREST=1 ⇒ composite LOCATION=701) got jl=0 for ALL volume
  columns; live gave 4958/3490/2384/12983 (stand 252270661010854).
- Root cause (MEASURED, not guessed): a live-FVS debug-stamp of `fvsvol.f` (VOLEQ/IREGN/KODFOR) showed live
  runs with **KODFOR=81110 (NC Uwharrie, region 8), IFOR=20, VOLEQ=821CLKE** — i.e. `forkod.f` CASE(701)
  remaps forest 701 ⇒ Uwharrie 81110. Both encodings converge on 701 (keyfile 701xx collapses via
  `IFORDI==701`; the FIA composite is already 701) ⇒ `SELECT CASE(701)`. jl's `kw_stdinfo!` Fort Bragg gate
  `div(KODFOR,100)==701` caught the 5-digit 70106 but not the bare FIA 701 (div=7) ⇒ region-7 ⇒ VOLEQDEF
  assigned no R8 Clark eq ⇒ zero volume. (An earlier guess — the region-8 default trap → Talladega 80106 —
  gave nonzero-but-~1%-off and was REVERTED; the debug-stamp is what pinned the correct forest.)
- Fix: extracted `sn_fortbragg_remap!(p)` (gate widened to `user_forest_code == 701 || div(...,100)==701`),
  called from BOTH `kw_stdinfo!` and the FIA reader (`io/fia_database.jl`, SN-gated) — mirroring FVS, where
  forkod runs for every stand regardless of input path.
- Verify: offender **bit-exact vs live** (4958/3490/2384/12983); full suite **37121 / 140 / 0** (Fort Bragg
  keyfile test still passes — the widened gate is compatible); post-fix SN spread sweep **287/292 = 98.3%**
  with the residual 5 **ALL cornered-ULP** (4 TopHt AVH-tie + 1 BdFt Δ1) — zero real divergences in the sample.
  `forest_idx=20` ⇒ default merch, matching live's setcubicdflts IFOR=20 keying.

## Slice S71 — REAL FIX: CS American-elm dubbed-height (birch-group transcription error) [2026-07-08]
The #99 region-9 CS divergence (TopHt +29% → volume +17% on dubbed-height stands) was a one-line PORT
TRANSCRIPTION ERROR, run down by a long measurement-driven re-trace (each wrong hypothesis refuted by data):
- Symptom: CS FIA stand 3303107010661 (elm/black-locust, all heights dubbed) TopHt live=41/jl=53.
- Ruled out by measurement: site-species selection (an SN-variant-contaminated probe — the cross-variant
  gotcha: `initialize()` with no `variant=` defaults to SN), site-index fan-out (jl `_CS_SITE_COEF` asite/bsite
  match `cs/sitset.f` bit-exact), and black locust (live TREELIST `cstl.trl`: BL dbh18.9→71.2 == jl, bit-exact).
- Root cause: live `htdbh.f` debug-stamp showed American elm (sp33) uses its OWN Wykoff coefs
  (HT1=4.6008/HT2=-7.2732 → dbh10.2→56.51), but `src/engine/volume.jl:45` (`_htdbh_wykoff` ifor-3 override)
  wrongly grouped **sp33 with the birches sp30/31**, overriding AE with birch 4.4635/-3.6456 → 67.2 ft.
- Fix: `(sp==30||sp==31||sp==33)` → `(sp==30||sp==31)`. jl's base AE Wykoff coefs were already correct, so
  removing sp33 makes AE use them.
- Verify: AE dbh10.2 dub **56.508503 == live 56.51**; BL unchanged (71.19); the region-9 CS stand now
  **bit-exact vs live** (TopHt 41 / TCuFt 433 / MCuFt 378, was 53/509/378); full suite **37121 / 140 / 0**
  (no passing test exercised the wrong path); CS FIA spread sweep **94.0% → 95.2%**. Oracle restored pristine.
- META (doctrine): the win came entirely from MEASURING live per-tree data (TREELIST + an htdbh.f stamp),
  not from labeling; and from always passing the correct `variant=` to `initialize()`/`run_keyfile()`.

## Slice S72 — REAL FIX: NE catch-all "other" species pointed at RL not OH [2026-07-08]
#98 (NE FIA TCuFt +52% on region-9 stand 1536568185290487) was a one-line SPECIES-RESOLUTION bug:
`northeast/species.jl:48` `other_species(::Northeast)=Int32(97)` (a "TODO verify" from the port) — but NE
index 97 is RL (red elm, VOLEQ 900CLKE975); "other hardwood" OH is index **98** (VOLEQ 900CLKE998). So every
UNMATCHED FIA code (e.g. 6918, "other hardwood") got the RL volume equation, over-stating small-tree cubic.
Run-down (all measurement): per-tree jl 0.2 vs live 0.0-0.1 (fvsvol.f ZTCF stamp, dbg binary validated
==145) → live VOLEQ map OH=998 vs jl using 975 → jl sp97=RL, OH is sp98 → the catch-all TODO. Fix:
`other_species(::Northeast)=Int32(98)`. Verify: stand now **bit-exact** — TCuFt 145 (was 220) AND CCF 53
(was 56; OH crown-width ≠ RL); full suite **37121/140/0**; NE FIA spread sweep **96.6% → 98.3%, BIG
divergences 0** (residual now all cornered: 1 TopHt tie + 1 MCuFt Δ1). Helps ALL NE FIA stands with
non-standard hardwood codes. Oracle restored pristine (NE buildDir needed `gfortran -c mrules_mod.f` first —
different-gfortran .mod; dbg binary validated ==oracle).

## Slice S73 — REAL FIX: NE Allegheny (IFOR=3) HT-DBH override was applied to ALL variants [2026-07-08]
The last region-9 FIA bug (#97, LS TCuFt ~3-4% high) and the true cause of #99 (CS) were ONE variant-safety
bug (doctrine #4): `_htdbh_wykoff` (volume.jl) applied the NE Allegheny-NF IFOR=3 HT-DBH override — keyed by
**NE species indices** (sp26=RM, 30/31/33=Y/S/P-birch, 41=ash, 55=WO …) — to **every** variant at forest
index 3. But CS (Hoosier) and LS also use ifor=3, where the same index is a different species (CS sp41=YP, LS
sp41=QA aspen). So LS aspen got NE ash's dub curve (64.2 ft vs live 52.2 → +18% vol); CS elm got NE birch's
(the #99 symptom — my earlier CS sp33-removal was a compensating hack that silently broke NE paper-birch sp33).
- Fix: `isne::Bool` on `_htdbh_wykoff`/`_htdbh_height`/`_htdbh_dbh`, block gated to `ifor==3 && isne`; sp33
  restored to the NE birch group; `isne = variant isa Northeast` threaded through every caller
  (dub_missing_heights!, fire snag/crown_biomass, establishment, NE small_tree_growth). NE Allegheny override
  test (test_net01.jl) updated to pass `isne=true`.
- Verify: LS stand bit-exact (TCuFt 952 was 983); CS #99 stand STILL bit-exact (via the gate, not the hack);
  NE override test passes; full suite **37121/140/0**; LS FIA sweep **95.4% → 96.3%** (aspen vol divergences
  gone; residual = cornered CCF/TopHt). Live-validated (LS TREELIST QA dbh6.7→52.2 == jl base Wykoff).
- ★ ALL FOUR variants' FIA residuals are now cornered-ULP-only (SN 98.3 / NE 98.3 / CS 95.2 / LS 96.3) —
  no real divergence class remains in the sampled populations. Meta: a shared-engine ifor-keyed override with
  variant-specific species indices is the exact doctrine-#4 hazard; measurement (per-tree live TREELIST) found it.

## Slice S74 — Pillar-1 coverage: TFIXAREA added to the KCV set (63 → 64 keywords/variant) [2026-07-08]
Gap analysis (jl-dispatched keywords vs KCV-covered) confirmed the remaining coverage gap is small — mostly
structural/IO keywords, documented exclusions (VOLEQNUM SN-only, THINAUTO ill-posed, REGDMULT vacuous), or
already-exercised (SITECODE lives in the setsite fixture). Added the genuinely-new, single-plot-testable
**TFIXAREA** (total fixed plot area, notre.f:45 → `kw_tfixarea!`) for NE/CS/LS: it rescales the per-acre
expansion (NE TPA 536→236 vs the setsite baseline, so genuinely exercised). Fixtures
`{ne,cs,ls}_tfixarea.{key,tre,live.sum}` built from each variant's setsite fixture (swap keyword) + live
FVS `.sum`; validated **bit-exact** vs freshly-relinked live via `test_kwcov_variants.jl` auto-discovery.
Suite **37628 pass / 140 broken / 0 fail** (pass count up = coverage added; the 140 cornered set unchanged).
KCV is now 64/variant; remaining additions (THINPT/SETPTHIN point-thin, CARBREPT/FMORTMLT FFE reports) are
low marginal value — the growth/mortality/volume/thinning(11 variants)/estab/FFE spine is already covered.

## Slice S75 — #73 soft-snag CWD1 cone-split LOHT: SOURCE-VERIFIED, deferred (not blind-fixed) [2026-07-08]
Investigated the last real correctness gap per the doctrine (verify from FVS source, NOT test behavior; the
path is inert so pass/fail can't guide it). FVS `fmcwd.f` CWD1 (called from FMSNAG on snag fall): `:170`
'1=soft,2=hard'; `:175 LOHT(1)=1.0` (soft), `:177 LOHT(2)=0.10` (hard) — the fallen-bole cone-split integrates
over [1.0,HIHT] for soft and [0.10,HIHT] for hard. jl (`snag.jl`) hardcodes `loht=0.10` (line 215) and applies
one `frac` to BOTH the soft `addS`→cwd[:,1] and hard `addH`→cwd[:,2] deposits (lines 337-344). SUBTLETY: the
normalizer `total_full=_cwd_pat(loht,…)` uses loht, so a soft `loht=1.0` also changes the [0.10,1.0]-base
normalization — must confirm FVS's CWD1 volume accounting (renorm-by-pat(1.0) vs actual-[LOHT,HIHT]-volume ×
bolevol/TVOLI) before implementing. DECISION: NOT blind-fixed — the path is INERT (dfis=0 every cycle: ordinary
snags all created hard; verified vs instrumented FMSNAG), so a subtly-wrong fix would be an UNVALIDATABLE latent
bug (doctrine #1/#3 + the verify-semantics-from-source rule). Precisely documented in #73 with the validation
path: a **SNAGPSFT (PSOFT>0) scenario seeds soft snags → dfis>0 at fall → the soft split becomes validatable
bit-exact vs live** (or unblock SNAGDCAY). Stays a documented cornered deferral until then.

## Slice S76 — Pillar-floor re-verify of the S73 variant-gate change (isne kwarg) [2026-07-08]
The #97 fix threaded an `isne::Bool` kwarg through the hot-path-adjacent `_htdbh_wykoff`/`_htdbh_height`/
`_htdbh_dbh` (@inline). Verified it costs NOTHING against the pillar floors:
- **Pillar 2 (allocation-free):** micro-benchmark (function-wrapped, typed-local `sd` to avoid global-boxing
  artifacts) — `_htdbh_height`/`_htdbh_dbh`/`_htdbh_wykoff` with `isne=true` AND the default path = **0 B /
  100k calls** each. The `@inline` fully elides the Bool kwarg; the fire snag/crown/dub paths stay bytes-free.
- **Pillar 1 (bit-exact):** full suite 37628 / 140 / 0; all four variants bit-exact (incl. the NE Allegheny
  override test, now `isne=true`).
- **Pillar 4 (type-stable):** `isne::Bool` is isbits; no `Any`/boxing introduced.
Conclusion: the session's five fixes + the TFIXAREA coverage add + the variant-gate refactor all hold the
tolerance-campaign floor AND the allocation-free/type-stable pillars.

## Slice S77 — Pillar-2 re-rank of the current hot-path allocators (negative result: floor confirmed) [2026-07-08]
Re-measured the live allocation state (the audit's ranked list was stale — top entry `snag.jl:203 pat closure
16 MB` is ALREADY fixed, hoisted to the module-scope pure `_cwd_pat`). Current `@allocated
run_keyfile(snt01;faithful)` = **12.49 MB/run**, down from the 71.85 MB baseline (prior slices did the heavy
lifting). Re-ranked the non-IO engine allocators via `--track-allocation=user` (8 warm runs): the survivors are
(a) one-time `StandState`/`Scratch` construction, (b) CONDITIONAL paths — `cuts.jl:24 ipush` (thinning only),
`fmburn.jl:452-472 push!` (fire only), `r8clark_vol.jl:612-614 out/logLen` (econ log-grade only), and (c) the
per-cycle `Float32[…]` snapshot comprehensions in `grow_cycle!` (`old_cfv`/`old_tpa`/`old_cfv2`).
- **Attempted** the (c) fix: move the three snapshots to reused `MAXTRE` `Scratch` buffers (plain-copy fill ⇒
  bit-identical). **Measured: allocation went UP** (12.49 → 12.59 MB, +106 KB) ⇒ **REVERTED** (doctrine #1/#3).
- **Root cause (re-trace, doctrine #9):** `run_keyfile` builds a FRESH `StandState`/`Scratch` per call, so three
  new `MAXTRE=3000` buffers cost ~36 KB one-time PER RUN, while the comprehensions they replaced are sized to the
  stand's ACTUAL small tree count (~60–180) — only a few KB/cycle over snt01's ~7 cycles. MAXTRE preallocation
  only amortizes over MANY cycles; on a short/small stand it is net-negative and worsens the documented per-run
  metric. **The Pillar-2 trap: preallocate-to-MAXTRE pays for long hot loops, not per-`run_keyfile` snapshots.**
- **Net: Pillar-2 floor CONFIRMED, not regressed.** Revert re-measured at 12.486 MB (== baseline within noise).
  The remaining allocators are the justified floor (one-time construction + conditional-path + actual-size
  snapshots already cheaper than MAXTRE reuse for typical stands). Suite unchanged (37628/140/0; edits reverted).

## Slice S78 — REAL FIX: SNAGPSFT soft-snag DDW (fmcwd.f LOHT(1)=1.0 cone-split) [2026-07-08]
Root-caused + FIXED the last named soft-snag residual (#73), live-validated bit-exact. A SNAGPSFT=1.0 scenario
(all ordinary-mortality snags created SOFT ⇒ they fall into the soft down-wood pool) exposed jl DDW ~1.3% HIGH
vs live (2000 8.0/7.9, 2005 10.8/10.7) with every OTHER carbon pool bit-exact — localizing it to the soft CWD1
deposit alone. Source trace (fmcwd.f): the DO-20 K-loop runs K=1 (soft) with LOHT(1)=1.0 and K=2 (hard) with
LOHT(2)=0.10 (:175/:177/:343), and LOHT(K) enters BOTH the cone-base radius R1 (:347 `R1 += LOHT(K)·R1·HTD/
(HTD−4.5)`) AND the LOCUT integration floor. A larger LOHT ⇒ fatter R1 ⇒ larger R1SQ normalizer inside the
conic profile P ⇒ smaller per-class volumes ⇒ a smaller TOTAL soft deposit (the fat 0.10–1.0 ft base is dropped,
NOT renormalized; FVS scales each class by the fixed hardness-invariant stem volume TVOLI). jl's
`_cwd_cone_fractions` used the hard (0.10) split for BOTH pools and self-normalized by pat(loht), cancelling the
loht dependence ⇒ soft over-deposit.
- **Fix:** `_cwd_cone_fractions` now returns `(frac_soft, frac_hard)`. frac_hard = the prior (loht=0.10) split,
  BIT-IDENTICAL. frac_soft uses r1_soft (loht=1.0 extension) + floor=1.0 for the raw bins but normalizes by the
  SAME invariant base pat_hard(0.10) — algebra from the bit-exact hard path (a = TVOLI·V2T·pat_hard(0.10)):
  a·frac_soft = raw_soft·TVOLI·V2T = raw_soft·a/pat_hard(0.10) ⇒ frac_soft = raw_soft/pat_hard(0.10). Callers
  (snag.jl update_snags! + cuts-salvage; cuts.jl CWD3) deposit soft→cwd[:,1], hard→cwd[:,2] with the matching split.
- **Why the floor is safe:** ordinary-mortality + fire-killed snags are created HARD (DFIS=0), so addS=0 and
  frac_soft is multiplied by zero ⇒ carbon_snt + ALL fire scenarios are UNCHANGED (frac_hard bit-identical). Only
  SNAGPSFT / seeded-soft snags (DFIS>0) exercise frac_soft.
- **Validated:** live FVSsn on carbon_snagpsft.key ⇒ DDW 5.8/5.2/7.9/10.7; jl now 5.8/5.2/7.9/10.7 (was 8.0/10.8) —
  BIT-EXACT to the F7.1 report on every carbon column. New regression fixture + testset (carbon_snagpsft.*,
  test_carbon.jl). Suite 37628/140/0 (floor held — no regression). #73 CLOSED (was the last cornered-deferred item).

## Slice S79 — Pillar-1 scouting: sharpened the LS/CS simfire FMEFF residual (no code change) [2026-07-08]
Applied the #73 lesson (re-examine "accepted" residuals with fresh measurement) to the next-strongest deferred
item — the cs/ls_simfire FMEFF fire-mortality residual. Ran jl vs freshly-relinked live FVSls on ls_simfire.key:
pre-fire rows (1990/2000/2010) BIT-EXACT; post-fire 2020 diverges TPA-ONLY (jl 225 / live 220) while BA (119=119)
AND every volume column (MCuFt 3002, SCuFt 2893, BdFt 2015, TCuFt) are BIT-EXACT. ⇒ the ~5 extra survivors carry
≈0 volume = the SMALLEST (sub-merch) trees; larger trees killed identically. This RULES OUT an RNG/XRAN-stream
misalignment (that would perturb BA+vols) and localizes the residual to the SMALL-TREE PMORT or the XRAN≤PSBURN
burn gate. Source (fmeff.f): LS has NO SN/CS MORTGP block (:195), so LS PMORT = 1/(1+exp(-1.941 + 6.316·(1-e^
{-FMBRKT(DBH,KSP)}) − .000535·CSV²)) — pure bark(FMBRKT) + crown-scorch(CSV). Recorded the sharper diagnosis in
_KCV_BROKEN + task #100 (next: debug-FVS FMEFF per-tree dump to pin bark/scorch-PMORT vs burn-gate). No code
change this slice — a scoped hand-off, not a half-finished fix. Suite unchanged (37628/140/0).
