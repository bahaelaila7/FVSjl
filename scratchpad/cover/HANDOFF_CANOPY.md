# COVER — CHUNK 2: CVIN reader + CANOPY report path (CVSHAP/CVCBMS/CVSUM/CVOUT)

Follows the merged beachhead (kt-variant-port HEAD **88522900**, `src/engine/cover.jl`
CoverState + cover_cvcw! + cover_init!).  This agent is git-isolated and CANNOT commit;
everything is STAGED under `scratchpad/cover/staged/`.  The main loop applies + wires +
commits from a real `kt-variant-port` checkout at 88522900.

================================================================================
## STATUS  —  PARTIAL (cover-canopy)
================================================================================
VALIDATED bit-exact-or-1ULP vs the LIVE oracle FVSem_g16 (measured, dump-replay):
- **CVIN reader** `kw_coverin!` — parses the COVER…END block, sets CoverState, schedules
  activity 900.  **INERT proven**: COVER-on vs no-COVER `.sum` DATA rows byte-identical
  (report-only — never perturbs growth/mortality/tree/.sum).
- **CVSHAP** (`vcovr/cvshap.f`) — crown-shape discriminant argmax → ISHAPE.  **351/351
  trees bit-exact** (integer argmax + EM species→11 map).
- **CVCBMS** (`vcovr/cvcbms.f`, COVOPT=2) — foliage biomass.  **337/351 hex-exact, 14 @
  exactly 1 ULP** (libm exp/log last-bit; cannot move the integer-rounded report).
- **CVSUM** (`covr/cvsum.f`) — 10-ft height-class canopy geometry (crown frustums →
  PROXHT/VOLXHT/CFBXHT + CRXHT/TXHT).  **239/240 hex-exact, 1 PROXHT value @ 1 ULP**.
- **CVOUT** (`covr/cvout.f`) — "CANOPY COVER STATISTICS" table emit, wired end-to-end
  (`run_keyfile`).  **TREES + COVER columns bit-exact vs oracle** at the inventory row.

All four kernels are covered by `test/unit/test_cover.jl` (15 assertions, passing).
`test_multicycle` 339/11 unchanged (edits inert for non-COVER runs).

### KNOWN GAP (does NOT block merge; report-only, honestly scoped)
The AREA/VOLUME/BIOMASS columns are NOT yet bit-exact end-to-end, for TWO reasons that
are **integration inputs, not COVER-logic errors** (the kernels themselves are proven):

1. **CRWDTH crown-width source (the blocker).**  CVCW reads FVS's `CRWDTH` array verbatim.
   For EM that array is filled by `base/cwidth.f → CWCALC` = **the WESTERN crown-width
   library `FVSem_buildDir/cwcalc.f`** (equations grouped by FIA code, IWHO=0 forest-grown).
   FVSjl has **only the eastern** crown-width library (`src/engine/crown_width.jl`,
   Bechtold/Bragg/Ek/Smith) — it returns the 0.5 floor for every EM species, and EM's
   growth path never needs a crown width (CCF uses the `em_tree_ccf` polynomial).  As a
   *placeholder* this chunk computes CRWDTH from **em/ccfcal.f MODE=2** (`_cover_crwdth_em`),
   which matches the true CRWDTH closely for merch trees (tree1 12.676 vs oracle 12.674)
   — enough that TREES + COVER round identically — but floors sub-inch trees differently
   (em/ccfcal gives cw≈1.2 where CWCALC floors at 0.5), which perturbs AREA/BIOMASS in the
   low height classes.  **Exact AREA/VOLUME/BIOMASS requires porting the western CWCALC
   crown-width library (IWHO=0), a separate subsystem** — this is the top next dependency.
2. **DDS timing at the inventory row.**  CVCBMS large-tree biomass uses
   `DDS=(2·D·DG+DG²)/FINT`.  At the cycle-0 (inventory) accumulation point some FVSjl
   `diam_growth` values are still 0 (dgf runs inside grow_cycle!), so DDS floors to 1e-4
   and their biomass diverges.  Fix = sample DG at the FVS GRADD phase (post-dgf) or seed
   the inventory DG from setup_growth!'s prediction.

Multi-cycle rows (years 2000/2010) additionally inherit the **known EM #206 OLDRN growth
straddle** (the `.sum` itself is cornered: MCuFt 1635 vs 1637), so COVER at later cycles is
cornered-by-inheritance, not by a COVER bug.

================================================================================
## APPLY  (from a clean kt-variant-port checkout at 88522900)
================================================================================
Staged files live under `scratchpad/cover/staged/`.

1. Full-file replacements:
   ```
   cp scratchpad/cover/staged/src/engine/cover.jl        src/engine/cover.jl
   cp scratchpad/cover/staged/test/unit/test_cover.jl    test/unit/test_cover.jl
   cp scratchpad/cover/staged/test/unit/cover_cvshap_g16_dump.txt test/unit/
   cp scratchpad/cover/staged/test/unit/cover_cvcbms_g16_dump.txt test/unit/
   cp scratchpad/cover/staged/test/unit/cover_cvsum_g16_dump.txt  test/unit/
   ```
   (`test/unit/cover_cvcw_g16_dump.txt` and the `include("unit/test_cover.jl")` line in
   `test/runtests.jl` are already present from the beachhead — unchanged.)

2. Three in-place wiring edits — apply `scratchpad/cover/staged/wiring.diff`
   (`git apply scratchpad/cover/staged/wiring.diff`) OR make them by hand:
   - **src/engine/keyword_dispatch.jl** (after the `WSBW` line, ~2447):
     ```julia
     elseif kw == "COVER";    kw_coverin!(s, rec, kr)    # COVER canopy/understory block (CVIN)
     ```
   - **src/io/summary.jl**: add `cover_year0 = 0` beside `prev_increment` (~211), then after
     the `cycle_hook` line (~225) in the `for c in 0:ncyc` loop:
     ```julia
     if s.cover !== nothing && s.cover.active
         compute_density!(s)
         c == 0 && (cover_year0 = Int(r.year))
         cover_fint = cycle_period_at(s.control, c == 0 ? 0 : c - 1)
         cover_accumulate!(s.cover, s, r.year, cover_year0, cover_fint)
     end
     ```
   - **src/engine/simulate.jl**: in `run_keyfile`, right after the `write_carbon_report_block`
     call (~873):
     ```julia
     (s.cover !== nothing && s.cover.active) &&
         cover_report(s.cover, out, String(sid), mid, strip(s.control.title))
     ```
   The state field `s.cover::Union{AbstractCoverState,Nothing}` (state.jl:985) and the
   `include("engine/cover.jl")` (FVSjl.jl:338) are already in place from the beachhead.

## VERIFY
```
JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project -e 'using FVSjl; using Test; include("test/unit/test_cover.jl")'   # 4 testsets pass
JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project -e 'using FVSjl; using Test; include("test/integration/test_multicycle.jl")'  # 339/11
```
INERT proof (report-only): run a keyfile with vs without the COVER block; the `.sum` DATA
rows are byte-identical (only the appended COVER report differs).

================================================================================
## ORACLE / REPRODUCTION  (all under scratchpad/cover/)
================================================================================
- `em_cov.key` / `em_cov.tre` — EM COVER keyfile (CANOPY+SHRUBS+COVER).  Run the shipped
  western binary (NO relink): `printf "em_cov.key\nem_cov.tre\n" | /workspace/.emwork/FVSem_g16`
  → `fort.16` (exit 10 = normal).  `em_cov.out` is the golden.  The CANOPY COVER STATISTICS
  table is identical with/without SHRUBS (LBROW does not touch the canopy arrays — verified).
- g16 instrumented dumps: `instr/{cvshap,cvcbms,cvsum}_instr.f` (copies of the pristine
  vcovr/covr sources + a single WRITE per routine) built by `instr/relink_canopy.sh`
  (3-way single-.o swap: g16obj minus cvshap/cvcbms/cvsum.o plus the instrumented .o +
  C objs + isoc23_shim.o).  Instrumented report proven byte-identical to clean (only the
  invocation-filename line differs).  Dumps: `cv{shap,cbms,sum}_g16_dump.txt`.
- Replay harnesses: `kernels.jl` (CVSHAP+CVCBMS), `cvsum_replay.jl` (CVSUM geometry).
- PRISTINE: `/workspace/ForestVegetationSimulator/{vcovr,covr}/*.f` untouched (verified 0
  instrumentation lines); scratch binaries removed.

================================================================================
## NEXT  (dep-ordered)
================================================================================
1. **Western CWCALC crown-width library** (`base/cwcalc.f` / `FVSem_buildDir/cwcalc.f`,
   IWHO=0 forest-grown, FIA-code grouped) → the real `CRWDTH`.  This UNBLOCKS bit-exact
   AREA/VOLUME/BIOMASS for the CANOPY table (and is reusable by all COVER variants).  Swap
   `_cover_crwdth` to call it; keep the em/ccfcal placeholder only as a fallback.  Dump
   CRWDTH per tree from g16 (instrument base/cwidth.f `CRWDTH(I)=CW`) and dump-replay.
2. **DDS/DG timing** at the inventory row — sample `diam_growth` at the FVS GRADD phase so
   large-tree CVCBMS biomass matches at cycle 0.  Then the full 5-row canopy table should be
   bit-exact at cycle 0 (later cycles cornered via the #206 growth straddle).
3. **CANOPY AND SHRUBS SUMMARY (case-3 / canopy-only)** table in CVOUT — the overstory
   summary row (STDHT/TPCTCV/TOTBMS/SDIAM/TRETOT/ICVAGE are already accumulated in CoverRow;
   only the emit formatting is unported).  Also post-thin (ITHN=2 / LTHIN) accumulation.
4. **Other variants' CVSHAP/CVCBMS species maps** — the MAP** DATA blocks for the remaining
   variants are in `vcovr/cvshap.f:157-284` (only EM + NI transcribed so far in `_cover_spmap`),
   and each variant's crown-width (dep #1).
5. **Shrub half** (CVBROW/CVSCON/CVBCAL/CVCLAS) + `pg/` serialization — the later chunks
   (unchanged from the beachhead NEXT list).
