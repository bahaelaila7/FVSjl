# COVER extension — SCOPE + beachhead apply/wire recipe

Staged in `scratchpad/cover/`. This agent ran isolated in an FVS-Fortran worktree
and CANNOT commit to FVSjl (git-isolation). The main loop applies + wires + commits
from a real `kt-variant-port` checkout (base `a7724d42`).

Source: `covr/` (11 f, canopy/report core), `vcovr/` (3 f, shrub browse/biomass/shape),
`pg/` (14 f, parallel-processing state serialization — see Q3, mostly out of scope).

================================================================================
## STEP-1 SCOPE (7 questions, file:line evidence)
================================================================================

### 1. What COVER computes — REPORT-ONLY (does NOT modify the tree record)
COVER produces an understory/canopy-cover text report; it never feeds back into tree
growth or mortality. Evidence:
- `grep -nE '^\s*(DBH|HT|PROB|DG|WK[0-9])\(.*=' covr/*.f vcovr/*.f` → **ZERO** writes
  to any tree array. The routines only READ the tree list.
- `covr/cvout.f:6` "PRINTS (PRE- AND POST-THIN) **COVER** STATISTICS BY CYCLE. CALLED
  FROM **MAIN** ONCE PER STAND, AT THE END OF PROJECTION." Writes to unit `JOSHRB`.
- Two halves:
  - **Canopy cover** (`covr/cvcnop.f` driver): CVCW (crown width/area) → CVSHAP (crown
    shape → profile area/volume) → CVCBMS (foliage biomass) → CVSUM (bin into 10-ft
    height classes) → CVCLAS (successional-stage code). Report = "CANOPY COVER
    STATISTICS" (TREES/COVER/AREA/VOLUME/BIOMASS per 10-ft class).
  - **Shrub / browse** (`vcovr/cvbrow.f` driver → `covr/cvscon.f` 35KB shrub-cover
    model + `covr/cvbcal.f` calibration): habitat-type-driven shrub cover, dormant
    twigs, successional stage. Report = "SHRUB STATISTICS", "CANOPY AND SHRUBS
    SUMMARY", "SHRUB-SMALL CONIFER COMPETITION".
⇒ "INERT" = no-COVER run byte-identical. "VALIDATED" = report text matches oracle.

### 2. Keywords
- **Base keyword** `COVER` = `base/keywds.f:32` TABLE index **46**. Dispatched at
  `FVS<v>_buildDir/initre.f:2265` "OPTION NUMBER 46: COVER" → `CALL CVIN` (:2268).
  `CALL CVINIT` at `initre.f:187`. (NOTE: the 36 "COVER" refs already in FVSjl src are
  the BASE crown-cover REPORT — unrelated; this is a separate extension.)
- **Sub-keywords** read by `covr/cvin.f` (TABLE at cvin.f:103), option→driver:
  1 END · 2 CANOPY (`LCNOP`, COVOPT) · 3 SHRUBS (`LBROW`; SAGE/IHTYPE/IPHYS/IDIST) ·
  4 SHRBLAYR / 5 SHRUBHT / 6 SHRUBPC (shrub calibration) · 7 DEBUG (no-op) ·
  8 NOCOVOUT (`LCOVER=F`) · 9 NOSHBOUT (`LSHRUB=F`) · 10 NOSUMOUT (`LCVSUM=F`) ·
  12 COVER (`OPNEW(KODE,IDT,900,...)` schedules activity **900** at cycle IDT; field 2
  = output unit JOSHRB) · 13 SHOWSHRB · 14 CVNOHEAD.
- **Activation gate** `covr/cvgo.f`: `LCVATV = LCOV .OR. OPSTUS(900…)>0 .OR.
  OPEVAC(900…)>0`. So COVER is "active" iff sub-keyword 12 COVER scheduled activity 900.

### 3. Linkage — ACTIVE in western binaries, STUBBED in eastern
- Real `cvin.o` compiled into **18** western variant buildDirs:
  ws ak bm ca cr ec op pn so tt wc ut ci em kt ie oc nc.
- `base/excov.f` stub (`excov.o`, all ENTRY CVINIT/CVIN/CVGO/CVCNOP/CVBROW/CVOUT/
  CVKEY/CVGET/CVPUT/CVACTV returning inert; CVGO sets LACTV=.FALSE.) linked instead in
  **6** eastern variants: ls ne on sn bc cs.
⇒ Oracle needs **NO relink** — the shipped western binary already fires COVER on the
  base keyword. (Contrast the insect models, which needed a stub→real relink.)

### 4. Driver + call sites
- `base/fvs.f:288` `CALL CVGO(LCVGO)`; `:294` `CALL CVBROW(.FALSE.)`; then `CALL CVCNOP`
  — initial (cycle-0) cover stats.
- `base/gradd.f:197` `CALL CVGO`; `:200` `IF(LCVATV) CALL CVBROW(.FALSE.)` — per-cycle
  post-growth/mortality.
- `base/grincr.f:303` `CALL CVGO`; post-thin `IF(ONTREM(7)>0 .AND. LCVATV) CALL
  CVBROW(.TRUE.)/CVCNOP(.TRUE.)`.
- `base/fvs.f:456` `CALL CVOUT` — end-of-stand report emit.
So COVER is BOTH a per-cycle accumulator (CVBROW/CVCNOP fill CVCOM arrays each cycle,
gated by LCVATV) AND an output-time emitter (CVOUT). Report-only throughout.

### 5. RNG — NONE
`grep -i 'RANN|RANDOM|ran(|URAND|iseed'` over covr/ vcovr/ pg/ → no COVER RNG. The
only hits are `pg/getstd.f:126 / putstd.f:128` serializing the ENGINE seed as stand
state (parallel-processing), not a COVER draw. **No LCG to port.**

### 6. Oracle — ACTIVE in shipped western binary (no relink)
- `/workspace/.emwork/FVSem_g16` (full g16 rebuild) and `/workspace/.ktwork/FVSkt_clean`
  both fire COVER. `scratchpad/cover/{cover_on.key, em_cov.key}` are working keyfiles;
  `cover_on.out` / `em_cov.out` are the golden reports. `cover_off.out ≡ stock` (no
  COVER section) confirms INERT at the oracle.
- g16 single-.o-swap dump recipe (used for the beachhead, keeps everything pristine):
  compile an instrumented `cvcw.f` to a SCRATCH `.o`, relink listing
  `/workspace/.emwork/g16obj/*.o` MINUS `cvcw.o` PLUS the scratch `.o` + C objs +
  `/workspace/.crwork/isoc23_shim.o` → scratch binary. See `relink.sh`. Instrumented
  COVER report proven byte-identical to clean except the run-timestamp line
  (`diff em_cov_instr.out em_cov.out`).

### 7. Variant coverage
- Canopy-cover half (CVCW/CVSHAP/CVCBMS): `covr/cvcw.f` now sets `TRECW(I)=CRWDTH(I)`
  (reads the base open-grown crown width — generic across all variants). `vcovr/cvshap.f`
  MAPS each variant's species onto the ORIGINAL 11 N-Rockies species for crown-shape
  coeffs (cvshap.f:153-320, per-variant blocks incl. CA/OC/OP mappings).
- Shrub half (CVSCON/CVBROW): driven by habitat type `IHTYPE`/physiography/disturbance
  — N-Rockies-centric coefficient tables. All 18 western variants LINK it, but the
  shrub-model coefficients are the Moeur/N-Rockies set regardless of variant.

================================================================================
## BEACHHEAD — VALIDATED (dump-replay bit-exact)
================================================================================
Ported the **CVCW crown-area first computation** (`covr/cvcw.f`):
`TRECW(I)=CRWDTH(I)`, `CRAREA=(Σ TRECW(i)²·PROB(i))·0.785398`.
- Kernel: `cover.jl :: cover_cvcw!` (Float32 running sum in tree order; literal
  `0.785398f0`). `cover_init!` = CVINIT; `CoverState` = CVCOM subset.
- Validation: `test_cover.jl` replays `cvcw_g16_dump.txt` (instrumented FVSem_g16, EM
  emt01 stand, 2 cycles → 3 CVCW invocations after tripling) → **7/7 hex-exact**:
  icyc0 26605.12 `46cfda3d`, icyc1 30998.373 `46f22cbf`, icyc2 34740.656 `4707b4a8`.
  Run: `JULIA_DEPOT_PATH=/workspace/.julia_depot julia test_cover.jl`.
- INERT argument: kernel only runs behind `active` (LCVATV); CoverState defaults inert;
  a run without the COVER keyword never constructs it — mirrors the insect gate.

NOT yet ported (do not claim done): the CVIN keyword reader, CVSHAP/CVCBMS/CVSUM/CVOUT
report path, and the entire shrub half.

================================================================================
## APPLY
================================================================================
1. `cp scratchpad/cover/cover.jl       src/engine/cover.jl`
2. `cp scratchpad/cover/test_cover.jl  test/unit/test_cover.jl`
   and `cp scratchpad/cover/cvcw_g16_dump.txt test/unit/fixtures/cover_cvcw_g16_dump.txt`
   (adjust the `include`/path in test_cover.jl to the fixtures location; drop the
   standalone `include("cover.jl")` shim once wired).

## WIRE — mirror the insect-model seams (WSBW is the nearest precedent)
### (1) src/core/state.jl
- After `abstract type AbstractWsbweState end` add:
  ```julia
  abstract type AbstractCoverState end   # COVER understory/canopy-cover REPORT extension; concrete CoverState in engine/cover.jl
  ```
- In `mutable struct StandState`, immediately AFTER the `wsbwe::Union{...}` field add:
  ```julia
      cover::Union{AbstractCoverState,Nothing}   # COVER report extension; nothing until the COVER keyword schedules activity 900 (report-only, INERT for growth/mort)
  ```
- Add ONE more trailing `nothing` (after the wsbwe one) in the `StandState(variant;
  faithful)` constructor call.
### (2) src/FVSjl.jl
- After `include("engine/wsbwe.jl")` add:
  ```julia
  include("engine/cover.jl")   # COVER understory/canopy-cover REPORT extension (covr/vcovr) — beachhead: CVCW crown-area dump-replay bit-exact; report path deferred
  ```
### (3) src/engine/keyword_dispatch.jl
- After the `elseif kw == "WSBW"; …` line add:
  ```julia
          elseif kw == "COVER";  kw_coverin!(s, rec, kr)   # keywds.f opt 46 → CVIN sub-keyword block; sets s.cover (report-only, INERT)
  ```
  (kw_coverin! is a NEXT chunk — until it exists, dispatch the COVER keyword to a
  skip-block reader so the keyfile parses; do NOT leave it unhandled.)
### (4) src/engine/simulate.jl
- NO growth/mortality seam (report-only). The per-cycle accumulator (CVBROW/CVCNOP)
  and the end-of-stand CVOUT emit are a separate report seam, added when the report
  path is ported (NEXT chunks). Until then COVER is byte-identical to stock.
### (5) test/runtests.jl
- After `include("unit/test_wsbwe.jl")` add `include("unit/test_cover.jl")`.

## VERIFY (from the real checkout)
`JULIA_DEPOT_PATH=/workspace/.julia_depot julia --project -e 'using Pkg; Pkg.test()'`
plus a no-COVER western .sum run byte-identical to baseline (INERT proof).

================================================================================
## NEXT (dep-ordered)
================================================================================
1. **CVIN keyword reader** (`covr/cvin.f`) → `kw_coverin!` + extend CoverState with the
   sub-keyword flags/calibration arrays (AVGBHT/AVGBPC/SHRBHT/SHRBPC). Validate INERT
   (no-COVER .sum byte-identical) + that a COVER block parses without perturbing growth.
2. **Canopy report path**: CVSHAP (`vcovr/cvshap.f`; variant→11-species map) → CVCBMS
   (`vcovr/cvcbms.f`) → CVSUM (`covr/cvsum.f`, 10-ft height-class binning + rounding) →
   CVOUT (`covr/cvout.f` "CANOPY COVER STATISTICS" table). Dump-replay each kernel at
   E22.13/hex via the g16 single-.o-swap recipe (relink.sh), then match the emitted
   integer table. This completes the canopy half (CVCW already done).
3. **Shrub half**: CVBROW (`vcovr/cvbrow.f`) → CVSCON (`covr/cvscon.f`, big) → CVBCAL
   (`covr/cvbcal.f`) → CVCLAS (`covr/cvclas.f`) → CVNOHD. Habitat-type driven; largest
   remaining piece.
4. **pg/** (get/put std, cvget/cvput): parallel-processing state serialization — needed
   only for multi-stand parallel/FVSSTAND stashing, LAST or skip for single-stand parity.

Golden/fixtures in scratchpad/cover/: cover.jl, test_cover.jl, cvcw_g16_dump.txt,
cvcw_instr.f (instrumentation diff), relink.sh (g16 recipe), replay_cvcw.jl,
{cover_on,cover_off,em_cov,em_cov_instr}.out (oracle goldens), *.key/*.tre.
