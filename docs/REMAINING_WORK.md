# FVSjl — remaining work and blockers

Status snapshot (suite 4223). The natural-process core + most management/disturbance keywords are
ported and validated bit-exact vs live Fortran. What remains, grouped by the **nature of the blocker**
(not just "todo") so the path for each is clear.

## A. FFE surface-fuel dynamics → grown-cycle Stand Carbon Report (8/9 COLUMNS BIT-EXACT)
Inventory cycle bit-exact (all columns). **Grown-cycle: 8 of 9 report columns now reconcile bit-exact**
vs the 4-cycle Fortran baseline (carbon_jenkins) — see test_carbon.jl. Ported + validated this session:
- ✅ `FMCWD` decay (fuel_decay.jl) · ✅ `FMCADD` litterfall + woody breakage (fuel_additions.jl) ·
  ✅ per-cycle driver `ffe_fuel_update!` · ✅ periodic-mortality → snags (mortality.jl) · ✅ dead coarse-
  root pool BIOROOT (Below-Dead column).
- ✅ THE BUG FIX: `crown_biomass` missed the V2T/=2000 rescale (fmvinit.f:1094) → woody cone weights
  ~2000× too large; root-caused via a Fortran XV dump (instrument fmcrow.f/fmcrowe.f, rebuild, dump,
  revert — the `DEBUG` keyword segfaults the stripped binary). One-line fix `sg = v2t * _FM_P2T`.
- ✅ RECONCILE: above/merch/below-live, **below-dead**, forest floor (every cycle), DDW (1990/95/2005),
  shrub — all bit-exact; Stand-Dead populated.
- ✅ **CWD2B crown-debris scheduling (FMSCRO) — DONE & bit-exact for Stand-Dead** (`cwd2b[4,6,60]`
  state, the TFALL table, the at-death crown spread, the snag = bole-only split). Landed the STAND-DEAD
  column bit-exact (5.2/4.5 = bole 3.72 + crown 1.46), validated vs an instrumented Fortran dump.
- ⛔ **REMAINING — DDW post-mortality column = the crown-LIFT term, FULLY CHARACTERIZED this session**
  (every piece validated; only the in-loop integration remains):
  - ✅ SOURCE: the crown-LIFT (`X·CROWNW·TPA·P2T`, fmcadd.f:95-102) — the lower crown shed as the crown
    base rises — is the dominant post-mortality down-wood addition (~0.39 t/ac/yr, instrumented FMCADD
    dump); breakage/cwd2b are minor. (Ruled out: crown-lift-negligible [buggy dump], CWD I/L structure,
    ordering swap, crown magnitude — all via dumps.)
  - ✅ FORMULA: `crown_lift_rate` (X = (NEWBOT−OLDBOT)/OLDCRL/CYCLEN, FMSDIT) — implemented + unit-tested
    (test_crown_lift.jl, 7 assertions).
  - ✅ PLUMBING: previous-cycle per-tree crown via `prev_height`/`prev_crown_len`/`prev_dbh` TreeList
    fields in `_TREE_VEC_FIELDS` (ride through compaction/tripling via `copy_tree!` — clean, no parallel
    arrays); cycle-start snapshot via `fmoldc!`. Validated to fire.
  - ✅ PLACEMENT: POST-grow (after `grow_cycle!`) — fires the lift at the right cycles (DDW jumps at
    2000 AND 2005), using this cycle's crown-base rise + the cycle-START crown (OLDCRW, from `prev_dbh`).
  - ⛔ ONLY REMAINING: apply it PER-YEAR inside a post-grow annual loop so each year's lift DECAYS over
    the rest of the cycle (a lump-sum `cyclen·X·CROWNW` overshoots ~1.5×; the per-year-with-decay
    magnitude hand-checks to the exact gap, ~1.76 mt/ha). = the FFE-update-AFTER-grow restructure: one
    post-grow annual loop carrying decay + litterfall + breakage + cwd2b-flow + crown-lift together
    (cycle-start crown for live terms, post-grow X for the lift), coupling to the cwd2b/Stand-Dead
    death-dating. A focused integration behind test_fire.jl, not an end-of-session rush. See
    FFE_FUEL_DYNAMICS_chunk_plan.md for the full validated trail.
- ✅ **the `.out` Stand-Carbon-Report writer — DONE** (`write_carbon_report`, carbon.jl): header block +
  per-row FORMAT byte-for-byte vs the Fortran `.out` (CARBREPT); inventory row bit-exact in every column,
  grown rows track within the LP growth tail. test_carbon.jl (18 assertions).
- **Blocker:** none external; the remaining item is the grow_cycle! integration (timing) + the .out writer.

## B. Validation-blocked by the stripped ground-truth binary
The rebuilt `/tmp/FVSsn_new` is a **stripped DBS build**: the DATABASE block accepts only
`DSNOUT/SUMMARY/COMPUTDB/TREELIDB`; every other DBS sub-keyword errors. So there is no ground-truth
SQLite to diff against (see `fvsjl-ground-truth-binary-limits`).
- **~13–14 DBS output tables** — `FVS_Carbon`, `FVS_StrClass`, `FVS_Fuels`, `FVS_PotFire`,
  `FVS_Mortality`, `FVS_Consumption`, `FVS_SnagSum`, `FVS_DWD`, `FVS_BurnReport`, `FVS_Hrv_Carbon`,
  `FVS_Summary2`, the 2 econ tables, … **Blocker:** need a fuller SN binary compiled with all DBS
  modules. NB the UNDERLYING data of several is validatable via the `.out` **text** reports (carbon via
  `CARBREPT` is done; `FUELREPT`/`POTFIRE`/`MORTREPT` text reports likely work too) — only the SQLite
  *table emission* is unverifiable here.

## C. `.sum`-inert + FFE-coupled (low value, hard to validate)
- **FINTM** (GROWTH mortality measurement period) — its only effect (cratet.f:486) scales INPUT
  dead-tree PROB into the FFE snag model. `.sum`-inert; needs the unported input-mortality→snag path
  (`FMSSEE`; FVSjl's `add_snag!` only fires from fire-kills); validatable only via FFE snag DBS (B).
- **FFE snag-input path (`FMSSEE`)** — seeds the initial snag list from input dead-tree records. Same
  snag-DBS validation blocker.

## D. Out-of-scope for the SN variant
- **8 insect/disease models** (MPB/DFB/DFTM/WSBW/BRUST/MISTOE/RDIN/ANIN/RRIN) — linked into SN as
  **`ex*.f` NO-OP stubs**; the real models (mpbmza/dfbmza/…) are **absent** from the SN build. Keyword-
  recognized but inert without infected/host input. **Blocker:** there is nothing to port for SN; a
  faithful port/validation would need a different variant binary + infected-stand scenarios.

## E. Long-tail accuracy (a residual, not a missing extension)
- **38-species single-species DG-calibration tail** (~1–2% TPA for non-snt01 species). snt01's species
  (22/27/33/65/89) are bit-exact; single-species stands exercise each species' calibration regression
  harder, exposing per-species COR/OLDRN residuals. **Blocker:** a long tail — each species needs a
  per-tree calib decomposition (ran/HTCALC/COR) vs oracle; partly transcendental-ulp. snt01 is the gate.
- **DGSCOR per-record cubic-volume ±0.03% drift** — documented as **irreducible** (transcendental-ulp
  through the bounded redraw), not a bug. See `DIVERGENCES.md`.

### Update — Stand-Dead bit-exact (8/9); DDW root cause = within-cycle deaths-timing
The snag stem-volume bole + CWD2B crown model landed the STAND-DEAD column bit-exact (5.2/4.5 =
bole 3.72 + crown 1.46), validated piece-by-piece against an instrumented Fortran dump. So 8 of 9
report columns now reconcile bit-exact. The ONE remaining: **post-mortality DDW** (Fortran 3.8/8.0
@2000/2005 vs FVSjl 2.1/6.1).

**ROOT CAUSE (pinned via a per-year Fortran dump of fmmain.f — supersedes the earlier crown-lift
guess, which an instrumented dump RULED OUT as negligible, ~0.0007 t/ac/yr).** The Fortran down-wood
grows GRADUALLY within each cycle — over 1995→2000 the per-year DDW is 2.56, 2.82, 3.04, 3.23, 3.40 —
because FVS spreads the cycle's mortality across the annual fuel loop, so EARLY-period deaths' crown
debris (CWD2B) flows into DDW *within that same cycle* (+0.84), while the still-in-waiting crown at the
boundary is the 1.46 already validated for Stand-Dead. FVSjl books ALL the cycle's mortality at the
boundary and runs the fuel flow UPD-before-GROW, so its crown debris flows a cycle LATE (the addition
lands at 2005, not 2000) — even though the TOTAL snags (Stand-Dead) reconcile. The report reads the
boundary value BEFORE the next mortality booking (carbon DDW 3.8 mt/ha @2000 = 3.39 t/ac = the year-
1999 dump value 3.40), confirming the read point.

So DDW bit-exact is NOT a missing addition term — it's the **FFE main-loop ordering**: the cycle's
mortality must be scheduled into CWD2B and then flowed through the annual loop WITHIN the same cycle
(FVS: GRINCR → FMSADD → annual FMSNAG/FMCWD/FMCADD), with the Stand-Dead crown then being the un-fallen
remainder. That is exactly the `grow_cycle!` hot-path wiring of `ffe_fuel_update!` (with the at-death
crown snapshot) — the integration step already on the list, now with the exact target behavior pinned.
A coupled refactor (touches the validated fire snag path via fmburn!); do it behind test_fire.jl, not
an end-of-session rush. Then the `.out` Stand-Carbon-Report writer.
