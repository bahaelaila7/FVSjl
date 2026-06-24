# FFE surface-fuel dynamics — chunk plan (the grown-cycle carbon remainder)

The **Stand Carbon Report** is now bit-exact at the **inventory cycle** (every column: Jenkins live
above/merch/below, snag, DDW ×0.5, forest floor ×0.37, shrub/herb via FULIV2, total — see
`test_carbon.jl` + `carbon.jl`). The remaining gap is the **grown-cycle** FFE pools: the Fortran
report decays the down-wood/litter/duff each cycle (e.g. DDW 3.8→2.5, Floor 9.1→6.6 from 1990→1995)
and accumulates new dead wood from tree mortality. That is a distinct FFE **surface-fuel-dynamics
subsystem**, not part of the carbon report itself — this plan scopes it from the Fortran flow.

## What FVSjl already has
- `fmcba!` (fire/fmcba.jl) — the per-cycle FFE cover-type + **initial** dead-fuel loading (FUINI) and
  live herb/shrub loading (FULIV / FULIV2 override). It loads `fire.cwd` once (`fuels_init`) and
  re-sets `fire.flive` each call. ✅ Correct for the inventory cycle.
- `fuel_model.jl` (FMCFMD3 dynamic fuel model), `fmburn.jl` (the fire-event behavior), `snag.jl`
  (FMSFALL snag falldown + hard→soft decay). These run **only on a fire event**.

## What's missing (the subsystem to port) — ⚠ DECAY + ADDITIONS ARE COUPLED, not separable
A first instinct is "just port the decay" — but the constants show that fails on its own (see below).
The grown-cycle DDW/Floor is `decay − additions` and BOTH must land together to validate.

1. ✅ **DONE (decay routine) — `FMCWD`** ported to `fire/fuel_decay.jl` (`fmcwd!`) + DKR/PRDUFF consts,
   unit-tested (`test_fuel_decay.jl`): duff persists (0.002/yr), litter crashes (0.65/yr, awaiting
   litterfall), DDW decays 3.42→1.91 t/ac per 5-yr cycle (the ~0.3 gap to the report's 2.5 is the
   chunk-2 woody-breakage additions). End-to-end report validation still pending chunk 2 (the coupling).
   For decay class L=1..4, size J:
   - duff (size 11): `cwd[11,1,L] *= (1−DKR[11,L]·1.1)^NYRS` (soft); `cwd[11,2,L] *= (1−DKR[11,L])^NYRS` (hard).
   - woody J=1..10: decayed amount `AMT = cwd[J,k,L]·(1−(1−DKR·{1.1 soft})^NYRS)`; `cwd[11,2,L] += AMT·PRDUFF[J,L]`
     (a fraction to duff); then `cwd[J,k,L] *= (1−DKR·{1.1 soft})^NYRS`; then hard→soft transfer (J<10)
     `TOSOFT = clamp(NYRS·ln(1−DKR[J,L])/ln(0.64),0,1)·cwd[J,2,L]` moved 2→1.
   - **DKR / PRDUFF constants (sn/fmvinit.f:70-115, → a CSV):** DKR[1:9,1]=0.11; DKR[·,2]= (0.11,0.11,0.09,
     then 0.07 for 4:9); DKR[1:9,3:4]=DKR[1:9,2]; **DKR[10,·]=0.65 (litter), DKR[11,·]=0.002 (duff)**;
     PRDUFF[·,·]=0.02. NYRS = cycle length.
   - ⚠ **Why decay alone DOESN'T reconcile:** litter (size 10) at DKR=0.65/yr ⇒ `(1−0.65)^5 ≈ 0.005` — it
     crashes to ~0 in one 5-yr cycle. Yet the Fortran report's Floor only goes 9.1→6.6. So litter MUST be
     replenished by annual **litterfall** (chunk 2). Porting FMCWD without it makes Floor far too low ⇒
     NOT independently validatable against the report. (Initial split for carbon_jenkins: FUINI 160s =
     litter 4.90 + duff 6.03 = 10.93 t/ac ⇒ ×0.37×2.2417 = 9.1; duff barely decays at 0.002/yr.)
   ✅ **2a DONE — litterfall** (`fmcadd_litterfall!`, fuel_additions.jl): `foliage·TPA/LEAFLF·P2T` →
   litter, per tree into its `dkr_cls`. The FFE update runs ANNUALLY (fmmain.f:226-259, NYRS=1) — the
   year loop `fmcwd!(1)+fmcadd_litterfall!` ×NYRS with the crown held at the cycle start reconciles the
   grown-cycle **Forest Floor BIT-EXACT** (carbon_jenkins 1990→1995: 9.1→6.6 = Fortran), which also
   implicitly validates crown_biomass foliage (FMCROWE). LEAFLF/dkr_cls are already in
   fire_species_props.csv. **⛔ 2b BLOCKED — woody crown breakage → DDW**: the FMCADD breakage logic
   (`LIMBRK·CROWNW(SIZE)·TPA·P2T`) is correct, but the `crown_biomass` WOODY components (`xv[2..6]`)
   are in FMCROWE "FFE-internal units, not literal tons" (an 8" loblolly gives 8.7/545/14394/20484 —
   not tons) → DDW ≈ 120 vs the report's 2.5. So **2b first needs the FMCROWE woody-component
   validation/fix (the F5/F6 crown-biomass chunk)**; foliage is already right, only the woody side is
   blocked. DDW at 1995 is 2.1 (decay only) vs Fortran 2.5 — the 0.4 gap is exactly these additions.
   - **Diagnosis (this session):** the bug is localized to the UMBTW **bole-tip cone weights** `u1..u4`
     in `crown_biomass` (the `sg·vol/P2T` terms, lines ~150-173) — they dominate `xv[3..5]` and are
     ~50× too large (an 8″ LP gives a 1-3" component of ~7 tons after ×P2T, vs a whole-tree biomass of
     ~0.5 t). Since `CROWNW(I,J)=XV(J)` directly (fmcrow.f:197) and Fortran's DDW is 2.5, Fortran's XV
     woody must be SMALL → it's a FVSjl cone-weight units bug. **To fix:** get the Fortran XV reference
     for the carbon_jenkins trees, then correct the cone weights. ⚠ The `DEBUG` keyword (which would
     dump `CROWNW=` via fmcrow.f:198) **SEGFAULTS the stripped binary** — so the reference needs an
     INSTRUMENTED REBUILD of fmcrow.f (an unconditional `WRITE` of XV), not the DEBUG keyword.

2. **Additions — `FMCADD` (fire/base/fmcadd.f:65-130) is the litter/wood input each cycle** (NOT fmsadd,
   which is salvage). Per live tree (FMPROB>0, decay class `DKRCLS(SP)`):
   - **Litterfall** = `CROWNW(I,0)·FMPROB(I)/LEAFLF(SP)·P2T` → size-10 litter. `CROWNW(I,0)` = foliage
     biomass (FVSjl has crown_biomass.jl / FMCROWE), `LEAFLF` = per-species leaf lifespan (→ CSV), P2T =
     lb→ton. This is the term that keeps litter from crashing — the crux of the decay/addition coupling.
   - **Woody crown breakage** = `LIMBRK·FMPROB·CROWNW(I,SIZE)·P2T` for SIZE 1..5 → woody CWD; plus
     crown-lift dead material `FMPROB·OLDCRW(I,SIZE)·P2T`.
   - **Snag debris falldown**: the year-1 pool `CWD2B(DKCL,·,1)` flows into CWD (couples to snag.jl
     falldown). Needs the `CWD2B` debris-in-waiting accumulator (new state).
   Dependencies to add: `LEAFLF` table, `LIMBRK` constant, the `CWD2B` pool, and the per-tree foliage +
   woody crown biomass (CROWNW) wired from crown_biomass.jl. Must land WITH chunk 1 to validate.
3. **Per-cycle FFE driver** — call the fuel update (`fmcba!` live/cover + `fmcwd!` decay + `fmsadd!` adds)
   every cycle for FFE-active stands from `grow_cycle!` (today only `fmburn!` runs, only on a fire year).
   Order vs growth/mortality must follow the FFE main (`fmmain.f`). ⚠ Fire-path regression gate is now in
   place (`test_fire.jl`) — fire_early/snt01 stand-4 post-fire `.sum` must stay within its residual.
4. **`stand_carbon_report` per-cycle emission + the `.out` report WRITER** — byte-exact like
   `write_structure_report` (SSTAGE). Only after 1-3; the inventory row is already bit-exact.

## ⚠ Mandatory regression gate (the lesson)
Wiring the fuel update into `grow_cycle!` changes `fire.cwd`/`fire.flive` at **fire time**, which
`fmburn!` reads — so it can move the validated fire results. **snt01 is ecounit 231Dd**, so the
already-committed FULIV2 override ALSO changes its `flive`; neither FULIV2 nor this chunk is exercised
by the current cycle-0/1 `test_snt01.jl` (the fire fires late). Before/after this chunk, run the full
snt01 stand-4 (and `fire_early`) to a fire cycle and diff the post-fire `.sum` (TPA/BA/mortality) vs
the Fortran baseline — the fire path must stay within its existing residual. Add that as a committed
fire-stand test (it is currently a coverage gap independent of this chunk).

## Validation target
The committed `carbon_jenkins.report.save` already has the 1995/2000 rows (DDW 2.5/…, Floor 6.6/…,
Shb 1.0/…, Total 129.9/…). Each chunk above is validated by another report column/cycle reconciling
bit-exact — same method as the inventory cycle. See [[fvsjl-ground-truth-binary-limits]] (the report
is a `.out` text report the stripped binary prints, so this is bit-exact-validatable).

## Status (after the chunk-1-3 + mortality-snag session)
✅ **FMCWD decay** (fuel_decay.jl) · ✅ **FMCADD litterfall + woody breakage** (fuel_additions.jl) ·
✅ **per-cycle driver `ffe_fuel_update!`** · ✅ **crown_biomass V2T/2000 fix** (validated vs a Fortran
XV dump) · ✅ **periodic-mortality → snags** (mortality.jl → add_snag!). The grown-cycle Stand Carbon
Report now reconciles vs the 4-cycle Fortran baseline (carbon_jenkins) for: **live pools every cycle,
Forest Floor every cycle, DDW at 1990/1995/2005, Stand-Dead populated** — see test_carbon.jl.

### Remaining (each a specific term, in priority order)
1. **CWD2B crown-debris scheduling (`FMSCRO`, fmsadd.f:306).** When a tree dies, only its BOLE should
   become the slow-falling snag; its CROWN components go to the `CWD2B` debris-in-waiting pool and fall
   over `TFALL` years (fast for foliage/fine, slow for coarse), flowed to DDW by FMCADD (fmcadd.f:113-125).
   FVSjl currently puts the whole-tree aboveground in the snag → DDW falls too slowly at the first
   post-mortality cycle (DDW 2000: 2.1 vs 3.8; it catches up by 2005). This split should also tighten
   the Stand-Dead column (bole-only, vs whole-tree now: 6.1 vs 5.2 at 2000). Needs the CWD2B state
   array (decay×size×year-to-fall) + FMSCRO + TFALL (a per-species crown-component fall-time table).
2. **Dead-root pool → the report's Below-Dead column** (FMCBIO RBIO at death + CRDCAY decay; currently 0).
3. **`grow_cycle!` hot-path wiring** of `ffe_fuel_update!` (fuel update BEFORE growth per the crown-timing
   finding; fire uses the evolved pools) — behind the `test_fire.jl` regression gate. ⚠ likely shifts the
   fire-stand fuels from initial→evolved, which may surface the known SIMFIRE fire-mortality residual.
4. **The `.out` Stand-Carbon-Report WRITER** (byte-exact like write_structure_report) — last, once 1-3 land.

## CWD2B crown-debris scheduling (FMSCRO) — the last reconciliation gap, fully specified
Status after this session: the grown-cycle Stand Carbon Report reconciles **7 of 9 columns** bit-exact
(above/merch/below-live, **below-dead**, forest floor, DDW@1990/95/2005, shrub) + Stand-Dead populated.
The two remaining (DDW timing @2000: 2.1 vs 3.8; Stand-Dead exact split: 6.1 vs 5.2 @2000) both come
from CWD2B. Everything needed to port it:

- **State:** `cwd2b[decay 1:4, size 0:5, year 1:TFMAX]`, TFMAX=60 (FMPARM.F77) — debris scheduled to
  fall N years out. Add to FireState (a `zeros(4,6,60)`); reset/persist per stand.
- **At death (FMSCRO, fmsadd.f:306 / fmscro.f:119-170):** the BOLE stays the (slow) snag; each crown
  component `CROWNW(SIZE)` (SIZE 0=foliage..5, from `crown_biomass`, now correctly scaled) is divided
  EQUALLY across years `next..min(TSOFT, TFALL(sp,SIZE))` into `cwd2b[dkcl, SIZE, fallyr]`.
- **TFALL (sn/fmvinit.f:1018-1058, keyed by `tfall_cls` ∈ 1..6, already a CSV column):**
  foliage TFALL(·,0)=1 (3 for redcedar); TFALL(·,1)=TFALL(·,2)= {5,3,2,1,1,1}[cls]; TFALL(·,3)=
  {10,6,5,4,3,2}[cls]; TFALL(·,4)=TFALL(·,5)= {25,12,10,8,6,4}[cls]. (cls6=pines fall fastest.)
- **Annual flow (FMCADD, fmcadd.f:113-125):** each year add `cwd2b[·,SIZE,1]` to the matching
  `fire.cwd` down-wood size class, then shift `cwd2b[·,·,k] = cwd2b[·,·,k+1]`. This is the fast
  initial falldown the current whole-tree snag lacks → fixes DDW@2000.
- **Snag = bole only:** with the crown moved to CWD2B, `standing_dead_carbon` should count the bole,
  not the whole-tree Jenkins aboveground (currently over-counts by the crown ⇒ Stand-Dead 6.1 vs 5.2).

After CWD2B: (2) grow_cycle! hot-path wiring (fuel update BEFORE growth, behind test_fire.jl); (3) the
`.out` Stand-Carbon-Report writer. These finish the FFE Stand Carbon Report extension.

### CWD2B — final spec details (size mapping, flow, extra dependency)
- **Size→cwd mapping (FMCADD flow, fmcadd.f:122-135):** CWD2B size 0 (foliage) → `cwd[10,2,dkcl]`
  (litter); CWD2B size 1-5 (woody) → `cwd[size,2,dkcl]`. Each year `DOWN = cwd2b[dkcl,size,1]`; add
  `DOWN·P2T` to the mapped cwd pool; zero it; then shift `cwd2b[·,·,k]=cwd2b[·,·,k+1]`. cwd2b is in
  CROWNW pounds (P2T converts on the way out).
- **Spread (FMSCRO):** each crown component `CROWNW(size)·density` spread EQUALLY over years
  `YNEXTY..min(TSOFT, TFALL(sp,size))`. ⚠ **extra dependency: `TSOFT` = FMSNGDK** (snag soft-transition
  time by species+dbh) — a routine that must also be ported (or approximated as ≥TFALL so TFALL bounds).
- **Bole-split:** at death `bole = jenkins_above − (foliage+Σwoody)·P2T`; store per snag cohort (add an
  `abv` field to SnagList — `add_snag!` has no cohort averaging, so it's a clean push). `update_snags!`
  falldown and `standing_dead_carbon` then use `sn.abv` (bole) instead of whole-tree Jenkins. Both the
  fire-kill (`fmburn!`) and mortality call sites pass the bole; SAFE for test_fire (it asserts the live
  `.sum`, not snag carbon).
⇒ ~10 touch points + the FMSNGDK dependency — a focused chunk, not a tail-end task.

### CWD2B attempt findings (the snag is DYNAMIC, not a static bole)
Attempted the full CWD2B port (cwd2b[4,6,60] state, FMSCRO spread, TFALL/TSOFT, the annual flow, and
a SnagList `abv` bole field with both add_snag! call sites) — then REVERTED it because it did not
reconcile and would have regressed the validated Stand-Dead. What it surfaced:
- **The snag biomass is DYNAMIC.** A static `abv = jenkins_above − all_crown` over-removes the crown:
  Stand-Dead went 6.1 → 4.6 (Fortran 5.2). The faithful model is the snag = whole tree MINUS the crown
  that has FALLEN so far; the standing biomass decreases each cycle as the coarse crown (long TFALL)
  comes off. So `standing_dead_carbon` must subtract the *cumulative* CWD2B outflow for each cohort,
  not a fixed bole — i.e. track per-cohort crown-remaining, or recompute from the cohort's CWD2B share.
- **Timing / crown snapshot.** The crown scheduled at a cycle's mortality must FLOW during that same
  cycle to land in that cycle's DDW; but litterfall needs the cycle-START crown. With the driver placed
  after mortality, the crown flows a cycle late (DDW@2000 stayed 2.1). Needs the fmmain.f:264 crown
  snapshot: record each tree's crown at cycle END for the NEXT cycle's litterfall, and run mortality →
  FMSCRO → the annual flow within the same cycle. The reverted attempt placed the flow in
  `ffe_fuel_update!` ahead of `grow_cycle!`, so the order was wrong.
⇒ CWD2B is more than the spec above: it needs the dynamic-snag accounting + the crown-snapshot
ordering. A focused port should do those two first, validating Stand-Dead AND DDW together.

### ROOT CAUSE of the Stand-Dead discrepancy (a different biomass model)
Checked how Fortran computes the Stand-Dead column: `BIOSNAG = TOTSNG(1)+TOTSNG(2)` (fmdout.f:280),
where TOTSNG = the standing snag **stem VOLUME × V2T** (SNVIS+SNVIH, fmdout.f:153) **PLUS the CWD2B
crown debris still in waiting** (fmdout.f:173: `+= P2T·(CWD2B+CWD2B2)`). So:
- The standing snag biomass is the **stem-volume bole** (a snag stem-volume model SNVIS/SNVIH × V2T),
  NOT the whole-tree Jenkins aboveground FVSjl uses → that alone is why FVSjl reads 6.1 vs 5.2 even
  with no crown split (Jenkins aboveground > merch stem volume).
- The crown is **still counted in Stand-Dead** (via CWD2B) until it falls; when it falls it moves to
  DDW. No double-count — Stand-Dead = bole + (crown not-yet-fallen); DDW = crown fallen + bole fallen.
⇒ Reconciling Stand-Dead bit-exact needs (a) the snag **stem-volume** model (SNVIS/SNVIH), not Jenkins,
and (b) Stand-Dead = stem-bole + Σ(CWD2B)·P2T. This is a deeper port than "subtract the crown"; it
replaces FVSjl's Jenkins-based `standing_dead_carbon` with the volume-based FFE snag model + CWD2B.
The 7-of-9 reconciled columns (incl. below-dead, floor, DDW@95/2005) are independent of this and stand.

### Stand-Dead needs the snag HEIGHT-LOSS + stem-volume model (final root cause)
`TOTSNG = (SNVIS+SNVIH)·V2T` where `SNVIH = FMSVOL(I, HTIH(I))·DENIH` (fmdout.f:140-155): the snag
biomass is the stem VOLUME computed at the snag's CURRENT height `HTIH/HTIS`, and that height
DECREASES each cycle as the snag loses its top (HTX height-loss rate). So Stand-Dead is dynamic on
two axes — falling density (already in `update_snags!`) AND shrinking height/volume (NOT in FVSjl).
To reconcile bit-exact, the FFE snag model must additionally:
1. track per-cohort snag HEIGHT (HTIH/HTIS) and apply the per-species HTX height-loss each cycle;
2. compute the standing bole as `FMSVOL(dbh, current_height)·V2T·density` (a stem-volume routine, not
   Jenkins) — FVSjl's per-tree cubic volume is the starting point but it's a full-height volume;
3. add the CWD2B crown-still-in-waiting to Stand-Dead (and move it to DDW as it falls).
⇒ The Stand-Dead + DDW-timing remainder is a focused FFE snag-volume sub-port (height-loss + FMSVOL +
CWD2B), not a line fix. The 7 reconciled columns (live, below-dead, floor, DDW@95/2005, shrub) stand.

### Stand-Dead validation targets (from an instrumented Fortran dump — fmdout.f SNGBOLE/SNGTOT)
Dumped TOTSNG split into the snag BOLE (stem-volume part) vs the BOLE+CWD2B total for carbon_jenkins
(biomass tons/ac; carbon metric t/ha = value·0.5·2.2417):
| year | bole (t/ac) | total (t/ac) | crown=total−bole | bole C (mt/ha) | total C (mt/ha) |
|------|-------------|--------------|------------------|----------------|------------------|
| 1995 | 0.031       | 0.044        | 0.013            | 0.03           | 0.05             |
| 2000 | 3.315       | 4.619        | 1.304            | **3.72**       | **5.18 ≈ 5.2** ✓ |
| 2005 | 2.925       | 3.989        | 1.064            | 3.28           | 4.47 ≈ 4.5 ✓     |
⇒ The faithful snag BOLE (3.72 mt/ha @2000) is FAR below FVSjl's whole-tree Jenkins snag (6.1) — it is
a small STEM VOLUME × V2T (with height loss), confirming the Jenkins basis is wrong. The CWD2B crown
(1.46 mt/ha @2000) is the separate piece. A future port now has independent targets: validate the
snag-volume bole vs 3.72/3.28 and the CWD2B crown vs 1.46/1.19, then the sum vs the 5.2/4.5 report.

### DDW residual — crown-lift RULED OUT (instrumented Fortran)
Instrumented fmcadd.f to dump the per-year DDW additions: the **crown-lift** (FMPROB·OLDCRW·P2T) and
the live-crown **breakage** (LIMBRK·CROWNW) are both NEGLIGIBLE — ~0.0007 and ~0.0005 t/ac/yr — not
the ~1.7 t/ac gap. So my crown-lift hypothesis was WRONG (good that the dump caught it). The DDW
addition must come from the snag-related path: the CWD2B woody-crown FALLDOWN to down wood (fmcadd.f:
122-135) and/or the snag-bole falldown timing. The next step is to fix the CFALL (cwd2b→DDW) dump
(the instrumentation had a match glitch) and the snag-falldown dump, to pin which one carries the ~1.7
— same instrument-then-port method that landed bole/crown/Stand-Dead. NB the Stand-Dead crown
reconciles (un-fallen cwd2b = 1.46), so the DDW gap is specifically the FALLEN crown/bole timing, not
the scheduled amount.

### DDW residual — structure understood, source still to pin (next: per-year per-source dump)
Two more things ruled out this turn (consulting the Fortran):
- **CWD structure**: Fortran `CWD(3,size,2,5)` is `(I=1:3 category, size, soft/hard, L=1:5)`. The DDW
  report sums `CWD(3,·,·,5)` = aggregate over I=1 (natural) + I=2 (fuel-treatment PILES) and decay
  L=1:4 (fmdout.f:101-109), then only **sizes 1-6** (`SMALL2+LARGE2`, fmdout.f:119-120), NOT 1-9. For
  carbon_jenkins (no piling, no >20" wood) this collapses to FVSjl's `cwd[1:9,:,1:4]` slice, so the
  I/L structure is NOT the gap here. (But FVSjl's `down_wood_carbon` should sum sizes 1-6, not 1-9, to
  be exact for stands with large wood — a minor correctness fix.)
- **crown-lift / breakage**: negligible (~0.0007 t/ac/yr, instrumented) — not the ~1.7 gap.
So the DDW gap is a per-size ADDITION-CONTENT difference in sizes 1-6, timing-related (Fortran adds at
2000, FVSjl at 2005). Next concrete step: instrument fmcadd.f to dump the per-YEAR per-source down-wood
additions (snag-bole falldown vs cwd2b-crown flow vs the FMSADD snag-init), which directly shows which
path carries the ~1.7 and when — the same dump method that landed bole/crown/Stand-Dead. The 8 bit-exact
columns stand; DDW is the only remaining one, now structurally understood.

### DDW residual — ROOT CAUSE: within-cycle deaths-timing (per-year Fortran dump)
Instrumented fmmain.f to dump the down-wood (sizes 1-6) per YEAR. The Fortran DDW grows GRADUALLY
within each cycle and jumps at the cycle boundary, e.g. over 1995→2000: 2.56, 2.82, 3.04, 3.23, 3.40
(year-by-year), then a boundary jump. The carbon-report DDW (3.8 mt/ha @2000 = 3.39 t/ac) matches the
year-1999 value (3.40), i.e. it's read at the cycle boundary BEFORE the next mortality booking.
⇒ The crown/snag debris flows into DDW GRADUALLY over the cycle (annual loop), which requires the
deaths to be spread across the cycle (FVS applies the cycle's mortality + schedules the debris before
the annual fuel loop). FVSjl books all the cycle's mortality at the boundary and flows the crown in the
NEXT cycle (UPD-before-GROW), so its DDW addition lands a cycle late (2005 not 2000) — even though the
TOTAL snags (Stand-Dead) reconcile. So DDW bit-exact needs the within-cycle flow ordering (mortality
schedules → annual flow same cycle) AND the Stand-Dead crown to then be the un-fallen REMAINDER (my
UPD-order full-crown = the Fortran's remainder was a coincidence at 1.46). This is the FFE main-loop
ordering (fmmain.f: GRINCR → FMSADD schedule → annual FMSNAG/FMCWD/FMCADD), i.e. the grow_cycle!
hot-path wiring with the crown snapshot — the integration step, now with the exact behavior pinned.

### DDW residual — ordering experiment DISPROVES the simple fix; it's crown MAGNITUDE
Ran both fuel/grow orderings on carbon_jenkins (experiment in tmp/ddw_exp.jl, NO production change):
- `flow_grow` (current/UPD-1st): DDW 2.06/6.10, Stand-Dead 5.18/4.45 (✓ vs 5.2/4.5), Floor ✓.
- `grow_flow` (GROW-1st):        DDW 2.86/4.82, Stand-Dead 3.72/3.26 (✗ — crown flowed out, bole only).
So grow_flow moves the crown Stand-Dead→DDW but (a) DDW is STILL short (2.86 vs 3.8) and (b) it breaks
Stand-Dead. Conservation: Fortran total mortality-dead @2000 ≈ Stand-Dead 5.2 + (DDW−FUINI ~2.1) +
Below-Dead 1.3 ≈ 8.6; FVSjl flow_grow ≈ 5.18 + ~0.36 + 1.3 ≈ 6.8 — short by ~1.8, ALL in the down-wood
the crown should contribute. ⇒ **The flow_grow Stand-Dead match (crown 1.46) was COINCIDENTAL** — my
full un-flowed crown equals the Fortran's POST-FLOW remaining. The Fortran's FULL scheduled crown is
larger (≈1.46 remaining + ≈2.1 that flows to DDW ≈ 3.5); mine is ≈1.46 and flowing ALL of it only adds
~0.8 to DDW. So the residual is a crown-MAGNITUDE shortfall (FMSCRO/crown_biomass), not just main-loop
ordering. NEXT: dump the Fortran's SCHEDULED crown at FMSCRO/FMSADD (before any annual flow) — if it's
~3.5 not ~1.46, fix crown_biomass magnitude, THEN apply grow_flow ordering with the fine→DDW /
coarse→remaining TFALL split so BOTH DDW and Stand-Dead reconcile. (This supersedes "ordering is the
fix"; the per-year dump showed WHEN, this experiment shows the amount is also short.)

### DDW residual — crown MAGNITUDE confirmed CORRECT; residual is flow-timing/report-point
Instrumented fmscro.f (the crown-scheduling routine) to dump the SCHEDULED crown per death year:
- death 1994: 0.013 · death 1999: **1.304** · death 2004: 1.064 · death 2009: 0.774 t/ac.
My `snag_crown_carbon` gives 1.3 t/ac for the 1995-2000 (death-1999) cohort ⇒ **crown magnitude is
bit-correct** (1.3 vs 1.304). Also confirmed `SCHT == CRWO` ⇒ the OLDCRW term (YRSCYC·OLDCRW·X,
fmscro.f:147) is ZERO here — nothing missing. So last turn's "crown magnitude shortfall" is ALSO ruled
out. The deaths are dated cycle_end−1 (1994/1999/2004/2009); the crown flows starting YNEXTY (~the
boundary year), and the per-year DDW dump's report value (3.39 t/ac @2000) equals the year-1999 value
(3.40), i.e. the report reads BEFORE the boundary-year crown flow. ⇒ The entire DDW residual is a
precise FLOW-TIMING + report-read-point detail in the FFE main loop (YNEXTY flow start vs the report
snapshot), with every magnitude (bole 3.72, crown 1.304, root/BIOROOT) already bit-correct. This is the
grow_cycle! integration wiring — the model pieces are validated; only the per-year flow/report
sequencing around the cycle boundary remains, which is exactly the hot-path wiring step (do it behind
test_fire.jl with the YNEXTY-aligned annual flow, not an end-of-session rush).
Ruled out across this investigation: crown-lift (negligible) · CWD I/L structure (collapses for this
stand) · simple grow/flow ordering swap (breaks Stand-Dead) · crown magnitude (confirmed correct).

### DDW residual — SOLVED (source identified + validated): the crown-LIFT term
The full decomposition, via an instrumented FMCADD (fmmain.f 3-point DDW dump + fmcadd.f per-source dump):
the within-cycle DDW additions (sizes 1-6) are **breakage 0.150 + crown-lift 0.392 + cwd2b-flow ~0 =
0.542 t/ac/yr**, against decay ~−0.30/yr ⇒ net +0.24/yr — exactly the Fortran's within-cycle DDW growth
(2.56→3.40 over 1995-99) that FVSjl was missing. So the entire post-mortality DDW residual is the
**crown-LIFT** term, which FVSjl had deferred with the (now-disproven) note "small for a closing canopy".

Crown-lift = `X · CROWNW(SIZE) · TPA · P2T` per year, where (fmsdit.f:103-117):
  OLDBOT = OLDHT − OLDCRL ; NEWBOT = HT − HT·ICR/100 ; X = (NEWBOT−OLDBOT)/OLDCRL/CYCLEN  (if >0 else 0)
i.e. X = the annual fraction of the OLD crown lifted into dead wood as the crown base rises. OLDCRW is
set to CROWNW at the cycle boundary (fmoldc.f:55) then scaled to X·CROWNW by FMSDIT.

IMPLEMENTATION BLOCKER (why it stays a focused task, not a one-liner): it needs the PREVIOUS-cycle
per-tree HT + crown length. A naive index snapshot FAILS here because the tree list changes every cycle
— carbon_jenkins regen grows it 6→18, and mortality compacts it — so OLD and NEW trees don't line up by
index (and crown_pct is 0 at the first snapshot point). FVS solves this with OLDCRW record-maintenance
in FMTDEL/FMTRIP/FMCMPR (the array is permuted in lock-step with every tree-list mutation). FVSjl needs
either a stable per-tree record id to match across cycles, or the same lock-step maintenance on
fs.oldht/oldcrl through regen/mortality/tripling. THEN: snapshot at cycle end (with valid crown_pct),
add `X·CROWNW·TPA·P2T` alongside breakage in fmcadd_woody!. Magnitude target: +0.39 t/ac/yr (validated).
This supersedes all earlier DDW hypotheses (crown-lift-negligible [buggy dump], ordering, magnitude) —
the source and formula are now PINNED and the only open part is the tree-record plumbing.

### Crown-lift rate X — IMPLEMENTED + tested (the upstream piece); plumbing remains
Landed `crown_lift_rate(oldht, oldcrl, ht, crown_pct, cyclen)` (fuel_additions.jl) — the pure FMSDIT
formula X = (NEWBOT−OLDBOT)/OLDCRL/CYCLEN with the fmsdit.f:106 guards — plus test_crown_lift.jl (7
assertions: worked example X=0.075, the no-rise/degenerate-OLDCRL → 0 guards, cyclen inverse-scaling,
and the OLDCRL-couples-OLDBOT subtlety). So the crown-lift SEMANTIC is now locked and validated
independently (principle 1+2: most-upstream least-dependent piece first, tested to the semantics). The
ONLY remaining part of the crown-lift (hence the last DDW column) is the tree-record PLUMBING: feed it
the previous-cycle per-tree oldht/oldcrl, maintained across the regen/mortality/tripling tree-list
mutations (compact_live!/tredel_compact!/establishment/tripling/sprout — the FVSjl analogues of FVS's
FMTDEL/FMTRIP/FMCMPR OLDCRW maintenance). Suite 4273.

### Crown-lift plumbing — PROVEN to work; remaining issue is POST-GROW placement (grow_cycle! wiring)
Implemented the full plumbing end-to-end (then reverted, see below): added `prev_height`/`prev_crown_len`
as TreeList fields in `_TREE_VEC_FIELDS` (so they ride through compaction/tripling via `copy_tree!` for
free — the clean architecture, no parallel-array maintenance needed), an `fmoldc!` cycle-end snapshot,
and the `(_FM_LIMBRK + X)·CROWNW·TPA·P2T` lift in `fmcadd_woody!`. Result on carbon_jenkins (flow_grow):
the lift FIRES with the right magnitude — DDW@2005 went 6.10 → 9.15. But it is MISTIMED:
- @2000 unchanged (2.06): its lift would come from the 1990 snapshot, but inventory crown_pct=0 there;
- @2005 OVERSHOOTS (9.15 vs 8.0): both cycles' lift piled onto 2005.
ROOT: in the flow-before-grow order the fuel update sees PRE-grow trees, so X measures the PREVIOUS
cycle's crown-base rise; FVS's FMSDIT runs POST-grow (after GRINCR), so X measures the CURRENT cycle's
rise and is applied in that cycle's annual loop. So the lift belongs AFTER grow_cycle!, using the
post-grow crown — exactly the grow_cycle! hot-path wiring. The cwd2b crown flow must STAY pre-grow
(it matched Stand-Dead), so the FFE update splits around grow_cycle!: pre-grow {cwd2b-fall, decay,
litterfall, breakage} and post-grow {crown-lift with the current-cycle rise}. Magnitude is validated
(~0.39 t/ac/yr); only this split placement remains. Reverted the wiring (it overshoots the DDW bound
test pre-split); kept `crown_lift_rate` (committed, tested). The fields/snapshot/lift are a ~30-line
re-add once the grow_cycle! split lands.
