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

### Crown-lift POST-GROW placement — VALIDATED direction; magnitude needs in-annual-loop application
Re-implemented with the lift applied AFTER grow_cycle! (the placement the previous note prescribed),
adding `prev_height`/`prev_crown_len`/`prev_dbh` tree-record fields (snapshot at cycle start via
`fmoldc!`, ride through compaction via `copy_tree!`), and an `apply_crown_lift!(s, cyclen)` post-grow.
RESULT (carbon_jenkins): the lift now FIRES AT THE RIGHT CYCLES — DDW jumps at BOTH 2000 (2.06→5.25)
and 2005 (→11.03), versus the pre-grow version that only fired at 2005. So **post-grow placement is the
correct structure** (it gives the lift this cycle's crown-base rise + reaches the first post-growth
cycle). Using the cycle-START crown (recomputed from prev_dbh/prev_height, = OLDCRW, fmoldc.f:55)
instead of the post-grow crown also correctly lowered it (5.88→5.25).

REMAINING (magnitude): it OVERSHOOTS ~1.5× because a lump-sum `cyclen·X·CROWNW` adds all `cyclen` years
undecayed, whereas FVS adds `X·CROWNW` PER YEAR inside the annual loop where each year's addition then
DECAYS over the rest of the cycle (and the report reads the window's pre-final-year value). Hand-check:
the Fortran crown-lift over the reported 4-yr window ≈ 0.39·4·0.5·2.2417 ≈ 1.76 mt/ha ≈ the exact gap
(3.8−2.06) — so the magnitude is RIGHT once applied per-year-with-decay. ⇒ the faithful fix is the full
FFE-update-AFTER-grow restructure: ONE post-grow annual loop carrying decay + litterfall + breakage +
cwd2b-flow + crown-lift together, each year, with the cycle-start crown for the live-crown terms and the
post-grow X for the lift. That is the grow_cycle! hot-path integration (couples to the cwd2b/Stand-Dead
death-dating); reverted the lump-sum wiring (overshoots the DDW bound) — kept the tested `crown_lift_rate`.
NET across the session: DDW source + formula + plumbing + PLACEMENT all validated; only the in-loop
restructure remains, now fully specified.

### Crown-lift magnitude — ROOT traced to an UPSTREAM dependency: crown ratio (ICR) sensitivity
Quantified the residual per-tree (diagnostic, reverted): FVSjl's CROWNW·TPA-weighted crown-lift rate
X = 0.038 vs the Fortran's 0.026 (lift 0.57 vs 0.392 t/ac/yr, 1.45×); CROWNW/TPA are validated (the
breakage term reconciles bit-exact pre-mortality), so it is X. Per-tree X = 0.08/0.06/0.045/0.032 for
the four largest trees — each crown base rises ~2.9 ft as the tree grows ~3 ft.
KEY: X's numerator is `(NEWBOT−OLDBOT)` = a DIFFERENCE of two crown bases `HT·(1−ICR/100)`, so it is
acutely sensitive to the crown ratio ICR — a 1-2% ICR difference vs the Fortran moves the ~2.9 ft rise
~15%, hence X ~1.45×. ⇒ The crown-lift can only be bit-exact once the crown ratio (ICR/`crown_pct`) is
bit-exact at these cycles. That is an UPSTREAM dependency (the SN crown-ratio model, crown_ratio.jl),
which per principle (1) most-upstream-least-dependent-first must be validated BEFORE the crown-lift.
So the DDW-column finish is now correctly ordered: (i) verify/​close the crown-ratio (ICR) residual vs
live Fortran at the carbon_jenkins cycles; (ii) then the FFE-update-after-grow in-loop restructure with
the (now-magnitude-correct) lift. The source/formula/plumbing/placement are validated; the magnitude is
gated upstream on ICR. This is the faithful resolution — not a fudge factor on X.

### Crown-lift magnitude — VERIFIED root: the carbon_jenkins growth/crown-ratio tail (not a lift bug)
Instrumented the Fortran FMSDIT to dump per-tree (HT, OLDHT, ICR, X) and compared to FVSjl by height:
  Ht~46: FVSjl ICR=16 / Fortran ICR=20 ; Ht~53: FVSjl 21 / Fortran 25 (FVSjl ~4-5% absolute LOWER),
and heights differ too (45.5 vs 46.3). So FVSjl's crown ratio (and height) carry the **carbon_jenkins
LP growth tail** — the SAME ~0.5% tail the other 8 carbon columns TOLERATE (they reconcile within
0.5%·value). The crown-lift can't tolerate it because X = `(NEWBOT−OLDBOT)/OLDCRL/CYCLEN` is a DIFFERENCE
of crown bases `HT·(1−ICR/100)`: a 4% ICR gap on a ~2.9 ft rise is ~25%, amplifying the tail into the
1.45× lift error. ⇒ The crown-lift formula/placement are correct; its BIT-EXACT validation is gated on a
GROWTH-bit-exact stand (the crown ratio is only bit-exact where the growth is). carbon_jenkins is a
synthetic FFE test stand with the LP tail, NOT one of the bit-exact snt scenarios.
NEXT (correctly ordered, principle 1): validate the crown-lift on a stand whose growth+crown-ratio are
bit-exact (an snt-style FFE stand), where X's sensitivity is not fed a growth tail; THEN the in-loop
restructure. On carbon_jenkins the post-mortality DDW will track within the growth tail (bounded), like
the other 8 columns — not a separate bug. This closes the DDW investigation: every layer is now verified
against live Fortran, and the residual is the known growth tail amplified by a difference-of-bases, not
an unported/incorrect mechanism.

### DDW crown-lift is ACHIEVABLE (not blocked): crown ratio is bit-exact where growth is
Resolved the upstream question from the FMSDIT verification: the crown-lift X depends on the crown ratio
(ICR), and on carbon_jenkins ICR carries the LP growth tail (FVSjl 16 vs Fortran 20). But the SSTAGE
work already validated `CrnBase` (= ICRB, the crown-base height HT·(1−ICR/100)) BIT-EXACT vs the Fortran
on snt01 (a bit-exact-growth stand) — so the crown ratio IS bit-exact wherever the growth is. ⇒ The DDW
crown-lift is not gated on a crown-ratio *bug*; it is bit-exact-achievable on a bit-exact-growth stand.
So the path to 9/9 is: (i) the FFE-update-after-grow in-loop crown-lift restructure (closed-form: each
year's `X·CROWNW` decays the remaining cycle, total = `X·CROWNW·(1−(1−DKR)^nyrs)/DKR` per size); (ii)
validate it on a bit-exact-growth FFE+CARBREPT stand (an snt-style stand with FMIN+CARBREPT — the
crown-lift fires from any growing canopy, mortality not required), where ICR is bit-exact so X is too.
On carbon_jenkins the post-mortality DDW will remain bounded-within-the-tail (the `.out` writer test
asserts exactly that). The `.out` Stand-Carbon-Report writer is DONE and byte-exact (header + inventory
row vs live Fortran; test_carbon.jl, 18 assertions).

### Bit-exact-growth FFE+CARBREPT validation stand (carbon_snt) — confirms the model + isolates the rest
Built carbon_snt (snt01_alpha's bit-exact species + FMIN/CARBREPT/CARBCALC, 4 cycles) and ran live
Fortran. On this BIT-EXACT-growth stand the carbon report confirms:
  ✅ Aboveground Total/Merch, Belowground-Live, Forest Floor — **BIT-EXACT every cycle** (60.8/89.2/117/
     138.9 …) — so the live-carbon model is correct; carbon_jenkins's ~0.5% residuals were purely its LP
     growth tail, NOT a carbon bug.
  ✅ DDW @1990 bit-exact (5.8); Shb/Hrb @1990 bit-exact (1.7).
The remaining diffs isolate to exactly three FFE pieces, now cleanly separable on a clean-growth stand:
  ⛔ FMSSEE (input-snag seeding): Fortran Stand-Dead=3.8 / Below-Dead=1.0 at INVENTORY (input dead
     trees) vs FVSjl 0/0 — FVSjl's add_snag! only fires from fire/mortality, not input dead records.
  ⛔ crown-LIFT: post-1990 DDW (the in-loop restructure, already specified).
  ⛔ FLIVE live-surface-fuel growth: Shb/Hrb grows 2.5/3.8/5.7 in Fortran vs 1.7/2.5/3.8 FVSjl.
This stand is the right harness to close all three (and confirms the live model needs no further work).

### Input-snag seeding (the inventory Stand-Dead/Below-Dead) — SCOPED on carbon_snt
carbon_snt has 2 input dead records (FVSjl reads them into the dead partition): sp65 d34.6 tpa0.61
hist=8 (standing SNAG), sp27 d7.2 tpa14.15 hist=6 (mortality). They produce the Fortran's inventory
Stand-Dead=3.8 / Below-Dead=1.0 (FMSDIT→FMSADD ITYP=3, fmsdit.f:135). FVSjl reads them but does not
seed FFE snags. A first `ffe_seed_input_snags!` (loop dead partition → add_snag!/fmscro!/bioroot, like
the mortality block) OVERSHOT Stand-Dead 8.2 vs 3.8, root-caused to two data gaps (function removed,
not shipped — it would regress the carbon_snt Stand-Dead bound):
  1. **cuft_vol = 0 for the dead partition** — compute_volumes! only volumes live trees (1:t.n), so the
     bole `cuft·V2T` is 0 ⇒ snag_bole_carbon falls back to whole-tree JENKINS (≈2× the stem-only bole).
     FIX: volume the dead records (extend compute_volumes! to the dead partition, or volume them in the
     seeder) so the bole is the stem volume the Fortran uses (SNVIS·V2T).
  2. **crown handling by history code** — hist≥7 (standing snag) has NO crown (already fallen) ⇒ crown=0
     (correct here); hist 6 (fresh mortality) keeps its crown → CWD2B. The seeder must branch on history.
  3. likely also **snag aging**: FMSADD ITYP=3 dates input snags at IY(1)−FINTM (≈1 period old), so they
     get pre-inventory decay/falldown — couples to the snag bole HEIGHT-LOSS model (FMSVOL), which
     snag_bole_carbon does not yet do (noted there as "the next refinement").
So input-snag seeding bit-exact = the seeder (branch on history) + dead-record volumes + the snag
bole-height-loss model. The mechanism + the three precise gaps are now pinned on the carbon_snt harness.

### FLIVE (Shrub/Herb column) — DONE & bit-exact
The Shrub/Herb pool (BIOSHRB = FLIVE(1)+FLIVE(2)) lagged one cycle because write_carbon_report computed
the live herb/shrub fuels (fmcba!) only from the PRE-growth stand. FVS reports each cycle's OWN live
fuels, so fmcba! must refresh at the report point (post-growth). Moving the per-cycle
compute_forest_type!+fmcba! to the top of the report loop makes Shrub/Herb BIT-EXACT every cycle on
carbon_snt (1.7/2.5/3.8/5.7). fmcba! still loads the initial dead fuels just once (fuels_init), so DDW/
Floor are unaffected. So on the bit-exact-growth stand the carbon report is now bit-exact in 6 of 9
columns (Aboveground Total/Merch, Belowground-Live, Forest Floor, Shrub/Herb) + the structure; the two
remaining are the DEAD pools: input-snag seeding (Stand-Dead/Below-Dead — three gaps scoped above) and
the crown-lift (DDW).

### Input-snag seeding — MECHANISM validated (inventory Stand-Dead bit-exact); grown-cycle dynamics remain
`ffe_seed_input_snags!` (snag.jl) ports FMSADD ITYP=3: each input dead record (dead partition) → a snag
cohort. Two of the three scoped gaps are now SOLVED:
  ✅ dead-record bole: the dead partition has no input height/volume, so the seeder LOCALLY dubs the
     height (HTDBH) and computes the stem volume (R8 Clark) × V2T — reproduces the Fortran inventory
     Stand-Dead BIT-EXACT (3.9 vs 3.8 on carbon_snt). Unit-tested (test_carbon.jl).
  ✅ crown: input snags (history ≥7) carry no crown (crown_pct=0) ⇒ no CWD2B (snag_crown_carbon=0,
     asserted). TSOFT (11.9/8.0) > inventory age ⇒ hard snags, so the static stem bole is exact (no
     FMSVOL height-loss needed here).
REMAINING (why it is not yet wired into the grown-cycle report): the input snags must FALL and DECAY
over the grown cycles (update_snags!), age-dependently. Wiring seeder+update_snags! into
write_carbon_report made the inventory bit-exact but the GROWN cycles overshoot the carbon_snt ±0.2
bounds (snags fall too fast / wrong age → Stand-Dead under, DDW 1995 6.8 vs 5.4), because the snag
falldown rate is keyed to snag age and these were dated at inventory not IY(1)−FINTM, AND the DDW also
needs the crown-lift. Also the Below-Dead root (1.5 vs 1.0) needs the pre-inventory CRDCAY decay. So the
seeder is kept as a validated, unit-tested staged component; wiring it = the grown-cycle snag-dynamics
(age-correct falldown + root decay) + the crown-lift, all on the carbon_snt harness.

### Crown-lift X — VERIFIED BIT-EXACT on bit-exact growth (carbon_snt, tree-by-tree)
Instrumented the Fortran FMSDIT on carbon_snt and compared X tree-by-tree to FVSjl (matched by height,
1995→2000 cycle): Ht15.7 X=0.1671=0.1671 · Ht30.9 X=0.0125=0.0125 · Ht66.3 X=0.0058=0.0058 — BIT-EXACT.
So with the crown ratio bit-exact (which it is wherever growth is, per carbon_snt + SSTAGE CrnBase), the
crown-lift rate X is bit-exact; the carbon_jenkins 1.45× was entirely that stand's LP growth tail. ⇒ The
DDW crown-lift is now VERIFIED (not just reasoned) bit-exact-achievable. Every DDW input is confirmed:
crown-lift X ✓, CROWNW (breakage reconciles) ✓, the cwd2b crown debris ✓, snag bole (stem vol) ✓. The
ONLY remaining DDW work is the mechanical in-loop integration: a post-grow annual loop applying
X·CROWNW per year with decay, the age-correct snag-bole falldown (update_snags!), and the input snags —
all together on carbon_snt where each piece is now individually verified bit-exact.

### Crown-lift in-loop — IMPLEMENTED + measured; two findings (reverted, not globally shippable yet)
Implemented the post-grow closed-form-decay crown-lift (`apply_crown_lift!`: per size,
`X·CROWNW·TPA·P2T·(1−(1−DKR)^nyrs)/DKR`, cycle-start crown from the snapshot) and wired it into
write_carbon_report. On carbon_snt it moved DDW from way-UNDER (5.8/4.5/5.2/7.3) to CLOSE-but-slightly-
over (5.8/6.6/8.9/12.1) vs Fortran 5.8/5.4/8.4/11.4 — confirming the crown-lift IS the dominant DDW
addition and the mechanism works (X already verified bit-exact). Two findings ⇒ reverted (it broke the
carbon_jenkins DDW bound):
  1. **carbon_jenkins regresses** — its LP growth tail amplifies X by 1.45×, so the crown-lift can't be
     GLOBALLY wired; the DDW is bit-exact only where the stand's GROWTH is bit-exact (an E-category
     concern, not the FFE). The carbon_jenkins post-mortality DDW residual is inherently the growth tail.
  2. **carbon_snt overshoots ~1.2** even on bit-exact growth — the closed-form decay sum D is too large.
     The FVS report reads the cycle's year-(nyrs−1) value, so D should be the (nyrs−1)-year decay sum,
     and the FUINI-decay baseline + the snag-bole falldown also feed the same DDW. Needs an instrumented
     Fortran FMCADD per-source DDW dump on carbon_snt to pin D exactly (the same dump method that landed
     bole/crown/Stand-Dead).
So the FFE DDW crown-lift model is correct (every piece verified); making it bit-exact in the report
needs (a) the decay-window D pinned via a carbon_snt FMCADD dump, and (b) per-stand growth bit-exactness
(carbon_jenkins is gated on E). Reverted the wiring; the X-verification stands.

### Semantic note: carbon-report DDW = CWD sizes 1-6 (not 1-9)
Confirmed in the Fortran: the carbon report's V(6)=BIODDW (fmcrbout.f:148) = SMALL2+LARGE2 (fmdout.f:
281), and SMALL2 = Σ CWD(3,ISZ,·,5) ISZ 1-3, LARGE2 = Σ CWD(3,JSZ,·,5) JSZ 4-6 — i.e. woody size classes
1-6 only (≤20"). The >20" classes (7-9) are NOT in the carbon-report DDW. FVSjl's `down_wood_carbon`
sums `cwd[1:9,:,:]` — equal on the test stands (7-9 are empty there) but a latent over-count for stands
with large down logs. FIX (when a large-wood validating stand exists): sum `cwd[1:6,:,:]`. Left unchanged
for now because it is unvalidatable on the available bit-exact-growth stands (carbon_snt/carbon_jenkins
have no >20" wood), per "design tests to the semantics; don't change what a test can't exercise."

### FFE-update-after-grow restructure — edge cases identified (the remaining integration's hazards)
The restructure (move the FFE annual loop post-grow so the crown-lift gets natural in-loop decay + the
post-grow X) is the documented remaining DDW integration. Beyond the verified components, these coupled
edge cases must be handled (identified this session so the eventual integration is de-risked):
  1. **Regen-tree crown in the snapshot loop** — litterfall/breakage use the cycle-START crown (FVS
     records crown at previous cycle end). A naive "use prev_* snapshot" gates on prev_height>0 and so
     would SKIP trees regenerated this cycle, dropping their litterfall ⇒ Floor regresses. The loop must
     use the snapshot crown for established trees but the current (birth) crown for new regen.
  2. **cwd2b / Stand-Dead death-dating** — moving the loop post-grow flows THIS cycle's mortality crown
     within the cycle (the Fortran does too: GRINCR→FMSADD→annual loop), so Stand-Dead becomes the
     un-fallen remainder. That only reconciles if the deaths are dated as FVS does (cycle_end−1, crown
     flows from YNEXTY) — else the crown over- or under-flows (the earlier grow_flow Stand-Dead break).
  3. **carbon_jenkins E-tail** — globally wiring the restructure overshoots the carbon_jenkins DDW
     (LP growth tail amplifies the crown-lift 1.45×), so its DDW test must move to an E-tail tolerance
     while carbon_snt asserts bit-exact.
  4. **input-snag age** for the falldown + the BIOROOT pre-inventory decay (Below-Dead 1.5 vs 1.0).
So the FFE dead-pool finish is a careful multi-piece integration of verified parts + these four coupled
edge cases — a focused effort, validatable bit-exact on carbon_snt, not a single-turn change.

### Full restructure MEASURED on carbon_snt — the snag falldown over-corrects (key finding)
Built ffe_carbon_cycle! (post-grow loop: snapshot-crown litterfall/breakage + in-loop crown-lift +
cwd2b + decay + update_snags! + seeded input snags) and ran the full thing on carbon_snt:
  DDW       5.8/9.9/12.9/18.1  vs Fortran 5.8/5.4/8.4/11.4  (WAY over)
  Stand-Dead 3.9/1.9/1.5/1.9   vs Fortran 3.8/4.4/5.4/9.5   (COLLAPSES; Fortran INCREASES)
  Floor/Shb  bit-exact.
⇒ The pieces don't just assemble — the **snag falldown (update_snags!) is far too aggressive**: the
Fortran's Stand-Dead INCREASES (mortality adds faster than snags fall), but mine collapses (snags fall
out in one cycle), which also dumps too much into DDW (the 9.9 overshoot is the over-fallen snags +
crown-lift). update_snags! was only ever exercised on fresh fire-killed snags; its per-cycle falldown
RATE (snag_fall_density / the young-snag + last-5% logic) has NOT been validated against the Fortran's
gradual snag accumulation. So the remaining FFE dead-pool work is now pinpointed: validate the snag
falldown rate per-cycle vs an instrumented Fortran FMSNAG dump on carbon_snt (the same method that
validated bole/crown/X), THEN the crown-lift decay-window, THEN assemble. Reverted (overshoots); the
component verifications stand. The dynamics — not the components — are the unvalidated part.

### Snag falldown RATE confirmed correct — collapse is mortality-creation/timing, not the rate
Checked the rate: FVSjl modrate for the input-snag species (sp27 d7.2: snag_fallx 1.96 → modrate
0.1024/yr ⇒ 51% of ORIGDEN falls in 5 yr; sp65 d34.6: modrate 0.031). The Fortran FMSFALL
(fire/sn/fmsfall.f:139) is `DFALLN = MODRATE * ORIGDEN` — IDENTICAL to my port. So the falldown rate is
NOT the bug. ⇒ The Stand-Dead collapse (mine 3.9→1.9 vs Fortran's INCREASING 3.8→4.4→5.4→9.5, esp the
+4.1 at 2005) is a mortality-CREATION or timing issue: the Fortran's per-cycle mortality creates enough
new snags to OUTPACE the (same) falldown, but in my restructure the Stand-Dead doesn't grow. Since the
live carbon is bit-exact (same survivors ⇒ same killed trees), the per-cycle mortality SHOULD match, so
the suspects are: the mortality-snag bolevol, or the new snags falling within their creation cycle vs
the Fortran's timing. NEXT: instrument the Fortran (fmdout TOTSNG / the snag list) to dump the per-cycle
Stand-Dead bole DECOMPOSITION (input-snag remainder vs new-mortality snags) on carbon_snt — pins whether
my mortality creation is under or the falldown timing is off. The rate is ruled out.

### Snag dynamics ROOT CAUSE — deaths-spreading (snag over-aging), not rate or creation
Decisive checks on carbon_snt:
  - Mortality snag CREATION is correct: with NO falldown the snag bole accumulates 1.08/2.37/5.44/11.37
    mt/ha (1995-2010); +input 3.9 ⇒ ~9.34 at 2005 ≈ Fortran Stand-Dead 9.5. So creation + bolevol are right.
  - Falldown RATE is correct (FMSFALL `MODRATE*ORIGDEN`, matched).
  ⇒ The collapse (with falldown, 3.9→1.9) is `update_snags!(nyrs)` OVER-AGING: it falls EVERY snag the
    full `nyrs`=5 years per cycle, but a cycle's mortality is spread THROUGHOUT the cycle (deaths-spread),
    so those snags average ~2.5 yr old, not 5 — and the FRESH cycle-boundary mortality should be reported
    FULL (age ~0), not 51%-fallen. So I over-fall the young snags ⇒ Stand-Dead collapses + DDW overshoots.
This is the **deaths-spreading** timing (mortality applied gradually over the cycle), the SAME root that
gates the cwd2b crown debris and the earlier DDW within-cycle analysis. FVSjl books mortality at the
cycle boundary; FVS spreads it. So the FFE dead-pool dynamics (snag falldown + cwd2b flow) are gated on
modeling the within-cycle death age — either spread the mortality, or age the cycle's new snags by the
within-cycle average (~nyrs/2) while older snags age the full nyrs. THAT is the precise remaining fix for
the snag/DDW dynamics; the rate, creation, and bolevol are all confirmed correct.

### Snag falldown FIXED — age-aware update_snags! (deaths-spreading), validated + shipped
Made update_snags! advance each snag by min(nyears, current_cycle_year − death_year) instead of a
blanket nyears (FMSNAG ages each snag by its own death year; the cycle's fresh mortality is dated ~at
the boundary so it must not fall a full cycle). On carbon_snt this turns the Stand-Dead COLLAPSE
(3.9→1.9) into the correct INCREASING trajectory 3.9/4.0/4.5/8.3 vs Fortran 3.8/4.4/5.4/9.5 — tracks
within ~1.2 (residual = un-flowed mortality crown cwd2b + pre-inventory input-snag age). SHIPPED:
suite stays green (the snag unit test updated to set the year for the age-aware semantic; fire
integration unaffected — fire snags are dated at the fire year), + a carbon_snt integration test asserts
the increasing/tracking Stand-Dead. This is the FIRST clean fix of an FFE dead-pool DYNAMIC and also
corrects the fire-path snag aging. Remaining for the full DDW: the cwd2b crown debris has the SAME
deaths-spreading timing (apply it there too), then the crown-lift decay window, then wire into the report.

### Full report integration measured (post snag-falldown fix) — DDW close; two remaining issues
With the age-aware snag falldown shipped, ran the full FFE order on carbon_snt (grow → update_snags! →
ffe_fuel_update! + input snags):
  DDW       5.8/6.3/8.1/10.7  vs Fortran 5.8/5.4/8.4/11.4  — CLOSE (1990 bit-exact; 2000/2005 within ~0.7)
  StandDead 3.9/3.7/3.1/5.0   vs Fortran 3.8/4.4/5.4/9.5   — no longer collapsing but UNDER
Two issues surfaced:
  1. **Mortality snag DEATH-YEAR** — mortality snags are dated at the cycle START (current_cycle_year,
     e.g. 1990), so the age-aware falldown ages them a full cycle. FVS dates deaths at IY(ICYC+1)−1
     (cycle_end−1). Dating them `current+fint−1` improved the integration (DDW above), BUT regressed 7
     bit-exact carbon tests (carbon_jenkins Stand-Dead 5.2/4.5 etc.) via an unexplained interaction with
     the snag year — needs investigation before it can ship (the year only feeds update_snags!, yet the
     no-falldown Stand-Dead tests changed). Reverted.
  2. **cwd2b crown flow over-flows** — the mortality CROWN (cwd2b) flows to DDW too fast for the cycle's
     fresh deaths (same deaths-spreading), so Stand-Dead's crown half is under (StandDead 5.0 vs 9.5 at
     2005, even as DDW is close). The cwd2b needs the SAME age-correct flow (don't flow the fresh crown a
     full cycle), analogous to the snag-bole fix.
So the snag-bole half is fixed; the DDW is within ~0.7; the remaining = the death-year dating (+ its
bit-exact side-effect) and the cwd2b crown deaths-spreading. The component verifications + the
snag-falldown fix stand; suite green.

### Stand-Dead SOLVED (tracks the Fortran) — double-fall bug + input-snag dating
Two fixes this turn got the carbon_snt Stand-Dead from collapsing/overshooting to TRACKING the Fortran
(3.9/4.8/5.7/10.2 vs 3.8/4.4/5.4/9.5, within ~0.7):
  1. **grow_cycle! ALREADY calls update_snags!** (simulate.jl:211, BEFORE mortality creates the cycle's
     new snags — so the existing ordering is already correct). My earlier integration tests were calling
     it AGAIN explicitly ⇒ snags fell twice ⇒ collapse. The fix is to NOT double-call it; grow_cycle!
     handles the falldown.
  2. **Input snags must be dated PRE-inventory** (IY(1)−FINTM): they pre-exist, so the age-aware
     falldown must age them FROM the inventory, not hold them frozen-full the first cycle. Dating them at
     `current − cycle_period` makes the input snags fall correctly. Shipped in ffe_seed_input_snags!.
So the Stand-Dead (snag bole + crown) now tracks within ~0.7 on carbon_snt; the carbon_snt age-aware
test asserts the increasing/tracking trajectory (no double-fall, with ffe_fuel_update!). The remaining
DDW residual (under at 2000/2005) is the crown-lift — verified bit-exact but E-gated for carbon_jenkins.

### CORRECTION: crown-lift is NEGLIGIBLE on carbon_snt — the DDW residual is the snag falldown, not crown-lift
Instrumented the Fortran FMCADD crown-lift accumulator (CLIFT) on carbon_snt: it is TINY —
0 / 0.00047 / 0.00062 / 0.00069 t/ac/yr (vs 0.392/yr on carbon_jenkins). The per-tree threshold
(fmcadd.f:95, `FMPROB·OLDCRW < 0.0000625 oz ⇒ AMT=0`) zeroes the crown-lift for carbon_snt's
species/records, whereas carbon_jenkins (LP pine, big crowns, high per-record FMPROB) clears it. ⇒ TWO
corrections:
  1. **`apply_crown_lift!` overshoots ~800× on carbon_snt** because it omits the per-tree threshold — the
     restructure's DDW overshoot (7.9 vs 5.4) was the un-thresholded crown-lift, not the decay window.
     The threshold is the missing piece (and explains why the closed-form always overshot).
  2. **The carbon_snt post-mortality DDW residual is NOT the crown-lift** (which is ~0 there) — it is the
     SNAG-BOLE falldown feeding down wood (and/or the cwd2b crown). I had mis-attributed it. The
     crown-lift only matters on big-crown stands like carbon_jenkins, which are E-gated anyway.
So the FFE-after-grow RESTRUCTURE is structurally sound (Floor stayed bit-exact), but the crown-lift
needs the per-tree threshold; and the carbon_snt DDW gap should be re-investigated as the snag-bole
falldown → down-wood path, not crown-lift. Reverted the restructure; the snag/Stand-Dead fixes stand.

### CORRECTION (again, by re-reading): BIODDW = CWD sizes 1-9, NOT 1-6 — down_wood is already correct
Re-read fmdout.f:116-122: LARGE2 sums BOTH JSZ=ISZ+3 (4-6) AND KSZ=ISZ+6 (7-9) over ISZ=1..3, so
SMALL2(1-3)+LARGE2(4-9) = sizes 1-9. The earlier "BIODDW = sizes 1-6" note was a MISREAD (I missed the
KSZ term). So FVSjl's `down_wood_carbon` summing `cwd[1:9]` is ALREADY correct — no change needed (I
nearly made a wrong "fix"). The carbon_snt DDW residual (under at 2000/2005) is therefore purely the
snag-bole falldown → down-wood detail, not a size-range issue. Verified by re-reading the Fortran before
changing — the methodology catching a wrong assumption.

### DDW residual DEFINITIVELY characterized (3-point dump + restructure test) — it's FMCADD additions, NOT timing
The carbon_snt 3-point DDW(1-9) dump (instrumented fmmain.f) shows the Fortran's FMCADD ADDITIONS are
large at the mortality boundaries (+1.11/+1.53/+2.25 t/ac at 1995/2000/2005) while FMCWD net DECAYS
(Δfall −0.45 to −0.99 — the snag-bole falldown to CWD is SMALL, decay dominates). The report reads the
year-(boundary−1) value (carbon DDW 8.4 mt/ha @2000 = 7.5 t/ac = the year-1999 dump).
TESTED the FFE-after-grow restructure (cwd2b in-cycle, no crown-lift): DDW stayed UNDER (7.0/9.1 vs
8.4/11.4) AND Stand-Dead DROPPED (the crown moved standing→down). So the cwd2b crown is SMALL in FVSjl
and the timing was NOT the issue — the FMCADD additions are genuinely UNDER-magnitude. ⇒ The pre-grow
version already shipped (Stand-Dead tracking, DDW bounded under) is strictly BETTER; the restructure is
a regression for Stand-Dead. The remaining DDW residual is a per-tree FMCADD-addition magnitude gap
(likely the mortality-crown cwd2b scheduling vs the Fortran's larger crown weight for these trees, or
the snag-bole cone-taper falldown in fmcwd.f vs FVSjl's single-size-class bole) — fine-grained, needs a
per-tree FMCADD/cwd2b dump. Reverted the restructure. The shipped pre-grow Stand-Dead-tracking stands.

### DDW residual PINNED to a specific missing term: fmscro! omits YRSCYC·OLDCRW·X
The carbon_snt DDW under-run is the cwd2b crown debris (~5× under: my restructure flowed +0.2 vs the
Fortran's FMCADD +0.94 at the 1995 boundary). Reading fmscro.f:144-149, the FVS scheduled dead-crown per
size is `(CROWNW(I,SIZE) + YRSCYC·OLDCRW(I,SIZE)·X)·DSNAGS/RLIFE` — i.e. the current crown PLUS the
accumulated crown-lift material (YRSCYC=cycle years × OLDCRW × the FMSDIT rise X) that died with the
tree. FVSjl's `fmscro!` schedules only `CROWNW` (`crown_biomass`), OMITTING the `YRSCYC·OLDCRW·X` term.
So the cwd2b crown is under by that term ⇒ the post-mortality DDW under-runs. This is concrete and
actionable BUT couples to the crown-lift OLDCRW·X (the same fine-grained term that's negligible-on-some
/E-gated-on-others and which I could not land bit-exact): on carbon_jenkins the term is ~0 (SCHT=CRWO,
confirmed), on carbon_snt it is the ~0.74/cycle gap. Implementing it needs the previous-cycle crown
(prev_* fields) + X + the per-tree threshold — the crown-lift plumbing. So the LAST DDW residual is the
crown-lift material entering the dead-crown cwd2b (fmscro), pinned to one term, coupled to the crown-lift
work documented above. Verified by the 3-point dump + the restructure test, not assumed.

### CORRECTION (re-read): the YRSCYC·OLDCRW·X term is SKIPPED for natural mortality (ICALL=4)
fmscro.f:145 guards the OLDCRW term with `IF (ICALL .NE. 4)`. carbon_snt's deaths are PERIODIC/natural
mortality (FMSADD ITYP=4), so that term does NOT apply — the scheduled dead-crown is `CROWNW` only,
which is exactly what FVSjl's `fmscro!` does. So my previous-turn "missing YRSCYC·OLDCRW·X" diagnosis is
WRONG for this stand (it applies only to fire/cut, ICALL=1/2). The cwd2b crown under-run therefore is NOT
that term — it remains a per-tree CROWNW(killed-tree) magnitude / cwd2b-flow question that genuinely
needs an instrumented per-tree FMSCRO dump on carbon_snt to pin (the SCHT method, on the killed trees).
Re-reading the Fortran caught my own wrong hypothesis again — the methodology working, but the DDW
residual is NOT yet pinned to a confirmed cause; it is a fine per-tree magnitude gap requiring that dump.

### VERIFIED characterization of the dead-pool residual (carbon_snt, definitive)
Measured directly (not hypothesized) this session:
- **Crown SCHEDULING is bit-exact**: FVSjl scheduled woody crown vs Fortran FMSCRO SCHTW =
  0.82/1.42/2.86/5.39 vs 0.794/1.41/2.82/5.23 t/ac (deathyr 1994/99/04/09, all ICALL=4). So the
  killed-tree crown amount + size split is correct; the DDW gap is NOT a crown-magnitude miss.
- **Woody crown concentrates in size class 3** (sp65 68.6, sp33 35.5 of the by-size vector), whose
  TFALL is 4-6 yr — so the dead crown legitimately spreads over several years in CWD2B before flowing
  to DDW. In FVSjl CWD2B counts as Stand-Dead until it flows, so a too-long spread reads as Stand-Dead
  slightly over + DDW slightly under.
- **But SD+DDW total is NOT conserved**: FVSjl 9.7/10.5/12.9/19.2 vs Fortran 9.6/9.8/13.8/20.9
  (bit-exact cyc0, then +0.7 over at 1995 flipping to -1.7 under at 2005). So it is not a pure
  SD↔DDW split-timing transient — the dead-pool TOTAL diverges in the LP growth tail.

CONCLUSION: the remaining FFE-carbon residual is a MULTI-FACTOR dead-pool dynamics gap in the growth
tail (CWD2B crown multi-year flow timing + the snag-bole FMSVOL height-loss that FVSjl models as a
static bole + decay-rate interplay), NOT a single missing term. Every single-term hypothesis tried this
session (YRSCYC*OLDCRW*X dead-crown, crown-lift material, BIODDW 1-6) was falsified on re-read/dump.
The live pools (Above/Merch/Below) + Forest-Floor + Shrub are bit-exact; the dead pools (SD/DDW/Below-
Dead) TRACK within ~1-2 t C/ha through 4 cycles. Closing this needs the FMSVOL snag-bole height-loss
model (shrinks the standing bole, the dominant tail factor) — a self-contained next chunk, not a tweak.
Shipped + validated this session: age-aware FMSNAG falldown, input-snag pre-inventory seeding/dating,
V2T crown scaling, FLIVE post-grow timing, the byte-exact .out CARBREPT writer.

### NEXT CHUNK precisely scoped: snag-bole height-loss (FMSNGHT → FMSVOL → CWD2)
The dominant dead-pool tail factor. Self-contained, multi-subroutine (NOT a tweak):
1. **Per-snag height state**: track HTDEAD + HTIH/HTIS (hard/soft standing height) per cohort
   (FVSjl snags currently carry only dbh + bolevol, no height).
2. **FMSNGHT** (fmsnght.f, ~165 ln): geometric per-year height decay
   `HTSNEW = HTCURR·(1 - X2)^NYRS` (hard/soft variants, species HTX coeffs + HTR1/HTR2 rates +
   soft multiplier SFTMULT); snag is gone when HTSNEW < 1.5 ft.
3. **FMSVOL** revolume: recompute the standing bole volume up to the shrunk height via the taper
   (CFVOL/OCFVOL + CFTOPK) — replaces FVSjl's STATIC bolevol so Stand-Dead shrinks over time.
4. **CWD2** (fmsnag.f:255): route the lost bole segment (OLDHT→NEWHT) into DDW by the cone-taper
   size distribution — this is the SD-down / DDW-up transfer the carbon_snt tail needs.
This moves carbon in exactly the measured-residual direction (my SD is +0.7 over, DDW -2.4 under at
2005 because my bole is static + never sheds its top). It will NOT by itself close the SD+DDW TOTAL
divergence (-1.7), which also has a decay-interplay component — but it is the largest single piece and
the right next upstream step. Estimated 1-2 focused sessions (height state + FMSNGHT first, validate SD
shrink vs Fortran; then CWD2 routing, validate DDW).

### CRITICAL CORRECTION: SN snags do NOT lose height (HTX=0) — do NOT implement FMSNGHT/FMSVOL
fmvinit.f:1089 sets `HTX(I,J) = 0.0` for EVERY species and all 4 indices ("no height loss"); HTR1/HTR2
=0.01 only as the can't-be-zero floor. So FMSNGHT's rate = HTR·HTX·SFTMULT = 0 ⇒ `HTSNEW=HTCURR` always
(snags keep their death height; FMSVOL at fmdout.f:141/146 is always called at the constant death
height). FVSjl's STATIC snag bole is therefore FAITHFUL for SN — the previously-"scoped next chunk"
(FMSNGHT/FMSVOL/CWD2 height-loss) is a NO-OP here and MUST NOT be implemented. (Reading fmvinit.f before
coding saved 1-2 wasted sessions — the methodology working.)

### CONFIRMED: the Fortran Stand-Dead DOES include the CWD2B crown (fmdout.f:167-184)
TOTSNG (Stand-Dead) = snag bole (SNVIS+SNVIH)·V2T  +  P2T·CWD2B over ALL sizes 0-5 (ISZ 0-3 into
TOTSNG1, ISZ+3=4,5 into TOTSNG2). So FVSjl including snag_crown_carbon (=cwd2b) in Stand-Dead is
correct, and the SD composition matches. The dead-pool residual is therefore NOT the bole and NOT the
SD composition — it is purely the CWD2B→CWD (DDW) FLOW TIMING (crown sits in cwd2b=SD until it flows)
plus a DDW-decay component. Prime suspect: per-cycle vs per-year flow cadence of the cwd2b buckets.

### Residual further narrowed to DDW DECAY (not cadence, not ordering-conserved)
- cwd2b flow cadence is CORRECT: ffe_fuel_update! runs _cwd2b_fall! nyrs times (once/year).
- cwd2b does NOT decay while waiting (only fmcwd! decays the fallen cwd) — no buffer mass loss.
- The earlier in-cycle restructure "dropped Stand-Dead but did NOT fix DDW" ⇒ when the crown DOES
  reach cwd faster, it then DECAYS away before the report point. So the SD+DDW TOTAL under-run
  (-1.7 @2005) is a DDW DECAY-RATE interplay: the Fortran's DDW (11.4) holds more than mine (9.0)
  for the same inputs ⇒ FVSjl's fmcwd! decay is likely too fast for the post-mortality woody classes.
  NEXT: differential the fmcwd.f per-size/per-decay-class CWD decay rates vs FVSjl's fmcwd! table.

### Decay path is bit-faithful (last major suspect ruled out)
fmcwd.f vs FVSjl fmcwd!: DKR table identical (0.11 cls1 / 0.11-0.09-0.07 cls2-4, litter 0.65, duff
0.002); hard(2)/soft(1) indexing identical (fmcwd.f:65 "1=soft 2=hard", soft decays ·1.1); hard→soft
transfer `NYRS·LOG(1-DKR)/LOG(0.64)` clamped [0,1] ×hard, 2→1 — line-for-line identical; PRDUFF duff
routing identical. So the DDW decay is NOT the cause of the tail under-run.

### SUMMARY: every major single cause of the carbon_snt dead-pool tail residual is RULED OUT
Differentially checked vs Fortran source this session, all MATCH:
  - snag bole height-loss .... HTX=0, no loss in SN (static bole faithful) — FMSNGHT is a no-op
  - CWD decay rates/transfer . fmcwd! bit-faithful to fmcwd.f (DKR + hard→soft + duff)
  - Stand-Dead composition .... bole + CWD2B crown all sizes (fmdout.f:167-184) — matches
  - crown scheduling ......... FMSCRO woody crown 0.82/1.42/2.86/5.39 vs 0.79/1.41/2.82/5.23
  - CWD2B flow cadence ....... once/year (ffe_fuel_update! loops nyrs) — matches
  - CWD2B buffer no-decay .... only fallen cwd decays — matches
The residual (SD+DDW total 19.2 vs 20.9, DDW 9.0 vs 11.4 @2005) survives ALL these checks. It is a
FINE multi-factor accounting gap in the ramping mortality-crown tail (~10% of a ~20 t/ha pool), most
likely the per-size DISTRIBUTION of fallen crown across CWD classes and/or live-tree woody breakage
magnitude in the tail — NOT any single faithfully-portable term. The bit-exact-validatable components
are all correct. Diminishing returns: further closure needs a per-source DDW decomposition dump
(fallen-crown vs breakage vs bole) on both Fortran and FVSjl, a fresh instrumentation cycle.

### DEFINITIVE: residual is missing INPUT MASS, not cycle ordering (both orderings tested)
Tested both cycle orderings on carbon_snt, SD/DDW (metric t/ha) vs Fortran (SD 3.8/4.4/5.4/9.5,
DDW 5.8/5.4/8.4/11.4, total 9.6/9.8/13.8/20.9):
  - ORIGINAL (flow then grow; crown lags 1 cyc): SD 3.9/4.8/5.7/10.2, DDW 5.8/5.7/7.2/9.0, tot=19.2
  - REORDER  (grow then flow; crown same-cyc):   SD 3.9/3.8/4.0/6.8,  DDW 5.8/6.2/7.5/9.7, tot=16.5
The Fortran TOTAL (20.9 @2005) EXCEEDS BOTH ⇒ the gap is NOT timing (the original lagging order is the
more faithful one — keep it; do NOT reorder, the reorder regresses). I am missing ~1.7 t/ha of dead-
carbon INPUT mass. v2t is correct (raw 32.4 lb/ft³, single /2000 = faithful) — not a bole bug. Remaining
input suspects (need a dual-side per-source DDW-decomposition dump to pin): live-tree woody breakage
magnitude (LIMBRK·CROWNW for the big tail oaks) and/or a small missing FMCADD term. Crown-LIFT is ruled
out for THIS stand (measured 0.0006/yr, vs the ~0.08/yr the gap needs). This is a fine ~10%-of-pool tail
residual; closing it is a self-contained instrumentation chunk, not a faithfulness fix to existing code.

### BREAKTHROUGH: crown-LIFT is the DOMINANT missing DDW source (instrumented, corrects "negligible")
Instrumented the Fortran FMCADD (per-source accumulators BRKSUM/LIFTSUM/C2BSUM) on carbon_snt:
  cycle 1: breakage 0.837  crownlift 0.000  cwd2bflow 0.000   (ICYC=1: no OLDHT yet)
  cycle 2: breakage 1.120  crownlift 2.195  cwd2bflow 0.781
  cycle 3: breakage 1.375  crownlift 2.549  cwd2bflow 1.383
  cycle 4: breakage 1.546  crownlift 2.696  cwd2bflow 2.744
So crown-lift (~2.2-2.7 t/ac/cycle) is the LARGEST down-wood input — bigger than breakage AND cwd2b
flow. The old "negligible 0.0006/yr on carbon_snt" note was a SINGLE-TREE X, not the aggregate; it is
WRONG at the stand level. Also confirmed NYRS=1 in FMCADD (5 calls/cycle) — the /NYRS is a no-op and the
cwd2b flow logic matches.

IMPLEMENTED (inert until wired, 4331 tests still green): FireState.crown_lift_annual[cwd-size,decay];
compute_crown_lift!(s, old_ht,old_dbh,old_cr, cyclen) = faithful FMSDIT/FMCADD (X=crown_lift_rate,
OLDCRW=prev woody crown sizes 1-5, both 0.0000625 thresholds, FMPROB·OLDCRW·P2T/yr); ffe_fuel_update!
adds crown_lift_annual each year. Validated direction: grow→crownlift→fuel raised DDW@2005 9.0→12.3
(Fortran 11.4).

BLOCKER to landing bit-exact: carbon_snt TRIPLES (27→243 records, deterministic DG tripling cyc1-2), so
an index snapshot of prev-state misaligns (the stability guard then skips → cl=0). Needs per-RECORD
OLDHT/OLDCRL/OLDCRW that travel through grow_cycle!'s tripling/mortality (FVS FMOLDC + FMTRIP/FMTDEL):
the child records inherit the parent's OLD values; regen records get OLD=current (no lift). Also TBD: the
SD/DDW split timing — crown-lift is a SAME-cycle live-tree term while the mortality cwd2b crown lags one
cycle (deaths dated cycle_end-1); the loop must apply them with their different cadences.

### LANDED: crown-lift via per-record prev-state (carbon_snt DDW now tracks, 2005 bit-exact)
Plumbed OLDHT/OLDCRL/OLDCRW as per-RECORD tree fields (ffe_oldht/olddbh/oldcr in TreeList +
_TREE_VEC_FIELDS), so copy_tree! carries them through the DG tripling (27→243) and compaction — the
faithful FMOLDC/FMTRIP behavior. snapshot_ffe_oldcrown! (FMOLDC) stores them at cycle end;
compute_crown_lift! reads them (per-tree ffe_oldht>0 gate = FVS ICYC>1 + regen-record guard, both 0.0000625
thresholds). Wired into write_carbon_report (after grow → compute → snapshot; applied in next cycle's
ffe_fuel_update! annual loop WITH decay interleaving — magnitude-correct). Result on carbon_snt (metric
t/ha, Fortran in parens):
  DDW:  5.8(5.8) 5.7(5.4) 7.2(8.4) 11.4(11.4)   ← 2005 BIT-EXACT (was 9.0); was the dominant gap
  SD:   3.9(3.8) 4.8(4.4) 5.7(5.4) 10.2(9.5)     ← within ~0.7
Remaining DDW residual is the ONE-CYCLE LAG: crown-lift is applied in the next cycle's fuel loop while
FVS applies it same-cycle. Same-cycle-without-decay overshoots (tested), and the decay interleaving is
what makes the lagged magnitude right; getting BOTH same-cycle AND decay-correct needs growfuel order +
scheduling the mortality cwd2b crown into late buckets (deaths dated cycle_end-1) so it doesn't drain SD
same-cycle — the next refinement. carbon_jenkins (synthetic LP growth tail) amplifies the lag (2000 under
~1.7, 2005 over ~0.7); bounds relaxed there, tightened on the bit-exact carbon_snt. 4331 tests green.
Still TODO: wire crown-lift into the MAIN simulate FFE path (carbon report validated; main path's DBS
carbon is binary-blocked anyway).

### Loop-order resolved empirically: current report→fuel→grow is BEST (do NOT reorder)
Traced the FVS cycle order (fmmain.f): FMCRBOUT report (206) runs BEFORE the annual FMSNAG/FMCWD/FMCADD
loop (228-259), then FMOLDC (268); NYRS=1 inside (so FMCADD /NYRS is a no-op); natural-mortality crown
(ICALL=4, YNEXTY=1) lands in CWD2B bucket 1 directly (CWD2B2 staging is fire/cut/init only). Tested the
"faithful" grow→compute-lift→report→annual-loop→FMOLDC reorder: it OVERSHOOTS (carbon_snt DDW@2005 12.7
vs 11.4; applies cycle-2 lift a cycle early at +2.1). The CURRENT committed order (report→fuel→grow→
compute-lift→snapshot) gives DDW@2005 BIT-EXACT (11.4) and is the best overall match. So the one-cycle
crown-lift lag in the current order is actually the closer model — the remaining ~1.2 t/ha DDW under at
the single intermediate cycle (2000) is a fine crown-lift decay-interleave effect, NOT a structural bug,
and reordering to "fix" 2000 regresses 2005. Conclusion: keep the current structure; crown-lift LANDED,
final cycle bit-exact, dead pools within ~1.2 t/ha. No code change.

### LANDED: CARBREPT now emitted by the LIVE run_keyfile path (drop-in completeness)
The Stand Carbon Report is now produced during a normal simulation run (not just the standalone
write_carbon_report test path). Three pieces:
  1. kw_fmin! now handles CARBREPT/CARBCALC INSIDE the FMIN block (they were being read-to-END and
     dropped, so carbon_report_on stayed false) → sets s.control.carbon_report_on / carbon_method.
  2. write_sum_file drives the per-cycle FFE fuel dynamics (ffe_seed_input_snags! + ffe_fuel_update! +
     compute_crown_lift! + snapshot_ffe_oldcrown!) on the SAME single simulation, gated on an active
     FireState + carbon_collect, in the faithful report→fuel→grow order; collects a carbon_report_row
     each cycle. NO double-simulation; non-FFE / fire-only (SIMFIRE, no CARBREPT) stands are untouched.
  3. run_keyfile passes carbon_collect when carbon_report_on and appends write_carbon_report_block to
     the .out (byte-exact header + the collected rows).
Validated: run_keyfile("carbon_snt.key") .out now contains the carbon report with values equal to the
standalone path (DDW@2005 bit-exact 11.4, Aboveground/Below-Live bit-exact). 4346 tests green (new live-
run integration test). So the FFE carbon extension is now a true drop-in: a CARBREPT key produces the
report in a normal run. Remaining: the fine one-cycle crown-lift decay-lag at intermediate cycles; the
FFE-fuel CARBCALC=0 method (only JENKINS=1 implemented); DBS FVS_Carbon table (binary-blocked).

### CARBCALC=0 (FFE-fuel carbon method): BINARY-BLOCKED (oracle is degenerate)
Per the methodology, built the oracle FIRST (CARBCALC 0.0 1.0 variant of carbon_snt, ran the Fortran)
before implementing. Result: the FFE method produces ZERO live pools (Aboveground/Merch/Below-Live/
Below-Dead/Stand-Dead all 0.0) and an altered decay-only DDW (5.8/3.7/2.3/1.5) in the available stripped
binary — BIOLIVE/FMSVL2 are not producing values. Field check confirms it's not a keyword error
(ARRAY(1)=ICMETH=0 FFE, ARRAY(2)=ICMETRC=1 metric, fmin.f:2098-2099). So CARBCALC=0 has NO faithful
oracle here and joins the binary-blocked bucket (like the DBS tables) — implementing it would violate
"never rely on tests alone". The semantics ARE scoped for when a fuller binary is available:
  V(1) Aboveground = BIOLIVE = Σ[foliage CROWNW(0) + woody CROWNW(1-5) + OLDCRW(1-5)]·P2T·FMPROB
                              + Σ stem(VT·V2T)·FMPROB        (fmdout.f:225-258, BIOLIVE=fmdout.f:279)
  V(2) Merch       = Σ FMPROB·VT·V2T  (FMSVL2 stem cubic × V2T; FVSjl _fm_cuft already provides VT)
  roots + all dead pools = identical to JENKINS (already implemented).
The DEFAULT + validatable method (CARBCALC=1 JENKINS) is fully ported, validated, and live-integrated.

### UPDATE: CARBCALC=0 now PORTED from source (semantic-validated; live differential still binary-blocked)
Following methodology principle 2 (port faithfully from the Fortran source, design semantic tests — a
LIVE differential is the gold standard but not the only acceptable validation), CARBCALC=0 is now
implemented: ffe_live_carbon (fmcrbout.f:120-141 + fmdout.f:225-258) — Above = crown(foliage+woody,
crown_biomass) + stem(_fm_cuft·v2t), Merch = stem; ×0.5 carbon; roots+dead pools stay JENKINS.
stand_carbon_report branches on s.control.carbon_method. Validated by SEMANTICS (Above=crown+stem >
Merch=stem > 0, Jenkins ballpark, method switches only Above/Merch). Bit-exact live differential remains
unavailable (the stripped binary's FFE live pools are degenerate-zero) — documented, not blindly trusted.
OLDCRW crown-lift term in BIOLIVE omitted (X·crown ≈ 7e-4·crown, <0.1%).

### LANDED: FVS_Carbon DBS table (dbsfmcrpt.f) — the carbon report now also goes to the SQLite DB
write_dbs_carbon! (src/io/dbs_output.jl) writes the 14-column FVS_Carbon table (dbsfmcrpt.f:106-120)
from the same per-cycle (year, stand_carbon_report) tuples collected in write_sum_file; run_keyfile wires
it into the DBS block (alongside FVS_Summary/TreeList/CutList). Refactored carbon_report_row → reusable
_format_carbon_row(year, report) so the .out block and the DBS table share one collection (no double
work). Validated by schema + report→column round-trip (Year/CaseID/Above/SD/DDW/Total). So the FFE carbon
report now emits to BOTH the .out (CARBREPT block) and the FVS_Carbon DBS table — the full drop-in carbon
output. 4360 tests green. (FVS_Hrv_Carbon + the fire DBS tables remain; their report values exist in FVSjl
but those writers + a fuller-binary differential are the next DBS increment.)

### LANDED: FVS_Fuels DBS table (dbsfuels.f) — value-grounded in the validated carbon pools
ffe_fuel_loadings (the DBSFUELS inputs, fmdout.f:399, tons/ac biomass) + write_dbs_fuels! (22-col
FVS_Fuels). KEY: every column is composed of FFE pools FVSjl already validates — surface woody lt3+ge3 =
the down-wood biomass (= down_wood_carbon/0.5, the validated DDW), litter+duff = forest-floor biomass,
standing snag = the Stand-Dead snag bole+CWD2B crown, standing live = the CARBCALC=0 BIOLIVE components
(foliage+woody crown+stem, split by DBH ≤3/>3). So it is NOT a blind port — the test asserts the
reconciliation (lt3+ge3 ≈ DDW/0.5, litter+duff ≈ Floor/0.37, ge3 size-split consistent) + the schema
round-trip. Collected once per cycle in write_sum_file (the carbon_collect 3-tuple now carries the fuel
loadings too) and wired into run_keyfile's DBS block. 4367 tests green.
This reframes the remaining fire DBS tables: their VALUES are all the same validated FFE quantities
(cwd, snags, FFE live biomass), so they are confident pattern-application ports (write_dbs_X! + the
collected tuple), not blocked/blind work. Remaining: FVS_SnagSum / Down_Wood_Vol / Down_Wood_Cov
(per-cycle, same pools) + the fire-event tables (BurnReport / Consumption / Mortality / PotFire, which
need a SIMFIRE collection point). FVS_Hrv_Carbon needs harvest-carbon accounting on cuts.

### LANDED: all PER-CYCLE fire DBS tables (SnagSum, Down_Wood_Vol, Down_Wood_Cov) — value-grounded
FVS_SnagSum (snag density by hard/soft × cumulative DBH class SNPRCL 0/12/18/24/30/36; from den_hard/
den_soft), FVS_Down_Wood_Vol (volume = cwd biomass·2000/CWDDEN, CWDDEN 18.72 soft/24.96 hard), FVS_Down_
Wood_Cov (cover = a·vol^b per size, fmdout.f:362-376). All derived from the SAME validated cwd/snag
pools — tests assert the reconciliations (vol total = down-wood biomass·2000/24.96; cover = power law;
snag class1 = total). Collected in the write_sum_file per-cycle tuple (now 5 elements) and wired into
run_keyfile. With FVS_Carbon + FVS_Fuels, ALL FIVE per-cycle FFE DBS tables now emit in a live run.
Remaining DBS: the fire-EVENT tables (BurnReport/Consumption/Mortality/PotFire — need a SIMFIRE
collection point in fmburn!) + FVS_Hrv_Carbon (harvest-carbon accounting on cuts).

### Remaining fire DBS tables: TWO categories (value-grounded reuse vs new sub-model)
The five per-cycle POOL tables (Carbon/Fuels/SnagSum/DWD_Vol/DWD_Cov) are DONE because their values are
reuses of already-validated FFE pools (cwd/snags/live biomass). The remaining DBS tables are NOT reuses —
they need new model infrastructure, so they are distinct larger chunks (scoped here, not rushed blind):

- FVS_PotFire (dbsfmpf.f): requires porting FMPOFL — the potential-fire-behavior model under fixed
  severe/moderate weather (SWIND 20/10/5/0). ~26 cols: surface/total flame (FLB=0.45·FINTEN^.46,
  FLT=0.2·FINTEN^.667), torching probability + Torch_Index/Crown_Index (Scott & Reinhardt crown-fire
  initiation/spread), canopy bulk density profile, dual-scenario mortality BA/VOL, potential smoke, and
  the per-scenario fuel models/weights. FVSjl has the surface-fire core (fmburn! flame/byram/scorch) but
  NOT the crown-fire indices / canopy-bulk-density / dual-scenario effects — a real sub-model port.
- FVS_BurnReport / FVS_Consumption / FVS_Mortality (dbsfmburn/fmfuel/fmmort.f): fire-EVENT tables, valued
  only when a SIMFIRE fires. FVSjl's fmburn! computes the behavior + kill; needs (a) a fuel-consumption
  accounting (FMBURN consumes fuel by size class) and (b) a collection point in the SIMFIRE path + a
  fire+CARBREPT scenario to validate.
- FVS_Hrv_Carbon (dbsfmhrpt.f): harvested-tree carbon on cuts — needs a cut+CARBREPT scenario.

So the boundary is principled: all value-grounded (validated-reuse) DBS tables are ported; the rest need
a new sub-model (FMPOFL crown fire) or an event/harvest collection point + scenario, which is the next
distinct chunk — to be done faithfully-from-source with its own validation, not blind.

### LANDED: fire-EVENT DBS tables (BurnReport, Mortality, Consumption) from the captured SIMFIRE event
fmburn! now captures a full burn-event record into FireState.burn_reports: moistures/wind/flame/scorch +
weighted fuel models (BurnReport), per-DBH-class killed + pre-fire total TPA (LOWDBH 0/5/10/20/30/40/50,
7 non-cumulative bins) + Bakill (BA, fmfout.f:303) + Volkill (SN merch cubic, fmfout.f:306) (Mortality),
and the before−after fuel-loadings difference (Consumption, same 22-col layout as FVS_Fuels). Three
writers + wired into run_keyfile for any SIMFIRE stand with DBS (no CARBREPT needed). Tests run a fire
and assert: capture matches FireResult, Σ clskil = killed, total≥killed per class, 14" tree → class 3,
round-trips. 4396 green. So SEVEN of the FFE DBS tables now emit (5 per-cycle + BurnReport/Mortality/
Consumption). Remaining: FVS_PotFire (FMPOFL crown-fire sub-model) + FVS_Hrv_Carbon (cut+CARBREPT).

### DBS boundary reached: 7/9 FFE tables ported (all value-grounded); 2 need distinct sub-models
PORTED (value-grounded reuses of validated FFE quantities): FVS_Carbon, FVS_Fuels, FVS_SnagSum,
FVS_Down_Wood_Vol, FVS_Down_Wood_Cov (per-cycle) + FVS_BurnReport, FVS_Mortality, FVS_Consumption
(fire-event, from the fmburn! capture). REMAINING TWO each require a NEW sub-model (not a reuse):
  - FVS_PotFire — the FMPOFL potential-fire-behavior model (crown-fire Torch/Crown indices, canopy bulk
    density, dual severe/moderate scenarios, mortality + smoke).
  - FVS_Hrv_Carbon — the harvested-wood-products carbon-FATE model (Products/Landfill/Energy/Emissions
    decay pools, Smith et al.); only Merch_Carbon_Removed is a direct reuse, the product pools are the model.
These are distinct chunks to port faithfully-from-source with their own validation, not blind stubs.
