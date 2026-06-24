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
