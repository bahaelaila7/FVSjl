# Confirm B3 — FFE crown-lift: inventory crown snapshot creates a cycle-1 lift FVS gates off

## Verdict: CONFIRMED

The Julia port snapshots the inventory crown state into `ffe_old*` BEFORE the first growing
cycle. That gives the first growing cycle a valid `OLDCRW`, so `compute_crown_lift!` produces a
NONZERO crown-lift down-wood input for cycle 1. FVS suppresses exactly this contribution two
independent ways (initialized `OLDHT=0` that is never set until the *end* of cycle 1, plus an
explicit `IF (ICYC.GT.1)` gate). Instrumented live FVS confirms cyc1 lift = 0.0, cyc2 = 947.

## FVS logic (what the Fortran actually computes)

- `fminit.f:970-971` — at FFE init, `OLDHT(I)=0.0`, `OLDCRL(I)=0.0` for all records.
- `fvs.f:203` sets `ICYC=0`; the inventory pass runs DISPLY/VOLS/etc. but **never calls TREGRO**
  (TREGRO is only reached at `fvs.f:376` after `ICYC=ICYC+1` at `:367`). TREGRO → GRADD → FMMAIN
  (`gradd.f:118`) → FMOLDC (`fmmain.f:268`). So **FMOLDC / FMMAIN do not run at inventory** and
  `OLDHT` stays 0 going into the first growth cycle.
- `fmsdit.f:93` — the entire crown-lift block is gated `IF (ICYC.GT.1) THEN`. At `ICYC=1` (the
  first growing cycle) the block is skipped and `OLDCRW(I,J)=0.0` (fmsdit.f:116-118 else-branch is
  not even reached; the values are left as the FMINIT/FMCROW state and contribute nothing). The
  per-year lift consumed by FMCADD (`fmmain.f:241`) is therefore 0 for all of cycle 1.
- First nonzero lift is at `ICYC=2`: `OLDBOT=OLDHT-OLDCRL` uses the end-of-cycle-1 `OLDHT` recorded
  by FMOLDC (`fmoldc.f:52-53`), `NEWBOT=HT-HT*ICR/100` uses the cycle-2 state.

Net: FVS produces ZERO crown-lift for the first growing cycle, by design.

## Julia logic (what the port actually computes)

- `src/io/summary.jl:100-101` (main `.sum` path) — before the cycle loop, for any FFE-active stand:
  `fill!(s.fire.crown_lift_annual, 0f0)` then `snapshot_ffe_oldcrown!(s)`. The snapshot
  (`fuel_additions.jl:163-172`) writes `ffe_oldht/olddbh/oldcr = current (inventory) tree state`.
- `src/engine/fire/carbon.jl:388` — the carbon-report test path does the identical inventory
  snapshot, with a comment (`:403-404`) claiming "Both no-op on the first grow (ffe_old* unset ⇒
  zero), matching FVS's ICYC>1 gate." That claim is **self-contradicted by line 388**, which sets
  `ffe_old*` precisely so they are NOT unset on the first grow.
- First loop iteration `c=0` (= FVS `ICYC=1`): `grow_cycle!` runs (summary.jl:183 / carbon.jl:400),
  THEN `compute_crown_lift!(s, per)` (summary.jl:187 / carbon.jl:406). Inside
  `compute_crown_lift!` (`fuel_additions.jl:132-151`) the guard is `oldht = t.ffe_oldht[i];
  oldht > 0f0 || continue`. Because of the inventory snapshot, `ffe_oldht > 0`, so the guard PASSES.
  `crown_lift_rate(oldht=inventory, oldcrl=inventory, ht=post-cyc1-growth, …)`
  (`fuel_additions.jl:29-35`) returns `rise/oldcrl/cyclen > 0` since the grown trees' crown base
  rose above the inventory base → NONZERO `x` → nonzero `cl[sz,dkcl]`.
- There is no `ICYC>1`/cycle-index gate anywhere in `compute_crown_lift!`; the only suppressor is
  the `oldht>0` guard, which the inventory snapshot defeats.

So Julia computes a cycle-1 crown-lift (the inventory→cyc1 growth delta) that FVS gates off.

## Precise semantic divergence

FVS: `OLDHT` initialized 0 (fminit.f:970) AND `IF(ICYC.GT.1)` (fmsdit.f:93) ⇒ first-grow lift = 0.
Julia: inventory `snapshot_ffe_oldcrown!` (summary.jl:101 / carbon.jl:388) sets `ffe_oldht>0` and
there is no `ICYC>1` gate ⇒ the `oldht>0` guard passes and the first-grow lift is nonzero. The
port mirrors neither of FVS's two suppressors. The "FVS calls FMOLDC in the inventory FMMAIN"
justification in the code comments is false — FMMAIN/FMOLDC never run at `ICYC=0`.

## Faithful fix (do NOT apply here)

Restore FVS's first-cycle-zero behavior:

1. Delete the inventory `snapshot_ffe_oldcrown!(s)` call at `src/io/summary.jl:101` (keep the
   `fill!(s.fire.crown_lift_annual, 0f0)` at :100) and the matching call at
   `src/engine/fire/carbon.jl:388`. With `ffe_oldht` left 0 (its constructed default), the existing
   `oldht > 0f0 || continue` guard makes the first growing cycle's lift exactly 0 — reproducing
   `fminit.f:970` (`OLDHT=0`).
2. (Belt-and-suspenders, mirroring `fmsdit.f:93` directly) add an explicit cycle gate so
   `compute_crown_lift!` is a no-op on the first growing cycle (Julia loop `c==0`).

Both reduce to: the inventory→cyc1 delta must contribute zero crown-lift.

## Upstream rank: MID (UPSTREAM for SIMFIRE stands)

`crown_lift_annual` is applied by `ffe_fuel_update!` into `fs.cwd` (the down-wood pools,
`fuel_additions.jl:190-191`). `fs.cwd` feeds: (a) the DWD/Carbon report columns — LEAF; and (b)
FMCFMD's (SMALL,LARGE) fuel-model selection (summary.jl:91-96 documents this is why cwd must evolve)
→ fire behaviour (byram/flame/scorch) → FMEFF fire mortality → the live tree list and downstream
cycles. So for a stand that actually burns in an early cycle this is UPSTREAM (changes mortality);
for non-fire FFE stands it is a LEAF feeding only the carbon/down-wood reports.

## Reachability

Exercised. The main path (`write_sum_file`, `ffe_on`) runs for every FFE-active stand in the
keyword-coverage suite, and the dedicated carbon-report path (`run_carbon_cycles!`) is the
carbon_snt / carbon_jenkins tests. Fixing it WILL move the DDW/Carbon down-wood columns in the
early cycles (it is a silent value path; the suite will register the change).

## Masked-bug warning

The inventory snapshot was added to "recover ~1.9 t/ac of fine down-wood" (comments at
summary.jl:101-102 / carbon.jl:385-387), and the "FFE crown-lift LANDED" memory note records
`carbon_snt DDW@2005 bit-exact` achieved WITH this snapshot. Removing the snapshot will therefore
likely REGRESS the early-cycle DDW. That regression is a SIGNAL, not a reason to keep the bandaid:
the underlying masked bug is the ONE-CYCLE LAG in how Julia applies crown-lift. Julia computes the
cyc-N lift AFTER `grow_cycle!` (summary.jl:187) and applies it in the NEXT cycle's pre-growth fuel
loop (summary.jl:173-182, next iteration). FVS computes the cyc-N lift in FMSDIT at the START of
cyc N (after that cycle's HT is grown in GRADD) and consumes it in the SAME cycle's FMMAIN annual
loop. So Julia is consistently one growth-delta behind: in the cyc-2 application window FVS applies
the (cyc1→cyc2) delta while Julia (with the snapshot) applies the spurious (inv→cyc1) delta. The
snapshot papers over the lag by injecting a leading term. The genuinely faithful fix is to (a) gate
the first-cycle lift to zero (this flag) AND (b) realign the lift computation/application to the
same cycle as FVS (a separate, more-upstream fix). Sequence: fix the lag first, then the cyc1 gate
falls out naturally; fixing only the gate without the lag will leave the cyc-2 window under-applying.


---
## ⚠️ RETRACTED — VERDICT CORRECTED TO FALSE_POSITIVE
This report mis-mapped the jl loop timing. FVS UPDATE (update.f:65/115) runs at gradd.f:180 AFTER FMMAIN
(gradd.f:118), so cyc1 FMOLDC captures the PRE-growth=inventory crown; the jl inventory snapshot is its faithful
analog (inv→postcyc1 lift added in CYC2, matching FVS ICYC=2=947). carbon_snt DDW is bit-exact vs LIVE FVS WITH
#6. #6 is FAITHFUL — keep it. See INDEX.md + FAITHFULNESS_AUDIT.md retraction.
