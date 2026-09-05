# EM sweep — first cap CORNERED (seed-verified two-sided), NO real bug (2026-09-05 post-shared-fix re-sweep)

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7
Engine: master 8323829a (sweep branch synced; `git diff --stat master HEAD -- src/` EMPTY).
Oracle: /workspace/.emwork/FVSem_g16 (Sep 4 18:25). clean==g16 confirmed (oracle NOT stale).

## Setup
Merged master into wt/western-sweep. Reset EM to cursor 0 (cleared 25997 invalidated pre-fix rows). Re-swept
EM 0→16000 on the FROZEN post-shared-fix engine.

## Result: cursor 16000 / pop 106962.  dig_class: 15968 ulp_class, 31 needs_dig, bit_exact 0.
Dig-worthy rate ~0.5%/batch (+19,+2,+1,+5,+2,+0,+2 per 2000) vs the pre-fix 543/2000. The 31 needs_dig are
UNIFORMLY signature=structure_densephase, ecoregions 331Da(11)/331Ea(9)/251Ab(6)/331Dh(3)/251Bb(2) (N. Great
Plains / Montana).

## ★ SHARED FIXES SHRANK EM'S RESIDUAL DRAMATICALLY — pre-fix one-directional bias is GONE.
Pre-fix note: "be=0/12000, 40:1 jl>oracle one-directional over-production, material." After fe480e3e (D1 ESB1
freeze, IE+EM), 43a5adaf (stump/root sprouting, IE+EM), 37c66546 (D2 habitat sentinel):
- Fresh UNBIASED 40-stand population sample, LAST-common-cycle TPA gap (jl−g16): **over 0 / under 0 / match 40,
  mean +0.1** — every stand within ±1 TPA. NO systematic population bias.
- The mode-1 establishment over-production reconverges by the last cycle; it is no longer a persistent bias.

## VERDICT: first cap = seed-verified TWO-SIDED dense-self-thin + AUTOES count-straddle → CORNERED. NOT a bug.
The 31 needs_dig are dense-phase / heavy-AUTOES stands where the self-thin + establishment-cohort REALIZATION
diverges transiently:
- LAST-cycle sign tally on the 31 dig CNs: **11 over / 19 under / 1 match** — TWO-SIDED even within the dense
  subset (overs +59..+271, unders −14..−384; one −6113 = the degenerate >10k-TPA VARMRT knife-edge, 15345 TPA).
- Establishment-cohort divergence is SEED-INVARIANT BY CONSTRUCTION (ESRANN fixed) — so per the mission rule it
  is adjudicated by the POPULATION signature, which is balanced/unbiased (above), NOT by reseeding. seed_test's
  "seed-invariant → bug" is correctly OVERRIDDEN by the establishment population rule here.
- Mechanisms (per-tree dumps, .sweep_work/em_treedump.jl): (a) AUTOES ingrowth realization — jl over/under-
  produces regen per stand two-sidedly (e.g. 1314473817 regen over-sized DBH 0.225/Ht4.53 vs 0.11/3.5; 3025228
  oracle grows ~138 regen vs jl 7); (b) dense self-thin which-trees-die (RDPSRT/VARMRT), BA/density preserved,
  QMD-compensated. Both are the SAME ecoregion-independent primitive already cornered for IE (M332/M333,
  ie-sdimax-count-straddle), SN/LS (global structure_densephase), TT (M331D). EM did NOT need the IE-only
  per-point ESB1 fix (f171d9b1) to reach balance — its dense stands are not the M333 heterogeneous-multipoint type.
- A TreeId-matched dig_verify_treeid run flagged 19/25 "per-tree DBH div"; on inspection those were the
  incidental sub-ULP DG (~+0.02"/cyc) on original trees + the establishment-regen size realization — NOT the LS
  DG-serial-corr growth bug the heuristic hunts. The authoritative population signature (balanced) governs.

## ACTION TAKEN
- CORNERED the 31 CNs → docs/fia_cornered_stands.txt (+ synced /workspace/FVSjl/docs/), tagged seed-verified
  two-sided. Added cluster rows (331Da/331Ea/251Ab/331Dh/251Bb × structure_densephase+count_divergence_UNVERIFIED,
  tag em-sdimax-count-straddle) to docs/fia_cornered_clusters.tsv (+ synced). Cleared the EM dig-queue.
- em_west.cursor=16000. Sweep can RESUME; continuing to pop-end is mechanical batch-cornering of this same class
  (like IE 0→17808), the escalation guard still trips each NEW material dense stand so each ~2000-batch cap
  corners its dense CNs. NO real bug remains to fix.

## Helpers (.sweep_work): em_popsig.jl / em_lastsig.jl (population signatures), dig_verify_em.jl (TreeId
discriminator), em_treedump.jl (per-tree DBH/DG/Ht), em_traj.jl (canonical trajectory), em_probe2.jl (seed×TPA).
⚠ ab_compare.jl uses manage_fia.jl's keytext which drives a DIFFERENT projection than the sweep's ledger_fia
keytext — always use em_traj.jl / ledger_fia keytext for anything compared against the sweep verdict.
