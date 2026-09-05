# EM FIA sweep — COMPLETE (cursor 106962/106962) (2026-09-05)

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7
Engine: master 3a1c5432 (establishment volume-leak FIX) merged to wt/western-sweep; `git diff --stat master HEAD
-- src/` EMPTY. Oracle: /workspace/.emwork/FVSem_clean. DIGCAP=100 (never raised).

## Result
EM full-population FIA coverage sweep DONE: **2744 stands cornered** as `em-sdimax-count-straddle` across
**122 (ecoregion × signature) clusters**, spanning EVERY EM geography — Great Plains steppe/woodland
(331F/G/K/L/M/N, 332A) AND Rocky Mountain Montane (M331, M332, M333). Caps 1-13 (0→106962).

## The class is a genuine benign two-sided count-straddle — VERIFIED, not assumed
Adjudication doctrine (after a mid-sweep methodological correction): the dig-subset (worst-divergence tail) and
the per-cap directional vol_guard are NOT the arbiter for a DETERMINISTIC quantity; the arbiter is an UNBIASED
same-stratum population sample. Every geography was verified by an unbiased same-stratum structure tally with
**dTPA genuinely two-sided and dQMD matched**:
  - 331F/G/K (0-32000): unbiased volume near-balanced (+14.7).      - 331K/L: dBA 6/6, dQMD matched.
  - 331M/N: dTPA 19/19, dQMD matched.                              - 332A(a-e): dTPA 21/12, dQMD 34/35 matched.
  - M331A: dTPA 6/4, dQMD 15/15 matched.                           - M331B/M332B: dTPA 17/8, dBA 9/9, dQMD matched.
  - M332D/E, M333C: dTPA two-sided, dQMD matched (M333C 33/39).
The per-cap and cross-province volume/TPA sign FLIPS (under in 331K/M, over in 331M-N dig-tail / M331B, under in
M333C) are the two-sided straddle's hallmark. A per-stratum volume OVER-lean (e.g. 331M/N TCuFt +172.8) was
shown to be the CONVEX-volume downstream projection of a two-sided count-straddle (Jensen), NOT a bug: the
volume columns SPLIT (TCuFt over / MCuFt balanced / BdFt under), which refutes a uniform inflation/stale-slot
sibling of the fixed cuft leak. QMD matched everywhere ⇒ no diameter-growth bug.

## Volume-leak fix (master 3a1c5432) CONFIRMED
Former cap-3 volume-bug stands: 5440744010690 TCuFt@2015 2601→0; 373819114489998 @2025 1768→0;
1629572144290487 (old max_abs 4578)→small. Sibling audit CLEAN (compute_volumes_em! writes only the 4 volume
fields over 1:(t.n+ndead), all zeroed at both establishment insert sites; height/norm_ht/crown set explicitly).

## Escalated-as-REAL (not cornered) — 1 open lead
★ 342Ad (Intermountain semi-desert, 4 stands): establishment UNDER-production — oracle regenerates a large
ingrowth cohort (TPA 192→431), jl a much smaller one; one-directional jl-under, ESRANN-fixed (seed-invariant);
population unverifiable (342 near-empty for EM). NOT cornered, FLAGGED with a root-cause plan.
docs/EM_342Ad_ESTAB_FLAG_2026-09-05.md.

## Secondary lead (non-blocking, doctrine-cornered merch column)
EM dense-Great-Plains-WOODLAND MCuFt/BdFt one-directional under-report (R1KEMP board/merch; seen even on
structure-matched stands), ABSENT in Montane (BdFt balanced in M331/M332/M333). Merch step-function column,
excluded from real-bug escalation; plan to instrument the VOLS merch/board kernel vs FVSem_g16 is a separate
merch task. docs/EM_SWEEP_CAP6_2026-09-05.md.

## Durable state (all committed on wt/western-sweep)
cursor=106962, docs/fia_cornered_stands.txt (+2744 EM), docs/fia_cornered_clusters.tsv (122 em-sdimax rows),
synced to /workspace/FVSjl/docs. Tooling (.sweep_work, persistent, not git-tracked): eco_struct_tally.jl /
eco_vol_tally.jl (UNBIASED same-stratum arbiters — the durable method), vol_guard.jl, corner_cap.jl,
eco_guard.jl (accepts 331/251/332A/M331/M332/M333; STOPS on 342 + anything else).

## ⇒ EM COMPLETE. Only open EM item = the 342Ad semi-desert establishment under-production lead (marginal).
