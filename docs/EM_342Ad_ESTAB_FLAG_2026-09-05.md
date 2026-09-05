# EM — FLAGGED (NOT cornered) potential real bug: 342Ad semi-desert establishment UNDER-production (2026-09-05)

Session: https://claude.ai/code/session_014RMPT9hw2GKinUa9hrorT7. Found at cap-7 of the EM FIA sweep.

## What
4 stands in ecoregion 342Ad (Intermountain Semi-Desert) tripped the eco_guard. Adjudication:
- seed_test (5 seeds): all 4 SEED-INVARIANT, one-directional jl-UNDER TPA (−77, −220, −184, −76).
- Unbiased 342 population read: ZERO stands with TPA>200 in EM's swept range (342 barely exists for EM —
  Eastern Montana is Great Plains + Rocky Mtn foothills; semi-desert is marginal). So two-sidedness CANNOT be
  established, and per doctrine these were NOT cornered — they are FLAGGED as an open lead and removed from the
  dig-queue so the sweep proceeds. eco_guard STILL stops on 342.

## Mechanism (trace of 303089960489998, O = FVSem_clean, J = FVSjl)
    yr    O[TPA] J[TPA]   note
    2014   114    114     identical
    2024   112    105
    2034   192    166     O regenerates ~+80, J ~+60
    2054   431    206     O regen cohort +239 (192->431!), J barely +43  => dTPA -225
    2064   427    202
The ORACLE produces large REGENERATION/INGROWTH cohorts (TPA 192->431 at 2054); jl produces a much smaller
regen cohort. QMD tracks this (O 9.0 = diluted by many new small trees; J 11.8). => jl UNDER-produces
establishment in the 342 semi-desert habitat. This is the OPPOSITE sign of the fixed IE M333 AUTOES
OVER-production; ESRANN-fixed establishment => seed-invariant (matches seed_test). BA/QMD/volume follow the
regen deficit. NOT a volume bug, NOT the count-straddle, NOT the establishment cuft leak.

## Verdict: REAL one-directional establishment under-production, but MARGINAL (4 stands, near-empty ecoregion
## for EM). Correctly NOT cornered. FLAGGED CNs: 12336192010690, 303089960489998, 4754520010690, 511377763126144.

## PLAN to fix (separate establishment task)
Root-cause the AUTOES/regen rate for the 342 semi-desert habitat vs FVSem_g16: instrument ie_autoes_establish!
(EM shares IE AUTOES) on 303089960489998 at the 2054 cycle — compare the per-point tally (ESB1/STOMLT/tally_pt)
and the habitat/plant-association gate that sets the regen count, since the deficit is a smaller regen COHORT
(not sizing). Likely a habitat-code or per-point establishment-count gate that differs for the 342 semi-desert
plant association. Low priority (marginal ecoregion) but real; fix belongs with the establishment-keyword work.
