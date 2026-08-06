# EM #137 dense-cohort self-thin — REPRODUCED + characterized (2026-08-06)

## Reproduction (constructed dense PLANT scenario)
em_plant.key plants only 400 TPA (non-self-thinning; jl TPA BIT-EXACT there, see #137 note). To exercise the
dense-cohort self-thin, cloned it to em_plant_dense.key with `PLANT 1990 3.0 6000.0` (6000 TPA). jl vs FVSem_clean:

  year   liveTPA  jlTPA   liveBA  jlBA
  2000   5455     5455    2       0
  2010   5350     5350    14      9
  2020   5246     5247    30      23
  2030   5145     5145    51      39
  2040   5045     5045    78      58     <- TPA still bit-exact, BA already −26%
  2050   4810     4947    108     81     <- TPA diverges: jl UNDER-thins
  2060   4354     4783    133     106
  2070   4001     4369    158     127
  2080   3699     3973    182     154
  2090   3167     3566    190     184    <- jl TPA +12.6%, BA −3% (converged)

## Characterization
- DIRECTION: jl UNDER-thins (TPA HIGH: 3566 vs 3167 = +12.6% by 2090). This CORRECTS the goal-file's "over-kill"
  label — it is UNDER-kill, the SAME direction as BM #140 (resolved).
- ROOT is UPSTREAM of mortality: TPA is bit-exact 2000-2040 while BA is already one-directional LOW (58 vs 78 at
  2040, −26%). So jl under-grows the dense cohort's BASAL AREA first; the self-thin under-kill (2050+) is the
  DOWNSTREAM consequence (smaller QMD → higher T85D10/BAMAX self-thin target → fewer trees killed → QMD-feedback).
- The BA lag is a TIMING lag (converges to −3% at 2090): the dense planted small-tree cohort grows diameter slower
  in jl early, so it crosses the small→large-tree threshold later. NOT the crown-init gap (planted trees get crowns
  at PLANT, crown≠0) — it is the small-tree DG/HTG under DENSE competition (SMHTGF/SMDGF) or the DGSCOR calibration
  on the dense cohort.

## Next (fresh chunk)
Instrument jl small-tree DG + HTG (em/regent.f SMHTGF/SMDGF) vs live on the dense planted cohort at an early cycle
(2010/2020) — measure per-tree height + diameter growth to see whether jl under-grows HEIGHT (delaying the 4.5'
crossing) or DIAMETER directly, under the high-competition (dense TPCCF/BAL) inputs. Whether it is a real bug or the
accepted dense-cohort small-tree timing tail (cf. the CORNERED EM/IE growth-only ~7% BA tail) depends on that
measurement — it converges by 2090, which leans cornered, but the +12.6% TPA under-thin is a notable persistent
divergence worth the instrument pass. Harness: em_plant_dense.key in /workspace/.emwork (live) + run_keyfile (jl).
