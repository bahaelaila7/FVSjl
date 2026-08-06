# BM real-FIA sweep (2026-08-06) — dense CONIFER-seedling DG UNDER-growth [#149]

6-stand BM sample (bm_test.db from the 70GB DB, VARIANT='BM') vs FVSbm_clean, 2 cyc. 3 dense-seedling + 3 treed.

| stand | type | cyc1 jl (TPA/BA) | cyc1 live | verdict |
|-------|------|------------------|-----------|---------|
| 504443988126144 | seed (ES/WL/LP, 44025 TPA) | 38775/38 | 38775/80 | ★ TPA EXACT, BA 2× LOW |
| 850488769290487 | seed | 29687/115 | 26136/136 | under-thin +14% TPA, BA −15% |
| 1544772750290487 | seed (AF/LP) | 32044/84 | 34473/96 | BA −12.5% |
| 1127530489290487 | treed | 739/154 | 757/160 | −2%/−4% (cornered) |
| 1127530768290487 | treed | 56/104 | 56/104 | ✓ bit-exact |
| 1127530848290487 | treed | 293/51 | 293/51 | ✓ bit-exact |

## Finding [#149]
Dense CONIFER-seedling stands (ES/WL/LP/AF at DBH=0.1", ~35-44k TPA) UNDER-grow BA ~2×. ★ On 504443988 the TPA is
EXACT (38775=live) ⇒ NO mortality/RNG divergence ⇒ a DETERMINISTIC small-tree DG under-growth (QMD jl 0.42 vs live
0.615). Distinct from IE WH #146 (ZZRAN tripled-record realization) and BM #140 (mortality). Treed stands bit-exact
⇒ specific to the dense conifer-seedling cohort. BM = Stage-SDI + SMHTGF/SMDGF small-tree. NOT cornered. NEXT:
instrument jl bm/regent.f SMDGF per-tree DG vs live for an ES/LP 0.1" seedling — likely density-feedback
over-suppression (dense ~38000 TPA) or an SMDGF coeff/form gap. META: the real-FIA sweep found a real bug in the
5th variant swept (EM/IE/UT/TT all had one) — dense sub-1" seedlings remain the systematically under-tested regime.
