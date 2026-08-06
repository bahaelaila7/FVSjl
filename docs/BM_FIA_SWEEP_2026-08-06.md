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

### #149 LOCALIZED (2026-08-06): height-growth 4.5' THRESHOLD effect
Instrumented jl bm/regent.f small-tree loop on 504443988: ALL small trees (d<3) stay hk=H+HTG ≤ 4.5' after cyc1
(measured hk=3.0-4.42 for 0.1" ES/WL/LP/AF seedlings; h_start=1.01, htgr=2.0-3.4). bm/regent.f:390-394: for hk≤4.5
DG(K)=0, DBH(K)=D+0.001·HK (tiny) — the small-tree DIAMETER only grows once HK crosses 4.5' (the HTDBH/AX-BX branch).
⇒ jl's seedlings land JUST below 4.5' (LP hk=4.42) → DG=0 → BA under-grows 2×; live's evidently cross 4.5' and grow
diameter (QMD 0.615 vs jl 0.42). So a SMALL height-growth under-prediction has an OUTSIZED effect via the 4.5'
threshold. NEXT: instrument live bm/regent.f htgr/HTG + PCTRED (density modifier) for a seedling on 504443988 vs jl —
is jl's small-tree HTG deterministically low (PCTRED density over-suppression on the dense ~44000-TPA stand, or the
bm_smhtgf POTHTG) or a ZZRAN realization (htgr includes zzran·0.1)? TPA is EXACT ⇒ the mortality-RNG is synced, so
a deterministic HTG gap is the leading hypothesis. Threshold-sensitive ⇒ could be partly cornered like IE WH #146.

### #149 ROOT-CAUSED (2026-08-06): crown_pct=0 ⇒ VIGOR floor ⇒ HTGR under → 4.5' threshold miss
Instrument-replay (bm/regent.f:350 vs jl regent.jl) on 504443988 sp2(WL) 0.1" seedling:
  live: POTHTG=11.05 PCTRED=0.8186 VIGOR=0.426/0.382 CON=1.0 → HTGR=3.86/3.46 (crosses 4.5' next subcycle → DG grows)
  jl:   POTHTG=12.06 PCTRED=0.8186 VIGOR=0.30       CON=1.0 → HTGR=2.96 (stays <4.5' → DG=0)
POTHTG/PCTRED/CON match (jl POTHTG even slightly higher). The ENTIRE gap is VIGOR: jl 0.30 (the FLOOR) vs live
0.426/0.382. VIGOR=150·CR³·exp(−6·CR)+0.3 ⇒ jl's crown_pct=0 (measured) forces the 0.3 floor; live's seedling has
crown ratio ~12% (VIGOR 0.426 ⇒ CR≈0.12). ⇒ ROOT: jl assigns ZERO crown ratio to the dense FIA seedlings (DBH=0.1,
HT=None) where live estimates ~12%. The suppressed HTGR keeps them below the 4.5' breast-height threshold ⇒ DG=0 ⇒
BA 2× low. NEXT (the fix): jl's crown-ratio initialization for read sub-1" trees with no measured crown must estimate
CR (the BM/shared crown model, not leave 0). Likely affects ALL variants' dense-seedling stands (crown init shared).
DETERMINISTIC (TPA exact) ⇒ fixable, NOT cornered. ⇒ #149 root = crown-ratio init, not the regent/DG itself.
