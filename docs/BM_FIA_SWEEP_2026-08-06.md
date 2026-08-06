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

### #149 FIX LOCATION (2026-08-06): crown-init skips sub-1" trees
crown_ratio_update! (bm/crown.jl:61) SKIPS d<1 trees at lstart (`(d < 1f0 && lstart) && continue`, comment "small
trees → REGENT bm/crown.f:237") ⇒ dense read seedlings keep crown_pct=0. But live's REGENT is "CALLED FROM CRATET
DURING CALIBRATION AND FROM TREGRO DURING CYCLING" (bm/regent.f header) — the CRATET-time REGENT call SETS the
small-tree crown BEFORE the growth-cycle regent's VIGOR reads it. jl never runs that calibration-time small-tree
crown-set ⇒ crown=0 at the first regent ⇒ VIGOR=0.30 floor ⇒ HTGR under ⇒ 4.5' miss ⇒ DG=0 ⇒ BA 2× low.
THE FIX: port the small-tree crown assignment from bm/regent.f (the JCR/ICR set for d<XMAX at LESTB/calibration) so
read sub-1" trees get a crown ratio (~12%, not 0) before the growth regent. Likely a shared pattern across variants'
regent (VIGOR = f(crown) everywhere). ⇒ #149 is a CROWN-INIT chunk (deterministic, fixable), not regent/DG. Verify
on 504443988 (BA 38→~80) + ttt01/emt01/etc. no-regress (they have measured crowns, so init-skip is inert there).

### #149 COMPLETE FIX SPEC (2026-08-06): BM missing CI-parity crown-init (bm_dubscr)
The BM setup dispatch (simulate.jl:98-100) is MISSING the crown-init that CI has (simulate.jl:102-104):
  compute_density!(s); crown_ratio_update!(s, s.variant; lstart=true)
CI's comment there LITERALLY describes #149: "Without it, 0.1" seedlings keep crown_pct=0 ⇒ ... seedlings never
reach breast height (4.5') ⇒ DBH growth skipped ⇒ small-tree DG low". BM's crown_ratio_update! (bm/crown.jl:61)
SKIPS d<1 at lstart; CI's instead has a d<1 lstart BRANCH (ci/crown.jl:128-136) that estimates the crown via
ci_dubscr → clamp[10,95] → t.crown_pct[i]. THE FIX (mirror CI, source live bm/crown.f:336/370 CALL DUBSCR):
(1) port bm_dubscr (BM has NO dubscr yet; KT/CI do — likely the shared DUBSCR crown model with BM coeffs);
(2) replace bm/crown.jl:61 skip with a d<1 lstart branch: cr=bm_dubscr(...); crown_pct=clamp(cr·100,10,95);
(3) add `compute_density!(s); crown_ratio_update!(s,s.variant;lstart=true)` to the BM setup branch.
VALIDATE: 504443988 BA 38→~80 (=live), bmt01 no-regress (its trees have measured crowns ⇒ line-60 skip keeps them
inert). ★ LIKELY CLUSTER-WIDE: check EM/TT/UT setup too — only CR+CI currently call the lstart crown-init; the
dense-seedling HTGR-via-VIGOR(crown) starvation is the SAME mechanism (BM #149, and the EM/IE/UT/TT dense-seedling
findings may share this crown-init root where the DG isn't the specific variant bug I fixed). This is the campaign's
likely UNIFYING root: read sub-1" seedlings need a dubscr crown estimate at inventory, which most variants skip.
