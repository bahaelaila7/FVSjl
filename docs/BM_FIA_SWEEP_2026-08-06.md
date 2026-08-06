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

### #149/#150 FIX ATTEMPTED (2026-08-06): bm_dubscr ported, over-corrected + broke bmt01 → REVERTED
Ported bm_dubscr (bm/dubscr.f logistic, BCR0-10+CRSD extracted, sp13/14/16/18 linear-rescale) + a d<1 lstart
branch + the setup crown-init call (mirror CI). RESULT on 504443988: BA 38→203 (OVER-shot live 80) + TPA 38775→
15977 (over-kill) ⇒ the crown estimate came out ~39% (VIGOR→1.0 max) not live's ~12%, so the DUBSCR INPUTS
(TPCCF/AVH/BA/RMAI) I passed differ from what live passes at inventory. ALSO broke bmt01 (establishment.jl:266 /
height_growth.jl:26 error — the lstart crown call + compute_density! disturbed bmt01's ESTAB path). REVERTED
(doctrine #4). ⇒ NEXT: (1) instrument live bm/dubscr.f to dump ISPC/D/H/BA/TPCCF/AVH/RMAI/CR for a 504443988
seedling ⇒ match jl's inputs exactly (the coeffs are verified-extracted; the inputs are the discrepancy — likely
AVH or TPCCF or RMAI≠0); (2) gate/guard the lstart crown call so it doesn't perturb the ESTAB path (bmt01 has
ESTAB). Root (crown=0→VIGOR floor) + fix DIRECTION (add DUBSCR crown-init) CONFIRMED; the port needs input-exactness.

### CORRECTION (2026-08-06): bmt01 error is PRE-EXISTING, not the crown fix
Re-checked: the REVERTED (original) code ALSO errors on bmt01 (MethodError getindex(::Nothing) in the ESTAB path)
— bmt01.key has ESTAB 1992, and the BM establishment path is unwired (like utt01's :essprt_fsp). So the crown-init
fix did NOT break bmt01; bmt01 was already erroring. ⇒ the ONLY real issue with the fix attempt is the crown
OVER-estimate (~39% vs live ~12%) = wrong DUBSCR INPUTS (TPCCF/AVH/BA/RMAI), coeffs are correct. The no-regress
check must use a BM stand WITHOUT ESTAB (bmt01-first-stand growth-only / emc2), NOT the full bmt01.key. NEXT: (1)
instrument live bm/dubscr.f inputs on 504443988; (2) match jl's; (3) re-apply + validate vs a no-ESTAB BM stand.
SEPARATE pre-existing bug logged: jl BM ESTAB path errors on bmt01.key (BM AUTOES/establishment unwired).

### #149 ROOT-CAUSE DEFINITIVELY MEASURED (2026-08-06 cont.): CRATET DUBSCR runs on a DEAD-INCLUSIVE density + site-species RMAI
Instrumented live bm/dubscr.f (D<0.5 dump, ICYC) + maical.f on 504443988126144. The lstart DUBSCR inputs are:
  ISPC=2 D=0.1 H=1.0  BA=55.560  TPCCF=84.3369  AVH=85.0748  RMAI=128.0  → CR .119/.076/.056/.051 (WL)
  ISPC=3/4 (DF/GF) → CR .950 (saturated). MAICAL: ISISP=4(GF) SSSI=61 ISICD=15 RMAI=128 (ADJMAI grp9, capped).
jl's crown-init (reverted attempt) passed BA=35.36 TPCCF=63.39 AVH=45.22 RMAI=71 → WL CR .28-.35 (2× high) →
seedlings over-grow → BA 200 (vs live 80). TWO measured discrepancies, BOTH now understood:

1. **RMAI 71 vs 128 — FIXED (commit 1fbdd92).** jl bm_sitset! reset ISISP=0 and re-derived it from the ecoclass
   (LP idx7, SITEAR 70) instead of KEEPING the DB SITE_SPECIES=17→GF(idx4, SITEAR 61). Live bm/sitset.f:121 only
   overrides ISISP when ≤0. Fix: seed isisp from p.site_species. ADJMAI(15,61) grp9 = 129→cap128. Bit-exact on
   8 BM sweep stands. (This is a general FIA-DB site-resolution bug: also mis-propagated SI to all species via HTCALC.)

2. **DEAD-INCLUSIVE CRATET density — the remaining #149 blocker (NOT YET FIXED).** Live's lstart DENSE (feeding
   DUBSCR) COUNTS the inventory standing-dead trees (HISTORY=6/8). This stand: 15 overstory records, 5 live (HIST=1)
   + 10 dead (HIST=6/8, D 6.8-12.8"). jl correctly partitions the 10 as dead (ndead=10) ⇒ its live-only top-40 AVH
   = 45.2 (5 live overstory 20 TPA + 20 TPA of h=1 seedlings). LIVE's AVH=85.07 = the DEAD-INCLUSIVE top-40:
     (127·.999+109·.999+96·6.018+88·6.018+89·6.018+57·6.018+118·6.018+55·6.018+74·1.9)/40 = 85.08 ✓ (hand-verified).
   Likewise BA 55.56 (live+dead) vs jl live-only 35, TPCCF 84.3 vs 63. So live's CRATET init density includes the
   HIST=6/8 dead trees; the projection density (the .sum BA=35) excludes them. jl has no dead-inclusive density pass.
   ⇒ The crown-init (bm_dubscr, coeffs+ADJMAI+RMAI all verified) CANNOT be bit-exact until the lstart DUBSCR sees a
   dead-inclusive BA/AVH/TPCCF. Implementation care: (a) compute dead-inclusive scalars into TEMPs, don't leave them
   in p.basal_area (the grow cycle recomputes live-only density before use, so it's safe if scoped to the crown call);
   (b) only the LIVE seedlings draw BACHLO in DUBSCR (dead-tree crown dubs come AFTER in bm/crown.f, so live-seedling
   RNG order is unperturbed — jl can skip dead-tree crowns); (c) no-regress gate = the 8 BM sweep stands + a no-ESTAB
   growth stand (NOT full bmt01.key, which errors on the unwired ESTAB path).

★ META (extends the campaign META): CR/CI HAVE the lstart crown-init and are "complete", but their synthetic test
stands carry NO HISTORY=6/8 standing-dead inventory trees — so the DEAD-INCLUSIVE-density requirement was never
exercised there either. Real FIA stands routinely carry inventory dead. This is a SECOND cluster-wide gap layered
under #150: the lstart crown-init density must be dead-inclusive. Verify CR/CI on a real-FIA stand with inventory dead.

### #151 EXACT dead-inclusion + backdating rule — SOURCE-VERIFIED (base/dense.f LBKDEN pass)
bm/cratet.f:194 `LBKDEN = IDG.LT.2` → CALL DENSE → :610 CALL CROWN (→DUBSCR). So the DUBSCR density is the LBKDEN
(backdated) DENSE. base/dense.f:76-128 defines it exactly:
 1. **Live trees (DO 5 I=1,IREC1): BACKDATED DBH** WK3(I)=SQRT(D²·R), R=1−(2·D·G−G²)/D², G=DG(I)/BRATIO(IS,D,HT)
    (past bark-adjusted DG); if G invalid, R=BAGR (the stand-avg BA-growth ratio over the DG-sample trees). ⇒ BA/PCCF
    reflect the START-of-period (smaller) diameters.
 2. **HISTORY 6/7 dead (IMC=7, died ≤5yr): INCLUDED at CURRENT DBH** (dense.f:85 WK3=DBH).
 3. **HISTORY 8/9 dead (IMC=9, older dead): EXCLUDED** (dense.f:86 WK3=0).
AVH (top-40, height-weighted) is unaffected by the DBH backdating (monotone ⇒ same IND rank ⇒ same trees, CURRENT
heights) — that is why jl's AVH=85.07 was already bit-exact while BA/PCCF were not.
jl's c29c1c7 crown-init instead uses CURRENT-DBH density over ALL dead (HIST 6 AND 8) ⇒ BA 70.8/PCCF 94.2 vs live
55.6/84.3. THE #151/#150 FIX: the lstart crown-init must run on the LBKDEN density — backdate live DBH (jl's
calibrate_diameter_growth! already computes per-tree backdated D internally; expose/reuse it), include only IMC=7
(HIST 6/7) dead at current DBH, exclude IMC=9 (HIST 8/9). This is the SHARED mechanism for the whole #150 cluster
rollout (every variant's CRATET feeds DUBSCR the LBKDEN density). Near-floor crown ⇒ #149 net (38→87) already good;
this tightens BA/PCCF to bit-exact and de-risks KT/IE/EM/TT/UT.

### #151 DEEPER THAN A BACKDATE — the LBKDEN density is a calibration interpolation (base/dense.f:257-297)
Tracing the control flow of the LBKDEN DENSE (the one feeding DUBSCR): it is NOT "current density on backdated dbh".
 - **BA (dense.f:257-264, LBKDEN block)**: `RAT=FINTH/FINT; TEMP2=(BA−OLDBA)·RAT+OLDBA; BA=OLDBA; OLDBA=TEMP2` —
   the BA that LEAVES DENSE (and that DUBSCR reads) is set to **OLDBA**, with RELDEN←RELDM1 restored. So the DUBSCR
   BA/RELDEN are the PREVIOUS-state calibration-density values, interpolated by FINTH/FINT — not the raw backdated BA.
   Nailing jl's 70.8→live 55.56 requires reproducing OLDBA/RELDM1 as they stand at the lstart CRATET call (measure them).
 - **PCTILE skipped** under LBKDEN (:272 `IF(LBKDEN) GO TO 55`).
 - **AVH (:285-297) computed REGARDLESS of LBKDEN** — uses IND (dbh-sort) + PROB(TPA) + HT(height). Only the SORT RANK
   reflects backdating (monotone ⇒ same top-40 trees), heights/TPA unchanged ⇒ AVH invariant. THIS is why jl's AVH=85.07
   was already bit-exact with the current-dbh dead-inclusive density, while BA/PCCF were not.
⇒ REVISED #151 plan: (a) AVH is already correct — do NOT disturb it. (b) BA/RELDEN need the LBKDEN OLDBA/RELDM1
interpolation — instrument dense.f to dump BA/OLDBA/RELDEN/RELDM1 at the CRATET LBKDEN call, then replicate. jl's
_backdate_dbh! (southern/diameter_growth.jl:266, already BM-aware) is the WK3 half; the OLDBA/RELDM1 interpolation is
the missing half. NOTE: the crown is near-floor for the #149 seedlings ⇒ the committed c29c1c7 (BA 38→87 vs 80, 0
regressions) already MEETS the bar; #151 is a bit-exactness refinement, best done with its own instrument-replay pass.
