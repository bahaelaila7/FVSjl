# UT real-FIA sweep (2026-08-06) — dense woodland-seedling mortality OVER-KILL [#147]

6-stand UT sample (ut_test.db from the 70GB DB, VARIANT='UT') vs FVSut_clean, 2 cyc. 3 dense-seedling + 3 treed.

| stand | type | cyc0 jl:live | cyc1 jl (TPA/BA) | cyc1 live | verdict |
|-------|------|-------------|------------------|-----------|---------|
| 12449454010690 | seed (sp814 oak, 44979 TPA) | 48470/116 = | 26688/88 | 43220/154 | ★ jl kills 45% vs 11% |
| 276412671489998 | seed | 45262/83 = | 27796/64 | 40012/98 | ★ −31% TPA |
| 31538752010690 | seed | 27332/146 = | 14457/89 | 24756/160 | ★ −42% TPA |
| 11751442010690 | treed | 2503/117 = | 2209/120 | 2492/134 | −11% (cornered-ish) |
| 11755085010690 | treed | 5520/63 = | 5164/72 | 5342/62 | close |
| 11936970010690 | treed | 537/106 = | 537/120 | 522/118 | ✓ close |

## Finding [#147]
The 3 dense-seedling stands (dominated by sp814 = Gambel oak woodland-hardwood, ~45000 TPA at DBH=0.1") show jl
OVER-KILLING ~4× (cyc1 TPA 30-45% below live). cyc0 bit-exact. Mortality (TPA over-kill) dominates; BA also low.
NOT cornered (30-45% is far beyond the tie-break bar). UT uses Zeide SDI self-thin (vs EM/BM Stage SDI). Likely the
Zeide self-thin over-kills the dense sub-1" woodland-seedling cohort, OR the woodland species' mortality coeffs.
NEXT: instrument jl mortality! vs ut/morts.f (Zeide DR10/self-thin target) on 12449454 at cyc1. The 3 treed stands
are bit-exact-or-cornered ⇒ the bug is SPECIFIC to the dense woodland-seedling cohort. META: the real-FIA sweep on a
NEW variant (UT) immediately surfaced a real bug the ttt01/utt01 synthetic stands (conifer-heavy) never exercised —
same lesson as the EM CRVAR find (#145): population sweeps expose what single test stands miss.

## FIX 1 (2026-08-06, 6e57347) — Zeide DR10 (partial resolution of #147)
ROOT (ut/morts.f:218-219,260-263): UT is Zeide-SDI but jl's utah/mortality.jl computed dq10/dq0 as QMD (sqrt),
not Reineke DR10=(Σp·(D+G)^1.605/T)^(1/1.605). On dense sub-1" cohorts QMD over-stated D10 (0.7785 vs 0.5187) ⇒
TMD10 uncapped ⇒ TN10 low ⇒ RN over-kill. FIXED. Post-fix vs live: 31538752 14457→24611 (=24756 ✓), 11751442
2209→2487 (=2492 ✓), 11936970 522 bit-exact; 12449454 26688→29750, 276412671 27796→29750 (both improved).
REMAINING (2 stands): jl kills to TN10=29750, live to 43220 — measured live self-thin T=35000 vs jl 48470 (RN
0.016 vs 0.047). The SDI loop sums all P (DBHZEIDE=0) ⇒ T should be 48470; live's 35000 is likely the MORTS
IPASS QMD-convergence iteration (kill→recompute DQ10/T→re-kill) that jl's single-pass utah/mortality.jl lacks.
NEXT: port the IPASS loop (ut/morts.f IPASS) to utah/mortality.jl. Also a pre-existing SEPARATE bug: jl utt01
ERRORS (KeyError :essprt_fsp) — UT establishment/sprout path not wired (utt01.key fires ESTAB; unrelated to mort).

## CROSS-VARIANT (2026-08-06, #148): TT + CI share the Zeide-QMD mortality bug
Found while fixing UT #147: TT and CI are BOTH Zeide-SDI (LZEIDE=.TRUE.) but their custom mortality! uses
dq10=sqrt(QMD), NOT Reineke DR10 — the identical bug. ttt01/cit01 (conifers) hide it; dense small-tree stands
expose it. TT (teton/mortality.jl:145) uses the TMD10/TN10 self-thin ⇒ DR10 applies directly (like UT). CI
(centralidaho/mortality.jl:44) is BAMAX-based (DELTBA/ba10/tb) ⇒ needs per-term analysis (DELTBA=0.005454·D10²·T
uses D10=DR10 for Zeide, but the BA→TPA conversion is QMD). ★ TT DR10 swap ATTEMPTED → CRASHED ttt01 (DomainError
−0.0545 at volume.jl:61: (D+G)^1.605 surfaced a latent D+G<0 tree the old QMD (D+G)² tolerated) → REVERTED
(doctrine #4). Tracked as #148. The UT #147 DR10 fix (6e57347) stands — validated, no crash on the UT sweep stands.

### #148 TT DR10 fix — surfaces MULTIPLE latent TT bugs (2026-08-06, deferred)
Re-attempted the TT DR10 mortality fix; it cascades through THREE latent jl bugs (all tolerated by the old QMD
d*d/(d+g)² but not the Zeide ^1.605, and all masked on ttt01's normal treelist — the DR10 kill change alters
which trees survive to volume): (1) fpow(negative,3f0) DomainError in the woodland volume cubic → FIXED generally
in fmath.jl (gfortran-match); (2) a tree with NEGATIVE DBH reaches d^1.605 (jl-specific — live D≥0; old d*d
tolerated); (3) after guarding both, a BoundsError in compute_volumes_tt!:61 (t.saw_cuft_vol[i] OOB) — the altered
kill desyncs the treelist array sizing. ⇒ TT DR10 is CORRECT (Zeide) but blocked on these latent bugs + the MSTEM/
FCLASS woodland-volume gap. REVERTED (doctrine #4). The general fpow fix is kept+committed. UT #147 stands.

### #148 CORRECTION: CI verified NOT affected (2026-08-06, doctrine #2)
ci/morts.f:262-263 uses DQ10=SQRT (QMD) in its BAMAX-based DELTBA=0.005454·DQ10²·T / BA10 / TB mortality — QMD is
CORRECT there (BA≡0.005454·QMD²·T). CI is Zeide only in the SDIMAX INPUT, not the mortality diameter metric. jl's
centralidaho/mortality.jl:44 (dq10=sqrt) MATCHES ci/morts.f ⇒ NO bug; CI #142 stays the cornered DGSCOR verdict, do
NOT swap CI to DR10. The Zeide-QMD bug is confined to the TMD10/TN10-self-thin variants: UT (fixed) + TT (blocked).
