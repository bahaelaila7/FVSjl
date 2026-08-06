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
