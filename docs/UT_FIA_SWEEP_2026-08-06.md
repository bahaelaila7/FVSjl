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
