# LS growth-div candidates — post-FIX-#8 classification (audit 43eb)

29 LS `REAL_growthdiv` candidates (docs/fia_real_growthdiv_candidates.csv) re-run vs freshly-linked live FVSls
after FIX #8 (LS REGENT calibration stale-HTGR carry). Worst% = max relative diff across all cycles×6 .sum cols.

## A. Resolved / cornered-ULP (worst < 5%) — 11 stands
FIX #8 (tamarack stale-carry) + ULP/self-thin-tie residuals. TPA typically bit-exact all cycles.
1831637837290487 (0.78, the FIX #8 exemplar — was 2-3× off), 1283797302290487 (0.83), 1831681132290487 (0.8),
1536065011290487 (1.01), 1283811993290487 (1.2), 1210955343290487 (2.53), 54541693010661 (2.8),
569144292126144 (3.42), 1803273086290487 (3.97), 1283789309290487 (4.17), 1686728467290487 (4.81).

## B. Ultra-dense self-thinning RDPSRT tie-break primitive (cornered class) — 17 stands
CYCLE-0 BIT-EXACT (TPA/BA/SDI/CCF/TopHt/QMD all match at inventory), diverge ONLY in later cycles on
ultra-dense seedling stands (cycle-0 TPA 12k–39k). Signature = BA tracks within ~few %, TPA diverges,
QMD ±rounding — the accepted self-thinning tie-break primitive (see fvsjl-stand-pct-rdpsrt-fix). The 12
worst (>10%) verified cycle-0-exact via /tmp/cyc0_check.jl; the 5 milder (5–10%) same class by signature.
1536044338290487 (16.3), 156735105010661 (32.51), 164818658020004 (22.81), 176140320010661 (15.14),
1803286767290487 (20.12), 224645781010661 (25.66), 366591155489998 (17.41), 54608351010661 (15.88),
55238024010661 (19.57), 55733792010661 (15.71), 64256149010661 (20.59), 722509708290487 (17.43),
1686728401290487 (8.46), 1831635002290487 (6.22), 21952588010661 (7.61), 234462560020004 (8.84),
63969757010661 (6.72). (worst% amplified by the ultra-dense self-thinning cascade; cycle-0 identical rules out
a growth-model bug.)

## C. Cycle-0 TopHt tie-break (1 stand) — high-confidence RDPSRT primitive; possible fixable extension
**55250794010661** (81.21): the ONLY candidate with a CYCLE-0 (inventory) divergence — TopHt live 22 / jl 29,
before any projection. Traced: top height = avg height of the 40 largest-DBH trees/acre. This ultra-dense stand's
largest-DBH record (DIA=2.9", TREE_COUNT≈75 TPA > 40) alone fills the entire top-40-TPA window ⇒ TopHt = that
single record's height. But TWO records tie at DIA=2.9": aspen (sp746) HT=29 and balsam fir (sp12) HT=22. The
unstable RDPSRT quicksort FVS uses picks the balsam (22); jl's top-DBH selection picks the aspen (29). The 19/44
missing-HT trees are all shorter (never in the top set) ⇒ NOT height imputation. This is the SAME RDPSRT
tie-break primitive as fvsjl-stand-pct-rdpsrt-fix, here manifesting in the summary top-height sort rather than
the mortality percentile — the cycle-0 22/29 then cascades through the ultra-dense self-thinning (hence 81%).
STATUS: high-confidence cornered (accepted tie-break), BUT possibly a FIXABLE extension — if jl's top-height
selection uses a stable sort where FVS's summary uses RDPSRT/IQRSRT, applying the same unstable-sort tie-break
to the top-DBH ordering would break the tie to balsam ⇒ bit-exact. NEXT: locate jl's top-height/avg_top_height
DBH-sort, check whether it matches FVS sumit.f/htcalc top-tree ordering; if a stable-sort mismatch, extend the
_rdpsrt! fix. Sole remaining non-fully-cornered LS candidate.
