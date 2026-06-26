# FVSjl Southern Drop-In — Task Tracker

Persistent progress tracker for the SN drop-in-replacement effort (survives container
restarts). Goal: **FVSjl (Southern) is a bit-exact drop-in for live FVSsn; the only
accepted divergences are ULP floating-point and the COMPRESS eigensolver.**

Canonical test runner: `julia --project=. test/runtests.jl`. Current suite:
**4467 pass + 18 broken** (broken = tracked, documented divergences, NOT failures).
Live-Fortran oracle: `/tmp/FVSsn_new` (rebuilt via `tests/fortran_baseline.sh`).

Keyword-coverage harness: `test/integration/test_keyword_coverage.jl` over
`test/keyword_coverage/scenarios/*` (37 scenarios). See
[memory: fvsjl-keyword-coverage-state] for the gate design.

---

## DONE (committed)
- [x] Wire the 37-scenario keyword-coverage harness into the main suite.
- [x] 10-yr-cycle **mortality** root-caused + fixed (linear FINT-extrapolated SDI trajectory, morts.f:225).
- [x] **s37_thinauto** MAI terminal-row quirk (evtstv.f:414/disply.f:392) — bit-exact.
- [x] Column-aware gate: structural cols strict ULP; merch/board cols 0.3% (Scribner/Behre FP) → **s14/s34** conform.
- [x] Structured-yaml **named schemas** for THINHT/THINRDEN/THINAUTO/VOLUME/BFVOLUME/MCDEFECT/BFDEFECT/SERLCORR/RESETAGE.
- [x] Robust yaml writer: lossless **raw passthrough** on parse failure / unknown keyword → yaml==key 36/37.

## OPEN — engine divergences (must reach ULP/compression-only)
- [ ] **s5/s9** — 10-yr board-foot tail (~4%): calibrated-species HTCALC/AGET + tripling-variance DG/HTG precision under a non-5 period. (TPA already within 3.) [memory: fvsjl-10yr-cycle-mortality]
- [ ] **s26_estab** — CORRECTED (the earlier "245.4 establishment-reduction" was a `.trl` column-parse ERROR — disregard it). Via a targeted GRADD-end LP-sum print in a recompiled FVS, the real per-cycle LP TPA is: FVSjl 2005=270.0/2010=63.31/2015=7.27 ; FVS 2005=270.0/2010=60.31/2015=8.75. **ESTABLISHMENT IS BIT-EXACT** (270 both at 2005). The divergence is the DENSE LP-cohort SELF-THINNING MORTALITY 2005→2010 (FVS 270→60.31, FVSjl 270→63.31, ~5% under-kill). LP are dbh 2.6→3.2 (near dbh_zeide≈3 SDI threshold) ⇒ likely small-tree SDI-inclusion or VARMRT distribution for a dense sub-merch cohort. NOT establishment, NOT the linear-G fix (5-yr=identity). TOOLING WIN: stripped FVS global DEBUG segfaults (fvsvol.f:530/fvs.f:376 uninit WRITEs), but targeted per-cycle prints (recompile one .f at -O0, relink) WORK. INSTRUMENTED (mortality.jl S26DBG dump): at 2005 cycle tt=753, sdimax=359, TOTkill=305.9; the LP cohort (PCT≈5, the smallest trees, VARMRT eff≈0.79) takes LPkill=206.69 → surv 63.3, vs FVS 209.69 → 60.31. So it's a ~3-tree (~1.5% of cohort) difference in the VARMRT kill DISTRIBUTION to the dense small cohort and/or the SDI tokill total — the deepest mortality-kernel layer. Next: targeted prints in FVS varmrt.f to compare the LP per-record PCT/efficiency (snt01 mortality stays bit-exact, so the fix must be small-tree/dense-cohort specific). ROOT FOUND (instrumented FVS morts.f vs FVSjl at the 2005 cycle): T=753.09 (match), DIA0=5.916 (match), SDIMAX=359.08 (match), but **D10=6.594 (FVS) vs 6.57 (FVSjl)**. DIA0 matches ⇒ start DBH identical ⇒ the DIFFERENCE IS THE DIAMETER GROWTH G (FVS grows slightly more) → lower tn10 → more tokill → ~3 more LP killed. dbh_zeide=0 so all trees incl the LP are in D10. So s26's VARMRT difference is DOWNSTREAM of a sub-percent diameter-growth precision diff — the SAME CLASS as s5/s9. UNIFYING INSIGHT: the remaining growth residuals (s5/s9 + s26) are all small DG/HTG precision in non-snt01 stand/species calibrations (s26 stand = forest type 231 + planted LP), amplified downstream (board-foot / mortality). Fix both via the per-tree DG of forest-type-231 / small-tree growth vs FVS dgdriv (snt01 forest type 117 stays bit-exact). **COMPLETE ROOT (7-layer trace):** s26 stand WITHOUT the PLANT is BIT-EXACT ⇒ entire divergence is the established LP. Chain: LP self-thin 270→63.31 vs FVS 60.31 ⇒ mortality D10 6.57 vs 6.594 ⇒ DG ⇒ established LP DBH 2.5409 (FVSjl) vs 2.5609 (FVS) at SAME height 20.0 + IDENTICAL htdbh coeffs (SNALL[*,13]=243.86/4.2846/-0.4713). FVSjl establishment.jl:98 uses PLAIN _htdbh_dbh; FVS derives it in REGENT (regent.f:317-321, LHTDRG-calibrated ht-dbh) → +0.02 dbh. FIX = derive established dbh via the calibrated REGENT path (snt01 has no establishment ⇒ stays bit-exact). ✅ PARTIAL FIX LANDED: established dbh now adds 0.001*HK (regent.f:334) → LP dbh 2.5609 matches FVS (establishment.jl). REMAINING: the LP SMALL-TREE HEIGHT GROWTH 2005→2010 = FVSjl 20.0→24.907 vs FVS 20.0→24.967 (0.06ft/1.2% less) → lower dbh growth → D10 → ~3-tree mortality diff. Next = REGENT small-tree height increment for the LP (small_tree_growth.jl) vs FVS regent/htgf. ✅ 2 FIXES LANDED this session: (1) established dbh +0.001*HK (regent.f:334); (2) NPTIDS=points_inv-nonstockable (esplt2.f:74) → regen record count 44→50 matches FVS, RNG draws re-aligned. s26 diff cells 16→9; LP 2010 63.31→58.89 (FVS 60.31). REMAINING: a smaller RNG-state residual (LP↔stand interaction shifts the RANN draws when the planted cohort is present; the no-plant stand is bit-exact). Both fixes keep BARE/snt01 bit-exact (4487+28).
- [ ] **s32_volume** — VOLUME zeroes SCF (cols past 80) → all trees prod=01. Per-tree .trl: FVS scuft=0 below dbh 10, mcuft=0 below ~6; FVSjl leaks small-tree merch (+19 each). mrules.f fixes DON'T apply (sawDib=6 overshoots; merchL=10 BREAKS 3 bit-exact defaults) → SN R8 uses merchL=8/sawDib=7-9; multi-threshold. Default path bit-exact.
- [ ] **s22_compress** — COMPRESS different eigensolver — **ACCEPTED per spec** (keep as documented-broken).

## OPEN — yaml / conversion (user-requested)
- [x] **Conversion tool operational (drop-in paths)** — `test/integration/test_translate.jl`: .key→.yaml→.key engine-equal + .tre→.csv preserves TreeRecords (the engine reads either form). Regression test in suite.
- [ ] **csv→tre legacy re-emission bug** — `write_tree_file` drops an F-field that shares a column with a packed nI1 field in the SN T-specifier layout (e.g. field `T60,F3.1` value 5.0 → 0, overlapping `T54,7I1`). Tracked @test_broken (10). Does NOT affect the drop-in (engine reads the .csv); fix the overlap resolution in write_tree_file for a byte-faithful legacy export.
- [ ] **s20_spgroup** — 2-record SPGROUP round-trip (group name + species-list record) — handle in the hierarchical redesign.
- [ ] **Task 8 — hierarchical, order-aware YAML redesign.** ⚠️ ORDER MATTERS: the current flat form is a `list` so order is preserved. The redesign must NOT become fully order-independent — some stand operations are genuinely order-dependent and must keep their relative order:
      • a SPGROUP / species-group definition must precede a THIN that references the group;
      • multiple same-cycle thinnings/operations apply in sequence;
      • COMPUTE / event-monitor variable definitions vs. their use.
      Make date-scheduled keywords order-free where safe, but preserve declared order within a stand for operations whose effect depends on it (model as ordered sub-lists / explicit sequence, not an unordered map).
- [ ] **Ensure the YML/CSV ↔ .key/.tre conversion tool is operational** (`bin/fvsjl-translate.jl`): verify both directions round-trip on all coverage scenarios + the examples; add a regression test.
- [ ] **Keyword documentation**: update `docs/KEYWORDS.md` to document each keyword in BOTH yml and .key form, how they translate to each other, a good explanation of the keyword + each named parameter, and usage examples.

## OPEN — ill-posed scenarios
- [ ] **s30_thinqfa**, **s36_readcord** — no FVSsn baseline (bad/edge param layout; SCF fields past col 80). Re-author as realistic.

## Deferred
- NE variant port (explicitly deferred while the SN drop-in is in progress).
