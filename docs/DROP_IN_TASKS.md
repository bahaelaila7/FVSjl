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
- [ ] **s5/s9** — 10-yr board-foot tail (~4%): calibrated-species HTCALC/AGET + tripling-variance DG/HTG precision under a non-5 period. (TPA already within 3.) [memory: fvsjl-10yr-cycle-mortality] GENUINE BUG FOUND (verified vs instrumented FVS, NOT applied): FVSjl's AUTCOR used the wrong 'old' period for the first 10-yr cycle (autcor(10,10) not autcor(10,5)=measurement base) AND didn't carry PVMLT=VMLT across cycles (dgdriv.f:116/121). Fixing both makes cov/corr match FVS EXACTLY (4.35029/0.24042 then 5.31562/0.18082) and the early cycles go bit-exact + board-foot shrinks. BUT it EXPOSES a separate 10-yr MORTALITY residual (one cycle TPA 360 vs 350). LANDED (committed) WITH its companion mortality fix: diameter_growth.jl: cyc==0 oldp=5; add c.vmlt=vmlt after corr (snt01 uniform-5yr unaffected). COMPANION FIX LANDED (mortality.jl:295 d10n recompute now uses the linear _mort_traj_g, not the raw sqrt fint-yr growth — the sqrt form understated d10n on a 10-yr cycle and forced a spurious extra QMD iteration that FVS does not do). RESULT: timeint10 FIRST 10-yr cycle (2000) now BIT-EXACT (TPA/BA/cuft/bdft), and cycle-2 tn10 matches FVS post-clamp (424.86/411.85 vs 425.0/412.0). EARLIER NOTE: with the autcor fix, cycle-1 tt/dia0/d10/sdimax MATCH FVS exactly (589.65/4.701/6.41/348.44), but the Pretzsch self-thinning TN10 differs (FVSjl≈522.4 vs FVS≈516.5) → cycle-1 under-kills ~6 TPA. So the companion fix is in _pretzsch_tn10 / the QMD-convergence iteration at the 10-yr cycle (the SDI inputs match; the tn10 line/iteration is the gap). Land BOTH together.
- [ ] **s26_estab** — CORRECTED (the earlier "245.4 establishment-reduction" was a `.trl` column-parse ERROR — disregard it). Via a targeted GRADD-end LP-sum print in a recompiled FVS, the real per-cycle LP TPA is: FVSjl 2005=270.0/2010=63.31/2015=7.27 ; FVS 2005=270.0/2010=60.31/2015=8.75. **ESTABLISHMENT IS BIT-EXACT** (270 both at 2005). The divergence is the DENSE LP-cohort SELF-THINNING MORTALITY 2005→2010 (FVS 270→60.31, FVSjl 270→63.31, ~5% under-kill). LP are dbh 2.6→3.2 (near dbh_zeide≈3 SDI threshold) ⇒ likely small-tree SDI-inclusion or VARMRT distribution for a dense sub-merch cohort. NOT establishment, NOT the linear-G fix (5-yr=identity). TOOLING WIN: stripped FVS global DEBUG segfaults (fvsvol.f:530/fvs.f:376 uninit WRITEs), but targeted per-cycle prints (recompile one .f at -O0, relink) WORK. INSTRUMENTED (mortality.jl S26DBG dump): at 2005 cycle tt=753, sdimax=359, TOTkill=305.9; the LP cohort (PCT≈5, the smallest trees, VARMRT eff≈0.79) takes LPkill=206.69 → surv 63.3, vs FVS 209.69 → 60.31. So it's a ~3-tree (~1.5% of cohort) difference in the VARMRT kill DISTRIBUTION to the dense small cohort and/or the SDI tokill total — the deepest mortality-kernel layer. Next: targeted prints in FVS varmrt.f to compare the LP per-record PCT/efficiency (snt01 mortality stays bit-exact, so the fix must be small-tree/dense-cohort specific). ROOT FOUND (instrumented FVS morts.f vs FVSjl at the 2005 cycle): T=753.09 (match), DIA0=5.916 (match), SDIMAX=359.08 (match), but **D10=6.594 (FVS) vs 6.57 (FVSjl)**. DIA0 matches ⇒ start DBH identical ⇒ the DIFFERENCE IS THE DIAMETER GROWTH G (FVS grows slightly more) → lower tn10 → more tokill → ~3 more LP killed. dbh_zeide=0 so all trees incl the LP are in D10. So s26's VARMRT difference is DOWNSTREAM of a sub-percent diameter-growth precision diff — the SAME CLASS as s5/s9. UNIFYING INSIGHT: the remaining growth residuals (s5/s9 + s26) are all small DG/HTG precision in non-snt01 stand/species calibrations (s26 stand = forest type 231 + planted LP), amplified downstream (board-foot / mortality). Fix both via the per-tree DG of forest-type-231 / small-tree growth vs FVS dgdriv (snt01 forest type 117 stays bit-exact). **COMPLETE ROOT (7-layer trace):** s26 stand WITHOUT the PLANT is BIT-EXACT ⇒ entire divergence is the established LP. Chain: LP self-thin 270→63.31 vs FVS 60.31 ⇒ mortality D10 6.57 vs 6.594 ⇒ DG ⇒ established LP DBH 2.5409 (FVSjl) vs 2.5609 (FVS) at SAME height 20.0 + IDENTICAL htdbh coeffs (SNALL[*,13]=243.86/4.2846/-0.4713). FVSjl establishment.jl:98 uses PLAIN _htdbh_dbh; FVS derives it in REGENT (regent.f:317-321, LHTDRG-calibrated ht-dbh) → +0.02 dbh. FIX = derive established dbh via the calibrated REGENT path (snt01 has no establishment ⇒ stays bit-exact). ✅ PARTIAL FIX LANDED: established dbh now adds 0.001*HK (regent.f:334) → LP dbh 2.5609 matches FVS (establishment.jl). REMAINING: the LP SMALL-TREE HEIGHT GROWTH 2005→2010 = FVSjl 20.0→24.907 vs FVS 20.0→24.967 (0.06ft/1.2% less) → lower dbh growth → D10 → ~3-tree mortality diff. Next = REGENT small-tree height increment for the LP (small_tree_growth.jl) vs FVS regent/htgf. ✅ 2 FIXES LANDED this session: (1) established dbh +0.001*HK (regent.f:334); (2) NPTIDS=points_inv-nonstockable (esplt2.f:74) → regen record count 44→50 matches FVS, RNG draws re-aligned. s26 diff cells 16→9; LP 2010 63.31→58.89 (FVS 60.31). REMAINING: a smaller RNG-state residual (LP↔stand interaction shifts the RANN draws when the planted cohort is present; the no-plant stand is bit-exact). Both fixes keep BARE/snt01 bit-exact (4487+28). PRECISE RNG RESIDUAL (instrumented FVS regent RAN vs FVSjl): both make exactly 50 LP RAN draws (count matches after NPTIDS fix), but FVSjl's values = FVS's OFFSET BY 2 (FVSjl[1]=FVS[3]=-0.61953 etc). So FVS makes 2 EXTRA NON-LP main-stream (RANN) draws before the LP small-tree growth that FVSjl doesn't — triggered by the planted cohort's presence (no-plant stand is bit-exact). Deterministic growth matches exactly. FIX = find the 2 missing RANN draws (likely the LP records' processing at the 2005 establishment-cycle REGENT(LSTART) call, esgent.f:48, which FVSjl's establish! doesn't replicate). LOCALIZED: FVS REGENT processes 3 SK(sp65) small trees BEFORE the LP(sp13) — FVSjl small_tree_growth! loops species ASCENDING (sp13 before sp65). The 2-draw RANN offset is the small-tree PROCESSING ORDER (and/or a 2-record diff) with the LP present. snt01 bit-exact (single dominant small-tree species / order coincides). FIX = match FVS's small-tree tree-processing order (IND1) so the per-record RANN draws align. Deep RNG-order alignment; the 2 committed fixes (dbh, NPTIDS) already match the deterministic growth + record count exactly. CORRECTION: FVS REGENT DOES loop ISPC ascending (regent.f:140) like FVSjl — the earlier 'sp65 before sp13' was a previous-cycle misread. The real 2-draw offset = FVSjl draws 2 EXTRA RANN in some species index < 13 (processed before the LP) within the 2005→2010 cycle ⇒ a 2-small-tree-record difference in another species' set (possibly tripling-related). Deterministic growth + LP record count match exactly. Next = diff the per-species small-tree record counts FVSjl vs FVS for sp<13. FINAL LOCALIZATION: small-tree species are only 13/22/33/65 (LP=13 processed FIRST, none below it), so the 2 extra RANN draws are NOT in the small-tree path — they're in the LARGE-tree/tripling path BEFORE small_tree_growth, triggered by the LP's presence (no-plant stand is bit-exact, so the original stand's RANN aligns WITHOUT the LP but draws 2 extra WITH it). FVSjl's LP RAN sequence = FVS's shifted by 2 (clean offset). Deterministic growth + record counts match exactly. The 2 draws are the deepest layer — a stochastic-DG/tripling interaction with the added cohort. DEEPEST LAYER: height_growth draws no RNG; the 2 extra main-stream draws are in the DGSCOR rejection loop (serial_correlation.jl:83, redraws bachlo until |frm|≤bound) during diameter_growth (which processes the LP before small_tree_growth overrides them). The LP's presence shifts the stand variance/ssig, flipping the |frm|≤bound rejection for ~2 borderline trees → ±2 redraws → the 2-draw offset. This is a stochastic rejection-boundary sensitivity from the added cohort — near the FP boundary. To close: instrument FVS dgscor.f redraw counts vs FVSjl for the s26 stand.
- [ ] **s32_volume** — VOLUME zeroes SCF (cols past 80) → all trees prod=01. Per-tree .trl: FVS scuft=0 below dbh 10, mcuft=0 below ~6; FVSjl leaks small-tree merch (+19 each). mrules.f fixes DON'T apply (sawDib=6 overshoots; merchL=10 BREAKS 3 bit-exact defaults) → SN R8 uses merchL=8/sawDib=7-9; multi-threshold. Default path bit-exact. CORRECTION (this session): the 'FVS scuft=0 below dbh 10' was ANOTHER bad .trl parse. FVS INCLUDES the small-tree sawtimber — zeroing scf for d<bfmin drops mcuft/scuft by ~570 (target diff is only +19), so that's wrong. s32 is a genuine ~0.7-1% MERCH-CUBIC PRECISION diff under the degenerate zeroed-SCF card (reliable .sum: mcuft 2765 vs 2746, scuft 1809 vs 1790), NOT a structural exclusion. The .trl column parsing is unreliable — use .sum or recompiled-FVS prints only.
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


## CONSOLIDATION (root-cause analysis, this session)
The three remaining non-ULP divergences reduce to TWO roots:

1. **DG PRECISION** (DGSCOR serial-correlation + WK3 calibrated-species COR evolution under a
   non-5 period) — drives BOTH:
   - **s5/s9**: the 10-yr board-foot/volume tail (first 10-yr cycle now bit-exact after the
     AUTCOR + d10n fixes; the tail is the COR evolution).
   - **s26**: PROVEN this session to be the SAME root, NOT a mortality bug. The 2005 dense-cohort
     mortality iterates 4x and matches FVS step-for-step (tn10 595.78/534.90/508.14/496.51 vs FVS
     595.52/534.84/508.30/496.83; d10 within 0.003). Every mortality intermediate (tt=753.09,
     dia0=5.916, sdimax, tn10, tokill, the line-reset) matches FVS to <0.07%. The ONLY seed is a
     ~0.03% iter-1 d10 difference, which (since the 2005 stand state is .sum-bit-exact) can only
     come from `g`=diam_growth. So s26 is the DG-precision root propagating through d10, amplified
     by the nonlinear self-thinning iteration + VARMRT distribution. RULED OUT as causes: SDI
     dbh-inclusion (DBHZEIDE=0 in BOTH, T=753.09 both), BA-percentile, all 90 shade_adj, the
     line-reset — all bit-match FVS.

2. **NVEL R8 CLARK TAPER GEOMETRY** (s32) — independent: emergent small-tree sawtimber-top floor.

Closing the DG-precision root (the YR-vs-FINT calibration split for DGSCOR/COR) would close BOTH
s5/s9 AND s26. Closing s32 needs the R8 taper sawtimber-boundary work.


## s26 DG-root REFINED (cohort split, this session)
Splitting the s26 2005 dense-cohort SDI by DBH (D<5 vs D>=5), inside-bark, vs an instrumented
morts.f confirms the stand STATE is bit-exact (n=105/524.41 TPA small, n=188/228.68 large — both
match FVS exactly) and the DG difference is CONCENTRATED in the D>=5 (large, original) cohort:
  small (D<5):  FVSjl DGib 0.577 / d10c 3.1969  vs FVS 0.5778 / 3.1973  (0.14% — matches)
  large (D>=5): FVSjl DGib 0.95  / d10c 11.8849 vs FVS 0.9626 / 11.8885 (~1% — the gap)
The large trees are the ORIGINAL (WK3-calibrated sp33/65) stand. So s26's seed is the COR-
calibration DG precision on the calibrated species — the SAME root as s5/s9 — NOT establishment
DG-init, NOT a competition table, NOT the LP cohort. (A prior split that omitted the bark
division falsely suggested an 8-16% DG gap; the bark-corrected split is the correct comparison.)


## s26 ROOT CORRECTED (all-cycle DG evidence — supersedes the 'COR root' claim above)
Dumping the D>=5 (large/original) cohort mean inside-bark DG for EVERY cycle, FVSjl vs an
instrumented morts.f, REFUTES the earlier 'same COR root as s5/s9' consolidation:
  cycle:   1990     1995     2000     2005      2010
  FVS:     1.23468  1.03142  0.99915  0.96263   0.97731
  FVSjl:   1.23468  1.03142  0.99914  0.95005   1.00811
The large-tree DG is BIT-EXACT through 2000 and diverges ~1.3% EXACTLY at 2005 — the cycle the
LP cohort is established. A COR-calibration root would diverge from cycle 1; this does not. So
s26's seed is the ESTABLISHMENT cohort's effect on the EXISTING trees' diameter-growth COMPETITION
(FVSjl OVER-suppresses the large-tree DG when the 270-TPA small LP are present — consistent with
including the small LP crowns in CCF/RELDEN where FVS effectively does not, OR an establishment-vs-
DG ordering difference). This is s26-SPECIFIC (establishment x DG-competition), DISTINCT from the
s5/s9 COR/AUTCOR root. Net: the three divergences are THREE roots, not two — correcting the prior
consolidation. Next probe: the DG competition terms (CCF, RELDEN, PCT, BAL) for a large tree at
the 2005 cycle with vs without the LP present, to pin which term over-counts the small cohort.


## s26 ROOT RE-CORRECTED (DDS-matches evidence is conclusive)
The 'establishment x DG-COMPETITION' claim above was based on a measurement artifact and is WRONG.
Re-measuring the FVS DDS *after* its forest-type terms (it adds FTUPHD*KUPHD at dgf.f:321-328, which
my first dump captured before) shows the FULL deterministic ln-DDS MATCHES:
   FVSjl meanDDS = 2.54193   vs   FVS meanDDS = 2.54182   (0.004%, D>=5 cohort, 2005)
Every deterministic DG input for the s26 large trees at 2005 matches FVS:
   - calibration  CONSPP (dg_const+dg_cor): 0.33673 both
   - competition  BA 171.875 both; AVH 64.39 both; PCT 49.5596 both; PBA 193.89 vs 194.36 (0.24%);
                  PBAL 96.62 vs 96.81 (0.19%, and PNTBL=-0.004 so its DDS effect is ~1e-3, negligible)
   - forest type  IFORTP 520 -> KUPHD=1 (FVSjl ftgrp=uphd) both
   - FULL ln-DDS  2.54193 vs 2.54182
So s26 is NOT a competition/forest-type/calibration bug. With DDS matching, the residual is the
STOCHASTIC DGSCOR (frm=exp(frmbase+corr*oldrn) serial-correlation: corr/oldrn) — the SAME family as
s5/s9 — surfacing when the LP cohort is added. The stand-level impact is sub-0.5% (.sum 2010: TPA 401
vs 403, BA 151/151, cuft 3294 vs 3292 = 0.06%), i.e. ULP-FP-level accumulation. (The apparent '1.3%
DGib' at the mortality point is inconsistent with the 0.06% cuft and is a post-DDS/measurement
nuance, not a deterministic-input gap.) NET: s26 folds back to the DGSCOR-precision family with s5/s9
— deterministic path fully excluded. Lesson logged: always measure Fortran intermediates at the FINAL
assignment point (FT terms are appended after the base DDS).


## s26 FINAL MECHANISM (the stochastic path identified)
At cycle 2005 tripling is long finished, so the large-tree frm comes from `dgscor!`
(diameter_growth.jl:632) — the STOCHASTIC serial-correlation rejection-loop RNG draw
(frm=exp(frmbase+corr*oldrn) with a bounded bachlo redraw). snt01 runs this identical
`dgscor!` path every cycle and is BIT-EXACT, so the kernel + main-RNG are correctly aligned
there. s26 differs because the ESTABLISHMENT consumes main-RNG draws (the LP crown-ratio
bachlo draws, establishment.jl:133-143) that OFFSET the subsequent `dgscor!` draw sequence for
the existing trees. Net: s26's residual is an RNG-DRAW-ALIGNMENT difference introduced by the
establishment, surfacing in the bounded dgscor! draws — sub-0.5% at the stand level because the
draws are bounded/averaged. HONESTLY: this is NOT ULP-FP (it is a different RNG sequence, just
small). The fix is to match FVS's exact main-RNG consumption order through ESTAB crown assignment
so the dgscor! sequence stays phase-locked (cf. the establishment ESRANN/main-RNG split already
done for heights). Scoped to establishment.jl crown draws x dgscor! ordering.

## HONEST CLASSIFICATION OF ALL THREE (re the drop-in spec)
None of the three remaining residuals is honestly 'ULP FP':
  - s5/s9 : COR/AUTCOR calibration evolution under a non-5 period (1st 10-yr cycle now bit-exact).
  - s26   : establishment main-RNG offset of the dgscor! serial-correlation draw sequence.
  - s32   : NVEL R8 Clark taper geometry at the sawtimber-top boundary (~0.7%).
They are small (sub-1%) and precisely scoped, with all logic/coefficient/table causes excluded,
but they are GENUINE divergences requiring real fixes (RNG phase-locking, the COR YR/FINT split,
the taper boundary) — NOT ULP-FP. The drop-in spec is not yet met; claiming ULP-FP for these
would be inaccurate.


## s26 PINPOINTED to a discrete +6 RNG-draw offset (RANN-count instrumentation)
Counting main-RNG (RANN) draws per cycle, FVSjl vs an instrumented rann.f, the streams are
PHASE-LOCKED bit-for-bit THROUGH 2005 (cum count 60/120/309/1710 identical at every cycle's dgf
entry — so the establishment draws matched exactly), then FVSjl makes EXACTLY 6 MORE draws during
2005->2010 (3267 vs FVS 3261). Per-phase counting (gated on 2005) localizes all of it to the GROWTH:
   after dgf=2811, after small_tree_growth=3267 (+456); htg/mort/esuckr/estab/crown add 0.
So the +6 is in the 2005->2010 dgf (dgscor!) + small-tree (REGENT) growth. The LP cohort sits right
at the dbh~=3 small/large boundary (2.6->3.2 over the cycle), so 1-2 LP trees are classified small-
vs-large differently (REGENT vs dgf path) — or a dgscor!/REGENT rejection redraw fires once — costing
the 6 extra draws and offsetting the remainder of the stream. NET: s26 is a DISCRETE RNG-count
divergence at the dbh-3 growth boundary, NOT a precision residual. NEXT: dump the LP dbh near 3.0 at
the 2005 growth on both sides; if a tree is within ~1 ULP of dbh 3.0 and classified oppositely, the
offset is ULP-rooted (acceptable per spec); else it is a small/large-classification or REGENT-draw-
count logic gap to align. This is the most actionable of the three and is well-scoped.


## s26 +6 SPLIT across the growth boundary (DGDRIV vs REGENT)
Instrumenting the RANN counter at FVS's growth-call boundaries (grincr.f DGDRIV@437 / REGENT@449)
vs FVSjl's phase counts at 2005:
                 entry   afterDGF/DGDRIV   afterSMALL/REGENT
   FVSjl :       1710        2811               3267
   FVS   :       1710        2802               3261
   delta :         0          +9                 +6 (i.e. -3 within REGENT)
So the +6 net = +9 in dgf (large-tree DGDRIV: dgscor + new-tree oldrn-init) and -3 in REGENT
(small-tree). The -3/+3 is ~3 trees shifting across the dbh-3 small/large boundary (FVSjl classifies
~3 as large where FVS keeps them small, so they leave REGENT and enter dgf's dgscor); the remaining
+6 in dgf is extra serial-correlation/oldrn-init draws (the new LP cohort's oldrn initialization or
dgscor rejection redraws). NEXT: split FVSjl's dgf draws into oldrn-init vs dgscor! to attribute the
+6, and dump the ~3 boundary trees' dbh (within ULP of 3.0 => ULP-rooted classification flip; else a
small/large-threshold comparison difference to align with REGENT XMAX handling).


## s26 ROOT MECHANISM (dgscor! rejection depends on per-tree oldrn)
dgscor! (diameter_growth.jl:632) rejection loop: `frm = bachlo*rhocp + rho*oldrn[it]; reject if
abs(frm) > bound`. CRUCIAL: the reject test depends on oldrn[it] — the TREE's persisted serial-
correlation residual — not only on the drawn bachlo value. So the number of redraws (and thus the
RANN draw count) is a function of WHICH tree sits at each draw position. Consequences for s26's +9
in dgf at 2005:
  - bachlo is bit-exact (snt01) and the stream is phase-locked at the 2005 dgf entry (1710 both),
    so a pure value-only rejection could NOT diverge.
  - It diverges because, AFTER the LP cohort is inserted, FVSjl's per-species tree PROCESSING ORDER
    (species_sort! / IND1) and/or the new LP's oldrn INITIALIZATION differ from FVS's, so a different
    tree (different oldrn) lands at some draw positions -> a few extra rejection redraws -> +9.
  - The -3 in REGENT is the ~3-tree dbh-3 small/large split (operator d>=3 matches FVS; so it is the
    dbh values of a few near-3 trees, or the same order/oldrn perturbation feeding which trees REGENT
    sees). 
FIX DIRECTION: phase-lock the post-establishment tree order (IND1 after regen insertion) and the LP
oldrn init to FVS, so dgscor!'s oldrn-dependent rejection sequence stays aligned. This is the precise,
well-scoped root; it is an RNG-ORDER alignment after regen, not ULP and not a coefficient error.
NEXT PROBE: dump FVSjl species_sort! order vs FVS IND1 for the 2005 stand (and the LP oldrn values)
to confirm the order/oldrn mismatch and align it.


## s26 ROOT CONFIRMED — regen trees never get oldrn initialized (LP oldrn=0)
Direct dump at the 2005 dgscor!: the LP cohort enters with oldrn=0.0 (dbh 2.561, all 0). FVSjl
initializes the serial-correlation residual (old_random) ONLY ONCE, in calibrate_diameter_growth!
called at setup (simulate.jl:42). Trees ESTABLISHED mid-run (the LP, added 2005) are therefore never
oldrn-initialized and enter dgscor! with oldrn=0. FVS re-initializes new-tree OLDRN each cycle in
DGDRIV (dgdriv.f:402-433/534-605: regression OLDRN=BNY+(EDDS-BNX)*SLOP, or the <5-GST bachlo path).
Because dgscor!'s rejection (abs(bachlo*rhocp + rho*oldrn) > bound) DEPENDS on oldrn, the LP's wrong
oldrn=0 shifts the rejection sequence -> the +9 dgf draws. snt01 has NO regen (all trees present at
setup, inited once), so it is bit-exact and never exercises this -- s26 is the first regen scenario in
the suite to hit it. THIS IS A REAL, REGEN-SPECIFIC BUG, not ULP.
FIX: initialize new (established) trees' old_random the way FVS does -- the deterministic regression
path (calibrate's bny+(exp(wk2)-bnx)*slop, already implemented at diameter_growth.jl:444) applied to
regen trees when they first enter growth, leaving existing trees' persisted oldrn untouched. Must
verify (a) it drives the s26 dgf draw count to FVS's (closing the +9) and (b) snt01/BARE stay
bit-exact (no oldrn=0 large trees there). The -3 REGENT piece likely follows once the order/oldrn is
aligned. This is the concrete fix for s26 and is well-scoped.


## s26 ROOT — CORRECTED PATH (LP species is UNCALIBRATED -> bachlo init, which DRAWS)
Refinement: sp13 (loblolly) is NOT present at setup (the plant is at 2005), so calibrate_diameter_
growth! never sees it and it has NO growth samples -> calibrated[13]=false. So FVS's per-cycle new-
tree oldrn init takes the <5-GST BACHLO path (dgdriv.f:599-605: Z=BACHLO(0,SIGMA); reject Z>DGSD*SIGMA;
OLDRN(I)=Z), which CONSUMES random draws for each new LP record. FVSjl runs its oldrn init ONLY once at
setup (simulate.jl:42), so the LP are never inited (oldrn=0) and FVSjl skips those bachlo draws. The
net 2005->2010 bookkeeping (+9 dgf / -3 REGENT) is the combination of (a) FVSjl skipping FVS's LP
oldrn-init bachlo draws and (b) the resulting oldrn=0 changing the dgscor! rejection sequence.
THE FIX (concrete, scoped): replicate FVS DGDRIV's per-cycle new-tree (oldrn==0) oldrn initialization
for trees ESTABLISHED mid-run — for an uncalibrated species that is the bachlo rejection path
(bachlo(0,sigma); reject z>DGSD*sigma; oldrn=z), drawing in the SAME order FVS does (species-sorted,
after the calibration but before the per-tree dgscor!). For a calibrated regen species it is the
deterministic regression path. Must match FVS's draw count/order exactly and keep snt01/BARE bit-exact
(they have no mid-run regen, so the once-at-setup init still covers them). This is the precise close
for s26's RNG offset.
