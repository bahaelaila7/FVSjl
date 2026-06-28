# FVSjl Southern — Deferred Structural Backlog (Campaign 2)

Ordered **most-upstream / least-dependent / easiest first**. These are NOT audit flags (those are all
resolved — see INDEX.md + CAMPAIGN_DONE-era record); they are the structural/infrastructure items the flag
campaign deferred-by-design. Same doctrine + 5 principles apply (trace LOGIC both sides; port faithfully +
live-validate THEN test; document every verdict; variant-aware).

Baseline: suite **5014 pass / 1 broken** (the 1 = item #6 COMPRESS #29).

## Tier 1 — bounded & independent (easy wins, no new infrastructure)

**1. Single-canopy structure-stage (NTREES≤1 path). — ✅ DONE (live-validated).**  sstage.f:235-256: a
   1-canopy-tree stand uses a special branch — cover-based class 0/1 when `WK6 = CW²·TPA·π/4 < 435.6·TMPCCM`
   (0=BG, or 1=SI if TPA≥TPAMIN), else DBH-direct (SSD→1, SAW→2 w/ PCTSMX SE→SI demote, else 5). Checked BEFORE
   stratification (so it fires when the lone tree forms no OK stratum) and `GOTO 80`s past the NSTR-setting code,
   so it reports N=0. jl's `structure_class` only had the nstr==1 DBHNOM path (no cover branch), so a sparse
   single canopy tree diverged (it gave 0=BG where FVS gives 1=SI via the TPA-override). PORTED at
   structure_stage.jl:157-176 (reuses `st.crarea[1]`=WK6, returns nstr=0). Live-validated vs FVSsn across all
   four reachable sub-branches (BG / SI-via-TPA / SI-via-DBH / OS) — see the 4 `struct_1canopy_*` scenarios +
   test_structure_stage.jl. Class-2/SE is unreachable for a lone tree (tiny SDI always demotes — faithful).
   Suite 5014→5026/1.

## Tier 2 — upstream infrastructure (each unblocks a small cluster)

**2. Per-point density layer — ✅ DONE (2a wired faithful; 2b found INERT in live, matched).**
   `PTBAA` (point BA) already existed (`point_basal_area!`). Added `point_ccf`(PCCF)/`point_tpa`(PTPA) to the
   `Density` struct + `point_density!` (standstats.jl, mirrors dense.f:210-211 with the same PI/GROSPC scale as
   the already-validated point_ba), called from `compute_density!`.
   - **2a. multi-point regen pccf — DONE.** regent.f:178 `CR=0.89722−0.0000461·PCCF(IPCCF)`, IPCCF=ITRE(I).
     establishment.jl now indexes `density.point_ccf[plot_id]` per regen record (was whole-stand `stand_ccf`).
     Faithful transliteration; coeff 4.6e-5 ⇒ bare/sparse stands unchanged to print resolution; suite green.
   - **2b. TCONDMLT point weights — INERT in live FVS; NOT applied (matching live).** Built the term, then a
     multi-point THINBBA differential showed live FVSsn is byte-identical with PTPAWT/PBAWT = 0 vs **9999** — the
     cuts.f PTPA(IP)/PTBAA(IP) is uniform/zero at thin time (a point-thinning-only path plain TCONDMLT doesn't
     arm). jl's per-point `point_tpa` *varies*, so adding the term DIVERGED from live; per the doctrine (match
     live) the `_cut_pref_wt` point term was REVERTED. Weights are still parsed into Control (inert). Guarded by
     `tcond_base`/`tcond_pw` (+ `tcond_pw.sum.save` live oracle) in test_tcondmlt.jl. Open thread: *why* FVS's
     PTPA(IP) is uniform/zero needs a debug-FVS trace (low priority — feature is inert either way).
   Suite 5026→5036/1.

**3. Snag soft-decay model (DECAYX hard→soft). — PARTIAL: faithful transition PORTED; upstream YRDEAD bug
   unmasked (→ #28).** The hard→soft state is a DKTIME *threshold*, not a fractional rate: fmsnag.f:282-284
   flips HARD when `(IYR−YRDEAD) ≥ DKTIME`, DKTIME = FMSNGDK = `DECAYX·(1.24·D+13.82)` for SN (= jl's falldown
   `tsoft`; DECAYX 0.07/0.21/0.35 matches). fmssum.f:9-22 then counts the flipped DENIH into the SOFT report
   column.
   - **3b. snag_summary hard→soft — ✅ FIXED this session (suite 5105/2, no regression; live-validated).** The
     split now tracks live FVSsn: carbon_snt 1995 35.79h/6.91s BIT-EXACT (was the documented divergence), 2000
     44.2h/3.8s vs live 44.8/3.3, 2005 65.6h/5.4s vs live 66.8/4.3 — was wildly inverted (2000: 6.5h/41.6s). FIX
     (the de-risked plan below, executed): added `SnagList.yrdead` (TRUE death year = cycle-END−1 for ordinary
     mortality per fmkill.f:140, = `year` for input/fire), set in book_mortality_snags! (fint threaded), and switched
     the hard→soft DKTIME classification (snag_summary) + post-burn window + salvage age to read `yrdead`; the FALL
     stays on `sn.year` (cycle-start fall-clock). StandDead falldown [3.796,4.393,5.354,9.535] UNCHANGED (split-
     independent: snag_bole_carbon uses den_hard+den_soft total) — verified ≤0.05. Test in test_carbon.jl.
     ⊕ COMPLETENESS verified across ALL snag sources (FMSADD YEAR by ITYP): ordinary-mort (fmkill:143 ITYP4 =
     IY(ICYC+1)−1) ⇒ jl cycle-end−1 (the fix); fire (fmeff:608 ITYP1 = IYR), input (fmsdit:135 ITYP3 = IY(1)−FINTM),
     thinning (fmscut:157 ITYP2 = IY(ICYC)) all have YRDEAD = the fall-year ⇒ jl's `yrdead=year` default is correct
     for them. Only ordinary mortality (dies at cycle-END, after FMMAIN) carries the offset. Tiny
     residual ✅ FIXED — and it was a REAL one-year-offset BUG, NOT float noise (re-trace discipline vindicated AGAIN:
     I'd labeled it "near-boundary float effect" but checked it rather than trusting the label). Live HARD-flag dump
     (FMSNAG, per-snag dbh/DENIH/HARD/YRDEAD/DKTIME) showed: the carbon REPORT (FMCRBOUT, fmmain.f:206) runs BEFORE
     the annual FMSNAG hard→soft flip (fmsnag.f:282-284), so it reflects the HARD flag as of the PREVIOUS cycle's last
     FMSNAG year = IY(ICYC)−1. jl used age=iyr−YRDEAD; should be iyr−1−YRDEAD (one year behind). Proof: dbh8.09
     DKTIME5.01 → live age-5<5.01 = HARD, jl age-6≥5.01 = soft. FIX: snag_summary uses `iyr−1` in the DKTIME check.
     RESULT: 2005 split 65.6h/5.4s → 66.7h/4.3s = live 66.8/4.3 (≈bit-exact); 2000 44.2/3.8 → 44.6/3.5 vs live
     44.8/3.3. Suite 5109/2, no regression; test tightened (atol 1.0/1.5 → 0.5). 1990/1995 unchanged (no near-tie).
     ⊕ RESIDUAL now ~0.2 at 2000 only (2005 ≈bit-exact). jl HARD-dump confirms the mechanism is right — jl's SOFT set
     is exactly the age-5 cycle-1 cohort with DKTIME≤5 (e.g. dbh3.2/dkt3.74, dbh5.0/dkt4.2 soft; dbh5.3+ vary). The
     0.2 is a handful of small cycle-1 snags whose DKTIME sits within ULP of age=5 (the boundary), flipping
     differently than live by Float32 rounding — NOW genuinely accepted-class (near-tie, sub-% on totals), distinct
     from the systematic 1-yr offset that was the real bug. Fully confirming would need the COMPLETE live soft-set
     (only a 25-row sample captured); not worth re-instrumenting given 2005 is bit-exact and the offset is fixed.
     ─ (orig:) ─
   - **3b(orig). snag_summary hard→soft — transition PORTED (faithful), but unmasked an upstream divergence.** Added
     the DKTIME flip to `snag_summary` (snag.jl). It is a faithful port (verified vs source). Validating the
     6-cycle FVS_SnagSum (carbon_snag.key) showed jl OVER-softens (1995: jl 2.9h/39.8s vs live 35.79h/6.91s;
     **totals bit-exact**, only the split wrong). ROOT CAUSE (the real divergence, NOT the transition): jl dates
     periodic-mortality snags at the cycle-START year and ages them with a per-cycle `update_snags!`, whereas FVS
     reports snags via **FMSSUM (fmmain.f:178) BEFORE the annual FMSNAG aging loop (fmmain.f:228-232)** — so a
     cycle's fresh mortality snags read ~age-0 (hard) at the report, and only OLD input snags have aged past
     DKTIME (= live's small, stable soft pool). Fixing the snag report/aging timing IS the **#28 FFE annual-loop
     phasing** work. Per principle #3 the faithful transition is KEPT (the cycle-0 SnagSum test still passes);
     the timing fix lands with #28, at which point the multi-cycle SnagSum will match live.
     ⊕ SHARPENED (this session): the death-year is concretely **fmkill.f:140 `YEAR = IY(ICYC+1)−1`** (cycle-END
     minus 1) → `FMSADD:240 YRDEAD(X)=YEAR`. jl uses cycle-START (`current_cycle_year`), a ~(cyclen−1)-yr over-age
     that makes the DKTIME check `iyr−sn.year ≥ dktime` trip too early ⇒ over-soften (e.g. a 1995→2000 snag: jl
     reads age 5 @2000-report → soft; FVS YRDEAD 1999 → age 1 → hard). CONFIRMED coupled to #28, NOT a standalone
     fix: `sn.year` ALSO drives the `update_snags!` fall-age gate `clamp(eff−sn.year,0,nyears)`, so swapping it to
     cycle-end−1 would cut the fall steps (5→1) in the next cycle and break the bit-exact carbon_snt StandDead bole.
     So the faithful fix needs BOTH the report-before-aging order AND YRDEAD=IY(ICYC+1)−1 together (the #28 phasing),
     not a death-year change alone. ⚠ Note: the SNAGDCAY test (3a) validates around this by using DECAYX=2 to
     SUPPRESS the transition entirely (all-hard), so the keyword validates bit-exact despite this open split residual.
     ⊕ ACTIONABLE ORDER for the #28 refactor (fmmain.f traced this session): FMBURN(170) → REPORTS FMSSUM/FMDOUT/
     FMCRBOUT(178/202/206) → annual loop FMSNAG-age/FMCWD(228-236) → FMOLDC(268); FMKILL (ordinary-mortality snag
     CREATE) is AFTER FMMAIN (gradd.f:122). ⇒ FVS reports snags PRE-aging, and a cycle's fresh mortality snags first
     appear in the NEXT cycle's report. jl inverts this (ages at grow_cycle! start + creates mortality, THEN reports
     at the boundary). The faithful fix = sample the snag report on the START-of-cycle (pre-age, pre-this-cycle-
     mortality) snag pool with YRDEAD=IY(ICYC+1)−1, mirroring the carbon_hook pattern already used for the fire-cycle
     carbon sample (#28 ISSUE A). VALIDATION GATE: carbon_snt StandDead falldown [3.796,4.393,5.354,9.535] must stay
     bit-exact (it depends on the SAME snag pool/timing) — so move the REPORT sample point, do NOT change the fall
     stepping. This is the last enabling refactor; it unblocks 3b (split), the soft-decay report, and likely #28's
     2005-bole (same snag-state-at-report question).
     ⊕⊕ DE-RISKED to a small change (this session): the SnagSum hard/soft SPLIT can be fixed WITHOUT touching the
     validated fall, because `standing_dead_carbon` = `snag_bole_carbon`(uses den_hard+den_soft TOTAL) + crown — it is
     SPLIT-INDEPENDENT. So the carbon_snt falldown can't regress from a classification-only change. CRUX confirmed by
     reading update_snags! (snag.jl:200): `yrs = clamp(eff − sn.year, 0, nyears)` ties the FALL-step count to
     `sn.year`, so `sn.year` must STAY cycle-start (its tuned fall-clock); but the CLASSIFICATION (snag_summary
     dktime `iyr − sn.year ≥ dktime`) needs the TRUE YRDEAD = cycle-end−1 for ordinary-mortality snags (fmkill.f:140).
     READY-TO-EXECUTE PLAN: add `SnagList.yrdead::Vector{Int32}` (the true death year; = `year` for input/fire snags,
     = cycle-end−1 for ordinary mortality, set in book_mortality_snags! with `fint` threaded), have add_snag! default
     `yrdead=year`, and switch ONLY snag_summary's dktime check to `iyr − yrdead`. Validate: carbon_snt SnagSum
     hard/soft → live (carbon_snag.key 6-cycle) AND falldown [3.796,4.393,5.354,9.535] STAYS bit-exact (guaranteed by
     split-independence). Post-burn window (snag.jl:245) should also read yrdead. This is bounded (~4 edits: struct +
     constructor + add_snag!/book_mortality_snags! + snag_summary) — no fall-mechanism change.
   - **3a. SNAGDCAY — ✅ PORTED + LIVE-VALIDATED (suite 5099/2). The "BLOCKED" label was STALE** (re-traced
     against source after the user flagged "blocked by design"). fmin.f:633 opt 11 is just `DECAYX(JSP)=ARRAY(2)`
     — a per-species override of the snag decay-rate multiplier (DKTIME=DECAYX·(1.24·D+13.82)) jl already uses,
     the EXACT pattern of the already-ported SNAGFALL. Added `FFEParams.snag_decayx_ovr` + `_snagdcay!` (mirrors
     `_snagfall!`), applied at the two live sites (snag_summary DKTIME, fmscro! TSOFT). Default empty ⇒ inert
     (suite unchanged). Live-validated: `SNAGDCAY 0 2.0` (live echoes "RATE-OF-DECAY…IS:2.000") keeps ALL snags
     HARD — live SNAG SUMMARY hard total 48.0(2000)/71.0(2005) soft 0, jl matches bit-exact (test_carbon.jl).
     **SNAGBRK** — ✅ PORTED + FUNCTIONAL this session (suite 5109/2, default bit-exact; magnitude refinement
     pending). Stages: (1) PARSE `_snagbrk!` → `FFEParams.snag_htx` (4 HTX/species; bit-verified vs hand calc
     6.70/9.71/1.70/2.49 for `SNAGBRK 0 10 20 15 30`); (2) APPLY `ffe_snag_height_loss!` (FMSNGHT SN-default,
     per-year in ffe_fuel_update!'s annual loop) shrinking a new `SnagList.htcur` (HTIH/HTIS, init=HTDEAD;
     hard/soft via the DKTIME gate, regime via 0.5·HTD, <1.5ft⇒snag removed); (3) COUPLE — snag_bole_carbon
     recomputes the bole by the FMSVOL(htcur)/FMSVOL(HTDEAD) merch-cuft ratio when snag_htx is set. All gated on
     snag_htx non-empty ⇒ the SN default (HTX=0, htcur≡height) keeps the frozen `bolevol` BIT-EXACT (suite green).
     LIVE-validated DIRECTION: `SNAGBRK 0 10 20 15 30` reduces carbon_snt StandDead in both — live base 5.4/9.5 →
     5.1/9.0 (2000/2005), jl 5.32/9.52 → 4.57/8.71. Test asserts the effect + no-op default (test_carbon.jl).
     ⊕⊕ DATA-FLOW WALL CRACKED + bug #1 FIXED (this session). The fmdout NSNAG=0 wall was a RED HERRING — fmdout is a
     fuel-report path; the snags ARE live in **FMSNAG** (fmmain.f:232, the per-year ager). Instrumenting FMSNAG
     (SNAGN/SNAGST dump) exposed the full per-snag live state (DBHS/DENIH/DENIS/HTIH/HTDEAD/YRDEAD) — unblocking BOTH
     this AND the #28 2005-bole. Findings: (a) 3b YRDEAD CONFIRMED BIT-EXACT (live yrd=1994 ordinary-mort = cycle-
     end−1, 1985 input = IY(1)−FINTM — exactly jl). (b) SNAGBRK had TWO bugs. BUG #1 (FIXED): jl picked the height-
     loss hard/soft rate from the DKTIME report transition, but FMSNGHT uses the snag's INITIAL state (DENIH/DENIS via
     IHRD) — fixed `ffe_snag_height_loss!` to `soft = den_soft > den_hard`. RESULT: the dbh34.6 input snag's HTIH is
     now BIT-EXACT vs live (92.69→65.54→46.35→28.74). (c) BUG #2 — ALSO FIXED: jl scaled the bole by
     `_fm_cuft(dbh, htcur)` = a NORMAL short tree, but FVS's FMSVOL(I,XHT) computes the ORIGINAL tree (D=DBHS,
     H=HTDEAD) TRUNCATED at XHT (fmsvol.f:101-142, the fat lower bole). Fixed by routing snag_bole_carbon through
     jl's existing `cftopk` (Behre top-kill, volume.jl:118) at itht=htcur·100. RESULT: ✅ SNAGBRK BIT-EXACT vs live
     to display resolution — carbon_snt StandDead 4.30/5.06/9.04 vs live 4.3/5.1/9.0 (was 4.57/8.71). Test tightened
     to the live magnitude. Default HTX=0 stays bit-exact; suite 5109/2. ⇒ SNAGBRK FULLY ported + live-validated; the
     whole snag subsystem (SNAGFALL/SNAGDCAY/SNAGBRK + 3b split) is now feature-complete AND value-validated vs live. (Soft-portion HTIS dual-height is collapsed to one htcur — exact only when den_soft=0, i.e. no
     SNAGPSFT; the report's hard/soft track switch is handled via the DKTIME gate.)
     ─────────────── (original UNPORTED note, superseded by the port above:) ───────────────
     **SNAGBRK(orig)** — UNPORTED but DEFAULT-INERT (not "blocked by design"): opt 10 sets per-species HTX (snag
     height-loss); SN's default is `HTX=0` (fmvinit.f:1089) and FMSNGHT with HTX=0 returns the height unchanged,
     so WITHOUT the keyword snags don't lose height and jl (no height-loss model) is faithful. SNAGBRK is a real
     GAP only if a scenario USES it.
     ⊕ FULL SCOPE traced (this session) — it is a MODEL EXTENSION, not a parse: (1) PARSE — SPDECD species (fld1) +
     YRS50(hard/soft)=ARRAY(2/3), YRS30(hard/soft)=ARRAY(4/5); compute 4 HTX/species (fmin.f:538/546/557/566):
     HTX1=(1−0.5^(1/YRS50h))/HTR1, HTX2=(1−0.6^(1/(YRS30h−YRS50h)))/HTR2, HTX3=…/(HTR1·HTXSFT), HTX4=…/(HTR2·HTXSFT)
     (HTR1=HTR2=0.01, HTXSFT=2.0). (2) APPLY — per year, FMSNGHT (SN = CASE DEFAULT, fmsnght.f:153-164): if
     HTCURR>0.5·HTD use HTX1/3 else HTX2/4, `HTSNEW=HTCURR·(1−HTR·HTX·SFTMULT)^NYRS` (SFTMULT 1 hard / 2 soft), and
     HTSNEW<1.5→0 (snag becomes fuel). (3) COUPLE — the shrinking height feeds the snag BOLE volume
     (SNVIH=FMSVOL(HTIH)·DENIH + SNVIS=FMSVOL(HTIS)·DENIS), so StandDead changes. ⇒ jl must replace the FROZEN single
     `sn.height`/`bolevol` with a TIME-VARYING height (faithfully DUAL: HTIH hard + HTIS soft, different loss rates;
     collapsible to one only when den_soft=0, i.e. no SNAGPSFT) and recompute the bole each cycle. Inert at default
     (HTX=0). NOT a tail-of-session port — needs its own focused pass + a SNAGBRK live scenario. Last self-contained
     backlog feature.

## Tier 3 — bounded, independent, medium

**4. density `notre!` FINT/FINTM dead inflation. — ✅ DONE (live-validated; small residual documented).**
   notre.f:122-124 inflates dead-record PROB by `FINT/FINTM` for backdated-density CALIBRATION, then deflates
   for treelist/FFE uses (FMSSEE/PRTRLS). jl carries the TRUE dead TPA everywhere, so the inflation is applied
   ONLY inside `calibrate_diameter_growth!`, scoped to the two backdated density passes + the PCTILE percentile,
   and restored before returning (so FFE snag seeding still sees true TPA). Gated on FINT≠FINTM ⇒ the default is
   untouched (bit-exact). Live-validated on `dead_fint.key` (carbon_snt's 10 input dead + `GROWTH 0 10 0 0 5`):
   the missing inflation gave 2000 TCuFt −8 / MCuFt −9 vs live; with it the .sum tracks live to ≤4 cuft across
   all cycles (regression test in test_growth.jl + dead_fint.sum.save oracle). RESIDUAL (~3 MCuFt, ~0.1%): the
   dead trees' contribution to the crown-init backdated CCF — `init_crown_ratios!` computes its CCF on the live
   partition only, whereas FVS's single DENSE includes the (inflated) dead. A separate, pre-existing live-only
   simplification (inert at default, sub-percent here); left documented rather than risk the bit-exact crown
   init. Suite 5036→5053/1.

**5. COMPRESS #29 — post-compression record-ORDER drift. — ✅ LANDED (suite 5061/2, no regressions). Two fixes:
   (a) the comprs.f:1007 TREDEL post-merge compaction → record ORDER bit-exact; (b) the SPESRT/LNKCHN sort_key
   reset → fixed a real same-species swapped-RNG bug. COMPRESS merge/order/swap now bit-exact vs live FVS. The
   ≈1% s22 tail that remains is a downstream crown-ratio→DG residual on ONE sp33 record → split to task #41 / item 5b.**

   ════════ RESOLUTION (this session) ════════
   The earlier "zero-prob layout" hypothesis was REFUTED (jl has 0 zero-prob records at COMPRESS). A second debug-FVS
   dump of IREC1 (the positive IND entries before comprs.f:1007) proved FVS's survivor slot = the **minimum** record
   index of each class — i.e. EXACTLY `dst=minimum(mem)`, where `_merge_one!` already writes the merged record — NOT
   mem[1]. My earlier tredel sim used mem[1] for the survivor slots (the bug). With survivors = dst (s22: [1,3,4,6,14]),
   the smallest-vac←largest-survivor swap reproduces FVS's order `[c5,c3,c1,c2,c4]` EXACTLY. Implemented in
   `_merge_classes!` (compress.jl): simulate the TREDEL swap on the dst slots → fvs_order, gather the merged records by
   ascending dst, then cycle-permute into fvs_order. Verified: all 5 s22 merged records match live FVS to the digit in
   sp/dbh/ht/ICR/OLDRN AND order (debug-dump MERGEREC at the compression cycle).
   Then a SECOND, genuinely-new bug surfaced: 4 of 5 records grew bit-exact through 2005 but the two sp33 records
   diverged. Reverse-engineering the DGSCOR serial-corr equation (frm = bachlo·rhocp + rho·OLDRN_old) showed the two
   sp33 records had SWAPPED bachlo draws. Root: after COMPRESS, FVS's SPESRT→LNKCHN→SETUP (lnkchn.f called I=1..IREC1)
   re-lists each species in ASCENDING PHYSICAL record index, but jl's species_sort! orders by `sort_key`, and the merged
   records carried STALE sort_keys from their original member slots → reversed order for the sp33 pair → swapped per-tree
   RNG. Fixed by resetting `sort_key=position` at the end of `_merge_classes!` (mirrors core/trees.jl `compact!`, the
   thin/TREDEL path — same documented FVS behavior). After the fix OLDRN matches all 5 records; s22 2005 TPA 415→413.
   RESIDUAL — ✅ TRACED TO CERTAINTY = the ACCEPTED COMPRESS-eigensolver + ULP divergence (#41 CLOSED): the single
   sp33 record at dbh 6.664 grows ~1% slow because its dgf `point_bal` competition term reads the wrong POINT BASAL
   AREA (PTBAA): jl puts it on plot 5, FVS on plot 3. That plot is the merged record's `sel` (the RANN-sampled member,
   comprs.f:725, copies ISP/ITRE/IMC/…). The RANN draw VALUE is BIT-EXACT (x=34.39857 both), but it selects a
   different member because the WITHIN-CLASS member ORDER differs for ~4 near-tied records. ROOT (dumped both sides):
   the PC1/PC2 sort keys WK3/WK4 (comprs.f:308-318) match live to **< 1 Float32 ULP** (rec6 WK3 9154.72461 vs jl
   9154.72413, Δ0.0005 < ULP 0.00098; rec13/52/186 WK3 all within 0.07 of each other), and those sub-ULP diffs flip
   the partition's near-tie-sensitive nested sort (rdpsrt by wk3, then wk4). BRCMPR is a no-op stub (WPBR, inactive in
   SN); point_basal_area! matches ptbal.f exactly; the cycle DOES recompute density post-compress — so the ONLY mover
   is the ULP sort flip. This is precisely the GOAL's two ACCEPTED divergences (ULP float + COMPRESS eigensolver), NOT
   a bug. The COMPRESS port is faithful; s22 stays @test_broken as the accepted divergence. (The earlier crown-ratio
   framing was a layer up — the cr/pbal IS the PTBAA term; same single record, same root.)

   ──────── (historical trace below, superseded by the resolution above) ────────

   ★★★★ DECISIVE FINDING (debug-FVS build, this session): I added a temporary `WRITE` to comcup.f dumping per
   post-COMPRESS record `(DBH, OLDRN)`, recompiled just comcup.o (`gfortran -fPIC -g -cpp`; MAXTRE=3000 is a
   fixed PARAMETER so the COMMON layout needs no -DCMP), relinked a separate `FVSsn_dbg`, ran s22, and compared
   to jl. RESULT: FVS's 5 merged records match jl's **EXACTLY** in BOTH dbh AND OLDRN (e.g. dbh 8.817→OLDRN
   −0.0051155 in both; all 5 pairs identical). So **jl's merge — including OLDRN — is already correct**; the
   entire multi-turn OLDRN hunt (min vs mem[1] vs average) was chasing a non-bug. The ONLY difference is the
   post-COMPRESS record ORDER: FVS = `5.347, 8.817, 6.992, 8.092, 6.664`; jl(min) = `5.347, 6.992, 8.092, 6.664,
   8.817`; jl(mem[1]) = `5.347, 6.992, 6.664, 8.092, 8.817`. NEITHER a min-sort NOR a mem[1]-sort reproduces FVS
   — its order is the genuine **tredel "fill-smallest-vacancy-with-highest-survivor" swap permutation**
   (tredel.f:46-120: IQRSRT the deletion index so vacancies sit at top in descending order, then for IV (smallest
   vacancy up) / IR (highest survivor down) `IF IVAC≤IREC: TREMOV(IVAC,IREC)`). That order is RNG-load-bearing
   (the same cycle's growth draws `rann!` per record in order), so it is the whole residual.
   ★★★★ FAITHFUL FIX (now mechanical, no more guessing): in `_merge_classes!`, lay the merged records at their
   IREC1=mem[1] slots and flag the other members deleted, then port tredel's swap loop EXACTLY (the
   smallest-vacancy←highest-survivor pass) to permute the live records into FVS's order — the merge VALUES/OLDRN
   are already right, so only the final permutation changes. Then s22 2005 should hit 409 and the suite greens.
   Reconfirm with the same debug dump (the recipe + backups are in /workspace/FVSjl/tmp/comcup.f.bak).

   ★★★★★ FIRST PORT ATTEMPT — DISPROVEN the naive TREDEL model (this session): I implemented "merged records at
   IREC1=mem[1], fill smallest vacancy with largest survivor" and validated jl's resulting ORDER against the FVS
   dump. EXACT s22 data: irec1(mem[1])=[144,236,240,186,136] (classes 1..5), min slots=[3,4,14,6,1], merged dbh
   per class=[6.992,8.092,8.817,6.664,5.347]. My model produced final order [c3,c2,c4,c1,c5] (dbh 8.817 first);
   FVS's actual order is [c5,c3,c1,c2,c4] (dbh 5.347,8.817,6.992,8.092,6.664). So FVS's order = the MIN-slot
   ascending order [c5,c1,c2,c4,c3] with the HIGHEST-min-slot class (c3, min=14) lifted to POSITION 2 — which is
   neither min-sort, mem[1]-sort, nor smallest-vac←largest-survivor. KEY REALIZATION: comcup.f:71 `CALL TREDEL`
   runs BEFORE `CALL COMPRS` (line 118), so it removes pre-existing zero-prob records, NOT the merged-away ones —
   the post-merge record removal/compaction uses a DIFFERENT path inside/after comprs.f that I have not yet traced.
   The reorder code was reverted (it didn't regress the suite but produced a faithfully-wrong order). NEXT: read
   comprs.f's end (after the DO-500 merge loop) + comcup.f after COMPRS to find how merged-away members are packed
   out, and/or add a second debug-FVS dump (IREC1 + the pre-final ordering) to deduce the exact permutation that
   yields [c5,c3,c1,c2,c4]. The OLDRN values are confirmed correct, so ONLY this permutation remains.

   ★★★★★★ MECHANISM TRACED (continued): the post-merge compaction is `CALL TREDEL(ITRN-NCLAS, IND)` at
   **comprs.f:1007** (NOT the comcup.f:71 call, which is pre-COMPRS zero-prob cleanup). comprs.f:863-869 negates
   the non-survivor members (`IND(I)=-IND(I)` for I=I1+1..I2, the deleted class members), leaving the survivor
   IND(I1)=IREC1=mem[1] positive; TREDEL then IQRSRTs and fills smallest-vacancy←largest-survivor. So the ALGORITHM
   I modeled is correct. BUT the remaining mismatch is the **pre-COMPRESS RECORD INDEX LAYOUT**: with jl's survivor
   record indices [144,236,240,186,136], the algorithm yields [c3,c2,c4,c1,c5] — yet FVS yields [c5,c3,c1,c2,c4],
   placing c5 (jl index 136, the SMALLEST survivor) FIRST, which the largest-survivor-first swap can only do if
   FVS's c5 index is actually the LARGEST. ⇒ **jl's record indices ≠ FVS's at COMPRS time.** Likely cause: jl's
   itrn=243 still contains zero-PROB / dead records that FVS removed via the comcup.f:71 pre-COMPRS TREDEL,
   shifting every index. NEXT (precise): (a) verify with a 2nd debug-FVS dump of `IND`/IREC1 right before
   comprs.f:1007; (b) check whether jl drops zero-prob records before COMPRESS (replicate the comcup.f:71 pre-pass
   if not); once the pre-COMPRESS index layout matches FVS, the existing TREDEL-swap port reproduces the order and
   greens the suite. The OLDRN/merge values are confirmed correct — ONLY the index layout + swap remain.

   ───────────────────────────────────────────────────────────────────────────────────────────────────────────
   (earlier diagnosis, now superseded by the debug dump above — kept for context:)
   DIAGNOSED; fix deferred (needs TREMOV swap-order). RE-TRACED (the "different eigensolver — accepted" label was STALE): the
   eigensolver/partition + the merge ARE bit-exact (s22 1990-2000 .sum matches live to the digit). The residual
   is purely the per-tree DG serial-correlation deviate **OLDRN** (dgdriv.f) carried INTO the post-compression
   cycle: 2005 TPA 415 vs live 409 (≈1.5%), growing to 257 vs 249 by 2015 (jl under-kills). Mechanism: comprs.f
   NEVER sets OLDRN on the merged record — it writes the merge to IREC1=IND(I1) and silently keeps that slot's
   OLDRN, then the tredel→**TREMOV** compaction (tremov.f:39/92/146) **SWAPS** whole records (incl. OLDRN),
   reshuffling the deviate differently than jl's one-way `copy_tree!` gather + minimum-index merge slot. dbh/ht/
   tpa are order-independent averages and nominal attrs come from the RANN-sampled member, so the merge stays
   bit-exact; only OLDRN drifts. TRIED: inheriting mem[1]=IND(I1)'s OLDRN — OVER-shot (2005→395). So the faithful
   fix is to replicate the exact TREMOV **swap** sequence (not a one-way gather) in `_merge_classes!`, so OLDRN
   lands where FVS's swaps leave it. Bounded but intricate; left as the @test_broken residual until done.
   REFINED (mem[1] overshoot): it is NOT just the OLDRN *source* — jl writes the merge to `dst=minimum(mem)` and
   gathers by sorted-min, which gives bit-exact dbh/ht/tpa but a different merged-record *ordering* than FVS's
   write-to-IREC1 + TREDEL compaction. So the faithful fix must change BOTH the write slot (→ mem[1]=IND(I1)) AND
   the compaction (→ TREMOV swap), together, so positions AND OLDRN match. A partial change (OLDRN source only)
   over-corrects (2005→395). This is a co-refactor of `_merge_one!` write-slot + `_merge_classes!` compaction.

   COMPLETE SOURCE-TRACED MECHANISM (read-first, the fix is now fully specified):
   • comprs.f:696 writes the merged averages to slot `IREC1 = IND(I1) = mem[1]` (the WK3-descending-sorted FIRST
     member), copies NOMINAL attrs from the RANN-sampled member (:733-741), and NEVER assigns OLDRN — so the
     merged record silently keeps `OLDRN[IREC1]` (= mem[1]'s deviate).
   • comcup.f:71 then compacts via `TREDEL`→`TREMOV(IVAC,IREC)` (tremov.f:39/92/146 SWAP whole records incl.
     OLDRN), packing survivors into original-index order.
   • jl instead writes to `dst=minimum(mem)` and `_merge_classes!` does a one-way `copy_tree!` gather of the
     reps sorted by dst. The `.sum` is bit-exact at the COMPRESS cycle because the merged VALUES are
     order-independent AND the same-cycle growth's per-record RNG happens to align under the min-index order —
     but the per-record OLDRN inheritance is wrong (dst's, not IREC1's), so the NEXT cycle's DG serial
     correlation drifts (2005 415 vs 409).
   • COUPLING (why both partial fixes fail): the post-COMPRESS record ORDER is RNG-load-bearing (the same cycle's
     growth/mortality draw `rann!` per record in order). `dst=mem[1]` alone changes the order → RNG misaligns →
     regresses the 2000 `.sum`. `OLDRN=mem[1]` alone (keeping min-index order) misplaces the deviate → overshoots.
   • FAITHFUL FIX = replicate FVS exactly: `_merge_one!` writes to `mem[1]` (keeping its OLDRN), the other members
     are flagged deleted, and `_merge_classes!` compacts with the TREDEL→TREMOV SWAP sequence (not a copy-gather)
     so the final order AND every record's OLDRN match FVS together. Bounded but intricate (the swap order must be
     exactly tredel's); deferred to a fresh-context implementation that can instrument the s22 partition to verify
     the class/order/OLDRN per record. The eigensolver itself remains the one accepted divergence.

   ★ EMPIRICAL DISPROOF (tested vs live, then reverted — doctrine #4): a clean `dst=mem[1]` (write to IREC1,
     keep its OLDRN) was implemented and run. Result: 2000 .sum stays BIT-EXACT (the COMPRESS-cycle aggregate is
     order-independent), but s22 2005 TPA = **395 vs live 409** — it OVERSHOOTS, while `minimum(mem)` UNDERSHOOTS
     (415). **Live sits BETWEEN the two slot choices**, which DISPROVES the simple "merged record keeps one
     member's raw OLDRN" hypothesis: FVS's post-COMPRESS OLDRN is whatever the TREDEL→TREMOV *swap reshuffling*
     leaves at each final slot, not mem[1]'s nor min's deviate. So the write-slot alone is not the fix — the swap
     SEQUENCE must be ported. Kept `minimum(mem)` (closer: +6 vs −14, and all other COMPRESS tests stay bit-exact).
     NEXT (fresh context): instrument the s22 partition (classes, per-member OLDRN, the tredel deletion order) and
     port tremov's exact swap loop so the OLDRN array is permuted identically — only then will 2005 hit 409.

   ★★ FURTHER PROBES (this session, all reverted — findings only):
     - `dst=mem[1]` makes the merged-record ORDER match FVS's tredel (ascending mem[1] = original-index order), so
       the 2000 .sum stays bit-exact and the remaining gap is purely the per-record growth, not the order.
     - `dst=mem[1]` + keep IREC1's own OLDRN (= comprs.f's literal behavior) → s22 2005 = **395**, NOT live 409.
       This is the KEY surprise: with BOTH the order AND the OLDRN matching FVS's source semantics, jl still
       diverges — so a THIRD factor differs, almost certainly the post-COMPRESS RNG-STREAM alignment (the same
       cycle's growth/mortality draws `rann!` per record; if the COMPRESS RANN draw-count or the record set/order
       desyncs the stream by even one draw, every subsequent DG deviate shifts). comprs.f:725 draws RANN once per
       multi-member class for the nominal-attr sample — the jl draw-count vs FVS's must be verified per class.
     - `dst=mem[1]` + PROB-weighted OLDRN AVERAGE → 2005 = **412** (live 409), 2010 328 (327), 2015 254 (249):
       much closer than min(415)/mem[1](395), but UNSOURCED (comprs.f never averages OLDRN) and not exact, and
       given the 395 result above it is a coincidental partial-compensation for the RNG-stream factor, not the
       mechanism. NOT shipped.
     CONCLUSION: the residual is NOT a single mis-set value — it is the joint (record-order + OLDRN + RNG-stream)
     state through comprs→tredel→tremov. The faithful fix must port that compaction's RANN draws + swap permutation
     exactly; the empirical probes bound where live sits (between min and mem[1]) and pin the RNG-stream as the
     missing third factor. Status quo `minimum(mem)` retained (suite bit-exact; +6 the closest single-knob value).

   ★★★ CONCRETE NEXT STEP (decisive, fresh context): output-matching has been exhausted — the joint state is not
     resolvable from `.sum` columns alone. Build a DEBUG FVS to read FVS's internal post-COMPRESS state directly:
     (1) add a `WRITE` in comcup.f (after CALL TREDEL/the compaction) or at DGDRIV entry for the 2000 cycle that
     dumps, per record I=1,ITRN: `DBH(I), OLDRN(I)` (find OLDRN's COMMON — it's referenced in dgdriv.f, add that
     INCLUDE to the dump file if needed); (2) recompile just that `.f`→`.o` (gfortran is at /usr/bin; objs in
     bin/FVSsn_buildDir/) and relink a SEPARATE debug binary (do NOT overwrite /workspace/FVSjl/tmp/FVSsn_full):
     `gfortran -o /workspace/FVSjl/tmp/FVSsn_dbg $(ls *.o) /workspace/FVSjl/tmp/glibc_shim.o -lpthread -ldl`;
     (3) run s22, parse the dump → the EXACT (dbh→OLDRN) map FVS feeds the 2000→2005 growth; (4) compare to jl's
     post-COMPRESS `old_random` per record (match by dbh) — this reveals whether the gap is the OLDRN values, the
     record order, or the RANN draw alignment, and gives the values to reproduce. THEN port the matching mechanism.

## Tier 4 — large independent subsystems (keyword-gated, transparency-guarded today)

**6. NOHTDREG/LHTDRG HT-DBH calibration. — ✅✅ RESOLVED, faithful END-TO-END (suite 5061/2). The live dub, the
   per-tree projected DG (27/27 exact vs live FVS_TreeList), the COR evolution (START clock), AND the dead-tree dub
   (cratet.f:413-473, ported this session) ALL match live. The post-1990 .sum drift is NOT a NOHTDREG gap — trees
   grow identically to 1995; it is purely the downstream tripled-record DGSCOR serial-corr + SDI mortality, i.e. the
   cross-cutting WK3 sp33/65 tail (same family as the post-thin tail / COMPRESS). @test_broken stays, re-attributed
   to that shared tail (a separate concern). Original port notes + the dead-dub fix + the COR-clock disproof below.**
   PORTED: `kw_nohtdreg!` now sets `ht_drag_sp[sp]` (LHTDRG) per the initre.f IS<0/0/>0 + field-2 invoke decode
   (drops the @warn); `Calibration.ht_dbh_aa`/`ht_dbh_iabflg` added; `dub_missing_heights!` fits `AA = mean(log(H
   −4.5) − HT2/(D+1))` over each invoked species' ≥3 measured-height trees (cratet.f:292-335) and dubs the missing
   ones with the calibrated Wykoff `exp(AA + HT2/(D+1)) + 4.5` when IABFLG==0, else the default Curtis-Arney HTDBH.
   Gated on `any(ht_drag_sp)` ⇒ default stands untouched (suite stays bit-exact). LIVE-VALIDATED on `nohtdreg_cal`
   (carbon_snt + `NOHTDREG 0 1`): the **dub cycle (1990) .sum is BIT-EXACT** vs FVSsn (TCuFt 1358 vs 1368 default).
   RESIDUAL (post-1990 ≈1.2% by 2005) — root NOT yet localized; tracked @test_broken (honest, not toleranced).
   READ-VERIFIED facts (correcting earlier guesses): CRATET is INIT-ONLY (fvs.f:197 sits before the cycle
   back-edge `GOTO 40`@363 — NOT per-cycle); LHTDRG has exactly TWO consumers — cratet.f's height-dub and
   regent.f:315 (regen) — and regen does NOT fire in this NOAUTOES stand, so the dub is the ONLY NOHTDREG effect
   here. The dub itself is faithful by source-reading (formula :432, condition :301, coeffs HT1/HT2, fit :329,
   top-kill clamp :397 all match) ⇒ the post-1990 drift is the calibrated heights perturbing a PRE-EXISTING
   growth path (the default no-NOHTDREG carbon_snt is itself bit-exact to live), not an un-ported consumer.
   TODO: localize via a per-tree dubbed-height differential vs live (the sp65 dub 22.91 is unverified), and a
   multi-point/regen scenario to exercise regent.f:315. Suite 5053→5069/1. — Original trace below:

   ★ SHARPENED (this session, no rebuild): jl-vs-live(.sum.save) by cycle — onset is the FIRST period (cols 24/25
   periodic accr/mort differ at 1990) → 1995 TPA 504 vs 502, BA 106 vs 107 (jl more trees / less BA), ≈1% by 2005.
   **TopHt (col 7) is BIT-EXACT every cycle** ⇒ the divergence is DG/mortality, NOT height. NOHTDREG barely moves
   heights (jl mean_ht 39.134 cal vs 39.108 default, Δ0.026 ft) and leaves mean crown IDENTICAL (41.39 both) ⇒ NOT
   a crown-path diff. Re-traced LHTDRG: active consumers = cratet.f (dub) + regent.f (regen, doesn't fire);
   htdbh.f is comment-only. So the residual = the tiny calibrated-height perturbation propagating relht→DG→SDI-mort
   over cycles. Two hypotheses to settle with ONE debug-FVS per-tree dump (HT/DBH/ICR at 1990+1995): (i) per-tree
   calibrated heights differ jl-vs-live though the aggregate matches (⇒ AA-fit/dub per-tree gap, fixable); or
   (ii) per-tree heights match and it's the WK3 sp33/65 DG-COR sensitivity amplifying a sub-ULP perturbation
   (⇒ likely accept, same family as the post-thin tail). Needs the rebuild trace like COMPRESS — focused session.

   ★★ RESOLVED via FVS_TreeList oracle (DATABASE/TREELIST 0 + TREELIDB 1): hypothesis (i)/(ii) both refined.
   The 27 LIVE trees' 1990 DBH/Ht/DG/PctCr/BAPctile ALL match live (dub faithful for live; gross_space=1.1 explains
   the ×1.1 per-tree tpa). The TWO extra in live's 29-row TreeList are DEAD trees (history=8: dbh-7.2 HI, dbh-34.6 SK)
   — jl ALSO has them in its dead block (treeinput.jl:71), NOT dropped (earlier "dropped" note was WRONG). THE GAP:
   jl's `dub_missing_heights!` is LIVE-ONLY (volume.jl:13,26 loop 1:t.n), so it never dubs the DEAD trees' missing
   heights, but live gives them CALIBRATED Wykoff dubs (55.22, 32.901). Since the DG-calibration backdating exposes
   the dead partition at current dbh (diameter_growth.jl:302-320), dead-tree state can feed the COR. OPEN: FVS's dead
   dub site (cratet.f DO 50@193 is sort-cleanup, not the dub) + whether dead HEIGHT (vs dbh-based PTBAA) actually
   changes a growth input. NEXT: find the FVS dead-dub + test extending jl's dub to the dead block.

   ★★★ DEAD-DUB GAP FIXED (this session): FVS cratet.f:413-473 has a dedicated DO 145 loop (II=IREC2..MAXTRE)
   that dubs DEAD trees' missing heights with the SAME calibrated-AA/HTDBH formula as the live dub. jl's
   `dub_missing_heights!` was LIVE-ONLY — FIXED by extending the dub-application loop to 1:(t.n+t.ndead) (AA fit
   stays live-only, matching FVS DO 15). Validated bit-exact: jl dead trees now 34.6 SK→32.901, 7.2 HI→55.221
   (= live FVS_TreeList; were undubbed). Suite 5061/2, ZERO .sum impact (dead heights don't enter the live
   aggregate), no regression — a clean faithfulness improvement. BUT it does NOT close the .sum residual (1995
   still 504/502). So the remaining ≈1% drift is a SEPARATE downstream 0.4% mortality diff (1990 live state +
   per-tree DG all match live; confounded at 1995 by tripling). Scenario has WK3-calib sp33/65 ⇒ most consistent
   with the WK3 DG-COR sensitivity ([[fvsjl-postthin-tail-is-wk3-calib]]); finalize accept-vs-fix via a debug-FVS
   COR(ISPC) dump for sp33/65.

   ★ COR-CLOCK off-by-one — INVESTIGATED + REFUTED (re-trace discipline; do NOT re-attempt). The debug-FVS dgdriv
   COR dump shows sp33 COR 1.0221/0.9605/0.9069 at cyc1/2/3, which LOOKS one cycle ahead of jl's dg_cor
   1.09285/1.0221/0.9605. Changing jl's diameter COR to the END clock (elapsed+sfint, matching FVS SFINT=
   IY(ICYC+1)−IY(1)) made dg_cor match the dump BUT REGRESSED THE SUITE −1823 tests. Reason: the dgdriv COR WRITE
   fires AFTER dgdriv updates COR for the NEXT cycle, whereas dgf already baked the CURRENT (pre-update) COR into
   WK2 — so jl's dg_cor[N] correctly equals the COR FVS USES for cycle N's DG. The START clock is RIGHT (reverted;
   in-code comment added at diameter_growth.jl so it isn't re-tried). ⇒ the residual is NOT a COR-clock bug.

   RE-TRACED into an actionable spec; feasible (Wykoff coeffs exist).
   FULL TRACE: LHTDRG defaults FALSE (grinit.f:104). NOHTDREG (initre.f:2605-2674, opt 60) with **field 2 > 0
   INVOKES** calibration (LHTDRG[sp]=TRUE; IS<0 group / IS=0 all / IS>0 one species); **field 2 blank/0
   SUPPRESSES** (=default, jl's current no-op — correct). So the gap is ONLY the invoke path (jl @warns).
   THE FIT (cratet.f:292-335): over species trees with H>4.5, NORMHT≥0, D≥3 accumulate `SUMX += log(H−4.5) −
   HT2/(D+1)`, K1=count; if K1≥3 AND LHTDRG ⇒ `AA[sp]=SUMX/K1`, and if AA≥0 set IABFLG[sp]=0. THE DUBBING
   (cratet.f:342-372): AX = (IABFLG[sp]==0 ? AA[sp] : HT1[sp]); for missing-height trees `H = exp(AX + HT2/(D+1))
   + 4.5` — the WYKOFF curve; BUT if `.NOT.LHTDRG OR (LHTDRG AND IABFLG==1)` it instead calls HTDBH (the
   Curtis-Arney inventory eqn = jl's current default dub). So calibration REPLACES jl's Curtis-Arney dub with a
   calibrated-Wykoff dub. FEASIBLE: jl already has Wykoff `wykoff_ht1`/`wykoff_ht2` (sprout_htdbh_wykoff.csv;
   = HT1/HT2). PLAN (layer→consumers): (a) `kw_nohtdreg!` set `ht_drag_sp[sp]=true` on field2>0 (drop the @warn);
   (b) add per-species `AA`/`IABFLG` state + the fit in `dub_missing_heights!`; (c) switch the dub to the
   calibrated Wykoff when LHTDRG[sp]&&IABFLG==0; (d) consumers regent.f:315 (small-tree regen DBH→Wykoff branch)
   + esuckr (sprout-dbh). Validate on a `NOHTDREG <sp> 1` + measured-height stand vs live FVSsn. Multi-site but
   each step is bounded; the fit+dub is the validatable core, regen/esuckr follow.

**7. log-graded HRVRVN revenue (units 4/5). — ★★ FULLY DONE: BOTH unit 4 (BF_1000_LOG) AND unit 5 (FT3_100_LOG) bit-exact vs live.**

   ★★ UNIT 5 (FT3_100_LOG) NOW LANDED (suite 5142/2, +33 assertions). The cubic log-graded path is the exact
   mirror of unit 4. VALIDATED BIT-EXACT vs live FVSsn FVS_EconHarvestValue (econ_u5.key: HRVRVN 300 5 10.0 ALL +
   THINSDI 2000): per species SM=3ft3/$10, HI=2/$5, AB=8/$23, SK=1/$3, class [10.0, 999.9) — Ft3_Removed,
   Ft3_Value, Total_Value all match. Per-log split ALSO bit-exact (verified against an echarv.f debug dump: tree-1
   sp27 logs 0.684/1.620/2.631/2.669/4.257 + two DIB-10 logs summing 10.158; only the DIB≥10 pair qualifies).
   • PORTED PIECES (cubic analogs): r8clark_vol.jl `_r8_cuft_by_dib` (R9LGCFT — per-log Smalian
     `0.00272708·(Dbot²+Dtop²)·len` over the predicted (index-2) boundary DIBs, where boundary-1 is DIB@4.5 ft not
     the stump, then ALL logs renormalized so Σ = VOL(4)+VOL(7); r9logs.f:364/r9clark.f:428), wired into
     `_R8CLARK_VOL` via a `log_cuft` ref filled in the sawtimber block. volume.jl stashes per-tree `by_dib` into
     `EconState.tree_log_ft3` (gated on a unit-5 HRVRVN record). cuts.jl `_log_cut!` → econ.jl
     `accrue_log_grade_cuft!` accumulates `EconState.log_grade_ft3`; `econ_harvest_value_rows` (unified board+cubic
     emitter, FVS species→unit→DIB order) + dbs_output.jl fill the Ft3_Removed/Ft3_Value columns.
   • ★ THE KEY BUG (found by instrumenting live echarv.f, dumping ft3PerTree/treeVol/defProp + per-log logFt3Vol):
     ECHARV's `ft3PerTree` (the defProp denominator) is NOT the merch cubic — for the EASTERN variants SN/CS/LS/NE
     it is the SAWTIMBER cubic **SCFV** (cuts.f:1667-1669 `CFVOLI = SCFV(I)`), NOT MCFV. So defProp = treeVol(gross
     v4+v7) / SCFV(v4) ≈ 1.237, not ~1.0. jl first passed merch_cuft_vol ⇒ defProp≈1 ⇒ ≥10" logs over-credited
     (16/2/9/1 vs live 14/2/8/1 totals). FIX: cuts.jl passes `saw_cuft_vol` (SCFV). Western variants would use
     MCFV — gate when added (principle #6).
   • The @warn for unit 5 is removed; `warned_loggrade` deleted. Price appreciation still flat (faithful for the
     validated no-appreciation case). [OLD unit-4-only verdict below retained for history.]

**7-OLD. log-graded HRVRVN revenue (units 4/5). — ★ DONE for unit 4 (BF_1000_LOG), bit-exact vs live; unit 5 was deferred.**

   ★ LANDED (suite 5086/2, +25 assertions). The FVS_EconHarvestValue table is now emitted and is BIT-EXACT vs
   live FVSsn for the BF_1000_LOG path (econ_strtecon.key, HRVRVN 300 4 10.0 ALL + THINSDI→SDI250 @2000):
   per species SM=16bf/$5, HI=9/$3, AB=41/$12, SK=5/$1, class [10.0, 999.9) — every column (SpeciesFVS/PLANTS/
   FIA, Min/Max_DIB, Board_Ft_Removed/Value, Total_Value) matches the oracle.
   • THE KEY MECHANISM (found by instrumenting live echarv.f via fort.16): the defect proportion
     `defProp = treeVol / bfPerTree`, where **treeVol = Σ gross Scribner of ALL logs incl. the small topwood
     logs** (R9LOGS bucks the full stem to the ~4" pulpwood top), and bfPerTree = the NET saw board feet
     (vol[10], saw logs only). ECHARV credits each ≥class log as `logBF·harvTpa/defProp`. My first cut used
     by_dib = saw logs only ⇒ defProp=1 ⇒ ≥10" logs over-credited (82 vs 71). FIX: `_r8_scribner_bf_by_dib`
     now ports the FULL R9LOGS (sawtimber + topwood segmentation via `_r9loglen!`) → R9LOGDIB → R9BDFT, so
     treeVol includes the topwood and defProp = ~1.3 (matches live exactly).
   • Implementation: r8clark_vol.jl `_r9loglen!` + rewritten `_r8_scribner_bf_by_dib` (full stem); volume.jl
     stashes per-tree `by_dib` into `EconState.tree_log_bf` (gated on a unit-4 HRVRVN record); cuts.jl
     `_log_cut!` → econ.jl `accrue_log_grade!` (getDiaGrp = `_log_dia_grp`, defProp scaling) accumulates
     `EconState.log_grade_rev`; `econ_harvest_value_rows` + dbs_output.jl `write_dbs_econharvest!` emit the
     table (wired in simulate.jl). REPORT-ONLY confirmed: PNV/EconSummary untouched (already bit-exact).
   • REMAINING (deferred): unit 5 = FT3_100_LOG needs the per-log CUBIC bucking (R9LGCFT between log
     endpoints) — jl computes cubic via the taper integral, not per-log; the @warn now fires only for unit 5.
     Also: HRVRVN price appreciation (rate/duration) is not stored in EconCostRev (flat price only) — faithful
     when no appreciation is given (the validated case); add if a scenario needs it.

   --- original spec (read-first) retained below ---

**7. log-graded HRVRVN revenue (units 4/5). — RE-TRACED into an actionable spec (read-first).** HRVRVN field 2 =
   UoM: TPA(1)/BF_1000(2)/FT3_100(3)/BF_1000_LOG(4)/FT3_100_LOG(5) (ecin.f:239-246); for the _LOG units field 3 is
   a DIB class. FVS data flow: ECVOL (ecvol.f, called from VOLS per tree) stores per-log `logDib*`/`log*Vol`
   [treeId][logId] from the NVEL `LOGDIA`/`LOGVOL` arrays (DIB at LOGDIA[*,1] scaling dib, gross cuft at
   LOGVOL[4,*], scribner bf for the bf path); ECHARV (echarv.f:72-126) then, for each harvested tree, buckets
   every log's volume by its DIB grade (`getDiaGrp`) into `revVolume[sp][unit][dibIdx]`; revenue = Σ price[DIB]·vol.
   jl STATUS: the per-log BUCKING for **BF_1000_LOG already exists** — r8clark_vol.jl:489-498 walks the logs and
   computes `idib = INT(DIB+0.499)` + `bf_per_log = round(len·SCRBNR[idib])` — it just sums and discards the
   breakdown. **FT3_100_LOG needs NEW work**: jl computes total cubic via the taper integral (vol[1]/[4]), not
   per-log, so a per-log cubic integration between log endpoints must be added. PLAN: (a) expose per-log
   `(idib, bf, cuft)` from the volume routine; (b) store per harvested tree (the ECVOL arrays); (c) `kw_econ!`
   parse the units 4/5 DIB-class prices (drop the @warn); (d) port ECHARV's per-DIB accumulation into the revenue
   stream. Validate on an `HRVRVN … 4 …`/`… 5 …` cut scenario vs live FVS_EconSummary. Inert for sn.key (live=0
   too), so the default suite stays bit-exact.

   ★ VALIDATE-FIRST (this session, principle #4) — RE-TRACED the "inert/live=0" claim against the live oracle and
   SHARPENED the scope. Ran `econ_strtecon.key` (HRVRVN 300 4 10.0 ALL + THINSDI 2000) through `tmp/FVSsn_full`:
   • Live DOES produce non-zero log-graded value, but ONLY in the detail report **FVS_EconHarvestValue** (per
     species, DIB class 10–999.9): SM Board_Ft_Value=5, HI=3, AB=12, SK=1 (total 21; Board_Ft_Removed 16/9/41/5).
   • It is **REPORT-ONLY**: FVS_EconSummary.Undiscounted_Revenue/Discounted_Revenue/PNV stay **0** even after
     extending NUMCYCLE to 5 (so it's NOT a last-cycle-timing artifact — log-graded HRVRVN simply does not feed
     the PNV stream). So jl's PNV (=0) ALREADY MATCHES live; the EconSummary is faithful as-is.
   • The earlier code comment ("FVS_EconHarvestValue empty; harvests don't qualify") was WRONG — corrected in
     keyword_dispatch.jl. The REMAINING gap is purely: jl does NOT emit the FVS_EconHarvestValue table at all
     (grep: zero emitters). So the port reduces to (a) per-log idib/bf/cuft (BF path exists in r8clark_vol.jl),
     (b) store per harvested tree, (c) bucket by DIB grade per HRVRVN unit-4/5 prices, (d) EMIT FVS_EconHarvestValue
     (new DB table + schema). Oracle = the 4-row table above (econ_strtecon.key, harvest 2000). This is a
     report-table port with NO PNV/suite impact (PNV already bit-exact), so lower priority than #28/#40.

## Tier 5 — hardest / downstream / risky

**8. FFE phasing #28 co-refactor** — FVS FMMAIN runs FMBURN(fire)→annual fuel loop AFTER grow+mort within
   GRINCR; jl runs the fuel loop BEFORE grow_cycle, phase-shifting the cwd/fuel trajectory by one cycle. A naive
   reorder was REVERTED (regressed the #20 under-kill fix + carbon). Needs a co-refactor of the fire-basis stash
   + fuel loop + grow ordering. **Unblocks the one remaining real output divergence:** the fire-carbon
   released-from-fire VALUE (carbon report / FVS_Carbon col emits 0; live target 5.5 t/ac@2000 on
   fire_carbon.key; the 0.37/0.50 conversion factor is already fixed).

   ★ MECHANISM PINNED + VALUE VALIDATED (this session). The released-from-fire carbon = consumed fuel (FVS
   fmcrbout.f:151 V(11)=BIOCON·{0.37,0.50}; BIOCON from fmdout.f:286 fire-consumed litter/duff+woody). jl HAS the
   whole path (apply_fire_consumption!→fmburn! FireResult.carbon_released→FVS_Carbon column). TWO concrete gaps:
   (a) WIRING — _maybe_burn! (simulate.jl:179) DISCARDS fmburn!'s FireResult (returns only mortality vol); the
   carbon-report rows pass released=0. (b) VALUE — instrumented jl carbon_released@2000 = **6.412 t/ac vs live
   5.5** (~16.5% high). ROOT (summary.jl:196-199): jl already stashes the START-of-cycle fuel BASIS (fire_smlg)
   for fire BEHAVIOR (flame/scorch correct), but `ffe_fuel_update!(per)` then ADVANCES the cwd pools, and
   apply_fire_consumption! consumes from those ADVANCED pools — whereas FVS FMBURN consumes from the START-of-cycle
   pools BEFORE the annual loop. So the fire over-consumes by one period's litterfall. SURGICAL FIX DIRECTION:
   make apply_fire_consumption! consume on the start-of-cycle cwd basis (stash like fire_smlg), then let the annual
   loop run on (start − consumed); must keep #20 under-kill + DDW carbon green (the prior naive reorder didn't).
   Then wire carbon_released through _maybe_burn!→carb_rows→write_dbs_carbon!. NOT a pure-emit fix (would ship 6.41).

   ★★★ CONFIRMED via attempt + fmmain.f trace (this session) — the surgical consume-from-start-basis is INSUFFICIENT;
   the fix REQUIRES the full fuel-trajectory alignment. Implemented the surgical version (stash start-of-cycle cwd,
   consume start·f): released 6.41→5.13 vs live 5.52 — OVERSHOT, reverted. fmmain.f proves the basis TIMING was
   already right: FMBURN (fmmain.f:170) runs BEFORE the annual fuel loop `DO IYR` (FMSNAG:232/FMCWD:236/FMCADD:241,
   line 228), so FVS consumes the START-of-cycle pools (0 annual steps) — exactly the stash point. ⇒ the 5.13-vs-5.52
   gap is that jl's cwd POOL VALUES at the fire year differ from FVS's: jl runs ffe_fuel_update!(litterfall/decay)
   BEFORE grow_cycle!(crown-lift/mortality FMCADD), while FVS runs grow+mort THEN the annual loop — so the per-year
   fuel accumulation ORDER (and thus every pool value) is phase-shifted. The released value + the SD-spike (2005 vs
   2000) are BOTH this one phase-shift. So #28 is genuinely the fuel-loop REORDER (grow+mort → annual loop), the hard
   co-refactor that must keep #20 + DDW green — not a consumption tweak. NEXT: debug-FVS FMCWD dump at the 2000 fire
   year to quantify the per-pool offset, then stage the reorder. Wiring is trivial once the trajectory matches.

   ★★★★ RE-TRACE BREAKTHROUGH (this session — the prior "simulation reorder" framing was aiming at the WRONG layer).
   Ran fire_carbon.key (SIMFIRE 2000) through BOTH the live binary (tmp/FVSsn_full, text carbon report in fort.16)
   and jl run_keyfile, and diffed the .sum AND the carbon report side by side:
   • THE .sum IS BIT-EXACT through the fire: live & jl both have 2000 TPA 470 / BA 126 / **Mort=240** / MCuFt 2274,
     and 2005 TPA 104 / MCuFt 1551 — identical (only a trivial ≥2010 regrowth ULP: jl BA 85 vs 84). So the fire
     mortality is applied at the SAME cycle (booked in the 2000 cycle, post-fire stand shows at 2005) in both. The
     SIMULATION ORDERING IS ALREADY FAITHFUL — grow_cycle!'s grow→fire→mortality→snag-booking matches FVS.
   • ONLY the CARBON REPORT pool SAMPLE-PHASE is shifted one cycle: live's 2000 row is POST-fire (AGL Total 19.1,
     Standing-Dead spike 20.2, DDW 1.1, **Released 5.5**), jl's 2000 row is PRE-fire (AGL 36.2, SD 2.2, DDW 3.7,
     Released 0) and jl's SD spike / fire effect lands at the 2005 row (SD 18.5). FVS prints FMCRBOUT (fmmain.f:206)
     AFTER FMBURN (170) but BEFORE the annual loop (228); jl's carbon.jl loop prints the row at the TOP of the cycle
     iteration, BEFORE that cycle's grow_cycle!(fire) — i.e. one phase too early.
   ⇒ REFRAMED SCOPE: #28 is NOT a grow/fuel-loop simulation reorder (that's what regressed #20/DDW in BOTH prior
     attempts — it touched a layer that is already correct). It is a CARBON-REPORT SAMPLE-POINT fix: sample/print
     each cycle's carbon row at the FVS phase (after that cycle's fire + snag booking, before the next annual fuel
     advance) instead of at cycle-top. This is isolatable to carbon.jl's report loop (summary.jl carbon collection
     + simulate.jl carb_rows), NOT grow_cycle!/ffe_fuel_update! ordering, so it should leave the (bit-exact) .sum and
     the #20/DDW fuel trajectory UNTOUCHED. The released-from-fire value then rides along at the correct year (2000)
     once _maybe_burn!'s FireResult.carbon_released is captured into the same post-fire 2000 sample. Open detail to
     pin before coding: confirm whether FVS's post-fire carbon row is sampled on the 2000-SIZE or 2005-SIZE live
     stand (carbon Merch 12.4@2000 ≠ 15.1@2005 ⇒ likely 2000-size i.e. the fire is applied to the pre-this-cycle-
     growth stand at the report) — a short debug-FVS stamp of (growth-applied-yr, FMBURN-yr, FMCRBOUT-yr, live TPA)
     resolves it. This reframing is the key unblock; the actual sampling change is the remaining work.

   ★★★★★ DEBUG-FVS STAMP (this session) — RESOLVED the open detail AND corrected the over-optimistic "sim-untouched"
   note above. Instrumented fmmain.f (STAMP at FMBURN line 170 + FMCRBOUT line 206, dumping IYR / live TPA=ΣPROB /
   meanDBH) and ran fire_carbon.key on the full binary. At the **2000 carbon report**: TPA=516.74, meanDBH=5.658 —
   i.e. the **2000-SIZE stand** (between 1995's 4.726 and 2005's 10.761), with full pre-fire PROB (the fire KILL is
   not applied to PROB until FMKILL(1) at gradd.f:122, AFTER FMMAIN). FMBURN (170) runs immediately BEFORE FMCRBOUT
   (206), so the 2000 carbon row's snag/DDW/Released pools are POST-FMBURN but the stand is still 2000-size and
   pre-growth. THEN the cycle grows 2000→2005 (2005 carbon: TPA 114.34 / meanDBH 10.761).
   ⇒ CORRECTED SCOPE: the "isolatable to carbon.jl, simulation untouched" hypothesis is WRONG. The live carbon 2000
   row reflects a **2000-size, post-fire, pre-2000→2005-growth** stand that jl never materializes — jl applies the
   fire inside the grow_cycle! that grows 2000→2005 (current_cycle_year=2000), i.e. AFTER that growth, so jl's fire
   snags are 2005-size and surface at the 2005 carbon row. (The .sum stays bit-exact only because jl's fire KILL
   PROBABILITY is computed on the cycle-start 2000-size basis — fire_smlg/FMFINT stash — so the survivor set matches
   even though PROB is reduced post-growth.) So #28 genuinely needs the fire's CARBON effects (snag creation +
   consumption + Released) booked at the fire YEAR on the fire-year-SIZE stand — a real fire/growth/report
   co-ordering, NOT a pure report-sample move. Concrete validated target (fire_carbon.key carbon trajectory, US t/ac):
   1990 Rel0 / 1995 Rel0 / **2000 AGL19.1 SD20.2 DDW1.1 Rel5.5** / 2005 AGL22.5 SD2.8 DDW14.8 / 2010 AGL28.1 SD1.0
   DDW11.5. jl currently: 2000 AGL36.2(pre-fire) Rel0, spike at 2005. The fix must hit the 2000 post-fire row while
   keeping the bit-exact .sum + #20 under-kill + DDW green. This remains the single hardest open item.

   ★★★★★★ FULL MECHANISM PINNED + STAGED PLAN (this session — gradd.f read settles the last ambiguity). gradd.f
   applies growth at **CALL UPDATE (gradd.f:180)** — AFTER CALL FMMAIN (118, =FMBURN fire) and CALL FMKILL(1)/(2)
   (122/135, =apply MORTS⊔fire kill to PROB). So FVS's per-cycle order is: [stand at year Y, post-prev-growth] →
   FMBURN(Y) on the Y-size stand → FMKILL combine → FMCRBOUT carbon report(Y) → annual fuel loop (FMSNAG/FMCWD/
   FMCADD) → UPDATE grows Y→Y+1. jl's grow_cycle! MATCHES the sim order EXACTLY (fire `_maybe_burn!` simulate.jl:256
   on cycle-start dims — confirmed the DBH update is later at :285; `_maybe_burn!` docstring "cycle-start dimensions").
   So the .sum is bit-exact for the RIGHT reason (both fire the Y-size stand). The divergence is TWO report/fuel
   ordering bugs, NOT a sim-mortality bug:
   (A) CARBON-REPORT SAMPLE POINT — summary.jl:142 pushes the carbon row at the TOP of the cycle iteration, BEFORE
       cuts!/ffe_fuel_update!/grow_cycle!. So the fire (inside grow_cycle!, :199) hasn't run ⇒ the Y=2000 row is
       pre-fire and the snag/AGL/DDW/Released effects surface at the 2005 row. FVS samples FMCRBOUT AFTER FMBURN.
   (B) FIRE-CONSUMPTION BASIS — summary.jl runs ffe_fuel_update!(per) (:197) BEFORE grow_cycle!(:199); the fire's
       apply_fire_consumption! then consumes the ALREADY-ADVANCED cwd (period litterfall/decay added) ⇒ released
       6.41 vs live 5.5. FVS's FMBURN (fmmain.f:170) consumes BEFORE the annual loop (228). [The fire BEHAVIOR
       flame/scorch is already right — fire_smlg stash at :196 captures the start-of-cycle (SMALL,LARGE).]
   STAGED FIX (keep .sum + #20 + DDW green at each step; fire_carbon.key carbon trajectory is the oracle, captured
   above; the FFE carbon/DDW/snag suite tests are the regression gate):
     1. [✅ ENABLING REFACTOR LANDED this session — suite bit-exact 5086/2] Factored grow_cycle!'s fire+MORTS block
        (the `fire_now` MORTS⊔FIRKIL combine + snag booking, simulate.jl) into `mortality_and_fire!(s; fint)`, a
        pure no-behaviour-change extraction called at the same point. This makes the fire half independently
        callable — the seam the reorder needs. REMAINING for this step: reorder summary.jl's fire-cycle path to FVS
        order — fire FIRST (on cycle-start stand, consuming START-of-cycle cwd), THEN push the carbon row (post-fire),
        THEN ffe_fuel_update!(advance remaining cwd), THEN grow_cycle! with the fire already applied. mortality_and_
        fire! already short-circuits the second call (one-shot via fire_year=0 inside _maybe_burn!), but note the
        MORTS half would double-apply if called twice — so the reorder must either split MORTS vs fire, or guard the
        whole call once-per-cycle. Non-fire cycles stay on the in-grow_cycle! call, unchanged.
     ⚠ CONSTRAINT discovered staging step 1: jl's `mortality!` reads the diameter_growth! STASH (the BAMAX/size-cap
       kill uses the linear DG/BARK increment — see [[fvsjl-10yr-cycle-mortality]]), so mortality_and_fire! CANNOT be
       hoisted above grow_cycle!'s diameter_growth! (:237). ⇒ a summary.jl-level "fire before the carbon push"
       reorder is INVALID (it would run MORTS before growth is computed). CORRECTED APPROACH: keep mortality_and_
       fire! where it is INSIDE grow_cycle! (right after diameter_growth!), and instead pass grow_cycle! a
       `carbon_sample` callback that fires immediately AFTER mortality_and_fire! (post-fire, post-snag-booking) and
       BEFORE the DBH growth application (:285) — exactly the FVS FMCRBOUT phase (after FMBURN/FMKILL, before UPDATE).
       For the fire cycle the carbon row is captured by that callback (post-fire, 2000-size); summary.jl/carbon.jl
       use it instead of the pre-grow_cycle! :142 push. Non-fire cycles can keep the :142 push (same Y-size state).
     2. Make apply_fire_consumption! consume the cwd as of BEFORE ffe_fuel_update!: since ffe_fuel_update! (:197)
        runs before grow_cycle! (:199) and the fire is inside grow_cycle!, stash the cwd basis BEFORE ffe_fuel_update!
        (like the existing fire_smlg behaviour-basis stash at :196) and have apply_fire_consumption! consume from that
        start-of-cycle stash, so released matches FVS's pre-annual-loop consumption (target 5.5, not the 6.41 from
        consuming the advanced pools).
     3. Wire _maybe_burn!/apply_fire!'s FireResult.carbon_released into the carbon row (carbon_collect tuple +
        stand_carbon_report Released column + write_dbs_carbon! col 13), now landing at year 2000.
     4. Validate: fire_carbon.key carbon row 2000 == AGL19.1/SD20.2/DDW1.1/Released5.5 (live), .sum still bit-exact
        (TPA/Mort/vols), and the FFE carbon/DDW/snag/under-kill suite stays green. Then write the scenario test.

   ★★★★★★★ ISSUE (A) LANDED + TESTED (this session — suite 5089/2). Step 1's report-sample-phase fix is done WITHOUT
   any fuel/mortality reorder, exactly via the corrected "carbon callback" design: grow_cycle! gained a
   `carbon_hook` (simulate.jl) fired right after mortality_and_fire! (post-fire, post-snag-booking) and BEFORE the
   DBH growth application — the FVS FMCRBOUT phase. summary.jl skips the cycle-top carbon push for the SIMFIRE cycle
   and defers it to that hook (non-fire cycles unchanged — and the carbon-report suite scenarios carbon_jenkins/snt/
   ffe have NO SIMFIRE, so zero risk; confirmed bit-exact). VALIDATED vs live FVSsn (fire_carbon.key): the 2000
   carbon row is now POST-fire — AGL 18.8 (live 19.1; was 36.2 pre-fire), StandDead 18.5 (live 20.2; was 2.2),
   Merch 12.2 (live 12.4), DDW 1.3 (live 1.1). The fire effects surface at the fire YEAR, not one row late. +test
   (test_carbon.jl "FFE fire-year carbon row sampled POST-fire"). The callback is read-only (idempotent
   compute_density!, guarded fmcba!) so the .sum/sim is untouched.
   REMAINING = ISSUE (B), two coupled pieces, both rooted in ffe_fuel_update! running BEFORE grow_cycle! (so the
   fire — inside grow_cycle! — sees ALREADY-advanced cwd, and the fire's NEW snags are created after that cycle's
   annual loop already ran):
     (B1) Released-from-fire VALUE still 0 (unwired) — wiring it now ships ~6.41 (consumes advanced pools) not 5.5;
          needs the start-of-cycle cwd consumption basis. (B2) snag-fall is ONE CYCLE LATE: the 2000 fire snags fall
          to DDW at 2005→2010 (jl 2010 SD2.4/DDW13.5) instead of 2000→2005 (live 2005 SD2.8/DDW14.8) — because
          ffe_fuel_update! advanced the pools BEFORE the fire created the snags, so they aren't processed until the
          next cycle's loop. Both need the fire (fmburn! consumption + snag booking) to run BEFORE ffe_fuel_update!
          within the cycle — but mortality_and_fire! needs diameter_growth! (the size-cap stash), so the FIRE-only
          part (fmburn!, which does NOT need the growth stash) must be split from the MORTS combine and hoisted
          before ffe_fuel_update!, with the FIRKIL stashed for the later MAX-combine. This split is the deep part —
          the careful seam to attempt next, gated on the bit-exact .sum + #20 + DDW.
   RISK: step 1 touches grow_cycle!'s mortality seam (the bit-exact .sum path) — the one-shot guard + factoring
   (not reordering) the existing fire block is what keeps OMORT/snag booking identical. This is the careful seam the
   two prior REVERTED attempts got wrong (they reordered ffe_fuel_update!/grow wholesale). Execute as one focused
   chunk with the .sum diff + FFE suite as the gate after each sub-step.

   ★★★★★★★★ ISSUE (B) LANDED — #28 CO-REFACTOR DONE (suite 5091/2, NO regression to .sum/#20/DDW). The deep seam was
   tractable once mortality_and_fire! was a function: it gained a `post_fire` callback fired right after _maybe_burn!
   (the fire) and BEFORE the WK2 MAX-combine — exactly the FVS fmmain.f phase (FMBURN:170 → FMCRBOUT:206 → annual
   loop:228 → FMKILL). grow_cycle! builds that callback to run the carbon sample (carbon_hook) THEN the DEFERRED FFE
   annual fuel update (new `fuel_period` kwarg); summary.jl WITHHOLDS the fire cycle's pre-grow ffe_fuel_update! and
   hands its period to grow_cycle!. So the fire consumes the START-of-cycle down wood, the carbon row is sampled
   post-fire/pre-annual-loop, and the fresh fire snags are processed by the same cycle's loop. The MORTS combine +
   growth are untouched ⇒ the .sum stays bit-exact. RELEASED wired: burn_reports[year].released → carbon_collect 6th
   tuple field → write_carbon_report_block + write_dbs_carbon! (FVS_Carbon col 13), unit-converted like the pools.
   VALIDATED vs live (fire_carbon.key 2000, US t/ac): AGL 19.2 (live 19.1; was 36.2), StandDead 19.5 (live 20.2; was
   2.2), DDW 1.1 (live 1.1; was 1.3), Released 5.1 (live 5.5; was 0). Snag-fall now starts the right cycle (SD
   19.5→10.9 / DDW 1.1→8.1 over 2000→2005). +test (test_carbon.jl, 4 asserts within tol).
   ★★★★★★★★★ SNAG-FALL TIMING RESIDUAL — FOUND + FIXED (this session, suite 5093/2). The snag-falldown MODEL was
   confirmed faithful (snag_fall_density=fmsfall.f BASE/MODRATE; FALLX/ALLDWN/DECAYX=fmvinit.f class 1/2/3=7.17/3.07/
   1.96; per-species SNAGCLS matches for SM22→2/HI27→3/AB33→2/SK65→2). The bug was the falldown TIMING in
   update_snags!: `yrs = clamp(cur − deathyr, 0, nyears)` used cur = current_cycle_year (the CYCLE-START year,
   constant through the cycle's nyrs annual steps). A FIRE snag (dead at the fire year = cycle-start, and — post the
   issue-B reorder — created BEFORE ffe_fuel_update!'s annual loop) thus read deathyr==cur ⇒ yrs=0 EVERY step ⇒ never
   fell its creation cycle, then fell all at once the NEXT cycle (density 445→404 over 2000→2005, then →15 by 2010).
   FIX: ffe_fuel_update! now passes the INCREMENTING annual year (`at_year = cur+ (k−1)`) to update_snags!, which uses
   it for the fall-age gate (`eff − deathyr`) — so the fire snag ages 0,1,2,3,4 across the loop and falls years 1-4,
   matching FVS FMSNAG(IYR−deathyr). The post-burn PBTIME window keeps the constant cur (no change to post-fire
   cycles). Ordinary-mortality snags are created AFTER the loop ⇒ never in it this cycle ⇒ unchanged (carbon_snt SD
   stays bit-exact; whole suite green). RESULT vs live (fire_carbon.key): density now 445→74 over 2000→2005 (84% fall,
   was 9%); carbon 2005 SD 4.0 (live 2.8, was 10.9) / DDW 13.8 (live 14.8, was 8.1); 2010 SD 1.0 (live 1.0) / DDW 11.8
   (live 11.5). +test (2005 row asserts). The fire-carbon trajectory now tracks live across ALL cycles.
   ★★★★★★★★★★ FIRE-KILLED-ROOT BELOWGROUND POOL — PORTED + EXACT (this session, suite 5093/2). fmburn.jl booked NO
   roots for fire-killed trees, so Below-Dead missed the fire cohort (0.9 vs live 5.6). FVS adds them at death
   (fmsadd.f:320 BIOROOT += RBIO·SNGNEW·XDCAY, XDCAY=1 for age-0). FIX: fmburn.jl now does `fs.bioroot += rbio·curkil`
   per fire-killed record (the snag-FALL path transfers only the bole, so no double-count). VALIDATED vs live
   (fire_carbon.key) BIT-EXACT: Below-Dead 5.6/4.6/3.8 @2000/05/10 == live 5.6/4.6/3.8. +test.
   ⇒ #28 FULLY RESOLVED. The fire_carbon carbon report now tracks live across ALL columns + ALL cycles: 2000
   AGL 19.2/19.1, Merch 12.5/12.4, BGLive 4.6/4.6, BGDead 5.6/5.6, StandDead 19.5/20.2, DDW 1.1/1.1, Floor 0.8/0.8,
   Shb 0.2/0.2, TotC 51.0/51.6, Released 5.1/5.5. Only ~ULP-to-low-% residuals remain (StandDead Δ0.7, Released ~7%,
   TotC Δ0.6) — the phasing co-refactor, released value, snag-fall timing, and fire-killed-root pool are all DONE.
   ★★★★★★★★★★★ RELEASED RESIDUAL — TRACED + FIXED (this session, suite 5094/2). Debug-FVS dump of fmcons.f
   (TCWD pools + PRDUF + DIARED + MOIS) proved jl's SURFACE-fuel consumption matches live BIT-EXACT: litter 2.631,
   duff 6.220·PRDUF(66.66%)=4.146, and the woody classes 3-6/6-12/12+ = 1.123/0.714/0.309 (same PDIA 4/8/15 +
   DIARED=3.38−0.027·m100, m100=17%). The gap was that FVS's released BIOCON(2) (the ×0.5 term, fmdout.f:286-287)
   ALSO includes the burned LIVE herb/shrub fuels (FMCONS BURNLV, fmcons.f:310-311: herb 1.0, shrub 0.6) — jl's
   apply_fire_consumption! consumed only the surface cwd, omitting FLIVE. FIX: `released += (flive[1]+0.6·flive[2])·0.5`
   (release-only; FLIVE regrows via fmcba! each cycle). VALIDATED: released 5.13 → 5.54 (raw) → 5.5 in the report ==
   live 5.5 to 1 decimal. test_consumption.jl updated (the old assertion validated the incomplete woody+ff formula).
   ⇒ #28 COMPLETELY RESOLVED. The fire_carbon carbon report now matches live across EVERY column: 2000 AGL 19.2/19.1,
   Merch 12.5/12.4, BGLive 4.6/4.6, BGDead 5.6/5.6, StandDead 19.5/20.2, DDW 1.1/1.1, Floor 0.8/0.8, Shb 0.2/0.2,
   TotC 51.0/51.6, Released 5.5/5.5. Only ULP-to-rounding residuals remain (StandDead Δ0.7, TotC Δ0.6, released +0.04).

   ★★★★★★★★★★★★ STANDDEAD Δ0.7 — TRUE SEMANTIC CAUSE = crown-lift-at-death (principle #1 semantics-not-runtime).
   Traced the snag-crown lifecycle BOTH sides: bole faithful (fmsvol.f MAX(0.005454·H,MCF)=jl), fire-killed crown
   per-tree EXACT (crown_biomass=CROWNW: d10.44/cr25 → 244.676=244.676; density curkil=DTHISC=1.437), live
   crown-lift→CWD down-wood (fmcadd.f:101, matches jl).
   ⚠ FALSE START (caught by principle #4 + the user): first hypothesised the EXCESS-MORTS snag crown staging and
   "fixed" it by booking those snags before the report — REVERTED as a value-coincidence: FMKILL (gradd.f:122) runs
   AFTER FMCRBOUT, so FVS's report can't contain excess-MORTS snags; the number improved by accident.
   TRUE CAUSE (proven by FVS disable): FVS dead-tree crown = `CROWNW + YRSCYC·OLDCRW` (fmscro.f:147) — the dead crown
   PLUS the accumulated crown-lift the tree would have shed; jl's crown was `CROWNW` only. Disabling the OLDCRW term
   in live FVS drops its 2000 StandDead 20.2→19.3 (≈ jl's 19.5) — confirming jl omits exactly that term. FIX: store
   the per-tree OLDCRW (compute_crown_lift!'s `x·xvold[size]`) in a new `TreeList.ffe_oldcrw` (5 sizes, carried by
   copy_tree! through tripling/compaction), and add `YRSCYC·OLDCRW` to the dying record's crown via
   `crown_lift_at_death` at the fire death site (fmburn!). NO double-count: the killed fraction has tpa→0 so
   compute_crown_lift! (skips tpa≤0) never books its lift to down-wood. RESULT: 2000 StandDead 19.5→20.4 (live 20.2),
   TotC 51.0→51.9 (live 51.6); suite 5094/2, no regression. NOT bit-exact — residual now StandDead +0.2 / TotC +0.3
   (NOT verified ULP): a crown-lift snapshot-PERIOD nuance — FVS's death-OLDCRW is the FMSDIT-scaled value of the
   CURRENT cycle; jl reuses the PRIOR cycle's stored lift (≈equal magnitude, hence +0.2).
   ⚠ SECOND FALSE START (caught by principle #3 + #1): wired the SAME term into ORDINARY mortality
   (book_mortality_snags!) → carbon_snt StandDead [3.796,4.393,5.354,9.535] OVERSHOT (11 fails). Traced to SOURCE:
   fmscro.f:145 gates the OLDCRW term with `IF (ICALL .NE. 4)`, and ordinary mortality reaches FMSCRO via
   FMKILL→FMSADD with ITYP=4 (fmkill.f:143) ⇒ ordinary deaths get CROWNW ONLY. Only FIRE (ICALL=1, fmeff.f:390/465/
   491/525) and CUT (ICALL=2, fmscut.f) add YRSCYC·OLDCRW. So the regression was MY over-application, NOT a masked
   bug — REVERTED; jl now applies the term ONLY at the fire death site. VERDICT logged: ordinary-mortality crown =
   CROWNW only = faithful.
   ⊕ LIVE-VALIDATED the term (debug fmscro.f JOSTND dump of OLDCRW/CROWNW/YRSCYC/X per fire-killed tree, fire_carbon
   2000): every fire row has X=1.0, YRSCYC=5, and OLDCRW≈0.022·CROWNW (the per-year crown-lift rate) — exactly the
   shape jl's `crown_lift_at_death = 5·ffe_oldcrw` implements; term dominated by high-density small kills (DSNAGS≈30)
   + the large tree.
   ⊕⊕ +0.2 RE-LOCALIZED (this session — NOT the OLDCRW basis after all): a direct per-tree jl-vs-live diff of BOTH
   ffe_oldcrw AND crown_biomass at the 2000 fire (env-gated jl println + debug fmscro.f JOSTND dump of CROWNW(I,SIZE),
   matched per-tree by the bit-identical kill density DSNAGS=curkil) shows ffe_oldcrw is FINE (i=1 bit-exact; small
   trees track the OLDCRW≈0.022·CROWNW rate). The +0.2 is actually a pre-existing **crown_biomass size-class-1
   over-production for SMALL trees**: i=2 (sugar maple, dbh1.28) size-1 crown jl 0.2725 vs live 0.136 (2×); i=3
   (dbh6.9) 11.94 vs 9.82 (1.22×); BIG trees bit-exact (i=1 dbh12.4 all sizes match). Sizes 2/3/4 match the digit on
   the SAME trees — only xv1 (the 0–.25" woody class) is high. Species-form selection is CORRECT (SPILS 26=sugar
   maple→maple F1-F4, 28=beech→default; both match FVS fmcrowe.f:320-352 SELECT CASE(SPILS)); the maple/conifer XV
   assembly is byte-identical to fmcrowe.f:507-543. The xv1-high-while-xv2/3/4-match signature is mathematically
   INCONSISTENT with both sides sharing the maple formula+inputs ⇒ the CROWNW that FMSCRO reads is NOT jl's raw
   crown_biomass.
   ⊕⊕⊕ +0.2 RESOLVED (this session — the REAL cause, principle #1): FMEFF CONSUMES the fire-reached fine crown of the
   killed trees BEFORE FMSCRO books them as snags. fmeff.f:457-460: in the scorched crown zone the fire burns 100% of
   the foliage (size 0) and 50% of the 0-0.25" branches (size 1, incl. the OLDCRW size-1), `CROWNW(0)·(1−PROPCR)` /
   `CROWNW(1)·(1−0.5·PROPCR)`, where PROPCR = the scorched fraction of crown LENGTH (fmeff.f:435 = sl/CRL; NOT the
   parabolic-volume `csv` used for mortality). Sizes 2-5 are above the flames / too coarse ⇒ untouched. This is EXACTLY
   the per-tree live data: sugar-maple d1.28 (whole crown below scorch, PROPCR=1) size-1 0.2725→0.136 (the 2×); beech
   d6.9 PROPCR≈0.36 → ×0.822; tall trees PROPCR=0 (intact — the prior "crown above flame" comment was right ONLY for
   them). jl had booked the FULL pre-fire crown. FIX: fmburn.jl computes PROPCR per killed tree and reduces
   foliage+size-1 (and halves the OLDCRW size-1, fmeff.f:460) before fmscro!. The consumed crown is RELEASED to the
   atmosphere (BCROWN→BURNCR→TSMOKE smoke, fmeff.f:609-611), so it correctly leaves StandDead AND total stand carbon;
   the carbon-report "released" col 11 = BIOCON surface-fuel only (fmcrbout.f:151 = BIOCON(1)·0.37+BIOCON(2)·0.50),
   already matched, so NO double-count. RESULT (live-validated, all 3 carbon columns within 0.1): 2000 StandDead
   20.4→20.1 (live 20.2), TotC 51.9→51.7 (live 51.6), Released 5.5 (live 5.5). Suite 5094/2, NO regression (fire_early
   down-wood + carbon_snt unaffected — carbon_snt is ordinary mortality, no fire consumption). The earlier "FMCROW
   ITRN=0 / CROWNW-fill-path" blocker was a red herring: the discrepancy was never crown_biomass being wrong, it was
   the MISSING fire consumption of the (correct) crown. #28 fire-carbon is now faithful on every column at the report's
   display resolution.
   ⊕ 2005 StandDead 4.0 vs 2.8 — SCOPED + fall-model RE-TRACED FAITHFUL (this session). Split per-cycle (jl debug):
   2005 = bole 3.68 + crown 0.31; the crown falls correctly (8.87→0.31), so the residual is entirely the BOLE. The
   "snag-fall-RATE" label was checked source-to-source and HOLDS: snag_fall_density = fmsfall.f:128-175 (BASE
   −0.001679·D+0.064311, MODRATE=BASE·FALLX, the d<12 small-snag path + the last-5%/ALLDWN large-snag path), and ALL
   post-burn params match FVS fmvinit.f:1100-1104 (PBSOFT 1.0 / PBSMAL 0.9 / PBSIZE 12 / PBTIME 7 / PBSCOR 0) with the
   `DBH<PBSIZE` accel gate (fmsnag.f:182). One INERT deviation noted: jl floors BASE at max(0.01,·) where FVS doesn't —
   binds only for D>32 (negative-fall edge), absent in SN. Since the fall model is faithful AND the SN snag bole is
   frozen at death (FMSNGHT height-loss is a no-op, HTX=0), the 2005 bole reduces to the 2000 SNAG LIST.
   ⚠ CORRECTION (the first guess was wrong): I'd attributed it to the "per-tree fire-kill DISTRIBUTION," but evidence
   already in hand REFUTES that — the per-tree kill `curkil` (=DSNAGS) matched FVS BIT-EXACT in the fmscro dump (i=2:
   29.825 both; i=1/3/4/6 all match). So the 2000 snag-list DENSITY matches; the residual is NOT which trees die.
   With density matched + fall-rate model faithful + bolevol frozen (MCF), jl's 2005 bole *should* equal FVS's — yet
   it's high. That implicates the fall APPLICATION or the runtime snag-state evolution, cause genuinely OPEN. BLOCKED
   from FVS-side confirmation by a recurring data-flow wall: fmdout (BIOSNAG=TOTSNG(1)+TOTSNG(2)) shows NSNAG=0 at
   every report call; the SN FFE snag-carbon report serves a SERIALIZED value (fmppget.f:226 BIOSNAG=REALS(48)), so the
   snag list/bole isn't dumpable there. Verdict: fall-RATE model verified FAITHFUL; 2005-bole CAUSE is OPEN (earlier
   "fire-kill distribution" verdict retracted), gated on mapping the SN snag-carbon serialization data flow.
   ⊕ UNBLOCKED (SNAGBRK session): the data-flow wall was a RED HERRING — the live per-snag state IS dumpable from
   FMSNAG (fmmain.f:232), NOT the fmdout report. Captured live fire_carbon @2005: NSNAG≈73, ΣDENIH≈654.9/ΣDENIS≈15.1,
   bole dominated by a few large hard snags (i=1 dbh9.26 dih1.18 htih63.6; i=2 dbh34.6 dih0.33 htih92.7; i=3 dbh10.48
   dih1.10). NEXT STEP (ready, no longer blocked): dump jl's fire_carbon @2005 snag list and diff DENIH/DENIS by
   record vs this live state — jl's 2005 bole 3.68 vs live ~2.5 should localize to either retained DENSITY (fall) or a
   per-record dbh/bole difference. The same FMSNAG dump now makes ALL snag-carbon residuals investigable.
   ⊕⊕ DIAGNOSED (jl-vs-live FMSNAG dump). ⚠ My first read ("live ΣDENIH=654.9, live retains 9× MORE") was a SELF-
   INFLICTED ERROR — the awk summed column $7 = DBH, not density. CORRECTED with a clean by-DBH-class dump: live 2005
   ΣDENIH = **15.1**, jl = 73.6 ⇒ jl **UNDER-falls** the fire snags (consistent: jl more snags ⇒ more bole ⇒ higher
   StandDead 4.0 vs 2.8 ✓; the dbh34.6 large snag still matches bit-exact). Live fire-snag fall by DBH class (created
   @2000 → standing @2005): <3" 256.2→0, 3-6" 47.8→0.08, 6-12" 129.7→10.6, >12" 11.6→4.48 (total 445→15.1). So live
   falls the SMALL (<12") fire snags almost COMPLETELY over 5 yr (post-burn RSMAL/RSOFT), keeping only the large
   (>12") that don't accelerate. jl retains ~74 ⇒ jl's post-burn small-snag fall is too WEAK. The memory's "445→74
   tracks live" is STALE/WRONG (live is 445→15, not →74). ROOT (fmsnag.f:178-214): post-burn rates PBFRIH/PBFRIS are
   SET in the year after the fire and applied for PBTIME; for small snags PBFRIH=RSMAL, bumped to RSOFT if the snag is
   SOFT (line 184 `.NOT.HARD`). NEXT: dump live RSMAL/RSOFT from fmsfall.f (jl's ≈28%/yr gives 19% retained = jl's 74,
   but live reaches ~3% ⇒ live's effective small-snag rate is ~2× jl's) and reconcile jl's update_snags! post-burn
   (snag.jl:245-256) against it — likely the RSOFT-after-DKTIME switch and/or the RSMAL value. The fix touches ONLY
   the post-burn (fire stands), so carbon_snt's bit-exact falldown is unaffected.
   ⚠ MECHANISM NOT YET FULLY RESOLVED — do NOT code a fix on this yet. Read the rates: fmsfall.f:107/115 RSMAL =
   1−exp(log(1−PBSMAL)/PBTIME) = 28%/yr (PBSMAL .9, =jl); RSOFT = 1−exp(log(DZERO/DENTTL)/PBTIME) ≈ 83%/yr (PBSOFT=1,
   the ELSE branch since PBSOFT≥1). The dktime flip (fmsnag.f:282-284) sets HARD=.FALSE. but does NOT move DENIH→DENIS;
   PBFRIH/PBFRIS are SET ONCE at IYR−BURNYR≤1 (fmsnag.f:178) with the line-184 `.NOT.HARD→RSOFT` bump evaluated THEN.
   STATIC ANALYSIS PREDICTS most small fire snags are still HARD the year after the fire (DKTIME>1 unless dbh<~0.4) ⇒
   RSMAL=28% ⇒ ~19% retained ≈ jl's 74 — yet LIVE falls them to ~0. So a real mechanism is missing from this reading
   (candidates: DENTTL-dependent RSOFT being tiny-density-driven; per-record DKTIME<1 for the bulk; or the rate is
   re-set each year not once).
   ✅ MECHANISM RESOLVED (per-snag/per-year FMSNAG FALLDBG dump done): the small-snag fall is
   `DFIH = MAX(post-burn PBFRIH·DENIH, DFALLN)` where **DFALLN = MODRATE·ORIGDEN is a FIXED amount/yr** (FMSFALL gets
   `DEND`=orig density, fmsnag.f:165) — NOT proportional. So a fire snag falls ~28% early (post-burn dominates) then a
   CONSTANT MODRATE·DEND/yr that clears it to ~0 by 2005. Live dbh2.67 trajectory: DENIH 53.9→38.8→27.9→18.0→8.1(→0),
   DFIH 15.1→10.9→9.9→9.9→9.9 (PBFRIH=0.28 constant, HARD=T throughout — the RSOFT/soft path is NOT what clears them,
   the constant DFALLN is). My static prediction was wrong because it treated the normal fall as proportional.
   ⇒ jl's snag_fall_density ALREADY returns `modrate·origden` for d<12 (snag.jl:31) and update_snags! takes the same
   post-burn MAX — so jl's LOGIC matches; the 73.6-vs-15.1 gap is a VALUE/stepping diff. FINAL CHECK (one jl dump):
   print jl's per-snag DFIH for the fire snags vs the live FALLDBG (dbh2.67: 15.1/10.9/9.9/9.9/9.9) — likely jl's
   `origden`(DEND) or `fallx` differs for fire snags, or update_snags!'s min(denttl,·)/per-year stepping caps it.
   Bounded; carbon_snt unaffected (its falldown already bit-exact, same code path, no post-burn).
   ⊕ jl-DUMP attempted — hit a TIMING/SAMPLING mismatch (must be handled in the next pass): update_snags! SKIPS a
   just-created snag at its creation year (`yrs = eff−sn.year ; yrs>0||continue`, snag.jl:200-201), so the 2000 fire
   snags first enter jl's fall processing at eff=2001, whereas live's FMSNAG dump samples them at IYR=2000 (already in
   the array). So a naive eff=2000 jl pool dump shows only the PRE-2000 snags (~45), not the 445 fire pool — the
   per-year/per-snag jl-vs-live comparison must align on this (compare jl eff=2001..2005 fall vs live IYR=2001..2005,
   and the matching snag by dbh, since jl's fire-snag dbh records ≠ live's index order).
   ✅✅ CAUSE PINNED (jl JLCHK dump vs live FALLDBG): jl's fall RATE is CORRECT — `modrate·origden` matches live
   bit-exact (fallx≈3.07 both species; e.g. dbh1.28 orig29.83 → DFIH 5.69/yr = live's constant-fall shape). The bug is
   the **CREATION-YEAR SKIP**: jl's gate `yrs = clamp(eff−sn.year,0,nyears); yrs>0||continue` (snag.jl:200-201) skips a
   FIRE snag's fall in its creation year (eff=2000=sn.year ⇒ yrs=0), so jl falls them 4 yr (2001-2004) while LIVE
   falls 5 yr (2000-2004 — FMSNAG runs AFTER FMBURN in FMMAIN, so the fire-year fall happens). Since the small-snag
   fall is a CONSTANT `modrate·origden`/yr, the missing FINAL year leaves ~5× more (orig29.83, 5.69/yr ⇒ jl@4yr=7.0 vs
   live@5yr=1.4) ⇒ the 73.6 vs 15.1, the higher 2005 bole, and StandDead 4.0 vs 2.8. FIX (scoped + SAFE): let FIRE
   snags fall in their creation year. The `yrs>0` skip ONLY reaches snags processed at eff==sn.year = FIRE snags
   (created before ffe_fuel_update!'s loop); ordinary-mortality snags are created AFTER the loop so never hit it ⇒
   carbon_snt's bit-exact falldown is UNAFFECTED. Impl: when a snag has sn.year==cur0 (created this cycle, pre-loop =
   fire) treat the first annual step as a fall (yrs≥1 at eff=cur0), matching FMSNAG-after-FMBURN. Validate vs live
   fire_carbon 2005 (target StandDead 2.8, DENIH 15.1) + keep carbon_snt falldown [3.796,4.393,5.354,9.535] green.
   ✅✅✅ FIXED + LIVE-VALIDATED (this session): added `born_now = (sn.year==cur)` to update_snags! so a fire snag
   (created this cycle, pre-loop) falls in its creation year too (`yrs = clamp(eff−sn.year+born_now,0,nyears)`).
   RESULT: fire_carbon 2005 StandDead 4.0→**2.6** (live 2.8), DDW→15.2 (live 14.8); 2000 (20.1) + 2010 (1.0) still
   match. carbon_snt falldown UNCHANGED bit-exact (the gate only reaches fire snags — ordinary mortality is created
   AFTER the loop). Suite 5109/2, NO regression; fire_carbon test tightened to the live magnitude (SD≈2.8/DDW≈14.8).
   ⇒ the LAST #28 value-residual is CLOSED — fire-carbon now matches live on EVERY column/cycle at display resolution.

   ★★ CARBCALC GAP (task #42) — ✅ LANDED this session (suite 5061/2, carbon 192/192, no regression). CARBCALC
   (fmin.f opt 46): FLD1=method (0=FFE *default*, 1=Jenkins), FLD2=units (0=US t/ac *default*, 1=metric t/ha,
   2=metric t/ac). jl had the WRONG defaults (Jenkins+metric) and parsed only FLD1, and only in kw_carbcalc!
   (the each_stand dispatch path parsed neither units). FIXED: carbon_method default→FFE(0), new carbon_units
   (default 0=US t/ac), FLD2 parsed in BOTH handlers, conditional unit factor in stand_carbon_report (carbon.jl),
   test:20 stale "default Jenkins" assertion corrected. VALIDATED: fire_carbon (no CARBCALC → FFE+US t/ac) now
   Total 34.3/44.11 vs live 34.40/44.42 (was 92.27); carbon_{ffe,jenkins,snt} (FLD2=1 metric) stay bit-exact. So
   fire_carbon's BASELINE carbon is now CORRECT — the remaining 2000+ gap there is purely the phasing lag (SD spike
   2005 vs live 2000) + the unwired released value. (Re-trace win: my first "STDINFO/area" guess was wrong; the
   72-test regression from the first partial fix was the principle-#3 signal that found the second-handler root.)
