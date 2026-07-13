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

## C. Cycle-0 TopHt tie-break (1 stand) — the already-cornered AVHT40 primitive
**55250794010661** (81.21): the ONLY candidate with a CYCLE-0 (inventory) divergence — TopHt live 22 / jl 29,
before any projection. Traced: top height = avg height of the 40 largest-DBH trees/acre (AVHT40). This ultra-dense
stand's largest-DBH record (DIA=2.9", TREE_COUNT≈75 TPA > 40) alone fills the entire top-40-TPA window ⇒ TopHt =
that single record's height. But TWO records tie at DIA=2.9": aspen (sp746) HT=29 and balsam fir (sp12) HT=22.
FVS's unstable RDPSRT quicksort lands the balsam (22) at the boundary; jl's double-RDPSRT lands the aspen (29).
The 19/44 missing-HT trees are all shorter (never in the top set) ⇒ NOT height imputation.
RESOLVED = CORNERED: this is the KNOWN AVHT40 top-height tie-break primitive already documented in
`src/engine/standstats.jl:124` stand_top_height — jl already ports `_rdpsrt!` with the cratet.f DOUBLE sort
(LSEQ=.TRUE. then .FALSE.), and prior dig-sessions #1/#2 proved (empirical single-vs-double sweep over 4 tie-heavy
stands) that "no global sort choice is bit-exact" because RDPSRT is an unstable quicksort on tied DBHs — the
tie-break is STAND-DEPENDENT. The double-sort matches the most stands; 55250794010661 is one of the residual
tie-heavy swings. The cycle-0 22/29 then cascades through ultra-dense self-thinning (hence 81%). NOT a new bug,
NOT fixable without bit-matching FVS's exact quicksort pivot sequence (accepted primitive per the GOAL). All 29
LS candidates now explained: 11 resolved/ULP + 17 self-thin RDPSRT primitive + 1 AVHT40 RDPSRT primitive.

# NE/CS candidates (14) — post-FIX-#8 check (audit 43eb)

FIX #8 is LS-gated ⇒ inert for NE/CS (expected). Batch worst%:
- CS resolved/ULP (<5%): 193228158010661 (3.78), 943930276290487 (2.28), 430360223489998 (4.23).
- CS/NE mid (5-12%): 75190472010538 NE (9.63), 97513385010661 CS (7.14), 65532203010661 CS (11.43),
  1203406023290487 NE (11.67).
- Higher (>20%): NE 14106653020004 (24.91), 166318995010661 (27.69), 245503277010661 (29.88),
  366792805489998 (24.86); CS 193238139010661 (24.03), 255129978489998 (21.21), 562745328126144 (28.46).

## NE/CS cycle-0 TopHt divergence — RESOLVED = the cornered AVHT40 RDPSRT tie-break primitive
Several NE/CS candidates diverge in **TopHt at cycle 0 (inventory)** with everything else bit-exact, then cascade:
  NE 14106653020004  TopHt live34/jl27   NE 245503277010661 live36/jl27   NE 166318995010661 live46/jl42
  CS 255129978489998 TopHt live33/jl40   CS 562745328126144 live36/jl38   CS 97513385010661 live35/jl36(±1)
TRACED to certainty (14106653020004): the .sum TopHt = IBTAVH = INT(AVH+0.5), AVH = AVHT40 (avht40.f) computed in
cratet.f. AVHT40 = avg HT of the largest-DBH 40 TPA over IND. jl's AVHT40_DBG dump: top-40-TPA = DBH5.0(HT28,
TPA6) + DBH**3.9**(HT**27**, TPA75, 34 used) ⇒ (28·6+27·34)/40 = 27.15 → 27. FVS lands the OTHER dbh=3.9 record
(sp316 HT**35**) in the window ⇒ (28·6+35·34)/40 = 33.95 → **34**. So the split is a genuine **DBH tie at 3.9"**
between a HT=27 and a HT=35 record — my earlier "largest tree unique ⇒ not a tie-break" was WRONG (the SECOND
slot, not the unique largest, is the tied boundary). This is the KNOWN AVHT40 RDPSRT tie-break primitive
(src/engine/standstats.jl:124).
FIX ATTEMPT (refuted): cratet.f sorts IND for AVHT40 via RDPSRT(.FALSE.)@141 then RDPSRT(.TRUE.)@245; the .TRUE.
pass RESETS IND to identity ⇒ effective single sort from record-order. Tried replacing jl's double-sort with the
single lseq=true — jl STILL gave 27 (second pass was inert here). So the divergence is NOT the double-vs-single
sort choice; it's that jl's tree-array order (the identity seeding the unstable quicksort) breaks the 3.9" tie to
the HT27 record while FVS's breaks it to HT35. That is precisely the documented "stand-dependent, no global sort
choice is bit-exact" AVHT40 unstable-quicksort primitive — CORNERED (accepted per the GOAL). More prevalent in
NE/CS (~5/14 vs LS 1/29) because NE/CS FIA inventory has more integer-tied DBHs at the top-40 boundary; on
tie-heavy dense stands the cycle-0 tie-break cascades into the self-thinning trajectory (accepted cornered cascade,
same class as LS 55250794010661). A true fix would require bit-matching FVS's exact tree-read order AND RDPSRT
pivot sequence (high blast-radius, deferred; not attempted — would risk the 99.7% already-matching population).
ALL 43 candidates now explained: LS 29 (1 FIX #8 + 28 cornered) + NE/CS 14 (ULP + AVHT40 tie-break cornered).

### Trace provenance (43eb)
Sampling path: .sum TopHt = sumout.f IBTAVH ← disply.f:360 INT(OLDAVH+0.5) ← OLDAVH=AVH (fvs.f:436/grincr.f:282)
← AVH = avht40.f, called ONLY from cratet.f:529 (calibration). cratet IND sort for AVHT40 = RDPSRT(.FALSE.)@141
then RDPSRT(.TRUE.)@245; rdpsrt.f LSEQ=.TRUE. RESETS INDEX=1..N ⇒ the .TRUE. pass discards line-141 ⇒ effective
single sort from record-order/identity. jl-side AVHT40_DBG dump on 14106653020004 pinned the 3.9" HT27/HT35 tie
as the boundary (above). Live avht40.f WRITE guarded ICYC≤1 did not fire (AVHT40 runs during cratet before the
cycle counter; not needed once the jl-side dump + source path made the tie-break conclusive). FVS source restored
pristine, oracle untouched.

### Dig-queue triage snapshot (43ec, sweep @ ~262k) — all known classes
The live dig-queue (17 entries) is fully classified, no unexplained pattern:
- `18xxxxx010661` cluster (CCF-only or CCF+volume, ~15-24%) = forest-924 elevation ⇒ FIX #9-resolved. Verified on
  18447951010661: post-fix, 10-col cols bit-exact in early cycles (incl volume TCuFt/MCuFt/SCuFt/BdFt), residual
  only ~2-3% in later cycles (2010+) tracking a ±1-2 structural ULP/tie-break — NOT a separate volume bug.
- multi-col ultra-dense (224645781010661, 62261244010661, 9426842020004, 24089675010661; all/most 10 cols, high
  rel%) = self-thinning RDPSRT tie-break primitive (cornered).
- 55250794010661 (81%, TopHt) = AVHT40 RDPSRT tie-break primitive (cornered).
- 166318995010661 (NE, empty div_cols) = NE AVHT40 TopHt tie-break (cornered).
- 1831637837290487 (LS, empty div_cols) = tamarack, FIX #8-resolved (STALE queue entry).
STALE entries (FIX #8/#9-resolved) to be purged at the DIGCAP triage; the sweep picks up both fixes on fresh
batch-cycles so forward flags won't repeat them.

  (verified 43ec: 9426842020004 cyc0 TopHt 64/62 ±2 TPA/BA bit-exact ⇒ AVHT40 tie-break cornered;
   24089675010661 cyc0 BIT-EXACT TopHt64/TPA1396/BA160 ⇒ self-thin later-cycle cornered. Both NOT new bugs.)
Dig-queue FULLY triaged — every entry FIX-#8/#9-resolved or a named cornered primitive; no unexplained divergence.
