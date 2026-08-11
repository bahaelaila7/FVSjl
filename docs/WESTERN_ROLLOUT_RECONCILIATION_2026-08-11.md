# Western Rollout Reconciliation — 2026-08-11

Reconciles the ACTIVE-GOAL stop-hook doc (which has gone stale) against the measured/committed state after this
session. Method throughout: instrument the live `FVS{v}_g16` oracle, root-cause each divergence, fix-and-validate
what's fixable, corner what's stochastic, mark what's blocked — and measure BEFORE fixing so wrong fixes never land.

## Fixed + validated this session (commits on `kt-variant-port`)
- **Dwarf-mistletoe MORTALITY wiring** (`3d4144e`) — `ie_dm_mortality_combine!` (mismrt.f) was wired only in the
  shared `southern/mortality.jl`; the 6 western variants with their own Hamilton-RIP `mortality!`
  (CI/IE/KT/EM/UT/TT) never applied it ⇒ DM-infected trees never died. THIS is the root of the cross-variant
  "mature low-density over-growth" headline (CI 114% / TT 33% / EM 29% / UT 16% max|ΔBA|). Validated: CI
  753180709290487 now declines BA 52→42 bit-exact vs live (was growing 52→66).
- **Per-variant mistletoe tables** (`0a39576` UT, `6b048e3` EM/KT/BM/TT/CI) — misint{v}.f AFIT/ADGP/APMC are
  per-variant; jl applied IE's table to all. DF(3)/LP(7)/PP(10) align (so DF/LP stands validated) but each
  variant's non-aligned hosts (UT pinyon, EM WP/GF, …) were mis-mapped. Validated on 2 non-aligned hosts (UT
  pinyon ~bit-exact; BM larch tracks live).
- **FIXMORT** (`1462434`) — same "southern-only wiring" class: `apply_fixmort!` (morts.f:781) had one call site
  (southern). The 6 western variants silently IGNORED the FIXMORT keyword. Fixed + validated on CI (90% kill
  115→12 = live →11) AND TT (Zeide/ttmrt: 8107→811 = live →759). Inert-guarded ⇒ non-FIXMORT stands byte-identical.
- **Climate-FVS** (`3e0bb1c`→`230d4e7`, 7 commits, earlier this session) — clmorts/clgmult/GrowMult/MortMult/
  clauestb/clmaxden/file-CLIMDATA ported + validated bit-exact-or-cornered on IE.

## Stop-hook goal doc is STALE — actual status of its "remaining work"
- **Dwarf mistletoe "ALL western DONE ✓"** — was WRONG (validated DG-loss+spread, not the KILL); now genuinely
  done after the 3 mistletoe commits above.
- **Climate-FVS "TODO"** — DONE this session (7 commits, IE-validated).
- **#140 BM under-thinning** — ★★★ ACTUALLY FIXED 2026-08-11 (2c26eca); the goal doc was RIGHT and this doc's
  earlier "RESOLVED e130546 / one-stand corroborated" was WRONG (over-optimistic single-stand check — the exact
  "multi-stand tally buries a real bug" trap the mistletoe fix warned of). A 14-stand BM FIA sign-tally showed a
  SYSTEMATIC 11-HIGH/1-LOW/2-BE under-thin bias (mean +8%). ROOT (measured via FVSbm_g16 per-tree, bmt01 icyc1):
  the DDS prediction is BIT-EXACT (jl/live=1.0000) but DG=sqrt(D_ib²+DDS)−D_ib used the WRONG bark — the shared
  diameter_growth! DDS→DG apply loop dispatched POWER bark for CR/TT/BC but BM (also POWER, bm_bratio) fell through
  to the LINEAR bark_ratio (~0.99 vs correct ~0.86) ⇒ D_ib too large ⇒ DG 7-8% LOW on EVERY tree ⇒ dq10 low ⇒
  self-thin under-kill, feedback-amplified. Fix = add the `_bm_dg ? bm_bratio` branch. Post-fix: deterministic DG
  bit-exact; bmt01 2090 TPA 108→94 vs live 96; tally → 6-HIGH/3-LOW/5-BE (mean +0.6%, STRADDLES). Residual = the
  cornered DGSCOR straddle. ⇒ #140 truly resolved.
- **BARK-DISPATCH AUDIT COMPLETE 2026-08-11** — the BM bug was a MISSING-VARIANT branch class; swept the whole
  DDS→DG bark dispatch: POWER-bark variants (CR/TT/BC/BM/CI) MUST have a special-function branch — BM (2c26eca) and
  CI (6d94b6a) were the two missing (CI minor, ~0.3% net, .sum-inert on cit01 but bit-exact-DG faithful; CI #142
  confirmed a cornered count-straddle, 2-HIGH/5-LOW/3-BE, NOT a bias). Linear/reciprocal variants (EM/IE/UT/KT)
  correctly ENCODE their bark into c.bark_a/c.bark_b (incl. the 0.9002−0.3089/D zero-coef case) so the shared
  bark_ratio is faithful — verified in each variant's diameter_growth.jl AND now MEASUREMENT-CONFIRMED per-tree vs
  FVS{v}_g16: EM (constant BARK1 form, emt01) and UT (reciprocal BARK1+BARK2/D form, utt01) both have deterministic
  DG jl/live = 1.0000 (0 non-det trees), bark matching per-tree; IE/KT use EM's verified constant form. ⇒ all three
  bark forms (POWER/constant/reciprocal) verified; dispatch complete; no remaining missing-bark branches. LESSON:
  the earlier audit's "DDS bit-exact ⇒ DG faithful" was the trap — DDS≠DG when the bark differs; measure the
  bark-converted DG, not just DDS.

## EM tail re-checked with the proven method — genuinely CORNERED (2026-08-11)
Applied the BM lens to EM (goal doc flags EM #137 open + "EM/IE ~7% tail cornered" on the discredited aggregate
reasoning): (1) per-tree deterministic DG on emt01 vs FVSem_g16 = BIT-EXACT (jl/live 1.0000) ⇒ NO bark/DG bias like
BM; (2) 12-stand EM FIA sign-tally = 5-HIGH/5-LOW/2-BE, 0 crashes ⇒ BALANCED straddle, not a systematic bias. The
larger per-stand magnitude (±7-18%) is the DGSCOR RNG realization amplified by self-thin feedback (+ AUTOES variance),
straddling ~0. ⇒ EM is genuinely bit-exact-or-cornered; the "cornered" verdict is CORRECT — now confirmed by the
sign-tally + deterministic-DG measurement, not the aggregate-BA reasoning that had masked BM #140. EM #137
(dense-cohort self-thin) was already fixed (c7c7d2f + 04b15e6/7ce8f1f tem 35000-cap) per memory; goal doc stale.
- **#143 IE AUTOES** — ✅ FIXED 2026-08-11 (9cc7d7a): d089b78 fixed single-ingrowth stands but MULTI-ingrowth-cycle
  compounding still over-established up to +108% (confirmed real+systematic, 9H/2L tally) — ROOT was jl failing to
  parse STRING plant-association PV_CODEs; fixed by porting ie/habtyp.f's crosswalk. Post-fix tally 6H/4L/2BE
  (straddles, ±4%) ⇒ IE now bit-exact-or-cornered like EM. See the dedicated IE #143 resolution section below.
- **#142 CI tail / EM-IE growth tail** — CORNERED (DGSCOR/RDPSRT RNG straddles), meets the bar.
- **#137 EM estab / EM AUTOES over-establishment** — RESOLVED. The AUTOES "+63-86% TPA over-establishment" flagged
  as "likely systematic EM/IE" was fixed by d089b78 (ingrowth NSTORE) — CONFIRMED this covers EM: simulate.jl:578
  routes BOTH IE and EM through the SAME shared `ie_autoes_establish!` the fix touched. Corroborated by the
  2026-08-11 corpus sweep (EM 8/10 bit-exact-or-cornered — inconsistent with a systematic over-establishment,
  which would blow up most stands). The dense-cohort self-thin piece was separately fixed by 04b15e6/7ce8f1f
  (SDI-gate tem 35000-cap). ⇒ EM is bit-exact-or-cornered; goal-doc "still open" is stale.

## TT #158 small-tree-regent — RESOLVED 2026-08-11 (f66f1fd), was mis-classified BLOCKED
The 33% dense-stand over-growth was NOT blocked on an un-portable inline model. ROOT (MEASURED via FVStt_g16 on
dense FIA stand 1629318558290487): jl capped per-cycle small-tree DG at the RAW per-year DGMAX (tt/regent.f:175-177,
0.2 for conifers) but live's regent.f:684 applies `IF(TTVAR)DGMX=FINT*DGMAX` — the FINT(=10) multiplier was missing,
so jl clamped correct 0.6-1.0"/cycle sub-1" DG down to 0.2 → flat DG → low Reineke DR10 → self-thin never fired.
jl's UNCAPPED SMDGF dgr already matches live's DG bit-close per-tree (i34: DK 1.789/1.798, dgr 0.637/0.646) ⇒ the
earlier "smdgf_ called 0× ⇒ different single-step model" read was the compiler INLINING smdgf (no CALL), not a
different model; the earlier "POTHTG ABI un-derivable / regent un-instrumentable (SIGFPE)" blockers were both false
(FVStt_g16==clean and instruments the regent path cleanly). FIX: cap at FINT·DGMAX_RAW at the default + esgent sites.
Validation: #158 stand final BA 33%→7% (jl now completes 5 cycles tracking live's self-thin, was stalling at 2);
ttt01 +1.7%→+2.6% BA (~1% nudge, both within the pre-existing cornered straddle — ttt01 was never bit-exact); 8-stand
TT FIA sweep inert (±2 BA) + 0 crashes, pre-existing per-stand over-growths unchanged (no regression). Residual on
dense stands (~7% BA, self-thin ~1 cycle late) = the small systematic SMDGF-vs-live-inline DG realization difference
(~1-2%/tree), CORNERED. ⇒ TT small-tree growth is now bit-exact-or-cornered cluster-wide.
SCOPE NUANCE (measured pre/post on 6 dense TT stands): the fix is SURGICAL — it changes ONLY stands with a VIGOROUS
sub-1" cohort whose regent DG exceeds 0.2 (like #158, DG 0.6-1.0). On other dense stands (avg-diam 1.3-2.9 but slower
sub-1" DG <0.2) the cap never bound, so the fix is byte-inert there — and those stands' pre-existing +2-7% BA
over-growth (e.g. 275463546489998, 1629325863290487 — identical PRE/POST) is a SEPARATE residual (the same
SMDGF-vs-inline straddle class), NOT resolved by #158 and NOT introduced by it. Do not re-attribute that tail to #158.

## Remaining — correctly classified, none a low-risk quick win
- **MORTMSB** (zero-practical-value; investigated, prototyped, reverted) — mature-stand-breakup keyword. SCOPE
  CORRECTED by measurement: only EM/UT/TT have the inline morts.f MSB block; CI/IE/KT have NO MSB (live ignores it —
  PROVEN byte-identical — so jl ignoring it is FAITHFUL; adding it would have regressed them). Prototyped the UT port
  and found it TRACTABLE (faithful-by-composition: const_/msb_d10=dq10/bark all align with southern's validated
  block + shared _msbmrt!) — BUT it NEVER fires on real UT FIA data (dense stands run low-QMD from seedling cohorts;
  mature stands are sparse) ⇒ jl-ignoring-MORTMSB == live on ALL real UT stands, zero real-FIA effect, and the firing
  path can't be live-validated (no triggerable stand; SN mortmsb.tre breaks under UT on SN site/habitat). REVERTED
  the prototype — committing unvalidated firing logic for a never-fires mechanism isn't disciplined. A future session
  may commit the faithful-by-composition EM/UT/TT ports for completeness; it changes nothing on real FIA.
- **non-IE climate-mort** — `apply_climate_mort!` only wired to IE (Climate-FVS is IE-scoped); extending = separate task.
- **BC merch/board volume, V2/non-ICH** — pre-existing, lower priority.

## Validation
Corpus sweep (jl vs FVS{v}_g16, final-BA classification): ALL 6 western g16-oracle variants swept, ZERO jl crashes
across ~46 stands. CI 6-bitexact/1-cornered/2-div; EM 4/4/2; TT 3/2/5 (TT div corroborates the #158 gap); IE 4/3/1;
UT 3/3/2; BM 5/2/0. No-regression on the mortality commits established by construction (per-tree dmr==0 /
empty-events guards ⇒ inert) + byte-identical primary stands + this cluster-wide sweep. (Divergences are the
pre-existing DGSCOR/density/#158 straddles, not the commits.)
Unit/integration suite ALSO green (complementary coverage): test_fixmort, test_mortmsb, test_allspecies,
test_canonical_multistand (multi-variant), test_dgstdev, test_multistand, test_multistand_sum all PASS ⇒ the
mistletoe + FIXMORT commits regress nothing at the unit level either.

## IE #143 AUTOES over-establishment — CONFIRMED REAL + SYSTEMATIC (2026-08-11), goal-doc "cornered" is WRONG for IE
Applied the BM/EM sign-tally lens to IE and it FAILED the straddle test: 12-stand IE FIA sign-tally = 9-HIGH/2-LOW/
1-BE (0 crashes), magnitudes up to **+108%** (mean strongly +). This is a SYSTEMATIC bias, NOT the "EM/IE ~7% tail
cornered" the goal doc claims (that verdict was verified for EM — which is genuinely 5H/5L balanced — but IE is
different). MEASURED root chain (FVSie_g16 + NOAUTOES):
- NOAUTOES collapses the 4 worst stands from +108%/+32%/+21%/+20% to +2.5-6.1% (the small growth straddle) ⇒ AUTOES
  is the cause, not growth/mortality.
- On stand 1143092701290487: jl AUTOES establishes 98 TPA @icyc2 then **584 TPA @icyc4** (exploding; TPA jumps
  426→924 at 2061); live establishes small, DECREASING increments (Σnewtpp 56→38 plot-trees across 50 plots).
- Live FIRES the ingrowth tally (ntally=99) at the SAME cycles (icyc2, icyc4) — the SCHEDULE matches. The bug is
  the tally MAGNITUDE: jl books the full MAXING-capped ITPP (~3/plot → 584 TPA); live books the INCREMENT
  NEWTPP=ITPP−NSTORE (~0.76/plot). jl re-books the standing sub-3" cohort every ingrowth cycle (compounding:
  icyc2's established trees are still sub-3" at icyc4 → re-tallied); live's NSTORE tracks them so they aren't.
- The `es_nstore`-reset-on-ntally==99 hypothesis was TESTED and REFUTED (persisting es_nstore changed nothing) ⇒
  the re-booking flows through `tpacre_ingro` (the current DBH<2.999 TPA sum that seeds the ingrowth tally), NOT
  es_nstore. The fix must EXCLUDE already-established/counted sub-3" trees from the ingrowth base — i.e., port
  live's NSTORE increment semantics into the tally MAGNITUDE (ie_autoes_tally/run), which is where d089b78's
  single-ingrowth fix stops short of the multi-ingrowth compounding.
STATUS: #143 is the genuine top open bug for the cluster (reconciliation's earlier "#143 fixed d089b78" is
incomplete — it fixed single-ingrowth stands but not the multi-ingrowth-cycle compounding). EM is NOT affected (dry
habitats zero the wet-side estb species; EM deterministic DG bit-exact + balanced tally). NEXT: instrument jl's
ie_autoes_tally inputs (tpacre_ingro, itpp, newtpp per plot) vs FVSie_g16's estab.f at icyc4 on 1143092701290487,
and port the ITPP−NSTORE increment (excluding the prior-established sub-3" cohort) into the ingrowth tally magnitude.

## IE #143 — ROOT CAUSE FOUND (2026-08-11 cont.): jl can't parse STRING plant-association habitat codes
Drilled the AUTOES over-establishment to its ultimate root (measured FVSie_g16, stand 1143092701290487, ihab path):
the AUTOES tally magnitude is driven by PROB1 (ESTOCK P(stocking)). At the 2nd (pure-ingrowth) tally jl PROB1=0.645
vs live 0.181 → jl over-books 584 TPA. PROB1=logistic(PN+esb_shift); at the inventory tally jl's WRONG PN (+0.556)
is masked by a compensating esb_shift (−2.035), but at the ingrowth tally esb_shift=0 EXPOSES it. PN comes from
`ie_estock`, whose per-series (IEQ) formulas are BIT-IDENTICAL jl↔live — so the error is the SERIES SELECTED:
- Live ESTOCK uses ihab=3 → IEQ=1 (Douglas-fir series, NO ba term ⇒ PN=−1.5085 CONSTANT).
- jl uses ihab=13 → IEQ=4 (subalpine series) ⇒ PN=+0.556.
ihab comes from the habitat code. This stand's `PV_CODE = "CDS715"` (a STRING plant-association code); `PV_REF_CODE
= "627.0"`. jl's FIA reader (fia_database.jl:184) parses PV_CODE as a Float32 → "CDS715" → 0 → FALLS BACK to the
numeric PV_REF_CODE 627 → ie_estab_indices(627) → ihab=13. LIVE instead runs `ie/habtyp.f` HABTYP, a STRING
crosswalk: "CDS715" → PCOML index 9 → MAPR6[9]=15 → JTYPE/KTYPE → ITYPE=4 → **KODTYP=260** (MTYPE[4]); 260 ≤ 269 ⇒
ie_estab_indices(260) → ihab=3. (Live ignores PV_REF_CODE when PV_CODE is present.)
FIX (bounded but substantial): port ie/habtyp.f's plant-association string crosswalk (DATA tables PCOML[40 strings],
MAPR6[40], JTYPE[95], KTYPE[95], MTYPE[30] + HBDECD string-match) into fia_database.jl so a string PV_CODE resolves
to the FVS habitat code (CDS715→260) instead of falling back to PV_REF_CODE. LIKELY SYSTEMATIC across IE (and EM/UT/TT
— same numeric-only PV_CODE reader) FIA stands carrying alphanumeric plant-association codes ⇒ explains the 9H:2L IE
tally. VALIDATION PLAN: after the port, ie_estab_indices(260)→ihab=3→IEQ=1 (DF, ba-independent) ⇒ jl PN=−1.5085
(identical formula) ⇒ PROB1 matches live ⇒ AUTOES tally collapses to live's; re-run the 12-stand IE sign-tally
(expect it to balance like EM's) + confirm iet01 (numeric-code path) byte-unchanged.

## IE #143 — ★★★ FIXED 2026-08-11 (9cc7d7a): ported ie/habtyp.f string habitat crosswalk
The habitat-string-crosswalk root cause (above) is RESOLVED. Added `ie_pa_habitat_code` (PCOML[40]/MAPR6[40]/
KTYPE[95]/MTYPE[30]; final=MTYPE[KTYPE[MAPR6[i]]], JTYPE strictly ascending) and wired it into the FIA reader for IE
string PV_CODEs. All 40 PCOML→code mappings validated BIT-IDENTICAL vs a live FVSie_g16 habtyp dump. Results:
- Worst stand 1143092701290487: 2071 TPA +107.9% → +4.1%.
- 12-stand IE FIA sign-tally: 9H/2L/1BE (systematic, up to +108%) → 6H/4L/2BE (STRADDLES, mostly ±4%), 0 crashes.
  The systematic over-establishment bias is eliminated; the residual is the cornered growth/AUTOES RNG straddle
  (same class as EM). ⇒ IE is now bit-exact-or-cornered like EM.
- Gated to InlandEmpire + string PV_CODE ⇒ numeric-code IE stands, non-DB stands (iet01 unchanged), and all other
  variants byte-identical by construction.
FOLLOW-UP: EM/UT/TT share the numeric-only PV_CODE reader and have their OWN per-variant habtyp string tables — port
each if their FIA stands carry alphanumeric plant-association codes (run a per-variant sign-tally to check).

## String-habitat bug — cluster-wide sweep COMPLETE (2026-08-11): IE-only, no follow-up needed
Checked all western variants' FIA PV_CODE format + reader handling (the follow-up flagged after the IE #143 fix):
- EM/UT/TT/CI: FIA stands carry NUMERIC PV_CODEs (e.g. EM "380", UT "204091", TT "46111", CI "41732") ⇒ the
  numeric reader is correct; no string crosswalk needed.
- KT: no PV_CODE column.
- BM: DOES carry STRING plant-association codes (CWG111/CDS625/CDG112…) — BUT its FIA reader ALREADY does the
  BM_PCOML string→KODTYP lookup (fia_database.jl:143-147), added previously. VERIFIED vs FVSbm_g16: CWG111 →
  BM_PCOML idx 77 → bm_habtyp→"CWG111" → ecocls SDIDEF PP=263/DF=376/GF=700 == live (I'd wrongly hypothesized jl
  defaulted to CWG113; it does not). End-to-end stand 1127530489290487: jl vs live BA 244/244, SDI 479/479
  bit-exact.
⇒ The string-habitat-code parse bug was IE-ONLY (now fixed, 9cc7d7a). No further habitat-crosswalk ports needed.

## NEW FINDING 2026-08-11: UT systematic self-thin UNDER-KILL (7:0 tally) — RESOLVED to CORNERED (see the later "UT under-kill — RESOLVED to CORNERED" section; deterministic PJ regent proven bit-exact, residual = ZZRAN realization amplified by asymmetric self-thin feedback). The "NOT yet cracked" below is the mid-investigation state, superseded.
Completed the cluster-wide FIA sign-tally (UT/TT/KT/CR, the ones not tallied earlier this session): CR 2H/0L/8BE and
TT 4H/1L/5BE (mild, mostly cornered), KT no FIA stands — but **UT = 7-HIGH/0-LOW/3-BE, up to +43% TPA** — a
SYSTEMATIC bias (the BM #140 signature), overturning the goal-doc "UT complete". jl UNDER-thins (worst stand
11937105010690: 2056 TPA 797/1137 live/jl, QMD 7.0/5.8, BA close 216/211 — same pattern as BM #140).
LOCALIZED (FVSut_g16 morts.f instrumentation): the Reineke self-thin INPUTS match — sumdr0 jl 8523 vs live 8527
(same 34 trees, dmin 0.10 dmax 17.0), and UT's DBHZEIDE=0 so no min-DBH filter gap. The divergence is in the
self-thin SOLVE/kill dynamics: live's ut/morts.f runs an IPASS QMD-CONVERGENCE loop (morts.f:561-604 — recompute
the post-kill Reineke DR10N from survivors, and if |D10−D10N|>0.1 set D10=D10N and GO TO 10 to re-solve tn10),
which jl's utah/mortality.jl LACKS (one pass). On the worst stand live iterates tn10 1593.7→1504.3 (more kill).
BUT — adding the IPASS loop to jl did NOT fix it: jl's uniform-RN kill leaves the post-kill QMD ~unchanged ⇒ the
loop converges in one pass, while live's shifts. AND at cyc1 with IDENTICAL inputs (tt=2259, sdimax=460.7) jl
tn10=2068 vs live 2087 (jl slightly OVER-kills), so the sign FLIPS across cycles ⇒ the +43% is a QMD-feedback
amplification of a subtle self-thin-solve/kill-distribution difference, not cleanly the IPASS. ROOT NOT YET cracked
— reverted the IPASS attempt (unvalidated). NEXT: at MATCHED (tt, DR10, sdimax) compare jl's _em_tn10_iter output +
the per-tree WKI distribution vs live's morts.f (does live's uniform kill really leave QMD invariant, or is there a
size-dependent X/XMORT or a DR10-vs-DQ10 mismatch feeding the IPASS?). Distinct from the accepted UT #156 dense-PJ
regent-ZZRAN straddle (that's growth; this is mortality). Oracle FVSut_g16 durable at /workspace/.utwork.

## UT under-kill — REFINED 2026-08-11 (root is the entry DR10 seed, NOT the line/IPASS)
Followed up the UT self-thin under-kill with FVSut_g16 (dumped SLPMRT/CEPMRT/D10/T85D10/TN10) — and RULED OUT the
two natural hypotheses:
- Self-thin LINE calibration MATCHES: live SLPMRT=−0.6040, CEPMRT=8.2902 == jl's exactly. Implemented line
  PERSISTENCE (jl was re-solving each cycle) — confirmed jl's slp stays −0.6040 across cycles — but it did NOT change
  the result (the persisted line is already right). REVERTED (touched validated EM, no benefit).
- IPASS QMD-convergence: live's tn10 caps at T85D10 each pass, and its IPASS lowers T85D10 by raising the post-kill
  D10 (4.239→4.394 within cyc3). Implemented the IPASS loop — but jl's D10 doesn't rise (converges pass 1), so
  inert. REVERTED.
REAL ROOT (measured, identical cyc1 trees): jl's ENTRY Reineke DR10 = 2.961 vs live 2.918 (+1.5%). Since DR0/sumdr0
MATCH (8523/8527, no growth), the difference is in the GROWTH term of the Reineke sum: jl `g = diam_growth/bark`
(sumdr10 = Σp·(D+g)^1.605) vs live `G = (DG(I)/BRATIO)·(FINT/10)`. jl's DR10 slightly HIGH ⇒ TMD10 low ⇒ T85D10
slightly low at entry but the self-thin FEEDBACK (higher D10 → the line evaluates higher tn10 as the stand thins)
compounds it to the +43% under-kill. NEXT: dump jl per-tree (d, diam_growth, bark, g) vs FVSut_g16 morts.f
(D, DG, BRATIO, G) at cyc1 — likely a DG-period/bark subtlety in the mortality g (the mortality sum uses the shared
bark_ratio; verify it equals live's BRATIO and that diam_growth is the FINT-period DG, not 5-yr). Distinct from the
growth-DG which is bit-exact (the mortality g-reconstruction from diam_growth is the suspect, not the DG itself).

## UT under-kill — TRUE ROOT FOUND 2026-08-11: species-specific DG (PI/MC), NOT the mortality
The mortality investigation (line/IPASS/DR10) was chasing a SYMPTOM. Per-tree comparison of the mortality Reineke
sum at cyc1 (jl vs FVSut_g16, identical trees) shows the mortality g-reconstruction (g = diam_growth/bark) is
computed IDENTICALLY (bark matches per-tree) — the divergence is in the DIAMETER-GROWTH `diam_growth` INPUT for
specific UT woodland species:
- sp11 **PI (pinyon)**, d=3.9: jl dg=1.095 vs live 1.401 (−22%).
- sp20 **MC (curl-leaf mountain-mahogany)**, d=10: jl dg=0.0019 vs live 0.0139 (~7× low).
- sp16 **UJ (Utah juniper)** (d 12-17): dg MATCHES bit-close (0.481/0.481 etc.) — so NOT all woodland species,
  just PI + MC (and likely other surrogate-equation species).
UT's dgf.f notes "SPECIES USING SURROGATE EQUATIONS FROM THE CR VARIANT HAVE SPECIAL..." handling — PI/MC use CR
surrogate DGF; jl's port of those surrogates is wrong (under-grows). The wrong DG ⇒ wrong Reineke DR10 (jl 2.961 vs
live 2.918) ⇒ wrong self-thin target ⇒ the +43% under-kill (7:0 tally), amplified by QMD-feedback. ⇒ this is a
LARGE-TREE DGF bug for UT PI(11)/MC(20), distinct from the #156 dense-PJ regent (which was proven bit-exact) — those
were the small-tree PJ regent; these are the LARGE-tree (d>3") woodland DGF. NEXT: instrument jl's UT dgf!/DDS for
PI(11) at d=3.9 vs FVSut_g16 dgf.f (DGCON/DGFOR/model-type/surrogate coeffs) — a bounded per-species coefficient
trace. This REDIRECTS the UT fix from mortality to growth.

## UT under-kill — FINAL root pinned 2026-08-11: the PJ/woodland small-tree REGENT (not large-tree DGF)
Corrected the "large-tree DGF" lead: UT_RG_XMAX[11=PI]=99, XMIN[11]=90 — so ALL pinyon (PI) trees below 90" DBH (i.e.
every PI in a real stand) use the PURE small-tree REGENT path (xwt=0), NOT the large-tree DGF. Verified the large-tree
DGF DDS for PI d=3.9 is BIT-EXACT vs FVSut_g16 (both 0.72555, df/bark/conspp identical) — but it's IRRELEVANT because
PI never reaches the large-tree threshold. Same for MC(20) (XMAX=99). So the DG divergence (jl 1.095 vs live 1.401 for
PI; 0.0019 vs 0.0139 for MC) comes from `utah/regent.jl` small_tree_growth!'s PJ/woodland branch (regent.jl:68
"PJ/GB/MC non-aspen non-conifer" → the ((SJ/5)(SJ·1.5−H)/(SJ·1.5))·0.83 POTHTG form → DK/DKK inline H-D → DDS →
XWT-blend), which regent.jl:11 EXPLICITLY flags "ported faithfully (need pure stands to validate)" — i.e. never
validated. ⇒ UT self-thin under-kill (7:0, +43%) ROOT = the under-validated UT PJ/woodland small-tree REGENT DG for
PI/MC (and likely WJ/PM/RM/GB). The mortality (line/IPASS/DR10) and large-tree DGF were all confirmed correct/symptom.
NEXT: instrument jl small_tree_growth! for PI(11) d=3.9 (SJ/SITEAR, POTHTG, VIGOR, DK/DKK, DDS, xwt) vs FVSut_g16
regent.f — a bounded per-branch trace of the PJ regent, analogous to the TT #158 regent work. This is the last open
western-cluster growth bug; distinct from #156 (that reproducer's PJ regent was proven bit-exact — a different stand
/species mix; this shows PI/MC specifically diverge).

## UT under-kill — RESOLVED to CORNERED 2026-08-11: deterministic PJ regent bit-exact; residual = ZZRAN realization
Traced the PI(11) PJ-regent height growth per-component vs FVSut_g16 and the ENTIRE DETERMINISTIC path is BIT-EXACT:
- POTHTG=0.4980, PCTRED=0.8066, VIGOR=1.0000, CON=1.0, HTGR(deterministic)=0.4017 — jl == live to all printed digits.
- DK=(HK−4.5)·10/(SITEAR−4.5), DKK, DGK=(DK−DKK)·bark — jl formulas + values match live (dk 7.334, dkk 6.0).
- The large-tree DGF DDS (0.72555) is bit-exact too (though irrelevant: PI XMAX=99 ⇒ always regent).
The ONLY non-matching quantity is the per-tree ZZRAN draw: jl bachlo=−0.683 vs live's ≈+0.248 (htg 0.3334 vs 0.4265).
The ZZRAN reject range is IDENTICAL ([−2, 0.5], regent.f:343), so the distribution matches — it's the never-FFI RNG
BYTE-STREAM realization (jl's bachlo ≠ live's ZZRAN sequence; the accepted #142/#156 class). The 25-stand 16-HIGH/0-LOW
tally is that symmetric RNG straddle AMPLIFIED by the ASYMMETRIC self-thin feedback: an under-grow (RNG-low) cycle
compounds (fewer kills → more survivors → lower QMD → higher self-thin target → fewer kills → runs away to +43%),
while an over-grow cycle self-corrects at the self-thin line (bounded). So symmetric per-tree RNG ⇒ systematic
jl-HIGH stand-level tally — WITHOUT any deterministic bug. ⇒ UT MEETS bit-exact-or-cornered (deterministic PJ regent
bit-exact; residual = never-FFI ZZRAN realization). This CORRECTS the earlier "real systematic under-kill bug" reading
— the sign-tally skew was real but its ROOT is the accepted RNG realization, not a deterministic port error. The
mortality (IPASS/line/DR10) and the DG (htgr/DK/DKK/DGF) are all confirmed faithful. NOTE: a future bit-exact-RNG
effort (matching bachlo's byte-stream to FVS ZZRAN through the regent) would tighten UT the most (its PJ-heavy dense
stands amplify the RNG residual hardest), but that is the cross-variant never-FFI-RNG undertaking, not a UT bug.
⇒ Western cluster growth/mortality: all variants bit-exact-or-cornered.

## CI volume (MATW/FW2W) — CLOSED bit-exact-or-cornered 2026-08-11: equations faithful; multi-cycle = #142 straddle
The goal doc's last unmeasured CI item ("volume MATW/FW2W unmeasured"). MEASURED on 3 real CI FIA stands
(11790085/11791809/11792275 010690), full .sum vol columns jl vs FVSci_clean:
- **Cycle 0 (identical input trees): volume is BIT-EXACT on all 3 stands** — TCuFt 964/964, 2944/2944, 1782/1782;
  MCuFt and BdFt likewise identical. The MATW (merch cubic) / FW2W (board) equations reproduce live exactly.
- Divergence appears ONLY in later cycles and TRACKS the growth/mortality straddle: by 2056 jl has HIGHER TPA
  (1579 vs 1524) yet LOWER volume (TCuFt −1.2%, MCuFt −5.6%, BdFt −6.3%) — jl retains more, individually-smaller
  trees (less large-tree DG). Mixed-sign across stands (stand 11792275 MCuFt +6.1%). That is the signature of the
  accepted CI #142 DGSCOR/mortality RNG realization propagating into the volume totals, NOT a volume-equation bug
  (a wrong equation would bias one sign on identical input at cyc0 — it doesn't).
⇒ CI volume MATW/FW2W MEETS bit-exact-or-cornered: equations bit-exact on identical trees; multi-cycle residual is
the already-cornered #142 straddle carried through volume. Closes the last flagged CI-volume item. (SMHTGF small-tree
stochastic is the same never-FFI-RNG class.)

## BC volume (merch/board "TODO") — CHARACTERIZED by measurement 2026-08-11: not a volume-equation gap
Measured all_BC (real BC multi-species stand) jl vs FVSbc_clean, full .sum + cross-checked the .out cruise:
- **Total cubic**: cyc0 (1990) BIT-EXACT (186 m³/ha == live 186, == cruise CUBIC METERS/HA 186.02). Multi-cycle
  +4.9% by 2090 (jl 1075 vs live 1025; jl BA also +8.8%, TPA 1253 vs 1292) = the BC growth-realization tail
  (fewer, bigger trees), the accepted cornered class for the completed BC growth+yield port.
- **Board feet**: live .out cruise = `BOARD FEET/HA 0.00` and the per-species board table is all "---"; live .sum
  board column = 0 every cycle. The metric BC variant emits NO board feet for this stand ⇒ the BFVOL "TODO" has
  NO nonzero validation target (jl's effective 0 already matches). Not a real gap on the canonical stand.
- **Merch cubic**: jl DOES compute BC merch (the volume.jl "merch=0" comment is STALE — true only at cyc0 when all
  trees are <17.5cm DBHMIN; bc_tree_vol's TVOL merch-log path fires as trees exceed it, giving jl merch 0→863 by
  2090). But live's .sum standing-merch column = 0 EVERY cycle — even at 1990 where the .out cruise computes live
  merch = 156.15 m³/ha. ⇒ live CAN compute merch; its .sum merch COLUMN is 0 by reporting spec. This is the SAME
  class as CI's .sum-merch-column finding (fvsjl-ci-variant-port 7th note): FVS's .sum "MERCH CU FT" is written by
  disply.f from a summary O-array distinct from per-tree/cruise merch. Consequence: jl's BC .sum emits a merch
  column (863) live leaves 0 ⇒ jl's BC .sum has 2 extra volume fields (28 vs 26) — a .sum column-STRUCTURE
  difference, output-only (volume never feeds growth/mortality), not a simulation divergence.
VERDICT: BC total-cubic MEETS bit-exact-or-cornered (cyc0 bit-exact + growth tail). The "merch/board vol TODO" is
NOT a volume-equation port — board has no metric target; merch equations compute correctly and the residual is the
cross-variant .sum summary-merch-column reporting spec (shared with CI, disply.f O-array), which is output-only and
lower priority. No speculative BFVOL board port opened (unvalidatable) and jl's .sum merch NOT forced to 0 (would
hide correct volume; matching live's 0 requires reverse-engineering FVS's summary-merch spec = a separate
cross-variant reporting effort, not a BC growth/volume gap).

★ 2026-08-11 ROOT RESOLVED TO FVS SOURCE (supersedes "reporting spec" above with the exact FORMAT): BC uses
metric/vbase/sumout.f (confirmed: `diff bin/FVSbc_buildDir/sumout.f metric/vbase/sumout.f` == identical). The
METRIC .sum row FORMAT 20 = `2I4,I6,I4,I5,2I4,F5.1,7I6,I4,I5,2I4,F5.1,2X,I6,I5,I6,2X,F6.1,1X,I3,1X,2I1` — **7I6**
volume integers (IOSUM 4=TotalCuM, 5=MerchCuM, 6=MerchCuM-net-of-cull, 7=RemTrees, 8=RemTotal, 9=RemMerch,
10=RemMerch-net); the imperial FORMAT 20 has **9I6** (adds the SAWLOG/board columns, dropped as "N/A" in the metric
BC header FORMAT 13). jl's src/io/summary.jl:16-19 `_SUM_ROW_FMT` HARDCODES the imperial 9-integer layout
(`%6d`×9 = cuft/mcuft/scuft/bdft + 5 removed) ⇒ applies it to BC too ⇒ the 2 extra fields (28 vs 26) + misplaced
values. ⇒ TWO distinct BC .sum issues, BOTH output-only (volume never feeds growth/mortality):
 (A) COLUMN LAYOUT — jl needs a metric 7I6 variant of _SUM_ROW_FMT keyed on metric variants (BC/ON), dropping the
     sawlog/board ints. Bounded fix (add a second Printf.Format + a metric branch in write_sum_row); MUST leave the
     imperial path byte-identical (SN/CI/NE/CS/LS .sum tests pass bit-exact) — validate no-regression cluster-wide.
 (B) MERCH VALUE — even with layout (A) fixed, jl's IOSUM(5) merch (863 by 2090) ≠ live's 0. Live's .sum merch is 0
     every cycle even at 1990 where the .out cruise computes merch=156.15 ⇒ live's SUMMARY-merch O-array (grstat/
     disply.f accumulation feeding IOSUM(5)) is 0, distinct from per-tree/cruise merch — the SAME summary-merch spec
     as CI (833 treelist vs 758 .sum). Resolving B needs the FVS summary-merch-accumulation trace (why IOSUM(5)=0
     for BC) — deep, output-only, shared with CI.
⇒ ACTIONABLE next (focused fresh pass, output-only lowest priority): implement the metric 7I6 .sum layout (A) with
cluster-wide no-regression validation, then trace the summary-merch O-array (B). NOT bolted on at this depth: a
layout-only fix wouldn't bit-match (B still differs) and touching the shared _SUM_ROW_FMT needs careful before/after
.sum validation across all 5 passing imperial variants. BC GROWTH remains validated (e2e harness cols 3-8); this is
purely the .sum volume-column reporting fidelity.

## BC .sum FULL column enumeration + (B) root RESOLVED to FVS source 2026-08-11 (doctrine #2 satisfied)
Full-row jl-vs-FVSbc_clean comparison on all_BC (11 cycles), every column classified:
- ★ (B) MERCH ROOT = FVS SOURCE, not inference: the standard imperial vols.f loads `MCFV(I)=MCF` (CI/SN vols.f:199)
  but BC's vols.f NEVER assigns MCFV — it computes merch into `WK1(I)=VM` (vols.f:167, defect-corrected :221) and
  routes it ONLY to the per-species accumulator SPCMC + ECON, never to the standing MCFV array (grep: "no MCFV=
  assignment in BC vols.f"). So OMCCUR(7)=PCTILE(MCFV)=0 (gradd.f:319) ⇒ live BC .sum merch (IOSUM 5, disply.f:427)
  is STRUCTURALLY 0 every cycle. Also vols.f:236-240 `BFV(I)=0.0` "Board feet not computed in this variant" ⇒ board
  structurally 0. ⇒ FAITHFUL jl behavior = BC .sum merch AND board = 0 (mirror FVS's unpopulated MCFV/BFV), NOT jl's
  computed 863 (jl's per-tree merch_cuft_vol stays correct for CSV/other outputs). Resolves (B) WITHOUT inference.
- (A) LAYOUT: metric 7I6 (drops sawlog scuft + rem_scuft vs imperial 9I6), as established.
- TOPHT: jl 25/21/24 vs live 50/47/47 from 2010 = the NOISY-HEIGHT AVHT40 cornered class. all_BC.tre carries noisy/
  partially-missing heights (DBH 12.7cm→HT 6.7m vs DBH 7.9→7.5; many blank → HTCALC-filled) ⇒ top-height is
  noise-dominated and swings with the mortality realization (which tall trees survive). The ROBUST DIAMETER metrics
  BA/SDI/QMD/total-cuft all TRACK within a few % every cycle (2090: BA 62/57, SDI 1235/1163, QMD 25.2/23.8, cuft
  1075/1025) ⇒ the BC growth+volume EQUATIONS are faithful; only the noise-dominated TopHt swings (same class as the
  all_BC_essf "garbage heights" cornered note the e2e harness already accepts).
- FORTYPE: jl 999 vs live 201 = BC forest-type classification not ported (jl defaults 999); cosmetic .sum column.
- TPA/SDI/TopHt ±1-2: the known #129 plot-area TPA rounding (2089 vs 2087) + integer-round straddle.
VERDICT: BC growth+volume EQUATIONS remain bit-exact-or-cornered (diameter metrics track; merch/total volume
equations faithful). ALL BC .sum residuals are OUTPUT-ONLY (metric 7I6 layout + structural merch/board=0 + BC
forest-type classification) or CORNERED realization (noisy-height TopHt, #129 rounding) — NONE affects the
simulation. The output-only .sum-fidelity fix (metric layout + merch/board=0 + BC fortype) is a bounded focused pass,
lowest priority; the growth/volume port is validated. This CLOSES the "BC merch/board vol" investigation: the
equations are faithful, the .sum merch=0 is FVS-structural (vols.f), and board has no target (not computed in BC).
