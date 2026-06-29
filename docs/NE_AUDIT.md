# FVS Northeast (NE) port — audit ledger

Analogue of the SN audit campaign (`docs/audit/INDEX.md`). Doctrine: TRACE LOGIC NOT RUNTIME;
UPSTREAM-FIRST; REGRESSION = MASKED-BUG SIGNAL; PORT FAITHFULLY + validate vs LIVE FVSne THEN test;
DOCUMENT EVERY verdict; VARIANT-AWARE (keep SN bit-exact). Every flag below records both-sides logic
and a verdict (OPEN / FAITHFUL / FALSE-POSITIVE / FIXED).

Oracle: `bash test/harness/ne_oracle.sh <key> <outdir>` (relinks live FVSne). net01.sum.save is stale
for cyc1+ — validate against the live binary.

## Methodology
Sweep the NE port subsystem-by-subsystem; for each, run a scenario that EXERCISES it, diff jl vs live,
and root-cause every divergence to either a real bug (FIX) or a documented faithful/ULP class. Do NOT
accept an aggregate match (BA/QMD within a few %) as proof — distribution-level errors hide under matching
aggregates (see A1). Flag with debug dumps from the live Fortran, not inferred from test pass/fail.

---

## A1 — Stand-2 thin divergence — SYSTEMIC BUG FIXED; residual pinned to OLDRN ★ flagged by user "−40 vs −22 BA not acceptable"

> **RESOLUTION (read the "A1 — FIX LANDED" section below first).** The root cause was NOT an unalignable
> RNG realization (the earlier conclusion in this section, kept for the investigation record). It was a
> CONCRETE coefficient bug: 20 species had the IFOR=3 Allegheny HT-DBH override values instead of the
> defaults. Fixing it closed the first/largest draw divergence (cyc2 +27→0). The "RNG not alignable" framing
> below was WRONG — the draw sequence WAS alignable once the upstream coefficient was corrected; this is why
> the draw-counter method (find first divergence → trace the term → fix) is the right approach, not accepting
> the divergence. Smaller cyc4+ divergences remain (same method).

**Symptom.** net01 stand-2 (repeated THINDBH every 3 cycles) `.sum`: through ~2110 jl tracks live within
~Δ4 BA, then at the **2130 thin (cycle 15)** jl removes ~2× the wood — after-cut BA **83 (jl) vs 99 (live)**,
removed TPA 20 vs 14, removed BdFt 6658 vs 3238. ~14% BA divergence that propagates to 2140+.

**Not the thin logic.** The cut-selection (cuts.f) is variant-agnostic and validated faithful in isolation
(THINBBA/THINABA/THINSDI on net01 stand-1 all match live within cyc-1 drift — test_net01.jl). At 2130 both
sides parse the THINDBH specs identically (live debug-dump of VALMIN/VALMAX/CTPA/CBA confirmed) and apply
the same per-class TPA targets.

**Root cause = the input DBH DISTRIBUTION, not the mean.** Per-class stocking at 2130 (live debug dump of
CLSSTK CSTOCK vs a jl pre-thin histogram):

| DBH class | LIVE CSTOCK | jl cstock | THINDBH target | consequence |
|---|---|---|---|---|
| 4–8″  | 75.5 | 94.0 | 125 | both < target, no cut |
| 8–12″ | 66.0 | 44.6 | 60  | live cuts 6, jl cuts 0 |
| 12–16″| 35.7 | 17.5 | 35  | ~at target |
| 16–20″| **9.5** | **27.0** | 15 | **live no-cut; jl cuts 12** |

jl's distribution is too **spread** — fat tails (excess 4–8″ AND 16–20″), hollow middle (deficit 8–12″,
12–16″). Aggregates hide it: pre-thin BA 123 (jl) vs 121 (live), QMD 11.3 vs 11.1 — ~2%.

**Evidence chain (each step traced to live FVSne / the NE Fortran, not inferred from tests):**

1. *Thin parse + targets identical.* Live debug dump of cuts.f VALMIN/VALMAX/CTPA/CBA per class at 2130
   matches jl's; both apply the same per-DBH-class TPA targets. The cut-selection is variant-agnostic and is
   independently faithful (THINBBA/ABA/SDI/HT/CC + FIXMORT, A3). So it's not the thin.
2. *It's the input DBH distribution, not the mean.* Per-class CSTOCK at 2130 — jl is mean-preservingly
   OVER-SPREAD (16–20″ jl 27.0 vs live 9.5; hollow 8–12″/12–16″), while pre-thin BA 123 vs 121 (~2%).
3. *It's the record TRIPLING.* NOTRIPLE-on-both at 2130: jl 16–20″ = 8.9, live 10.8 — AGREE (both < target,
   no cut). With TRIPLE, live stays 9.5 but jl balloons to 27.
4. *Not the tripling spread or coefficients.* FU/FM/FL (= ne/dgdriv.f:626-628), cadence ICL4=2
   (= ne/grinit.f:183), dgscor! bound (jl = dgscor.f `|frm|>DGSD·SIGMA` redraw + dds>4 taper), DGSD=2
   (= NE grinit), and comcup! pruning (every cycle, PROB≤1e-5 = base/comcup.f:50) are ALL faithful. Per-class
   16–20″ MATCHES live at cyc3 (2020) right after tripling — the initial spread is correct.
5. *It's the stochastic realization, not a systematic bias.* Per-class DBH histograms are over-spread in BOTH
   stands but in OPPOSITE directions:

   | @2090 unthinned stand-1 | 8–12 | 12–16 | 16–20 | 20+ | BA agg |
   |---|---|---|---|---|---|
   | LIVE | 17.6 | 6.1 | **63.5** | **23.8** | 194 |
   | jl   | 21.4 | 7.9 | **49.4** | **31.0** | 192 |

   Unthinned stand-1: jl's 16–20″ is LOW (49 vs 63), 20+″ HIGH (31 vs 24). Thinned stand-2: jl's 16–20″ is
   HIGH (27 vs 9.5). A systematic growth/mortality bias would push the SAME way in both; opposite-direction
   per-class diffs with faithful aggregates ⇒ a different-but-valid stochastic REALIZATION.

**ROOT CAUSE (final).** The NE BACHLO/dgscor stochastic RNG stream is **not bit-aligned to live FVSne**. The
model is faithful — aggregates, NOTRIPLE, every coefficient, cyc1-2 growth, and A2/A3/A4 all match — but the
per-tree stochastic draws are a different valid realization (the 9× tripled records draw a 9× BACHLO stream
that drifts from FVSne over the cycles). SN bit-aligned its RNG against the Oracle-A 1:1 transliteration; NE
has NO such reference, so its multi-cycle stochastic stream was never aligned. Repeated class-target thinning
(stand-2) AMPLIFIES the realization difference at the THINDBH boundaries (a mean-preserving spread → a discrete
over/under-cut), which is why the `.sum` divergence surfaces there and not in aggregates. This is the SAME
class as the accepted SN COMPRESS eigensolver divergence — but the user has (rightly) asked it be closed.

**Seed is aligned — it's the DRAW SEQUENCE.** jl's NE main-stream seed = 55329 (= `ne/blkdat.f` DATA
S0/SS; net01 has no RANNSEED override), identical to FVSne. So the two LCG streams START identically; they
diverge because jl makes a different NUMBER/ORDER of `rann!`/`bachlo` draws than FVSne somewhere across setup
+ cyc0-2 (deterministic tripling) before the first stochastic cycle (cyc3, 2020) — a single extra/missing
draw offsets the whole downstream stream. Cyc1-2 aggregates match (deterministic, no RNG), so the offset is
either in cyc0 setup, the mortality (VARMRT draws), or the cyc3 dgscor draw order.

**FIRST DIVERGENCE LOCALIZED (draw-counter instrumentation, both sides).** Added a RANN draw counter to jl
`rann!` and to FVSne `rann.f` (COMMON /RNDBG/, dumped per cycle from grincr.f). Cumulative main-stream draws,
unthinned stand-1 (jl cyc N ↔ FVS ICYC N+1):

| cycle | jl draws | FVS draws | Δ |
|---|---|---|---|
| 1 (1990→2000) | 87 | 87 | **0 — exact** |
| 2 (2000→2010) | 321 | 294 | **+27 (jl)** ← first divergence |
| 3 (2010→2020) | 1131 | 1131 | 0 (then cascades) |

Cycle 1 matches the draw count EXACTLY; the first divergence is entirely within the **2nd growth cycle**, and
phase-boundary counters show ALL of it is in **`small_tree_growth!` (REGENT)** — DG/HTG/mortality/tripling draw
ZERO. jl draws +27 rann = **+9 `bachlo`** there. Both jl (small_tree_growth.jl:76, `for l in 0:nrec-1`) and
FVS (regent.f:263, the per-triple loop) both draw one bachlo per tripled sub-record. CONFIRMED: the +9 bachlo
is in REGENT, cycle 2 (the 2nd tripling cycle), unthinned stand (thinning-independent). LEADING HYPOTHESIS
(not yet confirmed at the record level): jl processes ~3 extra small-tree records because their cycle-2-entry
DBH lands on a different side of the 5″ small-tree boundary (`d < NE_REGENT_XMAX`, small_tree_growth.jl:54)
than in FVS — seeded by the cycle-1 per-tree DBH not being BIT-exact (the ~0.5-1% A2 growth residual: cyc-1
deterministic tripling growth is faithful in aggregate but not per-tree). A few records near 5″ flip
membership ⇒ different REGENT draw count ⇒ rejection-sampling cascade. If true, A1 ⊂ "make NE per-tree growth
bit-exact" — the boundary flip is a symptom; the upstream cause is per-tree DG/HTG precision. NEEDS the
per-record dump (below) to confirm the mechanism before any fix.

**CONFIRMED (per-record dump, both sides) — A1's upstream root is a SPECIES-SPECIFIC REGENT small-tree DG
discrepancy.** Dumped every record REGENT draws for at the year-2000 entry (jl cyc1 / FVS regent.f:158 ICYC2;
both skip `D ≥ XMX=5″`, same threshold). The year-2000 small-tree DBHs differ BY SPECIES:

| species | jl DBHs | FVS DBHs | verdict |
|---|---|---|---|
| 9 (WP, conifer) | 1.354, 1.297, 1.323 | 1.357, 1.300, 1.326 | **MATCH (~0.3%)** |
| 27 (SM) | 2.100, 2.103, 2.144 | 2.202, 2.205, 2.253 | **jl ~4.7% LOW** |
| 27 (SM) near 5″ | 4.94, 4.94, 4.97 (3 recs <5) | 4.90, 4.99 (2 recs <5) | **boundary FLIP** |
| 30 (YB) | similar ~4-5% low | — | jl low |

So white pine (sp 9) matches bit-close but sugar maple (27) / yellow birch (30) are ~4-5% LOW in jl. This is a
real per-species REGENT (small-tree) growth error in CYCLE 0 (1990→2000) — NOT RNG phase, NOT ULP. It shifts
the year-2000 hardwood small-tree DBHs down, so ~1 hardwood parent's tripled set keeps 3 records < 5″ in jl
where FVS has only 2 (the upper crossed 5″) → jl draws for the extra record → the +9 bachlo / +27 rann → the
RNG stream cascade → the per-tree distribution divergence → the THINDBH over-cut. The whole A1 chain reduces to
this. (Aggregates hid it: small hardwoods are a small BA fraction; A2/A4 matched within ~1%.)

**Next (task #50).** Find the per-species REGENT DG error for sp 27/30 (hardwoods). Candidates in
small_tree_growth.jl (NE): the Wykoff HT→DBH inverse `_htdbh_dbh` coefficients, the `regent_min_diam` (FVS
DIAM budwidth floor, regent.f:103), the NC-128 height curve (`ne_htcalc_*`), or BALMOD `b3_dg`. Instrument one
sp-27 small tree's cycle-0 DG step-by-step (htgr → hk → dkk/dk → dgsm → dg) vs FVS regent.f and find the
diverging term. Fix it (NE small-tree path only — keep SN `<3″` bit-exact), then re-run the draw counter to
confirm cyc-2 Δ→0 and the stand-2 .sum 2130 thin converges. This is the head of the A1 chain.

**Status: OPEN** (model faithful, stochastic stream not yet aligned). Previously mis-closed as "faithful
within drift" — a lax verdict the user correctly rejected; re-opened.

---

## A1 — FIX LANDED (partial): IFOR=3 Allegheny HT-DBH overrides were applied unconditionally

ROOT CAUSE FOUND + FIXED. `data/northeast/htdbh_coeffs.csv` carried the FVS **IFOR=3 (Allegheny) HT-DBH
override** values (ne/sitset.f:428-490, `IF(IFOR.EQ.3)`) for **20 species** (26,27,30,31,33,40,41,42,44,54,
55,60,64,67,69,71,93,102,106,108) instead of the DEFAULT values (sitset.f:207-423). net01 is IFOR=2, so FVS
uses the defaults; jl used the Allegheny values → the small-tree (REGENT) HT→DBH inverse `HT2/(ln(H-4.5)-HT1)-1`
was wrong for these hardwoods (e.g. sp27 DK 0.7068 vs FVS 0.7396) → ~4-5% low small-tree DBH → boundary flip
→ +9 REGENT bachlo at cyc2 → RNG cascade → the stand-2 2130 thin over-cut. CONFIRMED by step-trace: HTGR
(height) matched; DK (HT→DBH) was low; back-solved FVS AX=4.5832 = default HT1=4.4834-path, not jl's 4.6354.

FIX: CSV ht1/ht2 for the 20 species set to the FVS defaults. RESULT: jl's per-cycle RANN draw count now matches
FVS **EXACTLY cyc0-3** (cyc2 474=474, was 501 / +27; cyc3 1605=1605, was 1632) — the first/largest divergence
is CLOSED. Suite 5214/2 (no regression; SN uses a different CSV).

CYC4 RESIDUAL — CORRECTED VERDICT (re-trace discipline; the earlier "ULP-class" call was WRONG). A ~0.5-1%
DBH difference is NOT ULP (Float32 ULP is ~1e-7). Traced the sp27 cyc0 small-tree DG term-by-term (RNG aligned
at cyc0): dk now MATCHES FVS (0.7396, the HT-DBH fix worked); but **dgk (the LARGE-tree DG blended in) is
~0.2-0.7% LOW in jl for every sp27 tree** (d=1.2: jl 0.99137 vs FVS 0.99819; d=4.0: 0.92715 vs 0.93255;
d=1.9: 0.98785 vs 0.99229; d=0.1: 0.9557 vs 0.95722). So the residual is a REAL systematic error in the NE
LARGE-TREE diameter growth (ne/dgf.f) for sp27 — the actual #50/A2 "drift" root, NOT ULP and NOT the small-tree
path. RULED OUT: the DG coefficients B1/B2 (FVS B1(27)=.0007439 / B2(27)=.0706905 = jl's dg_b1/dg_b2 EXACTLY).
⚠ CORRECTION (re-trace discipline): the "residual = OLDRN" verdict is NOT confirmed and may be wrong. Dumped
FVS OLDRN(I) for sp27 at cyc0 — for the matched d=1.2 tree, FVS OLDRN = -0.0229 = jl's -0.02289 EXACTLY (OLDRN
MATCHES). The sp27 dgf-dump DBHs also don't cleanly line up jl-vs-FVS (multi-stand/order mixing: FVS shows
d=1.002/1.587, jl shows 0.1/4.0/1.9), so the dgk comparison that drove the ~0.5% claim may have conflated
different trees/stands. So A1's residual source is UNRESOLVED: OLDRN matches for at least one tree; the apparent
dgk gap needs a CLEAN per-tree comparison (isolate stand-1, match trees by record index, dump dgk+OLDRN+ssigma
for the SAME tree both sides). It may be smaller than the ~0.5% the conflated comparison suggested. The
SYSTEMIC IFOR=3 fix stands (draw count exact cyc0-3, independently verified); only the fine residual's
attribution is now open. NEXT: clean stand-1-isolated per-record dgk comparison before any further claim.

BOTTOM-REACHED HYPOTHESIS (now under review): A1 residual = calibration-time RNG draws for sp27's OLDRN seed. jl dump: sp27 has
dg_cor[27]=0 (UNcalibrated) and varied per-tree oldrn (0.079/0.020/-0.023/0.036/...) ⇒ these are the random
BACHLO z draws (the `else` branch, diameter_growth.jl:487 `z = bachlo(s.rng, 0, sigma[sp])` with the DGSD·sigma
rejection bound), NOT a regression line. So sp27's serial-correlation seed OLDRN is a stochastic z drawn at
CALIBRATION (LSTART, before the cycle loop). It diverges from FVS because the calibration-time RNG stream —
consumed BEFORE the cyc0-3 window the IFOR=3 fix aligned — is not yet matched to FVS dgdriv.f (the species-sorted
order/count/bound of the OLDRN z draws). CONFIRMED sp27 is UNcalibrated on BOTH sides (jl dg_cor[27]=0 via the fn[sp]>=fnmin gate; FVS COR(27)=0) ⇒
OLDRN=random z on both, so it is NOT a calibrated-flag bug — it is genuinely the calibration-time z-draw
RNG alignment. FIX PATH: apply the same draw-counter method to the LSTART calibration
(instrument the bachlo z draws in calibrate_diameter_growth! vs ne/dgdriv.f), align the per-species draw
order/count, so sp27's OLDRN seed matches. This is the FULL A1 chain end-to-end: user-flagged thin → RNG
divergence → IFOR=3 HT-DBH (FIXED) → residual large-tree DG ~0.5% (sp27) → dgf predictor faithful → serial-corr
exp(frmt) → OLDRN seed → calibration-time bachlo z RNG alignment. Deepest layer; well-scoped.

FINAL PIN (back-solved from matched dgk/DDS, no new instrumentation): for d=1.2 (dib=1.104, DDS=3.2558 both),
exp(frmt) = (((dgk+dib)^2 - dib^2)/DDS): jl frmt=-0.0261, FVS frmt=-0.0167 (Δ~0.0094). Since ssigma≈0.093 and
FM=-0.14228, the deterministic FM·ssigma·rhocp term is bounded at -0.0132 (rhocp≤1) — too small to reach
jl's -0.0261 alone. So the per-tree OLDRN (DG calibration residual, frmt = FM·ssigma·rhocp + corr·OLDRN) MUST
contribute and is where jl/FVS diverge for sp27. ⇒ A1's final residual = the per-tree DG-calibration residual
OLDRN (seeded from net01.tre's measured DG) computed slightly differently for sp27 vs FVS dgdriv.f. This is the
deepest stochastic-DG layer; the fix is aligning the OLDRN measured-DG-residual computation. END OF TRACE.

WK2 SPLIT (decisive): dumped the dgf DDS predictor (exp(wk2)) for sp27 cyc0 both sides — MATCHES bit-close
(dbh 0.1→1.09096 both; 1.2→3.25577; 1.9→4.48258; 4.0→7.80665). So the dgf PREDICTOR is FAITHFUL. The ~0.5%
dgk gap is purely the SERIAL-CORRELATION factor exp(frmt) applied after (frmt = FM·ssigma·rhocp + corr·oldrn).
Further ruled out: SIGMAR(27)=.093 = jl dg_resid_sd; jl DOES have the SIGMA sample-variance calibration (line
455 = dgdriv.f:552); FVS COR(27)=0.0 ⇒ no calibration fired for sp27 ⇒ SIGMA=SIGMAR=.093 on both ⇒ ssigma
matches. So the residual is in the REMAINING serial-correlation components for sp27: the per-tree residual OLDRN
(the calibration residual seeded from measured DG) or the ARMA multipliers (vmlt/corr/rhocp from autcor). This
is the deepest stochastic-DG layer (dgdriv.f serial correlation), ordered-work #2/#3. EXHAUSTIVE TRACE COMPLETE:
DDS predictor + B1/B2/B3/SITEAR/COR/BKRAT/SIGMAR all match; residual isolated to exp(frmt)'s oldrn/ARMA terms.

DEEPER TRACE (all flat coefficients now ruled out): for sp27 cyc0, term-by-term — B1(27)=.0007439, B2(27)=
.0706905 (= jl), B3(27)=.016240 (= jl dg_b3, BALMOD), SITEAR(27)=71.456 (= jl), COR(27)=0.0 (= jl), BKRAT(27)=
.920 (= jl bark; NE BRATIO ignores D/H and returns the per-species constant BKRAT). ALL MATCH. Yet jl's dgk
(large-tree DG) is ~0.2-0.7% low. So the residual is NOT a flat coefficient — it's in the value dgk CARRIES
beyond the dgf predictor: jl's tripling-deterministic serial-correlation factor exp(frmt) (frmt = FM·ssigma·
rhocp + corr·oldrn; FM=-0.14228 matches; ssigma from VARDG; oldrn the per-tree calibration residual) and/or the
10x annual iteration arithmetic (Float32 op-order). This is the central NE DG variance/serial-correlation
(dgdriv.f), the deepest layer of ordered-work #2. Next: dump jl's wk2 (dgf DDS predictor) vs FVS WK2 at cyc0 —
if wk2 matches, the gap is purely the serial-correlation/iteration factor; trace VARDG(27)/ssigma there.

REMAINING source (the NE DG predictor `POTBAG=B1·SITEAR·(1-exp(-B2·D))·0.7` → BALMOD → 10x annual iterate → DDS
→ `WK2=log(DDS)+COR`): SITEAR (site index), COR (DG calibration), BALMOD (B3), or the iteration. This is the
central NE DG model (ordered-work #2) not yet bit-exact for sp27. Next: trace the dgf predictor terms (SITEAR
/COR/BALMOD) for one sp27 tree vs ne/dgf.f.

CYC4 PER-RECORD (REGENT records dumped both sides, stand-1 cyc3): both have 45 records but different
per-species survival — jl sp9=15/sp27=21/sp30=9 vs FVS sp9=17/sp27=19/sp30=9. The 2 extra jl sp27 records sit
at 4.95-4.99″ (just below the 5″ REGENT cutoff) where FVS pushed them just PAST 5″ (excluded). So the IFOR=3
fix reduced sp27's small-tree DBH error from ~4-5% to a small ~0.5-1% RESIDUAL at the boundary — flipping 2
boundary trees, which changes the bachlo rejection count (the -12). This residual is approaching ULP / minor-
coefficient territory (a remaining sp27 term in the small-tree path: NC-128 height curve `ne_htcalc`, BALMOD
b3, bark_ratio, or regent_min_diam). Much smaller than the systemic IFOR=3 bug. Continue the method if pushing
to bit-exact; the per-tree small-tree DG is now within the ~1% general NE growth tolerance (A2).

CYC4 LOCALIZED (FVS phase-split counter, both sides): FVS ICYC4 DG = 2553-1605 = 948 = jl's 948 EXACTLY ✓;
REGENT+post FVS 183 vs jl 171 → the -12 is in REGENT again (4 small-tree records), NOT dgscor. The sitset.f
IFOR=3 block overrides ONLY HT1/HT2 (20 each) and is the ONLY IFOR-conditional block — so it was the single
SYSTEMIC coefficient bug (fixed), not a class. The cyc4 REGENT -12 is therefore a subtler VALUE-DRIFT boundary
effect: a few records near the 5" REGENT threshold shifted by accumulated per-tree value differences from
cyc0-2 (draw COUNT matches cyc0-3, but per-tree VALUES can still differ from a coefficient that doesn't change
the draw count). Pinning it needs a per-record DBH dump at cyc3-REGENT-entry (jl vs FVS) to find the 4 records
+ their species + the diverging coefficient (NC-128 height curve / regent DIAM budwidth / BALMOD / another
HT-DBH species). Next iteration.

NEXT-ITERATION DATA (cyc4 = jl cyc3, the new first divergence after the fix): phase-counter breakdown
shows jl cyc3 = DG(dgscor) 948 draws + REGENT 171 = 1119, vs FVS ICYC4 1131 (Δ-12). It's now in a STOCHASTIC
cycle (post-tripling), so the -12 is likely downstream of a VALUE divergence in cyc0-2 (matching draw COUNT
cyc0-3 does NOT guarantee matching draw VALUES — a coefficient that shifts a value without changing the draw
count would drift the state, surfacing as a cyc3 rejection-count difference). Stand-1 .sum @2090 after the fix:
TPA 105 vs live 111, QMD 18.4 vs 17.9, TCuFt 7605 vs 7456 (+2%) — the .sum is slightly WORSE than pre-fix
(110/17.9/7456→...) because the now-correct coefficients removed the wrong-coefficient coincidence; the cyc4+
residuals dominate. Per doctrine #3 the faithful fix STAYS; iterate the draw-counter method on cyc4+ (find the
phase, then the term — dgscor value-divergence or another REGENT/coefficient table) to drive the residual to 0.

REMAINING (A1 not fully closed): stand-1 draw counts re-diverge by smaller amounts at cyc4+ (Δ-12 at cyc4,
wobbling to ~Δ40 by stand-1 end) — SEPARATE, smaller bugs to fix the SAME iterative way (find first draw
divergence → trace the term → fix). The stand-2 .sum 2130 thin has NOT yet converged (it's downstream of the
accumulated cyc4+ residuals + the thinning record-survival question). FOLLOW-UPS: (a) the IFOR=3 Allegheny
override must be re-added as a CONDITIONAL (forest_idx==3) so Allegheny stands stay faithful — currently jl
uses defaults for ALL IFOR (right for IFOR≠3, now wrong for IFOR=3; no IFOR=3 test exists so no regression);
(b) iterate the draw-counter method on the cyc4+ divergences.

## A2 — Full NE species-set cycle-0 volume/crown/density (VERDICT: FAITHFUL)

net01 only exercises ~6 of 108 NE species, so the per-species volume (R9 Clark cubic + R9LOGS board feet),
crown (CWCALC), and density coefficients were largely unvalidated. Built a multi-species differential: rewrote
net01.tre's 30 records' species codes (cols 34–35) to span ALL 108 species in 3 batches, ran cycle-0 jl vs
live. Result — TPA/BA/SDI/CCF/TopHt/QMD **exact** every batch; stand volume within ULP/~0.5%:

| batch (species) | LIVE TCuFt/MCuFt/SCuFt/BdFt | jl |
|---|---|---|
| 1 (30 incl. BF/WP/SM/HI/AB/oaks) | 1551/1286/186/1023 | 1550/1284/185/1026 |
| 2 (39 incl. ashes/YP/spruces/oaks) | 1549/1285/184/1037 | 1557/1290/185/1036 |
| 3 (39 incl. birches/elms/OH)       | 1520/1241/103/515  | 1515/1237/103/512  |

VERDICT: FAITHFUL — the CSV-driven per-species coefficients are loaded correctly for the full species set.
Minor residual: batch-2 TCuFt 0.5% (likely one ash/oak species' cubic coefficient at single precision over a
30-tree aggregate; below the bar to chase). Pinned in test_net01.jl. Coverage 6 → 108 species (cycle-0).

Multi-species GROWTH also checked (batch-1, 30 species, cyc1-2): 2000 TPA/BA/SDI/QMD exact (524/105/209/6.1),
TopHt Δ1 (70 vs 71), TCuFt 2368 vs 2362 (0.25%); 2010 within ~1% (TPA 471/469, BA 133/134, TCuFt 3328/3293).
⇒ the per-species diameter/height-growth coefficients track live for the early cycles; the 2010 +1% TCuFt is
the earliest trace of A1 (tripling over-dispersion) beginning to bite. Per-species DG/HTG = FAITHFUL.

## A3 — Mid-run keyword paths net01 omits (VERDICT: FAITHFUL)

Injected each into net01 stand-1 (unthinned control) at 2000, diffed jl vs live:
- THINBBA (resid BA, from below), THINABA (from above), THINSDI (resid SDI) — within cyc-1 drift (test).
- THINHT (thin by height < 40 ft): 2010 TPA 136 EXACT, BA Δ1, QMD Δ0.1 — faithful.
- THINCC (resid CCF): jl tracks live (2010 TPA 479 vs 475 = cyc-1 drift).
- FIXMORT (mortality override, all-species rate 0.20): 2010 TPA 420 / BA 116 **BIT-EXACT** vs live.

VERDICT: FAITHFUL — the cut-selection (cuts.f) and mortality-override paths drive correctly off the NE
growth/volume. The one cut-related divergence (A1) is upstream in the DG distribution, not in these paths.

## A4 — Establishment / regen over the full projection (VERDICT: FAITHFUL)

net01 stand-5 (BARE + ESTAB/PLANT) over 9 cycles (1992→2082), jl vs live: TPA EXACT (800/777/755/…),
BA within Δ2, TCuFt 6% early (2012, the planted seedlings' tiny-DBH volume) converging to <1% by 2082
(6481 vs 6462). The PLANT path, regen growth, and density all track live end-to-end. #49 had validated
only cycle-0; this extends it to the full run. VERDICT: FAITHFUL.

## A1 supporting verification (record pruning is faithful — ruled out as the lever)
jl's per-cycle record pruning matches FVS: `comcup!` is called every cycle at simulate.jl:336 (= FVS
grincr.f:391) with the PROB ≤ 1e-5 threshold (= base/comcup.f:50), and `tredel_compact!` runs after each
thin (= TREDEL swap-from-end). So the 230-vs-165 record-count gap at cyc15 is NOT a missing-prune bug; it
is the per-record MORTALITY ALLOCATION differing (matching total TPA, but FVS fully empties more records
→ removed, jl leaves more at low-but->1e-5 TPA), which feeds back into the A1 RNG-stream drift. Every
discrete component of A1 is faithful; the residual is the COMPOUND stochastic evolution under repeated thins.

## Validated-faithful so far (breadth, vs live FVSne)
- Cycle-0 stand state (TPA/BA/SDI/CCF/QMD/TopHt) — bit-exact.
- R9 Clark cubic + International-¼″ board-feet volume — bit-exact (#48).
- THINBBA / THINABA / THINSDI thin-selection — within cyc-1 drift (test_net01.jl).
- FFE crown-fire subsystem (CBD, fuel moisture, FMCFMD, crowning/torching, SNAGINIT) — bit-close (#47).
- `.sum` -999 header variant code — FIXED (was stamping "SN" for all variants).

## Audit queue (subsystems to sweep, exercised by a live differential)
- Diameter growth distribution (A1) — IN PROGRESS.
- Mortality per-DBH-class distribution (does the A1 over-dispersion also implicate mortality skew?).
- Height growth / crown ratio across the full NE species set (net01 has ~6 of 108 species).
- Volume across the full NE species set (only net01's species are validated).
- Establishment / regen (net01 stand-5 BARE; #49 cycle-0 only).
- Mid-run keyword paths not in net01: THINATA/THINCC/THINHT, FIXMORT, SETSITE, FERTILIZE, species multipliers.
