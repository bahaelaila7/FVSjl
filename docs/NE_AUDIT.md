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

## A1 — REFINED: the cyc3+ residual is PURELY small-tree REGENT (dgscor DG now bit-exact)

Phase counter (cyc3, after all 3 fixes): DG (dgscor) = 948 draws = FVS 948 EXACTLY — the three fixes aligned
the STOCHASTIC dgscor draw stream too, not just cyc0. The entire cyc3 Δ-18 is in REGENT (small-tree): jl 165
vs FVS 183 (~6 small records / bachlo rejections differ at the 5″ boundary). HTG/MORT draw 0. So A1's last
residual is ENTIRELY the small-tree REGENT growth at cyc3+ — the small-tree DBH near the 5″ boundary (REGENT's
own growth over cyc0-2, which determines the cyc3 small-tree set + the ±10% random rejection counts) — which
perturbs OLDRN and feeds the height (TopHt -2) + volume (+1.8%) divergences. NEXT: per-record REGENT dump at
cyc3 (small-tree DBH + draw count) jl vs FVS to find the ~6 boundary records and whether REGENT's small-tree
DBH is bit-exact through cyc0-2. The deterministic DG core (large + dgscor) is bit-exact; small-tree REGENT
is the last layer.

## A1 — UNIFIED: remaining height+volume residual is ONE thing — the cyc3+ multi-cycle RNG stream

The height residual and the volume residual are NOT separate subsystem bugs. NE height growth uses
`htg = htg1·(1+OLDRN)·gmod` (height_growth.jl:26) and the DG uses `exp(frmt)` with frmt = FM·ssigma·rhocp +
corr·OLDRN — BOTH driven by the SAME per-tree serial-correlation residual OLDRN. OLDRN is bit-exact at cyc0
(dgk bit-exact) but EVOLVES each cycle via dgscor (cyc3+ stochastic). The draw counter (after all 3 fixes)
still shows cyc0-3 EXACT then cyc3 Δ-18 (jl 2718 vs FVS 2736) — a small multi-cycle RNG-stream offset (~6
bachlo / ~2 records) at the first stochastic cycle. That offset perturbs OLDRN at cyc3+, which then feeds BOTH
the diameter (volume +1.8%) AND the height (TopHt -2) divergences. So A1's entire remaining residual = the
cyc3+ dgscor/REGENT multi-cycle RNG stream not bit-aligned to FVSne (the small-tree-boundary or mortality-
record-selection draw-count offset at cyc3). This is the SN-class multi-cycle RNG alignment (SN had Oracle-A;
NE needs the draw-counter method extended through the stochastic cycles). The three SYSTEMATIC DG fixes
(IFOR=3 + VMLT + oldp) made the deterministic core bit-exact; this stochastic stream is the last layer.

## A1 — LARGE-TREE DG NOW BIT-EXACT (3 fixes); residual is NOT DG ★ re-trace caught the z-draw over-claim

After all THREE fixes (IFOR=3 HT-DBH + calibration VMLT + cyc0 ARMA oldp), the sp27 large-tree dgk is
BIT-EXACT for EVERY tree (stand-1 isolated, cyc0): d=0.1→0.957219, 1.2→0.998194, 1.9→0.992287, 4.0→0.932551,
12.7→1.19269 — all = FVS exactly. So the "calibration z-draw realization" verdict was WRONG (re-trace
discipline, 6th catch): there is NO z-draw residual — the per-tree dgk scatter was the corr/oldp bug, now
fixed, and OLDRN matches everywhere (the fixes are GLOBAL across all species, not just sp27). ⇒ the NE large-
tree diameter growth is bit-exact per-tree at cyc0. The remaining stand-1 .sum Δ (2030 TPA 352 vs 356; 2090
109 vs 111, BA 193/194) is therefore a DIFFERENT subsystem — MORTALITY (morts, ordered-work #3) and/or the
small-tree/regen path — NOT the large-tree DG. Per-cycle decomposition (stand-1, accret/mort): mortality Δ is small (2010 14/15, 2020 37/40, 2030 77/78,
2040 75/77), but ACCRETION (volume growth) diverges MORE (2010 134/139, 2020 140/145, 2040 116/121 — Δ~5).
Since cyc0 DIAMETER growth is bit-exact, the accretion Δ is HEIGHT growth (htgf) and/or the VOLUME computation
and/or cyc1+ DG — the other growth subsystems (ordered-work #3), NOT the cyc0 diameter. So the NE residual has
moved DOWNSTREAM of the (now bit-exact) large-tree diameter growth into height/volume + a small mortality tail.
Next: trace the mortality (which trees die / the mortality RNG)
or the REGENT small-tree growth for the Δ-TPA; the large-tree DG is DONE.

## A1 — THIRD FIX LANDED: cyc0 ARMA `oldp` also hardcoded the SN 5-yr period

After the VMLT fix, ssig/vmlt matched FVS but CORR still differed (jl 0.14799 vs FVS 0.18082). Root: the
per-cycle ARMA multiplier (diameter_growth.jl) set the FIRST cycle's `oldp` (AUTCOR old period) to a hardcoded
5 — the SN measurement base — but NE's is YR=10. covmlt=AUTCOR(YR,YR).covar drives CORR; oldp=5 under-counts it
for NE. FIX: `oldp = cyc==0 ? htg_period(s.variant) : ...` (5 SN, 10 NE). RESULT: CORR=0.18082 and RHOCP=0.9834
now MATCH FVS EXACTLY, and dgk is BIT-EXACT where oldrn matches (sp27 d=4.0: jl 0.932551 = FVS 0.932551).
Variant-aware, SN bit-exact (5214/2). This is the THIRD instance of the same root (SN 5-yr measurement period
hardcoded in shared DG code; NE needs the variant YR) — alongside IFOR=3 HT-DBH and the calibration VMLT.

With ssigma/corr/rhocp now ALL bit-exact, A1's only remaining per-tree residual is the calibration OLDRN z-draw
realization (uncalibrated species), ~the stand-1 Δ2 TPA — the documented-stochastic class (NE calib RNG stream
not bit-aligned to dgdriv.f). The THREE systematic DG bugs are fixed; the model's large-tree DG is now
bit-exact per-tree wherever the stochastic z seed matches.

## A1 — SECOND FIX LANDED: NE DG calibration VMLT used the SN 5-yr measurement period

ROOT of the ~0.5% large-tree DG residual FOUND + FIXED (the re-traced verdict, after OLDRN was cleared).
calibrate_diameter_growth! hardcoded `autcor(5,5)` for the calibration VMLT — the SN measurement base period
YR=5. But NE's YR=10 (blkdat DATA YR/10.0/; SIGMAR is on a 10-yr basis; FVS dgdriv.f VMLT=VMLTYR from the
YR-period autcor). Direct FVS dump (sp27): SSIGMA=0.0930 (=SIGMA), VMLT=29.40; jl had SSIGMA=0.14954 because
its vardg used vrnext(5)≈11.15 instead of vrnext(10)=29.40 (vardg 2.6x high ⇒ ssigma 84% high ⇒ the serial-
correlation factor exp(frmt) off ⇒ ~0.5% dgk error on EVERY NE tree — the #50/A2 drift root). FIX: use
`autcor(htg_period(s.variant), htg_period(s.variant))` (5 SN, 10 NE). VARIANT-AWARE, SN bit-exact (suite
5214/2). RESULT: NE stand-1 moved toward live — 2090 TPA 105→109 (live 111), QMD 18.4→18.0 (live 17.9), SDI
274→277 (live 279). A real second NE bug fixed (after IFOR=3). Residual now much smaller; stand-2's repeated
thinning still amplifies the remaining realization difference.

FINAL RESIDUAL CHARACTERIZED (~0.1%): the post-vmlt-fix per-tree dgk differences are VARIED IN SIGN
(d=0.1: -0.15%, d=1.2: +0.05%, d=1.9: -0.10%). Both-sign per-tree scatter (not a uniform offset) ⇒ it is the
stochastic per-tree z-draw REALIZATION (the calibration OLDRN bachlo z for uncalibrated species like sp27),
NOT a deterministic op-order or coefficient error. So with the two SYSTEMATIC bugs fixed (IFOR=3 HT-DBH +
calibration VMLT period), A1's remaining ~0.1% is the calibration-time z-draw RNG realization — a
different-but-valid stochastic draw because NE's calibration RNG stream is not bit-aligned to FVS dgdriv.f
(the SN-class alignment SN got via Oracle-A; NE has no reference). This is the documented-stochastic class;
closing it to literal bit-exact needs aligning the LSTART calibration z-draw order/count (draw-counter method
on the calibration). The two systematic fixes are the substantive result; this tail is near-ULP scale.

VERIFIED the vmlt fix closes the dgk gap: post-fix sp27 cyc0 trace shows ssig=0.093 (= FVS exactly, was
0.14954) and dgk now within ~0.1% of FVS (d=1.2 jl 0.99871 vs FVS 0.998194 = +0.05%, was -0.68%; d=1.9
0.991263 vs 0.992287 = -0.10%; d=0.1 0.955828 vs 0.957219 = -0.15%). The ~0.5% large-tree DG residual is
CLOSED to near-ULP (~0.1%, Float32 op-order). Stand-1 .sum tracks live (2090 TPA 109/111). Remaining: the
near-ULP large-tree tail + the small-tree REGENT boundary (cyc4 draw Δ-18) + stand-2 repeated-thinning
amplification — all much smaller than the two fixed systemic bugs.

VERDICT TRAIL (for the record): "ULP" (wrong) → "OLDRN" (wrong, OLDRN matches) → ARMA ssigma/vmlt (RIGHT,
confirmed by direct FVS SSIGMA/VMLT dump). Two re-trace corrections; the doctrine's re-trace discipline caught
both before they stuck.

## A1 — Stand-2 thin divergence (investigation history) ★ flagged by user "−40 vs −22 BA not acceptable"

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


### A1 final localization (per-record, cyc3 stand-1 isolated)
FVS has 19 sp27 small REGENT records, jl has 18 - FVS keeps one at 4.979 in (just under 5 in) that jl grew PAST 5.
Shared DBHs ~0.5% off with VARIED sign (FVS 4.31481/jl 4.2942 low; FVS 4.46957/jl 4.4932 high; FVS 3.14229/jl
3.1442 high). Draw count matched cyc0-3 => the small-tree REGENT growth is ~0.5% off DETERMINISTICALLY (op-order
or a small REGENT-path term, not RNG), flipping one boundary tree at cyc3 and cascading the draw stream (-18).
A1 LAST residual = small-tree REGENT diameter/height growth (small_tree_growth.jl dgsm/htg) vs ne/regent.f,
~0.5%, analogous to the (fixed) large-tree DG bugs. Large-tree diameter core (dgf + dgscor) is BIT-EXACT.
NEXT: trace one sp27 small tree REGENT growth (htgr->htg->hk->dkk/dk->dgsm->dg) at cyc0 vs regent.f.


### A1 cyc0 fully bit-exact (verified all 3 growth paths)
After the 3 fixes, ALL cyc0 growth is bit-exact for sp27: large-tree DG (dgk 0.1->0.957219 etc.), large-tree
HEIGHT (htgf: d=10.0 htg 6.471574 = FVS, others within ~6e-6 ULP), AND small-tree REGENT (cyc0 htg+dg match
exactly: d=1.2 htg 10.41718/dg 0.924942 = FVS). So NO cyc0 growth term diverges. The ~0.5% cyc3 small-tree DBH
divergence (one boundary tree flip) is ENTIRELY in the cyc1-2 MULTI-CYCLE path: the tripled-record growth
(cyc0-1 tripling) or the mortality record-selection over cyc0-2, NOT any cyc0 growth model. The NE growth
MODELS (DG + HTG + REGENT) are bit-exact at cyc0; the last residual is the multi-cycle tripling/mortality
record evolution -- the SN-class record-management/RNG alignment (no NE reference). NEXT: trace cyc1 (tripled
small records) or the mortality record selection over cyc0-2.


### A1 cyc2-entry refinement (cyc0-1 result)
cyc2-entry DBHs (cyc0+cyc1 growth): the 1.9-2.2 in small-trees MATCH FVS to 4 digits (FVS 1.94516/jl 1.9457,
2.09953/2.1001) but the 3.2+ in small-trees diverge ~0.7% (FVS 3.20675/jl 3.2294). The larger small-trees blend
MORE with the large-tree dgk (xwt=(d-1.5)/3.5: 0.49 at d=3.2 vs 0.14 at d=2.0), so the residual is the
dgk-blend at cyc1-2 for HIGHER-xwt small-trees (or trees crossing 5 in mid-projection). cyc0 growth is fully
bit-exact (all models); this is a fine-grained cyc1-2 multi-cycle small-tree-blend residual that accumulates
to one tree flipping 5 in at cyc3, cascading the draw stream. NEXT: dump the cyc1 small-tree dg for a
high-xwt (d~3-4) sp27 tree vs regent.f -- check the dgk blend / the small tree crossing 5 in mid-cycle.
The three systematic DG fixes made cyc0 bit-exact; this last layer is the higher-xwt small-tree blend over
the tripling cycles.


### A1 cyc1 finding (NOT uniformly bit-exact at cyc1)
cyc1 sp27 dgk (large-tree DG blended into small trees): the d~4.98 tree MATCHES (FVS 0.90038/jl 0.90038), but
the d~3.0 tree DIVERGES ~4% (FVS dgk 0.976285 vs jl 0.93751). Same tree (~3.0), real ~4% gap. So the cyc1
large-tree DG diverges for SOME small trees, while cyc0 was fully bit-exact. dgk = bsc(sqrt(dib^2 +
DDS*exp(frmt)) - dib); DDS depends on the dgf COMPETITION term (point-BAL pbal = pba*(1-crown_ratio/100) +
stand BA), frmt on OLDRN. So the cyc1 divergence is the dgf competition (point density PTBAA / crown ratio PCT)
or OLDRN at cyc1 - downstream of the cyc0 stand state (mortality record selection / crown / tripling).
NEXT: dump cyc1 dgf inputs for the d~3.0 sp27 tree (pba/PTBAA, crown_ratio/PCT, oldrn, DDS, ba_v) jl vs
ne/dgf.f. cyc0 is bit-exact; cyc1 dgf competition for small trees is the next layer.


### A1 cyc1 wk2 split (op-order/ULP floor reached)
cyc1 wk2 (dgf DDS predictor) for d~3.0 sp27 trees: the d=2.9095 tree is BIT-EXACT (FVS wk2 1.81522226 / jl
1.81522), but the d~3.0 trees have slightly different cyc1-ENTRY DBHs (FVS 2.98971/3.00604 vs jl 2.9924/3.0021,
~0.1%) so their wk2 differs because the INPUT DBH differs. These are TRIPLED small-tree records (U/L sub-
records); their cyc0 tripling DG (dgU/dgL via the REGENT ±10% random for small trees, or the dgk-blend) diverges
~0.1% while the CENTRAL records are bit-exact. So the cyc1 dgf itself is faithful (wk2 bit-exact for a tree whose
input DBH matches) -- the divergence is the ~0.1% cyc0 tripled small-tree U/L DBH, at the Float32 op-order/ULP
floor. The earlier "cyc1 dgf ~4% dgk" was the COMPOUND of this ~0.1% entry-DBH diff through the iteration +
boundary; the per-input dgf is faithful. So A1's residual is now at the ULP/op-order floor for the tripled
small-tree records -- the documented-ULP class. The deterministic NE growth (dgf + dgscor + htgf + REGENT
central) is bit-exact; only the tripled small-tree U/L sub-records carry a ~0.1% op-order tail.


### A1 correction: the tripled U/L residual is ~0.1% (NOT confirmed ULP)
Re-trace discipline on my own prior note: a ~0.1% difference is ~1000x Float32 ULP (ULP~1e-6 relative), so it is
NOT pure ULP and I should not imply it might be. The tripled small-tree U/L sub-records carry a REAL ~0.1%
difference at cyc1-entry (FVS 2.98971/3.00604 vs jl 2.9924/3.0021). The central (l=0) records and the dgf are
bit-exact, so the divergence is specifically the U/L (l=1,2) sub-record growth. Since htg=htgr+ran*0.1*htgr,
a ~0.1% htg/DBH diff implies the ±10% random `ran` for the U/L sub-records differs ~1% -> the REGENT bachlo
draw VALUES for the tripled sub-records are misaligned (the draw count matched cyc0-3, so it is a draw-ORDER /
sequence difference, not a missing/extra draw). So A1's residual is the REGENT small-tree tripled-record bachlo
draw-ORDER vs ne/regent.f -- a real (small) RNG-sequence misalignment for the U/L sub-records, NOT ULP and NOT
yet fixed. VERDICT: OPEN, ~0.1% per-U/L-record, the SN-class RNG-order alignment for the small-tree tripling.
NEXT: instrument the REGENT bachlo draw sequence (which sub-record consumes which draw) jl vs regent.f for one
sp27 small tree at cyc0.
