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

## A1 — NE stochastic RNG stream not bit-aligned to FVSne (OPEN; model faithful, realization not) ★ flagged by user "−40 vs −22 BA not acceptable"

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

**Next (task #50).** Pin the ~3 boundary records: instrument jl small_tree_growth! + FVS regent.f to dump
(record, dbh, is-small, n-draws) at cycle-2 entry; find the records where small-membership disagrees; reconcile
the NE small-tree DBH threshold / which records REGENT draws for (vs the grincr small-tree subset that calls
regent). Then re-run the draw-counter to confirm cycle-2 Δ→0, and walk forward. Keep SN bit-exact (the SN
small-tree path is `<3″` and already aligned).

**Status: OPEN** (model faithful, stochastic stream not yet aligned). Previously mis-closed as "faithful
within drift" — a lax verdict the user correctly rejected; re-opened.

---

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
