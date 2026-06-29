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

## A1 — Tripled-record DBH over-dispersion under REPEATED thinning (OPEN, bounded) ★ flagged by user "−40 vs −22 BA not acceptable"

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

**Localized to the diameter-growth record TRIPLING.** Decisive `NOTRIPLE`-on-both test at 2130:

| 16–20″ stocking | TRIPLE | NOTRIPLE |
|---|---|---|
| LIVE | 9.5 | 10.8 |
| jl   | **27.0** | **8.9** |

Under NOTRIPLE jl and live agree (both below the 15 target → no cut). Under TRIPLE live stays at 9.5 but
jl balloons to 27. jl also carries 230 records at cyc15 vs live's 165.

**What's ruled out.** The tripling spread factors `FU=1.271 / FM=−0.14228 / FL=−1.549`
(src/.../diameter_growth.jl DG_FU/FM/FL) are BIT-identical to NE `ne/dgdriv.f:626-628`. The tripling
cadence (ICL4=2, first 2 growth cycles) matches NE `ne/grinit.f:183`. Per-record DG + mortality are fine
(NOTRIPLE matches). So the bug is the **stochastic evolution of the tripled records in cycles 3–15**: the
high-`old_random` upper records (seeded rnU = FU·ssigma·rhocp + corr·rnpar at tripling) run away too wide
through the serial-correlation `dgscor!` (AR1/BACHLO) path, where FVS keeps them tight and prunes/compresses
back to 165 records. Suspects: (a) jl not pruning/merging emptied tripled records the way FVS does over many
cycles; (b) the dgscor! AR1 persistence on the seeded rnU diverging from FVS's RNG stream.

**Localized further (16–20″ CSTOCK per thin cycle, live debug dump vs jl histogram):**

| thin yr (cyc) | LIVE | jl |
|---|---|---|
| 2020 (3)  | 0.0  | 1.3  |
| 2050 (6)  | 7.1  | 10.0 |
| 2080 (9)  | 7.9  | 21.8 |
| 2110 (12) | 35.2 | 15.7 |
| 2130 (15) | 9.5  | 27.0 |

They MATCH at cyc3 (2020), right after tripling — so the initial tripling spread is fine. They diverge from
cyc6+ and OSCILLATE OUT OF PHASE (the 16–20″ class fills then is thinned to 15 each cycle; jl's fill/empty
timing leads/lags live), with jl averaging ~28% higher (15.2 vs 11.9 over the 5 thins). ⇒ the bug is the
**post-tripling stochastic DG evolution of the 9× tripled records** (dgscor! AR1 on the seeded upper/lower
`old_random`, walking the BACHLO RNG stream), NOT the tripling spread and NOT per-record DG (NOTRIPLE matches).
The 9× record count makes the NE dgscor!/BACHLO RNG stream draw 9× and drift from live FVSne's stream across
15 cycles — SN aligned this via the Oracle-A transliteration; NE has no such reference, so the multi-cycle
tripled-record RNG stream was never bit-aligned. Mild systematic over-dispersion (+28% avg) rides on the
phase noise, so it is not purely RNG-phase.

**Candidate fixes ruled out (all faithful):** the tripling factors FU/FM/FL (= ne/dgdriv.f:626-628), the
tripling cadence ICL4=2 (= ne/grinit.f:183), per-record DG+mortality (NOTRIPLE matches), the dgscor! per-cycle
bound (jl replicates dgscor.f's `|frm|>DGSD·SIGMA` redraw + the dds>4 taper), and DGSD=2 (= NE grinit). The
OLDRN init clamp (dgdriv.f:634, ±DGSD·SIGMA) is a one-time calibration step, and an unclamped tripled-record
seed self-bounds after one dgscor cycle — not the systematic lever.

**Remaining cause = RNG-stream misalignment of the tripled records.** With 9× records jl's NE BACHLO/dgscor
draw stream diverges from live FVSne over 15 cycles; the count mismatch (jl carries 230 records at cyc15 vs
live's 165 — jl does not prune emptied records the way FVS record-management does) changes the per-cycle draw
count, guaranteeing stream divergence. SN bit-aligned this against the Oracle-A transliteration; NE has no such
reference, so the multi-cycle tripled-record RNG stream was never aligned.

**Scope NARROWED — it's the tripling × repeated-thinning interaction, not pure growth.** The UNTHINNED
stand-1 tracks live faithfully through all 10 cycles (2090: TPA 110/111, BA 192/194, SDI 274/279, QMD 17.9
exact, TCuFt 7538/7456 = 1.1%, BdFt 42978/43258 = 0.6%). So the tripled-record growth is fine in isolation;
the over-dispersion only appears in stand-2, which thins to per-DBH-class targets EVERY 3 cycles (5 thins over
16 cycles). Mechanism: each THINDBH removes a slightly different set of records in jl vs live (post-thin TREDEL
compaction order + which partial-cut records survive → a different surviving-record RNG stream), and over 5
repeated thins this compounds into the 16–20″ class divergence. Same class as the SN COMPRESS work (TREDEL
compaction + RDPSRT order, already done for SN). **Practical severity is bounded**: it requires repeated
class-target thinning of one stand; normal NE projection (stand-1 unthinned, single thins A3) is faithful.

**Next (task #50).** Trace the post-thin record survival + compaction order (cuts.jl TREDEL path) for the
THINDBH class-target case against FVS cuts.f/TREDEL, and bit-align which records survive a partial cut + their
order, so the surviving-record RNG stream tracks FVSne across repeated thins. Focused multi-step effort; the
one substantive OPEN item, now bounded to the repeated-thinning scenario.

**Status: OPEN.** This was previously (wrongly) closed as "faithful within drift" — a lax verdict the user
correctly rejected. Re-opened.

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
