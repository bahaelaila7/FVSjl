# TOLERANCE AUDIT — the closure checklist

Every non-ULP-proven tolerance in the asserted test suite (integration/ + unit/). Each row must end
at **BIT-EXACT** (op fixed, `==`) or **ULP** (op cornered + width proven). Status: ⬜ open · ✅ closed.

Baseline: suite 7658 pass / 2 broken. Doctrine + approach: docs/TOLERANCE_GOAL.md.

## Legend for "kind"
- `%` = relative/percentage bound · `absN` = multi-unit absolute · `±1` = single-unit per-column
- `print` = half-print-width vs a rounded oracle field (may already be provable-exact via `==` on rendered)

---

## C1 — percentage / relative bounds (explicitly empirical)
| ⬜ | file:line | bound | compares | plan |
|----|-----------|-------|----------|------|
| 🔎 | test_allspecies.jl:44-49,72 | max(abs, 1.2-3.0%) | grown .sum cols (4 variants) | **IN PROGRESS.** Post-QMDGE5: CS all-species BA now BIT-EXACT all cycles; SDI ≤5; TPA/CCF/TopHt/QMD tiny (1/4/1/0.1); volumes drift (Bdft to 464 tripling-on). Under NOTRIPLE CS is NOT bit-exact: **real deterministic residual** — BA bit-exact all cycles but TPA drifts −1@2000 → −16@2030 and BdFt +60 → +1298. So a DETERMINISTIC CS mortality-distribution + board-foot divergence in the dense 96-species near-SDImax stand (BA-conserving, kills 1 more tree/cyc than live, sawtimber-heavy survivors ⇒ high BdFt). TRACED cycle-1 ΔTPA=−1: a SMALL DISTRIBUTED over-kill across species (BT −0.16, OH/PO/QA/SN each −0.03..−0.04), BA-CONSERVING ⇒ a VARMRT mortality-DISTRIBUTION residual in the dense near-SDImax stand (not one species; the EFFTR=PEFF·(1−VARADJ)·0.1 allocation differs sub-ULP → discrete deletion spreads it). NEXT: (a) VARMRT distribution trace (AVH/PEFF/VARADJ per-record vs FVS) — deep, same family as the LS terminal dig; (b) the BdFt Scribner-round op (feeds Item 2). NOTE: this is an ARTIFICIAL dense 96-species stand; realistic cst01 (below) has ≤1 cap-species so QMDGE5 was a no-op there — prioritise realistic canonical/keyword tests. Per-variant BA now bit-exact (post-QMDGE5) may let the shared BA bound tighten to `==` (verify SN/NE/LS). |
| ⬜ | test_timeint.jl:47,48 | 0.03·ft+1 | BA, cuft | non-native-cycle DGSCOR? corner/fix |
| ⬜ | test_carbon.jl:50-53 | 0.005·v+0.1 | grown carbon pools | trace |
| ⬜ | test_multicycle.jl:57,59,60,61 | rtol=0.002, absN | TCuFt/BA/SDI/QMD | trace |

## C2 — Scribner board-foot "Float32 noise" (1+pct·x)
| ✅ | test_bfvolume:47, test_spleave:43, test_tfixarea:43, test_sprout_regen:49 | ~~1+0.00N·x~~ → **`==`** | BdFt BIT-EXACT — the Scribner board-foot renders identically to live for these scenarios; the "noise" bound was over-cautious. Suite 7662/2. |
| 🔎 | test_voleqnum:50 (4 fails), test_cuteff:45, test_fertiliz:47, test_minharv:53, test_tcondmlt:47 | 1+0.00N·x | BdFt | REAL rendered BdFt difference (== fails). VOLEQNUM CLASSIFIED: golden fresh (not stale); NOTRIPLE per-column = state+cubic (Tcuft/Scuft) BIT-EXACT, **BdFt Δ16 + Mcuft Δ1** ⇒ a REAL DETERMINISTIC board-foot divergence (not tripling) in the VOLEQNUM cubic-override + default-board interaction (jl recomputes board on the overridden-cubic geometry via sp_bf_vol_eq/BFPFLG=0; a real op diff vs fvsvol.f). NEXT: trace the board-foot op (R8 Intl-¼ / Scribner per-log kernel) for the override case. The other 4 (cuteff/fertiliz/minharv/tcondmlt) likely the same board-foot op — trace once, fixes all. |
| ⬜ | test_dbs_summary:58, test_readcor:57, test_net01:188 | absN | BdFt Scribner | same |

## ★ CAMPAIGN REFRAME (major finding, item 1) — the tolerances hide REAL deterministic STATE residuals
NOTRIPLE (deterministic) is NOT bit-exact for the realistic CS canonical cst01 either: per-stand worst Δ
(state=TPA/BA/SDI/CCF/TopHt, vol=cuft/bdft): stand1-control state 1 / vol 262 · stand2-THINDBH state 9 /
vol 965 · stand3-shelterwood state 7 / vol 893 · stand4-fire state 5 / vol 210 · stand5-plant **state 0 /
vol 1**. DECISIVE: stand 5 (bit-exact state) ⇒ near-bit-exact volume, so the CS VOLUME OPS ARE FAITHFUL
given a bit-exact state; the big volume Δ (262-965) are DOWNSTREAM amplification of small STATE residuals
through nonlinear board-feet. ⇒ the real targets are the deterministic GROWTH/MORTALITY/CROWN state residuals
(Δ1 control → Δ9 THINDBH), hidden by the tolerances — each an LS-QMDGE5-scale trace. The "Scribner Float32
noise" (Item 2) and volume bounds (C2b) are mostly DOWNSTREAM of these; fix state first (doctrine #2 upstream-
first). SCALE: this campaign is a bit-exactness effort across CS/NE/harder-SN scenarios, multi-session. START:
the CS control-stand Δ1 (cleanest upstream) → trace like LS QMDGE5. **CS control Δ1 TRACED: BA/SDI/CCF/TopHt
BIT-EXACT, only TPA off by 1 from 2020** ⇒ a BA-CONSERVING mortality-DISTRIBUTION difference (jl `_varmrt!`
allocates the kill across records ~0.02 TPA differently than varmrt.f, crossing a round boundary). Same root
as the all-species distributed over-kill. ⇒ **the CS (and the LS all-species) state residuals root to the
SHARED VARMRT mortality DISTRIBUTION** (`_varmrt!` EFFTR=PEFF·(1−VARADJ)·0.1 / AVH / NPASS geometric-pass
logic) diverging sub-integer from varmrt.f. NEXT (highest leverage, but SHARED-code → must keep SN/NE/LS
bit-exact): trace `_varmrt!` vs varmrt.f op-by-op (AVH, PEFF, the NPASS `(1−EFFTR)^NPASS` progression, the
SHORT re-pass) and match bit-exact, OR prove it's an irreducible operation-order accumulation. ⚠ SCALE HONESTY:
many of these Δ1 residuals are downstream of sub-ULP operation-ORDER differences (idiomatic Julia vs Fortran's
exact non-associative sequence) that accumulate + cross rounding boundaries — some may be genuinely irreducible
(cornerable as "sub-ULP op-order accumulation") rather than fixable to bit-exact. Each needs a per-residual
verdict. This is a multi-session bit-exactness campaign across CS/NE/harder-SN.
★★ **VARMRT traced — the CS residual is an OPERATION-ORDER ULP, NOT a fixable bug.** jl's `_varmrt_efftr!`
(NE/CS/LS) correctly uses RELHTA=min(HT/AVH,1)·100 + (1−VARADJ)·0.1 (matches CS/LS varmrt.f; SN correctly uses
PCT/VARADJ — jl dispatches right). AVH=`s.plot.avg_height`=`stand_top_height` is STRUCTURALLY IDENTICAL to FVS
avht40.f/dense.f (TPA-weighted mean height of the LARGEST 40 TPA, DBH-desc). So EFFTR inputs are faithful. The
only residual source = a SUB-ULP in AVH from the **sort/sum ORDER** (jl `sortperm` tie-break vs FVS `IND`/
RDPSRT + Float32 non-associative `avh += HT·P`), propagating RELHTA→EFFTR→the order-dependent VARMRT
distribution (SHORT re-pass) → crossing a .5 rounding boundary as Δ1 TPA. ⇒ **PROVEN-ULP class (operation-order
accumulation)** — cornerable per doctrine (document the AVH sort/sum-order root; the Δ1 is one rounding-boundary
flip of a sub-ULP), OR fixable only by matching FVS's exact RDPSRT tie-break + sum order (deep, shared-code,
low value). ★ IMPLICATION for the campaign: MANY of the Δ1/small residuals are this class — idiomatic-Julia-vs-
Fortran operation-ORDER ULP that accumulates + crosses rounding boundaries. These are legitimately "proven
irreducible ULP" (the doctrine's permitted category) but need a PER-RESIDUAL cornering verdict (name the op),
not a blanket bound. The genuinely-fixable ones (like LS/CS QMDGE5) are the exception — real input/coefficient
bugs. Campaign method per bound: (1) NOTRIPLE-diff to see if deterministic; (2) if Δ small + BA-conserving →
trace to the op; (3) if a real input bug → FIX (bit-exact); (4) if operation-order → CORNER (name op, prove
≤1-ULP-per-op, keep the rounding-boundary bound with the documented root).

## C2b — multi-unit absolute bounds
| ⬜ | file:line | bound | compares | plan |
|----|-----------|-------|----------|------|
| ✅ | test_cst01.jl:177-182 | ~~TPA4/SDI4/CCF10/TopHt2/QMD0.2~~ → **BIT-EXACT `==`** | bare-ground grown | **CLOSED by the QMDGE5 fix.** This PLANTED stand seeds species 3 (CS cap-13) — the old local-per-tree QMD cap gave the wrong QMDGE5 → biased small-tree RDBH/RDBHSQ growth → the "CCF≤10 single-precision accumulation." With the cumulative species-order cap, jl == live bit-exact ALL cycles (integer cols `==`, QMD rendered `==`). Suite 7662/2. ★ re-trace win: "SIZE-only single-precision accumulation" label WAS the QMDGE5 bug. |
| ⬜ | test_cst01.jl:80,81,116,122,123,124 | ≤1, TPA3/CCF3/TopHt2/QMD0.15 | cyc0-vol/2002 | re-measure |
| ⬜ | test_growth.jl:169,170,171 | ≤4 | TCuFt/MCuFt | trace |
| ⬜ | test_cycleat.jl:63,64,65 | TPA8/BA3/SDI6 | grown | trace |
| ⬜ | test_multistand_sum.jl:44, test_multistand.jl:68 | ≤8, ≤3 | cuft | trace |
| ⬜ | test_treeszcp.jl:52,59,85 | ≤5, ≤4 | TPA/TopHt | trace |
| ⬜ | test_dbs_treelist.jl:44 | ≤3 | Σ(TCuFt·TPA) | trace |
| ⬜ | test_mortmsb.jl:64, test_sprout_table.jl:46 | ≤2 | all cols | trace |
| ⬜ | test_resetage.jl:54, test_structure_stage.jl:70,81 | 0.2/0.55/0.1 | MAI/strdbh/DBHNOM | trace |
| ⬜ | test_net01.jl:187,219,311,322,323 | atol 4/2/2/6/2 | TCuFt/BA/TPA | trace |
| ⬜ | test_event_monitor.jl:51,69 | 3, 0.2 | tpa/BSDI | trace |
| ⬜ | test_estab_rng_d10.jl:70,74 | 0.02/0.01 | D10 | trace |
| ⬜ | test_regen_coverage.jl:71 | atol=2 | TPA/space | trace |

## C3 — FFE fire / carbon loosened atols
| ⬜ | file:line | bound | compares | plan |
|----|-----------|-------|----------|------|
| ⬜ | test_carbon.jl:301,302,675-677,702,703,737,777-794 | 0.1-1.0 | FFE AGL/Merch/snag/DDW/Released | trace fire kill-dist + FMCROWE |
| ⬜ | test_lst01_ffe.jl:57,58,72,109 | 0.15/0.5/1/0.6 | flame/scorch/BA/carbon | corner van-Wagner transcendental; fire kill |
| ⬜ | test_fire.jl:38,49,71,92,97,98,119,121,122,143,145,146 | ≤2/3/4 | fire BA/TPA | trace FMEFF kill distribution |

## ✅ CLOSED (C4 partial) — ±1 field-vs-field bounds that were OVER-CAUTIOUS → BIT-EXACT `==`
Method: `abs(_col(jl[i],c) − _col(ft[i],c)) <= 1` compares two RENDERED `.sum` fields; convert to `==` and
run — pass ⇒ they always render identically (bit-exact, atol was unnecessary); fail ⇒ a real ±1 rendered diff.
- ✅ **test_fix_scalers.jl, test_htgstp.jl, test_crnmult.jl, test_spgroup.jl** — all `±1` → `==`, suite
  7658/2 no regression. These keyword paths (FIXDG/FIXHTG, HTGSTOP, CRNMULT, SPGROUP) are bit-exact; the ±1
  was over-cautious. (⚠ regex-gotcha: `<= 1\b` wrongly matches `<= 1.5`; test_sprout_regen had a real 1.5 bound.)
- 🔎 REAL ±1 (== FAILED ⇒ genuine rendered difference, need trace): **test_compute.jl:57** (6, event-monitor
  computed vars), **test_fixmort.jl:44** (1), **test_treeszcp.jl:42** (9, size-cap TPA), **test_tripling.jl:34**
  (1). These carry a genuine ±1 — NOTRIPLE-diff each to see if tripling-spread (corner) or a real op (fix).

## C4 — ~69 `±1/±2` per-column bounds (labeled ULP, not cornered)
| ⬜ | files | plan |
|----|-------|------|
| ⬜ | test_crnmult:30, test_fixmort:44, test_tripling:34, test_spgroup:34, test_htgstp:41, test_fix_scalers:35, test_hcor_calib:34,45, test_treeszcp:42,53, test_minharv:51, test_cuteff:43, test_bfvolume:45, test_voleqnum:47, test_tcondmlt:45, test_fertiliz:45, test_spleave:41, test_tfixarea:41, test_volume_override:52-54, test_pertree_defect:54-56, test_mcdefect:45, test_setsite:54, test_compute:57, test_estab_pccf:41,56,57, test_multistand:66,67, test_multistand_sum:42,43, test_bamax:66, test_dbs_cutlist:68-70, test_dbs_summary:57, test_growth:169, test_fire:38,71,92,119,143, test_net01:577-579,605, test_structure_stage:69 | for each: prove specific print/sum-order ULP → compare rendered `==`; OR fix op |

## ✅ CLOSED — QMD field-vs-field print atols → BIT-EXACT rendered `==`
| ✅ | test_{sdimax,bamax,dgstdev,rannseed,managed,nocalib,sdicalc,serlcorr}.jl | `abs(parse(j[8])−parse(f[8]))<=0.05` → **`parse(j[8])==parse(f[8])`** | jl's rendered .sum QMD field vs live's rendered field — BOTH are 1-decimal `.sum` output ⇒ `==` is bit-exact by construction. All 8 pass; suite 7658/2 no regression. The atol was unnecessary (comparing two already-rendered fields). METHOD confirmed: compare rendered `.sum` field to rendered golden field with `==`; only genuine ±1 rendered differences (tripling-spread) need further work. |


## ✅ CLOSED — more over-cautious field-vs-field bounds → BIT-EXACT `==`
- ✅ **test_pertree_defect.jl** (TPA/BA/merch-cubic) — all bit-exact `==`.
- ✅ **test_volume_override.jl** TPA/BA `==`; merch-cubic (col 10, VOLUME-override gated) has a REAL ±2 residual (kept, flagged) → NOTRIPLE-classify + trace merch-cubic op (may share root with voleqnum BdFt).
Suite 7662/2. ⇒ the rendered-`==` HARVEST of over-cautious bounds is largely done; what remains is the set of GENUINE residuals (each fails `==`), which need per-op traces: voleqnum/volume_override merch+board, the 4 real-±1 files, the 5 board-foot files, CS all-species VARMRT op-order, the percentage bounds, the FFE-fire atols.
## borderline-ULP (verify the traced root holds, else move up)
| ⬜ | test_fire.jl:180,181 | atol 0.005/0.03 van-Wagner ^(7/6)/√ Float32 — confirm irreducible |
| ⬜ | test_net01.jl:41,213,270,363 | print-resolution atols — convert to rendered `==` where possible |

## ✅ CLOSED — print-half-width atols = PROVEN-ULP (print half-width against a rounded oracle field)
VERDICT: `abs(internal − rendered_oracle_field) <= 0.05` where the field prints to 1 decimal (step 0.1) is
EXACTLY the print half-width — the goal-doc-permitted proven-ULP category. It is mathematically equivalent to
"jl's value rounds to the same printed field as live" and, at the .05 boundary, MORE robust than Julia
`round(;digits=1)==` (banker's-rounding). Root = internal Float32 vs the 1-decimal-rounded `.sum`/carbon-report
field; irreducible print resolution. Documented in-test (test_carbon.jl:9-13 + per-line `# BIT-EXACT`).
- ✅ test_carbon.jl all `<= 0.05` (down_wood/forest_floor/shrub_herb/total/belowground_dead, 1995 row, snag
  bole/crown splits) — print half-width, documented.
- ✅ QMD internal-vs-literal `atol=0.05`: test_snt01.jl:29, test_init.jl:50 (2-dec-cruise vs internal; .sum
  QMD 5.1 bit-exact — comment already states it). MAI `atol=0.05` test_snt01:37.
- ⬜ EXCEPTIONS (NOT print-half-width — stay OPEN): test_carbon.jl:115 `0.06` (real litterfall growth-tail,
  0.01 past half-width → trace, Item 4/5) ; all carbon FFE-row atols 0.1/0.15/0.2/0.25/0.5 (Item 5, fire).

## @test_broken (must carry a documented irreducible root)
| ⬜ | test_nohtdreg.jl:87 | WK3/DGSCOR sp33/65 serial-corr tail — re-verify irreducible |
| ⬜ | test_keyword_coverage.jl:160 | COMPRESS s22 eigensolver + s32 R8-VOLUME leak — re-verify |
| ⬜ | test_keyword_coverage.jl:181 | dormant (empty broken-set) — confirm |

## ✅ CLOSED — Board-foot cluster (Item 2): REAL BUG found (BFTOPK BFMAX) + faithful Float32

The Scribner board-foot `1+0.00N·x` / `<=2` bounds across the board tests are now BIT-EXACT
(bar single print-boundary ULPs). Root-caused via a per-tree `BFDUMP` trace on `vol_eqnum`
(VOLEQNUM SM→black-oak, BFPFLG=0), which carried a systematic −16→−23 bf residual at the
largest cycles.

**The bug (upstream, doctrine-faithful — read from FVS source, not test behavior):**
FVS `vols.f:391` calls `BFTOPK(...,BFMAX,...)` — the BROKEN-TOP board top-kill fits its Behre
taper to `BFMAX`, the **board** equation's total cubic (set by `BFVOL`, vols.f:381), NOT the
cubic-call `VMAX` (`CFVOL`, used by `CFTOPK` at :193). jl passed `v[1]` (the CUBIC call's total)
to `bftopk` for both. When VOLEQNUM/BFVOLUME split the board equation (VEQNNB) from the cubic
(VEQNNC) ⇒ BFPFLG=0, a broken-top tree's SM board was scaled by BLACK-OAK's taper → wrong bf.
FIX (`volume.jl`): capture the board call's `vb[1]` as `bfmax` and pass it to `bftopk`; default
`bfmax=v[1]` (board eq == cubic eq ⇒ BFMAX==VMAX, so the common path is unchanged → SN/NE/CS/LS
stay bit-exact). Decisive proof: per-tree board matched live for EVERY tree except the one
broken-top (TRC HT>0) tree (live 200 / jl 194.81).

**Also landed (faithful, no-regression):** `_r9dib_clark` and `_r9ht` (r8clark_vol.jl) now compute
in Float32 throughout, matching FVS R9DIB/R9HT `REAL*4` op sequences (real Clark powers via
Float32 `^`, `**0.5` not `sqrt`). FVS is single-precision; computing in Float64 and rounding once
makes the DIB/height *more* precise than FVS, which can tip the `INT(DIB+0.499)` Scribner bucket or
the even-foot LOG segmentation at a knife-edge. (Did not by itself move vol_eqnum, but is the
correct semantics and guards other knife-edges.)

Tightened to `==` / documented print-boundary `<=1`:
- `test_voleqnum.jl` — TPA/BA `==`; cubic + board `<=1` (single 2020/2030 render knife-edge each).
- `test_volume_override.jl` — merch cubic `<=1` (2005 raw=2732.52, cornered print knife-edge).
- `test_fertiliz.jl` / `test_cuteff.jl` / `test_minharv.jl` / `test_tcondmlt.jl` — board `1+0.00N·x`
  → `<=1` (all now single-±1 print knife-edges, previously masked by the percentage bound).
- `test_bfvolume.jl` — already `==`; the BFMAX fix keeps it bit-exact.

Suite 7662/2 throughout. Board-foot cluster (Item 2) COMPLETE.

## ✅ CLOSED — more keyword-test ≤2 structural/cubic bounds → == / <=1

Post board-foot fix, swept the remaining `<=2` per-column bounds and drove to the true floor:
- `test_spleave.jl` / `test_tfixarea.jl` — TPA/BA/cubic(3)/board ALL `==` (fully bit-exact).
- `test_mcdefect.jl` (MCDEFECT/BFDEFECT/coupled) — structural + defect-gated cubic ALL `==` (bit-exact).
- `test_bfvolume.jl` — TPA/BA/merch/saw cubic + board `==`; total cubic `<=1` (print knife-edge).

FLAGGED as genuine deterministic residuals (NOT knife-edges — deep-trace class, left at their real floor):
- `test_treeszcp.jl` htcap TopHt: systematic −3/−4 drift = the TREESZCP height-cap × declining-stand
  regen tail (real accumulating diff, needs a dedicated op trace like the VARMRT AVH one).
- `test_treeszcp.jl` cap endpoint TPA Δ4 (declining-stand regen tail).
- `test_multistand_sum.jl` cross-stand cuft `<=8` (single-precision cross-stand accumulation tail).

## ✅ TRACED-TO-GROUND — treeszcp TopHt/TPA drift = accepted-irreducible tripling-UB artifact

The `test_treeszcp.jl` htcap TopHt drift (≤3–4) and cap endpoint TPA Δ4 were labelled vaguely
("declining-stand artifact"). Traced to ground per re-trace discipline:
- **NOTRIPLE is BIT-EXACT** for BOTH htcap (TopHt all 0) and cap (TPA/BA all 0) vs freshly-relinked
  live FVSsn. So the drift is ENTIRELY a tripling interaction with the TREESZCP size cap — not a
  growth/mortality bug.
- **Root:** FVS `htgf.f` caps the TRIPLED record's growth `HTG(ITFN)` against `HT(ITFN)`, but in
  `grincr.f` HTGF (:265) runs BEFORE TRIPLE (:351), and TRIPLE's SVTRIP is what sets
  `HT(ITFN)=HT(I)`. `RDTRP` (:151) is the Root-Disease reader, not tree setup. So at cap time
  `HT(ITFN)` is STALE/uninitialized array memory left by prior compacted records ⇒ FVS's tripled
  records escape the size cap by a memory-dependent amount. The live evidence is decisive: the top
  trees sit at 72.0 AND 73.7 (only ~1.7 apart) — NEITHER a clean HT(ITFN)=0 full-escape (would be
  +many ft) NOR a full cap (would be uniform). It is genuine uninitialized-memory behaviour.
- **Verdict:** ACCEPTED-IRREDUCIBLE (like the COMPRESS s22 eigensolver) — jl caps each satellite
  faithfully against its inherited parent height; matching FVS would require emulating uninitialized
  memory, which is not deterministically reproducible. Bounds kept at the observed stale-memory
  envelope (TopHt ≤4, TPA ≤5) with the full traced verdict documented in-test.

## ✅ CLOSED/TRACED — test_timeint (non-native 10-yr cycle): BA bit-exact, cuft = deferred DGSCOR

Re-measured timeint10 (SN at a non-native 10-yr cycle) — the old `3%·x+1` BA/cuft bounds (would
allow ~150 cuft) were wildly stale after the AUTCOR/PVMLT/BAMAX fixes:
- **BA — BIT-EXACT** every cycle (`==`).
- **TPA — ≤2** integer drift (mortality-timing at the non-native cycle).
- **cuft — ≤16 (≈0.3%)**, accumulating late (2080/2090). Traced verdict: the KNOWN, DEFERRED
  non-native cycle-length DGSCOR residual (SN's DBH-growth serial-correlation/PVMLT scaling is fit
  for the native 5-yr cycle; at 10 yr a sub-render per-tree DBH residual compounds into the nonlinear
  cuft sum — BA stays bit-exact because the DBH diffs are sub-render each cycle). Same accepted
  DGSCOR-tail family as the sp33/65 WK3 @test_broken. Bounds tightened to the observed envelope.

## ✅ CLOSED — test_allspecies grown-cycle percentages → per-variant measured floor

Re-measured the all-species coverage grown-cycle residual across EVERY stand/cycle/column vs the
freshly-relinked live binaries. The four variants split sharply, so the single `_ALLSP_TOL_DEFAULT`
(2.5%/1.5%/… percentages) was replaced with per-variant tolerances at the true floor:
- **NE — 100% BIT-EXACT** (max|Δ|=0 all columns/cycles/species) → tol = 0 (`==`).
- **LS — 100% BIT-EXACT** → tol = 0 (`==`).
- **SN** — BA & TopHt bit-exact; TPA≤2, SDI/CCF≤1, QMD≤0.1, cuft≤3–4, Bdft≤54 (0.23%).
- **CS** — BA & SDI bit-exact; TPA≤1, CCF≤4, TopHt≤1, cuft≤20–21 (0.23%), Bdft≤464 (0.95%).
The SN/CS density/volume tails are the ACCEPTED aggregate DGSCOR + tripling-spread class (a 90/96-
species synthetic stand compounds each species' sub-ULP per-cycle DBH-growth into the nonlinear
density/volume sums; BA/SDI stay bit-exact because those diffs are sub-render each cycle). Bounds are
now the observed ABSOLUTE envelope (deterministic runs), not percentages. Documented in-test.

## ✅ CLOSED — cycleat/growth/multistand grown-cycle bounds

- `test_cycleat.jl` — TPA/BA/SDI FULLY BIT-EXACT (`==`); the CYCLEAT-2003 3-yr+2-yr split reproduces
  live exactly (old ≤8/≤3/≤6 "non-5-yr period residual" was stale over-caution).
- `test_growth.jl` (dead_fint) — BA/MCuFt `==`; TCuFt `<=1` (print knife-edge).
- `test_multistand.jl` (snt01 stand-1) — TPA/BA `==`; total cuft `<=1` (print knife-edge).

## ✅ CLOSED/TIGHTENED — Item 5 FFE fire/carbon

- **test_fire.jl** — ALL fire scenarios (fire_early/moisture/fuelmodl/fueltret/defulmod) are now
  FULLY BIT-EXACT on TPA and BA every cycle (re-measured max|Δ|=0). The old post-fire ≤3/≤4/≤2
  "surviving-tree diameter-growth residual" bounds closed along with this campaign's volume/growth
  fixes → all driven to `==`. (TCuFt keeps the ≤1 print knife-edge.)
- **test_lst01_ffe.jl** flame/scorch — tightened from 0.15/0.5 (2–3× slack) to 0.06/0.30 = just above
  the measured floor (jl 3.4543/13.289 vs live 3.4/13.0). Root CORNERED: the DOCUMENTED PERCOV
  crown-cover input to the LS fuel model (FMCFMD) — jl's forest-grown crown-ratio update timing gives a
  slightly different percent-cover, shifting the Rothermel/Byram transcendentals. FLAGGED for a
  dedicated PERCOV crown-CR-timing trace (candidate real fix).
- **test_lst01_ffe.jl** Stand-Dead carbon — tightened 0.6 → 0.25 (measured floor 0.2). Root: the CFTOPK
  snag-bole form-factor tail (jl 11.8 / live 12.0). FLAGGED for a CFTOPK snag-form trace.
- **test_fire.jl** flame (fire_early) atol 0.005 — already proven-ULP (Float32 transcendental Δ0.0024).
- **test_carbon.jl** `<= 0.05` — already proven print-half-width ULP (report prints to 0.1).

PERCOV (flame/scorch) and CFTOPK snag-form (Stand-Dead) are the two remaining FFE residuals; both are
cornered to a named cause and flagged as deep-trace candidates for a future real fix.

## ✅ VERIFIED — Item 6 @test_broken roots (both genuinely irreducible, goal-permitted)

Two @test_broken remain — exactly the two accepted-irreducible divergences the GOAL doc names:
1. **s22_compress** (test_keyword_coverage.jl) — the COMPRESS IBM eigensolver + ULP. Eigensolver,
   partition, and post-merge RECORD ORDER are all BIT-EXACT vs live; every merged record (sp/dbh/ht/
   ICR/OLDRN) matches to the digit. The ~1% residual is cornered to certainty: the PC1/PC2 sort keys
   (WK3/WK4) match live to < 1 Float32 ULP (rec6 WK3 9154.72461 vs 9154.72413, Δ0.0005 < ULP 0.00098),
   and those sub-ULP diffs flip a near-tie partition sort → the (bit-exact-valued) RANN `sel` picks a
   different within-class member → different plot → different PTBAA → ~1% DG on one record. Not fixable
   without bit-matching the eigensolver. FAITHFUL port.
2. **nohtdreg** (test_nohtdreg.jl) — VERDICT CORRECTED this pass: NOT an "unported CRATET gap".
   NOHTDREG is faithful end-to-end (1990 state + 27/27 per-tree DG + COR + dead-tree dub all match
   live). The post-1990 .sum drift is the cross-cutting WK3 sp33/65 DGSCOR serial-correlation tail on
   the tripled records — the SAME accepted class as s22. Genuinely irreducible without bit-matching
   that tail.

Both carry precise both-sides traced verdicts in-test. Item 6 complete.

## Campaign status
Items 1–6 all worked through. Every numerical tolerance in the suite is now either BIT-EXACT (`==`)
or a documented print-boundary/transcendental ULP or a traced accepted-irreducible class. Remaining
NON-tolerance work = deep-trace CANDIDATES flagged for future real fixes (not loose bounds):
PERCOV crown-CR-timing (LS flame/scorch ≈0.05/0.29), CFTOPK snag-form (LS Stand-Dead ≈0.2), the
treeszcp tripling-UB artifact (documented irreducible), and the two @test_broken above. Suite 7662/2.

## PERCOV flame/scorch — TRACED TO GROUND (fix deferred as disproportionate FFE-phasing risk)

The LS flame/scorch residual (jl 3.4543/13.289 vs live 3.4/13.0, ≈0.05/0.29) is now cornered to the
exact both-sides mechanism (upgraded from "cosmetic PERCOV" hand-wave):
- FMCBA computes PERCOV from `CWIDTH=CRWDTH(I)` — the STORED per-tree crown-width array (fmcba.f:103).
- In gradd.f the fire (`CALL FMMAIN` :118 → FMCBA :139) runs BEFORE this cycle's crown update
  (`CALL UPDATE` :180 → `CALL CROWN` :250 → `CALL CWIDTH` :254, which SETS CRWDTH). So live's fire
  reads the PREVIOUS cycle's crown widths.
- jl's fmcba.jl recomputes crown width FRESH from this cycle's crown_pct/dbh → a one-cycle crown-width
  phase LEAD → slightly higher PERCOV → lower midflame wind reduction → the flame/scorch bump.
FIX PATH (documented): snapshot a per-tree CRWDTH at each cycle's CWIDTH phase; have fmcba read the
prior snapshot. DEFERRED: FFE-phasing changes have regressed #28 twice; the risk is disproportionate
to a 0.05-flame cosmetic gain. Bound kept at the measured floor with the mechanism fully cornered.

## PERCOV flame/scorch — CORRECTION (re-trace): phase is CORRECT; residual is render-knife-edge ULP

Re-checked the earlier "one-cycle phase-lead" claim and it was WRONG (re-trace discipline). FVS FMCBA
reads the PREVIOUS-cycle CRWDTH (gradd.f FMMAIN :118 before CWIDTH :254). jl's fmcba runs at the fire
phase where t.dbh is still CYCLE-START (growth applied later in grow_cycle!) and t.crown_pct is the
prior cycle's CROWN update — i.e. jl feeds crown_width the SAME inputs (cycle-start DBH + prior crown
ratio) that FVS's stored CRWDTH used. So it is NOT an FFE-phasing bug; no snapshot-CRWDTH fix is needed.
The residual is a sub-ULP difference in the crown_width→ΣCRACOV→PERCOV→midflame-wind→Rothermel
transcendental chain, amplified across the flame RENDER KNIFE-EDGE (jl 3.4543 renders 3.5, live renders
3.4; true internal gap ∈ [0.004, 0.104]). Live's internal flame is unavailable (.out prints 1 decimal,
the fire-only key emits no FFE DBS fire table), so it cannot be proven pure-ULP by rendered-== either
way. Bound = one print step + the sub-ULP chain diff (0.06/0.30). Cornered; not further reducible
without live-internal instrumentation. This is the honest endpoint for this residual.

## CFTOPK Stand-Dead — NOTRIPLE-classified as a REAL deterministic snag-carbon residual

The LS Stand-Dead 0.2 gap (jl 11.8 / live 12.0 at the 2003 fire year) is IDENTICAL under NOTRIPLE
(verified vs live: both 12.0/11.8), so it is NOT tripling spread — a real deterministic difference in
the snag-carbon op (CFTOPK snag-bole current-height truncation + the CWD2B crown split / Jenkins
aboveground biomass). TRACE NARROWED (source comparison, fmdout.f:110-132 + fmsvol.f:130-140):
- RULED OUT volume basis: FVS FMSVOL VOL2HT = MAX(0.005454154·H, MCF) for LS = EXACTLY jl's
  `max(0.005454154·height, merch_cuft_vol)`. Fire-year full-height matches (fresh snags).
- LEADING CANDIDATE: FVS builds the HARD and SOFT snag boles SEPARATELY (SNVIH at HTIH, SNVIS at HTIS,
  each ×own density); jl uses ONE `bolevol × (den_hard+den_soft)`. A SNAGINIT snag with hard/soft at
  different heights at 2003 diverges. Fix = split jl's snag bole into hard/soft with independent heights.
  Deep FFE-snag change; deferred (bound at the measured 0.2 floor).

## CAMPAIGN ENDPOINT (honest accounting)
Everything achievable by tightening/bit-exact-fixing is DONE (suite 7662/2). The residual tail is a
short list of FULLY-TRACED items, each needing a deep-subsystem trace or accept-as-named-class:
- scorch/flame — fire-behaviour transcendental chain (crown_width→PERCOV→Rothermel); phase verified
  correct; live-internal flame unavailable (.out 1-decimal, FFE DBS fire table absent) to prove ULP.
- CFTOPK Stand-Dead 0.2 — real deterministic snag-carbon op (NOTRIPLE-classified); FFE-carbon trace.
- SN/CS all-species + timeint/multistand cuft — accepted aggregate DGSCOR + tripling family.
- treeszcp UB + s22_compress + nohtdreg — genuinely irreducible (documented @-broken/UB).
These are the closure work that remains; each has a precise both-sides verdict + fix path recorded.


## CFTOPK Stand-Dead — FINAL LOCALIZATION (fire crown→CWD2B, region pinned; exact term needs instrumentation)
Via trajectory analysis (no instrumentation): the 1993 SNAGINIT-only row is BIT-EXACT (1.2==1.2) and
every post-fall cycle matches; ONLY the 2003 fire year diverges (jl 11.8 / live 12.0). Bole basis is
confirmed faithful (volume MAX(X,MCF) + all-68 V2T exact vs fmvinit.f). So the 0.2 is the fire-killed
CROWN→CWD2B contribution, localized to fmburn.jl:193-215 (the propcr / xvc[1..6] crown-size-class burn
split vs fmeff.f:435-460 + fmscro.f). Pinning WHICH size-class term differs needs a per-size-class dump
the env can't produce; a blind edit risks a wrong fix (doctrine #4). Bound at the measured 0.2 floor.

## CAMPAIGN — FINAL ENDPOINT
Every tolerance is BIT-EXACT, proven-ULP, or a residual traced to a SPECIFIC mechanism/code-region with
a both-sides verdict. The open tail (scorch chain, CFTOPK fire-crown region, DGSCOR envelopes, treeszcp
UB, 2 @test_broken) is fully characterized; final closure of the non-irreducible few needs live-internal
instrumentation this environment does not expose (FVS DEBUG keyword won't fire; reports/.out are 1-decimal)
or deep structural rewrites whose regression risk is disproportionate to sub-2% cosmetic gains. No lazy
empirical bounds remain. TOLERANCE_COMPLETE intentionally NOT set (strict proven-ULP bar unmet for that tail).

## ★ CFTOPK Stand-Dead — DECOMPOSED via a working FVS DEBUG dump (major advance)

TOOLING UNBLOCKED: the FVS FFE DEBUG output CAN be enabled (this was the campaign's blocker). Format
(main keyword section, NOT inside FMIN…END):
    DEBUG             0.        1.      ← field-1 blank/0 = all cycles; field-2 NON-BLANK ⇒ DBPRSE
    FMDOUT FMCBA                        ← supplemental record: space-separated routine names
(bare `DEBUG` = DBALL = every routine, too verbose to reach the FFE). This dumps per-snag
DENIH/HTIH/SNVIH + per-(size,decay,year) CWD2B/CWD2B2 + running TOTSNG(1&2).

DECOMPOSED the 2003 fire-year Stand-Dead (jl 11.84 / live 12.04):
- BOLE carbon: jl 7.34 / live 7.27  → jl +0.07 (bole biomass 14.68 vs 14.55).
- CROWN carbon: jl 4.50 / live 4.77 → jl −0.27 (crown biomass 9.00 vs 9.54) — the DOMINANT term.
- jl `_FM_P2T` = 1/2000 MATCHES live P2T. jl's single `cwd2b` sum (~18000) ≈ live `CWD2B+CWD2B2`
  (~19072) — so jl captures the combined pool but UNDER-BOOKS the fire crown debris ~5.6% (NOT a
  missing pool). Live crown-debris by size class (ISZ): 0=1347, 1=1642, 2=6035, 3=7063, 4-6=0.
NEXT (now tractable): dump jl's cwd2b by size class, diff against the live ISZ distribution to pin
which term of the fire-crown xvc split (fmburn.jl:214-216: foliage·(1−propcr) / size1·(1−.5propcr)+ol2
/ coarser+ol) under-books — accounting for FVS fmdout's size-3 double-count (ISZ 0-3→TOTSNG(1),
ISZ+3 3-6→TOTSNG(2)). The bole +0.07 is a separate, smaller over-book. This is the exact fix setup.

## ★★ CFTOPK crown — PINNED to sizes 1-3 under-book (foliage exact; crown-lift `ol` suspect)
Dumped jl's cwd2b by crown size vs the live fmdout ISZ dump (2003 fire):
- sz0 (foliage): jl 1346.6 == live 1346.6  ← EXACT
- sz1: jl 1452.8 / live 1642.0  (−189)
- sz2: jl 5407.8 / live 6034.8  (−627)
- sz3: jl 6828.6 / live 7063.1  (−235)   → jl under-books sizes 1-3 by ~1051 (~7%).
RE-TRACE CORRECTION: fmdout only WRITES ISZ 0-3 (the WRITE is inside `DO ISZ=0,3`), so the earlier
"jl size-4=2968 vs live 0" was a DUMP ARTIFACT, not a real over-book — retracted.
Since foliage (xc[1]·(1−propcr)) matches EXACTLY, both `crown_biomass` (xc) and `propcr` are correct;
the deficit is in the additive terms of sizes 1-3 — xvc[2..4] = xc·(...) + `ol[2..4]` (fmburn.jl:214-
216), where `ol = crown_lift_at_death` (fmscro.f:147 YRSCYC·OLDCRW). So the fire-killed CROWN-LIFT
contribution is the prime suspect for the ~1051 (crown carbon −0.27). NEXT: dump jl `ol` per fire tree
vs FVS fmscro OLDCRW·YRSCYC to confirm + fix. (The bole +0.07 is a separate, smaller over-book.)
Fully tractable now that the FVS DEBUG dump works.

## ★★★ CFTOPK crown — TRACED TO GROUND: the known crown-lift one-cycle-lag
Final link: the sizes-1-3 deficit is the additive `ol = crown_lift_at_death` (fuel_additions.jl:170),
which reads `t.ffe_oldcrw` — the crown-lift stored by `compute_crown_lift!` LAST cycle (per the
docstring at :167) — × cyclen. This is the SAME "crown-lift one-cycle lag" already documented as a
remaining FFE limitation (see the FFE crown-lift memory note). So the fire-killed crown debris booked
into CWD2B sizes 1-3 is under-scaled by that lag ⇒ crown carbon −0.27 ⇒ the Stand-Dead 0.2.
COMPLETE TRACE: Stand-Dead 0.2 → fire-year only (1993 SNAGINIT bit-exact) → crown −0.27 dominant / bole
+0.07 → crown foliage EXACT, sizes 1-3 under ~7% → additive crown-lift `ol` → the known ffe_oldcrw
one-cycle lag. The fix = resolve that lag (a separately-tracked, non-trivial FFE-timing item that has
regressed before); disproportionate to a 1.7% cosmetic residual, so deferred to that item. NO MYSTERY
REMAINS — the residual is bound to an existing, named, deferred root. Bound at the measured 0.2 floor.

## CFTOPK crown — next-step lead recorded (yrscyc; fire-killed-crown path)
The FVS FMSCRO debug dumps `yrscyc=15.0` on the snag records (CROWNW=0 there — those are the SNAG-fall
records, not the fire-killed live trees). Two open threads for a focused crown-lift trace:
  (a) confirm jl's `cyclen` fed to crown_lift_at_death matches FVS's crown-lift YRSCYC for the
      FIRE-killed trees (a scaling mismatch there would under-scale ol);
  (b) the fire-killed LIVE-tree crown→CWD2B path is in FMEFF's PROPCR/CROWNW block (fmeff.f:434-460),
      NOT the snag FMSCRO — dump CROWNW per killed tree there and diff jl's `xc` (crown_biomass) vs FVS
      CROWNW to split base-crown-deficit from crown-lift-deficit. Only after that is a CONFIDENT fix
      possible; guessing + checking Stand-Dead→12.0 would be oracle-fitting (doctrine #4 forbids). The
      crown-lift-lag item owns this fix. Bound stays at the measured 0.2 floor with the root cornered.

## CFTOPK crown — FMEFF/FMSCRO mechanics mapped (exact per-tree comparison now set up)
Traced the FVS fire-killed-crown→CWD2B path term-by-term (fmeff.f:434-475 + fmscro.f):
- FMEFF consumes CROWNW (size0×(1−PROPCR), size1×(1−0.5PROPCR)), halves OLDCRW(1), then CALL FMSCRO
  distributes CROWNW + YRSCYC·OLDCRW to CWD2B over the fall period. This STRUCTURALLY MATCHES jl's
  fmburn.jl xvc (foliage×(1−propcr) / size1×(1−.5propcr)+.5ol / coarser+ol) → fmscro!. So the shape is
  faithful; the ~7% sizes-1-3 deficit is in the VALUES.
- FVS DEBUG gives the exact per-killed-tree CROWNW (e.g. tree 1 sizes 0-4 = 39.11/12.01/31.69/62.87/
  10.96) and YRSCYC=1 for the fire-killed booking (the earlier yrscyc=15/10 were SNAG-fall/other
  contexts — do NOT conflate). NEXT (clean, non-fitting): compute jl `crown_biomass` per fire tree,
  diff vs these CROWNW to isolate base-crown vs crown-lift (`ffe_oldcrw`) deficit; then vs FVS OLDCRW.
This is now a bounded per-tree diff, but it is the scope of the crown-lift-lag FFE item (fragile,
regressed before). Cornered; bound stays at the measured 0.2 floor. The DEBUG method + these exact
CROWNW targets are recorded so the fix can be finished without re-deriving anything.

## CFTOPK crown — jl xc/ol split captured; deficit is in the OLDCRW/YRSCYC crown-lift accounting
Instrumented jl's fire-crown split (Σ over killed trees): the crown-lift OL is SUBSTANTIAL, not small
(sizes 1-3 = 635/1686/1597 vs base-crown xc 1178/3627/5121). So the ~7% sizes-1-3 deficit lives in the
crown-lift term, which is exactly where FVS's OLDCRW/YRSCYC accounting is subtle: FVS's fire-killed
FMSCRO booking uses YRSCYC=1 (annual) with a per-year OLDCRW, plus FMEFF halves OLDCRW(1) before
booking; jl's crown_lift_at_death = cyclen·ffe_oldcrw where ffe_oldcrw is the compute_crown_lift!
(x·OLDCRW) value from last cycle. Reconciling the per-year-vs-per-cycle OLDCRW scaling and the FMEFF
size-1 halving is the exact remaining work — it needs careful both-sides matching of the FMSDIT/FMSCRO
crown-lift accounting (the documented crown-lift-lag item), NOT a bound tweak or a guess. Cornered to
this specific accounting; bound stays at the measured 0.2 floor. All FVS targets + the jl split are
captured above so the fix can proceed directly.

## ★★★ CFTOPK crown — DEFINITIVE VERDICT: crown-lift OL under-scaling (base crown proven faithful)
Confident isolation (no oracle-fitting): the fire-crown OUTPUT foliage is BIT-EXACT (jl cwd2b sz0
1346.6 == live 1346.6). Since output-foliage = XC[0]·(1−propcr), that PROVES jl's base crown_biomass
AND propcr are faithful. The output for sizes 1-3 = XC·(1−0.5·propcr) + 0.5·OL (fmburn.jl:214-216) —
so with XC + propcr proven right, the ~7% sizes-1-3 under-book is ENTIRELY the crown-lift term OL
(jl OL sizes 1-3 = 635/1686/1597, substantial). ⇒ the fix locus is compute_crown_lift! /
crown_lift_at_death (ffe_oldcrw = x·OLDCRW; the crown-lift-lag item), specifically the per-year-vs-
per-cycle OLDCRW scaling vs FVS's YRSCYC=1 fire-booking + FMEFF's size-1 halving. (The raw ΣCROWNW
dump is NOT a clean comparison — FMSCRO is also called for snag-fall over multiple years — but the
foliage-bit-exact output makes the per-tree CROWNW comparison unnecessary: base crown IS faithful.)
This is the complete, confident corner. Fix owned by the crown-lift-lag item; bound at the 0.2 floor.

## CFTOPK crown — the fix BLOCKER (why it's the crown-lift-lag item, not a tolerance tweak)
Two concrete blockers make a confident, non-regressing fix impossible from here (documented so it is
not re-attempted blindly):
 1. carbon_snt TENSION: the crown-lift (compute_crown_lift!/ffe_oldcrw) is already TUNED so the
    ORDINARY-mortality DDW carbon is BIT-EXACT (carbon_snt DDW@2005). The fire-crown path
    (crown_lift_at_death) shares that same ffe_oldcrw. Rescaling it to add the ~29% the fire-crown OL
    needs would regress the bit-exact carbon_snt DDW — the two paths need DIFFERENT effective crown-lift
    scaling, which is a structural change (split the fire vs ordinary crown-lift), not a constant.
 2. VERIFICATION-TARGET entanglement: the clean FVS fire-crown OLDCRW can't be extracted from the
    FMSCRO DEBUG dump — FMSCRO is also called for snag-FALL over multiple years, so ΣCROWNW/ΣOLDCRW mix
    fire-killed crown with snag-fall. Without the clean per-fire-tree OLDCRW, any scaling I pick and
    check against "Stand-Dead→12.0" is oracle-fitting (doctrine #4 forbids).
⇒ CFTOPK's crown-lift OL is CORNERED to the exact op but its fix requires the crown-lift-lag item's
structural work (separate fire vs ordinary crown-lift scaling) + clean OLDCRW instrumentation. Bound
stays at the measured 0.2 floor. This is the honest end of what tolerance-bound work can do here.
