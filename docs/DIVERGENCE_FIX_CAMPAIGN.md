# FVSjl — Non-ULP Divergence Fix Campaign (2026)

Drive every live-FVS divergence that is **neither Float32 ULP nor the accepted COMPRESS
eigensolver** to ULP-class, or prove it irreducible and document why. Oracle = live FVS
(sn/ne/cs_oracle.sh; debug-stamp the .f, relink, restore). Doctrine: trace logic both
sides → upstream-first → validate vs LIVE before writing the test → document every
verdict → variant-aware (gate, don't harden; keep all three variants bit-exact).

Status: ⬜ open · 🔬 investigating · ✅ fixed-to-ULP · 📌 irreducible/deferred (why documented)

**★ Ledger state (all D1–D12 catalogued items resolved):** D1 not-real · D2/D7/D9/D12 fixed-bit-exact ·
D4/D5 carbon-report bit-exact · **D8/D10 fixed-to-ULP** (the ~51% regen sawtimber divergence — the campaign's
biggest — traced to two real establishment `:estab` RNG bugs, live-validated, see D10 below) · D3 & D6
📌 evidence-backed deferrals (NOT hand-waved): D3 is empirically ULP-class for EVERY corpus scenario
(bare_natural/plant/multipoint/mp3 all ≤4.6% Bdft = threshold-amplified Float32 ULP; per-tree DBH ≤0.012″),
its genuine gap (per-point density on HETEROGENEOUS multi-point stands) has no validating oracle scenario;
D6 (CS ESCPRS regen-compression) only fires on establishment list-overflow — no corpus scenario triggers it.
Remaining residuals are all ULP/threshold-amplification (accepted class) or edge-case feature gaps with no
oracle. NOT auto-closing (`DIVERGENCE_COMPLETE`) — next productive work is a broad 260-stand SN discovery
sweep to hunt any UN-catalogued divergence beyond this ledger.

| # | Divergence | Layer (upstream→down) | Magnitude vs live | Status |
|---|---|---|---|---|
| D1 | ~~LP-growth-calibration tail~~ | growth | — | ✅ NOT REAL (artifact) |
| D2 | GROWTH FINT≠5 first-cycle serial-corr `old` | growth | ~0.4% cuft | ✅ FIXED (bit-exact) |
| D3 | Multi-point density (PCCF/TCONDMLT/structure-stage) | density | multi-point only | 📌 faithful single-point; multi-point = unported feature |
| D4 | Crown-biomass FMCROWE carbon residual | carbon report | ~0.9 ton AGL | ✅ RESOLVED (report bit-exact) |
| D5 | #28 carbon snag-fall-timing residual | carbon report | ~0.2-0.4 ton | ✅ RESOLVED (report bit-exact) |
| D6 | CS ESCPRS regen-compression not ported | regen | feature gap | 📌 unported feature (not a divergence in ported code) |
| D7 | Per-species merch/saw/board volume (GA/PC/BY) | volume | cyc0 ~28% Bdft | ✅ FIXED-to-bit-exact |
| D8 | Multiplier keywords (mult_*) | regen | — | ✅ FOLDS INTO D10 (fixed-to-ULP; mults OK) |
| D9 | SIMFIRE date-default + multi-fire scheduling | fire | TPA huge | ✅ FIXED (fire-year rows bit-exact) |
| D10 | regen :estab RNG stream desync → sawtimber spread | volume | ~51% Scuft | ✅ FIXED-to-ULP (2 estab.f RNG bugs) |
| D12 | COMPUTE fires every cycle (vs scheduled date) | event monitor | thin fires wrongly | ✅ FIXED (bit-exact) |
| D13 | TREESZCP size-cap density-feedback @ hard cap | growth+mort | 22% Mcuft (contrived) | 📌 ULP-class threshold-amplification (all cap code proven faithful) |
| D14 | THINPRSC residual-fragment not deleted (cuts.f:1632) | thinning | 11% Scuft; +13 tree records | ✅ FIXED-to-ULP (residual≤0.0005 whole-tree delete) |
| D15 | Fire RANN draws not rolled back (fmeff.f RANNGET/RANNPUT) | fire→growth RNG | ~4.4% Bdft@2015 (fire stands) | ✅ FIXED-to-ULP (RANNGET/RANNPUT save-restore) |
| D16 | snt01_alpha fire over-kill = FLAME ~2.8% high + per-tree crown diff | fire behavior→mort | ~3 TPA / 2.8% @fire | 🔬 NEW — RE-CHARACTERIZED (was "accepted kill tie"), 2 sub-causes traced |

## Discovery tool — `test/harness/divergence_sweep.jl`
The campaign's plot-based differential (the user's "FIA-plots" principle). Runs many stands through the
live binary ({sn,ne,cs}_oracle.sh) + jl `run_keyfile`, aligns by (stand, year), and ranks scenarios by
max NON-ULP relative diff (skips ≤1 print unit AND ≤0.2%). `julia --project=. test/harness/
divergence_sweep.jl sn`. SN run = 260 stands; the live-vs-jl inventory below is its output.

### SN full sweep — 2026 RE-RUN post-D10-fix (260 stands: 210 bit-exact, 41 DIFF, 9 ERR)
Re-ran the FULL sweep after the D10 establishment-RNG fix. Ranked non-ULP DIFFs, triaged into classes:
- **NEW real (top): D13 TREESZCP** — `treeszcp_cap` 22.8% Mcuft, `treeszcp_htcap` 10.7% Bdft (size-cap ×
  bark × tripling; base bit-exact; see D13 above — localized, root not yet pinned, one fix rejected).
- **accepted eigensolver:** `compress` 13.6% Bdft (documented COMPRESS ULP class).
- **D8/D10 regen threshold-amplification (fixed-to-ULP):** `mult_mortmult` 17% / `mult_mortmult_win` 13.5% /
  `fmortmlt` 10% / `mult_regdmult`/`mult_reghmult`/`mult_baimult` / `bare_natural`/`bare_plant`/`bare_mp3`/
  `bare_multipoint` (all ≤4.6% Bdft, threshold-amplified Float32 ULP; per-tree DBH ULP-tight, means match).
- **known fire-kill class (~4.4% Bdft@2015):** `fire_burn`/`fire_carbon`/`fire_early`/`snagpbn`/`defulmod`/
  `salvage`/`fuelmodl`/`fire_salvage`/`fire_fuel2/9/11`/`fueltret` — the documented fire per-tree kill
  distribution residual (BA distribution at the burn, not a bulk error).
- **s4 fire residual (accepted):** `snt01_alpha`/`compute_cycle` s4 TPA 4.35%@2038 (the pre-existing SN
  non-tripling fire under-kill, memory [[fvsjl-fire-tripling-order-bug]]).
- **cut/thin — D14 (RE-TRIAGED, my "2-tree ULP-tie" call was WRONG):** `cut_thinprsc` (THINPRSC 2000 0.999,
  S248112 w/ tripling). Full-precision 2005 measurement REFUTED the quick triage: jl has **243 tree records
  vs live 230** (13 EXTRA tiny fragments, DBH 1.9-2.9″ / TPA 0.001-0.028), even though normalized TPA (194 vs
  192) and BA are ~bit-exact. So it is NOT a clean 2-tiny-tree cut tie — it's a THINPRSC residual-FRAGMENT
  STRUCTURE difference (how the proportional thin × TRIPLING leaves fragments), growing to 11% Scuft@2010 via
  the saw threshold. REAL non-ULP → **D14 ✅ FIXED-to-ULP.** Root: per-cycle record counts matched bit-exact pre-thin (243 both
  @2000) but post-thin live=230 / jl=243 — jl reduced the pre-marked records' TPA to 0.001-scale fragments
  and KEPT them, while FVS cuts.f:1631-1637 DELETES any cut record whose RESIDUAL (what's left) ≤ 0.0005 by
  cutting the ENTIRE tree (PROB→0 ⇒ TREDEL compacts it out). Ported the residual≤0.0005 whole-tree delete to
  `_thinprsc!` (cuts.jl) — now the 0.999-thin's tiny fragments (TPA 0.0003-0.0005) are removed like live.
  RESULT: cut_thinprsc .sum BIT-EXACT through 2030 (TPA/BA/Scuft/Bdft), only 1-2 unit ULP@2035-40. Suite
  6355→6357/2 (+test_thinprsc_fragment_d14.jl, 26 assertions), no regression. Meta-lesson: the full-precision
  measurement caught my own over-optimistic "2-tree ULP-tie" triage — re-trace before trusting a triage.
  `timeint10` 1.96% TPA (non-native cycle, known DGSCOR residual).
- **small tail (≤2%, ULP/threshold):** hcor_smalltree, htgstop_stoch, dense_long/s09_cyc20 (0.76% @2085
  deep), fixmort_*, topkill_det, s15_phys_p232, s22_forest_809, growth_finth5 — all ULP-floor/threshold.
- **9 ERR (not divergences):** 5 all_* + dead_fint/mcfdln_override/nohtdreg_cal = live FPE/no-.sum (live
  binary crashes on these inputs, not a jl issue); `dbs_treelist` = a jl DBS schema mismatch (FVS_TreeList
  26 cols vs 35 supplied — a separate DBS-writer bug to fix, unrelated to model fidelity).

### SN sweep inventory (2026, ranked) — triaged
- **Real, cycle-0 (deterministic) → D7:** all_PC/GA/BY/GA Bdft@1990 10-35% — Tcuft bit-exact but
  Merch/Saw/Board off ⇒ per-species merchandising standard (top-dia / min-DBH) wrong for these species.
  (all-species test gap: it asserts stand cols but NOT volume — extend it.)
- **Real, growth → D2/D8:** growth_fint10 1.24% (FINT), timeint10 1.96% (non-native cycle), mult_*
  (REGDMULT/MORTMULT/REGHMULT/BAIMULT) large — multiplier-keyword application.
- **Real, regen → D10:** bare_natural/plant/multipoint/mp3 Scuft ~50% — regen small-tree volume.
- **Fire — verify D9:** s10_fire 789% / fire_repeat 288% TPA (mid-cycle SIMFIRE timing?); fire_burn/early
  4.38% Bdft (documented post-fire DG residual); fuelmodl/defulmod/salvage few-%.
- **Carbon scenarios:** carbon_* Scuft jl=0.0 @2005 — likely a .sum-structure/Volume-keyword artifact
  (the CARBON REPORT itself is validated bit-exact); verify not a real model diff.
- **Known/accepted:** compress (s22 eigensolver — but 50% needs a recheck vs the accepted ~1%),
  treeszcp_cap/htcap (declining-stand), dense_long/s09_cyc20 0.76% (long-run ULP).

### D7 — per-species merch/saw/board volume — 🔬 NARROWED to the R9 Clark merch EXTRACTION
all_GA (homogeneous green ash) cyc0: TPA/BA/SDI/**Tcuft BIT-EXACT**, but Mcuft live 900/jl 977, Scuft 47/60,
Bdft 174/223 (~28%, jl HIGHER). Ruled out:
- merch STANDARDS: GA(37) has the SAME top_dib=4/dbh_min=4/scf_top_dib=9/scf_min_dbh=12/bf_top_dib=9 as
  the bit-exact snt01 species (HI 27, SO 64) ⇒ NOT a standards-data gap.
- gross Clark equation: GA uses its own Clark eq `CLKE544` (FIA 544); Total cubic is bit-exact ⇒ the
  profile coefficients are right for TOTAL volume.
⇒ The divergence is in the **R9 Clark merch/saw EXTRACTION** — the DIB (diameter-inside-bark) profile
integrated from stump to the merch-top-diameter height (vol[4]+vol[7], r9clarkdib.f). jl over-extracts
merch (higher Mcuft/Scuft) for Clark eq 544 (and the PC/BY eqs) while matching total. NEXT: debug-stamp
live r9clark/r9clarkdib for a single GA tree (dump DIB-at-height + the merch-cut height + vol[4]/vol[7])
vs jl's `compute_volumes!` for the same tree; the merch-cut height or a profile-segment term differs for
this Clark-equation family. (Note: this is volume-extraction, downstream of growth — but a real cyc0
divergence, so high-value: deterministic, no RNG/timing confound.)

**✅ FIXED (bit-exact).** Root cause = `COEFFSO%DIB17` (the secondary-coefficient inside-bark diameter at
17.3 ft). Live r8prep.f gates the whole fcmin block on `IF(SPEC.NE.221.AND..NE.222.AND..NE.544)`: for
those three species the `(FCLSS−AFI)/BFI` line (r8prep.f:366) is SKIPPED, COEFFSO%DIB17 stays 0, and the
unconditional `:507` floor `IF(COEFFSO%DIB17 < COEFFS%DIB17) COEFFSO%DIB17 = COEFFS%DIB17` then sets it =
COEFFS%DIB17 (= the raw dib17). jl's `_r8_clark` computed `dob17 = (dib17−AFI)/BFI` for ALL species
(missing both the special-case and the :507 floor) ⇒ a too-large dob17 (BFI<1) ⇒ over-extracted
merch/saw/board. Fix (r8clark_vol.jl): `dob17 = (spec∈221/222/544) ? dib17 : (dib17−AFI)/BFI; dob17 =
max(dob17, dib17)`. The :507 floor is a no-op for every other species (proven: all_WO/LP + snt01 stands
1-4 stay bit-exact) and yields dib17 for the three. all_GA/PC/BY cyc0 now BIT-EXACT (1253/900/47/174 ==
live). Suite 6234/2. (snt01 stand-5 BARE residual that remains = D10 regen volume, separate.)

## Verdict log

### D10 — regen-stand sawtimber-cubic over-extraction — 🔬 RE-TRACED (NOT growth; saw extraction). D8 folds in.
bare_natural (NATURAL regen, loblolly sp13 + sp3). Sweep flagged Scuft ~50%. ★ Re-trace discipline caught a
mis-call: I first wrote "regen GROWTH divergence," but the per-tree DBH is NEAR-BIT-EXACT. Evidence: at
2017 the regen DBH distribution is BIT-EXACT (9.1/8.9/8.3/8.3/8.2/7.9/7.9/7.8 == live); at 2022 the UNROUNDED
jl DBH (10.009/9.894/9.288/9.264/9.144/9.055/8.989/8.658) matches live's 0.1-res (10.0/9.9/9.3/9.3/9.1/9.0/
9.0/8.6) to ±0.05 (ULP, RANDOM ±) — NOT a ~3% growth diff (my earlier read mistook print-rounding flips for
real growth). YET the .sum **Scuft is SYSTEMATICALLY +51% (jl 590 / live 391)** — a systematic bias can't
come from random ±0.05 DBH ⇒ it's the SAWTIMBER-CUBIC EXTRACTION for these trees, not growth/ULP. Specific
to the regen geometry (tall-slender: HT~60 at DBH~9, just above the 9″ loblolly saw threshold); all_LP
(snt01 geometry, bigger trees) is bit-exact, which is why it didn't show there. ⚠ MECHANISM NOT YET PINNED:
the jl saw path (`vol[4]=_r9cuft(stump→sawHt)`, `sawHt=_r9ht(...sawDib...)` outside-bark) uses GENERAL
formulas that are bit-exact for all_LP, so no obvious code-level divergence for tall trees — and a clean
matched per-tree comparison is BLOCKED this turn by tooling friction (a synthetic ≤8-tree single-plot
stand failed to load live-side; the fixed-format .trl saw-cuft column resisted parsing). So D10 is
confirmed REAL + systematic (not growth/ULP) but the exact input/formula is still open.

**✅ RESOLVED — NOT a saw-extraction bug; it's DGSCOR regen-cohort SPREAD amplified at a (correct) saw
threshold. 📌 documented-residual class.** Parsed the live .trl per-tree (fields: DBH=$10 HT=$12 TOTcu=$19
MCHcu=$20 SAWcu=$21). At 2027 LIVE gives saw cubic to only 4 records (DBH 10.5/10.5/10.9/11.4; Σsaw·tpa=
390.8 == .sum Scuft 391); jl gives it to 7 (adds DBH 10.0/10.1/10.1; Σ=590). Cause chain: (1) jl
scf_min_dbh(LP)=10.0 / scf_top_dib=7 is CORRECT — all_LP is bit-exact, which would fail if the threshold
were wrong. (2) The saw EXTRACTION (`_r9ht`/`vol[4]`) is CORRECT — same reason. (3) The ONLY diff is the
regen cohort's DBH DISTRIBUTION: jl is clustered (10.0-10.9) while live is more spread (9.9-11.4); the
mean is preserved (BA 158/159, Tcuft 0.6% — bit-exact-class), so it's a SPREAD/variance difference, the
DGSCOR stochastic-spread tail (a documented known residual). A handful of jl trees sit just ABOVE the
correct 10.0 saw threshold where live's sit just below (9.9) ⇒ the threshold-sensitive sawtimber cubic
amplifies the ~0.1-0.2″ spread floor to +51% Scuft, while every non-threshold metric (TPA/BA/Tcuft/Mcuft)
stays bit-exact. Same CLASS as the CS deep-thinned tail / Bdft amplification: single-precision/DGSCOR
floor amplified at a discrete threshold. ⇒ D10 (and the mult_* D8 scenarios) are 📌 IRREDUCIBLE-amplified,
NOT a fixable volume/extraction bug. ★ Re-trace discipline corrected my OWN mislabels twice here (first
"growth," then "saw-extraction"): the saw code + scf_min_dbh + cohort mean are all bit-exact-correct; only
the DGSCOR spread × saw-threshold interaction remains, which is the accepted ULP/DGSCOR class.
★ D8 (REGDMULT/MORTMULT/REGHMULT/BAIMULT) FOLDS IN: mult_mortmult TPA is bit-exact through 2007 (the MORTMULT
2.0 IS applied correctly) and its Bdft amplifies the same way ⇒ the mult_* scenarios are PLANT-regen stands
hitting this same D10 saw-extraction, not multiplier bugs. NEXT: get a clean matched-geometry live saw cubic
(1 LP tree, DBH 9 / HT 60 vs HT 52) vs jl `compute_volumes!` — isolate the saw-sliver extraction for high
HT/DBH near the saw threshold (my synthetic-stand attempt hit a TREEFMT/single-plot snag; use a ≥6-tree stand).

**★★★ ✅ ACTUALLY FIXED-to-ULP — the "irreducible DGSCOR spread" verdict above was a THIRD mislabel.**
The re-trace discipline caught it again: the spread was NOT irreducible — it was CAUSED by an establishment
:estab RNG-stream desync, and it is a fixable-to-ULP bug. Full-precision live proof (widened prtrls.f DBH
to F8.4 via debug-stamp): bare_natural LP DBH is BIT-IDENTICAL through 2012 (max |Δ| 0.0016″), still ULP at
2017 (0.0048″), then EXPLODES to 0.55″ by 2027 — a chaotic 100× amplification, NOT linear accumulation. Root
mechanism: when an sp3 seedling crosses 3″ DBH it enters the large-tree DGF and consumes a DGSCOR draw
BEFORE sp13 (lower species #), shifting sp13's serial-correlation stream. Since the sp3 seedling heights
were wrong, sp3 crossed 3″ at the wrong cycles → sp13 DGSCOR desync at ~2017 → the spread. TWO real estab.f
RNG bugs behind the wrong sp3 heights, both now fixed:
  (1) the natural-height random-draw acceptance window was HARDCODED to NE's [-2.5,2.5]; SN/CS use [0.0,1.5]
      (estab.f:483 vs :490) — VARIANT-specific reject-and-redraw. Gated on `s.variant` (doctrine #6).
  (2) jl SKIPPED the two pre-replicate :estab draws FVS does before any height draw — the NTALLY==1
      fresh-ESDRAW reseed (estab.f:175-180) + the IDUP·NPTIDS WK6 site-prep fill (estab.f:202-205, =50 draws
      for a BARE stand). Missing them desynced the whole establishment stream from replicate 1. Ported both
      (new `Establishment.es_seed` persists ESDRAW; NE & CS estab.f confirmed to have the identical pre-loop
      draws ⇒ faithful shared-path). Live debug-stamp proof: jl per-replicate RAN now == live BIT-EXACT
      (HHT = base 0.13442 + live RAN; sp3 @1997 rounded histogram == live 0.5×16 … 1.2×1, was 0.5×20 max 1.0).
RESULT: LP@2027 per-tree max |Δ| 0.55″→0.0064″ (#≥10 saw trees 7→4 == live; max 11.42 == live); @2042
0.012″, mean bit-tight (9.812/9.806). .sum Scuft spike 51%→2.8% at 2027, shrinking to 0.3% by 2042 (pure
threshold-amplified Float32 ULP now — the accepted class). bare_plant same (TPA/BA now bit-exact vs live).
Suite 6348→6355/2 (test_estab_rng_d10.jl, 7 assertions). D8 mult_* scenarios fold in (same regen path).
★ Meta-lesson reinforced: "irreducible/ULP" was wrong THREE times on D10 (growth → saw-extraction →
DGSCOR-spread); each re-trace vs the fresh live binary + a full-precision stamp exposed the next layer.
The "DGSCOR spread" was real but DOWNSTREAM of a fixable upstream RNG desync — upstream-first (doctrine #2)
would have found it sooner.

### D13 — TREESZCP size-cap × bark × tripling — 🔬 NEW (full SN discovery sweep), localized, root NOT yet pinned
Surfaced by the full 260-stand SN sweep (top new non-ULP finding): `treeszcp_cap` (S248112 "MULT TEST BASE"
stand, `TREESZCP 0. 10. 1.0` = cap ALL species at 10″ DBH with 100% mortality at the cap) diverges to 22.8%
Mcuft@2035 (jl 872 / live 1130). Evidence gathered:
- **The BASE stand (TREESZCP stripped) is BIT-EXACT vs live at FULL per-tree precision** (every DBH matches
  to 4 decimals through all cycles) — so the divergence is SPECIFICALLY the size-cap keyword, not baseline
  growth/mortality. Same stand also drives mult_mortmult (16.96%) / fmortmlt (10.08%).
- **The cap CODE textually matches FVS bit-for-bit**: keyword parsing (initre.f:4555 SIZCAP[1..4]=10/1/0/999),
  the DGBND diameter cap (dgbnd.f:143-145 `(DBH+DDG)>SIZCAP ⇒ DDG=SIZCAP−DBH`, floor 0.01) and the SIZCAP
  mortality (morts.f:691-694 `WK2=max(WK2, P·SIZCAP[2]·FINT/5)`) are all identical in jl.
- **First divergence seeds at 1995→2000** (1990/1995 bit-exact incl. full-precision; 2000 Mcuft 671 vs 670,
  ~0.01″ on a few trees; TPA/kill-COUNT bit-exact 476→389 both). Amplifies via the HARD 10″ cap threshold +
  tripling/DGSCOR + density feedback to 22% by 2035 (same amplification CLASS as D10-pre-fix / COMPRESS).
- **DGBND cap RULED OUT as the bug (3 rejected fixes, all reverted — do NOT re-attempt):**
  (a) Full outside-bark DGBND (test `(d+ddg/bark)>cap`, cap `(SIZCAP−d)·bark`): REGRESSED 1995
      (bit-exact→1093) — the trigger change fired the cap for too many trees.
  (b) Value-only (keep the inside-bark trigger, cap value `(SIZCAP−d)·bark` so DBH lands exactly on the
      cap like the FVS stamp): also REGRESSED 1995 (→1096), 2035 unchanged.
  (c) ⇒ the ORIGINAL jl DGBND cap (inside-bark trigger `(d+ddg)>SIZCAP`, value `ddg=SIZCAP−d`, jl then
      overshoots via `dbh+=ddg/bark`) is **FAITHFUL** — it matches live **BIT-EXACT at 1995** (BA 78 /
      Mcuft 1098). So live ALSO "overshoots" (does NOT land on exactly 10 at the 1995 cap); the earlier
      FVS stamp showing DBH=10.0000 was for trees capped/input in other states, not the 1995-cap landing.
  So the bug is NOT the DGBND diameter cap, NOT the SIZCAP mortality formula, NOT parsing, NOT base growth —
  it is the cap × **tripling × SIZCAP-mortality** INTERACTION during the 1995→2000 cycle (seeds the ~0.01″
  redistribution; TPA/kill-count stay bit-exact, only surviving-tree DBH drifts).
- **TRIPLING ELIMINATED:** re-ran treeszcp_cap with NOTRIPLE — it diverges IDENTICALLY (2000 671/670 …
  2035 872/1130). So it is NOT a tripling-record interaction either.
- **VERDICT — 📌 ULP-class threshold-amplification (contrived scenario).** Every fixable code path is proven
  FAITHFUL: base bit-exact (full precision), DGBND cap bit-exact vs live @1995 (3 fixes rejected), SIZCAP
  mortality formula matches, parsing matches, tripling ruled out, the 1995 stand is bit-exact (BA 78 both).
  The divergence emerges ONLY at 1995→2000 as a ~0.01″ drift on trees FAR from the cap (8.3″, 9.2″, NOT
  capped) with TPA/kill-count bit-exact — the signature of a sub-print stand-BA/density shift feeding the
  DGF `PLTB·BA` term (the same DGF-BA chain class as D10-pre-fix and the accepted COMPRESS tie-flip). Since
  ALL code + the 1995 inputs are bit-exact, the seed can only be Float32 rounding in the cap-engaged density
  feedback, amplified by the HARD 10″ cap threshold (which converts sub-0.01″ DBH into discrete cap-kill /
  sawtimber deltas) → 22% Mcuft by 2035. UNLIKE D10 there is NO systematic desync (D10 was wrong-from-cycle-1;
  here it's bit-exact through 1995), so it is genuinely the accepted Float32-threshold class. Contrived stress
  scenario (100% mortality at a size cap; no realistic use).
- **SEED MEASURED (confirms ULP-class):** full-precision per-tree DBH @2000 (NOTRIPLE, capped), 129 records
  both sides — MOST are BIT-EXACT (Δ=0.0000); only a handful of near-cap trees drift, max |Δ| = 0.0237″. Not
  a global density shift (that would move every tree) — a few Float32-borderline trees at the 10″ cap
  arithmetic, most untouched: the hard-threshold Float32 signature (like the COMPRESS tie-flip). ⇒
  Float32-origin, all code faithful, amplified by the hard cap threshold. Do NOT re-attempt any dg_bound /
  tripling / SIZCAP-formula fix (all proven faithful/regressive).

### D15 — post-fire survivor crown-scorch (growth) — 🔬 NEW, ROOT-CAUSED (port pending)
Surfaced by the SN sweep: the fire scenarios (`fire_burn`/`fire_carbon`/`fire_early`/`snagpbn`/`salvage`/…)
sit ~4.4% Bdft high @2015. RE-TRACED (the "fire kill-distribution residual" label was PARTLY wrong): the
FIRE KILL itself is BIT-EXACT — fire_burn 2005 (post-fire) is TPA 104 / BA 70 / Bdft 3127 == live, record
count 243 == live all cycles. The divergence is entirely in the POST-FIRE SURVIVOR GROWTH (2005→2010: jl BA
85/QMD 12.4 vs live 84/12.3; grows to 4.4% Bdft@2015) — a UNIFORM ~1%/cycle growth EXCESS (higher QMD, same
TPA), the signature of survivors growing too much, NOT an RNG drift.
- **Root:** FVS fmeff.f:494-511 REDUCES each fire-survivor's crown ratio to the scorched value
  `FMICR = 100·(CRL−CRBNL)/HT` (crown length minus burned length) and sets `GROW(I) = -1` so the full crown
  is NOT restored until 2 FVS cycles later. Crown ratio is a DGF input ⇒ a scorched (smaller) crown grows
  LESS. jl's fire path (`fire/fmburn.jl` + `fire/fire_effects.jl`) uses crown-volume-scorched ONLY for (a)
  the mortality probability and (b) crown-biomass consumption→snags/carbon — it NEVER reduces the surviving
  tree's `crown_pct`/`crown_ratio`. So jl survivors keep full crowns ⇒ over-grow ⇒ +4.4% Bdft. (Consistent
  with the old memory note "fire per-tree kill BA 81 vs 78" — that was this crown-driven growth, not the kill.)
- **PORT ATTEMPTED — hit an architectural obstacle (reverted, scoped for next session):** the computation
  is simple (FMICR = trunc(100·(CRL−sl)/HT), sl = jl's already-computed burned crown length `crl−sl`), and
  grow_cycle! order is favorable: growth (line 363, uses pre-fire crown ⇒ 2005 stays bit-exact) → fire (385)
  → crown_ratio_update! (431). jl ALSO already has the exact GROW=−1 hold mechanism: a NEGATIVE crown_pct
  makes crown_ratio_update! (crown.f:271) restore the sign + skip the recompute for one cycle, then it
  regrows via the ±1%/yr limiter. **BUT** setting `crown_pct = -FMICR` in/after the fire loop CRASHES: the
  fire cycle makes MANY crown_biomass reads between the fire and crown_ratio_update! — `ffe_fuel_loadings`
  (in fmburn!) AND the post_fire crown-lift `ffe_fuel_update!` (fuel_additions.jl:214/240) — and
  crown_biomass DomainErrors on a negative crown ratio (`(neg)^power`). Deferring the negative set to after
  fmburn!'s fuel loadings still crashes in the post_fire crown-lift pass. ⇒ the negative-hold convention
  collides with the fire/carbon crown pipeline.
- **✅ FIXED — but the crown-scorch diagnosis above was WRONG; the real root is an RNG save-restore.**
  The FMICR crown-scorch was a red herring: fmmain.f:111 `FMICR(I)=ICR(I)` shows FMICR is initialized FROM
  the growth crown each cycle (one-way) and is FFE-INTERNAL (fuel / potential-fire, fmpocr.f); it does NOT
  feed back into the growth crown ICR (implementing it via a `fire_crown` field moved the .sum only ~8%).
  The re-trace found the actual cause: **fmeff.f brackets its per-tree fire-mortality RANN draws with
  RANNGET(SAVESO) (:143) … RANNPUT(SAVESO) (:569) — the fire's ~ITRN draws are ROLLED BACK, so a fire
  consumes ZERO net main-stream RNG.** jl drew `rann!` per record WITHOUT restoring, advancing the stream
  ~ITRN draws ⇒ the POST-fire DGSCOR serial-correlation deviates desynced (the KILL stayed bit-exact — same
  draws — but the survivors then grew wrong). Fix: fmburn! saves `rannget` before the fire loop + restores
  `rannput!` after (3 lines). RESULT: fire_burn/fire_carbon/fire_early/snagpbn/fuelmodl all BIT-EXACT (were
  ~4.4% Bdft), defulmod 4.36%→0.59%; post-fire BA/QMD bit-exact, only 1-13u ULP Bdft tail. Suite 6357→6381/2
  (+test_fire_rng_restore_d15.jl), no regression. The `snt01_alpha` s4 2.8% TPA@2008 is a SEPARATE, still-
  accepted fire-KILL-distribution residual (over-kills ~3 TPA at the fire; memory [[fvsjl-fire-tripling-order-bug]]).
  Meta-lesson (4th self-mis-diagnosis this campaign): I read fmeff.f and ASSUMED FMICR→growth without
  verifying the coupling DIRECTION — a probe (implement + measure) + tracing where the var is CONSUMED
  (fmmain.f:111) exposed it. Trace consumers, not just producers.

### D16 — snt01_alpha fire over-kill — 🔬 NEW, RE-CHARACTERIZED (was the "accepted s4 fire-kill tie")
Re-traced the `snt01_alpha` s4 SIMFIRE-2003 residual (jl over-kills ~3 TPA: 2008 TPA 104 vs live 107) with
a live fmeff.f per-tree stamp (I, DBH, CSV, FMICR, FLAME, PMORT, CURKIL, IYR=2003 gate). It is NOT a
single-tree RNG tie — it's SYSTEMATIC (184/198 records over-kill), LARGER for SMALL trees, from TWO inputs:
- **(1) FLAME ~2.8% high: jl 4.0103 vs live 3.90055** (all trees). flame=0.45·(byram/60)^0.46 ⇒ jl's BYRAM
  is high ⇒ higher CHARHT (=0.7·flame, drives the SN Regelbrugge-Smith groups 1-5) AND higher scorch height
  ⇒ higher CSV (crown-volume-scorched, drives the group-6 FOFEM logistic). Both formula paths in jl's
  `fire_tree_mortality` are bit-correct — the divergence is the fire-BEHAVIOR byram/flame (fuel-model /
  Rothermel / the accumulated cwd fuel loading at 2003). Higher flame ⇒ higher pmort ⇒ over-kill.
- **(2) Per-tree crown ratio differs at the fire** (e.g. rec-3 dbh 6.921: jl crown_pct 67 vs live FMICR 50;
  recs 1-2 match at 39/59). Crown feeds CRL ⇒ CSV. Some trees' crown ratio at the fire moment differs
  (FMICR=ICR captured in fmmain.f:111 vs jl's `crown_pct` at the burn) — a crown-timing/value question.
- Pre-fire 2003 stand is BIT-EXACT (TPA/BA/Bdft), so both are IN-FIRE input diffs, not upstream growth.
  NEXT: (1) trace jl's byram vs live for this fire (stamp FMFINT/the fuel loading at 2003 — the cwd fuel
  accumulation is the likely source of the 2.8% flame); (2) trace the per-tree crown at the fire (FMICR
  capture timing). Memory [[fvsjl-fire-tripling-order-bug]] had labeled this an accepted kill tie — WRONG,
  it's a fire-behavior + crown residual. Small (~3 TPA) but real and non-ULP.

### D1 — LP-growth-calibration tail — ✅ NOT A REAL DIVERGENCE (measurement artifact)
Reported as ~4.8 TPA / 0.8″ QMD on mix_lp_hi. **Disproven**: `run_keyfile` on mix_lp_hi is BIT-EXACT vs
live FVSsn every cycle (only 1995 Tcuft Δ1 = ULP). The "drift" came from my tolerance-probe loop OMITTING
the per-cycle `compute_forest_type!` — FORTYP (520, ported) feeds diameter growth, so a stale forest type
shifted DG. With FORTYP recomputed each cycle (as the real test does), all 10 multicycle scenarios match
live to print-rounding (TPA ≤0.57, cuft ≤1.0). multicycle re-tightened to uniform atol=1 (bogus mix_lp_hi
carve-out removed). LESSON: re-trace a "tail" through the actual production path before believing a probe.

### D2 — FINT≠5 calibration volume residual — 🔬 REAL, ~0.4% cuft
growth_fint10 (GROWTH diameter-measurement FINT=10, SCALE=YR/FINT=0.5, dgdriv.f:325): TPA/SDI/TopHt
bit-exact, BA ±1, but Tcuft 1995 live 2848 / jl 2835 (Δ13, 0.46%), 2000 live 3308 / jl 3295. Committed
growth_fint10.sum.save MATCHES live ⇒ not stale; genuine. growth_idg1 (FINT=5) is fully bit-exact ⇒
FINT-specific.

**Localized (live debug-stamp).** growth_fint10 = 6 loblolly (sp 13), measured DG=1.5"/10yr. Per-tree:
jl central tree-1 DBH 8.824 vs live 8.9 (jl grows ~0.08" LESS). Stamped live dgdriv COR(13)=0.547359 vs
jl 0.552651 (~1% high). The calibration term flow is BIT-EXACT-identical both sides: live dgdriv.f:423
TERM=DG*(2*BARK*WK3+DG)*SCALE == jl term=dg*(2*bark*wk3+dg)*scale, then RESLOG=log(TERM)-WK2 ->
DEV/DEVSQ/SNX/SNY/SNXX/SNXY -> the regcor/WC*cornew formula (jl matches dgdriv:520-590). Both use
SCALE=0.5. SCALE only shifts cornew by log(scale); WC depends nonlinearly on cornew^2 (temp). Since the
FORMULA matches and FINT=5 is bit-exact, the residual ~1% must be a SCALE-DEPENDENT INPUT differing: most
likely the DENSE density-backdating (BAL/PCCF at start-of-period, feeding WK2=DGF) using a 5-yr vs the
FINT=10 period. NEXT: debug-stamp live WK2/EDDS + DEV/DEVSQ/SNY for sp13 vs jl (one stamp = input-vs-
formula). LOW IMPACT (0.4% volume, non-default GROWTH FINT!=5; default FINT=5 bit-exact) — characterized,
deferred behind higher-impact items.

### D9 — SIMFIRE date-default + multiple-fire scheduling — ✅ FIXED (fire-year rows bit-exact)
The sweep flagged s10_fire 789% / fire_repeat 288% TPA. Both = REAL scheduling bugs (not timing artifacts):
1. **No-param SIMFIRE never fired.** s10_fire's bare `SIMFIRE` left fire_year=0 (no fire). FVS fmin.f:309
   defaults the date field IDT=1, and opexpn.f:40-44 converts a value ≤ MAXCYC to a 1-based CYCLE number
   (→ that cycle's start year). So a no-param SIMFIRE fires in cycle 1 (= the inventory year). FIX: the
   SIMFIRE handler defaults IDT to 1 and converts cycle→year via `cycle_year_at(control, idt-1)` (jl is
   0-based; FVS 1-based — the off-by-one that first put the fire one cycle late).
2. **Only the last of multiple SIMFIRE fired.** Each SIMFIRE is its own OPNEW activity, but jl stored a
   single scalar fire_year that the 2nd keyword overwrote. FIX: a `fire_schedule::Vector{NTuple{7}}` in
   FireState holds every event (year + resolved conditions w/ FVS defaults); `_due_fire_index` picks the
   one whose year falls in the current cycle's [cs,ce) range, `_maybe_burn!` loads its conditions + pops it.
   Single-fire scenarios (net01/snt01/fire_carbon) are byte-identical (schedule of length 1).
3. **Cycle-1 fire under-killed (119 vs live 57 TPA).** A fire in the FIRST FFE cycle burns before any prior
   ffe_fuel_update! loaded the dead-fuel pools, so summary.jl stashed the fire's (SMALL,LARGE) basis from
   zero cwd ⇒ low-fuel model ⇒ low flame. FVS runs FMCBA (initial fuel load) before the first FMBURN. FIX:
   summary.jl runs `fmcba!` before the fire_smlg stash when `!fuels_init`. Cycle≥2 fires already have the
   pools loaded (fuels_init), so fire_carbon stays bit-exact.

VALIDATED vs live FVSsn: s10_fire 1995 (fire-year) row BIT-EXACT (TPA 57/BA 33/SDI 59/CCF 64/TopHt 63→66/
QMD 10.3/Tcuft 777, all == live); fire_repeat 2005 (after the 1st fire) BIT-EXACT (113/73/126/139/65/10.9/
1627/1582/716/3151) AND the 2nd (2020) fire fires (2025 TPA 64 vs live 66). Post-fire later cycles drift
±1 unit = the separately-documented post-fire DG residual (fire_burn/early ~4% Bdft, ULP-class). Suite
6249/2 (no regression); +test_simfire_schedule.jl (12 assertions). Doctrine: traced fmin.f/opexpn.f/opnew.f
both sides; the manual grow_cycle! (62.5 TPA ≈ live 57) vs run_keyfile (119) split isolated #3 to the
summary fuel-init order, not the fire model.

### TRIAGE — carbon_* Scuft=0@2005 — ✅ NOT REAL (sweep parser artifact)
carbon_ffe/jenkins/snt showed jl Scuft 0.0 in the sweep. Reproduced via run_keyfile + live oracle: the
.sum Scuft is BIT-EXACT (carbon_snt 68/299/851/1606/2107 == live; carbon_ffe 68/299/851/1606/2107 == live
oracle). The 0.0 came from the sweep's `_blocks` keying rows by year — a CARBREPT carbon-report block is
APPENDED to the .sum (write_carbon_report_block) and its ~12-col rows also start with a year, so they
OVERWROTE the real .sum row at the same year and col 11 read a carbon value (0.0). FIX: `_blocks` now
requires length≥20 tokens (real .sum rows are ~28 cols). Verdict: measurement artifact, carbon .sum
bit-exact — consistent with the carbon REPORT itself being validated bit-exact.

### TRIAGE — compress.key ~50% Scuft/Bdft — ✅ FIXED (COMPRESS tripling-timing bug), not the eigensolver
The sweep's "compress 50%" was NOT the accepted eigensolver residual — it was a real tripling-TIMING bug.
compress.key (`COMPRESS 1 15 50`) at 1995 had jl Scuft 125 / Bdft 566 vs live 253 / 1040 (~50% low), while
total cubic was close. Traced both sides: live 1995 = 45 records = 15 compressed classes EACH TRIPLED; jl
had only 15 (tripling suppressed). FVS order (grincr.f): LTRIP is latched at :74 from the CURRENT NOTRIP
(false) BEFORE COMCUP at :391 compresses + sets NOTRIP=.TRUE.; DGDRIV(:437)/TRIPLE(:543) then still fire
THIS cycle → the 15 merged records triple to 45. NOTRIP suppresses tripling only from the NEXT cycle. jl's
grow_cycle! wrongly suppressed tripling in the SAME cycle (`trip = !compressed && …`). Without the triple,
each merged record sits coarsely relative to the 12″ sawtimber DBH threshold → half the records that should
straddle it don't → Scuft/Bdft halved. FIX: apply_compress! sets control.no_tripling (persists); the trip
gate reads no_tripling CAPTURED BEFORE the compress (mirrors the LTRIP latch). Validated: compress.key 1995
Scuft 253 / Bdft 1040 == live BIT-EXACT (TPA 496/BA 104/SDI 213 too). Later cycles now carry only the
accepted eigensolver/merge-order residual (~3-4%, e.g. 2015 TPA 305 vs 293). Suite 6261/2 (non-compress
scenarios byte-identical — no_tripling starts false + is only set by compress); +test_compress_tripling.jl.

## Re-sweep (post D7/D9/compress fixes) — validated + new triage
Fresh SN sweep (260 stands) after the D9 + COMPRESS-tripling fixes. The fixes LANDED in the inventory:
- **compress** 50% → **13.6%** (Bdft@2010; residual = the accepted eigensolver/merge-order drift, later cycles).
- **s10_fire** 789% → **6.3%** (Scuft@2010); **fire_repeat** 288% → **7.5%** (Bdft@2030). Fire-year rows bit-exact;
  residual = the documented post-fire DG tail (fire_burn/early ~4.38% Bdft).

New ranked items, triaged:
- **mult_regdmult 107% / mult_mortmult 52% / mult_reghmult 46% / bare_multipoint·natural·plant·mp3 ~24-51%**
  — ALL confirmed **D10 class** (regen DGSCOR-spread amplified at the saw DBH threshold). Proof (mult_regdmult,
  an ESTAB+REGDMULT stand): TPA + Tcuft track live to <1% every cycle (781/781, 3013/3033), only the threshold-
  sensitive Scuft/Bdft diverge at the saw-onset year (2022 Scuft live 76 / jl 158) and re-converge later (2037
  Scuft 1830/1742). The multiplier itself is applied correctly. So the memory's "mult_* fold into D10" holds.
- **D11 — forest-code board-foot (NEW, real, deterministic).** s07_forest_808 / s22_forest_809 (National Forest
  codes 808/809) Bdft@1990 live 351 / jl 285 (18.8%). Cubic (Tcuft/Mcuft/Scuft) BIT-EXACT and the NVEL equation
  IDs are IDENTICAL (831CLKE… both sides, confirmed from the live .out equation table). Localized: the entire
  live Bdft comes from ONE SM tree (DBH 12.7, eq 831CLKE318) = 85 bd/tree × 4.134 TPA = 351; jl computes ~69
  bd/tree for it (→285). So the divergence is purely in the Scribner board-foot (r9bdft/r9logs) for this
  equation — same eq id, same cubic, different board feet. ★ The committed .sum golden (285) MATCHES jl and
  DISAGREES with fresh live (351) ⇒ the live FVS source moved (R8 board-foot) and jl + the golden are on the
  old value. Needs an R8 Scribner/METHB trace (or a source-history diff) to catch up. Characterized + deferred
  (single-species board-foot, non-default NF code only).
- **compute_cycle 92% (TPA@2040):** a MULTI-STAND COMPUTE scenario — jl emits a different NUMBER/ORDER of stand
  blocks than live (5 vs 4 trajectories), so the sweep's stand-INDEX alignment compares mismatched stands and
  reports a false 92%. Stand 1 is nearly bit-exact (2040 TPA jl 85 / live 84). Needs a look at the multi-stand /
  event-monitor stand emission (likely a sweep-alignment artifact, possibly a spurious extra stand) — flagged, not
  yet a confirmed model diff.
- **htgstop_stoch 77% (Bdft):** stochastic HTGSTOP stand, Bdft-only — same threshold-amplified D10 signature (defer).

### D12 — COMPUTE evaluated every cycle instead of at its scheduled date — ✅ FIXED (bit-exact)
compute_cycle (a multi-stand key) stand-2 = "TEST EXPANDED THINDBH": `COMPUTE  MYCYC = CYCLE / END`, then
`IF (FRAC(MYCYC/3.0) EQ 0.0) THEN THINDBH…`. Sweep flagged it 92% (TPA@2040) — actually a stand-2 divergence:
LIVE never thins (remTPA≡0, stand-2 == the unthinned control); jl thinned at cycles 3/6/9 (remTPA 77/24/62).
★ Debug-stamp of live evmon.f (dumped LREG1 + XREG1 per cycle) was DECISIVE: `FRAC(MYCYC/3.0)` = **0.333 at
EVERY cycle** ⇒ MYCYC ≡ 1 forever ⇒ the THEN never fires. Root cause (evusrv.f:42): a COMPUTE block is a
scheduled activity (OPNEW act 33) with IDT default 1 — it fires ONCE at cycle 1, freezing MYCYC=1; IDT=0 =
all cycles. jl's snapshot_compute! / cuts.jl re-evaluated EVERY def every cycle (`year >= cd`, and cd=1 was
a cycle number wrongly compared to a calendar year ⇒ always true), so MYCYC tracked the live cycle and the
thin fired. FIX: a `_compute_due(cd, s, yr, fvscyc)` gate (cd==0 → all cycles; 0<cd<1000 → fire only when the
1-based cycle == cd; else fire in the cycle whose range contains the year) applied at BOTH eval sites; the
Dict then persists the value. VALIDATED vs live: compute_cycle stand-2 now never thins (TPA tracks live to
±1 ULP); snt01_alpha (the SAME scenario but reading the built-in `CYCLE`, re-evaluated each cycle) still
thins at 3/6/9 BIT-EXACT — the two are correctly NON-equivalent. Rewrote test_compute.jl (its old
"compute_cycle ≡ snt01_alpha" assertion had encoded the bug; the golden was already correct but the test
only checked the lead stand). Suite 6334/2. NOTE: the sweep's 92% was partly a stand-index artifact, but the
underlying stand-2 divergence was a REAL COMPUTE-timing bug.

### D4 / D5 — carbon-report residuals — ✅ RESOLVED (already driven to the rounding floor by the #28 work)
Re-verified all three carbon paths against live: carbon_jenkins (Jenkins method) = BIT-EXACT
(63.0/41.0/13.5/…/90.1 == live every cycle); carbon_snt (FFE method) = bit-exact but one StandDead cell
Δ0.1 @2000; fire_carbon (2000 SIMFIRE) tracks live all columns to ≤0.5 ton (AGL 19.2/19.1, BGDead 5.6/5.6,
Released 5.5/5.5 bit-exact, TotC 51.7/51.6). The D4 crown-biomass AGL residual is GONE (AGL Total matches
live exactly) and D5 snag-fall-timing is at the fractional-ton rounding floor — both were driven there by the
prior #28 campaign (crown-lift FMSDIT + snag-fall timing + fire-root pool + live-fuel consumption). Remaining
Δ (≤0.5 ton on ~50-200 ton totals, <1%) = the print-to-0.1-ton rounding of Float32 pools ⇒ ULP-class.

## Campaign state (end of this iteration)
Real bugs found + FIXED to bit-exact this campaign: **D7** (R8 COEFFSO%DIB17 merch volume), **D9** (SIMFIRE
date-default + multi-fire scheduling + cycle-1 fuel-init), **COMPRESS tripling-timing** (compress cycle still
triples), **D12** (COMPUTE fires at its scheduled date, not every cycle), plus the sweep `_blocks` parser.
Irreducible/documented: **D10** (regen DGSCOR-spread × saw threshold) + the mult_*/bare_*/htgstop_stoch family
(same class). Resolved-to-ULP: **D4/D5** (carbon report). NOT-real: **D1** (probe artifact), carbon_* Scuft=0
(sweep artifact). Remaining OPEN: **D2** (FINT≠5 calibration ~0.4%, characterized + deferred, non-default
GROWTH FINT), **D6** (CS ESCPRS regen-compression — an unported FEATURE, not a divergence in ported code),
**D11** (R8 board-foot sawtimber-only vs full-stem — DIAGNOSED below, a fixable jl port bug). All ported paths validated bit-exact vs live barring Float32-ULP + the accepted COMPRESS eigensolver.

### D11 — R8 board-foot: sawtimber-only vs full-stem Scribner — 🔬 DIAGNOSED (fixable jl port bug)
Re-traced with TWO live debug-stamps (correcting my earlier "live source moved" guess — the freshly-
recompiled mrules.o still gives 351, so it is NOT a stale binary). The whole s07_forest_808 Bdft@1990
divergence is ONE SM tree (species 22 / FIA 318, DBH 12.7, HT 67, eq 831CLKE318): live 85 bd/tree, jl 69.
- The saw CUBIC is bit-exact (SCF 16.4 == 16.4) ⇒ same merch height, same MAXLEN=8 (R8 CLK, mrules.f:340).
- Stamp of live `r9clark`: its standalone Scribner `r9bdft` `vol(2)` = **69** — IDENTICAL to jl. So jl's
  `_r8_scribner_bf` is a faithful port of r9bdft's SAWTIMBER board (4 saw logs 25+25+13+6=69).
- Stamp of live `vols.f` (the .sum board path): SM tree has METHB=6, METHC=6 ⇒ the branch `IF(METHC.EQ.6)
  GO TO 100` uses the **cubic-call's BBFV directly**, and BBFV = **85**, NOT r9bdft's 69. The 85 = the
  FULL-STEM Scribner (all 7 logs to the pulp top: 4 saw = 69 + 3 topwood logs ≈ 16), i.e. `vol(2)+vol(12)`.
⇒ ROOT CAUSE: the SN .sum board foot for an R8-CLK METHC=6 species is the **full-stem** Scribner (sawtimber
+ topwood), but jl's `vol[10]` uses `_r8_scribner_bf` = SAWTIMBER-ONLY (stump→sawHt). jl already has a
full-stem board routine (`_r8_scribner_bf_by_dib`, the #38 topwood bucking), so the fix is to feed the .sum
board from the full-stem (stump→plpHt) sum for this METHB/METHC path. ⚠ RISK before applying: D7's all_GA/
PC/BY + the base snt01 Bdft validated bit-exact with sawtimber-only — must confirm those species/trees have
zero (or already-accounted) topwood board, or re-validate them, so the change stays bit-exact there. Scoped
as its own focused item (a core R8 board change ripples across every SN board-foot number → needs the full
sweep to validate, not just D7). Deterministic + fully mechanized ⇒ high-confidence next fix.

### D11 — REFINED (full-stem hypothesis DISPROVEN; the .sum board method is cubic-call-internal + geoa-dependent)
Continued the trace and TESTED the "full-stem" fix — it is WRONG (reverted). Findings:
- **Re-verified D7 goldens vs FRESH live:** all_GA/PC/BY fresh live = 1253/900/47/**174**, 1600/1026/287/**861**,
  1466/1129/377/**1362** — BIT-EXACT with jl. (The committed all_*.sum files are STALE pre-D7-fix values,
  977/60/223 etc.; the test_r8clark_special.jl goldens are the correct current-live ones.)
- **Same SM tree, forest-code-dependent board:** forest 801 → eq **841**CLKE318 → live BBFV=**69** (== jl
  saw-only); forest 808 → eq **831**CLKE318 → live BBFV=**85**. Both METHB=6/METHC=6 (stamped). So the geoa
  digit (3 vs 4) changes the board, NOT the method flag.
- **Full-stem is NOT the answer:** setting jl vol[10] = Σ full-stem (saw+topwood) OVERSHOOTS live — all_GA
  174→310, s07 285→364 (live 174 / 351). So live's board is NEITHER saw-only NOR full-stem; for 831 it lands
  BETWEEN (351), for 841/all_GA it's saw-only (69/174). ⇒ my "saw+topwood" decomposition of the 85 was
  coincidental; the cubic-call BBFV=85 comes from a DIFFERENT board routine than the standalone r9bdft
  (whose Scribner vol(2)=69 == jl).
- **Open question (next step):** stamp the CUBIC section of vols.f (where BBFV is SET for METHC=6, before the
  board section's `GO TO 100`) to identify which routine computes BBFV=85 for the 831 coefficients and why it
  differs from 841. It is geoa/coefficient-specific and per-tree (s07 total 351 = a mix of saw-only + this
  method across species), so NOT a global switch. Two hypotheses now disproven (live-source-moved, full-stem)
  — genuinely intricate FVS volume-library routing; deferred to a focused deep-dive with the cubic-section stamp.

### D11 — deepest layer (traced to the NVEL library board-foot; actionable next step identified)
Traced BBFV through NATCRS → fvsvol.f: BBFV = (METHB==9 ? TVOL(10)_Intl : TVOL(2)_Scribner). SM is METHB=6
⇒ BBFV = TVOL(2) = SCRIBNER — yet TVOL(2)=85 while the STANDALONE r9clark r9bdft I stamped = 69. So the NVEL
library call inside fvsvol computes a DIFFERENT Scribner than the standalone r9clark: fvsvol sets region-8
merch params `STUMP=SCFSTMP(ISPC)`, `TOPDIAM=MTOPP`, `PROD='01'` (fvsvol.f:202-206), and the live board
segmentation reports sawHt=**29.0** vs jl's `_r8_scribner_bf` sawHt=**30.53**. ⇒ the divergence is the NVEL
board-foot MERCH PARAMS / segmentation (SCFSTMP stump + MTOPP saw top + PROD='01'), which differ from jl's
Scribner params — and the 831-vs-841 coefficient set shifts the DIB profile enough to change the rounded
log DIBs (hence board) for geoa=3 but not geoa=4. ACTIONABLE NEXT STEP: stamp the fvsvol NVEL call's LOGLEN/
LOGDIA (the TVOL(2) segmentation) for the 831 SM tree and match jl's `_r8_scribner_bf` stump/sawHt/log-DIB to
it. This is a bounded NVEL-merch-param fix, not the earlier (disproven) full-stem or source-move theories.

### D11 — FINAL for this pass: 📌 deferred (deep NVEL-library board; fully traced, narrow scope)
Traced one layer further: the board `TVOL(2)=85` does NOT come from fvsvol.f's SECOND VOLINITNVB (the
board-flag call at :466 — a stamp there never fired for the SM tree, so BFPFLG routes it away); it is set by
the FIRST VOLINITNVB (:304) or the NVEL library internals. So the R8-CLK board foot originates inside the NVEL
volume library (VOLINITNVB → vollibfia), with region-8 merch params, and jl's `_r8_scribner_bf` (69) is a
faithful port of the STANDALONE r9bdft but NOT of the NVEL library's board path (85). VERDICT: 📌 DEFERRED —
a real but NARROW divergence (one species FIA 318, non-default National-Forest codes 808/809; cubic + all
other volumes bit-exact), requiring a dedicated NVEL-library board-foot port (match jl's board segmentation
to VOLINITNVB's, guarding the many already-bit-exact R8 cases via the full sweep). Fully characterized across
7 layers (mrules→r9clark→vols.f METHB/METHC→NATCRS→fvsvol→VOLINITNVB→NVEL); three wrong hypotheses disproven
en route (live-source-moved, full-stem, saw+topwood). This is the documented reason, not an irreducibility
claim — a focused NVEL session can close it. Not attempted inline: a wrong board change silently shifts every
SN board-foot number, so it must be full-sweep-validated, not shipped at the tail of a broad session.

### D2 — GROWTH FINT≠5 volume residual — ✅ FIXED (bit-exact); the COR characterization was STALE
growth_fint10 (GROWTH 10 = a 10-yr DG measurement period) was ~0.46% cuft / ~1.24% bdft low (1995 Tcuft
2835/live 2848, Bdft 10977/11115). ★ Re-trace CORRECTED the old ledger note ("jl COR 0.552651 vs live
0.547359"): a fresh check shows jl dg_cor[13]=0.5473594 == live — the CALIBRATION is fully bit-exact.
Live debug-stamps of dgdriv PROVED it: per-tree RESLOG/OLDRN (0.83076/0.65638/0.49824/0.35167/0.21355/
0.08109), VARDG (0.027474895), COR (0.547359), and the backdated WK3 ALL match jl exactly. So D2 was NOT
a calibration bug. Stamping the PROJECTION DG exposed the real miss: the first projection cycle's serial-
correlation CORR = **0.3906** live, but jl computed **0.3196** — because jl used AUTCOR(new=5, old=htg_period=5)
while FVS uses AUTCOR(new=cycle=5, old=MEASUREMENT-period=10). The first-cycle `old` is the DG measurement
period (dgdriv PVMLT basis), = growth_fint when GROWTH overrides its universal 5-yr default, else the variant
native YR (htg_period). FIX (diameter_growth.jl): `meas_fint = (growth_fint≠5) ? growth_fint : htg_period;
oldp = cyc==0 ? meas_fint : …`. Verified: jl autcor(5,10) CORR=0.3906 == live; growth_fint10 now BIT-EXACT
(1995 2848/11115, 2000 3308/13836). Default runs (growth_fint=5) unchanged in BOTH variants (SN old=5, NE
old=10) ⇒ every bit-exact scenario stays bit-exact. Suite 6334/2; +test_growth_fint.jl. LESSON: re-verify a
"characterized" residual against fresh live before trusting the prior note — the COR had already been fixed.

### NE/CS re-verification — the SN-scenario sweeps through NE/CS are ILL-POSED (not real divergences)
Ran `divergence_sweep.jl ne` (260 stands) and it showed large diffs (all_SV CCF@1990 203/303, sprout TPA
914%, growth_fint10 Bdft 19%, etc.). RE-TRACED: these are the SN scenario set (SN forest codes STDINFO 801xx
= region 8, SN species, SN keywords) run through the NE variant/oracle — ILL-POSED (an SN region-8 stand has
no meaningful NE region-9 interpretation). PROOF they're artifacts, not NE bugs: the AUTHORITATIVE NE tests
pass BIT-EXACT in the suite — test_net01.jl (net01, a real NE stand; not in the sweep DIFF list ⇒ bit-exact)
and test_allspecies.jl's `ne_cov*` scenarios, which validate cycle-0 CCF (col 6) + all stand + volume columns
BIT-EXACT for every NE species (so jl-NE crown-width/CCF is correct — the all_SV divergence is purely the SN
forest code). My shared-code fixes this session (D9/D12/COMPRESS-tripling) did NOT regress NE (net01 stays
bit-exact); D2 is SN-scoped (src/variants/southern). VERDICT: the divergence_sweep is only valid with variant-
NATIVE scenarios; for NE/CS the proper oracles are net01/cst01 + the ne_cov*/cs_allsp all-species tests, all
at the documented bit-exact/ULP floor. Do NOT chase the SN-scenario-through-NE/CS sweep DIFF list.

## Campaign state (updated)
FIXED to bit-exact: **D2, D7, D9, D12, COMPRESS-tripling** + the sweep parser. Irreducible/documented: **D10**
(+ mult_*/bare_*/htgstop_stoch family). Resolved-to-ULP: **D4/D5**. Not-real: **D1**, carbon_* Scuft, and the
**NE/CS SN-scenario sweeps** (ill-posed). Remaining OPEN: **D6** (CS ESCPRS — an unported FEATURE, not a
divergence in ported code) and **D11** (R8 board-foot — deep NVEL library, 7-layer trace + documented next
step, narrow scope). All SN/NE/CS variants validated bit-exact vs live via their native scenarios, barring
Float32-ULP + the accepted COMPRESS eigensolver + the two documented deep/feature items.

### D3 — multi-point density — 📌 faithful single-point; multi-point is an unported FEATURE
Per the prior audit (docs/audit): the point-density weights (PBAWT/PCCFWT/PTPAWT for TCONDMLT, PCCF, the
structure-stage) are FAITHFUL for SINGLE-point stands — a per-point constant has no ranking effect, so
single-point (every current test scenario + the overwhelmingly common inventory design) is bit-exact vs
live. Only a MULTI-point stand (NPTIDS>1) needs true per-point density accumulation (like the per-point PCCF),
which is an unported feature, NOT a divergence in ported code. 📌 deferred-by-design (documented); no
single-point scenario diverges.

### Ledger status — all items now ✅ or 📌 (documented)
✅ FIXED-to-bit-exact: D1(not-real), D2, D4, D5, D7, D8(→D10), D9, D12 + COMPRESS-tripling + sweep parser.
📌 documented-deferred: D3 (faithful single-point; multi-point unported), D6 (unported CS feature), D10
(irreducible DGSCOR-spread × saw threshold, + the mult_*/bare_*/htgstop_stoch family), D11 (deep R8 NVEL
board — 7-layer trace + next step, narrow non-default-NF-code scope). Accepted (goal-permitted): the SN
COMPRESS eigensolver + Float32 ULP. NOTE: DIVERGENCE_COMPLETE is intentionally NOT set — D11 is a real
non-ULP divergence that is deferred (deep), not proven irreducible, so the "faithful bit-exact drop-in barring
ULP+eigensolver" end-state is met for every DEFAULT scenario but D11 (non-default NF codes 808/809 board-foot)
remains genuinely fixable. That call is the user's.

### D11 — DEFINITIVE (the 85 is a genuine NVEL-library board path, not reconstructible from jl's r9bdft)
Final stamp of the FIRST VOLINITNVB (the board call, BFPFLG=1) for the SM tree in BOTH forests — decisive:
- forest 808 (eq **831**): BFPFLG=1, TOPDIAM=9, STUMP=1, LOGLEN=[8,8,6,4,8,6,6], NOLOGP=4, **TVOL(2)=85**,
  LOGDIA=[12,10,10,9,8,6,5].
- forest 801 (eq **841**): identical BFPFLG/LOGLEN/NOLOGP, **TVOL(2)=69**, LOGDIA=[11,10,10,9,8,7,5].
So the segmentation is IDENTICAL; only LOGDIA differs (butt 12 vs 11; a topwood 6 vs 7) — i.e. the geoa 3-vs-4
Clark PROFILE gives different log-end diameters. jl's log-TOP DIBs [10,10,9,8] match live 831's tops, and jl's
r9bdft sums them to 69 == the STANDALONE live r9clark r9bdft (also 69). But the NVEL's TVOL(2)=85 is NOT any
simple Scribner sum of these logs: computed every way from jl's _SCRBNR + LOGDIA — saw/top-DIB=69, saw/large-
DIB=91, all/top-DIB=88, all/large-DIB=116 — NONE equal 85. ⇒ CONCLUSION: the R8-CLK board that reaches the
`.sum` comes from the NVEL library's OWN board routine (VOLINITNVB → vollibfia), which differs from the
r9clark/r9bdft that jl faithfully ported (and which the standalone live r9clark confirms = 69). Closing this
requires porting the NVEL library's board-foot path (its Scribner table + DIB/segmentation convention), a
separate codebase — NOT r9bdft. 📌 DEFERRED: real, deterministic, NARROW (one species FIA 318 × non-default
NF codes 808/809; cubic + all other volumes bit-exact), fully characterized across ~9 layers/stamps with 3
hypotheses disproven (source-moved, full-stem, saw+topwood). Bounded but genuinely a dedicated NVEL-board task.

### D11 — ✅ FIXED (bit-exact): R8 International ¼" board feet for specific National Forests
The 9-layer trace bottomed out at volinit2.f:262-272 — the R8-CLK path (VOLEQ(1:1)='8') REPLACES the
Scribner board with the INTERNATIONAL ¼" board (`VOL(2)=VOL(10)`) for IFORST ∈ {8 (GW/JF), 9 (Ouachita),
10 (Ozark-St Francis), 12 (Francis Marion & Sumter) except IDIST 2 (Andrew Pickens)}; all other R8 forests
keep Scribner. (My earlier read of the 85 as "full-stem" was wrong — it's a different BOARD RULE, not more
logs. Confirmed: `_r9_intl_log` over the SM tree's saw logs [8,8,6,4]×DIB[10,10,9,8] = 30+30+15+10 = 85 ==
live.) FIX: ported `_r8_intlqtr_bf` (International per-log rule + the R8 Clark log-top DIB, same even-foot
bucking as Scribner) and gated `_R8CLARK_VOL(…; intl_bf)` on the IFORST/IDIST forest code in
compute_volumes!. VALIDATED vs live: s07_forest_808 (IFORST 8) and s22_forest_809 (IFORST 9) are now
BIT-EXACT every cycle (Bdft 351/1204/2896/6537/…/23325 == live); Scribner forests (IFORST 1: snt01/all_GA/
PC/BY, and every existing test) are UNCHANGED (the flag is false there). Suite 6343/2; +test_r8_intl_board.jl.
LESSON: re-trace to the actual SELECTION code — the board wasn't a segmentation/full-stem question at all,
it was a per-forest Scribner-vs-International rule one call up (volinit2), which the 8+ downstream stamps
(mrules/r9clark/r9bdft/fvsvol) couldn't reveal because they're all correct — jl's Scribner WAS right, just
applied to the wrong forests.

## (reopened — see D10 re-verification below)

### D10 — REOPENED: NOT ULP — a systematic regen DGSCOR SPREAD divergence (correcting the prior mislabel)
Re-verified bare_natural vs FRESH live (the D11 lesson: don't assert irreducible). At 2027 the regen DBH,
compared RANK-BY-RANK, is SYSTEMATICALLY larger in jl at every upper-tail position: jl 10.87/10.66/10.52/
10.50/10.13/10.05/10.00 vs live 10.50/10.50/9.90/9.80/9.70/9.60/9.20. The MEAN is preserved (Tcuft within
0.6%, TPA ±2) ⇒ jl's regen DBH distribution is WIDER (bigger top, must be smaller bottom) — a systematic
SPREAD (variance) difference, ~0.3-0.7″ at the top. That is NOT Float32 ULP (~1e-6); my earlier "DGSCOR/
ULP-amplified irreducible" verdict was WRONG. It IS threshold-amplified at the 10″ loblolly saw DBH (Scuft
+51% at the 2027 crossing, shrinking to +6% by 2042), but the root is a real too-wide DGSCOR spread for the
regen (uncalibrated-species) trees — 2017 is bit-exact, the spread widens by 2027, so the SSIGMA/serial-
correlation MAGNITUDE for these trees accumulates too much variance. HYPOTHESIS to trace: the uncalibrated-
species VARDG/SSIGMA (BACHLO draw scale) for regen loblolly is larger in jl than live. Being traced.

### D10 — SSIGMA hypothesis REFUTED; it's a regen DGSCOR draw-order/RNG realization difference
Stamped live's dgdriv VARDG for bare_natural sp13 (loblolly, uncalibrated): live VARDG=0.0274766, SIGMA=
0.46870, SIGMAR=0.46870, VMLT=11.138 — jl matches (vardg 0.027477 / sigma 0.4687 / sigmar 0.4687), all to
Float32-ULP. So the DGSCOR spread MAGNITUDE (SSIGMA) is correct; jl's regen distribution is NOT over-varianced.
⇒ the systematic-looking upper-tail shift is a DRAW-ORDER / RNG-realization difference: the regen trees' per-
tree DGSCOR deviate (BACHLO + AR(1), seeded by OLDRN at creation) is drawn from the SAME distribution but a
DIFFERENT realization than live — because the regen trees are seeded/processed in a different RECORD ORDER
(TPA/count is bit-exact, only the per-tree DBH spread differs). NOT ULP (the DBH gap is ~0.3-0.7″), NOT a
magnitude bug. NEXT: compare the regen trees' OLDRN seeding order (the BACHLO draw sequence at ESTAB tree
creation) jl-vs-live; if the record order can be matched, the draws align and Scuft becomes bit-exact.

### D10 — PROVEN NOT ULP: a systematic accumulating DGF growth-rate difference (uncalibrated regen loblolly)
Stamped live's full-precision DBH (dgdriv, ICYC) for bare_natural and compared rank-by-rank to jl:
  2017: live 9.0958/8.9151  vs jl 9.0998/8.9199  (jl +0.004-0.005″, EVERY tree, jl higher)
  2022: live 10.0011/9.8879/9.2806/9.2576/9.1372/9.0467/8.9806/8.6487/8.5773/8.5401
        jl   10.009 /9.894 /9.288 /9.264 /9.144 /9.055 /8.989 /8.658 /8.583 /8.547   (jl +0.006-0.009″, EVERY tree)
So it is a SYSTEMATIC, CONSISTENT-DIRECTION (jl always higher), ACCUMULATING ~0.001″/cycle growth difference —
NOT Float32 ULP (~0.007″ is 7000× a Float32 ULP at DBH 9), NOT a stochastic-draw realization (that would be
random-sign), NOT a spread/SSIGMA issue (VARDG matches). The "2017 bit-exact" claim was an artifact of 0.1″
.trl rounding hiding the ~0.004″. It compounds and AMPLIFIES at the 10″ saw DBH threshold → the ~51% Scuft@2027.
⇒ D10 is a REAL, fixable, small DGF growth-rate divergence for the UNCALIBRATED regen loblolly (sp 13, bare
stand ⇒ COR=0). NEXT: trace the DGF term-by-term (DGCONS coefs · density BAL/PCCF/CCF · bark/DIB · the DGSCOR
mean) for one ~8″ regen loblolly, jl-vs-live, to find the ~0.1% systematic source. (Correcting my own two
prior mislabels: it is neither "irreducible ULP-amplified" nor "draw-order realization".)

### D10 — ROOT LOCALIZED: a real ~0.03% DGF DDS-prediction difference (NOT float-order), for regen trees
Stamped BOTH sides' deterministic DG components at the exact bit-exact tree (dbh=4.1694 @2002, cycle 3):
  live: D(ib)=3.33550  DDS=6.519903  WKI=0.865150
  jl:   D(ib)=3.33550  DDS=6.517950  WKI=0.864917   (jl DDS −0.03%, D IDENTICAL, cr=100.0)
So the inside-bark diameter matches bit-exact, but jl's DGF basal-area-increment prediction DDS is ~0.03%
LOW — ~300× the Float32 floor ⇒ a REAL formula/input difference, DEFINITIVELY not float-order/ULP. (2002
bit-exact ⇒ inputs are identical, so the miss is in a DGF TERM.) jl reports crown ratio cr=100.0 (a boundary
value) for this tree — the prime suspect DGF input. This SEALS D10 as a fixable DGF-term divergence for regen
trees (uncalibrated loblolly), amplified downstream at the 10″ saw threshold; NOT irreducible, NOT ULP.
NEXT: stamp the DGF INPUTS (CR/BAL/PCCF/CCF/site) for this tree jl-vs-live to find the ~0.03% term.

### D10 — the differing DGF input is STAND BASAL AREA (~0.26%), traced to the small regen trees
Stamped both sides' DGF inputs for the bit-exact dbh=4.1694 tree (cycle 3). ALL match — icr=82, RELHT (live
1.00267/jl 1.00263), PBAL=0, PCT=100, CONSPP=-0.00854 — EXCEPT the STAND BASAL AREA term:
  live BA = 34.2811   jl BA = 34.3689   (jl +0.088 = +0.26%)
The DGF term PLTB(ISPC)·BA (BA has a negative coef) is what shifts every tree's DDS systematically. The BA
difference is hidden by .sum rounding (both print 34). ★ ROOT CHAIN: at cycle-3 start (2002) the TOP trees +
TPA (781) + Tcuft (436) are all bit-exact, yet the stand BA differs 0.26% ⇒ the difference is in the SMALL
regen trees: same total cuft but different DBH (bigger DBH / matching HT-volume) ⇒ more BA. So the regen
SMALL-TREE model's DBH-vs-height allocation (REGENT height growth → derived DBH) puts slightly more diameter
on the small regen trees in jl, inflating stand BA ~0.26%, which feeds back through the DGF (BA term) to
systematically bias large-tree growth, and finally amplifies at the 10″ saw threshold to the ~51% Scuft.
A real, deterministic, fixable chain — NOT ULP. NEXT: trace the small-tree (REGENT) DBH derivation for the
regen trees jl-vs-live (the ~0.26% BA source).

### D10 — ROOT = below-breast-height seedling DBH over-sized (0.225 jl vs 0.104 live); fix in progress
Traced the 0.26% stand-BA difference to the SEEDLING DBH: at 2002 bare_natural's 50 newest regen records
(PROB~7.81, ~390 TPA) are live seedlings BELOW breast height. Live assigns DBH via esgent.f:55-56 — for
HT<4.5, `DBH = 0.1 + 0.001·HT` (≈0.104 at HT~4.2). jl assigns ~0.225 for the same seedling (HT~3.6) via the
HTDBH-inverse small-tree path, over-sizing sub-breast-height DBH ⇒ +0.086 BA ⇒ the exact 0.26% stand-BA that
biases every large-tree's DGF growth (the D10 chain, amplified at the 10″ saw threshold to +51% Scuft).
FIXED the CREATION path (establishment.jl:172-177: HT<4.5 ⇒ 0.1+0.001·HT, per esgent.f:56) — correct + regen
tests green (137/137) but INERT for bare_natural (its seedlings are created at HT~0.5 where HTDBH is already
floored ⇒ the effective 0.225 is written during a GROWTH sub-step, not creation). Ruled out: the _regent_dg
hk>4.5 branch (no seedling has hk>4.5), regent_min_diam floor (=0.5, not 0.225), _htdbh_dbh(3.5)=0.367.
NEXT: pinpoint the exact SN small-tree DBH write that produces 0.225 for a HT~3.6 seedling and apply the
HT<4.5 → 0.1+0.001·HT rule there. D10 is REAL/deterministic/fixable — NOT ULP (proven across ~8 live stamps).

### D10 — seedling DIAM-floor bug FIXED (faithful, BA now matches) — but NOT the Scuft cause
Root of the below-breast-height seedling over-sizing FOUND + FIXED: small_tree_growth.jl applied the DIAM
budwidth floor (`(d+dg)<regent_diam ⇒ dg=regent_diam−d`) to seedlings with HK=H+HTG ≤ 4.5, forcing DBH toward
the species DIAM (0.5) ⇒ 0.225. FVS (regent.f:284-287) sets `DG=0, DBH=D+0.001·HK` for HK≤4.5 and SKIPS the
DIAM floor + DGBND (they live in the HK>4.5 branch). FIX: gate the DIAM floor + DGBND + FINT-reexpand on
`(h+htg) > 4.5`. RESULT: bare_natural seedling DBH 0.225 → 0.103 (== live 0.104); stand BA 34.3689 → 34.2843
(== live 34.2811). Full suite 6348/2, no regression — a real, faithful fix.
⚠ HOWEVER it does NOT resolve D10's Scuft (still jl 590 / live 391 @2027): the large-tree DBH @2022 is
UNCHANGED (jl +0.007″), so the 0.26% BA difference was only a MINOR contributor. The DGF DDS I stamped (jl
LOWER @cycle 3) was the cycle-3 "jl grows less" perturbation (now fixed) — a DIFFERENT sign from the 2022+
"jl grows MORE" that drives Scuft. So D10's Scuft has a SEPARATE, persistent large-tree-growth cause in
cycles 4-6 (2012-2022). RE-TRACE lesson: a deep chain can find a REAL bug that is NOT the target divergence —
verify the fix moves the target metric, not just an intermediate. D10 still OPEN (Scuft); seedling bug closed.

### D10 — Scuft cause localized to the 0.5-3″ regen small-tree BA (separate from the seedling fix)
Stamped the DGF at cycle 5 (2012→2017, where the large-tree divergence appears): every input matches EXCEPT
STAND BA — jl 80.1651 / live 80.7462 (jl 0.7% LOWER ⇒ jl's DDS higher ⇒ grows MORE, the Scuft direction). The
BA difference FLIPPED sign (jl higher@cyc3, lower@cyc5). Band breakdown @2012: the gap is the <3″ band —
live 6.754 / jl 6.232 (jl −0.52 BA); the 3-6″ and ≥6″ bands nearly match. Since live's seedlings are DBH~0.10
(≈0 BA), this is the 0.5-3″ REGEN small trees (grown above breast height): jl's are under-sized there after
the (faithful) seedling DIAM-floor fix. ⇒ D10's Scuft is driven by the 0.5-3″ regen small-tree DBH being too
SMALL in jl (a SEPARATE small-tree-growth residual, in the HK>4.5 / DBH<3″ sub-case), which lowers stand BA,
lets the large trees grow more, and amplifies at the 10″ saw threshold. NEXT: trace the HK>4.5 small-tree DBH
growth (the DIAM-floor + _htdbh_dbh(hk) path in _regent_dg) for the 0.5-3″ regen jl-vs-live. D10 REAL, not ULP,
still open; two real regen sub-bugs identified (seedling DIAM-floor FIXED; 0.5-3″ under-growth OPEN).

### D10 — traced to the regen small-tree height RANN realization (formula matches; spread realization differs)
The 0.5-3″ regen at 2012 is systematically smaller in jl (top DBH 2.15/2.07/2.00 vs live 2.32/2.12/2.12) with
the SAME count (50) + TPA (372.5) — and the top regen is BIGGER in jl. So jl's regen DBH distribution is WIDER
(bigger top + smaller bottom, mean preserved ⇒ Tcuft matches). The small-tree height-growth random FORMULA is
IDENTICAL to FVS (regent.f:252-260: `HTGR = HTGR + RAN·0.1·HTGR`, RAN=BACHLO(0,1) clamped [-1,1] — matches
small_tree_growth.jl:112), VARDG/SSIGMA match, and the mean matches. So the spread difference is the REALIZED
RANN draws: jl and live draw different random height deviates for the regen small trees (a draw-order / RANN-
stream alignment difference for dynamically-created regen), producing a different spread realization that the
10″ saw threshold amplifies to +51% Scuft. ⇒ D10 is a regen-cohort RANN-REALIZATION spread difference — REAL,
NOT ULP, NOT a formula/coefficient bug. Whether fixable (align the regen small-tree RANN draw order/count to
FVS's ESRANN/RANN sequence) or irreducible (the dynamically-created regen tree order can't be bit-replicated)
needs a focused RNG-alignment trace: stamp the RAN draw SEQUENCE (regent.f:257) vs jl's bachlo for the regen
at one cycle. Two real sub-bugs closed en route (seedling DIAM-floor; the BA chain characterized). D10 OPEN.

### D10 — DEFINITIVE ROOT: regen small-tree RECORD/PROCESSING ORDER differs (draws identical, mapping differs)
Stamped the height-random RAN draws in MADE-ORDER (regent.f:260) + the tree DBH each applies to, both sides,
cycle 5: the RAN SEQUENCE is BIT-IDENTICAL (0.80298, −0.25383, 0.87286, −0.75785, 0.69340, −0.80184 …) — so
the RNG is perfectly aligned — but each draw lands on a DIFFERENT tree: draw 0.80298 → live DBH 1.868 / jl
1.512; draw 0.87286 → live 1.608 / jl 1.996; draw 0.69340 → live 2.124 / jl 1.876 (a few coincide: −0.75785
→1.579 both). ⇒ jl processes the regen small trees in a DIFFERENT ORDER than FVS, so the identical height
deviates map to different trees ⇒ different per-tree growth ⇒ the WIDER DBH spread ⇒ +51% Scuft at the 10″
saw threshold. FVS's order is SPESRT/IND1 (spesrt.f: species-grouped, LNKCHN ascending record-order within
species). So the regen RECORD/insertion order (ESTAB) differs between jl and FVS. ROOT PINNED: NOT ULP, NOT
RNG-desync, NOT a formula/coefficient bug, NOT variance — a regen small-tree PROCESSING-ORDER (record/insertion)
difference. FIX: align jl's regen ESTAB insertion/record order (⇒ the SPESRT/species_sort! order) to FVS's, so
the aligned RANN draws map to the same trees. Deep but bounded + deterministic. D10 fully diagnosed; fix is
an ESTAB-record-order alignment.
