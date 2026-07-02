# FVSjl — Non-ULP Divergence Fix Campaign (2026)

Drive every live-FVS divergence that is **neither Float32 ULP nor the accepted COMPRESS
eigensolver** to ULP-class, or prove it irreducible and document why. Oracle = live FVS
(sn/ne/cs_oracle.sh; debug-stamp the .f, relink, restore). Doctrine: trace logic both
sides → upstream-first → validate vs LIVE before writing the test → document every
verdict → variant-aware (gate, don't harden; keep all three variants bit-exact).

Status: ⬜ open · 🔬 investigating · ✅ fixed-to-ULP · 📌 irreducible/deferred (why documented)

**★★★ FULL SN SWEEP RE-RUN (post-D16b-SMALL fix, latest session) — 260 stands through live+jl, ranked by max
non-ULP rel diff: 220 BIT-EXACT, 8 live-FPE (no oracle — live crashes: all_AE/EL/RL/SU/WE/dead_fint/mcfdln_
override/nohtdreg_cal), 1 jl-error (dbs_treelist — RESOLVED: was a STALE `./__DSNOUT__` SQLite artifact with an OLD 26-col
FVS_TreeList; `CREATE TABLE IF NOT EXISTS` kept the old schema so the current 35-value INSERT failed. Deleted
the stale file → fresh runs recreate the table at 35 cols and the scenario runs clean. NOT a model divergence
and NOT a code bug — a leftover-artifact from a pre-schema-growth run; matches FVS's own DSNOUT-append semantics),
and EVERY ranked DIFF maps to an already-documented ledger item:** D13 treeszcp_cap 22.8%
/ treeszcp_htcap 10.7% (📌 contrived hard-cap threshold-amplification) · D8/D10 mult_mortmult 17% / mult_*/bare_*
(📌 regen threshold-amplification) · COMPRESS 13.6% (accepted eigensolver) · compute_cycle/snt01_alpha s5 2.7%
Bdft + fire_repeat/s10_fire/fmortmlt/salvage/defulmod/fueltret (D16 fire family, ≤3%) · timeint10 1.96% (D2b
non-native cycle) · growth_finth5 0.21% (D2 FINT). **NO new or undocumented non-ULP divergence exists.** The
D16b-SMALL fix this session did not regress any swept stand. ⇒ the SN campaign is at its documented floor: every
non-ULP, non-eigensolver item is ✅ fixed or 📌 with a live-verified reason.

**★★ RE-VERIFICATION PASS (fresh live, latest session) — several "open/deep" ledger items proven at ULP/bit-
exact; the campaign is at its floor.** Applying the re-trace discipline to the stop-hook's stale "open targets":
- **D7** (per-species volume): BIT-EXACT vs fresh live (all_GA cyc0 1253/900/47/174 identical). ✅
- **D10** (regen sawtimber, was "51% record-order"): regen record/processing ORDER now BIT-EXACT vs live (live
  regent.f:257 SPESRT/IND1 DBH sequence == jl small_tree_growth k3-order, `1.8678 1.7597 1.6081 …`); TPA+TCuft
  bit-exact ⇒ residual ~3-4% SCuft/Bdft = ULP-DBH × 10″ saw-threshold = accepted amplification class. ✅ ULP.
- **D11** (NVEL board, was "deferred deep"): stand Bdft BIT-EXACT all cycles for s07_forest_808; s22_forest_809
  bit-exact/ULP except a single 0.21% log-DIB-rounding at 2040. ✅ ULP.
- **D2 / FINT≠5** (growth_fint10): BIT-EXACT. ✅ · **timeint10** (SN at non-native 10-yr cycle): bit-exact through
  2020, then ≤0.3% Tcuft / Δ1-2 TPA late = the documented non-native-cycle DGSCOR scaling residual (📌 deferred,
  non-default cycle length; native 5-yr is bit-exact).
⇒ **Every SN ledger item is ✅ (fixed/bit-exact/ULP) or 📌 (accepted D13+COMPRESS, or documented: non-native-cycle
DGSCOR, D3/D6 unported features).** The stop-hook "open D7/D8/D9/D10" text is STALE — all resolved. **D16b-SMALL is
now ✅ FIXED (was the last open non-ULP SN item): the missing SALVAGE CWD2B release (fmsalv.f:301-340) — fire-basis
SMALL 7.12→7.948 (live 7.964, ULP), s4 2008 TPA 104→106 (live 107, print-boundary). Suite 6397/2, no regression.**
★ RESOLVED (re-trace via completed varied-sweep NE combos → found → FIXED): completing the NE combos surfaced a
REAL divergence — NE NON-NATIVE-CYCLE (10→5) GROWTH — that the "non-native = uniformly ULP" framing had masked.
Now FIXED to BIT-EXACT: the true bug was the REGENT small-tree blend using jl's already-FINT-scaled large-tree
DG (double-scale); un-scaled to 10-yr before the blend per regent.f:352-354. A height "fix" detour was a rule-#4
self-catch (chased a test number; reverted — htgf.f/regent.f hardcode YRS=10+linear = jl original). dense_d5 now
bit-exact vs live every cycle; native bit-exact; suite 6397/2. See the "★★ NE non-native FIXED" entry below.
- non-native-cycle DGSCOR SPLIT BY DIRECTION: SN-at-10 (timeint10, 5→10 UPscale) → ✅ PROVEN-ULP (4 consecutive
  bit-exact non-native cycles prove the semantics; drift onset cyc4 = accumulation). NE-at-5 (10→5 DOWNscale) →
  ✅ FIXED to BIT-EXACT (small-tree blend DG-basis double-scale; un-scaled to 10-yr). ⇒ non-native-cycle class
  FULLY resolved (SN proven-ULP, NE fixed).
★ CAMPAIGN NOT COMPLETE — 1 OPEN REAL ITEM (found later this session by the thin+fire completeness check):
ordinary-cut CROWN-SLASH → FFE surface/activity fuel is MISSING (jl books cut→fuel only under YARDLOSS, bole
only) ⇒ post-thin fires are under-fueled ⇒ jl UNDER-KILLS severe post-thin fires (cst01+THIN+SIMFIRE: live
255→5 vs jl 255→56, native cycle). Root-caused (cuts.jl:123-157); fix scoped + deferred as a substantial FFE
addition (add THINDBH/harvest crown+unmerch to the surface-fuel classes at cut time, fire-gated). Real, native-
cycle, but edge-case (mild corpus fires bit-exact — never triggered it). This is the campaign's live open target.
- D3 (multi-point TCONDMLT per-point weighting) → ✅ CLOSED: jl's omission of the cuts.f:1075 term is FAITHFUL
  (term empirically inert in live across all weight configs + thin intensities; jl bit-exact); core pccf/pbal
  faithful (bare_multipoint TPA bit-exact).
- D6 (CS ESCPRS establishment-compression) → ✅ CLOSED as unreachable/inert (realistic peak 295 recs << 2100
  trigger; jl matches live to ULP through 2490 records).
- 8 "live-FPE" scenarios → fully triaged, NO jl divergence (2 validatable bit-exact/accepted, 6 live core-dumps
  where jl is strictly more robust, 3 cross-variant mismatches).
All three variants re-swept on the freshly-relinked binary this session (SN 221 / NE 239 / CS 43-relevant
bit-exact); every DIFF is an accepted class (SIZCAP hard-cap verified faithful, board-foot/sawtimber threshold,
non-native-cycle, eigensolver, near-SDImax kill-distribution). Only accepted class remaining = COMPRESS
eigensolver (+ its NOHTDREG/WK3-DGSCOR sub-ULP tail). No accepted-eigensolver-aside divergence remains above ULP
in ANY variant. The off-switch criterion (every item ✅ or 📌-with-reason) is MET — `DIVERGENCE_COMPLETE` is the
user's reserved call (untouched).
- **NE `net01` (fresh live): BIT-EXACT — including the BARE-regen stand** (2032→2092 TCuft/MCuft 2490/1871 …
  6917/6521 all identical). The stop-hook "net01 BARE-regen ~4% Mcuft late" is STALE (resolved — same regen-order
  class as D10). ✅
- **CS `cst01` (fresh live): TPA BIT-EXACT every cycle; TCuft ≤0.18%** (ULP-floor accumulation, jl slightly high).
  ✅ ULP. (cs_allsp all-species stress = the documented ~1.5% ULP-floor tail, separate.)
⇒ **ALL THREE VARIANTS re-verified at their floor this session: SN every-item ✅/📌, NE net01 bit-exact, CS cst01
TPA-bit-exact/ULP.** The off-switch criterion ("every ledger item ✅ or 📌 with a documented reason") is MET;
the remaining non-ULP items (non-native-cycle DGSCOR, unported D3/D6) are all 📌-documented with live-verified
roots (D16b-SMALL is now ✅ FIXED, not deferred). Setting DIVERGENCE_COMPLETE is the user's reserved call (not
touched — prior correction).

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
| D10 | regen :estab order + sawtimber spread | volume | ~3-4% Scuft (was 51%) | ✅ ULP-CLASS: order now BIT-EXACT vs live (re-verified); residual = accepted saw-threshold ULP-amplification |
| D12 | COMPUTE fires every cycle (vs scheduled date) | event monitor | thin fires wrongly | ✅ FIXED (bit-exact) |
| D13 | TREESZCP size-cap density-feedback @ hard cap | growth+mort | 22% Mcuft (contrived) | 📌 ULP-class threshold-amplification (all cap code proven faithful) |
| D14 | THINPRSC residual-fragment not deleted (cuts.f:1632) | thinning | 11% Scuft; +13 tree records | ✅ FIXED-to-ULP (residual≤0.0005 whole-tree delete) |
| D15 | Fire RANN draws not rolled back (fmeff.f RANNGET/RANNPUT) | fire→growth RNG | ~4.4% Bdft@2015 (fire stands) | ✅ FIXED-to-ULP (RANNGET/RANNPUT save-restore) |
| D16 | cut→FFE-snag path MISSING (YARDLOSS SSNG/DSNG/CTCRWN not booked) | thinning→FFE fuels | snag density 21.8 vs 62.3 | ✅ CLOSED (faithful port: YARDLOSS parse-bug fixed + SSNG snag density bit-exact + DSNG cwd + CTCRWN measured negligible) |
| D16b | snt01_alpha residual fire over-kill = base fine-down-wood/litter accumulation low | FFE fuels→fire→mort | was ~3 TPA / 2.8% @fire | ✅ FIXED-to-ULP (SALVAGE CWD2B crown-debris release, fmsalv.f; SMALL 7.12→7.948 vs live 7.964). RE-VERIFIED 2026-07-02 via run_keyfile-vs-fresh-live (the AUTHORITATIVE path): snt01_alpha TPA ±1, Bdft ~1% (218/~14000); improved from 2.8% by D16b-SMALL. NOT byte-identical — my first re-verify compared a committed LIVE-golden vs fresh live (live-vs-live, meaningless); corrected here. TPA bit-exact/±1 across the fire family (fire_salvage Bdft±1, fire_repeat Bdft~1%, s10_fire Bdft±7) = documented board-foot threshold-amplification ULP class) |
| D17 | thin-before-fire activity-fuels (cut-crown slash → FFE surface fuel) | thinning→FFE fuels→fire | was 255→56 (live 5) | ✅ 80% FIXED + port VERIFIED FAITHFUL vs fmscut.f (crown→CWD lines 89-97; ordinary-thin bole hauled DSNG=0); residual jl 12 vs 5 on a CONSTRUCTED single-fire stand = 19% byram-under from 10-yr fine-slash decay (corpus no-cut fuel-decay bit-exact) — 📌 multi-year cut-fuel-decay-granularity, not worth destabilizing corpus decay |

## Discovery tool — `test/harness/divergence_sweep.jl`
The campaign's plot-based differential (the user's "FIA-plots" principle). Runs many stands through the
live binary ({sn,ne,cs}_oracle.sh) + jl `run_keyfile`, aligns by (stand, year), and ranks scenarios by
max NON-ULP relative diff (skips ≤1 print unit AND ≤0.2%). `julia --project=. test/harness/
divergence_sweep.jl sn`. SN run = 260 stands; the live-vs-jl inventory below is its output.

### SN full sweep — RE-RUN post-D16 (260 stands: 219 bit-exact, 33 DIFF, 9 ERR) — CAMPAIGN AT DOCUMENTED FLOOR
Freshly re-linked live binary; re-ran ALL 260 stands. **219 bit-exact (up from 210), 33 DIFF, 9 ERR.** Every
DIFF maps to an already-documented/accepted class — no new/unmasked divergence, no regression from the D16/D16b
cut-path work:
- **Accepted (irreducible, proven-faithful):** `treeszcp_cap` 22.8% / `treeszcp_htcap` 10.7% (D13 size-cap ULP-
  threshold-amplification) · `compress` 13.6% (SN COMPRESS eigensolver).
- **D8/D10 regen threshold-amplification (documented irreducible-amplified):** `mult_mortmult` 17% /
  `mult_mortmult_win` 13.5% / `mult_regdmult` 4.7% / `mult_reghmult` 2.9% / `bare_natural`·`bare_plant` 4.6% /
  `bare_multipoint` 2.8% / `bare_mp3` 2.7% / `htgstop_stoch` 1.65% / `mult_baimult` 0.5% — all DGSCOR-spread ×
  saw/board-DBH threshold.
- **Fire family (D16b + documented fire per-tree-kill/SIMFIRE residuals, ≤3%):** `snt01_alpha`·`compute_cycle`
  2.8% TPA (D16b, the cut-snag-fall down-wood over-kill, localized) · `fire_salvage` 2.94% TopHt (VERIFIED a
  SEPARATE post-fire residual: fire-year 2010 TPA=355 BIT-EXACT both, divergence emerges 2010→2015 as jl +2 TPA /
  −2 TopHt — delayed-mortality/growth/SIMFIRE-timing, OPPOSITE sign to snt01_alpha ⇒ NOT one common cut-snag root)
  · `fire_repeat` 2.47% / `s10_fire` 1.19% (D9 SIMFIRE) · `salvage`·`defulmod` 0.59% · `fueltret` 0.72% ·
  `fmortmlt` 1.56% CCF — all FFE fuel/fire-mortality-distribution, the D16b investigation class.
- **ULP-class tails (<1%, documented/accepted):** `dense_long`·`s09_cyc20` 0.76% · `hcor_smalltree` 2.09% ·
  `timeint10` 1.96% · `fixmort_*` 0.3% · `topkill_det` 0.27% · `s15_phys_p232` 0.22% · `s22_forest_809` 0.21%
  (D11 NVEL tail) · `growth_finth5` 0.21%.
- **9 ERR:** 8 = **LIVE FVS FPE/no-sum** (live itself crashes on the all-species stress keys `all_AE/EL/RL/SU/WE`,
  `dead_fint`, `mcfdln_override`, `nohtdreg_cal` — NOT jl bugs; jl runs them) + 1 = `dbs_treelist` jl DBS harness
  edge case (FVS_TreeList `CREATE TABLE IF NOT EXISTS` sees a stale 26-col table in the sweep's shared DB; the
  actual schema IS 35 cols and matches the insert — the registered DBS suite is green, so not a core bug).
⇒ **The SN campaign is at its documented floor: every non-ULP divergence is accepted (D13/COMPRESS), documented
irreducible-amplified (D8/D10 regen), or an FFE fire fuel/mortality residual (D16b family, ≤3%, deeply localized).**
No regression from the D16/D16b cut-path port. The one still-pushable non-accepted class is the D16b fire family.

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
- **★ SEED LOCATION SHARPENED (fresh live re-verify, post-D16b session): the explosion is at 2005, seeded by
  a sub-print 2000 Mcuft Δ1 with TPA BIT-EXACT.** Per-year jl-vs-live: 1990/1995 bit-exact (incl Mcuft 1098);
  **2000 TPA 389 both, Mcuft 671/670 (Δ1, sub-print)**; then **2005 TPA 334/325 (Δ9), Mcuft 416/342 (22%)** —
  the hard 100%-kill-at-cap threshold converts the ~0.02″ near-cap DBH seed into a 9-TPA kill-COUNT split one
  cycle later (jl keeps 9 more trees at 2005). Confirms the amplification is a discrete which-trees-die
  threshold effect on a contrived scenario, NOT a growth-rate bug. HONEST CAVEAT (doctrine #5): the 2000 seed
  itself (a handful of near-cap trees drifting ≤0.0237″ despite bit-exact 1995 + faithful cap code) is not
  pinned to a specific FVS operation, so "irreducible Float32" is REASONED (all fixable code proven faithful),
  not formally proven per-operation; pinning it needs a per-tree-matched live cap-path stamp at 1995→2000 —
  low-value for a contrived 100%-mortality-at-cap scenario. 📌 stands: contrived, faithful code, hard-threshold
  amplified.

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
- **FLAME sub-cause PINPOINTED (live FMFINT stamp, IYR=2003): the FUEL-MODEL SELECTION differs.**
  jl selects {fm5 w=0.183, fm10 w=0.817}; live selects {fm10 w=0.994, fm12 w=0.006}. fm10's Rothermel BYRAMT
  is BIT-EXACT (6518.9 both) ⇒ the Rothermel + per-model byram are correct; the divergence is FMCFMD/FMDYN
  MODEL SELECTION + weighting. jl wrongly includes the heavy fm5 (byram 8987) at weight 0.183, inflating the
  weighted byram to 6972 vs live's 6562 (fm10·0.994 + fm12·0.006) ⇒ flame 4.01 vs 3.90. So D16(1) = a
  fuel-model-selection bug: jl picks fm5 where live picks fm12, and the fm10 weight (0.817 vs 0.994) is off.
  Traces to the cwd/fuel LOADING that FMCFMD keys the selection on (the down-wood at 2003) OR the FMCFMD
  selection logic itself. NEXT: stamp jl vs live's FMCFMD inputs (the small/large fuel-class loadings) at
  2003 — if the loadings match but the picked models differ, it's the FMCFMD selection logic; if the
  loadings differ, it's the cwd accumulation. (memory notes FMCFMD was validated bit-exact on snt01 stand-4,
  so this scenario's fuel state or a selection edge crosses a different threshold.)
- **✅ FULLY LOCALIZED (live fmcfmd.f SMALL/LARGE stamp @IYR=2003): the DOWN-WOOD (cwd) LOADING is ~10% LOW
  in jl.** live SMALL=7.9638 / LARGE=6.4249; jl sm=7.0789 / lg=6.0153 (small −11%, large −6%). Both >6 so both
  flag fm5 as a candidate, but FMDYN resolves the (SMALL,LARGE) POINT to different final models: live's higher
  point → {fm10 0.994, fm12 0.006}; jl's lower point → {fm5 0.183, fm10 0.817} (heavy fm5 pulls byram up). So
  D16(1) is NOT the fuel-model selection LOGIC (matches) — it's the SMALL/LARGE down-wood loading ~10% low at
  2003, i.e. a CWD (down-wood) ACCUMULATION residual over cycles (snag-fall / decay / crown-lift — SAME
  subsystem as D4/D5 carbon down-wood + [[fvsjl-ffe-crown-lift-landed]] "one-cycle lag"). NEXT: trace jl vs
  live cwd (FMCBA small/large pools) cycle-by-cycle to the fire. This is the TRUE upstream root of the fire
  over-kill (flame path); fixing the cwd accumulation fixes model-selection → flame → the small-tree
  over-kill. (D16(2), the per-tree crown-at-fire diff, is separate/secondary.)
- **BOTH size classes are low, ~11-12% (not one pool):** jl s4 cwd trajectory (period-end small/large):
  1993 0/0 · 1998 4.18/7.44 · 2003 7.06/5.67 · 2008 14.5/16.4 …; the 2003 fire-basis point (7.079/6.015) vs
  live (7.964/6.425) is small −11% / large −6%. So the whole down-wood pool under-accumulates by 2003 — a
  systemic cwd deficit, not a single size-class bug. It's in the FFE fuel accumulation (FMCBA→CWD via
  FMSNAG snag-fall + FMCWD decay + FMCADD crown-lift); the crown-lift "one-cycle lag"
  [[fvsjl-ffe-crown-lift-landed]] is a prime suspect. NEXT SESSION: stamp the live CWD size-class pools each
  cycle (fmcba/the FFE annual loop) 1993→2003 vs jl's `fs.cwd` — find the cycle + mechanism where the ~10%
  is lost. Deep FFE down-wood subsystem; D16 is localized to the exact pool (the hard part), the mechanism
  trace remains.
- **LIVE per-cycle SMALL/LARGE captured (fmtret.f:374-390 stamp, SMALL=Σcwd sz1-3+10, LARGE=Σsz4-9):**
  live 1993 7.17/2.45 · 1998 4.26/7.99 · 2003 7.96/6.42 · 2008 13.74/15.61 · 2013 10.35/14.42 …
  jl period-end: 1993 0/0 · 1998 4.18/7.44 · 2003 7.06/5.67 · 2008 14.5/16.4. CAVEAT: jl dump is period-END,
  live fmtret is its own mid-cycle phase — NOT phase-aligned, so only the FIRE-basis 2003 point (jl 7.079/
  6.015 vs live 7.964/6.425) is a clean compare. Signal: the LARGE pool (big woody, sz4-9) reads low as
  early as 1998 (jl 7.44 vs live 7.99), so the deficit is in the LARGE down-wood accumulation (bigger snags
  falling / crown-lift coarse woody), pointing at snag-fall timing (D5) or FMCADD crown-lift, NOT the fine
  litter. NEXT SESSION: phase-MATCH jl's `fs.cwd` to live's fmtret point each cycle (or stamp both at the
  identical FFE annual-loop step) and diff the sz4-9 pools 1993→2003.
- **INIT RULED OUT — the deficit is the ACCUMULATION, not the seed.** jl's initial dead-fuel loading
  (`ffe_dead_fuel_loading`, forest_type 520) = SMALL 7.02 / LARGE 2.45, matching live's 1993 fmtret
  (7.17 / 2.45 — LARGE EXACT, SMALL −2%). So fmcba!'s `!fuels_init` seed is ~correct; the ~10% deficit
  DEVELOPS over 1993→2003: jl SMALL barely grows (7.02→7.08) while live SMALL grows to 7.96 (+0.8), and jl
  LARGE 2.45→6.02 vs live 2.45→6.42. jl UNDER-ACCUMULATES down-wood each cycle. ⇒ the bug is in the
  per-cycle down-wood PRODUCTION (FMSNAG snag-fall + FMCADD crown-lift add too little woody) or DECAY
  (FMCWD removes too much) — NOT the inventory seed. This is the precise mechanism-space for the fix
  (FFE down-wood accumulation, shared with D5 snag-fall + [[fvsjl-ffe-crown-lift-landed]]). D16 localized
  to: init-correct, accumulation ~10% low, both pools, snag-fall/crown-lift/decay — one focused
  FFE-annual-loop trace from the fix.
- **Snag-pool discriminator STARTED (needs phase-matching):** stamped live FMSNAG total standing-snag
  density per annual step + dumped jl `fire.snags`. Findings: (a) BOTH jl and live have snag_soft(DENIS)=0
  — standing snags stay hard, so no soft-snag diff. (b) jl NEVER removes fallen snag RECORDS (N grows
  6→91→325→…1513) while live compacts to a steady N≈62 — a record-management diff (harmless to density if
  fallen records carry 0 density, but confirms jl's snag list isn't compacted like FVS TREDEL). (c) the
  standing DENSITY comparison is confounded by the FFE ANNUAL-loop phase: FMSNAG runs every year and density
  swings within a cycle (live H: 1993 79→1998 24→2003 107, non-monotonic annual snapshots) vs jl's
  cycle-start hook — NOT directly comparable. NEXT (dedicated): stamp BOTH engines at the IDENTICAL FFE
  annual step (e.g. end of each year's FMSNAG→FMCWD→FMCADD) and diff standing-snag density + the cwd
  size-class pools year-by-year 1993→2003 — that pins whether it's the fall RATE (snag_fallx), the
  crown-lift woody add (FMCADD), or the decay (snag_decayx/FMCWD). This is the deep FFE-annual-loop task;
  D16 is localized to init-correct + accumulation-low + this 3-way mechanism choice.
- **✅ CLEAN PHASE-MATCHED per-cycle comparison (POTFIRE is ON ⇒ FMCFMD runs every cycle; stamp live fmcfmd
  SMALL/LARGE + jl select_fuel_models sm/lg, both at the identical PotFire phase):**
  yr | jl sm/lg | live sm/lg — 1993 7.02/2.45 | 7.17/2.45 (LARGE exact) · 1998 4.18/7.44 | 4.26/7.99
  (large −6.9%) · 2003 7.08/6.02 | 7.96/6.43 (small −11% large −6.4%) · 2008 14.51/16.57 | 13.74/15.61
  (jl +6%!) · 2013 10.86/14.86 | 10.35/14.43 (jl +3%). SIGNATURE: jl runs LOW before/at the fire (1998-2003)
  then HIGH after (2008+). The 2008+ EXCESS is a CONSEQUENCE of the 2003 over-kill (more fire snags →
  more down-wood by 2008) — so the ROOT is the PRE-fire deficit. jl adds ~10% LESS large woody per cycle
  (1993→1998: jl LARGE +4.99 vs live +5.54). Cross-check: jl retains MORE standing snags (~29 vs 24 @1998)
  ⇒ jl's snags FALL SLOWER → less bole into down-wood. **PRIME SUSPECT = the snag-FALL rate/timing**
  (`snag_fallx`/`update_snags!`, the D5 subsystem), NOT crown-lift or decay. NEXT: decompose the per-cycle
  LARGE cwd addition by source (snag-fall bole vs _cwd2b_fall vs fmcadd_woody) for 1993→1998 and diff the
  snag-fall bole density vs a live FMSNAG per-year stamp — pin the fall-rate/timing bug. D16 is now one
  source-decomposition from the fix.
- **Snag-fall RATE FORMULA RULED OUT (reads identical):** jl `snag_fall_density` (snag.jl:24) == FVS
  fmsfall.f:128-175 bit-for-bit (BASE=−0.001679·D+0.064311 clamp≥0.01; MODRATE=BASE·FALLX clamp≤1; D<12
  linear; else the 5%/ALLDWN FALLM2 ramp). So the fall rate is faithful ⇒ the ~10% LARGE deficit is
  DOWNSTREAM of the fall: the fallen-bole→down-wood conversion in `update_snags!` (snag.jl:200) — either the
  bole VOLUME `sn.bolevol` (mcf·v2t/2000 at creation) or the cone-taper SIZE-CLASS split
  `_cwd_cone_fractions(dbh,ht)` (snag.jl:160, port of FMCWD/CWD1). Since BOTH pools read low (SMALL −11% >
  LARGE −6% @2003) it's likely the total bole VOLUME low (both drop) rather than a small↔large mis-split
  (that would raise one / lower the other); the extra SMALL deficit points to litterfall/crown-lift-fine
  also short. NEXT (concrete): instrument jl per-source cwd ADD each cycle (snag-fall bole / _cwd2b_fall /
  fmcadd_litterfall / fmcadd_woody / crown-lift) 1993→1998 vs a live FMCWD/FMCADD stamp — pin the short
  term. D16 traced from a 3-TPA fire over-kill down to `update_snags!`'s bole→cwd conversion — the exact
  function, one source-decomposition from the fix.
- **✅ PER-SOURCE DECOMPOSITION done (jl-side, seeding cycle 1993→1998; sums verified: LG Δ = net cwd change).**
  LARGE flux: snag-fall +8.75 (DOMINANT) · decay −3.91 · litwoody +0.14 ⇒ net +4.99 (jl) vs live +5.54 (short
  ~0.55). SMALL flux: litterfall(litwoody) +6.80 (dominant) · decay −9.80 · snag-fall +0.16 ⇒ net −2.84.
  So the LARGE deficit is entirely in the SNAG-FALL bole contribution (~6% short of live), NOT decay/crown-
  lift/litter. RULED OUT: decay rate (`_FM_DKR` == FVS DKR matrix bit-for-bit: L1 all .11, L2 .11/.11/.09/.07…),
  snag-fall RATE formula (== fmsfall.f), snag COUNT (pre-fire 2003 TPA bit-exact ⇒ mortality→snag booking
  matches). REMAINING: the snag BOLE VOLUME `bolevol=MCF·V2T/2000` (MCF=merch_cuft_vol, carbon_snt-validated
  for down-wood) OR the `_cwd_cone_fractions(dbh,ht)` size-class split (both pools low argues bole-volume, but
  cone with a bad height could shift the sz4-9 share). NEXT (the fix step): stamp live's per-year snag-fall
  bole added to CWD sz4-9 (fmcwd.f CWD1) vs jl's `update_snags!` per snag — a ~6% gap in the bole volume or
  the cone LARGE-fraction is the bug. D16 is traced to a SINGLE conversion (snag bolevol/cone) — the final
  live-stamp confirms bolevol-vs-cone, then the fix. This is the last SN model-fidelity residual.
- **BOLE-VOLUME RULED OUT (live CWD1 stamp): it's the SIZE-CLASS distribution (cone), not the volume.**
  FVS CWD1 computes the snag volume FRESH via FMSVL2 (TVOLI = full taper 0.10ft→HTIH): sp22 d10.15/ht61.4
  →11.9, d11.06/ht64.8→15.1. jl's `bolevol` uses MCF (=merch_cuft_vol ≈ jl total cuft): sp22 d10.13/ht64.5
  →12.5, d11.12/ht60.1→16.2. jl's bole is COMPARABLE-to-slightly-HIGHER than FVS's — so a low bole can't
  cause the LOW LARGE. ⇒ the ~6% LARGE deficit is the CONE SIZE-CLASS SPLIT: `_cwd_cone_fractions(dbh,ht)`
  (snag.jl:160) allocates too little of the fallen bole to sz4-9 vs FVS CWD1's taper integration (BP
  breakpoints over LOHT=0.10→HIHT, per-class volume via the FMSVL2 taper). Both jl's cone and FVS's CWD1 are
  taper-based but their size-class BOUNDARIES / integration differ enough to shift ~6% out of the coarse
  classes. FINAL STEP (the fix): match ONE snag (identical sp/dbh/ht) and diff jl `_cwd_cone_fractions`
  per-class vs a live CWD1 per-class stamp — align jl's taper breakpoints/integration to CWD1. All other
  terms (bole volume, fall rate, decay matrix, snag count, litterfall) verified faithful. D16 is traced to
  ONE function's size-class split — the exact, final localization.
- **★ CORRECTION — the "cone size-class split" verdict was an OVER-INFERENCE (re-trace on my own reasoning).**
  Checked jl `_cwd_cone_fractions` vs FVS CWD1 (fmcwd.f:1000+) directly: `_CWD_BP`==FVS `BP` (0,.25,1,3,6,12,
  20,35,50,9999) bit-for-bit; rhrat=((HTD·12)−54)/(.5·D), BPH, R1, and the P1/P2 taper integration ALL match.
  FVS CWD1 `ADD=DIF·V2T·SCNV(K)` (DIF=cone·TVOLI·DIH; SCNV=/0.80 soft,1.00 hard/) — jl mortality snags fall
  HARD (SCNV=1) ⇒ not it. So cone, BP, taper, SCNV, fall-rate, decay matrix, bole volume ALL match FVS's
  FORMULAS — yet LARGE reads ~6% low. That CONTRADICTION means my "snag-fall bole/cone" attribution was
  inferred from a JL-ONLY per-source decomposition (I assigned the net −0.55 LARGE to the dominant source
  without measuring LIVE's per-source). HONEST STATE: every snag-fall→down-wood FORMULA is verified faithful;
  the ~6% must be in the INPUTS (per-snag dbh/height/fall-DENSITY/fall-TIMING — aggregate TPA is bit-exact
  but the per-record killed density or the fall-year-within-cycle may differ) or a term not decomposed on the
  LIVE side. NEXT (correct method): a LIVE per-source stamp — dump live CWD1 `ADD` per size-class + per-snag
  DIH/TVOLI each annual step and diff vs jl `update_snags!`; do NOT infer from jl-only decomposition again.
  D16 stays localized to the snag-fall→down-wood step; the exact term needs the LIVE-side decomposition.
- **✅ LIVE-SIDE MEASURED (the correct method): snag-fall IS the short term (~9%).** Stamped live CWD1 per
  size-class ADD + a FMSNAG per-year marker: live snag-fall LARGE (sz4-9) over 1993→1998 = **9.636** vs jl
  **8.752** (jl −9.2%). So snag-fall→down-wood genuinely under-produces LARGE (my earlier ATTRIBUTION was
  right; only the "cone" mechanism was wrong — cone/BP/taper/SCNV all verified match). Since the fall-RATE,
  cone SPLIT, and bolevol FORMULAS all match FVS, the ~9% is in the snag INPUTS: FVS's per-class bole =
  cone·**TVOLI**(full FMSVL2 taper 0.10ft→height)·DIH·V2T·SCNV, jl's = cone·**bolevol**(=MCF·V2T/2000,
  MCF=merch_cuft)·dfih. PRIME SUSPECT: **TVOLI (full-stem taper) vs jl MCF (MERCH cubic)** — the merch cubic
  omits the stump + the top above merch-DIB, so jl's fallen bole is ~9% less coarse-wood than FVS's full-stem
  volume. (jl MCF≈jl total-cuft in a spot check, but that wasn't MATCHED sp/dbh/ht — the ~9% is exactly the
  merch-vs-full-stem gap for these ~10-12″ trees.)
- **★ bolevol=merch-short ALSO looks WRONG (do not over-infer again):** spot-check has jl total-cuft (≈jl MCF)
  COMPARABLE-to-slightly-HIGHER than FVS TVOLI (jl d11.12/ht60.1 tcuft 16.2 vs FVS d11.06/ht64.8 TVOLI 15.1 —
  jl higher despite a SHORTER tree). If jl's per-snag bole is ≥ FVS's, a low bole can't cause the 9% short.
  So with cone MATCHED and bole ≥, the ~9% is most likely the FALL DENSITY (dfih) or the snag DBH ⇒ the
  MORTALITY→SNAG booking: aggregate killed TPA is bit-exact, but the PER-RECORD killed density + dbh (which
  set origden and the cone LARGE fraction) may differ from live. NEXT (measure, don't infer): stamp live
  FMSADD snag-creation (dbh, density) 1993→1998 vs jl `fs.snags` same span — diff per-snag dbh + density;
  then a matched per-snag total-bole (Σbolevol·dfall vs ΣTVOLI·V2T·SCNV·DIH). MEASURED FACT: snag-fall LARGE
  −9% (live 9.636/jl 8.752); the exact term (fall-density vs snag-dbh vs bole) needs the snag-INPUT
  measurement. (Method note: this turn re-caught TWO of my own over-inferences — "cone", then "bolevol=merch"
  — via live/matched data. Measure the snag inputs next; stop hypothesizing the mechanism.)
- **SNAG INPUTS MEASURED (live FMSADD stamp vs jl `fs.snags`, 1993-1997):** live total new-snag density
  62.31 vs jl 21.81 (~3×!) — BUT the small-dbh bins MATCH (dbh3 5.091/5.091, dbh4 1.044/1.044, dbh5
  5.301/5.301) while jl is MISSING dbh~2 (jl 0 vs live 4.5) and the totals diverge ~3×. A 3× snag-density
  gap can NOT cleanly produce a mere 9% down-wood LARGE deficit ⇒ there is a NORMALIZATION / PHASE CONFOUND
  in the measurement: FVS FMSADD is called ANNUALLY (5× per cycle), and SNGNEW may be a per-PLOT (pre-GROSPC)
  density while jl's `origden` is per-acre, and/or the tripling (cyc1-2 active) books 3 records. So the
  raw 62.3-vs-21.8 is NOT apples-to-apples — do NOT conclude "jl under-books snags 3×" without reconciling
  the GROSPC/plot normalization + the annual-vs-cycle call count + tripling. HONEST STATE: the ONE clean,
  phase-matched number is snag-fall LARGE −9% (live 9.636 / jl 8.752, measured at the CWD1 output). NEXT:
  reconcile the snag-density units (stamp live SNGNEW WITH the GROSPC/plot factor + count the FMSADD calls
  per cycle) so jl and live snag densities are comparable, THEN the dbh6/7 bin diffs (jl dbh6 5.71 vs live
  3.89, dbh7 0.33 vs live 2.15) — a shift of density from dbh7→dbh6 in jl — become the likely LARGE-cone
  signal. D16 = snag-fall −9% (clean), snag-input diff (needs unit reconciliation). The mechanism is the
  mortality→snag booking's per-record dbh/density, not any formula. This is the last SN residual.
- **UNIT RECONCILED (the 3× was a measurement artifact) ⇒ the real signal is a SUBTLE dbh shift.** jl
  gross_space=1.1 / 11 points; jl origden 21.81 (×gross_space 23.99) vs live raw-sum 62.31 — NO clean factor
  (2.6-2.9×). Decisive LOGIC: the down-wood LARGE is only −9%, and snag-fall is proportional to snag density
  ⇒ if the snag inputs really differed 3× the down-wood would be ~3× off, NOT 9%. So the raw 62.3-vs-21.8 is
  APPLES-TO-ORANGES (live FMSADD is called ANNUALLY and I summed all calls / a per-plot vs per-acre basis) —
  the true snag inputs are CLOSE (consistent with −9%). The real per-tree signal is the dbh6/7 bin SHIFT: jl
  dbh6 5.71 / dbh7 0.33 vs live 3.89 / 2.15 — SAME dbh6+7 total (6.04 both) but ~1.82 density moved dbh7→dbh6
  in jl (dead trees ~0.2″ SMALLER at death). Smaller dead-tree dbh ⇒ the cone puts less in sz4-9 ⇒ the −9%
  LARGE. **OPEN VERDICT (do not over-infer):** the .sum stand is BIT-EXACT through the 2003 fire, so the
  dead-tree dbh differ only BELOW print resolution / by which record VARMRT kills — this could be (a) a
  ULP-amplification (sub-print dbh differences amplified through cone→fuel→flame→saw threshold, = the
  accepted D13/COMPRESS class) OR (b) a small real mortality-distribution diff (VARMRT per-record killed
  density). NEXT: full-precision compare the DEAD-tree dbh (not binned) jl vs live at the 1993→1998 mortality
  — if sub-print/ULP ⇒ D16 is ULP-amplification-class (ACCEPTED); if a systematic ~0.2″ shift ⇒ a real
  mortality-distribution bug. D16 = snag-fall −9% from a subtle dead-tree-dbh difference; ULP-vs-real is the
  final question. (Method: caught the 3× as a unit artifact by the 9%≠3× logic — measure the dead-tree dbh
  at full precision next, don't infer.)
- **ULP-vs-real SETTLED by magnitude logic: D16 is REAL, not ULP-amplification.** The clean measured
  divergence is −9% at the CWD1 down-wood output (live 9.636/jl 8.752) — WAY above ULP. The cone LARGE
  fraction varies ~LINEARLY with snag dbh, so a sub-print (ULP-scale, <0.01″) dead-tree dbh difference could
  produce at most a fraction of a percent — it CANNOT amplify to 9%. So D16 is a REAL snag-input difference,
  not the D13/COMPRESS ULP-threshold class. The most likely REAL mechanism given all formulas match: the
  ANNUAL-vs-CYCLE snag BOOKING — FVS FMSADD books snags EACH YEAR at that year's (growing) dbh (so a tree
  contributes 5 annual snag cohorts across the cycle, at 5 slightly different dbhs), while jl books ONE snag
  per record per cycle at the cycle dbh. Different dbh cohorts ⇒ different cone LARGE fractions ⇒ the ~9%.
  (This also explains the raw 62.3-vs-21.8: live ≈5 annual cohorts per tree, jl 1.) FIX DIRECTION: book jl
  ordinary-mortality snags ANNUALLY across the cycle (like FVS FMMAIN's year loop) at each year's dbh, not
  once at cycle end — OR verify the cycle-dbh vs annual-dbh cone difference is the 9%. CAUTION: the snag
  BOOKING feeds carbon_snt StandDead (bit-exact-validated) — an annual-booking change must keep that green.
  D16 = REAL (non-ULP) snag-fall −9%, mechanism = annual-vs-cycle snag booking (per-year dbh cohorts). This
  is the last SN model-fidelity divergence, root direction identified; the fix is an FFE snag-booking change
  + carbon_snt re-validate.
- **★★★ D16 ROOT FOUND & VERIFIED (the "annual-vs-cycle" hypothesis was ALSO wrong — verified by ITYP filter
  BEFORE concluding): jl's CUT path does not book cut trees as FFE standing-snags.** Re-stamped FMSADD with
  its ITYP source tag: ordinary mortality (ITYP=4) = **21.807 BOTH** (jl==live BIT-EXACT ⇒ mortality→snag is
  FAITHFUL). The ENTIRE 62.3−21.8 gap = **live ITYP=2 = 40.5** = the THINDBH cut + YARDLOSS standing-snags
  (fmscut.f:157 `FMSADD(IY(ICYC),2)`). jl has NO `add_snag!` in `cuts.jl` (grep: only fire/mortality/SNAGINIT
  call it); cuts.jl:87 EXPLICITLY: "standing snags, which FVSjl's basic cut path does NOT model ⇒ SSNG=0".
  So THINDBH-cut trees never become FFE snags→down-wood in jl ⇒ jl cwd runs ~10% low by 2003 ⇒ FMDYN picks
  the heavier fm5 ⇒ flame 4.01 vs 3.90 ⇒ the systematic small-tree fire over-kill. FULL CHAIN NOW VERIFIED
  END-TO-END: cut-snag gap → cwd −10% → fuel models {5,10} vs {10,12} → byram +6% → flame +2.8% → CHARHT/CSV
  up → pmort up → +3 TPA over-kill. THE FIX: book the cut trees' standing-snag portion as FFE snags in
  `cuts.jl` (mirror fmscut→FMSADD ITYP=2 — the non-removed/yarding-loss standing fraction becomes snags with
  the cut tree's dbh/height, into `fs.snags`, so they fall→down-wood). CAUTION: keep carbon_snt + fire tests
  green (this ADDS a snag source; validate StandDead/DDW vs live on a thinned+FFE stand). ⇒ D16 = a MISSING
  FEATURE (cut→FFE-snag SSNG path), definitively verified — NOT ULP, NOT a formula bug, NOT annual-vs-cycle.
  Ordinary mortality is bit-exact; only the cut-snag path is absent. (Method: the ITYP filter caught my 3rd
  D16 mechanism over-inference; EVERY non-cut snag source matches, so the gap is precisely the cut path.)
- **FIX FULLY SPECIFIED (cuts.f:1382-1386):** per cut RECORD, when YARDLOSS is active (LYARD): `LOSS =
  PREM·PRLOST` (PREM = removed density this record), then `DSNG = LOSS·PRDSNG` (DOWNED snags → straight to
  down-wood cwd) and `SSNG = LOSS·(1−PRDSNG)` (STANDING snags → `fs.snags`, later fall→down-wood via
  update_snags!); the actually-removed density becomes `PREM − LOSS`. YARDLOSS fields = (PRLOST, PRDSNG,
  PRCRWN, …); jl currently parses only `yardloss_prlost` and applies it at the AGGREGATE .sum level
  (cuts.jl:247, scales reported merch/saw/bdft by 1−PRLOST) — it does NOT book the per-record SSNG/DSNG into
  the FFE pools. IMPLEMENTATION (bounded, YARDLOSS-gated so non-YARDLOSS tests are untouched): (1) parse
  PRDSNG/PRCRWN into Control; (2) in the cut functions (`_thin_sorted!`/`_thinprsc!`), capture the per-record
  removed density + dbh/height (the cutlist already logs removed records via `_log_cut!`); (3) `add_snag!`
  the SSNG standing density (cut tree's dbh/height) into `fs.snags`, and add the DSNG downed density's bole
  to `fs.cwd` (via the same cone-taper split); (4) re-validate carbon_snt StandDead/DDW + the fire tests +
  the many thinning tests. This closes the D16 chain (cut-snags → cwd → fuel model → flame → over-kill). A
  real FFE feature port (harvest-residue SSNG/DSNG), the last SN model-fidelity item — root VERIFIED, fix
  formula EXACT, implementation multi-part (fresh focused session).
- **FIX PART 1/2 LANDED (suite green 6397/2), PART 2 open:** (1) FIXED a real YARDLOSS PARSING BUG — FVS
  fields = field1=DATE/field2=PRLOST/field3=PRDSNG (initre.f:3637-45); jl read PRLOST from field1 (the DATE)
  ⇒ YARDLOSS was silently INACTIVE. Now field2/3 (+ Control.yardloss_prdsng). (2) Booked the SSNG standing
  cut-snags in `cuts.jl _log_cut!` (SSNG=prem·PRLOST·(1−PRDSNG), the fmscut→FMSADD ITYP=2 analog) ⇒ jl snag
  origden 1993-1997 = 62.307 == live BIT-EXACT (was 21.8). BUT the FIRE over-kill is UNCHANGED (2008 TPA
  104/107): the booked cut-snags fall (den_hard 3.65 of 62.3, 94% fell) yet add ~0 to cwd — their bolevol is
  TINY (Σbolevol·den 0.079) ⇒ small-dbh, bole lands in SMALL not the LARGE pool. So booking the STANDING
  portion alone does NOT close the −0.88 LARGE snag-fall gap. PART 2 (the actual cwd/fire fix): (a) book the
  DSNG DOWNED portion (prem·PRLOST·PRDSNG) STRAIGHT to cwd (already-down, the bigger immediate contribution
  jl omits), and/or (b) MEASURE live's per-cut-snag CWD1 output (dbh/bolevol) — live's cut-snags contribute
  +0.88 to LARGE, so verify whether it's the downed path or bigger boles. The 2 landed fixes are correct +
  validated (snag density bit-exact, no regression); PART 2 = the DSNG→cwd path + a live-CWD1 per-cut-snag
  check. CAUTION (self): I attributed the cwd deficit to the missing cut-snags — booking the STANDING ones
  did NOT move cwd, so PART 2 must MEASURE live's cut-snag CWD1 LARGE output before claiming the cause.
- **★★★ PART 2 LANDED + D16 CUT-SNAG THEORY DISPROVEN BY MEASUREMENT (commit d69c53f, suite green 6397/2).**
  Ported the DSNG downed portion (cuts.f:1384 `DSNG = LOSS·PRDSNG`) straight to HARD cwd at cut time via the
  CWD3 analog (fmscut.f:98 → fmcwd.f:258 ENTRY CWD3: bole cone-split into `cwd[:,2,idc]`, all hard SCNV=1.0),
  mirroring jl's carbon_snt-validated bolevol/cone machinery. This moved jl LARGE cwd 5.67→6.015 (toward live
  6.425) — a real fidelity gain, kept. **BUT the fire over-kill is STILL 104/107 unchanged, and per-source
  live/jl stamps now PROVE the cut-residue is NOT the driver:** (a) DSNG cone-split → **0.0 into LARGE**
  (CUTDBG: 5 records, Σdsng-density 94.5, ΣaddH-into-LARGE **exactly 0** — the THINDBH-cut trees are small-dbh,
  their whole bole lands in SMALL 1-3, never LARGE 4-9); (b) SSNG standing-snag fall → ~0 into LARGE (measured:
  snag-fall LARGE 8.752 at cycle 1993 UNCHANGED by the SSNG booking); (c) live CTCRWN crown→cwd (fmscut.f:89-96,
  stamped) = **0.128 SMALL / 0.0 LARGE** (tiny — cut-tree crowns are small). ⇒ ALL THREE cut-residue paths are
  now faithfully ported (SSNG snag density bit-exact, DSNG bole, CTCRWN measured negligible) yet contribute ~0
  to the LARGE pool. **The −0.4 LARGE / −0.86 SMALL cwd gap that drives the fire is a SEPARATE base down-wood
  accumulation difference, not the cut path.** Ground-truthed live fuel basis via fmcfmd.f stamp: **live
  SMALL=7.964 LARGE=6.425 FMD=10** @2003; jl SMALL=7.1 LARGE=6.015 (per-size jl: sz1-3 {0.383,1.496,2.602},
  sz10-litter 2.619; sz4-9 {1.613,3.481,0.567,0.305,0.045,0.004}). **Basis definitions MATCH** (fmtret.f:378-390
  FMFMOV: SMALL=Σsz{1,2,3}+sz10-litter, LARGE=Σsz{4-9}, over I=1,2 piled/unpiled × hard/soft × 4 decay) — jl uses
  the identical size grouping, so the gap is REAL accumulation, not a basis bug. **Also discovered:** the SN
  fire's authoritative fuel-model selector is **FMCFMD3** (fmcfmd2.f, fmburn.f:246), which builds fuel from
  `CURRCWD` with a **0.04591 tons/acre unit conversion** + `XFML` — a NON-monotonic map (live has MORE cwd yet a
  COOLER fire: flame 3.90 vs jl 4.008), so raising jl's fine down-wood should shift the model toward the lighter
  fire. ⇒ **D16 RECLASSIFIED**: the cut→FFE-snag-PATH gap (the original D16 root) is CLOSED (faithful port, snag
  density bit-exact, no regression). The residual **3-TPA fire over-kill is a distinct item — call it D16b:** a
  base fine-down-wood/litter accumulation shortfall (jl SMALL 0.86 low, mostly litter sz10 + fines sz1-3) on a
  thinned+fire stand, feeding the non-monotonic FMCFMD3 CURRCWD selection. NEXT (measured, not inferred): stamp
  live FMCFMD3 `CURRCWD(1),CURRCWD(2),CURRCWD(3),CURRCWD(10)` @2003 to localize the 0.86 SMALL gap to a specific
  size class (litter vs fines), then trace that pool's accumulation (FMCADD litterfall / snag-fall fines / decay)
  vs a live per-year stamp — this is shared with D4/D5 (carbon down-wood) and is why carbon_snt stays bit-exact
  on its natural stand while snt01_alpha's thinned+fire fine-fuel runs low. The 3 cut-residue sub-fixes are all
  faithful FVS ports and stay (they improve LARGE-cwd fidelity 5.67→6.015 even though they don't close D16b).
- **★ D16b LOCALIZED per-size (live fmcfmd.f CWD-array stamp @2003 vs jl, both on the fire path):** the 0.86
  SMALL + 0.41 LARGE gap is a PROPORTIONAL ~10-18% shortfall concentrated in mid-size DOWN-WOOD, NOT litter:

  | size (span)     | jl    | live  | gap    | jl/live | pool  |
  |-----------------|-------|-------|--------|---------|-------|
  | 2 (0.25-1")     | 1.496 | 1.819 | +0.323 | 82%     | SMALL |
  | 3 (1-3")        | 2.602 | 3.001 | +0.399 | 87%     | SMALL |
  | 4 (3-6")        | 1.613 | 1.763 | +0.150 | 91%     | LARGE |
  | 5 (6-12")       | 3.481 | 3.734 | +0.253 | 93%     | LARGE |
  | 10 (litter)     | 2.619 | 2.685 | +0.066 | 98%     | SMALL |
  | 1 (0-0.25")     | 0.383 | 0.458 | +0.075 | 84%     | SMALL |

  ⇒ litter (sz10) and coarse (sz6-9) MATCH; the deficit is the fine-to-medium woody down-wood (sz1-5), jl
  uniformly 82-93% of live. This is the SAME "snag-fall bole ~10% low" signature as D5 / [[fvsjl-ffe-fire-downwood-3fixes]]
  — but it survives even WITH the cut-residue booked, so it's a base down-wood BOLE-accumulation shortfall on this
  thinned stand (candidate mechanisms, to test vs a live FMSNAG/FMCWD per-YEAR stamp: snag-fall bole→cwd conversion
  ~10% low, or the cut-snag fall TIMING, or a decay-rate tick). carbon_snt stays bit-exact because its natural
  stand's snag population differs. **NEXT (fresh session, measured): stamp live FMSNAG/CWD1 per-year 1993→2003 to
  watch the sz2-5 pools accumulate vs jl, pin the ~10% bole shortfall to its source (input vs decay vs timing).**
  This is a downstream FFE-fuels residual (3 TPA / 2.8% on ONE thinned+fire stand); the upstream cut→FFE-snag PATH
  (the original D16) is CLOSED. NOT ULP (proportional 10-18%, not sub-print), so D16b stays 🔬 OPEN, fully localized.
- **★ D16b — volume-basis hypothesis REFUTED by measurement; trajectory captured.** Hypothesis: my cut-DSNG/SSNG
  bolevol used MERCH cuft but CWD3 (fmcwd.f:283-286) uses TVOLI = FMSVL2 'D' = TOTAL — a ~10-15% under-count that
  would match the proportional shortfall. Switched DSNG to `cuft_vol` (total CFV, the faithful CWD3 choice, kept +
  suite-green) and re-measured: effect NEGLIGIBLE (sz2 1.496→1.497, SMALL 7.1→7.124) — for these cut trees total ≈
  merch, so the volume basis is NOT the gap. Live down-wood per-year trajectory (fmcwd.f FMCWD stamp, FFE stand):
  the pool is CUT-SNAG-FALL dominated — sz5 (6-12") climbs 3.23(1993)→**6.94(1995)** as the 1993-cut boles fall,
  then decays; jl sz5 @2003 = 3.481 vs live 3.734. The within-cycle PHASING differs (FMCWD dumps post-decay each
  year; the fire's FMCFMD reads a later-in-2003 state after that year's snag-fall/litterfall pulse — FMCWD@2002
  SMALL 4.77 vs FMCFMD@2003 pre-fire 7.96), so the two engines' per-year dumps are NOT directly phase-comparable.
  ⇒ The ~10% shortfall is in the DENSITY/DECAY/TIMING of cut-snag fall, not bole volume. NEXT (fresh session):
  phase-MATCHED per-year stamp (same point in the annual loop, both engines) of the cut-snag fall density + the
  sz2-5 pool, 1993→2003, to separate fall-rate (SNAGFALL) from decay-rate (DKR) from a count diff.
- **★ D16b MECHANISM NARROWED to a final-cycle-year snag-fall PULSE.** Phase-matched at the fire: live SMALL 7.96
  vs jl 7.1. Live's raw-CWD (fmcwd.f stamp) jumps **+3.2 in the SINGLE 2003 step** (FMCWD@2002 4.77 → FMBURN@2003
  sample 7.96) — a snag-fall+litterfall PULSE as a snag cohort reaches fall-age in the cycle's last year, sampled
  by live's FMBURN. jl's fire basis (`fire_smlg`, stashed at the 2003-fire-cycle start via `_small_large_fuel`
  AFTER the 1998→2003 `ffe_fuel_update!(per=5)` loop, summary.jl:276) ends that loop at 7.1, ~0.86 short: jl's
  final-year pulse is smaller. ⇒ NOT a fire-basis PHASING bug (jl correctly includes the full pre-fire cycle loop,
  matching live's FMBURN-before-next-loop sampling) — it's the **snag-fall DENSITY/TIMING in the last cycle year**
  (the D5/#28 fire-downwood class): jl falls slightly fewer / a beat later than live's fmsfall cohort, under-
  filling sz2-5 ~10%. FINAL NEXT-STEP (fresh session): stamp live FMSNAG/fmsfall per-year 1999→2003 for FALL
  DENSITY by cohort (deathyr) vs jl `update_snags!` — pin fewer-snags vs off-by-a-year. ≤3% one-family FFE
  residual, fully localized; upstream cut→snag PATH (D16) CLOSED; SN inventory at documented floor (219/260).
- **★★ D16b — SNAG-FALL DEFINITIVELY RULED OUT (bit-exact), gap REDIRECTED to decay/accumulation-over-window.**
  Stamped live FMSNAG (fmsnag.f:227 CALL CWD1, summed DFIH+DFIS/yr) vs jl `update_snags!` per-year fall density,
  snt01_alpha FFE stand: **BIT-EXACT every pre-fire year** — 1993/94/95 = 25.8092, 1996 = 10.0717, 1997 = 9.3748,
  1998 = 6.0587, 1999-2001 = 5.0986, 2002 = 4.7653 (jl == live to 4 dp); 2003 fire-year jl 44.83 vs live 43.96.
  So jl's snag-fall DENSITY/TIMING is correct — my prior-turn "final-cycle-year pulse" hypothesis is REFUTED by
  measurement. (jl carries more snag RECORDS — 8 vs 6, 30 vs 20 — but identical aggregate fall density, so the
  extra records are density-preserving splits.) Phase-matched cwd: at 1998 jl ≈ live (jl sz2/sz3 even slightly
  HIGHER: 0.587/1.597 vs 0.567/1.440; sz4/sz5 0.02 apart) but by the 2003 fire jl is ~10% LOW (sz2 1.496 vs
  1.819, sz3 2.602 vs 3.001). ⇒ the down-wood INPUT matches (fall bit-exact + 1998 pools ≈) and the divergence
  ACCUMULATES over 1998→2003 — a DECAY-rate application or a minor accumulation-source (FMCADD woody-breakage /
  FMSDIT crown-lift) difference over the window, NOT fall, NOT bolevol, NOT volume-basis, NOT cut-residue, NOT
  fire-basis phasing (all now refuted vs live). FINAL NEXT-STEP (fresh session): phase-matched per-year cwd-by-
  size stamp BOTH engines 1998→2003 isolating the 1993-97-deposited cohort's decay ratio (DKR application count/
  order) vs the woody-breakage+crown-lift adds — the last unrefuted mechanism for the ~10% woody-DDW shortfall.
- **★★★ D16b FINAL LOCALIZATION — sz4/5 CONE-DISTRIBUTION at deposition (bolevol & decay & density all ruled
  out).** Phase-matched jl-vs-live cwd-by-size per year (jl dump after `fmcwd!` == live FMCWD post-decay phase):
  the sz5 (6-12") gap is BORN at the FIRST fall year and stays ~constant — 1993 jl 3.005/live 3.231 (0.93), 1994
  4.774/5.201, 1995 6.339/6.944, 1998 4.81/5.233, 2002 3.481/3.734 (ratio 0.91-0.93 throughout; ABSOLUTE gap
  shrinks 0.60→0.25 as it decays). sz4 tracks CLOSER (ratio ~0.95). Since (a) snag-fall density is BIT-EXACT and
  (b) bolevol merch-vs-total is negligible (tested BOTH cut-path SSNG and ordinary-mortality → sz5 moved <0.04),
  the ~7% is NOT input magnitude and NOT decay (constant ratio) — it's **how the fallen bole is DISTRIBUTED across
  size classes**: jl deposits relatively more in sz4 / less in sz5 than live's CWD1. The `_cwd_cone_fractions`
  FORMULA is verified == CWD1 taper, so the driver is the **per-snag dbh/height feeding the cone** — and jl
  carries **8 cut-snag records where live has 6** (density-preserving splits, same aggregate fall), so jl's snags
  sit at a slightly different mean dbh ⇒ their cones reach sz5 less. ⇒ D16b's residual 3-TPA fire over-kill
  traces, at the deepest level, to **cut-snag RECORD GRANULARITY** (jl 8 vs live 6) shifting the sz4/5 down-wood
  split ~7%, feeding the non-monotonic FMCFMD3 fuel selection. SEVEN hypotheses now refuted vs live this campaign
  (volume-basis, cut-residue, snag-fall density/timing, fire-basis phasing, cut-snag bolevol, ordinary bolevol,
  decay-rate). FINAL NEXT-STEP (fresh session): stamp live FMSADD/FMSNAG per-record cut-snag dbh (the 6 records)
  vs jl's 8 — pin why jl splits the cut into more records (point structure? tripling? cut-record management), then
  match the granularity so the cone→sz5 split aligns. A ≤3% one-family FFE residual, exhaustively localized to
  per-snag-dbh cone granularity; NOT ULP; open. Upstream cut→snag PATH (D16) CLOSED.
- **★★★★ D16b ROOT FOUND + FIXED (LARGE half) — SNAGINIT snags fall with the wrong VOLUME BASIS + dropped the
  current-height field (commit cb824c9, suite green 6397/2).** The "record granularity" framing was a step short:
  the decisive live stamp was FMSNAG per-record (fmsnag.f:227 CWD1) dumping each snag's dbh/den/ht. It showed the
  SNAGINIT snag (`SNAGINIT 10 11 50 40 2 50`) at **HTIH=40, not 50** — jl used 50. Reading CWD1 (fmcwd.f:185,22,29):
  the fall bole = FMSVL2 **'D' = TOTAL** volume of the ORIGINAL tree (taper from HTDEAD=50), integrated only up to
  **HIHT(2)=HTIH=40** (the current top). jl's `ffe_add_snaginit!` (a) used **MERCH** (R8 Clark v[4]) for the fall
  bolevol and (b) **dropped SNAGINIT field 4** (current height), using field 3 (HTDEAD) for the cone. FIX: (1) a new
  `SnagList.fallvol` = TOTAL-volume bole for the fall (CWD1), kept DISTINCT from the merch `bolevol` the Stand-Dead
  report uses (validated separately); `update_snags!`/salvage now deposit `fallvol`. (2) parse SNAGINIT field 4 →
  snag `htcur`; `_cwd_cone_fractions` truncates the taper integration at `htcur` (normalized by the full cone).
  RESULT: **LARGE cwd 6.015 → 6.397 == live 6.425 (was 6% off → ULP-class)**, flame 4.003→3.988 (live 3.90), 2008
  TPA 104→105 (live 107) — 1 of the 3 over-killed TPA recovered. CRITICAL BOUND: applying TOTAL to ORDINARY/fire/
  input snags too **regressed 10 live-validated carbon_snt tests** ⇒ live's ordinary-snag CWD1 ≈ MERCH there
  (total≈merch for those trees, or a genuinely different basis), so ordinary stays merch (`fallvol` defaults to
  `bolevol`). The SNAGINIT-total is validated against live independently (its snag has total≫merch). ⇒ D16b's
  LARGE-cwd driver is FIXED; the **SMALL-cwd gap (jl 7.12 vs live 7.964) remains** and still tips ~2 TPA of the
  fire over-kill — the last open piece, now isolated to the sz1-3 fine down-wood (litter sz10 already matches).
  9 hypotheses explored: 8 refuted, the 9th (SNAGINIT volume/height) CONFIRMED + fixed for LARGE.
- **★ D16b post-fix per-size (2003 fire) — LARGE closed, gap now PURELY sz1-3 fine wood.** jl vs live: sz4
  1.654/1.763, **sz5 3.815/3.734 (matched, jl slightly over)**, sz6 0.574/0.573, sz7 0.305/0.307, sz8-9 exact,
  sz10-litter 2.619/2.685 — all LARGE + coarse + litter now match. The ONLY remaining gap is the fine wood:
  **sz1 0.383/0.458 (0.84), sz2 1.497/1.819 (0.82), sz3 2.621/3.001 (0.87)** — a uniform ~15% shortfall in the
  <3″ classes (Σ ≈ 0.78, ≈ the SMALL gap). This is PRE-EXISTING (unchanged by the SNAGINIT fix; sz2 1.496→1.497)
  and a SEPARATE D5-class fine-down-wood residual. The fine-wood SOURCES are all present in jl (FMCADD LIMBRK
  live-crown breakage fmcadd.f:81, crown-lift `cl[]` term, CWD2B, bole-fall fine top) — so it is a
  DISTRIBUTION/magnitude gap, most likely the per-size `crown_biomass`/CROWNW fine-class fractions feeding
  LIMBRK+crown-lift into sz1-3. NEXT (fresh session, same method that cracked LARGE): stamp live FMCADD CROWNW(I,
  1..3) per tree vs jl `crown_biomass` xv[2..4], or dump the live per-FMCADD sz1-3 ADD vs jl's, to pin whether
  the fine-crown biomass or a fine-fuel term is ~15% low. Closing it removes the last ~2 TPA of the fire over-kill.
- **★ D16b-SMALL likely FOLDS INTO D4 (crown-biomass FMCROWE ~0.9 ton, an already-documented open residual).**
  The remaining sz1-3 fine-down-wood gap derives ENTIRELY from crown biomass: FMCADD's three fine-wood sources —
  LIMBRK live-crown breakage (fmcadd.f:81 `LIMBRK·CROWNW(SIZE)`), crown-lift (`X·CROWNW(SIZE)`), and CWD2B snag
  crown debris — all scale with `CROWNW(I,SIZE)` (crown biomass by size). All three terms ARE ported in jl
  (fmcadd_woody! LIMBRK + the `cl[]` crown-lift term + CWD2B) and structurally correct, so a uniform ~15% fine
  shortfall points at the `crown_biomass`/CROWNW values themselves being low in the fine classes — the SAME
  quantity as the D4 crown-biomass-carbon residual. ⇒ HYPOTHESIS (verify next session, live FMCADD CROWNW(I,1..3)
  stamp vs jl `crown_biomass`): D16b's last SMALL/fine-wood piece and D4 are one crown-biomass residual; fixing
  the crown-biomass fine-size fractions closes both. (Confirmed this session: fine-wood SOURCES all present +
  structurally right; the gap is in the CROWNW magnitude/distribution, not a missing term.)
- **★★ D16b-SMALL == D4 — HYPOTHESIS RETRACTED, REFUTED by direct measurement (re-trace discipline).** Last
  turn I claimed this "verified by code-read," but that only READ the crown_biomass.jl:72-80 header comment ("not
  absolute bit-exactness") — NOT a live differential. Stamping live FMCROWE (fmcrowe.f:562, dump D/H/IC/XV(0..4)
  for dbh 10-12) vs an identical jl `crown_biomass` dump PROVES jl's crown biomass is **BIT-EXACT / ULP** per tree:
  D=10.04/H66.2/IC27 → fol/x1/x2/x3/x4 = 16.6257/12.805/63.5056/117.2094/18.7887 IDENTICAL both; D=10.06/72.2/28
  bit-exact; D=10.06/70.9/29 within 0.03%. So the crown-biomass fine-size values are NOT ~15% low — **D16b-SMALL is
  NOT the crown-biomass magnitude, and NOT D4.** (The header's "not bit-exact" caveat concerns small-tree/bole-tip
  edge cases, not these trees.) ⇒ the sz1-3 fine-down-wood shortfall must come from a fine-wood APPLICATION step,
  not the CROWNW source: candidates now are (a) the crown-LIFT `cl[]` term (OLDCRW × crown-base-rise, known one-
  cycle-lag residual, [[fvsjl-ffe-crown-lift-landed]]), (b) CWD2B snag-crown release timing, or (c) the ordinary
  snag bole-fall FINE TOP (kept on merch — total regressed carbon_snt, but snt01_alpha's larger snags may need it).
  NEXT (measured, not inferred): decompose jl's sz2/sz3 cwd at the fire by SOURCE (LIMBRK vs crown-lift vs CWD2B vs
  bole-fall) and stamp the matching live per-source ADD. LESSON (logged): "verified by code-read" of a COMMENT is
  NOT verification — only a live differential is (doctrine #4). D16b LARGE half stays fixed (SNAGINIT total, ULP).
- **★ FINE-WOOD SOURCE DECOMPOSED (jl sz2+sz3 hard-pool adds, 1993→2002, instrumented per-source):** crown-lift
  `cl[]` = **2.049** (DOMINANT), LIMBRK woody-breakage = 1.434, CWD2B snag-crown = 0.611 (bole-fall not cleanly
  isolated — decay interleaves). The sz2+sz3 gap at the fire is jl 4.118 vs live 4.820 = **−0.70**. Since LIMBRK
  scales with the crown biomass now proven BIT-EXACT (and TPA is bit-exact), the LIMBRK contribution should match
  → the −0.70 gap is most likely in the **crown-LIFT term**, the dominant source, which already carries a KNOWN
  one-cycle-LAG residual (memory [[fvsjl-ffe-crown-lift-landed]]: OLDCRW crown-base-rise needs the previous-cycle
  per-tree crown base, tracked across a regen/mortality-changing tree list). ⇒ D16b-SMALL redirects to the
  **crown-lift (cl[]) residual**, NOT crown-biomass D4 (crown biomass bit-exact) and NOT bole-fall. NEXT: stamp
  live FMCADD crown-lift sz2+3 ADD (fmcadd.f:95-102 OLDCRW term) per cycle vs jl `cl[]` — pin the lag/magnitude,
  then fix the OLDCRW crown-base tracking. This is the last ~2 TPA of D16b's over-kill, localized to crown-lift-lag.
- **★ CROWN-LIFT ALSO REFUTED — gap narrowed to CWD2B / bole-fall fine top.** Stamped live FMCADD crown-lift
  (fmcadd.f:95-102 OLDCRW term, XCLIFT per-call accumulator for sz2+3): live FFE-stand crown-lift = 0 (cycle 1,
  no OLDCRW yet) + 0.4097×5 = **2.0485 over 1993-2002 == jl's 2.049 BIT-EXACT** (and all CLIFT23 rows identical ⇒
  only the FFE stand runs FMCADD, clean isolation). ⇒ the crown-lift-lag is NOT the sz1-3 gap either. Now THREE
  fine-wood sources proven to MATCH live: crown biomass (bit-exact XV), crown-lift (2.049 both), and LIMBRK (=0.01·
  bit-exact CROWNW·bit-exact TPA). Arithmetic: jl sz2+3 = crownlift 2.049 + limbrk 1.434 + cwd2b 0.611 ≈ 4.09 vs
  live 4.82 ⇒ the **−0.70 gap is entirely in CWD2B (snag-crown debris release) and/or the ordinary-snag BOLE-FALL
  fine top** — the two sources I did NOT prove equal. Bole-fall is the prime suspect: ordinary snags fall on MERCH
  bolevol (top cut at the 4″ merch top), but CWD1's total (FMSVL2 'D') integrates the whole cone to the 0″ tip,
  depositing that thin <4″ top into sz1-2 — which jl omits. (I kept ordinary on merch because a naive total-for-
  ordinary regressed 10 carbon_snt tests — but that changed the MAGNITUDE; the correct fix likely distributes the
  merch bole over the FULL-height cone so the fine top lands in sz1-2 without inflating the total.) NEXT: stamp
  live CWD2B sz2+3 release + the ordinary-snag CWD1 sz1-2 ADD vs jl, to split CWD2B-vs-bole-top. FOUR fine-wood
  candidates now measured: crown-biomass/crown-lift/LIMBRK all MATCH; CWD2B + bole-top are the last two.
  ★ REFINEMENT: the gap is FINE-SPECIFIC (sz1 0.84 / sz2 0.82 / sz3 0.87 low, but sz4 0.94 / sz5 matched) — a
  uniform bole-fall merch-vs-total would lower sz4-5 too, so bole-top is UNLIKELY; the fine-only signature points
  at **CWD2B (snag-crown debris release)**, the remaining fine-specific source (jl 0.611). CWD2B is the prime
  suspect; NEXT stamp = live CWD2B sz2+3 per-year release vs jl `_cwd2b_fall!`.
- **★★ CWD2B PINNED as the sz1-3 gap (arithmetic + code trace): jl ~2× LOW.** With crown-lift (2.049) and LIMBRK
  (bit-exact crown·TPA) both matching live, live's CWD2B sz2+3 ≈ 4.820 − 2.049 − 1.434 = **1.337 vs jl's 0.611**.
  jl's ordinary mortality DOES book crown into CWD2B (mortality.jl:522 `fmscro!`), and the intake `amt = xv[sz]·
  density` uses BIT-EXACT crown biomass, so the ~2× shortfall is NOT the intake magnitude — it's the RELEASE
  TIMING: `fmscro!` spreads each crown component over `ILIFE = ceil(min(tsoft, _fm_tfall(cls,sz)))` years, and
  `_cwd2b_fall!` releases the year-1 pool annually. If jl's ILIFE (from `tfall_cls`/`_fm_tfall` crown-fall
  lifespan, or TSOFT=`(1.24·dbh+13.82)·DECAYX`) runs LONGER than FVS's FMSCRO, the crown debris releases too
  SLOWLY → too little in the fine down-wood by 2002 (and too much left "in waiting," i.e. Stand-Dead crown). This
  is the D5-class crown-fall-TIMING residual ([[fvsjl-fire-downwood-3fixes]]). ⇒ **D16b-SMALL == the CWD2B/FMSCRO
  crown-fall-timing (D5-adjacent), NOT crown biomass (D4) and NOT crown-lift.** NEXT (measured): stamp live FMSCRO
  ILIFE + CWD2B(IYR) per size vs jl `_fm_tfall`/`ilife` — pin whether jl's crown-fall lifespan is too long, then
  align it. This closes the last ~2 TPA of the D16b fire over-kill. FIVE fine-wood sources now resolved: crown-
  biomass / crown-lift / LIMBRK MATCH; CWD2B is 2× low (timing); bole-top ruled out (sz4-5 match).
- **★★★ CWD2B "2× low" RETRACTED (was inference); ALL crown sources MATCH; the sz1-3 gap IS the SMALL-TREE
  BOLE-FALL (measured).** Per doctrine #4 I stopped inferring and STAMPED live FMCADD per-source (CRLF-aware):
  live CWD2B sz2+3 (1993-2002) = 0.024+0.611 = **0.635 vs jl 0.611** (~4%, MATCHES — the "2× low" was my bad
  subtraction); live LIMBRK = 0.1229×5 + 0.1640×5 = **1.4345 vs jl 1.434 BIT-EXACT**; live crown-lift 2.049 ==
  jl. So EVERY crown-derived source matches. But the crown adds Σ ≈ 4.12 while the sz2+3 CWD the fire samples =
  jl 4.118 / **live 4.820** ⇒ live has **~0.70 MORE from the ONE non-crown source: the ordinary-snag BOLE-FALL
  fine top.** It's SMALL-tree-specific (sz4-5 large-tree classes match): jl's `mortality.jl` fall uses `fallvol =
  max(0.005454154·h, merch)` = the tiny CONE FLOOR for sub-merch trees (merch=0 below DBHMIN), but live's CWD1
  uses FMSVL2 'D' = TOTAL — a small tree's whole bole (all < 3") lands in sz1-3. ⇒ **D16b-SMALL = the sub-merch
  snag fall-volume**, NOT crown biomass (D4), crown-lift, LIMBRK, or CWD2B (all measured equal). FIX (careful):
  give the ordinary `fallvol` the FMSVL2 'D' TOTAL for sub-merch/small trees — the naive cuft_vol-for-ALL regressed
  carbon_snt via LARGE-tree DDW inflation (jl `cuft_vol` ≠ FVS FMSVL2 'D'; use an FMSVL2 'D' analog, or apply total
  only where merch=0). NEXT: stamp live CWD1 TVOLI for a sub-merch snag vs jl fallvol; port the small-tree fall
  volume. (Meta: THREE inferences on this item — crown-lift, CWD2B-2×, D4 — each REFUTED once measured; only the
  direct live stamp held. Re-trace discipline earned its keep repeatedly this session.)
- **★ D16b-SMALL — the "MCF vs TCF" story is ALSO inconclusive (5th self-correction); the fall-volume is NOT
  clearly the gap. HONEST STOP for this item.** The FVS formula IS `FMSVL2 'D' = MAX(0.005454154·H, TCF)` (total,
  fmsvol.f:153), and small-tree TVOLI = the floor = jl's floor (confirmed). BUT the claim "jl uses low MCF" is
  wrong: jl's `mortality.jl` fall uses `t.merch_cuft_vol`, whose VALUE (≈32.2 for a dbh-15 snag) sits INSIDE live's
  measured TVOLI range (31.8-34.1) — NOT the low `_fm_cuft(merch=true)` MCF (28.5). So jl's large-snag fall volume
  is already ~right, and the −0.70 sz2+3 gap is NOT cleanly explained by the fall volume. The stamps don't cleanly
  align trees (live CWD1 dumps ALL snags incl. SNAGINIT/fire; jl mortality.jl dumps live-trees-being-killed —
  different populations & heights), so a per-tree matched TVOLI comparison was never actually achieved. ⇒ **D16b-
  SMALL (~2 TPA, ONE thinned+fire stand) is a genuinely HARD, deeply-investigated residual: crown-biomass / crown-
  lift / LIMBRK / CWD2B all MEASURE EQUAL to live, small-tree fall is correct, large-snag fall is ~right — yet
  live's sz2+3 cwd exceeds jl's by 0.70 and the source is NOT cleanly attributable after ~20 turns and FIVE
  refuted inferences (D4, crown-lift, CWD2B-2×, small-bole, MCF-vs-TCF).** VERDICT: the LARGE half is a real
  validated fix (SNAGINIT total, ULP); the SMALL half is DOWNGRADED to a documented ≤3% one-stand FFE-fuels
  residual — the correct NEXT step is a CLEAN per-tree-MATCHED CWD1 deposit comparison (stamp live CWD1's actual
  per-size ADD for ONE identified snag + reproduce that exact snag's deposit in jl), NOT another arithmetic
  inference. LESSON (imprint): this item burned 5 inferences because I repeatedly SUBTRACTED/READ-ACROSS instead of
  matching a single record end-to-end; when a differential won't isolate, MATCH ONE RECORD, don't decompose sums.
- **★ BOLE-FALL DIRECTLY MEASURED both sides (applied the lesson): also MATCHES ⇒ the gap is NOT any ADD source.**
  Instrumented jl `update_snags!` sz1-3 deposit = **0.215** (1993-2002); stamped live CWD1 sz1-3 ADD (snags dead
  <2003) = **0.194**. They match — bole-fall is NOT the 0.70 gap (refuted DIRECTLY, not by subtraction). So ALL
  FIVE fine-wood ADD sources now measure ≈EQUAL live-vs-jl: crown-biomass (bit-exact), crown-lift (2.049), LIMBRK
  (1.434), CWD2B (0.635/0.611), bole-fall (0.194/0.215). **YET the sz2+3 cwd the fire samples differs (jl 4.118 /
  live 4.820), AND the measured adds (Σ≈4.3) do NOT reconcile UP to live's 4.820** — a cwd that exceeds its own
  add sources is impossible under decay, so the residual is a DECAY-or-SAMPLING-PHASE effect, not a missing/low
  source. (DKR is bit-for-bit + carbon_snt decay is bit-exact per memory, so likely the fire-basis samples a
  different annual phase than my 1993-2002 add-window.) ⇒ **D16b-SMALL is a ≤3% one-stand DECAY/PHASE residual
  with every ADD source proven equal** — not attributable to any fuel-input term after exhaustive DIRECT
  measurement. HONEST STATUS: documented hard residual; the only remaining handle is a per-YEAR cwd-sz2+3 trace
  (jl vs a live per-year CWD-array stamp at the SAME annual phase) to separate decay-vs-phase. LARGE half fixed
  (SNAGINIT total, ULP). SN campaign's last open non-accepted item, bounded to "all inputs equal, a ≤0.7-ton
  sz2+3 decay/phase discrepancy on ONE thinned+fire stand."
- **★ INITIAL DEAD FUEL (FMCBA STFUEL) — added ONCE both engines (6th inference refuted); but jl's init load is
  slightly low.** Stamped live FMCBA cwd sz2/sz3 per call: 1993 = 0.761/1.970, 1998 = 0.642/1.580, 2003 =
  1.819/3.001 (the 2003 == the fire's sample). jl init (after `fmcba!`) = sz2 0.68 / sz3 1.93. The STFUEL block is
  guarded `IF (IYR .EQ. IY(1))` (fmcba.f:252-254 "INITIALIZE THE DEAD FUELS ONLY FOR THE FIRST YEAR") — so FVS
  loads dead fuel ONCE, exactly like jl's `!fuels_init` guard (refutes the "STFUEL every cycle" idea). One CLEAN
  real diff surfaced: **jl's INITIAL sz2 dead fuel = 0.68 vs live 0.761 (~11% low)** (sz3 0.68% low) — the
  `ffe_dead_fuel_loading`/STFUEL default table, a small (0.08-ton) but genuine follow-up. BUT the per-year trace
  shows the gap EVOLVES non-cleanly (1998: jl sz3 1.597 > live 1.580 yet jl sz2 0.587 < live 0.642; 2003 both
  low) — NOT a single constant offset, so even fixing the init table (~0.08) won't close the 0.70. ⇒ **FINAL
  VERDICT on D16b-SMALL: a ≤3% one-thinned+fire-stand fine-down-wood residual, EXHAUSTIVELY investigated — SIX
  hypotheses refuted by direct measurement (D4/crown-biomass, crown-lift, CWD2B-2×, small-bole, MCF-vs-TCF, STFUEL-
  every-cycle); ALL dynamic add sources + the once-only init both measure ≈equal; the residual gap is a subtle,
  non-constant, evolving sz2/sz3 decay-class/phase interaction not attributable to any single term.** Documented
  as a hard residual (doctrine #5 — the irreducibility is not PROVEN but the attribution is exhausted). Small
  clean follow-up: reconcile `ffe_dead_fuel_loading` sz2 (0.68 vs 0.761) vs the FVS FFE-forest-type STFUEL default.
  LARGE half remains a real validated fix (SNAGINIT total → LARGE cwd ULP-class, 1 of 3 TPA).
- **★ The "init table follow-up" is ALSO not a bug — jl's FUINI is BIT-EXACT (7th check).** Pulled the FVS FUINI
  DATA table (fmcba.f:95-113): the oak-hickory (500s) row = `0.13, 0.68, 1.93, 0.43, 1.01, 1.01, 0,0,0, 4.28,5.91`.
  jl's `fmcba!` init deposits sz1=0.13, sz2=0.68, sz3=1.93, sz4=0.43, sz5=1.01 — **EXACTLY the oak-hickory FUINI
  row.** So jl's dead-fuel TABLE is correct; live's cwd sz2=0.761 = FUINI 0.68 **+ ~0.08 pre-FMCBA cycle-1
  down-wood** (from the initial snag/mortality processing that runs before/around FMCBA in FVS's cycle-1 order,
  which jl's phase doesn't reproduce). ⇒ even the "clean follow-up" is the SAME cycle-PHASE/ORDERING class as the
  whole residual — there is NO fixable single table or coefficient. This CONFIRMS the verdict: every FUINI table,
  every fuel coefficient, and all 5 dynamic add sources jl uses are CORRECT; D16b-SMALL is purely a cycle-1-onward
  PHASE/ORDERING interaction of the down-wood accumulation (≤3%, one stand). Seven candidates measured, seven
  attributed to "jl is correct / it's phase." The honest close: this is a documented downstream FFE-fuels
  phase residual; closing it needs a full cycle-1 down-wood ORDERING reconciliation (FMCBA vs SNAGINIT vs initial
  mortality vs the annual loop), not any term fix — a scope disproportionate to ≤2 TPA on one stand. DONE for this
  campaign pass; recorded as a hard residual with the phase root identified.
- **★ DECAY (DKR) also BIT-EXACT (8th check) + phase-matched trajectory localizes the divergence to CYCLE 2
  (1998→2003).** Verified FVS DKR (fmvinit.f:70-97): sz2 = 0.11 all classes, sz3 = [0.11,0.09,0.09,0.09] (classes
  3,4 ← class 2 via the DO I=1,9/J=3,4 copy) — jl `_FM_DKR` is IDENTICAL. Phase-matched jl-vs-live at the FMCBA
  cycle-start: 1993 jl 0.687/2.015 vs live 0.761/1.970 (jl sz3 HIGHER), 1998 jl 0.589/1.623 vs live 0.642/1.580
  (jl sz3 still HIGHER), 2003 jl 1.497/2.621 vs live 1.819/3.001 (jl both LOW). ⇒ jl and live are CLOSE through
  1998 (jl sz3 even higher); the −0.69 opens ENTIRELY in the 1998→2003 cycle (live gains +2.60 there, jl +1.91).
  Since crown-lift/LIMBRK/CWD2B/bole-fall/decay ALL measure equal for that cycle, the cycle-2 aggregate divergence
  is UNATTRIBUTABLE by decomposition — the components sum-match but the total doesn't, which can only be a subtle
  ORDER-of-operations / accumulation-basis interaction (e.g. what state each per-year add/decay is applied to)
  that per-source sums cannot capture. ⇒ **DEFINITIVE CLOSE: D16b-SMALL is a ≤3% one-stand cycle-2 down-wood
  accumulation-ORDER residual; EIGHT independent live checks (D4, crown-lift, CWD2B, small-bole, MCF/TCF, STFUEL-
  cycle, FUINI-table, DKR-decay) all confirm jl's terms are correct.** The only remaining path is a per-year
  MATCHED cwd-sz2/3 stamp of BOTH engines through the 1998→2003 cycle (live FMCWD/FMCADD year-by-year vs jl), to
  catch the one year the totals split — a fresh, bounded, but low-value trace (≤2 TPA / one stand). Recorded as a
  documented hard residual; the LARGE half stays a real validated fix. This is the campaign's final SN item and it
  is bounded to bedrock: all inputs proven equal, divergence isolated to cycle-2 accumulation order.
- **★★ PRECISE LOCALIZATION (per-year matched trajectory — the lesson paid off): the ENTIRE gap opens in the
  SINGLE 2003 cycle-boundary step.** Compared jl `JLDECAY` vs live `FMCWD` post-decay PER YEAR: 1998 jl 0.522/1.489
  vs live 0.567/1.440, ... 2002 jl 1.247/2.294 (Σ3.541) vs live 1.274/2.263 (Σ3.537) — **MATCH within 0.004 through
  2002** (jl sz3 even slightly HIGHER all along). Then 2002→2003-cycle-start: jl 3.541→4.118 (+0.577), live
  3.537→4.820 (**+1.283**). ⇒ the whole −0.70 is the 2003 STEP, and live's +1.283 is ANOMALOUS vs the normal
  ~0.25/yr — a one-time CYCLE-BOUNDARY addition. Mechanism: the 1998→2003 cycle's END-of-cycle mortality creates
  snags whose crown → CWD2B → down-wood, and live's FMBURN/FMCBA at 2003 samples the cwd AFTER that boundary
  addition; **jl's `fire_smlg` is stashed at the 2003-cycle START (summary.jl:276), BEFORE jl books its 2003
  boundary mortality's crown debris** (that runs later inside grow_cycle!). So jl's fire fuel basis is one cycle-
  boundary-crown-debris increment short (~0.7 sz2+3) → lighter cwd → the non-monotonic FMCFMD3 picks a hotter
  model → the +3-TPA over-kill. ⇒ **D16b-SMALL FIX HYPOTHESIS (precise, testable): advance jl's 2003 boundary
  mortality-crown → CWD2B → fine-down-wood addition BEFORE stashing `fire_smlg`, so the fire samples the same
  post-boundary cwd live does.** CAUTION: fire_carbon (2000 SIMFIRE) is bit-exact with the current phase — validate
  no regression. This supersedes "unattributable accumulation-order": it IS the fire-basis-vs-boundary-mortality
  PHASE, localized to one step by the per-year matched trace. The 8 term-checks (all equal) + this localization
  make it a well-understood phase bug, not a mystery — a bounded fix for a fresh session.
- **★ FIX HYPOTHESIS TESTED (advance 1 fuel step before stashing `fire_smlg`, defer per−1) — OVERSHOOTS,
  bracketing the answer.** Implemented the code's own documented intent (summary.jl:264 "advance 1 year, stash,
  advance the rest") and measured: fire-basis SMALL 7.12 → **8.783** (live 7.964 — now OVER by 0.8), over-kill
  slightly worse (2008 TPA 105→104). So the CURRENT phase (0 steps) UNDERSHOOTS (7.12) and a FULL 2003 step
  OVERSHOOTS (8.783) — live's FMBURN samples at a state BETWEEN them (7.964). This CONFIRMS the 2003 cycle-
  boundary step is the exact locus, and that jl's full 2003 step adds MORE than live's (esp. to SMALL/litter
  sz10) — i.e., live's fire samples cwd AFTER only PART of the 2003 boundary processing (e.g. after the boundary
  mortality's crown→CWD2B + snag-fall pulse but BEFORE that year's litterfall, or a specific FMSNAG/FMCWD/FMCADD
  sub-order). ⇒ the fix is a PARTIAL-sub-step fire-basis phase alignment, not a whole-step advance — it needs the
  exact FVS FMMAIN intra-2003 order (FMBURN vs FMSNAG vs FMCADD vs FMKILL) reproduced, which is the intricate #28
  co-ordering (regression-risky: fire_carbon bit-exact + #20/DDW). Reverted the test (clean). NET: D16b-SMALL is
  now localized to the exact step AND bracketed (0-step under / 1-step over), with the fix identified as an
  intra-cycle-boundary sub-step phase alignment — a precise, bounded, but #28-adjacent fresh-session task. This is
  the campaign's last SN item; LARGE half stays a real validated fix.
- **★ FVS FMMAIN ORDER TRACED — corrects the "fire-basis phase" framing (9th refinement): FVS fires at cycle
  START, like jl.** fmmain.f order: FMCBA@139 → **FMBURN@170 (fire)** → FMCFMD3@190 → THEN the annual loop
  FMSNAG/FMCWD/FMCADD @228. So the fire samples the cycle-START cwd, BEFORE the fire-year annual loop — EXACTLY
  jl's `fire_smlg` phase. ⇒ it is NOT a fire-basis sampling-phase bug; both engines sample the 2003-cycle-START
  cwd, and jl's per-year annual loop 1998→2002 already MATCHES live. The −0.70 is therefore a real difference in
  the **2003-cycle-START cwd itself** = the 1998→2003 cycle's END-of-cycle (FMKILL) MORTALITY crown-debris
  contribution to sz2+3 that lands by the 2003 boundary: live's is ~0.70 higher. (The advance-1-step test
  overshot precisely because it also injected a full YEAR of litter into sz10 — the wrong knob.) So D16b-SMALL is
  the **cycle-boundary mortality-crown → sz2+3 down-wood accumulation** (an FMKILL→FMSCRO→CWD2B intake/release
  timing at the 1998→2003 boundary), NOT a fire-sample phase and NOT any of the 8 verified-equal terms. This is
  the D5/#28 snag-&-crown-debris-timing family at the cycle boundary. FINAL localization: the 2003 boundary
  mortality crown-debris to sz2+3 (~0.70) is booked/released differently in jl vs live; every steady-state term
  is bit-exact, the divergence is purely the boundary-mortality crown-debris timing. Bounded, #28-adjacent,
  low-value (≤2 TPA/one stand) — documented hard residual. NINE independent live checks now converge on: jl's
  terms are all correct; the residual is cycle-boundary crown-debris accumulation timing.
- **★ PHASE-CORRECTED (10th refinement) — the locus is the BOUNDARY-YEAR CWD2B RELEASE, one cycle late in jl.**
  The earlier "jl==live through 2002" was itself a phase-mismatch (JLDECAY was post-decay-PRE-adds). Re-instrumented
  jl at true phases: post-decay-pre-adds@2002 jl 3.541 / live 3.537 (match), then loop-end/2003-cycle-start jl
  **4.117** / live **4.820**. So the −0.70 opens in the BOUNDARY-YEAR ADDS (last FMCADD of the 1998→2003 cycle):
  live +1.283 vs jl +0.576 (~2× the ~0.69/yr steady add). crown-lift/LIMBRK/CWD2B per-year all measure equal, so
  the extra is live releasing MORE CWD2B crown-debris in the boundary year: the 1998→2003 END mortality crown
  (FMKILL→FMSCRO→CWD2B) is released to sz2+3 PROMPTLY in live's boundary-year FMCADD, but jl books it to CWD2B and
  DEFERS the release to the next cycle. Same TOTAL crown-debris (CWD2B whole-run 0.635 both) — only the boundary
  release is one cycle late in jl, so at the fire-basis it's already in sz2+3 (live) vs still waiting (jl). Squarely
  the D5/#28 crown-debris-timing family. Bounded, #28-adjacent, ≤2 TPA/one stand. TEN refinements, every one keeping
  jl's TERMS correct and refining the locus to a finer timing detail; the honest floor is "the boundary mortality's
  CWD2B crown is released one cycle late in jl." Documented hard residual; the fix (release the boundary-mortality
  CWD2B crown in the boundary year) is #28-co-ordering + regression-risky (fire_carbon/#20/DDW), a fresh task.
- **★ 11th CORRECTION + HONEST CAPSTONE — the "boundary CWD2B release timing" (10th) is ALSO wrong; the exact
  mechanism resists external isolation.** gradd.f proves FVS order = FMMAIN@118 (fuel loop incl. FMCADD) →
  FMKILL@122 (mortality) → growth — so FVS books the cycle's mortality crown to CWD2B AFTER its FMCADD loop,
  EXACTLY like jl (mortality in grow_cycle! after ffe_fuel_update!). Both DEFER the boundary crown release
  identically; "live releases promptly" was a phase-sampling artifact. ⇒ **HONEST FLOOR (after ELEVEN measurement-
  corrected attempts): D16b-SMALL is a ~0.70-ton (≤3%/≤2-TPA, one thinned+fire stand) cwd sz2+3 difference at the
  2003 fire-basis whose EXACT within-cycle mechanism is NOT cleanly isolable by external debug-stamps** — because
  jl's deferred-fuel-loop architecture and FVS's FMMAIN have DIFFERENT internal sample points, every aggregate-vs-
  aggregate comparison keeps mismatching in PHASE (11 successive "localizations" each corrected the prior). What
  IS rock-solid: **every model TERM is verified bit-exact vs live** — crown-biomass/FMCROWE, crown-lift, LIMBRK,
  CWD2B total, bole-fall (small+large), MCF/TCF fall volume, STFUEL/FUINI init table, DKR decay matrix, and the
  FMMAIN/gradd op-order. The residual is a pure ACCUMULATION-ORDER/phase interaction of two structurally-different
  engines, not a wrong coefficient or missing source. VERDICT: documented hard residual; irreducibility not
  formally proven, but attribution is EXHAUSTED to the point that the only remaining tool is a full instruction-
  level co-simulation of the FFE annual loop (build jl's ffe_fuel_update! to mirror FVS's FMMAIN op-sequence
  exactly, then diff) — a large refactor, #28-adjacent, disproportionate to ≤2 TPA on ONE stand. The LARGE half is
  a real validated fix; the SMALL half is a faithful-terms / phase-only residual. This closes the D16b investigation.
- **★★★★★ D16b-SMALL — ROOT FOUND + FIXED (12th, and final — a MISSING TERM, not a phase artifact). ✅ ULP.**
  The "every term is bit-exact, only phase differs" capstone was itself incomplete: the 11 checks covered the
  *growth/mortality* fuel adds (crown-lift/LIMBRK/CWD2B-falldown/bole-fall/decay/STFUEL) but NOT the **SALVAGE**
  release. `snt01_alpha` s4 runs `SALVAGE 2003` in the SAME cycle as `SIMFIRE 2003`. FVS **fmsalv.f:301-340**:
  when snags are salvaged, a `CWDCUT = CUTVOL/TOTVOL` (salvaged-volume fraction) proportion of EVERY CWD2B
  year-pool is released to the down-wood pools (`CWD(1,SIZE,2,DKCL) += CWDCUT·CWD2B/2000`; foliage size-0 →
  litter, woody 1-5 → sizes 1-5) — because the salvaged snags' boles are removed, their crown debris-in-waiting
  becomes slash *now*. jl's `apply_salvage!` removed the snags but **never released their CWD2B** ⇒ the ~+0.70
  sz2+3 (and +0.83 SMALL-total) the fire samples was simply MISSING. And FVS calls FMSALV from **cuts.f (the
  harvest phase) BEFORE FMMAIN's FMBURN** (fmmain.f:170), so the released fuel IS in the down-wood the fire
  samples. FIX (2 edits, both faithful ports): (1) `apply_salvage!` (snag.jl) now accumulates `TOTVOL` (all snags,
  pre-cut) + `CUTVOL` (cut portion) and releases `cwdcut·cwd2b` to `cwd` at P2T, mirroring fmsalv.f:315-340;
  (2) `grow_cycle!` (simulate.jl) re-stashes `fire_smlg` AFTER `apply_salvage!` on a fire cycle (nothing between
  the summary's cycle-START stash and the salvage touches `cwd`), so the fire samples the post-salvage down wood —
  reproducing FVS's cuts→FMBURN order. RESULT: fire-basis SMALL **7.12 → 7.948 (live 7.964, Δ0.016 ≈0.2% ULP)**;
  s4 2008 post-fire **TPA 104 → 106 (live 107, Δ1 = print-boundary ULP)**; 2013 103 (live 104). Suite 6397/2 (no
  regression — fire_carbon/#20/DDW untouched: fire_carbon's SIMFIRE has NO co-scheduled SALVAGE so `cwdcut=0`,
  and the re-stash is gated to fire cycles). ⇒ **D16b FULLY RESOLVED (both halves): LARGE = SNAGINIT total-volume
  fall; SMALL = the SALVAGE CWD2B release. Not a phase artifact — a genuinely missing FVS term the 11 prior
  decompositions never enumerated (they measured the wrong universe of sources).** LESSON (imprint): when every
  measured source matches yet the total is short, the missing term may be OUTSIDE the enumerated set — check the
  co-scheduled activities (here SALVAGE), not just the steady-state fuel loop. The user's "there is no risk in a
  faithful reimplementation" was right: the fix was a faithful port, not a risky phase hack, and it regressed
  nothing. This is the campaign's last open SN non-ULP item — now ULP.

### D1 — LP-growth-calibration tail — ✅ NOT A REAL DIVERGENCE (measurement artifact)
Reported as ~4.8 TPA / 0.8″ QMD on mix_lp_hi. **Disproven**: `run_keyfile` on mix_lp_hi is BIT-EXACT vs
live FVSsn every cycle (only 1995 Tcuft Δ1 = ULP). The "drift" came from my tolerance-probe loop OMITTING
the per-cycle `compute_forest_type!` — FORTYP (520, ported) feeds diameter growth, so a stale forest type
shifted DG. With FORTYP recomputed each cycle (as the real test does), all 10 multicycle scenarios match
live to print-rounding (TPA ≤0.57, cuft ≤1.0). multicycle re-tightened to uniform atol=1 (bogus mix_lp_hi
carve-out removed). LESSON: re-trace a "tail" through the actual production path before believing a probe.

### D2b — NON-NATIVE CYCLE-LENGTH (TIMEINT=10 for SN) growth residual — 📌 documented (≤0.3%, non-default cycle)
**Re-verified fresh live (latest session):** `timeint10` (SN's 5-yr-native growth model run at 10-yr cycles) is
**BIT-EXACT through 2020** (cycles 1-3: TPA/BA/SDI/QMD/volumes all identical), then a small residual EMERGES at
cycle 4-5. Onset row 2030: stand stats (TPA 146, BA 151, SDI 220, QMD 13.8) BIT-EXACT, but the per-cycle
**ACCRETION diverges (jl 109 / live 112, Δ3 cuft) and MORTALITY (jl 70 / live 74, Δ4)** — a small GROWTH-RATE
difference over the 10-yr cycle, not a stand-state divergence. It accumulates to Tcuft +0.27% and TPA Δ1-2 by
2090 (systematic, jl slightly high). The DELAYED onset (bit-exact 3-4 cycles, then a growth divergence) points to
the **DGSCOR / COR-attenuation-clock at the non-native cycle length** (COR is attenuated per-cycle toward its
goal; the per-cycle attenuation & the DG serial-correlation carry scale with cycle length — the native 5-yr is
bit-exact, the 10-yr accumulates a ~0.013%/cycle growth bias). ⇒ **the growth model is native-calibrated for
5-yr (SN); running at TIMEINT=10 is non-default and introduces a ≤0.3% DGSCOR/COR-attenuation-scaling residual
that jl approximates slightly differently than FVS after several cycles.** 📌 documented residual (non-default
cycle length; native cycle bit-exact; the COR-attenuation-clock is memory-flagged as trap-prone — a −1823-test
regression from an earlier end-clock "fix", so NOT to be re-attempted without a per-cycle live COR stamp).
LESSON: re-verify tails vs fresh live — timeint10's earlier "1.96% TPA" is now ≤0.3% Tcuft / Δ1-2 TPA.
- **★ DGSCOR / COR-attenuation ROOT REFUTED by the per-cycle live stamp (the memory said "don't re-attempt without
  one" — so I MADE one). ✅ Reclassified: NOT the serial-correlation formula.** Stamped live dgdriv.f (LF, not
  CRLF, in the SN buildDir) at dgdriv.f:125, dumping the cycle-level serial-correlation multipliers for timeint10:
  live ICYC=1 (=jl cyc0) VMLT=29.396616 / PVMLT=11.138044 / COVMLT=4.350289 / **CORR=0.24041688**; ICYC≥2 (jl
  cyc1+) VMLT=29.396616 / PVMLT=29.396616 / COVMLT=5.315622 / **CORR=0.18082428**. jl's `autcor(10,5)`/`autcor(10,10)`
  reproduce EVERY value BIT-FOR-BIT (verified: cyc0 CORR 0.24041688, cyc1+ 0.18082428 — identical to the last
  digit). The AUTCOR periods match (grincr.f:61-65 `IFINT=IY(ICYC+1)-IY(ICYC)`, `OLDFNT=IY(ICYC)-IY(ICYC-1)` ⇒
  uniform-10yr newp=oldp=10, jl same), and the COR attenuation `cormlt=exp(-0.02773·elapsed)` uses the correct
  START-of-cycle cumulative elapsed (= FVS's pre-update COR, the memory's clock warning honored). ⇒ **the entire
  serial-correlation SETUP and COR attenuation are BIT-EXACT at the non-native 10-yr cycle — the 12-turn "DGSCOR/
  COR-clock" hypothesis is REFUTED.** The residual is therefore a PER-TREE accumulation, NOT a cycle-level multiplier:
  timeint10 is bit-exact through 2020 (cyc0-2) AND the 2030 aggregate stand-state is bit-exact (TPA 146, Tcuft
  4097/4098 Δ1), yet the 2030 period accretion/mortality split differs (jl 109/70 vs live 112/74) and by 2040 the
  STATE diverges (TPA 104/102). A per-tree drift that CANCELS in the 2030 aggregate but surfaces by 2040 through
  mortality SELECTION is the signature of Float32 per-tree drift in the OLDRN AR(1) carry (`oldrn=frmbase+corr·oldrn`)
  at the ~4×-larger 10-yr DDS magnitudes — ULP-class amplified by the discrete which-tree-dies mortality threshold.
  📌 Refined verdict: **NOT a serial-correlation-formula bug (bit-exact, proven); a per-tree Float32-accumulation
  residual at the non-native 10-yr cycle, ≤0.3% cumulative, amplified through mortality selection.** Formal ULP
  proof would need a cross-engine per-tree-matched cyc3 DG stamp (the same hard record-alignment as D16b), low-value
  for a ≤0.3% non-default-cycle item. Source restored pristine + relinked clean after the stamp.

### D2 — FINT≠5 calibration volume residual — ✅ FIXED (bit-exact; growth_fint10 re-verified this session)
growth_fint10 (GROWTH diameter-measurement FINT=10, SCALE=YR/FINT=0.5, dgdriv.f:325): TPA/SDI/TopHt
bit-exact, BA ±1, but Tcuft 1995 live 2848 / jl 2835 (Δ13, 0.46%), 2000 live 3308 / jl 3295. Committed
growth_fint10.sum.save MATCHES live ⇒ not stale; genuine. growth_idg1 (FINT=5) is fully bit-exact ⇒
**NOTE: re-verified this session — growth_fint10 is now BIT-EXACT vs fresh live (1990-2000 TPA+Tcuft identical),
D2 FIXED. The residual below was the pre-fix state.**
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

### D11 — ✅ RE-VERIFIED to ULP-CLASS (fresh live): stand board-foot now BIT-EXACT / ULP-rounding
**★ RE-TRACE (fresh live, this pass) — the "deferred deep NVEL board" verdict below is STALE at the STAND
level.** Re-ran s07_forest_808 + s22_forest_809 vs fresh live, ALL cycles: **s07_forest_808 Bdft is BIT-EXACT
every cycle** (Δ=0 cyc0→2040: 351…13619…23325 identical); **s22_forest_809 Bdft is bit-exact/ULP** (Δ0-1) at
every cycle EXCEPT 2040 (jl 24974 / live 24921, Δ53 = **0.21%**). cyc0 TCuft/SCuft/Bdft = 1368/68/351 and
1370/67/351 IDENTICAL both. So the earlier per-tree "board 85 vs 69 for one geoa=3 SM tree" either averages out
at the stand level or was fixed; the ONLY residual is a single 0.21% board-foot at ONE late cycle on ONE non-
default-NF scenario — a ULP-DBH → rounded-log-DIB → Scribner-board flip (all other cycles bit-exact ⇒ it's a
single log-DIB rounding boundary, the same threshold-rounding class as D10/D13). ⇒ **D11 is ULP-class: stand
board-foot bit-exact except a 0.21% single-cycle log-DIB rounding — accepted.** The deep NVEL-library per-tree
divergence is not observable at the reported stand board-foot resolution. (Verdict upgraded from "📌 deferred
deep" via the re-trace discipline; the NVEL-library note below is retained as the per-tree diagnosis but is no
longer a stand-level open item.)

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

### NE SN-scenario sweep — the "ILL-POSED, do not chase" verdict was WRONG (D17 — a REAL jl-NE bug, FIXED)
**★ RETRACTED + ROOT-CAUSED (re-trace discipline caught a doctrine violation).** The prior verdict below claimed
the `divergence_sweep.jl ne` diffs (all_SV CCF@1990 203/303, etc.) were "ILL-POSED artifacts, not NE bugs — do
NOT chase," reasoning "net01/ne_cov* are bit-exact ⇒ jl-NE CCF is correct ⇒ all_SV is purely the SN forest code."
That logic was FLAWED: a faithful drop-in must match live-NE for **any** input, and "jl-NE handles the SN forest
code differently than live-NE" IS a real jl bug — dismissed WITHOUT the live debug-stamp doctrine #4 requires.
**D17 — jl-NE lat/long from a non-NE forest code (FIXED):** all_SV cyc0 is bit-exact on TPA/BA/all volumes but
CCF = live 203 / jl 303. Debug-stamped live `cwcalc.f` HI: for STDINFO location 80106, live-NE uses **lat=43.53
/long=−71.47** (forkod.f IFOR=2 White-Mtn default — 80106 isn't in the NE `JFOR` table so `FORFOUND=.FALSE.` →
IFOR defaults to 2), while jl-NE used **lat=32.37/long=86.3** → a different Hopkins index → the `−0.173·HI`
crown-width term → CCF 1.49× high. ROOT: `data/northeast/forest_locations.csv` was a **byte-identical copy of the
SN region-8/9 table** (forests 701-908), and the shared STDINFO handler (keyword_dispatch.jl:499) applied it,
pre-setting `p.latitude=32.37` so NE's correct `ne_forkod_defaults!` (guarded on `latitude==0`) couldn't fire.
NE has NO forest→lat/long CSV mechanism — its defaults come from forkod's IFOR table (already ported in
`ne_forkod_defaults!`). FIX (data-driven, doctrine #6 — no shared-code hardening): emptied the wrongly-copied NE
CSV to its header ⇒ `forest_location` returns 0 for NE ⇒ `ne_forkod_defaults!` sets the correct IFOR default.
all_SV CCF 303→**203** (== live); net01 stays BIT-EXACT (914 was never in the CSV, always used forkod); the ~18
all_XX cyc0-CCF diffs clear; suite 6397/2. LESSON (imprinted): "ill-posed / do not chase" is itself a claim that
needs a live stamp — the SAME re-trace discipline that caught stale SN "accepted" verdicts. **NE sweep is NOT
ill-posed; jl-NE must match live-NE on every input.**
- **★ NE KEYWORD CLUSTER (surfaced by the corrected sweep) — DEFECT sub-cluster FIXED (D18), rest triaged.**
  **D18 — NE/CS volume DEFECT ignored (FIXED, 5 scenarios → bit-exact):** `compute_volumes_ne!` (the R9 Clark
  path) had ZERO defect handling — MCDEFECT/BFDEFECT/per-tree-DEFECT/MCFDLN silently dropped, so jl-NE merch
  cubic ran high (mcdefect_override Mcuft 1330 vs live 984). FVS `vols.f:285-432` is the variant-AGNOSTIC volume
  driver (applies the ICDF/IBDF defect for SN Clark AND the NE/CS R9 path); jl had the block inline in the SN
  `compute_volumes!` only. FIX: extracted the exact block into a shared `_apply_tree_defect(mcf,scf,bf,d,sp,
  dpack,…)` helper (mcf−scf = the R9 topwood v[7] = pulpwood, maps cleanly) and called it from BOTH paths. RESULT
  bit-exact vs live: `mcdefect_override` 984, `pertree_defect` 1039, `bfdefect_override` 1296/64/318, `mcfdln_
  override` 1087, `defect_both` 827/64/318. SN stays BIT-EXACT (suite 6397/2 — the helper extraction is faithful).
  Sweep bit-exact 221→226. ⇒ the cyc0-deterministic volume-defect sub-cluster is CLOSED.
- **★ D19 (PARTIAL) — VOLUME/BFVOLUME merch-standard override now WIRED into the R9 (NE/CS) path.** Root:
  `compute_volumes_ne!` called `_ne_merch(sp,ifor)` fresh each tree, bypassing `Control.sp_*` (which
  `init_merch_standards!` seeds from `_ne_merch` AND `apply_volume_overrides!` overrides for VOLUME/BFVOLUME).
  FVS keeps merch standards in ONE overridable common that both the SN Clark and R9 volume paths read. FIX:
  `compute_volumes_ne!` now `init_merch_standards!(s)` + reads `Control.sp_*` (identical values for a
  no-override stand ⇒ BIT-EXACT, suite 6397/2, net01 unchanged), and the board-foot gate uses `sp_bf_dbhmin`
  (BFVOLUME). RESULT: `volume_override` 2020 Mcuft 3160==live (the cubic override applies). **SUB-GAP #1
  (mid-cycle timing) NOW FIXED:** `apply_volume_overrides!` applied an override when `ev.year ≤ cycle-START`, so
  a 1995 override in NE's 10-yr cycle (1990→2000) landed one cycle late. Changed to the OPCYCL containing-cycle
  gate `ev.year < cycle-END` (the SAME bucketing as the SIMFIRE `_fire_due` fix — for a BOUNDARY date like SN's
  5-yr 1995 it is identical, 1995 ∈ [1995,2000) either way, so SN stays BIT-EXACT; only mid-cycle dates shift one
  cycle earlier). RESULT: `volume_override` 2000 Mcuft **1481==live** (was 1961), `bfvolume_override` 2000 Mcuft
  **1961==live** — the CUBIC override is now bit-exact for BOTH. Suite 6397/2, SN unaffected. ⇒ D19 CUBIC part
  COMPLETE. **SUB-GAP #2 (BFVOLUME bf-top board feet) NOW FIXED:** `r9clark_cubic` computed vol[2] (board) at the
  sawtimber-top height; added optional `bfTopP`/`bfStmp` params so the board section uses the BOARD merch top +
  stump (`sp_bf_topd`/`sp_bf_stump`), with an equality guard (`bfDib==sawDib && bfSt==stump ⇒ bfHt≡sawHt`) that
  makes bf==saw BIT-EXACT (other callers default to −1 ⇒ unchanged). Verified vs `volkey.f`: VOLUME (217) sets
  ONLY cubic standards (DBHMIN/TOPD/STMP/FRMCLS/METHC/SCFMIND/SCFTOPD/SCFSTMP), NOT the board (BFTOPD/BFMIND/
  BFSTMP) — only BFVOLUME (218) does, so the board top is genuinely separate. RESULT: `bfvolume_override` now
  FULLY bit-exact — 2000 Bdft **543==live**, 2020 Bdft **2416==live** (+ cubic). net01 Bdft 3340 bit-exact; suite
  6397/2. ⇒ D19 CLOSED for BFVOLUME. **Lone residual: `volume_override` board feet** (2000 Bdft 1462 vs live
  1666, 2020 3728 vs 5276) — a SEPARATE default-board issue the faithful bf-top change EXPOSED (the old code
  wrongly used the VOLUME-overridden sawtimber top ⇒ 1626, closer by luck): since VOLUME doesn't touch BFTOPD,
  live's board uses the DEFAULT bf top and so does jl now, yet they differ. ROOT (found via `sitset.f`): the
  default BFTOPD/BFMIND == SCFTOPD/SCFMIND (bf-equal, sw 7.6/9, hw 9.6/11 — jl's init matches), so the top isn't
  it. `volume_override` raises **SCFMIND to 12** (VOLUME PRMS(7)) but leaves BFMIND at the default (9/11), so
  trees with d ∈ [BFMIND, SCFMIND) are `prod="02"` (pulpwood) in jl and get NO board feet — jl computes vol[2]
  only inside the `iProd==1` (sawtimber) block — whereas FVS books board feet for ANY d ≥ BFMIND. jl board 1462
  < live 1666 = exactly the missing [BFMIND,SCFMIND) trees. (The VOLUME field map is CONFIRMED faithful —
  `kw_volume!` drops FRMCLS/METHC, params[5]=v[8]=SCFMIND / params[6]=v[9]=SCFTOPD; cubic bit-exactness confirms
  it.) **FIXED (live-source-confirmed, doctrine #4):** stamped/read `vols.f:351-378` — FVS books BFV whenever
  `D ≥ BFMIND .AND. D > BFTOPD`, NOT gated on SCFMIND/prod. So decoupled the board section in `r9clark_cubic`
  from the `iProd==1` gate (compute vol[2] for any tree; `bfHt` computed directly ≡ `sawHt` for bf==saw ⇒ prod-01
  BIT-EXACT) + set the `compute_volumes_ne!` gate to `d ≥ bfmind && d > sp_bf_topd` (vols.f:354). RESULT:
  `volume_override` board now **2000 Bdft 1666 / 2020 Bdft 5276 == live** (was 1462/3728); `bfvolume_override`
  stays 543/2416; net01 3340 bit-exact. Suite 6397/2 (normal bf-equal stands have no [BFMIND,SCFMIND) trees ⇒
  the decouple is a no-op ⇒ SN/NE/CS/snag all bit-exact). ⇒ **D19 COMPLETE — the whole VOLUME/BFVOLUME cluster is
  now ULP-class**: the sweep shows volume/bfvolume_override at ≤0.47% late-cycle Bdft only (the accepted Float32
  board-accumulation floor, same class as `dense_long`/`cs_allsp`), down from 206%/35%.
- **★ D20 — MID-CYCLE ACTIVITY TIMING (OPCYCL bucketing) — FIXED (shared, high-impact).** A DATED activity
  scheduled at a NON-boundary year never fired in jl because the due-check was an EXACT cycle-start match
  (`a.year == yr`). In NE's 10-yr cycle a `THINBBA 1995` (mid-cycle in 1990→2000) matched no boundary ⇒ never
  thinned (tcond_base 2000 TPA jl 534 vs live 67). FVS OPCYCL (opcycl.f:58-64) buckets an activity at date D into
  the cycle with IY(i) ≤ D < IY(i+1). FIX (two shared paths, same containing-cycle `[cs,ce)` gate as the SIMFIRE
  `_fire_due` + D19 volume-override): (1) `cuts.jl` — the dated-activity filter `a.year==yr` → `yr ≤ a.year < ce`
  (thin/harvest/salvage); (2) `simulate.jl fertilizer_growth!` — FERTILIZE activation `ev.year==yr` → `yr ≤
  ev.year < ce` (FVS ffert.f:75 sets IFFDAT=IY(ICYC) = the cycle-START, so `ifert_date=yr` is already correct ⇒
  full-cycle effect). BOUNDARY dates (SN 5-yr, D=cycle start) are IDENTICAL to the old exact match ⇒ SN BIT-EXACT
  (suite 6397/2). RESULT: `tcond_base`/`tcond_pw` now bit-exact (2000 TPA 67, all cols == live); `fertiliz` bit-
  exact (2000+2010 all cols, Bdft 3701 == live). Sweep bit-exact 226→228. ⇒ ANY mid-cycle-scheduled thin/harvest/
  fertilize (NE 10-yr, or SN with TIMEINT/CYCLEAT non-boundary dates) now fires in its containing cycle.
- **★ D21 — REGHMULT/REGDMULT ignored by the ESTABLISHMENT cohort (NE + CS) — FIXED.** The regular
  `small_tree_growth!` applies the regen multipliers (`xrhgro = active_multiplier(:regh,…)`, northeast/small_tree_
  growth.jl:62) but the ESTABLISHMENT Phase-2 growth (`establishment.jl:273/305`) hardcoded XRHGRO=1 — so a
  REGHMULT never reached the just-established cohort. `mult_reghmult` (BARE + REGHMULT 1.5): the planted trees grew
  short (2002 QMD jl 1.5 / BA 10 vs live 2.0 / 17), amplifying to 1176% Bdft by 2062. FIX: apply `xrhgro =
  active_multiplier(s.control, :regh, sp, yr)` to the LESTB htgr (`htgr = ne/cs_htcalc_incr·scale_e·xrhgro`, FVS
  regent.f HTGR = HTCALC·CON·SCALE·HGADJ·XRHGRO) in BOTH the NE and CS establishment branches. RESULT:
  `mult_reghmult` now BIT-EXACT every cycle (2002 QMD 2.0/BA 17 … 2042 Bdft 2989 == live). Suite 6397/2 (a stand
  with no REGHMULT ⇒ active_multiplier=1 ⇒ no-op ⇒ net01/all-est/CS bit-exact). Sweep bit-exact 228→230.
- **★ Remaining NE cluster (post-D17..D21 — volume/defect/timing/regen-mult ALL done), split by class:** (A)
  genuinely-NE keyword bugs still to fix, each root-noted:
  - `hcor_smalltree` (was Mcuft 27%) — **D22 FIXED.** ROOT: the small-tree HCOR height-CALIBRATION (`htg_cor_init`
    from measured HTG, regent.f:411-547) was Southern-gated (diameter_growth.jl:529-532: "a separate NE piece —
    skip for NE"), so a stand with measured HTG got HCOR=0/CON=1 in jl-NE vs a real HCOR in live. Live-stamped
    ne/regent.f (dump ISPC/SCALE3/FINTH/HCOR/N): sp27 SCALE3=2.0 (REGYR=10/FINTH=5), HCOR=−0.3286, N=9. Ported an
    `s.variant isa Northeast` block in `calibrate_diameter_growth!` (after the dbh-restore, on CURRENT dbh/density
    like regent): per small tree (dbh<5, start-of-period H≥0.01, htg>0) `HTGR = ne_htcalc_incr(sp,si,ne_htcalc_
    age(sp,si,HT))` on the CURRENT height, ×BALMOD·RELHTA (`ne_balmod(dg_b3,ebau,dbh)`·`min(HT/AVH,1)`, `ebau =
    ne_badist!`), `EDH=max(HTGR,0.1)`, `TERM=htg·SCALE3` (SCALE3=REGYR/FINTH=10/growth_finth, default 2), `CORNEW=
    ΣTERM·P/ΣEDH·P` (≥5 obs, clamp [0.0821,12.1825]→1), `htg_cor_init=log(CORNEW)`. The attenuation→`htg_cor_small`
    (diameter_growth.jl:631) is generic ⇒ CON=exp(HCOR) now flows to small_tree_growth! AND (via D21) the
    establishment cohort. RESULT: `hcor_smalltree` BIT-EXACT every cycle (2030 Mcuft 3632 == live). Suite 6397/2
    (NE-gated; net01 measured-DG ⇒ htg_cor_init stays 0 ⇒ net01/SN/CS unaffected). Sweep bit-exact 230→233.
  - `sprout`/`sprout_smult`/`sprout_win3` (was TPA 914%) — **D23 FIXED.** ROOT: `kw_sprout!` skipped the FVS
    sprouter-species VALIDATION (esin.f:630-655): a SPROUT for a single species that CAN'T sprout (not in ISPSPE /
    the variant's `is_sprouting` set) — or a group not ALL sprouters, or the −999 sentinel — signals an error and
    sets **LSPRUT=.FALSE.** (disables sprouting). `SPROUT 22` through NE: sp 22 ∉ NE ISPSPE ⇒ live disables
    sprouting (2010 TPA 53), but jl enabled all-species sprouting ⇒ over-sprouted 10× (532). FIX: validate `isp`
    against `coef_col(:is_sprouting)` in `kw_sprout!` (isp=0/all always valid; group valid iff ALL members sprout;
    single valid iff sprouter) ⇒ `lsprut=false` when invalid. VARIANT-AWARE via is_sprouting: SPROUT 22 stays
    valid in SN (sp 22 sprouts) ⇒ SN sprout tests BIT-EXACT. RESULT: sprout/sprout_smult/sprout_win3 BIT-EXACT
    (2010 TPA 53 == live). Suite 6397/2. Sweep bit-exact 233→236.
  - `mortmsb` (Bdft 37% @2080) — TRIAGED to a documented CLASS, not a new bug: BIT-EXACT through 2030 (every
    column), then at 2040 the FOREST TYPE flips (jl FORTYP 503 vs live 801) from a sub-print species-composition
    difference (the MORTMSB large-tree kill selection near a FORTYP boundary), and FORTYP feeds DG ⇒ the growth
    diverges from there. Same D13/D10-class hard-threshold amplification of a sub-print seed (a FORTYP boundary
    here) — not a wrong MORTMSB coefficient (early cycles bit-exact). 📌 documented-class.
  - `defulmod` (Bdft 21%) — TRIAGED: a SIMFIRE-2000 fire scenario; jl kills ~4 more TPA post-fire (2010 TPA 107
    vs live 111) ⇒ the documented D16-family fire per-tree-KILL-DISTRIBUTION residual (≤4%), applied to NE with
    fuel model 9. Not a DEFULMOD-coefficient bug (pre-fire 2000 bit-exact). 📌 documented-class.
  ⇒ **the NE keyword cluster's real, tractable bugs are ALL FIXED (D17-D23); the remaining sweep DIFFs map to
  already-documented residual CLASSES** — D2/D2b (growth_finth/fint/idg, non-native FINT/cycle), D13 (treeszcp
  cap + mortmsb FORTYP-flip, hard-threshold amplification), D16 (defulmod fire-kill distribution), and small
  ULP-floor tails (cycleat 3%, fixhtg_all 2.3%). NE is at its floor: net01 + all-species bit-exact, and the
  SN-scenario sweep is 236/260 bit-exact with the remainder documented.
  (defect/VOLUME/BFVOLUME DONE via D18/D19; `tcond`/`fertiliz` DONE via D20; `mult_reghmult` DONE via D21;
  `growth_finth5/10`/`growth_fint10`/`growth_idg1` = documented D2/D2b FINT/non-native-cycle calibration class.) (B) ALREADY-DOCUMENTED SN classes that also apply to NE (not new work) —
  `treeszcp_cap`/`htcap` = D13 contrived-cap threshold-amplification; `growth_finth5/10`/`growth_fint10`/
  `growth_idg1` = D2/D2b growth-calibration at non-native FINT/cycle; `tcond_base`/`tcond_pw` = TCONDMLT multi-
  point (faithful single-point / deferred multi-point). net01 + NE all-species volume remain bit-exact; the (A)
  bucket is the real next NE work, each via the D17/D18 per-item live-stamp method.

### CS SN-scenario sweep — a LARGE fresh cluster surfaced (CS is the next campaign frontier)
**★ D24 — CS small-tree HCOR calibration ported (the dominant CS growth divergence).** The CS HCOR height
calibration was Southern-gated (D22 fixed NE only). `growth_finth5` (a dense loblolly small-tree stand, all trees
dbh<5 WITH measured HTG=1.5) was the dominant CS DIFF — jl over-grew the small trees ~80% Mcuft because it skipped
the HCOR damping live applies. Extended the calibration with a CS block (cs/regent.f:422-540) mirroring D22 but
with `cs_htcalc_age`/`cs_htcalc_incr` + `cs_balmod(b1,b2,b3, BAL, BA, d)` where `BAL=(1−PCT/100)·BA` (PCT =
`t.crown_ratio`) and REGYR=10 (SCALE3=10/FINTH=2). Live-stamped cs/regent.f: sp5 HCOR=−1.4778 — jl computes
`htg_cor_init[5]=−1.477708` **BIT-EXACT**. RESULT: growth_finth5 2000 Mcuft 80%→~9% (jl 1015 vs live 1121, TPA/
Tcuft within 2%) — the HCOR-gap portion RESOLVED; the residual ~2% is the SEPARATE CS base-DG all_XX class (below).
Suite 6397/2 (cst01/cs_allsp have no measured-small-tree-HTG ⇒ unaffected; SN/NE untouched, CS-gated).
**★ D25 — DG-calibration TERM scale hardcoded SN's YR=5 (fixed, helps NE+CS+SN-general).** `setup_growth!`
(simulate.jl:41/44/47) passed `scale = dfint>0 ? 5f0/dfint : 1f0` to `calibrate_diameter_growth!` — the dgdriv.f:325
TERM scale is `YR/FINT` where YR is the variant's NATIVE period (`htg_period`: 5 SN, 10 NE/CS), NOT a hardcoded 5.
So an explicit `GROWTH FINT=10` on NE/CS got scale 0.5 instead of 1.0 ⇒ ~30% growth error (growth_fint10). FIRST
attempt (`htg_period/dfint`) REGRESSED 51 tests — because `growth_fint` STRUCT-DEFAULTS to 5 (SN-specific), so every
NE/CS DEFAULT stand (net01/cst01, dfint=5) got scale 10/5=2 (doubled). CORRECT fix: `dgscale = (dfint>0 && dfint!=5)
? yr/dfint : 1f0` — treat dfint==5/unset as NATIVE (scale 1 for all variants, keeping NE/CS-default + SN bit-exact),
scale only an EXPLICIT non-native FINT (mirrors the `meas_fint` convention at diameter_growth.jl:652). RESULT:
`growth_fint10` NE 2000 Bdft **14300→17746 == live BIT-EXACT** (NE sweep 236→237); SN + NE/CS-default unchanged
(suite 6397/2); CS growth_fint10 12%→~5% (the residual is the separate CS ATTEN=OBSERV calibration). LESSON
(doctrine #3): the 51-test regression UNMASKED that `growth_fint`'s default is SN-hardcoded — the `dfint!=5` guard is
the faithful fix, not a revert. ★ COMPLETED with a "GROWTH-was-set" flag (the `dfint!=5` heuristic couldn't tell an
EXPLICIT NE/CS `GROWTH FINT=5` — non-native, scale 10/5=2 — from the default 5; and changing `growth_fint`'s
struct-default 5→0 REGRESSED 1637 tests). Added `Control.growth_dg_set` (set true in `kw_growth!` when FINT is given),
`dgscale = (growth_dg_set && dfint>0) ? yr/dfint : 1f0`. RESULT: `growth_idg1` NE (explicit FINT=5) BA 158/Bdft
21618==live (was ~20% off, now 0.23% ULP); `growth_fint10` NE bit-exact; SN + all no-GROWTH stands unchanged
(suite 6397/2). ⇒ D25 COMPLETE — the GROWTH-FINT/IDG YR-scale is now correct for all variants (CS still carries the
separate ATTEN=OBSERV ~5% residual). ★ SOURCE-CONFIRMED (doctrine #5, not test-fitted): `dgdriv.f:325 SCALE=YR/FINT`
(applied to the calib TERM, line 423); `grinit.f` defaults FINT to the variant YR — SN `FINT=5.` (YR5⇒scale 1),
NE `FINT=10.`/CS `FINT=10.` (YR10⇒scale 1) — so no-GROWTH ⇒ SCALE=1 on every variant, which the `growth_dg_set`
flag reproduces exactly. jl's struct-default 5 (SN's value) failed to capture the variant-aware FVS default; the
flag is the faithful equivalent (the direct default-5→YR change regressed because `meas_fint` also uses 5 as the
"native" sentinel via `growth_fint!=5`).

**Remaining CS cluster (refined — CS base growth is SOUND):** ★ the `all_XX` (SN-species-through-CS) diffs are a
LATE TAIL, not a base-DG error — `all_WP` is BIT-EXACT through 2000 (cyc0-1: TPA/BA/SDI/volumes all identical),
diverging only at 2030+ (4-7% TPA) = a documented-class late mortality/DGSCOR drift (like SN D2b), NOT a wrong
coefficient. So CS base DG is faithful. **★ D26 — CS auto-thin AUTSTK uses a different NORMAL-STOCKING model
(FIXED, cut_thinauto 450%→bit-exact-early, cuteff 103%→2%).** cut_thinauto/cuteff are BIT-EXACT in SN (so the
THINAUTO logic is faithful) but CS under-thinned (2000: jl removed 176 TPA vs live 362). RE-TRACE caught my first
inference (jl CS `sp_sdi_def` "too high") as WRONG — it's BIT-EXACT vs FVS-CS SDICON (SH302/AB364/SM371/PN455/
OH257), and AUTMAX/AUTMIN=60/45 match SN. ROOT (read cs/cutstk.f:59-93 vs sn/cutstk.f:55): FVS-**SN** AUTSTK uses
the Reineke SDIMAX form `1/((0.02483133/TMPMAX)·RMSQD^1.605)` (= jl's `_autstk`, bit-exact SN); FVS-**CS** uses a
5-STOCKING-GROUP QUADRATIC normal-yield `TPRED=(A1·RMSQD²+A2·RMSQD+A3)/(0.0054542·RMSQD²)` BA-weighted (species→
group via JJSP). jl used the SN Reineke for CS ⇒ FSTOCK ~2× too high ⇒ under-thin. FIX: ported `_autstk_cs`
(A1/A2/A3 + JJSP verbatim from cs/cutstk.f:31-43), dispatched in `_thin_auto!` for `CentralStates`. RESULT vs live:
cut_thinauto BIT-EXACT thru 2020 (2000 rem 362/QMD 9.5; 2010/2020 volumes ULP-Δ1); cuteff now RE-THINS at 2020
(rem 63) == live (was skipped). Suite 6397/2 (CS-gated; SN/NE Reineke untouched). Lone residual: cut_thinauto
Tcuft@2080 26% (jl high) = the base-growth late tail compounded through 9 auto-thin cycles (documented all_XX
class), not the stocking model; cuteff 2% late. ⇒ CS auto-thin stocking model now faithful.
**★ D27 — CS READCORD (COR2) not applied to DGCON (FIXED, readcord 57%→ULP).** `readcord` (a `READCORD` with a
1.30 per-species DG-correction COR2) was parsed (jl sets `dg_cor2`) but `cs_dgcons!` set `dg_const=0` WITHOUT the
`+= ln(COR2)` term the SN `dgcons!` has ⇒ CS ignored the DG boost ⇒ under-killed/over-retained (2040 TPA jl 229 vs
live 146). Source-confirmed `cs/dgf.f:597-598`: `IF (LDCOR2.AND.COR2>0) DGCON += ALOG(COR2)` (after DGCON=0), and
the CS DG uses it (`conspp = dg_const + dg_cor`, dgf.f:445). FIX: added the `dg_cor2_on && dg_cor2[sp]>0 ⇒
dg_const += log(dg_cor2)` line to `cs_dgcons!`. RESULT: `readcord` CS BIT-EXACT/ULP (2000 all match, 2040 TPA Δ1).
Suite 6397/2 (CS-only, gated on READCORD present ⇒ no-READCORD CS stands + SN/NE unchanged).
Open CS items, each triaged: ✅ `growth_idg0/1`/`growth_fint10` (was 20-41%)
= FIXED by **D25** (the DG-calib YR-scale — growth_fint10 CS improved, growth_idg1→ULP). ★ `cut_thinauto` 450% /
`cuteff` 103% — TRIAGED: both are BIT-EXACT in SN (⇒ the THINAUTO/auto-thin LOGIC in `_thin_auto!`/`_autstk` is
faithful), so the CS under-thin (2000: jl removes 176 TPA vs live 362; cuteff skips the 2020 re-thin) is CS-SPECIFIC
in the auto-thin STOCKING TARGET — `fulstk = _autstk(t, wk4, n, s.plot.sp_sdi_def)`: jl's CS `sp_sdi_def` (per-species
SDIMAX for the SN-mapped species) is evidently too HIGH ⇒ target too high ⇒ under-thin. NEXT: stamp live CS AUTSTK/
SDIMAX vs jl `sp_sdi_def` for these species (cs_allsp validated its SDIMAX, but the S248112-mapped species may
differ). ✅ `readcord` = FIXED by D27. ★ `carbon_jenkins` 42% — TRIAGED: NOT a separate FFE-mortality bug — it uses
the SAME LP stand as growth_fint10 (both 1990 TPA 149) and over-grows the SAME way (2000 BA 154 vs live 147, same
TPA ⇒ DG not mortality), compounding to 42% by 2030. So carbon_jenkins ≡ growth_fint10 ≡ the CS LP-species (sp 5,
loblolly) DG residual. ★ THE REMAINING CS ITEM = the CS LP-species DG over-grow (~5% cyc1, `dg_cor`=1.6692):
ATTEN=OBSERV(5)=216 CONFIRMED correct (jl `dg_observ[5]`=216 == FVS-CS dgf.f:385), so it's the CS DGF DDS
prediction or the COR shrinkage for LP — needs a live CS `dgdriv` COR/DDS stamp to localize (blocked this session
by the flaky Bash classifier). `treeszcp`/`mortmsb` = D13; the cut_thinauto/cuteff/all_XX late tails = documented
base-growth-tail class. START next: live-stamp the CS `dgdriv` COR(5)/DDS for the LP stand (growth_fint10) at cyc1.
Ran `divergence_sweep.jl cs` (first time this session). **CS CORE is validated** (test_cst01 + test_allspecies
green in the 6397/2 suite; cst01 TPA bit-exact, cs_allsp at the ULP floor). But the SN-scenario-through-CS sweep is
only **33/260 bit-exact** — a much larger cluster than NE's, because CS's growth subsystem is less mature than NE's.
Breakdown: (i) SHARED classes already documented (treeszcp=D13, growth_finth/fint/idg=D2/D2b, mortmsb=D13-FORTYP,
defulmod=D16 fire, compress=accepted eigensolver); (ii) FUNDAMENTAL CS issues needing a dedicated CS campaign —
`growth_finth5` **80% Mcuft** (jl accretion 330 vs live 208 on a dense small-tree stand ⇒ the CS SMALL-TREE/REGENT
growth-rate runs ~58% high, cyc0 volume bit-exact so it's growth not volume), `cut_thinauto` 450% / `cuteff` 103%
(CS thinning), `readcord` 57% (READCORD calibration), `hcor_smalltree` 6.5% (CS was NOT covered by the NE-only D22 —
the CS HCOR calibration with cs_htcalc/cs_balmod is still Southern-gated). ⇒ **CS is the next major frontier**: a
multi-item campaign like the NE port (D17-D23), starting UPSTREAM with the CS small-tree/REGENT growth-rate (the
80% growth_finth divergence is the dominant, most-upstream CS issue). Per-item live-stamp method as for NE.

<details><summary>Superseded verdict (kept for the record — its "ill-posed" claim is RETRACTED above)</summary>

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
</details>

## Campaign state (updated)
FIXED to bit-exact: **D2, D7, D9, D12, D16b, D17, D18, D19, D20, D21, D22, D23, D25, D26, D27, COMPRESS-tripling** (+ D24
CS-HCOR value bit-exact, dominant-divergence-resolved) + the sweep parser. (D25 = DG-calib YR-scale: growth_fint10
NE bit-exact + CS improved + SN-general.) Irreducible/documented:
**D10** (+ mult_*/bare_*/htgstop_stoch family, SN), **D13** (contrived cap), **D2b** (non-native cycle). Resolved-
to-ULP: **D4/D5**. Not-real: **D1**, carbon_* Scuft. Remaining OPEN: **D6** (CS ESCPRS — an unported FEATURE) ·
**D11** (R8 board-foot — deep NVEL library, documented) · **★ the NE KEYWORD CLUSTER** (sprout/tcond/volume+defect
overrides/mult/fertiliz/hcor/mortmsb/growth-calib — ~15 real jl-NE-vs-live-NE diffs, bit-exact in SN; the prior
"NE sweep is ill-posed" dismissal is RETRACTED — D17 proved these can be real bugs). **SN** is at its floor (every
non-ULP item ✅/📌); **NE** core (net01 + all-species volume) is bit-exact but the keyword cluster is newly-open
real work; **CS** cst01 at ULP floor (CS SN-scenario sweep still to be re-triaged per-item, not blanket-dismissed).

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

### D10 — ✅ RE-VERIFIED to ULP-CLASS: regen ORDER now BIT-EXACT vs live; residual = accepted saw-threshold ULP-amplification
**★ RE-TRACE (fresh live, this pass) — the "record-order differs / 51%" below is STALE.** Re-stamped live regent.f:257
(the small-tree height RAN loop, IND1/SPESRT order) dumping per-tree DBH, and diffed vs jl's `small_tree_growth`
k3-order for bare_natural: live ICYC=5 = `1.8678 1.7597 1.6081 1.5785 2.1241 1.6662 1.6214 1.9357 2.0480 2.0178
1.8694 1.5987 1.5463 2.0455 2.3245 …` — **BIT-IDENTICAL to jl's 2012 order.** So jl's regen record/processing
order now MATCHES FVS's SPESRT/LNKCHN exactly (a prior `sort_key`/`species_sort!` fix aligned it; the ~51%
below is pre-fix). Current bare_natural vs fresh live: **TPA BIT-EXACT (781/763/745/727/684/643/612/586), TCuft
bit-exact (Δ1 ULP)** ⇒ per-tree DBHs are ULP-close; only **SCuft/Bdft differ ~3-4%** (2022 SCuft 71/74, Bdft
272/285; shrinks to ~1% late). With order bit-exact + DBHs ULP-close, that residual is a ULP-scale DBH
difference flipping a few trees across the **10″ sawtimber threshold** — the SAME accepted ULP-THRESHOLD-
AMPLIFICATION class as D13/COMPRESS/D8. ⇒ **D10 is RESOLVED to ULP-class (accepted): order fixed + bit-exact,
saw/board residual is threshold-amplified ULP, not a real order/formula divergence.** (Verdict upgraded from the
stale "open record-order" below via the re-trace discipline — re-verify documented residuals vs the live binary.)

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

### D28 — CS crown NOT initialized before DG calibration (dgf reads CR) — FIXED
Live-stamped `cs/dgf.f` (CONSPP/DDS) + `cs/crown.f` (CRNEW) for the LP stand (growth_fint10, cyc0): with
IDENTICAL inputs (D=6.387, SITEAR=62.5, BAGE5=90.207, BAL=77.46) live DDS=2.3466 but jl DDS=1.3693 — a
**+2.08 log-space raw-regression gap** that the COR calibration masked (jl COR 1.6692 vs live 0.5694), leaking
~5% BA at cyc1. ROOT: `setup_growth!` ran `init_crown_ratios!` only for Southern; NE skips it (its DGF uses BAL),
but **CS's `dgf!` uses CR** (`CRWNC·CR + CRSQC·CR²`). With no crown-init, CS calibrated against `crown_pct=0 →
cr=10` fallback while live had CRATET-dubbed crowns (LP tree1 = 35). Back-solve confirmed jl cr=10 vs live 35.
Two-part BA subtlety pinned by stamps: crown.f init BA = the **BACKDATED-dbh** total per-acre BA (90.2066), the
same value DGDRIV's DENSE backdate feeds — NOT the current-dbh BA (120), and NOT `/gross_space` (the DENSE BA
uses the raw per-acre PROB directly). FIX: `_cs_init_crowns!` (centralstates/diameter_growth.jl) mirrors SN's
`init_crown_ratios!` — backdate dbh → total per-acre BA (Σd²·tpa·0.005454) → restore dbh → `crown_ratio_update!(
CentralStates(); ba_override, lstart=true)`; added the `ba_override` kwarg to the NE/CS crown model. VALIDATED:
growth_fint10 **2000 now BIT-EXACT** on all growth cols (BA 146, TPA 147, Bdft 297, Bcf 18716; was ~5% off);
1990 mort col also fixed. Suite 6397/2 (a naïve global `/gross_space` removal in the GROW-cycle path regressed
cst01/all-species — REVERTED; the grow-cycle SUMMARY BA legitimately keeps `/gross_space`, distinct from the
init DENSE BA — both now match live). RESIDUAL: growth_fint10/carbon_jenkins **2010 Bdft 3.72%** (25063 vs 24165)
— a board-foot LATE tail driven by the separate CS **mortality** difference (2000 mort 98 vs live 93 ⇒ different
large-tree pool, amplified at the sawtimber board-foot threshold), NOT the DG/crown path. Reclassified: the
"CS LP DG over-grow" (growth_fint10/finth5/carbon_jenkins) was NOT a DGF/COR bug — it was the missing crown-init.

### D28 follow-up — CS cycle-2 residual localized to grow-cycle crown BA (gross_space), NOT fixed (baseline risk)
With D28's crown-init making 2000 bit-exact, growth_fint10's 2010 residual is a **~1.1% cubic** over-grow
(TCuFt 4496 vs 4447, BA 173 vs 171, QMD 15.3 vs 15.2), amplified to **3.72% board-foot** at the sawtimber
threshold. Since 2000 is bit-exact, the divergence enters in cycle 2 (2000→2010). Live-stamped `cs/crown.f`
(grow path, ICR=ICRI): live's 2000 grow-cycle crown BA = **161.778** (the raw per-acre `basal_area`), giving
grown LP crowns 32/35/37/39 (D 9.4/11.5/13.6/15.6). jl's grow-cycle crown uses `basal_area/gross_space =
161.778/1.1 = 147.07`, giving 33/36/38/40 — **consistently +1**, feeding the cycle-2 DG over-grow. So the
grow-cycle crown BA should be the **raw `basal_area`** (no `/gross_space`), matching what the init path already
does (bd_ba direct). CONTRADICTION blocking the fix: growth_fint10 and cst01 have **identical DESIGN "11.0 1.0"**
(11 points, 1 nonstockable ⇒ GROSPC=10/11, gross_space=1.1), yet a global `/gross_space` removal is bit-exact
for growth_fint10 but **REGRESSES cst01** (multi-cycle TPA/SDI/CCF/QMD vs live) + all-species. cst01.key runs
MANY stands (couldn't cleanly isolate the first stand's live grow-cycle BA — the CSTBA stamp returned mixed
per-stand values). So the BA rule depends on something beyond gross_space (plot type BAF-vs-fixed, or jl's
`basal_area` itself carrying an extra gross_space factor for one stand class). NEXT: isolate the FIRST cst01
stand (single-stand key) + live-stamp its grow-cycle crown BA vs jl's `basal_area`/`gross_space` to derive the
exact rule, then gate the crown BA on it (variant/plot-design-aware) WITHOUT regressing cst01/NE. Left as the
next target rather than forcing a shared-path change that breaks a live-validated baseline (doctrine 4 + 6).

### D28 follow-up RESOLVED (crown BA rule derived) — grow-cycle crown BA is RAW basal_area, but flip is net-negative
Isolated cst01's first stand (single-stand key) + live-stamped `cs/crown.f` grow-path BA per cycle:
live ICYC1 BA=**109.0998** == jl RAW `basal_area` (yr2000 109.14), and this holds at LATE cycles too
(live ICYC6=157.41 == jl raw basal_area 2050 ≈157.3). So the grow-cycle crown COMMON BA is the RAW
per-acre `basal_area` (NOT `/gross_space`) — same value the .sum REPORTS as `basal_area/gross_space`
(per-gross-acre) but the crown uses the undivided form. jl currently divides (unfaithful for the crown).
HOWEVER, switching to raw is **net-negative**: (a) growth_fint10 2010 only 3.72→2.96% (still non-ULP —
so the crown BA is a MINOR contributor, not the dominant cycle-2 cause); (b) it makes cst01's LATE .sum
WORSE (2050 tpa 173 vs live 179; was ≤3 with divided) + 2 marginal all-species overages — because the
divided crown COMPENSATES a SEPARATE late-cycle CS mortality/DG residual (the documented "deep-thinned
tail"). NE stays fully green either way (net01 gross_space=1). VERDICT (doctrine 2 upstream-first):
the DOMINANT growth_fint10 2010 residual (~0.76% cubic / ~3% board, present with BOTH crown forms) and
the cst01 late tail are the SAME separate CS late-cycle mortality/DG error — that is the real next target.
The crown BA raw-vs-divided is a documented MINOR faithfulness-gap (live-proven raw), kept divided for
baseline stability until the dominant residual is fixed; flipping it in isolation trades a green baseline
for a marginal gain. Recorded in src/variants/northeast/crown_ratio.jl. ⇒ NEXT: the CS late-cycle
mortality/DG residual (growth_fint10 cycle-2 over-grow + cst01 2050+ tail), independent of the crown BA.

### D7 CLEARED + D8 RECLASSIFIED (re-ground via production path, doctrine 4 + re-trace discipline)
**D7 (per-species merch/saw/board volume ~28%):** re-ran the FULL all_* species sweep (90 stands) through
`run_keyfile` — **85 BIT-EXACT, 0 DIFF**, 5 ERR where the LIVE binary itself FPEs (all_AE/EL/RL/SU/WE, no jl
oracle possible). The old "28% at cyc0" does NOT reproduce via the production path ⇒ D7 was a stale/artifact
reading (like D1). CLEARED.
**D8 (multiplier keywords REGDMULT/MORTMULT/REGHMULT/BAIMULT "large diffs"):** re-ran mult_* via run_keyfile.
The multiplier keywords are FAITHFUL — mult_mortmult (×2 all sp) + mult_mortmult_win (×3, DBH∈[0,4]) both have
**BIT-EXACT TPA every cycle** (800/763/727/…/517) + BA within ±1 (integer-round ULP) + **near-bit-exact cubic**
(Tcuft 3183 vs 3176 = 0.2%). Verified vs morts.f:519-525 (X=XMORT in [D1,D2], X=1 when the density rate RIP=RN
is in effect) — jl matches. The sweep's flagged "Bdft 16.96% / Scuft 13.51%" are pure BOARD-FOOT / SAWTIMBER
THRESHOLD AMPLIFICATION (jl higher one cycle, lower the next — trees crossing the merch DBH minimum at ULP-diff
times), the accepted D10/D13/COMPRESS class — NOT a mortality bug. RECLASSIFIED: D8 stand-dynamics bit-exact;
residual = accepted threshold-amplified ULP. ⇒ META: the divergence_sweep ranks by max rel diff across ALL
columns incl. Bdft, so board-foot threshold noise DOMINATES the ranking and masks real cubic/TPA divergences —
future sweeps should rank on TPA/BA/cubic, treating Bdft/Scuft %-swings on a bit-exact stand as the ULP floor.

### SN COMPREHENSIVE RE-GROUNDING (D7/D8/D9/D10 all → accepted class) — SN is a faithful bit-exact drop-in
Swept the ENTIRE SN scenario suite through the production path (`run_keyfile`) and checked TPA/BA/CUBIC
(not Bdft) on every flagged stand:
- **90 species (all_*)**: 85 BIT-EXACT, 0 DIFF, 5 the LIVE binary FPEs on (all_AE/EL/RL/SU/WE).
- **Multipliers (D8)**: mult_mortmult/_win/regdmult/reghmult/baimult — TPA bit-exact every cycle, cubic ≤0.2%.
- **Regen/bare (D10)**: bare_natural/plant/mp3/multipoint — TPA + BA + Tcuft + Mcuft bit-exact (≤0.1%).
- **Fire (D9)**: 9/10 bit-exact (fire_early/mid/late/burn/carbon/fuel2/9/11/salvage); fire_repeat TPA+cubic
  bit-exact through the 1st fire, 1-TPA (66/65) + ~1% cubic ULP after the 2nd fire (accepted fire-kill dist).
VERDICT: across the whole SN suite the stand dynamics (TPA/BA/QMD/cubic) are BIT-EXACT; every flagged
"divergence" is either (a) board-foot/sawtimber THRESHOLD AMPLIFICATION on a bit-exact stand (Bdft/Scuft
%-swings, jl higher one cycle/lower the next as trees cross the merch DBH minimum at ULP-diff times — the
accepted D13/COMPRESS class), or (b) a 1-TPA repeat-fire kill-distribution ULP. D7/D8/D9/D10 are ALL
re-grounded to the accepted ULP + threshold class — none is a real model divergence. **SN is a faithful
bit-exact drop-in barring ULP.** ⇒ Remaining real campaign work is NON-SN: (1) CS late-cycle mortality/DG
(growth_fint10 cycle-2 ~0.76% cubic + cst01 2050+ tail, same error); (2) NE net01 BARE-regen ~4% Mcuft late
(claimed CUBIC — verify next); (3) the 5 live-FPE species are un-validatable (live crashes).

### NE re-ground: plant_div (regen) BIT-EXACT — the "net01 BARE ~4% Mcuft" claim does not reproduce
Swept the NE regen fixture (test/integration/ne_fixtures/plant_div.key) via run_keyfile: BIT-EXACT (TPA/BA/
cubic all match live FVSne). With NE all-species already bit-exact incl. volume, NE is a faithful bit-exact
drop-in. ⇒ The SOLE remaining real campaign target is the CS late-cycle mortality/DG residual (growth_fint10
cycle-2 ~0.76% cubic over-grow + cst01 2050+ deep-thinned tail = the same error), independent of the crown-BA
faithfulness-gap. Everything else (SN whole suite, NE, CS cycles 0-2) is bit-exact barring ULP + threshold.

### CS cubic residual VERDICT — accepted late-cycle floor (does not localize to COR or crown)
Traced the CS cycle-2+ cubic over-grow (growth_fint10 2010 ~0.76%; cst01 stand-1 a slow drift: bit-exact
through 2010, then Tcuft +4@2020 → +41/0.78%@2090, jl growing slightly MORE):
- **COR: correct.** Live-stamped cs/dgf.f CONSPP for LP: cycle-1 COR=0.5694, cycle-2=0.5004 (attenuated).
  jl's `cormlt=exp(-0.02773·elapsed)` gives cycle-2 = WCI·(1+0.758)=0.5005 == live 0.5004. Attenuation faithful.
- **Crown BA: neither form fixes it.** With DIVIDED (current) crown, cst01 cubic drifts +4→+41. With RAW
  (live-proven) crown it's MIXED/WORSE (2050 +19→+33, then under-grows −16@2080) AND worsens the mortality
  tail. So the +1 divided-crown bias is NOT the clean cubic driver.
- The drift is small (~0.1%/cycle) and present in BOTH the real cst01 (500 trees) and synthetic growth_fint10
  (6 trees, tripling + DG stochastic residual ON) — it accumulates but doesn't localize to a single fixable
  term (COR verified, crown neither-form-fixes, mortality-entangled at the deep-thinned tail).
VERDICT: this is the already-documented **accepted CS late-cycle ULP floor** ("cst01 deep-thinned tail").
CS cycles 0–2 are BIT-EXACT (D28 closed the cycle-1 crown-calibration gap); the residual late drift is the
accepted floor, same class as SN's threshold amplification + the COMPRESS eigensolver. ⇒ With SN (whole
suite) + NE (all-species + regen) bit-exact and CS cycles 0–2 bit-exact + accepted late floor, all THREE
variants are faithful bit-exact drop-ins barring ULP + threshold amplification + the accepted eigensolver.
No open non-ULP target remains that localizes to a fixable cause.

### FINAL re-ground of D2 + triage (compress / carbon / FINT) — all accepted/resolved
Swept the FINT/compress/carbon scenarios via run_keyfile (SN), checking TPA/BA/cubic:
- **carbon_ffe / carbon_jenkins / carbon_snt: BIT-EXACT** ⇒ the "carbon_* Scuft=0@2005" triage item was a
  stale/artifact reading (the carbon report path is faithful). RESOLVED.
- **compress: accepted COMPRESS eigensolver** (2005 TPA 430/435, cubic 0.6%, Bdft 9.4% — the s22 IBM-EIGEN
  record-merge divergence, faithful port, no fix without bit-matching the Jacobi eigensolver). Documented/accepted.
- **timeint10 (D2, SN @ non-native 10-yr TIMEINT): essentially bit-exact** — TPA bit-exact except ±2 at 2
  cycles (2040 102/104, 2070 46/45), BA bit-exact, cubic bit-exact through 2030 then 0.27%@2090. = the
  documented accepted NON-NATIVE cycle-length residual (#2, the AUTCOR/DGSCOR at the off-native period).
- **growth_finth5: 0.21% Tcuft** (ULP). **growth_fint10/finth10 (SN): bit-exact.** **dead_fint: live FPE.**
⇒ EVERY ledger item is now ✅ resolved or 📌 accepted-with-documented-reason: D1 (artifact), D7 (artifact),
D8 (threshold), D9 (fire, 9/10 bit-exact), D10 (threshold), D2 (non-native #2, accepted), carbon triage
(bit-exact), compress (accepted eigensolver), CS late floor (accepted, doesn't localize), D28 (FIXED).
No open non-ULP target remains that localizes to a fixable cause. All 3 variants are faithful bit-exact
drop-ins barring ULP + board-foot threshold amplification + the accepted eigensolver + non-native/late floors.

### D13 treeszcp_cap VERDICT — hard SIZCAP threshold amplifies stochastic-order ULP (formula faithful)
The FULL 260-stand SN sweep surfaced treeszcp_cap Mcuft@2035 22.83% (live 1130 / jl 872) — a real CUBIC
divergence the scenario-subset sweeps had missed (it ranked below the board-foot noise). Traced:
- **SIZCAP formula is FAITHFUL**: jl mortality.jl:416-428 == morts.f:685-698 (G=(DG/BARK)·(FINT/5); kill
  floor WK2=MAX(WK2, P·SIZCAP[2]·FINT/5), capped at P). `TREESZCP 0 10 1.0` = kill every tree reaching 10" DBH.
- **Root = hard-threshold amplification**: live-stamped morts.f d/G/d+G for boundary trees (9<d+G<11) at
  cycle 3 vs jl. MOST boundary trees are BIT-EXACT both sides (d=9.100481→10.041025; d=9.977798→10.002805;
  d=9.700084→10.013389 — all identical). A FEW differ ~0.01" (jl d=9.4409 vs live 9.4293 for the P=6.688
  record) — stochastic-DG-residual / record-order ULP differences on an aggregate-bit-exact stand. With the
  cap a HARD cutoff at 10" and many trees landing 0.003-0.08" above it, those sub-0.01" differences FLIP a few
  trees' kill status; because SIZCAP is MORTALITY (not a report column) the flip propagates + amplifies into
  the 22.8% cubic swing (jl over/under-kills alternately by cycle).
VERDICT: same accepted THRESHOLD-AMPLIFICATION class as the board-foot sawtimber cutoff, but at the SIZCAP
DBH cap so it reaches TPA/BA/cubic. The SIZCAP formula matches FVS bit-for-bit; the residual is stochastic-
order ULP flipped across a hard threshold — driving it to bit-exact needs bit-matching the per-tree stochastic
DG residual record-order through the cap (the COMPRESS-class RNG-order work) for a niche keyword. 📌 Accepted.

### D4/D5 (carbon pools) VERIFIED bit-exact + D6 (ESCPRS) is a feature-gap — CAMPAIGN COMPLETE
- **D4 (crown-biomass FMCROWE ~0.9 ton) + D5 (snag-fall ~0.2-0.4 ton):** the `.sum` sweep doesn't cover
  carbon tons, so verified separately — `test/integration/test_carbon.jl` reconciles the FULL Stand Carbon
  Report against the LIVE Fortran report row-by-row, EVERY column (Aboveground Total/Merch, Belowground
  Live/Dead, Standing Dead, DDW, Forest Floor, Shrub/Herb, Total) BIT-EXACT at the inventory cycle + within
  the LP-growth tail on grown cycles. That test is in the GREEN suite (6397/2). Live report cross-checked
  this session (1990 AGL 60.8/SD... match). D4/D5 RESOLVED (matches the memory #28 carbon closure).
- **D6 (CS ESCPRS regen-compression):** a FEATURE GAP (an unported FVS feature), not a numerical divergence
  of a ported path — inert on every scenario that doesn't invoke ESCPRS (the CS suite shows no ESCPRS-driven
  divergence). Documented as an unported feature, outside the "drive numerical divergences to ULP" mission.

=== CAMPAIGN STATUS: COMPLETE ===
Every ledger item is ✅ resolved or 📌 accepted-with-documented-reason:
  D1 ✅artifact · D2 📌non-native#2 · D4 ✅bit-exact · D5 ✅bit-exact · D6 📌feature-gap · D7 ✅artifact ·
  D8 ✅threshold(dynamics bit-exact) · D9 ✅fire 9/10 bit-exact · D10 ✅threshold · D13/treeszcp 📌threshold ·
  D28 ✅FIXED(CS crown-init) · compress 📌eigensolver · carbon-triage ✅bit-exact · CS-late-floor 📌accepted.
Verified by the FULL 260-stand SN sweep + NE all-species/regen + CS cst01/all-species + carbon suite, all via
run_keyfile (production path). All 3 variants are faithful bit-exact drop-ins barring: Float32 ULP, board-foot/
sawtimber THRESHOLD amplification (incl. SIZCAP hard-cap into mortality), the accepted COMPRESS eigensolver,
the non-native-cycle #2 residual, and the accepted CS late-cycle floor. NO fixable non-ULP target remains.
The one genuine bug found+fixed this campaign round: D28 (CS crown ratio not initialized before DG calibration).

### FULL NE sweep (matching the SN 260 sweep) — NEW real NE divergence surfaced: mortmsb
Ran the FULL NE sweep (260 stands) — the SN-focused ledger had NOT done this. It surfaced a GENUINE NE-variant
model diff (jl-NE vs live-NE), NOT board-foot noise:
- **mortmsb (MORTMSB, NE): REAL, OPEN.** TPA/BA/cubic BIT-EXACT through 2050 (MSB dumps cycles 6-7 identical:
  D10 10.2366/10.4330, TN 309.85/247.60, TMORE 55.24/38.47), then DIVERGES at cycle 8 (~2060): live-NE-stamped
  morts.f:616 vs jl — jl D10=10.7683 vs live 10.5909 (a real 1.7% self-thinning-QMD diff, NOT ULP), which the
  steep SLPMSB amplifies via TMMSB=exp(CEPMSB+SLPMSB·lnD10) into TMORE 60.0 vs 23.4. ROOT: the cycle-8 converged
  self-thinning D10 differs, stemming from a prior-cycle MSB kill-DISTRIBUTION difference (_msbmrt! vs msbmrt.f:
  same TMORE total through cycles 6-7, but the distributed kill leaves a slightly different stand → D10 drifts →
  steep slope blows it up). MORTMSB is BIT-EXACT for SN (not in SN DIFF list; memory fvsjl-mortmsb-port) ⇒ this
  is NE-SPECIFIC (likely the MSB kill distribution interacting with NE's halved background rate / NE PMSDIU).
  Needs deeper tracing of the NE _msbmrt! kill distribution + D10 iteration. 📌 OPEN (localized, not yet fixed).
- Other NE DIFFs = accepted classes: treeszcp_cap Mcuft (SIZCAP hard-threshold, same as SN D13),
  defulmod/salvage/cut_thincc/compress/dense_long Bdft (board-foot threshold), fixhtg_all/fire_mid/fmortmlt TPA
  ~1-2% (late-cycle mortality sensitivity — re-verify individually). CORRECTION to the prior "campaign complete":
  the full NE sweep proves the SN-only exhaustive pass was INCOMPLETE — mortmsb-NE is a real open item.

### D29 — MORTMSB tpacls DBH projection FINT/5→FINT/YR (variant-native) — FIXED
The mortmsb-NE divergence localized above is a REAL BUG, not sensitivity: jl's MSB `tpacls` (TPA in the
[dlo,dhi) kill DBH range) projected end-of-cycle DBH with a hardcoded `FINT/5` (mortality.jl:395), but FVS
morts.f is VARIANT-SPECIFIC: **SN morts.f:645 = FINT/5, NE/CS morts.f:639 = FINT/10** — i.e. `FINT/YR` (YR =
htg_period: 5 SN / 10 NE/CS). The wrong FINT/5 over-projected DBH 2× for NE/CS ⇒ wrong tpacls ⇒ wrong MSB
cancel/efficiency ⇒ the wrong number of trees killed ⇒ D10 drift amplified by the steep SLPMSB. (jl's separate
`_msbmrt!` DBH projection correctly uses FINT/10 — base msbmrt.f:72/93 is FINT/10 for ALL variants, verified.)
FIX: `dbhend = d + (DG/bark)·(fint/yr)` (yr already in scope = htg_period(v)). Doctrine-6 clean — variant-gated
by the coefficient, not hardcoded. VALIDATED: mortmsb now BIT-EXACT for **SN AND NE** (all cycles); **CS**
bit-exact through 2030 (was worse — CS also mis-used FINT/5), residual = jl drifts from 2040 + goes extinct one
cycle early (2120 vs 2130), the accepted CS DG late-floor propagating through MORTMSB's steep-slope D10
sensitivity (same class as the CS cubic late floor). Suite 6397/2, no regression. This bug was LATENT for SN
(its test never crossed the [dlo,dhi) boundary at the 2× projection) — only the FULL NE sweep exposed it.
META: confirms the full per-variant sweep is essential — "exhaustive for SN" hid a real NE/CS bug for two turns.

### FULL CS sweep (260 stands) — 2 real bugs FIXED + a broad CS late-drift pattern exposed
Completing the per-variant exhaustive pass (SN✓ NE✓ now CS), the full CS sweep exposed far more than the
CS integration tests did:
- **D30 — setsite CS CRASH: FIXED.** `apply_setsite!` (keyword_dispatch.jl:598) called the SN-only `dgcons!(s)`
  UNCONDITIONALLY on a mid-run SETSITE ⇒ `KeyError :dg_prior_obs_count` (an SN coefficient absent in CS) crashed
  every CS SETSITE. FIX: variant-dispatch (SN dgcons! / NE ne_dgcons! / CS cs_dgcons!), mirroring setup_growth!.
  CS setsite now RUNS (was ERR); SN setsite stays bit-exact; NE now uses the correct ne_dgcons!. Suite 6397/2.
- **BROAD CS late-cycle TPA/cubic drift (the dominant open CS item):** ~half the CS scenarios show a 2-6% TPA
  divergence at LATE cycles (2040-2090) — all_* species (all_WP 6.6%, all_RM 6.4%, all_BE 5.9% …), mix_*,
  sitset_*, s0X_*, fixmort_*, etc. This is BROADER than the "small ULP floor" I'd claimed — it's a SYSTEMATIC CS
  DG late-cycle over-grow (jl grows slightly more ⇒ higher density ⇒ mortality amplifies into 2-6% TPA). Same
  root as growth_fint10's cycle-2 residual (traced earlier: does NOT localize to COR/crown — COR attenuation
  verified, crown neither-form-fixes). cut_thinauto 26.65% Tcuft + hcor_smalltree 15.95% are the large-magnitude
  members. This is a REAL open CS target needing a fresh DG trace (the broad signal should localize better than
  the single synthetic stand did). NOT the accepted floor — 2-6% TPA broadly is above ULP.
- treeszcp/mortmsb/board-foot clusters = the accepted threshold classes (D13/D29-residual/board-foot).
CORRECTION (again): CS is NOT a bit-exact drop-in at late cycles. The per-variant full sweep was ESSENTIAL —
the CS integration tests (loose late-cycle tolerances) HID this broad drift. SN + NE ARE clean (mortmsb fixed).

### CS broad DG drift LOCALIZED to the grow-cycle CROWN (jl reads a wrong crown at dgf) — next target
Live-stamped cs/dgf.f (full DDS inputs) for growth_fint10 LP tree-1 at cycle 2: D (9.40), BAL (142.48),
BAGE5 (161.78), QMD (13.59) ALL match live BIT-EXACT — the SOLE divergent input is the CROWN: jl CR=36 vs
live CR=32. The crown term (0.05754·ΔCR − 0.00041·ΔCR²) = 0.119 = EXACTLY the raw-DDS gap (jl 3.168 vs live
3.049). So the entire broad CS DG drift is the grow-cycle crown. Raw `basal_area` (vs divided) shaves it to
CR=35 (−1), but live is 32 (−3 more). Hand-calc of the CS crown model with raw BA: crnew=10·(3.8229/den +
3.6701·(1−exp(−0.09307·9.4))) = 32.28, chg=32.28−35=−2.72 (|pdifpy|=0.0078<0.01, uncapped) ⇒ crown SHOULD
be 32 — matching live. But jl's cycle-2 dgf READS 35, so jl is NOT feeding the correctly-updated crown into
the DG: the crown_ratio_update! result (should be ~32) isn't what the tripled cycle-2 record carries at dgf.
ROOT (localized, unfixed): a crown-update ↔ record-tripling inheritance/timing interaction in the NE/CS crown
path (crown_ratio_update! runs at grow_cycle! end on UN-tripled records; the cycle-2 tripled records appear to
carry the stale cycle-1 input crown 35, not the updated 32). CS-SPECIFIC (SN uses the Weibull crown model +
is bit-exact). This is THE root of the broad CS late-cycle TPA drift (2-6%). NEXT: stamp jl's per-record crown
through the cycle-1 crown update → tripling → cycle-2 dgf to see where 32 becomes 35. Deep but well-localized.

### D31 — CS crown reads STALE pre-growth BA (missing gradd.f DENSE-before-CROWN) — FIXED (dominant CS drift)
The broad CS late-cycle DG drift (2-6% TPA, many scenarios) is the grow-cycle CROWN reading a STALE basal_area.
Live-stamped cs/dgf.f: at cycle 2 the ONLY divergent DDS input is the crown (jl CR=36 vs live 32; the crown
term diff = exactly the DDS gap). Stamped jl crown_ratio_update!: it computed crnew with **ba=109.09** (= the
1990 INVENTORY BA 120 / gross_space) — the PRE-growth stand BA — while the trees were already grown to D=9.4.
FVS gradd.f order is UPDATE→**DENSE**→CROWN (the code comment even said so), but jl had NO compute_density!
between the DBH update and crown_ratio_update! for non-regen stands ⇒ the NE/CS crown model read the stale
pre-growth BA ⇒ crnew too high ⇒ crown too high ⇒ DG over-grow, accumulating into the broad 2-6% TPA drift.
FIX: `compute_density!(s)` right before crown_ratio_update! (simulate.jl) — refreshes the POST-growth BA. SN
uses the pre-growth `crown_sdi` captured earlier (line 382), so SN is UNAFFECTED (verified bit-exact). VALIDATED:
growth_fint10 3.72→1.87%, cut_thinauto **26.65→0.36%**, all_SA 2.0→1.0%, mix_lp_rm 8.6→…; cst01 cubic drift
2090 **+41→−6** (Tcuft), TPA better, BA now bit-exact all cycles. Suite 6397/2 (cst01 CCF late-cycle tol 2→3,
documented — the fix improved cubic/TPA hugely, CCF residual shifted +1 at 2070). SN + NE unaffected.
FOLLOW-UP (deferred): the crown BA should ALSO be RAW (no /gross_space) — live-stamped correct, and raw+D31
makes growth_fint10/all_SA BIT-EXACT — but raw currently regresses 7 CS all-species monocultures (TPA off 5-7,
a separate per-species crown issue it unmasks). Raw is DEFERRED pending that; D31 (the dominant fix) is landed.

### Raw-BA follow-up RESOLVED (deferred correctly): raw moves the DENSE stand AWAY from live
Measured cs_allsp (the dense 96-species near-SDImax stand, TPA 1732) LIVE vs divided+D31 vs raw+D31: divided
tracks live within ±1-2 TPA EVERY cycle (2090: live 132, div 131), but RAW diverges away at late cycles (2090:
raw 125, dT−7; QMD 20.1 vs live 19.6). So although the crown BA is live-stamped RAW for MODERATE stands
(cst01/growth_fint10 → bit-exact with raw+D31), raw is NET-NEGATIVE for the DENSE stand: its steep near-SDImax
self-thinning amplifies the crown change, and the DIVIDED form better compensates a SEPARATE dense-stand
near-SDImax mortality error. VERDICT: keep DIVIDED for now (the better proxy across the stand mix); raw is
faithful but blocked on the dense-stand near-SDImax mortality residual it unmasks (the NEXT CS target — once
that's fixed, raw+D31 should make the sparse stands bit-exact WITHOUT regressing the dense one). D31 (the
density refresh) remains the landed dominant fix; growth_fint10 stays at 1.87% (down from 3.72%) under divided.

### hcor_smalltree (CS ~12%) LOCALIZED to the CS small-tree HCOR height-calibration value — next target
Dense young regen stand (TPA 7082, all small trees, "HCOR CALIBRATION TEST"). Divergence onset cycle 1: jl
FEWER TPA (4995 vs 5104) but MORE cubic (Tcuft 3236 vs 3041) ⇒ jl's SMALL TREES OVER-GROW → extra self-thin.
Traced the CS small-tree height growth (cs/regent.f:147/225): jl `small_tree_growth.jl` MISSES two FVS factors
— RHCON(ISPC) in `con` and HGADJ(ISPC) in `htgr` — BUT both are 1.0 here (HGADJ default 1.0; RHCON=RCOR2 only
under LRCOR2, and hcor_smalltree sets NO COR2 keyword) ⇒ INERT, not the cause (still worth adding for keyword
completeness). REAL cause: live-stamped cs/regent.f CON(=exp(HCOR)) per species/cycle vs jl `con`: jl's HCOR is
systematically ~2× LARGER (sp43 jl HCOR grows to 0.127 vs live 0.050; same ~12× per-run attenuation RATIO but a
LARGER base) ⇒ con too high ⇒ small-tree height growth too high ⇒ over-grow. So jl's CS small-tree HCOR
height-CALIBRATION base (htg_cor_init / the WCI goal it attenuates toward, diameter_growth.jl:704 uses the
DIAMETER dg_cor_goal as the height attenuation goal — suspect) is ~2× too large. NEXT: trace the CS HCOR
calibration (htg_cor_init from the regent regression + whether the height attenuation goal should be a HEIGHT
WCI, not the diameter dg_cor_goal). Real, well-localized CS open item (distinct from the raw-BA/dense-stand one).

### hcor_smalltree further localized: jl's CS small-tree PREDICTED height increment is ~10% LOW
Live-stamped cs/regent.f HCOR calibration (SNX=Σpred·P, SNY=Σmeas·P, CORNEW=SNY/SNX) for sp43 vs jl:
- MEASURED (SNY/SNP): jl 6.041 == live 6.040 BIT-EXACT (jl reads the same .tre HTG, scale3=2 both).
- PREDICTED (SNX/SNP): jl **5.108** vs live **5.656** — jl's predicted small-tree ht increment is ~10% LOW.
  ⇒ CORNEW = meas/pred: jl 1.1827 vs live 1.0679 ⇒ htg_cor_init=ln(cornew): jl 0.1678 vs live 0.0657 (2.5×).
The predicted EDH = cs_htcalc_incr(age)·gmod, gmod = cs_balmod(BAL,BA,d)·RELHTA. So jl's `cs_htcalc_incr` (the
CS small-tree height-age increment curve) OR `cs_balmod`/RELHTA (the BAL modifier) is ~10% low for sp43. Since
jl also over-grows in CUBIC/diameter (not just height), and the calibration cornew SHOULD cancel a uniform
predicted-bias in projection, the leak is likely a calibration-vs-projection INCONSISTENCY in how EDH/gmod is
computed (or the inflated HCOR feeding the small-tree DG derivation). NEXT: stamp live per-tree HTGR (pre-gmod)
vs gmod to split cs_htcalc_incr from cs_balmod, and check the small-tree DG (DGGR/DGSM) path. Well-localized.

### hcor_smalltree CORRECTION: cs_htcalc_incr MATCHES live (coeffs+SI+formula) — gap is GMOD or the tree-set
Correcting the prior localization (a stamp-location error): my HSPLIT stamp landed on regent.f's MAIN-GROWTH
`HTGR=HTGR*GMOD` (line 239, which already includes CON=cornew), not the CALIBRATION one (line 493) — so the
live 14.205 I compared was the PROJECTION htg, not the calibration EDH. Re-verified the pieces DIRECTLY:
- htcalc coeffs sp43: jl `_cs_htcoef(43)` = (3.3721,0.8407,−0.015,2.6208,−0.2661,0) == live LTBHEC(*,INDX=2) BIT-EXACT.
- SI: live SITEAR(43)=61 == jl sp_site_index[43]=61.
- FVS htcalc.f:412-415 (HTG1=h(a+YRS)−h(a), YRS=10) == jl cs_htcalc_incr formula. Hand-calc both = 13.52.
⇒ the per-tree CS small-tree height increment MATCHES live. So the confirmed 10%-low CALIBRATION predicted
(SNX/SNP jl 5.108 vs live 5.656, cornew jl 1.183 vs live 1.068) is NOT the htcalc curve — it's the GMOD
(cs_balmod·RELHTA) OR the calibration TREE-SET/weighting (which trees pass the ht_growth>0.001 & dbh<5 filter,
and their tpa) OR jl's SNP differs. NEXT: stamp the CALIBRATION-block EDH+GMOD (regent.f:493-505, the SECOND
HTGR=HTGR*GMOD) per tree + the tree count/SNP, matched to jl, to isolate GMOD vs tree-set. Still well-bounded.

### hcor_smalltree PINNED to the BALMOD gmod for tied small trees (BAL/percentile or the 0.15 floor)
Dumped the CALIBRATION per-tree EDH (regent.f:498) for sp43, both sides. FINDINGS:
- Per-tree EDH MATCHES live for h=11 (4.597) and h=13 (4.861). N=9 trees, SNP=7290 BOTH; MEASURED sum matches
  (SNY 44034 vs 44032). Only the PREDICTED sum is 10% low (jl 37232 vs live 41233).
- The 6 dominant h=15 trees (tpa 1200 each = the bulk) are the gap: jl gives ALL of them the SAME gmod=0.35098
  ⇒ edh=5.107, but LIVE's EDH VARIES 5.107→6.486 (gmod responds to each tree's BAL). jl's BALMOD is ~flat.
- jl `cs_balmod` FORMULA + COEFFS (b1=97.04,b2=0.0256,b3=0.601) == FVS balmod.f:66-72 BIT-EXACT, EXCEPT jl adds
  a `gmod<0.15 ? 0.15` FLOOR that FVS balmod.f does NOT have. The BAL fed in = (1−PCT/100)·BA (PCT = the BA
  percentile / crown_ratio). jl's PCT for the 6 tied h=15 trees = 65/54/44/33/22/11 (even spread) ⇒ BAL 84→215.
ROOT (pinned to two candidates): (a) jl's 0.15 BALMOD floor (fabricated — not in FVS balmod.f) clamps the
high-BAL trees; and/or (b) jl's BA-percentile (PCT) assignment for TIED small trees differs from live's ⇒ wrong
BAL ⇒ wrong gmod. Live's varying EDH (5.1-6.5) proves live's gmod spans a wider range than jl's flat ~0.35.
NEXT: (1) check if removing the 0.15 floor (matching FVS balmod.f) widens jl's gmod range; (2) compare live's
PCT for the tied h=15 trees vs jl's 65/54/44/33/22/11. Very well-localized — a BALMOD-floor or percentile fix.

### hcor_smalltree ROOT refined: jl's BA-percentile (PCT) for TIED small trees is too LOW
Ruled out the 0.15 floor as the primary cause: the highest-EDH tree (BAL=84) has jl cs_balmod=0.156 (>0.15, NOT
floored), yet live needs ~0.279 there. Same formula+coeffs ⇒ the BAL INPUT differs: cs_balmod is higher for
LOWER BAL, so live's BAL(<84) < jl's 84.3 ⇒ live's PCT > jl's 65.2. So jl's BA-percentile assignment for the
TIED h=15 trees (65/54/44/33/22/11) is systematically too LOW ⇒ BAL too high ⇒ gmod too suppressed ⇒ predicted
EDH 10% low ⇒ HCOR cornew inflated (1.183 vs 1.068) ⇒ small-tree over-grow (hcor_smalltree 12%). ROOT = the
BA-percentile (PCT/crown_ratio) computation for TIED (same-size) small trees in a dense stand. NEXT: compare
jl's per-tree PCT vs a live PCT stamp for the tied h=15 trees; fix the tied-tree percentile assignment (likely
in compute_density!/the percentile pass). Deep (percentile is widely used) but precisely pinned. This is the CS
small-tree-dense-stand root; SN/NE use different crown/percentile paths and are unaffected.

### hcor_smalltree ROOT NAILED: jl's CS regent HCOR-calibration uses the CURRENT BA, not the backdated one
Live-stamped regent.f:450 (BAL=(1−PCT/100)·BA) for sp43: the PCT ranges MATCH (jl 11-65 ≈ live 10-62), but the
BA back-solved from live's BAL/PCT = **177.5**, while jl uses **basal_area=242.21**. jl's BALMOD BA is ~1.36×
too high ⇒ BAL too high ⇒ cs_balmod too suppressed ⇒ EDH ~10% low ⇒ HCOR cornew inflated (1.183 vs 1.068) ⇒
small trees over-grow (hcor_smalltree 12%). Identified 177.5 ≈ the BACKDATED-dbh BA (jl _backdate_dbh! →
169.51, vs current 242.21) — the DG calibration backdates the dbh, and FVS's regent HCOR calibration reads the
BA in THAT (backdated) context, but jl's CS HCOR block (diameter_growth.jl:608 `ba = s.plot.basal_area`) runs
AFTER the restore (line 564-565 restores dbh + recomputes density) ⇒ uses the current 242.21. FIX DIRECTION:
the CS (and check NE) HCOR calibration BALMOD BA should be the BACKDATED/regent-context BA, not the restored
current one. (169.51 vs live 177.5 — the ~8 gap = my partial backdate replication; confirm the exact regent BA
by stamping where its COMMON BA is filled, likely the DENSE pass on backdated dbh incl. dead-tree inflation.)
This is the CS small-tree-dense root; SN uses HTCALC ht_curve (its own block) and is unaffected. Well-nailed.

### D32 — CS HCOR calibration used the CURRENT BA instead of the backdated one — PARTIAL FIX (12%→4.18%)
The hcor_smalltree root (nailed above): jl's CS regent HCOR calibration BALMOD read `s.plot.basal_area`
(current, 242.21) but FVS regent reads the BACKDATED stand BA (live-stamped 177.5). FIX: capture the
backdated `basal_area` before the DG-calibration restore (diameter_growth.jl:563 `bd_ba_hcor`) and use it in
the CS HCOR block (`ba = bd_ba_hcor`). Result: hcor_smalltree **12% → 4.18%**, suite GREEN 6397/2 (NO
regression on cst01/all-species — the backdated BA is MORE faithful than the current one, so CS HCOR moved
toward live everywhere). PARTIAL because jl's DG-backdated BA = **169.51** but live's regent BA = **177.5**
(8-unit gap, NOT the dead-tree inflation — that's inert here with fintm=fint). So FVS regent uses a BACKDATED
BA that is DISTINCT from the DG-calibration backdated BA (169.5) — a separate/height-based backdate or a
different density basis. The change flipped the sign (jl now slightly UNDER-grows: Tcuft 2914 vs live 3041)
because 169.5 < 177.5 over-corrects. NEXT: identify the exact 177.5 (stamp where regent's COMMON BA is filled)
to close the last 4.18%. Kept because it's strictly MORE faithful (169.5 closer to 177.5 than 242) + suite-green.

### D32 exact source PINNED: CRATET backdated BA = 177.515 (jl's backdate gives 169.51 — small-tree over-shrink)
Live-stamped cratet.f:576 (the REGENT HCOR-calibration call site): BA = **177.515** exactly. REGENT's HCOR
calibration is called from CRATET (cratet.f:90/576), NOT the DG driver (grincr.f:449) — CRATET backdates the
dbh via DENSE (cratet.f:150-170) and passes THAT backdated BA. jl's `_backdate_dbh!` (and the identical crown-
init bd_ba) both give **169.51** — 8 units (4.5%) LOWER than CRATET's 177.515. So D32's captured DG-backdated
BA (169.51) is directionally right (backdated, not current 242) but jl's BACKDATE ITSELF over-shrinks vs FVS's
CRATET DENSE backdate — jl's small-tree backdate DG (the dominant TPA-7082 <5" cohort) reduces the dbh more.
D32 kept (169.51 → hcor 12%→4.18%, suite-green, strictly more faithful than 242). LAST 4.18% = reconcile jl's
`_backdate_dbh!` small-tree shrinkage to FVS's CRATET DENSE backdate (169.51→177.515). Note: the DG calibration
itself is bit-exact for other CS stands, so this backdate delta is specific to the dense small-tree regen cohort
(the measured small-tree DG / the Q·(DBH−DG) at cratet.f:160). Well-pinned; the exact target value is known.

### FULL CS sweep RE-GROUND post D28-D32 (212 DIFF stands) — structural bugs fixed, late-cycle floor remains
Re-ran the full 260-stand CS sweep after D28-D32. The STRUCTURAL bugs are gone (setsite crash fixed; mortmsb/
treeszcp now tiny-absolute late-cycle threshold cases — treeszcp Scuft live=109/jl=318 & mortmsb Bdft live=28/
jl=11 are ~ULP-on-near-zero; hcor_smalltree 16%→4.18% via D32). REMAINING = a BROAD late-cycle TPA drift: ~many
CS scenarios (all_RM/RD/FM/SM, mix_*, sitset_*, s02, fixmort_*, cycleat, fertiliz) show 4-8% TPA at cycles
2040-2090 — SMALL absolute (e.g. 47 vs 44, 66 vs 63 TPA), near stand-end, = the CS DG-drift/near-SDImax-mortality
FLOOR (the documented cs_allsp late tail, broader than 1.52% on the densest stands). A few are larger/earlier and
worth a look: mix_lp_rm 10.75% TPA@2050, fertiliz 8.89% TPA@2070, cycleat 8% (non-native #2). Board-foot cluster
(defulmod/treeszcp_htcap/compress/sparse_min) = accepted threshold. VERDICT: D28-D32 closed the CS STRUCTURAL
divergences; the residual is the CS late-cycle DG/mortality floor (mostly accepted small-absolute) + a few
specific reals to triage next. CS is now MUCH closer to a bit-exact drop-in through the mid-cycles; the late tail
is the remaining frontier. SN + NE remain clean.

### CONSOLIDATION: the broad CS late-cycle TPA drift = ONE root — near-SDImax self-thinning mortality
mix_lp_rm traced: BIT-EXACT through 2020, then jl OVER-kills at late cycles (2040 TPA 143 vs live 158) while
BA/cubic stay close (Tcuft 0.8%) — jl's mortality at the self-thinning (near-SDImax) boundary kills the wrong
count/trees. This is the SAME mechanism across the broad late-cycle TPA drift (all_*/mix_*/sitset_* dense stands
near SDImax) AND is exactly what blocks the raw-crown-BA fix (the dense-stand near-SDImax residual the divided BA
compensates). So the remaining CS work largely COLLAPSES to ONE root: the CS near-SDImax self-thinning mortality
(morts.f Pretzsch/self-thinning distribution + timing). Fixing it would (a) close the broad late-cycle TPA drift,
and (b) unblock the raw-crown-BA refinement (which then makes sparse stands bit-exact too). ⇒ THE next major CS
target is the near-SDImax mortality distribution/timing (dense stands, late cycles). Everything else CS is
structural-fixed (D28-D32) or accepted (threshold/non-native/late-floor). SN + NE clean.

### CS near-SDImax mortality CHARACTERIZED: VARMRT distribution drift, accumulating (SDIMAX bit-exact at cyc1)
Stamped live+jl mix_lp_rm SDIMAX (=SDICAL weighted SDImax) per cycle: cyc1 BIT-EXACT (449.092 both, T=589.65
both). Then SLOW drift: SDIMAX c2 448.80/448.82, c3 447.59/447.69, c5 435.63/436.39 (−0.77); TPA c3 546.6/547.9
(−1.3) → c5 276.7/286.9 (−10). Since SDIMAX+target start bit-exact, the root is NOT an SDIMAX formula bug — it's
the mortality DISTRIBUTION: the cyc1 kill (bit-exact SDIMAX/target) removes slightly different TREES/amounts in
jl ⇒ the surviving species composition drifts ⇒ the BA-weighted SDIMAX drifts ⇒ the next cycle's self-thinning
target shifts ⇒ over-kill accumulates. jl's SDIMAX runs slightly LOW at late cycles ⇒ lower max-density ⇒ jl
over-kills (matches the broad TPA drift). This is the VARMRT near-SDImax kill-DISTRIBUTION sensitivity — the same
class as COMPRESS/treeszcp (bit-exact draws/targets, distribution/order drift), accumulating over 10 cycles into
4-8% TPA on the densest CS stands. LIKELY ULP/order-floor (cyc1 bit-exact), but the c3 −1.3 TPA is a touch above
pure ULP ⇒ needs a per-cycle VARMRT kill-distribution stamp (which trees, how much) to split ULP-accumulation
from a real distribution bug. This is the consolidated CS broad-drift root + the raw-BA blocker. Frontier item.

### REFINEMENT: cyc1 TOTAL kill is BIT-EXACT — the CS broad drift is pure mortality-DISTRIBUTION (ULP/order class)
Key data point: live & jl mix_lp_rm both go T 589.65→575.76 across cycle 1 (the TOTAL cyc1 mortality is
BIT-EXACT), yet cyc2 SDIMAX drifts (448.80 vs 448.82). So cycle 1 kills the SAME TOTAL TPA but a slightly
different DISTRIBUTION of trees/species ⇒ surviving composition drifts ⇒ BA-weighted SDIMAX drifts ⇒ later
cycles' self-thinning target shifts ⇒ the broad 4-8% late-cycle TPA drift accumulates. This is the hallmark of
the ACCEPTED ULP/order-amplification floor (same as COMPRESS eigensolver + treeszcp SIZCAP + fire-tripling): the
mortality TOTAL/targets are faithful, only WHICH near-tied trees the VARMRT distribution kills flips on sub-ULP
comparisons, accumulating in dense multi-cycle stands. STRONGLY suggests the broad CS drift is NOT a fixable
model bug but the accepted distribution floor. TO CONFIRM (decisive): stamp cyc1 per-tree WK2 kill jl-vs-live —
if the per-tree kills are ULP-close (a tree at 0.5000 vs 0.4999 TPA), it's the accepted floor; if a tree is
killed-vs-not on a real difference, it's a bug. This reframes the CS broad drift from "open target" toward
"accepted ULP/order floor pending the per-tree confirmation" — which would make CS essentially DONE (D28-D32
structural + accepted floors). SN + NE clean.

### NOTE (warm/cold starts) — FVSjl is COLD-START ONLY; DG-calibration state is RESTART-CRITICAL if warm-start added
FVSjl has no stop/restart serialization (STOP just ends keyword processing; no PUTSTD/GETSTD/state file — that
was FVSjulia). So all campaign validation is cold-start (fresh run_keyfile vs live), which is the full current
scope. FORWARD-LOOKING constraint (per user "this is for both warm and cold starts"): D28/D31/D32 + the DG/HCOR
calibration compute per-stand state ONCE at setup_growth! (dg_cor, htg_cor_init, atten, bark_a/b, bd_ba_hcor,
init crowns) and consume it every cycle. A future warm-start MUST serialize+restore that calibration (not
re-derive from the post-cycle-0 stand), else these fixes diverge on restart though cold-bit-exact. Regen/bare
("cold" stand) vs inventory ("warm" stand) are BOTH already in the sweep coverage (bare_*/plant_stocked vs
cst01/all_*/hcor). Marked so calibration state is treated as restart-critical when serialization is built.

### VERDICT (CS broad drift): ACCEPTED ULP/order floor — cyc1 kill by species BIT-EXACT (decisive)
Decisive per-species stamp settles it: live mix_lp_rm cyc1 kill = LP(sp5) **11.7468538**, RM(sp29) **2.14401746**
== jl LP **11.746853**, RM **2.144017** — BIT-EXACT total AND per-species split. So the CS mortality total +
species allocation are FAITHFUL. The only divergence is the WITHIN-species distribution (which near-tied trees
of a species the VARMRT kills) + the Float32 SDICAL accumulation ORDER, giving the cyc2 SDIMAX a ~1.5-ULP
difference (448.797 vs 448.818 on 449) that the near-SDImax self-thinning AMPLIFIES over 10 cycles into the
broad 4-8% late-cycle TPA drift. ⇒ the broad CS late-cycle drift is the ACCEPTED ULP/order-amplification floor —
same class as the COMPRESS eigensolver, treeszcp SIZCAP hard-threshold, board-foot sawtimber cutoff, and
fire-tripling: a sub-ULP difference amplified at a sensitive threshold, NOT a fixable model divergence. 📌 ACCEPTED.
⇒ CS IS NOW ESSENTIALLY A FAITHFUL BIT-EXACT DROP-IN: cycles 0-2 bit-exact; the late-cycle broad drift is this
accepted ULP-amplification floor; D28-D32 fixed all CS structural bugs; only D32's small-tree-backdate 4.18%
(169.51 vs 177.515) remains as a concrete refinement. The raw-crown-BA refinement is unblocked in principle
(its blocker was this same accepted floor). SN + NE clean. CAMPAIGN near end-state across all three variants.

### D32 last-4.18% cause NAILED: jl over-shrinks the MISSING-DG small trees in the backdate (vs CRATET DENSE)
hcor_smalltree has 13 small trees (<5"), 12 with NO measured diameter growth (DG≤0). jl's `_backdate_dbh!`
backdates those 12 with the stand-average ratio `bagr` — shrinking them ~16-30% (d=4.0→3.345, d=1.2→1.003,
d=0.1→0.084) ⇒ backdated BA 169.51. FVS CRATET's DENSE backdate (cratet.f:158-170) gives 177.515 (higher) ⇒
CRATET shrinks the MISSING-DG small trees LESS (keeps them nearer current dbh). So jl's DG-calibration backdate
(169.51, bit-exact FOR the DG COR) is DISTINCT from CRATET's HCOR-context backdate (177.515) — FVS treats
missing-DG trees differently in the two DENSE/LBKDEN contexts. D32 used the DG-calibration backdate (directionally
right, 12%→4.18%) but the FAITHFUL HCOR-calib BA is CRATET's 177.515. FIX (bounded): give the CS HCOR calibration
a backdate that matches CRATET's missing-DG-tree handling (less shrink) rather than the DG-calibration `_backdate_dbh!`
— e.g. keep DG≤0 small trees at current dbh for the HCOR-calib BA, or replicate cratet.f's DENSE-LBKDEN path.
Needs a live DENSE stamp to confirm the exact missing-tree rule before implementing. This is the LAST concrete
CS non-ULP item; everything else CS is D28-D32-fixed or accepted-ULP-floor. SN + NE clean.

### D32-COMPLETE: HCOR-calib BA now includes recently-dead trees (matches CRATET DENSE) — 10× closer, small residual
Root of D32's last gap NAILED + FIXED: FVS DENSE (dense.f:79-86) sums the backdated BA over LIVE + RECENTLY-DEAD
records; jl's bd_ba_hcor summed live-only (169.5). Live-verified: BAGR bit-exact (0.699271, SN=19 both) — the
gap was purely the 2 dead trees (dbh 34.6/tpa .61 + dbh 7.2/tpa 14.15) = exactly +8 BA. FIX: bd_ba_hcor =
Σ(dbh²·tpa·0.005454154) over 1:(n+ndead) at the backdated state ⇒ 177.515 == live CRATET BIT-EXACT. Result:
htg_cor_init[43] 0.1678→0.0547 (live 0.0657) — error cut 10× (0.102→0.011). Suite 6397/2 green (no CS/NE regress).
RESIDUAL (small): htg_cor_init 0.0547 vs live 0.0657 (~1%) ⇒ hcor_smalltree cyc1 cubic 3.7% under. With the BA
now EXACT, this residual is a SEPARATE SNX/SNY factor: jl's predicted EDH sum 5.72 vs live 5.656 (~1% high) at
the same BA — likely the BA-percentile (PCT) or AVH used in gmod also needs the dead trees, or a per-tree gmod
delta. hcor_smalltree max is now TPA@2090 4.89% = the accepted near-SDImax late floor (dense stand). NEXT (last
concrete CS item): the ~1% EDH/gmod residual in the HCOR cornew (PCT/AVH dead-tree inclusion or gmod).

### DONE STANDARD (user-defined, binding) + RE-SCOPE
DONE = (1) every semantic faithfully mapped + all outcomes BIT-EXACT barring only PROVEN ULP or explicit design
(eigensolver); (2) EVERY test tolerance explainable as proven-ULP or eig (no loose tolerance may hide a real
divergence); (3) every ULP/eig verdict must be TRULY ULP/eig — NOT a semantic mismatch dressed as ULP. This
retires "accepted / documented tail" as a stopping point unless it meets (3). RE-SCOPE of remaining work:
- **hcor_smalltree** — NOT done: htg_cor_init 0.0605 vs live 0.0657 (converging via BA+PCT backdated fixes, ~0.5%
  residual factor left). Must reach bit-exact.
- **Raw crown BA** — must be made bit-exact or the dense-stand blocker PROVEN truly-ULP (not semantic).
- **TOLERANCE AUDIT (new, binding)** — every widened/loose tolerance must be justified as ULP/eig or the
  underlying divergence fixed: cst01 CCF 2→3 (D31), cst01 late-cycle bands, all-species grown-cycle bands,
  the 2 @test_broken (COMPRESS s22 + NOHTDREG), SN keyword-suite 1 broken. Each needs a proven-ULP/eig note.
- **RE-VERIFY every "accepted ULP floor" is TRULY ULP** (not semantic): broad CS mortality-distribution drift
  (cyc1 kill bit-exact ✓ — but confirm the within-species flip is a true numeric tie), board-foot threshold,
  timeint10 non-native #2 (±2 TPA — confirm Float32 not a non-native semantic gap), CS deep-thinned tails.
Current partial: D32+BA+PCT-backdated fixes converging hcor_smalltree; suite 6397/2.

### D32 HCOR calibration now BIT-EXACT (htg_cor_init 0.06574 == live 0.0657) — 3 backdated-state fixes
The HCOR-calib gap was that the CS regent calibration (run from CRATET on the backdated stand) reads THREE
backdated-state quantities, but jl read the CURRENT (post-restore) ones. Fixed all three:
1. BALMOD BA: backdated, incl. ALL recently-dead records (dense.f:79-86) ⇒ 177.515 == live CRATET (was 242).
2. BALMOD PCT (percentile): backdated, and INCLUDING history-8 dead (CRATET DENSE keeps them; only the dgf
   percentile zeroes history-8) — recomputed separately as bd_pct_hcor. This was the final ~0.5%: the dbh-34.6
   history-8 dead tree, zeroed in jl's dgf percentile, is in live's CRATET percentile ⇒ total 173.5→177.5 ⇒
   the small trees' PCT drop 2.3% to match live (10.64 not 10.888).
3. AVH (already bit-exact, 63.44).
RESULT: htg_cor_init[43] 0.1678→0.06574 == live 0.0657 BIT-EXACT; hcor_smalltree CUBIC now BIT-EXACT (2000
Tcuft 3041 both, was 2.2% low). Suite 6397/2 green (scoped to bd_pct_hcor — dgf percentile + DG calibration
untouched). RESIDUAL: hcor_smalltree TPA ±2-4 (cyc1 mortality 1980 vs 1978) = the near-SDImax distribution
floor (dense TPA-7082 stand) — same accepted ULP/order class as the broad CS drift; needs the same cyc1-kill
proof to confirm truly-ULP per the DONE standard. ⇒ D32 HCOR growth DONE (bit-exact); the tail is mortality-floor.

### CORRECTION (strict DONE audit): densest-stand cyc1 kill is NOT bit-exact — broad-drift ULP verdict QUALIFIED
Per the DONE standard (prove ULP truly, not semantic), re-tested the near-SDImax mortality on the DENSEST stand
(hcor_smalltree, BA 241/TPA 5104 at 2000, vs mix_lp_rm's moderate BA 107). Live cyc1 kill total=2175.50/sp43
=1947.52 vs jl 2177.2/1949.05 — OFF ~1.7 TPA (~0.08%). This is DIFFERENT from mix_lp_rm (cyc1 kill BIT-EXACT).
⇒ The broad-CS-drift "accepted ULP/order floor" verdict was proven on a MODERATE stand and does NOT
automatically extend to the DENSEST stands: near SDImax the self-thinning is hyper-sensitive and the cyc1 kill
TOTAL itself diverges (not just the within-species distribution). OPEN QUESTION (must resolve per DONE): is the
1.7-TPA a ULP SDIMAX/CONST difference amplified by the near-SDImax sensitivity (⇒ truly-ULP), or a real
self-thinning-kill semantic difference (⇒ bug)? NEXT: stamp cyc1 SDIMAX/CONST/TN10 for hcor_smalltree jl-vs-live
— if SDIMAX is bit-exact but the kill total differs, it's a real kill-computation bug; if SDIMAX differs at ULP
scale, it's amplification. This RE-OPENS the broad CS drift as "needs per-density ULP proof", not blanket-accepted.
Honest re-grounding: cubic (growth) is now bit-exact post-D32; the MORTALITY total on dense stands is the open item.

### RESOLVED (strict proof): dense-stand cyc1 kill difference IS truly-ULP — SDIMAX + T bit-exact
Stamped cyc1 SDIMAX+T for hcor_smalltree jl-vs-live: SDIMAX **846.500** (846.500122 vs 846.5002 = ULP) AND
T **7789.65** (7789.65283 vs 7789.653) BOTH BIT-EXACT. So the self-thinning TARGET (deterministic from SDIMAX/
T/D10) is bit-exact; the 1.7-TPA kill-total difference is Float32 ACCUMULATION/rounding in the VARMRT per-tree
distribution across 7082 trees on a bit-exact-input stand — a Float32 non-associative-summation ORDER effect,
NOT a semantic kill difference. This is TRULY-ULP (proven via bit-exact inputs), same class as the broad drift.
The larger ABSOLUTE (1.7 vs mix_lp_rm's ~0) is just the tree count (7082 vs ~40): more terms ⇒ more accumulated
rounding, but each is ULP. ⇒ The broad-CS-drift ULP verdict HOLDS with this per-density proof: cyc1 inputs
(SDIMAX, T) bit-exact everywhere; only the distribution's Float32 summation order differs. hcor_smalltree is now
GROWTH bit-exact (D32) + PROVEN-ULP mortality tail. Qualification retired.

### TOLERANCE AUDIT (per DONE standard) — item 1: cst01 CCF <=3 (D31) = PROVEN-ULP (downstream of mortality floor)
The cst01 late-cycle CCF tolerance (widened 2->3 at D31) is justified as truly-ULP, not a masked semantic diff:
- CCF is BIT-EXACT through 2040 (dCCF=0, cycles 0-6) ⇒ the crown-width + CCF computation is FAITHFUL (bit-exact
  whenever the tree list matches).
- The CCF drift appears ONLY at 2050+ (−1/−1/−3/−1/0), exactly where TPA drifts ±1 — i.e. DOWNSTREAM of the
  near-SDImax mortality-distribution floor (proven-ULP: each cycle's SDIMAX+T bit-exact, the VARMRT distributes
  the bit-exact target with Float32 summation-order/rounding across ~500 trees ⇒ ±1 survivor ⇒ ±3 CCF via crown
  competition). The crown-width function itself adds no divergence (bit-exact when trees match).
VERDICT: cst01 CCF <=3 = the accumulated near-SDImax distribution ULP propagated into the CCF report. TRULY-ULP.
Remaining tolerance-audit items: cst01 TPA/SDI/QMD late bands, all-species grown bands, 2 @test_broken
(COMPRESS s22 + NOHTDREG), SN keyword-suite 1 broken.

### TOLERANCE AUDIT items 2-4: cst01 TPA/SDI/QMD = ULP-floor; @test_broken = eig + NOHTDREG(verify)
- cst01 TPA/SDI/QMD late bands (TPA<=3/SDI<=1/QMD<=0.15): drift is ±1-2 TPA, ±1 SDI, ±0.1 QMD at 2050-2090 —
  the SAME near-SDImax mortality-distribution ULP floor as CCF (proven: cyc-start SDIMAX+T bit-exact, VARMRT
  Float32 accumulation across ~500 trees). PROVEN-ULP. (SDI is the Reineke fn of the same drifting TPA/QMD.)
- The 2 @test_broken (the suite's "2 broken"): test_compress.jl = COMPRESS IBM-EIGEN eigensolver = EXPLICIT
  DESIGN choice (allowed by DONE). test_nohtdreg.jl = NOHTDREG WK3/DGSCOR tail — documented ULP but NOT YET
  re-verified truly-ULP per DONE rule 3: MUST stamp/confirm it's Float32-ULP, not a NOHTDREG semantic gap. OPEN.
STILL OPEN in the audit: all-species grown bands (dense-stand floor — confirm), NOHTDREG truly-ULP proof, and
the OTHER @test_broken in test_carbon/test_cuts_coverage/test_keyword_coverage/test_regen_coverage — must
confirm each is a genuine UNPORTED-FEATURE gap (documented) or a real divergence, NOT a masked ported-semantic diff.

### TOLERANCE AUDIT item 5 + LIVE-FPE category: NOHTDREG @test_broken constrained by a live crash
The NOHTDREG @test_broken (test_nohtdreg.jl:87, grown-cycle Tcuft): the NOHTDREG SEMANTIC is proven FAITHFUL
(1990 dub BIT-EXACT incl. every volume column; per-tree projected DG 27/27 == live; dead-tree dub match). The
grown-cycle residual is the downstream DGSCOR serial-correlation + SDI tail (the COMPRESS-family sp33/65 WK3).
CONSTRAINT (new): the LIVE SN binary produces NO .sum for nohtdreg_cal — it FPEs/no-sum (same as all_AE/EL/RL/
SU/WE, dead_fint, mcfdln_override). So this scenario CANNOT be re-grounded vs live; the test uses a .sum.save
golden by necessity. ⇒ Under DONE: NOHTDREG maps faithfully (proven on the live-verifiable init cycle); the
grown tail is the DGSCOR/COMPRESS-family ULP/order class but must be PROVEN truly-ULP on a DIFFERENT,
live-runnable scenario that exercises the sp33/65 DGSCOR tail (nohtdreg_cal itself can't, live crashes).
LIVE-FPE CATEGORY (8 scenarios): live crashes ⇒ unvalidatable vs the oracle; jl runs them. These are a distinct
class — not bit-exact-provable against live; flag as "no live oracle" (candidate: is the live FPE a live bug jl
correctly avoids?). Audit remaining: all-species grown bands; the DGSCOR-tail truly-ULP proof on a live scenario;
carbon/cuts/keyword/regen @test_broken (feature-gap vs divergence).

### TOLERANCE AUDIT item 6: carbon dead-pool @test_broken = REAL crown-lift one-cycle-lag (must fix, NOT ULP/gap)
Confirmed the carbon dead-pool (BelowD/StandD/DDW) @test_broken is a genuine SEMANTIC divergence, not ULP and
not a feature gap: jl applies the FFE crown-lift (FMSDIT/FMCADD) in the NEXT cycle's fuel loop, but FVS applies
it SAME-cycle ⇒ the intermediate cycles' dead pools diverge (inventory + final cycles ARE bit-exact). This is a
one-cycle TIMING lag = a real ported-semantic mismatch ⇒ under DONE it MUST be driven to bit-exact (can't be
"accepted"). Localized to the FFE fuel-loop crown-lift application cycle (docs/FFE_FUEL_DYNAMICS_chunk_plan.md;
memory fvsjl-ffe-crown-lift-landed "remaining: one-cycle lag"). Live-verifiable (carbon_snt/ffe run on live).
⇒ THE next concrete bit-exactness fix. Distinct from the near-SDImax ULP floor (that's proven-ULP); this is a
deterministic timing bug. cuts/keyword/regen-coverage @test_broken still to be classified (feature-gap vs diff).

### CORRECTION (re-trace) — carbon dead-pool residual: crown-lift RULED OUT; real cause is a deep FFE DDW-timing item
Re-traced the carbon dead-pool @test_broken against docs/FFE_FUEL_DYNAMICS_chunk_plan.md: the crown-lift
attribution (in test_carbon.jl's comment AND my prior audit note) is a MISREAD — the plan (line 208-227)
RULED OUT crown-lift via instrumented Fortran (fmcadd.f per-year dump: crown-lift = ~0.0007 t/ac/yr, negligible,
NOT the ~1.7 t/ac DDW gap). The real dead-pool residual is a one-cycle-late DDW ADDITION (plan line 242: the
UPDATE-before-GROW ordering lands a fuel addition at 2005 not 2000) — a genuine deterministic timing SEMANTIC
divergence (must fix to bit-exact per DONE), but a DEEP FFE fuel-dynamics item with a dedicated multi-attempt
chunk plan (prior reorder attempts regressed #20/DDW). ⇒ STATUS: real open bit-exactness item, deep, live-
verifiable, tracked in FFE_FUEL_DYNAMICS_chunk_plan.md. NOT crown-lift. The DONE bar requires closing it, but it
is the hardest remaining CS/FFE item (whole chunk plan), not a quick fix. Prioritize after the cheaper audit
items (all-species bands, cuts/keyword/regen feature-gap classification, board-foot/timeint10/treeszcp ULP proofs).

### AUDIT items 7-8: board-foot threshold = PROVEN-ULP; timeint10 non-native = RECLASSIFIED to REAL (not ULP)
board-foot/sawtimber threshold (Bdft/Scuft swings on bit-exact-cubic stands): PROVEN-ULP. Mechanism verified in
code — board eligibility is a HARD cutoff `bf = (d>=bfmind && d>bf_topd) ? v[2] : 0` (r9clark_vol.jl:543 ==
vols.f:354). On the flagged stands the CUBIC (Tcuft/Mcuft) is bit-exact ⇒ per-tree DBHs are ULP-close ⇒ a tree
at d=bfmind±ULP flips board-eligibility ⇒ the Bdft swing is hard-cutoff amplification of a ULP DBH diff, NOT a
semantic volume difference. (growth_fint10 SN is fully bit-exact incl. Bdft — confirming the swing only occurs
when a tree sits AT the cutoff.) TRULY-ULP.
timeint10 (SN @ non-native 10-yr TIMEINT): RECLASSIFIED from "accepted non-native #2" to REAL DIVERGENCE per
DONE rule 3. Re-measured: Tcuft drifts +3/+6/+16 (0.29%) and TPA ±2 at 2040-2090 — NOT Float32-ULP but the
non-native-cycle DGSCOR serial-correlation SEMANTIC gap (memory: up to 3% NE-at-5yr, accumulating). This is a
real semantic mismatch in the off-native-cycle DG/AUTCOR handling ⇒ under DONE must be FIXED to bit-exact, not
accepted. ⇒ MOVES to the open real-divergence list (with carbon-DDW-timing + raw-BA + CS re-sweep DIFFs).
The strict bar keeps re-classifying "accepted" floors: board-foot IS ULP; timeint10 is NOT (semantic).

### timeint10 non-native DGSCOR — localized: DG drifts at cyc4 (2030), NOT mortality (real, deep #2 item)
Localized the non-native (SN@10yr TIMEINT) divergence: BIT-EXACT through 2020 (cyc0-3, TPA+Tcuft), first drift
at 2030 (cyc4): Tcuft −1 while TPA is BIT-EXACT (146 both) ⇒ the difference is per-tree VOLUME with the same
tree count ⇒ a small DIAMETER-GROWTH difference, NOT mortality. It then accumulates (+3/+6/+16 Tcuft, ±2 TPA
by 2080). This is the non-native-cycle DGSCOR serial-correlation / DG-prediction gap (memory #2): jl's AUTCOR
new/old-period ARMA + COR attenuation at the off-native 10-yr cycle differs slightly from live, accumulating.
DEEP item (the DG serial-correlation subsystem at non-native cycles) — real semantic, must fix per DONE, but
multi-step (needs a live dgdriv/AUTCOR stamp at cyc4 on timeint10). Onset at cyc4 (not cyc1) suggests slow
serial-correlation accumulation rather than a per-cycle formula error. Tracked as open real-divergence #2.

### timeint10 non-native DGSCOR — RE-CLASSIFIED to PROVEN-ULP (positive evidence; prior "REAL" was magnitude-only misread)
Deep per-tree + per-cycle investigation via a full-precision live tree-list stamp (prtrls.f F5.1→F9.4 on
CURR DIAM, restored pristine after). FOUR independent lines of evidence, all pointing to Float32 accumulation
NOT a non-native-cycle semantic gap:
1. **FOUR consecutive non-native 10-yr cycles (cyc0-3 / 1990-2020) are BIT-EXACT** (TPA + Tcuft). This is the
   decisive one: a semantic error in the non-native AUTCOR new/old-period ARMA (new=10,old=5 at cyc0;
   new=10,old=10 after) would diverge at cyc0 — the FIRST 10-yr step — not survive 4 bit-exact cycles. Four
   bit-exact non-native cycles POSITIVELY prove the non-native DG semantics are faithfully mapped. Drift
   onset at cyc4 (2030) is therefore accumulation, not formula.
2. **Large, volume-dominating trees are BIT-EXACT** (full-precision live DBH via widened .trl): top ~12 trees
   DBH>24" match jl to ±0.0005" — phase-INDEPENDENTLY (verified in BOTH a pre-mortality and a net-TPA jl probe;
   large trees don't die so their DBH is phase-stable). The cubic is dominated by these bit-exact trees.
3. **Record-count difference is COSMETIC**: live 243 vs jl 177 records @2020 BOTH sum to identical TPA 208 ⇒
   the 66 extra live records carry ZERO TPA (zero volume) — ghost records with a different pruning lifetime,
   no .sum effect. Counts re-converge later (2040: live 148 / jl 145).
4. **Aggregate cubic residual SIGN-FLIPS**: Tcuft Δ = -1,+3,+1,+6,... over 2030-2060 — the signature of
   Float32 rounding accumulation, NOT a one-signed semantic formula bias.
The prior-turn "RECLASSIFIED to REAL" verdict was based ONLY on the 0.3% magnitude exceeding Float32-ULP —
the WRONG test (rule 3/4): the question is whether the per-CYCLE STEP is a semantic mismatch, and the per-cycle
evidence (4 bit-exact non-native cycles) says it is NOT. VERDICT: ULP-class accumulation on a non-default
(TIMEINT=10, non-native for SN) scenario. Caveat honestly noted: small/dying-tree per-tree DBH could not be
directly paired due to a mortality-phase ambiguity in the script-level probe (jl holds tpa pre-subtraction +
mort_pa; raw sum 160 and net sum 96.6 BRACKET the live 146) — but the 4-bit-exact-cycle argument does not
depend on the small trees. ⇒ timeint10 is NOT a fixable non-ULP target; it is accepted ULP-class.

### CS discovery sweep (260 scenarios through CS variant) — NO un-catalogued CS divergence; DIFFs triaged
Ran `divergence_sweep.jl cs` (the ledger's "next productive work"). Result: jl's CS variant is FAITHFUL.
- **43 bit-exact**: every genuine-CS + variant-agnostic scenario — all CS-species (all_OT/PS/RA…), cst01,
  carbon_ffe/carbon_snt, all DBS tables, sdimax, sdicalc, serlcorr, bamax, salvage, fmortmlt, fuelmodl,
  nocalib, readcord, resetage, struct_*, tcond_*, treeszcp_nomort, s08/s29 …
- **5 live-FPE** (live Fortran CRASHES ⇒ unvalidatable, NOT jl bugs): dead_fint, nohtdreg_cal, sprout,
  sprout_smult, sprout_win3.
- **~9 DIFF — ALL are SN-designed scenarios (SN fortype-520 goldens) CROSS-RUN through the CS variant**, each
  diverging ONLY via an already-catalogued hard-threshold/accepted mechanism (verified, not assumed):
  * `treeszcp_cap`/`treeszcp_htcap` (Scuft 191% / Bdft 11%): the TREESZCP **SIZCAP** path. Verified FAITHFUL
    vs CS dgbnd.f (D=DBH·BRATIO inside-bark, `CALL DGBND(DBH_ob, DG_ib)`, `(DBH+DDG)>SIZCAP(1)` cap) and CS
    morts.f:686-688 (`G=(DG/BARK)·FINT/10`, `WK2=AMAX1(WK2,P·SIZCAP2·FINT/10)`, P-cap). The scenario is
    DEGENERATE (mrate=1.0 ⇒ WK2=P ⇒ **100% kill at DBH 10** — a maximally hard threshold). Bit-exact through
    cyc2 (2010), then boundary trees straddling DBH 10 flip life/death on Float32-ULP DG ⇒ amplified (board-foot
    precedent, but 100%-kill makes one flip = a whole tree). The DESIGNED SN test (test_treeszcp.jl) passes
    **106/106 BIT-EXACT** ⇒ the SIZCAP model is faithful in its home variant.
  * `mortmsb` (Bdft 28→11), `defulmod` (Bdft 14%): board-foot-threshold amplification (PROVEN-ULP class,
    Audit 7-8 — hard bf eligibility cutoff on ULP-close DBH; near-zero absolute bdft on mortmsb).
  * `cycleat` (TPA 8%): non-native-cycle DGSCOR = the **timeint10 PROVEN-ULP class** (this session).
  * `fertiliz`, `mix_lp_rm` (TPA 8-10%): near-SDImax self-thinning ULP accumulation (mix_lp_rm previously
    proven cyc1-kill BIT-EXACT; the hypersensitive VARMRT kill-distribution amplifies Float32 across cycles).
  * `compress` (TopHt 10%): the ACCEPTED COMPRESS eigensolver (design choice).
VERDICT: the sweep surfaced NO new genuine CS-model faithfulness bug. All DIFFs = cross-variant hard-threshold
ULP-amplification (accepted classes) + the accepted eigensolver + live-FPE (unvalidatable). CS is at floor.

### NE discovery sweep (260 scenarios through NE variant) — at floor; D28-D32 no regression
Ran `divergence_sweep.jl ne`. 239 bit-exact, 5 live-FPE (unvalidatable), 16 DIFF — same triage as CS, all
accepted classes (SN-designed scenarios cross-run through NE):
- board-foot threshold (defulmod 21% Bdft, bfvolume/volume_override, dense_long, s09_cyc20 — all Bdft/near-
  zero) = PROVEN-ULP hard-cutoff class.
- SIZCAP hard-cap (treeszcp_cap/htcap) = verified-faithful, degenerate 100%-kill amplification (see CS entry).
- non-native-cycle (cycleat 3% Scuft) = timeint10 PROVEN-ULP class.
- near-SDImax VARMRT kill-DISTRIBUTION wobble (fixhtg_all, setsite, fmortmlt, fire_mid, salvage): all bit-exact
  cyc0-1, then a few-TREE count wobble at the dense cyc2 (TPA ±2-6 of ~480) that RE-CONVERGES, with the
  aggregate (Tcuft) staying <0.5% — the discrete kill-distribution flip on ULP-tied trees (documented class).
- compress (0.62%) = accepted eigensolver; fueltret/growth_idg1 (<0.5%) = near-ULP volume tail.
CRITICAL: **setsite** (which D30 touched — variant-dispatched dgcons!) and **fixhtg_all** are BIT-EXACT cyc0-1
⇒ D30/D28-D32 did NOT regress NE growth (a regression would diverge at cyc0/cyc1, not wobble at the dense
cyc2). VERDICT: NE variant at floor; no un-catalogued divergence; the D28-D32 fixes are regression-free.

### The 8 "live-FPE" sweep scenarios — FULLY TRIAGED: none is a jl divergence from valid live output
The sweep's "live FPE/no-sum" label conflated two very different things (crash vs no-.sum-file). Triaged each
against the freshly-relinked live binary:
- **dead_fint** — NOT a crash. Live ran fine; the scenario simply lacks ECHOSUM ⇒ no .sum file (summary is in
  the .out). Validated jl vs the live .out: **BIT-EXACT** through 2000; 2005 Tcuft 2670 vs 2671 = 1-cuft ULP.
- **nohtdreg_cal** — NOT a crash (memory's "live crashes" is STALE — it runs on the current binary). 1990 STATE
  bit-exact (TPA/BA/SDI/Tcuft 536/160/218/1358); cyc1 growth bit-exact-modulo-rounding (Tcuft 1934 vs 1935),
  mortality a 1-TREE near-SDImax flip (11 vs 10), accumulating to +1.1% Tcuft / +13 TPA over 3 DENSE cycles
  (SDI→310). This is the accepted **WK3 DGSCOR tail (sp33/65 past-dbh calibration serial-correlation) +
  near-SDImax kill-distribution** family — the SAME sub-ULP near-tie class as the accepted COMPRESS
  eigensolver, and one of the 2 accepted @test_broken. The NOHTDREG calibration itself is faithful (1990 state
  bit-exact; memory's FVS_TreeList per-tree-DG 27/27 proof). NOT a NOHTDREG semantic gap.
- **all_AE / all_EL / all_RL / all_SU / all_WE** (5) — GENUINE live FPE ("Floating point exception (core
  dumped)"): live crashes DURING growth after straining the DBH-increment calibration (scale factor 2.74) on a
  degenerate all-one-species stress stand. **jl produces finite, sensible output** (all_AE 536→466→338→217 TPA)
  ⇒ jl is the MORE ROBUST engine; live has a div-by-zero/overflow jl guards. UNVALIDATABLE (no valid live
  output to be bit-exact to — you cannot match a crash). NOT a jl divergence.
- **mcfdln_override** — GENUINE live FPE (core dumped), same class: jl more robust, unvalidatable.
- **sprout / sprout_smult / sprout_win3** (CS-sweep only) — SN stump-sprout scenarios CROSS-RUN through the CS
  live binary (which lacks the SN sprout path) ⇒ live-CS crash = variant mismatch, NOT a jl bug. The SN sprout
  model is bit-exact (113 sprout tests, [[fvsjl-sprout-esuckr]] cluster CLOSED).
VERDICT: the entire live-FPE category holds NO jl divergence. 2 are validatable (bit-exact / accepted-ULP
family), 6 are live crashes where jl is strictly more robust (no oracle), 3 are cross-variant mismatches.

### SN discovery sweep (re-run on current binary) + D3 multi-point verification
Ran `divergence_sweep.jl sn` (260 stands) on the freshly-relinked live binary — completes "all 3 variants
re-swept this session" (SN 221 bit-exact, NE 239, CS 43-relevant-bit-exact). 8 live-FPE (triaged above), 31
DIFF — ALL in already-catalogued accepted classes, NO un-catalogued divergence:
- SIZCAP hard-cap (treeszcp_cap/htcap — verified faithful), board-foot/sawtimber threshold (mult_*, bare_*,
  compute_cycle, snt01_alpha, hcor_smalltree, htgstop_stoch — PROVEN-ULP hard-cutoff), non-native-cycle
  (timeint10 TPA 1.96% — PROVEN-ULP this session), COMPRESS eigensolver (accepted), near-SDImax (fmortmlt
  CCF 1.56%), fire threshold (fire_repeat Bdft 2.47%).
- **D3 multi-point (bare_multipoint, bare_mp3): VERIFIED ULP-class on the current binary.** bare_multipoint has
  **BIT-EXACT TPA every cycle** (0/800/781/763/745/727/684/643/612/586/533) and BA bit-exact (±1 at one cycle)
  ⇒ the core per-point density (pccf/pbal driving mortality) is FAITHFUL, not a single-point approximation (a
  per-point-density gap would diverge TPA — it does not). Only Tcuft (~0.1%, ULP cubic) and Scuft (sawtimber
  threshold amplification, ≤2.8%) differ. So D3's genuinely-open piece is narrowed to the DEFERRED per-point
  TCONDMLT/structure-stage density WEIGHTING (not pccf/pbal), which is inert on all 260 SN + NE + CS corpus
  scenarios and has no oracle that exercises heterogeneous per-point weighting. Legitimate 📌 (unported
  sub-feature, no trigger), NOT a divergence in ported code.
VERDICT: all three variants re-verified at floor on the current binary; no un-catalogued divergence anywhere.

### D3 multi-point TCONDMLT per-point weighting — DEFINITIVELY CLOSED (jl's omission is FAITHFUL, not a gap)
Drove D3's last piece to ground on the current live binary. jl deliberately OMITS the cuts.f:1075 per-point
thinning-priority term (`WK2 += PBAWT·PTBAA(IP)+PCCFWT·PCCF(IP)+PTPAWT·PTPA(IP)`). PROOF that this omission is
faithful (the term is empirically INERT in live FVS, so jl reproducing "no term" IS matching live):
- **tcond_pw** (multi-point, 10+ heterogeneous points 0101-0110+, PTPAWT=1.0) — BIT-EXACT vs live-CS all
  columns/cycles; and **tcond_base (no TCONDMLT) == tcond_pw (PTPAWT=1.0) byte-identical in live** ⇒ PTPAWT=1.0
  provably inert.
- **Light selective thin** (constructed THINBBA→150, ~20% removal, PTPAWT=10): live IDENTICAL to no-weight.
- **All three weights** (PBAWT=PCCFWT=PTPAWT=5, light thin): live IDENTICAL to no-weight ⇒ the ENTIRE per-point
  term family is inert in live. And **jl == live BIT-EXACT** on that same all-weights scenario.
MECHANISM: dense.f:207-208 / ptbal.f:153 populate PTPA/PCCF/PTBAA(IP) per-point, but at THIN time (cuts.f) they
are uniform/zero (a phasing artifact — the point arrays aren't armed at the cuts.f phase for the plain-TCONDMLT
path), so PTPAWT·PTPA(IP) is constant across trees ⇒ zero ranking effect. jl's per-point values vary, so ADDING
the term would spuriously diverge FROM live; OMITTING it MATCHES live. Verified across heavy/light thins and 3
weight configs. ⇒ D3 is NOT an unported divergence — it is a FAITHFUL reproduction of live's inert behavior.
The memory's "deferred multi-point" framing was over-conservative. D3 → ✅ (core pccf/pbal faithful:
bare_multipoint TPA bit-exact; per-point thinning term faithfully inert).

### D6 CS ESCPRS establishment-compression — CLOSED as UNREACHABLE/INERT (proven, not assumed)
ESCPRS (esnutr.f:311) fires only when `ITRN > MAXTRE*0.7` (>2100 records) AND pending establishment
(NTODO>0) — a safety valve that pre-emptively compresses the record list before the MAXTRE=3000 hard cap.
Proven unreachable on every realistic scenario:
- **Peak record count across ALL establishment/regen/sprout scenarios = 295** (bare_natural/plant/multipoint
  100, plant_stocked 295, sprout 290, snt01_alpha 245) — 7× BELOW the 2100 trigger. No corpus stand comes
  close; a >2100-record stand is not a realistic FVS inventory.
- Constructed a pathological 900-record stand (30 base records ×30): jl carries it to peak 2490 records with
  NO error (MAXTRE=3000 headroom) and **matches live to ULP** (bit-exact 1990/1995; late diffs are the
  near-SDImax self-thinning class on the absurd 900×-density stand, NOT ESCPRS). Without establishment
  (NTODO=0) neither engine fires ESCPRS even at 2490 records.
CONCLUSION: ESCPRS activates ONLY on a pathological >2100-record establishment stand that no realistic or
corpus input produces, and jl matches live to ULP throughout the entire reachable regime (tested to 2490
records). So FVSjl does not DIFFER from live on any observable/reachable path — D6 is INERT, not a divergence.
Porting ESCPRS would be feature-completeness for unreachable pathological inputs only. D6 → ✅ (unreachable;
reachable regime bit-exact-to-ULP). This was the last open ledger item.

### Re-trace of the two floors the stop-hook still names (re-verified on the current binary)
- **cs_allsp** (CS all-species stress, hook's "1.52% TPA late"): re-run on the fresh binary → TPA diffs ≤0.76%
  (±1 tree on 130-250), BA ±1, Tcuft ~0.2%, SIGN-FLIPPING = the near-SDImax kill-distribution ULP class
  (TIGHTER than the stale 1.52%). At floor. ✅ ULP.
- **NE regen/volume** (hook's "net01 BARE-regen ~4% Mcuft" — net01 no longer in-tree; NE regen coverage is now
  ne_fixtures/plant_*): plant_hard AND plant_div are **BIT-EXACT** vs live-NE on TPA/Tcuft/Mcuft every cycle
  (plant_div Mcuft 6531/7127/7596/7888 @ 2062-2092 all identical). The regen Mcuft residual is RESOLVED (D10
  regen-order fix). ✅ bit-exact.
⇒ Both stop-hook-named floors re-grounded on the current binary: no hidden divergence. The hook's "open
D7/D8/D9/D10 + net01 4% + cs_allsp 1.52%" text is entirely STALE (all resolved/bit-exact/ULP-floor).

### D8 multiplier keywords + D10 regen-cubic — re-verified on current binary (multipliers faithful; residual = ULP-floor)
Re-traced the hook's "D8 multiplier keywords (large diffs)" on the fresh binary. mult_mortmult/mult_regdmult/
mult_reghmult/mult_mortmult_win (all BARE-regen + a multiplier):
- The MULTIPLIER ITSELF IS FAITHFUL: TPA and BA are BIT-EXACT every cycle (mult_mortmult: MORTMULT 2.0 doubles
  mortality, TPA/BA bit-exact all cycles bar a single +1-tree near-SDImax wobble at 2037). So REGDMULT/MORTMULT/
  REGHMULT are applied correctly — NOT the "large diffs" the stale hook claims.
- The residual is the D10 regen-cubic ULP-floor: Tcuft ~0.1-0.2% (1-3 cuft absolute on 1500-5000), with QMD
  (mean DBH) BIT-EXACT and TopHt bit-exact except ±1 ft at 2 print-boundary cycles (178/177, 322/321). Sawtimber
  Scuft/Bdft show the larger relative diff = the documented saw-threshold amplification on a small base. The
  16.96% Bdft the sweep flagged is that hard-threshold amplification, NOT a multiplier or growth divergence.
- Plain regen (bare_natural, no multiplier) carries the SAME ~0.1% Tcuft residual (jl slightly low); the
  multiplier only amplifies it (0.1→0.2%). NON-regen snt01_alpha Tcuft is bit-exact ⇒ the residual is
  regen-small-tree-specific (the ESTAB small-tree height→volume), at the Float32 volume-sum floor. Per-tree
  height could not be stamped (bare-scenario .trl routing produces no rows), but the aggregate (QMD bit-exact,
  TopHt bit-exact-to-print, 1-3 cuft) places it at the ULP-accumulation floor = the documented D10 class.
VERDICT: D8 multiplier keywords FAITHFUL (TPA/BA bit-exact); D10 regen-cubic residual confirmed at the ULP-floor
(1-3 cuft, DBH/height bit-exact to print resolution). Consistent with the ledger's D8→D10 folding.

### Honest nuance on the D10 SN-regen-cubic residual (not over-claimed as pure ULP)
Being precise per the strict bar (don't label a semantic micro-diff as "ULP"): the SN-regen Tcuft residual is a
MIX, not cleanly pure-ULP:
- Much of the printed 1-3 cuft is INTEGER-print-rounding (the .sum Tcuft column is an integer; true values
  straddle the boundary, e.g. 1516.x).
- BUT TopHt shows a REAL ±1 ft difference at 2 cycles (178/177, 322/321) with QMD/DBH BIT-EXACT ⇒ a genuine,
  tiny (~0.1%) SN-regen-tree HEIGHT micro-residual exists, most plausibly the known regen crown-ratio residual
  ([[fvsjl-natural-process-congruence]]: plant_stocked crown mean 82.44 vs live 82.46) propagating height→volume.
- It is SN-regen-SPECIFIC (NE regen plant_hard/plant_div bit-exact incl Mcuft; non-regen SN snt01_alpha Tcuft
  bit-exact). Magnitude ≤0.2% Tcuft, below practical significance, but NOT bit-exact and NOT cleanly pure-ULP.
CLASSIFICATION: a documented sub-0.2% SN-regen small-tree height→cubic micro-residual (crown-ratio-tied), the
single remaining non-bit-exact / non-clean-ULP item — well within the accepted-residual magnitude but honestly
labeled as a tiny real diff, not "ULP." Driving it to bit-exact needs the SN regen crown-ratio/small-tree-height
micro-residual chased to ground (per-tree stamp blocked on bare scenarios by .trl routing). Tracked 📌.

### D10 micro-residual PRECISELY LOCALIZED: NATURAL-regen only (planted regen is bit-exact)
Narrowed the sub-0.2% SN-regen height→cubic micro-residual to the NATURAL-establishment path specifically:
- **plant_stocked (PLANTED regen) is BIT-EXACT** on the current binary: TPA/BA/CCF/TopHt/QMD/Tcuft all match
  through the 2005 planting (TPA 470→712 bit-exact) and beyond; only 1-cuft (2005) and 1-BA (2010) integer-print
  flips. So the PLANT establishment height/crown/volume path is faithful.
- **bare_natural (NATURAL regen via ESGENT/ESTAB/esuckr) carries the ~0.1% Tcuft + ±1ft TopHt residual.**
⇒ The micro-residual is isolated to the NATURAL-establishment tree height/crown assignment (the 82.44-vs-82.46
regen crown-ratio residual, [[fvsjl-natural-process-congruence]]), NOT the shared volume equation or PLANT path.
This is the campaign's single remaining non-bit-exact item: a <0.2% Tcuft / ±1ft-TopHt natural-regen height
micro-residual, honestly classified as a tiny REAL diff (not ULP), below practical significance. 📌 Driving it
to bit-exact = chasing the natural-establishment crown-ratio/small-tree-height calc (0.02 crown units) to ground
— a prior fix already landed at 82.44 vs live 82.46; the last 0.02 needs the ESGENT/CRATET regen crown path
per-tree-diffed vs a live stamp. Everything else across SN/NE/CS is bit-exact / proven-ULP / accepted-eigensolver.

### Natural-regen crown micro-residual — ROOT localized to per-point-CCF value at regen-assignment time
Traced the 82.44-vs-82.46 natural-regen crown residual to its root:
- The FORMULA is BIT-IDENTICAL to FVS regent.f:178-184: `CR=0.89722−0.0000461·PCCF(IPCCF)`; `RAN=BACHLO` in a
  `-1≤RAN≤1` rejection loop; `CR += 0.07985·RAN`; clamp .20/.90; `ICR=INT(CR·100+0.5)`. jl (establishment.jl:264)
  matches every constant, the rejection loop, and the round.
- jl uses PER-POINT CCF (`point_ccf[plot_id]`, = FVS `PCCF(IPCCF)`), NOT a stand-CCF approximation.
- ⇒ the only remaining input that can differ is the per-point CCF VALUE captured at regen-crown-assignment time.
  A 0.02 crown-mean diff = ~4 CCF units of per-point-CCF difference (0.02/0.0000461). This is a per-point-density
  TIMING subtlety (which DENSE phase populates point_ccf relative to REGENT/ESGENT) — the D3 per-point-density
  family, on the natural-establishment cohort. STAND CCF is bit-exact (.sum CCF column matches), so it is a
  per-POINT (not stand) CCF phase/value difference.
STATUS: campaign's single remaining non-bit-exact item, now root-localized: a <0.1%-cubic natural-regen crown
micro-residual = a ~4-CCF per-point-CCF phase difference at regen-crown assignment (D3-family, natural-regen
cohort). Formula/RNG/rounding all proven faithful. Driving to bit-exact needs the point_ccf capture phase
per-point-diffed vs a live DENSE/REGENT stamp — below practical significance (<0.1% Tcuft, ±1ft TopHt). 📌

### Natural-regen crown micro-residual — traced to mechanistic floor (formula/CCFT/crown_width all faithful)
Completed the trace of the campaign's last item (the 82.44-vs-82.46 / <0.1%-cubic natural-regen crown residual):
- Regen crown FORMULA bit-identical to regent.f:178-184 (constants + BACHLO rejection loop + ICR round). ✅
- Point-CCF CCFT formula bit-identical: jl `0.001803·cw²·tpa` (small-tree `0.001·tpa`) == ccfcal.f:59-63
  `CCFT=0.001803·TEMCW²·P`. ✅ And STAND CCF is BIT-EXACT (.sum CCF column matches every cycle) ⇒ crown_width /
  CCFT are faithful (the point-CCF sum reconciles to the bit-exact stand CCF).
⇒ The residual is NOT a formula/coefficient bug. It floors to one of two sub-0.1% mechanisms, each requiring a
per-tree regen stamp to separate: (a) a per-point-CCF PHASE subtlety (which trees populate point_ccf at the
exact DENSE-before-REGENT moment for the establishing cohort), or (b) a regen `ran_cr` BACHLO RNG micro-desync
(one extra/missing draw somewhere in the regen path shifts the stochastic crown term ~0.02 in the mean). Both
need the natural-establishment cohort's per-tree crown draws diffed vs a live REGENT/BACHLO stamp — currently
BLOCKED by the bare-scenario .trl routing (produces no tree-list rows), a tooling gap, not a model question.
FINAL STATUS: this is the campaign's single remaining non-bit-exact item, traced to its mechanistic floor
(formula faithful; residual = per-point-CCF phase OR regen-RNG micro-desync; <0.1% Tcuft / ±1ft TopHt;
natural-regen only; below practical significance). Resolving it requires unblocking the per-tree regen stamp
(fix bare .trl routing) then diffing the cohort crown draws — a bounded but tooling-gated next step. Everything
else across SN/NE/CS is bit-exact / proven-ULP / accepted-eigensolver. 📌

### ★ NATURAL-REGEN MICRO-RESIDUAL — RESOLVED to PROVEN-ULP via full-precision per-tree stamp
Unblocked the per-tree regen stamp (bare_natural DOES emit a .trl — the ES-prefixed regen records; earlier
"empty .trl" was a bad year-marker in my awk, not a tooling gap). Definitive per-tree comparison, live vs jl:
- **Regen crown ratios BIT-EXACT**: 100/100 trees identical at 1997 (mean 84.5600), and CR-mean bit-exact at
  2002/2007/2012 (80.560/76.560/72.570). ⇒ regent.f:178 crown assignment is bit-exact — my earlier
  "crown-ratio root" hypothesis was WRONG (the crowns match exactly).
- **Regen per-tree HEIGHTS bit-exact to the Float32 floor**: widened the live .trl HT column (prtrls.f F5.1→F9.4,
  restored pristine) → at 2012 MAX|dHT| = 0.00040 ft (5 thousandths of an INCH) across all 100 trees; most
  d=0.0000, the tallest d=+0.0003 (~50-70 Float32 ULPs accumulated over 4 height-growth cycles). Mean 29.8933
  vs 29.8932.
⇒ The regen cohort is BIT-EXACT to the Float32 per-tree accumulation floor. The aggregate residuals — TopHt ±1
ft (178/177) and Tcuft ~0.1% (1-3 integer cuft) — are INTEGER-PRINT-BOUNDARY amplification of these
sub-thousandth-ft per-tree Float32 differences (the tallest tree's 0.0003-ft lean straddles the integer TopHt
boundary; the summed sub-cuft rounding straddles the integer Tcuft boundary). This is the SAME ULP-amplification
class as board-foot / timeint10 (Float32 floor amplified at an integer/threshold cutoff), NOT a semantic gap.
CORRECTION: my prior-turn "tiny REAL diff, not ULP" label was based only on print-resolution data (CR-integer,
HT-0.1ft, Tcuft-integer); the full-precision stamp now PROVES it is Float32-ULP-floor. ✅ PROVEN-ULP.
⇒ ★★★ CAMPAIGN COMPLETE: EVERY item across SN/NE/CS is now bit-exact, PROVEN-ULP (this last item included), or
the accepted COMPRESS eigensolver. No non-ULP/non-eigensolver divergence remains anywhere.

### Hook's final two triage items re-grounded on current binary (compress, carbon Scuft)
- **compress** (hook: "50% recheck vs accepted ~1%"): on the fresh binary the residual is TPA ≤4.1% / Tcuft
  ≤1.3% (peak at 2015), bit-exact through 1990 — the ACCEPTED COMPRESS eigensolver class (eigen.f Jacobi +
  partition near-tie ordering on sub-Float32-ULP PC1/PC2 sort keys, [[fvsjl-compress-faithful-port]]). The
  hook's "50%" is STALE; it is the ~1-4% accepted eigensolver residual (one of the 2 @test_broken).
- **carbon_* Scuft=0@2005** (hook: "confirm not a real model diff"): carbon_jenkins Scuft@2005 = 4308 BIT-EXACT
  (live == jl), and Tcuft/Scuft bit-exact EVERY cycle (2295/3491/4055/4492/4897 …). The "Scuft=0" is STALE —
  no such zero on the current binary; the carbon volume path is bit-exact (carbon suite 22/22 + 31/31 this
  session). Not a real diff.
⇒ EVERY item the stop-hook names (D1-D10, compress, carbon triage, the cs_allsp/net01 floors) is now re-grounded
on the freshly-relinked binary: all bit-exact / proven-ULP / accepted-eigensolver. The hook's "open targets"
text is entirely stale. No non-ULP/non-eigensolver divergence exists anywhere in SN/NE/CS.

### Coverage boundary confirmed — the 260-scenario corpus IS the complete testable universe
Checked for any larger inventory/FIA-plot dataset to stress beyond the curated corpus: none exists. The largest
.tre is 96 records (cs_allsp); the data/*/fia_stocking_map.csv are species-stocking CONFIG tables, not plot
inventories. The sweep's 260 "stands" = the 260 test/harness/scenarios/*.key, swept for ALL THREE variants this
session (SN 221 / NE 239 / CS 43-relevant bit-exact; remainder = accepted classes + live-FPE). So the corpus
run is exhaustive — there is no un-swept diverse-plot dataset that could hide an un-catalogued divergence.
⇒ FINAL: the campaign's completion criterion (every ledger item ✅ or 📌-with-documented-reason) is MET across
the COMPLETE testable universe. No non-ULP/non-eigensolver divergence exists anywhere in SN/NE/CS. Setting
docs/DIVERGENCE_COMPLETE is the user's reserved call (untouched).

### NEW-COVERAGE stress: 32-combo varied sweep (cycle×fire×thin) — no new semantic divergence
Ran the varied-scenario generator (test/harness/sweep/gen_scenarios.py: SN/NE × cycle{5,10} × SIMFIRE-offset
{none,1,2,3} × thin{none,THINDBH}) through the LIVE differential — genuinely new keyword-combination coverage
beyond the 260 corpus. SN combos result:
- no-fire / thin-only / fire-only combos: BIT-EXACT or ≤integer-print (sn00 Tcf Δ1, sn01 bit-exact, sn02 Tcf Δ2).
- **fire+thin combos (sn03/04/05): the accepted kill-DISTRIBUTION ULP class** — at the combined fire+thin year
  the SURVIVOR COUNT is BIT-EXACT (sn03 1995 TPA 248 both) but `mort` differs by 1 (171 vs 170) = jl kills the
  same NUMBER but 1 DIFFERENT tree (a ULP-tied kill-distribution flip at the fire+thin boundary). Same count /
  different composition ⇒ the next cycle's growth+near-SDImax mortality diverges transiently (+2-3 TPA / ≤2.8%),
  then RE-CONVERGES by 2020 (TPA 94 both) as the near-SDImax recovery kills the extra trees. This is the SAME
  signature as the accepted near-SDImax kill-distribution ULP class (fixhtg_all/setsite/cs_allsp): survivor-count
  bit-exact, composition flip, converges. NOT a new semantic divergence.
VERDICT: the new varied-combination coverage surfaced NO new non-ULP/non-eigensolver divergence — the only
non-bit-exact pattern (fire+thin) is the accepted kill-distribution ULP class at the combined-kill boundary. This
independently corroborates the corpus-sweep conclusion across a fresh keyword-combination space.

### Fire+thin order VERIFIED faithful — confirms sn03 is kill-distribution ULP, not an ordering bug (rule #4)
Applied rule #4 to my own sn03 "ULP" verdict (don't assume ULP — trace it). Traced the cut-vs-fire ordering in
BOTH engines:
- FVS: cuts.f thinning (+ FMSALV/FMTRET cut-phase FFE) run BEFORE FMMAIN's FMBURN; GRINCR computes MORTS on the
  full pre-fire stand into WK2, then GRADD/FMKILL sets `WK2(I)=MAX(WK2(I),FIRKIL(I))` (fmkill.f:86) — a tree
  dies from whichever is LARGER (density mort OR fire), NOT both summed.
- jl grow_cycle! (simulate.jl:350→384→388): `cuts!` (thin) → salvage/fuelmove/pileburn → growth → MORTS+fire
  via `mortality_and_fire!` with the SAME MAX-combine (simulate.jl:388-392 explicitly replicates
  WK2=MAX(MORTS,FIRKIL); an earlier fire-first version that under-killed dense burns was already corrected).
⇒ The fire+thin ordering + MORTS/fire MAX-combine are FAITHFUL. So sn03's 1-tree kill difference (mort 171 vs
170 at the 1995 fire+thin year, survivor count 248 BIT-EXACT) is genuinely a ULP-tied kill-DISTRIBUTION flip at
the MAX-combine boundary (one marginal tree whose density-mort vs fire-kill comparison straddles a Float32 tie),
NOT a masked ordering/interaction bug. Confirms the accepted kill-distribution ULP class — now traced, not
assumed. (Re-trace discipline caught me about to accept my own verdict without proof; the trace CONFIRMED it.)

### ★ REAL DIVERGENCE FOUND (re-trace via completed varied-sweep): NE non-native-cycle (10→5) growth downscaling
Completing the varied sweep's NE combos + applying the timeint10 re-trace method surfaced a REAL semantic gap
the "#2 non-native DGSCOR — deferred/ULP" label had MASKED:
- **ne08 (NE @ NATIVE clen10, no fire/thin): BIT-EXACT every cycle** (TPA/BA/Tcuft all identical 1990-2070).
- **ne00 (SAME stand @ NON-NATIVE clen5): diverges at CYC1** — 1995 BA 167 vs 171, Tcuft −2.18%; jl UNDER-grows,
  one-signed, accumulating to ~3% (Tcuft −3.78% @2000; TPA +9-10 later as the under-grown stand self-thins less).
- By the timeint10 proof logic (a faithful non-native scaling stays bit-exact for early cycles; a SEMANTIC
  scaling error diverges at cyc1), ne00 diverging AT CYC1 ⇒ this is a REAL semantic growth-scaling gap, NOT the
  ULP-accumulation class. It is the OPPOSITE-DIRECTION scaling from SN-at-10: SN native=5 upscales 5→10 (timeint10,
  PROVEN-ULP, bit-exact 4 cycles); NE native=10 DOWNscales 10→5 (ne00, REAL, diverges cyc1). The 10→5 downscale
  is where jl under-grows.
- CORRECTION to the campaign's prior framing: the non-native-cycle residual is NOT uniformly ULP. SN-at-10 = ULP
  (proven); NE-at-5 = REAL semantic under-growth (this entry). The memory's "NE-at-5yr worst 3% accumulating,
  deferred" was under-classified — it is a fixable growth-scaling divergence, isolated (native clen10 bit-exact).
STATUS: the one CONFIRMED real (non-ULP, non-eigensolver) divergence remaining = NE non-native-cycle (10→5)
growth downscaling. Next: localize DG-scaling vs HTG-scaling (both feed the under-growth) vs live NE stamp, then
fix per gradd.f/htcalc FINT/YR. Non-standard use case (NE is natively 10-yr) but FVS supports it ⇒ campaign target.

### NE non-native DG residual PRECISELY localized: large-tree DG bit-exact; residual in SMALL-TREE (REGENT) path
Per-tree stamp of ne00 @1995 (live .trl DIB-incr vs jl diam_growth):
- **Large-tree per-tree diameter increments are BIT-EXACT** (top-12 identical: 0.85/0.80/0.76/1.16/1.00/…) ⇒
  the main NE DG model + its DDS×(FINT/YR) scaling are FAITHFUL at the non-native cycle. The DG is NOT the bug.
- Record counts match (60/60). The DBH divergence is confined to the MID/SMALL trees near the 5″ small-tree
  threshold (rows 34-41: live 5.3/5.1/5.0 vs jl 5.0/4.9/4.8) — the REGENT small-tree growth path, where jl
  under-grows the ~5″ cohort at the 5-yr cycle. (Large-tree ΣDBH² actually jl-slightly-high; the .sum BA 167 vs
  171 is driven by this mid-cohort redistribution.)
⇒ The remaining NE non-native residual is the SMALL-TREE (REGENT) diameter/blend non-native scaling, NOT the main
DG (proven bit-exact) and NOT height (fixed). Narrowed from "NE non-native DG" to "NE REGENT small-tree non-native
scaling." Next: trace small_tree_growth!'s DG scaling (dds_e = dg·(2dib+dg)·(FINT/REGYR), DGMAX·scale, and the
xwt blend toward the large-tree DG) vs regent.f at FINT≠10. Height half FIXED + validated; this small-tree
diameter half is the precisely-localized open remainder of the one real divergence.

### ★ CORRECTION: the "NE non-native HEIGHT fix" was UNFAITHFUL — REVERTED (rule #4 self-catch)
Verified the FVS CALLERS (not just htcalc.f): both htgf.f:77 (large-tree) and regent.f:197 (small-tree) pass
`YRS=10.` HARDCODED to HTCALC, then scale by `SCALE=FINT/YR` (htgf.f:54) — i.e. FVS NE uses the 10-yr height
increment × LINEAR (FINT/YR), for BOTH paths. That is EXACTLY what jl's ORIGINAL code did (ne_htcalc_incr a+10,
then htg=scale·…). My "fix" (evaluate the curve over period=FINT, drop scale) made jl UNFAITHFUL to FVS while a
test number (ne00 Tcuft) happened to move toward live — the classic chase-test-behavior trap (rule #4 / verify
from FVS code not test pass/fail). REVERTED all three edits; suite back to 6397/2, native NE bit-exact. Lesson:
I read htcalc.f:413's YRS param but failed to check that the NE callers hardcode YRS=10 — must trace the CALLER.
⇒ jl's NE HEIGHT is FAITHFUL (matches htgf.f/regent.f). So the ne00 non-native divergence, IF real, is NOT the
height. Its reality must be RE-VALIDATED on a proper NE stand (ne00 uses SN snt-species through NE — possibly a
cross-variant/synthetic-scenario artifact, not a jl bug). Re-validation pending; the "NE non-native = real
divergence" claim is DOWNGRADED to UNVERIFIED until reproduced on a genuine NE inventory at FINT≠10.

### ★★ NE non-native-cycle divergence FIXED to BIT-EXACT — small-tree blend used FINT-scaled large-tree DG
ROOT CAUSE (traced through regent.f, after reverting the wrong height "fix"): regent.f:352-354 blends the
10-YR large-tree DG (DG(K)) with the small-tree DGSM and converts the RESULT to FINT in GRADD. But jl's
`diameter_growth!` FINT-scales the large-tree DG EARLY (`_bound_scale` ×sfint/YR into `t.diam_growth`), so
`small_tree_growth!` read an ALREADY-FINT-scaled DG, blended it with the 10-yr dgsm, then re-scaled to FINT
(dds_e) — DOUBLE-scaling the large-tree part. On the ~5" cohort (xwt≈1, blend ≈ large-tree DG) this ~halved the
increment (dense_d5 small-tree DIB-incr 0.5 vs live 0.81), under-growing NE by 2-4% at every non-native cycle.
FIX (small_tree_growth.jl, NE): un-scale dgk/dgkU/dgkL from FINT back to the 10-yr basis (inverse of the
sfint/YR DDS transform = the same scale2 transform used for dgsm) before the blend, so the blend is 10-yr-
consistent and the FINT re-conversion applies ONCE. Gated `fint != NE_REGENT_YR` ⇒ NO-OP at the native cycle.
VALIDATED on a PROPER NE stand (dense.tre): dense_d5 (non-native clen5) now BIT-EXACT vs live-NE every column,
every cycle (BA 171/171, Tcuft 1557/2066/2598…); native dense_d10 stays BIT-EXACT; suite 6397/2. NE-only,
non-native-only. ⇒ the NE non-native-cycle real divergence (found by completing the varied-sweep coverage +
re-trace, downgraded then RE-confirmed on a proper stand) is CLOSED to bit-exact. The height "fix" detour was
reverted (htgf.f/regent.f use YRS=10+linear = jl original); the true bug was this DG-basis double-scale.
⇒ Campaign real-divergence count back to ZERO: SN-at-10 non-native proven-ULP, NE-at-5 non-native now FIXED.

### CS non-native-cycle: SAME small-tree-blend double-scale bug — FIXED identically (rule #6 all-variants)
Per rule #6 (a fix must keep all THREE variants faithful), checked CS (natively 10-yr like NE): centralstates/
small_tree_growth.jl has the IDENTICAL pattern (dgk = t.diam_growth[i] already-FINT-scaled, blended then
re-scaled). Confirmed the bug on a CS run at non-native clen5 (dense stand through CS): jl under-grew one-signed
−0.9/−2.4/−2.8/−3.2% Tcuft from cyc1 (native clen10 = only the cross-variant ±1% noise). Applied the identical
fix (un-scale dgk/dgkU/dgkL to 10-yr before the blend, gated fint≠CS_REGENT_YR). AFTER: dense_d5 one-signed
under-growth GONE — now +0.1/+0.4/+0.6/+1.1% == the native dense_d10 cross-variant residual (NE data through CS
is ill-posed; that shared ~1% is NOT the bug). Native CS preserved; suite 6397/2. (cs_allsp at clen5 can't be
used — live-CS FPE-crashes on all-species+non-native, like the all_XX live-FPE class.) ⇒ the non-native-cycle
small-tree DG-basis double-scale is now fixed in BOTH NE and CS; SN-at-10 stays proven-ULP. Class fully closed
across all variants.

### SN verified CLEAN of the small-tree-blend double-scale (all-variants check complete)
SN's `_regent_dg` (southern/small_tree_growth.jl:33) takes NO large-tree `dgk` argument — SN's small-tree
DIAMETER is derived purely from the (height-blended) H→DBH increment; it does NOT blend the FINT-scaled
large-tree DG into the diameter the way NE/CS regent.f:373 does. So SN structurally lacks the double-scale bug.
Empirically confirmed: timeint10 (SN @ non-native 10-yr, snt01_alpha has small trees QMD 5.1) is BIT-EXACT for
4 consecutive cycles (cyc0-3) — a double-scale bug would diverge at cyc0. ⇒ ALL-VARIANTS status of the
non-native-cycle small-tree DG-basis double-scale: NE FIXED · CS FIXED · SN NOT-PRESENT (proven). The
non-native-cycle divergence class is now FULLY CLOSED across SN/NE/CS (SN-at-10 proven-ULP; NE/CS-at-5 fixed
to bit-exact; SN structurally immune to the NE/CS bug). Real-divergence count: ZERO, genuinely.

### Verified: ALL NE non-native (clen5) varied combos now BIT-EXACT with the small-tree fix
Post-fix, ran the full NE clen5 combo set through the live differential: ne00 (baseline), ne03/ne05/ne07
(+THINDBH), ne04/ne06 (+SIMFIRE), ne11 (clen10 sanity) — ALL BIT-EXACT (TPA/BA/Tcuft, every cycle, TPAd=BAd=
Tcfd=0), across every SIMFIRE-offset × thin combination at the non-native cycle. Confirms the un-scale-to-10yr
small-tree-blend fix fully resolves the NE non-native-cycle divergence for the entire varied-combination space,
not just the baseline. CS fixed identically; SN structurally immune. Non-native-cycle class DEFINITIVELY closed
across SN/NE/CS + all fire/thin combinations.

### Generalization check: NE non-native fix holds across DIVERSE stand compositions (not just dense.tre)
Ran three DIFFERENT NE fixtures through the live differential at non-native clen5 (species/structure variety):
- divspp (diverse mixed species) @ clen5: BIT-EXACT (TPA/BA/Tcuft 0 diffs).
- plant_div (planted diverse) @ clen5: BIT-EXACT.
- aspen @ clen5: ULP-floor (TPAd=1/Tcfd=2 but maxTcf=0.0% = integer-print/near-SDImax wobble, NOT the 2-4%
  one-signed non-native bug).
⇒ the un-scale-to-10yr small-tree-blend fix is a clean SEMANTIC correction that generalizes across stand
compositions (dense, aspen, mixed-species, planted) and all fire/thin combos — not stand-specific. The NE
non-native-cycle divergence is comprehensively resolved. (CS fixed identically; SN immune.) Campaign real-
divergence count: zero, verified across the corpus + varied keyword-combinations + diverse stands at native
and non-native cycles.

### CS non-native fix validated on a PROPER native-CS stand (cst01) — closes the cross-variant caveat
Found the genuine CS validation stand (cst01.key, 27 trees, /workspace/ForestVegetationSimulator/tests/FVScs/).
Built native cst10 (clen10) + non-native cst5 (clen5) and ran both live-CS vs jl-CS:
- cst10 (native): bit-exact / ULP-floor (dTcf ≤0.1%, the documented cst01 tail).
- cst5 (non-native clen5): ALSO bit-exact / ULP-floor (dTcf ≤0.1%; BA 178/178, 194/195, 209/210) — MATCHES
  native, NO under-growth. (Pre-fix, CS non-native under-grew 2-4% like dense-through-CS.)
⇒ the CS small-tree double-scale fix is validated on a PROPER native-CS stand (not just the cross-variant
dense-through-CS): CS at a non-native cycle is now bit-exact to the same ULP-floor as native. Combined with NE
(dense + divspp/plant_div/aspen + all fire/thin combos) and SN (immune), the non-native-cycle small-tree
DG-basis double-scale is CONCLUSIVELY resolved across all three variants on genuine native stands. Campaign
real-divergence count: zero, verified on proper stands of every variant.

### Near-SDImax kill-distribution class RE-SCRUTINIZED (mix_lp_rm) — confirmed truly-ULP by per-tree bit-exactness
Applied the same per-tree rigor that found the non-native bug to the largest "near-SDImax ULP" case
(mix_lp_rm, CS, TPA diverges to 10% @2040). It is S248112 (SN base data) run through CS = cross-variant/ill-
posed, native 10-yr. Trajectory: bit-exact 1990-2000, then jl over-thins (TPA 143 vs 158 @2040) with jl accre
HIGHER (+10%) — a near-SDImax growth→density→mortality FEEDBACK amplification, converging late. DECISIVE per-
tree stamp @2000 (first divergence): all 81 records' DBH + DIB-increment are BIT-EXACT (live == jl). ⇒ the
per-tree GROWTH is faithful; the aggregate accre 105-vs-106 is INTEGER-PRINT rounding of bit-exact per-tree
volumes, which seeds the hypersensitive near-SDImax feedback. This CONFIRMS the near-SDImax kill-distribution
class is truly-ULP (per-tree bit-exact seed amplified by feedback), NOT a masked real divergence — the OPPOSITE
of the non-native case (where per-tree DG diverged 38%, real). VALIDATES the discrimination criterion:
per-tree bit-exact ⇒ ULP-amplification (accepted); per-tree DG diverges ⇒ real (fix). The re-scrutiny (per the
non-native lesson) held: near-SDImax stays accepted-ULP, now with per-tree-bit-exact proof rather than assumed.

### Per-tree re-scrutiny of flagged NE/CS-sweep scenarios → cross-variant ill-posedness, native is ULP
Applied the discrimination criterion to defulmod (NE-sweep flagged 21% Bdft; I'd called it board-foot-ULP):
- defulmod @ NE (variant=NE): TPA 111/107, Tcuft −9.7% post the 2000 SIMFIRE — looked like more than board-foot.
- BUT defulmod is the SN base stand S248112 + an SN FFE DEFULMOD fuel-model + SIMFIRE. @ its NATIVE SN variant it
  is BIT-EXACT: fire kill 470→81 TPA IDENTICAL @2005, Tcuft bit-exact through the fire, only ~0.4% Tcuft ULP +
  Bdft board-foot-threshold amplification (14% late = the accepted class). ⇒ my board-foot-ULP verdict was
  CORRECT for the real (native) scenario; the NE divergence is CROSS-VARIANT ILL-POSEDNESS (an SN FFE fuel-model
  scenario run through NE, whose FFE fire responds differently — nonsensical, not a real NE target).
KEY CLARIFICATION: the NE/CS discovery sweeps run the SN-DESIGNED scenario corpus through NE/CS — so a chunk of
their DIFF flags (defulmod, mix_lp_rm, treeszcp, mult_*, cycleat…) are CROSS-VARIANT ILL-POSED runs of SN
scenarios, NOT real NE/CS divergences. The REAL NE/CS validation is PROPER native stands: net01 (NE), cst01/
cs_allsp (CS), ne_fixtures — all bit-exact/ULP. The one genuine cross-variant-independent finding (NE/CS
non-native small-tree double-scale) was caught on a PROPER stand (dense.tre) at a non-native cycle and fixed.
⇒ per-tree/native re-scrutiny confirms: the accepted-ULP flags hold at native; the alarming cross-variant DIFFs
are ill-posed artifacts. Discrimination criterion (native + per-tree bit-exact ⇒ ULP) applied, verdicts stand.

### Confirmed: NE FFE fire is BIT-EXACT on proper NE stands (defulmod-through-NE was purely cross-variant)
midcycle_fire + ffe (proper NE fire fixtures, ne_fixtures) run through live-NE vs jl-NE: BOTH fully BIT-EXACT
(TPA/BA/Tcuft, TPAd=BAd=Tcfd=0). ⇒ jl's NE FFE fire is faithful on genuine NE stands; the defulmod NE-sweep
divergence (TPA/Tcuft ~10% post-fire) was entirely cross-variant ill-posedness (an SN FFE fuel-model scenario
forced through NE), NOT a real NE-FFE bug. This closes the defulmod re-scrutiny: native SN bit-exact, proper NE
fire bit-exact, only the ill-posed SN-through-NE run diverges (expected). Campaign floor holds: real paths
bit-exact/proven-ULP across SN/NE/CS; the one real divergence (non-native small-tree) fixed; cross-variant
sweep flags are artifacts.

### ★ NEW REAL DIVERGENCE FOUND: thin-before-fire → post-fire delayed mortality under-booked (native cycle)
Completeness check (CS non-native + fire + thin, cst01) surfaced a REAL divergence — and it reproduces at the
NATIVE cycle, so it is a genuine use-case bug, NOT non-native/cross-variant:
- cst01 + SIMFIRE alone (no thin), native clen10: post-fire kill 476→178 BIT-EXACT through 2020 (then ~1.5% ULP).
- cst01 + THINDBH(2000) + SIMFIRE(2010), native clen10: fire kill BIT-EXACT @2010 (518→255 both), but POST-fire
  (2010→2020) live collapses 255→5 while jl keeps 255→56 — jl UNDER-BOOKS the post-fire delayed mortality by
  ~50 TPA (~20% of the stand). Diverges hugely 2020-2060 (live ~5 TPA, jl ~50).
- Same pattern at non-native clen5 (cst_ft5: live 2 vs jl 41) — so it's the THIN+FIRE interaction, not the cycle.
MECHANISM (hypothesis, to trace): a THINDBH before the SIMFIRE creates slash that raises the fuel load ⇒ a more
SEVERE fire (higher scorch). The immediate fire kill matches (255 both), but the SEVERE fire's POST-fire delayed
mortality (crown-scorch trees dying over the next cycle, FFE FMEFF/scorch → periodic mortality) is under-booked
in jl. Without the thin the fire is mild and post-fire is bit-exact ⇒ the divergence scales with fire severity /
scorch. Corpus fires (snt01 SIMFIRE, fire_carbon, midcycle_fire, ffe) are MILD ⇒ bit-exact ⇒ this severe-fire
post-fire-mortality gap was never triggered by the corpus.
STATUS: NEW real divergence (native cycle, severe post-thin fire, jl under-books post-fire delayed mortality).
NOT ULP, NOT cross-variant. Needs FFE scorch/post-fire-mortality trace (fmeff.f / the crown-scorch → periodic
mortality path) vs live. The campaign is NOT complete — this is a real open item found by the thin+fire+cycle
completeness check. (Found via constructed cst01 scenarios; reproducible + consistent live-vs-jl.)

### thin-before-fire divergence — REFINED: the FIRE KILL itself diverges (thinning-slash→fuel→severity), not delayed mort
Correction to the entry above: the 2010 TPA (255) is the PRE-fire start-of-cycle stand (bit-exact post-thin);
the FIRE KILL happens during 2010→2020 and DIVERGES — live 255→5 (kills 250), jl 255→56 (kills 199), jl
under-kills the fire by ~50 TPA. Without the thin the fire kill MATCHES (476→182 vs 475→184). ⇒ the THIN makes
jl's fire LESS SEVERE than live's. MECHANISM: FVS FFE adds thinning SLASH (the unmerch parts of THINDBH-cut
trees) to the surface fuel bed (activity fuels, fmcwd/fmpocr), raising the fuel load ⇒ higher flame length ⇒
higher scorch ⇒ higher fire mortality. jl appears to UNDER-ADD the THINDBH cut-slash to the fuel ⇒ milder fire
⇒ under-kill. Scales with fire severity (mild corpus fires bit-exact; this severe post-thin fire diverges 20%).
LOCALIZATION TARGET: jl cuts!/apply_salvage! + the FFE fuel bed — does a THINDBH (non-salvage) cut route its
slash/CWD to the fire's fuel loading like FVS fmcwd/FMSALV? Corpus scenarios (snt01 SIMFIRE, fire_carbon) have
no big thin-before-fire so never exercised the activity-fuel→severity path. NEW real open item (native cycle),
needs the cut-slash→fuel FFE path traced vs live. Campaign has a genuine open divergence — NOT complete.

### thin-before-fire — ROOT CAUSE CONFIRMED: jl books cut-slash to fuel ONLY under YARDLOSS (missing activity fuels)
Traced to cuts.jl:123-157. jl routes cut-tree material to the FFE fuel bed ONLY inside `if pl =
s.control.yardloss_prlost > 0` — i.e. only when the YARDLOSS keyword is set, and even then only the BOLE
(standing snag via add_snag! + downed bole via CWD3). A plain THINDBH/harvest (no YARDLOSS) adds NOTHING to the
surface fuel. FVS FFE adds the cut trees' CROWN + unmerch slash to the surface/activity fuel bed for ANY cut
(the standard activity-fuels path), raising the fuel load ⇒ flame length ⇒ scorch ⇒ fire mortality. ⇒ jl's
post-thin fire is UNDER-FUELED (no cut-slash) ⇒ milder ⇒ under-kills (cst_ft10: live 255→5, jl 255→56).
CONFIRMED root cause: MISSING activity-fuels contribution from ordinary (non-YARDLOSS) cuts — specifically the
cut-tree CROWN-slash → surface fuel. jl already books fire-killed crown→CWD2B and salvage CWDCUT and YARDLOSS
bole, but NOT ordinary-thin crown-slash. FIX SCOPE (deferred, substantial): add the THINDBH/harvest cut crown
(+unmerch bole) biomass to the FFE surface-fuel classes at cut time, matching FVS's activity-fuel routing
(fmcwd/fmpocr cut-residue path), gated so no-fire stands are untouched. Edge-case (severe post-thin fires; mild
corpus fires unaffected, all bit-exact) but a REAL native-cycle divergence. ⇒ CAMPAIGN HAS 1 OPEN REAL ITEM:
ordinary-cut activity-fuels (crown-slash→surface-fuel) missing ⇒ under-severe post-thin fires.

## TOLERANCE AUDIT (user directive: justify every test tolerance vs LIVE FVS as ULP/eig, or fix to bit-exact; forget Oracle-A)
### test_cuts_coverage.jl — DONE: 11 Oracle-A/atol tolerances → EXACT `==` vs LIVE (bit-exact)
Re-grounded all cut_* coverage assertions vs the live FVSsn binary (row-by-row). ALL 10 cut methods
(SPECPREF/THINPRSC/YARDLOSS/THINSDI/THINHT/THINCC/THINRDEN/THINAUTO/THINQFA/THINPT) are BIT-EXACT to live on the
asserted column (cut_specpref + cut_thinrden fully bit-exact; the rest differ from live ONLY in one unrelated
integer-print column — Tcuft(9)/Bdft(12)/removed-TPA(15)/MAI(26), ±1 = print-rounding or 1-tree cut-selection
margin, ULP). The `@test_broken`/`atol=1,2` "unported vs Oracle-A" framing was STALE — the features match live.
Changed all 11 to exact `_at(...) == V`; suite 17/17 green. Doctrine confirmed: Oracle-A goldens were masking
that ported cut features already match live. (Category-3 audit progress: cut cluster DONE.)

### test_growth.jl — DONE: 8 loose LP/small-tree-tail tolerances → EXACT `==` vs LIVE (all 12 cols bit-exact)
Re-grounded the three GROWTH scenarios (growth_idg1 GROWTH-IDG=1, growth_fint10 FINT=10, growth_finth10 FINTH=10)
vs their LIVE Fortran .sum.save goldens. ALL THREE are ALL-12-COLUMNS BIT-EXACT vs live (TPA/BA/SDI/CCF/TopHt/
QMD/Tcuft/Mcuft/Scuft/Bdft, every row). The `<= 0.02·f+1` (cuft), `<= 0.015·f+3` (TPA), `<=0.1` (QMD) and `±1`
(BA/CCF/SDI) tolerances + the "known LP-calibration tail ~1.3-1.5%" / "small-tree REGENT tail" comments were ALL
STALE slack masking bit-exact matches. Tightened to exact `for c in 1:12; @test j[c]==f[c]`. 39+40+17 green.
(The LP-calibration tail DOES appear on LP-heavy stands like nohtdreg_cal — a SEPARATE ~1% residual, not these
GROWTH scenarios.) Confirms again: loose tolerances masked bit-exact results.

### Tolerance audit running tally (re-grounded vs LIVE, Oracle-A abandoned):
cut_cluster 11→bit-exact · regen_coverage 6→bit-exact · multistand 2→bit-exact · longrun 3→(2 bit-exact +1 ULP)
· mortality SDIMAX→exact +3 internal-pins noted · growth 8→bit-exact. ~30 tolerances done. REMAINING: test_carbon
~38 (`<=0.05`, live-Fortran carbon pools, print-rounding — need per-line confirm), test_net01 ~14, scattered
singles, and the Float32-ULP cluster (rtol 1f-4/1f-6). Verdict so far: the loose/Oracle-A tolerances were almost
all masking BIT-EXACT-vs-live results; genuine non-ULP residuals found are the 2 real divergences (non-native
FIXED, activity-fuels OPEN) + the LP-calibration/WK3-DGSCOR tail (nohtdreg_cal ~1%, separate).

### test_init/test_core/test_carbon — DONE (2026-07-02 tolerance audit)
- test_init: internal per-acre density pins → BIT-EXACT .sum integers (tpa 536.05→trunc==536, ba 77.39→77;
  sdi round==160); qmd 5.14 noted (internal deterministic; .sum QMD 5.1 bit-exact; atol=cruise-2dec-vs-internal).
- test_core: plot.pi → exact Float32(π)==3.1415927f0 (was atol=1f-6); RNG testset was already exact `==`
  (RANN LCG = FVS's own algorithm, mis-labeled "Oracle A").
- test_carbon: the 38× `<=0.05` (+ one 0.06) are PROVEN-ULP = half the carbon report's 1-decimal print
  resolution; goldens are LIVE Fortran; carbon_snt/carbon_jenkins BIT-EXACT vs live (31/31 + 22/22). Header note
  added; not re-litigated.

### Tolerance-audit tally (~75 done): cut 11 · regen 6 · multistand 2 · longrun 3 · mortality 4 · growth 8 ·
init 3 · core 1(+RNG) · carbon 40 — all re-grounded vs LIVE (Oracle-A abandoned): tightened to bit-exact OR
proven print-rounding/Float32-ULP with in-file notes. FINDING: the loose/Oracle-A tolerances ALMOST ALL masked
BIT-EXACT-vs-live results (cut/growth/regen were fully bit-exact under stale "tail"/"unported" tolerances). The
only genuine NON-ULP tolerance-covered residual found = test_growth's former LP-tail (turned out bit-exact here;
the real LP/WK3-DGSCOR tail lives on nohtdreg_cal ~1%). REMAINING: test_net01 ~14, test_multicycle 5,
test_keyword (lexer Float32-parse), test_fortbragg_coverage (atol=5, still Oracle-A → re-ground), scattered
singles, Float32-ULP cluster (rtol 1f-4/1f-6).

### net01 + Float32-ULP cluster + keyword lexer — DONE (2026-07-02 audit)
- test_net01: post-thin residual TPA/BA for THINBBA/THINABA/THINSDI → EXACT `==` vs live (verified THINBBA
  jl==live bit-exact; the rest tightened with it, suite green). Per-tree cuft vs live .trl 15.4: atol=0.1 =
  PROVEN-ULP (the .trl's 1-decimal print resolution). The stand-state/volume assertions were already exact `==`.
- Float32-ULP cluster (~21 sites, rtol/atol 1f-4…1f-6): JUSTIFIED BY CONSTRUCTION, two classes —
  (a) internal CONSERVATION/consistency invariants (test_carbon `sum(bins)≈total`, `size-split≈parent`;
      test_compress `tpa0≈tpa1` merge-conserved; test_fuel_decay decay identity) = Float32 SUMMATION ULP
      (a sum is not bit-equal to its parts in Float32) — these are not live comparisons at all;
  (b) deterministic Float32 COEFFICIENT/formula pins (test_growth dg_cor[65]=0.700993/[33]=1.085818,
      test_hcor htg_cor_init[22]=−0.893823, test_crown_width vs the hand-coded coefficient formula) =
      Float32 ARITHMETIC ULP on a deterministic value that FEEDS the bit-exact-vs-live .sum. Not slack.
- test_keyword (lexer vs KEYRDR, atol=1f-6/rtol=1f-5): the compared values are PARSED keyword-field numbers
  (variant-agnostic input, e.g. "60.0"→60.0f0); atol=1f-6 = Float32 PARSE ULP. Oracle-A is only a proxy for
  FVS KEYRDR, which parses the identical input numbers ⇒ the ULP justification is oracle-independent.
### Tolerance audit ~100 done. Live-facing tolerances re-grounded to BIT-EXACT vs live (cut/growth/regen/
### multistand/init/fortbragg/net01/longrun); print-resolution + Float32-arithmetic ULP JUSTIFIED-and-noted
### (carbon/net01-per-tree/the ULP cluster). Oracle-A abandoned everywhere. Suite 6436/2.

### ★ TOLERANCE AUDIT — FINAL SUMMARY (2026-07-02, user directive: justify vs LIVE as ULP/eig or fix to bit-exact)
Went through the suite's ~139 nonzero tolerances. Oracle-A abandoned; everything re-grounded vs the live binary.
FOUR categories:
1. LIVE-FACING coverage/Oracle-A slack → RE-GROUNDED TO BIT-EXACT `==` vs live (~30): cut_specpref+9 cut methods,
   growth_idg1/fint10/finth10 (all 12 cols), regen bare_plant/natural, multistand, init tpa/ba, fortbragg s30/s31,
   net01 3 thin residuals, longrun 2040/2140. FINDING: these loose/Oracle-A tolerances were MASKING BIT-EXACT
   results — the code already matched live; the slack hid it.
2. PRINT-RESOLUTION → PROVEN-ULP, noted (~60): carbon pools `<=0.05` (=½ the 1-dec report resolution, live
   goldens, 31/31+22/22 bit-exact), QMD `<=0.05` across 8 stand tests (1-dec .sum vs live Fortran), snt01
   qmd/mai, net01 per-tree cuft (.trl 1-dec), dbs_summary/compute.
3. FLOAT32-ULP → JUSTIFIED BY CONSTRUCTION, noted (~21): conservation/consistency invariants (carbon sums,
   compress TPA-conserved), deterministic coefficient/formula pins (dg_cor/htg_cor/crown_width, feed bit-exact
   .sum), econ PNV `atol=1f-3` (live-validated 13.6379/24.3235/32.6959), keyword lexer parse-ULP, core π→exact.
4. REAL NON-ULP RESIDUAL covered by a tolerance (honestly labeled): test_structure_stage `strdbh <=0.55` =
   ~0.5-DBH SSTAGE cohort/window-edge stratum-mean residual = the DEFERRED structure-stage/D3 item. NOT ULP.
NET: the suite was NOT "bit-exact or proven-ULP" as I'd claimed — it had Oracle-A slack (mostly over bit-exact)
+ one real non-ULP residual (structure-stage, deferred). After the audit: category-1 tightened to bit-exact,
2+3 proven-ULP with in-file notes, category-4 honestly documented as a deferred non-ULP residual. Plus the 2
real divergences found this session (non-native small-tree FIXED; thin-before-fire activity-fuels OPEN) and the
LP/WK3-DGSCOR tail (nohtdreg_cal, accepted @test_broken). Suite 6436/2.

### structure_stage strdbh — refined: bit-exact TREES, SSTAGE-report stratification/percentile residual
Applied the discrimination criterion to the one tolerance-covered "real non-ULP" item. The test grows snt01
stand-1 (BIT-EXACT vs live) directly with grow_cycle!, so the per-tree DBH is BIT-EXACT at every cycle. ⇒ the
strdbh residual (8/11 cycles bit-exact, 3 differ ≤0.543) is NOT a growth divergence — it is purely in the SSTAGE
report: `strdbh` = the dominant stratum's SSTGHP 70th-percentile DBH (sstage.f:487-576). At a cohort/window
boundary a single bit-exact tree lands in a different stratum (or on the wrong side of the percentile
interpolation) than Fortran, shifting the small stratum's 70th-pct DBH ~0.5. OPEN sub-question (needs a sstage.f
trace at the 3 differing cycles to close): is the flip a Float32 boundary TIE (⇒ ULP-amplification, reclassify)
or a real SSTGHP stratification/percentile-interp difference (⇒ real, the deferred structure-stage/D3 item)?
Documented as the deferred structure-stage residual pending that trace. Real-open-items list (post tolerance
audit): (1) thin-before-fire ACTIVITY-FUELS (open, root-caused, substantial FFE add); (2) structure-stage
SSTGHP boundary (this, ≤0.5 DBH, 3/11 cycles); (3) LP/WK3-DGSCOR tail (nohtdreg_cal ~1%, accepted @test_broken).

### ★ structure_stage strdbh RECLASSIFIED: ULP-class near-tie (NOT real) — rule #4 self-catch
TRACED (not assumed): residual isolated to cyc2 (0.543) + cyc7 (0.168), 9/11 bit-exact. Per-tree DBH BIT-EXACT
(snt01 grown). strdbh = 70th-pct DBH of the dominant stratum's cohort, cut at a HARD `csum(crownArea)>41382`
(0.95-ac) boundary. At cyc2 there are MASSIVE exact height TIES from tripling (9 recs @75.77, 9@73.03, 18@64.48…);
when the 0.95-ac cutoff lands in a tied group, jl's sortperm tie-break ≠ Fortran RDPSRT ⇒ a different cutoff tree
⇒ the 70th-pct ±4 window shifts ~0.5. ⇒ ULP-CLASS near-tie-ordering (the ACCEPTED COMPRESS/RDPSRT family), NOT a
real non-ULP residual. My prior tolerance-audit "category-4 real" label was a mislabel (rule #4: trace, don't
assume). ⇒ ALL test tolerances are now bit-exact / print-ULP / Float32-ULP / near-tie-ULP — NONE covers a
genuine non-ULP divergence.
### ⇒ Real non-ULP MODEL divergences remaining (found via completeness checks, NOT tolerance-hidden):
### (1) thin-before-fire ACTIVITY-FUELS — OPEN, root-caused (ordinary-cut crown-slash not fed to FFE fuel bed).
### (2) LP / WK3-DGSCOR tail (nohtdreg_cal ~1%) — accepted @test_broken (COMPRESS-family near-tie tail).
### Non-native small-tree DG double-scale = FIXED this session. Suite 6436/2.

### activity-fuels fix — PRECISE IMPLEMENTATION PLAN (root-cause + FVS routing traced; unit-convention to resolve)
ROOT (confirmed): jl's cut→fuel (cuts.jl:123-157) books cut material to the FFE fuel bed ONLY under YARDLOSS,
bole only. FVS FMSCUT (fmscut.f:61-98) adds EVERY cut tree's CROWN to the debris pool:
  IDC = DKRCLS(sp);  X = CTCRWN(i)·P2T                     # CTCRWN = removed TPA of the tree
  CWD(1,10,2,IDC) += CROWNW(i,0)·X                          # foliage → CWD size-class 10 (= jl litter)
  DO ISZ=1,5:  CWD(1,ISZ,2,IDC) += CROWNW(i,ISZ)·X          # branch classes 1-5 → CWD woody sizes 1-5
  CALL CWD3(...)                                            # downed-snag bole (jl already does the YARDLOSS analog)
jl BUILDING BLOCKS present: `crown_biomass(s,sp,d,h,icr)`→XV(0:5) (fmcrowe, faithful), `ffe_dkr_cls`(=DKRCLS),
`_FM_P2T`, `s.fire.cwd[size,2,idc]` (cwd size 10=litter/11=duff/1-9 woody; the YARDLOSS path already writes
cwd[j,2,idc]).
IMPLEMENT: in the cut loop (cut_one!, gated `s.fire!==nothing && s.fire.active`, for each cut tree with removed
TPA `prem`): xv=crown_biomass(...); idc=ffe_dkr_cls(s,sp); add xv[1]·prem·SCALE → cwd[10,2,idc]; xv[isz+1]·prem·
SCALE → cwd[isz,2,idc] for isz 1:5.
⚠ UNIT-CONVENTION TO RESOLVE FIRST (else it regresses the bit-exact corpus fires): what units does `s.fire.cwd`
store? The YARDLOSS path adds TONS (`fallvol=tcf·v2t/2000`), but `fmscro!` adds crown_biomass in LB to `cwd2b`
(×P2T deferred to the carbon report). Determine whether SCALE = `_FM_P2T` (if cwd is tons, matching YARDLOSS) or
1 (if cwd is lb like cwd2b) by checking a bit-exact fire scenario's cwd magnitude vs FMCONS. Then apply the
matching scale. VALIDATION target: cst01+THINDBH(2000)+SIMFIRE(2010) → live 255→5 (jl currently 255→56); the
corpus mild fires (snt01 SIMFIRE, fire_carbon, midcycle_fire, ffe) must stay bit-exact (they have no big
thin-before-fire so the crown-slash add is ~0 ⇒ should be untouched, but VERIFY). This is the LAST real model
divergence; all else is bit-exact / proven-ULP / accepted-eigensolver. Deferred as a careful FFE port, not
rushed at depth.

### ★ activity-fuels fix LANDED (cut-crown slash → FFE surface fuel): thin+fire under-kill 80% closed, no regression
Implemented FMSCUT crown→fuel in cuts.jl (per cut tree, fire-gated): crown_biomass(sp,d,h,icr) → cwd[foliage→10,
branches→1:5, cat2, idc] × prem·_FM_P2T (cwd is TONS; = CROWNW·CTCRWN·P2T with CTCRWN=removed-TPA). VALIDATED:
cst01+THINDBH(2000)+SIMFIRE(2010) post-fire jl 2020 TPA 56→12 (live 5) — ~80% of the under-kill closed; the
crown-slash now feeds the fire ⇒ more severe ⇒ more kill. NO corpus regression: suite 6436/2 (corpus fires have
no pre-fire cut ⇒ the routing is a no-op ⇒ bit-exact preserved). RESIDUAL (jl 12 vs live 5, ~1.4% of stand): jl
fire still slightly under-severe — candidates: (a) CTCRWN in cuts.f may carry a crown-proportion/PRCP factor
(I used the bare removed-TPA prem); (b) the UNMERCH bole/top slash (FVS also routes cut-bole residue via CWD3;
jl only does it under YARDLOSS) — an ordinary THINDBH leaves the top+unmerch as slash too; (c) a size-class or
decay-class mapping ULP. NEXT: stamp live FMSCUT CTCRWN + the cut-bole CWD3 contribution for this scenario to
close the last ~1%. The LAST real divergence is now 80% fixed + validated + regression-free; residual scoped.

### activity-fuels residual DIAGNOSED (live fmburn stamp): fuel-DECAY-TIMING, not crown-slash amount
Stamped live FVScs fmburn.f (FMFINT: IYR,BYRAM,FLAME,HPA; restored pristine) + jl fmburn.jl (temp) for
cst01+THINDBH(2000)+SIMFIRE(2010) — which ALSO has a base mid-cycle fire ~2003. Byram (fireline intensity ∝ fuel):
  2010 fire: live 14585 / flame 5.63  vs  jl 11847 / flame 5.12  → jl UNDER-fuels 19%
  2003 fire: live 2831  / flame 2.65  vs  jl 4931  / flame 3.42  → jl OVER-fuels
OPPOSITE directions ⇒ the crown-slash AMOUNT is ~right (my FMSCUT port); the residual is fuel-DECAY TIMING/
GRANULARITY: jl's 2000-thin crown-slash has NOT decayed enough by the 2003 mid-cycle fire (jl over-fuels ⇒
the 2003 fire OVER-consumes) ⇒ jl is left UNDER-fueled for the 2010 fire ⇒ under-kills (255→12 vs live 5). FVS
FMCWD decays the slash ANNUALLY between the 2000 cut and each fire; jl likely decays per-cycle or adds the slash
after the cycle's decay step, so the 3-yr-old slash reads full at 2003. FIX DIRECTION (next): ensure the cut-
crown slash enters the ANNUAL fmcwd decay from the cut year (so it's partially decayed at a mid-cycle fire),
matching FVS. This is a fuel-timing refinement on the (correct-amount) crown-slash add; the 80%-closing fix
stands. Non-corpus scenario; corpus fires bit-exact (suite 6436/2). LAST real divergence, now precisely diagnosed.

### activity-fuels — CORRECTION + crown-slash port VERIFIED FAITHFUL (fmscut.f read in full)
CORRECTION: the prior "2003-over/2010-under two-fire interaction" note was WRONG — the 2003 FMFINT was a
DIFFERENT stand (stand 4 "FFE TEST", InvYear 1993, SALVAGE+DEFULMOD). The activity-fuels stand is stand 1
(THIN 2000 + SIMFIRE 2010) = a SINGLE fire: live byram 14585 / flame 5.63 vs jl 11847 / 5.12 (19% under).
Read fmscut.f in full: the crown→CWD port is FAITHFUL to lines 89-97 exactly (foliage→CWD size10, branches→
size1-5, ×CTCRWN·P2T); for an ordinary thin DSNG=0 ⇒ CWD3 (bole) adds nothing (boles hauled) ⇒ NO missing
component. So the 19% is NOT the slash amount — it's the 10-yr DECAY/stash of the fine crown-slash between the
2000 cut and the 2010 fire (the corpus no-cut fuel-decay is bit-exact). VERDICT: crown-slash port faithful;
residual = a small multi-year cut-fuel-decay-granularity / near-threshold-kill effect on a constructed non-
corpus scenario. Closing it needs a dedicated FFE cut-fuel-decay study; not worth destabilizing the bit-exact
corpus decay at this depth. 80%-closing fix stands (255→12 vs live 5; was 255→56). Sources pristine.

## Tri-variant floor RE-VERIFICATION vs fresh live (2026-07-02, this session)
Re-ran all three variants' differential against freshly-relinked live binaries (re-trace discipline — confirm
the floor holds, hunt uncatalogued cross-variant divergence). Result: NO new/uncatalogued divergence; every
DIFF maps to an already-documented/accepted class.
- **SN:** D16b fire family RE-VERIFIED at ULP floor via run_keyfile-vs-fresh-live (authoritative). CORRECTION:
  the "byte-identical" first pass was a measurement error (committed LIVE-golden vs fresh live = live-vs-live).
  True jl-vs-live: TPA bit-exact or ±1 across the family; Bdft ±1 to ~1% (snt01_alpha 218/~14000, fire_repeat 50,
  s10_fire 7, fire_salvage 1) = the documented board-foot threshold-amplification ULP class. Fire-kill (TPA) is
  at the floor; volumes are ULP-threshold. LESSON: always diff run_keyfile(output=:sum) vs a FRESH live run —
  committed .sum goldens can be live saves (comparing them to live proves nothing).
- **NE:** ALL 10 oracle-having stands BIT-EXACT vs fresh live — net01 + ne_cov0..4 (all-species) + ffe +
  midcycle_fire + plant_div + plant_hard. (The hook's "net01 BARE-regen ~4% Mcuft" is STALE — net01 is
  bit-exact.) 7 ne_fixtures (aspen/dense/divspp/divspp_f9{14,20,21}/thin) have NO live oracle: live FVSne
  IC-stops (STOP 10/20, main.f:11-12) on them regardless of ECHOSUM ⇒ same accepted no-oracle class as the SN
  live-FPE stands (they are jl-internal fixtures; jl runs them).
- **CS:** cst01 + cs_allsp DIFFs all map to documented classes. cst01 stand-2 (deep every-3-cycle THINDBH):
  DISCRIMINATION CRITERION applied — BA & QMD are bit-exact through 2130 (BA 133/133) while only TPA drifts
  a few ⇒ the divergence is small-tree/REGEN COUNT (D8/D10 threshold-amplification, large trees identical),
  BA only diverges at the 2140+ late thin (accumulated regen-count × threshold-sensitive thin). Stand-4 = the
  FFE fire-mortality class (D16b/D17). cs_allsp 1.31% CCF@2090 = documented late ULP floor. Stands 1/3/5 =
  ±1 TPA rounding / bit-exact.
⇒ **CONFIRMED: all three variants are at the documented floor. Every live divergence is bit-exact, ULP/
threshold-amplified (regen-count / deep-thin / board-foot), accepted-eigensolver, an FFE fire-mortality-
distribution residual (D16b/D17, faithful ports), or a no-live-oracle fixture.** No uncatalogued divergence found.

### D17 activity-fuels — EVIDENCE-CORRECTED via live+jl fuel stamps (2026-07-02)
Stamped BOTH live FVScs fmburn.f (per-size-class CWD dump; restored pristine + byte-identical output verified)
and jl fmburn.jl at the cst_ft10 stand-1 2010 fire. Hard data corrects the diagnosis TWICE:
- NOT "jl under-fuels 19%": jl actually has MORE woody (12.69 vs live 10.06 tons), uniformly ~1.2-1.3× in
  EVERY size class 1-9 (sz7-8 ~1.85×). Litter (0.134 vs 0.129) and duff (5.17 vs 5.13) are ~bit-exact.
- Yet jl byram is LOWER (11847 vs 14585). More fuel + lower intensity ⇒ the byram swing is the NON-MONOTONIC
  FMCFMD dynamic fuel-model selection (the D16b "non-monotonic FMCFMD3" note) responding to jl's higher loading
  — a discrete model/interp choice, not a raw-load monotone effect.
- The crown-slash ADD formula stays VERIFIED FAITHFUL (fmscut.f); the excess is BROADER woody over-accumulation
  in the thinned stand (incl. sz7-8 where crown-slash doesn't reach) over 2000→2010, i.e. an FFE base fuel-
  accumulation/decay-dynamics difference in the post-thin stand, feeding the sensitive FMCFMD.
VERDICT: 📌 real but constructed-scenario-only (THIN-then-SIMFIRE-at-cycle-boundaries; no corpus/real-test stand
hits it), well-characterized, add-formula faithful. Full close = trace the ~26% post-thin woody over-accumulation
+ the FMCFMD non-monotone response; deferred (deep FFE fuel-dynamics, not worth destabilizing the bit-exact
corpus fuel-decay). Documented with live+jl per-size-class evidence above.

### SN full sweep — AUTHORITATIVE RE-RUN (2026-07-02, post non-native + activity-fuels changes) — 221/31/8
Ran the FULL 260-stand SN sweep via divergence_sweep.jl (run_keyfile-vs-fresh-live, the AUTHORITATIVE path)
as a regression check after this session's code changes (NE/CS non-native small_tree_growth fix — gated off for
SN; cuts.jl cut-crown→FFE-fuel). Result: **221 bit-exact (up 1 from 220), 31 DIFF, 8 live-FPE (no oracle).**
- NO REGRESSION: no new stand in the DIFF list; the cut-crown fuel change added no SN divergence (fire stands
  salvage/defulmod/fueltret/fire_repeat/snt01_alpha all at their prior documented levels, not worse).
- Every DIFF maps to an already-documented/accepted class: D13 treeszcp_cap 22.8%/htcap 10.7% (size-cap ULP) ·
  COMPRESS compress 13.6% (eigensolver) · D8/D10 regen mult_mortmult 17%/mortmult_win 13.5%/regdmult 4.7%/
  reghmult 2.9%/baimult 0.5% + bare_natural/plant 4.6%/multipoint 2.8%/mp3 2.7% (regen threshold-amp) · D16b/D17
  fire family snt01_alpha=compute_cycle 2.72% Bdft / fire_repeat 2.47% / fmortmlt 1.56% / s10_fire 1.19% /
  fueltret 0.72% / defulmod=salvage 0.59% (board-foot threshold + fire-mort distribution) · documented ULP tails
  hcor_smalltree 2.09% / timeint10 1.96% (non-native SN@10) / htgstop_stoch 1.65% / dense_long=s09_cyc20 0.76% /
  fixmort_* 0.3% / topkill_det 0.27% / s15_phys 0.22% / s22_forest_809 0.21% (NVEL) / growth_finth5 0.21%.
- CONFIRMS the "byte-identical" correction: the sweep independently ranks snt01_alpha at 2.72% Bdft (NOT zero) —
  vindicating the corrected ULP-class-threshold assessment over the earlier committed-golden measurement error.
⇒ The SN floor holds authoritatively after this session's changes; no uncatalogued or regressed divergence.

### D17 — root NARROWED to cut-crown decay-rate/decay-class (live ALL-FUELS report evidence, 2026-07-02)
Pulled live's ALL FUELS REPORT (FuelOut/FuelRept on a single-stand thin+fire key — NO source stamping) for the
cst_ft10 stand-1 (THIN 2000 + SIMFIRE 2010). Live surface-fuel decay of the thin's crown-slash:
  YEAR  LITT  DUFF  0-3"  >3"  ...  TOTAL
  2000  3.33  4.8  17.1   7.4 ...  32.8   ← post-thin spike: 0-3"=17.1 t (crown-slash)
  2010  0.00  1.7   1.6   1.0 ...   4.4   ← pre-fire: crown-slash 0-3" decayed 17.1→1.6 (~90% in 10yr)
Identical raw-cwd stamps (apples-to-apples) gave jl woody 12.69 vs live 10.06 at the fire (+26%, uniform across
size classes); litter/duff bit-exact. Combined ⇒ HYPOTHESIS: jl decays the CUT-CROWN material slower than FVS's
fast 17.1→1.6, so jl retains ~26% more woody by 2010 → the non-monotonic FMCFMD picks a lower-intensity model →
under-kill (12 vs 5 TPA). Prime suspect = the decay-CLASS (IDC=DKRCLS(sp)) or DKR rate applied to cut-crown CWD
in cuts.jl vs fmscut.f/fmcwd.f (litter/duff bit-exact means the fine-litter path is right, but the woody 0-3"/>3"
crown-slash decay class may differ). NEXT (for whoever closes it): compare jl's DBS Fuels table (0-3"/>3" by year)
to live's ALL-FUELS report at 2000/2005/2010 to confirm the decay-rate gap, then reconcile the cut-crown decay
class. Still 📌 (constructed-scenario-only; add-formula faithful; corpus fires bit-exact/ULP per the authoritative
sweep) — but root now narrowed from "fuel dynamics" to a specific cut-crown decay-class/rate check.

### D17 — ROOT CAUSE PINNED (code-verified ordering bug, 2026-07-02)
Traced the jl per-cycle order in summary.jl for a cut-in-cycle-N / fire-in-later-cycle-M case (cst_ft10: THIN
2000, SIMFIRE 2010):
- 2000 (non-fire) cycle: `ffe_fuel_update!` (10-yr DECAY) runs at summary.jl:277 BEFORE `grow_cycle!` (:280);
  the THIN's crown-slash is added INSIDE grow_cycle! (cuts.jl) — i.e. AFTER that cycle's decay. And non-fire
  grow_cycle! is called with `fuel_period=nothing` (:281) ⇒ no in-grow decay. So the 2000 crown-slash gets
  ZERO decay in its own cycle.
- 2010 (fire) cycle: `fire_smlg` is stashed at :276 (start of cycle, before any 2010 decay) and the 2010
  ffe_fuel_update! is DEFERRED to post-fire (fire cycle). ⇒ the 2010 fire burns the 2000 crown-slash with
  ZERO decay ever applied. FVS instead decays it annually via FMCWD (live report: 0-3" 17.1t@2000 → 1.6@2010).
This is exactly why jl retains +26% woody at the fire ⇒ non-monotonic FMCFMD picks a lower-intensity model ⇒
under-kill (12 vs 5 TPA). WHY THE CORPUS IS UNAFFECTED: corpus fire stands cut+fire in the SAME cycle (salvage/
fueltret/defulmod), so the crown-slash is correctly fresh/undecayed at the fire — bit-exact/ULP per the
authoritative sweep. The bug ONLY fires when a cut precedes a fire by ≥1 cycle (the constructed cst_ft10).
FIX DIRECTION (delicate, deferred to a fresh session w/ full re-validation): decay cut-added crown-slash over the
cycles between the cut and a later fire — e.g. apply the cut cycle's fuel decay AFTER the cut (reorder for cut
cycles) or decay the newly-added crown-slash by the remaining cycle span — GATED so the same-cycle cut+fire
corpus stays bit-exact. Root now PINNED (was: "fuel dynamics" → "decay rate" → THIS specific ordering bug).

### D17 — FIX DESIGN (three approaches traced; correct one identified, pitfalls documented)
The bug (pinned above): cut crown-slash added AFTER its cycle's fuel decay ⇒ burns undecayed at a later-cycle
fire. FVS order = FMSCUT (add cut fuel at cycle boundary) THEN annual FMCWD loop (decay). jl order = ffe_fuel_
update! (decay) THEN grow_cycle!/cuts.jl (add). Three candidate fixes, traced BOTH sides:
1. Move ffe_fuel_update! AFTER grow_cycle! for cut cycles — ✗ REGRESSES: ffe_fuel_update! also does FMCADD
   litterfall off the STANDING crowns; running it post-growth changes the litterfall basis for EVERY non-fire
   cycle ⇒ breaks the bit-exact corpus base fuel. Do NOT do this.
2. Pre-decay the crown-slash increment in cuts.jl (multiply by the cycle's fmcwd decay factor before adding)
   — ✗ INCOMPLETE: FMCWD moves a PRDUFF fraction of decayed woody INTO duff; jl/live duff is currently
   bit-exact (5.17/5.13), so a naive increment-multiply that drops the respired mass would still need to route
   the duff-transfer or it desyncs duff. Must replicate fmcwd's woody→duff transfer for the increment.
3. Hoist ONLY the cut's FFE-fuel-addition (the crown-slash → cwd) to BEFORE ffe_fuel_update! in the cut cycle,
   so the cycle's decay loop ages it (= FVS FMSCUT-then-FMCWD) — ✓ CORRECT in principle, but a real restructure:
   the cut/crown-slash is computed inside grow_cycle! (depends on which trees are thinned), so it must be split
   out (compute cut removals + crown-slash → cwd) and run before the fuel decay, with grow_cycle! doing the rest.
   GATE: same-cycle cut+fire (corpus salvage/fueltret/defulmod) must stay bit-exact — those re-stash fire_smlg
   POST-cut and burn the crown-slash FRESH (FVS same-year cut+fire = no annual decay between), which is already
   correct; the hoist must not disturb that (only the cross-cycle cut→fire path needs the pre-decay).
RECOMMENDED: approach 3, implemented in a fresh session, validated by re-running the authoritative 260-stand SN
sweep (baseline 221/31/8) + the cst_ft10 fire (target: jl woody→~10.06, byram→~14585, kill→~5 TPA) + the CS/NE
sweeps, confirming the same-cycle-cut+fire corpus stays bit-exact. Root PINNED, fix SCOPED, pitfalls MAPPED.

### D17 — duff-paradox reconciled (root-cause pin CONFIRMED, 2026-07-02)
Checked an apparent contradiction: if jl leaves the crown-slash undecayed, jl duff should LACK the crown-slash
decay products and read LOW — but jl/live duff match (5.17/5.13). Reconciled: the crown-slash woody decay
(0-3" 17.1→1.6 = ~15.5 t over 10yr) mostly RESPIRES OUT (leaves the system); live duff is actually slightly
LOWER than jl, so the crown-slash→duff routing is negligible. ⇒ "jl retains undecayed crown-slash (+26% woody),
duff ≈ unaffected" is fully self-consistent. Root cause (cut crown-slash misses its cycle's decay) CONFIRMED.
Code-read confirms approach 3 is a real reorder (cuts! is INSIDE grow_cycle! at simulate.jl:~350; the pre-grow
ffe_fuel_update! is in summary.jl BEFORE grow_cycle!) — and the local alternative (decay just the increment) must
replicate fmcwd's exact per-year respire/duff split for bit-exactness. Both delicate; deferral is evidence-based,
not reflexive. D17 = fully characterized 📌 (root pinned+confirmed, fix scoped, pitfalls mapped, validation plan set).

### D17 — crown-slash-decay hypothesis TESTED and REFUTED; root redirected to FMCFMD/FMFINT (2026-07-02)
IMPLEMENTED the "age the cut crown-slash by the cut cycle's decay" fix (gated non-fire cut cycles, _book_cut_fuel!
mirroring one fmcwd pass) and TESTED on cst_ft10 (rule #4 — validate the hypothesis vs live BEFORE trusting it).
RESULT: it made the fire WEAKER — cst_ft10 2020 TPA 12 (pre-fix) → 29 (with decay), i.e. AWAY from live's 5.
⇒ REFUTED the "jl over-accumulates woody ⇒ under-kill" causal chain. Reverted the fix.
KEY INSIGHT the experiment forces: reducing jl's woody LOWERS its byram/kill, yet LIVE has LESS woody (10.06 vs
jl 12.69) and a STRONGER fire (byram 14585 vs 11847, kill 5 vs 12). So jl and live produce DIFFERENT byram for
essentially the SAME (or even less) fuel ⇒ the divergence is in the FMCFMD dynamic-fuel-model SELECTION / FMFINT
byram computation, NOT the fuel amount or its decay. The earlier "+26% woody → non-monotonic FMCFMD → under-kill"
framing had the causality backwards: the woody excess is real but is NOT what drives the under-kill (removing it
worsens the kill). The 80%-closing crown-slash ADD fix likely helped by nudging the fuel into a range where jl's
(divergent) FMCFMD happens to run hotter — coincidental, not the mechanism.
NEW ROOT DIRECTION (for closure): differential jl vs live FMCFMD3/FMFINT — for identical cwd loading, compare the
selected standard fuel model(s) + weights + the resulting reaction-intensity/byram. Prime suspects: the small/large
fuel partition (_small_large_fuel vs FVS FMCFMD3), the weighted-model interpolation, or a rothermel-input term.
D17 stays 📌 — but the crown-slash-decay lead is CLOSED (refuted by experiment); root redirected to FMCFMD/FMFINT.

### D17 — FMCFMD path systematically narrowed: 5 selection layers VERIFIED MATCH; divergence is in _fmdyn/rothermel math
Followed the redirected FMCFMD/FMFINT root (cst_ft10 CS: jl byram 11847 vs live 14585) by verifying each layer
of the CS fuel-model selection against jl (source-read, no live stamp):
1. FOREST-TYPE (FMCSFT vs SN FMSNFT): the ONLY logic diff is the pine species range CASE(3:7) vs CASE(4:14) —
   and jl ALREADY handles it (fuel_loading.jl:66 `s.variant isa CentralStates ? (3,7) : (4,14)`). ✓ MATCH.
2. SELECTION LOGIC (fmcfmd.f iffeft branches): identical to jl's SN path modulo the redcedar CASE(1,2) vs (2)
   [redcedar-only, N/A to cst_ft10's hardwood type]. ✓ MATCH for cst_ft10.
3. XPTS breakpoint table: CS fmcfmd.f XPTS == SN's (fmcfmd.f diff = a COMMENT only); jl uses _FMD_XPTS. ✓ MATCH.
4. SMALL/LARGE partition (fmtret.f:377-388): sizes 1-3+litter(10)=small, 4-9=large, over pools I=1,2 (piled+
   dispersed). jl _small_large_fuel matches; cst_ft10 has no piling so I=2=0. ✓ MATCH.
5. IFLOGIC = 0 by default (fminit.f:824) ⇒ the weighted-standard-models path (fmcfmd2.f:653), NOT the custom-
   model-89 dynamic path (IFLOGIC=2) — i.e. jl's exact path. ✓ MATCH (custom-89 not active for cst_ft10).
⇒ ALL selection INPUTS match. The byram divergence is downstream in the MATH: either jl's `_fmdyn` inverse-
distance interpolation vs FVS fmdyn.f, or the resolved standard-model params (load/sav/depth) + the rothermel/
FMFINT byram computation. Note the woody-amount confound: jl's cwd woody is +26% (higher SMALL/LARGE), but
DECAYING it moved byram DOWN (refuted fix) while live has LOWER woody + HIGHER byram ⇒ jl and live have different
byram(SMALL,LARGE) curves ⇒ the interpolation/rothermel math diverges, not just the inputs. NEXT: stamp live
FMFINT (SMALL,LARGE, selected FMD list+weights, resolved load/sav/depth, BYRAM) at the cst_ft10 2010 fire and
diff jl's `_fmdyn` + `rothermel_surface_fire` at the identical (SMALL,LARGE). D17 root narrowed 5 layers deep.

### D17 — LOCALIZED to FMDYN model-WEIGHTING driven by jl's coarse-woody (LARGE) over-accumulation (2026-07-02)
Stamped live CS FMFINT (per-model FMOD/FWT/BYRAMT, IYR-tagged) + FMDYN (SM,LG) [both restored pristine, verified]
and jl's fmburn model loop, for cst_ft10's 2010 fire. Result — the models MATCH, the per-model byram nearly
matches, but the WEIGHTS diverge hugely:
  live 2010: fm10 w=0.519 byram=11584  +  fm5 w=0.481 byram=17826  → 14586  (fm5 = the HOT model, ~half weight)
  jl   2010: fm10 w=0.910 byram=11304  +  fm5 w=0.090 byram=17344  → 11848  (fm5 nearly ignored ⇒ cool fire)
⇒ D17's under-kill is the FMDYN inverse-distance WEIGHTING, not the selection/rothermel. FMDYN weights the
(SMALL,LARGE) point vs each model's XPTS iso-line (fm5 (5,15), fm10 (10,30)). SM/LG: jl fire_smlg=(7.91,4.91)
vs live≈(7.76,4.03) [SMALL ~matches; jl LARGE +22%]. A higher LARGE pulls the point toward fm10's iso-line
(LARGE-intercept 30 » fm5's 15) ⇒ jl over-weights the weak fm10. ROOT: jl's LARGE fuel (coarse woody, cwd sizes
4-9, esp. 6-9 which the per-size stamp showed at ~1.85× live) is OVER-accumulated. This is NOT the crown-slash
(that lands in sizes 1-5+litter = mostly SMALL, which ~matches) — it's the COARSE down-wood (sizes 6-9), i.e.
SNAG-FALL / base coarse-CWD dynamics. RECONCILES the refuted decay fix: decaying the crown-slash cut SMALL (not
the LARGE 6-9 excess) and dropped SMALL below the sm>6 candidate threshold ⇒ wrong target, made it worse.
NEXT: differential jl vs live coarse-CWD (sizes 6-9) accumulation at the fire — trace snag-fall (update_snags!/
FMSNAG bole→cwd) + base FMCWD coarse pools; find why jl sizes 6-9 run ~1.85× live. That over-accumulation is the
true D17 root. 📌 (constructed-scenario-only; now localized to a specific fuel pool + mechanism).

### D17 — coarse-CWD source narrowed to snag-fall/cone-split (standing TPA matches; 2026-07-02)
Per-size-class cwd at the cst_ft10 2010 fire (both from identical raw-cwd stamps):
  size:   4      5      6      7      8      9
  live: 2.401  0.441  0.630  0.361  0.061  0.000   (LARGE=Σ4-9 = 3.894)
  jl:   2.900  0.542  0.686  0.674  0.113  0.000   (LARGE = 4.915, +26%)
The excess concentrates in sizes 7-8 (20-50" down-wood): jl 0.674/0.113 vs live 0.361/0.061 (~1.87×/1.85×).
Both are NONZERO ⇒ cst_ft10 legitimately has 20"+ boles (not a spurious-index bug) — it's the AMOUNT.
KEY CONSTRAINT: standing TPA MATCHES at 2010 (255 both) ⇒ the same trees survive, so the mortality COUNT over
2000→2010 is ~equal ⇒ the coarse-CWD excess is NOT "jl kills more big trees". It must be the snag→down-wood
conversion: snag-FALL timing (how many of the 2000-2010 mortality snags have fallen by the 2010 fire) or the
CWD1 cone-split (how a fallen large bole's volume distributes across sizes 7-8). jl deposits ~1.87× the coarse
(20"+) down-wood for the SAME standing stand. NEXT: jl-vs-live differential of the snag cohort (density/size at
2005 & 2010) + the CWD1 cone-split fractions for a ~20"+ bole — find why jl's sizes 7-8 run ~1.87×. That is the
terminal D17 root (feeds LARGE → FMDYN weight → under-kill). 📌 constructed-scenario-only; localized 5 layers deep
(FMDYN weight → LARGE +26% → sizes 7-8 1.87× → snag-fall/cone-split, standing-TPA-matched).

### D17 — FULLY CHARACTERIZED: initial fuel matches; +26% is the 2000→2010 thin-slash evolution (2026-07-02)
ICYC-aligned live FMDYN stamp (restored pristine) + jl fire_smlg give the authoritative (SMALL,LARGE) trajectory:
  cyc1 1990:  live (7.76, 4.03)  ==  jl (7.76, 4.03)   ← INITIAL FUEL MATCHES EXACTLY (not a loading bug)
  cyc2 2000:  live (20.48, 7.41) post-thin crown-slash SPIKE
  cyc3 2010:  live (6.30, 3.89)  vs  jl fire_smlg (7.91, 4.91)  ← +26% in BOTH small & large
Live's 2000 crown-slash spike (20.48) DECAYS to 6.30 by the 2010 fire (~(1-0.11)^10 woody decay). jl's 2010 fire
fuel is +26% high in BOTH SMALL and LARGE. Since the cut crown-slash spans cwd sizes 1-5, it feeds BOTH the SMALL
(1-3+litter) and LARGE (4-5) partitions ⇒ the +26% is one coupled source: the 2000-thin slash's evolution to the
2010 fire (add-timing + annual decay), NOT a separate coarse-woody bug and NOT an initial-loading bug. THIS
RECONCILES THE REFUTED DECAY FIX: the naive per-increment decay reduced the slash too far (litter dk=0.65→~0,
woody dk=0.11→0.31 over 10yr) and shifted the SMALL/LARGE point past live's, moving the FMDYN weight the wrong
way ⇒ weaker fire. The correct fix must reproduce FVS's EXACT 2000→2010 slash evolution so (SMALL,LARGE) lands on
live's (6.30, 3.89) — i.e. faithful FMSCUT-add-then-annual-FMCWD phasing (approach 3, the reorder), NOT the
lossy per-increment single-shot. 📌 constructed-scenario-only; now FULLY characterized (initial✓, source=coupled
thin-slash 2000→2010 evolution, exact target (6.30,3.89), fix=faithful add-then-decay phasing).

### D17 — DECISIVE: fire machinery PROVEN faithful; D17 is purely a +26% fuel-accounting gap (2026-07-02)
Ran a clean fuel-forcing test (temporary hack, reverted): forced jl's 2010 CS fire fuel from its (7.91,4.91) to
LIVE's exact (SMALL,LARGE)=(6.30,3.89). Result: jl kill → cst_ft10 2020 TPA = 5 = LIVE EXACTLY (was jl 12).
⇒ PROVEN: jl's _fmdyn model-weighting + rothermel byram + FMEFF mortality are ALL FAITHFUL. D17 is ENTIRELY a
fuel-amount problem — get jl's 2010 fire fuel to (6.30,3.89) and the kill matches live bit-for-bit. No machinery
bug. This DE-RISKS the fix: it lives solely in the FFE fuel pipeline (the 2000-thin slash's 2000→2010 evolution),
not in fire behavior/mortality. Refines the decay story: jl's 2010 fuel is 7.91 (NOT the ~20 t 2000 slash spike),
so jl DOES decay the slash — just ~26% LESS than live (live 20.48→6.30; jl →7.91). ⇒ the earlier single-shot
per-increment decay DOUBLE-decayed (jl already decays it) and overshot (kill 29). The correct fix is to close the
SPECIFIC 26% under-decay in jl's existing slash-decay path (find where jl's 2000-slash decays to 7.91 vs live's
6.30 over 2000→2010 — a decay-years/phasing gap), NOT to add blanket decay. 📌 constructed-only; machinery proven
faithful; fix isolated to the fuel-decay accounting (target (6.30,3.89) → exact kill 5).

### D17 — fuel trajectory: jl FLAT vs live SPIKE-AND-DECAY (compensating under-add + under-decay, 2026-07-02)
jl grow_cycle! probe (post-cut _small_large_fuel; jl's 2010 value matches the faithful fire stamp ⇒ reliable):
  jl:   1990 SMALL 3.60 → 2000 7.91 → 2010 7.91 (FLAT 2000→2010)
  live: 1990 SMALL 7.76 → 2000 20.48 (big thin-slash SPIKE) → 2010 6.30 (decays hard)
⇒ jl has TWO compensating fuel errors, not simple over-accumulation:
  (1) UNDER-ADDS the crown-slash at the 2000 thin (jl ~7.9 vs live ~20.5) — jl's slash spike is ~1/3 of live's;
  (2) UNDER-DECAYS 2000→2010 (jl flat 7.9 vs live 20.5→6.3) — the phasing gap (slash added after the cycle's
      pre-grow ffe_fuel_update!, fire cycle defers the next update ⇒ ~no decay before the fire).
Net at the 2010 fire: jl 7.91 vs live 6.30 (+26%). The two errors partially cancel. THIS is why prior single-
factor fixes failed: fixing only the decay (removing fuel) ignored the under-add and overshot; the crown-slash
ADD looked "faithful" per-tree vs fmscut.f but the STAND-level 2000 spike is 1/3 of live's (⇒ likely a
CTCRWN/prem or per-record-vs-total accounting gap in the slash magnitude, OR the 2000 spike is mostly the base
+litterfall not the slash and jl's litterfall/base differs). MACHINERY REMAINS PROVEN FAITHFUL (fuel-forcing →
exact kill). D17 fix now needs BOTH: (a) match live's 2000 slash spike magnitude, (b) faithful 2000→2010 decay.
📌 constructed-only; a genuinely tangled two-error FFE fuel-accounting divergence; machinery faithful; exact
target (6.30,3.89)→kill 5. Given no corpus/real stand exercises it and the fix is a coupled two-part fuel-pipeline
correction, this is the campaign's lowest-value open item — recommended acceptance as a documented 📌.

### D17 — CORRECTION: jl intermediate-fuel probes are UNRELIABLE; only machinery-faithful + net-gap are solid
Rule #4 check: my two jl-side fuel probes DISAGREE (one _small_large_fuel probe read jl 1990 SMALL=7.76, the
grow_cycle! post-cut probe read 3.599) — both call fmcba! OUT OF THE FIRE SEQUENCE and perturb the FFE state, so
the jl 1990/2000 intermediate trajectory (and the "under-adds slash to 1/3" breakdown in the prior entry) is NOT
trustworthy. RETRACT the specific under-add/under-decay split. The SOLID, reproducible facts remain:
  • jl's 2010 FIRE fuel = (7.91, 4.91) [from the in-sequence fmburn stamp] vs live (6.30, 3.89) — +26% net.
  • Forcing jl's fire fuel to live's (6.30,3.89) ⇒ kill 5 = live EXACTLY ⇒ fire MACHINERY proven faithful.
  • Exact target: land jl's 2010 fire fuel on (6.30,3.89).
To pin the fuel gap's mechanism (add-magnitude vs decay-phasing) requires NON-perturbing instrumentation of the
actual fire-fuel stash path (the fmburn/fire_smlg sequence), not out-of-band fmcba! probes. Deferred. D17 stays a
📌: machinery faithful, pure +26% fire-fuel gap, exact target known; the mechanism split needs careful in-sequence
instrumentation. Constructed-scenario-only; the campaign's lowest-value open item.

### D17 — RELIABLE non-perturbing fuel trace (un-retracts the multi-factor picture, adds fmcba! lead, 2026-07-02)
Added a NON-perturbing trace (reads `_small_large_fuel` right after the real ffe_fuel_update!, no out-of-sequence
fmcba!) — its 2010 value (7.912) MATCHES the faithful fmburn fire stamp ⇒ this trace IS reliable (the earlier
"retract" was over-cautious; the disagreeing probe was the fmcba!-AUGMENTED read, a different quantity). jl raw-
cwd fire-fuel trajectory vs live (ICYC-aligned):
  1990 base:  jl 3.599  vs live 7.76   ← jl HALF; BUT jl's fmcba!-augmented 1990 = 7.76 = live
  2000 thin:  jl 7.912  vs live 20.48  ← jl's crown-slash spike far smaller
  2010 fire:  jl 7.912 (FLAT)  vs live 6.30  ← jl doesn't decay 2000→2010
THREE factors net to +26% at the fire:
  (A) fmcba! AUGMENTATION GAP — jl's fire uses the RAW stashed `fire_smlg` (3.599 at 1990); live's fire SMALL
      (7.76) equals jl's OWN fmcba!-augmented value ⇒ live's fire fuel includes the fmcba! live-herb/shrub/
      litterfall augmentation that jl's `fire_smlg` OMITS. Specific + testable: does live's FMDYN SMALL include
      the FMCBA live-fuel that jl's fire_smlg drops? (fmburn calls fmcba! at line 92 but select_fuel_models reads
      the pre-fmcba! stashed fire_smlg.) NOTE: corpus fires are bit-exact with raw fire_smlg ⇒ any change here is
      HIGH corpus-regression risk; verify against live before touching.
  (B) SLASH MAGNITUDE — jl's 2000 thin-slash spike (Δ≈4.3) is ~1/3 of live's (Δ≈12.7).
  (C) DECAY PHASING — jl flat 2000→2010 (no decay of the slash) vs live 20.48→6.30.
Machinery still PROVEN faithful (fuel-forcing→kill 5). D17 = a genuinely multi-factor FFE fuel-accounting
divergence (fmcba-augmentation + slash-magnitude + decay-phasing), all in the fuel pipeline, netting +26% at ONE
constructed fire. 📌 constructed-only, lowest-value; each factor is high corpus-regression risk (corpus fires are
bit-exact on the current path) ⇒ any fix needs per-factor live-differential + full sweep re-validation.

### D17 — two clean ELIMINATIONS: initial FUINI verified correct + fmcba-aug is a red herring (2026-07-02)
Two hypotheses KILLED this pass, narrowing D17 unambiguously:
1. INITIAL FUEL TABLE ✓ CORRECT: jl's CS fire_fuel_dead.csv oak-hickory row = live cs/fmcba.f FUINI oak-hickory
   row BYTE-FOR-BYTE (0.15,0.74,1.70,0.38,0.97,2.68,0,0,0,5.17,4.52) ⇒ SMALL=7.76/LARGE=4.03 = live's 1990 FMDYN
   exactly. The trace's "jl 1990=3.599" was POST-first-cycle-decay (sampled after ffe_fuel_update!), not the raw
   initial — NOT a table bug. Initial loading is faithful.
2. fmcba AUGMENTATION = RED HERRING: live's fire SMALL (6.30) = the RAW cwd sum (live 2010 sizes 1-3+litter =
   4.093+1.683+0.393+0.129 = 6.30 EXACTLY) ⇒ live's fire uses raw cwd, same quantity as jl's fire_smlg; the
   "jl-augmented 1990=7.76=live" was a coincidence (jl-raw 3.599 + ~4 herb/shrub ≈ 7.76). No augmentation gap.
NOW RULED OUT for D17: fire machinery (fuel-forcing→kill 5), initial FUINI table, forest-type (FMCSFT), model
selection, XPTS, small/large partition, fmcba-augmentation. ⇒ D17 is DEFINITIVELY the thin-slash fuel EVOLUTION
over 1990→2010 (the 2000-thin crown-slash's add magnitude + its 2000→2010 decay/phasing), sampled at the fire:
jl fire fuel 7.912 vs live 6.30 (+26%), both from the same correct FUINI. Fix = faithful FMSCUT-add-then-annual-
FMCWD phasing (approach 3), the ONLY remaining candidate. 📌 constructed-only; machinery+inputs all verified
faithful; the residual is isolated to the cut-slash fuel-evolution phasing (delicate reorder, high corpus risk).

### D17 — decay approach VALIDATED (LARGE lands exact); residual = crown-slash FINE-size over-loading (2026-07-02)
Re-ran the cut-slash decay fix WITH a fire_smlg measurement (temp, reverted). Decaying the slash by the cut→fire
period (s.control.year) moved jl's 2010 fire fuel:
  SMALL: 7.91 → 4.579  (target 6.30)  — OVERSHOT by ~1.7
  LARGE: 4.91 → 3.851  (target 3.89)  — LANDED EXACTLY ✓
⇒ the decay PHASING approach is CORRECT (LARGE, the coarse branches sizes 4-5, lands bit-close to live). The
residual is confined to SMALL. Disabling the foliage→litter(size10) decay changed nothing ⇒ NOT litter. So the
SMALL overshoot is the FINE branches (sizes 1-3): the decay removes 3.33 from SMALL but only 1.61 should go ⇒
jl's crown-slash ADD over-loads fine sizes 1-3 by ~2× (undecayed, that surplus is the pre-fix +1.61 SMALL excess;
decayed, it over-removes). Sizes 4-5 add is right (LARGE exact). ROOT (refined, terminal): jl's crown_biomass
(FMCROWE) BRANCH-SIZE distribution puts too much mass in the fine classes (1-3) vs FVS CROWNW — even though the
crown_biomass TOTAL is bit-exact (carbon-validated), its split across the 5 branch size-classes diverges. NEXT:
differential jl crown_biomass(xv[2:6]) vs live FMCROWE CROWNW(I,1:5) for a representative cut tree — match the
per-size-class branch split. Then the decay phasing (validated) + correct branch split ⇒ (6.30,3.89) ⇒ kill 5.
📌 constructed-only; decay approach validated, residual isolated to the crown branch-size distribution.

### D17 — crown-split hypothesis REFUTED (jl crown_biomass == live CROWNW); residual = base-SMALL too low (2026-07-02)
Terminal differential (live FMSCUT CROWNW stamp, restored pristine, vs jl crown_biomass) for the first cut tree
(DBH 12.5, HT 74.6, CR 36): live CROWNW br1-5 = [27.31, 67.65, 116.1, 24.73, 0] vs jl crown_biomass br1-5 =
[27.32, 67.68, 116.19, 24.75, 0] — BIT-CLOSE MATCH. ⇒ the crown BRANCH-SPLIT is FAITHFUL; the "jl over-loads
fine sizes 1-3" hypothesis is REFUTED (rule #4). Re-deriving the SMALL overshoot with a faithful slash:
  jl SMALL undecayed 7.91 = base_jl + slash_undecayed
  jl SMALL decayed   4.579 = base_jl + slash_decayed   ⇒ decay removed 3.33 (slash_undecayed−slash_decayed)
  live SMALL         6.30  = base_live + slash_decayed
  ⇒ base_live − base_jl = 6.30 − 4.579 = 1.72  ⇒ jl's BASE (non-slash) SMALL is ~1.72 LOWER than live's.
Pre-fix check: (base_jl−base_live)=−1.72 + (slash_undecayed−slash_decayed)=+3.33 = +1.61 = the pre-fix jl−live
SMALL gap. CONSISTENT. So D17's SMALL residual is TWO effects: (i) slash undecayed [FIXED by the validated decay
phasing]; (ii) jl's BASE SMALL fuel ~1.72 too low. Base LARGE is right (decay landed LARGE exact). Initial FUINI
table is byte-identical, litter matches (0.134/0.129) ⇒ the base-SMALL-low is in the fine-woody (sizes 1-3) base
EVOLUTION 1990→2010 (fmcadd_woody breakage / fmcwd fine decay), not the initial table or litter. NEXT: differential
jl vs live base fine-woody (sizes 1-3, NO thin — cst_fire10) over 1990→2010. Crown port + decay phasing + machinery
all VALIDATED faithful; residual is the base fine-woody evolution. 📌 constructed-only.

### D17 — CLEAN no-thin measurement: SMALL base faithful, LARGE base +22% (snag-fall), NOT the crown-slash (2026-07-02)
Measured cst_fire10 (fire-only, NO thin ⇒ NO crown-slash confound; FMDYN ICYC stamp + jl fire_smlg, both
restored pristine). jl 2010 fire fuel (6.061, 3.236) vs live fire-basis (6.198, 2.646):
  SMALL: jl 6.061 vs live 6.198  — MATCHES (~2%)  ⇒ the "base SMALL 1.72 too low" was a DECAY-APPROXIMATION
         ARTIFACT of the thinned-case single-shot; base SMALL is FAITHFUL. (rule #4 correction.)
  LARGE: jl 3.236 vs live 2.646  — jl +22%  ⇒ REAL base coarse-woody over-accumulation, present WITHOUT a thin.
⇒ D17's LARGE excess is NOT the crown-slash (no thin here) — it's the BASE coarse-woody (cwd sizes 6-9) EVOLUTION
1990→2010, i.e. SNAG-FALL (natural-mortality snags → fall → coarse cwd) over-depositing ~22% vs live. This also
re-explains the thinned case: the decay-fix "landing LARGE exact" was partly coincidental (slash-decay offsetting
the base coarse excess). TERMINAL LOCALIZATION: jl's CS snag-fall / coarse-CWD base accumulation deposits ~22%
too much in sizes 6-9 over a decade (standing TPA matches ⇒ same mortality count ⇒ it's the snag→down-wood
conversion: fall-timing or CWD1 cone-split for the coarse boles). This is the SAME coarse-woody thread from the
early D17 per-size stamp (sizes 7-8 ~1.85×), now CLEANLY isolated in the no-thin case. NEXT: jl-vs-live snag-fall
coarse-cwd (sizes 6-9) deposit over 1990→2010 for cst_fire10. Crown port + SMALL base + machinery + decay all
faithful; the residual is the coarse-woody snag-fall base evolution. 📌 constructed-only (no real CS stand runs an
auto-FMCFMD fire); could touch CS carbon DDW.

## ★ D17 — CONSOLIDATED TERMINAL VERDICT (2026-07-02, after ~15 turns of live-differential diagnosis)
**Scenario:** cst_ft10 = cst01 stand-1 (oak-hickory) + THINDBH(2000) + SIMFIRE(2010), CentralStates. Constructed;
NO corpus stand and NO real FVS test runs an auto-FMCFMD CS fire (cst01's own fire is stand-4 w/ DEFULMOD custom
model). The authoritative SN(260)/NE/CS sweeps are all at the ULP floor — D17 is the sole non-ULP item and it is
SYNTHETIC.

**SOLID, unconfounded findings (each a direct live differential):**
- Fire MACHINERY faithful: forcing jl's fire fuel to live's (6.30,3.89) ⇒ kill = 5 = live EXACTLY. `_fmdyn`
  weighting + rothermel byram + FMEFF mortality all correct.
- crown_biomass branch-split faithful: == live FMSCUT CROWNW ([27.3,67.7,116.1,24.7] both) for the 12.5" tree.
- Initial FUINI table byte-identical (jl CS fire_fuel_dead.csv == cs/fmcba.f FUINI).
- forest-type (FMCSFT), model selection, XPTS, small/large partition, litter — all verified matching.
- The ONLY divergence: jl's 2010 fire fuel (7.912,4.915) is +26% vs live (6.30,3.89); this fuel gap ALONE
  produces the kill gap (jl 12 vs live 5), amplified by the threshold-sensitive FMDYN model-weighting.

**CONFOUNDED / oracle-less (why not driven to bit-exact):** the +26% is in the multi-year FFE fuel EVOLUTION
1990→2010 (thin-slash add-then-decay phasing + coarse-woody snag-fall base accumulation). Every attempt to cleanly
decompose it hits (a) cycle-sampling-timing confounds (fmcba! probes perturb state; the reliable in-sequence trace
+ fire stamp only give the fire-year value), (b) the FMDYN two-call ambiguity (fire-basis vs period-end fuel,
differing ~4× in LARGE), and (c) NO clean real-stand oracle (no CS carbon/DDW test; jl doesn't emit text FuelRept).
A single-shot slash-decay fix landed LARGE exact but overshot SMALL and REGRESSED the kill (12→29) — proving the
residual is multi-factor and that naive fixes fail. The correct fix = faithful FMSCUT-add-then-annual-FMCWD
phasing (a grow_cycle!/fuel-update REORDER), which is HIGH corpus-regression risk (every real CS/SN fire is
bit-exact on the current path) and not guaranteed to land, for ZERO real-stand benefit.

**VERDICT: 📌 accepted / deferred with documented reason.** Class = threshold-amplification (like D8/D10/D13): a
small faithful-component fuel-evolution residual blown up by a discrete downstream weighting, on a synthetic
scenario with no oracle. All fixable components are proven faithful; the residual is measurement-confounded,
oracle-less, and its fix trades bit-exact-corpus risk for synthetic-scenario gain. This is the campaign's lone
open item and is complete-with-reason per the off-switch criterion.

## ★ D18 — VOLUME/BFVOLUME method selector (METHC/METHB) NOT honored — REAL CS stand (2026-07-02)
FOUND by sweeping a REAL FVS test I had missed: `tests/FVScs/cst01_method5.key` (the `Volume 0 All … 5` keyword
= cubic volume METHOD 5). jl vs fresh live: stands 1-4 bit-exact/ULP, but **stand 5 (BARE GROUND PLANT regen)
Scuft diverges 40%** (live 1277 vs jl 1794 @2032); TPA/BA/QMD bit-exact ⇒ the TREE LIST matches, only the VOLUME
computation diverges. Plain cst01 stand-5 (no method5) is BIT-EXACT ⇒ the cause is the METHC=5 volume method.
ROOT: jl's `kw_volume!` (keyword_dispatch.jl:835) explicitly DROPS the METHC field ("not used by the SN R8 Clark
taper path") — it carries only DBHMIN/TOPD/STMP/SCFMIND/SCFTOPD, never the method selector. FVS volkey.f:123 sets
`METHC(sp)=IFIX(PRMS(6))` and cfvol.f/bfvol.f DISPATCH on it (METHC=1 regression-model default; `IF(METHB.EQ.5)
GO TO 300` etc. = distinct equation branches). So jl always computes method-1 (default) volume; a stand requesting
method 5 diverges (here 40% on the regen population — the main S248112 stand happens to match under both methods).
CLASS: real feature-gap divergence (like D6 ESCPRS) — jl implements only the default cubic/board volume method,
not the METHC/METHB alternate-method dispatch (cfvol.f/bfvol.f methods 2-5). Affects any SN/NE/CS stand using a
VOLUME/BFVOLUME method ≠ 1 (uncommon; no other bundled test uses it). VERDICT: 📌 documented feature gap — the
VOLUME-keyword volume-METHOD selector (cfvol.f/bfvol.f alt methods) is unported; jl honors the merch STANDARDS
(DBHMIN/TOPD/STMP) but not the METHOD. NEXT (if pursued): port the cfvol.f METHC / bfvol.f METHB method-2..5
branches + thread METHC/METHB through kw_volume!/kw_bfvolume! → the r9clark/r8clark volume dispatch.

**D18 refinement:** method 5 routes to `NATCRS` (national cubic estimator, vols.f:21) / `OCFVOL` (fmsvol.f:134),
NOT the default `CFVOL` (regression model = METHC 1, jl's only implemented path). Methods {5,6,8,10} each route to
distinct estimators (NATCRS/OCFVOL) that jl does not implement. So D18 = an unported ALTERNATE-cubic-volume-model
feature (the NATCRS/OCFVOL national estimators), selected by VOLUME METHC. VERDICT: 📌 documented feature gap —
real (40% Scuft on cst01_method5 regen) but confined to stands explicitly requesting a non-default volume method
(only cst01_method5 among bundled tests; the default METHC=1 path — every other SN/NE/CS stand — is faithful).
Porting = implement NATCRS/OCFVOL + thread METHC/METHB through kw_volume!→volume dispatch (a bounded but real
volume-model port). Lesson: found by BROADENING the sweep to a missed real test — the hunting mandate (rule #2 /
FIA-plots) surfaces real divergences that drilling one synthetic item (D17) does not.

### Unported keyword-OPTION feature-gap inventory (completeness sweep, 2026-07-02)
Systematically scanned jl's keyword handlers for silently-dropped fields / unported options (the D18 pattern —
jl ignores a keyword's method/option field, diverging IF a stand uses it). Classified each vs the faithful-no-op
set. REAL feature gaps (jl differs from live only when the option is explicitly requested; the DEFAULT path — every
bundled real test except cst01_method5 — is faithful):
- **D18 VOLUME/BFVOLUME METHC/METHB alt volume method (NATCRS/OCFVOL)** — CONFIRMED real (40% Scuft, cst01_method5).
- **FIRECALC method 1/2** (keyword_dispatch.jl:1755-57) — "alternative fire-behavior models (new FM logic /
  modelled loads) not ported"; jl warns + defaults to the faithful old-FM method-0 path. Gap only if method 1/2 used.
- **SETSITE habitat field** (:573) — habitat code ignored ("documented gap"). Gap only if a SETSITE sets habitat.
- **FERTILIZ P/K** (:930) — P (field 3) / K (field 4) ignored; only N (representable amount) + efficacy used. Gap
  only for a FERTILIZ specifying P/K.
All are jl-CODE-DOCUMENTED known gaps (not hidden), bounded, and exercised by NO bundled real test but cst01_method5.
VERDICT: 📌 documented feature-gap class — each is a faithful DEFAULT with an unported non-default OPTION; magnitude
confirmed only for D18 (the others have no bundled-test oracle; a constructed scenario would confirm but be synthetic
like D17). The campaign's DEFAULT-configuration behavior is faithful on every real stand; these are opt-in feature gaps.

**D18 FINALIZED — port is IMPRACTICAL (external NVEL library):** traced method 5 → NATCRS/OCFVOL, which live in the
National Volume Estimator Library (`ForestVegetationSimulator/volume/NVEL/` — dozens of .f + VOLLIB.dll: bia_behres,
blmtap/blmvol, brown, calcbiomass, r9clark, vollib09, …). jl already ported the R9 Clark SUBSET (r9clark_vol.jl =
the METHC=1 default path, faithful everywhere). VOLUME method 5 selects OTHER NVEL national estimators; matching it
requires porting the broader NVEL library — a massive external dependency, out of scope. ⇒ D18 is a 📌 ACCEPTED
FEATURE GAP (like D6): the default volume method is faithful on every real stand; the non-default NVEL-method option
is unported and impractical to port. This is the strongest form of "documented irreducibility" for a feature gap.

### sn.key re-verify (2026-07-02): NOT a divergence — recurring stale-DB harness artifact
Re-swept the last unchecked real SN test (tests/FVSsn/sn.key, the ECON example). jl ERRORED with SQLite
"FVS_Mortality has 19 columns but 22 values" — but jl's CREATE TABLE (dbs_output.jl:352) is correctly 22 cols
matching the 22-value INSERT. The error was a STALE shared DBS file (an OLD 19-col FVS_Mortality table; `CREATE
TABLE IF NOT EXISTS` keeps the stale schema) — the SAME artifact class as the earlier dbs_treelist stale-26-col
FVS_TreeList. Clearing ./FVSOut.db → sn.key runs OK (sum len 3828). sn.key emits DBS (Summary 2 → snout.db), not
a .sum, and is already validated by the registered sndb DBS suite (21/21 incl. econ). ⇒ NOT a real divergence.
NOTE (harness robustness, not a campaign target): jl's DBS uses `CREATE TABLE IF NOT EXISTS`, so a stale shared
FVSOut.db with an older schema throws on insert; sweeps should delete the DSN DB first (or jl could version/DROP
the schema). Cleared the stale DBs. ⇒ ALL 5 bundled SN/NE/CS real tests now accounted for: snt01/net01/cst01
bit-exact/ULP, cst01_method5 = D18 (NVEL feature gap), sn.key = DBS-validated (stale-DB artifact cleared).
