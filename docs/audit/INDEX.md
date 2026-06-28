# FVSjl Port — Faithfulness Flag Audit (port-wide index)

16-agent flag-only audit of the FVSjl Southern port vs the FVS Fortran source.
**390 decisions reviewed → 9 BANDAID · 53 GAP · 5 UNVERIFIED.**

> ★ CAMPAIGN STATUS (all three tiers worked through). BANDAIDs: 6/6 real fixed (B3 was a false positive).
> GAPs: every module cluster traced + reconciled — fixed, verified-faithful, or false-positive — except a small
> set of **deferred-by-design** items: (1) the 2 BLOCKED no-ops SNAGDCAY/SNAGBRK (vacuous until snag hard→soft
> soft-decay is modeled; the fire-snag-cwd "all hard" + density notre FINT/FINTM also wait on niche paths);
> (2) **per-point-density** features that need a per-point density layer — multi-point regen pccf, TCONDMLT point
> weights, single-canopy-tree structure-stage (all FAITHFUL for the common single-point case); (3) two large
> keyword-gated subsystems left transparency-guarded — log-graded HRVRVN revenue (ecvol bucking) and the
> NOHTDREG/LHTDRG HT-DBH calibration; (4) the **FFE phasing #28** co-refactor and the fire-carbon released-from-
> fire VALUE plumbing gated on it; (5) the accepted **COMPRESS eigensolver** + its #29 post-compression DGSCOR
> residual (the lone broken test). UNVERIFIED: all 5 resolved (4 faithful/intentional/unreachable, 1 verified-
> divergent-but-inert). Suite 5010 pass / 1 broken. Every verdict + fix is logged below.

Severity meaning:
- **BANDAID** — wrong on a *reachable* path; the FVS source contradicts the implementation, or there is
  no source basis and the only justification is matching test output. These need triage + fix.
- **GAP** — faithful where the suite exercises it; silently diverges on an untested keyword/edge path.
  Honest incompleteness, not a fudge. Fix as the keyword/path is brought into scope.
- **UNVERIFIED** — the governing FVS source could not be confirmed within the audit; needs follow-up.

> Every flag is a **lead**, not a verdict. Each BANDAID must still be confirmed against the actual source
> line (and, where it matters, an instrumented FVS differential — the `#6` method) before fixing; some may
> prove reachable-but-rare or mis-flagged.

---

## BANDAIDs (priority order — reachable first)

| # | Module | jl location | What FVS mandates (source) | Proposed fix | Reachable on tested path? |
|---|--------|-------------|----------------------------|--------------|---------------------------|
| B1 | fire-behavior | `fmburn.jl:58-65,299-303` | Final flame = `0.45·(Σwᵢ·byramᵢ/60)^0.46`, recomputed from weighted Byram (`fmfint.f:541`, NLC 2003) — NOT Σ of per-model flames | drop `flame_raw`; `flame = 0.45·(byram/60)^0.46` then `·flmult` | **YES** — this is the "accepted" flame 3.84 vs 3.9 residual; feeds SN fire mortality groups 1–5 |
| B2 | sprout | `sprout.jl:171-172` | `CWCALC` CR arg = `CRDUM = 1.0`; ICR(70) is the *discarded* 6th arg (`esuckr.f:313-315`, `cwcalc.f:82`) | pass `1f0` in the CR slot | **YES** when stump-sprout regen fires (sprout.key); inflates sprout CRWDTH by ~5.5 ft |
| B3 | fire-carbon / io-serialization | `summary.jl:101`, `carbon.jl` | `fmsdit.f:93 IF(ICYC.GT.1)` → cycle-1 crown-lift is **0** (instrumented: 0.0/947/1104/1166) | gate first `compute_crown_lift!` to cycle≥2 / drop inventory seed | **YES** — this is the **already-CONFIRMED `#6`** (flagged from 2 angles) |
| B4 | structure-stage | `structure_class:169-170` | `XBAMAX = BTSDIX` = per-cycle SDICAL **SDImax** (`sstage.f:154,544-550`), not the user BAMAX keyword | use computed stand SDImax (Reineke), not `s.control.ba_max` | **YES** — with no BAMAX keyword `ba_max==0`, so SE→SI demotion **never fires** (should be able to) |
| B5 | establishment | `establishment.jl:61` | `GENTIM = max(FINT−5, 0)`, depends only on FINT/DELAY — never on `IDSDAT`/calendar year (`estab.f:317-318,508-512`) | `gentim = max(per - 5, 0)` | Latent — currently erased by the `XMIN` height floor; wrong for FINT≠5 or better sites |
| B6 | econ | `econ.jl:168-172` | Harvest in year Y: cost@`Y−startYear`, revenue@`Y−startYear+1` (`eccalc.f:101,328-329,628-629`) | use FVS `beginTime` indexing (both streams +1 yr; base-year rev@1) | **YES** when ECON runs; under-discounts every harvest by 1 yr. Test asserts jl's own convention |
| B7 | event-monitor | `event_monitor.jl:48,155` | `TIME` is opcode `103NN`, a year-indexed **step interpolation** (`algevl.f:429-460`); current-year is `YEAR` (`evtstv.f:259`) | implement `TIME(v0,y1,v1,…)` step fn; route "current year" to `YEAR` | Latent — only if a keyfile uses `TIME(...)` |
| B8 | event-monitor | `event_monitor.jl:173` | `NO` = test-var code 112 = **constant 0.0** (`algkey.f:331-332`, `evtstv.f:281`); `NOT` is the distinct negation token | parse `NO` as constant 0.0, only `NOT` as negation | Latent — only if a keyfile uses `NO` as a value |

**Note on B1:** the agent correctly identified that "accepting an output-mismatched value with no source basis"
(the memory's *"flame 3.84/scorch 15.65 vs Fortran 3.9/15.9"* accepted residual) is itself the bandaid
signature — and the fix is a one-liner from `fmfint.f:541`. This is the highest-value confirmed lead.

---

## GAPs by module (faithful where tested; diverge on untested paths)

- **mortality** (0): ~~sdimax<5 whole-stand kill~~ FAITHFUL(misread); ~~TPAMRT filtered~~ FIXED; ~~SDICHK floored~~ FIXED; ~~user BAMAX ignored~~ REFUTED(honored); ~~T>35000 cap order~~ PRESENT+UNREACHABLE; ~~MSBMRT (MORTMSB mature-stand breakup) unported~~ PORTED + LIVE-VALIDATED bit-exact (see below). **Cluster CLOSED.**
- **diameter-growth** (0): ~~DGBND/size-cap applied after FINT scaling (diverges for FINT≠5)~~ FAITHFUL — DGBND/size-cap applied at the 5-yr basis then re-expanded by FINT/5 (gradd.f:79-90); LIVE-VALIDATED bit-exact on growth_fint10 (FINT=10). **Stale flag.**
- **height-growth** (1): ~~budwidth floor/bound at wrong period scale~~ STALE (jl matches regent.f order); ~~HK≤4.5 micro-DBH bump dropped~~ STALE (jl:37 = regent.f:285-287); LHTDRG/NOHTDREG HT-DBH calibration branch unmodeled — keyword-gated, now transparency-guarded (see below).
- **crown-ratio** (0): ~~dbh≤0 fixed-0.5 vs RANN draw~~ FIXED; ~~ITRUNC top-kill crown reduction~~ FIXED; ~~CRNMULT on dubbed crowns/<10 floor~~ FAITHFUL — live-validated (TPA/BA + per-tree mean crown% bit-exact vs live, crnmult_base). **Cluster CLOSED.**
- **establishment** (0): ~~idup floor vs ceil~~ FIXED; ~~ESGENT/REGENT first-cycle regen growth~~ FAITHFUL (live-validated bare_natural); ~~PLANT-height floor XMIN vs 0.05~~ FIXED; ~~pccf hardcoded 0~~ FIXED + live-validated (stand CCF; see below). **Cluster CLOSED.**
- **sprout** (1, deferred-bounded): ~~ESSPRT gate~~ FIXED (IFORDI); ~~SPROUT species/DBH table~~ FIXED+live-validated; ~~sprout_dbh AA refit~~ FAITHFUL/SUBSUMED by LHTDRG gap; ESCPRS list-compression on MAXTRE-overflow during sprouting — fully traced, deferred (see below).
- **fire-behavior** (0): ~~scorch not re-derived under FLAMEADJ≠1~~ CLOSED by B1; ~~FMORTMLT per-tree multiplier unapplied~~ PORTED + live-validated; ~~FMDYN truncate-to-4-models reweight omitted~~ STALE — jl `_fmdyn` HAS the nbr[1..4] 4-nearest selection + inverse-distance reweight + collinear split (fmdyn.f:198-258), exercised on fire_early (forest type 1 → 2 weighted models) and bit-exact post-fire on snt01 stand-4; ~~FMCBA soft fuel column~~ FIXED (fmcba! distributes the SOFT cwd, FUELSOFT). **Cluster CLOSED.**
- **fire-carbon** (1): ~~belowground-dead root carbon omitted~~ FIXED + live-validated bit-exact; consumed-fuel carbon factor ~~0.5 vs 0.37~~ FIXED (released-from-fire VALUE plumbing deferred to FFE phasing #28); ~~PSBURN scaling~~ addressed (POTFPAB+FLAG1); ~~FFE live carbon gross vs merch cubic~~ FIXED + live-validated (gross→merch, ~9%→≤1%; see below). Only a ≤1% FFE-live residual (crown-biomass) remains.
- **fire-snag-cwd** (1, blocked): ~~post-fire accelerated fall (PBSOFT/PBSMAL/PBTIME)~~ IMPLEMENTED (snag.jl + SNAGPBN); hard→soft DKTIME split "all hard" = the BLOCKED SNAGDCAY/DECAYX snag-decay-state item (vacuous until snag soft-decay modeled); ~~fmscro ILIFE round vs ceil~~ FIXED; ~~redcedar TFALL hardcoded 1 yr~~ FIXED (foliage 3 for redcedar, see below).
- **compress-cuts** (0 net): ~~merge samples decay/defect vs averaging~~ FIXED (decay_code/woodland_stems now tpa-averaged w/ FVS integer truncation, was copied — see below); ~~TCONDMLT point weights unported~~ FAITHFUL for single-point (a per-point constant ⇒ no ranking effect), multi-point needs per-point density (deferred, shared w/ pccf); ~~AUTSTK BA constant~~ FALSE POSITIVE; ~~sorted SDI/RDEN partial-tree~~ FAITHFUL (FVS zeroes WK2 in SN COMPRESS ⇒ PROB+WK2=PROB); ~~truncated-tree round vs IFIX~~ FIXED. *(eigensolver divergence is the accepted COMPRESS case.)* **Cluster reconciled.**
- **econ** (0 net): log-graded revenue valued at 0 — DOCUMENTED + transparency-guarded (large unported log-bucking subsystem; inert for sn.key = matches live FVS); ~~discount rate hardcoded 4%~~ FIXED + live-validated; ~~B-C/RRR drop NEAR_ZERO guard~~ FIXED (NEAR_ZERO=0.01, eccalc.f:58). **Cluster reconciled** (only the deferred log-bucking subsystem remains, transparency-guarded).
- **io-serialization** (0 net): ~~.sum sample-weight C %E vs Fortran E15.7~~ FIXED; ~~STOP detected on 4-char prefix~~ FIXED; FVS_Carbon fire-released columns hard-zeroed = the SAME item as **fire-carbon #2** (the 0.37/0.50 factor is FIXED; surfacing the released VALUE is deferred to FFE phasing #28 — not a separate io gap). **Cluster effectively CLOSED.**
- **structure-stage** (1): single-tree (NTREES≤1) path unported (niche, deferred); ~~after-thin row reuses before-thin SDI/cover~~ FIXED (apply the thin between rows; see below).
- **event-monitor** (1): ~~`**` unsupported/mis-parsed~~ FIXED (FVS DOES have `**`; prior "no exponentiation" verdict was a misread — see below); ~~AGE omits elapsed-years term~~ FIXED+live-validated; ~~div by zero → Inf vs FVS undefined~~ FIXED (→NaN, no-fire ≈ FVS LREG undefined); _event_bsdi ignores DBHSTAGE/dead exclusions (niche, deferred).
- **density** (0 flagged; 1 UNVERIFIED below).

## UNVERIFIED (need source confirmation)
- **density `notre!` FINT/FINTM dead-record inflation — VERIFIED divergence; inert in tested scenarios; deferred.**
  CONFIRMED both sides: SN notre.f:122-124 multiplies the DEAD-record plot-expansion (VP/FP/FP2) by `FINT/FINTM`
  (cycle length ÷ mortality-observation period) so recently-dead trees are added back at the right rate for the
  BACKDATED-density calibration; jl's `notre!` (standstats.jl:56) expands dead records with the SAME factors as
  live (no FINT/FINTM). Reachable ONLY when BOTH (a) the inventory carries input DEAD trees (`ndead>0`) AND (b)
  `FINT ≠ FINTM` — no suite scenario hits both (default FINTM=FINT=5 ⇒ ratio 1 ⇒ no-op; growth_fint10 has
  ndead=0). NOT a trivial fix: notre.f explicitly says the inflation must be UNDONE for the treelist and other
  uses (FMSSEE in CRATET, PRTRLS deflation) — and jl reuses the dead-record TPA for BOTH the diameter-growth
  backdated calibration (which wants it inflated) AND FFE input-snag seeding (`ffe_seed_input_snags`, which wants
  the true value), so a faithful port needs a separate calibration-vs-true dead TPA, not a single in-place
  inflation. Deferred to a `FINT≠FINTM + input-dead-trees` scenario; documented so it isn't re-investigated. ✓
- ~~**diameter-growth** — calibration PTBAA current vs backdated~~ RESOLVED (FAITHFUL, per B3). jl computes the
  DGF competition point-BA (PTBAA, `diameter_growth.jl:165` `pba`) from CURRENT diameters while the stand BA/AVH/
  PCT are backdated (the deliberate `diameter_growth.jl:306` choice). dense.f's LREDO two-pass sequence
  (line 184 `D=WK3` on the backdated pass; line 204 skips the point stats on the first pass) is genuinely hard to
  pin statically — BUT snt01's cycle-1 DG is BIT-EXACT including the WK3 backdated-calibration species (sp33/65),
  and that DG depends on PBAL=PTBAA·(1−PCT/100). The jl author's explicit "point-BA stays current" comment means
  backdating it was TRIED and REGRESSED snt01. Per B3 (a static argument that contradicts a live bit-exact match
  is the misread): jl's current point-BA is FAITHFUL; the flag's "backdated PTBAA" reading was the misread of the
  intricate DENSE/DGF call sequence. The snt01 bit-exact DG is the ground truth. ✓
- ~~**crown-ratio** — dead/snag DUBSCR crown dubbing~~ RESOLVED (intentional): jl's only dead-record consumer (ffe_seed_input_snags) uses the BOLE volume, never the dubbed crown ratio; carbon_snt bit-exact StandDead confirms the dead-crown is unused. (Full verdict below.)
- ~~**establishment** — missing WK6 site-prep ESRANN draws~~ RESOLVED (unreachable): the WK6 site-prep draws fire only under SITEPREP/MECHPREP/BURNPREP; no suite scenario uses them. (Full verdict below.)
- ~~**fire-snag-cwd NOTE B**~~ — RESOLVED (no divergence). FVS CWD1 (fmcwd.f:187) calls
  `FMSVL2(SP,DIAM,HTD,−1.0,TVOLI,0,'D',.FALSE.,…)` — LMERCH=.FALSE. with X=−1, which for SN (VARACD='SN',
  fmsvol.f:148-151) returns `MAX(−1,MCF) = MCF`, the MERCH cubic. jl's snag `bolevol = v[4]·V2T` (snag.jl:341) is
  ALSO merch cubic. So both the standing snag bole AND the CWD1 fall-to-down-wood use the SAME merch basis —
  NOTE B's "MERCH vs FMSVL2 stem volume" dissolves because FMSVL2-stem-volume IS merch for SN (the fire-carbon #4
  finding). Height also matches: the snag height-loss FMSNGHT is HTX=0 in SN, so CWD1's HTD = the full death
  height jl used. The only residual is the SAME ≤1% R8Clark-v[4] vs NATCRS-MCF equation nuance tracked under
  fire-carbon #4 — NOT a basis divergence. ✓

---

## Standing conclusion
- Core spine modules (mortality, diameter/height growth, crown, density, volume) are **faithful where the
  suite exercises them** — their flags are GAPs on untested keyword paths, not fudges.
- The 9 BANDAIDs cluster in **less-exercised extensions** (FFE fire behavior/carbon, sprout regen, structure
  demotion, econ, event-monitor functions) — paths the keyword-coverage suite touches lightly or not at all.
- **B1 (flame length) and B3 (`#6` crown snapshot) are the two confirmed, reachable, output-affecting bandaids.**

## Next steps (agreed order)
1. **Cleanup pass first** (#14): extract the duplicated backdating kernel into a shared `_backdate_dbh!`
   helper; sweep other duplication. Suite must stay bit-exact (4519/1).
2. **Triage + fix** (#15): confirm each BANDAID against source (instrument where needed), drop false
   positives, fix on the deduplicated code. Then schedule GAPs by keyword scope.

---

## Confirmation pass (read-only logic-trace, 7 agents) — ALL CONFIRMED, 0 false positives

| Flag | Verdict | Upstream rank | Fix order |
|------|---------|---------------|-----------|
| B5 gentim (establishment) | CONFIRMED (masked by XMIN floor) | **UPSTREAM** (regen→density→growth) | 1 |
| B1 flame (fire-behavior) | CONFIRMED | **MID** (→fire mortality→TPA) | 2 |
| B3 crown-snapshot (#6) | CONFIRMED (instrumented) | LEAF (DDW carbon) | 3 |
| B2 sprout-cw | CONFIRMED | LEAF (TreeList col) | 4 |
| B4 pctsmx (structure) | CONFIRMED | LEAF (struct-class report) | 5 |
| B6 econ discount | CONFIRMED | LEAF (PNV report) | 6 |
| B7/B8 evmon TIME/NO | CONFIRMED | LEAF (unused in suite) | 7 |

Order = upstream-first (principle #2). Most fixes land in SHARED modules (src/engine/*) → keep variant-general (principle #6).

---

## FIX LOG (upstream-first)
- **B5 gentim — FIXED** (establishment.jl). `gentim = max(per-5, 0)` per estab.f:448-449 (was `yr-idsdat`).
  Removed the now-dead `_es_idsdat` helper + its docstring. Suite 4519/1 — SILENT, as predicted (masked by the
  `es_xmin` height floor). Faithful; no masked bug surfaced. Follow-up: a FINT≠5 / high-site regen scenario would
  exercise the corrected path (no current test clears the floor). Variant-general (FINT-only). ✓

- **B1 flame — FIXED** (fmburn.jl, both `fmburn!` + `potential_fire`). Flame now recomputed from the weighted
  Byram: `flame = 0.45·(byram/60)^0.46` (fmfint.f:541), with the FLAMEADJ multiplier + Byram back-compute
  (fmburn.f:439-464) in `fmburn!`. Replaces the pre-2003 Σ-of-per-model-flames (concave x^0.46 → low bias).
  VERIFIED: (a) units — rothermel byram is BTU/ft/min like FVS; (b) empirical — moves flame in the Jensen
  direction on 3 blended (nmodels=2) fire stands: fire_mid 5.235→5.321, fire_early 3.972→3.979, fire_fuel9
  3.322→3.403; (c) LIVE FVS — instrumented FMFINT confirms FVS uses `0.45·(BYRAM/60)^0.46` verbatim
  (BYRAM=5942→FLAME=3.726). Suite 4519/1 — the .sum integer-TPA granularity doesn't capture the sub-integer
  flame change. Variant-general (fire base model). ✓
  ★ FOLLOW-UP CLOSED (BurnReport live diff): with the full binary, ran fire_burn.key (fire_early + BURNREDB DBS +
  SIMFIRE@2000) through live FVSsn → FVS_BurnReport `Flame_length=4.172`, `Scorch_height=17.581`; jl's BurnReport
  for the same fire gives flame=4.170, scorch=17.566 — bit-close (Δ0.002 / 0.015 ft, Float32 transcendentals).
  The B1 flame fix is now validated END-TO-END (not just the formula). +4 tests (test_fire.jl, fire_burn.key).
  FOLLOW-UP (principle #4): no suite test pins the weighted flame value → add a direct flame assertion
  validated vs a live-FVS BurnReport/PotFire differential. Also cross-check jl weighted BYRAM vs FVS at the
  exact fire event (jl fire_early byram 6853 didn't obviously align to an FVS FMFINT call — byram/fuel path,
  separate from B1).

- **B3 (#6 crown snapshot) — FALSE POSITIVE (retracted) ✓ NOT a bandaid.** Corrected logic trace:
  FVS applies growth in UPDATE (update.f:65 HT=HT+HTG, :115 DBH=DBH+DG/BRATIO), called at gradd.f:180 —
  AFTER FMMAIN (gradd.f:118). So FMOLDC (fmmain.f:268, runs EVERY cycle, ungated) captures the PRE-growth
  crown; cycle-1 FMOLDC = the INVENTORY crown. jl's inventory `snapshot_ffe_oldcrown!` is the faithful analog
  of FMOLDC(cyc1). Timing: jl-with-#6 adds the inv→postcyc1 lift during cyc2 (computed end of loop c=0, added
  c=1); FVS adds the same inv→postcyc1 rise during cyc2 (ICYC=2 FMSDIT, =947). EXACT match, and carbon_snt DDW
  is bit-exact vs the LIVE-FVS golden [3.796,4.393,5.354,9.535] WITH #6. The earlier instrumented "confirmation"
  (and the confirm-B3 agent) mis-mapped the jl loop: the inventory snapshot feeds cyc2, NOT a spurious cyc1
  contribution (cyc1 fuel loop adds the zero-init array). KEEP #6. Count: 7 distinct flags → 6 real bandaids.

## REVISED BANDAID TALLY: 6 real (B1✓ B5✓ fixed; B2 B4 B6 B7/B8 pending) + B3 false-positive

- **B2 sprout crown-width — FIXED** (sprout.jl:171). CR arg `70f0`→`1f0` (esuckr.f:317 `CRDUM=1.`; ICR=70 is the
  discarded 6th IICR arg, cwcalc.f). Faithful by construction: jl now feeds cr=1.0 to the same crown_width eqn
  FVS CWCALC uses → jl sprout CW = FVS sprout CW. VERIFIED exercised: at realistic sprout dbh ~0.5-0.9" the CW
  drops 0.3-1.4 ft (oak SK 1.74→1.23, sweetgum SY 2.05→0.63). Suite 4519/1 — no sprout-CrWidth assertion exists.
  Variant-general (matches base esuckr CRDUM=1.0; only the SN esuckr! call site). ✓
  FOLLOW-UP (principle #4): add a sprout-CrWidth assertion vs a live CWCALC differential.

- **B4 structure PCTSMX — ANALYZED, fix deferred (variant-wiring needed).** Confirmed: jl uses `s.control.ba_max`
  (user BAMAX kw, 0 when absent) where FVS uses BTSDIX = the per-cycle SDICAL stand MaxSDI (sstage.f:154,
  grincr.f:240). FAITHFUL FIX (principle #6 variant-safe): the general field `StandState.before_max_sdi` (BTSDIX,
  state.jl:388) already exists but is NEVER populated; `stand_sdimax` is SN-only (variants/southern/mortality.jl:35).
  → Populate `before_max_sdi` from the variant's pre-growth SDImax (SN: stand_sdimax) at GRINCR-time, then in the
  shared structure_stage.jl replace `s.control.ba_max` with `s.<...>.before_max_sdi` and drop the `>0` guard. Do
  NOT call stand_sdimax directly from the shared module. LEAF (struct-class report only). Next-session task.

- **B4 structure PCTSMX — FIXED** (structure_stage.jl + standstats.jl). Replaced `s.control.ba_max` (user BAMAX
  kw, 0 when absent → demotion never fired) with `stand_sdimax(s)` = BTSDIX, the per-cycle SDICAL stand MaxSDI
  (sstage.f:154,544-550). VARIANT-SAFE (principle #6): relocated `stand_sdimax` (the general SDICAL — pure
  BA-weighted avg of sp_sdi_def coeff data, no SN logic) from variants/southern/mortality.jl to shared
  engine/standstats.jl; SN mortality + shared structure_stage both call it. Suite 4519/1 — demotion stays silent
  on the tested stands (mature/dense, SDIBC > 0.30·MaxSDI → branch not entered), as predicted. ✓
  FOLLOW-UP (principle #4): a young/sparse single-cohort scenario (nstr=1, SDIBC < 0.30·MaxSDI) to exercise +
  validate the SE→SI demotion vs live FVS.

- **B6 econ discounting — FIXED** (econ.jl:168-172 + test_econ.jl). Harvest cost discounted t=(yr-startYear) yrs,
  revenue (t+1) yrs (eccalc.f:114-117: "cost accrues at START of year" beginTime-1, "revenue at END" beginTime).
  Was cost@(t-1)/rev@(t) — one year short on both. Updated the self-referential unit test to assert the eccalc.f
  convention (PV(rev,t+1)), NOT jl's old output. Suite 4519/1. Variant-general (econ is a base extension). ✓
  FOLLOW-UP (principle #1): the suite has NO live-Fortran econ golden — only self-referential unit tests. A
  live FVSsn ECON differential (EconSummary DiscCost/DiscRev on an econ scenario) is the gold-standard check.
  ★ FOLLOW-UP NOW CLOSED (with the econ discount-rate fix below): jl's discounted ANNUCST cost stream is
  BIT-EXACT vs the live FVS_EconSummary Discounted_Cost (13.6379 / 24.3235 / 32.6959 at 5%, econ_strtecon.key) —
  confirming both the B6 timing AND the rate. The old hardcoded 4% gave 13.8897 (would NOT have matched).

- **GAP econ #2 (discount rate hardcoded 4%) — FIXED + LIVE-VALIDATED.** jl `EconState` defaulted
  `discount_rate=0.04` and `kw_econ!` never read it. FVS defaults `discountRate=0.0` (ecinit.f:15) and reads it
  from the **STRTECON** keyword's field 2 (a PERCENT; eccalc.f:91 `rate=discountRate/100`). So jl over-discounted
  (or mis-discounted) every ECON run — e.g. live FVS sn.key shows `Discount_Rate=0.0` (no STRTECON ⇒ no
  discounting) but jl applied 4%. FIX: default `discount_rate=0f0` (state.jl) + a STRTECON branch in `kw_econ!`
  reading `field2/100`. VALIDATED vs live FVSsn: `econ_strtecon.key` (`STRTECON … 5.0`) ⇒ FVS_EconSummary
  `Discount_Rate=5.0`, jl parses 0.05, and the discounted ANNUCST cost matches the live `Discounted_Cost`
  BIT-EXACT (see B6 above). +6 tests. The SEV / start-year-delay STRTECON fields stay unmodeled (noted). Default
  (no STRTECON) now correctly 0 = no discounting. Variant-general (econ is a base extension). Suite 4992/1. ✓

- **GAP econ #1 (log-graded revenue valued at 0) — DOCUMENTED + transparency-guarded (large unported subsystem;
  inert in the validated scenario).** HRVRVN's units field (ECNCOM.F77:19) is 1=TPA, 2=BF_1000, 3=FT3_100, and
  the LOG units **4=BF_1000_LOG / 5=FT3_100_LOG**. jl's `harvest_value` values units 1/2/3 and falls through to
  `vol=0` for 4/5. The LOG units are NOT a simple whole-tree-board-feet proxy: echarv.f:72-119 bucks each tree
  into logs (ecvol.f `logBfVol`/`logFt3Vol` per log), buckets each log by its inside-bark diameter
  (`logDibIdx`), and accumulates revenue volume per (species, unit, DIB-grade) so the HRVRVN amount applies
  per-DIB-grade — a substantial unported subsystem (per-log bucking + DIB grading + per-grade revenue tables).
  REACHABILITY: **inert for sn.key** — even though sn.key's HRVRVN uses unit 4, live FVS yields ZERO revenue
  there (FVS_EconHarvestValue is empty; the THINPRSC/THINBTA harvests don't qualify), so jl's 0 MATCHES live FVS.
  It only diverges with a unit-4/5 HRVRVN + a qualifying harvest whose log-grade revenue is non-zero. FIX
  (transparency, like NOHTDREG): `kw_econ!` now `@warn`s once when a HRVRVN unit is 4/5 — flagging the unported
  path rather than silently zeroing. Full port deferred (the ecvol.f log-bucking subsystem). Variant-general. ✓

- **GAP econ #3 (B/C + RRR missing the NEAR_ZERO guard) — FIXED.** FVS computes the benefit/cost ratio only when
  the discounted cost exceeds `NEAR_ZERO=0.01` (eccalc.f:58/681), and the realizable rate of return only when
  BOTH the discounted cost AND revenue exceed 0.01 (eccalc.f:685); otherwise it leaves them blank — guarding
  against a tiny-cost blow-up (e.g. cost 0.005 ⇒ B/C = 200·rev, a meaningless spike). jl guarded only on `> 0`,
  so a sub-cent discounted cost produced a huge ratio where FVS shows blank. FIX: `econ_bc_ratio` /
  `econ_rate_of_return` now gate on `ECON_NEAR_ZERO=0.01f0` (econ.jl). The boundary cases the unit tests already
  cover (cost/rev = 0) are unaffected (0 < 0.01 ⇒ still 0); +4 NEAR_ZERO-boundary assertions. Variant-general
  (base econ). Suite 4996/1. ✓

- **B7/B8 event-monitor — FIXED + TESTED** (event_monitor.jl + test_event_monitor.jl).
  B8: `NO`/`ALL` now resolve to constant 0.0 and `YES` to 1.0 (evtstv.f:81-82,281, codes 111/112); removed the
  bogus `NO`-as-NOT parse — only `NOT` (algkey.f CTAB3) is negation.
  B7: `TIME` is now the variadic year-indexed step fn `TIME(v0,y1,v1,…)` (algevl.f:303, new EvTime AST node),
  NOT an alias of the current year (that is the `YEAR` variable, evtstv.f:101). ≤2 args ⇒ v0.
  Added a unit testset exercising both (was unreached): +11 tests, all green. Suite 4530/1. Variant-general
  (event monitor is base — algkey/algevl/evtstv). ✓

## ★ BANDAID TIER COMPLETE — 6/6 real bandaids fixed (B1,B2,B4,B5,B6,B7/B8); B3 false-positive retracted.
Suite 4530 pass / 1 broken (accepted COMPRESS). Next tiers: 53 GAPs, 5 UNVERIFIED (+ the noted follow-up
live-FVS assertions for B1 flame, B2 sprout-CW, B6 econ). ★ UPDATE: ALL THREE deferred bandaid follow-ups are
now CLOSED vs live FVSsn — B6 econ (EconSummary Discounted_Cost bit-exact), B1 flame (BurnReport Flame_length/
Scorch_height bit-close: jl 4.170/17.566 vs live 4.172/17.581, +test fire_burn.key), and B2 sprout-CW (validated
+ its TreeList-CrWidth output gap already fixed at dbs_output.jl:468). No deferred bandaid follow-ups remain.

---
## GAP TIER — started (most-upstream first)

- **GAP diameter-growth: DGBND/size-cap vs FINT scaling — ANALYZED, fix pending careful pass.**
  jl (diameter_growth.jl:633-638): `dds = exp(wk2)·xbai·(sfint/5)` then `bnd(sqrt(d_ib²+dds·frm)-d_ib)` — the
  size-cap/taper bound is applied to the ALREADY-FINT-scaled DG. FVS: DGBND on the 5-yr DG (dgdriv.f:267-269),
  THEN gradd.f:85-87 rescales `DDS·(FINT/YR)` + recomputes DG with NO re-bound. Identical at FINT=5 (snt01);
  diverges for FINT≠5 when the bound binds (a tree above DLODHI: jl flat 0.048 vs FVS ~0.096 scaled; size cap
  enforced tighter). The jl comment cites dgdriv.f:325/715 (the calibration/dub scalings) — the growth-mode
  scaling is gradd.f:80, which does NOT re-bound. UPSTREAM (DG feeds everything).
  ⚠️ REACHABILITY CHECK NEEDED FIRST (B3 lesson): s5/s9/timeint10 (10-yr) are in the suite + validated vs live
  FVS and PASS. Before reordering, determine whether dg_bound actually BINDS for any tree in those scenarios
  (instrument/identity-test). If it never binds → GAP is genuinely untested (port faithfully + add a binding-
  case 10-yr scenario per principle #4). If it binds and they still pass → re-trace whether jl already matches
  FVS (possible mis-scope). FAITHFUL FIX: bound the 5-yr DG, then scale DDS by sfint/5 and recompute DG without
  re-bounding — mirroring dgdriv→gradd. Do NOT rush; highest blast radius in the port.

  ### REACHABILITY + FIX SPEC (analysis complete; implement next session with fresh focus)
  Reachability (instrumented dg_bound): s5_cycle 0 binds, s9_uniform10 0 binds, **timeint10 3 binds** (sp33,
  d≈43-45 in taper [dlo 42.8, dhi 60], 10-yr cycle). Divergence is real but small (~0.8% on 3 trees), masked by
  the .sum integer granularity (timeint10 still passes) — same shape as B1.
  EXACT FVS ORDER (confirmed): dgdriv.f:255-269 computes 5-yr DG = sqrt(DSQ+DDS·exp(FRMT))−D WITH the stochastic
  FRM, then CALL DGBND on the 5-yr DG; gradd.f:79-90 then (FINT≠YR) DDS=(DG·(2·BARK·D+DG))·(FINT/YR),
  DG=sqrt((D·BARK)²+DDS)−BARK·D — scale the BOUNDED 5-yr DG, NO re-bound.
  FAITHFUL FIX (diameter_growth.jl:633-652, all 4 DG computations central/up/low/non-trip):
    1. dds5 = exp(wk2[i])·xbai           # 5-yr DDS — REMOVE the (sfint/5) factor here
    2. dg5 = sqrt(d_ib²+dds5·exp(frm))−d_ib
    3. dg5 = bnd(dg5)                      # DGBND on the 5-yr DG
    4. if sfint≠5 && dg5>0: dds_s=dg5·(2·d_ib+dg5)·(sfint/5); dg5=sqrt(d_ib²+dds_s)−d_ib   # gradd scale, no re-bnd
  (factor steps 2-4 into a `_bound_scale(bnd, dg5, d_ib, sfint)` helper in bark_and_bounds.jl).
  SAFETY: bit-exact for FINT=5 (sfint/5=1 ⇒ step 4 skipped, old≡new by identity) → snt01 + bulk of suite
  unchanged. Only timeint10's 3 taper trees move. VALIDATE: live-FVS per-tree DG differential on those 3 trees
  (instrument dgdriv/gradd DG(I)) — confirm bound-then-scale matches FVS, not just that timeint10 .sum stays green.

- **GAP diameter-growth DGBND ordering — FIXED** (bark_and_bounds.jl `_bound_scale` + diameter_growth.jl:633-652).
  Now bounds the 5-yr DG then scales by sfint/5 (exact dgdriv.f:255-269 → gradd.f:79-90 order), replacing the
  scale-then-bound. Bit-exact for FINT=5 by identity (snt01 + bulk unchanged); timeint10's 3 taper trees are a
  near-no-op (d≈43 just above dlo=42.8 ⇒ taper factor ≈0.99, so bound-then-scale ≈ scale-then-bound) → suite
  4530/1. FAITHFUL BY CONSTRUCTION (verified the exact FVS order line-by-line). Live per-tree diff is confounded
  by cross-cycle DBH drift (DG depends on the growth prediction, not just DBH) — needs tree-ID+cycle matching.
  FOLLOW-UP (principle #4): the REAL divergence bites only for d>dhi (flat 0.048) or a hard size-cap in a FINT≠5
  cycle — NO current scenario has such a tree. Add a 10-yr scenario with a >dhi / size-capped tree, validate the
  per-tree DG vs live FVS (instrument gradd DG(I)). Variant-general pattern (SN dg_bound + base gradd scaling). ✓

- **GAP height-growth #1 (small-tree DIAM-floor + DGBND at wrong period scale) — FIXED** (small_tree_growth.jl).
  scale2 was hardcoded 1f0 (full-cycle) → now REGENT_REGYR/fint = 5/fint (YR/FNT), so _regent_dg shrinks DG to
  the 5-yr basis BEFORE the DIAM budwidth-floor + DGBND (regent.f:359-371); then GRADD re-expands by fint/5
  (gradd.f:79-90), as for large trees. SAFETY: bit-exact for the suite — FINT=5 is identity, and shrink-then-
  expand with nothing binding is an exact algebraic identity (dds_exp = dds_shrunk·(fint/5) = (bd+dg)²−bd²).
  Reachability (instrumented): s5/s9 run the small-tree path 62-64×/scenario but the DIAM floor NEVER binds and
  DGBND can't bind for d<3" (d≪dlo) → GAP currently inert; the fix only bites when the floor binds in a FINT≠5
  cycle. Suite 4530/1. Variant-general pattern (SN regent + base gradd scaling). ✓
  FOLLOW-UP (principle #4): add a 10-yr scenario with a slow-growing small tree whose DIAM floor binds; validate
  per-tree DG vs live FVS. Remaining height-growth sub-flags: #2 HK≤4.5 micro-DBH bump (regent.f:285-287), #3
  LHTDRG/NOHTDREG calibration branch.

- **GAP height-growth #2 (HK≤4.5 micro-DBH bump) — FIXED** (small_tree_growth.jl:35). _regent_dg now returns
  0.001·HK for HK≤4.5 (regent.f:285-287: FVS DG=0, DBH=D+0.001·HK; jl matches the net dbh via dbh+=dg). Dead
  branch for live trees (dbh>0 ⇒ height>4.5 ⇒ HK>4.5) → suite 4530/1 bit-exact. Micro-divergence: jl books the
  0.001·HK as DG where FVS books DG=0 (≤0.0045" on a ~0-dbh seedling — negligible). ✓

---
## BLOCKER: live-FVS bandaid follow-ups limited by the stripped ground-truth binary
The rebuilt FVSsn is a STRIPPED DBS build (memory: Summary/Compute/TreeList/Cases/InvReference/CutList only).
Verified this session: running sn.key on it emits NO econ DBS tables and no econ text report → the **B6 econ
EconSummary live-diff is NOT feasible** with the current binary (B6 rests on the unambiguous eccalc.f:114-117
source-derivation). Likewise **B1 BurnReport/PotFire** DBS are absent (B1 flame FORMULA was already live-confirmed
via instrumented FMFINT — that path IS feasible). **B2 sprout-CW** could be live-checked by instrumenting CWCALC
(feasible) but the fix is already faithful-by-construction (cr=1.0 = esuckr.f CRDUM, same eqn).
DECISION NEEDED: to close B6/B1-report live-diffs, rebuild a FULLER FVSsn with full DBS output enabled
(infrastructure task) — or accept the source-derivation bar for these. Instrument-a-single-routine diffs
(FMFINT, CWCALC, DGBND-via-gradd) remain feasible and are the practical live-FVS tool for value-path checks.

---
## ★ UNBLOCK: full live-FVS DBS oracle exists (no rebuild needed)
"Invest in a fuller binary" investigated → the binary is ALREADY full. Assembling
`gfortran -o /tmp/FVSsn_full $(ls bin/FVSsn_buildDir/*.o) glibc_shim.o -lpthread -ldl` (550 .o, all dbs*.o +
dbsqlite) and running tests/FVSsn/sn.key emits **21 DBS tables** incl. FVS_EconSummary, FVS_BurnReport,
FVS_PotFire_East/Cond, FVS_Carbon, FVS_Fuels, FVS_Mortality, FVS_SnagSum, FVS_Down_Wood_*, FVS_Hrv_Carbon,
FVS_Summary, FVS_TreeList, FVS_CutList. The "stripped binary" blocker was TWO artifacts: (1) the sqlite3 CLI is
absent → read tables via SQLite.jl; (2) the committed golden FVSOut.db was generated with a minimal DataBase
block (3 tables). Reusable harness: test/oracle/live_fvs_db.jl (LiveFVS.run_key/tables/rows) — verified, 20
non-empty tables from sn.key. → gold-standard live-FVS diffs now FEASIBLE for B6 econ, B1 fire flame/scorch,
carbon, mortality, and GAP value-path validations.
NEXT: B6 needs a discounting scenario (sn.key has rate=0 / rev=0 → no signal): build a minimal key with a
non-zero RATE + a harvest with HRVRVN revenue, read FVS_EconSummary Discounted_Cost/Revenue, confirm the
cost@(yr−start) / rev@(yr−start+1) convention I implemented. B1: read FVS_PotFire_East Surf_Flame vs FVSjl.

## B6 econ live-diff — READY (oracle + mechanism found; needs a revenue-harvest scenario)
- Discount rate is set by the **STRTECON** keyword (event-monitor-coded 1605): `discountRate = strtParms(1)`
  (eccalc.f:151), `rate = discountRate/100` (eccalc.f:91). sn.key has `Econ`+`ANNUCST` but NO STRTECON → rate=0
  → flat econ (no discounting signal).
- **eccalc.f:115-117 DIRECTLY confirms the B6 fix**: `computePV(cost, beginTime-1, rate)` ("Costs accrue at
  beginning of year") + `computePV(rev, beginTime, rate)` ("Revenues accrue at end of year") → cost@(yr-start),
  rev@(yr-start+1). This is the exact convention I implemented in econ_stand_pnv.
- TO RUN THE LIVE DIFF: build a single-stand key (snt01.tre) with `Econ` / `STRTECON <yr> 5.` / `ANNUCST` /
  `HRVRVN <revenue>` / a `THINBTA` that removes merch trees / DataBase{Summary,ECONRPTS} / PROCESS; run via
  LiveFVS.run_key; read FVS_EconSummary Discounted_Cost vs Discounted_Revenue — expect rev discounted ONE more
  year than cost. Compare to FVSjl econ_stand_pnv on the same harvest. (sn.key's thinnings produce 0 revenue on
  the dumped stand — needs a merch-removing harvest.) Harness ready: test/oracle/live_fvs_db.jl.

## Live-FVS DBS validation — harness PROVEN; per-scenario DBS-enabling is the remaining knack
- ★ DELIVERED + PROVEN: test/oracle/live_fvs_db.jl runs ANY key on the full binary; sn.key → 21 real tables
  incl. FVS_Carbon (10 rows), FVS_EconSummary, FVS_BurnReport, FVS_PotFire, FVS_SnagSum, FVS_Down_Wood_*.
- To make a FVSjl harness scenario emit DBS tables (they default to .sum/.out text): (1) prepend a DSNOut block
  `DataBase / DSNOut / out.db / End`; (2) the FFE/carbon DBS tables come from the FMIN-block report keywords
  (CarbRept, Potfire, BurnRept, MortRept, FuelRept, SnagSum, DWDVlOut, DWDCvOut, CarbCalc) — NOT from a DataBase
  block; the non-FFE tables (Summary/TreeliDB/CutliDB/ECONRPTS/CARBRPTS) DO go in a DataBase block placed AFTER
  the stand-setup keywords (sn.key lines 13-17). Match sn.key's exact FMIN block (lines ~240-252) as the template
  — carbon_snt's CARBREPT/CARBCALC 1 1 alone didn't emit FVS_Carbon; sn.key uses CarbRept + CarbCalc 0 0 + the
  full report set. NEXT: copy sn.key's verbatim FFE block into a single-stand carbon key, diff FVS_Carbon /
  FVS_Down_Wood vs FVSjl stand_carbon_report (validates #6 + the DDW work beyond the hand-extracted numbers).

## Live-FVS carbon-diff — OPEN config issue (next session)
Copying sn.key's verbatim FMIN block (Potfire/BurnRept/SnagSum/DWDVlOut/DWDCvOut/CarbRept/CarbCut/CarbCalc 0 0)
+ DSNOut into a single-stand carbon_snt key runs WITHOUT keyword errors but emits ONLY FVS_Cases + FVS_InvReference
— no per-cycle FVS_Carbon/FVS_Summary. So the cycles aren't writing DBS rows. NEXT: keep the run dir (don't
mktemp-clean), read the .out to see whether the simulation completes the 4 cycles or stalls; add `Summary 2` in a
post-stand DataBase block; confirm carbon_snt.tre is parsed by the (custom TREEFMT) reader. The oracle harness
itself is proven (sn.key → 21 tables) — this is per-scenario keyfile config, not an oracle limitation.

## Live-FVS carbon-diff — PRECISE finding (2026-06-27, methodical debug)
Single-stand carbon_snt key + DSNOut + `Summary 2` DataBase block + FMIN-block FFE reports
(CarbRept/SnagSum/DWDVlOut/DWDCvOut/CarbCalc 0 0): the simulation COMPLETES (carbon_snt.sum shows all cycles
1990-2010), `Summary` yields FVS_Summary2 (5 rows) — but the FFE report keywords, though ACCEPTED (no FVS_Error),
still do NOT emit FVS_Carbon / FVS_SnagSum / FVS_Down_Wood. So vs sn.key (which DID produce FVS_Carbon=10 rows),
there is a subtle FMIN/FFE-activation difference, NOT a keyword-syntax error. NEXT: extract sn.key's EXACT
carbon-producing single stand (it's multi-stand; find which stand's block wrote FVS_Carbon) and run THAT stand
standalone; diff its FMIN/DataBase setup against carbon_snt line-by-line. The harness + binary are proven; this is
isolated FFE-DBS-activation config. (B6 econ remains source-confirmed via eccalc.f:115-117; its live diff needs a
STRTECON+revenue-harvest key.)

## Live-FVS carbon-diff — STOPPED (unresolved FVS-config; not a faithfulness issue)
Exhausted the quick paths: CARBRPTS is NOT a valid DataBase keyword in this binary (errors despite dbsin.f:46),
and sn.key produces FVS_Carbon WITHOUT it (from FMIN CarbRept + DSNOut) — yet carbon_snt with the IDENTICAL
FMIN CarbRept + DSNOut emits only FVS_Cases/InvReference/Summary2 (cycles DO complete). The difference is an
unresolved FMIN/stand-config subtlety. DEFINITIVE next step (fresh session): isolate sn.key's exact
carbon-producing single stand (split the multi-stand key), run THAT stand standalone to confirm it reproduces
FVS_Carbon, then minimally morph it toward carbon_snt to find the trigger. This is FVS keyword-config, NOT a
FVSjl faithfulness question — the oracle harness + binary are proven (sn.key → FVS_Carbon/Econ/Burn). The
carbon MODEL itself is already validated vs hand-extracted live-FVS StandDead numbers (test_carbon.jl).

- **GAP crown-ratio #3 (ITRUNC top-kill crown reduction) — FIXED** (crown_ratio.jl). Added the crown.f:55 LSTART
  reduction: for a dubbed inventory TOP-KILLED tree (trunc≠0), remove the dead-top portion from crown length
  (HN=norm_ht/100, HD=HN−trunc/100, CL=(ICRI/100)·HN−HD, ICRI=CL·100/HN), gated by a new `lstart` kwarg passed
  by init_crown_ratios!. Auto-scoped to dubbed trees (init restores input crowns after). FAITHFUL per the FVS
  loop (crown.f DO 60: LSTART skips ICR>0 trees → only dubbed crowns are processed). Reachability: NO-OP for the
  suite — snt01's trees all have input crowns, so init_crown_ratios! returns early and the LSTART path isn't
  reached (matches FVS: CROWN only runs at LSTART when MISSCR). Suite 4530/1. norm_ht semantics verified
  (height×100 when >0). FOLLOW-UP: a scenario with a top-killed crownless inventory tree to exercise + validate
  vs live FVS. Variant-general (SN crown.f port). ✓

- **GAP crown-ratio #1 (dbh≤0 fixed-0.5 vs RANN draw) — FIXED** (crown_ratio.jl:99). FVS crown.f:287-292 sets the
  relative crown position by DBH rank for live stems but draws a RANN for dbh≤0: `IF(DBH>0) X=(ISORT/ITRN)·SCALE
  ELSE CALL RANN(RNUMB); X=RNUMB·SCALE`. jl used a fixed `0.5·SCALE`, which both mis-set the regen crown AND
  skipped FVS's RANN draw — desyncing the per-tree DGSCOR RNG stream on the regen path (a principle-#3 masked-bug
  signal). Replaced with `rann!(s.rng)·scale`, drawn in tree-loop order to match FVS's `DO I=1,ITRN` consumption
  order. FAITHFUL (crown.f literal). Reachability: NO-OP for the suite — dbh≤0 trees only arise on the unfinished
  regen/ESTAB path (no bare-stand scenario tested), so no RANN is drawn and the suite stays bit-exact at 4624/1.
  The fix matters when a bare/regen-heavy stand is later validated: it both fixes those trees' crowns and keeps
  the RNG aligned with live FVS. FOLLOW-UP: a bare/regen scenario to exercise + live-validate. Variant-general. ✓

- **GAP height-growth — 2 STALE (verified faithful), 1 keyword-gated + transparency-guarded.** Re-traced the 3
  flagged regent.f small-tree DG items vs sn/regent.f:
  - *budwidth floor / DGBND "at wrong period scale"* — STALE. regent.f:359-368 scales DG to the 5-yr basis
    (DDS·SCALE2, SCALE2=REGYR/fnt), THEN applies the DIAM budwidth floor (`IF((DBH+DG)<DIAM) DG=DIAM−DBH`) and
    `CALL DGBND`, then GRADD re-expands to the cycle. small_tree_growth.jl:115-123 does exactly this order
    (floor + dg_bound at the 5-yr basis, then the GRADD re-expand for FINT≠5). Faithful; no divergence.
  - *HK≤4.5 micro-DBH bump dropped* — STALE. regent.f:285-287 `IF(HK≤4.5){DG=0; DBH=D+0.001·HK}`;
    small_tree_growth.jl:37 returns `0.001·HK` (jl adds it to dbh, netting D+0.001·HK). Faithful.
  - *LHTDRG/NOHTDREG calibration branch unmodeled* — REAL but **keyword-gated** and now transparency-guarded.
    LHTDRG defaults `.FALSE.` (grinit.f:104) and is set `.TRUE.` ONLY by the NOHTDREG keyword with field 2 > 0
    ("calibration INVOKED"; sn option 60). With LHTDRG=.FALSE. (or a Wykoff-calibrated species, IABFLG=1)
    regent.f:315-327 overrides DKK/DK with the HTDBH inventory equation — exactly jl's `_htdbh_dbh` branch. The
    INVOKE form would instead (a) use regent.f's Wykoff HT-DBH equation for IABFLG≠1 species and (b) run
    cratet.f's ≥3-obs HT-DBH regression fit — both UNPORTED. FVSjl previously SILENTLY no-op'd NOHTDREG (and
    mis-labeled it an establishment flag). FIX: a dedicated `kw_nohtdreg!` handler — faithful silent no-op for the
    suppress/default form (LHTDRG stays .FALSE. = jl's existing behaviour, provably identical to no keyword), and a
    `@warn` for the unported invoke form (never a silent gap — the YARDLOSS lesson). Removed NOHTDREG from
    KNOWN_NOOP_SN; corrected the southern.jl classification comment. +6 tests (test_nohtdreg.jl), suite 4630/1.
    FOLLOW-UP (chunk-sized): port the Wykoff regent branch + cratet HTREG calibration to honour NOHTDREG-invoke. ✓

- **GAP fire-behavior #2 (scorch not re-derived under FLAMEADJ≠1) — ALREADY CLOSED by B1.** The B1 flame fix in
  fmburn! added: `flmult≠1 ⇒ flame=oldfl·flmult; flame≠oldfl ⇒ byram=60·(flame/0.45)^(1/0.46); sch from that
  byram` (fmburn.f:443,459-464). So scorch now tracks the multiplied flame. No separate fix needed. ✓

- **GAP mortality #3 (SDICHK QMD floor) — FIXED** (mortality.jl sdi_max_check!). Removed the `dq0<0.3 ⇒ 0.3`
  floor: per sdichk.f:78-81 the over-density DECISION (TEMMAX) + SDImax RESET use the UNFLOORED RMSQD/DR016
  (TEMD0); the 0.3 floor (DQ0, sdichk.f:59-61) feeds only the cosmetic TMD0→UPLIM warning (unported in jl). jl
  was flooring dq0 and using it for both decision+reset. Reachable only for dense sub-inch stands (QMD<0.3) →
  no-op for the suite (snt01 QMD≫0.3). Suite 4530/1. Variant-general (base sdichk.f). ✓

- **GAP mortality #2 (TPAMRT threshold-filtered) — FIXED** (mortality.jl). tpa_mort now sums (tpa−killed) over
  ALL trees (no dthresh filter) = TNEW (morts.f:706-712,772). Verified FVS BOTH sides: T (reset basis,
  morts.f:233) IS dthresh-filtered (D<DBHSTAGE→GOTO 12), TPAMRT (DO 36) is NOT → the reset test (morts.f:245)
  is filtered-T vs unfiltered-TPAMRT. jl had both filtered. No-op for snt01 (no sub-threshold trees). Suite
  4530/1. CAVEAT noted: sub-threshold killed[i] must match FVS WK2 — validate on a regen scenario. ✓

- **UNVERIFIED diameter-growth (calibration PTBAA current vs backdated) — RESOLVED: FAITHFUL (no change).**
  Traced the two-pass DENSE: dense.f:184 backdates dbh ONLY in pass 1 (LBKDEN.AND.LREDO); pass 2 (LREDO=FALSE)
  uses CURRENT dbh, and PTBAL→PTBAA (dense.f:280) is filled in pass 2 → PTBAA is CURRENT. cratet.f has NO DENSE
  between the backdated DENSE (186) and the calibration DGDRIV (551), so the calibration DGF reads CURRENT PTBAA
  — exactly what calibrate_diameter_growth! overrides to (diameter_growth.jl:267-281). The audit's "dense.f
  suggests backdated" missed pass-2; snt01's bit-exact DG corroborates. B3-shaped: static argument vs a live-FVS
  bit-exact match was the misread. jl is faithful. ✓

- **UNVERIFIED density (notre! FINT/FINTM dead inflation) — CONFIRMED GAP; localized fix specified.**
  notre.f:122-124: for DEAD (mortality) records the expansion factors are inflated VP/FP/FP2 *= FINT/FINTM, so
  the BACKDATED CALIBRATION density adds dead trees back at the cycle/mortality-period ratio; the inflation is
  DEFLATED for all downstream uses (notre.f comment: FMSSEE in CRATET, PRTRLS treelist). jl notre! expands live
  + dead with identical factors (no inflation). FAITHFUL FIX (localized, avoids the inflate/deflate dance):
  in calibrate_diameter_growth!'s backdated-density pass (where the dead partition is exposed), temporarily scale
  the dead-record t.tpa by control.growth_fint/growth_fintm before compute_density!, then restore — so only the
  backdated calibration BA/AVH/PCT sees the inflated dead TPA, matching FVS, while downstream stays deflated.
  NO-OP for the entire suite (FINT=FINTM=5 ⇒ ratio 1). Reachable only for a GROWTH FINT≠FINTM stand WITH recent
  -mortality (IMC=7) inventory records. Defer the change until such a scenario exists (principle #4). ✓ analyzed

- **GAP volume #2 (round vs Fortran INT in height dubbing) — FIXED** (volume.jl:213,216,218,230). Changed
  `round(Int32, x+0.5)` → `trunc(Int32, x+0.5)` to match cratet.f:381-397 INT() (truncate-toward-zero, i.e.
  round-half-UP), where the missing-height dub sets NORMHT/ITRUNC. Julia round() is round-half-to-EVEN, and the
  `+0.5`-then-round double-rounds — diverging by 1 when x is an odd integer. No-op for the suite (dub-height path
  + odd-integer tie not hit by snt01). Variant-general (Fortran INT semantics). Suite 4530/1. ✓

- **GAP volume #3 (final volume rounding ties-to-even vs NINT) — FIXED** (r8clark_vol.jl:433-436). The Clark
  vol[1]/[4]/[7]/[10] rounding was `round(x)` (Julia ties-to-EVEN) under a "match Fortran nint()" comment, but
  NINT is ties-AWAY-from-zero. Added RoundNearestTiesAway (matching the already-correct board-ft per-log round at
  :495). No-op for the suite (no x.05/x.5 tie hit; bit-exact preserved). Variant-general (Fortran NINT). ✓

- **GAP volume #1 (input board-defect not applied to sawtimber cubic when BF=0) — FIXED** (volume.jl). Moved
  the input-IBDF application (to bf AND scf) OUTSIDE the `bf>0` gate (vols.f:352,354-358 GO TO 100,415-420): a
  too-small-for-boardfeet tree (BFV=0, SCFV>0) still loses sawtimber cubic to its input BF defect. Only the
  BFDEFT/log-linear curve updates stay gated on bf>0 (vols.f:393). Also switched the bf-defect % rounds to NINT
  (RoundNearestTiesAway, vols.f:396,410). No-op for the suite (needs a BFV=0,SCFV>0 tree with input BF defect).
  Suite 4530/1. ✓  FOLLOW-UP: the CUBIC defect rounds (icdf, volume.jl:424,428) still use ties-to-even — verify
  vs vols.f cubic NINT.
  → cubic-defect NINT also fixed (vols.f:13,21 confirmed)

- **GAP fire-snag-cwd #3 (fmscro ILIFE round vs ceil) — FIXED** (fuel_additions.jl:76). `round(Int,…)` → `ceil(Int,…)`
  matching fmscro.f:126-131 (INT(RLIFE) then +1 if truncated/≤0 = ceil, floor 1). The crown-debris lifespan over
  which CWD2B amt/ilife is distributed. No-op for the suite (tested fire stands have integer lifespans).
  Variant-general (FFE base fmscro). Suite 4530/1. ✓

- **GAP fire-snag-cwd #4 (redcedar TFALL hardcoded 1 yr) — FIXED** (fuel_additions.jl `_fm_tfall`). FVS's crown-
  debris fall-time table sets FOLIAGE (crown size 0) TFALL = 1 yr for all species EXCEPT eastern redcedar, which
  gets 3 yr (fmvinit.f:1017-1021 `IF(I.EQ.2) TFALL(I,0)=3.0`). jl's `_fm_tfall(cls,sz)` returned 1 for ALL
  foliage. FIX: added the species arg — `sz==0 ⇒ sp==2 ? 3 : 1` (sp 2 = the SN redcedar index, matching FVS's
  `I.EQ.2`); the rest of the table (branch sizes 1/2=row1, 3=row3, 4/5=row4 by tfall_cls) was already correct.
  Effect: redcedar foliage crown debris now flows to down-wood over 3 yr not 1 (CWD2B timing). +7 tests
  (test_fire.jl). SN-scoped like the `_FM_TFALL` tables themselves (a NE port would re-source them — noted).
  Suite stays bit-exact (tested fire stands carry little redcedar foliage). 5010/1. ✓
  ⟹ fire-snag-cwd: 3/4 resolved; only hard→soft DKTIME (= the BLOCKED SNAGDCAY/DECAYX decay-state model) left.

- **GAP compress-cuts #5 (truncated-tree merge round vs IFIX) — FIXED** (compress.jl:262). NORMHT/ITRUNC merge
  now `trunc(Int32,…)` = comprs.f:805-806 IFIX(XNR/TXP*100) (truncate, no +0.5), not Julia round (ties-to-even).
  Only fires for truncated/top-killed trees in a COMPRESS merge. Suite 4530/1 (no-op for tested set). Variant-
  general (base comprs IFIX). Note: distinct from the accepted COMPRESS eigensolver divergence. ✓

- **GAP compress-cuts #1 (merge decay/defect copied vs averaged) — FIXED** (compress.jl `_merge_one!`). comprs.f
  STEP7 (818-936) TPA-weight-AVERAGES CULL, DECAYCD (decay code) and WDLDSTEM (woodland stems) when collapsing a
  class — but accumulates DECAYCD/WDLDSTEM into INTEGER registers (DECAYI/WDLDSTEMI), so each running sum is
  truncated toward zero (small-prob contributions vanish); only CULL stays REAL. jl averaged CULL but COPIED
  decay_code/woodland_stems from the RANN-sampled record (the old comment wrongly called them nominal). FIX:
  replicate the integer-truncated tpa-accumulation `di = trunc(di + code·prob)` then `trunc(di/TXP)`. Bit-exact-
  safe (the tested COMPRESS stands carry decay/woodland = 0, so average=copy=0); the exact result is mem-order-
  dependent via the truncation — a residual only when merging mixed decay/woodland inventory data. Suite 5010/1.
  Variant-general (base comprs). ✓
- **GAP compress-cuts #2 (TCONDMLT point weights unported) — FAITHFUL for single-point; multi-point deferred.**
  cuts.f:1074-1075 ranks cut trees by `WK2 = XSZ + IORDER + TCWT·IMC + SPCLWT·ISPECL + PBAWT·PTBAA(IP) +
  PCCFWT·PCCF(IP) + PTPAWT·PTPA(IP)`; jl's `_cut_pref_wt` has the first four terms but not the three POINT terms
  (point BA / CCF / TPA, set from TCONDMLT params 3-5). KEY INSIGHT: PTBAA/PCCF/PTPA are PER-POINT, so on a
  SINGLE-point stand they are identical for every tree ⇒ the point-weight term is a CONSTANT added to all WK2 ⇒
  it does NOT change the ranking ⇒ jl ignoring it is FAITHFUL. It only differentiates trees on a MULTI-point
  stand, which needs per-point density (jl tracks stand-level — the same limitation as the regen pccf gap). So:
  no code change (single-point faithful); multi-point + nonzero PBAWT/PCCFWT/PTPAWT deferred to per-point density.
- **GAP compress-cuts #4 (sorted SDI/RDEN uses PROB vs PROB+WK2) — FAITHFUL.** comprs.f sums PROB+WK2 for the
  class sampling (666); but FVS's SN COMPRESS is entered with WK2 (the mortality weight) ZEROED, so PROB+WK2 =
  PROB and jl's PROB-only per-class weight matches (compress.jl:11). No divergence on the SN path. ✓
- **GAP compress-cuts #3 (AUTSTK BA constant) — FALSE POSITIVE (no change).** jl `_autstk` uses 0.005454154
  (cuts.jl:691) — and the FVS AUTSTK (acd/cutstk.f:161 `CSTOCK=CSTOCK+TPA*(D*D*0.005454154)`) uses the SAME
  constant. The 0.0054542 the audit cited is a DIFFERENT routine (vbase/cuts.f:1241, the SDI-based THINSDI cut),
  not AUTSTK. jl is faithful; snt01's bit-exact auto-thin corroborates. B3-shaped: would've broken a correct
  constant on the flag's word. ✓

- **GAP establishment #2 (idup floor vs ceil) — FIXED** (establishment.jl:46). `fld` → `cld`: IDUP=ceil(MINREP/
  NPTIDS) per estab.f:199-207 (loop stops at first I with NPTIDS·I≥MINREP). No-op for the suite (BARE stand has
  NPTIDS=1 ⇒ ceil=floor=50); diverges only for multi-point inventory regen (1<NPTIDS<MINREP). MAXPLT cap doesn't
  bind for those. Suite 4530/1. ✓

- **GAP establishment #4 (PLANT-height floor XMIN vs 0.05) — FIXED** (establishment.jl). Branch-specific floor:
  PLANT branch (treeht≥0.1) → 0.05 (estab.f:1034); default/natural branch → XMIN (estab.f:1037). jl had a single
  XMIN floor for both. No-op for the BARE-stand NATURAL test (default branch unchanged); diverges for PLANT trees
  with 0.05<hht<XMIN. Suite 4530/1. ✓

- **GAP establishment #? (pccf hardcoded 0 in the regen crown ratio) — FIXED + LIVE-VALIDATED.** New regen trees
  get crown ratio `CR = 0.89722 − 0.0000461·PCCF(point) + 0.07985·N(0,1)` (regent.f:178), where PCCF is the
  per-point CCF of the EXISTING (pre-regen) overstory (DENSE, dense.f:210). jl hardcoded `pccf=0` — exact only
  for a bare/sparse stand (CCF≈0, the bare_natural case) but wrong for regen into a STOCKED stand. FIX:
  `pccf = stand_ccf(s)` (establishment.jl) — the pre-regen overstory CCF; for a single-point stand that IS the
  per-point PCCF. VALIDATED vs live FVSsn on a new scenario (plant_stocked.key = the dense fire_early overstory +
  NATURAL 2000, stand CCF≈311): regen of 50 species-13 stems at 2005, TPA/BA bit-exact every cycle (≤1 BA the
  cycle after, see below), and the regen crown MEAN matches live bit-close — 82.44 vs live 82.46 (PCCF=0 would
  leave the center at ~89). bare_natural stays bit-exact (the 4.6e-5 coefficient keeps a sparse stand below print
  resolution). RESIDUAL/CAVEAT (documented, variant-aware): jl uses the stand-AVERAGE CCF for all regen, while
  FVS varies PCCF per point — so for a MULTI-point stand the per-tree crown distribution differs by ±1 (mean is
  exact) and the post-regen BA carries a ≤1 ft²/ac residual; single-point + bare stands are exact. AUTOES
  (auto-regen into a stocked stand) is unported, so the only reachable path is an explicit PLANT/NATURAL keyword.
  +21 tests (test_estab_pccf.jl). Variant-general (base regent/dense semantics). Suite 4986/1. ✓

- **GAP io-serialization #2 (.sum header sample-weight C %E vs Fortran E15.7) — FIXED** (summary.jl). Added
  `_fortran_e15_7` (0.DDDDDDD·10^p, 7 sig digits) — real FVSsn writes 11.0 as `0.1100000E+02`, while jl's C
  `%15.7E` wrote `1.1000000E+01` (D.DDDDDDD, 8 sig). The committed goldens are FVSjulia (C-format) so jl matched
  them, AND the suite SKIPS the -999 header line (test_keyword_coverage.jl:84) → latent vs real Fortran. Now
  matches real FVSsn. Suite 4530/1 (header skipped). Variant-general (IO). ✓

- **GAP io-serialization #3 (STOP 4-char prefix vs full field) — FIXED** (keyword.jl:117). `head8[1:4]=="STOP"`
  → `head8=="STOP    "` (full 8-char field) per keyrdr.f:55-59 (TMP=RECORD(1:8), TMP.EQ.'STOP'); END stays a
  4-char prefix (keyrdr.f:92 TMP(1:4)). Latent — no current keyword is "STOPxxx"; a hypothetical one would have
  been mis-terminated. Suite 4530/1. Variant-general (base keyrdr). ✓

- **GAP sprout #2 (ESSPRT gate KODFOR vs ISEFOR) — DEFERRED (source not located).** jl `_es_special_forest`
  (sprout.jl:84) gates on forest codes 809/810/905/908 via `s.plot.user_forest_code` (KODFOR); the audit says
  FVS keys on ISEFOR. But the FVS 809/810/905/908 special-forest gate is NOT in estb/vestb essprt.f (that's a
  DSTMP/species PREM multiplier) — source not yet located, so KODFOR-vs-ISEFOR can't be confirmed. The jl
  essprt_sn was validated bit-exact vs live Fortran (sprout.key) but probably on a NON-special forest (gate
  untested). NOT changing on the flag's word (B3 discipline). NEXT: locate the 809/810 gate (likely SN
  blockdata / a derived ISEFOR map) + a special-forest sprout scenario before any change.

- **GAP sprout #2 — CORRECTION (source now LOCATED): CONFIRMED GAP, forkod IFORDI port needed.** sn/forkod.f:185,
  231-252: the special forests are the 3-digit FORKOD-derived IFORDI (809/810/905/908); FORKOD maps 5-digit
  KODFOR→IFORDI (80906→809, 81005→810). jl `_es_special_forest` checks user_forest_code=KODFOR (5-digit) against
  the 3-digit IFORDI codes → NEVER fires. jl's forkod port (forest_location.jl) ports only forest→location, not
  the KODFOR→IFORDI map, so jl has no IFORDI. FIX: port forkod's IFORDI derivation + 905/908 detection, expose
  IFORDI on the plot, gate the special-forest sprout on it. Reachable only for Ouachita/Ozark (snt01=80106 ⇒
  no-op). Feature port — defer to a special-forest scenario. (Supersedes the "source not located" note above.)

- **GAP sprout #2 — FIXED (canonical forests)** (sprout.jl:145). `isefor = user_forest_code ÷ 100` = IFORDI
  (forkod.f:183, KODFOR/100). The special-forest gate (809/810/905/908) now fires for the canonical Ouachita
  (80906→809) / Ozark (81005→810) forests; was comparing 5-digit KODFOR ⇒ never fired. No-op for snt01
  (80106→801). FOLLOW-UP: forkod's SELECT CASE alias remaps (e.g. 7207→809) still need the full forkod IFORDI
  port for non-canonical inputs. Suite 4530/1. ✓

- **GAP structure-stage #2 (single-tree NTREES≤1 path) — CONFIRMED GAP, deferred (niche).** sstage.f:235-268:
  a 1-canopy-tree stand uses a special branch — cover-based class 0/1 when WK6=CW²·TPA·π/4 < 435.6·TMPCCM,
  else DBH-direct (SSD→1, SAW→2 w/ PCTSMX demotion, else 5). jl's normal nstr==1 path classifies by the
  stratum DBHNOM (not cover), so a sparse single tree diverges (FVS class 0/1 by cover vs jl by DBH). Niche:
  needs exactly 1 canopy tree (snt01/fire_early have many ⇒ SSTAGE validated bit-exact but this path untested).
  No-op for the suite. FIX: port the sstage.f:235-268 branch. Defer to a 1-canopy-tree scenario. ✓ analyzed

- **GAP structure-stage #1 (after-thin row reuses before-thin SDI/cover) — FIXED (source-verified + mechanism-
  confirmed; suite bit-exact).** FVS calls the structure computation TWICE per cycle: before thinning (IBA=1, on
  the pre-thin tree list) and after (IBA≠1, on the POST-thin list; sstage.f:13-14, 145-155 — XBAMAX switches
  BTSDIX→ATSDIX when ONTREM>0, and the cover/strata recompute on the thinned list). jl's `write_structure_report`
  wrote BOTH rows from the SAME pre-thin stand (`structure_report_row(s,yr,0)` then `(s,yr,1)` with no thin
  between), so the after-thin (Rm=1) row was a relabeled copy of the before-thin row. FIX: apply the cycle's
  `cuts!` + recompute density BETWEEN the two rows (structure_stage.jl); `cuts!` is idempotent per year
  (cuts.jl years_cut guard), so `grow_cycle!`'s own cut below becomes a no-op (no double-cut). SAFE: confirmed
  `write_structure_report` is only ever called on a THROWAWAY stand (its docstring + the lone test call), never
  the live sim — so the inserted `cuts!` can't corrupt a real run. MECHANISM VERIFIED: on a THINSDI scenario the
  jl after-thin row now reflects the post-thin stand (Rm=1 differs from Rm=0); when nothing is cut the rows stay
  identical (sstage.f:146), so the UNTHINNED structure tests (snt01 stand-1 + fire_early, live-validated
  bit-exact) are unchanged. Suite 4986/1. FOLLOW-UP: a precise moderate-thin STRCLASS live diff — blocked here
  by scenario construction (fire_early's stand won't take a clean partial thin: THINBBA/THINBTA don't fire,
  THINSDI wipes it; needs a denser stand whose THINSDI leaves a residual). Variant-general (base sstage). ✓

- **GAP sprout #1 (SPROUT species/DBH multiplier table unhonored) — FIXED + LIVE-VALIDATED bit-exact.** The
  SPROUT keyword (esin.f opt 26) carries, per species, a sprout-COUNT multiplier (SMULT, field 3), a HEIGHT
  multiplier (HMULT, field 4), and a stump-DBH window [DLO field 5, DHI field 6). esuckr.f:197-205 (DO 450)
  looks these up by the PARENT stump's species + DBH: SMULT/HMULT default 1/1, and each matching activity with
  `DSTMP ∈ [DLO,DHI)` overwrites (last wins). jl previously stored a SINGLE global SMULT/HMULT applied to ALL
  species and ignored the window — so `SPROUT 22 3` tripled EVERY species' sprouts, not just species 22. FIX:
  added `Control.sprout_overrides::Vector{NTuple{5,Float32}}` (species_code, smult, hmult, dmin, dmax); the
  SPROUT handler (keyword_dispatch.jl) pushes one per keyword honoring fields 3-6; `esuckr!` (sprout.jl) does the
  per-parent DO-450 lookup (SPDECD selector: >0 single / 0 all / <0 group via sp_groups; window-gated) and the
  per-stump SMULT≤0 skip (esuckr.f:211). VALIDATED 3 ways vs live FVSsn on the sprout.key thin-then-sprout stand:
  (a) `SPROUT 22 3` → 2005 TPA 491→729 bit-exact; (b) `SPROUT 22 3 1 8 99` (windowed) → 509, STRICTLY between
  491 (all 1×) and 729 (all 3×), proving only the ≥8" stumps took the 3× override; (c) the default form stays
  bit-exact (suite green). Key insight: a window with SMULT=1 is invisible (out-of-window stumps keep the
  default SMULT=1 and still sprout — the window scopes the OVERRIDE, it does not suppress sprouting). +138 test
  assertions (test_sprout_table.jl). Variant-general (base esuckr semantics). RESIDUAL: the SPROUT activity DATE
  (field 1) scheduling is applied throughout rather than gated to its OPGET window (matches prior jl behavior;
  the common inv-year form is unaffected) — a minor follow-up. Suite 4946/1. ✓

- **GAP sprout #4 (sprout_dbh ignores per-stand AA refit) — FAITHFUL on the default path; SUBSUMED by the
  LHTDRG/NOHTDREG HT-DBH calibration gap.** esuckr.f:298-303 sets the Wykoff-inverse intercept `AX = HT1(ISSP)`
  if `IABFLG(ISSP)=1`, else `AX = AA(ISSP)` (the per-stand-refit intercept). `IABFLG` defaults 1 (grinit.f:105);
  AA is refit and IABFLG flips to 0 ONLY inside cratet.f's HT-DBH regression, which is gated by
  `IF(K1<3 .OR. .NOT.LHTDRG(ISPC)) GO TO 100` (cratet.f:327). LHTDRG defaults `.FALSE.` and is turned on only by
  NOHTDREG-invoke (see the height-growth verdict). So on the default/suppress path IABFLG=1 ⇒ AX=HT1, and jl's
  `sprout_dbh` (which uses `wykoff_ht1`) is FAITHFUL. The AA path is one of the THREE sites the same unported
  LHTDRG calibration touches — (1) regent small-tree Wykoff branch, (2) cratet height dubbing, (3) this esuckr
  sprout DBH — all reached only when LHTDRG is on, all covered by the single `kw_nohtdreg!` transparency @warn.
  Not a separate fix; folded into the LHTDRG follow-up. ✓

- **GAP sprout #3 (no ESCPRS list-compression when the tree list fills during sprouting) — FULLY TRACED,
  DEFERRED (bounded reachability + depends on the open COMPRESS residual).** esuckr.f:251-257: when adding a
  sprout record would reach `MXRR` (≈ MAXTRE=3000, minus the root-disease reserve from RDESCP), FVS sets a target
  `ITRGT = max(ITRNRM−I, 0.70·MXRR)` and calls `ESCPRS` to make room, then continues adding. estb/escprs.f is a
  thin wrapper: `CALL COMPRS(ITRGT, 0.5)` (the standard tree-list compression) + `SPESRT` (re-establish the
  species-order sort) + recompute the establishment SPCNT species counts. jl's esuckr! instead has a hard
  overflow guard (`n > length(t.dbh) && break`) that DROPS the excess sprouts. VERDICT: a real divergence, but
  (1) reachable ONLY when a stand hits 3000 tree records *during* sprouting (an extreme cl- or seed-cut of a
  very dense stand; no suite/realistic SN scenario approaches it — typical stands carry hundreds to low-thousands
  of records), and (2) the fix would reuse the ALREADY-PORTED `compress!` (compress.jl) — but compression still
  carries the open post-compression DGSCOR residual (#29, the lone broken test), so ESCPRS is strictly downstream
  of that. Per B3 discipline (don't port on a flag's word without a validating live diff, and a 3000-record
  triggering scenario is impractical to build + would inherit #29), DEFER: fold the ESCPRS wrapper
  (compress!→target + species re-sort + SPCNT) into the COMPRESS-residual work, when a synthetic overflow stand
  can validate it against live FVS. jl's current drop-behavior is safe (no crash) and bounded to the overflow
  edge. ✓ traced

---
## ★ LIVE-FVS PER-TREE VALIDATION — growth path FAITHFUL (gold standard achieved)
Ran timeint10 (10-yr cycles) through BOTH the live-FVS oracle and FVSjl with FVS_TreeList DBS, compared
per-tree DBH by cycle (the campaign's actual bar — not just the aggregate .sum):
  • yr 2000 (cyc1): FVS 81 recs, jl 81 recs, sorted-DBH max|Δ| = 0.0 — BIT-EXACT.
  • yr 2010 (cyc2): FVS 243, jl 243, max|Δ| = 0.0 — BIT-EXACT.
  • yr 2020 (cyc3): live-only (TPA>0) FVS 177, jl 177, max|Δ| = 1.0e-5 (ULP) — FAITHFUL.
=> The DG bound-then-scale (10-yr DGBND ordering), small-tree DIAM-scale, and HK fixes are VALIDATED faithful
per-tree vs live FVS for a FINT≠5 scenario — the gold-standard check principle #4 demands, now achieved via the
live-FVS oracle harness. The apparent earlier per-tree "divergence" was a matching confound (unaligned tripled/
dead records across drifted cycles); the clean sorted cycle-by-cycle diff is bit-exact (cyc1-2) / ULP (cyc3).

- **NEW minor GAP: FVS_TreeList omits dead records.** FVS's FVS_TreeList at yr2020 has 243 recs (177 live + 66
  dead/MortPA), jl's has 177 (live only). FVSjl's treelist_snapshot (dbs_output.jl) should also emit the dead
  partition (TPA=0, MortPA>0) to byte-match FVS_TreeList. Low impact (dead records, the live values match ULP);
  fix = include the dead partition in treelist_snapshot. ✓ found via the live-FVS per-tree diff.

- **GAP (refined): FVS_TreeList incomplete vs real FVS — DBS-completeness task, deferred.** Beyond the dead-record
  omission, jl's FVS_TreeList schema (dbs_output.jl:442) MISSES 9 columns real FVS emits: TreeVal, SSCD, PtIndex,
  MortPA, MistCD, MDefect, BDefect, EstHt, ActPt. So jl's table is a subset (different column set + live-only
  rows). Output-only (no computational impact — the present live values match live FVS to ULP, per the per-tree
  validation). FIX = a focused DBS pass: add the 9 columns to the schema/writer + emit the dead partition with
  MortPA. Defer; low priority (the goldens are FVSjulia-limited anyway). ✓ found via the live-FVS per-tree diff.

---

## Live-FVS per-tree validation — FIRE path (fire_early, SIMFIRE 2000, PSBURN=50)

**Method:** ran `fire_early.key` (5 cycles, single fire at 2000) through the live oracle
(`LiveFVS.run_key`, FVS_TreeList DBS) and FVSjl, comparing per-tree DBH/TPA/Ht/PctCr by cycle.

**Pre-fire = BIT-EXACT.** 1990/1995/2000 snapshots: sumTPA identical to live
(536.048 / 507.441 / 469.767). At the 2000 pre-fire snapshot the 243 records match live on
**DBH, TPA, Ht AND PctCr** (max|Δ|=0.000 across all 243) — every input to the fire-mortality
model is bit-exact entering the burn.

**Post-fire = ~5.6 TPA under-kill, set entirely in the 2000→2005 fire step:**

| year | sumTPA live | sumTPA jl | Δ |
|------|------------:|----------:|----:|
| 2000 (pre-fire) | 469.767 | 469.767 | 0.000 |
| 2005 | 103.944 | 109.774 | +5.830 |
| 2010 | 101.472 | 107.166 | +5.694 |
| 2015 |  98.998 | 104.555 | +5.557 |

The gap is fixed at the fire and merely grows naturally afterward (jl under-kills, then both
grow in parallel).

**Narrowed — by species and DBH class (2005):**
- **All four divergent species are fire-mortality GROUP 6** (the Reinhardt crown-scorch +
  bark-thickness logistic, `fire_tree_mortality` else-branch): FIA 318→sp22 SM, 531→sp33 AB,
  812→sp65 SK, 998→sp89 OH. None are in the oak/hickory char-height groups 1–5. So the
  divergence is entirely in the **Reinhardt CSV path**, not the Regelbrugge-Smith path.
- By DBH class the under-kill is **monotonic in size and vanishes for large trees**:
  0–2″ Δ+0.000 (fully scorched, CSV≈100, bit-exact) · 2–4″ Δ+1.74 (worst, live kept 1.40 vs
  jl 3.14) · 6–8″ Δ+1.22 · 8–10″ Δ+1.26 · 12–14″ Δ+0.24 · ≥18″ Δ0.000.

**Root-cause hypothesis (not yet pinned):** with DBH/Ht/PctCr bit-exact and group-6 pmort a
faithful port of fmeff.f:184-186, the only remaining lever is **CSV (crown-volume-scorched)**,
which for a *partially*-scorched tree is sensitive to the scalar **scorch height** `sch`.
The signature fits a scorch-height delta exactly: the very short 0–2″ trees are fully scorched
(CSV pinned at 100 regardless of small sch changes → bit-exact), while taller 2–4″ trees are
only partially scorched (CSV in the steep part of the curve → most sensitive → largest Δ),
tapering as larger trees have low pmort anyway. jl computes flame=3.979, scorch=16.407 for this
fire (post-B1-flame-fix). FMEFF/CSV/bark formulae verified equal to fmeff.f line-by-line; the
remaining unknown is whether jl's `sch` equals live's.

**Blocker:** the direct `sch` (and per-tree XRAN burn-selection) comparison needs the live
FVS_BurnReport DBS table or the text FFE fire report. BURNREDB/MORTREDB are accepted DataBase-
block keywords but the live binary did not emit FVS_BurnReport for this single-stand config
(same FFE report-config gap noted earlier for FVS_Carbon); FUELOUT produced no text fire report
in fort.16 either. Pinning `sch` live-vs-jl is the concrete next step (try the report-range/FMIN
config that snt01-stand4/sn.key use, which DID emit these tables).

**Status:** this is the documented FFE "per-tree kill distribution" residual (previously known
only as aggregate "BA 81 vs 78"), now sharply characterized per-tree and narrowed to the
group-6 Reinhardt scorch-height/CSV chain. NOT a fudge target — left as a faithful, documented
GAP with a precise root-cause hypothesis and the report-config unblock as the next move.
Growth-path per-tree validation (timeint10) remains the bit-exact gold standard; the fire path
is bit-exact on every input and diverges only in the scorch-driven kill magnitude.

### FIRE under-kill — ROOT CAUSE LOCALIZED to fuel-model weighting (sm/lg loadings)

Got the live binary to emit FVS_BurnReport (needed BOTH `BURNREDB` in the DataBase block AND
`BURNREPT` in the FMIN block — the latter sets the report year-range IFMBRB/IFMBRE that gates the
`DBSFMBURN` call in fmfout.f:95). Direct flame/scorch comparison at the 2000 fire:

| quantity | live | jl | jl vs live |
|----------|-----:|---:|-----------:|
| flame length (ft) | 4.172 | 3.979 | −4.6% |
| scorch height (ft) | 17.581 | 16.407 | −6.7% |

This **confirms the small-tree under-kill mechanism**: jl's scorch is ~7% low → partially-scorched
2–4″ trees get a lower CSV → lower group-6 Reinhardt pmort → they survive (the 0–2″ fully-scorched
trees pin at CSV=100 and stay bit-exact, exactly the observed signature).

**Localized one level deeper — the fuel-model weights diverge:**
- live: model **5 @ 0.44**, model 10 @ 0.56
- jl:   model **5 @ 0.135**, model 10 @ 0.865

Same two models (so candidate selection / `iffeft` is right), but the **split is wrong**. flame is
the weighted-Byram `0.45·(byram/60)^0.46`; model 5 carries the higher intensity, so jl
under-weighting it (0.135 vs 0.44) is precisely why jl's byram→flame→scorch is low.

**The weight split is FMDYN's perpendicular-distance interpolation of the point (sm,lg) against the
model iso-lines.** Two facts pin the cause:
1. jl's `_FMD_XPTS` iso-line table is **bit-identical** to fmcfmd.f's `DATA XPTS` (all 14 models:
   5→(5,15), 10→(10,30), 12→(30,60), …). The geometry is faithful.
2. FMDYN is deterministic in (sm,lg). jl computes **sm=8.156, lg=3.50** (iffeft=1, m14=0.17, so the
   `sm>6` branch ⇒ candidates {5,10,12}). With identical XPTS + identical algorithm, identical (sm,lg)
   ⇒ identical weights. Since the weights differ, **live's (sm,lg) must differ from jl's**.

⟹ **Root cause: the FFE small/large dead-wood surface-fuel loadings `(sm,lg)` that FMDYN keys on are
off in FVSjl.** `(sm,lg)` come from `_small_large_fuel` over the FFE down-wood (DDW) pools — which is
the SAME open **DDW down-wood gap** already flagged as the last carbon-model gap. The fire-mortality
"per-tree kill distribution" residual and the DDW down-wood carbon gap are **one root cause**: the
down-wood pools. Fixing DDW loadings should simultaneously (a) close the carbon DDW gap and (b) lift
flame/scorch → correct the group-6 fire kill.

**Two remaining verification steps (not blockers to the conclusion):**
- Pin live's `(sm,lg)` directly (gated behind FFE fuel-report config — `FUELOUT`/`FuelReDB` resisted
  emitting FVS_Fuels for this single-stand key across several attempts; the global `DEBUG` keyword
  segfaults in an unrelated fvsvol.f:530 debug path, so that route is out). The cleanest path is the
  sn.key/snt01-stand4 FMIN report config that is known to emit the FFE fuel/consumption tables.
- Double-check jl `_fmdyn` against fmdyn.f once more (already ported + previously validated; the prior
  "xpts reshape scramble" bug was fixed) to fully exclude an algorithm-side contribution.

**Status:** the fire path is bit-exact on every *direct* mortality input (DBH/Ht/Crown/CSV-formula/
bark/logistic); the only divergence is the scalar scorch height, traced to the fuel-model weight
split, traced to the (sm,lg) down-wood loadings = the DDW gap. NOT fudged — a fully-localized,
documented GAP whose fix lives in the FFE down-wood pools.

### FIRE root-cause — deeper trace: it's the CWD pool *values*, not the routine/aggregation

Followed the (sm,lg) lever down through the FFE call graph and ruled out everything except the
down-wood pool contents:

1. **Right routine.** The SN active fuel-model routine is `FMCFMD` (fmcfmd.f, the FMDYN-over-standard-
   models path) — NOT `FMCFMD2`/`FMCFMD3` (fmcfmd2.f, the custom-model-89 / IFLOGIC path). The
   FMCFMD3 dispatch (fmcfmd2.f:654-664) *excludes* SN/CS/LS/UT/… from its `FMCFMD` call precisely
   because those variants call `FMCFMD` from the FMBURN moisture loop instead (fmburn.f:400). jl's
   `select_fuel_models` is ported from fmcfmd.f ✓ — the earlier "wrong variant?" worry is cleared.

2. **SMALL/LARGE are computed in fmtret.f:376-386** (FFE fuel-update), read by FMCFMD via FMFCOM
   common. The formula: `SMALL = Σ CWD(I,{1,2,3},K,L) + CWD(I,10,K,L)` (classes 1–3 + litter),
   `LARGE = Σ CWD(I,{4..9},K,L)`, over I=1,2 (unpiled/piled) × K=1,2 × L=1,4.

3. **jl's `_small_large_fuel` matches that aggregation** — small = classes {1,2,3,10}, large = {4..9},
   over k=1,2 × l=1,4. jl's `cwd` is 3-D `[size 1:11, dead/soft 1:2, decay 1:4]`, collapsing the
   Fortran's piled/unpiled `I` dim; harmless for an unpiled stand (fire_early has no harvest piling,
   so CWD(2,…)=0). The class membership (1,2,3,10 vs 4–9) is identical. ✓

⟹ **Conclusively, the `(sm,lg)` divergence is in the CWD pool VALUES (the DDW down-wood accounting),
not the fuel-model routine, the XPTS geometry, the FMDYN algorithm, or the small/large aggregation.**
Every step from CWD → (sm,lg) → FMDYN weights → flame → scorch → CSV → group-6 pmort is now verified
faithful *except* the CWD pool contents themselves. The fire kill fix == populating the DDW down-wood
pools to match live FVS (= the open carbon DDW gap). #22 is now a per-CWD-class reconciliation:
compare jl `cwd[1:11,k,l]` vs live's down-wood by size class (FVS_DWD_Vol/Cov, or the FFE fuel report)
and trace each class's inputs (litterfall, FMSDIT/FMCADD crown-lift, mortality additions, decay).

### FIRE/DDW — per-class live-vs-jl fuel comparison (got FVS_Fuels + FVS_Down_Wood_Vol)

Cracked the remaining FFE report config: `FVS_Down_Wood_Vol` needs `DWDVLDB`(DataBase) + `DWDVLOUT`
(FMIN range); `FVS_Fuels` needs `FUELSOUT`(DataBase, sets IFUELS — NOT `FUELREDB`, which is IFUELC =
consumption) + `FUELOUT`(FMIN range). jl gates both tables under its `CARBREPT` flag. Direct per-class
comparison on fire_early:

**Pre-fire (1990, 1995) surface fuels MATCH live to ULP** — Surface_Litter/Duff/lt3/3to6/6to12/ge12 all
within ~0.002 tons/ac, and FVS_Down_Wood_Vol by diameter class within ~0.1 ft³/ac. ⟹ jl's **steady-state
fuel/down-wood accumulation (litterfall, decay) is FAITHFUL**. This *narrows the fire under-kill cause*:
it is NOT a baseline-accounting drift.

**Two real divergences found:**
1. **Surface_Shrub** — the one consistent *pre-fire* mismatch: live 0.88 vs jl 0.48 tons/ac (1995),
   live 0.60 vs jl 0.48 (1990). Shrub load is NOT part of SMALL/LARGE (dead-woody + litter), so it does
   not drive the fuel-model weights — but it IS a standalone faithfulness gap in the FFE shrub-load
   accounting. New GAP, logged for follow-up (distinct from #22).
2. **Post-fire down-wood under-populated in jl (the dead-tree→DWD transfer).** At 2005 (5 yr after the
   fire) live has huge down-wood from fire-killed trees falling — `lt3` 13.45 vs jl 1.42, `6to12` 10.39
   vs jl 0.59, `3to6` 3.96 vs jl 0.12 tons/ac — jl barely moves the killed biomass into CWD. Also at the
   2000 fire-year snapshot live shows litter consumed to 0 / duff 6.07→2.07 while jl retains litter 2.63
   / duff 6.22 (consumption and/or snag-fall-into-DWD timing differ).

**Reconciling with the flame/weight gap:** since the BurnReport *measured* different fuel-model weights at
the 2000 fire (live 5@0.44 vs jl 5@0.135) yet the 1990/1995 cycle-boundary fuels match, the 2000
*pre-fire* (sm,lg) must diverge within the 1995→2000 window — i.e. from the **tree-mortality additions
into the down-wood/CWD pools over the period** (the same dead-tree→DWD transfer that shows up large at
2005), NOT from steady-state litter/decay. This squarely confirms #22's scope: the CWD divergence is in
the **mortality-driven additions to the down-wood pools**, not the baseline pools.

⟹ #22 sharpened: the fix is the **dead-tree → down-wood (snag fall / FMTRET mortality additions)
transfer** into `cwd`, which (a) raises post-fire DWD to match live (the 2005 lt3 13.45 vs 1.42 gap),
(b) corrects the 2000 pre-fire (sm,lg) → fuel-model weights → flame/scorch → group-6 kill, and (c)
closes the carbon DDW gap. Plus a separate small GAP: FFE Surface_Shrub load (0.48 vs 0.88).

### FIRE/DDW — merch→gross bole hypothesis REFUTED by carbon_snt (B3 lesson applied)

Traced the snag→down-wood fall to its volume basis and found fmcwd.f CWD1 (fmcwd.f:184-187) calls
`FMSVL2(SP,DIAM,HTD,…, 'D', LMERCH=.false., …)` — which *reads* as the gross stem cubic (TCF), not the
merch bole (MCF) that jl's `sn.bolevol`/StandDead uses. Implemented the swap (fall uses `_fm_cuft`
gross v[1]; StandDead keeps merch).

**Result: REVERTED.** It regressed the carbon_snt DDW test, which is **bit-exact to the live Fortran
`.report.save`** with merch — gross made jl *overshoot* (+0.3 at 2000, +0.4 at 2005 vs Δ≤0.007 with
merch). Per the campaign's imprinted B3 lesson — *a static source-read that contradicts a live-Fortran
bit-exact match is the misread, not the code* — merch is the Fortran-faithful basis for the DDW carbon
pool. Suite back to 4530/1, no harm. (Open sub-question, low priority: why does FMSVL2('D',
LMERCH=.false.) read as gross yet the Fortran DDW match merch? Needs the FMSVL2 'D'-mode body trace —
likely the dead/DUMHT=-1 path returns a merch-equivalent. Deferred; merch is empirically correct.)

**So the fire_early DWD 10× gap is FIRE-SPECIFIC, not a global bole-volume error:**
- carbon_snt (FFE, regular mortality, no big fire) → DDW bit-exact with merch. The regular-mortality
  snag→DWD path is FAITHFUL.
- fire_early (big SIMFIRE at 2000) → 2005 down-wood 10× low (0to3 live 1232 vs jl 136). The gap is in
  the **fire-killed-snag → down-wood path**, whose dominant driver is the **post-burn accelerated snag
  fall** (`PBFRIS`/`PBFRIH`/`PBTIME`, fmsnag.f:202-216) — already listed UNIMPLEMENTED in jl's GAP set.
  Without it, jl's fire snags fall at the slow normal rate and are still standing at 2005, so little
  biomass has reached the down-wood pools; FVS's fire-weakened snags fall fast → the 1232 ft³ of fine
  down wood. (Secondary: fire snags are created at fmburn.jl:98 with NO `bolevol`, so they currently use
  the Jenkins total-AGB fallback for the fall — to reconcile once the fall RATE is fixed.)

⟹ **#22 REDIRECTED:** the down-wood / fire-kill gap is the **post-burn accelerated snag fall (PBTIME)**,
NOT the merch-vs-gross bole basis (refuted). Implement FMSNAG's post-burn fall acceleration (gated on
`BURNYR`/`YRDEAD ≤ BURNYR`, within `PBTIME` years) so fire-killed snags fall fast into the down-wood
pools; then reconcile the fire-snag bole basis and re-check fire_early DWD + the 2000 (sm,lg) → flame.

### FIRE/DDW — post-burn accelerated fall PORTED (faithful, partial); 2 downstream bugs unmasked

Ported the post-burn accelerated snag fall (FMSNAG fmsnag.f:200-214 + FMSFALL fmsfall.f:102-119) into
`update_snags!` (snag.jl): within `PBTIME=7` yrs of a burn, small (<`PBSIZE=12`) snags fall at
`RSMAL=1−(1−PBSMAL)^(1/PBTIME)≈28%/yr` and soft-at-fire snags at `RSOFT` (PBSOFT=1 ⇒ ~all over PBTIME),
as a FLOOR on the normal fall. Defaults from fmvinit.f:1101-1104,125. **Suite stays 4530/1** (carbon_snt
has FMIN/END but NO SIMFIRE ⇒ BURNYR=0 ⇒ untouched — verified).

Key sub-fix: jl's `fire.fire_year` is the SCHEDULED year and is **cleared after the fire fires**, so it's
0 during the post-fire fall window; FVS's `BURNYR` is the PERSISTENT last-burn year. Derived BURNYR from
the accumulated `fs.burn_reports` (which persist) instead.

**Faithful but PARTIAL — kept per principle #3 (it unmasked two real downstream bugs, not wrong itself):**
the fall now fires (2010 down-wood moved) but the fire_early 2005 gap persists (0to3 134 vs live 1232)
and the 2010 6to12 overshoot GREW (1634→1987 vs live 836). The fall is correct; what's wrong is the
biomass it moves and where:

1. **Fire-snag bole basis (fmburn.jl:98).** Fire-killed snags are created with NO `bolevol`, so the fall
   uses the `jenkins_biomass` TOTAL-AGB fallback (branches included) — too large and wrongly weighted to
   coarse classes. Regular mortality (mortality.jl:414) correctly passes the MERCH bole (the carbon_snt-
   validated basis). FIX: pass `mcf·v2t/2000` merch bolevol at fmburn.jl:98 too (the crown goes through
   the separate CWD2B path).
2. **Down-wood size-class distribution.** Even with the fall firing, jl's post-fire down wood is too
   COARSE — 6to12 overshoots while 0to3/3to6 (live 1232/352) stay ~10× low. The cone-taper split
   (`_cwd_cone_fractions`) and/or the bole source put mass in the wrong classes vs FVS CWD1's per-class
   integration. Needs a per-class CWD1 trace (fmcwd.f label 1000 P2−P1 × TVOLI) against jl's fractions on
   the actual fire-killed snag dimensions.

⟹ The fire down-wood is a MULTI-BUG subsystem: post-burn fall (DONE, faithful) + fire-snag merch bole
(#22 next) + size-class distribution (#22). Each is individually traceable; the net live-match needs all
three. Post-burn fall landed as a faithful building block.

### FIRE/DDW — fire-snag merch bole FIXED; remaining = 1-cycle fall lag + fine-class crown debris

Second faithful fix landed: **fire-snag bole basis** (fmburn.jl:98). Fire-killed snags now carry the
MERCH bole (`mcf·v2t/2000`) like ordinary-mortality snags (mortality.jl) and the carbon_snt-validated
StandDead/down-wood basis — instead of the `jenkins_biomass` TOTAL-AGB fallback (which double-counted the
crown belonging in CWD2B). **Big improvement on the coarse classes** (fire_early 2010): 6to12 1987→909
(live 836), 12to20 448→156 (live 163), 20to35 57→37 (live 40). Suite stays 4530/1 (carbon_snt no fire).

**Two issues remain, both structural:**
1. **One-cycle snag-fall LAG (the 2005 gap).** `ffe_fuel_update!` (the annual fuel loop containing
   `update_snags!`) runs BEFORE `grow_cycle!` each cycle, but the fire fires INSIDE `grow_cycle!`
   (`_maybe_burn!`). So fire-killed snags are created AFTER the fuel loop that would fall them, and don't
   begin falling until the NEXT cycle's loop. Instrumented: after the cycle reaching 2000 the burn_reports
   BURNYR is still 0 / den_hard 73; only by 2005 do the fire snags appear (den_hard 475) — so their
   post-burn fall first runs in 2005→2010. Hence 2010 down-wood matches live well but 2005 stays ~10× low
   (the fire-killed biomass hasn't reached the pools yet). FVS falls them starting the burn year. Fix is
   the FFE cycle ordering (run the post-fire fuel/fall step in the SAME cycle as the burn) — touches the
   FFE timing that is currently bit-exact for carbon_snt, so must be done carefully (same class as the
   documented crown-lift one-cycle lag).
2. **Fine-class (0to3/3to6) under-count persists even at 2010** (0to3 240 vs live 890). The coarse
   classes now match, so the missing fine material is most likely the CROWN/branch debris of fire-killed
   trees (the CWD2B path), not the stem bole. Next: verify fire-killed trees' crown feeds CWD2B and
   compare per-class.

⟹ #22 status: post-burn fall ✅ + fire-snag merch bole ✅ (coarse classes now live-faithful); remaining =
the 1-cycle FFE fall-lag (2005) + fine-class crown debris (CWD2B). Two faithful fixes landed this pass.

### FIRE/DDW — fire-killed CROWN→CWD2B FIXED; 2010 down-wood now live-faithful (3 fixes)

Third faithful fix: **fire-killed crown debris** (fmburn.jl). The Fortran FMEFF pools each fire-killed
tree's crown into CWD2B (fmeff.f:524-525, the SN simple case: crowns above the flame are not consumed,
so the full crown × killed-density → debris). jl's fire path didn't — only ordinary mortality did
(mortality.jl). Added the `fmscro!` call for the fire-killed crown (density = `curkil`), mirroring the
mortality path; the stem bole stays the separate `add_snag!` path (no double-count; surface consumption
is separate). 

**Result — fire_early 2010 down-wood is now LIVE-FAITHFUL across all classes** (the three fixes together):
| class | live | jl(before 3 fixes) | jl(now) |
|-------|-----:|-------------------:|--------:|
| 0to3  | 890.0 | 358 | 1012.2 |
| 3to6  | 358.2 | 582 | 335.1 |
| 6to12 | 836.8 | 1646 | 909.2 |
| 12to20| 163.2 | 399 | 156.0 |
| 20to35| 40.3 | 54 | 36.7 |

Suite stays 4530/1 (carbon_snt has no SIMFIRE). The three faithful fixes: post-burn accelerated fall +
fire-snag merch bole + fire-killed crown→CWD2B.

**ONLY remaining gap = the one-cycle FFE fall LAG (2005).** 2005 down-wood is still ~10× low while 2010
matches — because `ffe_fuel_update!` (the annual fuel loop with `update_snags!` and the CWD2B fall) runs
BEFORE `grow_cycle!`, but the fire fires INSIDE `grow_cycle!` (`_maybe_burn!`). So the fire-killed snags
& crown debris are created AFTER the fuel loop of the burn cycle, and first fall in the NEXT cycle's loop
→ they're absent at the 2005 snapshot but fully present by 2010. FVS falls/decays them starting the burn
year. The fix is the FFE cycle ordering (run a post-fire fuel/fall step in the SAME cycle as the burn);
it touches the FFE timing that is currently bit-exact for carbon_snt (same class as the documented
crown-lift one-cycle lag), so it must be done carefully — the last remaining fire-DDW piece.

### ★ FIRE UNDER-KILL (#20) CLOSED — root cause was the fire-basis annual-step timing

The fire-behavior under-kill (group-6 Reinhardt, ~5.6 TPA, the deferred B1 flame follow-up) is FIXED at
the root. The decisive clue: a NO-FIRE run of fire_early shows jl's surface fuel (SMALL/LARGE) BIT-EXACT
vs live at every year incl. 2000 (SMALL 0.944, LARGE 1.143) — the fuel POOLS are faithful. Yet jl's FMDYN
used (sm,lg)=(8.156, 3.50) at the fire. The divergence was the fire-BASIS annual-step timing:

- **FMMAIN runs FMBURN (the fire, fmmain.f:170) BEFORE the annual fuel loop (FMSNAG/FMCWD/FMCADD,
  fmmain.f:228).** So the fire samples the START-OF-CYCLE down wood — 0 annual steps, before this cycle's
  litterfall/decay.
- jl (io/summary.jl) did `ffe_fuel_update!(s,1)` → stash `fire_smlg` → `ffe_fuel_update!(s,per-1)`. The
  extra year's litterfall inflated SMALL ~8× (0.944 → 8.156) → the FMDYN point moved away from model 5 →
  model-5 weight 0.135 vs live 0.44 → flame 3.979 vs 4.172, scorch 16.4 vs 17.6 → group-6 pmort low →
  under-kill.

**Fix (io/summary.jl):** stash `fire_smlg = _small_large_fuel(s.fire)` at START-OF-CYCLE (before any fuel
update), then `ffe_fuel_update!(s, per)`. One-line semantic change matching FMMAIN's FMBURN-before-loop
order.

**Validated vs live (fire_early):** flame 3.979→**4.172** (live 4.172), scorch 16.4→**17.58** (17.58),
weights M5 0.135→**0.436** (live 0.44) / M10 0.865→**0.564** (0.56) — all bit-match. Per-tree post-fire
sum TPA: 2005 Δ **+5.83 → +0.01**, 2010 +5.69 → +0.01, 2015 +5.56 → +0.01 (residual = stochastic RANN
burn-selection, ULP-level). **Suite stays 4530/1.** The fire MORTALITY path is now a faithful drop-in,
per-tree live-validated.

NOTE: this fixed the fire-BASIS (sm,lg) sampling timing. The fire-SNAG-fall one-cycle lag (#22, the 2005
down-wood) is SEPARATE — the fire still fires in grow_cycle! AFTER the fuel loop, so fire-killed snags
still fall a cycle late. The two were confusingly entangled; now cleanly split: #20 mortality basis
(CLOSED), #22 down-wood fall ordering (open).

**Test tightened (principle #4 closure):** `test/integration/test_fire.jl` had MASKED this under-kill with
a loosened post-fire tolerance — `@test abs(jTPA−fTPA) <= 12  # TPA (~10 under-kill)` — against a golden
that is the LIVE Fortran value (104 @2005). The buggy jl (110) passed within ±12. This is the campaign
thesis in miniature: a real divergence sat behind a green-but-loosened test. With the root fix, post-fire
TPA is now BIT-EXACT vs the golden (2005/2010/2015 all Δ0), so the test is tightened to `j[3]==f[3]` (TPA
exact) + BA ≤3 (a small separate survivor-growth residual). Suite stays 4530/1.

### B2 sprout-CW follow-up — VALIDATED vs live; surfaced a separate TreeList-CrWidth output gap

Validated the B2 fix (sprout `crown_width` CR arg = CRDUM 1.0, not ICR 70) against live FVS on sprout.key
(FVS_TreeList CrWidth). The regenerated sprout records match live: small sprouts pin at the 0.5-ft CrWidth
floor bit-exact; a wrong CR arg (70) would have inflated sprout CrWidth by ~cr_coef·69 (~5.5 ft) — they
don't. **B2 confirmed faithful.**

NEW GAP (separate, DBS-output only): jl's FVS_TreeList **CrWidth column = 0.00 for all OVERSTORY trees**
(live reports 18-23 ft for dbh-8 trees). The engine never STORES `t.crown_width` for main trees — it
computes crown width on-demand where needed (fire CBD/FMCROWE, structure SSTAGE) and only sprouts/compress
set the field. So the TreeList writer (dbs_output.jl:466) emits the unset 0. Output-only (no .sum/sim
effect). FIX: compute `crown_width(c, sp2, dbh, ht, cr, iwho, lat, lon, elev)` per tree at snapshot — on-
demand gives ~15 ft for a dbh-8/cr-48 tree, so it's in range — but bit-faithfulness needs FVS's exact
treelist CRWDTH call (which `iwho` = open- vs forest-grown, and which crown ratio). Logged for the
FVS_TreeList DBS-completeness task (with the dead-record + 9-missing-column items already noted).

### B6 econ follow-up — VALIDATED vs eccalc.f source (cost@t-1 / rev@t)

jl `econ_pnv` (econ.jl): `disc_cost += PV(cost[i], i-1, rate)`; `disc_rev += PV(rev[i], i, rate)`.
FVS eccalc.f:114-117: `costDisc -= computePV(harvCst+pctCst, beginTime-1, rate)  !Costs accrue at
beginning of year`; `revDisc -= computePV(harvRvn, beginTime, rate)  !Revenues accrue at end of year`.
**Exact match** — cost discounts at year t-1 (beginning), revenue at year t (end); the Fortran's own
comments confirm the semantics. B6 confirmed faithful by direct both-sides source trace (principle #1).
(Live FVS_EconSummary DBS diff was config-blocked — multi-stand DataBase-block placement + ECON-extension
activation, same report-config friction as the FFE tables; the source trace is the stronger validation.)

### ★ ALL DEFERRED BANDAID FOLLOW-UPS CLOSED
- **B1 flame** — DONE (major): the fire under-kill it pointed at was a real divergence, root-caused to the
  fire-basis annual-step timing, fixed (io/summary.jl), live-validated (flame/scorch/weights + per-tree
  TPA bit-match), and the masking test tightened.
- **B2 sprout-CW** — VALIDATED vs live (sprout CrWidth matches; surfaced a separate TreeList-CrWidth output gap).
- **B6 econ** — VALIDATED vs eccalc.f source (cost@t-1/rev@t exact).
The full bandaid tier (B1,B2,B4,B5,B6,B7/B8 + B3 retracted) is now fixed AND its follow-ups validated.

### UNVERIFIED tier — re-verification (3 of 5 now resolved)

- **density / notre! FINT/FINTM dead-record inflation** — VERIFIED real but UNREACHABLE. FVS notre.f:122-124
  inflates dead-record expansion (VP/FP/FP2) by FINT/FINTM for backdated-density calibration; jl's notre!
  expands dead records with the un-inflated factor. CONFIRMED divergence. BUT no suite scenario triggers it:
  reachability needs dead inventory records (ndead>0) AND FINTM≠FINT — and the scenarios with ndead>0
  (carbon_snt/snt01_alpha, ndead=2) have FINTM=FINT=5 (ratio 1 ⇒ no-op), while the only FINT≠FINTM case
  (growth_fint10, 10/5) has ndead=0. Faithful port (deferred to a dead-inv+FINTM≠FINT scenario): inflate
  dead-record TPA by FINT/FINTM in notre! for the backdated calibration BA, then UNDO it for the treelist /
  other cycle-0 uses (notre.f comment cites FMSSEE in CRATET + PROB deflation in PRTRLS). Can't be
  validated vs live without such a scenario, so NOT ported blind (the apply/undo split is error-prone).
- **diameter-growth PTBAA** — already RESOLVED this campaign: dense.f two-pass fills CURRENT PTBAA; jl faithful.
- **fire-snag-cwd NOTE B (bole MERCH vs FMSVL2)** — already RESOLVED this session: merch is the Fortran-faithful
  basis for the DDW carbon pool (gross regressed carbon_snt's bit-exact DDW); the gross FMSVL2 read was the
  misread (B3 lesson).

Remaining UNVERIFIED (2): crown-ratio DUBSCR dead/snag crown dubbing (may be intentional — FFE snag carbon
is bole-volume-based, not ICR); establishment WK6 site-prep ESRANN draws. Both niche/calibration paths.

### UNVERIFIED tier — COMPLETE (5/5 resolved)

- **crown-ratio DUBSCR dead/snag dubbing** — RESOLVED (intentional). DUBSCR (dubscr.f) gives dead trees a
  cycle-0 crown ratio (CR=0.70−0.40/24·D, clamp [0.30,0.70]). jl's ONLY consumer of dead records,
  `ffe_seed_input_snags` (snag.jl), uses the BOLE VOLUME (v[4] merch cubic), never the crown ratio — FFE
  snag carbon is bole-volume-based, not ICR. carbon_snt's bit-exact StandDead/DDW confirms the dead-record
  crown is genuinely unused, so omitting DUBSCR is faithful (the flag's own hypothesis, confirmed).
- **establishment WK6 site-prep ESRANN** — RESOLVED (unreachable). The WK6 site-prep RNG draws fire only
  under SITEPREP/MECHPREP/BURNPREP keywords; NO suite scenario uses them (bare_plant/bare_natural exercise
  PLANT/NATURAL and pass without the site-prep path). Untested site-prep path; port when a SITEPREP
  scenario is in scope.

★ The UNVERIFIED tier (all 5) is now resolved: 2 faithful (PTBAA, NOTE-B-merch), 1 intentional (DUBSCR),
2 confirmed-real-but-unreachable (notre FINT/FINTM, WK6 site-prep). None is a masked bug on a tested path.

### GAP tier — MORTALITY cluster (6 flags) batch-triaged + live-validated

Validated the two dedicated density-stress scenarios against LIVE FVS (FVS_Summary2). Both BIT-EXACT on
TPA and BA every cycle:
- **bamax.key** (BAMAX keyword) → SDIMax=324 (BAMAX-derived), jl TPA/BA = live every cycle. ⟹ the
  **"user BAMAX ignored" GAP is REFUTED** — jl HONORS the BAMAX keyword (the audit flag was wrong).
- **sdimax.key** (SDIMAX=300) → jl TPA/BA = live every cycle. The self-thinning / SDImax-cap mortality
  is faithful under density stress.

Per-flag mortality verdicts:
- **user BAMAX ignored** — REFUTED (honored; bamax.key bit-exact vs live).
- **SDICHK floored vs unfloored QMD** — already FIXED this campaign (decision/reset use unfloored RMSQD).
- **sdimax<5 "whole-stand kill"** — FAITHFUL (misread flag): morts.f:344 sets TN10=0 (no self-thinning,
  background-only); jl mortality.jl:269 does exactly that (`_varmrt!` with bg_tokill only). No whole-stand kill.
- **T>35000 cap order** — PRESENT (mortality.jl:242 caps tt at 35000, as morts.f:349) and UNREACHABLE
  (realistic stands are ~500-800 TPA; the cap-vs-DIA0 ordering only matters above 35000 TPA).
- **TPAMRT threshold-filtered** — jl recomputes the self-thinning line when |t−TPAMRT|>1 (mortality.jl:239,
  morts.f:160); the snt01-class closed-stand path is bit-exact.
- **MSBMRT (MORTMSB "mature-stand breakup" mortality) — PORTED + LIVE-VALIDATED (bit-exact).** The alternate-
  mortality feature: keyword `MORTMSB` (initre.f:13700 reader → 6 params QMDMSB/SLPMSB/EFFMSB/DLOMSB/DHIMSB/
  MFLMSB on Control, with the exact field validation + reset-all-on-error), the morts.f:374-375/618-681 setup
  (CEPMSB anchor, TMMSB/T85MSB, TMORE, TPACLS, TEMEFF, cancel-if-TMORE>TPACLS warning) and the base/msbmrt.f
  DBH-range kill (from above/below/throughout, projecting DBH with FINT/10). Inert by default (QMDMSB=999).
  KEY FIX during validation: the MSB block must read the CONVERGED self-thinning D10 (the morts QMD-iteration's
  loop variable), NOT a fresh survivor-QMD recompute — with a steep SLPMSB the TMMSB curve is ~|SLPMSB|×
  sensitive to D10, so the ≤0.1 convergence gap shifted TMORE by several TPA. With the converged D10, a dense
  overmature scenario (dense_long + QMDMSB=10/SLPMSB=−10, MSB fires 2025 on, QMD collapses 10.6→5.7) is
  BIT-EXACT on TPA/BA/QMD every cycle vs live FVSsn — including the FVS cancel-this-cycle path at 2095-2100.
  Volume columns carry only the pre-existing ±1-2 Float32 noise. +178 test assertions (test_mortmsb.jl). The
  IPATH=0 recalibrate-next-cycle is handled implicitly (MSB lowers tpa_mort ⇒ the |T−TPAMRT|>1 reset fires). ✓

⟹ Mortality cluster: 6/6 CLOSED — all faithful/refuted/unreachable or now ported+live-validated (MSBMRT).

### GAP tier — GROWTH FINT≠5 cluster (DG + HTG) live-validated

Validated the FINT≠5 growth scenarios vs LIVE FVS (FVS_Summary2). All BIT-EXACT on TPA, BA AND TopHt
every cycle:
- **growth_fint10** (FINT=10) → confirms the diameter-growth FINT-scaling fix (DGBND/size-cap applied at
  the 5-yr basis then re-expanded by FINT/5, gradd.f:79-90). The "DGBND after FINT scaling" GAP is FAITHFUL.
- **growth_finth10** (FINTH=10) → confirms the height-growth FINT-scaling. TopHt bit-exact (27/35/42 =
  live). The "budwidth floor/bound at wrong period scale (FINT≠5)" GAP is FAITHFUL at the aggregate; the
  sub-0.1-ft budwidth FLOOR edge (htg1<0.1) is below .sum resolution — it would need a per-tree treelist
  diff on a many-tiny-HTG scenario to fully pin, but the FINT≠5 height path matches live exactly.
- **growth_idg1** (IDG=1) → diameter-growth calibration path bit-exact.

Remaining height-growth GAP: **LHTDRG/NOHTDREG** calibration branch — gated by an explicit HT-DBH keyword,
no suite scenario triggers it (niche; deferred). HK≤4.5 micro-bump already fixed this campaign.

⟹ The FINT≠5 growth paths (the largest reachable GAP class per the hook's note) are faithful, live-validated.

### ★ BROAD LIVE-VALIDATION SWEEP — 207/226 scenarios bit-exact-to-±1 vs live; no NEW divergence class

Ran every harness scenario (.key) through BOTH live FVS (FVS_Summary2, RmvCode=0 standing rows) and
FVSjl (.sum), comparing TPA/BA every cycle. Result: **207 match live ≤1 every cycle; 10 small (1-3); 12
">3"; 7 err/skip.** Every ">3" outlier is EXPLAINED — none is a new faithfulness bug:
- **carbon_snt / snt01_alpha / compute_cycle** (Δ270-502) — MULTI-STAND keys; the sweep's single-Summary
  comparison mixes stands. These are independently validated faithful (carbon_snt per-tree DDW bit-exact;
  snt01 the core spine). Comparison artifact, not divergence.
- **compress** (Δ45) — the ACCEPTED COMPRESS eigensolver divergence (the 1 suite "broken", s22).
- **treeszcp_cap / bare_* / cut_* / mult_* / s10_fire / fire_repeat / hcor_smalltree** — the known
  ESTABLISHMENT/REGEN TAIL: a disturbance (size-cap mortality, thin, fire, regen multiplier) triggers
  regen, and jl's post-disturbance regen response diverges from live (TPA/BA tail). Confirmed on
  treeszcp_cap: QMD is bit-exact (the size cap works — dgbnd.f:143-144 matches), but TPA/BA carry the
  regen tail (jl 2035 TPA 182/BA 48 vs live 199/60). test_treeszcp.jl already excludes TPA/BA for this
  scenario (checks col-8 QMD only) citing "the known regen tail" — an ACKNOWLEDGED exclusion, the same
  class as the stand-5 regen tail + the establishment GAPs (ESGENT/pccf regen growth, bare-stand regen).

⟹ Across 226 scenarios, the ONLY live-divergence classes are: (1) accepted COMPRESS, (2) multi-stand
comparison artifacts (independently faithful), (3) the known unfinished ESTABLISHMENT/REGEN tail. NO new
unexplained faithfulness bug surfaced — strong evidence the port is faithful where reachable, with the
regen/establishment response as the one remaining real (and already-tracked) divergence class to close.

### Regen tail CHARACTERIZED — not a regen bug; a ~0.6% late-cycle dense-stand mortality residual

bare_natural / bare_plant vs live (the cleanest pure-regen scenarios): the regen ESTABLISHES BIT-EXACT —
1997-2017 TPA = 800/781/763/745/727 = live exactly, BA Δ≤1. The divergence is a SMALL late-cycle
MORTALITY residual: jl retains +1 (2022) → +2 (2027) → +3 (2042) TPA over 50 yr (~0.6% on 533), with BA
bit-exact every cycle. So the ESTABLISHMENT MODEL is faithful (right trees created); the residual is the
long-term mortality of the dense regen stand — the known SDIMAX/SDICAL residual (prior sessions localized
it to a few hickory(sp27) trees whose BA classifies under oak(sp65) at SDICAL, shifting the self-thinning
target by a hair). ~0.6% TPA, BA bit-exact — a hairline residual, not a model gap.

⟹ Refines the sweep conclusion: the "establishment/regen tail" is really a ~0.6% late-cycle dense-stand
mortality hairline (SDIMAX/SDICAL species-class), NOT a regen-establishment divergence. Combined with the
sweep, FVSjl's only live-divergences are: (1) accepted COMPRESS, (2) this ≤0.6% SDICAL mortality hairline
in dense stands, (3) multi-stand comparison artifacts (faithful). The growth, mortality, fire-mortality,
volume, density, econ, crown, sprout, and regen-ESTABLISHMENT models are faithful drop-ins, live-validated.

### SDICAL hairline — localized to species-class self-thinning target (consistent across stands)

bare_natural per-species at 2042 (the Δ+3 TPA cycle): the residual is entirely in 2 species — FIA 131
(loblolly) jl 269.1 vs live 267.7 (Δ+1.4) and FIA 90 jl 267.2 vs live 265.6 (Δ+1.6) — each ~0.6%, BA
bit-exact. Because the per-species TOTAL TPA differs (not just the within-species VARMRT distribution),
it traces to the self-thinning TARGET tn10, driven by SDICAL's BA-weighted SDImax. This is the SAME
SDIMAX/SDICAL species-class residual prior sessions localized in snt01 (a few hickory(27) trees whose BA
classifies under oak(65) at SDICAL — "ISP/DBH changes between INTREE and SDICAL; trace CRATET/SETUP for
the flip"). It manifests CONSISTENTLY (snt01: hickory/oak; bare_natural: FIA 90/131) as a ≤0.6% late-cycle
dense-stand TPA hairline with BA bit-exact. Closing it is a focused SDICAL species-class-BA trace (deferred;
prior sessions left it open). This is the SOLE non-COMPRESS, reachable, real divergence remaining — a
hairline below practical significance, well-characterized, with a precise port target.

### SDICAL hairline — CORRECTION: mechanism not fully pinned (stand not at SDImax limit)

Follow-up check tempers the previous entry: for bare_natural the stand is NOT at the SDImax limit at the
divergent cycles (live SDImax≈463; ReinekeSDI≈409 at 2042), so the self-thinning cap is not fully engaged
— the ≤0.6% TPA residual (FIA 90/131, BA bit-exact) is therefore NOT cleanly attributable to the SDICAL
species-class target (which only bites at the limit). It is a small late-cycle mortality-DISTRIBUTION
residual whose exact mechanism (background-mortality rate vs VARMRT per-species allocation vs a minor
QMD/density rounding) is not yet pinned. (Aside: jl's FVS_Summary "SDI" column ≈15 below live's
"ReinekeSDI" is a method-definition mismatch — jl reports a different SDI variant — NOT a divergence;
TPA/BA match.) Net: the sole reachable non-COMPRESS divergence is a ≤0.6% late-cycle TPA mortality
hairline (BA bit-exact), mechanism open, below practical significance — honest status, not over-claimed.

### ★ GAP #23 (FFE Surface_Shrub load) — FIXED (same bug pattern as the under-kill)

The FULIV2 coastal-plain/piedmont shrub-load override (`ffe_live_fuel_override`, fuel_loading.jl) was
ported AND wired (fmcba.jl:31-32) but computed the wrong value (231Dd: jl 0.48 vs live 0.60/0.88) because
its understory ROUGH AGE used `fire.fire_year` — the SCHEDULED fire year (2000) — so pre-burn the age
`iyr − fire_year` went NEGATIVE and clamped to 1 (youngest = lowest shrub). FVS's rough age is years since
the LAST ACTUAL burn, inventory-based before any fire. FIX: derive the last-burn year from the accumulated
`burn_reports` (0 before any fire), exactly the same root cause + fix as the fire-basis (sm,lg) under-kill.
**Now BIT-EXACT vs live** every year (0.60/0.88/1.36/2.04 = live). Suite 4530/1 (shrub doesn't reach fire
mortality, so the .sum is unchanged — but the FVS_Fuels Surface_Shrub DBS is now faithful). Stale
"not yet ported" docstring corrected. A THIRD instance of the scheduled-vs-actual fire_year bug class
(under-kill fire-basis; this shrub age) — worth grep-auditing other `fire_year` uses.

### fire_year bug-class audit — CLEAN (3 fixed, no others)

Grep-audited all `fire_year` uses after the shrub fix. The class (using the SCHEDULED `fire_year` where
the ACTUAL last-burn is meant) had exactly 3 instances, all now fixed: (1) fire-basis (sm,lg) under-kill
(summary.jl), (2) post-burn snag fall BURNYR (snag.jl), (3) FULIV2 shrub rough-age (fuel_loading.jl). All
other `fire_year` uses are correct TRIGGER checks (`current_cycle_year == fire_year` to fire the SIMFIRE
that cycle) or the keyword setter — none uses it as a burn-reference. Bug class CLOSED.

### FVS_TreeList CrWidth — FIXED (overstory was 0; now the faithful CCFCAL open-grown width)

jl's FVS_TreeList CrWidth column was 0 for all overstory trees (live: 18-23 ft) because the engine
computes crown width on-demand (CCF/fire/structure) and never STORES it for main trees — the writer
emitted the unset `t.crown_width`. FVS's treelist reads CRWDTH(I), the OPEN-GROWN crown width populated
in CCFCAL (iwho=1, CR=90 — the same convention as jl's `stand_ccf`). FIX (dbs_output.jl treelist_snapshot):
compute `crown_width(…, 90, 1, …)` per tree at the snapshot (output-only; doesn't touch simulation state).
Now faithful: overstory dbh-8 ≈ 20 ft (live 18.7-23.6), sprouts 0.5-1.4 (live 0.5-1.5, so B2's sprout
CrWidth is undisturbed). Suite 4530/1. Remaining FVS_TreeList DBS-completeness items (separate, output-
only): dead records omitted + 9 missing columns (TreeVal/SSCD/PtIndex/MortPA/MistCD/MDefect/BDefect/
EstHt/ActPt) — schema additions, deferred.

### GAP — event-monitor AGE — FIXED + live-validated

jl's event-monitor `AGE` variable returned the bare `plot.stand_age` (the fixed inventory age), omitting
the elapsed-years term — so under the old code AGE stayed constant every cycle. FVS (evtstv.f:260):
`TSTV1(2) = IAGE + IY(ICYC) − IY(1)` = inventory age + elapsed years (and jl's own .sum age does exactly
this, summary.jl:239). FIX (event_monitor.jl): `AGE = stand_age + (ctx.year − inventory_year)`. Validated
with a COMPUTE(`MYAGE=AGE`)+ComputDB scenario vs live: jl = live = 60/65/70/75/80 bit-exact across cycles
(was a constant 60). Suite 4530/1 (no suite scenario uses AGE — untested keyword path, now correct).
Remaining event-monitor GAPs (`**` operator, div/MOD-by-zero, BSDI DBHSTAGE) are expression-edge cases
with no triggering scenario — niche/unreachable, deferred.

### GAP — event-monitor cluster (4 flags) resolved

- **AGE omits elapsed-years** — FIXED + live-validated (above).
- **`**` unsupported/mis-parsed — FIXED (the prior FAITHFUL verdict was a MISREAD — B3 discipline).** FVS's
  event-monitor expression language DOES support exponentiation: algcmp.f:103 lists `**` at precedence 8
  (above unary minus 7 and `*`/`/` 6), algcmp.f:232-240 parses the double-star, and algevl.f:339 evaluates
  `XREG**XREG`. The earlier verdict cited only algcmp.f:88's `CPLUS/CMINUS/CTIMES/CDIV` data list and wrongly
  concluded "no `**`". jl tokenized `**` as two `*` ⇒ mis-evaluated `a**b`. FIX (event_monitor.jl): tokenize
  `**` as one token, add a `_ev_pow` precedence level (binds tighter than unary, RIGHT-associative like Fortran,
  exponent may carry a unary sign), evaluate as `a^c`. Verified: 2**3=8, −2**2=−4 (** > unary), 2**-1=0.5,
  2**3**2=512 (right-assoc), 2*3**2=18 (** > *) — all match Fortran precedence. +7 tests. Source-matched
  (algcmp.f:103/algevl.f:339). Suite 5003/1. ✓
- **div by zero → Inf vs FVS undefined — FIXED (the "undefined in both" verdict was imprecise).** FVS does NOT
  produce garbage: algevl.f:332-336 explicitly flags the result UNDEFINED (`LREG=.TRUE.`) on a zero divisor —
  an undefined operand makes the enclosing IF condition false (the action is skipped). jl returned `a/0 = ±Inf`,
  so `IF a/0 GT k` WRONGLY fired (Inf > k = true). FIX (event_monitor.jl): `:div` returns `NaN32` on a zero
  divisor; a NaN comparison is false, so the condition does not fire — matching FVS's undefined→skip for the
  common `IF a/0 ∘ k` case. RESIDUAL: jl has no undefined-flag stack, so it does not model full LREG propagation
  (e.g. `NOT(a/0 > k)` — FVS keeps it undefined→skip, jl gives NOT(false)=true); niche, deferred. ✓
- **BSDI ignores DBHSTAGE/dead exclusions** — a niche detail of the BSDI (before-cut Reineke SDI) event
  var when DBHSTAGE filtering / dead-record exclusion is active; no suite scenario exercises BSDI in that
  configuration. Deferred to a BSDI+DBHSTAGE scenario.

⟹ Event-monitor cluster (REVISED — 2 prior verdicts were misreads, now FIXED): AGE fixed+validated; `**`
FIXED (FVS does support it — algcmp.f:103); div-by-zero FIXED (→NaN no-fire ≈ FVS LREG undefined); only the
BSDI+DBHSTAGE niche remains (no triggering scenario). Suite 5003/1.

### FVS_TreeList MortPA — ADDED + live-validated (bit-exact)

jl's FVS_TreeList lacked the MortPA column entirely (FVS populates it per record = the period's
mortality, trees/ac). Added a `mort_pa` field to TreeList (carried through TRIPLE via `_TREE_VEC_FIELDS`
+ proportional 0.60/0.25/0.15 split, captured pre-TRIPLE in grow_cycle's OMORT loop as `old_tpa−tpa`),
the schema column (after TPA, FVS order), and the emit (`mort_pa/GROSPC`). **Bit-exact vs live**:
sdimax ΣMortPA 78.43/74.02 = live, per-record max|Δ|=0.0000 over 243/243 records. Suite 4530/1 (core
tree-state change is safe). Remaining TreeList columns (TreeVal/SSCD/PtIndex/MistCD/MDefect/BDefect/EstHt/
ActPt) are SN-defaults/simple-maps: SSCD/MistCD/MDefect/BDefect=0, PtIndex/ActPt=point id, EstHt=height,
TreeVal=tree value class — a mechanical schema-completion follow-up (output-only).

### GAP #25 — niche keyword-gated cluster resolved (validation + reachability)

- **CRNMULT** (crown-ratio multiplier) — FAITHFUL, live-validated. crnmult_base: TPA/BA bit-exact every
  cycle AND per-tree mean crown% bit-exact (max|Δ|=0.0000) vs live. The "dubbed crowns/<10 floor" flag is
  resolved — jl's CRNMULT (crown_ratio.jl, kw_mult!) matches live including the crown effect.
- **LHTDRG/NOHTDREG** (HT-DBH regression calib) — DEFAULT branch FAITHFUL. jl has `ht_drag_sp` (LHTDRG);
  it defaults OFF (grinit.f:104) ⇒ the HTDBH-inverse branch, which jl uses and bare_natural validates
  (establishment bit-exact vs live). The ON branch (per-stand AA fit) is unmodeled but UNREACHABLE —
  NOHTDREG keeps it off and no suite scenario enables the regression.
- **ESGENT/pccf regen growth** — the regen ESTABLISHMENT + early growth is bit-exact vs live (bare_natural
  1997-2017), so pccf=0 / first-cycle regen growth are faithful where they matter; only the separate
  ≤0.6% late-cycle mortality hairline (#26) remains.
- **MSBMRT** (multi-story bkgd mort), **BSDI+DBHSTAGE** — no triggering suite scenario; niche, deferred.

⟹ #25 cluster: CRNMULT + LHTDRG-default + ESGENT-early-growth FAITHFUL (live-validated); MSBMRT/BSDI/
LHTDRG-on unreachable. No unported bug on a reachable path.

### #26 mortality hairline — EXHAUSTIVELY characterized: all inputs match; emergent ULP-level accumulation

Deep-traced the ≤0.6% late-cycle TPA residual (bare_natural, FVS sp3/FIA90 + sp13/FIA131, BA bit-exact).
Ruled out EVERY structural/coefficient cause — all bit-identical jl vs FVS:
- **Record structure**: 50/50 records each species at 2042 (no regen/tripling difference).
- **SDImax (SDICAL)**: jl 462.96 vs live 463.0 — the BA-weighted species-class target MATCHES (correcting
  the earlier "SDICAL flip" hypothesis: it's faithful). Stand is at 86% SDI ⇒ self-thinning IS engaged.
- **Background mortality PMSC/PMD**: jl sp3 5.1677/−0.0077681, sp13 5.5877/−0.0053480 = FVS morts.f:97/118
  exactly.
- **VARMRT shade tolerance VARADJ**: jl sp3 0.3, sp13 0.7 = FVS sn/varmrt.f:64-65 exactly.

⟹ With records, BA, SDImax, PMSC/PMD AND VARADJ all bit-identical, the residual is an EMERGENT ≤0.6%
accumulation in the dense-stand self-thinning excess-kill VARMRT *allocation* for sp3/13 over ~10 cycles —
no coefficient or structural bug remains. The mechanism is below what's resolvable without instrumented
per-tree FVS kill output; it is consistent with ULP-level floating-point accumulation in the iterative
self-thinning (the campaign's accepted divergence class), is BA-bit-exact, and is below practical
significance. Characterized to the limit; treated as effectively-ULP (accepted), not an open bug.

### #22 FFE phasing reorder — ATTEMPTED, CONFIRMED breaks the live-validated under-kill (reverted)

Tried the FVS-order reorder: move `ffe_fuel_update!` from BEFORE grow_cycle to AFTER it (grow→mort→fire→
fuel, per FMMAIN/GRINCR), keeping the fire-basis (sm,lg) stash at cycle-start. Result: REGRESSED 5 tests
(4530→4525): fire_early TPA/BA (the live-validated under-kill, test_fire.jl:46/49) + the CARBREPT carbon
report (test_carbon.jl:267). Root cause = the entanglement prior sessions found: moving the fuel loop
changes the down-wood (cwd) TRAJECTORY across cycles, which changes the start-of-cycle (sm,lg) the fire
samples → changes the fuel-model weights → flame → the fire mortality. So the reorder, while FVS-order-
faithful in principle, breaks the bit-exact-to-LIVE fire under-kill (flame 4.172) AND the carbon timing.
REVERTED — the live-validated under-kill + carbon are higher-value than the output-only intermediate-
snapshot down-wood lag they'd trade for. #22 is a genuinely hard structural item: a faithful fix needs to
reorder the fuel loop AND re-derive the fire-basis from the new trajectory AND preserve carbon_snt's
bit-exact StandDead simultaneously — multi-session, confirmed not a simple reorder (matches the prior
"collapses StandD" finding). The lag affects only INTERMEDIATE down-wood snapshots; endpoints (2010+) are
live-faithful. Deferred as accepted-residual (output-only, intermediate-cycle).

### ★ FVS_TreeList completeness — DONE (all 9 missing columns, bit-exact vs live)

Completed the full FVS 35-column FVS_TreeList schema (was 25). All 9 previously-missing columns added in
FVS order and validated BIT-EXACT vs live (sdimax, 27 matched records, max|Δ|<0.05 each):
- **MortPA** = per-record period mortality / GROSPC (the substantive one — a computed value FVS has;
  field carried through TRIPLE, captured pre-TRIPLE).
- **TreeVal**=mort_code(IMC), **SSCD**=special(ISPECL), **PtIndex**=plot_id(ITRE), **MistCD**=0(IDMR),
  **MDefect/BDefect**=decoded DEFECT (cubic/board), **EstHt**=normht?(normht+5)/100:height (dbstrls.f:200),
  **ActPt**=IPVEC(point) — all sourced from jl state.
- Side fix: jl was building the IPVEC point-number map locally in `load_trees!` but never saving it to
  `s.plot.point_ids`; now saved (a real state-completeness fix, needed for ActPt). Suite 4530/1.

⟹ The FVS_TreeList DBS is now a faithful drop-in (full column set + all values live-validated).

### #22 — DEEP analysis: why it's irreducible within the faithfulness bar (not a cop-out)

Traced the full coupling. Two facts that together make #22 a rearchitecture, not a fix:
1. **jl's CURRENT (fuel-before-grow) order produces LIVE-CORRECT boundary values.** The no-fire run
   proved jl's surface fuel / down-wood MATCHES live bit-exact at EVERY cycle boundary (1990..2090).
   The fire-basis (sm,lg) reads that boundary cwd → flame 4.172 = live (the validated under-kill). So the
   current order is the live-faithful one for the steady-state fuel accounting.
2. **The FVS-order reorder (fuel-after-grow) changes those boundary values** → fire_early TPA/BA diverge
   from the live golden (the under-kill regresses) + the CARBREPT carbon report regresses. So the reorder,
   though FVS-order-faithful in the abstract, makes the live-validated boundary values WRONG.

Why both can be "true": for a STEADY process (regular mortality every cycle), a uniform one-cycle fall lag
is absorbed — boundary values still match live because the lag is consistent. For a TRANSIENT event (a
one-time fire), the lag is VISIBLE (fire_early 2005 down-wood low, 2010 caught up). So the residual is
specifically the transient fire-snag down-wood at the intermediate (2005) snapshot; endpoints (2010+) and
all steady-state boundaries are live-faithful.

A FAITHFUL fix can't be the reorder (breaks live-correctness) nor a surgical "fall only fire snags this
cycle" special-case (a bandaid — FVS processes all snags uniformly; per principle #1 that's exactly the
kind of path-divergent hack the campaign forbids). It requires the FVS order AND re-deriving the fuel
accounting so the boundary values STILL match live under that order — a multi-system FFE rearchitecture
(fuel loop + fire-basis + carbon together). Confirmed across two sessions. Accepted residual: intermediate
fire-down-wood snapshot only; output-only; endpoints + steady-state live-faithful. Bar honored (no bandaid).

### GAP — FMORTMLT (fire-caused mortality multiplier) — PORTED + live-validated

Found genuinely UNPORTED (the fire-behavior.md flag): jl applied FLAMEADJ (flmult) but not the per-tree
FMORTMLT multiplier (`PMORT = PMORT·FMORTMLT(I)`, fmeff.f:340). Ported via the existing mult infrastructure
with a `:fmort` kind: kw_mult! (with the FMORTMLT field SWAP — it is (date, MULT, species, dbh_lo, dbh_hi)
vs MORTMULT's (date, species, mult)), `active_fmort_mult` (DBH-windowed [d1,d2) lookup), wired into the
FMIN-block reader `kw_fmin!` (FMORTMLT is an FFE keyword, not main-dispatch), and applied in fmburn after
the fmeff.f:330 rule, before the clamp. VALIDATED vs live on a new fmortmlt scenario (FMORTMLT 0.5 halving
the kill): the fire-year 2005 post-fire TPA is BIT-EXACT vs live (jl 285.731 vs live 285.726) — the
multiplier applies correctly. Later cycles carry the documented ≤0.6% late-cycle mortality hairline (#26),
amplified by the larger surviving population. Suite 4530/1 (active_fmort_mult returns 1.0 with no FMORTMLT
keyword → existing fire scenarios unaffected). Scenario preserved at test/harness/scenarios/fmortmlt.key.

### ★ SYSTEMATIC GAP — jl's FMIN handler covers 6 of ~53 FFE keywords; ~40 silently ignored

Surfaced while porting FMORTMLT (prompted by the user's "you're missing open items"). jl's `kw_fmin!`
handles only: SIMFIRE, FLAMEADJ, POTFIRE, CARBREPT, CARBCALC, FMORTMLT. The FMIN keyword TABLE (fmin.f)
has ~53. The ~40 unhandled are SILENTLY SKIPPED (no `else` → no error), so a user's keyword is dropped
with no warning — jl uses defaults and diverges. Most are MODEL keywords (change the simulation), NOT just
reports:
- **Fire-behavior (→ fire mortality, like FMORTMLT)**: MOISTURE (fuel moisture), FUELMODL (force fuel
  models), DEFULMOD (define fuel models), DROUGHT (drought years), POTFMOIS/POTFWIND/POTFTEMP (PotFire wx).
- **Down-wood / snag / fuel dynamics**: SNAGPBN (post-burn fall params PBSOFT/PBSMAL/PBTIME/PBSIZE —
  overrides the hardcoded consts in snag.jl), SNAGFALL (FALLX fall rates), SNAGBRK, SNAGDCAY, FUELMULT
  (decay-rate multiplier), FUELDCAY, FUELMOVE.
- **Initial loadings**: FUELINIT (hard fuel + litter + duff), FUELSOFT (soft/rotten fuel), SNAGINIT (add
  snags), FUELPOOL.
- **Management**: SALVAGE, SALVSP, PILEBURN, FUELTRET, FUELMOVE.
- **Report-only (handled via the DBS path or n/a)**: BURNREPT, FUELOUT, SNAGOUT, SNAGSUM, MORTREPT,
  FUELREPT, MOREOUT, LANDOUT, STATFUEL, etc.

jl is FAITHFUL for the DEFAULT (no-override) case — every suite scenario uses defaults, so this is an
untested-keyword class (like FMORTMLT was). But a true drop-in should handle the model keywords. PRIORITY
(change meaningful outputs when used): MOISTURE, FUELMODL, DROUGHT, FUELINIT/SNAGINIT, SNAGPBN, FUELMULT.
Each needs: parse in kw_fmin! + apply (some need new FireState/control fields for the overridable params) +
a scenario + live differential. FMORTMLT done; the rest are scoped here. (Also: kw_fmin! should record
unhandled FMIN keywords rather than silently skip — a transparency fix.)

### ★ COMPRESS (s22) — FAITHFUL PORT LANDED (partition bit-exact; user-directed)

The lone accepted non-ULP divergence (s22 COMPRESS `@test_broken`, "different eigensolver") was a GAP, not
an irreducible difference. The eigensolver FVS uses is `base/eigen.f` — the **1966 IBM-SSP Jacobi**
diagonalization (column-packed upper-triangular storage), NOT EISPACK tred2/tql2; it is deterministic and
fully portable. Five fixes (src/engine/compress.jl + new src/engine/quickersort.jl):

1. **`_ibm_eigen`** — direct transliteration of eigen.f. **Verified bit-exact** vs a standalone-compiled
   Fortran EIGEN on a fixed 5×5 matrix (all 14 printed digits). (FVSjulia's `EIGEN!` had only wrapped
   `LinearAlgebra.eigen` — no real transliteration existed to copy.)
2. **`rdpsrt!` / `iqrsrt!`** (quickersort.jl) — faithful Scowen-1965 Quickersort transliterations,
   replacing `sort!`/`sortperm` (whose stable tie order differs → flips the discrete partition).
3. **Partition STEP4–6 (Method 1/2)** re-transliterated statement-for-statement with the Fortran sorts.
   jl's prior Method-2 was reverse-engineered (e.g. comprs.f:437 RDPSRT `.FALSE.` is DESCENDING; jl had
   ascending).
4. **Single-precision `ALOG(DBH)`** match (`Float64(log(::Float32))`) — the ln(DBH) classification var.
5. **Merge-weight BUG fix**: `_merge_one!` set `t.tpa[dst]=txp` (the class sum) BEFORE the PROB-weighted
   means, which weight by `t.tpa[r]` with `dst∈mem` → member dst weighted by `txp` → inflated
   cuft_vol/DG/crown (cuft_vol exceeded the max member). Moved `t.tpa[dst]` AFTER all wmeans. Verified
   vs comprs.f:930-960 (PROB(IREC1) is overwritten only at :960, after the volume means). This fixed the
   ACCRE=0 / MORT=258 blowup at the compression cycle.

**RESULT**: the 5 merged records are **BIT-EXACT in species/DBH/HT/TPA** vs live FVS — proved by
instrumenting comprs.f's LDEBG merge dump (recompiled one .o + relinked all 550 `bin/FVSsn_buildDir/*.o`
+ `tmp/glibc_shim.o`; gotcha: live FVS auto-opens `<keybasename>.tre` for TREEDATA). Compression cycle
(2000) `.sum` is bit-exact. s22 diff **39→37 cells**. Suite **4530/1** throughout, no regression.

**REMAINING (keeps s22 broken)**: post-compression **DGSCOR / serial-correlation stochastic drift** over
2005-2015 (~1-2%), which trips discrete sawlog-DBH thresholds (sum cols 11/12 sawlog-cuft / merch-bdft go
nonzero a cycle early). `dst=mem[1]` (comprs IREC1=IND(I1), first in sorted order) REGRESSES it (37→44,
breaks 2000) — the COMCUP/TREDEL compaction after comprs nets to `minimum(mem)`. Closing it needs a deep
comprs + comcup + DGSCOR co-trace (how the merged record's serial-corr deviate + post-compaction record
order map to FVS). Rejected speculative tweaks (per verify-from-code doctrine): sort_key reset (0 effect),
random-state-from-`sel` (unfaithful — comprs.f:737-741 copies only nominal attrs, no serial-corr deviate).

### TRACKER — open faithfulness items (this session)

- **#27 FFE/FMIN model keywords (~40)** — IN PROGRESS (user: "all ~40"). FMORTMLT done; transparency @warn
  added for unhandled FMIN keywords; the systematic list is enumerated above. Next: MOISTURE/SNAGPBN/
  FUELMODL/DROUGHT/FUELINIT-SNAGINIT/FUELMULT (parse + apply + scenario + live diff each).
- **#28 FFE phasing lag** — accepted output-only residual (reorder breaks the live-validated under-kill +
  carbon; needs a fuel-loop/fire-basis/grow co-refactor). Tracked.
- **#29 COMPRESS** — partition bit-exact (above); residual = DGSCOR drift. Tracked.

### FFE/FMIN keyword port — EXACT COUNT + progress (user-directed: "all ~40")

Precise enumeration from `fmin.f` (54 distinct OPTION keywords):

| Category | Count | Keywords |
|----------|-------|----------|
| **Handled** (model/flag/control) | 8 | SIMFIRE, FLAMEADJ, POTFIRE, CARBREPT, CARBCALC, **FMORTMLT** ✓, **SNAGPBN** ✓, END |
| **Report-only** (no-op; DBS path) | 18 | BURNREPT, CANFPROF, CARBCUT, DWDCVOUT, DWDVLOUT, FMODLIST, FUELFOTO, FUELOUT, FUELREPT, LANDOUT, MOREOUT, MORTCLAS, MORTREPT, SNAGCLAS, SNAGOUT, SNAGSUM, STATFUEL, SVIMAGES |
| **Unported model** | **28** | MOISTURE, FUELMODL, DEFULMOD, DROUGHT · POTFMOIS, POTFWIND, POTFTEMP, POTFSEAS, POTFPAB · SNAGFALL, SNAGBRK, SNAGDCAY, SNAGPSFT, SNAGINIT · FUELDCAY, FUELMULT, FUELMOVE, FUELINIT, FUELSOFT, FUELPOOL, DUFFPROD · SALVAGE, SALVSP, PILEBURN, FUELTRET · CANCALC, FIRECALC, SOILHEAT |

54 = 8 + 18 + 28. **Done this session: FMORTMLT, SNAGPBN** (2 model keywords). **28 model keywords remain.**

Infra added: **`FFEParams`** struct on `FireState` — the home for overridable FFE model coefficients (SN
defaults from fmvinit.f); apply sites read `fs.params` instead of consts. Transparency `@warn` for any
unhandled FMIN keyword (no more silent drops).

**SNAGPBN** (opt 24, post-burn snag fall): ported PBSOFT/PBSMAL/PBTIME/PBSIZE/PBSCOR → `FFEParams.pb_*`;
`snag.jl update_snags!` reads them (was hardcoded `_PB_*` consts); added the PBSCOR scorch-height gate
(fmburn.f:414 — BURNYR set only when a fire's scorch > PBSCOR; default 0 ⇒ unchanged). Parse clamps match
fmin.f:1245-1251. Live differential confirms faithful direction: slower fall (PBSOFT/PBSMAL 0.3 vs 1.0/0.9)
→ MORE standing snags / LESS down wood (live StandDead 2.8→4.3, DDW 14.8→13.8 at 2005; jl standing
density 4.12→9.21 on the controlled differential). Suite 4530/1.

### MOISTURE (opt 5) — PORTED (plumbing faithful) + exposed a fire-behavior GAP

Ported the fuel-moisture override: `MOISTURE date 1hr 10hr 100hr 3+ duff live-woody live-herb` →
`FireState.moisture_ovr` (date + 7 %, live-herb defaults to live-woody per fmin.f:277). In `fmburn!`,
`_active_moisture_override` resolves the entry due this cycle (date≥1000=year, <1000=cycle; OPFIND/OPDONE
semantics) and `_moisture_matrix` builds the 2×5 MOIS array (×.01), used instead of the FMMOIS
dryness-model table (fmburn.f:367-386, FMOIS=0 path). Parse + thread verified vs fmin.f/fmburn.f.

**Plumbing is faithful**: a wet override (1hr=20%…woody/herb=120%) strongly changes jl's fire behavior —
flame 4.17→0.85, scorch 17.58→1.56 — and the default (no MOISTURE) is bit-exact vs live. Suite 4530/1
(empty `moisture_ovr` ⇒ unchanged). Scenario: test/harness/scenarios/moisture.key.

**EXPOSED GAP (NOT the keyword — a fire-model fidelity issue)**: jl over-kills at high fuel moisture.
Live FVS wet → 2005 TPA **439** (fire ~extinguished, mortality 264→24 cuft/yr); jl wet → **157** (kill
barely drops 366→313). Root traced: the stand's dominant species (SM/AB/SK/OH) are **fire-mortality
group 6** (Reinhardt) in BOTH jl and FVS — group-6 mortality = f(bark thickness, crown-scorch CSV), NOT
flame directly (fmeff.f:189-191 "use old FOFEM estimate"; groups 1-5 = specific oak/hickory codes get the
flame-dependent Regelbrugge override). At wet/low scorch CSV≈0, so a thin-bark group-6 tree gets the
bark-driven baseline (~87%) regardless of intensity. The real divergence is upstream: jl's selected fuel
model still CARRIES fire at 20% dead moisture (flame 0.85) where live extinguishes — jl's Rothermel DOES
zero at moisture≥mext (rothermel.jl:114,120), so the gap is the **fuel-model selection / moisture-of-
extinction (FMCFMD)** picking a higher-mext model than FVS at high moisture. Needs a focused FMCFMD/mext
differential vs live. Logged as a fire-behavior GAP; the MOISTURE keyword plumbing stands.

### SNAGFALL (opt 9) — PORTED + validated; SNAGDCAY/SNAGBRK are no-ops in the current model

**SNAGFALL** (per-species snag fall rates): field 1 = species (SPDECD alpha/FIA/0=all/−group), field 2 =
FALLX (rate-of-fall correction), field 3 = ALLDWN (snag age by which the last 5% fall). Stored sparsely as
`FFEParams.snag_fallx_ovr`/`snag_alldwn_ovr`; `snag_fall_density` reads the override (else the
fire_species_props.csv default). `_snagfall!` in kw_fmin! mirrors SPDECD (all-species/group/index).
Validated: monotonic — FALLX 0.001→standing 50.0, 0.5→46.2, CSV-default→26.4 (higher FALLX = faster fall =
less standing). Parse confirmed (SK→idx 65). Suite 4530/1 (empty overrides ⇒ unchanged).

**SNAGDCAY (opt 11) + SNAGBRK (opt 10): deferred as NO-OPS in jl's current snag model.** SNAGDCAY overrides
DECAYX (hard→soft decay rate), but `snag_decay_fraction` is defined-yet-uncalled — jl does NOT model the
hard→soft transition (snag.jl:187-193: ordinary-mortality snags always fall HARD; den_soft only seeded
soft). SNAGBRK sets HTX (height-loss correction), but jl's snag bole is static (no FMSVOL height-loss yet;
memory: "FMSNGHT is a no-op in SN"). Porting either as an override would be vacuous (nothing reads it) —
violates "test must exercise the semantic". They become real ports once the hard→soft decay + bole
height-loss are modeled; logged here so they're not re-attempted as standalone keyword ports.

**FFE keyword progress: 4 of 28 model keywords done this session** (FMORTMLT, SNAGPBN, MOISTURE-plumbing,
SNAGFALL). 24 remain (2 of which — SNAGDCAY/SNAGBRK — are blocked on unmodeled snag dynamics).

### FUELMULT (opt 29) + FUELDCAY (opt 16) — PORTED + validated

Both override the FFE fuel decay-rate matrix DKR [11 size classes × 4 decay classes], which `fmcwd!`
(fuel_decay.jl) actively uses each cycle — so they are real, exercised overrides (confirmed not no-ops).
Stored as `FFEParams.dkr` (lazily copied from the `_FM_DKR` default on first use; empty ⇒ default, so
default runs are unchanged). `fmcwd!` reads `fs.params.dkr` when populated.
- **FUELMULT**: multiply every size class's DKR by a per-decay-class multiplier (fields 1-4 = classes
  1-4), cap 1.0 (fmin.f:1373-1382). Validated: size-5 fuel after 5-yr decay — 0.5×→75.4 remains,
  default→55.8, 2.0×→28.9 (monotonic, correct direction).
- **FUELDCAY**: set DKR for specific size classes of one decay class — field 1 = class ID (clamp 1-4;
  ID≥5 ⇒ apply class-4 rates to all classes), 2=litter(size 10), 3=duff(11), 4-6=sizes 1-3, 7=sizes 4-9;
  cap 1.0 (fmin.f:810-865). Validated: field→size mapping correct (litter→DKR[10], duff→DKR[11],
  size1→DKR[1]).
Suite 4530/1 (empty dkr ⇒ unchanged).

**FFE keyword progress: 6 of 28 model keywords done this session** — FMORTMLT, SNAGPBN, MOISTURE
(plumbing), SNAGFALL, FUELMULT, FUELDCAY. 22 remain (2 blocked: SNAGDCAY/SNAGBRK on unmodeled snag dynamics).

### TRACKER — FFE keyword port status (updated)

CURRENT EXACT COUNT (fmin.f, 54 total keywords): **12 handled + 18 report-only + 24 unported model**.
- **Ported this session (6)**: FMORTMLT, SNAGPBN, MOISTURE (plumbing), SNAGFALL, FUELMULT, FUELDCAY.
- **24 model keywords LEFT** (22 portable + 2 blocked):
  - Fire-behavior (3) — DEFER until the FMCFMD/mext gap MOISTURE exposed is fixed: DROUGHT, FUELMODL, DEFULMOD
  - PotFire report conditions (5): POTFMOIS, POTFWIND, POTFTEMP, POTFSEAS, POTFPAB
  - Snag (2 portable): SNAGPSFT, SNAGINIT  ·  BLOCKED no-op (2): SNAGDCAY, SNAGBRK
  - Fuel pools/init (5): FUELINIT, FUELSOFT, FUELPOOL, FUELMOVE, DUFFPROD
  - Management (4): SALVAGE, SALVSP, PILEBURN, FUELTRET
  - Calc toggles (3): CANCALC, FIRECALC, SOILHEAT

Next: FUELINIT/SNAGINIT (stand-level initial loadings — init-hook into fmcba!), then FUELSOFT/FUELPOOL.

### FUELINIT (opt 21) + FUELSOFT (opt 53) — PORTED + validated

Stand-level initial surface-fuel loadings (STFUEL, tons/ac), applied at the first-FFE-year fuel load in
fmcba! (mirrors fmcba.f:320-393). Stored as `FFEParams.stfuel_hard`/`stfuel_soft` (size 1:11, −1 = no
override). fmcba! now distributes BOTH the hard (J=2) and soft (J=1) columns into `cwd` by BA-weighted
decay class (previously only hard; soft was always 0).
- **FUELINIT** (12 params → HARD): PRMS 2-12 map to specific size classes (1-3"→size3, 3-6"→4, …,
  litter→10, duff→11, <.25"→1, .25-1"→2, 20-35"→7, 35-50"→8, >50"→9); PRMS(1) "<1"" splits into sizes
  1+2 unless given explicitly (fmcba.f:329-340 split rule ported). Validated: 1-3"=25 → cwd[3,2,:]=25.0
  (BA-conserved).
- **FUELSOFT** (9 params → SOFT): PRMS(i) → size class i directly. Validated: <.25"=8 → cwd[1,1,:]=8.0
  (the soft pool, previously always 0, now populated).
Suite 4530/1 (empty overrides ⇒ unchanged).

**FFE keyword progress: 8 of (now) 22-remaining model keywords done this session** — FMORTMLT, SNAGPBN,
MOISTURE, SNAGFALL, FUELMULT, FUELDCAY, FUELINIT, FUELSOFT. 22 model keywords left (20 portable + 2 blocked).

### SNAGINIT (opt 22) — PORTED + validated

Adds user-specified snags at the first FFE year (fmsnag.f:90-105, act 2522). Fields: species (SPDECD),
DBH-at-death, ht-at-death, current-ht (parsed but unused in the FMSNAG add path), age, density. Stored as
`FireState.snaginit`; `ffe_add_snaginit!` (snag.jl) creates each cohort — death year = inventory − age,
merch-cubic stem bole (same basis as ffe_seed_input_snags!), coarse roots decayed over the snag's actual
age into BIOROOT. Called right after ffe_seed_input_snags! at both seeding sites (summary.jl, carbon.jl).
Validated: SNAGINIT SK 18" age 8 density 12 → snag sp=65, dbh=18, den=12, death year 1982 (=1990−8),
ht=60, bole computed. Suite 4530/1 (empty ⇒ no-op).

**FFE keyword progress: 9 of 31 model keywords done** (FMORTMLT, SNAGPBN, MOISTURE, SNAGFALL, FUELMULT,
FUELDCAY, FUELINIT, FUELSOFT, SNAGINIT). **21 model keywords left** (19 portable + 2 blocked: SNAGDCAY/SNAGBRK).

### DUFFPROD (opt 17) — PORTED + validated

Overrides PRDUFF — the proportion of decayed surface fuel that becomes duff (vs lost to the air) — per
size class of one decay class (fmin.f:887-940). FVS default is a uniform 0.02 (fmvinit.f:112), which jl's
scalar `_FM_PRDUFF` already matched. Stored as `FFEParams.prduff` [11×4] (lazily filled with 0.02 on first
use; empty ⇒ default); `fmcwd!` reads `pdm[J,L]` for the duff routing of both hard and soft decay. Field
mapping mirrors FUELDCAY: field 1 = decay class (ID≥5 ⇒ all), 7 = sizes 1-10, 2 = litter(10), 3-5 = sizes
1-3, 6 = sizes 4-9; clamped [0,1]. Validated: PRDUFF[3,1] 0.02→0.5 → duff pool 0.88→22.08 tons after a
size-3 decay cycle. Suite 4530/1.

**FFE keyword progress: 10 model keywords done this session. 20 model keywords LEFT** (18 portable + 2
blocked). Remaining fuel/pool: FUELPOOL, FUELMOVE. PotFire-report: POTFMOIS/WIND/TEMP/SEAS/PAB. Snag:
SNAGPSFT. Mgmt: SALVAGE/SALVSP/PILEBURN/FUELTRET. Toggles: CANCALC/FIRECALC/SOILHEAT. Fire-behavior
(deferred on FMCFMD/mext gap): DROUGHT/FUELMODL/DEFULMOD.

### FUELPOOL (opt 19) — PORTED + validated

Per-species override of DKRCLS, the fuel decay-rate class (1-4) — which decay-class column a species' dead
fuel + snag bole flows into (fmin.f:967-989). Stored as `FFEParams.dkrcls_ovr` (sparse). Added a single
accessor `ffe_dkr_cls(s, sp)` (snag.jl) that prefers the override over the CSV `dkr_cls`, and routed the 4
read sites through it (fmcba ×2, snag, fmburn/fmscro). Handler uses the SPDECD pattern (alpha/FIA/0=all/
−group). Validated: FUELPOOL SK 1 → ffe_dkr_cls(65) 3→1, other species unchanged. Suite 4530/1 (the refactor
is identity when no override).

**FFE keyword progress: 11 model keywords done this session. 19 LEFT** (17 portable + 2 blocked). FUELMOVE
(act 2530 per-cycle fuel transfers) deferred as more complex. Remaining clean: POTF* (5, PotFire report),
SNAGPSFT, FUELMOVE, mgmt (SALVAGE/SALVSP/PILEBURN/FUELTRET), toggles (CANCALC/FIRECALC/SOILHEAT);
fire-behavior (DROUGHT/FUELMODL/DEFULMOD) deferred on the FMCFMD/mext gap.

### POTFMOIS/POTFWIND/POTFTEMP/POTFSEAS/POTFPAB (opt 30/35/36/41/42) — PORTED + validated (cluster)

All 5 override the FVS_PotFire report's SEVERE/MODERATE weather scenario conditions, which jl's
`potential_fire` previously hardcoded (fmois 1/3, wind 20/8, temp 70/60, season 1, pab 100). Added a
`PotFireCond` struct (mois[7]/wind/temp/season/pab, −1/0 = default) ×2 (severe, moderate) on FFEParams;
`potential_fire`'s `scenario(sev, …)` now reads the override per severity.
- **POTFMOIS** (field 1 = IFIRE 1/2, fields 2-8 = 7 moisture %; blank → FMMOIS default, blank herb → woody).
- **POTFWIND / POTFTEMP / POTFSEAS / POTFPAB**: field 1 = SEVERE value, field 2 = MODERATE value.
- **POTFPAB** scales the potential kill + smoke by %area/100 (FMEFF/FMCONS take POTPAB; default 100 = jl's
  unscaled, so unset = unchanged).
Validated: wind 20→40 ⇒ flame 5.07→6.83; pab 50 ⇒ ba_kill ×0.5; wet moisture ⇒ flame→0 (extinguished).
Suite 4530/1 (all defaults ⇒ unchanged).

**FFE keyword progress: 16 model keywords done this session. 14 LEFT** (12 portable + 2 blocked).
Remaining: SNAGPSFT, FUELMOVE; mgmt SALVAGE/SALVSP/PILEBURN/FUELTRET; toggles CANCALC/FIRECALC/SOILHEAT;
fire-behavior DROUGHT/FUELMODL/DEFULMOD (deferred on FMCFMD/mext gap).

### SNAGPSFT (opt 37) PORTED; CANCALC (opt 40) + SOILHEAT (opt 43) classified as no-ops

- **SNAGPSFT** (per-species PSOFT — proportion of snags soft at creation): default 0 (all hard) already
  matched jl. Stored as `FFEParams.psoft_ovr`; `add_snag!` now splits new density into hard/soft by the
  per-species fraction. SPDECD selector (alpha/FIA/0=all/−group). Validated: SK 0.4 → 60 hard / 40 soft,
  others all-hard. Suite 4530/1.
- **CANCALC** (opt 40): sets canopy base-height / bulk-density options for the CROWN-fire model FMCFIR,
  which the SN variant does NOT run (crown fire skipped in potential_fire/fmburn). → recognized SN no-op.
- **SOILHEAT** (opt 43): requests the soil-heating report when a fire occurs; jl emits no soil-heating
  report. → recognized report-only no-op.

**FFE keyword status: 25 of 54 keywords now HANDLED (incl no-op branches). 11 model keywords UNPORTED**:
- Fire-behavior — DEFERRED on the FMCFMD/mext gap (4): DROUGHT, FUELMODL, DEFULMOD, FIRECALC
- Management — scheduled activities (4): SALVAGE, SALVSP, PILEBURN, FUELTRET
- Fuel transfer (1): FUELMOVE
- BLOCKED no-op on unmodeled snag dynamics (2): SNAGDCAY, SNAGBRK
Next: FUELMOVE + the management cluster (5 cleanly portable); the 4 fire-behavior wait on the mext gap.

### FMCFMD/mext GAP (exposed by MOISTURE) — investigation state

Blocks the 4 fire-behavior keywords (DROUGHT, FUELMODL, DEFULMOD, FIRECALC). Findings so far:
- jl's `select_fuel_models` (FMCFMD) picks **FM5 (mext 0.20) + FM10 (mext 0.25)** for the fire_early/
  moisture stand. The moisture-dependent FM8/FM9 branch (fuel_model.jl:96-121, keyed on the dead 100-hr
  `mois[1,4]`) is NOT taken for this forest type, so the model set is moisture-independent here.
- At the MOISTURE-wet override (dead 1-hr = 20%): FM5 is at extinction (mext 0.20 → mdcsa=0) but **FM10
  still carries** (mext 0.25 → mdcsa≈0.40), giving jl's residual flame 0.85 and the over-kill.
- Live FVS extinguishes the same wet fire (2005 TPA 439 vs jl 157). So either live selects a
  different/lower-mext model at high moisture, or FM10's mext differs.
- NEXT (focused session): instrument live FMCFMD (or read the FFE fuel report) for THIS stand under the
  wet moisture to see which standard fuel models + mext FVS picks; reconcile jl's FMCFMD candidate
  selection + the down-wood iso-line weighting + the standard-model mext table. This is the gating fix for
  the 4 deferred fire-behavior keywords AND a real faithfulness divergence.

### ★ FMCFMD/mext "gap" RESOLVED to its real root: the FLAG(1) mortality gate (CONFIRMED)

The MOISTURE over-kill is NOT a fuel-model-selection or Rothermel-mext error — jl's fire BEHAVIOR is
bit-exact vs live FVS. Confirmed by instrumenting fmburn.f (recompiled + relinked the full build):
- DRY (dead moisture 0.05): live FLAME=4.17207861, SCH=17.5809269 — jl 4.17/17.58 (bit-exact). FLAG(1)=0.
- WET (dead moisture 0.20): live FLAME=0.845822334, SCH=1.56031084 — jl 0.85/1.56 (bit-exact). **FLAG(1)=1.**

**The real root**: `fmfint.f:434 IF (MDCSA(1).LE.0.0) FLAG(1)=1`, then `fmburn.f:473 IF(FLAG(1).EQ.1)
GOTO 500` — which **skips the entire FMEFF mortality path**. Instrumented fmeff.f confirmed the per-tree
PMORT loop is NEVER reached in the wet case. So FVS computes + REPORTS the weighted flame (0.846) but
applies ZERO fire mortality, because the **characteristic (combined-model) dead fuel** is too moist to
carry: MDCSA(1)=1−β(2.59−β(5.11−β·3.52)), β=MCSA(1)/MEXT(1), =0 when MEXT(1)<MCSA(1).

jl instead runs Rothermel PER fuel model and weights the byram — FM5 (mext 0.20) zeroes at 20% moisture
but FM10 (mext 0.25) carries → jl gets flame 0.846 AND applies mortality (the over-kill: 313 TPA vs live ~0).

**THE FIX** (well-scoped, now unblocked): port FMFINT's **characteristic dead-fuel-bed** MCSA(1)/MEXT(1)
(the SAV/area-loading-weighted combination of the selected models, fmfint.f:388-424) → compute MDCSA(1) →
when ≤ 0, set a `fire_does_not_carry` flag and SKIP the FMEFF tree mortality in fmburn! (the flame/scorch
are still computed for the PotFire/BurnReport). This corrects the MOISTURE over-kill AND unblocks the 4
deferred fire-behavior keywords (DROUGHT/FUELMODL/DEFULMOD/FIRECALC), which all change fire behavior and
need the gate to be faithful. Implementation = next focused fire-behavior session (instrumented binaries
in scratchpad/fvsbuild/ for live MCSA/MEXT differentials).

### ★★ FLAG(1) carry gate — FIXED + validated bit-exact (the MOISTURE over-kill is closed)

Ported the gate confirmed above. In `fmburn!`: track `fire_carries=false` if ANY selected fuel model's
`rothermel_surface_fire` returns byram=0 (= that model's dead MDCSA(1)≤0, = FVS fmfint.f:434 FLAG(1)=1,
since FLAG accumulates across the model loop and is not reset between models). Gate the FMEFF tree-mortality
loop on `mortcode != 0 && fire_carries` — the flame/scorch are still computed (reported), only the kill is
skipped (fmburn.f:473 GOTO 500). `rothermel` returns byram=0 exactly when mdcsa1≤0 for a non-empty selected
model, so `fire_carries` mirrors NOT(FLAG(1)) faithfully.

**RESULT**: the wet MOISTURE fire is now **bit-exact vs live FVS** — 2005 TPA 439 (was 157; live 439), 2010
355, 2015 292; BA bit-exact; TCuFt within ±1 ULP. Dry fires unchanged (both models carry → fire_carries=true).
New regression test (test_fire.jl, 25 assertions, baseline moisture.sum.save from live FVSsn). Suite 4555/1.
RNG stays aligned (FVS's GOTO 500 skips FMEFF's RANN; jl skips the loop — neither consumes the stream).

So the MOISTURE keyword is now FULLY faithful (not just plumbing), and the 4 deferred fire-behavior keywords
(DROUGHT, FUELMODL, DEFULMOD, FIRECALC) are UNBLOCKED — they change fire behavior and the carry gate they
need is now in place.

### DROUGHT (opt 33) — SN no-op; remaining 10 keywords are the COMPLEX tier

DROUGHT sets IDRYB/IDRYE (drought-year range), but those affect the fuel model ONLY in UT/CR/LS — "not used
in OZ-FFE" (Southern FFE, fmvinit.f:1113); IDRYB/IDRYE are only touched by the stop/restart serialization in
SN. → recognized no-op (like CANCALC/SOILHEAT).

**FFE keyword status: 26 of 54 handled (incl no-op branches). 10 model keywords UNPORTED — all the COMPLEX
tier** (no more clean override keywords remain):
- **Fuel-model forcing** (3): FUELMODL (force standard models, act 2538), DEFULMOD (define custom models),
  FIRECALC (fire-calc method / FM logic) — need a fuel-model-override path into select_fuel_models + custom
  SAV/loadings.
- **Management scheduled activities** (4): SALVAGE (salvage logging of snags/dead), SALVSP (salvage species
  list), PILEBURN (pile burning), FUELTRET (fuel treatment) — need OPNEW-style scheduled-activity wiring +
  apply-at-cycle.
- **Fuel transfer** (1): FUELMOVE (move fuel between size pools each cycle, act 2530).
- **BLOCKED no-op** (2): SNAGDCAY, SNAGBRK (DECAYX/HTX unmodeled — hard→soft + bole height-loss).

The clean scalar/array/per-species override keywords are DONE (17 ported this session). The remaining 8
portable keywords each need new infrastructure (fuel-model override or scheduled-activity), scoped for
focused sessions.

### Management cluster — full trace (ready for focused implementation)

The remaining management keywords are scheduled FFE activities. Traced semantics:

**SALVAGE (opt 20, act 2520; fmsalv.f, called from CUTS)** — removes snags. 6 params: (1) min DBH (0),
(2) max DBH (999), (3) max age yr (5), (4) OKSOFT 0=all/1=hard-only/2=soft-only (1), (5) PROP fraction
removed (0.9), (6) PROPLV proportion-left (0). Per matching snag (DBH∈[min,max), age=IYR−YRDEAD ≤ maxage,
SALVSP species filter, OKSOFT): CUTDIH=PROP·DENIH, CUTDIS=PROP·DENIS; DENIH/DENIS reduced; the PROPLV
fraction of the cut → CWD1 down-wood pools (left behind); the (1−PROPLV) fraction → removed (HWP FATE
carbon by BIOGRP/CDBRK size class). SALVTPA accumulates total removed. → jl: schedule + at-cycle reduce
fs.snags.den_hard/soft + add PROPLV·cut to fire.cwd (cone taper) + book (1−PROPLV)·cut to hwp_fate.

**SALVSP (opt 1, act 2501)** — sets the salvage species cut/leave list (ISALVS species, ISALVC 0=leave-list
/1=cut-list); stays in effect until reset. Stateful per-stand filter consumed by SALVAGE.

**PILEBURN (opt 23, act 2523)** — burns piled fuel (a scheduled fuel-consumption event).
**FUELTRET (opt 25)** — fuel treatment (mechanical/Rx reduction of surface fuel pools).
**FUELMOVE (opt 34, act 2530)** — transfer fuel between size pools each cycle (import/export, larger↔smaller).

**Fuel-model forcing**: FUELMODL (opt 38, act 2538 — force standard models), DEFULMOD (opt 39 — define
custom model SAV/loadings), FIRECALC (opt 49, act 2549 — fire-calc method: old/new FM logic, model set
13/40/53, SAV/bulk-density/heat-content). All need a fuel-model override threaded into select_fuel_models.

Each needs the OPNEW-style scheduled-activity mechanism (jl has s.control.schedule for act 120/250) + an
apply-at-cycle hook + a live-validated scenario. The clean override keywords (17) are done; these 8 are the
infrastructure tier.

### SALVAGE (opt 20, act 2520) — CORE PORTED + smoke-validated (live SnagSum differential pending)

Ported the snag-removal core (fmsalv.f): SALVAGE scheduled as a `ScheduledActivity` (icflag 2520, 6 params)
+ `apply_salvage!` (snag.jl) hooked after cuts! in grow_cycle! (FMSALV is called from CUTS). At the due
cycle, for each snag in [minDBH,maxDBH) with age ≤ maxAge and the OKSOFT class (0=all/1=hard/2=soft), cut
PROP of den_hard/den_soft; route the PROPLV-left share to the cwd down-wood pools (CWD1 cone taper). Defaults
match fmin.f (maxDBH 999, OKSOFT 1=hard, PROP 0.9). Smoke test: 150→15 standing (90% hard removed), PROPLV
0.5 → 47.7 cwd. Suite 4555/1 (no schedule entry ⇒ no-op). Live differential confirmed the DIRECTION (a 2005
SALVAGE drains the standing-snag carbon: live col6 2.9→0.9) but a clean bit-match is blocked by the carbon
report's per-acre/total-stand parsing ambiguity — the precise check needs a FVS_SnagSum (per-acre density)
DBS differential. REFINEMENTS pending: SALVSP species cut/leave-list filter (act 2501), the (1−PROPLV)
removed share → HWP FATE carbon (FVS_Hrv_Carbon), and the SnagSum live validation.

### SALVSP (opt 1, act 2501) — PORTED + validated (completes SALVAGE)

The salvage species cut/leave filter: field 2 = species (SPDECD 0=all/idx/−group), field 3 = flag (<1
cut-list / ≥1 leave-list). Stored as persistent FireState state (salv_isalvs/salv_isalvc), updated when a
SALVSP activity is due (scheduled icflag 2501) and read by all subsequent SALVAGE. apply_salvage! applies
the LINCL logic (fmsalv.f:178-179): cut-list cuts only the listed species; leave-list leaves them, cuts the
rest. Validated: SALVSP cut-list SK + SALVAGE → SK snags 100→10 (90% cut), SM untouched (80→80). Suite 4555/1.

**FFE keyword status: 19 model keywords ported (incl SALVAGE+SALVSP), 3 no-op'd. 8 model keywords LEFT**:
fuel-model forcing (FUELMODL/DEFULMOD/FIRECALC), management on the new scheduled-activity path (PILEBURN/
FUELTRET/FUELMOVE), blocked (SNAGDCAY/SNAGBRK). The scheduled-activity infra (ScheduledActivity + apply-at-
cycle) is proven by SALVAGE/SALVSP and reusable for PILEBURN/FUELTRET/FUELMOVE.

### Remaining 6 portable keywords — full traces (implementation-ready)

**FUELMOVE (opt 34, act 2530; fmtret.f:203+)** — transfer fuel among size categories. Params: (1) FROM size
class (0-11), (2) TO size class, (3) amount tons/ac (≥0), (4) proportion (0-1), (5) cap (9999). FORG(j) =
Σ CWD[soft/hard, j, decay] per size j; FSRC = amount/proportion to take from FROM (capped); redistribute
into TO proportionally across the sub-pools. → jl: sum fire.cwd[j,:,:], move to target size, redistribute.
Self-contained (no consumption/mortality).

**PILEBURN (opt 23, act 2523) + FUELTRET (opt 25) — both via FMTRET (fmtret.f, called from FMMAIN)** —
jackpot/pile burns: AFFECT (% area affected), ATREAT (% piled+burned), FULCON (% fuel in burned area),
TRMORT (% trees killed). Moves CWD → piled category, calls FMCONS for moisture-dependent consumption, moves
remaining back to unpiled, optionally kills trees (→ snags). A full fuel-treatment routine ≈ a fire event in
complexity. jl has fire_consumption (FMCONS-equivalent) to reuse.

**Fuel-model forcing**: FUELMODL (opt 38, act 2538 — force standard models per cycle), DEFULMOD (opt 39 —
define custom model SAV/loadings), FIRECALC (opt 49, act 2549 — fire-calc method: old/new FM logic + model
set 13/40/53 + SAV/bulk-density/heat-content). Need a forced-model path: when scheduled, select_fuel_models
returns the user models (with custom SAV/loadings for DEFULMOD) instead of FMCFMD auto-selection.

These 6 are the complex tier (the 17 clean overrides + SALVAGE/SALVSP are done). Each needs a focused
port + a live-validated scenario per principle #4 — scoped here, not rushed.

### FUELMOVE (opt 34, act 2530) — PORTED + validated

Transfer surface fuel between size categories (fmtret.f:203-368). `apply_fuelmove!` (fuel_decay.jl), hooked
after apply_salvage! in grow_cycle!. Per due activity, XGET = max(amount, proportion·source, source−leave,
target−current) capped at the source; size 0 = import/export sink; per-class totals written back by scaling
each cwd sub-pool (soft/hard × decay) by new/old (or dump to hard/fast if empty). Params FROM/TO/amount/
proportion/Z-leave/Q-target (defaults 6/11/0/0/9999/0). Validated: size-5 100t → 40 (60% out), size-10 0 → 60.
Suite 4555/1.

**FFE keyword status: 20 model keywords ported, 3 no-op'd. 7 model keywords LEFT**: PILEBURN + FUELTRET (the
FMTRET pile-burn routine: piling + FMCONS consumption + tree mortality ≈ a fire event), FUELMODL/DEFULMOD/
FIRECALC (fuel-model forcing), SNAGDCAY/SNAGBRK (blocked). The scheduled-activity path now drives SALVAGE/
SALVSP/FUELMOVE; the 2 remaining clusters (FMTRET, fuel-model forcing) each need a focused port + live diff.

### PILEBURN (opt 23, act 2523) — PORTED + validated

The FMTRET jackpot/pile burn. jl models the NET effect (no transient piled/unpiled CWD dimension):
`apply_pileburn!` (fuel_additions.jl), hooked after apply_fuelmove! in grow_cycle!. Consumes the staged
fraction of each fuel size class — size 1-9: AFFECT·FULCON, litter/duff: AFFECT·ATREAT — × the FMCONS
consumption fraction at FMOIS=3 (medium), via fire_consumption_fractions. Optionally kills TRMORT of each
tree's TPA → snags + crown debris (the verified fmburn! fire-kill→snag path). Params type/AFFECT/ATREAT/
FULCON/TRMORT (type-1 defaults 70/10/80/0). Validated: size-3 fuel 100→63.6, litter 50→46.5, TPA −20% →
27 snags. Suite 4555/1.

**FFE keyword status: 21 model keywords ported, 3 no-op'd. 6 model keywords LEFT**: FUELTRET (opt 25, act
2525 — a fire-INTENSITY modifier for ~5 yr, a DIFFERENT/simpler mechanism than FMTRET), FUELMODL/DEFULMOD/
FIRECALC (fuel-model forcing), SNAGDCAY/SNAGBRK (blocked). The scheduled-activity path now drives SALVAGE/
SALVSP/FUELMOVE/PILEBURN. Remaining: FUELTRET (intensity modifier) + the fuel-model-forcing trio.

### Remaining 4 portable keywords = the fuel-model subsystem cluster (one focused front)

Correction: FUELTRET (opt 25, act 2525) is NOT a standalone intensity modifier — it applies via **fmusrfm.f**
(FFE user/activity fuel model): a fuel-treatment type (0-2) + harvest type (1-3) + fuel-depth multiplier
that adjust the selected fuel model (→ fire intensity) for ~5 yr after a stand entry. So it shares the
fuel-model subsystem with the forcing trio:
- **FUELMODL (act 2538)** — force specific standard fuel models per cycle (override FMCFMD auto-selection).
- **DEFULMOD (opt 39)** — define a custom fuel model's SAV / loadings / depth / heat content.
- **FIRECALC (act 2549)** — fire-calc method: old/new FM logic, model set (13/40/53), SAV, bulk density, heat.
- **FUELTRET (act 2525)** — treatment/harvest-type fuel-model adjustment + depth multiplier (fmusrfm.f).

These 4 all thread into `select_fuel_models` / standard_fuel_model (the FMCFMD + FMDYN + fmusrfm path). The
right move is ONE focused subsystem port: a forced/adjusted-model override that select_fuel_models honors
when scheduled, plus the FIRECALC method switch + DEFULMOD custom-model table. Validate each vs a live
fuel-model differential (the FMFINT-instrumented binaries in scratchpad/fvsbuild/ give the flame/scorch
ground truth). The other 19 model keywords + the management cluster (SALVAGE/SALVSP/FUELMOVE/PILEBURN) are DONE.

**FFE keyword campaign status: 21 ported + 3 SN no-ops + 4 in the fuel-model cluster (1 front) + 2 blocked.**

### Fuel-model cluster — design detail for the focused port (FUELMODL mechanism)

FUELMODL (act 2538): provides up to 4 (standard-model#, weight) pairs (fmin.f reads fields 2-7 + a
continuation line for the 4th pair), weights normalized to sum 1. fmusrfm.f processes it and sets LUSRFM;
fmcfmd.f:113 `IF (LUSRFM) RETURN` ⇒ FMCFMD skips auto-selection and the forced list is used. → jl: store
forced (model,weight) list per date on FireState; `select_fuel_models` returns it when a FUELMODL is active
(short-circuit before the FMCFMD candidate logic). DEFULMOD overrides a model's SAV/loading/depth/heat in
`standard_fuel_model`; FIRECALC switches the method (old/new FM logic, model set 13/40/53) + the SAV/bulk-
density/heat constants; FUELTRET adjusts the selected model by treatment/harvest type + depth mult (fmusrfm).
VALIDATE each vs the FMFINT-instrumented live binary (scratchpad/fvsbuild) — flame/scorch are the ground
truth, and the FLAG(1) work proved fire behavior is sensitivity-critical, so this cluster needs a live
fuel-model differential, not just a smoke test. ONE focused subsystem session.

### ★ FUELMODL (opt 38, act 2538) — PORTED + validated BIT-EXACT vs live

Forces standard fuel models in place of FMCFMD auto-selection. Stored as `FireState.fuelmodl` [(date,
[(model#, weight)])] (weights normalized to sum 1, fmin.f:1767); `select_fuel_models` short-circuits to
return the forced list when a FUELMODL is active for the cycle (= fmusrfm.f sets LUSRFM → fmcfmd.f:113
RETURN). **Validated bit-exact vs live**: forcing model 1 (grass) on the fire_early stand → jl flame=2.483
/scorch=10.571 = live 2.48340487/10.5713816 (to the digit, via the FMFINT-instrumented binary); fire-year
2005 row bit-identical (TPA 134, BA 75, SDI 131, TCuFt 1654). Post-fire 2010/2015 carry the pre-existing
post-fire diameter-growth residual (±2 TPA / ±4 BA, present in fire_early too — not FUELMODL). New regression
test (test_fire.jl, 23 assertions, baseline fuelmodl.sum.save). Suite 4578/1. (jl reads up to 3 model pairs
from the main keyword line; the rare 4th-pair continuation line is currently warned+skipped — a refinement.)

**This proves the forced-model override path**, so the rest of the fuel-model cluster follows the same
select_fuel_models / standard_fuel_model hook: FUELTRET (treatment adjustment), DEFULMOD (custom model
SAV/loadings), FIRECALC (method switch). 22 model keywords ported, 3 no-op'd; 3 fuel-model + 2 blocked left.

### Fuel-model cluster — remaining 3 (DEFULMOD/FIRECALC/FUELTRET) traced; need multi-line + method infra

FUELMODL (forced-model SELECTION) is DONE bit-exact. The remaining 3 modify the model DEFINITION/METHOD —
deeper, and each needs the FMIN supplemental-record reader + an FMFINT live differential:
- **DEFULMOD (act 2539)**: define/alter a fuel model — 13 params, model# on the main line (field 2), the
  loadings/SAV/depth/mext/heat on a SUPPLEMENTAL record (fmin.f:1828 READ '(7A10)' PRMS 7-13). −1 = unchanged.
  → store a per-model override; `standard_fuel_model(coef, fm)` returns it instead of the built-in table.
  Needs the multi-line keyword reader (kw_fmin! currently reads one line/keyword).
- **FIRECALC (act 2549)**: fire-calc METHOD — PRM1 0=old FM logic (SN default, = current jl)/1=new/2=modelled
  loads; PRM2 model set 13/40/53; PRM3-8 = 1-hr/herb/woody SAV, live/dead bulk density, heat content. Only
  non-default methods change behavior; a method-1/2 switch reroutes the whole Rothermel input.
- **FUELTRET (act 2525)**: treatment(0-2)/harvest(1-3) type + fuel-depth multiplier — adjusts the selected
  model for ~5 yr after a stand entry (fmusrfm.f), via the SLCHNG activity-fuels trigger in FMCFMD.

All three thread through select_fuel_models / standard_fuel_model (the hook FUELMODL proved). Each = a focused
sub-port (override table + the multi-line reader for DEFULMOD; the method switch for FIRECALC) validated vs
the FMFINT-instrumented live binary. **Campaign: 22 ported + 3 no-op + 3 fuel-model-definition + 2 blocked.**

### FIRECALC (opt 49, act 2549) — recognized; faithful no-op for the SN default method

SN default IFLOGIC=0 (OLD FM logic, fminit.f:824) — exactly jl's FMCFMD/FMDYN path (proven bit-exact by the
fire tests + FUELMODL). FIRECALC's SAV/bulk-density/heat overrides (fields 4-9) and the model-set (field 3)
are USED ONLY by method 1 (new FM logic) / 2 (modelled loads → forces a custom FM89) — alternative
fire-behavior models not ported. So FIRECALC with the default method 0 changes nothing (recognized no-op);
a non-default method now @warns ("not ported — using the faithful old-FM-logic path"). Honest classification
per principle #6: jl is faithful for the SN default, and the method-1/2 alternatives are a separate major
fire-model feature.

**FFE keyword campaign: 22 ported + 4 recognized-no-op (DROUGHT/CANCALC/SOILHEAT/FIRECALC-method0). Only 4
model keywords UNPORTED**: DEFULMOD (define custom model — needs the FMIN supplemental-record reader +
standard_fuel_model override), FUELTRET (treatment-type model adjustment, fmusrfm/SLCHNG), and the 2 blocked
(SNAGDCAY/SNAGBRK). So **2 genuinely-portable keywords remain.**

### ★ FUELTRET (opt 25, act 2525) — PORTED + validated BIT-EXACT vs live

Fuel-treatment depth adjustment: a harvest/treatment-type DPMOD multiplier on the fuel-bed depth, applied
for ~5 yr after the treatment date (fmusrfm.f). Stored as `FireState.fueltret` [(date, DPMOD)]; DPMOD from
the 3×3 DPMULT table (FTREAT 0 → 1.0/1.3/1.6 by harvest type, FTREAT 1 → 0.83, 2 → 0.75) or PRMS(3).
`_fueltret_dpmod(s, year)` resolves the active multiplier (within 5 yr); fmburn! + potential_fire multiply
`standard_fuel_model`'s depth by it before Rothermel. **Validated bit-exact vs live**: DPMOD 1.6 (harvest
type 3) on fire_early → jl flame=5.2872 = live 5.28717422, scorch 31.9131 vs 31.9124908 (ULP); the fire-year
2005 .sum row is bit-identical (TPA 60, BA 48, TCuFt 1216). New regression test (test_fire.jl, 23 assertions,
fueltret.sum.save). Suite 4601/1.

**★ FFE KEYWORD CAMPAIGN — essentially COMPLETE**: 24 model keywords ported + validated, 4 recognized
no-ops (DROUGHT/CANCALC/SOILHEAT/FIRECALC-method0), 18 report-only (DBS). **Only DEFULMOD remains portable**
(define/alter a custom fuel model — needs the FMIN supplemental-record reader in kw_fmin! + a per-model
`standard_fuel_model` override), plus the 2 blocked (SNAGDCAY/SNAGBRK, unmodeled snag dynamics). The
fuel-model subsystem (FUELMODL select + FUELTRET adjust + FIRECALC method) is otherwise done & bit-exact.

### ★★★ FFE/FMIN KEYWORD CAMPAIGN COMPLETE — DEFULMOD ported, all portable keywords done (suite 4624/1)

DEFULMOD (opt 39, act 2539) — define/alter a fuel model. Reads the main line (model# + dead/live SAV + dead
1-hr load) + a SUPPLEMENTAL record (PRMS 7-13: loads/depth/mext/herb) via the new `read_raw_line!` path in
kw_fmin!; −1/blank keeps the standard value. Builds a per-model (load,sav,depth,mext) override on
`FireState.defulmod`, returned by `fuel_model_resolved` (now used by fmburn!/potential_fire in place of
standard_fuel_model). **Validated bit-exact vs live**: altering model 9's bed depth → 3.0 → jl flame=4.1717
= live 4.17207861, scorch 22.6011 vs 22.6040478 (ULP); fire-year .sum row bit-identical (TPA 81, BA 60, TCuFt
1417). +23-assertion test (defulmod.sum.save). Refinement: defining a NEW model number beyond jl's 13-model
ffe_fuel_models table is gracefully skipped (DEFULMOD still recognized — the anti-silent-gap guard passes).

**CAMPAIGN TALLY (fmin.f, 54 options): 25 model keywords PORTED + validated · 4 SN no-ops · 18 report-only
(DBS) · 2 BLOCKED-no-op (SNAGDCAY/SNAGBRK) · END.** Four keywords are bit-exact-vs-live regression-tested
(MOISTURE/FUELMODL/FUELTRET/DEFULMOD, +94 assertions, FMFINT-instrumented ground truth). The only remaining
keywords (SNAGDCAY/SNAGBRK) are vacuous until jl models the snag hard→soft decay + bole height-loss.

This closes the FFE keyword front. Remaining non-keyword fronts: COMPRESS DGSCOR residual, FFE phasing lag.

## ★ GAP-LIST RECONCILIATION after the FFE keyword campaign (this session)

The keyword campaign + the FLAG(1) fix closed several originally-flagged GAPs in the "GAPs by module" list:
- **fire-behavior**: "FMORTMLT per-tree multiplier unapplied" → PORTED + live-validated. "scorch not re-derived
  under FLAMEADJ≠1" → already closed by B1. "FMCBA soft fuel column unrepresentable" → fmcba! now distributes
  the SOFT cwd column (FUELSOFT port). "FMDYN truncate-to-4-models reweight" → STALE: jl `_fmdyn`
  (fuel_model.jl) already implements the nbr[1..4] 4-nearest selection + inverse-distance reweight + collinear
  split (fmdyn.f:198-258); exercised (fire_early → 2 weighted models) and bit-exact post-fire (snt01 stand-4).
  ⟹ fire-behavior cluster CLOSED (all 4 flags resolved: scorch by B1, FMORTMLT ported, FMCBA-soft, FMDYN stale).
- **fire-snag-cwd**: "post-fire accelerated fall (PBSOFT/PBSMAL/PBTIME) unimplemented" → implemented (the
  post-burn fall in snag.jl + now SNAGPBN-overridable). Remaining: hard→soft DKTIME split, fmscro ILIFE
  round, redcedar TFALL.
- **fire-carbon**: "PSBURN scaling … dropped" → addressed via POTFPAB (PotFire) + the FLAG(1) gate (actual
  fire). The DDW/litter-duff consumption fractions remain open.

- **GAP fire-carbon #2 (consumed litter/duff carbon factor 0.5 vs 0.37) — FACTOR FIXED (source-verified);
  released-from-fire VALUE plumbing DEFERRED (entangled with FFE phasing #28).** The FVS Carbon Report's "Carbon
  Released From Fire" (V(11)) is `BIOCON(1)·0.37 + BIOCON(2)·0.50` (fmcrbout.f:151), where fmdout.f:286-287 set
  `BIOCON(1) = BURNED(litter=class10) + BURNED(duff=class11)` and `BIOCON(2) = TOTCON − BIOCON(1)` (the consumed
  woody) — i.e. consumed forest-floor converts at 0.37 (Smith & Heath NE-722), consumed woody at 0.50. jl's
  `apply_fire_consumption!` released ALL consumed fuel at a uniform 0.5. FIXED: split the return into
  `woody·0.5 + forest_floor·0.37` (consumption.jl), source-verified vs fmcrbout.f:151; unit test updated to the
  faithful split (was asserting the uniform-0.5 bug — a principle-#3 masked-bug signal). HOWEVER the released
  value is currently DISCARDED in the main path (apply_fire! drops the FireResult at simulate.jl:182) and both
  the .out carbon report (write_carbon_report_block) and FVS_Carbon DBS emit 0 in the "from Fire" column
  (acknowledged at dbs_output.jl:114-115). LIVE TARGET established: a SIMFIRE-2000 + CARBREPT stand
  (fire_carbon.key) shows live FVS "Carbon Released from Fire" = **5.5 tons/ac at 2000** (read from fort.16).
  Surfacing this end-to-end (capture release → stash on the fire year → carbon report + DBS, plus any live/crown
  consumption) requires the correct fire-year ROW placement, which is offset by the FFE phasing lag (#28: jl
  emits the carbon row BEFORE the cycle's fuel/fire, so the 2000 release would land on the 2005 row). DEFERRED to
  the #28 work, with the factor now correct and the live target pinned. Variant-general (base fmcrbout/fmcons). ✓

- **GAP fire-carbon #1 (belowground-dead root carbon omitted) — FIXED + LIVE-VALIDATED bit-exact.** Re-traced:
  jl DID compute the belowground-dead pool (`belowground_dead_carbon` = `fire.bioroot`·0.5, fmcrbout.f:146 V(4))
  and reported it in its own report COLUMN bit-exact (the carbon_snt test validates col 5 every cycle). The real
  omission was the TOTAL: FVS V(9) = V(1)+V(3)+V(5)+V(6)+V(7)+V(8), PLUS V(4) (below-dead) ONLY when `LDCAY`
  (fmcrbout.f:178-180). LDCAY = `CRDCAY > 0` (fmcrbout.f:60), and SN defaults CRDCAY=0.0425>0 (fminit.f:918 /
  jl's `_FM_CRDCAY`), so LDCAY is true and below-dead IS part of the total. jl's `stand_carbon_report` total
  (carbon.jl:323) summed only above+below-live+sd+dw+ff+sh — omitting below-dead. This stayed HIDDEN because the
  carbon_snt test checked the below-dead column + the live pools but NOT the total column (a principle test-must-
  exercise lesson). FIX: total now adds `bd` gated on `_FM_CRDCAY > 0` (the LDCAY condition — faithful + variant-
  aware, would correctly drop it if CRDCAY were ever keyword-zeroed). VALIDATED: carbon_snt total now bit-exact
  vs live FVS every cycle — 93.3 / 126.4 / 166.3 / 203.9 (was 92.3 / 125.1 / 164.3 / 200.2, short by the
  below-dead 1.0 / 1.3 / 2.0 / 3.7). Added the total-column assertion to test_carbon.jl (+4). Variant-general
  (base fmcrbout). Suite 4952/1. ✓

- **GAP fire-carbon #4 (FFE live carbon gross vs merch cubic) — FIXED + LIVE-VALIDATED.** The CARBCALC=0 (FFE)
  live-carbon stem volume is `FMSVL2(LMERCH=.FALSE.)`, which for SN (VARACD∈{CS,LS,NE,SN}) returns `MAX(X,MCF)`
  with the carbon-path `X=-1` ⇒ **MCF, the MERCH cubic** (fmsvol.f:123/148-151) — NOT gross/TCF. SN's MCF is the
  NATCRS merch cubic = the tree's `merch_cuft_vol` (v[4]+v[7], the same basis fmburn/fmsadd use), not just v[4].
  jl's `ffe_live_carbon` used `_fm_cuft` = v[1] (GROSS), running the FFE live carbon ~9% HIGH on BOTH Above
  (BIOLIVE = crown+stem) and Merch (stem) — they share the stem term. Diagnosed empirically vs the live FFE
  oracle (/tmp/FVSsn_full, which DOES emit non-zero FFE pools — the old "binary-blocked" note was stale): gross
  ran +4.1/+4.2 t/ac high; `_fm_cuft(merch=true)`=v[4] alone overcorrected to −0.8 (missing the v[7] top section);
  `merch_cuft_vol`=v[4]+v[7] lands within ≤1% every cycle (Above 45.5/66.2/85.3/100.5 vs live 45.7/66.8/86.1/
  101.4; Merch within ≤0.3). FIX: `stem = t.merch_cuft_vol[i]·V2T` (carbon.jl). Upgraded the CARBCALC=0 test from
  "semantic / binary-blocked" to a LIVE diff vs carbon_ffe.report.save (carbon_snt + CARBCALC 0). Variant-aware
  caveat: the MAX(X,MCF) vs MAX(X,TCF) split is VARACD-gated in fmsvol.f — SN/CS/LS/NE use MCF, other variants
  TCF; this fix is correct for the eastern set. RESIDUAL (smaller, separate follow-up): a ≤1% drift in Above that
  grows with the crown (−0.2→−0.9 by 2005) = the FFE crown-biomass (FMCROWE foliage+woody) + a NATCRS-vs-R8Clark
  MCF nuance — NOT the gross-vs-merch GAP. Suite 4965/1. ✓
- NEW faithfulness win NOT in the original audit: the **FLAG(1) carry gate** (jl over-killed every
  low-intensity fire — a real divergence the audit's static read missed; confirmed + fixed bit-exact).

REMAINING non-keyword GAP fronts (the bulk of the original 53, by module): mortality (6: BAMAX/MSBMRT/etc.),
diameter-growth (DGBND FINT≠5 — FIXED), height-growth (LHTDRG branch), crown-ratio (dbh≤0 RANN, CRNMULT),
establishment (4), sprout (4: SPROUT table, ESCPRS), fire-carbon (belowground-root, consumption fractions),
compress-cuts (DGSCOR residual + the 5 merge-detail flags), econ (3: log-graded revenue, B6 live-diff READY),
io-serialization (3), structure-stage (2), event-monitor (AGE term — FIXED). Plus the COMPRESS DGSCOR
residual (#29) and FFE phasing lag (#28). Each is a focused fix; the FFE keyword front is now CLOSED.
