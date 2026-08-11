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
- **#143 IE AUTOES** — FIXED+VALIDATED (d089b78) per memory; goal doc stale.
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
