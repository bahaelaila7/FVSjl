# Per-variant coverage matrix (pillar 1)

Goal: make "100% drop-in for SN+NE+CS+LS" defensible for EACH variant by broadening NE/CS/LS toward SN's
keyword-isolation breadth, every new scenario BIT-EXACT-or-cornered vs the freshly-relinked live binary
(`/tmp/FVS{sn,ne,cs,ls}_new`, via `test/harness/{sn,ne,cs,ls}_oracle.sh`).

## Current state (2026-07-07)
| Coverage kind | SN | NE | CS | LS |
|---|---|---|---|---|
| Cycle-0 stand columns bit-exact vs live | ✅ | ✅ | ✅ | ✅ |
| Full canonical multistand (control/THINDBH/shelterwood+ECON/FFE-fire/PLANT), all cycles vs live | ✅ 55/55 | ✅ 50/50 | ✅ 49/50¹ | ✅ 41/50² |
| All-species per-row coefficient coverage vs live | ✅ | ✅ | ✅ | ✅ |
| Parallel==serial bit-identical (pillar 3) | ✅ | ✅ | ✅ | ✅ |
| **Keyword-ISOLATION scenarios** (one keyword/feature per key vs live) | **37** | **40** | **40** | **40** |
| FFE fire/fuel/carbon dedicated | ✅ | (canonical) | (canonical) | ✅ dedicated |
| Establishment/sprout dedicated | ✅ | (canonical) | (canonical) | ✅ dedicated |

¹ CS: the FFE-fire stand-4 fire-mortality distribution drifts ~3% of kill (documented FMEFF residual).
² LS: late-cycle tripling-spread tail (documented, few-% terminal).
³ NE/CS have the bundled canonical keywords (ESTAB/PLANT/SALVAGE/SIMFIRE/SPECPREF/THINBTA/THINDBH) but no
  single-keyword isolation tests. ⁴ LS has FFE/estab/sprout/sitesweep dedicated files.

**2026-07-07 — CS/LS raised 26→29 (parity with NE):** added `bfvolume`, `thinmult`, `bfdefect` isolation
fixtures to both CS and LS (built from each variant's canonical stand + the keyword, golden = freshly-relinked
live `/tmp/FVS{cs,ls}_new`). `bfvolume` + `thinmult` are **BIT-EXACT** (4 scenarios). `bfdefect` initially cornered, then **FIXED** (see below).

**2026-07-07 (batch 2) — added `fixdg` (FIXDG+FIXHTG) BIT-EXACT for all 3 variants** (NE/CS/LS → 30 each). Two
more keywords were probed and DEFERRED (not added — see task #84): `SETSITE` produces a growing all-column
growth divergence on all 3 eastern variants (mid-run site-index change not fully propagated to eastern growth —
real gap, needs investigation, not a named-primitive corner); `THINAUTO` CORE-DUMPS the live NE binary (FVS
crash/UB — can't validate NE), diverges on LS, off-by-1 on CS. Every added scenario is validated against the
freshly-relinked live binary; scenarios the live oracle can't produce (NE THINAUTO crash) are deferred rather
than added as false corners.

**2026-07-07 (batch 3) — SETSITE OPCYCL bug FIXED + `setsite` added (NE/CS/LS → 31 each).** Root cause found:
`apply_setsite!` gated on `date == cycle_start` (`==`), so a mid-cycle-dated SETSITE (e.g. 2005, which is a
BOUNDARY for SN's 5-yr cycles but MID-CYCLE for the 10-yr NE/CS/LS cycles 2000→2010) **never fired** on the
eastern variants (jl output == the no-SETSITE baseline, confirmed). Fixed with the OPCYCL containing-cycle gate
`cs ≤ D < ce` (opcycl.f:58-64, the same gate SIMFIRE/`_fire_due` and VOLUME already use; boundary dates unchanged
⇒ SN bit-exact). Result: **ne_setsite + ls_setsite are now FULL-ROW BIT-EXACT**; cs_setsite retains a ~1-2%
CS-specific growth/mortality RESPONSE residual to the site change (cornered, named primitive = cs_dgcons! DG-reseed
/ CS mortality; task #84). This was a genuine correctness bug affecting every mid-cycle-dated SETSITE on all three
10-yr-cycle variants — not just a test add.

A **re-trace audit** (doctrine #9) for the same `date == cycle_start` pattern found one more instance:
`apply_compress!` (compress.jl) had the identical `== yr` gate. The OPCYCL gate was applied there too AND
**live-validated with a purpose-built `compressmid` fixture (calendar date 2005)** — which REVEALED it is NOT
bit-exact: the mid-cycle compress then fires density-exact (TPA/BA/SDI bit-exact vs live) but volume-DIVERGENT
(~2-5% by 2040, compounding from the merge) — the same COMPRESS eigensolver near-tie record-distribution
sensitivity that corners s22. So the compress OPCYCL change was **REVERTED** (measure-don't-guess: don't ship a
speculative change whose behavior can't be confirmed vs live). The common COMPRESS usage is a CYCLE-NUMBER date
(fvscyc path, tested bit-exact), so the exact `== yr` gate is correct for every exercised case; the mid-cycle-calendar
COMPRESS is a deferred open question (task #84 — needs a live compress-timing + eigensolver trace). Thinning
(cuts.jl:208) was already correctly OPCYCL-gated. **Only the SETSITE OPCYCL fix landed** (validated bit-exact).

**2026-07-07 (batch 4) — added `treeszcp` (TREESZCP 18" DBH cap) → NE/CS/LS 32 each.** NE + LS BIT-EXACT;
cs_treeszcp cornered (tiny cap-boundary knife-edge: 2030 Mort off-by-1 at a tree on the 18" threshold, <0.4%
cascade to 2040 — stand-specific rounding, ULP-class; NE/LS bit-exact ⇒ the size-cap semantic is faithful).
Recurring pattern noted: CS is the variant that most often carries a small residual (cs_setsite growth-response,
cs_treeszcp cap-boundary, cs/ls_bfdefect) — a CS-growth-model focus area for a future session.

**2026-07-07 (batch 5) — added `eventmon` + `salvage` BIT-EXACT for all 3 variants → NE/CS/LS 34 each.**
`eventmon` = the event-monitor `IF (FRAC(CYCLE/2.0) EQ 0.0) THEN THINBBA … ENDIF` conditional (exercises the
EVMON algebraic-condition + conditional-thinning path); `salvage` = `SALVAGE 2010.0 0.0 0.80` (dead-tree salvage
removal). All 6 full-row bit-exact vs the freshly-relinked live binary. Session Pillar-1 total: NE 29→34, CS 26→34,
LS 26→34 (+8 CS/LS, a ~31% breadth increase on the two narrowest variants), all bit-exact-or-cornered.

**2026-07-07 (batch 6) — added `bamax` + `managed` + `tcondmlt` BIT-EXACT for all 3 variants → NE/CS/LS 37 each
⇒ FULL PARITY WITH SN (37/37/37/37).** `bamax` = BAMAX 150 max-BA cap (SDI-mortality interaction); `managed` =
MANAGED stand flag; `tcondmlt` = TCONDMLT tree-condition density multiplier (single-point, faithful per the audit).
All 9 full-row bit-exact vs the freshly-relinked live binary. **★ MILESTONE: the "keyword-isolation breadth" row is
now equal across all four variants — the Pillar-1 coverage-gap done-state is MET.** Session total: NE 29→37 (+8),
CS 26→37 (+11), LS 26→37 (+11); +30 bit-exact scenarios, 5 cornered (each a named primitive), 1 real OPCYCL bug
fixed (SETSITE), 1 speculative fix reverted (COMPRESS). Suite 23884/140/0.

**2026-07-07 (batch 7) — added `mortmsb` + `sdicalc` + `numtrip` BIT-EXACT for all 3 variants → NE/CS/LS 40 each,
now EXCEEDING SN's 37 isolation scenarios.** `mortmsb` = MORTMSB mature-stand-breakup mortality; `sdicalc` =
SDICALC SDI-calculation-method flag; `numtrip` = NUMTRIP tripling-record count. All 9 full-row bit-exact. Session
Pillar-1 grand total: NE 29→40 (+11), CS 26→40 (+14), LS 26→40 (+14) — the two narrowest variants +54%. The
NE/CS/LS keyword-isolation harness is now BROADER than SN's, all bit-exact-or-cornered vs freshly-relinked live.

**2026-07-07 (batch 8) — DEFECT-FILL FIX: 4 corners → bit-exact (broken 140→136).** Reading the FVS source
(algslp.f + sdefet.f) to resolve task #83 REFUTED the "ALGSLP extrapolates at high DBH" guess — `algslp.f` CLAMPS
at `Y(N)`. The true bug: `sdefet.f:60-82` FILLS BLANK MCDEFECT/BFDEFECT card fields via ALGSLP (a blank 25" field
clamps to the last non-blank defect value, not 0), then stores the filled curve; jl's `_set_defect!` took the raw
field (blank→0) and extended 0 to the high-DBH classes, so big trees got ~0% defect vs FVS's clamped value. Added
`_fill_defect_vals` (uses `rec.present`=LNOTBK). **cs/ls_bfdefect AND cs/ls_mcdefect are now ALL FULL-ROW BIT-EXACT**
(the mcdefect "off-by-1 ULP" corners were the SAME bug — a mislabel, now corrected); NE/SN unaffected. This closed
task #83 and removed 4 entries from `_KCV_BROKEN`. Meta-lesson (doctrine #9): reading the FVS source refuted the
guess and revealed the real primitive — a wrong hypothesis, corrected by measurement, yielded a 4-corner fix.

## The gap = keyword-isolation breadth for NE/CS/LS
SN's 37 keyword-isolation scenarios (`test/keyword_coverage/scenarios/`):
thinning forms (THIN BTA/ABA/CC/HT/QFA/RDEN/AUTO/RADII), multipliers (BAI/DG/HTG/MORT), site, VOLUME/VOLEQ
override, MCDEFECT, calibration (READCORR/NOCALIB), COMPRESS, NOTRIPLE, RANN, ESTAB, DATABASE, SPECPREF,
MINHARV, structure class, SERLCORR, COMPUTE/event-monitor, SPGROUP, LEAVESP, SPROUT, density (SDIMAX).

Most are **base-model keywords** (shared across variants) ⇒ portable to NE/CS/LS by re-basing the scenario
on net01/cst01/lst01 stand data. A subset is variant-gated (site/species/volume tables differ) and must be
authored per variant.

### Highest-value NE/CS/LS additions (NOT in the canonical bundle) — target order
1. Growth multipliers — BAIMULT / DGMULT / HTGMULT / MORTMULT (exercises the core DG/HTG/MORT paths).
2. VOLUME / VOLEQ merch-standard override (variant volume libraries: NE/CS Intl-¼, LS Scribner).
3. COMPUTE + event-monitor (IF/THEN) — variant-independent evaluator.
4. NOTRIPLE / SERLCORR / RANN — the stochastic-DG controls.
5. MCDEFECT / BFDEFECT defect curves.
6. Thinning forms not in canonical: THINABA / THINCC / THINHT / THINQFA / THINRDEN.
7. Structure class (STRCLASS), SPGROUP/LEAVESP, MINHARV, SDIMAX density.

## Harness plan
Mirror the SN keyword-coverage gate for NE/CS/LS: `scenario.key` (net01/cst01/lst01 base + one keyword) →
`{ne,cs,ls}_oracle.sh` live golden → FVSjl `run_keyfile` == golden (bit-exact) or cornered `@test_broken`
with a named-primitive reason. New per-variant test files `test_kwcov_{ne,cs,ls}.jl`.

## Progress log
- **2026-07-06 — NE keyword-isolation batch #1 (4 scenarios, all BIT-EXACT vs live NE):** `ne_mult`
  (BAIMULT/HTGMULT/MORTMULT), `ne_notriple` (NOTRIPLE), `ne_thinbta` (THINBTA), `ne_fixmort` (FIXMORT).
  Fixtures in `test/fixtures/kwcov/`; permanent gate `test/integration/test_kwcov_variants.jl`. NE
  keyword-isolation count 0 → 4. ⇒ FVSjl NE is a faithful drop-in for multipliers/notriple/thinbta/fixmort.
- ⚠ **DIAGNOSTIC TRAP (re-confirmed the memory note "run_keyfile probes need variant="):** running
  `run_keyfile(nekey)` WITHOUT `variant=Northeast()` defaults to SN ⇒ 5-yr cycles + SN growth model on NE
  data ⇒ *fabricated* gross divergence (looked like a huge NE mortality/height bug; was purely the wrong
  variant). Always pass `variant=` for NE/CS/LS probes AND in the test harness. Cost ~an hour of false alarm.
- **2026-07-06 — CS + LS keyword-isolation batch #1 (4 each, all BIT-EXACT vs live):** cs_/ls_ ×
  {mult, notriple, thinbta, fixmort}, validated with `variant=CentralStates()/LakeStates()`. NE/CS/LS now
  4 keyword-isolation scenarios each in `test_kwcov_variants.jl`. ⇒ multipliers/notriple/thinbta/fixmort
  are faithful drop-ins for ALL FOUR variants.
- **2026-07-06 — batch #2 (SERLCORR / RANNSEED / THINABA × NE/CS/LS, all BIT-EXACT vs live):** notably the
  stochastic-DG paths (SERLCORR ARMA serial correlation, RANNSEED RNG seed) are bit-exact for these stands.
  `test_kwcov_variants.jl` is now AUTO-DISCOVERING (`test/fixtures/kwcov/<prefix>_<kw>.*`, prefix→variant) so
  adding a scenario is just dropping the 3 fixture files in. NE/CS/LS keyword-isolation count now **7 each**.
- **2026-07-06 — batch #3 (THINCC / THINHT × NE/CS/LS, all BIT-EXACT):** count now **9 each** (27 total).
- **2026-07-06 — batch #4 (NOCALIB / DGSTDEV × NE/CS/LS, all BIT-EXACT):** the calibration-disable + DG-
  stddev-bound paths faithful; count now **11 each** (33 total). (THINAUTO dropped — bare keyword is a no-op
  without a companion target; live NE produced no .sum.)
- **2026-07-06 — comparison STRENGTHENED to FULL `.sum` ROW** (was TPA/BA/SDI/CCF/TopHt/QMD only). All 33
  existing scenarios still bit-exact on the volume + mortality columns too. This is now the semantic bar.
- **2026-07-06 — batch #5 (VOLUME / MCDEFECT) SURFACED A REAL GAP:** full-row compare exposed that the
  **R9 volume-override path (NE/CS/LS) diverges on the BdFt/MCuFt columns by ~1–2%** — 5 of 6 scenarios
  (ne_volume, cs_volume, ls_volume, cs_mcdefect, ls_mcdefect) diverge; ne_mcdefect bit-exact. SN (R8) +
  base stand are bit-exact. Tracked `@test_broken` + documented: `docs/MODERNIZATION_R9_VOLUME_OVERRIDE_GAP.md`.
  ⇒ a genuine "100% drop-in" gap for NE/CS/LS — a pillar-1 FIX target (correctness before optimization).
## Semantic map (2026-07-06, full-row vs live) — where NE/CS/LS is bit-exact vs where it diverges
**BIT-EXACT (full-row):** growth/mortality multipliers (BAIMULT/HTGMULT/MORTMULT), ALL thin forms
(THINBTA/ABA/CC/HT/RDEN), NOTRIPLE, FIXMORT, SERLCORR, RANNSEED, NOCALIB, DGSTDEV, STRCLASS, MINHARV,
MCDEFECT(NE), BFDEFECT, and the **BAIMULT+THINBTA combination** (interactions faithful).
**DIVERGENT (tracked @test_broken):**
1. **R9 merch-override board/merch cols** — VOLUME, BFVOLUME (BdFt col 12), MCDEFECT-CS/LS (MCuFt col 10).
   The R9 board/saw CALC under an active override. Task #78 (needs live fvsvol.f stamp). ← main gap.
2. **SDIMAX** — a REAL small self-thin bug (jl over-thins ~2-3% under a lowered max-SDI). **ROOT CAUSE
   CONFIRMED via live SDICAL stamp** (2026-07-06): jl `stand_sdimax` runs ~8 units LOW vs live SDICAL XMAX
   (e.g. jl 437.7 / live 446.0) because it BA-weights over LIVE trees only, while live SDICAL also includes
   the DEAD/being-killed records (the current-cycle BACKGROUND kills, DPROB-scaled) — which jl computes into a
   scalar buffer AFTER stand_sdimax. Lower SDImax → lower self-thin target → over-thin. Fix = match FVS
   mortality ORDERING (book background kills as dead records before the self-thin SDImax); base stays bit-exact
   (inert for uniform SDIDEF). Task #79. `ne_sdimax` @test_broken.

## Net semantic verdict (2026-07-06)
NE/CS/LS are **faithful full-row drop-ins** across 52 keyword-isolation scenarios. Two tracked non-bit-exact
items, both in UNCOMMON keyword paths, both root-investigated:
- **R9 merch-override volume** (VOLUME/BFVOLUME/MCDEFECT): the R9 board CALC is VERIFIED PER-TREE FAITHFUL vs
  a live fvsvol.f stamp; residual ~1% = permitted ULP/boundary class, NOT a bug. (Task #78.)
- **SDIMAX self-thin**: a REAL small bug (~2-3% over-thin), needs a mortality live-stamp. (Task #79.)
The "100% drop-in barring ULP/FVS-bug" claim holds for NE/CS/LS except the SDIMAX-override self-thin (a real,
small, tracked bug) — everything else is bit-exact or verified-faithful-with-ULP-residual.

## 2026-07-06 (later) — MORE coverage + a REAL BUG FOUND & FIXED
Added COMPRESS (eigensolver record compression) + CYCLEAT — both BIT-EXACT vs live NE. And RESETAGE surfaced
a **real reporting bug** (only caught by the full-row bar): post-RESETAGE-to-0, jl MAI = merch/age vs live 0.
FIXED (summary.jl): FVS zeroes MAI when the age is reset to ZERO (evtstv.f `ZERO==0` ⟺ `age_reset_age==0`);
gated so reset-to-nonzero (s17_managed→40) + bare-ground + non-RESETAGE are unchanged. ne_resetage now full-row
bit-exact. ⇒ open real bugs down to ONE (SDIMAX #79); R9 volume verified-faithful (#78); ~57 kwcov scenarios.
Everything else exercised so far is a faithful drop-in.

## 2026-07-07 — R9 volume-override BOARD divergence FIXED (bftopk board-top; #78 resolved)
The last substantive Pillar-1 divergence — the R9 VOLUME/BFVOLUME-override BdFt (~1%) — is FIXED. Root (found via
a live r9logs.f/r9clark.f/r9bdft stamp campaign + a decisive r9clark_cubic unit test that ruled out the mtopp/taper
red-herring): the BROKEN-TOP (tkill) board top-kill in compute_volumes_ne! used the SAWTIMBER scftopd/scfstmp for
`bftopd/bfstmp`, but a blank-SCFTOPD VOLUME card zeroes sp_scf_topd (keeping sp_bf_topd) ⇒ a broken-top board tree's
bftopk shifted. FIX: bftopk uses the board's OWN sp_bf_topd/sp_bf_stump (bf-equal ⇒ inert in the base). Result:
ne_bfvolume FULL-ROW BIT-EXACT (out of _KCV_BROKEN); ne/cs/ls_volume BdFt bit-exact; suite 19486/136/0, no regression.
Remaining volume residual = SawCuFt col off-by-1-2 (0.03-0.05%, genuine ULP-class, present in the base R9 too) —
cornered. ⇒ every substantive Pillar-1 divergence is now FIXED; only ULP-class residuals remain.

## 2026-07-06 (Pillar-1 coverage batch) — CS/LS broadened + recent fixes proven variant-safe
Added 8 NE/CS/LS keyword-isolation scenarios, ALL FULL-ROW BIT-EXACT vs the freshly-relinked live binary
(`{cs,ls}_oracle.sh`), and validated non-vacuous (each card's live output ≠ the no-card run, so the keyword
path is actually exercised — the s3-vacuous-pass lesson):
- **cs_sdimax / ls_sdimax** — SDIMAX on a PRESENT dominant species by NUMERIC sequence index (SM, maxSDI 250
  < 371 default ⇒ self-thin binds). PROVES the #79 SPDECD-sequence-index fix generalizes to CS + LS (bit-exact).
- **cs_resetage / ls_resetage** — RESETAGE-to-ZERO (2010→0), which exercises the #80 MAI-zeroing branch (the one
  that was buggy). PROVES the RESETAGE MAI fix is variant-safe (bit-exact).
- **cs_cycleat / ls_cycleat** and **cs_compress / ls_compress** — close the CS/LS-vs-NE gap (NE already had them).
CS/LS keyword-isolation count 16 → **20** each; NE 23. kwcov gate: **9724 pass / 6 broken / 0 fail** (the 6
broken = the documented R9 volume/mcdefect/bfvolume override set; unchanged). ⇒ the two Pillar-1 fixes landed
this session (SDIMAX, RESETAGE) are now proven faithful drop-ins across all THREE non-SN variants, not just NE.

### Batch 2 (same day) — THINQFA + COMPUTE/event-monitor × NE/CS/LS (6 scenarios, all BIT-EXACT)
- **{ne,cs,ls}_thinqfa** — Q-factor-of-area thinning (THINQFA + its `1` continuation record). Strongly
  NON-VACUOUS (e.g. ne TPA 524→71 at the 2005 thin) and jl reproduces the thin BIT-EXACT ⇒ the THINQFA
  thinning form is a faithful drop-in for all three variants.
- **{ne,cs,ls}_compute** — COMPUTE/END event-monitor block (`MYBA = BBA`). Exercises the variant-independent
  event-monitor parser+evaluator; output-neutral by design (defines a variable, doesn't alter the .sum), so
  it's a "evaluator runs + output unchanged" coverage test. BIT-EXACT.
NE 23→**25**, CS/LS 20→**22** each. kwcov gate: **10738 pass / 6 broken / 0 fail**.

### Batch 3 (same day) — LEAVESP + SPGROUP × NE/CS/LS (6 scenarios; 5 bit-exact, 1 ULP-cornered)
- **{ne,cs,ls}_leavesp** — LEAVESP SM + THINBBA 80 (leave the dominant species while thinning to 80 BA).
  cs/ls BIT-EXACT; **ne_leavesp cornered ULP-class** — ONE cell (col 26, a derived growth/MAI column) at 2010
  differs 36.4 vs 36.3, every density + volume column bit-exact ⇒ sub-integer volume rounding surfacing in a
  finer-rendered column; cs/ls prove the LEAVESP semantic faithful. (Named in `_KCV_BROKEN`.)
- **{ne,cs,ls}_spgroup** — SPGROUP (2-species group SM+WP/WO) + THINDBH group thin (−1 = group 1). BIT-EXACT
  for all three ⇒ species-group definition + group-referencing thin is a faithful drop-in.
NE 25→**27**, CS/LS 22→**24** each. kwcov gate: **11584 pass / 7 broken / 0 fail** (broken 6→7 = ne_leavesp ULP).

### Batch 4 (same day) — FERTILIZ × NE/CS/LS (3 scenarios, all BIT-EXACT)
- **{ne,cs,ls}_fertiliz** — FERTILIZ 2000 200 (N application). Verified NON-VACUOUS (fertilization boosts the
  eastern DG model — live BA differs 2 cycles vs no-fert) and jl reproduces it BIT-EXACT ⇒ the fertilization
  growth response is a faithful drop-in for all three variants. (SPECPREF was attempted but DROPPED: it is
  vacuous with an explicit THINBBA — SPECPREF only steers auto-thin/estab selection, not proportional BA thin —
  so it would not exercise its semantic; per "test must exercise the semantic".)
NE 27→**28**, CS/LS 24→**25** each.

### Batch 5 (same day) — ESTAB/PLANT × NE/CS/LS (3 scenarios; cornered, REAL gap surfaced → task #81)
- **{ne,cs,ls}_estab** — ESTAB + PLANT (300 TPA, 90% survival). Non-vacuous (adds a planted cohort). The
  cohort ESTABLISHES with the correct COUNT (TPA at first post-plant cycle 2010 BIT-EXACT, cs 722/722) but its
  early DIMENSIONS diverge (cs 2010 SDI 236 vs 234), cascading into slightly different self-thinning so TPA
  drifts by 2020+ (cs 609 vs 613). CORNERED @test_broken with the NAMED PRIMITIVE = PLANT/regen planted-tree
  SIZE init (seedling dbh/height + first-cycle growth) for NE/CS; ls_estab is ~2 ULP volume cells (LS
  regen-dimension already done). **Task #81** to fix (mirror the LS establishment-dimension work). This is a
  genuine drop-in gap the isolated harness surfaced — the canonical PLANT bundle passed, but the isolated
  planted-cohort dimensions weren't validated before.
NE 28→**29**, CS/LS 25→**26** each. kwcov gate: **12094 pass / 10 broken / 0 fail** (broken 7→10 = 3 estab).

### Batch 5b (2026-07-07) — cs_estab REAL BUG FIXED (bit-exact) via live cs/regent.f stamps
The cs_estab gap was ROOT-CAUSED to a real jl bug and FIXED. The CS establishment-height REGENT/BALMOD read a
STALE PRE-growth overstory (plot.basal_area/avg_height are refreshed by compute_density! only AFTER establish!,
simulate.jl:463) and used crown_ratio as a wrong proxy for FVS's BA-percentile PCT. Live cs/regent.f
factor-decomposition stamps proved essubh/CON/SCALE/raw-HTCALC all bit-exact — the whole ~1% SDI divergence was
the BALMOD height-modifier inputs (live BA 134.1/AVH 70.2/BAL 134.1 vs jl 109.1/67.6/~76). FIX (CS-localized):
(A) snapshot the POST-growth overstory BA/AVH at establish! entry PRE-seedling (stand_ba/stand_top_height,
mirroring ebau_pre/rmsqd_pre) — a first attempt reading them at phase-2 over-counted the new seedlings and
regressed BARE-GROUND-PLANT by 62 tests (doctrine-#1 revert + refine); (B) bal=ba_e (PCT≈0 for the smallest new
seedling). cs_estab now FULL-ROW BIT-EXACT; bare-ground + all CS canonical bundles green. Suite **19150/138/0**
(cs_estab out of _KCV_BROKEN). Remaining: ne_estab (18 cells, analogous fix for NE's ebau/ne_balmod structure —
task #82), ls_estab (~2 ULP vol cells). ⇒ REFUTED 5 wrong hypotheses en route via live stamps (doctrine #9).

### Batch 5c (2026-07-07) — ne_estab REAL BUG FIXED (bit-exact) via live ne/regent.f stamps
ne_estab root-caused + FIXED — a DIFFERENT bug from CS. Live ne/regent.f decomposition proved the NE establishment
BALMOD (GMOD/RELHT/AVH) all match; the divergence was the RAW height increment: the NE REGENT/LESTB path OMITTED
the CON = exp(htg_cor_small) (=RHCON·exp(HCOR)) factor that regent.f:224 + NE small_tree_growth.jl:48 apply. For WP
CON=0.914 (live HTCALC 17.46·0.914·SCALE 0.5 = 7.98 vs jl's 8.73 without CON). FIX: multiply the NE establishment
htgr by exp(s.calib.htg_cor_small[sp]). Ruled out the avh-timing hypothesis first (inert for NE — doctrine #9).
ne_estab now FULL-ROW BIT-EXACT; net01 canonical bundles GREEN; suite **19318/137/0**. ⇒ BOTH CS + NE establishment
now bit-exact drop-ins; only ls_estab remains (~2 ULP volume cells, not a REGENT-height issue). Live NE oracle intact.
Also applied the SAME CON fix to the LS establishment branch (regent.f:224 / LS small_tree_growth.jl:40) — inert for
JP (CON≈1) so ls_estab stays at its 2 ULP volume cells (off-by-1 at 2030, all density/height cols bit-exact = genuine
ULP), but it closes the latent CON≠1-species bug so all THREE eastern establishment paths (CS/NE/LS) are now faithful.
⇒ **the ESTAB/PLANT establishment-height gap is CLOSED**: CS + NE full-row bit-exact, LS genuine-ULP-cornered. Suite
19318/137/0 throughout; live CS+NE oracles restored+verified bit-exact; no regression.

### Stale-golden audit (same day) — the committed kwcov goldens are TRUSTWORTHY
Spot-checked the 6 OLDEST kwcov goldens (ne_mult, cs_notriple, ls_thinbta, ne_thinaba, cs_serlcorr, ls_nocalib)
by regenerating each vs the FRESHLY-relinked live binary and diffing the committed `.live.sum`: ALL 6 identical
⇒ no golden drift; the harness validates against a faithful baseline (satisfies the goal's "not a stale golden").

## 2026-07-06 (later still) — SDIMAX species-field RESOLVED (real bug FIXED, faithful; #79 closed)
The SDIMAX self-thin "over-thin" (#79) was NOT dead-record mortality ordering (red herring) — it was
WRONG-SPECIES RESOLUTION. FVS SPDECD (spdecd.f:97 `ISP=IFIX(ARRAY)`) reads a NUMERIC species field as the
species SEQUENCE INDEX; jl's kw_sdimax! used `tryparse(Int)` which missed the float format "19.0"/"27.0" and
fell through to resolve_species → mis-mapped to an FIA/OT species. FIX: float-parse then round → push as the
sequence index (keeps −N SPGROUP, 0/blank all-species, alpha-code branches). ne_sdimax now FULL-ROW BIT-EXACT
vs live (moved out of `_KCV_BROKEN`). Confirmed faithful by a live FVS crown.f stamp — FVS itself sets
SDIDEF(27)=500 for HI under the card; jl idx27=HI=FVS NSP#27, SDICON defaults match (276 both).
NEWLY CORNERED: s3_density (SN) — the correct resolution now targets HI (present), exposing a narrow
crown-ACR→growth knife-edge (jl RELSDI bit-exact cyc1-3; SDIDEF=500 inert in live but ~2 TPA/0.4% drift in jl
by cyc4). Tracked `@test_broken` in `_KC_FT_BROKEN` with the full trace (the old s3 "pass" was VACUOUS — the
mis-resolved card was a silent no-op). Suite green 15204/135/0. ⇒ Pillar-1 real open jl bugs now: R9 volume
(#78, verified-faithful ULP) + the s3 crown-coupling residual (small, traced). SDIMAX resolution is FAITHFUL.

- (next) VOLUME/VOLEQ override (needs volume-column comparison, not just TPA/BA), COMPUTE+event-monitor,
  MCDEFECT, THINQFA, STRCLASS, SDIMAX (species-field), NOCALIB/DGSTDEV.
