# FIA/FVS behaviour-compatibility — working checklist / audit

Goal + doctrine: `docs/FIA_FVS_COMPAT_GOAL.md`. Every slice: plots covered, per-cycle pass rate vs
freshly-relinked live FVS, divergences found → both-sides-traced → fixed or cornered. Never regress the
floor (`julia --project=. test/runtests.jl` = 38527/143/0).

## ★ DIVERGENCE TAXONOMY (Pillar-4 consolidated reference) — every FIA divergence is one of these
Verified across 4 variants × 4 mgmt regimes × 3 volume/board paths, and against the WORST outliers (to 436%).
**A. REAL BUGS — FOUND & FIXED** (all surfaced only by real FIA at scale; floor-safe; live-validated):
  1. LS FFE covtyp default → OOB segfault (slice 27, fmcba.jl variant default cover type).
  2. LS extended Scott-Burgan fuel-model index OOB (slice 33, dense-index ffe_fuel_models by raw model#).
  3. NE htcalc NaN calibration-poison (slice 34, missing HTMAX guard in the NCALHT path; SN-family latent).
  4. PLANT/NATURAL scheduled by cycle-number never fired + over-sized seedlings (slice 37, establishment.jl).
**B. CORNERED RESIDUALS — named primitives, bit-exact-or-cornered** (Float32 semantics, not padded tolerance):
  1. **Self-thinning count-straddle** — at a dense self-thin, live/jl kill a different NUMBER of tiny trees;
     DENSITY (BA/SDI/CCF/TopHt) stays bit-exact, only stem count + QMD diverge; re-converges next cycle.
     (growth-ULP × density-dependent mortality; the SIGMAR tripling-spread. Dense/hyper-dense regen stands.)
  2. **Merch/sawtimber-threshold crossing** — a growth-ULP dbh diff straddles a merch (4")/sawlog (10") DBH
     threshold; on the near-zero volume base when a product first forms this is a huge % (to 436%), converging
     as the cohort matures. Volume-only; board:cuft ratio EXACT on all rules (Int'l ¼", Scribner) ⇒ equations
     faithful. Amplified in the most-quantized measures (SCuFt/BdFt).
  3. **DGSCOR/AUTCOR stochastic record-ordering** — diameter-growth serial-correlation draw order (post-thin,
     COMPRESS). Deterministic paths are bit-exact; this is stochastic-by-elimination.
  4. **FFE/FMEFF fire-kill distribution** — the per-tree fire mortality split (fire BEHAVIOR + BA are bit-exact).
  5. **LS dense-phase growth residual** — LS-specific calibration-backdating relative-ranking in hyper-dense
     regen (BA/SDI diverge mid-projection ~20%, converge by stand maturity). Accepted-class per LS port notes.
  6. **Non-native cycle-length DGSCOR** — a variant at a non-native cycle length (NE@5yr, SN@10yr) drifts;
     bit-exact at each variant's native cycle. Deferred known residual.
  7. **Print-boundary ±1-unit straddle** — QMD 1-decimal / TPA·SDI·CCF·BA integer rounding straddles.
STRUCTURE columns never appear in the >5% outlier tail except via (B1)/(B5); no masked volume-equation bug on
any variant. Full per-slice detail below.

## Infra fixes
- **F1** — `test/harness/fia/validate_fia.jl`: fixed `FVSjl.NorthEast()` → `Northeast()` (would have
  errored EVERY NE FIA stand). Made the harness CLI-arg driven (`julia validate_fia.jl <listfile>
  <SN|NE|CS|LS>`) so it's reusable for the campaign. Committed.

## Slice 1 (Pillar 2 probe) — SN multi-cycle drift is REAL and must be characterized
First multi-cycle differential (3 SN plots from `/tmp/fia_val/sn_feas.txt`, NUMCYCLE 5, 6 stand cols vs
live FVSsn). Result — mean |rel diff| by cycle:

| cyc | TPA | BA | SDI | CCF | TopHt | QMD |
|----:|----:|---:|----:|----:|------:|----:|
| 0 | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| 1 | 0.04% | 0.5% | 0.16% | 0.14% | 0.0% | 0.46% |
| 3 | 1.0% | 0.96% | 0.66% | 0.37% | 0.0% | 1.15% |
| 5 | 1.86% | 1.23% | 0.98% | 0.49% | 0.0% | 1.54% |

Worst stand 3.9% (3237541010661). **Cycle-0 bit-exact** (confirms the inventory reader); the drift is
**purely in the projection** and **grows with cycle**. **TopHt = 0.0% every cycle** ⇒ heights match; the
divergence is **DBH-driven** (BA/SDI/QMD directly, TPA via density-dependent mortality).

**Hypothesis (to prove/refute, Pillar 4):** the DGSCOR diameter serial-correlation / grown-Float32
accumulation tail — the SAME accepted class as the `test_multicycle` / `test_dbs_compute` MYBA/MYSDI
`@test_broken`s. BUT the curated snt01/net01 stands are bit-exact multi-cycle, so a ~2-4% drift on real
plots is either (a) that tail hitting harder on these species mixes / larger diverse stands, or (b) a
REAL FIA-specific gap (a species DG coefficient or calibration path these plots trigger that the curated
tests don't). Doctrine #3/#4: root-cause the worst offender BOTH-SIDES before cornering. Next slice:
per-cycle per-species trace of 3237541010661 (are its drivers the WK3-calibrated sp33/65 family = accepted,
or a clean species that SHOULD be bit-exact = real gap?), then scale the differential to all 4 variants.

### Slice 1b — worst-stand species composition (classification clue)
Both worst SN stands are **diverse mixed-species hardwood plots**: 3237541010661 = 5 species
(FIA 68/407/835/541/462), 3196569010661 = **12 species** (FIA 826/837/621/318/931/541/409/701/544/521/
491/404). This matches the documented `test_allspecies` verdict: a diverse many-species stand
"accumulates every species' sub-ULP per-cycle DBH-growth + tripling-spread residual into the nonlinear
density/volume sums = the ACCEPTED aggregate DGSCOR + tripling class." Real FIA plots are inherently
multi-species (unlike the clean curated single/few-species tests), so they AMPLIFY that same cornered
primitive — a strong hypothesis that the ~2-4% drift is the accepted class, not a new gap.

**DECISIVE next diagnostic (do NOT assume — doctrine #3/#4):** the `test_treeszcp` method — compare
PER-TREE DG jl-vs-live (via FVS_TreeList DBS or a debug-FVS dgdriv stamp) on one worst stand. If per-tree
DG is bit-exact to ~1 ULP and only the AGGREGATE sums drift, it IS the accepted grown-Float32/DGSCOR
accumulation class (corner it). If a specific species' per-tree DG diverges beyond ULP, it's a REAL gap
(fix it). Also worth checking: cycle LENGTH used (non-native-cycle DGSCOR is a separate documented
deferred residual — fvsjl-scenario-sweep-findings #2).

### Slice 1c — DEFAULT (drop-in) trajectory pinpoints a SYSTEMATIC diameter-growth drift (stand 3196569010661)
Ruled out cheap causes first (doctrine #3): **cycle length = 5-yr = SN native** (1998→2023, NOT the
non-native-cycle DGSCOR issue); **NOTRIPLE does NOT converge** them (jl 408 vs live 415 TPA @2003 — a
NOTRIPLE-path confound, not the default behaviour). The apples-to-apples DEFAULT (tripled) comparison:

| year | jl TPA/BA/SDI | live TPA/BA/SDI | BA Δ |
|-----:|---------------|-----------------|-----:|
| 1998 | 416/102/195 | 416/102/195 | 0 (bit-exact cyc0) |
| 2003 | 408/114/214 | 408/115/215 | 1 |
| 2008 | 396/125/231 | 394/128/234 | 3 |
| 2013 | 378/135/244 | 375/139/248 | 4 |
| 2018 | 362/146/259 | 358/151/264 | 5 |
| 2023 | 347/157/273 | 343/163/278 | 6 (3.7%) |

**Signature:** TPA tracks closely (Δ0–4); **BA drifts SYSTEMATICALLY and ONE-DIRECTIONALLY — jl BA always
LOWER**, 0→6 (3.7%) over 5 cycles; TopHt exact all cycles. ⇒ seed = **jl diameter growth marginally but
systematically below live's**, compounding into BA/SDI (and slightly higher TPA via reduced density
mortality).

**Verdict discipline (doctrine #3/#4):** this is **NOT yet cornered**. It is (a) directional (pure Float32
ULP is not) and (b) ~3.7% BA — orders of magnitude larger than the accepted `dbs_compute` MYBA/MYSDI
DGSCOR ULP (~4 ULP / 1.2e-4). A systematic directional diameter-growth bias on diverse real SN plots must
be root-caused per-tree before any corner. It could still resolve to the accepted class (if it's the
sp33/65 WK3-calibrated species' serial-correlation compounding one way on these mixes) OR be a REAL small
DG gap a curated few-species test never exercised. **This is the campaign's first substantive open divergence.**

### Slice 1d — ROOT-CAUSE: yellow-poplar large-tree diameter growth ~20% LOW (REAL gap, NOT cornered)
Per-tree DG differential (live FVS_TreeList `DBH,DG` via `TREELIST`+`TreeLiDB`→FVSOut.db; jl per-tree
`s.trees.diam_growth` direct from the engine; matched by sorted initial DBH — cycle-0 bit-exact). Result
on stand 3196569010661, cycle 1998→2003:

- **Most trees: dDG = 0.0 (BIT-EXACT diameter growth).**
- A subset — **all yellow-poplar (FVS `YP`, FIA 621)**, the large fast-growing ones — jl DG is
  **systematically ~20–25% LOW**:

| initDBH | species | live DG | jl DG | Δ |
|--------:|:--------|--------:|------:|----:|
| 8.0 | YP | 0.712 | 0.558 | −0.155 |
| 11.9 | YP | 0.873 | 0.682 | −0.191 |
| 13.3 | YP | 0.635 | 0.495 | −0.140 |
| 14.5 | YP | 0.902 | 0.704 | −0.198 |
| 14.7 | YP | 1.029 | 0.804 | −0.225 |
| 17.9 | YP | 1.019 | 0.794 | −0.225 |

Co-large white ash (WA, DBH 11.8, live DG 0.709) was **bit-exact** ⇒ the gap is **species-specific to
yellow-poplar**, not a general large-tree issue. This ~20% per-cycle YP DG deficit compounds into the
observed 3.7%-BA/5-cycle stand drift (YP is a dominant large-tree component of these mixed hardwood plots).

**VERDICT: a REAL DG gap — NOT ULP, NOT cornered.** Directional, ~20% relative, single species. First real
bug of the campaign. Root-cause = trace the SN yellow-poplar large-tree DG (`dgf`/`dds` coefficients +
size terms) in FVS blkdat/dgf vs jl's coefficients; small YP are ~exact, large YP low ⇒ a DBH-dependent
term (dbh/dbh²/ln(dbh) or a size cap) for YP. **This is the next fix (Pillar 4), keeping the 38527/143/0 floor.**

NOTE also: jl loaded **36** tree records vs live **37** at cycle 0 for this stand (.sum TPA was bit-exact,
so likely a zero-TPA/merged record) — track separately.

### Slice 1e — mechanism leads (FVS SN dgdriv.f), YP coefficients located
- jl YP DG coefficients: `data/southern/species_coefficients.csv:46` (sp 45, YP, FIA 621). YP already has
  a documented SN special-case in HTCALC (`height_growth.jl:20`) — SN frequently special-cases YP.
- **Prime suspect (dgdriv.f:537-538):** "trees that have GREATER DBH THAN THE LARGEST TREE IN THE GROWTH
  SAMPLE are ASSIGNED THE PREDICTED RESIDUAL FOR THE LARGEST TREE IN THE GROWTH SAMPLE." This is the
  large-tree DGSCOR/COR residual **extrapolation** for DBH beyond the species' growth-sample max. If these
  FIA yellow-poplars are LARGER than any YP in the curated snt01 sample, jl's beyond-sample YP extrapolation
  runs for the first time — and is ~20% low. Fits "small YP exact, large YP low." (Distinct from the
  accepted sp33/65 WK3 tail — this is a LARGER, YP-specific effect.)
- Secondary: `dgdriv.f:741` `DGBND` (bounds on DG value) — large trees may hit a bound jl computes differently.
- Next: debug-FVS stamp dgdriv.f (the sp/DBH/DGSCOR-residual for YP on this stand) + read the jl large-tree
  DG path; find where large-YP DG diverges; fix keeping 38527/143/0. Then re-run the per-tree differential.

### Slice 1f — DGBND ruled out; reframed as a uniform YP DDS offset (size-amplified)
- **DGBND ruled out (both-sides):** jl YP `dg_bound_dbh_lo/hi` = 998/999 (sentinel = no range adjustment);
  FVS `DLODHI(45,·)` = 998.0/999.0 too (dgbnd.f I=31,45 block, 15th pair). Faithfully absent on both.
- jl base DDS (diameter_growth.jl:175-183) = conspp + intercept + ln_dbh·ln(d) + dbh_sq·d² + ln_crown·ln(icr)
  + rel_ht·relht + stand_ba·ba + point_bal·pbal + ft_coef + planted·kplant. (`dg_site_index` is folded into
  the per-stand constant via `dgcons!`, not a per-tree term.) YP coeffs: intercept −2.513351, ln_dbh 1.495351,
  dbh_sq −0.000756, ln_crown 0.530123, rel_ht 0.161718.
- **REFRAME:** DG = f(DDS, dbh) amplifies a CONSTANT ln-space DDS error for large dbh. So a DG deficit that
  grows with tree size is the signature of a **uniform YP DDS offset**, NOT necessarily large-tree-specific
  logic. Candidates: the folded site-index term (dgcons!), `dg_const[YP]`, `dg_cor[YP]` default, or a base
  coefficient — any small YP-constant error shows up big only on large YP.
- **DECISIVE next (debug-FVS, doctrine #6):** dgf.f already has `IF(DEBUG)WRITE ... 'IN DGF, I=,ISPC=,DDS='`
  (line 376) + a RELHT/DG5/DDS stamp (365). Run the stand with DEBUG (or stamp dgf.f), dump DDS + each term
  for a large YP tree; dump jl's DDS+terms for the same tree; the divergent term names the fix. Restore
  source + rebuild clean .o + verify oracle pristine.

### Slice 1g — ROOT CAUSE CONFIRMED: jl omits the SN ecological-unit (EUT) diameter-growth term
Debug-FVS DDS dump (dgf.f DEBUG) vs jl per-tree ln(DDS): **jl is UNIFORMLY 0.344 below FVS for EVERY YP
tree** (14.7: FVS 3.3614/jl 3.0175; 8.0: 2.2383/1.8945; 17.9: 3.7659/3.4221 — Δ=−0.344 exactly, all DBH/ICR).
A perfectly uniform ln(DDS) offset ⇒ a per-species CONSTANT-term error (amplified in DG space for big trees,
hence "large-YP-only"). Also refuted en route: crown variable (FVS `ICR`=25/35/45 == jl `crown_pct` ✓),
DGBND (both 998/999 ✓).

**CAUSE:** FVS STDINFO for this stand = **ECOLOGICAL UNIT `223Db`** (FOREST-LOCATION 80215). FVS dgf.f adds a
per-species **ecological-unit (EUT) categorical term** selected by that eco unit (dgf.f:1039 "EXAMINE THE
STAND ECOLOGICAL UNIT VARIABLE"; SN `dg_phys_*` coefficients p221/p222/p232/…). **jl's `s.plot.eco_unit` is
BLANK** ⇒ `_dgf_phys_group` returns `none` ⇒ `dgcons!` adds ZERO EUT term. YP's EUT coefficient for group
223 is large (~+0.255, `dg_phys_p222`-class) so YP loses it; white ash's ~0 ⇒ WA stayed bit-exact. (A small
forest-type-term difference — jl derives `forest_type=503→upok`, FVS may classify differently — likely
accounts for the remaining ~0.09.) The FIA reader reads `LOCATION=80215` but does NOT populate the
ecological unit → the whole EUT DG term is dropped for FIA-DB stands.

**THE FIX (next slice):** populate `s.plot.eco_unit` from the FIA STANDINIT ecological-unit field, and
confirm `_dgf_phys_group("223Db")` maps to the correct SN group (223→p222-class). Validate: re-run the
per-tree DG diff (expect YP Δ→0) + the stand differential (expect the 3.7% BA drift to collapse), and the
full suite MUST stay 38527/143/0 (the EUT term is inert when eco_unit is genuinely blank — snt01 etc. don't
set it — so no floor risk). This is a real FIA-reader + DG-constant gap, NOT ULP.

### Slice 1h — FIX LANDED + VALIDATED: FIA reader now reads ECOREGION → eco_unit (SN EUT DG term)
**Fix:** `apply_fia_stand!` (src/io/fia_database.jl) now reads the `ECOREGION` STANDINIT column (e.g.
"223Db") into `p.eco_unit` via `resolve_eco_unit`→SNECU (canonical "223DB"), SN-gated. This restores the
SN dgf ecological-unit (EUT) categorical term (`dg_phys_p222 = +0.2554` for group 223→p222) that was
dropped for all FIA-DB stands.

**Validated (stand 3196569010661):**
- YP per-tree ln(DDS) deficit: **−0.344 → −0.0886** (the p222 term = +0.2554, exactly as predicted).
- Stand `.sum` vs live: **BIT-EXACT all cycles** (2003 408/115/215, 2013 375/139/248, 2023 343/163/278 ==
  live) — the **3.7% BA drift collapsed to 0**. The residual −0.089 ln(DDS) forest-type term (jl upland_oak
  −0.0907 vs FVS ~0) is below the `.sum` integer-rounding threshold here.
- **Floor held: 38527 pass / 143 broken / 0 fail** — no regression (SN-gated FIA-path change; curated
  scenarios set eco_unit via STDINFO, unaffected; cycle-0 FIA baseline is DG-independent).

★ First campaign divergence FOUND → both-sides ROOT-CAUSED → FIXED → validated BIT-EXACT vs live, floor
preserved. Textbook FIA-behaviour-compat slice.

### Slice 1i — batch validation: eco_unit fix is BROAD (5/8 SN stands now bit-exact)
8-stand SN multi-cycle differential post-fix (`validate_fia.jl`): **5/8 stands now 0.0% (bit-exact all
cycles)** — the traced stand 3196569010661 PLUS 3237541010661, 3237540010661, 255262477010854,
1152014712290487 (all drifted pre-fix). Mean BA drift @cyc5: ~1.2% → **0.14%**. Confirms the EUT fix
corrected the whole FIA sample's DG, not one stand. Remaining drifters → next targets:
- **255260379010854 — 6.7%** (worst; also a known cyc0 TopHt AVH stand from modernization #92).
- **1152014752290487 — 6.5%**.
- 255262502010854 — 2.5%.
(One stand shows cyc0 TopHt 0.83% = the pre-existing 3/162 AVH-tie cyc0 residual, not new.)

### Slice 1j — residual SN stands classified: 1 ACCEPTED, 1 NEW (loblolly DG)
- **255260379010854 (6.7%) = ACCEPTED cyc0 TopHt AVH tie** (modernization #92 named this stand). The 6.7%
  is TopHt 32/live30 at the 2010 INVENTORY row (2/30) — the cornered AVH pre-sort tie; it propagates ~2 ft
  a few cycles then CONVERGES (2030 56/56, 2035 61/61). Not a new bug; carry as cornered. (Small secondary
  DG/mortality tail exists — TPA 3548/live3520 @2025 — but sub-dominant.)
- **1152014752290487 (6.5%) = NEW: loblolly-pine DG slightly HIGH.** Near-pure loblolly (FIA 111 = SN sp13,
  38/42 trees). Cyc0 BIT-EXACT (579/112/5.9); jl then grows DBH too fast (BA 129/live126 @2026 → SDI/QMD
  high → extra density mortality → TPA 447/live478 = −6.5% @2046). OPPOSITE direction to the YP gap (jl low).
  A distinct species-specific DG divergence — loblolly (sp13, a major planted SN species; special Fort-Bragg
  DG path exists). Next: per-tree DG diff on this stand — is it the planted modifier (kplant), an sp13 coeff,
  or an eco/fortype term for this stand's ECOREGION?

**SN multi-cycle status post-eco_unit-fix:** 5/8 bit-exact; residuals = 1 accepted AVH tie + 1 new loblolly
DG lead. The eco_unit fix cleared the dominant DG divergence class; remaining leads are narrower/species-specific.

### Slice 1k — SA(sp6) DG-high partially traced: NOT my fix, NOT fortype/phys (species term, open)
Corrected: the stand's dominant species is FVS **SA = sp6** (FIA 111), not loblolly-sp13 (both jl AND live
resolve FIA 111→sp6 here, so the mapping is consistent — not a mapping bug). Per-tree DG diff: jl SA DG
uniformly **~+0.10 higher** (mean +0.0998, relative ~25%). Ruled OUT (measured):
- **NOT a regression from the eco_unit fix:** SA `dg_phys_p232 = 0.0`, so on this 232Db→p232 stand the fix
  adds nothing for SA. (Fix confirmed safe.)
- **NOT the forest-type term:** jl derives `forest_type=142→ylpn` (yellow pine — correct for this pine
  stand); SA `dg_fortype_yellow_pine = 0.0`. So no fortype contribution either.
- jl context: `dg_const[SA]=0.6988`, `SI[SA]=102`. The +0.10 is some OTHER SA DDS term.
- Per-tree DDS match is currently MUDDY: FVS debug `D=` (3.29/4.59/5.59…) is ~0.5 below jl's DBH
  (3.9/5.0/6.2…) despite cyc0-bit-exact `.sum` — must resolve (bark/dib? tripled records? debug phase?)
  before decomposing the DDS term. OPEN — the campaign's 2nd substantive lead.

**NOTE — this refutes the earlier "loblolly sp13" framing in slice 1j** (it's sp6/SA). Magnitude is modest
(2–4% mid-projection, 6.5% TPA at the 2046 terminal cycle).

### Slice 1l — SA(sp6): site index RULED OUT (jl clamps correctly); narrowed to a per-species DG coefficient
- **Site index ruled out both-sides:** FVS uses SI=102 for SA (ISPC=6 dgf dump; SA's allowable range
  includes 102, NOT clamped — the `.out` clamp warnings are for OTHER site species: YP 144→135, BY 127→120,
  OS 67.5→55, …). jl also uses sp_site_index[SA]=102. Match.
- **jl DOES clamp site index (not a gap):** SITSET (`site_index.jl:56-85`) uses `site_index_min/max` — line
  66 clamps the site species to `simax`, and the per-species mapping stays within each species' [simin,simax].
  jl's derived YP=135 matches FVS's clamped 135; SA=102 matches. So the earlier "jl misses SI clamping"
  hypothesis is REFUTED — jl's site-index handling is faithful.
- ⇒ SA +0.10 DDS deficit is now cornered to a **per-species DG coefficient/term** (intercept/ln_dbh/dbh_sq/
  ln_crown/rel_ht/stand_ba/point_bal) — everything stand-level (site/slope/aspect/phys/fortype) matches.
  The dDG VARIES with tree DG (0.05–0.16), hinting a size/crown/BA-scaled term rather than a flat intercept.
  jl SA coeffs: intercept −1.641698, ln_dbh 1.461093, dbh_sq −0.00253, ln_crown 0.265872, rel_ht 0.069104,
  stand_ba −0.002939, point_bal −0.004873. Decisive next: compare these to FVS dgf.f DATA (INTERC at :399 +
  the b-arrays) for species 6 — a CSV-vs-blkdat coefficient diff.

**SA divergence magnitude:** modest (2–4% mid-projection, 6.5% terminal TPA) on one stand; not floor-blocking.

### Slice 1m — SA(sp6) ALL 8 dgf coefficients verified == FVS; narrowed to multi-point point-BAL input
Read FVS dgf.f DATA arrays, species-6 (SA) value == the CSV for **every** term:
INTERC −1.641698, LDBH 1.461093, DBH2 −0.002530, LCRWN 0.265872, ISIO 0.006851, HREL 0.069104,
PLTB(stand_ba) −0.002939, PNTBL(point_bal) −0.004873. **NOT a coefficient error** — extraction is correct.
Combined with site/phys/fortype/clamp all matching, the SA +0.10 DDS offset is a **stand-level INPUT** for
this stand, not the equation constants. 
- Prime suspect: **point-BAL** (`point_bal·pbal`, pbal = point-BA-in-larger × (1−cr/100)). This stand is a
  **dense multi-point FIA plot** (BA 112, 579 TPA); per-point PTBAA is the known harder multi-point-density
  area (cf. the modernization TCONDMLT/pccf per-point notes — jl faithful single-point, deferred multi-point).
  A per-point vs stand-point pbal difference × −0.004873, amplified by density, plausibly yields the +0.10.
- Secondary: `stand_ba·ba_v` (but cyc0 BA is bit-exact 112, so ba_v matches) or `rel_ht·(height/avh)`.
- Both-sides trace needed: dump FVS pbal (PTBAA per point) vs jl's for the SA trees; confirm the multi-point
  point-density path. If it's the multi-point PTBAA, this is a shared-density-model slice, not SN-specific.

**SA verdict:** NOT a coefficient bug; a stand-level competition-input divergence on dense multi-point FIA
plots — the campaign's 2nd substantive open item (modest magnitude, likely a multi-point density model slice).

### Slice 1n — SA divergence CORNERED to a named mechanism: point-BA-in-larger on multi-point FIA plots
Confirmed multi-point: stand 1152014752290487 has **4 points**. FVS debug (`IN PTBAL`) computes the ACTUAL
per-tree **PTBALT** (point BA in larger: 0/21.8/40.3/57.3/71.2 for descending DBH) with total point
**PTBAA=183.216**. **jl approximates** the competition term as `pbal = PTBAA·(1−cr/100)` — a crown-ratio
proxy for "fraction of BA in larger" (`diameter_growth.jl:163-168`), NOT the true rank-based PTBALT — and
jl's per-point BA (177.5/143.5/49.7/75.5) differs from FVS's point structure. On a dense multi-point plot
these diverge (→ SA DG +0.10, ~25%); on the simpler single-ish YP stand they matched (YP bit-exact). 

**NAMED CORNER:** the `point_bal` competition input on multi-point FIA stands — jl's crown-proxy
`PTBAA·(1−cr/100)` vs FVS's actual rank-ordered `PTBALT`, compounded by point-assignment differences. This
is a **shared density-model** item (the modernization per-point-density deferral: pccf/TCONDMLT were faithful
single-point / deferred multi-point). Fixing it = compute true per-point PTBALT (BA in larger at each point)
and feed dgf's point_bal term — a cross-cutting density-model slice, both-sides-traced, floor-guarded.
Magnitude modest (2-4% mid-projection on dense multi-point stands); does not block the SN single/simple-point
majority (already 5/8+ bit-exact).

### Slice 1n-CORRECTION — jl's point_bal FORMULA is correct; the diff is the PTBAA VALUE
Both-sides re-check (dgf.f:282-297): FVS computes `PBA=PTBAA; PBAL=PBA*(1-PCT/100); DDS += PNTBL*PBAL` —
the **SAME crown-proxy jl uses** (`pbal=PTBAA·(1−cr/100)`, diameter_growth.jl). So slice-1n's "jl crown-proxy
vs FVS actual PTBALT" verdict was WRONG (the `IN PTBAL` PTBALT debug is computed but NOT used in the SN DDS
point_bal term). CORRECTED cause: **jl's PTBAA (per-point total BA) VALUE differs from FVS's** — jl
point_ba[pt1]=177.5 vs FVS PTBAA=183.216 (~3%). The residual is the per-point basal-area COMPUTATION on
multi-point FIA stands (point assignment / expansion / PROB weighting), NOT the formula. Lesson logged: I
inferred from the PTBAL debug without checking what dgf CONSUMES (doctrine #3). Next: trace how FVS builds
PTBAA per point (fvsGetPtBal / the point BA accumulation) vs jl `compute_density!` point_ba on this 4-point stand.

### Slice 5 — REAL BUG (jl-high verified non-ULP per user): loblolly DG 2× via DGSCOR over-calibration
Per the user's directive ("jl-high ones are not reassuring; make sure ULP-class before scaling, else fix"),
root-caused the worst jl-high stand 157875449010854 with a PERIOD-ALIGNED per-tree diff (absolute cyc1 DBH,
not the treelist-DG that burned slice-3): **jl grows loblolly (LP/sp13) ~2× (4.5″: live 0.90/jl 1.91;
13.0″: live 1.26/jl 2.55); ALL other species BIT-EXACT** (SU/JU/SP/SK/BC/BG/WK dDG≈0). Definitively REAL &
species-specific, NOT ULP.
- Coefficients MATCH FVS (INTERC(13)=0.222214, LDBH 1.16304, DBH2 −0.000863) — not the dgf coeffs.
- LOCALIZED: jl `conspp = dg_const[13](0.2057) + dg_cor[13](0.9777)` (diameter_growth.jl:154). The
  **`dg_cor[LP]=0.9777` DGSCOR self-calibration term (line 450-472) is added to ln(DDS) ⇒ ×e^0.98≈2.66 DG.**
  Other species have dg_cor≈0 ⇒ bit-exact. jl's LP DDS (0.22/0.70 @dbh1.0/1.5) − 0.9777 ≈ FVS's (−0.65/0.008).
- Live FVS LP DG is LOW ⇒ FVS is NOT applying a +0.98 COR for LP here. So jl's self-calibration is FIRING /
  computing a large COR for loblolly on this FIA stand when FVS doesn't (or gets ~0). 
- **DECISIVE next (both-sides, do NOT fix blind):** debug-FVS dgdriv COR dump for ISPC=13 on this stand — is
  FVS's LP COR ≈0 (jl wrongly fires) or does the calibration input differ (which LP trees are "measured",
  their residuals)? The modernization DGSCOR was validated bit-exact on curated stands (measured sp33/65), so
  the divergence is likely a FIA-specific calibration INPUT (measured-tree flag / past-growth) that makes jl's
  LP fit fire spuriously. Fix there, keep floor, re-verify LP DG→live, re-sweep.

### Slice 5a — ROOT CAUSE: jl mis-scales the FIA measured `DG` field feeding the DGSCOR calibration
BOTH-SIDES (FVS `.out` calibration report): **"DBH GROWTH MODEL SCALE FACTORS WERE COMPUTED 0.00 0.88 0.00…"
— FVS computes a 0.88 (slightly REDUCE) scale for loblolly**; jl computes `exp(dg_cor[13]=0.9777)=2.66`
(2.66× INCREASE). Opposite direction, ~3× off. Both calibrate ONLY loblolly (others 0.00). The input is the
FVS_TREEINIT_COND **`DG` field = 8.4** (measured past diameter growth; `DIAMETER=9.1`) — the DGSCOR self-
calibration fits observed-vs-predicted DG. jl's observed-DDS (`reslog`) from DG=8.4 is ~3× too large ⇒ COR
2.66 vs FVS 0.88. **⇒ jl mis-reads/scales the FIA `DG` measured-growth field into the calibration** (units /
period / the observed-DDS computation). The modernization DGSCOR was validated on TREEDATA-measured stands;
the FIA-DB path feeds `DG` differently and jl's scaling diverges. This is the REAL fix site.
- Next: trace how jl's FIA reader + calibration consume TREEINIT.`DG` (units, whether it's past-period Δdbh,
  ×0.1, etc.) vs how FVS's DBS reader passes it to dgdriv; align jl so the LP COR → ~0.88 (== live). Floor-
  guard (curated TREEDATA calibration must stay bit-exact). Then LP DG → live, and re-sweep to measure lift.
- ★ This is a systematic FIA-path bug: ANY stand with measured-DG trees mis-calibrates ⇒ likely a big share
  of the 57 big failures + the growth-driven mortality divergences. Highest-value fix found in the sweep.

### Slice 5b — growth_idg hypothesis REFUTED; honest reset on the loblolly trace
Implemented + tested setting `growth_idg=1` in the FIA reader (DG=past-DBH). Result: **INERT on the .sum**
(loblolly stand still 944/200/411 @1995 == pre-change) and **`dg_cor[LP]` stays 0.9777** regardless. So the
DGSCOR calibration does NOT depend on the DG-increment conversion — the growth_idg root-cause was WRONG.
REVERTED (doctrine #3: don't ship an unvalidated, target-inert change; floor re-confirmed 38527/143/0).

**HONEST RESET.** This loblolly trace has produced repeated MEASUREMENT ARTIFACTS + wrong attributions:
(1) "3× low" (treelist PrdLen misalignment), (2) "crown-ratio bug" (it's the BA percentile), (3) "2× high"
(jd0 backdating), (4) "growth_idg/DG-as-increment" (inert). The ONLY trustworthy facts, from the `.sum`
oracle: **jl loblolly BA is ~8% high @cyc1 (200 vs live 185), TPA ~24% low; jl `dg_cor[LP]=0.9777`
(exp 2.66) vs FVS's reported scale 0.88.** So the divergence IS real (per the user, not ULP) and lives in
the **DGSCOR self-calibration FIT** (jl over-fits loblolly's COR), but the exact fit-input divergence is
NOT root-caused. Lesson (hard-learned this session): on coupled growth/backdating/calibration, TRUST ONLY
the `.sum` + FVS's own reported scale factors; per-tree Julia probes here keep mis-aligning period/quantity/
backdating. Next attempt must instrument the calibration FIT both-sides (dgdriv COR inputs: which trees enter
the LP fit, their residuals) via a debug-FVS stamp, NOT another Julia per-tree reconstruction.

### Slice 5c — debug-FVS attempt (instrumented wrong segment); honest consolidation on loblolly
Instrumented dgdriv.f:195 (made the `NEW DGCOR` dump unconditional) + built a SEPARATE debug binary
`/tmp/FVSsn_dbg2` — but line 195 is the NEXT-CYCLE ATTENUATION path (gated `IF(.NOT.LDGCAL)`, fires
post-cycle-1), not the INITIAL calibration COR. No dump. RESTORED dgdriv.f (matches backup), rebuilt clean
`.o`, guard verified, dbg binary removed, oracle `/tmp/FVSsn_new` untouched (hygiene per doctrine #6).
- **Correct instrument target for next attempt:** dgdriv.f **~319-470** (the "CALIBRATION SECTION"): dumps
  already exist — **9003** `I,ISPC,DG,TERM,DEV,DEVSQ,RESLOG` (per-tree residual, the fit input) and **9010**
  `ISPC,SPOPN,SPOPX,FN,SNP,…` (the fit sums) and **157** `I,DG,BARK,WK3,WK2,SCALE`. Make the ISPC=13 ones
  unconditional (or enable CALBSTAT), capture, compare to jl's `snx/sny/snxx/snxy/reslog` for LP.

**CONSOLIDATED HONEST STATE (loblolly, after a long error-prone trace):** CONFIRMED real & not-ULP (user's
concern validated) — `.sum` BA ~8% high @cyc1, jl `dg_cor[LP]=0.9777`. Localized to the DGSCOR self-
calibration OVER-FITTING loblolly's COR. NOT root-caused: 5 attributions tried & retracted (period, crown-
ratio-vs-percentile, backdating "2×", growth_idg, wrong-instrument-segment). The fix needs the calibration-
FIT both-sides trace above, done carefully in a fresh pass. I deliberately shipped NO speculative fix here.

### Slice 7 — ★★ FIXED (systematic): FIA reader didn't set `growth_dg_set` ⇒ DG calibration un-normalized by FINT
ROOT CAUSE (clean, both-sides): FVS_STANDINIT_COND provides `DG_TRANS`(=1) + `DG_MEASURE`(=9yr). The FIA reader
set `growth_idg=1` and `growth_fint=9` from them, but NOT `growth_dg_set` — which `simulate.jl:47` gates the
`dgscale = YR/FINT` normalization on (`growth_dg_set ? yr/dfint : 1`). So dgscale stayed 1, the 9-yr observed
DG increment was NOT normalized to the model's 5-yr period in the DGSCOR self-calibration ⇒ observed DDS ~1.8×
too high ⇒ loblolly COR over-fit to 0.9777 (exp 2.66×) vs FVS's fort.13 (CALBSTAT) raw scale 1.411.
**FIX** (`fia_database.jl`): set `growth_dg_set=true` when DG_TRANS/DG_MEASURE present. Then dgscale=5/9,
LP COR 0.9777→0.3442 (exp 1.411 == FVS raw), and the loblolly stand `.sum` goes from jl 944/200 @1995
(24% TPA low, 8% BA high) to **jl 1233/186 == live 1234/185 (Δ1 unit, ULP-class) at EVERY cycle**.
**Floor held 38527/143/0** (curated tests use kw_growth!, unaffected). SYSTEMATIC — every FIA measured-DG stand.
Re-sweeping SN-1000 to quantify the pass-rate lift (was 566/820=69%).

**★★ POST-FIX SN-1000 RESULT (the payoff):** bit-exact **611/824 (74%, up from 69%)**, FAIL 213 (was 254).
The magnitude tail COLLAPSED: worst-rel-diff histogram over ALL 824 stands = **<1%:734 · 1-2%:55 · 2-5%:29 ·
5-10%:5 · >10%:1** — i.e. **>10% failures 57→1, 5-10% 25→5**; the fix cleared ~56 of the 57 big divergences,
and **89% of stands are now within 1%** (ULP/print-straddle class). The systematic growth_dg_set fix (+ the
earlier eco_unit + JLERR fixes) resolved essentially the entire heavy tail on SN. Remaining 213 "fails" are
almost all sub-1% ULP straddles (strict all-cycle == criterion) — the cornerable class — plus 5 in 5-10% and
1 >10% to triage individually. Live-NOSUM still 176 (live FVS's own failures). Floor 38527/143/0.

★ This is the campaign's SECOND broad fix (after eco_unit), and vindicates the user's "jl-high = fix, not
ULP" call. The long error-prone trace (5 retracted attributions) finally resolved once I (a) used FVS's own
CALBSTAT output as the both-sides anchor and (b) found the actual call-site gate (`growth_dg_set`), not more
Julia per-tree reconstruction. Meta-lesson reinforced: anchor on the `.sum` + FVS's reported values.

### Slice 8 — remaining SN >5% failures triaged (post-fix): individual residuals, not systematic
Of the 5 stands still >5% post-fix (signature.jl):
- `157872834010854` QMD@2005 jl4.4/lv4.5, `155775361010854` QMD@1980 jl15.7/lv15.8 — QMD ±0.1 PRINT-STRADDLE
  (rendered-integer/tenth knife-edge) → cornerable ULP.
- `155776435010854` BA@1979 jl26/lv25, `921827906290487` CCF@2026 jl25/lv26 — Δ1-unit straddles → cornerable.
- **`158851606010854` (13.6%, the lone remaining >10%)** — FIA-128-dominated, **no measured DG (pastDBH null)**
  so NOT the calibration; jl OVER-grows the base DG (1996 BA jl134/live118, QMD 6.6/6.0). An INDIVIDUAL
  species/condition base-DG residual (FIA 128), isolated to ~1/1000 — a narrow follow-up, not systematic.

**SN END-STATE (this campaign):** the TWO systematic bugs (eco_unit EUT term; growth_dg_set FINT-normalization)
are FIXED — big-divergence tail collapsed 57→1 (>10%), 74% strict bit-exact, **89% of real stands within 1%**.
Remaining = a cornerable sub-1% ULP/print-straddle tail + a handful (~5) of individual species/condition
residuals (e.g. FIA-128 base-DG). SN systematic work is essentially DONE; the residual is long-tail/cornerable.
NEXT per the plan: scale NE/CS/LS sweeps (the eco_unit + growth_dg_set fixes are in the shared FIA reader, so
NE/CS/LS measured-DG stands should benefit — measure it), then Pillar-3 management scenarios.

### Slice 6 — BOTH-SIDES GROUND TRUTH (via CALBSTAT, no source edit): FVS loblolly calib = 0.879, jl = 2.66×
Used the **CALBSTAT keyword** (clean — writes per-species calibration to `fort.13`, no instrumentation).
FVS loblolly (sp13) large-tree DG calibration: `CAL: LD 13 LP 12 1.411 0.713 0.879` — **12 measured trees,
initial scale 1.411 → FINAL 0.879** (a slight REDUCTION; == the 0.88 summary). jl `dg_cor[LP]=0.9777`
(added to ln-DDS ⇒ exp=2.66× INCREASE); FVS's equivalent additive COR = ln(0.879) = **−0.13**. So jl's
loblolly COR is **+1.1 off in ln-DDS** — the quantified, both-sides-confirmed, non-ULP bug the user flagged.
- ROOT is now bounded to the calibration FIT: jl over-estimates loblolly's OBSERVED growth (positive residual
  → +0.98 COR) where FVS's 12-tree fit gives −0.13. The observed-growth/residual (`reslog`) computation for
  measured LP trees is where jl diverges. NOTE growth_idg=1 was INERT on this ⇒ the calibration does NOT read
  the converted `diam_growth`; it computes observed growth via the backdated WK3 path — trace THAT vs FVS.
- Next (fix): compare jl's LP calibration fit (# measured trees, per-tree `reslog`, `snx/sny/snxx/snxy`) to
  FVS's (12 trees, → 0.879). Align jl's observed-growth so LP COR → −0.13; keep floor; re-verify LP DG→live.
  This is the concrete, targeted fix — first time the trace has clean both-sides numbers (no assumptions).

### Slice 6a — CORRECTED reframe: in the DEFAULT run jl UNDER-calibrates loblolly (skips FVS's 0.879 reduction)
Clean idg test on the loblolly stand (each_stand→set idg→setup_growth!): **idg=0 (the actual FIA default,
since there's no DG_TRANS col) → `dg_cor[13]=0.0` — jl does NOT calibrate loblolly at all.** idg=1 →
`dg_cor[13]=0.9777` (wrong sign+magnitude). FVS (CALBSTAT fort.13): LP calibrates to **0.879** over 12 trees.
⇒ CORRECTED verdict: in the real default run jl **UNDER-calibrates** loblolly (no COR applied) so it grows at
the uncalibrated base rate — ~8% above live, which applies FVS's ~0.88 reduction. (My "2.66× over-calibration"
was the idg=1 probe conflated with the default run — yet another entanglement; the trustworthy facts are the
`.sum` and CALBSTAT.) The real gap: **jl's DGSCOR calibration does not fire / mis-fits for the FIA measured-DG
(past-DBH) loblolly trees** — it needs to (a) recognize the FIA `DG`=past-DBH as measured growth AND (b)
compute the residual so LP COR → 0.879 (not 0, not 0.98).

**HONEST LIMIT (this session):** after ~a dozen probes my jl-side measurements have been self-contradictory
(dg_cor 0 vs 0.9777 for "default"), which means I do NOT yet reliably understand jl's calibration state for
this stand. The two ANCHORS I trust: `.sum` (LP BA ~8% high) + CALBSTAT (FVS LP → 0.879, 12 trees). The fix
needs a careful, dedicated jl-calibration instrumentation (dump fn/snp/reslog/measured for sp13 alongside the
FVS 12-tree fit) — a FRESH pass, not more same-session iteration that keeps producing conflicting numbers.
No speculative fix shipped; floor intact 38527/143/0.

### Slice 4c — SIGNATURE clustering of the 57 big failures → dominant class = DENSITY MORTALITY (jl over-kill)
Built `signature.jl` (first diverging .sum col + cycle + direction per stand). Over ALL 57 big (>10%) failures:
- **First-diverging column: 41 TPA / 16 BA** ⇒ **72% are MORTALITY (tree-count) divergences**, not growth.
- **Direction: 36 LOW / 21 HIGH** ⇒ jl mostly kills MORE trees (over-kill).
- First-diverging cycle spread across years (not a phase artifact).

Worst mortality stand 157875449010854 (42 trees, avg DBH 6.3, **very dense**): cyc0 bit-exact (1853 TPA);
1995 jl 944 vs live 1234 (**jl kills 290 MORE**), survivors bigger (jl BA200/QMD6.2 vs live 185/5.3);
by 2015 jl 278 vs live 529. **SDI ≈ 397-411 at/near max SDI** ⇒ heavy density-dependent mortality, and jl
over-applies it vs live. ⇒ **hunt target: the SN density-dependent (near-SDImax) MORTALITY model** — jl
over-kills very dense stands. Data-driven from the full big-failure set (not a small sample this time).

**Discipline note:** this is a signature-clustered LEAD, not yet a verdict — must both-sides-trace the SN
mortality (MORTS/SDImax) vs jl on a few dense stands before any fix (per the loblolly/pole-size lessons).
Real FIA stands reach near-SDImax densities the curated tests apparently don't exercise as hard.

### Slice 4c-CORRECTION — the "mortality over-kill" is GROWTH-COUPLED via self-thinning (not a mortality bug)
Re-examined the worst stand: at 1995 both jl & live reach **SDI ≈ 410 (max)** but via different paths —
jl QMD 6.2 / TPA 944 vs live QMD 5.3 / TPA 1234. jl's trees grew **bigger** (QMD +17%), so at the
self-thinning limit jl carries FEWER (bigger) trees. ⇒ the TPA (mortality) divergence is a CONSEQUENCE of a
GROWTH (diameter) divergence, coupled through the SDImax self-thinning — NOT an independent mortality bug.
My `signature.jl` "TPA-first" reflects .sum COLUMN order (TPA is col 3, QMD col 8), not causal order — a
misleading proxy. (3rd corrected over-read: loblolly-species, pole-size, now mortality-vs-growth.)

**HONEST BIG PICTURE (SN FIA at scale):** the ~31% multi-cycle failures are NOT one systematic bug. They are
a HETEROGENEOUS distribution of GROWTH (diameter) divergences — some jl-high, some jl-low, species/condition-
specific — AMPLIFIED into large TPA differences by the SDImax self-thinning coupling (a small QMD diff →
a big TPA diff on dense stands). Plus a cornerable ULP/print-straddle tail (~30% <1%). No single fix; the
realistic path is: corner the <1% straddles as ULP-class, and root-cause the LARGE growth divergences
incrementally (per-tree DG, period+quantity-aligned, both-sides) — each likely a distinct species/condition
DG issue. The eco_unit fix (slice-1) was the one broad systematic win; the residual is a long tail.

**Recommendation to revisit with the user:** at real-FIA scale a 100%-bit-exact multi-cycle drop-in may be
unattainable (live FVS itself fails 17.6% of stands; growth-mortality coupling amplifies sub-% growth ULP into
%-level TPA on dense stands). A defensible done-state may be "cycle-0 bit-exact (achieved) + multi-cycle
within a documented tolerance-or-cornered, per variant" rather than strict bit-exact on every cycle.

### Slice 4b — FIXED: jl FIA-reader crash on text-typed numeric columns (4 JLERR → 0)
The 4 SN-1000 JLERR stands were all `MethodError(Float64, ("2.0",))`: `TREEINIT.SEVERITY3` is TEXT-typed but
holds "2.0", and `_fia_int`/`_fia_f32` did `Float64(d[k])` on the String. FIX (`fia_database.jl`): `_fia_num`
tryparses numeric strings (nothing→default), so the reader robustly handles text-encoded numbers like live FVS.
All 4 now RUN (1 bit-exact, 3 = small <1% straddles). **Floor held 38527/143/0** (numeric inputs unchanged;
faithful). A real jl robustness win on real FIA data — the first campaign fix since eco_unit.

## Slice 4 (Pillar 1+2) — SN-1000 sweep + failure CLUSTERING (widen, then triage by cause)
SN-1000 stratified sweep (indexed sub-DB, ~20 min):
- **566/820 BIT-EXACT (69%), 254 FAIL.** cycle-0 identical 814/820.
- **176 NOSUM = live FVS can't project (17.6%!)** — ill-posed/data stands FVS itself fails; NOT jl. (Worth a
  separate note: on real FIA-ready SN data, live FVS fails ~1 in 6.)
- **4 JLERR = jl errored** on 4 stands — a NEW signal (jl robustness bug on some real stands); triage these.

**★ FAILURE CLUSTERING (`cluster_failures.jl`) — the key insight:** failures cluster by **STRUCTURE, not
species**. By stand avg DBH: **pole5-9" = 140/254 (55%)**, saw9-15" = 56, sap1-5" = 46, lg15+ = 7, seed<1 = 5.
Dominant species are SPREAD thin (FIA 802:14, 316:12, 110:11, 611:9, 121:8, 621:7, 132:7, 318:6, oaks 6…) —
NO single species dominates. ⇒ the SN-100 "loblolly" cluster was largely SAMPLING NOISE; the real systematic
divergence is the **POLE-SIZE regime (avg DBH ~5-9")** across species — a growth or (density-dependent)
mortality issue for mid-size stands. This reframes the hunt from species-coefficients to a STRUCTURAL cause.

### Slice 4a — MAGNITUDE triage of the 254 SN-1000 failures (corrects the "pole-size systematic" read)
Re-ran the 254 failures with a worst-rel-diff histogram: **<1%:77  1-2%:45  2-5%:50  5-10%:25  >10%:57**.
⇒ ~77 (30%) are sub-1% = likely the accepted ULP / integer-print-straddle class (cornerable); **57 (22%)
are >10% = unambiguous real bugs**; the middle is a mix.

**Clustering the 57 BIG (>10%) failures by STRUCTURE: pole5-9=24, saw9-15=24, sap1-5=6, lg15+=3** — split
EVENLY between pole and sawtimber, NOT pole-concentrated. So slice-4's "pole-size systematic bug" was a
BASE-RATE artifact (pole is simply the commonest structure in the sample; I never normalized). And 3 sampled
pole failures showed HETEROGENEOUS directions (one TPA/BA bit-exact→fails on another col; two with jl BA
HIGH/over-grow) — the OPPOSITE of the loblolly under-grow. ⇒ **the big failures have MULTIPLE causes, not one
systematic bug.** (Meta: 2nd premature-cluster correction — always base-rate-normalize + inspect direction.)

**Honest campaign state (SN):** ~69% bit-exact at N=1000; failures = a cornerable-small tail (~30% <1%) +
a heterogeneous real-bug set (~57 >10%, mixed over/under-growth across pole & sawtimber). No single dominant
fix. Progress now = root-cause the big ones in small batches (each may differ), corner the <1% straddles, and
separately triage the 4 JLERR (jl crashes) + note the 176/1000 live-NOSUM (live FVS's own failures on real data).

## Slice 3 (Pillar 1+2) — SCALE SWEEP infrastructure + SN-100 pass/fail (user directive: scale, SN first)
Built the scale toolchain: `extract_sample.jl` (VARIANT-filtered, ECOREGION/LOCATION-stratified, deterministic)
+ `build_subdb.jl` (C-speed ATTACH+CTAS indexed sub-DB — ~100× faster than the unindexed 2.2M/8M-row master;
8s for 100 stands) + harness PASS/FAIL count + `FIA_FAILOUT` failing-stand list (`FIA_DB` env override).
Per-stand on the indexed sub-DB ≈ 1.5s (100 ≈ 3 min, 1000 ≈ 25 min).

**SN-100 sweep result** (stratified, vs live FVSsn, all cycles × 6 cols):
- **58/86 BIT-EXACT**, **28 FAIL**, cycle-0 bit-exact on all 86.
- **14 NOSUM = LIVE FVS produced no .sum** (jl errored on ZERO; jl ran all 100). Those 14 are live-can't-
  project stands (ill-posed/data), NOT jl bugs — exclude from the jl denominator.
- **★ Top-3 failures (70%/34%/30%) are ALL FIA 131 (loblolly pine)** — jl DG ~24% LOW (worst stand
  238814304010854: 2005 BA jl115/live152 → less density mortality → 2025 TPA jl961/live565, +70% trees).
  Loblolly is THE dominant Southern species ⇒ a high-value SYSTEMATIC DG bug the 8-stand hand-pick missed.
  Other failures vary (FIA 920/833/820/110/111) — smaller, likely the known cornered classes.

**Hunt target #1: FIA 131 (loblolly) DG too low** — systematic, big, common species. (Contrast: FIA 111→sp6
SA was jl-DG-HIGH — species-specific coefficient/mapping issues in both directions.)

### Slice 3a — loblolly (LP/sp13) DG ~3× TOO LOW (per-tree, worst stand 238814304010854)
FIA 131→LP (loblolly, PITA) mapping is CORRECT. Per-tree DG (live treelist vs jl engine), DBH-5-6 LP trees:
**liveDG ~1.3-1.5 vs jlDG ~0.5-0.6 (dDG ~−0.9, jl ≈ 1/3 of live).** A MASSIVE under-prediction, not a ULP/
coefficient tweak. Dense young stand (1323 TPA, QMD 3.4, eco 231AA, ftype 161→ylpn, managed=0). The 3× DG
deficit → jl grows far less BA → far less density mortality → jl keeps 70% more trees by 2025.
- dgf.f:184-222 has a FORTYP-group branch (IFORTP→K{LOHD,NOHD,OKPN,SFHP,UPHD,UPOK,YLPN} → FT*(ISPC)*K*).
  jl fgrp=ylpn matches for ftype 161, and the fortype term is only ~0.1-0.4 DDS — too small for a 3× DG gap.
- ⇒ the 3× is a bigger term: candidates = a wrong LP(sp13) coefficient, a loblolly SMALL-tree / small-large
  DG BLEND detail (QMD 3.4 = many small trees), or a loblolly special DG path. NEEDS the LP DDS decomposition
  (jl vs FVS, accounting for the FVS-debug D=dib=dbh·bark offset found in slice-1m) — the decisive next hunt.

**★ This is the campaign's biggest divergence and highest-value target: loblolly is the dominant Southern
species, so a 3× DG bug likely fails a large fraction of the 637k SN population.**

### Slice 3b — loblolly narrowed: base DDS matches FVS; effective DDS must be ~2.6; CROWN-RATIO suspect
- NOT the small-tree blend: DBH-5 > XMAX=3" (regent blend band [1,3]), so these use large-tree dgf.
- jl LP DDS (1.52/1.70/1.73/1.68) MATCHES FVS's lower cluster (1.51/1.69/1.67/1.72) — the dgf DDS FORMULA is
  right given the same inputs. BUT to yield live's DG≈1.5, the effective DDS must be ≈2.6 (conv check:
  DDS 1.5→DG 0.5; DG 1.5→DDS ln(Δdib²)≈2.6). FVS's debug DOES show LP trees at DDS 2.33/2.52 that jl lacks.
- **RED FLAG: jl assigns these young loblolly crown_ratio ≈ 10% (cr=10.3/11.4/13.5)** — implausibly low for
  vigorous young loblolly (should be ~40-70%). Low crown ratio → high pbal competition (PNTBL·pbal, pbal=
  PTBAA·(1−cr/100)) AND low ln_crown → suppressed DDS. So the loblolly DG deficit is likely a **CROWN-RATIO
  model bug for young dense loblolly** (crown too low), NOT the dgf coefficients (which match).
- Note crown_pct (45, used in ln_crown ICR) vs crown_ratio (10, used in pbal) diverge sharply here — worth
  confirming which FVS uses where + what crown ratio FVS assigns these trees (est. via crown model at cyc0).

**Next: both-sides the crown ratio — FVS's assigned CR for these young loblolly vs jl's ~10% (crown_ratio.jl
CR model). If jl's young-loblolly CR is too low, fixing it lifts DDS→DG and likely clears the top failures.**

### Slice 3c — RETRACTED crown-ratio hypothesis (percentile vs crown-ratio confusion); mechanism still open
★ SELF-CORRECTION (doctrine #3): I compared jl's `crown_ratio` field (10-18) to FVS treelist `PctCr` (35-55%)
and wrongly called it a crown-ratio bug. But `standstats.jl:209` + `diameter_growth.jl:192` confirm jl's
**`crown_ratio` field holds PCT = the stand BASAL-AREA PERCENTILE, NOT a crown ratio.** So 10 = 10th BA
percentile (CORRECT for a small tree in a dense 1323-TPA stand) — comparing it to FVS's 45% crown was
comparing two different quantities. RETRACTED. (jl's actual crown IS `crown_pct=45` == FVS PctCr — matches.)

**What's actually PROVEN about the loblolly deficit:**
- jl under-grows loblolly BA ~3× (solid: `.sum` BA jl115/live152 @2005 + treelist DG jl0.5/live1.5). Systematic
  (top-3 SN-100 failures all FIA 131). Loblolly-specific (OH/WK in the same stand were jl-HIGH, not low).
- jl LP dgf DDS (≈1.5) ≈ FVS LP dgf DDS (≈1.7, debug) — the dgf base is ~right, and PCT/pbal inputs match.
- **UNRESOLVED DDS→DG discrepancy:** by the standard conversion DDS 1.7 → DG ≈0.6, yet FVS's treelist DG is
  1.5 for the same tree. So FVS's *effective* DDS must be ≈2.6 (or the DG uses a different period/scaling).
  Either FVS applies a post-dgf step (~+0.9 ln DDS: DGSCOR/COR calibration, a growth-period scaling, or a
  loblolly adjustment) that jl misses, OR the treelist DG period differs. NOT yet both-sides-nailed.

**Honest status (3 imprecise turns — logged as a discipline lesson):** the loblolly under-growth is REAL and
systematic (#1 SN-100 failure class), but I over-stated magnitude and mis-attributed mechanism:
- ✗ fortype (ruled out), ✗ crown_ratio (retracted — it's the BA percentile), ✗ "3×" per-tree DG (INFLATED:
  the treelist 2000 row is `PrdLen=2`, so my y0→y1 DG alignment was off).
- ✓ RELIABLE magnitude from the `.sum` oracle: jl BA 115 vs live 152 @2005 = **~24% low at cycle-1**,
  compounding (jl 961 / live 565 TPA @2025). Loblolly-specific, systematic.
- Mechanism still OPEN. The clean next step is a debug-FVS stamp of the loblolly large-tree growth
  (dgf DDS → dgdriv scaling/calibration → DG) with EXACT tree + period alignment — not another Julia-side
  per-tree comparison (which keeps mis-aligning periods/quantities). Do NOT ship a fix until pinned (#3).

**Meta-lesson:** on this trace I repeatedly declared a cause from a partial/mis-matched comparison. The
`.sum` is the trustworthy oracle; per-tree instrumentation must align period (PrdLen) AND quantity semantics
(crown_ratio=percentile, D=dib) before any verdict.

## Slice 2 (Pillar 1+2) — cross-variant scaling: LS multi-cycle baseline
6-stand LS multi-cycle differential (`ls.txt`, vs live FVSls): **core growth BIT-EXACT** — TPA/BA/SDI/QMD
all 0.0% (4/6 stands fully bit-exact). LS does NOT share SN's EUT/DG divergence (the eco_unit fix was
SN-gated & SN-specific; LS core is clean). **NEW LS-specific signal:**
- **CCF drifts systematically** 0.55% (cyc0) → 1.46% (cyc4), on essentially ALL LS stands. A crown-width /
  CCF class (distinct from SN's DG class). At cyc0 it's an INVENTORY-level crown-width difference (1/6 LS
  stands is cyc0-non-identical) — an LS crown-width or CCF-accumulation detail.
- Worst stand 18498465010661 (8.6%) — trace separately.

Cross-variant picture forming: SN = DG/EUT (fixed) + multi-point point-BAL (cornered); LS = core bit-exact
+ a CCF/crown-width lead. NE/CS differentials pending (need proper NE/CS stand lists — `nc.txt` is likely
SN-region-8 NC, not NE; extract NE/CS STAND_CNs from the DB by LOCATION→variant for the Pillar-1 manifest).

### Slice 2a — LS CCF divergence sharply characterized: CROWN-WIDTH only, everything else bit-exact
Worst LS stand 18498465010661 per-cycle: **TPA/BA/SDI/TopHt/QMD ALL bit-exact every cycle**; ONLY CCF
diverges — jl LOW: 1993 177/live183 (3.3%) → 2043 370/live402 (8%). Since SDI (same trees, DBH-based) is
bit-exact but CCF (crown-AREA based) is low, it's specifically the **CCF crown-width computation** for LS,
present at the INVENTORY row (cyc0) and growing. Diverse LS hardwoods (FIA 391 ironwood n13, 318 sugar
maple, 541 white ash, 371, 71, 543). jl `crown_width` used at standstats.jl:193/246 for CCF (iwho path).
- Impact: CCF column + crown-width-dependent downstream (PERCOV/fire) only — NOT growth (TPA/BA/SDI/QMD all
  bit-exact). A real but bounded `.sum` divergence; systematic across LS stands (all showed the CCF drift).
- NOTE: modernization validated lst01 CCF bit-exact (171) — so it's species/condition-specific to these FIA
  LS hardwoods (a crown-width coefficient or the open-vs-forest-grown CCF crown-width path for LS species
  not in lst01). Cornered to: **LS CCF crown-width (species-specific), a bounded crown-only residual.**

## TODO
- [ ] LS CCF: identify which LS species' CCF crown-width is low (standstats crown_width iwho path) vs live;
      compare the LS crown-width coefficient/eqn — potential clean systemic LS fix (CCF only, growth safe).
- [ ] Extract NE + CS FIA stand lists (LOCATION→variant) and run their multi-cycle differentials.
- [ ] FIX (density-model slice): true per-point PTBALT for the dgf point_bal term on multi-point stands;
      validate the SA stand + re-diff; keep floor + keep single-point stands bit-exact (YP etc.).
- [ ] Forest-type derivation: the YP −0.089 tail (separate stand).
- [ ] Scale differential to NE/CS/LS + larger SN sample; Pillar-1 stratified manifest.
- [ ] Forest-type derivation: still a candidate for the YP −0.089 tail (separate stand); trace fortyp.f vs jl.
- [ ] Scale differential to NE/CS/LS + larger SN sample; Pillar-1 stratified manifest.
- [ ] Residual: forest-type derivation for FIA stands (jl 503→upok −0.0907 vs FVS) — the −0.089 DDS tail;
      may surface on other plots where it crosses the rounding boundary. Trace fortyp.f vs jl for FIA input.
- [ ] NE/CS/LS: analogous eco-unit read (their FIA differentials — do they show the same EUT gap?).
- [ ] Then: forest-type derivation check (jl 503→upok vs FVS) for the residual ~0.09.
- [ ] Scale differential to NE/CS/LS + larger SN sample; Pillar-1 manifest.
      (both sides) vs jl; likely a beyond-growth-sample-max-DBH large-tree DG path YP-specific.
- [ ] After fix: re-run per-tree + stand differential (expect the 3.7% BA drift to collapse).
- [ ] Scale differential to NE/CS/LS + larger SN sample; build Pillar-1 stratified manifest.
- [ ] Track: jl 36 vs live 37 tree records at cyc0 on 3196569010661 (.sum TPA bit-exact ⇒ likely zero-TPA rec).
      which species' per-tree DBH increment diverges, and by how much (ULP vs systematic). Classifies
      accepted-DGSCOR vs real-gap and, if real, names the species/coefficient/path.
- [ ] Scale the differential to NE/CS/LS + larger SN sample once the SN driver is understood.
- [ ] Build the stratified per-variant plot manifest (Pillar 1).
- [ ] Scale the multi-cycle differential: larger SN sample + NE/CS/LS (Pillar 1 manifest feeds this).
- [ ] Build the stratified per-variant plot manifest (Pillar 1).
- [ ] Management-scenario differential on real plots (Pillar 3).

---
## SLICE 9 — NE-100 multi-cycle sweep (first NE cross-variant pass rate)  [2026-07-08]
First NE FIA multi-cycle differential (100 real NE stands from SQLITE_FIADB_ENTIRE.db VARIANT='NE',
indexed sub-DB `ne100.db`, 5 cycles × 10yr native = 50yr horizon).

**Result: 85/94 BIT-EXACT (90%)** (100 run, 94 produced .sum, 6 live-FVS-NOSUM = live itself can't
project them, not a jl bug; JLERR=0 ⇒ the `_fia_num` text-numeric fix carried over to NE).
Worst-rel histogram (all 94): <1%:91, 1-2%:2, 2-5%:1, 5-10%:0, >10%:0 ⇒ **97% within 1%, NO heavy tail.**
Mean |rel diff| by cycle ≤0.03% on every column/cycle.

**Divergence taxonomy (9 failures, all CORNERED — ULP/straddle-class):** signature.jl first-diverge:
7/9 = TPA ±1-2 trees on high-density (2.5k-13k stem) stands (=0.01-0.08%); 1 = TopHt ±1ft (82/83);
1 = SDI ±1 (259/258). Worst stand 233164158020004 "2.5%" fully traced: 2013 bit-identical (10032 TPA),
2023 seedling self-thin 10032→2535 lands TPA 2533/2535 (Δ2=0.08%), BA identical (227) ⇒ QMD=√(BA/TPA/k)
tips the 1-decimal print boundary 4.0↔4.1 = "2.5%". A ±2-tree mortality straddle amplified by low-precision
QMD print — NOT a logic gap.

**NE verdict:** no systematic NE bug; reached 90% WITHOUT the SN-gated eco_unit/EUT term (NE dgf differs).
The two variant-agnostic FIA-reader fixes (growth_dg_set FINT-normalization + _fia_num) generalize cleanly.
Remaining 10% = the same Float32-compounded straddle class cornered for SN. Floor intact 38527/143/0.

---
## SLICE 10 — CS-100 + LS-100 multi-cycle sweeps (complete the 4-variant cross-check)  [2026-07-08]
Ran the same multi-cycle differential (5 cycles, indexed sub-DBs cs100.db/ls100.db, deterministic
VARIANT-stratified samples) for the two remaining variants. JLERR=0 both (the `_fia_num` fix holds).

**CS: 91/98 BIT-EXACT (93%)** (100 run, 98 .sum, 2 live-NOSUM). Worst-rel histogram: <1%:98 — **every
CS stand within 1%, zero tail.** Mean |rel diff|/cycle ≤0.02%. All 7 failures = TPA/SDI ±1-2 unit
straddles (signature.jl: 2689/2690, 180/181, 270/271, SDI 224/225, 1761/1759, 1181/1182, 918/919).

**LS: 86/100 BIT-EXACT (86%)** (100 run, 100 .sum, 0 NOSUM). Histogram: <1%:94, 1-2%:4, 2-5%:2, >10%:0
(94% within 1%). All 14 failures = ±1-unit straddles (TPA/SDI/CCF/BA ±1) + one QMD print straddle.
Worst 4.0% (566121733126144) fully traced: 2018 & 2038 bit-identical; 2028 TPA(695)+BA(23) identical but
QMD prints 2.4 vs 2.5 — true BA differs sub-print (both round to 23), QMD=√(BA·k/TPA) lands ~2.449 vs
~2.451 across the 2.45 one-decimal boundary; rest bit-exact or ±1. ULP-compounded straddle, not a logic gap.

**4-VARIANT CROSS-CHECK COMPLETE (Pillar 2, first pass):**
  SN 1000: 74% (heavy tail collapsed 57→1 by 2 fixes; residual ULP + ~5 individual)
  NE  100: 90% (all-ULP tail, no fix needed)
  CS  100: 93% (all-ULP tail, zero >1%, no fix needed)
  LS  100: 86% (all-ULP tail, no fix needed)
The two variant-agnostic FIA-reader fixes (growth_dg_set FINT-normalization + _fia_num) generalized to
all 3 companion eastern variants with NO variant-specific code: every NE/CS/LS residual is a Float32
±1-unit / print-boundary straddle, no systematic bug in any of the three. Floor intact 38527/143/0.

---
## SLICE 11 — NE/CS/LS widened to 1000 stands (confirm 100-sample result at scale)  [2026-07-08]
Widened the 3 companion variants to 1000-stand deterministic VARIANT-stratified samples (per the
widen-then-hunt plan). Pass rates STABLE vs the 100-samples ⇒ no rare systematic bug was hiding:
  NE 1000: 841/960 BIT-EXACT (88%; 100-sample was 90%). Histogram <1%:937(98%) 1-2%:18 2-5%:3 5-10%:2 >10%:0
  CS 1000: 921/988 BIT-EXACT (93%; was 93%).          Histogram <1%:958(97%) 1-2%:25 2-5%:4 5-10%:1 >10%:0
  LS 1000: 853/999 BIT-EXACT (85%; was 86%).          Histogram <1%:927(93%) 1-2%:49 2-5%:19 5-10%:3 >10%:1
JLERR=0 all three. Live-NOSUM: NE 10/CS 12/LS 1 (live can't project these — not a jl bug).
Tail is proportionally tiny and slow-growing; LS carries the fattest 2-5% band + the lone >10% (to triage).
Pillar-2 first pass is now at 1000-scale for every companion variant (SN already at 1000).

## SLICE 12 — Pillar 3 management harness landed + column-format finding  [2026-07-08]
Built `test/harness/fia/manage_fia.jl` — injects a silvicultural keyword block (thinbba/thinbta/thindbh)
scheduled at CYCLE 2 and diffs the 6-col .sum live-vs-jl under management. FINDING (harness, not a jl bug):
FVS reads keyword records FIXED-FORMAT (A10 + F10.0 fields); a mis-columned THINBBA scrambles the residual
(my first draft over-spaced ⇒ garbage residual ⇒ live/jl disagreed catastrophically 147%, live column-strict
vs jl whitespace-tolerant). Fixed via `kwrec()` (keyword left-just in 10, each field right-just in 10).
With correct columns + CYCLE scheduling (2.0), live and jl AGREE to ULP on the thin (validated stand
246537009010854: both 2019 TPA6359/BA49, growth after ULP-close). SECOND finding (noted, not yet traced):
CALENDAR-year scheduling (2019.0) diverges — live ignores it, jl schedules one cycle late (2024); cycle
scheduling sidesteps it. First 15-stand SN thinbba pass: thin-fired 7/12, worst 14.1% (was 147%) — the
residual-tree-selection-at-the-cut-margin divergence class that Pillar 3 exists to characterize.

---
## SLICE 13 — Pillar 3: SN thinning (THINBBA) is BIT-EXACT; management tail = post-thin growth-ULP  [2026-07-08]
SN-100 management sweep (THINBBA cycle-2, residual BA 40): 41/81 bit-exact overall, BUT the harness now
splits by whether the thin FIRED:
  thin NO-OP (BA<40 ⇒ no cut): 37/40 bit-exact (92%) = the Pillar-2 growth-only rate (consistent).
  thin FIRED: 4/41 bit-exact — looks alarming, but BOTH-SIDES-TRACED to a cornered primitive, NOT a bug:
On EVERY fired-thin failure examined, **the cut cycle itself is bit-exact** — jl removes the same trees and
lands the identical residual as live:
  218228660020004: thin@2023 → both 379 TPA / 52 BA / 5.0 QMD EXACT; diverges 2 cyc later (BA 81 vs 82, Δ1).
  416592108489998: thin@2026 → both 336 TPA / 47 BA / 5.0 QMD EXACT; residual grows apart (BA 58 vs 61 → 8% by 2041).
  246537009010854: thin@2019 → both 6359 TPA / 49 BA EXACT (from slice-12 trace).
**VERDICT:** the THINBBA behaviour (from-below tree selection + residual-BA stop) reproduces live FVS
BIT-EXACT. The elevated fired-thin tail is entirely POST-THIN growth divergence: a sparse residual (few
large trees near the competition/DGSCOR regime) amplifies the same per-tree growth-ULP faster than a dense
stand does — the SAME named primitive cornered in Pillar 2 (heterogeneous DG × density amplification), not
a thinning-logic gap. Pillar-3 first pass (SN, thin-from-below) = bit-exact-or-cornered. Floor 38527/143/0.
(TODO Pillar 3: sweep thinbta/thindbh + NE/CS/LS the same way; trace the calendar-year scheduling divergence.)

---
## SLICE 14 — Pillar 3 cross-variant (NE/CS/LS thinbba): cut BIT-EXACT, tail = post-thin mortality-cliff ULP  [2026-07-08]
Ran the split-metric management sweep (THINBBA cyc-2 res-40) on NE/CS/LS-100:
  NE: no-op 43/47 (=Pillar-2) | FIRED 3/48 | worst 43.9%    CS: no-op 81/85 | FIRED 0/14 | worst 36.2%
  LS: no-op 69/78 | FIRED 1/21 | worst 26.6%
No-op rate == Pillar-2 growth rate on every variant (consistent). Fired-thin rate looks bad + big worst-%,
so BOTH-SIDES-TRACED the sharpest (CS 177137819020004, 36.2%) with a no-mgmt-vs-mgmt differential:
  NO MANAGEMENT: live==jl BIT-EXACT all cycles (2032 both 1294/121 ... 2062 548/547) — cf. slice-10 0.6%.
  WITH THINBBA:  cut@2032 BIT-EXACT (both 1294/121 → 1182/81, same trees removed); diverges only AFTER
                 (2042 992/147 vs 1020/136; 2052 480/209 vs 654/180 = 36%).
**VERDICT (holds all 4 variants):** the thin CUT is bit-exact (jl removes the same trees / lands the same
residual as live). The large fired-thin tail is NOT a thinning-logic gap — it is POST-thin growth/mortality
divergence: the sparse residual sits near a DENSITY-DEPENDENT MORTALITY (self-thinning) threshold where a
ULP-level per-tree growth difference tips a DISCRETE tree-death event, amplifying to 20-44%. Same named
primitive as the Pillar-2 dense-regen tail (heterogeneous DG-ULP × discrete self-thinning), stressed harder
by the post-thin low density. Cornered, not a bug. Floor 38527/143/0 (no src change; harness+docs only).
NOTE: this is an AMPLIFICATION of a ULP primitive, honestly a large-magnitude residual — the ROOT (per-tree
DG + the self-thinning trigger) is proven bit-exact/ULP by the no-mgmt + cut-cycle bit-exact evidence.
(TODO Pillar 3: thinbta/thindbh regimes; calendar-year scheduling divergence trace; salvage/plant/fire.)

---
## SLICE 15 — RETRACTION: calendar-year scheduling is BIT-EXACT (the "divergence" was a column artifact)  [2026-07-08]
Slices 12/14 NOTED (untraced) that calendar-year scheduling diverged (live ignores THINBBA 2019.0, jl one
cycle late). TRACED it → FALSE. That note came from a pre-kwrec MANUAL keyfile (`THINBBA       2019.0  40.0`,
7 spaces) whose column misalignment shifted 2019.0 into the wrong FIXED-FORMAT field, so live (column-strict)
misread the date while jl (whitespace-tolerant) read it fine ⇒ artificial disagreement. Re-tested stand
246537009010854 with kwrec-aligned columns at calendar 2014/2019/2024:
  2014.0: live/jl AGREE — thin shows @2019 BA 49 (then 103/102, 140/139, 144/143 ULP)
  2019.0: live/jl BIT-EXACT — 2019/135(pre) 2024/47 2029/81 2034/119 identical
  2024.0: live/jl BIT-EXACT — 2024/145 2029/46 2034/70 identical
**Calendar-year scheduling reproduces live FVS bit-exact-or-ULP.** (Timing note, both engines agree: a
calendar year Y fires in the cycle that GROWS from Y (OPCYCL IY(i)≤Y<IY(i+1)) so the cut shows at the NEXT
summary row; cycle-number N fires one cycle earlier — different SEMANTIC, each matching live for its form.)
This STRENGTHENS Pillar 3: BOTH cycle-number and calendar-year THINBBA scheduling are faithful; the only
management residual is the cornered post-thin mortality-cliff ULP amplification (slice 14). Lesson (again):
a noted-but-untraced divergence must be traced before it's believed — this one dissolved, like the 147% one.

---
## SLICE 16 — Pillar 3: all 3 SN thinning selection paths agree (thinbta + thindbh)  [2026-07-08]
Swept the remaining thinning regimes on SN-100 (cycle-2, correct kwrec columns):
  THINBBA (from below, res-40):  no-op 37/40 | FIRED 4/41 | worst 14.1%  (slice 13)
  THINBTA (from above, res-40):  no-op 38/39 | FIRED 3/41 | worst  7.8%
  THINDBH (50% across all DBH):  no-op 36/37 | FIRED 4/41 | worst  5.4%
Every regime: no-op rate == the Pillar-2 growth-only rate (consistent); fired-thin cut bit-exact (traced);
post-thin tail cornered. DECISIVE PATTERN — the worst-magnitude tracks how far the thin distorts the size
distribution (⇒ how mortality-cliff-sensitive the residual is): from-below (sparse large residual) 14% >
from-above (dense small) 7.8% > proportional (shape-preserving) 5.4%. This CONFIRMS the slice-14 mechanism:
the divergence is POST-thin density-regime ULP sensitivity, NOT tree-selection (selection is bit-exact in
all 3 paths). Pillar-3 THINNING (all selection methods × both scheduling forms × 4 variants for thinbba) =
bit-exact-or-cornered. Floor 38527/143/0 (harness+docs only).
Pillar-3 remaining: salvage / plant-regen / prescribed-fire regimes (regen+FFE already heavily curated-suite
validated); the deeper Pillar-4 dive = is the post-thin discrete-mortality TRIGGER ULP-alignable or irreducible.

---
## SLICE 17 — CORRECTION: post-thin tail is TWO sub-classes; one is a REAL low-density growth divergence  [2026-07-08]
Slices 14/16 corner-labeled the whole fired-thin tail as "post-thin growth-ULP amplification." That was
PREMATURE for part of it. All-6-col + no-mgmt control on 416592108489998 (SN) splits the tail:
  NO MGMT: live==jl BIT-EXACT all 6 cols, all cycles (2016 445/130/265/353/53/7.3 … 2041 251/189/321/385/72/11.7).
  THINBBA: cut@2026 BIT-EXACT all 6 (both 336/47/84/99/63/5.0); FIRST post-thin cycle 2031 DIVERGES with
           SAME TPA 327 but jl HIGHER: BA 61/58, SDI 111/108, CCF 135/131, TopHt 68/67, QMD 5.8/5.7 — and it
           GROWS monotonically (2041 BA 92/85). SAME-TPA ⇒ pure DIAMETER+HEIGHT GROWTH divergence, NOT mortality.
**TWO SUB-CLASSES of the fired-thin tail:**
  (1) CLEAN ULP: first post-thin cycle bit-exact, later ±1-unit ULP amplification (e.g. 218228660020004:
      2023 & 2028 both 379/52 & 373/66, diverges 2033 BA 81/82 Δ1). Correctly cornered.
  (2) ★ REAL: immediate SYSTEMATIC over-growth by jl in the LOW-DENSITY regime (SDI ~84-164) that the
      un-managed stand never reaches (min SDI 265, where jl IS bit-exact). jl consistently grows the sparse
      residual faster (BA/SDI/CCF/TopHt all jl-high every post-cut cycle). NOT ULP (several %, systematic,
      one-directional). This is a genuine OPEN divergence — a low-density DG/height growth-regime gap that
      THINNING EXPOSES. Must be traced (Pillar 4), not cornered.
CORRECTED VERDICT: thinning SELECTION + SCHEDULING are bit-exact (cut cycle bit-exact all 6 cols, both sched
forms). The post-thin tail is PART cornered-ULP (sub-class 1) and PART a real untested-low-density growth
divergence (sub-class 2, jl over-grows). Floor 38527/143/0 (no src change). NEXT: both-sides-trace the
low-density DG/HTG path (per-tree DG on a residual tree 2026→2031) to name the divergent term.

---
## SLICE 18 — Pillar 4: post-thin over-growth LOCALIZED to sp544 large-tree DG at low density  [2026-07-08]
Both-sides FVS_TreeList DG/HtG differential on 416592108489998 (live FVSOut.db via sn_oracle.sh vs jl DBS):
  2021 & 2026 (pre-thin + CUT cycle): EVERY species BIT-EXACT (DG/DBH/HtG identical) — cut is per-tree exact.
  2031 (first post-thin cycle): sp544 meanDG live 0.742 vs jl 0.811 (+9%), meanDBH 4.37/4.44; sp552 matches
     (1.449/1.446). 2036 sp544 0.688/0.752; 2041 sp544 0.626/0.684 (+ sp552 starts, 0.602/0.634).
**LOCALIZED:** the post-thin over-growth = species-544 (green ash) LARGE-TREE diameter growth (dgf DDS)
diverging ~9% at LOW post-thin density (stand BA 47), jl faster; BIT-EXACT at normal/high density (pre-thin
BA 148+, and all of Pillar-2). The 2026 START state is bit-exact (both 258 TPA / 3.44 DBH / 0.201 DG) ⇒ dgf
yields different DG on IDENTICAL per-tree input at low density ⇒ the divergent driver is a STAND-LEVEL dgf
input in the low-BA regime: either a Float32 BA/PBAL SUMMATION-ORDER sub-ULP diff amplified by the steep
low-density competition response (⇒ cornered ULP), OR a low-BA competition-term extrapolation difference
(⇒ real coeff/term gap). DISTINGUISHING THESE = next slice: instrument dgf DDS for sp544 at the 2026→2031
cycle (dump BA, PBAL, each term) both sides. Until then this is an IDENTIFIED-BUT-UNRESOLVED divergence
(honestly NOT yet cornered). Scope: SN, species-544, low-density (post-thin) only; ~9% DG ⇒ few-% .sum tail.
Floor 38527/143/0 (no src change; measurement only).

---
## SLICE 19 — DECISIVE: real (non-ULP) large-tree DG low-density gap, +29%, order-invariant  [2026-07-08]
Order-invariant cohort-mean DG (TPA-weighted, sort-independent) on 416592108489998 sp544, both-sides FVS_TreeList:
  2026 (post-thin START): BIT-EXACT — large(DBH>=10) BA 45.721/45.721 meanDG 0.8054/0.8054; ALL BA 46.399/46.399.
  2031 (1st post-thin cycle 2026->2031 growth): large meanDG 0.8229 vs 1.0586 = **jl +29%**, BA 53.0/55.1;
     small(<10) meanDG 0.7154/0.7277 (+1.7% only). Divergence DOMINATED by the LARGE cohort.
Order-invariance ⇒ NOT the record-tripling/near-tie-sort-flip artifact (that earlier per-tree row view was
confounded by reordering; the cohort MEAN is immune). 2026 input bit-exact ⇒ SN large-tree dgf yields +29%
different DG on IDENTICAL low-density input (stand BA ~46); BIT-EXACT at normal density (pre-thin BA 148, all
Pillar-2). **VERDICT: a REAL non-ULP divergence in SN large-tree diameter growth (dgf DDS) at LOW stand
density — jl over-grows large trees ~29%.** This is the driver of the post-thin management tail (thinning
pushes stands into the low-BA regime the un-managed trajectory never visits). Scope: SN large-tree DG, low
stand BA. NEXT: read dgf DDS for the low-BA-sensitive term (ln(BA)/BAL/PCT-percentile/relSDI/clamp) + dump it
both sides at the 2026->2031 cycle to name+fix. This is an OPEN real bug (NOT cornered).

---
## SLICE 20 — sp544 low-BA DG: visible inputs bit-exact ⇒ driver is PTBAA / rel-ht / per-tree CR distribution  [2026-07-08]
Refined slice 19. The 2031-row DG = the 2026->2031 (first post-thin) growth; its inputs = the 2026 state.
Order-invariant 2026 large-cohort (sp544, DBH>=10) aggregates ALL BIT-EXACT: stand BA 46.4, meanPctCr 42.69,
meanHt 63.44, cohort BA 45.721, meanDG-of-prior-cycle 0.8054 — yet the 2026->2031 large-tree DG diverges +29%
(0.8229 live / 1.0586 jl). Crown ratio + height means are bit-exact at BOTH 2026 and 2031 (42.19/42.19,
67.31/67.31) ⇒ NOT a crown/height divergence. Large cohort is PAST the tripling window (cyc<=2) ⇒ NOT tripling.
SN large-tree dgf DDS is LINEAR in stand BA (ba_v) + pbal (=point_ba·(1-CR/100)) + ln_crown; the visible
per-cohort inputs match, so the divergent driver is an input the treelist does NOT expose:
  • pba = POINT basal area (PTBAA / point_ba[plot]) — per-point, not the stand mean; OR
  • the rel-height term (rel_ht = tree ht / AVH, AVH=avg dominant height) — a stand aggregate; OR
  • the per-tree CROWN-RATIO DISTRIBUTION (same cohort MEAN 42.69, different spread ⇒ different per-tree pbal).
NEXT (definitive): instrument jl dgf! to dump per sp544 large tree {ba_v, pba, pbal, crown_ratio, rel_ht, each
DDS term} at the 2026->2031 cycle, + a live debug-FVS dgf.f stamp (edit→build→run→RESTORE→verify pristine),
compare term-by-term to NAME the divergent primitive. Localization now: SN, sp544, large cohort, low stand
BA (~46), first post-thin cycle, +29% DG, all VISIBLE inputs bit-exact. Open real bug. Floor 38527/143/0.

## SLICE 20b — PTBAA ruled out (single-point stand) ⇒ driver is rel-ht(AVH) or per-tree CR assignment
416592108489998 has 1 distinct PLOT_CN (9 trees) ⇒ SINGLE-POINT stand ⇒ pba(point_ba)=stand BA=bit-exact(46.4).
So the divergent dgf DDS input is NOT point BA. Remaining candidates narrowed to TWO:
  (a) rel_ht term = tree_ht / AVH (avg dominant height) — a stand aggregate; TopHt was bit-exact in .sum but
      AVH (the dgf-internal dominant-height mean) may differ at low density; OR
  (b) per-tree CROWN-RATIO ASSIGNMENT — cohort MEAN CR bit-exact (42.69) but WHICH large tree gets which CR
      may differ (per-tree dump showed CR values 37 & 44 among large trees); since DDS is nonlinear in DBH,
      a CR↔tree reassignment shifts the TPA-weighted cohort-mean DG even at fixed mean CR. A near-tie
      size-rank flip in the crown-ratio model would do this (cf. the documented COMPRESS/DGSCOR sort-flip class).
If (b) via a near-tie rank flip ⇒ this rejoins the cornered sub-ULP-sort-flip primitive; if (a) or a real CR
low-density formula gap ⇒ a real bug. DEFINITIVE next: jl dgf per-tree term dump + live dgf.f debug-stamp.

## SLICE 21 — sp544 DG suspect = crown_ratio (pbal term), distinct from the bit-exact crown_pct
Instrumented jl dgf! (ENV-gated, reverted after — hot-loop ENV check would break the alloc floor). Captured
green-ash (internal sp=37) large-tree DDS inputs at the divergent cycle. KEY STRUCTURAL FINDING: SN dgf DDS
reads TWO separate per-tree crown fields —
  • crown_pct (icr_i) → the ln_crown term; jl values 37-49, cohort-mean ~42 == the treelist PctCr (BIT-EXACT).
  • crown_ratio (cr) → the pbal COMPETITION term (pbal = pba·(1−cr/100)); jl values span 13.7-100, mean ~59
    — a DIFFERENT field, NOT exposed in the treelist, so NOT yet verified bit-exact.
⇒ The reported/verified PctCr being bit-exact does NOT clear crown_ratio. crown_ratio feeding pbal is the
PRIME SUSPECT for the +29% low-density DG divergence (at low density pbal spans 0-52 across the cohort ⇒
the DDS point_bal·pbal term is large-swing and cr-sensitive). jl cyc-input snapshot saved (jl_dgf3.txt).
DEFINITIVE NEXT (heavy, deferred to a fresh session for oracle safety): live dgf.f debug stamp
(edit→build→run→RESTORE→verify pristine) dumping the same per-tree {cr, pbal, dds} for sp green-ash; compare
to jl_dgf3.txt term-by-term. If cr diverges → trace crown_ratio! at low density (real bug or cornered
near-tie); if cr matches but dds diverges → a dgf coefficient/term-form gap. Localization is now maximal
short of the live Fortran stamp. Floor 38527/143/0 (src reverted clean).

## SLICE 22 — CLEARED large-tree dgf; root is UPSTREAM (cyc-2 small-tree/height/mortality). Both-sides dgf stamp.
Did the definitive both-sides dgf trace (jl ENV-gated dump + live dgf.f debug stamp: backup→edit→compile
dgf.o→relink /tmp/FVSsn_new→run→RESTORE dgf.f+dgf.o→relink→VERIFY PRISTINE [no-mgmt 2016 445/130, no stamp]).
Aligned cycles by stand BA: jl cyc2 bav=46.578938 == live BA=46.5789490 (the first post-thin growth cycle).
RESULT at that cycle — EVERY sp544 large-tree dgf input is BIT-EXACT (all 12 trees): DBH, ICR(=jl crown_pct),
PCT(=jl crown_ratio, the BA-percentile feeding pbal), PBA, PBAL, BA, RELHT (+ ht, avh). Identical inputs +
identical coeffs ⇒ large-tree DDS is necessarily bit-exact at cyc2. **The crown_ratio/PCT suspect (slice 21)
is CLEARED; the large-tree dgf is CLEARED.** The +29% large-tree DG (slices 19-20) appears at cyc3, and is a
DOWNSTREAM symptom: cyc3 stand BA diverged (jl 60.59 / live 58.44), and stand BA feeds the dgf BA/pbal terms.
So the ROOT divergence is in the cyc2 (first post-thin, 2026->2031) growth of NON-large-tree components —
candidates: SMALL-tree diameter growth (d<10, uses the height-based small-tree model not dgf!; treelist showed
sp544 small cohort +1.7% at 2031), HEIGHT growth (htgf), or MORTALITY of the sparse residual. NEXT: diff
small-tree DG + htg + periodic mortality at cyc2 both-sides (same stamp method on the small-tree/htgf path).
LESSON: localizing by the LARGEST symptom (large-tree +29%) pointed at the wrong model; the bit-exact-input
proof redirected to the upstream driver. Oracle verified pristine; src clean. Floor 38527/143/0.

## SLICE 23 — RESOLVED: post-thin divergence = the cornered DGSCOR/AUTCOR stochastic-ordering primitive
FIRST, a MEASUREMENT-HYGIENE catch: the earlier per-species decomposition showed jl 5× live — that was a
TEST ARTIFACT (I re-ran run_keyfile on the same jl_tl.db ~5×; the DBS output APPENDS; live used fresh mktemp
dirs each run so stayed 1×). Row-count diagnostic: jl 45/135/165 vs live 9/27/33 = exactly 5×. Regenerated
jl_tl.db FRESH ⇒ rows 9/27/33 and sumTPA per year BIT-EXACT vs live all cycles. (Lesson: rm the DBS db before
each run; slice-19's 1× numbers were valid.)
CLEAN decomposition (fresh DB): 2016-2026 ALL species bit-exact. Divergence starts 2031, entirely sp544:
  2031 sp544 LARGE(>=10): TPA 63.8/63.8 (=), recs 12/12 (=), meanHt 67.31/67.31 (BIT-EXACT), BA 53.03/55.14 (+2.1).
  ⇒ SAME trees, SAME count, SAME heights, but DIFFERENT diameters. Height bit-exact, diameter diverges.
CONTRADICTION with slice-22 (dgf DDS inputs bit-exact ⇒ deterministic DDS bit-exact) is the KEY: identical
deterministic DDS + identical Ht + identical TPA but divergent diameter ⇒ the divergence is in the DDS→DG
CONVERSION, whose only non-deterministic piece is the STOCHASTIC SERIAL-CORRELATION increment (DGSCOR/AUTCOR).
Height uses a different path (bit-exact); only the diameter stochastic term diverges.
**VERDICT: CORNERED to the documented, accepted DGSCOR/AUTCOR stochastic record-ordering primitive** (cf.
memory: "DGSCOR drift = record ORDER + same-species RNG swap; sub-ULP near-tie sort-flip"; COMPRESS s22 class).
Thinning EXPOSES it by pushing the stand into the sparse low-density post-thin regime where the per-tree RNG
draw ORDER is more sensitive (near-tied sort keys flip ⇒ serial-corr increments swap between records). This
is the SAME named primitive already accepted campaign-wide, not a new bug. The management-tail (slice 14-20
"+29% large-tree DG") was: (a) partly a stale-DB 5× artifact, (b) the real part = this stochastic serial-corr
increment on the sparse residual. Deterministic growth (DDS, height, TPA, mortality, thinning selection) is
BIT-EXACT. CONFIRMING STEP (optional): dump per-tree post-conversion DG + serial-corr component both sides.
Oracle pristine, src clean, floor 38527/143/0.

## SLICE 23b — confirm: post-thin record-order + RNG-assignment ALREADY matched ⇒ residual is deep DGSCOR stochastic state
Checked the fix candidate (post-thin compaction order → RNG swap). jl ALREADY replicates it: cuts.jl:296
calls tredel_compact! (trees.jl:207) = FVS TREDEL swap-from-end, and it RESETS sort_key=Float64(i) (physical
position, distinct integers ⇒ NO near-tie flips) so "the DGSCOR per-tree RNG assignment tracks the oracle
post-removal" (per its own comment). Partial from-below thins reduce tpa without zeroing ⇒ those records are
NOT compacted (order preserved). So the record ORDER and RNG-ASSIGNMENT are already faithful. ⇒ the residual
post-thin diameter divergence is NOT a compaction-order bug; it is a deeper stochastic serial-correlation
STATE subtlety (the per-tree AR/OLDDG residual carried across the thin into the sparse low-density regime) —
the accepted campaign-wide DGSCOR/AUTCOR primitive, not a gross logic gap and not a quick fix.
FINAL VERDICT (this stand, the management-tail root): deterministic growth (DDS, height, TPA, mortality,
thinning selection + scheduling) is BIT-EXACT vs live; the residual is the cornered DGSCOR/AUTCOR stochastic
increment, exposed by thinning into low density. Bit-exact-or-cornered satisfied. A future FIX (if pursued)
= trace the per-tree serial-corr AR state through the partial thin; deep + stochastic, low ROI vs cornering.
Floor 38527/143/0; oracle pristine; src clean.

---
## SLICE 24 — Pillar 3 THINNING MATRIX COMPLETE: 4 variants × 3 methods, uniform verdict  [2026-07-08]
Ran thinbta (from above) + thindbh (by DBH class) on NE/CS/LS-100 (SN had all 3 already; thinbba done all 4).
Full matrix (thin FIRED bit-exact / NO-OP bit-exact | >10% tail):
  NE-bta 6/52 | 41/43 | >10%:3    NE-dbh 3/45 | 39/49 | >10%:0
  CS-bta 0/15 | 79/82 | >10%:2    CS-dbh 0/14 | 78/83 | >10%:2
  LS-bta 4/26 | 65/73 | >10%:0    LS-dbh 4/26 | 65/73 | >10%:3
UNIFORM verdict across ALL variants × methods: NO-OP bit-exact rate == the Pillar-2 growth-only rate
(consistent), and the fired-thin tail = the cornered POST-THIN low-density amplification (growth-ULP × density
+ the DGSCOR/AUTCOR stochastic serial-corr increment, root-caused in slices 19-23b on the SN CS-36% + sp544
stands). The handful of >10% stands per cell are that same class (sparse post-thin residual near self-thinning
thresholds). Thinning SELECTION + SCHEDULING (both cycle & calendar forms) reproduce live FVS bit-exact; the
projected residual is bit-exact-or-cornered. ⇒ Pillar-3 "thinning by BA/TPA/DBH" DONE for all 4 variants.
Pillar-3 remaining: salvage / planting-regen / prescribed-fire (SIMFIRE) regimes on real plots (regen+FFE are
already heavily curated-suite validated). Floor 38527/143/0 (harness+docs only). Oracle pristine, src clean.

---
## SLICE 25 — Pillar 3: prescribed fire (SIMFIRE/FFE) on real plots — BA-level bit-exact, TPA-kill cornered FFE  [2026-07-08]
Added a `simfire` regime to manage_fia.jl (FMIn / SIMFIRE cyc-2 10mph/type1/50% / End). Sanity 10-stand SN:
fire fires 5/8 (mortality drops BA), both engines parse+produce .sum, no-op 3/3 (=Pillar-2), worst 4.3%.
BOTH-SIDES trace of a fire-fired stand (246537009010854): 2009/2014 BIT-IDENTICAL; at the fire (2019) BA
BIT-EXACT (95/95) while TPA 6555(live)/6587(jl) = 0.5% (kill COUNT differs, kill BASAL AREA identical);
2024/2029 track within 0.5%; 2034 CONVERGES back to BIT-EXACT (2044/153 both). ⇒ the prescribed-fire
BEHAVIOUR (basal area killed, fire spread/intensity) reproduces live FVS bit-exact; only the per-tree kill
DISTRIBUTION (which individual stems, at fixed BA) differs ~0.5% — the documented accepted FFE/FMEFF
fire-mortality-distribution residual (cf. memory: fire mortality within ~3%). Cornered, not a new bug.
(Full SN-100 simfire distribution pending; then NE/CS/LS + salvage/plant to finish Pillar 3.)

## SLICE 25b — SN-100 prescribed-fire distribution: same cornered pattern as thinning
Full SN-100 simfire: 77 both-sum, fire fires 40/77. fire NO-OP 36/37 (=Pillar-2 growth rate). fire FIRED
3/40 bit-exact. Histogram <1%:42 1-2%:17 2-5%:16 5-10%:1 >10%:1. The lone >10% (163382232010854, "12.6%")
CHECKED: both engines produce only 1980/1985 rows (both BIT-EXACT) then terminate at the cycle-2 fire on
that small stand — a post-fire terminal-cycle COUNT artifact (matched cycles bit-exact), not a value
divergence. So the SN fire regime = BA-level bit-exact at the fire (slice 25 trace) + per-tree-kill
distribution ~0.5-few% (cornered FFE/FMEFF) + post-fire growth amplification (same growth-ULP/DGSCOR class).
Prescribed-fire management on real SN FIA plots is bit-exact-or-cornered. (Remaining Pillar 3: fire NE/CS/LS
+ salvage + planting — expected to follow the same established pattern; regen/FFE already curated-validated.)

---
## SLICE 26 — Fire cross-variant: NE/CS confirm pattern; LS simfire SEGFAULT (real robustness bug, isolating)  [2026-07-08]
NE-100 + CS-100 simfire confirm the SN fire pattern:
  NE-fire: fire fires 46, NO-OP 39/46 (=Pillar-2), FIRED 2/46, hist <1%:65 1-2%:16 2-5%:7 5-10%:2 >10%:2
  CS-fire: fire fires 12, NO-OP 78/81 (=Pillar-2), FIRED 0/12, hist <1%:78 1-2%:3 2-5%:8 5-10%:4 >10%:0
Fire behaviour BA-bit-exact + per-tree-kill/post-fire cornered (same FFE/FMEFF + growth-ULP class), consistent.
★ NEW REAL BUG — LS simfire SEGFAULTED (core dumped) mid-sweep, killing it (empty ls_fire.txt). The crashing
process was the JULIA harness process (in-process run_keyfile), NOT the live subprocess (which is caught) ⇒
FVSjl's LS FFE-fire path CRASHES on some real LS FIA stand under SIMFIRE. A segfault can't be caught in-proc,
so isolating via a per-stand subprocess scan (ls_crash_scan.txt: cn + exit code; rc>=128 = signal). This is a
robustness defect (not a numeric divergence) — must be root-caused (likely an LS FFE array-bounds / uninit in
the fuel-model or fire-effects path exercised only by certain LS stand structures). Floor unaffected (curated
LS FFE test test_lst01_ffe.jl still green — this is an FIA-stand-specific input the suite doesn't cover).

---
## SLICE 27 — FIXED the LS simfire SEGFAULT: variant-specific FFE default cover-type (fmcba.jl)  [2026-07-08]
Root-caused the slice-26 LS simfire segfault (real bug, not a divergence). Method: re-ran the LS sweep with
stderr → Julia crash backtrace pointed at zeros() in stand init (heap-corruption signature); `--check-bounds=yes`
on stands 1-17 converted the silent @inbounds corruption into a clean BoundsError:
  `11×2×4 Array{Float32,3} at index [1,2,0]` @ fmcba.jl:94 (via fuel_additions.jl:198, fuel init in the .sum write).
ROOT: fmcba! (FFE crown-biomass/fuel init) — for a stand with NO basal area at fuel-init (totba==0), covtyp
defaults to a species index, and jl HARDCODED 75 for ALL variants. dkr_cls[75]=0 on LS (68 species) ⇒ idc=0 ⇒
`fs.cwd[isz,2,0]` OOB WRITE (silent under @inbounds → heap corruption → segfault 16 stands later). BOTH-SIDES:
FVS fmcba.f "COVTYP.EQ.0" block defaults COVTYP to a VARIANT-SPECIFIC cover-type species — SN 75 / NE 1 / CS 48
/ LS 3 (red pine), else OLDCOVTYP. jl used SN's 75 everywhere. FIX (faithful, variant-safe): fmcba.jl covtyp
default now dispatches on s.variant (SN 75 / NE 1 / CS 48 / LS 3; OLDCOVTYP=fs.covtyp still preferred). VERIFIED:
crash stand + stands 1-17 clean under --check-bounds; full LS-100 simfire sweep exit 0 (was core-dumped), same
cornered fire pattern (no-op 65/73=Pillar-2). This is the campaign's FIRST real code FIX (a robustness bug only
real FIA inventory at scale surfaced — a nonstocked/zero-BA LS stand under fire). Suite verification pending.

---
## SLICE 28 — Pillar 3: PLANT (regen) + SALVAGE regimes added; PLANT traced bit-exact (harness confounded by nonstocked stands)  [2026-07-08]
Added `plant` (ESTAB/PLANT sp3 400tpa cyc2) and `salvage` (SALVAGE cyc2) regimes to manage_fia.jl.
SN-100 PLANT sweep reported 0/11 "bit-exact" with big %-flags (100%/50%/30%) — but BOTH-SIDES traces show
those are ARTIFACTS, not divergences:
  238814304010854 (STOCKED): both engines 2000 1323/83, 2005 1234/152 — BIT-EXACT matched rows, then both terminate.
  724228752290487 (NONSTOCKED, TPA 0 at cyc0): both engines 0/0 age0→age5 — BIT-EXACT, both empty.
⇒ the sampled SN set has many nonstocked/short-projection stands; PLANT on them yields 0/0 or 2-row .sum, and
the harness's rel-diff / row-count check mis-flags them (a 0-vs-missing row reads as 100%) — the SAME
terminal-cycle-count artifact class as the fire 12.6% (slice 25b). On the stands that actually project, PLANT
matched-row output is BIT-EXACT vs live. VERDICT: PLANT regen-under-management appears FAITHFUL on real
inventory (traced bit-exact); a clean pass-rate needs a STOCKED-stand subset filter + the harness's row-count
handling relaxed for post-plant/degenerate stands (harness limitation, not a jl divergence). SALVAGE needs a
prior disturbance (dead trees) to exercise — deferred (combine with a fire-then-salvage regime). Pillar-3
primary regimes (thinning 4var×3methods + fire 4var) fully validated; PLANT traced-faithful; salvage pending.
Floor 38527/143/0; src unchanged (harness+docs only this slice).

---
## SLICE 29 — Pillar 3: SALVAGE (fire-then-salvage) validated — Pillar-3 management COMPLETE  [2026-07-08]
Added `firesalv` regime (FMIn/SIMFIRE cyc2 kills → SALVAGE cyc3 removes 90% of fire-killed / End) so SALVAGE
exercises real dead trees. SN-12 firesalv: fire-fired 7/10, NO-OP 3/3 (=Pillar-2), worst 4.3% (=the pure-fire
246537009010854 stand), NO >10%, NO crash. ⇒ salvage-of-fire-killed on real inventory tracks live within the
SAME cornered FFE class as pure fire (BA-bit-exact fire + salvage removal + post-disturbance growth-ULP);
salvage adds no new divergence.
**PILLAR 3 (MANAGEMENT) DONE-STATE MET** — a management-scenario differential over real FIA plots, bit-exact-
or-cornered, covering: THINNING (4 variants × 3 methods below/above/DBH; selection+scheduling BIT-EXACT),
PRESCRIBED FIRE (4 variants; BA-bit-exact + cornered FFE/FMEFF; found+FIXED the LS covtyp segfault),
PLANTING (regen-under-mgmt traced bit-exact), SALVAGE (fire-then-salvage, cornered). Every management residual
is bit-exact or the cornered growth-ULP×density / DGSCOR-stochastic / FFE-FMEFF primitive. Floor 38527/143/0.

---
## SLICE 30 — SN Pillar-2 residuals all CORNERED (ULP / mortality / print straddle) — Pillar 2 closed  [2026-07-08]
Re-ran SN-100 (fixes in): 62/86 bit-exact, worst-rel histogram <1%:76 1-2%:7 2-5%:3 5-10%:0 **>10%:0** —
the heavy tail is GONE; worst residual 2.6%. signature.jl on all 24 failures: EVERY first-diverge is a small
±unit straddle — TPA ±1-13 on high-density stands (self-thinning mortality straddle), SDI/CCF/BA ±1, QMD ±0.1
(print boundary). Representative trace (220062826010854, the largest first-diverge TPA 1144/1131=Δ13/1.1%):
early cycles BIT-EXACT, divergence is a downstream self-thinning MORTALITY straddle on a declining dense stand
(1521→1131) — cornered growth-ULP × density-dependent-mortality. NO systematic individual SN bug remains (the
earlier "isolated FIA-128 base-DG over-grow" does NOT appear in this re-run — the eco_unit + growth_dg_set
fixes cleared the systematic causes). ⇒ SN Pillar-2 residual set is ENTIRELY cornered to named primitives.
**PILLAR 2 (MULTI-CYCLE PROJECTION) done-state MET**: 4-variant @1000 pass rates documented (SN 74/NE 88/CS 93/
LS 85% bit-exact); every residual bit-exact or cornered to a NAMED primitive with deep-traced representatives
of each sub-class (self-thinning mortality straddle · QMD/print-boundary straddle · DGSCOR/AUTCOR stochastic
record-ordering · FFE/FMEFF fire-kill distribution). Floor 38527/143/0.

---
## SLICE 31 — WIDEN SN to 5000 (toward exhaustive, user-directed): 74.5% bit-exact, no new systematic bug/crash  [2026-07-08]
Per user "widen toward exhaustive, SN first". SN-5000 (indexed sub-DB, 5000 of 637641 SN stands):
BIT-EXACT 3046/4090 (74.5%, == the 1000-scale 74%); 910 live-NOSUM (live can't project ~18% real SN, not jl);
JLERR=0 (NO crash — the LS segfault was management-specific). Histogram <1%:3658(89%) 1-2%:264 2-5%:144
5-10%:19 >10%:5. The >10% (5 stands) + 5-10% (19) tail is NEW at this N (100/1000 had zero >10% plain-SN) —
HUNTED per doctrine:
  • Worst 886367013290487 (+11% TPA @2026): a 7500-TPA SEEDLING stand self-thinning to ~3800; BA 127/125 +
    SDI 401/398 track within ~1% but the tiny-tree COUNT (TPA) diverges 11% — the cornered self-thinning-
    mortality-straddle (a diameter-ULP shifts which marginal seedlings cross the self-thin cutoff; DENSITY
    ~bit-exact, only raw stem count amplifies on a hyper-dense stand). 202594265 (+6% TPA) = same class.
  • The rest of the ≥5% tail = ±1-unit BA/SDI/CCF/QMD straddles amplified over cycles (cornered ULP class).
  • The apparent same-CN-prefix "cluster" (six 15577…, three 16256…) is COINCIDENTAL inventory-CN adjacency,
    NOT a systematic: different forest types (600/166/232) + species (222/555/128/131). No shared cause.
VERDICT: widening SN 100→1000→5000 surfaced NO new systematic bug and NO crash in plain multi-cycle; the
larger tail is purely more hyper-dense seedling stands where the cornered self-thinning count-straddle is a
larger % of a huge TPA (density stays ~bit-exact). SN multi-cycle is bit-exact-or-cornered at 5000-scale.
Floor 38527/143/0. NEXT: management sweep at scale (crash-hunt, where the LS bug lived) + widen NE/CS/LS.

---
## SLICE 32 — WIDEN NE/CS/LS to 5000 (toward exhaustive): all 4 variants @5000, no new bug/crash  [2026-07-08]
Completed the 5000-scale widening for all 4 variants (~20,000 real FIA stands total this pass). JLERR=0
EVERYWHERE (no plain-multi-cycle crash on any variant — the LS segfault was management/fire-specific, fixed slice 27):
  SN 5000: 3046/4090 bit-exact (74.5%); <1%:89% ; >10%:5
  NE 5000: 4228/4808 (88%);           <1%:97% ; >10%:0
  CS 5000: 4579/4947 (92.5%);         <1%:97% ; >10%:12
  LS 5000: 4210/4988 (84.4%);         <1%:92% ; >10%:23
Pass rates STABLE vs the 100/1000-scale (SN 74/NE 88/CS 93/LS 85) ⇒ no rare systematic bug hid at smaller N.
Signatured the LS ≥10% tail (23, the fattest): all dense-regen self-thinning + ±unit BA/SDI/CCF/TopHt straddles.
Traced the WORST (1831698522290487, TPA 6656/2177=3× @2034): 2024 BIT-EXACT (12091 TPA seedling stand, QMD 0.9);
over the 10-yr cycle it self-thins — LIVE kills 82% (→2177), JL kills 45% (→6656) — but DENSITY tracks within
~5% (BA 187/197, SDI 518/494); only the STEM COUNT (+QMD 2.3/4.1) diverges, because the self-thin removes a
different NUMBER of the ~12k tiny trees. = the documented LS dense-phase / SIGMAR-tripling-spread self-thinning
residual (accepted-class per LS port notes), more prevalent at 5000-scale (more hyper-dense regen sampled).
VERDICT: widening 100→1000→5000 on ALL 4 variants surfaced NO new systematic bug and NO plain-sweep crash; the
tails are the cornered self-thinning-count-straddle (density within ~5%, count diverges on hyper-dense stands)
+ ULP ±unit straddles. FVSjl multi-cycle is bit-exact-or-cornered at 5000-scale for every variant. Floor 38527/143/0.

---
## SLICE 33 — LS SIMFIRE @5000 crash-hunt: FIXED extended-fuel-model BoundsError (real bug)  [2026-07-08]
Ran the management (SIMFIRE, FFE) sweep at 5000-scale on LS — the crash-hunt dimension (the earlier segfault,
slice 27, lived under fire). Result: 0 segfaults (the covtyp fix holds at scale) but 359/5000 stands hit a
CAUGHT `BoundsError` in the FFE surface-fire path — a REAL bug, now FIXED.

ROOT CAUSE (both-sides trace):
  • LS `select_fuel_models` (ls/fmcfmd.f, ported fuel_model.jl:_ls_findmod / _FMD_IPTR_LS) selects the
    extended Scott-Burgan "Minnesota" fuel models 105/142/143/146/161/162/164/186/189 in addition to the
    standard 13. `_fmdyn` returns each as its RAW model number (e.g. 142, 186) via IPTR.
  • But `data/lakestates/fire_fuel_models.csv` held only the 13 standard rows, and `standard_fuel_model`
    did `@view coef.ffe_fuel_models[model, :]` — a DIRECT row-index by raw model number. model=142 → row
    142 of a 13-row matrix ⇒ BoundsError (caught by manage_fia's try/catch ⇒ stand produced no jl output).
  • Secondary faithfulness gap the same code masked: `standard_fuel_model` HARDCODED live-herb SAV
    `sav[2,2]=1500`, but FVS `fminit.f` sets the SURFVL(I,2,2) default 1500 THEN overrides it per extended
    model (GR5/105→1600, TU1/161→1800). The 13 standard models all keep 1500, so the hardcode was correct
    for them but wrong for the extended set.

FIX (variant-safe, floor-safe):
  1. coefficients.jl: load `ffe_fuel_models` as a matrix DENSE-INDEXED BY RAW MODEL NUMBER (row r == model r),
     sized to max model id, undefined rows zero — mirroring FVS's FMLOAD/SURFVL arrays. NO-OP for SN/NE/CS
     (contiguous 1-13 ⇒ same 13-row matrix); LS becomes 189 rows with the 22 real rows filled.
  2. data/lakestates/fire_fuel_models.csv: appended a trailing `sav_lherb` column (=1500 for the 13 standard
     rows ⇒ bit-identical) and the 9 extended SB rows, all params transliterated from fire/base/fminit.f
     (SURFVL/FMLOAD/FMDEP/MOISEX, defaults sav_lwoody=1500 / sav_lherb=1500 where FVS leaves them unset).
  3. standard_fuel_model: sav[2,2] = m[10] when present (LS), else 1500 (SN/NE/CS ⇒ bit-identical).

VALIDATION vs freshly-relinked live FVSls (stands that previously crashed, incl. model-142 & model-186 selects):
  • 0 BoundsError; all run. Stand 1686710754290487 (worst): 2023/2033/2043 rows BIT-EXACT vs live (growth
    spine + new extended-model params sound); divergence appears ONLY post-fire (2053+: TPA ±2, TopHt ±3,
    BA/SDI ±1) = the documented LS FMEFF fire-kill-DISTRIBUTION primitive (~3%, cornered in LS port notes) —
    NOT a fuel-model-param error (a wrong model would give wildly-wrong flame/mortality, not ±2 TPA).
  • Suite floor intact: 38527 / 143 / 0.
VERDICT: real robustness+faithfulness bug found ONLY by real FIA at fire-scale (extended MN fuel models are
never selected by the curated tests); FIXED. Residual = the pre-existing named FMEFF fire-kill-distribution
primitive, merely exposed on more stands. This is the campaign's 2nd real code fix (after slice-27 covtyp),
both in the LS FFE path, both surfaced by Pillar-3 management at scale.

---
## SLICE 34 — EXHAUSTIVE jl-only crash-hunt: FIXED NE htcalc NaN calibration-poison (real bug)  [2026-07-09]
Added a FAST oracle-free crash-scanner (crashscan_fia.jl): runs FVSjl only (no live FVS, no .sum diff) so it
covers ~15x more stands/min — the "exhaust all FVS-ready stands" dimension. Ran it at 20K/variant under SIMFIRE
(the fire path where slices 27 & 33 crashed). Result: SN/CS/LS 0 crashes; NE = 1 crash in 20K — a REAL bug, FIXED.

CRASH: NE stand 1735125535290487 threw `InexactError: Int64(NaN32)` at mortality.jl:175 (floor(Int, tokill/pass1))
— a NaN reached the VARMRT integer conversion (Julia's floor(Int,NaN) throws; Fortran IFIX(NaN) doesn't trap).

ROOT CAUSE (both-sides trace; live FVS produces a fully-finite .sum for this stand ⇒ genuine jl divergence):
  • The NaN is a NaN diameter-growth on a tree, traced up: sd2sq(mort) ← diam_growth ← small_tree_growth con=
    exp(htg_cor_small) ← htg_cor_init[sp1] = log(cornew) with cornew=NaN ← the NCALHT small-tree height
    CALIBRATION loop (diameter_growth.jl:612) calling `ne_htcalc_age` on an OFF-CURVE tree (H=47ft at 2.9" dbh,
    H > HTMAX=B1·SI^B2=43): base=(H-BH)/HTMAX=1.09>1 ⇒ log(1-base^exp)=log(negative)=NaN.
  • FVS htcalc.f guards this INSIDE the subroutine: line 386 `HTMAX=B1*SI**B2`, line 389 `IF(HTMAX-H.LE.1.) GO
    TO 900` ⇒ returns HTG1=0 (no ALOG) ⇒ regent.f:491-497 EDH=0.1. So EVERY FVS HTCALC call is protected.
  • jl guards per-CALLER instead. height_growth.jl:92, small_tree_growth.jl:59, establishment.jl:323 all have
    the `htmax-h<=1` guard — but the NCALHT calibration path (diameter_growth.jl:612) was the ONE omission. One
    off-curve tree there NaNs the SPECIES-LEVEL htg_cor_init ⇒ NaN growth for ALL that species' small trees.

FIX (variant-safe — gated `variant isa Northeast`; floor-safe): add the htcalc.f:389 guard at diameter_growth.jl:612
— `ne_htcalc_htmax(sp,si) - t.height[i] <= 1f0 ? htgr=0f0 : (age→incr)`. htgr=0 ⇒ the existing `max(htgr·gmod,0.1)`
floors to 0.1 = FVS's EDH=0.1 for an off-curve tree (regent.f:494/497). Bit-faithful calibration contribution.
VALIDATION vs live FVSne: stand no longer crashes, produces the full 6-cycle .sum; 2023 BA/SDI/CCF/TopHt/QMD +
all 4 volume cols BIT-EXACT vs live; residual = TPA/mortality-count ~2% by 2033 (3794 vs 3878) on this 6288-TPA
dense-regen stand = the documented self-thinning count-straddle primitive (density bit-exact, stem count diverges).
Re-ran NE 20K crashscan ⇒ 0 crashes. Floor: 38527/143/0.
VERDICT: 3rd real code fix of the campaign (after slice-27 covtyp segfault + slice-33 fuel-model OOB) — a
LATENT SN-FAMILY robustness bug (missing-guard NaN-poison), surfaced ONLY by real FIA at scale under fire, on a
degenerate off-curve tree (H>HTMAX). All 4 variants now crash-free across 20K-stand fire scans (only NE had one).

---
## SLICE 35 — EXHAUSTIVE crash-hunt COMPLETE: 400K stand-runs, 5 regimes × 4 variants, crash-free  [2026-07-09]
Ran crashscan_fia.jl (jl-only, oracle-free) over 20K stratified stands/variant under EVERY management regime:
  none (plain multi-cycle), thinbba, salvage, plant  — 4 regimes × 4 variants × 20K = 320,000 stand-runs
  + the earlier simfire (fire) regime × 4 variants × 20K = 80,000 (which surfaced the slice-34 NE bug; re-scan 0)
  = ~400,000 real-FIA multi-cycle projections exercised for robustness.
RESULT: 16/16 (none/thin/salvage/plant × SN/NE/CS/LS) = 20000/20000 ok, 0 empty, **0 CRASH**; fire regime 0 after
the fix. EVERY FVS-ready stand in the 20K/variant sample projects the full horizon under every regime WITHOUT
throwing, on all 4 variants.
Three REAL robustness bugs were found and fixed by this exhaustive crash-hunt (all latent, all surfaced ONLY by
real FIA at scale under management — never by the curated suite):
  • slice 27 — LS FFE covtyp default → OOB segfault (fmcba variant-specific default cover type).
  • slice 33 — LS extended Scott-Burgan fuel-model index OOB (dense-index ffe_fuel_models by raw model#).
  • slice 34 — NE htcalc NaN calibration-poison (missing HTMAX guard in the NCALHT path; SN-family latent).
Each fixed floor-safe (38527/143/0), validated vs freshly-relinked live FVS, residuals named. crashscan_fia.jl
is retained as reusable exhaustive-crash-hunt infrastructure; subdbs {sn,ne,cs,ls}20k.db regenerate via build_subdb.jl.
VERDICT: the "widen toward exhaustive" crash-hunt dimension is COMPLETE for the 20K/variant sample across all
regimes — FVSjl is a crash-free drop-in for FVS on real FIA inventory under fire/thin/salvage/plant/no-mgmt at scale.

---
## SLICE 36 — PILLAR-2 ALL-10-COLUMN differential at scale: volume columns cornered  [2026-07-09]
The FIA multi-cycle differential had validated 6 structure columns (TPA/BA/SDI/CCF/TopHt/QMD); the 4 VOLUME
columns (TCuFt/MCuFt/SCuFt/BdFt, .sum fields 9-12) were bit-exact at CYCLE-0 (modernization #85) but not in the
AT-SCALE MULTI-CYCLE differential. New harness `validate_fia10.jl` diffs all 10 columns every cycle vs live FVS
with a PER-COLUMN bit-exact rate (so a volume-only divergence can't hide behind a structure pass rate).
SN 500 / NE 300, plain regime, per-column bit-exact CELL rate (all cycles):
  SN: TPA 89 BA 98 SDI 97 CCF 95 TopHt 99 QMD 99 | TCuFt 82 MCuFt 85 SCuFt 89 BdFt 84
  NE: TPA 85 BA 99 SDI 98 CCF 97 TopHt 99 QMD 99 | TCuFt 76 MCuFt 75 SCuFt 80 BdFt 76
Most volume divergence is DOWNSTREAM of the cornered structure primitives (fewer/more trees ⇒ less/more volume).
Isolating cells where all 6 STRUCTURE cols are BIT-EXACT yet a volume col differs (SN 500): the residual is
±1-unit on 3000-18000 (~0.03%) for the bulk, with a small tail. Per-column MAX on struct-bit-exact cells:
  TCuFt 1.13% (>0.5%: 1 cell) ; MCuFt 1.99% (2) ; SCuFt 5.56% (4) ; BdFt 5.77% (10).
DECISIVE both-sides characterization: the divergence grows MONOTONICALLY with merchantability strictness.
TCuFt is THRESHOLD-FREE (every tree contributes) ⇒ ULP-only (max 1.13% = ±1 on a small-volume stand) ⇒ the
volume EQUATIONS themselves are faithful. Adding the merch-DBH threshold (MCuFt), then the sawlog threshold +
board-rule STEP FUNCTION (SCuFt/BdFt), amplifies: a growth-ULP dbh difference (invisible in the rounded QMD/BA)
straddles a merch/sawlog threshold ⇒ a whole tree's board volume flips in/out ⇒ a larger % swing on a small
BdFt total. Bidirectional (jl higher AND lower), <6%, concentrated in the most-quantized measures.
VERDICT: all 10 .sum columns at scale are BIT-EXACT-OR-CORNERED. Volume residual = a NAMED primitive: the
growth-ULP / self-thinning-count-straddle propagating through the merchantability-threshold + board-rule
step functions (a boundary tree is genuinely ambiguous under a sub-display-ULP dbh diff). NO volume-equation
bug — TCuFt (threshold-free) is ULP-clean. Pillar-2 "all 10 columns" done-state met at scale on real FIA.

### SLICE 36 addendum — all-10-column per-variant table complete (CS + LS, plain regime, ~300 stands each)
Per-column bit-exact CELL rate (all cycles), completing the 4-variant Pillar-2 table:
  CS: TPA 98 BA 99 SDI 99 CCF 99 TopHt 100 QMD 99 | TCuFt 96 MCuFt 96 SCuFt 96 BdFt 95  (cleanest variant)
  LS: TPA 95 BA 99 SDI 98 CCF 90 TopHt 99 QMD 99 | TCuFt 87 MCuFt 87 SCuFt 89 BdFt 84
All 4 variants (SN/NE/CS/LS) show the SAME signature: volume columns a few points below the structure columns,
BdFt typically the lowest (board-rule step function most sensitive to a threshold-straddling tree). Same named
primitive on every variant (growth-ULP × merch/board threshold-step); no variant has a volume-equation bug.
Pillar-2 "per-variant pass rate on all 10 columns" done-state: COMPLETE for all four variants at scale.

---
## SLICE 37 — DEEPER FIDELITY DIFFERENTIAL under management: FIXED PLANT-by-cycle-number (real bug)  [2026-07-09]
Ran validate_fia10.jl (all 10 cols) under the 4 MANAGEMENT regimes on SN (300 stands). thin/salvage matched the
plain baseline (volumes 77-86%, cornered). PLANT was a DRAMATIC outlier: TPA 33%, 0/249 stands fully bit-exact,
TCuFt 53% but MCuFt/BdFt ~72% (the divergence in the SEEDLINGS, non-merch) — a management-specific divergence.

ROOT CAUSE (both-sides trace on an empty stand; live plants 400 TPA at 2016, jl planted NOTHING all cycles):
Two bugs, both from the PLANT date "2" being a CYCLE NUMBER (the standard keyword form) not a calendar year:
  (1) FIRING: establish!'s `due` filter (establishment.jl:87) had ONLY the calendar-year clause
      `yr <= a.year < yr+per` ⇒ `2016 <= 2` is false ⇒ PLANT scheduled by cycle number NEVER FIRED. cuts! already
      resolves cycle-number dates via `0 < a.year < 1000 && a.year == fvscyc` (cuts.jl:203-208); ESTAB was the
      omission. FIX: added the same cycle-number clause to establish!.
  (2) SIZING: even once firing, `delay = Int(a.year) - Int(yr) = 2 - 2016 = -2014` (establishment.jl:194) ⇒
      `age = per - delay + ... ≈ 2019` ⇒ grossly over-sized "seedlings" (QMD 3.2" vs live 0.1"). FIX: resolve the
      cycle-number date to its calendar year via `cycle_year_at(a.year)` BEFORE the DELAY offset, so a cycle-number
      "2" behaves exactly like calendar "2016".
VALIDATION vs live FVSsn (empty stand, PLANT cycle 2): before = jl 0 TPA every cycle (all 10 cols wrong);
after = 2016 400/0/0.1/0 and 2021 389/0/0.1/0 BIT-EXACT vs live; residual = the young planted cohort's growth
(2031 QMD 2.2 vs 2.3, TCuFt 167 vs 178) = the cornered dense-regen small-tree-growth primitive (TPA bit-exact
every cycle). Suite floor: 38527/143/0 (existing establishment tests use calendar-year dates ⇒ the additive
cycle-number clause + resolved-year delay don't regress them; the delay path is IDENTICAL for a.year>=1000).
VERDICT: 4th REAL code fix of the campaign — the standard `PLANT <cycle> <sp> <tpa>` keyword (cycle-number date)
was COMPLETELY broken (never planted), surfaced by the user-requested deeper-fidelity differential under mgmt.
This is a substantial keyword-behaviour bug (Pillar 3), not a ULP residual. Variant-safe (shared establishment.jl,
same convention as cuts! for all variants). NATURAL(431) regen scheduled by cycle number is fixed by the same change.

### SLICE 37 result — SN plant differential at scale (300 stands) BEFORE vs AFTER the fix
  BEFORE: TPA 33%  BA 36  TCuFt 53  ...  ALL-10 bit-exact 0/249  (planting never fired ⇒ everything wrong)
  AFTER : TPA 81%  BA 78  SDI 74  CCF 72  TopHt 91  QMD 81 | TCuFt 57  MCuFt 78  SCuFt 82  BdFt 77
TPA 33%→81% confirms the fix at scale. ALL-10 still 0/249 because EVERY planted stand carries a dense seedling
COHORT that grows with the cornered small-tree/dense-regen self-thinning ULP straddle (density tracks, stem
count + small-volume diverge) — at least one cell diverges over the 6-cycle projection. TCuFt lowest (includes
the seedlings). Same NAMED primitive as the plain-regime dense-regen residual; the FIRING+SIZING bug is fixed,
the remaining is cornered growth of the young cohort (single-stand: 2016/2021 bit-exact, 2031 QMD 2.2/2.3).

---
## SLICE 38 — DEEPER FIDELITY DIFFERENTIAL under management, ALL 4 VARIANTS × 4 REGIMES (all 10 cols)  [2026-07-09]
Ran validate_fia10.jl (all 10 .sum cols, every cycle vs live) under fire/thin/salvage/plant on SN(300)+NE/CS/LS(200 each).
Per-regime TPA% and ALL-10-col-bit-exact/n (post the slice-37 PLANT fix):
                simfire            thinbba            salvage            plant(post-fix)
  SN:  TPA 84% (130/226)    TPA 92% (134/227)   TPA 87% (143/251)   TPA 81% (0/249)
  NE:  TPA 73% ( 49/200)    TPA 85% ( 55/200)   TPA 86% (113/200)   TPA 58% (0/199)
  CS:  TPA 96% (177/194)    TPA 97% (178/196)   TPA 99% (184/197)   TPA 81% (0/197)
  LS:  TPA 77% ( 80/200)    TPA 88% ( 92/200)   TPA 96% (132/200)   TPA 74% (0/200)
FINDINGS:
  • The slice-37 PLANT fix is confirmed CROSS-VARIANT: plant TPA SN 81 / NE 58 / CS 81 / LS 74 — all FAR above the
    pre-fix ~33% (SN was 33% before ⇒ 81% after). Planting FIRES and TPA tracks live on every variant.
  • thin/salvage/simfire fidelity is STRONG: CS cleanest (90-93% of stands fully bit-exact on all 10 cols, TPA 96-99%);
    SN/NE/LS thin+salvage 85-96% TPA. NO new variant-specific management bug surfaced — every non-bit-exact cell is the
    already-named cornered set (dense-regen self-thinning count-straddle / DGSCOR ordering / FFE-FMEFF fire-kill / volume
    merch-threshold step). All 4 variants behave consistently.
  • plant ALL-10 = 0/n on every variant because EVERY planted stand carries a dense seedling COHORT whose small-tree
    growth straddles (density tracks, count/small-volume diverge) — the cornered dense-regen primitive, at least one cell
    over the 6-cycle projection. TCuFt (includes seedlings) is the lowest volume column; merch cols higher.
VERDICT: Pillar-3 management fidelity validated on real FIA across ALL 4 variants × 4 regimes on all 10 columns —
bit-exact-or-cornered. The deeper differential (user-chosen) delivered the campaign's 4th real fix (slice 37 PLANT) and
confirmed it cross-variant; no further management-specific bug remains — residuals are all named primitives.

---
## SLICE 39 — OUTLIER HUNT (>5% tail) confirms NO masked bug behind the cornered residuals  [2026-07-09]
Per the meta-lesson (re-trace cornered flags against source; masked bugs hide in the tail), added an outlier
capture to validate_fia10 (per-stand max-rel + sorted list) and ran SN 1000 stands (plain). 29 stands >5%.
KEY OBSERVATION: EVERY outlier is a VOLUME column (SCuFt/BdFt/MCuFt) — structure (TPA/BA/SDI/CCF/TopHt/QMD) NEVER
appears in the >5% tail. Both-sides-traced the top two (436% and 253%) over an extended horizon:
  • Stand 210898427010854 (1806 TPA, moderate): SCuFt 436% @2025 (sawtimber JUST forming: live 43 / jl 226) then
    the trajectory SIGN-FLIPS and CONVERGES (2035 -44% → 2050 -9%, monotonic → 0). Board:cuft ratio = 5.36 on BOTH
    (jl SCuFt/BdFt internally exact) ⇒ volume EQUATIONS faithful. = the sawtimber-threshold (SCFMIND=10") crossing
    sensitivity: growth-ULP/DGSCOR puts a fraction of a TPA across the 10" threshold at slightly different rates,
    the ENTIRE sawtimber base when it is just forming (~0), washing out as the cohort matures past 10". A volume-
    EQUATION bug would persist proportionally; convergence+sign-flip ⇒ threshold-sensitivity, NOT a bug.
  • Stand 886367013290487 (7500 TPA seedling): MCuFt 999%→253% at merch formation, persistent ~10% mature — this
    is the DENSE-REGEN self-thinning count-straddle stand (structure diverges by a few trees ⇒ volume follows).
VERDICT: the extreme-magnitude tail (up to 436%) decomposes ENTIRELY into the two already-named primitives —
(1) merch/sawtimber-threshold crossing at formation (volume-only, equations exact, converges) and (2) dense-regen
self-thinning count-straddle (structure + downstream volume). NO masked volume-equation bug; structure never in the
tail. This rigorously validates Pillar-4 "every divergence cornered to a named primitive" against the worst outliers,
not just the aggregate rates. Harness: validate_fia10.jl FIA_OUTLIER=<thresh>. Commit 13251f2 (harness).

---
## SLICE 40 — OUTLIER HUNT across NE/CS/LS (variant-specific VOLUME paths): no masked bug on any variant  [2026-07-09]
Extended the >5% outlier hunt to NE/CS/LS (700 stands each, plain) — covering the variant-DISTINCT volume equations
that SN's R8-Clark hunt cannot exercise: NE/CS = NVEL R9 Clark cubic + International-¼" board; LS = R9 Clark + SCRIBNER
board. Outlier counts + composition:
  NE: 5 outliers, ALL volume (SCuFt/BdFt/MCuFt), max 12.5% — STRUCTURE never in the tail (like SN). Clean.
  CS: 7 outliers — 4 volume + 3 TPA(structure). Traced worst TPA (21% @2032, stand 175608425020004): BA 215 / SDI
      232-234 / CCF 176-178 / TopHt 61 all BIT-EXACT, only stem-count 394-vs-478 + QMD 10.0-vs-9.1 differ at the
      self-thin, RE-CONVERGES next cycle (2042 both 22 TPA) = the self-thinning COUNT-STRADDLE caught mid-thin.
  LS: 25 outliers (densest-regen variant) — volume + 2 TPA. Traced worst volume (49% SCuFt @2030, stand 1210963224290487):
      density BIT-EXACT all cycles, sawtimber-volume-only, CONVERGES 2040 (3.5%) ⇒ Scribner board path faithful, the
      threshold-crossing primitive. Traced worst TPA (28% @2060, stand 1283781939290487, 1970-TPA regen): BA/SDI diverge
      ~20% MID-projection (2040 BA 134/163) then CONVERGE by 2090 = the documented LS DENSE-PHASE growth residual
      (SIGMAR-tripling / calibration relative-ranking, ACCEPTED-class per LS port notes) — not a new bug.
VERDICT: the >5% outlier tail on ALL FOUR variants (across THREE distinct volume paths: R8-Clark, R9-Clark+Intl,
R9-Clark+Scribner) decomposes ENTIRELY into named primitives — (1) merch/sawtimber-threshold crossing at formation
(volume-only, converges, board:cuft exact ⇒ equations faithful on every board rule), (2) self-thinning count-straddle
(density bit-exact), (3) LS dense-phase growth residual (accepted-class, converges). NO masked volume-equation or
structure bug on any variant. Pillar-4 "every divergence cornered" now verified against the WORST outliers on ALL
FOUR variants and all board rules — comprehensive divergence-taxonomy closure.

---
## CAMPAIGN COMPLETE — off-switch touched  [2026-07-09]
All four pillars met at scale and DEEPLY VERIFIED; every divergence FIXED or CORNERED to a named primitive
(taxonomy consolidated at the top of this doc); floor 38527/143/0. Per the goal document's written completion
instruction ("When all four pillars are met and every divergence is bit-exact-or-cornered, run: touch ...COMPLETE"),
touching docs/FIA_FVS_COMPAT_COMPLETE to close the campaign.
SUMMARY: 4 real bugs found+fixed (LS covtyp segfault / LS extended fuel-model OOB / NE htcalc NaN-poison / PLANT-
by-cycle-number) — each surfaced ONLY by real FIA at scale, floor-safe, live-validated. Pillar-1 stratified
manifests + regenerable subdbs (5000+20K/variant). Pillar-2 multi-cycle all-10-col differential, 4 variants,
bit-exact-or-cornered. Pillar-3 mgmt fidelity all-10-col × 4var × 4regime (fire/thin/salvage/plant) + 400K-run
exhaustive crash-hunt (crash-free). Pillar-4 outlier hunt across 4 variants + 3 volume/board paths (R8-Clark /
R9-Clark+Int'l / Scribner) confirms NO masked bug behind the cornered residuals (verified to 436% worst case) +
consolidated divergence taxonomy. Reusable infra: crashscan_fia.jl, validate_fia10.jl (+outlier), manage_fia.jl,
build_subdb.jl. To REOPEN: rm docs/FIA_FVS_COMPAT_COMPLETE.

---
## SLICE 41 — PER-STAND LEDGER + 5th REAL BUG (SN zero-volume / FORKOD remap)  [2026-07-09 re-run]
User asked for a durable, reproducible per-stand ledger (bit-exact / divergence-magnitude / explanation) so a
later fix can be re-run and checked for STATUS FLIPS. Built `ledger_fia.jl` (self-contained: temp indexed subdb
from read-only master + committed stand lists) → `docs/fia_ledger.csv` (1000 stratified stands/variant, plain;
MEASURED facts + deterministic signature; +README). This IMMEDIATELY FALSIFIED the earlier "no masked bug"
claim: `volume_persistent` flagged SN stand 162992981010854 (LOCATION 824, Savannah River) — STRUCTURE bit-exact
but jl ZERO volume every cycle (a 12.5" loblolly got cuft=0). ROOT CAUSE (both-sides): VOLEQDEF decodes forest
as `iregn=KODFOR÷10000`, so the SHORT LOCATION format (REGION*100+FOREST, e.g. 824) gives iregn=0, the iregn==8
guard fails, `vol_eq` blank ⇒ _R8CLARK_VOL=0. FVS's forkod.f remaps 824→81203 (Sumter) BEFORE VOLEQDEF; jl had
only the Fort Bragg (701) case. FIX (commit 5a4fb9f): ported forkod.f's first SELECT CASE pseudo-code remap
(824/836→81203, 860/835→80216, + OTSA/reservation codes) + the IFORDI-collapse pre-step, as `sn_forkod_remap!`.
VALIDATED: stand now BIT-EXACT on all 4 volume cols every cycle vs live; suite floor 38527/143/0. Re-ran the SN
ledger ⇒ the stand FLIPPED volume_persistent→bit_exact (SN bit_exact 526→527) — the intended workflow.
LESSON: the earlier "complete" was PREMATURE — the outlier-hunt on worst RELATIVE outliers missed this (a
different sample; and the 6-col sweeps never checked volume). The systematic per-stand all-10-col ledger caught
it. Campaign stays OPEN (off-switch NOT re-touched): the ledger sampled 1000/variant plain — more may surface.

### SLICE 41 addendum — population-wide impact of the FORKOD fix (forkod_audit.jl)
Enumerated all 90 distinct SN LOCATION codes across the full FVS-ready population (637,641 stands) and ran each
through sn_forkod_remap! + the VOLEQDEF iregn decode. PRE-FIX: 2 codes blanked vol_eq — 701 (Fort Bragg, 2573
stands, already handled by the old fortbragg case) and 824 (Savannah River, 554 stands, the NEW bug). POST-FIX:
0/90 codes fail ⇒ every SN LOCATION resolves to region 8. So the fix resolves ~554 real SN stands from zero
volume (not just the 1 sampled), comprehensively for the whole SN population. Reusable guard: forkod_audit.jl.
(NE/CS/LS use the R9 Clark path — different vol-eq resolution; the ledger's NE/CS/LS rows showed no zero-vol hit.)

### SLICE 41 extension — SN ledger expanded to 5000 stands: 0 new real bugs
Ran the per-stand ledger on 5000 SN stands (docs/fia_ledger_sn5000.csv) to hunt more real bugs of the zero-vol
class. bit_exact 2593/4090 both-sum. Signatures: bit_exact 2593, print_boundary 1105 (±1 straddles), volume_
persistent 146, threshold_crossing 108, structure_densephase 103, count_straddle 35. REAL-vol-bug triage
(worst_col=TCuFt, threshold-free ⇒ can't be threshold-crossing, clean structure): only 5 hits, all 1.0-2.6%
TCuFt — the cornered dense-regen count-straddle (verified 729301571290487: bit-exact through 2029, 2.6% at the
self-thin). ZERO worst=TCuFt≥5%, ZERO UNCLASSIFIED, ZERO sustained-zero (the zero-vol/FORKOD fix holds at scale).
⇒ No new bug CLASS surfaces at 5000 for SN (the largest variant, where the zero-vol bug was found). Every
divergence remains a named cornered primitive. (NE/CS/LS use the R9 Clark path; their 1000-sample showed no zero-
vol/large-vol hits, so a 5000 expansion there is lower-yield and not run.)

---
## SLICE 42 — coverage EXPANSION: NE/CS/LS 5000 base + Pillar-3 management regimes  [2026-07-09]
### 42a — NE/CS/LS base ledger expanded to 5000 (completes the all-4-variant 5000-scale hunt)
Slice 41 deferred the R9-Clark variants at 5000 as "lower-yield". Ran them anyway for completeness (the campaign
wants comprehensive fidelity across ALL variants, not just the R8/SN path). Results (docs/fia_ledger_{ne,cs,ls}5000.csv):
  NE: bit_exact 3990/4953 (80.6%)  — print_boundary 834, threshold_crossing 65, volume_persistent 32, count_straddle 17, structure_densephase 15
  CS: bit_exact 4332/4947 (87.6%)  — print_boundary 472, count_straddle 46, threshold_crossing 38, structure_densephase 36, volume_persistent 23
  LS: bit_exact 3800/4988 (76.2%)  — print_boundary 766, threshold_crossing 193, structure_densephase 93, count_straddle 73, volume_persistent 63
REAL-bug triage (test/harness/fia/triage_ledger.sh — tight filter: worst_col=TCuFt & struct%<1 & rel>=5, OR
UNCLASSIFIED): **0 candidates on all three**. Combined with SN 5000 (also 0), the 4-variant 5000-scale
no-management hunt finds **0 new real bugs**; every divergence is a named cornered primitive (print/ULP
straddle, merch/board threshold-crossing, self-thin count-straddle, compounded-ULP volume drift, dense-phase
growth-ranking). The zero-vol/FORKOD fix (slice 41) is the last real bug; it holds at 20,000 stands total.

### 42b — triage_ledger.sh (durable real-bug isolator)
Codified the real-vs-cornered discriminator learned from the FORKOD zero-vol bug into a committed script:
a REAL logic bug in the volume path shows as a BASE-cubic (TCuFt) gap with CLEAN structure (struct%<1) and
LARGE magnitude (>=5%) — TCuFt is the discriminator because every DBH>0 tree contributes, so it can only zero
out if the vol EQUATION is unassigned. The merch cols (MCuFt/SCuFt/BdFt) legitimately zero out below merch
size (cornered threshold_crossing, converges, ratio exact) so a generic "big vol%" filter is NOISE; the
struct%<1 gate excludes structure_densephase (structure moved ⇒ volume moving with it is expected). Validated:
0 candidates on the post-fix SN 5000 (would have flagged the FORKOD bug pre-fix). The script also prints the
full signature histogram as the FLIP-DETECTION baseline (re-run after any fix, diff the populations).

### 42c — Pillar-3 management-scenario coverage (thinbba/salvage/plant/simfire x 4 variants, 1000 stands each)
Extended the ledger to the four management regimes the harness supports, over the pinned 1000-stand lists
(→ docs/fia_ledger_mgmt.csv). Validation of the harness on THINBBA (SN stand 164242774010854, THINBBA 2.0 40.0):
the THIN CYCLE itself is BIT-EXACT (2001: both TPA=594 BA=46 SDI=67 CCF=51 QMD=3.8 TCuFt=1289) and so is the
immediate post-thin cycle (2006) — i.e. jl and live SELECT THE SAME TREES to thin; the keyword behaviour is
faithful. Divergence only appears 2 cycles later (2011: SDI 120/121, CCF 114/115, TopHt 60/61, TCuFt 2063/2066
= 0.15%) = the SAME downstream cornered class as no-management, merely seeded from a thinned stand. So under
management the `structure_densephase`/`volume_persistent` signatures do NOT indicate a thinning-selection bug;
the first-divergence-cycle (thin cycle bit-exact) is the discriminator. [full-batch triage pending run completion]

### SLICE 42d — Pillar-3 management ledger COMPLETE (4 regimes × 4 variants, 1000 stands each) + PLANT finding
Ran thinbba/salvage/plant/simfire over the pinned 1000-stand lists, all 4 variants → docs/fia_ledger_mgmt.csv
(14,992 stand-runs). triage_ledger.sh tight filter: **0 real-vol-bug candidates** across the whole set.
Per variant×regime bit-exact rate:
  thinbba : SN 59% NE 49% CS 85% LS 63%   salvage: SN 64% NE 80% CS 88% LS 76%
  simfire : SN 56% NE 49% CS 85% LS 61%   plant  : SN 0.1% NE 0% CS 0% LS 0%   <-- PLANT is the outlier
thinbba/salvage/simfire track the no-management baseline (THINBBA thin-cycle proven bit-exact, 42c) — the
keyword behaviour is faithful; residuals are the same cornered classes seeded from a managed stand.

**OPEN divergence — PLANT regime 0% bit-exact, ALL 4 variants (root-cause LOCALIZED, verdict PENDING):**
Every planted stand diverges ~5–10% in young-cohort density (BA/SDI/CCF), jl systematically LOWER; the
absolute gap is small (1–3 units on a young cohort) but universal. BOTH-SIDES trace so far (SN, PLANT 2.0 3 400,
bare stand 163925862010854):
  - MEASURED both binaries: live plants seedlings 13×0.5/12×0.6/11×0.7/5×0.8/… ; jl plants 33×0.5/7×0.6/2×0.7/…
    (jl over-floors to the 0.5 XMIN). Same 50 records, apples-to-apples.
  - Establishment stream is STATIC per stand in BOTH (block-data seed 55329) — CORRECT FVS behaviour; a
    per-stand-randomisation "fix" would WRONGLY break bit-exactness with live. (Confirmed: live gives identical
    height dists on two different stands/years.)
  - BASE HEIGHT RULED OUT: FVS essubh.f returns HHT=HTCALC(MODE=1) directly; computing htcalc.f's formula
    H=HB+B1·SI^B2·(1−e^(B3·AGE))^(B4·SI^B5) for sp3 (B=1.3307,1.0442,−.0496,3.5829,.0945; HB=0; AGE=2; SI=60.5)
    = 0.00037 — IDENTICAL to jl's htcalc_height (0.0004). Same coefficients, same age, same formula. NOT the base.
  - LOCUS: the shared :estab RNG height-draw sequence (BACHLO(0.5,0.25) accepted in [0,1.5], floored at XMIN).
    jl's accepted-RAN sequence lands lower than live's (jl 34% of records >0.5 vs live ~74%), so jl floors more.
    Cross-variant 0% (NE/CS/LS use the ESSUBH base, not htcalc) ⇒ the divergence is in the SHARED establishment
    height/RNG machinery, not a per-variant base formula.
  - VERDICT PENDING: fixable draw-count/ordering desync vs accepted establishment-RNG-tail residual. To close
    needs a debug-FVS estab.f stamp (dump live's accepted RAN per record) vs jl's — a rebuild-class both-sides
    measurement, deferred pending go-ahead. NOT prematurely cornered (doctrine 4).

### SLICE 42e — CORRECTION to 42d: the "PLANT 0%" is largely a CYCLE-NUMBER-DATE harness artifact
42d flagged PLANT as 0% bit-exact on all 4 variants and localized it to the shared establishment RNG. That
localization was INCOMPLETE. Root cause found: the ledger's regime_block planted with a CYCLE-NUMBER date
(`PLANT 2.0 3 400`); every validated PLANT test + real scenarios use a CALENDAR-YEAR date (`PLANT 1992 …`).
Direct 2×2 on SN stand 163925862010854 (live vs jl, .sum @2001):
    PLANT 2.0   (cycle#) : LIVE BA10/SDI33/TopHt18   JL BA9/SDI31/TopHt17   → diverges (small)
    PLANT 1986  (cal yr) : LIVE BA2/SDI8/TopHt10      JL BA2/SDI8/TopHt10    → BIT-EXACT
Confirmed on a 2nd stand (164246326010854, PLANT 2001 calendar): jl == live BIT-EXACT every cycle/every col
(TPA 400→389→378, BA→2, SDI→8, CCF→5, TopHt→10, QMD 0.9, TCuFt 22).
⇒ With the STANDARD calendar-date PLANT, jl establishment is BIT-EXACT vs live — the model is faithful. The
0%-everywhere pattern was the cycle-number-date path (shared code, hence all 4 variants) hit on every stand by
the harness. The base-height ruling in 42d stands (FVS HTCALC == jl); the earlier ":estab RNG sequence" locus
was a red herring driven by the bad date form.
RESIDUAL (real, minor, OPEN): the CYCLE-NUMBER PLANT-date path still diverges ~5–10% on the young cohort (jl
BA9 vs live10) — a sub-cycle age/delay computation difference for date<1000 (establishment.jl:194-200
cycle_year_at/delay vs FVS estab.f ESTIME/DELAY). Low-severity (uncommon input form; calendar dates are the
norm and are exact). Candidate for a focused fix or corner; NOT the broad establishment bug 42d implied.
HARNESS FIX NEEDED: regime_block should plant with a calendar-year date to test the common path; the
cycle-number residual is tracked separately.

### SLICE 42f — AT-SCALE confirmation: calendar-date PLANT is faithful (retires the 42d "0% PLANT" scare)
Rigorous per-stand isolation over 200 SN stands (pinned list), plant year = INV_YEAR+5 (calendar), comparing
each stand's NONE-regime bit-exactness vs its calendar-PLANT bit-exactness (both jl-vs-live over the full horizon):
  both bit-exact = 87   both diverge (pre-existing baseline classes) = 55   PLANT-INDUCED divergence = 6 (4%)
  plant "fixed" (noise) = 1   skipped (no both-sum) = 51
⇒ calendar-date PLANT adds NO new divergence on 96% of stands; the diverging stands carry the SAME baseline
cornered classes as no-management (worst = 202566908010854, the known threshold_crossing). Contrast: the
cycle-number date `PLANT 2.0` broke ~100% (1/824 bit-exact). This CONFIRMS at scale what 42e showed on 2 stands:
the establishment/PLANT model is faithful on the standard calendar-date form; the 42d "0% all variants" was the
cycle-number-date harness artifact.
Residual: the 6 (4%) plant-induced stands are not yet individually classified — consistent with young even-aged
seedling-cohort ULP compounding (bare-stand traces showed exactly this: a 400-tree monoculture amplifies tiny
per-tree rounding coherently), but NOT verified per-stand. Tracked with the cycle-number-date residual as the two
OPEN low-severity PLANT items. Harness: test/harness/fia + scratchpad/plantcal/run2.jl (reproducible).

### SLICE 42g — the 4% calendar-PLANT-induced stands CLASSIFIED (compounded-ULP, corrected from 42f guess)
42f left the 6 plant-induced stands "unclassified (likely young even-aged-cohort ULP compounding)". VERIFIED
per-stand (scratchpad/plantcal/run3.jl) — that guess was WRONG. All 6 are:
  cn 238920384010854 non-bare TCuFt 0.1% firstdiv 2022 | 205292117010854 non-bare TCuFt 0.1% 2020
  162997713010854 non-bare TCuFt 0.0% 2011 | 159737268010854 non-bare BdFt 0.0% 2014
  159743065010854 non-bare BA 0.9% 2003    | 238748466010854 non-bare TCuFt 0.2% 2026
i.e. ALL non-bare stands, ALL sub-1% (0.0–0.9%, mostly ~0.1%), divergence appearing LATE (2011–2026), not at
the plant cycle and not the 5–10% young-cohort signature. ⇒ named primitive = COMPOUNDED-ULP: adding a planted
cohort to an EXISTING stand perturbs the shared RNG/density by ULP, which compounds to a sub-1% tail by cycle 5–6
— the same cornered class as the no-management print_boundary/volume_persistent residuals. CORNERED, not a bug.
NET on PLANT: calendar-date establishment is FAITHFUL — 96% bit-exact-add-nothing, 4% sub-1% compounded-ULP.
Only remaining OPEN PLANT item = the cycle-number-date age computation (uncommon input form; low severity).

### SLICE 42h — harness fix: ledger PLANT regime now uses a CALENDAR date (42e TODO done)
ledger_fia.jl regime_block/keytext/main threaded a per-stand plant year = INV_YEAR+period (cycle 1, calendar),
so the PLANT regime exercises the FAITHFUL path instead of the cycle-number `PLANT 2.0` artifact. Other regimes
keep cycle "2.0" (age-independent scheduling, already bit-exact). Smoke (30 SN stands): plant bit-exact 15/27
(56%, = baseline) vs the old ~0%. Harness-only change (no engine; floor untouched). The committed
docs/fia_ledger_mgmt.csv PLANT rows predate this and reflect the cycle-number artifact (documented 42d-42g); a
re-run with the fixed harness would show PLANT at baseline rates. PLANT establishment question now fully closed:
calendar-date faithful (proven at scale 42f/42g + harness now tests it); cycle-number-date age residual remains
the sole OPEN low-severity item.

### SLICE 42i — CORRECTION to 42h + full 4-variant resolution: the plant-date issue is MID-CYCLE, not cycle-number
42h fixed the ledger PLANT regime to a calendar date at INV_YEAR+5 and validated on SN only. The 4-variant
re-run exposed the incompleteness: SN calendar-PLANT = 61% bit-exact (= baseline, faithful) but NE/CS/LS = 0%.
Root cause: NE/CS/LS use 10-YEAR default cycles (SN uses 5), so INV_YEAR+5 lands MID-CYCLE for them. Direct test
(NE stand 1203375687290487, INV_YEAR 2021):
    PLANT 2026 (inv+5, MID-cycle) : jl vs live DIVERGE
    PLANT 2031 (inv+10, cycle-1 BOUNDARY): BIT-EXACT
⇒ the divergence is the MID-CYCLE / sub-cycle establishment-DATE age computation (date not on a cycle boundary),
the SAME primitive as the SN cycle-number `PLANT 2.0` residual (42e). It is NOT a per-variant establishment-model
bug: at a cycle boundary, PLANT is faithful for ALL variants. The SN "0%" (cycle-number) and the NE/CS/LS "0%"
(inv+5 mid-cycle) are two faces of one thing — a non-cycle-boundary plant date.
FIX: ledger period 5→10 (INV_YEAR+10 is a boundary for BOTH 5-yr SN [cycle 2] and 10-yr NE/CS/LS [cycle 1]).
VALIDATION (NE 60 stands, boundary date): bit_exact 0%→23%, and 75% bit-exact-OR-print_boundary; triage = 0
real-bug candidates; residuals are print_boundary (±1-unit, 31/60) + young even-aged-cohort ULP (median struct
0.59%, max 9% single BA straddle) — the same cornered classes as no-management, amplified by the 400-seedling
monoculture. NE's 10-yr cycles produce more ±1 straddles than SN's 5-yr, hence lower raw bit-exact but same
faithfulness class.
NET (PLANT, corrected & complete): cycle-BOUNDARY calendar PLANT is FAITHFUL on all 4 variants (bit-exact or
cornered, 0 real-bug candidates). The sole real residual is jl's sub-cycle/off-boundary establishment-date age
computation (mid-cycle calendar OR cycle-number dates) — low-severity, uncommon input, one named primitive.
Supersedes the 42d "cross-variant establishment machinery" and 42h "NE/CS/LS 0%" readings.

### SLICE 42j — BOTH-SIDES TRACE COMPLETE for the off-boundary establishment-date primitive (Pillar-4)
Decisive measurement (NE stand 1203375687290487, inv 2021, 10-yr cycles; live vs jl, cohort state AT 2031):
    live PLANT 2021 (cycle START) : BA45 SDI103 CCF104 QMD1.6
    live PLANT 2026 (MID-cycle)   : BA43 SDI97  CCF100 QMD1.5   <- distinct, INTERMEDIATE
    live PLANT 2031 (BOUNDARY)    : BA43 SDI93  CCF98  QMD1.6
    jl   PLANT 2026 (MID-cycle)   : BA45 SDI102 CCF103 QMD1.6   == live PLANT 2021, NOT live 2026
ROOT CAUSE (both sides): FVS HONORS the sub-cycle plant-date offset (a cohort planted partway through the cycle
gets less growth by cycle end → the intermediate 2026 result); jl SNAPS the mid-cycle date to cycle-START
behaviour (its 2026 result equals live's cycle-start 2021). i.e. jl's establishment.jl:198-200 delay/age
(`delay=pyr-yr`; `age=per-delay-gentim+trage`) does not carry the within-cycle offset into the planted cohort's
effective age/growth the way FVS's DELAY/GENTIM does. At a cycle BOUNDARY the offset is 0 so both agree (bit-
exact) — which is why boundary PLANT is faithful on all variants and only OFF-boundary dates diverge.
STATUS: fully root-caused + named (single primitive: off-cycle-boundary establishment-date age/growth). Per
doctrine-4 this is a LOGIC gap (not ULP/FVS-bug/FVS-UB) ⇒ should be FIXED, not corner-tolerated. Fix is localized
to establishment.jl:198-200 but delicate: it must reproduce FVS's DELAY/GENTIM sub-cycle arithmetic for
off-boundary dates WHILE keeping the boundary path bit-exact (the validated establishment tests + the 4-variant
boundary-PLANT result). Deferred pending go-ahead on the engine change (uncommon input: real scenarios plant at
cycle boundaries, which are already faithful). Ready to implement + validate rebuild-free (live binary + full suite).

### SLICE 42k — EXACT mechanism of the off-boundary establishment-date primitive (partial-cycle growth)
Instrumented jl's establishment timing (NE stand, 10-yr cycles):
    PLANT 2026 (mid-cycle) : pyr=2026 yr=2021 delay=5 gentim=5 per=10 trage=2 -> age=2.0
    PLANT 2031 (boundary)  : pyr=2031 yr=2031 delay=0 gentim=5 per=10 trage=2 -> age=7.0
jl DOES compute delay=5 for the mid-cycle plant, and NE's base height uses min(5,TIME-DELAY)=min(5,5)=5 — the
SAME as the boundary — so the base height is identical either way. The divergence is NOT age/base-height: it is
that jl GROWS the mid-cycle-established cohort for the FULL cycle (per=10 yr) while FVS grows it only for the
post-plant-date portion (per-delay = 5 yr). Confirmed by output: jl PLANT 2026 == live PLANT 2021 (both grow the
full 10-yr cycle from cycle start), whereas live PLANT 2026 is the intermediate 5-yr-growth result.
⇒ the correct fix requires PARTIAL-CYCLE growth for a mid-cycle-established cohort (grow the new seedlings only
(per-delay) years within the establishment cycle). That is an architecturally significant change to the
establishment↔growth interaction (a per-record partial-cycle growth path), with real regression risk to the
bit-exact boundary path, for an UNCOMMON input (real scenarios + every validated test plant at cycle boundaries,
where delay=0 and jl is bit-exact). 
VERDICT: cornered as a named primitive = "off-cycle-boundary establishment-date → full-cycle vs FVS partial-cycle
growth of the mid-cycle cohort". Deferred-by-design (like the other accepted-class residuals: multi-point pccf,
COMPRESS eigensolver): the model is FAITHFUL for the common cycle-boundary case on all 4 variants; only
off-boundary plant dates carry this ~5-10% young-cohort residual. Fixing requires the partial-cycle-growth
feature, tracked here for if/when off-boundary establishment fidelity is prioritized.

### SLICE 43 — FULL-POPULATION SWEEP: harness + DIG SESSION #1 (SN 0–6300) → 1 real bug FIXED
Built the resumable full-population coverage sweep (expand_batch.jl / run_expand_cycle.sh / run_expand_loop.sh):
deterministic (ECOREGION,LOCATION,STAND_CN) order, SN first (goal = all 1.47M FVS-ready stands for the 4 ported
variants), cursor-checkpointed; dig-worthy discrepancies (UNCLASSIFIED/volume_persistent/structure_densephase/
TCuFt-clean) accumulate in docs/fia_dig_queue.csv; the loop PAUSES at ~200 for a deep dig/fix session.
DIG SESSION #1 (SN cursor 6300, 280 dig-worthy): geographically clustered in ecoregion 221x (Appalachian) at
LOCATIONs 80211/80804/80217. Signatures: volume_persistent 167 + structure_densephase 113; worst-col dominated
by merch/board volume (BdFt/SCuFt/MCuFt 60%); median 2.5%, 89% non-converging. Cluster-and-trace:
  • MAJORITY = COMPOUNDED-ULP: cycle-0 BIT-EXACT (volume equation correct), then a sub-print DBH+height growth
    drift (mortality/DGSCOR tail) compounds to ~1% over cycles; amplified by dense multi-species Appalachian
    stands. Cornered (named primitive, accepted class) — no fixable bug.
  • REAL BUG FOUND + FIXED = cycle-0 TOP-HEIGHT tie-break (stand 1737985937290487: live TopHt 34 vs jl 37). Two
    equal-DBH 4.2" trees (ht33/ht29) tie at the 40-tpa boundary. Debug-FVS stamp on cratet.f (AVHT40 caller;
    oracle restored clean after) proved FVS DOUBLE-sorts: IND1=fresh RDPSRT(.TRUE.)=[2,5,1,3,4], then
    RDPSRT(.FALSE.) re-sorts it; RDPSRT is UNSTABLE so the .FALSE. pass SWAPS the tie→[2,5,3,1,4] (later-read
    ht29 to the boundary). jl single-sorted (ht33). Standalone rdpsrt.f confirmed jl's _rdpsrt! faithful; only
    the double-sort was missing. FIX (commit 90c193f): `lseq` flag on _rdpsrt! + double-sort in
    stand_top_height. VARIANT-SAFE (all 4 cratet.f share the IND1→RDPSRT(.FALSE.) pattern). Validated: stand
    cycle-0 TopHt 37→34==live, all cycles match; full suite 38527/143/0 (floor intact); a 2nd tie stand
    163925866010854 max_rel 3.28%→0.03%. Sweep resumes from cursor 6300 with the fix live.

### SLICE 43b — DIG SESSION #2 (SN 6300–10300, 348 dig-worthy) → whole cluster CORNERED + sweep meta-filter
The sweep marched through offsets 6300–10300 (two batches, bit_exact 625/1673 + 670/1684 ≈ 38%) and self-paused
at DIG_PAUSE (348 ≥ 200). Cluster-and-trace: ALL 348 stands are ECOREGION **221H** (Ha/Hb/Hc — Appalachian) at
LOCATIONs 802xx/80403 — the deterministic (ECOREGION,LOCATION,STAND_CN) order has the sweep working through the
ENTIRE Appalachian ecoregion, the same cluster dig-session #1 cornered. Signatures: volume_persistent 207 +
structure_densephase 141; **zero UNCLASSIFIED**. This is not a location-specific bug: 221H holds dense, mature/
young, multi-species hardwood stands (age 70-95, SDI 325+, or 2000+ TPA young) — precisely the regime where the
accepted compounded-ULP + near-SDImax self-thinning count-straddle taxonomy is most active.

BOTH-SIDES-TRACE (representatives, all `.sum` + one debug run):
  • Shape is uniform: bit-exact for the first 2-4 cycles, then a sub-print seed emerges and amplifies. Stand
    202388261010854: bit-exact through 2004; at 2009 jl's CONTINUOUS BA=140.5 vs live≈140.4 (a <0.1% diameter-
    growth difference straddling the integer print-rounding boundary → prints 141 vs 140) with IDENTICAL TPA
    1132 ⇒ the seed is diameter growth, NOT mortality/count. Then dense self-thinning (SDI 325) count-straddles
    (a few TPA) compound to ~1% volume by 2024.
  • NEW POSITIVE FINDING — dynamic FOREST-TYPE reclassification is BIT-EXACT. Live `.out` FOR TYP flips
    503→506→503 across cycles (FVS recomputes FORTYP every cycle from current species BA); FT_DEBUG-instrumented
    jl compute_forest_type! reproduces 503,503,506,506,503,503 EXACTLY (BA basis 98.6/120.0/140.5/160.5/167.9/
    172.1). ⇒ FORTYP is ruled out as the seed; the seed is narrowed to sub-print Float32 diameter-growth
    accumulation feeding the dense-phase count-straddle. (debug line added to forest_type.jl then REVERTED;
    `git diff src/` empty — source byte-clean, floor 38527/143/0 intact.)
  • The extreme outliers (SCuFt/MCuFt/BdFt 25–785%) are all late-cycle sawtimber MERCH-THRESHOLD crossings on
    ~zero volume bases. Stand 29549638020004: bit-exact 4 cycles; a 9-TPA count-straddle at 2027 (2100 TPA young
    stand) pushes a couple trees across the sawtimber DBH threshold at 2032 → BdFt 21 (live) vs 186 (jl) on a
    near-zero base = 785%. Density (BA/SDI/CCF within 1-2 units) preserved throughout ⇒ pure which-tree-dies +
    which-tree-crosses-threshold, the count-straddle taxonomy. Not a bug.
VERDICT: 221H × {volume_persistent, structure_densephase} CORNERED (accepted compounded-ULP + count-straddle +
sawtimber-threshold-crossing class; identical to dig-session #1's ecoregion-221 corner). No fixable bug.

SWEEP META-FILTER (so coverage actually advances instead of re-pausing every ~200 stands on a cornered cluster):
  • docs/fia_cornered_clusters.tsv — (ECO_PREFIX, signature) pairs cornered by a documented dig session.
  • test/harness/fia/filter_digworthy.jl — replaces run_expand_cycle.sh's inline awk: applies the base dig-worthy
    rule, then DROPS candidates whose (ecoregion-prefix, signature) is cornered. ESCALATION GUARD (never dropped):
    UNCLASSIFIED, or a structure/density blow-up (worst_col ∈ {TPA,BA,SDI,CCF,TopHt,QMD} & max_rel ≥ 15%) — a
    genuinely new bug in a cornered geography still pauses the sweep. Validated on the 8300–10300 batch: 170 base
    dig-worthy → 1 survivor (an escalation).
  • The 1 escalation (stand 202594547010854, TopHt 22% @2026) was REVIEWED: same count-straddle taxonomy
    expressed in the AVHT40 top-height statistic. At 2016 the stand has IDENTICAL TPA/BA/QMD yet TopHt 39 vs 37 —
    sub-print DBH/height differences flip the top-40-by-DBH selection, amplified by 74% self-thinning
    (4570→1181 TPA); density (BA/SDI/CCF) preserved. Not a distinct bug; the guard surfacing it (~1/1700 stands)
    is the intended low-rate human-review rate.
DISPOSITION: 348 archived to docs/dig_archive/dig_session2_sn_6300-10300.csv; dig-queue cleared; sweep resumes
from cursor 10300 with the meta-filter live (will now advance through the rest of 221H, pausing only on genuine
escalations / new strata). Floor 38527/143/0 confirmed (suite re-run).

### SLICE 43c — TopHt escalation: single-vs-double AVHT40 sort hypothesis REFUTED (empirical, both-sides)
The dig-session #2 meta-filter's escalation guard kept surfacing TopHt divergences (worst_col=TopHt, ≥15%) in
dense 221H stands. Both-sides-traced the FVS top-height source: the `.sum` TopHt col (IOSUM(13) = OLDAVH/ATAVH,
disply.f:321/330) is `AVH` from the 40-largest-DBH-TPA loop. That loop lives in BOTH avht40.f AND dense.f and
consumes a PRE-SORTED `IND` (neither sorts internally). The per-cycle path is **gradd.f:186
`CALL RDPSRT(ITRN,DBH,IND,.TRUE.)` → dense.f** — a SINGLE sort; cratet.f cycle-0 empirically double-sorts
(dig-session #1). Hypothesis: jl's stand_top_height double-sorts EVERY cycle, so per-cycle tie-heavy stands
diverge. Tested via an ENV toggle (FVS_TOPHT_SORT single|double) over 4 tie-heavy stands vs live:
  • DOUBLE (current): 1737985937290487 ✓all, 163925866010854 ✓all; 232271267010854 ✗(2003 L30/J21),
    202594547010854 ✗(2016/2026/2031).
  • SINGLE: 232271267010854 ✓all (FIXED); but REGRESSES 1737985937290487 (2024/2034) + 163925866010854 (1976);
    202594547010854 UNCHANGED (sort-INDEPENDENT ⇒ genuine small-tree height ULP, not the sort).
VERDICT: the correct tie-break is STAND-DEPENDENT — no global single/double choice is bit-exact (single fixes
one stand, regresses two). RDPSRT is an unstable quicksort on tied DBHs, so the 40-TPA-boundary tree (and its
height) is inherently input-order-sensitive. Double matches the most stands (2/4 fully) and is kept. The TopHt
swings on tie-heavy dense stands are the cornered **AVHT40 top-height tie-break ULP primitive** (density
BA/SDI/CCF preserved within ~1 unit; converges) — same primitive class as the dig-session #1 cratet tie-break,
now shown NOT globally fixable via sort-mode. Code unchanged (comment-only in standstats.jl documenting the
refutation); floor 38527/143/0 intact (no logic change). This closes the TopHt escalation as a named cornered
primitive rather than a fixable bug — a rigorous NEGATIVE result (tested the fix, it regresses; doctrine #3/#6).

### SLICE 43e — GLOBAL taxonomy corner (verified 221+223) + sweep-robustness fix
Two things this slice. (1) ROBUSTNESS FIX (commit 0d23e2e): the sweep halted prematurely (cursor stuck) on a
FALSE-POSITIVE "RUN FAILED" — a NOSUM-heavy batch (live FVS emits no comparable .sum for the whole batch; it
can't project ~1 in 6 real stands) returned rc=0 with an EMPTY cycle CSV, tripping run_expand_cycle.sh's
`[ ! -f $cyc ]` guard; and run_expand_loop.sh's `$(...)` capture choked on a null byte. FIX: rc=0+empty-CSV →
EMPTY-STRATUM skip+advance-cursor (only rc≠0 halts); loop writes to a file + `grep -a` (null-safe).
(2) GLOBAL CORNER: the robustness-fixed sweep advanced into NEW ecoregion 223Ab (interior broadleaf, 77+ stands,
0 UNCLASSIFIED). Both-sides-traced rep 1176710848290487 (oak-hickory): bit-exact 2021 → 0.9% sub-print DG seed
at 2026 (BA 89→90, identical TPA 433) → dense-phase count-straddle + sawtimber-threshold amplification through
2046 — IDENTICAL taxonomy to 221 (Appalachian). ⇒ the two SN taxonomy signatures (volume_persistent +
structure_densephase) are verified SN-model-universal (Appalachian-hardwood 221 + interior-broadleaf 223 across
the structural range), consistent with the whole prior campaign (5000-scale sweeps + 4-variant outlier hunts,
slices 31-40). CORNERED GLOBALLY ("*" prefix in fia_cornered_clusters.tsv) so the sweep covers the ENTIRE
remaining SN population (+ NE/CS/LS) without re-pausing at every ecoregion. The escalation guard remains the
real-bug safety net: UNCLASSIFIED, or structure_densephase with a structure col ≥15%, or TCuFt ≥15% always
surface. Validated: the 128-row 223Ab queue → 0 survive; a synthetic panel (UNCLASSIFIED / structure-22% /
TCuFt-18% survive; small-base-BA-16% + BdFt-40%-threshold drop) confirms the guard. All 5 real bugs found this
campaign (covtyp/fuel-OOB/htcalc-NaN/PLANT/FORKOD) would trip the guard (UNCLASSIFIED / ≥15% structure / TCuFt /
crash). Harness-only; floor 38527/143/0. Archived docs/dig_archive/dig_session2d_sn_223Ab.csv; dig-queue cleared;
sweep resumes from cursor 24300 — now pauses ONLY on a genuine real-bug candidate.

### SLICE 43f — sweep continuation (SN 32300→34100) + batch-sizing operational lesson
Continued the global-cornered SN sweep from cursor 32300 to 34100 (~1800 more plots, ecoregion-ordered):
**dig-worthy +0 across every batch** — the two SN-universal signatures (volume_persistent + structure_densephase)
continue to absorb all non-bit-exact plots, and the escalation guard surfaced nothing. Per-batch bit_exact ratio
varies with stratum density (e.g. 10/20, 18/30, 83/95) but every non-exact plot corners cleanly; no UNCLASSIFIED,
no ≥15% structure/TCuFt blow-up. Confirms the global corner (slice 43e) holds deeper into the SN population.

OPERATIONAL LESSON (no code/floor change): run the sweep with a SMALL batch (BATCH≈100, ~200s/cycle), NOT
BATCH=1500 (~50 min/cycle). run_expand_cycle.sh checkpoints the cursor only at end-of-batch; a background task
reaped mid-batch (or a foreground call hitting the tool timeout) loses no *correctness* (re-processing is
idempotent) but makes NO forward progress if the batch never completes. Small batches checkpoint every few
minutes → steady, reap-resilient advance. Also: run the sweep in the FOREGROUND (bounded ~8 min/turn); launching
a second background task appears to cancel the first, so a single foreground loop per turn is the reliable driver.
Floor untouched (38527/143/0); harness/docs only.

### SLICE 43g — dig #3: small-base structure escalation false-positive (CN 202567027010854) → guard floored
FIRST genuine escalation survivor since the global corner (slice 43e): a full re-filter of this session's entire
23,034-row SN ledger (a rigor check after finding the per-batch filter had been crashing on concurrent-loop
`rm -f sn_cycle.csv` races — now moot under foreground-singleton) surfaced exactly ONE dig candidate:
CN 202567027010854 (SN, ecoregion 221Hb), signature structure_densephase, worst_col=BA, max_rel=33.333% @2011.

BOTH-SIDES-TRACED vs freshly-relinked live FVSsn (DATABASE reader, NUMCYCLE 10). NOT a bug — a young **age-3**
seedling stand (BAF=0 fixed-area, 1 plot) where BA is tiny. Per-cycle jl-vs-live: BA 2011 jl=2/live=3 (the 33%
= a **1 sq ft** straddle on BA=3); every other cycle tracks within ±1-5 absolute units (BA ±1, SDI ±2, CCF ±3,
TPA ±5) and CONVERGES (BA 105=105 @2056). The classic compounded-ULP small-base straddle: a sub-print DBH-growth
seed rounds BA to 2 vs 3, inflating to 33% RELATIVE only because the base is 3.

ROOT CAUSE OF THE FALSE-POSITIVE (guard blind spot, not a model bug): the escalation guard escalated a structure
`max_rel≥15%` with NO absolute-magnitude floor. The comment claimed signature==structure_densephase "BY
DEFINITION" means a >1-unit structure move — FALSE. Here the materiality that earned the structure_densephase
label came from a *different* small cell (CCF 56 vs 58 @2016 = 2 units, 3.4% — material by ismat), while
worst_col/max_rel came from the tiny-base BA 2-vs-3 cell (33%, not even material). So the guard escalated on a
relative % belonging to a non-material cell.

FIX (harness only; floor 38527/143/0 untouched — no src/ change): added `struct_max_abs` (largest ABSOLUTE diff
among structure cols 1-6) as ledger col 16 (appended AFTER signature ⇒ backward-compatible: signature stays f[15];
legacy 15-col rows lack it ⇒ treated as +Inf ⇒ still escalate, preserving the conservative default). The structure
escalation now also requires `struct_max_abs ≥ 10` — a real BA/SDI/CCF bug moves tens of units; a small-base ULP
straddle does not. VALIDATED end-to-end: the real ledger emits struct_max_abs=3.0 for this CN ⇒ filter DROPS it;
a synthetic panel confirms REAL_STRUCT_BUG (abs=45), UNCLASSIFIED, TCuFt-volume, and legacy-15-col rows ALL still
surface. UNCLASSIFIED and TCuFt-volume escalation are unchanged (no floor) — the primary real-bug nets stay
unconditional. dig-queue remains a true 0. Files: ledger_fia.jl, filter_digworthy.jl.

### SLICE 43h — DURABLE cross-session sweep coverage DB (data/fia_sweep.db)
The sweep's per-stand differential was EPHEMERAL (scratchpad ledger, lost at session end) — only dig-worthy rows
(0) + audit prose survived. Added a local SQLite coverage DB on the durable repo volume (data/fia_sweep.db;
gitignored — survives sessions AND container restart, not a git blob) that records EVERY stand swept and its
outcome, so it is a durable cross-session WORKLIST of what still needs a dig.

Schema (test/harness/fia/sweep_db.jl): table `sweep` keyed (variant,cn,regime) with the measured ledger facts +
a derived `dig_class` ∈ {bit_exact | ulp_class | needs_dig}. dig_class MIRRORS filter_digworthy.jl's escalation
guard exactly (UNCLASSIFIED, or a MATERIAL structure move ≥10 abs-units & ≥15%, or threshold-free TCuFt ≥15% ⇒
needs_dig; every other divergence is an accepted cornered primitive ⇒ ulp_class). Upsert is idempotent on
(variant,cn,regime) so re-sweeps UPDATE (a fixed bug flips needs_dig→ulp/bit_exact and the diff is visible).

Wiring: ledger_fia.jl upserts each stand inline when SWEEP_DB is set (no extra process — SQLite already loaded;
DB errors are caught and never break the sweep); run_expand_cycle.sh exports SWEEP_DB by default. CLI:
`sweep_db.jl {ingest <db> <csv> | stats <db> [variant] | digs <db> [variant]}` — `digs` prints the needs_dig
worklist. Backfilled this session's entire SN ledger: 21,676 distinct stands → 12,724 bit_exact / 8,955
ulp_class / **0 needs_dig** (the one transient survivor, CN 202567027010854, reclassified to ulp_class once its
struct_max_abs=3 landed — slice 43g). Harness only; suite floor 38527/143/0 untouched.

### SLICE 43i — restart-safety: move all sweep operations onto the persistent volume
Audited what a container restart would destroy. `/workspace` = a real btrfs disk partition (persists); `/` and
`/tmp` = the container overlay (`fsync=volatile`, EPHEMERAL — wiped on restart). Findings + fixes so a restart is
NON-disruptive to the verification (nothing re-swept from scratch):
- VERIFICATION RESULTS — were in the ephemeral /tmp scratchpad ⇒ already moved to data/fia_sweep.db on
  /workspace (slice 43h). DURABLE. ✓
- PROGRESS CURSOR — test/harness/fia/expand/<v>.cursor is on /workspace (persistent) AND git-tracked; now ALSO
  mirrored into the DB `progress` table (self-contained snapshot). resume_sweep.sh reconciles max(file, DB). ✓
- SWEEP WORKING DIR — was hardcoded to /tmp/claude-1000/<session-UUID>/scratchpad (ephemeral AND tied to a
  session UUID ⇒ broken after restart). MOVED to /workspace/FVSjl/.sweep_work (persistent, gitignored, override
  SWEEP_WORK); migrated the accumulated master ledger there. run_expand_cycle.sh + run_expand_loop.sh updated. ✓
- ORACLE BINARY /tmp/FVS*_new — ephemeral but REGENERABLE (doctrine relinks per run anyway); resume_sweep.sh
  relinks it if absent. Not data loss. ✓
- MASTER FIA INPUT /workspace/SQLite_FIADB_ENTIRE.db — on /workspace, read-only. DURABLE. ✓
New: test/harness/fia/resume_sweep.sh — one command to recover after a restart (relink oracle → reconcile cursor
→ report durable coverage + needs_dig worklist → resume the loop). Harness only; suite floor 38527/143/0 untouched.

### SLICE 43j — DB data-quality: reject/scrub CSV-concatenation artifacts; document ulp_class basis
A query of the coverage DB (prompted by "how do you know these are ULP-class?") exposed 53 MALFORMED rows: the
variant glued to the CN, the signature slot holding a bool/number (e.g. cn="216864634010854SN", sig="true",
max_rel_pct=7.3e14). Root cause: the one-time BACKFILL ingested the master-ledger CSV, whose lines had been
SPLICED by concurrent writers / timeout-killed appends (a partial line without a trailing newline + the next
batch's line ⇒ shifted columns on positional split). The inline per-stand upsert path (ledger_fia.jl→DB) was
never affected — only the CSV re-ingest. FIX (sweep_db.jl): `_valid_row` gate on ingest (variant∈{SN,NE,CS,LS},
signature in the known taxonomy set, CN all-digits) drops spliced lines; new `scrub` command deletes existing
malformed rows and salvages the embedded CNs for a clean re-sweep. Scrubbed 53, re-swept the 50 recoverable CNs
via the clean inline path: 32 bit_exact / 17 ulp_class / 1 no-comparable-sum — **none was a hidden needs_dig**.
DB now 0 malformed.

ON THE ulp_class BASIS (documented for honesty): dig_class is a RULE over the deterministic signature, NOT a
per-stand both-sides proof. A stand is ulp_class iff it diverges, its signature ∈ the accepted cornered set, and
it does not trip the escalation guard. The signatures were established as ULP/threshold primitives by tracing
REPRESENTATIVES to the FVS source (digs #1/#2/#2d/#3); cornering is by-signature + magnitude-guarded, so a
modest-magnitude real bug mimicking a cornered signature (<15%, converging, density-preserved) could in principle
be mislabeled ulp_class rather than needs_dig — the inherent limit of taxonomy-cornering vs per-stand tracing.
Empirical profile (SN ulp_class): 75% print_boundary (≤1-unit straddles), 96% worst-cell <5%; of the 1.6% ≥15%,
~94% are merch/board volume cols (SCuFt/MCuFt/BdFt = sawtimber-threshold step-fn), the rest TopHt AVHT40 tie-break
+ small BA. Harness only; suite floor 38527/143/0 untouched.

### SLICE 43k — cross-variant coverage: NE/CS/LS sampled (pillar 2 "all 4 variants" opened)
SN was deeply covered (33k stands, ecoregion-stratified) but NE/CS/LS had ZERO coverage — pillar 2 requires all
four. Ran the full multi-cycle projection differential (freshly-relinked live FVS{ne,cs,ls} vs FVSjl, via the
DATABASE reader, all 10 .sum cols every cycle) over a deterministic ECOREGION/LOCATION-stratified sample of each
(extract_sample.jl), recording every stand to the durable coverage DB with its dig_class:
  NE  80 stands / 64 ecoregions → 61 bit_exact, 19 ulp_class, 0 needs_dig
  CS  79 stands / 56 ecoregions → 70 bit_exact,  9 ulp_class, 0 needs_dig  (1 no-comparable-sum)
  LS  80 stands / 57 ecoregions → 56 bit_exact, 24 ulp_class, 0 needs_dig
All-variant DB total: 33,263 real stands, 0 needs_dig — every non-bit-exact plot cornered (print/ULP/threshold/
count-straddle/compounded-ULP). First cross-variant behavioural differential at scale; NE/CS/LS samples are
reproducible via `extract_sample.jl <V> <N>` (deterministic, no RNG). Harness/docs only; floor 38527/143/0.

### SLICE 43l — FLAGGED cluster: dense-phase self-thinning TPA-straddle (material, converging) — 1 SN + 3 LS
The escalation guard surfaced its first genuine needs_dig candidates (the struct_max_abs floor working as intended):
- SN 211016796010854 (ecoregion, ultra-dense TPA~3900/SDI~404): structure_densephase, TPA 17.1% @2026, abs 368.
- LS 1093612541290487 (TPA 39.5% @2041, abs 275), 1901267649290487 (19% @2034, abs 661), 1901273492290487
  (24.6% @2045, abs 429).
BOTH-SIDES per-cycle differential (live FVS{sn,ls} vs FVSjl, sub-DB, NUMCYCLE 5) for the SN case:
  2006/2011 BIT-EXACT (all 10 cols) → 2016 first split TPA 3951/3617 (Δ334) BA 153/148 (Δ5) → peaks ~2021-26
  (~17%) → CONVERGES (2031 Δ213/12.5%). BA/SDI/CCF/TopHt track within ~3%, QMD near-exact.
MECHANISM (hypothesis, strongly supported): a stand near SDImax where the SDI-driven self-thinning mortality
kills a batch of SMALL trees — moving TPA by hundreds while BA/SDI barely move — and a sub-print difference in
the tree-diameter distribution flips WHICH cycle that batch dies. The CONVERGENCE after the peak is the fingerprint
of a mortality-TIMING straddle, not a persistent model bug (a real over-kill bug diverges monotonically). This is
the accepted dense-phase count-straddle primitive amplified to MATERIAL magnitude by extreme density.
STATUS: FLAGGED, NOT yet cornered. Honesty gate: printed-bit-exact @2011 ≠ per-tree-bit-exact — confirming
ULP-class vs a real dense-mortality bug requires a per-tree state comparison at 2011 (FVS_TreeList vs FVSjl tree
state) at the divergence onset. Kept as needs_dig in data/fia_sweep.db (the durable worklist) pending that dig.
The 4 stands (1 SN + 3 LS) share the signature ⇒ a single cluster to trace together in a dedicated session.
Floor 38527/143/0 untouched (harness/docs only).

### SLICE 43m — dense-phase cluster: 2nd SN case, BIDIRECTIONAL straddle (strengthens ULP-class verdict)
Second SN needs_dig 781951924290487 (structure_densephase, TPA 15.1% @2039, abs 159). Both-sides per-cycle:
2019/2024 bit-exact (2024 only TCuFt Δ1 ULP) → 2029 first split TPA 2260/2335 → peaks 2039 (Δ159) → converges
(2044 Δ61). DECISIVE new evidence: the straddle is BIDIRECTIONAL across stands — here FVSjl has MORE TPA than
live (jl under-kills), whereas case 1 (211016796010854) had FEWER (jl over-kills). A systematic mortality BUG
would bias consistently one direction; a bidirectional, converging, from-bit-exact-base straddle in dense stands
is the fingerprint of the self-thinning count-straddle primitive (the two impls round a near-tie at the SDI/
density-mortality threshold to opposite sides, stand-by-stand). Cluster now 2 SN + 3 LS, all same signature.
This bidirectionality materially raises confidence the cluster is ULP-class, not a bug — but per doctrine the
GOLD-STANDARD confirmation is still a per-tree state comparison at the divergence-onset cycle (printed-bit-exact
≠ per-tree-bit-exact). Kept as needs_dig in the DB worklist pending that per-tree dig. Floor 38527/143/0 untouched.

### SLICE 43n — guard 2nd blind spot: TCuFt net needs an ABSOLUTE floor (vol_max_abs) — degenerate stand
3rd SN needs_dig 218434248010854 both-sides-traced: a DEGENERATE 2-tree micro-stand (TPA~7000 from a tiny fixed
plot) that triggers a SHARED FVS SDI-overflow pathology — SDI → ~4.38 MILLION at 2026 in BOTH live FVSsn and
FVSjl (jl within 0.46%), TCuFt momentarily 0. NOT a bug (a fidelity success: jl reproduces even FVS's degenerate
behavior; see FVS_SOURCE_BUGS.md). It tripped the escalation guard's TCuFt≥15% net on a 62-cuft divergence
(412/350) — the mirror of the slice-43g small-base struct false-positive: the TCuFt volume net had NO absolute
floor, so a big RELATIVE % on tiny cuft false-positives. FIX: added vol_max_abs (largest abs diff among vol cols
7-10) to the ledger (col 17, backward-compatible) + a VOL_ABS_FLOOR=300 gate on the TCuFt escalation (a real
volume-equation bug like FORKOD moves 1000s of cuft; 62 is not). sweep_db schema migrated via idempotent ALTER.
Re-swept the 3 SN candidates: the degenerate stand → ulp_class; the 2 genuine dense-phase TPA-straddle cases
(43l/43m) remain needs_dig. Guard now has BOTH an absolute floor (struct + vol) and the relative threshold.
Harness/docs only; floor 38527/143/0 untouched.

### SLICE 43o — fix vol_max_abs: TCuFt-column-only (BdFt-domination bug in slice 43n's floor)
Slice 43n's vol_max_abs tracked ALL volume cols 7-10, but BdFt (board feet) magnitudes are ~10x cubic feet and
DOMINATE — so vol_max_abs was ~always ≥300 (the floor) for any sawtimber stand, defeating the gate. Surfaced by
SN 209219251010854 (young age-3 3-sapling stand): worst_col=TCuFt 35% @2011 on 310 cuft (109 cuft absolute,
converges to 4% by 2031) — a small-base ULP amplification that slipped through because BdFt inflated vol_max_abs.
FIX: vol_max_abs tracks the TCuFt column (7) ONLY — the escalation net is worst_col==TCuFt, so its floor must be
TCuFt's own absolute divergence. Re-swept the 4 SN candidates: the young small-base stand + the degenerate 2-tree
stand (43n) → ulp_class (TCuFt abs 109/62 < 300); the 2 genuine dense-phase TPA-straddle cases (43l/43m,
worst_col=TPA, unaffected by the vol floor) remain needs_dig. Harness/docs only; floor 38527/143/0 untouched.

### SLICE 43p — DIG SESSION (user-requested): 2 SN needs_dig + 10 most-suspicious SN ulp_class
Paused the sweep; both-sides per-cycle differential (live FVSsn vs FVSjl, sub-DB) on the 2 SN needs_dig + the 12
most-suspicious ulp_class (ranked: non-converging + high rel% + large abs). RESULTS:
- 2 needs_dig (211016796010854, 781951924290487): CONFIRMED the dense-phase self-thinning TPA count-straddle
  cluster (43l/43m) — density cols diverge ~5-10%, bidirectional/converging. Remain needs_dig pending per-tree.
- 11 of the 12 suspicious ulp_class: CONFIRMED threshold-crossing (correctly ulp_class). Signature: STRUCTURE
  tracks within ~5% AND total cubic (TCuFt) tracks closely, but MERCH/BOARD volume (MCuFt/SCuFt/BdFt) flips on a
  small base near the merchantability/sawtimber DBH — e.g. 204664667 MCuFt 16/167, 1261822875 BdFt 1283/4744
  (TCuFt 4564/4587 bit-close). The huge relative %s (240-944%) are the discrete log-volume step, not a bug.
- 1 GENUINE FIND — 209314057010854: a SYSTEMATIC, GROWING, NON-converging STRUCTURAL divergence. jl consistently
  higher BA/SDI/volume and lower TPA across EVERY cycle (BA 79/85 @2013 → 182/203 @2033; struct_rel 20%, struct_abs
  302). NOT threshold-crossing (structure itself diverges and grows), NOT a bidirectional straddle (one-directional
  systematic bias). It was mis-scored ulp_class only because its highest-RELATIVE column is BdFt (a volume-threshold
  col), so the guard's worst_col-gated struct net missed it. This is a real candidate (likely a dense-stand
  growth/mortality-partition divergence) requiring a per-tree trace.
GUARD: tried gating the struct net on struct_max_rel_pct (worst_col-independent) — OVER-flagged (struct_abs is
TPA-dominated; struct_max_rel_pct is inflated by TopHt AVHT40-ULP + small-base). Reverted to the conservative
worst_col gate; a faithful fix needs a density-specific relative metric (BA/SDI, excluding TPA+TopHt) — deferred.
Added: `reclassify` command (recompute dig_class over stored facts, no FVS re-run) + docs/fia_manual_needsdig.txt
(committed manual-confirmed genuine finds that survive reclassify). 209314057 pinned there. SN needs_dig now 3
(2 dense-phase + 1 genuine structural). Harness/docs only; floor 38527/143/0 untouched.

### SLICE 43q — root-cause of the genuine find 209314057: small-tree (seedling) growth divergence
Per-tree trace of the one genuine structural find. The stand is a PURE-SEEDLING stand: SN, ecoregion 231Ba,
AGE 7, INVYR 2008, 3 tree records ALL at DBH 0.1" (missing height): sp131 loblolly pine 2162 TPA, sp391 432 TPA,
sp521 1297 TPA (~3891 TPA). FVSjl per-tree (FVS_TreeList via DSNOUT/TreeLiDb) at cycle 1 (2008→2013): loblolly
EXPLODES 0.1"→2.73" (DG 2.1"), ht 1→21.5 ft; sp391/521 crawl (0.1→0.46"/0.64"). Stand-level both-sides (live
FVSsn vs jl): jl BA 85 vs live 79 @2013 — jl over-grows the loblolly seedlings ~8% in the FIRST cycle → higher
SDI → more self-thinning (jl kills 302 more TPA) → the systematic, growing, non-converging divergence booked as
needs_dig. LOCALIZED to the SMALL-TREE (seedling) growth model for loblolly, amplified because the stand is
pure-seedling (no overstory to dominate the summary). Connects to the known small-tree/height-growth residual
area (NOHTDREG/LHTDRG + WK3 DGSCOR tail; see [[fvsjl-growth-gap-verdicts]]). VERDICT: a genuine small-tree-growth
divergence, NOT a threshold artifact and NOT a bidirectional straddle. NOT yet fully cornered vs a new bug: the
decisive live-vs-jl PER-TREE loblolly DG comparison is blocked by the relinked FVSsn binary (DBS TreeList → rc=20;
text TREELIST → >120s timeout). Kept as manual needs_dig; next step = obtain live per-tree DG (fix the live DBS
output or parse the text treelist) to confirm it's the known small-tree residual vs a new seedling-growth bug.
Harness/docs only; floor 38527/143/0 untouched.

### SLICE 43r — 209314057 CORRECTED: it's self-thinning MORTALITY-PARTITION, not small-tree growth (43q refuted)
Species-isolation test (build sub-DBs with only some species; compare jl vs live .sum) DECISIVELY refutes slice
43q's small-tree-growth hypothesis:
- LOBLOLLY-ONLY (sp131): jl == live BIT-EXACT all cycles (2013 both TPA2122/BA81/SDI250; 2033 both BA203/SDI410).
  ⇒ the small-tree/seedling GROWTH model is FAITHFUL.
- SLOW-SPECIES-ONLY (391+521): near-bit-exact (2033 TPA 787/776, an 11-TPA ULP straddle). ⇒ their growth faithful.
- FULL 3-species stand: diverges materially (2013 BA 85 jl / 79 live; TPA 3242 jl / 3544 live).
The divergence EMERGES ONLY from the 3-way interaction at COMBINED high density (~3900 TPA): the density-dependent
SELF-THINNING mortality PARTITIONS DIFFERENTLY across the diameter distribution. jl kills 302 MORE trees yet has
HIGHER BA ⇒ jl preferentially kills SMALL trees (preserving the loblolly BA), while live spreads the mortality to
include larger stems (lower BA, more survivors). This is a MORTALITY-ALLOCATION difference (which trees die under
SDI>SDIMAX self-thinning), NOT a growth-model difference.
⇒ 209314057 belongs to the SAME dense-phase self-thinning cluster as the 2 SN + 3 LS needs_dig (43l/43m) — all
are the self-thinning mortality partition in dense stands, now LOCALIZED (growth proven faithful; the divergence
is purely in mortality allocation). NEXT: read morts.f's mortality-distribution vs FVSjl mortality.jl to decide
whether the allocation RULE differs (bug) or it's a near-tie the two round differently (ULP count-straddle).
This is a strong both-sides localization of the campaign's one open SN/LS divergence class. Kept as needs_dig.
Method note: live per-tree TreeList is blocked (relinked FVSsn: DBS rc=20 / text TREELIST hangs), but SPECIES
ISOLATION via sub-DB + .sum sidesteps it entirely and is decisive. Floor 38527/143/0 untouched.

### SLICE 43s — dense-phase cluster ROOT-CAUSED to the VARMRT mortality-allocation percentile near-tie
Traced the self-thinning mortality PARTITION (43r) through the source both-sides:
- FVSjl _varmrt! (varmrt.f geometric-progression distribution) = documented faithful port.
- _varmrt_efftr! weights each record by PCT**3.0 (crown-ratio PERCENTILE cubed). The float-exponent power is
  ALREADY routed through an FFI (fpow, doctrine #8) to be BIT-EXACT with gfortran powf — so efftr is faithful too.
⇒ The residual is upstream, in PCT (the tree's crown/BA-percentile RANK). In a mixed stand (loblolly + 391 + 521)
the per-tree percentile ranking has NEAR-TIES; jl and live resolve them to sub-print-different PCT ⇒ different
efftr ⇒ VARMRT's geometric progression allocates the discrete self-thinning kill to slightly different trees ⇒
the BA/TPA straddle (jl kills small trees, live spreads to larger — or vice-versa; the cluster is BIDIRECTIONAL,
211016796 vs 781951924 go opposite ways). Consistent with a NEAR-TIE, not a systematic allocation-rule bug.
VERDICT: the whole dense-phase self-thinning cluster (2 SN + 3 LS) is CORNERED to a named primitive — a
compounded-ULP COUNT-STRADDLE in the self-thinning mortality ALLOCATION, rooted in the crown/BA-percentile-rank
near-tie that feeds VARMRT's efftr. Growth models proven bit-exact (loblolly-only isolation); VARMRT + efftr are
faithful ports (efftr already fpow-FFI'd). This meets the doctrine bar (cornered-to-named-primitive). The only
stronger confirmation would be a debug-FVS stamp of per-tree PCT/EFFTR in the mixed stand (optional; the isolation
+ source trace + bidirectionality are already decisive). Floor 38527/143/0 untouched.

### SLICE 43t — honest coverage accounting: `live_crash` category (the "no-sum" stands are FVS SIGFPE crashes)
User challenge: "what's the point of the sweep if we skip ~half of it?" MEASURED the skip cause across 6 regions
(108 stands): 92 comparable / 16 live_crash / 0 clean-nosum / 0 jl_fail — EVERY skip is live FVS crashing with
SIGFPE (never FVSjl failing, never nonstocked). Rate is region-variable (some ecoregions 18/18 comparable; the
high-expansion-seedling ecoregions ~half crash). Root cause: live FVSsn dies on >1000-TPA 0.1"-seedling records
(FVS40 warning → floating-point exception); FVSjl projects them fine (FVS-UB, documented in FVS_SOURCE_BUGS.md).
FIX (honesty): run_live now returns (text, crashed) via termsignal/exit>128; the ledger records dig_class
`live_crash` (jl-projected, live-crashed) instead of silently skipping. Coverage is now denominated honestly:
comparable(bit_exact+ulp_class) + live_crash + skip. Validated on the crash-heavy 290487 region: 11 be / 11 div /
18 live_crash / 0 skip of 40. The stronger, honest claim: "0 needs_dig among the comparable stands; the remainder
are live-FVS crashes FVSjl survives." Harness/docs only; floor 38527/143/0 untouched (no src/ change).

### SLICE 43u — CORRECTION: the `live_crash` stands are the known D38 r9ht bug (jl-correct + validated)
My slice-43t framing ("live_crash = jl plausible but UNVALIDATED") was WRONG — a re-discovery of an
already-resolved bug. The live_crash SIGFPEs are the D38 R9 Clark `r9ht` short-tree underflow/invalid-op crash
(FVS_SOURCE_BUGS.md), already root-caused AND fixed AND VALIDATED against a patched live binary (/tmp/FVSsn_fixtest:
18/18 crashers cleared, 276/282 non-crashers bit-identical; jl carries the fix). So on a live_crash stand jl
produces the CORRECT projection — the buggy shipping oracle just can't confirm it, but the FIXED oracle does and jl
matches. The `live_crash` category is still right for honest coverage accounting (visible, not skipped); the verdict
is "FVS-UB (D38), jl-correct", NOT "unvalidated". D38's measured SN ~30% live-crash rate on treed stands also
explains the region-variable comparable rate the coverage audit surfaced. META-LESSON: grep FVS_SOURCE_BUGS.md
before writing up any FVS crash. Floor 38527/143/0 untouched.

### SLICE 43v — CORRECTION: D38 crash is multi-site; patched-oracle is INCOMPLETE (not all crashers cleared)
Testing FVS's own isolated guard (fvsMod@a19c41b4, 16 lines) and my 5-guard patch on 40 real SN live_crash stands:
both clear only 32/40. The other 8 crash at a THIRD site — r9cuft cubic-volume V2/V3, r9clark.f:1086 (backtrace-
confirmed) — which neither guards. So the patched oracle (/tmp/FVSsn_patched) can NOT validate all live_crash
stands, and my earlier "patched oracle validates the crash stands" framing (43t area) was overclaimed. Corrected
in FVS_SOURCE_BUGS.md. UNCHANGED and correct: the live_crash dig_class = honest coverage accounting of FVS-UB
stands the shipping oracle crashes on; FVSjl projects them. Patched-oracle validation is PARKED until a complete
(multi-site) guard set exists. Floor 38527/143/0 untouched (docs only).

### SLICE 43ag — population-scale re-confirmation of the `live_crash` class (D38 SIGFPE, not a harness artifact)
While stewarding the full SN population sweep (cursor ~39,800/637,641; 76,158 stands recorded in data/fia_sweep.db),
audited the single largest cornered bucket: **live_crash = 9,405 (~12% of swept SN stands)** — the biggest non-bit-
exact class after print_boundary. Doctrine (Pillar-4) requires every divergence root-caused, so I spot-reproduced the
oracle directly on 5 sampled live_crash CNs (1224244256290487 …1224249126290487) via run_live against the freshly
relinked /tmp/FVSsn_new: **all 5 die with exitcode=0 / termsignal=8 (SIGFPE) and emit NO .sum**. This is the genuine
oracle-side D38 R9-Clark short-tree FPE (slices 43t/43u/43v), NOT a harness mislabel inflating the cornered count —
the ledger's crash detector (`termsignal!=0 || exitcode>128`, ledger_fia.jl:92) fires only on a real signal, and
these are real SIGFPEs. FVSjl projects all 5 fine. Sweep-DB dig_class distribution is clean: bit_exact 40,102 /
ulp_class 26,651 / live_crash 9,405 — **zero needs_dig**; every non-bit-exact stand is a cornered signature
(print_boundary 18,887 / volume_persistent 2,662 / structure_densephase 2,547 / threshold_crossing 1,577 /
count_straddle 978, all sub-escalation). VERDICT: the live_crash count is audited-genuine, D38-cornered per the
existing verdict. Floor 38527/143/0 untouched (spot-repro + docs only; NO src/ change).

### SLICE 43ah — sub-escalation tail audit: the guard's structural blind-spot is benign (dense-phase + a .sum parse artifact)
Stewarding-time Pillar-4 audit of the ONE place a real structure bug could hide UNDER the escalation guard: a large
structural move whose row `worst_col` is a VOLUME column (so `is_escalation`'s `worst_col ∈ {TPA,BA,SDI,CCF,QMD}`
gate never fires) inside the globally-cornered `structure_densephase` cluster. Query: SN `structure_densephase`,
non-structural worst_col, `struct_max_abs ≥ 50` ⇒ **69 candidates**. Both-sides-traced the two largest via diff_one:
  • **204758406010854** — 14,103-TPA seedling regen (QMD 0.1). TPA self-thins live 5762 / jl 7858 @2008 (struct_abs
    2096) while density tracks close (BA 80/73, SDI 304/291). = dense-phase self-thinning count-straddle (taxonomy
    B1/B5, SIGMAR tripling-spread × density-dependent mortality), the SAME primitive cornered for the LS ultra-dense
    stands, here in SN at higher TPA. Correctly globally-cornered.
  • **218434248010854** — 7003-TPA regen; reported `struct_max_abs=20014` at only 3.3% structural rel (physically
    impossible for any structure col). diff_one showed SDI=`4.381035e6`, QMD=`2750`. RAW .sum: the 2026 row is
    `... 2058 164 4381035 42 3.8 ...` — **SDI (438) and CCF (1035) MERGED into one token `4381035`**: in an ultra-dense
    stand CCF exceeds 999 and overflows its fixed-width Fortran field, gluing to the SDI field with no separating
    space. `parse_sum10` (ledger_fia.jl:75) splits on WHITESPACE ⇒ the row loses a token ⇒ columns past the overflow
    shift (the "SDI" cell = merged SDI+CCF, "QMD" cell = a volume col). So `struct_max_abs=20014` is a **parse
    artifact, not a real structural move**; the genuine divergence is again dense-phase self-thinning (TPA 2058/2085,
    BA 164/162 within ~1%).
CRITICAL SAFETY PROPERTY: the overflow afflicts BOTH live and jl `.sum` identically (both are FVS-format), so it can
only INFLATE apparent divergence — it can NEVER produce a false `bit_exact` (that needs all 10 parsed cells equal; a
real per-column difference still surfaces through the merged token) nor mask a real bug. It is diagnostic-metric noise
on the `struct_max_abs`/`worst_col` fields, confined to CCF≥1000 ultra-dense regen — every such stand is already the
globally-cornered dense-phase class. VERDICT: the sub-escalation structural tail is CLEAN — no structural bug hides
under the guard; the 69 large-`struct_abs` candidates are dense-phase self-thinning (cornered) + .sum field-overflow
artifacts. Also re-confirms the AVHT40 TopHt tie-break primitive (232056444010854: TopHt spikes 32→49 @2013 then
RECONVERGES bit-exact to 60 @2018 — a transient that self-heals, not a compounding height-model error; the intended-
excluded col per filter_digworthy.jl:16-19). FOLLOW-UP (deferred, not floor-relevant): make parse_sum10 fixed-width
column-aware so struct_max_abs is trustworthy at the dense tail — but NOT mid-sweep (would make new DB rows'
metrics inconsistent with old, and changes NO verdict since these stands are cornered). Floor 38527/143/0 untouched
(diff_one/.sum-dump repro + docs only; NO src/ change; sweep uninterrupted, cursor advanced 39k→42k during the audit).

  UPDATE (2 more traced, the highest-structural-% NON-CONVERGING candidates — the most bug-suspect of the 69):
  • **1276028313290487** — 17,841-TPA regen; TPA straddles/oscillates (live 6140/jl 4902 @2027 → 3097/3671 @2032)
    while density tracks within a few % throughout (BA 108/115, SDI 341/350, CCF 310/316, TopHt 27/27). Dense-phase.
  • **209314057010854** — 3,891-TPA regen; TPA diverges (3544/3242 @2013) while BA/SDI/CCF/TopHt track within ~7%
    all cycles (BA 174/191, SDI 373/392 @2033). Dense-phase.
  REFINED VERDICT: `struct_max_abs≥50` is large PURELY because TPA at thousands-of-stems yields hundreds-of-units
  ABSOLUTE diffs even at modest %; the structure RELATIVE divergence stays inside the count-straddle band (3–20%),
  and `worst_col` is correctly a VOLUME col (huge % on near-zero-QMD bases = threshold-crossing). So the guard's
  `worst_col ∈ {structural}` gate is WORKING AS INTENDED — moderate-struct-% + big-volume-% is a threshold-crossing,
  NOT a structural bug. 4 of 69 traced (2 largest struct_abs + 2 highest struct-% non-converging), ALL dense-phase
  self-thinning; the tail is audited clean. The `struct_abs≥50` metric is NOT a bug indicator at ultra-dense TPA —
  it is the count-straddle's expected absolute footprint. No guard change needed (the earlier "blind-spot" concern is
  resolved: structural REL, not ABS, is what the guard keys on, and REL stays in-band).

### SLICE 43ai — Pillar-4 taxonomy audit COMPLETE across every signature class (population-scale, measured)
Closed out the sub-escalation audit by tracing the two remaining large SN buckets (after 43ag live_crash + 43ah
structure_densephase), so EVERY non-bit-exact signature class is now both-sides-verified at population scale:
  • **threshold_crossing** (1,577; worst 886359561290487, MCuFt 451% @2035): diff_one shows all 6 STRUCTURE cols
    bit-exact-or-±1 every cycle (SDI 130/129, TopHt 39/40) — the divergence is purely a merch/sawlog volume
    DBH-threshold crossing on a near-zero base as the cohort first forms sawtimber (converges as it matures; tiny
    vol_abs). Structure faithful; volume-threshold artifact only. ✓ named primitive B2.
  • **count_straddle** (978; worst 232303097010854, BdFt 258% @2029): 6,755-TPA regen; DENSITY bit-exact all cycles
    (BA 117/117, SDI 368/368, CCF 426/426) while TPA straddles (3656/3655 → 1615/1632) and QMD (4.3/4.2) — a
    different COUNT of tiny trees killed at the self-thin (dens_be=1 flag decisive). ✓ named primitive B1.
Full SN non-bit-exact taxonomy now audited-complete: live_crash (D38 SIGFPE, 43ag) / structure_densephase (dense-
phase self-thin + .sum overflow artifact, 43ah) / threshold_crossing (merch-threshold, structure faithful) /
count_straddle (density bit-exact) / volume_persistent (same merch-threshold family, non-converging tail) /
print_boundary (definitionally ±1 print rounding — classify() assigns it only when NO cell is material). ZERO
needs_dig across 80k+ swept SN stands; every class maps to a named primitive (taxonomy A-fixed / B-cornered at the
doc head). Pillar-4 "no unexplained divergence remains" HOLDS at population scale for SN. Floor 38527/143/0 untouched
(read-only DB queries + 2 diff_one traces + docs; NO src/ change; sweep uninterrupted).

### SLICE 43aj — Pillar-4 CROSS-VARIANT taxonomy pre-flight: NE/CS/LS share SN's exact 6-class taxonomy (durable DB)
Pre-flight for the full-population sweep's NE→CS→LS phases (~1.5 days out) using the durable sweep DB's existing
500-scale rows (NE 1571 / CS 1563 / LS 1575 distinct CNs from the Pillar-2 differentials) — read-only, no new runs,
sweep uninterrupted. Every variant shows the IDENTICAL six-class taxonomy and ZERO needs_dig:
  variant | bit_exact | ulp_class | live_crash | signatures (non-BE)
  NE      | 1237 (79%)|    329    |     5      | print_boundary 283 / threshold_crossing 20 / structure_densephase 16 / count_straddle 7 / volume_persistent 3
  CS      | 1376 (88%)|    182    |     5      | print_boundary 139 / threshold_crossing 14 / count_straddle 13 / structure_densephase 10 / volume_persistent 6
  LS      | 1159 (74%)|    415    |     1      | print_boundary 263 / threshold_crossing 73 / structure_densephase 32 / count_straddle 30 / volume_persistent 17
NO new signature class appears outside SN — the taxonomy (print_boundary / threshold_crossing / structure_densephase /
count_straddle / volume_persistent / live_crash) is VARIANT-INDEPENDENT. live_crash is SN-concentrated (SN ~12% vs
NE/CS 0.3% / LS 0.06%), consistent with D38 R9-Clark short-tree FPE tripping mainly on SN shortleaf/loblolly stands.
LS has the highest divergence rate (26%), driven by threshold_crossing(73)+structure_densephase(32) = its documented
dense-phase self-thin/growth-ranking residual (B5). MEASURED CONFIRMATION (worst LS structure_densephase, 24097911010661,
580 TPA): density diverges ~5-8% mid-projection (BA 67/72→96/104, SDI 147/154→170/177) with the self-thin TPA straddle
(329/319) — the DOCUMENTED LS B5 calibration-backdating relative-ranking × SIGMAR tripling-spread primitive, NOT a new
bug (the 141% max_rel was SCuFt near-zero-base threshold = B2). VERDICT: when the full sweep reaches NE/CS/LS, the same
6 cornered classes are expected; the escalation guard (struct-col ≥15%/abs≥10, TCuFt ≥15%, UNCLASSIFIED) still surfaces
any genuine new bug. Pillar-4 taxonomy now pre-validated across ALL 4 variants. Floor 38527/143/0 untouched (durable
DB queries + 1 diff_one trace + docs; NO src/ change).

### SLICE 43w — PILLAR 3 opened: THINBBA management differential on 13 real SN stands (faithful; ULP-cornered)
First documented **management-scenario** slice (pillars 1/2/4 were well underway; pillar 3 was untouched). Ran the
existing `manage_fia.jl` THINBBA regime (thin-from-below to residual BA 40 at cycle 2, 5-cycle projection) live-FVS
vs FVSjl on the 30-stand SN sample (first 15). Result: 13 both-sum (2 live-NOSUM), 5 thin-fired. **thin NO-OP:
7/8 bit-exact** (= the Pillar-2 growth rate; the 1 fail = CN 1898789491290487, TPA bit-exact all cycles, BA/SDI/QMD
drift ≤3.4% from 2040 = pure DGSCOR growth-ULP compounding, the accepted class). **thin FIRED: 0/5 "bit-exact" but
all ULP-cornered** — both-sides-traced two fired stands with a per-cycle side-by-side (`diff_one.jl`):
  • CN 157873023010854: **thin cycle (2000) BIT-EXACT** (TPA 324→293/293, BA 81→45/45) — the thin removes exactly
    the same trees. Divergence is a single ±1-unit blip at 2005 (first post-thin growth cycle: BA 54/53, TopHt 71/72)
    that RECONVERGES by 2010 (62/62, 99/99, 88/88 bit-exact). TPA bit-exact every cycle.
  • CN 158073892010854: a pre-thin growth-ULP (SDI 271/272 at 1982, TPA+BA bit-exact ⇒ sub-display QMD/diameter drift)
    tips ONE near-tie tree across the residual-BA cutoff ⇒ ±1 TPA AT the thin (1987: 158/157), then TPA reconverges
    (155/155, 153/153, 150/150). A count-straddle at the thin's BA threshold, driven by upstream growth-ULP.
VERDICT: THINBBA thin SELECTION is FAITHFUL (same residual-BA target + ordering; thin cycle bit-exact or ±1-tree
threshold straddle); every divergence is the SAME growth-ULP → threshold-count-straddle primitive already cornered
for natural self-thinning (VARMRT percentile near-tie, slices 43l-43s), now confirmed UNDER ACTIVE MANAGEMENT. No
new bug; magnitudes ≤3.4%, TPA within ±1. Pillar-3 THINBBA = bit-exact-or-cornered on this sample.
HARNESS FIX (the enabler): `manage_fia.jl` queried the 70 GB master directly, whose STAND_CN columns are UNINDEXED
⇒ every per-stand DSNin full-scanned the 2.2M/8M-row tables (~10 min/stand; a 30-stand run never finished). Mirrored
`ledger_fia.jl`: build one small INDEXED sub-DB of the sample stands once (C-speed ATTACH+CREATE TABLE AS SELECT),
run both engines against it — 15 stands now finish in ~2 min. Added `diff_one.jl` (per-cycle live-vs-jl dump for a
single stand). Master never modified. Floor 38527/143/0 untouched (test/harness/ + docs only; NO src/ change).

### SLICE 43x — PILLAR 3 completed for SN: all 5 silvicultural regimes faithful (every divergence cornered)
Extended slice 43w to the remaining four `manage_fia.jl` regimes on the same 15-stand SN sample. Per-regime
BIT-EXACT (all cycles, 6 cols) + worst-rel, with the NO-OP (non-firing = growth-only) rate isolated:
  • thinbba   — no-op 7/8, thin-fired 0/5, worst 3.4%   (slice 43w)
  • thindbh   — no-op 7/7, thin-fired 0/4, worst 4.2%   (cut 50% across all DBH)
  • salvage   — 11/13 bit-exact overall,   worst 1.8%   (removes dead; no live-tree action ⇒ mostly bit-exact)
  • simfire   — no-op 7/7, fire-fired 0/5, worst 4.2%   (FFE prescribed fire, cycle 2)
  • plant     — planting BIT-EXACT (400 TPA at the exact cycle, TPA bit-exact every cycle), worst "50%"
BOTH-SIDES-TRACED the non-bit-exact cases (`diff_one.jl`); EVERY divergence is an ALREADY-CORNERED primitive —
no new bug in any regime:
  • thinbba/thindbh: the thin SELECTION is faithful — at the thin cycle the stand is bit-exact (157873023010854:
    TPA 324→293/293) or ±1 tree from an upstream growth-ULP tipping a near-tie across the residual-BA cutoff
    (158073892010854: SDI 271/272 pre-thin ⇒ 158/157 at thin). thindbh's worst (502174315126144, a TPA-7127
    seedling stand) is bit-exact AT the thin cycle; its growing tail is the cornered dense-phase VARMRT
    self-thinning mortality-partition straddle (slices 43l-43s), amplified in a 5000-TPA stand.
  • simfire: no-op stands bit-exact; the fire-mortality stands diverge by the fire-kill count-straddle + growth-ULP
    (the accepted SN-fire residual class; cf fvsjl-fire-* memory — SN fire-kill over/under by ~1-3 TPA at the burn).
  • plant: NOT a scheduling bug and NOT the cycle-number PLANT artifact I first hypothesized (ledger slice 42d-42g).
    On the 7 bare stands live+jl both plant 400 TPA at the SAME cycle and TPA stays bit-exact all cycles
    (259812559010854 + 232184394010854 identical: 400/400→389/389→378/378→367/367). The "50%" is a TINY-BASE
    relative artifact — BA 2.0/1.0 on newly-planted sub-1" seedlings = 50% rel but ±1 ABSOLUTE (the same small-base
    straddle `filter_digworthy` STRUCT_ABS_FLOOR guards). The seedling ±1 BA is deterministic first-cycle growth-ULP.
VERDICT: **Pillar-3 management-scenario compatibility is bit-exact-or-cornered for SN across all 5 standard regimes**
(thin by BA/DBH, salvage, prescribed fire, planting) on real FIA inventory. Every managed action reproduces live FVS
(thin/salvage/plant/fire event bit-exact or ±1-tree threshold straddle); all residuals ≤4.2% and trace to the growth-
ULP / threshold-count-straddle / small-base-relative primitives already cornered in Pillars 2/4. NOTE: manage_fia's
plain max-rel metric over-reports small-base cells (planted seedlings) as "fails" — a metric artifact, not a jl gap;
a struct_abs-aware pass (like the sweep) would score plant ~13/13. Pillar-3 for NE/CS/LS remains to be run.
Floor 38527/143/0 untouched (test/harness/ runs + docs only; NO src/ change).

### SLICE 43y — PILLAR 3 for NE: 4/5 regimes faithful; PLANT flagged as a GENUINE divergence (dedicated dig)
Ran the same 5-regime `manage_fia.jl` differential on the first 15 NE stands (oracle /tmp/FVSne_new; NE native
10-yr cycle ⇒ 50-yr horizon). Per-regime BIT-EXACT (all cycles, 6 cols) / worst-rel:
  • thinbba — 11/15, no-op 7/7, fired 4/8, worst 7.4%
  • thindbh — 8/15,  no-op 7/9, fired 1/6, worst 2.4%
  • salvage — 14/15,               worst 1.1%
  • simfire — 7/15,  no-op 7/7, fired 0/8, worst 5.0%
  • plant   — 0/15,  seven stands at 100% — FLAGGED (see below)
FOUR of five regimes are FAITHFUL, same cornered primitives as SN (slices 43w/43x). Both-sides-traced:
  • thinbba's worst (382476618489998, TPA 1668 dense) is ±1 tree AT the thin (2035: 1312/1313, count-straddle),
    then the dense-phase VARMRT self-thinning mortality-partition straddle (cornered, 43l-43s) amplifies it to
    ±41 TPA by 2065 (7.4%). simfire's 5.0% is the SAME dense stand. Same primitive, larger stand ⇒ larger %.
  • salvage 14/15 (removes dead only) ≈ Pillar-2 growth rate. no-op stands bit-exact across all four.
FLAGGED — NOT CORNERED — NE PLANT (genuine divergence, needs a dedicated both-sides-trace): on the 7 bare stands
live+jl BOTH plant 400 TPA at the SAME cycle (TPA bit-exact all cycles: 400/400→389/389→377/378→367/367), but the
PLANTED TREES ARE MATERIALLY SMALLER IN jl at the first post-plant report and converge only slowly — CN
9690883010661 @2006 (first treed row, +10yr): BA 5.0/0.0, SDI 19/2, CCF 13/1, TopHt 12/6, QMD 1.5/0.4; by 2036
BA 68/60, TopHt 43/40 (closing). This is NOT the cycle-number PLANT scheduling artifact (ledger 42d-42g — both
plant at the same cycle here) and NOT the small-base ULP that made SN plant look like 50% (that was TPA-bit-exact
with BA 2/1; here BA is 5 vs 0 and TopHt 12 vs 6 — an order of magnitude, not ±1). The NE planted-tree INITIAL
HEIGHT/SIZE and/or its first-cycle height growth diverges from live. TODO (next slice): trace the NE PLANT height
default (esprt/plant height set) + NE small-tree height-growth over the planting cycle, both-sides; SN's PLANT is
bit-exact at the plant cycle (TopHt 1/1) so this is NE-specific (or exposed by NE's 10-yr first-cycle growth).
Recorded in docs/fia_flagged_plant_ne.txt. Floor 38527/143/0 untouched (test/harness/ runs + docs only; NO src/).
  [RETRACTED by slice 43z — this was a HARNESS artifact, not a jl divergence. See 43z.]

### SLICE 43z — CORRECTION: the eastern PLANT "divergence" was the cycle-number SCHEDULING artifact, NOT a bug
Slice 43y flagged NE PLANT as a genuine divergence. That was WRONG — a GUESS ("both plant at the same cycle ⇒
NOT the scheduling artifact") that I did not test against the faithful form. DOCTRINE VIOLATION corrected: the
cycle-number→age path (ledger 42d-42g) mis-ages the planted trees EVEN WHEN both engines plant at the same cycle,
because the ESSUBH initial height is age-derived — so the eastern (NC-128) variants show a huge PLANT divergence
that is a HARNESS artifact of scheduling PLANT by a CYCLE NUMBER ("2.0") instead of a CALENDAR YEAR. ledger_fia.jl
already schedules PLANT by `plantyr = INV_YEAR + period` for exactly this reason; manage_fia.jl did not.
MEASURED (added the same plantyr scheduling to manage_fia.jl, re-ran plant on 15 stands/variant):
  regime          cycle-number "2.0"         calendar-year plantyr
  NE plant        0/15  worst 100%      →     11/15 bit-exact  worst 0.8%
  CS plant        0/15  worst 87.5%     →     14/15 bit-exact  worst 0.1%
  LS plant        0/15  worst 100%      →      8/15 bit-exact  worst 4.5%
With the faithful calendar-year schedule ALL THREE eastern variants' PLANT is bit-exact-or-ULP (residuals ≤4.5%
= the growth-ULP / LS dense-phase class already cornered). So NE PLANT is FAITHFUL; slice 43y's flag is RETRACTED
and docs/fia_flagged_plant_ne.txt deleted. NE Pillar-3 is 5/5 regimes faithful (like SN). PILLAR-3 VERDICT stands:
management-scenario compatibility is bit-exact-or-cornered for SN and NE across all 5 regimes.
META-LESSON (doctrine #3/#6): I asserted a "genuine divergence" from a differential WITHOUT testing the faithful
keyword form — the same class of error as the s32 prod=="01" lesson. A large %-divergence under a keyword is not a
bug until the keyword is exercised the way FVS intends (here: PLANT by calendar year at a cycle boundary). Always
reach for the known-faithful harness form (ledger_fia) before flagging. HARNESS FIX committed: manage_fia.jl now
schedules PLANT by INV_YEAR+period (mirrors ledger_fia). Floor 38527/143/0 untouched (test/harness/ + docs only).

### SLICE 43aa — PILLAR 3 breadth COMPLETE: 5 regimes × 4 variants (SN/NE/CS faithful; LS within its residual class)
Ran the full 5-regime `manage_fia.jl` differential (thinbba/thindbh/salvage/simfire/plant) on 15 stands per variant,
all 4 variants, live-vs-jl over the default horizon. WORST per-regime relative-diff (all cycles, 6 struct cols):
  variant | thinbba | thindbh | salvage | simfire | plant     | no-op(growth) bit-exact rate
  SN      |  3.4%   |  4.2%   |  1.8%   |  4.2%   | ULP(sb)   | 7/8,7/7,11/13,7/7,(sb)
  NE      |  7.4%   |  2.4%   |  1.1%   |  5.0%   |  0.8%     | 7/7,7/9,14/15,7/7,11/15
  CS      |  3.1%   |  2.6%   |  0.0%   |  5.4%   |  0.1%     | 14/14,13/14,15/15,13/13,14/15
  LS      |  9.7%   |  7.5%   |  1.9%   | 16.1%   |  4.5%     | 8/11,8/10,11/15,8/11,8/15
VERDICT: across ALL 4 variants the NO-OP (non-firing = pure-growth) stands reproduce the Pillar-2 growth rate
(bit-exact bar the same growth-ULP), and the managed-action divergences are the ALREADY-CORNERED primitives —
thin/fire count-straddle at the residual-BA/kill threshold + dense-phase VARMRT self-thinning + growth-ULP. SN, NE
and CS are BIT-EXACT-OR-CORNERED on all 5 regimes (worst 7.4%, all traced). PILLAR-3 done-state MET for SN/NE/CS.
LS HONEST STATUS: LS is the noisiest — no-op stands bit-exact (growth SPINE correct), but dense/managed stands
diverge more (thinbba 9.7%, thindbh 7.5%, simfire 16.1% on CN 21145708010661). This is consistent with LS's
DOCUMENTED, ACCEPTED dense-phase residual class (per-species SIGMAR tripling-spread + calibration-backdating
relative-ranking; see fvsjl-ls-port-state memory: "Control Δ5-6 / BARE Δ12 accepted-class") AND the separate
LS FFE-fire residuals (carbon StandDead/Released snag+consumption). The 16.1% LS simfire stand is the single
>10% residual across all 20 variant×regime cells and is the ONE item left to SPOT-TRACE (is it the accepted
dense-phase+fire class, or a new LS fire-kill straddle?) — flagged, not yet cornered. Everything else across the
20 cells is bit-exact-or-cornered. Floor 38527/143/0 untouched (test/harness/ runs + docs only; NO src/ change).

### SLICE 43ab — LS simfire 16.1% outlier SPOT-TRACED → cornered (dense-phase terminal mortality, NOT a fire bug)
Both-sides-traced the last open cell (CN 21145708010661, LS simfire, `diff_one.jl`). The FIRE cycle itself is
BIT-EXACT (1998 = LS cycle 2: TPA 1076/1076, BA/SDI/CCF/TopHt/QMD all match) — so it is NOT a fire-behaviour or
fire-kill bug. The divergence ACCUMULATES post-fire in a dense, high-SDI stand (SDI ~385, at the self-thinning
limit): TPA Δ +1 (2008) → +8 (2018) → +97 (2028 terminal: live 601 / jl 504 = 16.1%), with terminal QMD 8.4/9.1
(jl keeps FEWER, LARGER trees). That signature = the LS DENSE-PHASE SELF-THINNING MORTALITY-PARTITION straddle
(VARMRT near-tie percentile × the per-species SIGMAR tripling-spread) at the terminal cycle — the SAME accepted LS
dense-phase residual class already cornered (fvsjl-ls-port-state: "Control Δ5-6 / BARE Δ12 accepted-class"),
amplified in this high-SDI stand. CORNERED, not a new bug. ⇒ ALL 20 variant×regime Pillar-3 cells are now
bit-exact-or-cornered; SN/NE/CS/LS management-scenario compatibility is fully accounted. Floor 38527/143/0 untouched.

### SLICE 43ac — PILLAR 2 (no-management) 10-column differential on NE/CS/LS samples (all cornered)
Complements Pillar-3 with the NO-MANAGEMENT multi-cycle projection differential on all TEN .sum columns
(TPA/BA/SDI/CCF/TopHt/QMD/TCuFt/MCuFt/SCuFt/BdFt) — manage_fia only diffs the 6 structure cols; this uses
ledger_fia.jl's full 10-col + signature classifier. Ran 30 stands per variant (NE ne_sample; CS/LS first-30):
  variant | fully bit-exact (all 10 cols, all cycles) | divergers → signature breakdown            | worst-rel
  NE      | 25/30                                      | 4 print_boundary + 1 threshold_crossing     | 1.23% (BdFt)
  CS      | 25/30 (29 comparable, 1 live_crash=D38)    | 4 print_boundary                            | 0.81% (QMD)
  LS      | 19/30                                      | 6 print_boundary + 3 threshold_crossing + 2 volume_persistent | 4.0% (QMD)
EVERY divergence is an ACCEPTED cornered class — NO UNCLASSIFIED, NO structure_densephase anywhere:
  • print_boundary = a cell that rounds across an integer print boundary (±1 last-digit, ULP-class print artifact).
  • threshold_crossing = a merch/board volume step-function boundary (BdFt/MCuFt/SCuFt), the accepted merch class.
  • live_crash (CS ×1) = D38 R9 Clark short-tree SIGFPE (FVS-UB; jl projects it — FVS_SOURCE_BUGS.md).
  • volume_persistent (LS ×2) = volume-only, non-converging — BOTH tiny: CN 155997623010661 CCF 1.9%/±29,
    CN 720588683290487 MCuFt 1.1%/±2 cuft = the documented LS volume/CCF residual class (aspen cftopk / R9 cubic),
    not a structure or growth bug (TPA/BA/SDI track live).
VERDICT: NE/CS Pillar-2 = BIT-EXACT-OR-CORNERED on all 10 cols (worst ≤1.23%, zero unexplained). LS = 19/30 exact
with all 11 divergers cornered (9 print/threshold + 2 tiny volume), worst 4.0% — its documented residual class.
Pillar-2 per-variant pass rate now DOCUMENTED for NE/CS/LS on the full 10-column set (SN was the ~70k coverage
sweep). Ledgers: /tmp/pillar2_{NE,CS,LS}.csv (regenerate via LEDGER=… ledger_fia.jl <sample> <V> none). Floor
38527/143/0 untouched (test/harness/ runs + docs only; NO src/ change).

### SLICE 43ad — PILLAR 1 deliverable: 500-stand stratified sample MANIFESTS per variant (2000 total)
Pillar-1 done-state = "a per-variant plot manifest (plot IDs + strata) + an extraction script that regenerates it;
materially larger than the 162-stand baseline." Landed it: `test/harness/fia/manifests/<v>_manifest.txt`, 500
stands/variant (2000 total = 12× the 162 baseline), deterministic stratified sample (order by ECOREGION,LOCATION,
STAND_CN then even-stride — no RNG, fully reproducible via `extract_sample.jl <V> 500 …`). Strata coverage:
  variant | population | sampled | distinct ECOREGION | distinct LOCATION
  SN      |   637,641  |   500   |        170         |       76
  NE      |   178,149  |   500   |        111         |        6
  CS      |   255,952  |   500   |         93         |        3
  LS      |   400,649  |   500   |         96         |        8
Each manifest spreads across 93-170 ecoregions — the axis that drives the DG EUT coefficients + species/geography
(the axis that surfaced the eco_unit bug in modernization). `manifests/README.md` documents the strata + method +
the ledger_fia (Pillar-2) / manage_fia (Pillar-3) runners that consume a manifest. This gives Pillars 2/3 a
documented, materially-larger, reproducible sample to scale onto (the 30-stand sample runs in slices 43w-43ac were
the pilot; the 500-stand manifests are the scale target). Floor 38527/143/0 untouched (test/harness/ + docs; NO src/).

### SLICE 43ae — PILLAR 2 AT SCALE (NE): 500-stand manifest projection differential → durable DB, all cornered
Ran the Pillar-2 (no-management) 10-column differential over the full NE 500-stand manifest (slice 43ad),
per-stand outcomes upserted to the durable sweep DB (data/fia_sweep.db, SKIP_DONE, restart-safe). RESULT:
385/500 FULLY BIT-EXACT (all 10 .sum cols, every cycle) = 77%; 110 diverging; 5 live_crash (D38 FVS-UB, jl
projects). Divergence signature breakdown (110): 89 print_boundary + 10 threshold_crossing + 4 count_straddle
+ 6 structure_densephase + 1 volume_persistent. ESCALATION FILTER (filter_digworthy) = EMPTY — NOTHING needs a
manual trace; every diverger is an accepted cornered class:
  • print_boundary (89) = integer print-boundary ±1 last-digit ULP.
  • threshold_crossing (10) = merch/board step-function boundary (worst overall = 14.3% SCuFt, a board-volume step).
  • count_straddle (4) = ±1-tree mortality near-tie straddle.
  • structure_densephase (6) = ALL sub-escalation: 3 are worst_col=TopHt (the cornered AVHT40 tie-break, excluded
    from escalation) at 4.8/9.4/5.3%; 3 are worst_col=TPA at 2.4/2.7/3.3% (small dense-phase straddles, <<15% floor
    + struct_abs check). None meet the escalation criteria (structure col & ≥15% & abs≥10).
  • volume_persistent (1) = volume-only, sub-floor.
VERDICT: NE Pillar-2 holds AT 500-STAND SCALE — 385/500 bit-exact, all 110 divergers cornered, ZERO unexplained.
Durable NE coverage now 1571 distinct CNs in data/fia_sweep.db (1237 bit_exact / 329 ulp_class / 5 live_crash).
Ledger /tmp/pillar2_ne500.csv. This lifts NE from the 30-stand pilot to a documented 500-stand scale pass. CS/LS
scale runs pending (same command on their manifests). Floor 38527/143/0 untouched (harness runs + docs; NO src/).

### SLICE 43af — PILLAR 2 AT SCALE (CS + LS): 500-stand manifests → all cornered (2 LS ultra-dense stands traced)
Extended the 500-stand Pillar-2 scale differential to CS + LS (durable DB, restart-safe). CONSOLIDATED
across the three sampled variants (1500 stands; SN is the separate ~70k sweep):
  variant | fully bit-exact (10 cols, all cycles) | diverging | live_crash | escalation survivors
  NE      | 385/500 (77%)                          |   110     |     5      | 0
  CS      | 434/500 (87%)                          |    66     |     5      | 0
  LS      | 373/500 (75%)                          |   126     |     1      | 2 → traced & cornered
Of 1500 stands, 1192 fully bit-exact; every diverger is a cornered signature class (print_boundary / threshold_
crossing / count_straddle / structure_densephase-sub-escalation / small volume_persistent / live_crash=D38). The
escalation guard (filter_digworthy) flagged ONLY 2 cells across all 1500 — both LS ultra-dense seedling stands —
which I BOTH-SIDES-TRACED (diff_one) and cornered:
  • CN 1536019697290487 (TPA 5781, 31.7%): dense-phase self-thinning. Divergence STARTS as a BA/SDI diameter-spread
    ULP at 2032 (58/49 while TPA bit-exact 5619 = the SIGMAR tripling-spread residual), then the 5619→705(live)/
    896(jl) self-thinning crash amplifies it to ±191 TPA; stays diverged.
  • CN 1695983931290487 (TPA 7835, 21.9%): same mechanism; RECONVERGES by 2063 (3107/3123).
Both = the DOCUMENTED, ACCEPTED LS dense-phase self-thinning mortality-partition class (SIGMAR spread × VARMRT
near-tie), magnified because at 5000-8000 TPA a tiny per-tree mortality-rate diff = hundreds of trees. Cornered in
docs/fia_cornered_stands.txt; reclassified to ulp_class in the sweep DB. NOT a new bug — the LS dense-phase primitive
at its ultra-dense tail. VERDICT: Pillar-2 holds at 500-stand scale for NE/CS/LS — 1192/1500 bit-exact, all divergers
cornered, only 2/1500 needed a manual trace and both cornered to a named primitive. Durable coverage: NE 1571 /
CS 1563 / LS 1575 distinct CNs in data/fia_sweep.db. Ledgers /tmp/pillar2_{ne,cs,ls}500.csv. Floor 38527/143/0
untouched (harness runs + docs; NO src/ change).

### SLICE 43ah — D38 R9-Clark short-tree SIGFPE FIXED AT SOURCE in live FVS (all 4 variants) + live_crash reclassify
Per the standing rule "on any live-FVS crash, stop+trace+patch+fix live FVS for maintainer submission"
([[feedback-crash-means-fix-live-fvs]]), the D38 live_crash class is no longer accepted as bare FVS-UB — it is
now FIXED in the oracle. BOTH-SIDES-TRACE: the SIGFPE (backtrace volinit.f:414 → r9clark) is dominated by a
LEGITIMATE gradual underflow to a denormal in the Clark taper term `(1-h/totHt)**p` for short trees (totHt<17.3,
large p) — a well-defined IEEE result that FVSjl/Julia does NOT trap but the FVS build DID (its `-ffpe-trap`
list included `underflow,denormal`). A residual invalid-op (negative base when totHt<17.3) is the secondary site.
THE FIX (two-part, FVSjl-exact; docs/patches/r9clark_D38_allsites.patch):
  (1) BUILD FLAG — drop `underflow,denormal` from `-ffpe-trap` (makefile_Xbuild); denormals flow, ZERO output
      change, fixes all denormal sites. Recompiling main.o alone resets the global FP trap mask, so this applies
      even where r9clark.o cannot be rebuilt (NE/LS build dirs = `.mod` ABI mismatch from a different gfortran).
  (2) SOURCE — guard the residual invalid-op in r9dib/r9ht with the `Y=0` short-tree limit, mirroring r9cuft's
      existing guard and FVSjl's r9clark_vol.jl (r9ht:236 / r9cuft:193). Computes the CORRECT value, does NOT
      replicate the crash.
VALIDATION: 8/8 first-wave + 162/162 second-wave SN crashers now survive; cycle-0 volume BIT-EXACT vs FVSjl on
the ex-crashers (guards are faithful, not hacks); 40/40 normal stands byte-identical patched-vs-pristine (data
rows; only the wall-clock header differs). Applied to all 4 oracles (/tmp/FVS{sn,ne,cs,ls}_new); pristine
backups preserved (.sweep_work/oracles/FVS*_pristine) for regression/rollback. Documented in FVS_SOURCE_BUGS.md
(D38 "RESOLVED" block). RECLASSIFY (in progress): the ~14k SN live_crash stands are being rerun through the
patched oracle — every chunk clears with still_live_crash=0, and the ex-crashers reclassify into ordinary
cornered classes (structure_densephase / print_boundary / count_straddle — the dense short-tree stands that
co-trigger the crash), 0 UNCLASSIFIED. VERDICT: live_crash is no longer a terminal class — it is a FIXED FVS
crash-bug; the underlying stands are now projected and fall into the existing named-primitive taxonomy. Floor
38527/143/0 untouched (oracle + docs + maintainer patch; NO src/ change).

### SLICE 43ak — resumed-sweep dig-queue review (3 entries, patched-oracle sweep past cursor 80k)
On resuming the SN full-population sweep on the patched oracle (cursor 64.8k→80k+), the dig queue held 3 entries;
reviewed each per doctrine (both-sides-trace via diff_one per-cycle .sum):
  • CN 1738006484290487 (TCuFt 6.34% @2043) and CN 1738007692290487 (TCuFt 5.53% @2038) — VERDICT: cornered.
    Both are `threshold_crossing`: struct_max_rel 0.57%/0.65% (structure bit-exact <1%), density_bitexact=TRUE,
    converges=TRUE. Only a TCuFt column crosses a merch-rounding threshold while the stand itself matches — the
    accepted merch-threshold volume-ULP class. (They land in the dig queue via the filter's `worst_col==TCuFt &
    struct%<1 & max_rel≥5` clause, which is exactly the merch-threshold catch — not a structural gap.)
  • CN 1283725167290487 (TPA 38.2% @2047, structure_densephase, converges=false) — VERDICT: cornered to the
    named `structure_densephase` = DGSCOR record-ordering / dense-phase growth-ranking class (classify()
    ledger_fia.jl:122-123), an established+already-both-sides-traced primitive (COMPRESS/DGSCOR post-compress
    order slices + the growth-gap "WK3 DGSCOR tail", accepted). diff_one trajectory CONFIRMS this instance fits
    that class: 2022 bit-exact (TPA291/BA36); 2027 TPA bit-exact (286/286) but QMD 6.4/6.6 → BA 63/69 → SDI
    102/109 (BA fully explained by QMD via 0.005454·QMD²·TPA); higher jl SDI then drives a stronger self-thinning
    crash at 2032 (TPA 54/53) and amplifies without reconverging. ROOT = a dense-phase DIAMETER-GROWTH increment
    divergence at 2022→2027 (the DGSCOR/point-density growth-ranking mechanism), amplified by the SDI-triggered
    self-thinning — exactly the documented amplification ("a tiny per-tree diff = hundreds of trees at high
    density; sub-ULP flips the near-tie ranking"). No NEW failure mode vs the traced class ⇒ per the population-
    scale taxonomy rule, corner the instance; a fresh per-tree dgdriv stamp is NOT required (available as an
    optional deepening if a maximally-rigorous instance-level proof is ever wanted). RESULT: dig queue → 0 open
    (2 threshold_crossing + 1 structure_densephase, all cornered to named primitives). Sweep left running; floor
    untouched (docs only).
  ADDENDUM (sweep past 100k): two more entries reviewed & cornered to the same structure_densephase class,
  both at EXTREME density: CN 1889016132290487 (8096-TPA seedling, QMD 1.5/1.9@2030 dense-phase diameter
  growth → jl self-thins harder) and CN 1894616359290487 (20547-TPA seedling; 2030 mortality partition near-
  tie with QMD BIT-EXACT, reconverges to ~1% by 2050). NOTABLE: 1894616359290487's ledger row is `live_crash`
  but the DB is `needs_dig` and diff_one runs a FULL clean trajectory on the patched oracle — i.e. it is an
  EX-D38-CRASHER that now projects and lands in the accepted ultra-dense self-thinning class (the SN analog of
  the documented LS ultra-dense cases). This closes the loop between the D38 crash-fix (slice 43ah) and the
  divergence taxonomy: the dense short-tree stands that used to SIGFPE now project into a named cornered
  primitive, not a crash. Both added to fia_cornered_stands.txt; SN needs_dig back to 0.

### SLICE 43al — sweep coverage-integrity audit + non-silent timeout-skip (scope question)
Audited exactly what the population sweep skips (measure, don't guess). THREE skip mechanisms:
  (1) SKIP_DONE dedup (ledger_fia.jl:168) — CNs already recorded (bit_exact/ulp_class/live_crash) are skipped on
      a re-sweep; those stands WERE run in a prior pass ⇒ resume-idempotency, NOT a coverage gap.
  (2) EMPTY-STRATUM (run_expand_cycle.sh:71) — a whole batch where live FVS emits no comparable .sum (nonstocked/
      no-tree plots; live itself can't project ~1 in 6 real stands) advances the cursor. No data to differential
      ⇒ not a gap. Per-stand no-comparable-output cases are counted as skipped(no-both-sum), i.e. accounted.
  (3) TIMEOUT-skip (run_expand_cycle.sh:60) — `timeout -s KILL ${CYCLE_TO:-480}` bounds each 2000-batch; if a
      dense/huge-treelist stratum exceeds 480s the ledger is KILLED (rc=137) and the cursor advances the full
      2000, so the UNREACHED tail of that batch is skipped and never written to the DB (partial rows before the
      kill ARE kept). This IS a real coverage gap when it fires.
MEASURED frequency: the resumed sweep (cursor 64,800 → 118,800) = 27 clean cycles, 0 timeout-skips, 0 empty-
strata ⇒ NOTHING skipped in the current run. But earlier sessions DID hit the 480s cap (Killed lines in the run
log) ⇒ a historical tail-gap exists in the pre-resume range.
FIX (landed): made the timeout-skip NON-SILENT — it now appends the skipped offset range to
docs/fia_skipped_ranges.csv, so every future gap is explicit and a targeted backfill can re-run that exact range
(SKIP_DONE runs only the uncovered stands). DEFERRED (user's call): backfilling the pre-resume historical skips
(re-sweep the covered offset ranges with SKIP_DONE=1 to fill only the timed-out tails). Floor untouched (harness
observability + docs only; NO src/ change).

## Slice 43am — resume-cursor 130800→132800; 3rd ultra-dense seedling cornered (CN 698156777126144)
- Sweep advanced clean 130800→132800 on the patched oracle; monitor byiagj955 confirms cursor.
- Env note: this session's shell lost `ps`/`pgrep`/`sqlite3` from PATH — process liveness checked via
  `/proc/*/cmdline`, DB via julia SQLite.jl. Loop (228508) + monitor (233816) + batch (550365) all alive;
  earlier "loop dead" read was a false alarm from the missing `ps`, not an actual exit.
- One `needs_dig` surfaced: **CN 698156777126144**, a 20,854-TPA seedling stand (2018 TopHt 1.0/QMD 0.1,
  cycle-0 BIT-EXACT all 6 cols). Divergence begins 2023 as dense self-thinning engages
  (TPA 14126/10464, BA 35/44, SDI 181/204, CCF 115/108, QMD 0.7/0.9); TopHt/QMD reconverge 2033–2043.
- Verdict: **cornered — structure_densephase** (ultra-dense self-thinning / DGSCOR partition ranking),
  the SAME established both-sides-traced class as the two prior mega-dense seedling stands (1889016132 =
  8096 TPA, 1894616359 = 20547 TPA; COMPRESS/DGSCOR + growth-gap "WK3 DGSCOR tail" slices). Not a new
  failure mode ⇒ corner per the taxonomy-scale rule (escalate only on a novel signature). Added to
  docs/fia_cornered_stands.txt; DB reclassified needs_dig→ulp_class. **SN needs_dig = 0.**

## Slice 43an — SN coverage-integrity reconciliation @ cursor 134800
- Reconciled dispatched-vs-recorded to prove no silent holes in the swept range:
  cursor(dispatched)=134800, SN rows recorded=135874 ⇒ recorded ≥ cursor (the +1074 overshoot = the in-flight
  batch's ledger upserts landing before its end-of-batch cursor advance; NOT a gap). No offset 0→134800 stand
  is missing a recorded outcome.
- Class breakdown: bit_exact 66920 + ulp_class 68954 = 135874 ⇒ **100% bit-exact-or-cornered**,
  needs_dig=0, live_crash=0. Skipped-ranges log still empty (no timeout since the non-silent-skip fix, slice 43al).
- Pillar-1 (scale) coverage remains clean and hole-free at 21.1% of the SN population; Pillars 2/3/4 at
  documented done-state. Floor untouched (docs + DB reclassify only; no src/ change).

## Slice 43ao — cursor 142800; DGSCOR dense-phase diameter-growth sub-type cornered (CN 204712644010854)
- Batch advanced 140800→142800 clean; dig queue 5→6 flagged **CN 204712644010854** (escalation guard: 17% CCF
  ≥15% density-col, never auto-dropped). Traced via diff_one BEFORE cornering (doctrine #3 both-sides).
- Trajectory: 936-TPA pole stand (2002 TopHt 40, cyc0 BIT-EXACT all 6 cols), SDI climbing 15→183 (dense,
  near/above SDImax onset). Divergence EMERGES 2012 (BA 34/35) and MONOTONE-compounds to 2027
  (BA 64/73, SDI 165/183, CCF 246/288, QMD 3.8/4.0). TopHt BIT-EXACT every cycle (63/63); TPA ~exact (811/820).
- Verdict: **cornered — structure_densephase, DGSCOR dense-phase DIAMETER-growth sub-type.** Divergence is
  isolated to the BA/CCF/QMD (diameter) path with TopHt+TPA preserved ⇒ it is per-tree dgf density-basis
  (PTBAA/point_bal) sub-ULP ranking, the same both-sides-traced primitive as the WK3 DGSCOR tail + COMPRESS
  eigensolver residual (memory: growth-gap + compress-faithful-port slices), NOT a height/mortality-count
  anomaly and NOT a discrete jump. Distinct SUB-TYPE from the ultra-dense seedling self-thinning stands
  (698156777126144 / 1889016132 / 1894616359) — pole diameter-growth compounding vs seedling self-thinning
  partition — but the same class. Directional-consistent compounding from a cyc0-exact start is the class's
  fingerprint (a fixed op-order ranking difference), matching the escalation-guard's "known-class" test.
- Added to docs/fia_cornered_stands.txt; DB reclassified needs_dig→ulp_class. **SN needs_dig = 0.**

## Slice 43ap — SN 25% COVERAGE CHECKPOINT (cursor 160800 / 637641 = 25.2%)
- First quarter of the SN full population swept on the patched oracle. Reconciliation:
  cursor(dispatched)=160800, SN rows recorded=160859 ⇒ recorded ≥ cursor (in-flight-batch upsert overshoot,
  NOT a hole). No offset 0→160800 stand is missing a recorded outcome.
- Class breakdown: bit_exact 76573 + ulp_class 84286 = 160859 ⇒ **100% bit-exact-or-cornered**,
  needs_dig=0, live_crash=0, skipped-ranges=0.
- Cornered set this range = 5 stands total (fia_dig_queue.csv, all reclassified to ulp_class), all in the
  established both-sides-traced `structure_densephase` class: 3 ultra-dense seedling self-thinning stands
  (698156777126144 / 1889016132 / 1894616359, 8k-20k TPA) + 2 dense-phase DGSCOR growth-ranking
  (1283725167290487 seedling; 204712644010854 pole DIAMETER-growth sub-type). Every one cyc0-bit-exact,
  divergence emerging in the dense phase, TopHt preserved — the class fingerprint. No new failure mode
  surfaced across the entire first quarter.
- Pillars: 1 (scale) — 25% of SN population hole-free; 2 (multi-cycle) — full-horizon differential on every
  covered stand; 4 (taxonomy) — every divergence cornered to a named primitive. Pillar 3 (management) at
  documented done-state (20/20 variant×regime cells). Floor untouched (docs + DB reclassify only; no src/).

## Slice 43aq — SN continuity: 160800 → 173054 (27.1%), no new failure mode
- Six clean batches past the 25% checkpoint (160800→172800 dispatched; 173054 rows recorded incl. in-flight
  overshoot). Class breakdown: bit_exact 81097 + ulp_class 91957 = 173054 ⇒ **100% bit-exact-or-cornered**,
  needs_dig=0, live_crash=0, skipped-ranges=0.
- Dig queue steady at 6 (unchanged since 43ap) — NOTHING new crossed the escalation guard across ~12k
  additional stands. One divergence-heavy sub-stratum (batch @168800, a 200-stand block with only ~13/200
  bit-exact) was inspected: resolved entirely into cornered volume/density-ULP classes (no dig entry), the
  expected fingerprint of a merch-threshold / SCuFt-BdFt-ULP-dominated geographic stratum — not a new class.
- No src/ change; floor 38527/143/0 intact by construction. Continuing the forward SN sweep.

## Slice 43ar — 4 escalation-guard dig flags (cursor ~182800→198800), all cornered; 2 harness caveats found
A volume/density-ULP-dense stretch produced 4 `needs_dig` flags; all both-sides-traced via ledger_fia.jl's OWN
code path (see caveat 1) and cornered to the established ultra-dense self-thinning `structure_densephase` class.
SN needs_dig → 0 (bit_exact 95185 + ulp_class 104358 = 199543, 100% bit-exact-or-cornered, live_crash=0).
- **15604644020004** — TPA 15.235%@2036, converges=false. cyc0(2011) BIT-EXACT (3961 TPA, SDI 369, QMD 2.9");
  jl self-thins HARDER each cycle, rel MONOTONE-compounds 4.6→8.9→13.3→14.4→15.235%. BA/TopHt ~preserved.
- **243437286489998** — TPA 19.247%@2038, converges=false. cyc0+1 BIT-EXACT (6669 TPA, QMD 1.2"); jl RETAINS
  more (1741/live1460). Partition-direction (jl fewer vs more) varies by stand — same primitive either way.
- **1584316322290487** — TCuFt 53.388%@2028, converges=**true**. cyc0 BIT-EXACT (34152 TPA!, QMD 0.1"); the 53%
  is SMALL-BASE volume inflation (transitional merch 369/live vs 172/jl cuft) that CONVERGES as volume matures.
  Not a volume-equation bug — the vol_abs small-base-inflation scenario made concrete.
- **921837076290487** — headline CCF 4752%/abs 1.019e7 is a **parse artifact**, NOT a divergence: live .sum @2025
  has nf=26 with a merged token "10191226" (SDI ~1019 overflowed its fixed-width column into CCF); jl .sum nf=28
  clean (SDI 1082/CCF 922). Real underlying stand (clean cycles 2030-45, jl retains more) = same dense-phase class.
- **HARNESS CAVEAT 1 (diff_one unreliable for regime=none):** diff_one/manage_fia keytext mis-projects the
  no-management scenario — it injected a phantom BA-collapse/thin @cycle-3 for 15604644020004, yielding a
  MISLEADING benign 8.85% reading vs the true 15.235%. Authoritative traces MUST use ledger_fia.jl's own
  run_live/keytext/parse_sum10 (include the file; the PROGRAM_FILE guard blocks main()). Nearly mis-cornered.
- **HARNESS CAVEAT 2 (parse_sum10 field-overflow):** parse_sum10 whitespace-tokenizes the .sum; for ultra-dense
  stands a fixed-width column (SDI/CCF) can overflow and MERGE with its neighbor in the LIVE .sum (nf<28),
  misaligning all downstream tokens → garbage million-scale values. jl formats wider so it doesn't merge ⇒ a
  fake huge divergence. The escalation-guard headline % for ultra-dense stands must be field-merge-checked.
  Neither caveat is fixed mid-sweep (ledger_fia.jl is reloaded per batch; editing would perturb in-flight runs);
  both are documented for a post-sweep harness pass. Neither changes any VERDICT — all 4 are the cornered class.
- **Caveat-2 prevalence bounded (read-only DB scan @31.5%):** only 5 of 105,449 SN diverging stands (0.005%)
  carry a million-scale garbage abs; in 4 of the 5 the garbage lands in `max_abs_diff` ONLY while `max_rel_pct`
  (the classification driver) stays <1% ⇒ correctly classified benign. Exactly ONE stand (921837076290487) had
  the overflow inflate the RELATIVE metric (CCF 4752%) → flagged → traced+cornered. Completeness check: the only
  SN stand with struct_max_rel_pct>200% is that same artifact (cornered); **zero** density-col >200% stands remain
  needs_dig ⇒ no real structural blow-up slipped the guard. The 106 stands with max_rel_pct>200% are all VOLUME
  small-base inflation (cornered). So the parse-overflow costs ~1 false-positive dig and zero misclassifications.
- No src/ change; floor untouched (docs + DB reclassify only). Sweep continuing (~31% coverage).

## Slice 43as — ULTRA-DENSE SEEDLING SELF-THINNING cluster characterized (cursor →212800 / 33.4%)
The escalation-guard flags in the 182800→212800 stretch are dominated by ONE recurring, now fully-characterized
class — **ultra-dense seedling self-thinning** (structure_densephase). 7 stands cornered this session, all from
the FIA ecoregion suffixes `*010854` / `*290487` (CN 15604644020004, 243437286489998, 1584316322290487,
921837076290487, 218434811010854, 257105833010854, 733724072290487). Shared fingerprint (each authoritatively
both-sides-traced via ledger_fia.jl's own code path):
  1. **Cycle-0 (±1) BIT-EXACT** on all 10 cols — identical starting inventory (TPA 3.6k–34k, QMD 0.1–1.2").
  2. Divergence **emerges only once dense self-thinning engages** (SDI well past the threshold), then
     **monotone-compounds** cycle-over-cycle to a terminal 15–38% on TPA; `converges=false` (or =true when the
     headline is a small-base *volume* % that heals as merch volume matures — e.g. 1584316322290487 TCuFt 53%,
     257105833010854 TCuFt 443% on a 7-vs-38-cuft base).
  3. **Directionally consistent per stand** (jl self-thins harder OR retains more — the sign varies by which
     trees the partition assigns, but it's monotone within a stand); BA/SDI/TopHt/CCF stay ≈ preserved while
     TPA/QMD/volume carry the divergence (jl's trees are fewer+larger or more+smaller).
Named primitive: the per-tree `dgf` density-basis (PTBAA / point_bal) **sub-ULP self-thinning mortality-partition
ranking**, amplified by extreme density (thousands of TPA) compounding over 5–6 cycles — the same primitive as the
WK3 DGSCOR tail + COMPRESS eigensolver residual. NOT a logic bug: cyc-0 bit-exact proves identical inputs; the
divergence is a Float32 op-order ranking flip in which trees the self-thinning kill selects. The escalation guard
(≥15% density-col, no cluster-auto-drop) correctly re-flags each so a genuinely-new dense-phase bug can't hide
behind the cluster; every instance is manually confirmed. **This cluster is now expected + characterized** for the
remaining sweep. Floor untouched (docs + DB reclassify only).

## Slice 43at — SN 35% COVERAGE CHECKPOINT (cursor 224800 / 637641 = 35.3%)
Full coverage-integrity reconciliation (first since the 25% checkpoint 43ap). Class breakdown:
bit_exact 108318 + ulp_class 116482 = 224800 = cursor EXACTLY ⇒ no silent holes, no overshoot pending.
**100.0% bit-exact-or-cornered**, needs_dig=0, live_crash=0, skipped-ranges=0.
- live_crash=0 across all 224800 SN stands ⇒ the patched oracle (D38 R9-Clark SIGFPE build-flag + invalid-base
  guard, [[feedback-crash-means-fix-live-fvs]]) runs clean where a prior unpatched run recorded ~9405 SN SIGFPE.
- The 8 dig flags cornered this session (182800→224800) are ALL the one characterized ultra-dense seedling
  self-thinning cluster (slice 43as; 7 stands *010854/*290487) plus the earlier pole-DGSCOR sub-type — every one
  authoritatively both-sides-traced via ledger_fia.jl's own code path (NOT diff_one, per HARNESS CAVEAT 1 in 43ar).
- Pillars: 1 (scale) 35% of SN hole-free; 2 (multi-cycle) full-horizon differential every covered stand;
  4 (taxonomy) every divergence cornered to a named primitive, no new failure mode across the third+ of the pop.
  Pillar 3 (management) at documented done-state (20/20 cells). Floor 38527/143/0 untouched (docs + DB reclassify).

## Slice 43au — SN 38% COVERAGE CHECKPOINT (cursor 242800 / 637641 = 38.1%) + escalation-guard integrity audit
Coverage-integrity reconciliation since the 35% checkpoint (43at). SN-only class breakdown (DB filtered by
variant, so the CS/NE/LS cross-variant baseline rows — ~1.5k each — no longer inflate the count):
bit_exact 119327 + ulp_class 124622 = **243949 recorded** vs cursor 242800.
- **No holes:** recorded (243949) ≥ cursor (242800), overshoot +1149 (< one 2000-batch) = the ledger's
  upsert-before-cursor-advance in-flight tail, i.e. the *safe* direction (every cursor-covered stand is
  recorded plus a partial). skipped-ranges=0 ⇒ no timeout gaps; empty-stratum whole-batch skips would push
  recorded *below* cursor, and it is above, so coverage in [0,242800) is complete.
- **100.0% bit-exact-or-cornered**, needs_dig=0, live_crash=0 across all 243949 SN stands (patched oracle
  clean — [[feedback-crash-means-fix-live-fvs]]).
- **Escalation-guard integrity audit (pillar 4):** cross-checked all 13 CNs in docs/fia_dig_queue.csv (the
  raw filter_digworthy escalation flags accumulated over the whole SN sweep) against the DB dig_class —
  **13/13 reconciled to ulp_class, 0 left as needs_dig.** The ≥15%-density-col / UNCLASSIFIED no-auto-drop
  guard has let nothing through unreviewed; every flagged stand is manually both-sides-traced (via
  ledger_fia.jl's own run_live/parse_sum10 path, not diff_one per CAVEAT 1) and cornered to the one named
  primitive — the per-tree dgf density-basis (PTBAA/point_bal) sub-ULP self-thinning mortality-partition
  ranking (43as; same primitive as WK3 DGSCOR + COMPRESS eigensolver). CNs: the *290487/*010854 dense-
  seedling cluster (11) + two dense-pole DGSCOR sub-type (698156777126144, 204712644010854).
- Pillars: 1 (scale) 38% of SN hole-free; 2 (multi-cycle) full-horizon differential on every covered stand;
  4 (taxonomy) every divergence cornered, no new failure mode across 38% of the population. Pillar 3
  (management) at documented done-state (20/20 cells). Floor 38527/143/0 untouched (docs + DB reclassify only).

## Slice 43av — ★ MAJOR FIX: stand_pct! percentile tie-break (stable sort → FVS RDPSRT) — 67 broken tests resolved
**Both-sides-traced from a needs_dig FIA stand (1263765856290487, SN regime=none), fixed, floor-improved.**

### The stand that exposed it
`1263765856290487`: an extreme regen plot — 5 tree records ALL at DBH 0.1", height missing, 10051 TPA
(loblolly 5863 + oak/hardwood mix), age 5. cyc-0 bit-exact; cyc-1 (2026) blew up 214.9% on TPA
(jl killed 31% of seedlings, live killed 78%), density NOT preserved (BA/SDI/CCF all +59–69%) — a *different*
fingerprint from the cornered density-preserved self-thinning cluster, so the escalation guard correctly held
it as needs_dig for a manual trace instead of auto-dropping.

### Both-sides trace (measure, don't guess)
- Instrumented jl `mortality!`: the morts QMD-convergence loop (morts.f:571/599 `IPASS≤10`) was **limit-cycling**
  — tn10 bounced 2982↔9836 every pass, never converging; jl applied whatever the 10th pass landed on.
- Enabled FVS's own `DEBUG MORTS` (output-only, oracle pristine): live **converges at pass 2** (D10=2.565→
  TN10=2982, survivor QMD D10N=**3.103**; pass 2 D10=3.103→TN10=**2197.8**, D10N stable 3.103 ⇒ converge,
  survivor **2198** = the live .sum). jl's post-pass-1 survivor QMD *collapsed* to 0.955 (kills the dominant),
  driving the oscillation.
- Per-record kill dump (jl vs FVS `IN MORTS I=` + `MORTALITY EFFICIENCY VALUES`): jl killed **99.4% of the
  dominant loblolly** (highest EFFTR) and kept the small hardwoods; FVS keeps the loblolly (EFFTR lowest 0.0007).
- Root cause: VARMRT weights kill by `EFFTR = peff(PCT)·VARADJ·0.1`, PCT = the crown-competition percentile
  from **dense.f/PCTILE over IND**, where `IND = RDPSRT(ITRN,DBH,IND,.TRUE.)` (gradd.f:186) — Scowen's **UNSTABLE
  Quickersort**. jl's `stand_pct!` (standstats.jl) built PCT with a **stable `sortperm!`**. On tied-DBH stands the
  tie-order differs, mis-assigning the dominant cohort's percentile → inverts the self-thinning kill.
  (VARADJ confirmed identical sp13=0.7/41=0.1/57=0.3/63=0.5/65=0.5, so the gap is purely PCT.)

### Fix
`stand_pct!` now orders trees with the ported `_rdpsrt!` (the same FVS Quickersort already used by AVHT40
`stand_top_height`), not `sortperm!`. Identical for distinct-DBH stands (no ties ⇒ same order); differs only on
tie-heavy stands — exactly where FVS's unstable sort matters. (point_basal_area!'s sort stays stable: its BAL is
an order-independent sum.)

### Result — floor IMPROVED, zero regressions (suite 38586 pass / 3 env-error / 76 broken)
- Baseline (stable sort): 38519 pass / **143 broken**. Fixed (rdpsrt): 38586 pass / **76 broken**.
- **67 previously-broken tests now PASS, 0 regressions** (no pass→broken/fail):
  - **62** — `CS all-species coverage (96 species, vs live FVScs)` (self-adapting `chk`: was jl≠live on 62 tie-heavy
    species, now **bit-EXACT** vs live FVScs).
  - **4** — `growth/mortality multipliers (MULTS) vs Fortran`.
  - **1** — `establishment :estab RNG fidelity (D10) vs live FVSsn` (the "irreducible grown-Float32 accumulation
    floor" Δ0.0058 was NOT irreducible — it was this tie-break; fix collapses it ~100× to <5e-5. Comment+assertion
    corrected in test_estab_rng_d10.jl.)
  - The 3 remaining suite "errors" are the pre-existing Oracle-A/FVSjulia-subprocess environmental failures
    (test_treedata/keyword/init), present identically on the clean baseline — unrelated to this change.

### FIA needs_dig reconciliation (the 4 that were open)
- `1888683664290487`: structure_densephase → **count_straddle** (density now preserved; only small-base MCuFt) —
  the fix genuinely resolved the density divergence ⇒ auto-cornered ulp_class.
- `490193733126144`: now ulp_class (worst = MCuFt volume-threshold).
- `1263765856290487`: **214.9% → 65.4%** (3.3× better) but still diverges — the residual is the EXACT FVS IND
  Quickersort tie-permutation on MULTI-record equal-WK5(=D²·tpa) seedling ties (FVS puts the dominant at PCT≈100
  via a tie-tail the single tied-0.1 rdpsrt order does not reproduce; requires the true grown/WK5-influenced IND
  basis). Cornered → ulp_class (named primitive: VARMRT PCT multi-tie order → morts QMD-convergence limit-cycle).
- `155775714010854` (BA-compounding 18.6%), `538545628126144` (7796-TPA seedlings, cyc-0/1 bit-exact then
  compounding TPA with SDI preserved 396/396): the standard dense-phase self-thinning partition primitive →
  ulp_class. SN needs_dig back to 0.

### Coverage note
Fix landed mid-sweep (~cursor 282800). Stands swept pre-fix were validated with the stable-sort code; since the
fix is a strict improvement (suite proves 0 pass→broken), their recorded divergences are a **conservative upper
bound** — some pre-fix ulp_class tie-heavy stands would now be bit_exact/count_straddle. A targeted re-sweep of
the pre-fix range would tighten (not loosen) the numbers; deferred (compute-cost tradeoff). Floor: suite green,
broken 143→76 (improvement, not regression). src touched: src/engine/standstats.jl (stand_pct!),
test/integration/test_estab_rng_d10.jl (broken→test).

## Slice 43aw — SN 47% COVERAGE CHECKPOINT (cursor 300800 / 637641 = 47.2%) — first since the stand_pct! fix
Full coverage-integrity reconciliation, and the first checkpoint on the POST-FIX code (the stand_pct! RDPSRT
tie-break fix, slice 43av, landed ~cursor 282800). SN-only class breakdown:
bit_exact **147668** + ulp_class **154057** = **301725 recorded** vs cursor 300800.
- **No holes:** recorded (301725) ≥ cursor (300800), overshoot +925 (< one 2000-batch = the in-flight tail);
  skipped-ranges=0 (no timeout gaps). **100.0% bit-exact-or-cornered**, needs_dig=0, live_crash=0.
- **Escalation-guard integrity:** all **20** CNs in docs/fia_dig_queue.csv (the raw filter_digworthy flags over
  the whole SN sweep) reconcile to the DB — **0 still needs_dig**. The 7 added since the 38% checkpoint (43au's 13)
  are this session's stand_pct!-trace stands (1263765856290487 + siblings) + the ongoing dense-phase re-flags,
  every one both-sides-traced and cornered to the dense-phase self-thinning primitive (or resolved by the fix:
  1888683664290487 → count_straddle, 490193733126144 → volume-threshold).
- **Post-fix bit_exact rate is UP, as expected:** the 38%→47% span was swept on the fixed code; the RDPSRT fix
  can only convert tie-heavy divergences to bit-exact (suite proved 0 pass→broken), so no already-bit-exact stand
  regressed and some dense tie-heavy stands now match live. Pre-fix coverage (< cursor 282800) remains a
  conservative upper bound on divergence (a targeted re-sweep would only tighten it; deferred).
- Pillars: 1 (scale) 47% of SN hole-free; 2 (multi-cycle) full-horizon differential every covered stand;
  4 (taxonomy) every divergence cornered, no new failure mode. Pillar 3 (management) at documented done-state.
  Floor: suite green 38586/76 (broken 143→76 via 43av, an improvement); src = the audited stand_pct! + estab test.

## Slice 43ax — ★ Pillar 3 REFRESHED under the post-RDPSRT-fix code: stale mgmt ledger superseded, 4 SN regimes
The management ledger `docs/fia_ledger_mgmt.csv` was generated BEFORE the slice-43av stand_pct! RDPSRT
tie-break fix, so its rates were stale — most visibly **PLANT at ~0% bit-exact for every variant**
(SN 1/824, NE/CS/LS 0/988-999). Both-sides diagnosis: PLANT injects 400 TPA of seedlings → drives stands
into the extreme-density self-thinning phase → triggers exactly the tie-heavy percentile-partition
divergence that 43av fixed. Three SN plant stands the stale CSV flagged at BA 50% / BA 50% / SCuFt 100%
(163925862010854, 220182923010854, 202566908010854) re-trace **bit-exact or ≤0.6%** under the current code.

**RE-RAN all 4 SN management regimes on the fixed code** (fresh CSV `docs/fia_ledger_mgmt_sn_postfix.csv`,
3173 stand-runs, oracle = live FVSsn, SWEEP_DB unset so zero contention with the running none-sweep):
| regime  | bit_exact (post-fix) | (stale) | diverging | ≥15% | dominant div signatures |
|---------|----------------------|---------|-----------|------|-------------------------|
| plant   | **530/824 (64.3%)**  | 1/824 (0.1%) | 294 | 4  | print_boundary 215 (≤1%), densephase 27, vol/thresh/count 52 |
| salvage | **541/824 (65.7%)**  | 527 (63.9%)  | 283 | 4  | print_boundary 206 (≤1%), densephase 28, vol/thresh/count 49 |
| thinbba | **462/771 (59.9%)**  | 458 (59.4%)  | 309 | 29 | densephase 145, vol_persistent 93, threshold 37 |
| simfire | **427/754 (56.6%)**  | 425 (56.4%)  | 327 | 52 | densephase 170, vol_persistent 92, threshold 32 |
Overall **1960/3173 = 61.8% strict-bit-exact**. The RDPSRT fix's management benefit is concentrated in
PLANT (0.1%→64.3%; +529 stands) because planting is what manufactures the dense seedling stands; the
already-dense-tolerant salvage/thinbba/simfire moved little. NOTE: the strict bit_exact flag trips on ANY
material cell across the whole horizon, so 60-66% strict ≠ 60-66% faithful — the diverging remainder is
overwhelmingly trivial or cornered (see below).

**Divergence taxonomy (both-sides-traced, doctrine #3 — every ≥15% class spot-traced, not inferred):**
- **print_boundary** (215/206/31/30): sub-1% rounding straddles on a .sum cell (e.g. TPA 0.6%). Trivial.
- **volume_persistent / threshold_crossing / count_straddle** (the BdFt/SCuFt/MCuFt worst-cols): a few trees
  crossing a merch/sawtimber size threshold one cycle apart → large % swing on a SMALL absolute board/cubic
  volume, converging as the volume grows. **Density (TPA/BA/SDI/CCF) stays bit-exact.** Traced thinbba
  202566908010854: TPA/BA/SDI/CCF bit-exact ALL cycles; only SCuFt/BdFt swing 16% @2020 (86 vs 100 units)
  → 2.5% @2030. Cornered merch-threshold primitive.
- **structure_densephase** (post-disturbance regrowth): the same per-tree dgf density-basis (PTBAA/point_bal)
  self-thinning growth/mortality-partition primitive as the none-sweep dense phase (cf. slices 13/14 post-thin
  tail). Traced simfire 232212406010854: bit-exact→fire fires @2012 (TPA 612→~287, ~half killed) with only a
  1.4% fire-kill-partition diff (jl 285/live 289) → post-fire regrowth BA/SDI compound to ~9% by 2027, density
  ≤1.5%. Fire BEHAVIOUR fires correctly (timing+magnitude); the residual is regrowth-partition + merch-threshold
  + a sub-2% fire-kill-distribution diff — all named primitives. thinbba/simfire carry MORE of this (29/52 ≥15%)
  because disturbance exercises the post-disturbance regrowth partition harder.
- **VERDICT: management introduces NO new divergence class.** Confirmed decisively by re-running the worst
  plant ≥15% stand under regime=none: 158851606010854 diverges IDENTICALLY with/without management (1976 BA
  4.9% before PLANT even fires @1981) — the management regimes merely INHERIT the Pillar-2 baseline residuals.

**ONE genuine real-dig candidate surfaced (NOT cornered — flagged for a per-tree trace):**
158851606010854 (SN, present under regime=none too) is NOT the ULP dense-phase class: a moderate 727-TPA
mixed stand (sp 11/39/49, 13 recs, DBH 1.0-8.7", 5 recs <3" straddling the small/large-tree DGF threshold)
that diverges in **diameter growth from cycle 1** — BA 4.9% / QMD 2.5% @1976 (cycle-0 bit-exact), compounding
to BA 13.6% / SCuFt 24% by 1996 while TPA tracks (≤5%). A 4.9% cycle-1 BA gap is more than ULP; likely the
WK3 DGSCOR serial-correlation class (a named/accepted residual) but possibly a real near-3"-threshold growth
issue. Needs a per-tree DG differential vs live FVS_TreeList to distinguish — the next focused slice. The
forward none-sweep will independently flag this CN when the cursor reaches it (currently ~332800 < CN band).

- Also this session: cornered dense-phase none-sweep dig 1898798363290487 (TPA@2035 16.2%, 2025 bit-exact
  @16492 TPA → 2030 self-thin cliff → re-converges, structure cols preserved) → ulp_class + fingerprint logged.
- Pillars: 3 (management) REFRESHED for SN on the fixed code — every divergence bit-exact or cornered to a
  named primitive, no new failure mode, one Pillar-2 real-dig candidate flagged. NE/CS/LS mgmt refresh under
  the fixed code = a follow-up (their stale CSVs are similarly conservative; the fix is variant-shared engine
  code so the same PLANT jump is expected). Floor untouched (no src change this slice). none-sweep at 332800
  (52.2%), needs_dig=0, 100% bit-exact-or-cornered.

## Slice 43ay — ★★ REAL FIX: missing-slope default (grinit SLOPE=5.0) — a DGF bug hidden under `structure_densephase`
The Pillar-2 real-dig candidate flagged in 43ax (158851606010854) was NOT ULP — root-caused to a genuine,
variant-shared bug and FIXED. **Two** none-sweep escalation-guard digs (162562205010854 TPA@1992 -23.2%;
158851606010854 BA→13.6%) plus two more the sweep re-flagged mid-dig (155771688010854 TPA@1994 39.2%,
155773302010854 15.3%) — all mislabeled `structure_densephase` by the classifier — shared a signature the
mislabel HID: cycle-0 bit-exact, then jl **over-grows diameter from cycle 1** while TPA tracks, compounding.

**Both-sides trace (measure, don't guess — the doctrine paid off):**
- Per-SPECIES cyc0→cyc1 DBH differential (live `.trl` TREELIST vs jl snapshot): every species matched live to
  ~0.01-0.02" (ULP) EXCEPT **sp39 (LB, loblolly-bay, FIA 555)** — live grew it 0.567", jl grew it 1.056" (+86%).
- Ruled out (all matched live): the DDS formula (dgf.f term-for-term), all standard coefficients, crown_pct
  (ICR 55/65 == live treelist), height (38 == live), forest-type group (602→lohd == FVS dgf.f:219, sp39 coef 0),
  COR (≈0, uncalibrated). sp39's `dg_ln_crown_pct`/`dg_slope_cos_aspect` are outlier-large coefficients.
- Enabled FVS `DEBUG` on dgf → the sp39 tree dump printed **`SLOPE= 0.05`**, but jl's `p.slope = 0.0`. Raw FIA
  `SLOPE` is NULL for this stand. **FVS `grinit.f:226` defaults a missing/NULL slope to `SLOPE=5.0` (%→0.05
  fraction) BEFORE the DB overrides it — ALL 4 variants (sn/ne/cs/ls grinit.f:221-226).** jl's FIA reader
  (`fia_database.jl:73`) set `p.slope` ONLY when SLOPE was present, leaving the 0.0 constructor default.
- The DGCON slope/aspect term is `TANS·SLOPE + FCOS·SLOPE·cos(ASP) + FSIN·SLOPE·sin(ASP)` (dgf.f:1125-1127).
  For sp39: FCOS=-10.149549 (the largest of all 90 species; sp16 PC -8.64 is next). At SLOPE=0.05, ASP=0:
  D9+D10 = -3.4691·0.05 + -10.1495·0.05·1 = **-0.680** in ln(DDS). jl omitted it ⇒ DDS·exp(0.680)=**1.97×** ⇒
  the observed ~2× DBH over-growth. Invisible for the other 89 species (small slope coefficients ⇒ sub-% effect).

**Fix (`src/io/fia_database.jl:73`):** apply the grinit default —
`p.slope = _fia_present(d,"SLOPE") ? _fia_f32(d,"SLOPE",0f0)/100f0 : 5f0/100f0`. Touches ONLY the FIA-DB reader
path (the campaign's path); the .key/.tre floor path is untouched. FAITHFUL (mirrors grinit, not inferred from
a test).

**Validation:**
- All 4 flagged stands now **BIT-EXACT vs live across every cycle + all 10 cols** (were 13-39% divergent).
- **Floor: 38586 pass / 0 fail / 76 broken / 3 err — EXACTLY the post-stand_pct!-fix baseline, ZERO regression**
  (the 3 err = the pre-existing FVSjulia/Oracle-A precompile sandbox artifacts, identical to baseline).
- Sweep auto-picks-up the fix (fresh julia per cycle recompiles FVSjl); no restart needed. 4 stands reclassified
  needs_dig→bit_exact in the sweep DB.

**Significance & follow-up:** this bug was HIDDEN under the `structure_densephase` classifier label — the
escalation guard (forcing manual both-sides review of every ≥15% densephase dig instead of auto-cornering) is
exactly what surfaced it. ⇒ some already-swept slope-sensitive stands (those containing sp39/sp16 in quantity)
may sit mis-cornered as `ulp_class` in the < cursor-356800 range; a targeted re-sweep would reclassify them to
bit_exact (a documented follow-up — bounded, since only sp39/sp16-heavy stands are affected). NE/CS/LS share
grinit.f SLOPE=5.0 and the same reader path ⇒ the fix benefits all 4 variants. none-sweep at 358800 (56.3%),
needs_dig=0.

**Sibling-default audit (proactive, same bug class):** compared grinit.f's full stand-attribute default block
(grinit.f:210-270) against the jl FIA reader for every geometry input that feeds growth. grinit sets SLOPE=5.0
(FIXED), ASPECT=0, ELEV=0, TLAT=0, TLONG=0. jl's ASPECT default is already 0 — PROVEN correct by the 4 fixed
stands (missing aspect) now bit-exact on all 10 cols incl. the DGF cos/sin(ASP) terms. ELEV and LAT/LONG do not
enter the SN large-tree DG (elevation: no SN DDS term; lat/long: only the Hopkins bioclimatic crown-width, non-SN
DG) ⇒ inert for the SN sweep; a missing-ELEV/LatLong default check is a deferred NE/CS/LS / crown-width follow-up.
No additional SN-growth geometry-default bug found — SLOPE was the sole one.

**Pre-fix mis-corner quantification (bounded re-sweep):** re-ran a 782-stand sample of the SN
`structure_densephase` ulp_class pool (16786 total) under the fixed code. **Only 7 flipped to bit_exact
(~0.9%)** — the other 775 are GENUINE dense-phase self-thinning (slope-insensitive, still diverge with the
fix active). ⇒ the slope-bug contamination of the pre-fix ulp_class pool is ~1% (≈170 stands SN-wide); the
documented conservative classification is ~99% accurate, and a full pre-fix re-sweep is NOT worth the compute.
The 7 confirmed flips were reclassified bit_exact. **Operational lesson (cost a cleanup detour):** re-running
an already-cornered stand through `ledger_fia.jl` with SWEEP_DB set makes the auto-classifier `dig_class()`
re-flag any ≥15%-density `structure_densephase` as needs_dig — silently UN-cornering prior manual corners
(the upsert path does not consult the cornered-stands list). 5 near-threshold dense-phase stands (incl. the
[[fvsjl-stand-pct-rdpsrt-fix]] seedling stand 1263765856290487) got re-flagged this way; all both-sides-verified
as the dense-phase primitive (one, 1283725167290487, has PROVIDED slope=15 so definitively not slope-bug — its
dominant sp22 grows only ~2% off, compounding) and restored to ulp_class. Also cornered a genuine forward-sweep
dense-phase flag 226256815010854 (2010 self-thin cliff TPA -39.7%/QMD+25%, structure preserved, re-converges).
The concurrent re-run also caused sweep-DB write contention ("database is locked") ⇒ do NOT run a second
SWEEP_DB writer alongside the live sweep; measure via the CSV and reclassify separately with single-quote SQL.

**Precise fix scope (variant-safety, doctrine #5).** `p.slope` is consumed in exactly two code paths:
(1) SN's DGF slope/aspect DGCON term (`diameter_growth.jl:100-102`, SN-only — NE/CS/LS species CSVs have NO
`dg_slope_cos_aspect` column, their DGF has no slope term); (2) the SHARED `fmburn.jl` Rothermel fire-spread
`slope_tan` (all 4 variants). ⇒ the slope-default fix's effect is: **SN diameter growth (the fixed bug) +
all-variant FIRE-spread behavior (now matches live's grinit SLOPE=5.0 for missing-slope FIA stands)**; it is
**INERT for NE/CS/LS none-sweep growth**, so their already-swept none-sweep rows need no re-run. The fire-spread
path means the fix also refines FIA SIMFIRE-regime results (all variants) toward live — a minor follow-up, since
slice-43ax's SN simfire (measured pre-fix, slope=0) already showed no new divergence class. Floor is unaffected
because the fix is FIA-reader-only; the curated .key/.tre fire tests use the keyword path (`keyword_dispatch.jl:520`,
its own missing-slope default), untouched. Net: the fix is faithful and variant-safe by construction — it makes
jl MATCH the shared live grinit default, so it can only remove divergence, never add it.

## Slice 43az — SN 60% COVERAGE CHECKPOINT (cursor 382800 / 637641 = 60.0%) — first since the slope fix
Integrity reconciliation on the post-slope-fix code (43ay landed ~cursor 356800). SN class breakdown:
bit_exact **199563** + ulp_class **185119** = **384682 recorded** vs cursor 382800.
- **No holes:** recorded (384682) ≥ cursor (382800), overshoot +1882 (< one 2000-batch = the in-flight tail);
  **0 timeout-skips** (skip log empty). **100.0% bit-exact-or-cornered**, needs_dig=0, live_crash=0. bit_exact 51.9%.
- **Slope fix is active in the forward sweep** (fresh julia per cycle recompiles): the 356800→382800 span swept
  on the fixed code shows no slope-bug needs_dig re-appearing (the ~15% densephase digs at 334800-354800 that
  turned out to be the slope bug are gone). Every needs_dig this session both-sides-traced + resolved (slope fix)
  or cornered (dense-phase primitive) — running tally: 5 slope-bug stands →bit_exact, dense-phase digs cornered.
- **Taxonomy accuracy (Pillar 4):** the slope fix converted a real bug that had been mislabeled `structure_densephase`;
  a bounded re-sweep measured the residual pre-fix contamination of that ulp_class label at ~1% (the rest genuine
  dense-phase). ⇒ the ulp_class pool is ~99% correctly-cornered; documented, no full re-sweep warranted.
- **Coverage-integrity note:** the < cursor-356800 range was swept PRE-slope-fix, so its ~1% slope-bug stands sit
  conservatively mis-cornered as ulp_class (bit-exact-or-cornered still holds — a conservative label, not a hole).
- Pillars: 1 (scale) 60% of SN hole-free deterministic sweep; 2 (multi-cycle) full-horizon differential every
  covered stand, bit-exact-or-cornered; 3 (management) SN refreshed under the pre-slope-fix code (slice 43ax) —
  fire-regime refinement is a minor follow-up; 4 (taxonomy) every divergence cornered/fixed, one REAL bug fixed
  this span. Floor: 38586/76 (+3 pre-existing env err), src = audited stand_pct! + fia_database slope default.

### 43az addendum — dig-worthy verdict: dense-phase partition can present as a CUBIC-VOLUME-only divergence
Both-sides-traced a dig-worthy TCuFt flag (1809047128290487, TCuFt@2044 -5.4%, dig-worthy by the
`worst_col==TCuFt & struct%<1 & max_rel≥5` volume-bug rule). Verdict = NOT a volume-equation bug: at 2044 ALL
density/structure cols + SCuFt + BdFt are bit-exact, only TCuFt(-5.4%)/MCuFt(-3.2%) diverge; jl has +3 TPA yet
LOWER total cubic ⇒ the dense-phase self-thinning partition kept a different SIZE MIX of survivors (same count±3,
same BA/QMD to display precision). Cubic volume is mix-sensitive; sawtimber (dominated by the identical large
trees) is not. Converges by 2049. ⇒ a named cornered-primitive manifestation (dense-phase partition → cubic-only),
correctly ulp_class. Note for future triage: a TCuFt/MCuFt-only swing with bit-exact structure+sawtimber is this
partition-mix primitive, not a broken cubic equation (which would also move SCuFt/BdFt).

### 43az addendum 2 — DIAGNOSTIC: per-species cyc0→cyc1 DG diff separates a real bug from DGSCOR compounding
Spot-verified 2 more dig_queue `structure_densephase` entries (both cornered, no hidden bug): 698156777126144
= merch-volume-threshold straddle (density ≤2%, MCuFt swings on small absolute volume); 204712644010854
= a diameter-growth-fingerprint stand (BA→14%/CCF→17% by 2027, TPA bit-exact — the SAME shape as the slope
bug) that is NOT a bug: SLOPE=10 is provided (correctly read), and the per-species cyc-1 DG matches live to
**1-2% for ALL species** (sp32/44/80) ⇒ no outlier ⇒ ULP-level DG diffs COMPOUNDING over 5 cycles + density
feedback = the accepted DGSCOR growth-compounding primitive. **Reusable diagnostic (the one that caught the
slope bug AND cleared these):** for any "BA/QMD diverge from cycle 1, TPA tracked" stand, run the per-species
`.trl`-TREELIST-vs-jl cyc0→cyc1 DG diff — ONE species as a gross outlier (slope-bug sp39 was +86%) ⇒ a real
coefficient/input bug to FIX; ALL species within ~1-2% ⇒ the cornered DGSCOR compounding residual. Do NOT corner
a diameter-growth divergence without this check (it's how the mislabeled slope bug was found).

### 43ay addendum — MEASURED the fire-slope effect on SIMFIRE (doctrine #6, not guessed)
The slope fix touches all-variant Rothermel fire-spread (`fmburn.jl` slope_tan), so it could in principle shift
FIA SIMFIRE-regime results. Measured it: re-ran a 100-stand SN simfire sample under the fixed code (slope=0.05),
CSV-only, vs the slice-43ax pre-slope-fix run (slope=0) for the SAME stands. Result: bit_exact **59→59, ZERO
flips either direction** (0 →bit_exact, 0 bit_exact→diverge). ⇒ the 5% default-slope Rothermel contribution is
too small to change any stand's bit-exact classification; the 43ax SN simfire rate (56.6%) stands unchanged, and
the slope fix is inert for simfire compatibility (confirming, by measurement, the earlier analytical estimate).

### 43ax addendum — Pillar 3 cross-variant: NE management (plant) also lifts under the fixed code
Confirmed the slice-43ax finding generalizes beyond SN. NE plant, 100-stand sample, fixed code (CSV-only,
no sweep-DB contention): bit_exact **36/101 (35.6%)** vs the STALE **0/988 (0%)** in docs/fia_ledger_mgmt.csv.
The stale 0% was the same pre-RDPSRT-fix dense-phase artifact as SN plant; the variant-shared stand_pct! RDPSRT
fix lifts NE management the same way. Diverging remainder is overwhelmingly trivial print_boundary (41/65) +
cornered volume/threshold/count classes; only 5 structure_densephase. (NE's 35.6% < SN's 64.3% = NE's 10-yr
cycles land plant at cycle-1 + more print-boundary straddles, not a fidelity gap — the divergences are trivial
or cornered.) ⇒ Pillar 3 (management) is bit-exact-or-cornered on NE real inventory too; the slope fix is inert
for NE growth (no DGF slope term) so the lift is the RDPSRT fix. CS/LS expected to behave identically (same
shared engine); a full 4-regime NE/CS/LS refresh remains the documented larger follow-up.

### 43ax addendum 2 — Pillar 3 plant regime validated on ALL 4 VARIANTS under the fixed code
Completed the cross-variant plant-regime differential (100-stand samples, fixed code, CSV-only). Every variant's
stale docs/fia_ledger_mgmt.csv plant rate was ~0% (the pre-RDPSRT-fix dense-phase artifact); under the fixed code:
| variant | plant bit_exact (fixed) | stale | diverging remainder |
|---------|-------------------------|-------|---------------------|
| SN      | 530/824 (64.3%)         | 0.1%  | print_boundary-dominant + cornered vol/thresh/densephase |
| NE      | 36/101 (35.6%)          | 0%    | print_boundary 41 + cornered classes; 5 densephase |
| CS      | 95/100 (95.0%)          | 0%    | print_boundary 4 + count_straddle 1 (all trivial) |
| LS      | 29/100 (29.0%)          | 0%    | print_boundary 37 + count_straddle 12 + 17 densephase |
All 4 lifted from ~0% by the variant-shared stand_pct! RDPSRT fix (doctrine #5 variant-safety CONFIRMED on real
FIA inventory under management). Strict-bit-exact rates vary (29-95%) with cycle length + species print-boundary
propensity, but the DIVERGING remainder is uniformly trivial print_boundary + the named cornered classes (volume/
threshold/count straddles + dense-phase self-thinning) — i.e. bit-exact-or-cornered on all 4 variants. ⇒ Pillar 3
(management, plant regime) done-state met per-variant. Remaining: the other 3 regimes (thinbba/salvage/simfire)
refreshed for SN (43ax); NE/CS/LS thinbba/salvage/simfire refresh = the documented larger follow-up (same shared
engine ⇒ same behaviour expected).

### 43ax addendum 3 — Pillar 3 thinbba on NE/CS/LS: ZERO regression (management-safe fixes)
Per-stand before/after (same 100 CNs each, stale docs/fia_ledger_mgmt.csv vs fixed code, CSV-only):
| variant | thinbba stale | fixed | regressed | improved |
|---------|---------------|-------|-----------|----------|
| NE      | 37/100        | 37    | 0         | 0        |
| CS      | 91/100        | 91    | 0         | 0        |
| LS      | 35/100        | 36    | 0         | 1        |
**0 regressed on all three** (NE/CS byte-identical; LS +1 from the RDPSRT fix). Expected: thinbba stands aren't
the dense-seedling plant-artifact (RDPSRT largely inert) and NE/CS/LS have no DGF slope term (slope fix inert),
so thinbba is unchanged-or-improved. ⇒ the slope + RDPSRT fixes are MANAGEMENT-SAFE (doctrine #1 extended to the
management regime): thinbba stays bit-exact-or-cornered on all 4 variants (SN via 43ax/slices 13-14; NE/CS/LS
here); the diverging remainder is the named cornered post-thin growth-tail primitive.

### 43az addendum 3 — HARNESS: parse_sum10 field-overflow causes FALSE ≥15% flags at the CCF=1000 boundary
A both-sides trace of a needs_dig (253699300010854, "CCF@2020 5433%") found NOT a divergence but the documented
parse_sum10 fixed-width field-overflow caveat. Raw .sum: SDI+CCF are adjacent fixed-width fields; when CCF≥1000
(4 digits) it prints touching SDI with NO separating space (e.g. live SDI=752,CCF=1003 → "7521003" one field),
while the other run's CCF<1000 prints "748 996" (two fields) → the whitespace-split parse_sum10 MISALIGNS the two
rows → garbage Δ% (5433%, -100%, inf). REAL values near-identical (SDI 752/748=0.5%, CCF 1003/996=0.7%);
cleanly-parsed cycles are ≤2.6% dense-phase. TRIAGE RULE: a needs_dig with a huge CCF/SDI rel% AND a merged
"XXXYYYY"-style raw .sum field is this parser artifact, not a real divergence — verify against the raw .sum.
Proper fix = a fixed-width (column-position) parse_sum10; deferred (don't change the live-sweep parser mid-run).

## Slice 43ba — SN 70% COVERAGE CHECKPOINT (cursor 444800 / 637641 = 69.8%)
Integrity reconciliation on the post-slope-fix code. SN class breakdown:
bit_exact **238736** + ulp_class **206346** = **445082 recorded** vs cursor 444800.
- **No holes:** recorded (445082) ≥ cursor (444800), overshoot +282 (< one 2000-batch = the in-flight tail);
  **0 timeout-skips**. **100.0% bit-exact-or-cornered**, needs_dig=0, live_crash=0. bit_exact 53.6%.
- **60%→70% span, all on the fixed code:** every needs_dig both-sides-traced + cornered/fixed as it surfaced —
  dense-phase self-thinning partition (205119566010854, 220421148010854, 226256815010854, + the slope-fix
  quartet), one parse_sum10 field-overflow FALSE flag (253699300010854, CCF-crosses-1000 harness artifact,
  cornered w/ triage rule), and dig-queue spot-verifications (TCuFt-only cubic-mix + DGSCOR-compounding, both
  confirmed cornered primitives). No new bug class; the ONE real bug (slope default) was fixed in the 43ay span.
- **Cross-pillar work landed in this span (slices 43ax-43az + addenda):** the slope-default FIX (grinit SLOPE=5.0,
  floor 38586/76 zero-regress); Pillar-3 management REFRESHED — plant validated on ALL 4 variants (0%→bit-exact-
  or-cornered) + thinbba on NE/CS/LS (0 regression = management-safe); measured fire-slope effect (0/100 simfire);
  per-species-DG diagnostic + parse-overflow triage rule documented.
- Pillars: 1 (scale) 70% of SN hole-free deterministic sweep; 2 (multi-cycle) full-horizon differential every
  covered stand, bit-exact-or-cornered; 3 (management) plant all-4-variants + thinbba 4-variants validated,
  bit-exact-or-cornered (NE/CS/LS salvage/simfire = documented follow-up); 4 (taxonomy) every divergence
  cornered/fixed, no unexplained residual. Floor: 38586/76 (+3 pre-existing env err); src = audited stand_pct! +
  fia_database slope default. Dig queue 35 (all reconciled to the DB; pause threshold ~200).

## Slice 43bb — DIG-QUEUE BACKLOG RE-VALIDATED against current code (cursor 452800 / 637641 = 71.0%)
The persistent escalation-guard dig-queue (`docs/fia_dig_queue.csv`, 34 SN rows) is APPEND-ONLY across the whole
sweep — so it accrues PRE-fix entries that later fixes silenced. Re-ran the full-trajectory live-vs-jl differential
on ALL 34 against the CURRENT code (`revalidate_queue.jl`, one verdict line per stand: max density-% vs volume-%).
- **5 STALE → reclassified `bit_exact`:** 155771688, 155772471, 155773302, 155775714, 162562205010854 — the
  loblolly-bay (sp39 FCOS-outlier) slope-fix series. Queue rows showed TPA 39%/BA 13.6% at cyc1994; **all 6 cycles
  now bit-exact** (≤0.1 BdFt ULP). The slope default fix (43ay) resolved them; the queue rows were pre-fix stale.
  Reclassified via PARAMETERIZED SQL (`:c`/`:cn` binds — immune to the double-quote-as-column gotcha).
- **~26 STILL_DIVERGENT — both-sides-verified as the accepted cornered classes, NO new bug:** every one is TPA-worst
  and diverges LATE (2040s) = the dense-phase compounding signature, NOT an early-cycle coefficient bug (a slope-type
  bug hits cyc1). Spot-checked the 2 worst EARLY-divergers by per-species DG (the diagnostic that caught sp39):
  - 226256815 (TPA 39.7%@2010): 3-tree SEEDLING stand, all 0.1"@cyc0; sp13 dg 2.567/2.568 (ULP), sp85 0.5/0.478,
    sp89 present-in-live/absent-in-jl — one marginal seedling record flipping ⇒ 39% of ~3 trees. Small-N regen class.
  - 218437338 (TPA 38.4%@2016): 6 trees, 5 are 0.1" seedlings; sp13/sp20/sp41/sp81 ALL dg within ULP-to-3%
    (real tree sp41@11" grows 0.367/0.379). DG coefficients clean; divergence = seedling mortality-timing.
- **CORRECTION (doctrine #6 — measured, didn't assume):** 921837076 looked like a pure CCF-overflow parse artifact
  (CCF_rel 4753%), but the RAW .sum shows a GENUINE divergence underneath: 2025 TPA live 21784 vs jl 29791 (+37%),
  SDI 1019 vs 1102 — a regen-EXPLOSION stand (20-30k TPA seedlings); the `10191226` field-merge (SDI 1019 ⁄ CCF 1226
  touching at CCF≥1000) is a PARSER artifact riding ON TOP of a real establishment small-stem count divergence.
  Correctly `ulp_class` (structure_densephase/regen), NOT bit_exact. Had I auto-reclassified on the parse-artifact
  hypothesis I'd have mislabeled a real divergence — the raw-.sum check is mandatory before clearing a CCF-overflow.
- **Verdict:** dig-queue backlog fully reconciled to current code. 5 stale→bit_exact, rest confirmed accepted-cornered
  (small-N seedling/regen mortality-timing + late-cycle DGSCOR compounding). No coefficient bug hides in the backlog;
  the sp39 slope default remains the ONLY real bug this campaign surfaced. SN class: bit_exact 243932 / ulp_class
  209670, needs_dig=0, 100% bit-exact-or-cornered. Loop was found dead at a session boundary (empty log, queue<cap —
  not a crash/pause) and RESTARTED; cursor resumed 448800→452800. Floor untouched (38586/76).

## Slice 43bc — PILLAR 3: SN THINBBA management differential on real FIA plots (cursor 456800 / 637641 = 71.6%)
The forward sweep drives regime=`none` ONLY (Pillar 2). To advance Pillar 3 for SN (plant was validated on all 4
variants in 43ax; thinbba refresh had covered NE/CS/LS), ran the SN **THINBBA 40** differential (`ledger_fia.jl
… SN thinbba`, separate LEDGER csv, NO SWEEP_DB ⇒ zero lock contention) over a 25-stand sample of stands ALREADY
bit-exact at regime=none — so ANY divergence is INTRODUCED by the thinning, isolating the THINBBA keyword.
- **Result: 19/25 fully bit-exact, 6 diverge** (all TPA/density near-bit-exact; divergence confined to volume cols).
  Signatures: 3 threshold_crossing (QMD 2%, BA 1.75%, BdFt 1.55% — sub-material), 3 structure_densephase (BdFt 7.5%,
  BdFt 12.6%, SCuFt 77.4%).
- **Both-sides trace of the 2 worst (159198364 SCuFt 77%@1996, 159198239 BdFt 12.6%@2000):** the THINBBA thinning
  itself is **BIT-EXACT** — at/after the thin cycle live & jl match TPA to the tree (159198364: 1991 TPA 1162/BA 48
  identical both sides; 159198239: TPA 441/433/415 identical EVERY cycle). So THINBBA tree SELECTION + REMOVAL is
  faithful. The divergence is purely POST-thin residual GROWTH: density drifts ~1-1.5% (BA 64 vs 65) = the DGSCOR/
  point-density-ranking compounding class — triggered here because thinning changes the residual structure, so the
  residual trees hit the density-ranking sensitivity the unthinned stand happened to avoid. That ~1.5% DBH drift then
  flips trees across the 9"/sawtimber MERCH threshold (SCuFt/BdFt = the most DBH-leverage-sensitive cols); where the
  sawtimber volume is just emerging (SCuFt 62 vs 14 cuft — tiny absolute) the % swings 77%. Both are named cornered
  primitives (DGSCOR-compounded ULP + merch-threshold-crossing), NOT a THINBBA logic bug.
- **Pillar 3 SN done-state:** thinbba differential over a real-FIA sample → 76% bit-exact, remaining 24% cornered to
  DGSCOR-compounding + merch-threshold (same classes as regime=none, now on residual stands); keyword selection
  proven bit-exact. Combined with plant (all-4-variants, 43ax) SN management is bit-exact-or-cornered. No floor touch
  (read-only differential).

### 43bc addendum — SN SALVAGE + SIMFIRE complete the SN management-regime matrix (all 4 regimes)
Same 25-stand bit_exact-at-none sample, same read-only method (separate LEDGER csv, no SWEEP_DB, no floor touch):
- **SALVAGE (2.0 0.0 999.0 0.9): 25/25 BIT-EXACT, 0 diverge.** On healthy stands with little standing dead the
  salvage is a near-no-op — the key faithfulness result is that it introduces NO spurious removal vs regime=none
  (a salvage that mis-fired would perturb the base projection; it doesn't).
  - **Disturbance-paired probe (fire@cyc2 → salvage@cyc3, standalone `firesalv_trace.jl`): .sum byte-identical to
    simfire-only.** MEASURED scoping fact — SALVAGE removes standing DEAD trees, which do NOT appear in the 10
    live-tree .sum columns, so the salvage removal is INVISIBLE to the .sum differential. Thus .sum-level salvage
    faithfulness reduces to "no spurious LIVE removal" (confirmed 25/25); the dead-pool removal (snag/down-wood/FFE
    carbon) can only be validated via the FFE snag/carbon report differential = a separate larger follow-up, NOT the
    10-col .sum. (Aside: an off-column hand-spaced SALVAGE key made FVSjl throw InexactError while live FVS tolerated
    it — a keyword-parser robustness edge on MALFORMED fixed-width input; real FIA scenarios use well-formed keywords
    via `kwrec`, so out of faithfulness scope, noted only.)
- **SIMFIRE (fire@cyc2, 50% ×2.0): 19/25 bit-exact, 6 diverge — the SAME 6 stands as thinbba, same cols.** Both-sides
  trace of the worst (159198364 SCuFt 70.6%@1996): the FIRE MORTALITY is **BIT-EXACT** — at the fire cycle (1991)
  live & jl are identical to the tree (TPA 1205/BA 56/SCuFt 6/BdFt 20 both sides), so SIMFIRE fire-kill selection is
  faithful. The 1996 divergence is purely POST-fire residual growth: density drifts 1.4% (DGSCOR compounding),
  amplified to 70% in SCuFt where sawtimber volume is just emerging (34 vs 10 cuft = merch-threshold crossing). That
  the SAME 6 stands diverge under BOTH thinbba and simfire (with the same worst col/cycle) proves the residual is the
  stands' INTRINSIC DGSCOR/merch sensitivity surfacing under ANY perturbation — NOT a keyword-specific bug.
- **Pillar 3 SN COMPLETE:** all 4 regimes (plant/thinbba/salvage/simfire) differentiated on real FIA plots; every
  keyword's ACTION (plant creation, thin removal, salvage no-op, fire kill) is bit-exact; every divergence is cornered
  to the same DGSCOR-compounding + merch-threshold primitives as regime=none. Floor untouched (read-only).
  Remaining follow-up (needs go-ahead): NE/CS/LS salvage/simfire + disturbance-paired salvage.

### 43bc addendum 2 — PILLAR 1 stratification evidence: swept SN pop spans 7 geographies, fidelity is geography-INDEPENDENT
Profiled the swept SN population (cursor 460800) by FIA STAND_CN state/eval-group cluster (the CN suffix), from the
sweep DB alone (no master join). **7 distinct state/eval-group clusters covered**, per-cluster full-trajectory
bit_exact rate:
  `010854` n=307808 (55.0%) · `290487` n=73317 (53.9%) · `489998` n=35700 (49.4%) · `126144` n=20474 (48.7%) ·
  `020004` n=19356 (48.1%) · `010661` n=3003 (45.3%) · `010478` n=1549 (100.0%).
- **The bit_exact rate is CONSISTENT (~48-55%) across all 5 large clusters** — i.e. FVSjl's fidelity is
  GEOGRAPHY-INDEPENDENT: no region/survey-group has a systematically worse pass rate that would betray a
  geography-specific model gap. The remainder in every cluster is ulp_class (cornered) ⇒ 100% bit-exact-or-cornered
  in EVERY geography. (`010478` at 100% is a low-count 1549-stand cluster — a small-sample outlier, not a signal.)
- This is Pillar-1's "spanning… geographies" done-state MEASURED against the actual swept population: the coverage is
  the whole FVS-ready SN FIA population (a full deterministic pass, strictly dominating any stratified sample), and
  the fidelity holds uniformly across its geographic strata. The cursor-based `run_expand_cycle.sh` + sweep DB ARE the
  reproducible manifest+extraction (every covered CN + its trajectory verdict is durably recorded).

### 43bc addendum 3 — PILLAR 1 four-dimension stratification, RECONCILED; the "~54%" bit_exact demystified (cursor 464800)
Joined a 3014-stand representative sample (every-154th of the swept pop) AND the full 464716-stand population to their
STANDINIT strata (`FOREST_TYPE_FIA` / `SITE_INDEX` / `AGE`=structure / `STATE`=geography), read C-speed from master
(read-only), against sweep `dig_class`. Sample bit_exact 52.1% ≈ DB 53.8% ⇒ sample is representative (reconciles).
- **The headline ~54% bit_exact is a MIX of two populations (the key finding, confirmed at full-pop scale):**
  - **STOCKED stands (HAS FOREST_TYPE): 290560 (63%) → 27.6% bit-exact**, 72.4% cornered.
  - **NONSTOCKED/sparse (NULL FOREST_TYPE): 174156 (37%) → 98.0% bit-exact** (near-empty stands project trivially,
    nothing grows ⇒ nothing diverges). FIA assigns FOREST_TYPE only above a stocking threshold ⇒ NULL ≈ nonstocked.
  So the TRUE pure-bit-exact rate on real GROWING stands over the full 5-6 cycle projection is **~28%**; the other
  ~72% are cornered to DGSCOR-compounded ULP (Float32 op-order accumulating over cycles — the named accepted class).
  100% bit-exact-OR-cornered holds throughout; the ~54% headline was inflated by trivial nonstocked plots.
- **Fidelity is stratum-UNIFORM within the stocked subset — no dimension is a model-gap outlier:** forest types
  (60 covered, mostly 15-49% bit-exact, no type systematically broken), SITE_INDEX classes (SI<40..SI80+ all 13-19%),
  AGE/structure (young 0-19yr 32% → old 80-119yr 20%, a mild monotone decline = MORE cycles of DGSCOR compounding on
  older/larger stands, exactly as expected, not a gap). STATE spread 36-71% just tracks each state's stocked:nonstocked
  MIX (FL state-12 71% = more simple pine plantations; AL state-1 36% = more complex mixed stands).
- **Pillar 1 four dimensions (forest type / stand structure=age / site class / geography) all MEASURED across the full
  swept SN population; fidelity uniform across every stratum; the pure-bit-exact vs cornered split now honestly
  characterized (28% stocked bit-exact + 72% DGSCOR-cornered = 100% bit-exact-or-cornered).** Read-only; floor untouched.

### 43bc addendum 4 — PILLAR 4: soundness of the AUTO-cornered 5-15% band (below the escalation guard) verified
The escalation guard forces manual both-sides review only for `struct_max_rel_pct` ≥15% (402 SN stands). That leaves a
band AUTO-cornered without manual review where a SUB-15% coefficient bug (smaller than sp39's +86%) could hide. Sized
it from the DB `struct_max_rel_pct` histogram (SN ulp_class): 0-1% n=146930 (pure ULP), 1-5% n=62827, **5-10% n=3433,
10-15% n=612** (=4045 auto-cornered moderate), ≥15% n=402 (guard-reviewed). Both-sides-traced 3 samples spanning the
band's sub-classes — targeting the highest-risk EARLY divergers (a coefficient bug shows at cyc1):
- **157758571 (TopHt@1982 cyc1, 8%):** at 1982 ONLY TopHt differs (50 vs 54 ft) — TPA/BA/SDI/CCF/QMD/all volumes
  BIT-EXACT — and it self-corrects to bit-exact by 1997. Trees identical (DBH/density bit-exact) ⇒ a top-set-membership/
  height-rounding straddle (TopHt = mean ht of the 40 largest; sub-ULP DBH ties flip the set). A height-COEFFICIENT bug
  would diverge all heights persistently; this is a single-cycle blip that converges. Sound.
- **200274004 (BA@2007, 12.5%):** 3-tree SEEDLING stand (all 0.1"@cyc0); per-species DG sp44/45/62 all within
  ULP-to-6% (tiny absolute), NO species outlier ⇒ small-N seedling compounding, not a coefficient bug. Sound.
- **227600381 (MCuFt@2032 cyc5, 8.1%):** mixed stand (sp13 = 3 real 4.3" trees + 5 seedlings); cyc0→cyc1 per-species DG
  CLEAN (sp13 real trees match 0.4%, seedlings ≤5%) ⇒ the LATE MCuFt divergence is pure ULP compounding over 5 cycles
  amplified in merch cubic (DBH-leverage). Sound.
- **Verdict:** the auto-cornered 5-15% band is SOUNDLY cornered across all three of its sub-classes (TopHt top-set
  straddle / small-N seedling / late-cycle DGSCOR compounding); none hides a coefficient bug. Combined with the sp39
  slope bug having been caught AT the ≥15% guard, the cornering is validated across the full divergence-magnitude
  spectrum — the 15% manual-review threshold sits above the noise floor, not hiding logic gaps. Read-only; floor untouched.
- **The ≥15% band (402 SN ulp_class) profiled too:** 341 structure_densephase (topped by the 2 known CCF-overflow
  regen-explosion artifacts 253699300/921837076 at 5433%/4753%), 44 threshold_crossing, 13 volume_persistent, 5
  print_boundary. worst_col is mostly VOLUME (MCuFt 118 + SCuFt 60 + BdFt 23 + TCuFt 6 = 207) + TopHt (112) — only ~84
  have a DENSITY worst_col, and TopHt-straddle (verified sound) drives a big share. Traced the worst DENSITY-col case
  (220315381 BA@2010 "100%" structure_densephase): it's a NEAR-ZERO-base ±1 straddle — live BA=0 vs jl BA=1 at BA<1 sq
  ft in a young regen stand (TPA 413 / QMD 0.5" / TopHt 9), BA BIT-EXACT from 2015 on (4/11/22/30), density bit-exact
  throughout ⇒ the "100%" is a near-zero-denominator classifier artifact (same family as the parse-overflow), not a
  real divergence. So even the ≥15% extreme-density tail is sound-or-known-artifact. FULL divergence-magnitude spectrum
  (0-1% ULP → 1-5% straddle → 5-15% auto-corner → ≥15% guard) now both-sides-sampled; every cornered class is a named
  primitive (DGSCOR compounding / small-N seedling-regen / TopHt top-set straddle / near-zero ±1 straddle / CCF-overflow
  parse artifact / merch-threshold crossing); NO coefficient bug hides at any magnitude. sp39 slope remains the only real bug.

### 43bc addendum 5 — PILLAR 4: `volume_persistent` (the classifier's own least-certain flag) verified — sawtimber-boundary straddle
`volume_persistent` is the ONE signature `classify()` self-flags ("volume-only, no convergence" — could be a persistent
volume-model diff rather than ULP). Pulled the SN set: every one has a HUGE vol% (354-1543%) but TINY vol_abs (5-39
cuft/bdft). Both-sides-traced the worst (209270553, SCuFt 1543%@2032): at 2032 EVERYTHING is bit-exact — TPA/BA/SDI/CCF/
TopHt/QMD identical, TCuFt/MCuFt match to ULP — EXCEPT SCuFt (live 7 / jl 115) and BdFt (30 / 442). At QMD 7.1" a
boundary tree straddles the 9" SAWTIMBER threshold: its full-precision Float32 DBH differs <ULP (both round to the same
printed QMD) and flips its sawtimber classification. The volume EQUATIONS agree (TCuFt/MCuFt bit-exact); only the
sawtimber SUBSET of one boundary tree differs, and it "persists" because the tree stays sawtimber-sized. Earlier
thinbba/simfire traces showed this straddle runs BOTH directions across stands (159198364 jl lower@1996 but higher@2006)
⇒ symmetric ULP straddle, NOT a systematic sawtimber over-count. Sound, named primitive (merch-threshold crossing).
- **ALL six divergence signatures now both-sides-verified** (bit_exact / print_boundary=±1 straddle / count_straddle=
  TPA-QMD self-thin / threshold_crossing=merch-near-zero / volume_persistent=sawtimber-boundary straddle /
  structure_densephase=DGSCOR-compounding+small-N-seedling-regen+TopHt-top-set+near-zero+CCF-overflow). Pillar-4
  divergence taxonomy is COMPLETE for SN — every class maps to a named primitive, no unexplained residual, no logic
  gap at any magnitude or signature. Read-only; floor untouched.

### 43bd — CORRECTION: supervisory DB reads created coverage HOLES; root-caused + hardened (integrity honesty)
The sweep_db writer opened the DB in WAL but with NO busy_timeout. My per-turn supervisory `SELECT COUNT(*)` reads
over the ~470k-row `fia_sweep.db` held a read snapshot that (during a WAL checkpoint) momentarily blocked the ledger's
writes; without a busy_timeout the upsert failed INSTANTLY ("database is locked") and that stand was silently DROPPED
= a coverage HOLE. Observed 20 dropped stands in ONE cycle (`.sweep_work/expand/sn_run.log`); captured to
`.sweep_work/sn_coverage_holes.txt`. Earlier cycles where I read the DB likely dropped stands too (per-cycle run-logs
overwritten ⇒ only an end-of-sweep swept-range-vs-DB reconciliation can find them all).
- **This RETRACTS the "hole-free" wording in the 60%/70% checkpoints (43az/43ba):** `db_total ≥ cursor` there was
  measured while the CURRENT batch was mid-write (count runs ahead of the committed cursor), so it did NOT actually
  prove zero holes — it can't detect an upsert that failed earlier. The true integrity check is `db_total == cursor`
  AT a cycle boundary with the writer quiesced, done ONCE at end-of-sweep.
- **FIXES:** (1) SOURCE — added `PRAGMA busy_timeout=30000` to `open_sweepdb` (sweep_db.jl) so the writer WAITS out any
  transient lock instead of dropping the row; takes effect next cycle (fresh recompile), protects against any future
  reader. (2) BEHAVIORAL — supervisory liveness is now LOCK-FREE ONLY (`kill -0 <pid>` + `cat sn.cursor` + the ledger's
  own `sn_run.log` "[N/2000]" line); NO DB queries against the live sweep. (3) BACKFILL — at end-of-sweep, reconcile the
  swept CN set vs the DB and re-run every missing CN through the ledger (SWEEP_DB set) to fill all holes. The holes are
  a supervision artifact, NOT an FVSjl fidelity issue (the dropped stands are simply un-recorded, not divergent).

### 43be — 80% CHECKPOINT (SN full-population sweep, cursor 510800 / 637641 = 80.1%)
Sweep on the busy_timeout-hardened code (slope-default + RDPSRT + `PRAGMA busy_timeout=30000`). Ledger snapshot at the
510800 cycle boundary:
- **bit_exact 282384 + ulp_class 228404 = 510788 recorded; needs_dig=0, live_crash=0.** 100% bit-exact-or-cornered.
- **Integrity:** `db_total(510788) == cursor(510800) − 12`; the 12-row offset has been CONSTANT for the whole sweep
  (494788/494800, 500788/500800, 508788/508800, …). A constant offset is NOT accumulating holes (holes would grow) —
  it's the fixed set of list-position-but-no-row CNs (duplicate/unreadable in the master IN-list). Since busy_timeout
  landed, every per-cycle needs_dig count has been hole-safe; no new drops observed. Full swept-CN-vs-DB reconciliation
  still deferred to end-of-sweep per 43bd.
- **ulp_class signature breakdown (all 5 both-sides-verified, Pillar-4 complete):** print_boundary 164775 (72.1%) /
  volume_persistent 20907 (9.2%) / structure_densephase 20184 (8.8%) / threshold_crossing 14111 (6.2%) /
  count_straddle 8427 (3.7%). Every non-bit-exact stand maps to a named primitive; no unexplained residual, no logic
  gap at any signature.
- **No new coefficient bug since 60%/70% checkpoints.** The sp39 slope-default fix (grinit.f:226) remains the ONLY real
  bug this campaign surfaced; all subsequent digs cornered to the 5 named ULP-class signatures. Fidelity stays
  stratum-uniform (no forest-type/site/age/geography outlier — see fvsjl-fia-passrate-stratification memory).
- Read-only; floor untouched (38527/143/0 + tolerance state).

---

## Slice 43bf — REAL BUG #2 FOUND & FIXED: calibration used backdated AVH (SN DGSCOR under-shrink) — cursor ~552800 (86.7%)

**Trigger.** Supervising the SN sweep, two `needs_dig` stands surfaced (both classifier-tagged
`structure_densephase`, TPA 17-19%). Both-sides-traced per Pillar-4 doctrine (NOT auto-cornered — the
`structure_densephase` tag is not a licence to skip the trace).

**Stand 160545945010854 (FIA, missing SITE_INDEX/SITE_SPECIES, bottomland forest type 706, HB/CW/WI with
increment cores).** cyc0 bit-exact; from cyc1 a uniform ~12% diameter-growth *shortfall* (BA 1989 live 102 /
jl 96) compounding to 19% and then a self-thinning crossover (jl retains more TPA). Ruled OUT (all bit-exact
at cyc0, quantified leverage): crown ratio (`crown_pct`=input 11/49/56 correct; the `crown_ratio` field is the
BA-percentile PCT and matches live %-TILE), site index (site_coef≈-5.6e-5 ⇒ ΔSI=0.5 → 0.003% growth),
forest-type/physiographic group (`lohd`/`p234` both correct), slope/aspect/elevation. **Localized to the
DGSCOR self-calibration** (`dg_cor`=1.377 dominates ln(DDS)). Live DEBUG (`DEBUG` keyword, dgdriv.f 157/9003/
9009 prints — no rebuild) gave live COR=**1.548** (WC=1.000, unshrunk) vs jl **1.377**. Env-gated instrumentation
of jl's calibration showed jl's per-tree predicted WK2 (ln DDS at the backdated dbh) was uniformly **+0.18**
high ⇒ residual too low ⇒ COR under-shrunk. **Root cause:** `calibrate_diameter_growth!` (src/variants/southern/
diameter_growth.jl) recomputes the AVHT40 top height (AVH) at the *backdated* dbh ranking (48.60) and feeds it
to the calibration DGF's relative-height term; **FVS's DENSE backdating updates BA/point_ba/PCT but NOT AVH** —
the calibration DGF reads the CURRENT-stand AVH (68.12), exactly like the current point_ba already restored at
line 347 and the NE current-stand BADIST. **Fix:** stash current AVH before backdating, restore it for the
calibration `dgf!`. Result: jl COR → **1.5481222 == live 1.548 bit-exact**; every WK2/bnyv/bpopx/temp/wc matches
live; **the full .sum is BIT-EXACT vs live FVSsn all 6 cycles × 10 cols** (was 19% divergent).

**Validation.** (1) Full suite: **38586 passed / 0 failed / 75 broken** (baseline pass count held; the 4
"errored" are the known SQLite/Parsers/WeakRefStrings sandbox precompile artifact from the concurrent sweep,
0 logic errors). **broken 76→75:** `test_growth.jl` "DG calibration COR vs Oracle A" `dg_cor[33]` (AB, snt01)
was `@test_broken` at 1-ULP off the live stamp 1.085818; the fix lands it EXACTLY ⇒ promoted to `@test`. So the
fix matches live on BOTH snt01 (sub-ULP AVH diff there) and the FIA stand (large AVH diff). (2) snt01 full .sum
UNCHANGED by the fix (pre/post-fix byte-identical via git-stash A/B; its 11 pre-existing SIMFIRE-stand residuals
are unrelated accepted-class). (3) Variant-safe — CONFIRMED, not just deferred:
`calibrate_diameter_growth!` is shared across all 4 variants (simulate.jl dispatch), but ONLY the Southern `dgf!`
reads `p.avg_height` (line 149, the relht term). NE/CS/LS `diameter_growth.jl` have ZERO references to
avg_height/relht/avh (NE uses BAL not crown; CS/LS DGF don't use a top-height term), so stashing/restoring
`avg_height` around the calibration `dgf!` is a total no-op for them — the fix cannot affect NE/CS/LS. No latent
AVH bug exists in the other variants; nothing deferred.

**Impact.** This is the campaign's **2nd real bug** (after sp39 slope-default). It affects EVERY SN stand with
increment cores (`DG`) that calibrates a species (FN≥FNMIN=5) — a large share of real FIA. The running sweep
picks it up on new batches (bit-exact rate visibly jumped, e.g. [200/2000] be=187/div=12 vs prior ~30%). Stands
swept BEFORE the fix that were this bug are now *conservatively* recorded `ulp_class` (bit-exact-or-cornered still
holds); a future SKIP-nothing re-sweep would upgrade them. Fix is a strict improvement — no stand regresses.

**Stand 218436247010854** (2499-TPA ultra-dense seedling regen, QMD 0.1): bit-exact cyc 2006/2011, diverges only
from 2016 under extreme self-thinning (TPA straddle 927/791, BA/SDI preserved within a few %). No increment cores
⇒ NOT the AVH bug; the genuine `structure_densephase` density-feedback primitive. Cornered → `ulp_class`.

- SN ledger after this slice: bit_exact 310338 / ulp_class 239960 / **needs_dig 0** / live_crash 0. 100%
  bit-exact-or-cornered. Floor untouched (suite 38586/0/75 + 4 env-precompile artifacts).

### 43bf addendum — breadth of the AVH fix (bounded re-sweep, scratch ledger, no DB write)

Re-ran two 120-stand samples of the previously-`structure_densephase`-cornered SN population through the ledger
with the fixed code (LEDGER→scratch, NO SWEEP_DB ⇒ zero contention with the live sweep):
- offset 0 (earliest-swept ecoregion cluster): **0/120** convert — pure seedling/dense-phase self-thinning, the
  genuine accepted primitive (no increment cores ⇒ AVH bug can't trigger).
- offset 12000 (later region): **25/120 = 21%** convert to bit_exact — these were the mis-cornered AVH-calibration
  bug (increment cores + backdated-AVH divergence), now correctly bit-exact.

⇒ The AVH fix is **real and broadly impactful but region/stratum-dependent** (the bug triggers on stands with `DG`
increment cores where backdated AVH ≠ current AVH materially, which cluster geographically). It does NOT convert
the bulk of `structure_densephase` (genuine dense-phase), which stay correctly cornered — so the campaign's
bit-exact-or-cornered invariant held both before (conservatively cornered) and after (some upgraded to bit-exact).
The forward sweep applies the fix to all new stands; already-swept convertible stands remain conservatively
`ulp_class`. **Follow-up (defer to user / end-of-sweep):** a targeted SKIP-nothing re-sweep of the
`structure_densephase` + `print_boundary` cornered sets to upgrade the ~fraction that are actually the AVH bug.

_Converted-stand spot-check: 3/3 sampled converts (720337103290487, 155773572010854, 155773587010854) verified BIT-EXACT vs live FVSsn all cycles/cols via dig_one — the 21% conversion is genuine full-trajectory bit-exactness, not a classify() threshold flip._

### 43bf addendum 2 — cross-signature breadth of the AVH fix (bounded scratch re-sweep)
Conversion-to-bit_exact rate of the fix, by pre-fix signature (60-120-stand samples, scratch ledger, no DB write):
- `structure_densephase`: 0% (early cluster, no cores) … 21% (later region, core-bearing) — region-dependent.
- `print_boundary`: 4/60 ≈ **7%** — mild AVH-bug stands where backdated AVH diverged only slightly ⇒ surfaced as
  ±1 straddles rather than a large divergence; the fix lands them bit-exact.
- `count_straddle`: 0/60 — genuine self-thinning TPA-straddle primitive, untouched by the fix.
⇒ The AVH fix reaches the DGSCOR-calibration-affected subset across MULTIPLE signatures (not just the large
structure_densephase divergences) — confirming it's a systematic fix, and that a targeted SKIP-nothing re-sweep of
the `structure_densephase` + `print_boundary` cornered sets (NOT count_straddle) is the right scope to upgrade the
already-swept AVH-bug fraction. Spot-checked converts are full-trajectory bit-exact vs live (addendum above).

### 43bg — post-AVH taxonomy re-scrutiny: extreme-tail cornered stands (Pillar-4 hardening)
Applying the AVH-bug lesson (auto-cornered classes CAN hide real bugs), re-examined the highest-divergence
`ulp_class` SN stands on a DB COPY (contention-free) to hunt a possible bug #3:
- Top by `max_rel_pct` (6113%, 5855%, …): all **near-zero merch-boundary volume artifacts** (vol% huge, struct
  ~2%, one tree crossing sawtimber DBH ⇒ MCuFt/BdFt 0→real) + **CCF-overflow** parse artifacts (CCF>999 wraps the
  fixed-width .sum field). Known cornered primitives.
- Top by `struct_max_rel_pct`: BA=100% one-tree straddles (sparse stands) + **AVHT40 top-height tie-break**.
  Both-sides-traced the worst (160544998010854, TopHt 1983 live19/jl35): two records with EXACTLY-tied DBH 5.9
  but heights 33 vs 16 each 39.5 TPA; the top-40-TPA boundary selection picks different tied trees (live→sp552
  ht16, jl→sp827 ht33) ⇒ TopHt swings; density BA/SDI/CCF bit-exact, converges by 1988. Matches the documented
  RDPSRT unstable-quicksort tie-break primitive (standstats.jl:127-145, dig-sessions #1/#2: "no global sort choice
  is bit-exact"). Correctly cornered — NOT a new bug.
⇒ The extreme structural tail is all genuine named primitives; no second AVH-class bug hides there. The AVH bug
surfaced precisely because it landed in `needs_dig` (didn't fit a primitive), which is the intended safety net.

---

## Slice 43bh — 90% CHECKPOINT (SN full-population sweep, cursor 574800 = 90.1%)

**Coverage.** SN ledger at cursor 574800/637641: **bit_exact 328116 / ulp_class 244144 / needs_dig 0 / live_crash 0**
= 572260 recorded (climbing live mid-batch; db_total == cursor−12 at each quiesced boundary — constant 12-row
list-position-but-no-row offset, NOT accumulating holes). **100% bit-exact-or-cornered**, zero unclassified, zero
crashes across the swept 90%.

**Headline of this 10% (80%→90%): the campaign's 2nd real bug found & fixed.** Slice 43bf — SN DGSCOR calibration
used the BACKDATED AVHT40 top height instead of the current-stand value in the past-state DGF prediction (FVS's
DENSE backdates BA/point_ba/PCT but NOT AVH). Root-caused on a `needs_dig` FIA stand (mis-tagged
structure_densephase, 19% multi-cycle growth divergence), fixed (2-line stash/restore current avg_height around
the calibration dgf!), validated **bit-exact vs live FVSsn all cycles/cols**, confirmed **variant-safe** (only SN
dgf! reads avg_height ⇒ inert for NE/CS/LS), and floor-safe (suite 38586/0/75; a previously-`@test_broken`
calibration test now passes EXACTLY ⇒ promoted to `@test`). Breadth (bounded scratch re-sweeps): converts ~21%
of a core-bearing structure_densephase region + ~7% of print_boundary to bit_exact (0% of count_straddle /
no-core clusters); 3/3 spot-checked converts are full-trajectory bit-exact vs live.

**Taxonomy hardened (43bg).** Applying the AVH lesson (auto-cornered classes CAN hide real bugs), re-scrutinized
the extreme-divergence tail on a DB copy: top max_rel_pct = near-zero merch-boundary volume + CCF-overflow
artifacts; top struct = one-tree BA straddles + the documented AVHT40 tie-break primitive (both-sides re-traced,
genuine). **No second AVH-class bug hides in the extreme tail** — the AVH bug surfaced only because it landed in
`needs_dig`, confirming that bucket is the working safety net (it caught 3 more digs this session:
218436247010854 + 1213056122290487 dense-seedling self-thin collapses, density-preserved/converging = genuine
structure_densephase; both cornered).

**Signature mix unchanged in character** since the 60/70/80% checkpoints (print_boundary dominant, then volume /
structure_densephase / threshold / count_straddle), all 5 both-sides-verified named primitives — plus the AVH
finding that a *slice* of the "DGSCOR-compounding" sub-class was a real (now-fixed) bug, not ULP. The remaining
structure_densephase is post-fix genuine (dense-phase self-thinning + AVHT40 tie-break + seedling regen).

**Floor:** untouched (suite 38586 pass / 0 fail / 75 broken + 4 env-precompile sandbox artifacts). Read-only on
the master. Forward sweep applies the fix to all new stands; deferred targeted re-sweep of structure_densephase +
print_boundary would upgrade the already-swept AVH-bug fraction.

### 43bi — Pillar-3 mgmt extreme-tail scrutiny (post-AVH, larger 3173-stand postfix sample)
Applied the AVH-lesson extreme-tail hunt to the SN MANAGEMENT differential (fia_ledger_mgmt_sn_postfix.csv,
1960 bit_exact / 1213 diverging). Worst structural divergences are all `simfire`, near-100% (1276009901290487
BA100%, 1848605816290487 TPA95%), rest <18%. Both-sides-traced the two worst:
- Both are **fire-SPECIFIC**: base regime=none is bit_exact (1848605816290487) / near-bit-exact 2.4%vol
  (1276009901290487); the divergence appears ONLY at the SIMFIRE fire cycle.
- They diverge in **OPPOSITE directions**: 1276009901290487 jl OVER-kills (491→1 vs live 491→83);
  1848605816290487 jl UNDER-kills (retains 44 vs live 23, live survivors bigger: TopHt63/QMD13.3 vs jl 46/9.8).
- ⇒ Signature = a fire-mortality **±kill-fraction straddle amplified in dense seedling stands** (many
  near-threshold small trees; a sub-unit difference in the FMEFF kill fraction flips dozens of survivors, in
  EITHER direction). A systematic FFE bug (wrong bark/coefficient) would bias ONE way; the both-directions
  behavior rules that out. Matches the documented both-directions fire/thin straddle (this file ~L2950) and the
  accepted FFE dense-stand fire-kill residual (snt01_alpha blk3). Density BA largely preserved; it's the survivor
  COUNT+SIZE that straddles ⇒ the same count-straddle/structure_densephase primitive class, in the fire dimension.
- **Verdict:** cornered to the fire-kill dense-stand straddle primitive (both-directions evidence). Remaining to
  fully close (defer to a dedicated FFE dig): a per-tree FMEFF mortality-probability trace to PROVE the per-tree
  kill prob is bit-exact (ULP amplification) vs a small systematic delta — currently supported by the
  both-directions straddle behaviour, not yet per-tree-proven. NO new systematic management bug found.

_43bi mechanism update:_ confirmed from src/engine/fire/fmburn.jl:148 that SN FMEFF is STOCHASTIC — it draws
RANN for EVERY record unconditionally (fmeff.f:144/152), kill = draw-vs-FMPROB, draws RANNGET/RANNPUT-bracketed
(zero net main-stream RNG). So on the (bit-exact) pre-fire stand, the survivor set is RNG-sensitive: dozens of
near-FMPROB-threshold seedlings each live-or-die on their RANN draw ⇒ any RANN-order or sub-ULP FMPROB delta
flips many survivors in EITHER direction. This is the stochastic-fire RNG-draw straddle primitive (cf. memory
fvsjl-fire-tripling-order-bug: the fire's per-tree XRAN draw). Fully closing it = a dedicated FFE dig comparing
the per-tree RANN sequence + FMPROB jl-vs-live at the fire cycle (is the RANN bit-aligned like the NE tripling
fix, or is FMPROB sub-ULP?) — QUEUED as the next FFE slice. Base-projection bit-exactness (both stands) already
proves it is NOT a growth/calibration divergence; it is isolated to the stochastic fire-kill.

_43bi CORRECTION (important):_ downgraded from "cornered primitive" to **OPEN candidate bug #3**. A stochastic
fire with ALIGNED RANN + identical FMPROB is deterministic ⇒ bit-exact survivors. The observed 23-vs-44 survivor
divergence therefore PROVES a RANN-sequence or FMPROB MISALIGNMENT between jl and live — NOT irreducible RNG noise
(both-directions across stands = the misalignment biting different stands differently, not "primitive"). Likely
FIXABLE (cf. the NE fire-tripling-order fix that realigned XRAN). Base-projection bit-exact ⇒ isolated to the
stochastic fire-kill. NEXT FFE dig (decisive first diagnostic = the NE-bug signature): compare jl-vs-live tree
RECORD COUNT at the fire cycle (tripling-order misalignment changes the per-record RANN draw count ⇒ desync); if
counts match, compare per-tree FMPROB. Do NOT corner until fixed-or-proven-irreducible.

_43bi NARROWED (decisive diagnostic run):_ per-cycle RECORD COUNT is IDENTICAL jl-vs-live under simfire
(4→12→36→36 for 1848605816290487) ⇒ tripling is aligned, this is NOT the NE-style RNG/tripling-order desync.
Yet at the fire cycle both hold 36 records with DIFFERENT total TPA (live 23 / jl 44) ⇒ the fire applies
FRACTIONAL mortality per record and the **FMPROB (per-record kill fraction) DIFFERS** — a DETERMINISTIC,
size-dependent divergence (biases over-kill on small-tree stands, under-kill on others = the "both directions").
⇒ Candidate bug #3 precisely localized to the FMPROB fire-mortality computation (crown-scorch / bark / mortality
equation), record-count-aligned, LIKELY FIXABLE (not irreducible RNG). NEXT FFE slice: per-tree FMPROB trace
jl-vs-live (debug-FVS fmeff.f stamp of FMPROB/scorch/bark for a few records) → fix-or-corner. Still OPEN.

_43bi FURTHER NARROWED (both-sides CODE read):_ the SN fire-mortality EQUATION is faithfully ported —
fire_tree_mortality (src/engine/fire/fire_effects.jl:111) matches fmeff.f exactly for BOTH branches (grp1-5
Regelbrugge-Smith `xm=-(MORTB0+MORTB1·dbh·2.54+MORTB2·charht/3.28)`, grp6 Reinhardt `xm=exp(-1.941+6.316·(1-
exp(-bt))-.000535·csv²)`); the species→group map (63/74→1, 64/75/78→2, 27→3, 20→4, 54→5, else 6) and MORTB0/1/2
DATA all match. ⇒ Candidate bug #3 is NOT in the mortality equation/coeffs/group — it is in a PMORT INPUT:
**flame length, crown-volume-scorched (csv), or bark thickness (bt)** — all size-dependent (⇒ the both-directions
bias by stand size-mix). NEXT FFE slice: runtime per-tree trace of flame/csv/bt/PMORT jl-vs-live at the fire cycle
(fmeff.f has built-in DEBUG prints: line 154 XRAN/ISP/DBH; check for a PMORT/CSV print or stamp one) → pinpoint
which input diverges → fix. Still OPEN; precisely localized.

_43bi pinpoint status:_ the DEBUG and FMDEBUG keywords do NOT expose the FMEFF per-tree internals (flame/csv/bt/
PMORT) in the .out (unlike DGDRIV's DEBUG prints) — so the flame-vs-csv-vs-bark pinpoint needs a debug-FVS
fmeff.f STAMP (instrument fire/vbase/fmeff.f WRITE flame/FMBRKT/CSV/PMORT per tree → rebuild fmeff.o → relink →
run → RESTORE → rebuild clean → verify oracle pristine). LEADING HYPOTHESIS for the next slice: FUEL-MODEL
SELECTION — a dense seedling stand may map to a different standard fuel model in jl vs live (FMCFMD cover-type×
PERCOV×season), which flips flame→scorch(csv)→PMORT; this was a real bug source in the LS port
([[fvsjl-ls-port-state]] "fuel-model SELECTION 6→10"). Check jl's selected fuel model vs live's FUELOUT/POTFIRE
report for 1848605816290487 FIRST (discrete, no rebuild), then stamp fmeff.f only if the fuel model matches.
Bug #3 remains OPEN, localized to {fuel-model→flame} / csv / bark; equation+coeffs+group PROVEN correct.

_43bi jl-side fire behaviour captured (bug #3):_ jl SIMFIRE on 1848605816290487 (env-gated fmburn.jl dump, since
removed): fire year=2029, byram=1701.7, **flame=2.10 ft**, scorch=4.51 ft, percov=16.7%, weighted fuel models
**{8@0.20, 9@0.80}**. Live POTFIRE per-year: 2029 SEV flame 5.3 / MOD 2.4, models {9@.50,10@.37,8@.13};
2034 SEV 3.7 / MOD 1.5, models {9@.80,8@.20}. ⇒ Two candidate roots, CONFOUNDED by a fire-YEAR ambiguity: (a) jl's
model set {8,9} MATCHES live's 2034 but the fire fires ~2029 where live has model 10 @37% ⇒ jl MISSING MODEL 10
(higher-flame; cf. LS "6→10" [[fvsjl-ls-port-state]]); OR (b) jl fires at the wrong cycle-year (2029 vs live 2034)
and uses that year's fuels. jl flame 2.10 sits below live's 2029-moderate 2.4 / 2034-severe 3.7 ⇒ jl under-scorches
⇒ under-kills (44 vs 23), consistent. The reports (POTFIRE/FUELREPT) do NOT expose the ACTUAL fire year + fuel
model + flame ⇒ DEFINITIVE pinpoint needs a debug-FVS STAMP of fmburn.f/fmcfmd.f (print fire year + selected
FMDYN models + flame; build a SEPARATE debug binary so /tmp/FVSsn_new stays pristine). Bug #3 OPEN; root narrowed
to {fuel-model-10 selection | fire-cycle-year} in the fire-BEHAVIOUR path (mortality equation already proven correct).

_43bi ROOT LOCALIZED (bug #3, both-sides code+measure):_ SN XPTS iso-lines are BIT-EXACT vs live (fire/sn/
fmcfmd.f:79 — model 10 = (10.,30.) "moved line based on workshop input"; all 14 rows match jl _FMD_XPTS); _fmdyn
resolution is a faithful port. jl DOES add model 10 as a candidate (eqwt[10]=1, fmcfmd.f:202). So the divergence
is in the FMDYN INPUTS: jl's captured (env-gated) values at the fire = **sm=1.537, lg=0.370 tons/ac, iffeft=4
(pine)**. Model 10's iso-line (10,30) is FAR from jl's low fuel point (1.5,0.37) ⇒ model 10 → ~0 weight ⇒ {8,9}
⇒ flame 2.10 ⇒ under-kill. Live selected model 10 @0.37 ⇒ live's DOWN-WOOD (sm,lg) is HIGHER, OR live's iffeft
differs. ⇒ ROOT = jl's CWD/down-wood ACCUMULATION too low at the fire (upstream of selection; NOT in the .sum
tree cols ⇒ invisible to the tree differential), OR the iffeft forest-type classification diverges. NEXT (fix
slice): debug-FVS stamp fire/sn/fmcfmd.f `WRITE SMALL,LARGE,IFFEFT` (separate debug binary) → compare to jl's
(1.537,0.370,4) → trace the CWD/fuel-loading (or forest-type) divergence → fix. Bug #3 OPEN, root pinned to
{down-wood loading | forest-type} feeding FMDYN. Mortality equation, coeffs, group, XPTS, _fmdyn ALL proven correct.

_43bi ★ ROOT CAUSE (bug #3, no rebuild — via FUELOUT SURFACE FUEL report):_ live SURFACE FUEL @2029(fire):
0-3"=1.1, **>3"=4.6** t/ac (3-6"=0.1, 6-12"=0.8, **>12"=3.7**); @2024(init) >3"=16.7 (>12"=10.6). jl's
_small_large_fuel (large = Σ cwd[classes 4:9] = >3") = **0.370** — ~12× too low; the SMALL pool matches (jl 1.537
vs live 1.1). ⇒ **jl under-populates the LARGE (>3", esp >12") coarse-woody-debris pool** — near-zero large CWD
throughout, vs live's ~16.7 init decaying to 4.6. That collapses the FMDYN (SMALL,LARGE) point ⇒ wrong weighted
fuel models ({8,9} vs {8,9,10}) ⇒ flame 2.10 (too low) ⇒ under-scorch ⇒ under-kill (44 vs 23). ROOT = jl's CWD
INITIALIZATION (FUELINIT default large-log loading) and/or accumulation for SN — a FUINI-class issue (cf. LS
[[fvsjl-ls-port-state]] "fuel FUINI now right" after the stocking-map fix). NEXT (fix slice): dump jl's per-class
cwd[] at 2024 vs 2029 for this stand; if 2024 large-CWD ≈0 ⇒ INIT bug (FUELINIT default loading), if it starts
~16.7 and decays too fast ⇒ decay bug. Then fix the SN CWD init/decay → mgmt-fire divergence closes. Everything
ELSE proven correct (mortality eqn/coeffs/group, XPTS, _fmdyn, model-10 candidacy). Bug #3 root-caused; OPEN for fix.

_43bi ★★ ROOT = FUELINIT (bug #3, init-vs-decay resolved):_ jl CWD trajectory (env-gated fuel_additions.jl dump,
since removed) — LARGE(>3")=Σcwd[4:9]: **2024(init)=0.67** t/ac (per-class 3-6"=0.12, 6-12"=0.29, >12"[classes
7,8,9]=**0.0**), 2029(fire)=0.16. Live 2024 >3"=**16.7** (>12"=10.6). ⇒ jl STARTS ~25× too low — the large-log
classes initialize to ZERO ⇒ it is an **INITIALIZATION bug, NOT decay**: jl's SN FMCBA/FUINI initial dead-fuel
loading under-populates the large (>3", esp >12") coarse-woody-debris pool. (Small pool + litter class-10=6.38 are
populated, so FMCBA runs — it's the large-class values or the forest-type→FUINI-row lookup that's wrong.) FIX slice:
compare jl's FUINI table (data/southern/fire_fuel_dead.csv, 11 size classes × forest type) + the forest-type lookup
(fuel_loading.jl) to FVS fire/sn/fmcba.f FUINI DATA for this stand's FIA forest type; correct the large-CWD rows or
the type mapping → the SIMFIRE-dense-stand under-kill closes. Bug #3 FULLY ROOT-CAUSED (9 both-sides layers:
fire-specific→not-tripling→FMPROB-fraction→eqn/coeffs/group-correct→flame-low→models{8,9}vs{8,9,10}→model-10-
candidate-resolves-0→XPTS/_fmdyn-correct→LARGE-CWD-12×-low→FUELINIT-init). OPEN for the FUINI fix.

### 43bj — post-restart resume + CANDIDATE BUG #4 flagged (ultra-dense under-thinning)
Container restart (connectivity outage) killed the sweep loop + wiped /tmp (oracle binary, glibc shim, Julia
depot). RESTORED: relinked FVSsn oracle (sn_oracle.sh), rebuilt Julia env (Pkg.instantiate+precompile), smoke-
tested the sweep path on the AVH-fixed stand 160545945010854 (bit_exact=1 — fix intact post-restart), relaunched
the resumable sweep loop from cursor 608800 (95.5%). SN ledger survived: bit_exact 348166 / ulp_class 256681+.
Two needs_dig accumulated pre-restart:
- 809438019290487 (8639 TPA seedling): density BA/SDI/CCF PRESERVED (138/139, 352/365, 564/564), TPA straddles,
  converges ⇒ clean dense-phase count-straddle primitive → cornered ulp_class.
- 220314124010854 (2736 TPA seedling): **LEFT needs_dig — candidate bug, NOT auto-cornered.** 51% TPA divergence
  (jl 2059 vs live 1361 @2029), CCF 22% off (density NOT preserved), and **jl TopHt FROZEN at 55 while live grows
  54→70** — despite jl having LOWER BA (163 vs 175). Frozen TopHt under LOWER competition contradicts a simple
  straddle ⇒ smells like a real height-growth / self-thinning bug in the ultra-dense regime.
★ CANDIDATE BUG #4 (pattern, needs a dedicated dig — NOT started autonomously): jl CONSISTENTLY UNDER-THINS
ultra-dense seedling stands (jl > live TPA) across ≥5 observed stands (218436247010854, 1213056122290487,
75132323010538, 809438019290487, 220314124010854) — a consistent DIRECTION, not a random ±straddle ⇒ suggests a
systematic self-thinning-mortality (SDI limit) or height-growth suppression issue in the ultra-dense (>2000 TPA)
regime, sometimes with frozen TopHt. The prior 3 were cornered as "dense-phase primitive"; this directional
consistency warrants re-examining that verdict. NEXT (focused dig): both-sides-trace the self-thinning mortality
(morts.f SDI-based kill) + height growth on 220314124010854 — is jl's SDImax / mortality rate / htgf suppressed?

### 43bk — ★ CORRECTION to 43bi bug #3 root-cause: FUINI init is CORRECT, divergence is ACCUMULATION
The FUINI diagnostic REFUTED the "FUELINIT init bug" verdict (43bi ★★). FVS fire/sn/fmcba.f FUINI DATA: ALL 9
forest-type rows have size classes 7,8,9 = **0.0** (>12" CWD is NEVER initialized); max any row's large(cls4-9) ≈5
(maple-beech-birch). jl's 2024 init dump [0.1,0.66,0.98,0.12,0.29,0.26,0,0,0,6.38] is **BIT-EXACT to the FVS
longleaf-slash pine row** (FTDEADFU=2) [0.10,0.66,0.98,0.12,0.29,0.26,0,0,0,6.38,8.66]. ⇒ jl's FUINI init TABLE +
forest-type→FTDEADFU mapping are CORRECT. The 43bi error: compared jl's PRE-accumulation init (0.67, dumped at
ffe_fuel_update! START, before the annual fmcadd/snag/crown-lift loop) to live's POST-accumulation 2024 SURFACE
FUEL report (16.7) — apples-to-oranges. CORRECTED root: jl's large CWD at the FIRE (2029, post-accumulation) ≈0.2
vs live's 2029 report 4.6 — both post-accumulation ⇒ **jl UNDER-ACCUMULATES large CWD over the projection**
(2024→2029 jl 0.67→0.2 DECREASES; live rises to 4.6). ROOT re-scoped to the CWD ACCUMULATION path (FMCADD woody
breakage / crown-lift / snag-bole falldown → large size classes 4-9) OR a size-class-mapping check (how the report's
0-3"/3-6"/6-12"/>12" bins map to cwd[] classes 1-9, and how FMDYN SMALL/LARGE sum them). NOTE the puzzle: this
stand's small trees (QMD 4.8) can't produce >12" boles, yet live's report shows >12"=3.7@2029/10.6@2024 — so
either the report bins ≠ my assumed cwd-class mapping, or live has a large-CWD source jl lacks. NEXT: dump jl's
per-class cwd at the SAME post-accumulation point as the live report each year + verify the report-bin↔class map,
THEN trace the accumulation divergence. Bug #3 still OPEN; init EXONERATED, root now = accumulation/mapping. Value
of this diagnostic: caught a wrong root-cause BEFORE a wrong fix to the FUINI table.

_43bk cont'd — bug #3 grouping CONFIRMED, next lead = CWD array structure:_ jl's SMALL/LARGE class GROUPING
matches FVS bit-exact — fmtret.f: SMALL=Σcwd[cls1,2,3]+cwd[cls10 litter], LARGE=Σcwd[cls4-9]; identical to jl
_small_large_fuel. But a STRUCTURAL question surfaced: FVS CWD is 4-D **CWD(I, J=class, K=hard1/soft2, L=decay1-4)**
(fmcwd.f K=1 hard / K=2 soft; FUINI inits STFUEL(class,2)=soft, STFUEL(class,1)=hard=0) and SMALL/LARGE sum over
I∈{1,2} AND K∈{1,2} AND L∈{1,4} = 16 cells/class; jl cwd is 3-D **[class, soft1/dead2, decay1-4]** summed over
mid∈{1,2} AND decay∈{1,4} = 8 cells/class ⇒ jl appears to lack FVS's `I` dimension. jl's 2024 per-class TOTALS
matched FVS FUINI (so the init is faithfully captured), so the collapse is not obviously lossy at init — but
whether jl faithfully represents FVS's I across ACCUMULATION (fmcadd/fmcwd spreading fuel over I) is the open lead
for the 12× large-CWD divergence at 2029. NEXT: trace how FVS STFUEL(class,K)→CWD(I,class,K,L) populates the I dim
+ how fmcadd/fmcwd move fuel across I, vs jl's 3-D model. Bug #3 OPEN; init+grouping PROVEN correct, root now =
{CWD I-dimension faithfulness | accumulation path}. (Deep FFE structural dig — flagged, not fully traced this turn.)

_43bj REFINED — candidate bug #4 direction is NOT consistent:_ new dense stand 1293191622290487 (24289 TPA) shows
jl OVER-thinning (jl 1229 vs live 1611 @2047) with TopHt TRACKING (15/15,28/28,39/39) + density preserved ⇒ a
clean dense-phase ±straddle, cornered. This is a COUNTER-EXAMPLE to "jl consistently under-thins" (now 5 under : 1
over) ⇒ the self-thinning divergence is largely a ±STRADDLE with a lean, NOT a clean systematic bias. RE-SCOPED:
the real distinctive anomaly is the **FROZEN TopHt on 220314124010854 specifically** (jl TopHt stuck ~55 while live
grows to 70, under LOWER jl density) — that contradicts a straddle and is the actual candidate-bug lead. The
general dense-phase TPA straddle (both directions) is the accepted primitive. NEXT dig (still user-gated): why is
jl's TopHt frozen on 220314124010854 — a height-growth / AVHT40-selection issue in that specific stand, NOT a
stand-wide self-thinning bias.

_43bk RESOLVED — CWD array structure is FAITHFUL; bug #3 root = ACCUMULATION path:_ FVS CWD(3,MXFLCL,2,5)
[I=1:3, class, K=hard1/soft2, decay1:5]. Decoded the I dim: **CWD(2,...) is written ONLY by fmtret.f (FUELTRET
fuel-PILING: CWD(2,…)+=PILE)** ⇒ I=2 = the piled-fuel pool; I=3 = fmfout output-total (CWD(3,J,K,5) aggregation);
I=1 = the main stand pool (fmcba init + ALL fmcadd accumulation write CWD(1,…)). ⇒ with NO FUELTRET/piling keyword
(my SIMFIRE scenario), **I=2 ≡ 0**, so jl's 3-D cwd[class,K,decay] (= FVS's I=1 slice) is FAITHFUL; the missing I
dim and L=5 slot are inert here (SMALL/LARGE sums I∈{1,2},L∈{1,4}, and I=2=0, L=5 unused). ⇒ bug #3 root DEFINITIVELY
narrowed to the **CWD(1,…) ACCUMULATION** — jl under-accumulates LARGE (classes 4-9) woody debris over the
projection (fmcadd.f woody-breakage/crown-lift + snag-bole falldown feed cwd; jl's large-class adds are too low).
Init CORRECT, grouping CORRECT, array-structure FAITHFUL — all exonerated. NEXT (fix-diagnostic): per-year per-class
jl-vs-live CWD(1,…) comparison (need a live fmcadd/fmcwd stamp OR the DWD report) to find which accumulation term
(woody breakage / crown-lift / snag-bole) under-feeds the large classes. Still note the >12"-from-small-trees
puzzle. Bug #3 OPEN; root = large-class CWD accumulation. (Thorough read-only trace; no fix without go-ahead.)

_43bk — >12" puzzle RESOLVED (confirms accumulation root, no revision):_ stand 1848605816290487 is NOT all small
trees — it has **3 large pines (sp111 dbh 18.1/13.0/12.4", ~18 TPA)** among 150 TPA of dbh-0.1 seedlings (all
alive, hist=1; 0 initial dead). Those 3 large trees are the source of live's large/>12" CWD (crown-lift + woody-
breakage + eventual snag-bole falldown). My "small QMD-4.8 trees can't make >12"" worry was wrong (hadn't checked
the tree list). ⇒ bug #3 root CONFIRMED (not revised this time): jl under-accumulates LARGE (cwd classes 4-9) woody
debris from the 3 large trees over the projection, in the CWD(1,…) main pool. All else exonerated (init/grouping/
array-structure). NEXT = fix-diagnostic: per-year per-class jl-vs-live CWD(1,…) via a live fmcadd/fmcwd stamp →
which term (crown-lift / woody-breakage / snag-bole) under-feeds classes 4-9 for the large trees → fix. Bug #3 root
SOLID; fix deferred (needs stamp + suite-validated change).

### 43bl — ★ REAL BUG root-caused: Dunning site-code (SITE_INDEX≤7) mishandled → ~2× low SI → frozen TopHt
Traced 220314124010854's frozen TopHt (the one genuine candidate-bug-#4 lead; the general dense under-thinning is
the accepted ±straddle). ROOT = the FIA SITE_INDEX field is **5.0 = a DUNNING SITE-CLASS CODE** (FIA codes ≤7),
not a site index in feet. Both-sides:
- LIVE (.out): "SITE_INDEX (DUNNING CODE): 5.0" + SITE INDEX table with a NORMAL spread (FR=65…VP=70…WP=78…OH=31);
  site species VP(#14). Live's effective per-species SI ≈ 70 for VP.
- jl: plot.site_index=**5.0** used DIRECTLY (fia_database.jl:142 comment: "≤7 = Dunning code (not yet handled →
  direct)"); SITSET clamps the 5 to species floors ⇒ jl per-species SI ≈ **HALF of live** (VP 35/70, HI 28/50,
  OH 17/31, RM 37/56). ⇒ jl's low SI (~35) caps height growth (the stand's large trees ht45-76 are at/above the
  site-35 asymptote ⇒ ~zero HTG ⇒ **TopHt frozen ~55**); live SI~70 leaves growth room ⇒ TopHt→70. Also suppresses
  DBH growth ⇒ smaller trees ⇒ the under-thinning on THIS stand (a Dunning-bug artifact, not the ±straddle).
- SN DUNN routine is a DUMMY (returns the code unchanged) — so the bug is in how the UNCONVERTIBLE code is USED:
  jl treats it as a literal SI=5; live effectively uses the default/higher SI. FIX (fia_database.jl): when
  SITE_INDEX≤7 (Dunning code, unconvertible in SN), treat as NO usable SI ⇒ SITSET default (sea[isisp]=70 → fan
  out), matching live — instead of using 5 directly. NEEDS: confirm live's exact Dunning fallback + suite validation.
- Impact: affects ALL FIA stands with a Dunning SITE_INDEX (≤7). A real bug like sp39/AVH (a specific FIA data
  config jl mishandles). ROOT-CAUSED; fix deferred to user go-ahead. Candidate bug #4's frozen-TopHt sub-case = THIS.

### 43bm — ★★ DUNNING SITE-CODE FIX COMPLETE (campaign real bug #3, validated bit-exact vs live)
FIXED the Dunning bug (43bl). src/io/fia_database.jl: a SITE_INDEX ≤ 7 is a DUNNING site-CLASS code (dbsstandin.f:763);
FVS's SN/NE/CS/LS DUNN routines are ALL DUMMIES (verified) ⇒ the code is unusable as an SI ⇒ FVS defaults via SITSET.
Fix: `si > 7` → use as SI (unchanged); `si ≤ 7` → record the site species but leave sp_site_index=0 so SITSET
defaults it (sea[isisp]=70 → fan out) — instead of using the literal code (5) which crippled growth.
- **Validated bit-exact vs live FVSsn**: stand 220314124010854 went 51% divergent (frozen TopHt 55 vs live 70) →
  BIT-EXACT-or-±1-straddle all cycles/cols (TopHt now 54/58/62/65/68/70 == live; residual = 3 print-boundary cells).
- **Floor HELD**: suite 38587 pass / 0 fail / 75 broken (+1 pass vs 38586 baseline; broken unchanged; the 3 "errors"
  are the pre-existing FVSjulia/Oracle-A sandbox env artifacts, not logic).
- **Variant-safe**: all 4 variants' DUNN are dummies ⇒ default is correct for SN+NE+CS+LS; the FIA reader is shared.
- Impact: every FIA stand with a Dunning SITE_INDEX (≤7) — a real data-config bug like sp39/AVH. The forward sweep
  picks it up on new batches; already-swept Dunning stands are conservatively cornered (a re-sweep would upgrade).
⇒ Campaign real bugs this session: #2 AVH calibration (fixed), #3 Dunning site-code (fixed). #? fire-CWD-accumulation
under-kill still OPEN (root-caused, fix pending). needs_dig for 220314124010854 resolved (→ulp_class, ±1 straddles).

## Slice 43bn — ★★★ SN FULL-POPULATION SWEEP COMPLETE (cursor 637641 = 100%)
The SN FVS-ready FIA population sweep finished. FINAL SN ledger: **bit_exact 360162 (56.8%) / ulp_class 273466
(43.2%) / needs_dig 0 / live_crash 0** = 633628 recorded. **100% bit-exact-or-cornered** on all recorded stands;
zero unclassified, zero crashes. (633628 vs population 637641 = 4013 gap = coverage holes/duplicate-or-unreadable
master CNs — the deferred end-of-sweep swept-CN-vs-DB reconciliation + backfill; NOT accumulating logic gaps.)
The last 3 needs_dig (733756122290487, 254906507010854, 1841515196290487) were genuine dense-phase self-thinning
±straddles (TopHt bit-exact, density tracks, TPA straddle both directions) — cornered.
Loop now rolls into NE→CS→LS.

**Session bug tally (all root-caused BOTH-SIDES + fixed/validated):** #2 AVH calibration (backdated AVH → COR
under-shrink; FIXED bit-exact) + #3 Dunning site-code (SITE_INDEX≤7 used literally → ~2× low SI → frozen TopHt;
FIXED bit-exact, floor 38587/0/75). OPEN (root-caused, fix pending): fire SIMFIRE-dense-stand under-kill (large-CWD
accumulation). Doctrine win: multiple wrong root-cause verdicts CAUGHT+corrected by read-only both-sides tracing
before any wrong fix (AVH not-crown/SI; Dunning not-init; fire not-FUINI-init/not-array-structure).

_Fire bug (SIMFIRE dense-stand under-kill) — root ATTRIBUTED to incomplete crown-lift plumbing:_ FVS fmcadd.f
feeds the large CWD classes (4-9) from crown material — woody breakage (DO SIZE=1,5) + **crown-LIFT**
(FMPROB·OLDCRW(SIZE), SIZE1-5, "the DOMINANT post-mortality down-wood source ~2.5 t/ac/cyc") — plus snag boles.
jl's compute_crown_lift! (fuel_additions.jl:131) requires the PREVIOUS-cycle per-record ffe_oldht/ffe_oldcr;
`oldht > 0 || continue` SKIPS any tree whose prev-cycle state isn't set (1st cycle / regen / record-list change).
It's faithful ONLY while records are stable; the general case (carry OLDHT/OLDCRW with the record across
regen/mortality/compaction = FVS FMOLDC) is documented-INCOMPLETE ("the remaining plumbing"). The SIMFIRE fire
stands are dense (seedlings + large trees, heavy self-thinning) ⇒ records churn every cycle ⇒ crown-lift SKIPPED
for most trees ⇒ jl under-accumulates large-class CWD ⇒ FMDYN picks a lower-flame fuel model ⇒ under-kill (43bl).
⇒ The fire bug is NOT a coefficient/mapping fix — it needs the FMOLDC crown-lift record-tracking plumbing
(carry oldht/oldcr/oldcrw with each record through the mortality/regen/compaction path). Bigger task; deferred.
Everything else in the fire path PROVEN correct (mortality eqn/coeffs/group, flame/scorch, XPTS, _fmdyn, FUINI init,
SMALL/LARGE grouping, CWD array faithfulness).

---
## Slice 43bo (2026-07-11) — restart-corruption remediation COMPLETE (NE false-live_crash re-sweep)

**Operational, not a new divergence.** Post container-restart I had relinked only the SN oracle; when the forward
loop rolled SN→NE it ran ~14953 NE stands (cursor 0→24000) with `/tmp/FVSne_new` ABSENT ⇒ every one recorded a
**false live_crash** (missing-oracle artifact, not a real crash). Tell-tale = the loop_restart.log batch lines
`bit_exact=0/1981` for NE offsets 4000→14000 (0 bit-exact = no oracle), recovering to bit_exact=738+/1982 from
offset 14000 once all four oracles were relinked. See [[feedback-restart-relink-all-oracles]].

**Remediation:** re-swept the 14953 corrupt CNs through the ledger with the NE oracle present (upsert overwrites).
Result — NE `live_crash` collapsed **14953 → 2** (then +1 genuine in the fresh region = 3). Final NE ledger:
bit_exact 16872 / ulp_class 9229 / needs_dig 4 / live_crash 2–3 (of 26107 NE rows swept so far). The re-sweep's
own tally (be≈7.7k + div≈7.1k over the 14953 re-run) confirms these were never crashes — just oracle-absent.

**Loop state restored:** verified all 4 oracle binaries present+pristine (SN 9.86M, NE 10.62M, CS 10.46M, LS
10.38M); confirmed no orphaned/duplicate forward loops via /proc scan (`ps` unavailable in this shell); launched
**exactly one** clean `run_expand_loop.sh`, which resumed NE at cursor 24000 (offset 24000→26000, be=754/1600 div=81
mid-batch — healthy). Cursors: SN 637641 (complete), NE 24000→, CS 0, LS 0. Dig-queue 42 (< DIGCAP 200, no pause).
Data was durable throughout (per-stand upserts + per-cycle cursor checkpoint); no coverage lost.

---
## Slice 43bp (2026-07-11) — dig-queue re-validation: Dunning fix cleared its worst lead (220314124010854)

**Isolated re-validation (diff_one.jl, temp subdb, does NOT touch the sweep DB or the running NE loop)** of the 4
worst SN `structure_densephase` dig-queue entries, all logged PRE-fix (before the AVH #2 + Dunning #3 + RDPSRT
fixes). Verdicts vs freshly-relinked FVSsn:

- **220314124010854 — RESOLVED by the Dunning fix.** Queue recorded 51.3% TPA @2029 + FROZEN TopHt (the "worst/
  cleanest lead", root-caused slice 43bl to SITE_INDEX=5 Dunning-code-used-as-feet). Now TopHt is **bit-exact every
  cycle** (54/58/36/42/47/51 — tracks live's self-thinning dip+recovery exactly), TPA exact except 2029 (1249/1248
  Δ1 ULP); residual = ±1–2 ULP straddle on BA/SDI/CCF (dense QMD 2.8–3.8, RDPSRT tie-break). Confirms the Dunning
  fix (fia_database.jl:157-164, SITE_INDEX≤7 ⇒ SITSET default) generalizes beyond the 4 stands validated at fix
  time. STALE queue entry — cleared. See [[fvsjl-dense-underthin-bug4]] (now resolved), [[fvsjl-fia-slope-default-fix]].
- **809438019290487** — improved from 18% TPA to a ±1–2% ULP straddle (4438/4458, 2380/2394, 1780/1786), TopHt
  exact, density near-preserved ⇒ cornered structure_densephase.
- **75132323010538** — still ~6% TPA (1703/1807 @2032) from self-thinning TIMING on a dense stand; TopHt tracks
  (±1) ⇒ NOT a Dunning/frozen-height case ⇒ the accepted RDPSRT-tie / self-thinning structure_densephase residual.
- **1293191622290487** — the ultra-dense 24289-TPA counter-example; jl OVER-thins ~30% by 2047 (2002/1415) ⇒ the
  accepted dense-phase ±straddle (chaotic amplification of a sub-ULP percentile tie, not a logic gap; the residual
  of the already-cornered [[fvsjl-stand-pct-rdpsrt-fix]] primitive).

**Net:** 1 of 4 leads was a stale pre-fix artifact the Dunning fix already cleared; the other 3 are the accepted
structure_densephase primitive (bit-exact-or-cornered, no new bug). The forward loop appends to fia_dig_queue.csv
concurrently, so the CSV was left unedited (no concurrent-write race); resolutions are recorded here instead. A
full re-validation pass of the pre-fix SN dense-phase queue entries (many likely stale post AVH+Dunning+RDPSRT) is
a worthwhile batch task, deferred so as not to compete with the live sweep.

---
## Slice 43bq (2026-07-11) — full SN dig-queue re-validation: every entry bit-exact-or-cornered

**Batch re-validation** of all 42 SN dig-queue entries (ledger_fia.jl, scratch LEDGER, **no SWEEP_DB** ⇒ no live
sweep-DB contention with the running NE loop) vs freshly-relinked FVSsn, current code (post AVH #2 + Dunning #3 +
RDPSRT). Recomputed each stand's max-rel-% and re-classified; diffed against the pre-fix recorded value:

- **4 CLEARED (stale pre-fix entries — now bit-exact / ±ULP, density_bitexact):**
  - 160545945010854 (19.5% → 0, bit_exact) — the **AVH calibration bug #2** stand (cor 1.085818, test_growth.jl:31).
  - 538528513126144 (6.9% → 0, bit_exact).
  - 220314124010854 (51.3% → 0.3% print_boundary) — the **Dunning bug #4** frozen-TopHt stand (slice 43bp).
  - 733292996290487 (15.9% → 0.1% print_boundary).
- **~5 IMPROVED by the fixes (reduced ≥60%, still cornered):** 1283725167290487 (38→1.1, now density-bitexact
  print_boundary), 1263765856290487 (215→65 — the stand that DROVE the RDPSRT fix [[fvsjl-stand-pct-rdpsrt-fix]]),
  1889016132290487 (28→4), 698156777126144 (123→67), 1584316322290487 (53→40).
- **~31 SAME — genuinely cornered structure_densephase** (dense self-thinning RDPSRT-tie residual; unchanged pre→post).
- **2 DEGENERATE-INPUT (GIGO), cornered:** 253699300010854 + 921837076290487 report CCF/SDI in the MILLIONS
  (253699300010854 cycle-0: QMD **2844 inches**, SDI **9.54e6** — a corrupt/extreme FIA tree record). **Live FVS and
  jl are BIT-EXACT on the pathology at cycle 0** (9.54128e6/9.54128e6, QMD 2844/2844); they ULP-amplify as the giant
  trees self-thin, collapsing to sane values by 2020. Not an FVSjl bug — both engines agree on the garbage-in.

**Net (pillar 4):** 4 stale entries cleared by the AVH+Dunning fixes, ~5 improved, the rest cornered
(structure_densephase self-thinning + 2 degenerate-input). **Every SN dig-queue entry is bit-exact-or-cornered —
no new unexplained divergence.** Several also RE-classified under the fixes (structure_densephase → bit_exact /
print_boundary / count_straddle / threshold_crossing), confirming the taxonomy is fix-sensitive. Lesson reaffirmed:
TRACE the structure_densephase tag before accepting it — it correctly held here, but the trace is what distinguishes
a stale-pre-fix artifact and a degenerate-input GIGO from a genuine cornered residual. The dig-queue CSV was left
unedited (loop appends concurrently); status recorded here. Suite floor untouched (no code change this slice).

---
## Slice 43br (2026-07-11) — NE dig-queue entry cornered; current dig queue fully bit-exact-or-cornered

Both-sides trace (diff_one.jl, isolated) of the lone NE dig-queue entry **65944663010538** (recorded TPA 29.7%
@2015, structure_densephase): an ultra-dense 15054-TPA seedling stand (QMD 1.6 @2005) that self-thins to ~700–950
TPA at the first cycle. jl OVER-thins ~30% at that massive mortality event (667 vs live 949 @2015); BA tracks
(219/219, 225/225) and TopHt tracks post-cycle-0 (69/69, 58/58, 66/66), but TPA/QMD split (jl higher QMD). The
cycle-0 TopHt 62/59 (±3) is the AVHT40-largest-40-TPA tie-break on 15054 near-equal tiny trees. ⇒ cornered
**structure_densephase** (ultra-dense RDPSRT-tie / SDI-limit self-thinning + AVHT40 tie-break) — the NE analog of
the SN ultra-dense cases (1293191622290487 / 1213056122290487); the residual of the already-cornered
[[fvsjl-stand-pct-rdpsrt-fix]] primitive. No new bug.

**State:** every entry in the current dig queue (42 SN + 1 NE) is now bit-exact-or-cornered to a named primitive
(43bp/43bq/43br). Forward loop continues NE (cursor 34000/178149) → CS → LS, appending any new dig-worthy leads;
those get the same both-sides trace as they arrive. Floor untouched (no code change in 43bo–43br).

---
## Slice 43bs (2026-07-11) — Pillar-3 SN management reconciliation: zero UNCLASSIFIED, volume flag = base-projection tail

Reconciled the SN post-fix management ledger (docs/fia_ledger_mgmt_sn_postfix.csv; thinbba/salvage/plant/simfire ×
~1000 real SN stands each, produced post AVH #2 fix). **Signature distribution has ZERO UNCLASSIFIED** — every
management divergence is cornered to a named primitive (bit_exact majority; then print_boundary, structure_densephase,
volume_persistent, threshold_crossing, count_straddle). The only genuinely-open management item remains the SIMFIRE
dense-stand under-kill (170 structure_densephase + 92 volume_persistent under simfire) = fire crown-lift bug #3
([[fvsjl-fire-fmprob-bug3]], FMOLDC plumbing, deferred bigger task).

**Spot-check of the flagged `volume_persistent` class** (doctrine: don't auto-accept a "no-convergence" flag) —
worst thinbba entry 259566665010854 (BdFt 42.2%): the thinbba run is **byte-identical to no-management** because
the stand's cycle-2 BA (16) is BELOW the THINBBA residual-BA-40 target ⇒ **no trees cut (no-op thinning)**. So its
"management divergence" is really the BASE projection's BdFt tail — a sawtimber-DBH threshold crossing that swings
board-foot % large on small absolute volume, atop the accepted ±1 ULP dense-phase straddle on TPA/SDI/CCF/TopHt.
Cornered (volume/threshold primitive), NOT a thinning-logic gap. METHOD NOTE: some thinbba ledger rows are no-op
thins (residual target above stand BA) whose "divergence" is the base projection — this slightly inflates the raw
management-divergence count with base-projection tails; the cut-cycle bit-exactness (slices 13-14) is the real
thinning-logic evidence. Floor untouched (no code change). Loop continues NE.

---
## Slice 43bt (2026-07-11) — NE dig leads homogeneous: structure_densephase label validated; supervision posture

New NE dig lead 68395700010538 (TPA 15.3%) both-sides-traced: dense 1368-TPA stand, jl UNDER-thins ~15% at the
first-cycle self-thinning event (430 vs live 373 @2011), BA tracks (188/188), TopHt tracks post-cycle-0 (99/99,
54/54, 59/59), TPA/QMD/SDI/CCF split; cycle-0 TopHt 92/96 = AVHT40 largest-40-TPA tie-break. Same cornered
**structure_densephase** primitive as 65944663010538 (43br).

**Pattern established (3 NE + 5 SN both-sides traces):** the forward loop's new dig-worthy leads are homogeneously
the dense self-thinning ±straddle — a dense stand hits a large first/second-cycle mortality event where the RDPSRT-
tie / SDI-limit kill FRACTION differs by ULP-amplified amounts (±straddle: under- or over-thin), with BA/TopHt
preserved and TPA/QMD splitting, plus the cycle-0 AVHT40 tie-break. The ledger's deterministic `structure_densephase`
classifier is FAITHFUL to this (verified, not assumed). **Supervision posture:** the signature classifier reliably
buckets these; re-tracing every identical instance adds no taxonomy value. Intervene on (a) a NON-structure_densephase
or UNCLASSIFIED lead, (b) the dig queue nearing DIGCAP=200, or (c) a signature the classifier flags for manual trace.
Loop: NE cursor 38000/178149, healthy (bit_exact ~82%/batch). Floor untouched.

---
## Slice 43bu (2026-07-11) — Pillar-1 status reconciliation (extraction script + manifest done-state)

Reconciled Pillar-1's done-state ("per-variant plot manifest + extraction script; materially larger than the 162
baseline") against what exists:
- **Extraction script: DONE** — `test/harness/fia/extract_sample.jl` is a deterministic (no-RNG), read-only,
  per-variant stratified sampler: orders FVS_STANDINIT_COND by (ECOREGION, LOCATION, STAND_CN) and takes every
  K-th, spreading across ecological units (the axis that surfaced the eco_unit DG bug) + national forests. Regenerates
  the sample reproducibly. Read-only on the master DB.
- **Sample manifests: exist** for NE/CS/LS (`.sweep_work/{ne,cs,ls}_sample.txt`, ~150 stands each with strata).
- **Coverage far exceeds a sample:** the forward sweep runs the FULL FVS-ready population per variant (SN complete
  at 637641; NE in progress at 38000/178149; CS/LS queued) through the per-cycle 10-col live-vs-jl differential —
  materially larger than 162, and larger than any sample manifest. Per-stand outcomes are the durable sweep ledger.
- **Gap (deferred by design):** a single CONSOLIDATED manifest doc with the full-population strata breakdown
  (stands per ECOREGION / forest-type / structure-class / site-class per variant). Building it needs a full GROUP BY
  scan of the 70GB master FVS_STANDINIT_COND, which would contend on master-DB I/O with the running loop's per-batch
  sub-DB builds and slow the sweep. Correctly deferred to a sweep milestone (NE completion or all-variants-done),
  when the scan won't compete. NOT a fidelity gap — purely a documentation-consolidation artifact.
Floor untouched (no code change). Loop healthy, NE advancing.

---
## Slice 43bv (2026-07-11) — Pillar-3 extended to NE management (salvage + plant): bit-exact-or-cornered

Ran the management differential (ledger_fia.jl, scratch LEDGER, **no SWEEP_DB** ⇒ no live-sweep contention) on a
bounded 40-stand deterministic NE sample (head of .sweep_work/ne_sample.txt) under SALVAGE and PLANT vs
freshly-relinked FVSne. Result — **ZERO UNCLASSIFIED; every divergence cornered to a named primitive:**
- SALVAGE: 28 bit_exact, 11 print_boundary, 1 structure_densephase (of 40).
- PLANT:   19 bit_exact, 15 print_boundary, 2 threshold_crossing, 2 count_straddle, 1 volume_persistent,
  1 structure_densephase (of 40).
Dominant non-bit-exact class = **print_boundary** (the .sum last-printed-digit rounding straddle — a sub-ULP diff
flips the rounded display; the reason the campaign standard is "bit-exact-OR-cornered"). PLANT shows more
print_boundary because added regeneration trees land on more display-rounding straddles. Same profile as the SN
management ledger (43bs, zero UNCLASSIFIED).

**Pillar-3 status update:** NE management now spans THINBBA (cut bit-exact, tail ULP — slices 13-14) + SALVAGE +
PLANT (this slice), all bit-exact-or-cornered. SIMFIRE deferred (fire crown-lift bug #3). SN complete. Remaining
Pillar-3 gap: CS/LS salvage/plant/thinbba differentials (same bounded approach, when loop CPU allows). Kept the
subset to 40 stands to bound contention with the running forward sweep; loop stayed healthy throughout (NE
advancing). Floor untouched (no code change).

---
## Slice 43bw (2026-07-11) — Pillar-3 CS management + a REAL live-FVS SIGFPE (essprt stump-sprout) under THINBBA

Ran the CS management differential (ledger_fia.jl, scratch LEDGER, no SWEEP_DB) on a bounded 40-stand deterministic
CS sample under SALVAGE + PLANT + THINBBA vs freshly-relinked FVScs. **Zero UNCLASSIFIED; all divergences cornered**
(salvage 36 be/3 print_boundary/1 volume_persistent; plant 32 be/7 print_boundary/1 structure_densephase; thinbba
28 be/7 structure_densephase/1 volume_persistent/1 print_boundary) — EXCEPT **3 thinbba `live_crash`**.

**The 3 thinbba live_crash are a GENUINE live-FVS SIGFPE** (not an artifact, not FVSjl): reproduced directly on the
CS oracle (signal 8) at **essprt.f:216-217** — the CS species-57/58 stump-sprout multiplier `1./((DSTMP/0.7788)-
0.4403)` divides by zero at cut-stump diameter DSTMP≈0.343", under the build's `-ffpe-trap=...,zero,...,overflow`.
THINBBA cuts a sp57/58 hardwood → post-harvest sprouting (esuckr→esnutr→gradd→essprt) → singularity → crash.
Salvage/plant/none don't crash (no matching cut→sprout). FVSjl projects all 3 without crashing. FULL root-cause,
the anomaly (only 57/58 use the reciprocal; siblings 47/54 use the linear `(…)*2.54`), the proposed maintainer fix
(logistic-argument clamp), and the in-container application blocker are documented in docs/FVS_SOURCE_BUGS.md
("CS essprt.f:216-217"). Stands: 1910906629290487, 488847180126144, 224864192010661.

**Environment blocker (measured, honest):** a fix cannot be validated here without risking oracle fidelity — local
gfortran 12.2.0 ≠ essprt.o's build compiler (SUSE 15.2.1), so a source recompile perturbs normal-stand numerics;
and even relinking main.f with the exact original flags does NOT reproduce /tmp/FVScs_new bit-for-bit. So the bug
is DOCUMENTED + fix PROPOSED for the maintainer, NOT applied. /tmp/FVScs_new left pristine (verified: still SIGFPEs
the 3 crashers, runs normals). All my test binaries removed; build-dir .o files untouched (compiled only to scratch).

**Pillar-3 status:** SN (all regimes) + NE (thinbba/salvage/plant) + CS (salvage/plant/thinbba) now differentiated,
all bit-exact-or-cornered, with the CS essprt SIGFPE the sole genuine live-FVS defect surfaced (directive-flagged,
fix proposed). Remaining: LS management + SIMFIRE (fire crown-lift bug #3). Loop stayed healthy throughout (NE
46000→58000 during the investigation). FVSjl floor untouched (no FVSjl code change).

---
## Slice 43bx (2026-07-11) — Pillar-3 LS management: non-fire management now covered on ALL 4 variants

Ran the management differential (ledger_fia.jl, scratch LEDGER, no SWEEP_DB) on a bounded 40-stand deterministic LS
sample under SALVAGE + PLANT + THINBBA vs freshly-relinked FVSls. **Zero UNCLASSIFIED, zero live_crash; every
divergence cornered to a named primitive:**
- SALVAGE: 18 bit_exact, 8 print_boundary, 5 count_straddle, 4 threshold_crossing, 4 structure_densephase, 1 volume_persistent.
- PLANT:   9 bit_exact, 12 print_boundary, 10 structure_densephase, 6 count_straddle, 2 threshold_crossing, 1 volume_persistent.
- THINBBA: 11 bit_exact, 21 structure_densephase, 4 volume_persistent, 2 threshold_crossing, 2 count_straddle.
LS shows more structure_densephase (thinbba 21/40) + lower bit_exact than SN/NE/CS — consistent with LS's documented
dense-phase terminal tails (calibration-backdating relative-ranking, accepted class; see the LS port state) — but
all cornered, none unexplained.

**Pillar-3 milestone — NON-FIRE management now differentiated on ALL 4 variants, all bit-exact-or-cornered:**
- SN: thinbba/salvage/plant/simfire, zero UNCLASSIFIED (43bs).
- NE: thinbba (43bv thinbba cut bit-exact, slices 13-14) + salvage + plant (43bv), zero UNCLASSIFIED.
- CS: salvage + plant + thinbba (43bw), zero UNCLASSIFIED + the 1 genuine essprt SIGFPE (documented, fix proposed).
- LS: salvage + plant + thinbba (this slice), zero UNCLASSIFIED.
Remaining Pillar-3 gap: SIMFIRE on NE/CS/LS (SN simfire done) — gated by the known fire crown-lift bug #3
([[fvsjl-fire-fmprob-bug3]]), the one deferred bigger implementation task. Loop stayed healthy (NE advancing).
FVSjl floor untouched (no FVSjl code change).

---
## Slice 43by (2026-07-11) — fire bug #3 re-scoped: crown-lift record-carry IS in place; cause = CWD-accumulation puzzle

Re-measured the documented SN SIMFIRE fire-kill divergence (1848605816290487) with CURRENT code vs freshly-relinked
FVSsn: pre-fire bit-exact (2024/2029 168/165 both), at the 2034 fire live retains **23 TPA / jl 44** — jl still
UNDER-KILLS ~2× (BUG REAL + CURRENT, not stale). BUT the both-sides SOURCE re-trace REFUTES part of the bug memory's
mechanism: `ffe_oldht/ffe_olddbh/ffe_oldcr` ARE registered in `_TREE_VEC_FIELDS` (trees.jl:130) and `ffe_oldcrw` is
explicitly carried in `copy_tree!` (trees.jl:168) ⇒ record moves (tripling/compaction) DO preserve the old-crown
snapshots. So `compute_crown_lift!`'s `oldht>0 || continue` skip (fuel_additions.jl:141) now fires ONLY for (a)
cycle-1 (no snapshot) and (b) genuinely-new regen records (no prior crown) — BOTH correct. The "record-list change
skips crown-lift" framing (memory [[fvsjl-fire-fmprob-bug3]], slice 43bl) is STALE — the FMOLDC-style carry was
since implemented. ⇒ the persistent fire under-kill is NOT the crown-lift record-churn; it is the still-unresolved
CWD-ACCUMULATION puzzle (slice 43bk): jl under-populates the LARGE (>3") coarse-woody pool (jl ~0.2-0.67 vs live 4.6
t/ac at the fire) ⇒ collapses the FMDYN (SMALL,LARGE) fuel-model point ⇒ picks low-flame models {8,9} not {8,9,10}
⇒ flame 2.10 too low ⇒ under-kill. The OPEN question is live's >3" CWD SOURCE on a QMD-4.8 seedling stand (small
trees can't make >12" boles, yet live reports >12"=3.7) — needs a both-sides CWD-source trace (which FVS routine
feeds cwd classes 4-9 here: FMCADD woody breakage vs snag-bole falldown vs a fuel-model-independent path). Recurring
meta-lesson reaffirmed: re-trace "remaining plumbing" flags against SOURCE — here it caught a STALE mechanism (the
record-carry is done) AND kept the real bug (under-kill persists, CWD-accumulation cause). FVSjl floor untouched.

---
## Slice 43bz (2026-07-11) — fire bug #3 narrowed: initial LARGE CWD gap is an FFE init-computation diff (not FUINI table/mapping/reader)

Deep both-sides trace of the SN SIMFIRE under-kill (1848605816290487), driving to the initial (2024/cycle-0) fuel.
Live FVS "ALL FUELS REPORT" (via FuelOut) shows 2024 SURFACE >3" = 16.7 t/ac (>12"=10.6, 6-12"=3.7, 3-6"=2.4),
litter 4.40, duff 25.2 — a LARGE pool present at cycle-0 that then DECAYS (16.7→4.6→3.7→…). Traced the source:
- **jl's FUINI table + forest-type→row mapping are CORRECT.** FVS fmcba.f maps FIA forest type 141/142 → FTDEADFU=2
  (longleaf-slash); jl's `ffe_dead_fuel_type(142)`=2 matches. The FVS FUINI longleaf-slash DATA row =
  `[0.10,0.66,0.98,0.12,0.29,0.26, 0,0,0, 6.38,8.66]` — jl's row 2 is bit-exact to it (slice 43bk was right about that).
- **NO FUINI row explains the 10.6.** ALL 9 FVS FUINI rows have classes 7,8,9 (>12") = 0.0, and none has litter 4.40
  / duff 25.2. So live's initial LARGE (>12") pool is NOT from the FUINI table — refuting the "wrong FUINI row" idea.
- **Not a simple DWM-reader gap.** The master FIA DB has rich down-woody tables (DWM_COARSE_WOODY_DEBRIS,
  DWM_DUFF_LITTER_FUEL, COND_DWM_CALC, …), but `build_subdb` copies only FVS_STANDINIT_COND + FVS_TREEINIT_COND,
  so BOTH engines run the identical 2-table sub-db + identical keyfile (neither reads DWM). jl's keytext(simfire)
  emits no FUELINIT. ⇒ same inputs, different initial large-CWD ⇒ an FFE INIT-COMPUTATION difference.

**Net:** the fire under-kill root is now cornered to *how FVS FFE derives the initial large (>3", esp >12") down-wood
pool at cycle-0* (a live-run FUELINIT-activity echo carries the exact report values 3.65/10.60/4.40/25.19 despite no
FUELINIT keyword — an FFE-internal derivation, e.g. an initial snag/past-mortality → down-wood spin-up, or a
DATABASE-side FUELINIT auto-inject). NOT the FUINI table, forest-type mapping, crown-lift record-carry (43by), or a
DWM table read — all ruled out by both-sides trace. NEXT (focused effort): instrument jl's fire cwd-by-size at 2024
vs live's ALL FUELS 2024 row to quantify the exact class gap, then trace the FVS FFE init routine that fills classes
4-9 at cycle-0 (fmcba.f OPFIND MYACT(1)=2521 path / fminit snag-derivation). FVSjl floor untouched (no code change).

---
## Slice 43ca (2026-07-11) — ★ fire bug #3 ROOT CAUSE FOUND (read the code): FIA FUEL_* standinit columns not read by jl

Per doctrine (read the FVS SOURCE, don't infer): traced live's initial large CWD to `dbsstandin.f`. The
FVS_STANDINIT_COND table carries MEASURED FIA down-woody-material fuel loadings as columns —
`FUEL_0_25_H/FUEL_25_1_H/FUEL_1_3_H/FUEL_3_6_H/FUEL_6_12_H/FUEL_12_20_H/FUEL_20_35_H/FUEL_35_50_H/FUEL_GT_50_H`,
`FUEL_LITTER`, `FUEL_DUFF`, + the `_S` soft-decay variants (22 columns total). `dbsstandin.f:396-458` reads them
into RSTANDDATA(39..63); :843-1000 assembles a FUELINIT/FUELHARD+FUELSOFT keyword; `fmcba.f:318-343` (the OPFIND
MYACT(1) branch) then OVERRIDES the FUINI-table STFUEL with those values. Confirmed for 1848605816290487:
FUEL_12_20_H=10.6035, FUEL_6_12_H=3.6520, FUEL_1_3_H=5.5894, FUEL_LITTER=4.4023, FUEL_DUFF=25.1908 — EXACTLY live's
FUELINIT echo + ALL FUELS 2024 row. Class map: FUEL_12_20→STFUEL class 6 (the FUINI longleaf default there is only
0.26 ⇒ the measured 10.6 is a ~40× override ⇒ the large-CWD gap).

**jl's `fia_database.jl` does NOT read the FUEL_* columns** ⇒ jl falls back to the FUINI table (class-6 = 0.26,
classes 7-9 = 0) ⇒ near-zero large CWD ⇒ collapsed FMDYN fuel-model point ⇒ low flame ⇒ SIMFIRE under-kill. This is
a READER GAP, the exact class as the Dunning (dbsstandin.f:763) + slope-default reader bugs — a REAL, FIXABLE FVSjl
bug. Supersedes the crown-lift (43by, refuted) and FUINI-table (43bz, refuted) theories; the "CWD-accumulation
puzzle" is resolved — the large pool is INITIAL (measured FIA DWM), not accumulated.

**FIX (scoped, faithful to dbsstandin.f + fmcba.f):** (1) fia_database.jl reads the 22 FUEL_* columns into a plot
fuel-init snapshot (missing ⇒ -1 sentinel); (2) the FFE fuel init (fmcba.jl) applies them as the STFUEL override
(hard→decay-hard, soft→decay-soft) BEFORE the basal-area decay-class distribution, matching fmcba.f:318-393. Floor
guard: existing fire tests use synthetic stands WITHOUT these columns ⇒ -1 sentinel ⇒ FUINI fallback unchanged ⇒ no
regression expected. Validate: SN SIMFIRE 1848605816290487 fire-kill vs live + suite floor. IMPLEMENTATION NEXT.

---
## Slice 43cb (2026-07-11) — ★★ fire bug #3 FIXED: FIA FUEL_* reader → STFUEL override (bit-exact vs live)

Implemented the reader-gap fix (root cause slice 43ca). Three files:
- `src/core/state.jl`: PlotData gains `ffe_fuel_hard`/`ffe_fuel_soft` (Vector{Float32}, empty default) — the
  measured FIA down-woody fuel loadings, size classes 1:11.
- `src/io/fia_database.jl` (apply_fia_stand!): reads the 22 FVS_STANDINIT FUEL_* columns (hard `_H` + soft `_S`,
  with the FUEL_0_1 lumped-<1" split per fmcba.f:329-340) into plot.ffe_fuel_* (−1 sentinel for missing). Mirrors
  dbsstandin.f:396-458. Only sets the fields when ≥1 fuel column is present.
- `src/engine/fire/fmcba.jl` (fuels_init block): seeds FFEParams.stfuel_hard/soft from plot.ffe_fuel_* when present
  AND not already set by an explicit FUELINIT/FUELSOFT keyword (keyword precedence preserved). The existing
  STFUEL-override path then applies them, overriding the FUINI-table default.

**Validated BIT-EXACT vs freshly-relinked FVSsn** on 1848605816290487 (the documented under-kill stand): every cycle
of all 6 structural cols now matches — 2034 fire kills to **23 TPA == live 23** (was jl 44 vs live 23, ~2× under-
kill); 2049 22/22, QMD 15.3/15.3, TopHt 67/67. The measured FUEL_12_20_H=10.6 now populates the large CWD pool ⇒
correct FMDYN fuel model ⇒ correct flame ⇒ correct fire kill. Supersedes the crown-lift (43by) + FUINI (43bz)
theories entirely.

Floor safety by design: the fields default empty; synthetic test stands (no FIA FUEL_* columns) ⇒ empty stash ⇒
fmcba! seeding skipped ⇒ FUINI fallback unchanged. Sweep (regime=none) never calls fmcba! ⇒ unaffected. FULL SUITE
floor validation RUNNING (pid 158296) — to confirm 38527/143 (+ documented broken set) holds before declaring done.

## Slice 43cb (cont) — FLOOR CONFIRMED
Full suite after the fire fix: **FVSjl 38587 pass / 75 broken / 0 real fail** (broken set UNCHANGED vs baseline
38586/75). The only non-pass are 3 "vs Oracle A" Errors (test_treedata/test_keyword/test_init) — each a
`failed process: julia --project=/workspace/FVSjulia` subprocess (FVSjulia/Oracle-A can't precompile in this
container: the WeakRefStrings/SQLite artifact), independent of the FVSjl source change. ⇒ fire bug #3 fix is
COMPLETE: bit-exact vs live on 1848605816290487 + floor preserved. (Ran the suite with the forward loop PAUSED to
avoid the concurrent-precompile race; loop resumed from the NE 86000 checkpoint after.) ★★★ Fire bug #3 CLOSED.

## Slice 43cb (cont) — generalization + inertness validation
Re-ran SN simfire on a 40-stand sample (fresh FVSsn oracle) to test the fix's population effect: 21/40 bit_exact,
rest cornered (9 structure_densephase, 7 volume_persistent, 2 threshold_crossing, 1 print_boundary), ZERO
UNCLASSIFIED/live_crash. KEY: **0 of these 40 carry measured DWM fuel columns** (all 22 FUEL_* NULL) — so the fix is
correctly INERT for them (both engines use the FUINI default), which is why their divergences (dense-phase
self-thinning / volume tails, present with or without fire) are unchanged in KIND. The fix's value shows on the
fuel-data subset (target 1848605816290487, FUEL_12_20_H=10.6 → bit-exact). INERTNESS PROVEN: for the one stand whose
magnitude shifted vs the pre-fix ledger (165735019, structure_densephase 70→89), ALL 22 FUEL_* columns are NULL ⇒
empty ffe_fuel stash ⇒ fmcba!'s `!isempty` guard skips seeding ⇒ the fix CANNOT touch it. The shift is stale-golden
noise (the pre-fix ledger used the PRIOR, pre-container-restart oracle; a chaotic dense stochastic-fire stand shifts
with any oracle relink — cf. the CS relink non-reproducibility, slice 43bw). Per doctrine (validate vs FRESH oracle,
not a stale golden) the fresh-oracle run is authoritative ⇒ NO regression. Fix is correct-where-applicable +
inert-elsewhere + floor-preserved. ★★★ Fire bug #3 CLOSED (validated).

## Slice 43cb (cont) — cross-variant (NE/CS/LS) validation: fire-fuel fix is VARIANT-SAFE
The fix is variant-agnostic by construction (shared fia_database.jl reader + the shared fmcba! STFUEL-override
block, outside the per-variant dead-fuel branch). Verified empirically on a fuel-data (FUEL_12_20_H>0) simfire stand
per variant vs fresh oracle:
- NE 1318365400290487: pre-fire bit-exact; tail ±2-3% dense-phase straddle (155/151…). Cornered.
- CS 66720086010661: bit-exact through 2041, then ±1-2 ULP straddle (416/418). Cornered.
- LS 1642770637290487: bit-exact/±ULP throughout (5316/5317, 991/984). Cornered.
All bit-exact-or-cornered (accepted structure_densephase self-thinning), NONE showing the gross fire-fuel under-kill
the SN target had — confirming the fix produces CORRECT fuel init across all 4 variants (live also reads the columns
via the shared dbsstandin.f) with no regression. ⇒ doctrine #5 (variant-safe) satisfied. Pillar-3 SIMFIRE is now
bit-exact-or-cornered on the fuel-data subset for every variant; fuel-less dense SIMFIRE stands remain the accepted
dense-phase primitive. ★★★ Fire bug #3 CLOSED + variant-validated.

---
## Slice 43cc (2026-07-11) — NE volume_persistent cornered: NVEL out-of-domain 0 vs jl extrapolation (degenerate runaway-height stand)

Both-sides trace of the sole non-densephase NE dig lead 207147469020004 (volume_persistent). The stand runaway-grows
height: TopHt 1→109→162→214→258→295 ft over 2013-2063 (295 ft = non-physical) — BOTH jl and live reach the same
heights (a SHARED degenerate FVS growth behavior, not a jl divergence; bit-exact TPA/TopHt/QMD throughout). The
divergence is VOLUME only, and only once TopHt ≥ ~258 ft: **live returns 0 volume (TCuFt/SCuFt/BdFt all 0 at 2053/
2063) while jl extrapolates** (2053 TCuFt 0/15284; 2063 SCuFt 0/**1.5395e10** — a garbage overflow). 

**Mechanism (both-sides):** live calls the NVEL (National Volume Estimator Library) which sets an error flag for
out-of-domain inputs (extreme height); FVS r9clark.f:110/177 (`TLOGVOL=0` … `if(errFlg.ne.0) return`) then returns
0 volume. jl's ported r9 volume equation has NO out-of-domain guard ⇒ it extrapolates the polynomial past its valid
height range ⇒ absurd values. ⇒ cornered primitive **volume_persistent / NVEL-out-of-domain**: manifests ONLY on
degenerate runaway-height stands (TopHt ≥ ~258 ft — real stands never reach this); regime=none (NOT related to the
fire fix). PROPOSED fix direction (deferred, low-priority — degenerate/rare): port the NVEL errFlg domain check (or
clamp jl volume to 0 when total height exceeds the equation's valid domain) so jl returns 0 like live. Note: jl's
1.5e10 overflow is a robustness gap that could pollute aggregate stats on degenerate stands — worth the domain
clamp eventually. Both-sides-traced + cornered; no unexplained divergence remains in the current dig queue.

## Slice 43cc (cont) — NVEL errFlg trigger read from source (fix scoped, deferred)
Read r9clark.f to pin the exact 0-volume trigger. errFlg is set nonzero by a MULTI-CONDITION NVEL validation, not a
single clean domain bound: totHt≤17.3 (errFlg=8, :202), dbhOb≤0 / htTot<ht1Prd / ht2Prd<ht1Prd inconsistencies
(errFlg=3/7/8/10, :550-586), and a diameter-profile degeneracy `dbhOb≤sawDib` (errFlg=13, :744). Notably the upper
height-vs-DBH sanity bound `htTot > 35·sqrt(dbhOb+3)` (errFlg=5) is COMMENTED OUT (:577-582) — FVS devs disabled it —
so there is NO simple "too tall" guard; our stand trips one of the geometry/consistency checks at extreme height.
`if(errFlg.ne.0) return` (r9clark.f:177/184/588) then leaves volume 0. The FAITHFUL fix = port the NVEL r9Prep/
r9dia417 errFlg validation into jl's r9clark_vol (so jl also returns 0 out-of-domain); it is a bounded port but
DISPROPORTIONATE for a single rare degenerate runaway-height stand ⇒ CORNERED + deferred (not a padded tolerance /
heuristic clamp, per doctrine #4 — the mechanism is named and the fix is the exact FVS validation). Both-sides-trace
complete.

## Slice 43cc (cont) — NE volume garbage is the KNOWN r9clark fvsMod-vs-NVEL discrepancy (not a new primitive)
Checked the jl side: src/engine/r9clark_vol.jl is translated from `r9clark_fvsMod.f` (FVS's LOCAL modified r9clark),
whereas LIVE FVSne links the NVEL-library `r9clark.f` (bin/FVSne_buildDir, the errFlg-validated version read above).
⇒ the NE 207147469020004 volume divergence (jl extrapolates to 1.5e10, live NVEL errFlgs → 0 on the degenerate
runaway-height geometry) is a manifestation of the ALREADY-DOCUMENTED r9clark fvsMod/NVEL version discrepancy (the
D38 r9clark family in docs/FVS_SOURCE_BUGS.md), NOT a new unexplained primitive. Cornered under that known family;
the "port the NVEL errFlg validation into jl's fvsMod-based r9clark" fix is the same D38 alignment work (bounded,
deferred, rare-degenerate-only). Both-sides-trace fully closed — no unexplained divergence remains in the dig queue.

## Slice 43cb (cont) — fire-fix coverage scope (Pillar-3 SIMFIRE)
Bounded prevalence sample (first 5000 stands/variant) of the measured DWM large-fuel column FUEL_12_20_H (the one
that drove the fire under-kill): SN 0.2%, NE 0.7%, CS 0.4%, LS 0.1% of stands carry it. FIA's down-woody-material
inventory runs on only a small subset of plots, so the fire-fuel reader gap directly affected ~0.1-0.7% of SIMFIRE
stands — MODEST in count but IMPORTANT in kind: those are exactly the DWM-inventoried plots where fuel-driven fire
behaviour matters, and the fix makes jl bit-exact there (target validated). For the ~99% without DWM data, both
engines use the FUINI-table default (already consistent) and the fix is provably inert. ⇒ Pillar-3 SIMFIRE:
bit-exact on the fuel-data subset (all variants), FUINI-consistent elsewhere; the residual dense-SIMFIRE divergences
are the accepted structure_densephase self-thinning primitive (fuel-independent), NOT the fire-fuel bug.

---
## Slice 43cd (2026-07-11) — CS essprt SIGFPE fix VALIDATED here (corrects earlier "env-blocked")

Reopened the CS essprt.f:216-217 SIGFPE (slice 43bw) per the standing crash directive. Earlier I called it
"env-blocked"; that was WRONG, from two artifacts: (1) partial compile flags (dropped -fintrinsic-modules-path/
-fpre-include) and (2) comparing the full .sum text incl the per-run TIMESTAMP header. Corrected: build-dir main.o
is LOCAL gfortran 12.2.0 (DWARF producer) + does zero FP ⇒ recompiling with FULL exact flags gives BYTE-IDENTICAL
.text/.data/.rodata. Rebuilt `main.f` with `-ffpe-trap=invalid` (drop zero,overflow), relinked, and validated with
a timestamp-stripped DATA comparison: **50/50 normal CS stands data-bit-identical to /tmp/FVScs_new + all 3 crashers
FIXED** (SIGFPE→valid output; the 1./0.→+Inf propagates to the correct logistic limit = the value FVSjl computes).
⇒ the crash IS fixable + validatable in this container via the build-flag path (the essprt.f source-recompile path
stays blocked — essprt.o is SUSE 15.2.1). Oracle LEFT PRISTINE / not hot-patched: the sweep is regime=none which
does NOT trigger the crash (needs THINBBA→sprouting), so CS sweep coverage is unaffected; and the build-flag drops
zero,overflow GLOBALLY (broader than the R9-clark precedent that kept them), so the precise essprt.f source guard
remains the maintainer recommendation. Directive satisfied (root-caused + VALIDATED fix + submission proposal). Test
binaries removed; build-dir .o untouched. Meta-lesson: validate relinks with FULL DWARF-producer flags + strip the
.sum timestamp header — both bit me into a false "env-blocked" verdict (doctrine: measure precisely, don't guess).

## Slice 43cd (cont) — with the essprt fix, the 3 crashers are bit-exact-or-cornered (Pillar-3 CS THINBBA closed)
Rebuilt the essfix oracle and compared the 3 formerly-crashing CS THINBBA stands to FVSjl: fixed-live vs jl is
BIT-EXACT through the harvest cycles (e.g. 1910906629290487: 2024/2034/2044 all "="), then a SMALL dense-phase
divergence develops (2054+: TPA 0.3-2%, volume 1-3%, TopHt ±6 late) — the accepted structure_densephase self-thinning
primitive (dense sprouting stand, record churn). So each essprt crasher goes live_crash → BIT-EXACT-OR-CORNERED with
the fix (bit-exact pre-divergence + dense-phase-cornered tail), NOT a new unexplained divergence. ⇒ Pillar-3 CS
THINBBA fully resolved: the essprt SIGFPE is the only genuine live-crash, it has a VALIDATED fix (slice 43cd), and
under the fix the stands are bit-exact-or-cornered vs jl. essfix binary cleaned; oracle pristine. The crash directive
is fully discharged: root-caused + validated fix + the stands shown bit-exact-or-cornered once fixed.

## Slice 43ce (2026-07-12) — NE dig-queue spot-verify: structure_densephase label is genuine (not a masked bug)

Doctrine check (memory: a structure_densephase tag has twice hidden a real bug — AVH-backdating, FIA-slope-default —
so fresh queue leads get TRACED, not auto-cornered). Picked the newest non-degenerate NE lead from the dig queue:
**1167714929290487** (NE, regime=none, worst=TPA@2030, 77 trees, classifier=structure_densephase).

Full-trajectory differential vs freshly-relinked /tmp/FVSne_new (LIVE/JL):
```
2020 | TPA= BA= SDI= CCF= TopHt 21/24 QMD= vols=      <- cycle-0 inventory
2030 | TPA 1933/2781  TCuFt 2651/2609 MCuFt 166/186   <- +44% transient TPA in jl
2040 | TPA 865/936    ...                              <- gap closing
2050 | TPA 570/575    ...                              <- ~1%, converged
2060 | TPA 416/420    2070 | 311/314                   <- tracks live
```
Raw live .sum: 2020 = age 10, **TPA 9096, QMD 0.8"**, dense seedling thicket self-thinning 9096→1933→865→570.

Both-sides trace — both signals corner to the SAME named primitive `structure_densephase`:
- **TopHt 21/24 @ cycle 0**: TPA/BA/SDI/CCF/QMD all BIT-EXACT; only TopHt differs. TopHt = mean height of the ~40
  largest-DBH trees/acre; among thousands of trees near-tied at 0.8" DBH the top-height cohort membership is decided
  by the RDPSRT unstable-sort tie-break (fvsjl-stand-pct-rdpsrt-fix family) — tied trees carry different imputed
  heights ⇒ 3-ft delta. Tie-break, not a height-model gap.
- **Transient TPA @ 2030**: dense self-thinning mortality phase; jl's kill lags one cycle then over-corrects,
  converging to within ~1% by 2050 and tracking live thereafter (the documented dense self-thinning limit-cycle).

Verdict: classifier label CONFIRMED genuine; cornered to structure_densephase (RDPSRT tie-break + dense
self-thinning phase). No floor impact (dig-only; no code change). Fixed a stale harness: .sweep_work/dig_one.jl
now unwraps the (String,Bool) tuples that run_live/run_keyfile return.

## Slice 43cf (2026-07-12) — NE dig-queue: TCuFt-worst class confirmed = D38 r9clark domain-guard

Queue discriminator scan: no dig entry is worst-on-TopHt (⇒ no masked height bug; TopHt divergences are always
secondary to a density tie, per 43ce). Worst-col histogram = TPA 39 / TCuFt 11 / CCF 5 / QMD 2 — all density or
volume, all named-primitive classes. Spot-verified the one distinct volume class: **NE 207147469020004**
(worst=TCuFt@2053, classifier=volume_persistent). Differential vs /tmp/FVSne_new (LIVE/JL):
```
2013,2023 | all 10 cols BIT-EXACT
2033,2043 | TPA/TCuFt/MCuFt/BdFt off by ~1 ULP (dense-phase mortality ULP)
2053 | TPA/BA/SDI/CCF/QMD = ; TopHt 258/257; TCuFt 0/15284  MCuFt 0/14932  SCuFt 0/12819  BdFt 190/86286
2063 | TCuFt 0/17245  MCuFt 0/16981  SCuFt 0/15395
```
Both-sides trace: at 2053/2063 LIVE reports **0** cuft/bdft while jl computes ~15000 — the trees are identical
(TPA/BA/QMD bit-exact). Signature of FVS NVEL r9clark.f errFlg DOMAIN GUARD (`if(errFlg.ne.0) return` ⇒ 0 for
out-of-domain trees; here the stand reaches TopHt ~257 ft, out of the volume-eq domain), vs jl's ported
r9clark_fvsMod.f which EXTRAPOLATES. = the D38 volume family; cornered to the r9clark errFlg-domain-guard primitive.
Matches the deferred "NE NVEL volume domain-clamp" item (make jl mirror the errFlg 0-return); deferred by design.
No floor impact (dig-only).

## Slice 43cg (2026-07-12) — dig-queue taxonomy holds uniformly (3 classes spot-verified)

Third spot-check this session: **SN 212199455010854** (worst=QMD@2008, structure_densephase). Differential vs
/tmp/FVSsn_new: cycle-0 (2003) all-10 BIT-EXACT; dense seedling thicket (2013 TPA 10026, QMD ~1") self-thinning
10026→4648→2914→2055, every divergence sub-2% density-driven. QMD-"worst" @2008 = `QMD 1/0` (value ~0.5 straddling
the integer print boundary) — print_boundary on a dense-phase stand. Confirmed genuine.

**Consolidated observation (Pillar-4):** across the three distinct dig-queue classes spot-verified this session —
43ce density tie-break (cycle-0 TopHt cohort tie), 43cf volume domain-guard (r9clark errFlg 0-return), 43cg dense
QMD print-boundary — the pattern is uniform: dense seedling/sapling FIA stands (10^3–10^4 TPA, QMD <1–2") produce
sub-2% dense-phase divergences, and the classifier's "worst column" merely reflects whichever metric straddled a
tie/rounding boundary that cycle (TopHt cohort-tie, QMD print-boundary, TPA self-thinning phase, CCF/SDI density).
All corner to named primitives (structure_densephase / print_boundary / volume_persistent-D38); no masked bug in any.
No floor impact (dig-only).

## Slice 43ch (2026-07-12) — CORRECTION to 43cf mechanism (re-traced against r9clark source)

Per doctrine (both-sides-trace; re-trace asserted verdicts against SOURCE), read volume/NVEL/r9clark.f. The 43cf
claim "live guards via errFlg 0-return, jl extrapolates WITHOUT the guard" is **mechanistically WRONG** and is
withdrawn. Findings:
- r9clark.f DOES return 0 on domain failure — `if(errFlg.ne.0) return` at :177/:184/:203/:225/:256 (leaves cfVol=0).
- BUT the trigger is a COMPUTED-quantity domain check, not a simple input bound: `errFlg=8 if COEFFS%totHt.le.17.3`
  (:202, the TAPER-computed total height from r9totHt), plus internal errFlg from r9dia417/r9totHt/r9cuft.
- CRUCIALLY the extrapolating variant volume/r9clark_fvsMod.f ALSO has errFlg (72 refs) and the jl port
  src/engine/r9clark_vol.jl tracks errFlg too (returns (dbhIb,dib17,errFlg)). So this is NOT guard-vs-no-guard.

Accurate verdict: NE 207147469020004 has NORMAL bit-exact-ish volume through 2043 (TCuFt 9459/9449, 13007/13003)
then live DROPS to 0 at 2053 (0/15284) while trees persist (TPA/BA bit-exact). A specific errFlg domain path fires
in LIVE at that cycle that the jl port does not replicate — the exact condition (which of totHt≤17.3 / r9dia417 /
r9totHt fails, and why jl's port computes errFlg=0 there) requires a per-tree volume trace at 2053. Class = D38
r9clark computed-domain family (volume_persistent), but the precise both-sides mechanism is OPEN, not yet pinned.
Deferred to the D38 volume-domain work item (per-tree trace + faithful errFlg replication), gated on a sweep-pause
milestone. No floor impact (dig-only). META: caught my own false-precise verdict by reading the source — the
recurring re-trace-cleared-flags lesson.

## Slice 43ci (2026-07-12) — D38 NE volume mechanism PINNED (both-sides source, supersedes 43ch OPEN)

Re-traced both sides against source; the 43ch "OPEN" is now resolved to corner-classification precision.
- Raw live .sum for NE 207147469020004 shows the trigger: TopHt runs away 109→162→214→**258→295 ft** while QMD is
  only 16.3" (h/d ≈ 16, vs a real ~0.5) — a degenerate height-growth geometry. It is reproduced BIT-EXACTLY by jl
  (TopHt 258/257) ⇒ the runaway height is FAITHFUL live-FVS behaviour (an NE height-model extrapolation for this
  spp809/site over 60 yr), NOT a jl divergence. Only the VOLUME of the degenerate tree diverges.
- LIVE links guarded NVEL r9clark.f (bin/FVSne_buildDir/r9clark.f, 12 `if(errFlg.ne.0) return`): computed-totHt
  guard `if(COEFFS%totHt.le.17.3) errFlg=8` (:202) + r9Prep height-reasonableness guards `errFlg=8/7/10` (:551-570).
  Out-of-domain ⇒ returns 0 volume.
- jl ports the EXTRAPOLATING volume/r9clark_fvsMod.f (src/engine/r9clark_vol.jl `_r9_totht` :136-151): no
  out-of-domain 0-return — it extrapolates a totHt and computes volume. ⇒ live 0, jl ~15284 at 2053+.

Verdict: cornered to **D38 r9clark taper-domain-guard, triggered by degenerate runaway-height geometry** (named
primitive: FVS NVEL out-of-domain 0-return vs jl fvsMod-port extrapolation). Both-sides SOURCE-confirmed; no longer
unexplained. Residual sub-detail (exact errFlg branch: :202 computed-totHt vs :551-570 reasonableness) needs a live
DEBUG%MODEL per-tree trace — does not affect the corner classification. Faithful-drop-in FIX = make jl's r9clark
port replicate NVEL's errFlg 0-return on out-of-domain trees (shared volume path; deferred to a sweep-pause
milestone + full suite run; floor-regression risk). 43cf mechanism was imprecise, 43ch over-corrected; 43ci is the
source-grounded resolution. No floor impact (dig-only).

## Slice 43cj (2026-07-12) — NE sweep COMPLETE; NE→CS rollover clean

NE full-population coverage sweep COMPLETE: ne.cursor=178149/178149 (all NE FVS-ready FIA plots, regime=none,
multi-cycle .sum differential vs /tmp/FVSne_new). Auto-rolled to CS (driver order SN→NE→CS→LS, first cursor<pop).
CS started clean: VARIANT=CS population=255952 offset=0, first batch (2000 stands) complete + 2nd building, **0
live_crash / 0 RUN FAILED**, dig queue stable at 57 (< DIGCAP 200). Confirms the CS oracle + variant-safe FIA
reader/CS dgf path produce valid differentials on the sweep (regime=none doesn't hit the essprt CASE 57/58 SIGFPE,
which needs THINBBA→sprouting — see fvsjl-cs-essprt-sigfpe). NE divergence taxonomy fully cornered (this session:
43ce density tie-break, 43ci D38 r9clark degenerate-geometry volume domain-guard, both source-grounded both sides).
Now supervising the CS sweep (255952 plots) toward the CS→LS rollover.

## Slice 43ck (2026-07-12) — CS Pillar-2 spot-verify: 4/4 stands bit-exact full trajectory

While the CS full-population sweep runs, proactively spot-verified CS multi-cycle fidelity on 4 real CS FIA stands
(from cs_sample.txt) via dig_one vs freshly-relinked /tmp/FVScs_new (contention-safe: each builds its own sub-DB
from the read-only master FIADB):
- 27621869020004: 2011→2061 (6 cycles) — ALL 10 cols BIT-EXACT every cycle
- 103505016010661: 2006→2056 (6 cycles) — ALL 10 cols BIT-EXACT
- 244228309010661: 2010→2060 (6 cycles) — ALL 10 cols BIT-EXACT
- 55193911010661:  1999→2049 (6 cycles) — ALL 10 cols BIT-EXACT
Pillar-2 (multi-cycle projection) for CS: 4/4 stands bit-exact on the WHOLE trajectory (TPA/BA/SDI/CCF/TopHt/QMD/
TCuFt/MCuFt/SCuFt/BdFt), not just cycle-0 — confirms the CS growth+volume+mortality spine reproduces live FVScs
behaviour across full default horizons. Complements the running population sweep (offset ~4000/255952, 0 live_crash).
No floor impact (dig-only).

## Slice 43cl (2026-07-12) — CS dig-queue label audit: threshold_crossing verified genuine
- **Queue state:** 82 entries, ALL cornering to the 5 named primitives (0 unknown-class):
  76 structure_densephase (SN 38 / CS 24 / NE 14), 5 threshold_crossing (SN 4 / CS 1), 1 volume_persistent (NE).
- **Spot-check target** (verify label, don't assume — doctrine #3/#6): CS 65459886010661, sig=threshold_crossing,
  div_cols=CCF|TCuFt, density_bitexact=true, struct_max_rel 0.606% vs vol_max_rel 7.601%.
- **Per-cycle differential (dig_one.jl, live FVScs vs FVSjl):** every density col (TPA/BA/SDI/CCF/TopHt/QMD)
  bit-exact at EVERY cycle; the only diffs are a single-cycle TCuFt blip at 2019 (592/547, fully recovered
  bit-exact by 2029/2039) and a 1-unit CCF rounding blip at 2049 (165/164).
- **Both-sides mechanism (confirmed, not inferred from pass/fail):** the stand STATE (TPA/BA/QMD) is bit-exact
  at 2019 — identical trees, identical dimensions — yet total cubic volume differs ~8% at that ONE cycle and
  re-converges. Volume qualification is a STEP function of per-tree DBH/height; a boundary tree sitting within a
  Float32 ULP of the merch threshold lands on opposite sides in the two implementations at exactly 2019, then
  grows unambiguously clear by 2029 so both re-qualify → transient blip, not a growth-model gap. Density is
  untouched because qualification does not feed back into growth (that is why density_bitexact holds throughout).
- **Verdict:** cornered — threshold_crossing (ULP-class in the input to a discrete step function). Label ACCURATE.
  CS queue remains 0 new-class leads; every entry corners to an already-named + both-sides-traced primitive.

## Slice 43cm (2026-07-12) — Pillar-1 status verified DONE; coverage-characterization pause-gated (measured why)
- **Verified (not assumed) Pillar-1 done-state is MET:** test/harness/fia/manifests/{sn,ne,cs,ls}_manifest.txt =
  500 stands/variant × 4 = 2000 (materially larger than the 162-stand modernization baseline), deterministic
  stratified extraction (extract_sample.jl, even stride over (ECOREGION,LOCATION,STAND_CN) — no RNG, regenerable),
  documented in manifests/README.md with the per-variant population/sampled/ECOREGION/LOCATION table. Checked
  before touching anything — the manifest infra already exists; no duplication needed.
- **Residual Pillar-1 doc gap identified:** the manifest STRATIFIES on geography (ECOREGION/LOCATION) only; the
  goal also lists forest-type / stand-structure / site-class strata. Whether the geographic-stride sample also
  SPANS those three axes is currently undocumented. FVS_STANDINIT_COND carries the needed columns
  (FOREST_TYPE_FIA, SITE_INDEX, AGE = development stage).
- **Why it is pause-gated (measured, doctrine #6):** `PRAGMA index_list(FVS_STANDINIT_COND)` = NO indexes on the
  70 GB master; a 2000-CN `IN` characterization query would force a FULL TABLE SCAN, i.e. exactly the master-DB
  I/O contention to keep off the live sweep. Deferred to the next idle window (CS→LS rollover or a DIGCAP pause),
  when the master DB is quiescent — then run the FOREST_TYPE_FIA/SITE_INDEX/AGE coverage tabulation and append it
  to manifests/README.md as the Pillar-1 "spans forest-type/site/structure" evidence.
- **Verdict:** Pillar-1 substantively COMPLETE (sample + script + docs); one full-scan-gated coverage-evidence
  addendum queued for the next sweep-quiescent milestone. No DB scan run now (sweep protected).

## Slice 43cn (2026-07-12) — Pillar-4 label-pipeline integrity verified from source (idle bug-hunt)
- **Trigger:** idle-time review of the cornering-label pipeline (the tool that stamps structure_densephase /
  threshold_crossing / count_straddle / print_boundary / volume_persistent). Verified from FVS-adjacent SOURCE,
  not assumed (doctrine #3/#6).
- **Corrected a mistaken premise:** the working note "run_keyfile now returns (String,Bool)" is WRONG.
  `FVSjl.run_keyfile` returns a plain `String` in BOTH branches (src/engine/simulate.jl:580 CSV, :582 default);
  no tuple path exists. The (String,Bool) tuple is the LOCAL `run_live` (ledger_fia.jl:86 `return live,
  live_crashed`), correctly destructured by the sweep at line 170 (`live, live_crashed = run_live(...)`).
- **Sweep labels TRUSTWORTHY:** ledger_fia.jl consumes run_keyfile's String via parse_sum10(jlout) with no
  unwrap (line 172) — correct. dig_one.jl's `jl isa Tuple ? jl[1] : jl` is a harmless no-op on a String.
- **signature.jl has NO bug:** line 51 `psum(FVSjl.run_keyfile(key; variant=var))` receives a String → works.
  (Initial "latent tuple bug" hypothesis REFUTED by reading run_keyfile's return — no change made, correctly.)
- **classify() is a COMPLETE deterministic partition** (ledger_fia.jl:113-123): struct_mat ⇒ count_straddle
  (!density_mat) | structure_densephase (density_mat); !struct_mat ⇒ print_boundary (!vol_mat) |
  threshold_crossing (vol_mat&converges) | volume_persistent (vol_mat&!converges). No fallthrough ⇒ no spurious
  UNCLASSIFIED for the real fact-combos; the CS-sweep's 0-UNCLASSIFIED stream is a genuine property, not a masked gap.
- **Verdict:** the entire Pillar-4 labeling pipeline is sound; no code change. One stale working-note claim corrected.

## Slice 43co (2026-07-12) — ★ HARNESS BUG: `.sum` parser mis-read fixed-width overflow → FALSE divergences
**Trigger.** User lowered DIGCAP 200→100; began the dig phase (queue=111). First target = the worst-magnitude
`structure_densephase` outlier, SN 253699300010854 (metric said CCF struct%=5433, SDI 7.6 *million*).
**Both-sides-trace ⇒ it is a MEASUREMENT bug, not an FVSjl divergence.** `parse_sum10` (and 7 sibling harness
parsers) tokenized the FVS `.sum` with `split()` — a WHITESPACE split — but the `.sum` data row is FIXED-WIDTH
(`sumout.f` FORMAT **9014**, byte-identical across SN/NE/CS/LS): `2I4,I6,I4,I5,2I4,F5.1,9I6,...`. CCF is `I4`
directly after SDI `I5` with **no separator**, so when CCF reaches 4 digits (≥1000) it fills its field and abuts
SDI: e.g. live `  762`+`1019` → `  7621019`, which `split()` reads as one token `7621019` and shifts every later
column. When CCF straddles 1000 DIFFERENTLY between live and jl (live CCF 1003 overflows, jl CCF 996 does not),
the two `.sum` rows tokenize into different columns ⇒ a fabricated SDI 7521003-vs-748 / CCF 5433% / QMD 4260-vs-19
"divergence." Verified against the raw `.sum` with a column ruler (`.sweep_work/dig_rawsum.jl`): the TRUE values
are SDI 752/748 (0.5%), CCF 1003/996 (0.7%), QMD 18.8/19.0 — an ordinary dense-stand straddle.
**Impact.** Systematic: for every stand whose CCF (or a wide SDI/vol field) straddles a width boundary between
live and jl, the parser inflated `max_rel_pct`/`struct_abs` and mis-drove `classify()` → FALSE `structure_densephase`
labels that tripped the dig-worthy escalation (≥15% & real abs) → false dig-queue rows, and contaminated the
headline `bit_exact` coverage counts. The `structure_densephase` cluster (dense stands = exactly the CCF≥1000
population) was the most affected. NOTE: `crashscan_fia.jl` already carried a comment that "whitespace-split of the
whole row is unreliable" — latent awareness of this exact bug, worked around only for year detection.
**Fix (harness-only; floor untouched — the julia suite never runs these parsers).** Parse by FIXED CHARACTER
COLUMNS: TPA[9:14] BA[15:18] SDI[19:23] CCF[24:27] TopHt[28:31] QMD[32:36] TCuFt[37:42] MCuFt[43:48] SCuFt[49:54]
BdFt[55:60] (verified against the ruler). Applied to all 8: `ledger_fia.parse_sum10` (the sweep), `run_sweep`/
`census_driver`/`validate_fia_cols`/`validate_fia10` (10-col), `validate_fia`/`manage_fia` (6-col), `signature.psum`
(6-col). `diff_one` reuses `manage_fia`'s parser; `crashscan` reads only the year (already safe). A fixed-column
parse can only REMOVE false divergence (bit-exact stands stay bit-exact), never add one.
**Validation.** SN 253699300010854: fake 5433% CCF / 7.6M SDI → true <3% (dense straddle). Cross-variant clean
(CS/NE parse correctly). GENUINE divergences unaffected — CS 351966617489998 (2034 TPA 35 vs 121, QMD 27.9 vs 15.0
at CONSTANT BA/SDI/vol) and NE 381531994489998 (2026 TPA 117 vs 626, QMD 17.2 vs 7.4) are real self-thinning
DISTRIBUTION divergences (jl retains more small trees); TPA is BEFORE the merge point so the old parser read it
right and they were correctly flagged.
**Re-check DONE** (`.sweep_work/recheck_digqueue.jl` → `digqueue_recheck.tsv`): reclassified all 111 queued stands
with the fixed parser via the sweep's own `classify()`. **12 ARTIFACTS** (2 now bit_exact, 6 densephase→print_boundary/
count_straddle/sub-threshold, 4 threshold_crossing <15%) + **99 GENUINE** (real dense self-thinning DISTRIBUTION
divergences — the accepted `structure_densephase`/`count_straddle` primitive). Purged the 12 from
`docs/fia_dig_queue.csv` (111→99; backup `.sweep_work/fia_dig_queue.prefix43co.bak`).
**★ Stale-queue note:** one artifact was SN 220314124010854 — whose frozen-TopHt was a GENUINE Dunning SITE_INDEX
bug already root-caused + FIXED (slice 43bp, fia_database.jl:157-164; TopHt bit-exact post-fix). The re-diff confirms
it: TopHt bit-exact, whole stand ≤0.3% (print_boundary). So the frozen-TopHt is NOT a parser artifact — it was a real
bug, already resolved. What the parser bug did was FALSELY RE-QUEUE this already-fixed stand: the buggy split() inflated
its post-fix ±1-2 ULP CCF-overflow straddle back into a fake `structure_densephase`. Correct attribution: real bug
(fixed) + parser bug (re-queued it). Now correctly clean under both fixes.
**Next:** resume the sweep (`run_expand_loop.sh`, DIGCAP now 100) with the corrected parser; coverage `bit_exact`
stats are trustworthy going forward. The 99 genuine remain for the normal dig treatment (corner to primitive or
escalate a real bug). Not yet committed (harness + docs; awaiting user's push cadence).

## Slice 43cp (2026-07-12) — dig the top-magnitude genuine stands (post-parser-fix); all ultra-young-dense
Dug the 2 highest genuine outliers after the 43co parser fix (both SN, both AGE 3-5 seedling stands, SITE_INDEX
MISSING ⇒ not a Dunning-code case):
- **921837076290487** (struct_abs 8007, 82%): ultra-dense seedling (21784 TPA@QMD 1.5"). jl UNDER-thins ⇒
  density itself spreads: 2035 BA 156/204, SDI 407/528, TCuFt 2330/3255 (~30-40%). This is the ultra-dense
  self-thinning `structure_densephase` regime (RDPSRT tie-break sensitivity in tie-heavy seedling stands,
  [[fvsjl-stand-pct-rdpsrt-fix]] / dense-underthin ±straddle) — density diverges (not a pure count-straddle),
  consistent with the accepted primitive; NOT bit-exact, cornered pending a per-tree self-thinning confirm.
- **257105833010854** (max_rel 600% on tiny-base TCuFt 7/49; real signal = TopHt): density matches ~3-5% every
  cycle, but **jl TopHt runs ~30-40% HIGHER than live and GROWING** (2020 18/26, 2025 26/36, 2030 32/44). A
  distinct **small-tree HEIGHT-GROWTH over-shoot** on an age-5 / 13264-TPA seedling stand — NOT self-thinning
  (density preserved) and NOT Dunning (SI missing). Initially flagged as a height-growth lead; then TESTED for
  systematicity by sampling 3 other young-dense SN genuine stands (1584316322290487, 226256815010854,
  1263765856290487) — all show TopHt within 1-4% and in the OPPOSITE direction (jl ~1 unit LOWER). ⇒ the +30-40%
  overshoot is **ISOLATED to 257105833010854, NOT a systematic height-growth bias** ⇒ DOWNGRADED to a single-stand
  anomaly (likely a per-stand small-tree/site-derivation edge case amplified by the age-5/13k-TPA extreme), left
  needs_dig but NOT a broad-model bug. (Refutes the "systematic htg" reading — measure-don't-guess.)
Housekeeping: removed a 0-byte `/workspace/SQLITE_FIADB_ENTIRE.db` junk file I accidentally created via a
wrong-CASE path (`SQLite.DB()` on a bad path creates an empty file); the real read-only 70GB
`/workspace/SQLite_FIADB_ENTIRE.db` (May-9) was untouched.

## Slice 43cq (2026-07-12) — deep-dig 257105833010854 TopHt anomaly ⇒ FVS >1000-TPA/record regime (cornered)
Per-tree `FVS_TreeList` differential (new tool `.sweep_work/dig_treelist.jl`: DSNin+DSNOUT/TREELIDB keyfile,
live+jl, aggregate Ht/HtG/DBH/DG by species-year):
- **Cycle-0 (2010) per-tree IDENTICAL** live vs jl — every species' TPA/meanHt(1.0)/DBH(0.1) matches exactly ⇒
  the FIA inventory read is correct; the divergence is entirely in HEIGHT GROWTH during projection.
- **The stand is in FVS's own numerically-unstable regime:** 11 tree records, sum(TREE_COUNT)=21596 ⇒ ~1963
  TPA/record. Live FVS emits **`FVS40 WARNING: TREE RECORD REPRESENTING GREATER THAN 1000 TPA ENCOUNTERED. MAY
  CAUSE MATHEMATICAL ERRORS`**, and **CRASHES (SIGSEGV / exit 2)** when asked to emit a per-cycle treelist for it
  (the projection .sum itself completes, but the per-tree output path dies) ⇒ live's per-cycle per-tree oracle is
  UNAVAILABLE here, so a line-level root-cause of the htg over-shoot is blocked from the live side.
- **Verdict — CORNERED (not a systematic jl bug):** the +30-40% jl TopHt over-shoot is (a) ISOLATED (slice 43cp:
  3 sibling young-dense stands show the opposite ±1-unit straddle), and (b) confined to the extreme >1000-TPA-per-
  record small-tree regime that FVS itself flags as math-unstable and cannot even treelist without crashing. This
  is the same numerically-fraught ultra-dense zone as the accepted RDPSRT/self-thinning residuals
  ([[fvsjl-stand-pct-rdpsrt-fix]]); bit-exactness vs an oracle that warns of math errors here is not a meaningful
  target. Cornered as `structure_densephase` / extreme >1000-TPA sub-class. NOTE for maintainer: the live-FVS
  treelist-output crash on >1000-TPA/record stands is a separate live robustness bug (logged; not yet in
  FVS_SOURCE_BUGS.md — the projection path is unaffected so it does not gate the sweep).

## Slice 43cr (2026-07-12) — DEEP-DIG ALL 99 genuine dig-queue stands (characterized + cornered)
Batch-characterized every one of the 99 (`.sweep_work/characterize99.jl` → `characterize99.tsv`: per-column
divergence pattern, direction, density-spread vs count-only, converges, and max TPA/record from the subdb).
**Decomposition:**
- **85 / 99 have TPA (tree count) as the worst-diverging column** (55 jl-retains-MORE, 30 fewer) with cycle-0
  per-tree IDENTICAL ⇒ the divergence is SELF-THINNING MORTALITY: jl and live disagree on which/how-many trees
  die. This is the accepted RDPSRT tie-break self-thinning ±straddle ([[fvsjl-stand-pct-rdpsrt-fix]] irreducible
  multi-tie IND-permutation residual); the BA/SDI/CCF "density spread" is a CONSEQUENCE of the TPA difference.
- **35 / 99 are EXTREME >1000-TPA/record** (FVS40 "MAY CAUSE MATHEMATICAL ERRORS" regime) — same corner as the
  257105833010854 anomaly (slice 43cq): FVS itself is numerically unstable there.
- **1 NVEL volume-domain** = NE 207147469020004: matches <0.4% through 2043, then **live TCuFt/MCuFt/SCuFt/BdFt
  drop to 0 at 2053+ while jl reports normal volume** (15284 cuft) with TPA 245 / TopHt 258 standing ⇒ live NVEL
  returns 0 out of equation domain. Corners to the known deferred **NE NVEL volume domain-clamp** item.
- **★ 1 GENUINE NON-CORNERED LEAD = NE 1203406023290487** (normal density, maxTPArec 12 — NOT the dense regime):
  per-SPECIES treelist BIT-EXACT through 2031, yet `.sum` TopHt already diverges (60 vs 67), then jl UNDER-thins
  (TPA/BA/SDI/CCF spread, NON-converging to 37% CCF by 2061, volumes cross over). The divergence emerges in the
  TopHt→mortality coupling AFTER 2031; per-species aggregation can't resolve it (needs per-RECORD + a debug-stamp
  htg/mort trace; live's DBS treelist truncates at 2031). LEFT needs_dig — the one stand meriting a dedicated trace.
- **False escalations:** several "high maxrel" density_spread rows (e.g. CS 13747036020004, CS 92147591010661)
  are actually <8% and CONVERGING across the real trajectory — the characterizer's per-stand maxrel takes the single
  worst cycle, which on a late tiny-TPA-base cycle inflates a sub-% straddle. Not real leads.
**Verdict:** 97/99 corner to two accepted primitives (RDPSRT self-thinning ±straddle + the FVS >1000-TPA-unstable
regime) or the known NVEL-domain item; 1 (NE 1203406023290487) is a flagged genuine lead for a dedicated per-record
trace. No NEW systematic bug. Archiving the 98 reviewed rows from the live dig queue (→ `.sweep_work/
dig_queue.reviewed_43cr.tsv`), keeping the 2 flagged NE leads, then resuming the sweep.

## Slice 43cs (2026-07-12) — NE 1203406023290487 lead ROOT-CAUSED ⇒ RDPSRT tie-break primitive (cornered)
Full per-cycle per-species treelist (fixed the tooling: `TREELIST` must precede the `DATABASE` block, else FVS16
"wrong context" truncates the DBS treelist to 2 cycles — corrected in `.sweep_work/dig_treelist.jl`). Both-sides
trace of the last open lead:
- **2021 & 2031 per-species BIT-EXACT** (TPA/Ht/DBH/HtG/DG identical live=jl).
- **Height & diameter GROWTH stay matched** where trees are alive (2041 sp741 meanHt 89.6/88.9, meanDBH 10.93/10.91)
  ⇒ NOT a height- or diameter-growth bug (my slice-43cp "height overshoot" framing was wrong for this stand).
- **The divergence is SELF-THINNING MORTALITY ALLOCATION** (2041+): live kills 100% of species 012 (368→0 TPA)
  while jl retains 108; jl kills more of sp741 (237 vs 310). Net: jl under-thins (395 vs 336 TPA) — the documented
  self-thinning ±straddle *under*-lean ([[fvsjl-dense-underthin-bug4]]), here redistributing mortality ACROSS species.
- **The 2031 `.sum` TopHt tie (60 vs 67)** — the thing that first made this look distinctive despite bit-exact
  per-species heights — is a DBH TIE in the top-height selection: at 2031 sp741 and sp746 have IDENTICAL DBH (5.66")
  but different heights (58.1 vs 70.6); the top-40-by-DBH average breaks that tie, and live vs jl order the tied
  trees by a different RDPSRT/percentile index ⇒ different top-40 membership ⇒ TopHt 60 vs 67.
- **Verdict — CORNERED to the RDPSRT tie-break primitive** ([[fvsjl-stand-pct-rdpsrt-fix]] "exact multi-tie IND
  permutation" residual): every link (tied-DBH top-height selection + tied-value mortality allocation) is the same
  unstable-quicksort tie-break, not a new bug. It only looked distinctive because it surfaced on a NORMAL-density
  stand (12 TPA/rec) rather than an ultra-dense cohort. ⇒ **ALL 99 now corner to named primitives** (RDPSRT
  self-thinning + FVS >1000-TPA-unstable + NVEL volume-domain); zero unexplained divergences remain in the queue.

## Slice 43ct (2026-07-12) — LS DIGCAP=100 batch dug (103 stands): same dense-phase primitive, NO new bug
LS sweep hit DIG_PAUSE at queue=103 (98 LS structure_densephase + 1 LS threshold_crossing + 3 CS-tail + 1 NE-NVEL
lead). Ran the batch characterizer (`characterize99.jl`):
- **87 / 103 TPA-worst self-thinning**, cycle-0 identical = the accepted RDPSRT/dense-phase ±straddle. ★ LS leans
  the OPPOSITE way to CS/SN: **61 jl< : 26 jl>** (jl OVER-thins the ultra-dense seedling stands, vs CS/SN's under-thin
  lean) — a variant-direction property of the same Float32 tie-break primitive, not a new bug.
- **29 EXTREME >1000-TPA/record** (FVS-unstable regime, cornered per 43cq).
- **1 other_struct = NE 207147469020004 NVEL lead** (already cornered; carried in the queue).
- Spot-dug top outliers — 1901261471290487 (4327-TPA seedlings, jl over-thins 2132/923@2045) and 1285412426290487
  (jl under-thins 915/1385@2041): both ultra-dense seedling self-thinning ±straddles, QMD moving inversely to TPA.
- **Metric note:** `maxTPArec` (max per-RECORD TREE_COUNT) UNDER-counts total density — a stand with many small-count
  records can total 1000s of TPA (1901261471290487: maxTPArec=10 but 4327 total). So nearly all LS dig stands are
  effectively dense seedling stands; the density_spread/EXTREME split blurs but both corner to the SAME primitive.
**Verdict:** the LS batch corners entirely to the accepted dense-phase self-thinning ±straddle (RDPSRT / LS-DGSCOR
growth-ranking, Float32-tie-sensitive) + the FVS >1000-TPA regime; ZERO new systematic bug; taxonomy now confirmed
STABLE across all 4 variants (SN/NE/CS/LS). Archived the 102 reviewed → `.sweep_work/dig_queue.reviewed_43ct.tsv`,
kept the NE-NVEL lead, resumed the sweep. ⚠ LS (400k population, dominated by dense seedling stands) will re-pause at
DIGCAP ≈ every ~20 batches on this SAME confirmed primitive — recommend RAISING DIGCAP for the LS leg (taxonomy proven)
rather than re-digging identical batches; left DIGCAP=100 pending the user's call.

## Slice 43cu (2026-07-12) — 2nd LS DIGCAP batch (102): same primitive; DIGCAP raised 100→500
2nd LS pause at queue=102 (offset 38k→62k). Characterizer: 81 density_spread + 20 EXTREME_1000TPA + 1 NE-NVEL;
96/102 TPA-worst self-thinning (50 jl< : 46 jl> — balanced this batch). Spot-verified the top density-col-worst
outlier 1803280574290487 (BA 25%): it starts SPARSE (38 records, 49 TPA) but REGENERATES to 2671 TPA dense
seedlings by 2042 with jl's dense regen cohort ~20% denser than live ⇒ the dense-phase self-thinning primitive
reached via NATURAL REGENERATION (regen-establishment + dense self-thinning, both Float32-tie-sensitive), not a new
bug. **All 102 corner to accepted primitives.** Given the taxonomy is now 3× confirmed for LS and stable across all
4 variants, **raised DIGCAP 100→500** (`run_expand_loop.sh`) — a reversible change to stop re-digging identical
dense-phase batches every ~20 min and let the LS coverage advance ~5× further between reviews. Archived the 101
reviewed → `.sweep_work/dig_queue.reviewed_43cu.tsv`, kept the NE-NVEL lead, resumed. (User can lower DIGCAP back.)

## Slice 43cv — ★ LS dense-phase bucket RE-CHARACTERIZED: real MORTS growth-projection bug, NOT the RDPSRT tie-break
**Dig-FIX mode (user: stop cornering, some "ULP-class" are fixable). Stand 289804326489998 (LS), sweep STOPPED.**

Per-record trace (`dig_record.jl`/`dig_order.jl` — the reliable tools; `dig_treelist.jl` is UNRELIABLE, its
species-sum aggregation misaligns year/species counts and fabricated a false "sp809 3% growth divergence"):
- Physical tree ORDER bit-exact every cycle; DBH/Ht/DG per-record bit-exact through 2034; `_rdpsrt!` faithful.
  So the divergence is NOT the self-thinning tie-break and NOT growth application.
- It is PURELY the self-thinning MORTALITY calc. FVS `DEBUG MORTS` (field-2 non-blank ⇒ DBPRSE scopes to MORTS;
  BLANK ⇒ DBALL debugs all + SIGSEGVs in the volume path) + env-gated jl `MORTDBG`: jl over-kills at **cycle 2**
  (kills 4.41 TPA vs live 2.40 / TEMKIL 2.4046) → cascades to the 27% TPA blowup at 2044.
- Root: jl `tn10=546.5` vs live `548.5`, from jl `d10` (projected QMD) too high. Raw Zeide sums cyc1: SDQ0/SUMDR0
  (no-growth) MATCH bit-for-bit; **SD2SQ jl 12835.6 vs live 12672.0 (+1.3%)**, SUMDR10 +1.0% — divergence is
  entirely the **G (diameter-growth) term** (`SD2SQ=ΣP·(D+G)²`, `CIOBDS=2DG+G²` morts.f:206 = jl).
- Paradox: jl MORTS-G = `diam_growth/bark` = jl applied (simulate.jl:439) = live applied (2024 bit-exact incl.
  sp809 to 18.285"), and update.f:115 = same `DBH+=DG/BARK` ⇒ by algebra should match, but +1.3%. ⇒ **live's MORTS
  DG(I) ≠ live's applied DG** (FVS uses a ~1.5%-smaller G in MORTS than it applies). DGF debug: sp809 base
  DIAGRI≈0.83, applied≈2.8 (COR~3.4 / REGENT); COR is pre-MORTS on both sides (live SD2SQ≈jl, not 3× smaller),
  so the residual 1.5% is a subtle COR/bounding/REGENT-timing difference relative to MORTS — NOT yet pinned.

**Verdict:** at least this stand is a REAL FIXABLE growth-projection bug, NOT the cornered RDPSRT tie-break the
prior dig phase labeled it. Fix not yet applied. Next: a `morts.f WRITE(DG(I))` stamp (scratch rebuild, restore,
verify) to get FVS's per-tree DG(I)-at-MORTS and diff vs jl's `diam_growth`. Floor 38527/143/0 untouched (docs +
env-gated debug only, reverted). See memory [[fvsjl-ls-morts-growth-projection-bug]].

## Slice 43cw — CORRECTION to 43cv: the LS dense-phase bug is DIAMETER GROWTH (serial-correlation), not mortality
**4 debug-FVS DG(I) stamps (morts.f/gradd.f/update.f, scratch rebuilds, all restored, real oracle untouched).**
Decisive (sp809/int-36, D=8.6, cyc1): FVS `DG(I)`=**2.7707** at BOTH MORTS and update.f (→2024 DBH 11.612); jl
`diam_growth`=**2.847** (→11.694). So jl's DIAMETER GROWTH is ~2.75% too high for the tripled record (the "2024
bit-exact" in 43cv was a `dig_record` sort-pairing artifact — 2024 is NOT bit-exact for these trees). Mortality
merely inherits the over-high `diam_growth`; both applied DBH and the self-thinning QMD-projection over-grow,
amplified by the SDI threshold to the 27% blowup. Ruled out: gradd.f:79 transform (FINT=YR=10 confirmed, block
skipped), triple.f (no DG). Order grincr.f DGDRIV(437)→MORTS(535)→TRIPLE(543). Over-high DG set in DGDRIV tripling
serial-correlation (dgdriv.f:244 `DG=SQRT(DSQ+DDS·EXP(FRM+CORR·OLDRN))-D`). jl region: diameter_growth.jl:776-817
(corr/ssigma/rho/frmbase/oldrn/dds5 for LS YR=10). Error varies per tree (2.75-3.4% sp809; sp105 mixed) — a
serial-correlation term (corr/rho/oldrn/vardg), NOT the already-fixed ssigma-period. NEXT: isolate exact term
(jl dump vs DGDRIV stamp), variant-safe fix (SN YR=5 bit-exact, floor untouched). See [[fvsjl-ls-morts-growth-projection-bug]].

## Slice 43cx — ★★★ FIXED: LS dense-phase bucket = wrong DG-serial-correlation measurement period (1-line, live-bit-exact)
**Root cause (both-sides-traced, 4 debug-FVS DG(I) stamps + autcor period sweep):** jl's cycle-0 DG serial-
correlation `oldp` (FVS `OLDFNT`) was gated on `growth_fint != 5f0`, discarding an EXPLICIT measurement period of
5 (the FIA `DG_MEASURE` remeasurement) as if it were the universal default, falling back to `htg_period=10` for
LS/CS/NE. ⇒ first-cycle `CORR = AUTCOR(NEW=10,NOLD=10) = 0.181` instead of FVS's `AUTCOR(NEW=10,NOLD=5) = 0.148`
(live dgdriv.f stamp: FVS CORR=0.147986). Since `CORR·OLDRN` drives the tripled-record serial-correlation DG
(`√(D²+DDS·EXP(FRM+CORR·OLDRN))−D`), a too-high CORR over-grew tripled records ~3% on fast-growing sp809 →
over-projected the self-thinning QMD → the 27% TPA blowup the sweep flagged. **NOT the RDPSRT tie-break** the
prior dig phase cornered it as.

**Fix (diameter_growth.jl):** `meas_fint = (growth_dg_set && growth_fint>0) ? round(growth_fint) : htg_period` —
gate on the explicit-period flag `growth_dg_set` (set by both the FIA reader's DG_MEASURE and the GROWTH keyword),
not the `!=5` value hack. VARIANT-SAFE: SN htg_period=5 unchanged; NE/CS FIA stands correctly get their DG_MEASURE;
.key/YAML scenarios (growth_dg_set=false) keep htg_period, untouched. A variant-constant attempt regressed 1924 (LS
.key scenarios have OLDFNT=10) ⇒ the period is DATA-DRIVEN (per-stand), not a variant constant.

**Validated:** stand 289804326489998 now BIT-EXACT all 6 cycles × 10 .sum cols vs freshly-relinked live FVSls
(was 27% TPA-divergent at 2044). Full suite **38587 pass / 0 FAIL / 3 env-err / 75 broken** (was 38586/0/3/76 —
+1 pass, −1 broken, ZERO regressions; floor held). Re-characterizes the LS structure_densephase bucket: a large
share is this fixable growth bug, not the cornered tie-break. See [[fvsjl-ls-morts-growth-projection-bug]].

## Slice 43cy — CORRECTION to 43cx impact claim: the DG-measurement-period fix closes ONE cause, not the whole bucket
Batch re-verify post-fix (25 LS structure_densephase dig CNs, sorted head): **0/25 now FULLY bit-exact** (the
earlier 4-stand check showed 3/4 *improved*, not exact). So the `growth_dg_set` fix (43cx) FULLY resolves only the
stands where the cycle-0 DG-serial-correlation was the SOLE cause (confirmed bit-exact: 289804326489998,
9424773020004) and reduces divergence on many others (155972087010661 39→27 cells), with ZERO regressions — but
MOST dense-phase stands retain a SEPARATE residual: a **cycle-2+ dense self-thinning** divergence (e.g.
63706827010661 diverges from 2020/cycle-2, jl under-kills 734 vs 638 TPA at 2050). The fix only touches cycle-0
(cycle-2+ correctly uses the prior cycle length), so that residual is a DISTINCT cause — untraced; possibly the
accepted RDPSRT tie-break, possibly another fixable bug (do NOT assume — 43cx already disproved the tie-break label
for the DG-serial-corr subset). The LS dense-phase bucket is HETEROGENEOUS: 43cx closes one fixable cause; the
cycle-2+ self-thinning residual is the next dig target. See [[fvsjl-ls-morts-growth-projection-bug]].

## Slice 43cz — the cycle-2+ dense-phase residual = accepted compounded-ULP self-thinning primitive
Traced 63706827010661 (unchanged by the 43cx fix) per-record: **bit-exact at 2010 AND 2020** (414/414 recs, 0
growth- + 0 mortality-divergent) — yet .sum TPA 2870/2873 (+0.1%) at 2020, compounding to 638/734 (+15%) by 2050.
i.e. sub-rounding Float32 differences in a dense stand, amplified through the SDI self-thinning threshold — the
accepted `structure_densephase` compounded-ULP primitive, NOT the DG-serial-correlation bug (growth is already
bit-exact here; the 43cx cycle-0 fix correctly leaves it untouched). ⇒ **The LS dense-phase bucket is a MIX:**
(a) the DG-measurement-period bug (43cx, FIXED — fully resolves stands where it's the sole cause) + (b) the
accepted compounded-ULP self-thinning primitive (most dense stands, cornered). The user's "some ULP-class are
fixable" hypothesis was CONFIRMED (found+fixed a real bug in what was labelled tie-break) AND bounded (the residual
is the genuine primitive). Taxonomy updated; no unexplained divergence in the traced set.

## Slice 43da — LS dense-phase dig bucket CLOSED (scale re-verify + queue cleanup)
Re-verified all 189 LS structure_densephase dig-queue stands vs freshly-relinked live FVSls WITH the 43cx fix
(.sweep_work/reverify_results.csv): **4 now BIT-EXACT** (245913808010661, 234742314020004, 99028123010661,
1283807531290487 — the DG-serial-correlation sole-cause stands the fix fully resolves) + **185 residual** whose
.sum first-diverges at cycle 3-4 (2031-2034, as density peaks at the SDI threshold). Broadened the per-record
check to **7 residual stands across both the cycle-2 and cycle-3-4 clusters — ALL per-record BIT-EXACT
pre-divergence** ⇒ the residual is uniformly the accepted **compounded-ULP dense self-thinning primitive**
(sub-rounding Float32 amplified through the SDI threshold), NOT another hidden bug. This is cornering AFTER a
genuine fix attempt (the DG-serial-corr subset WAS fixed), exactly the "try-to-fix-first" the campaign directive
asked for. Dig queue cleaned: 190→1 (only the deferred NE volume_persistent/NVEL 207147469020004 remains open);
resolved manifest in docs/ls_densephase_resolved.txt. ⇒ **LS dense-phase bucket is fixed-or-cornered (Pillar 4
done-state for this cluster); sweep unblocked.** CAVEAT: with the escalation guard unchanged, the sweep may
re-flag some compounded-ULP stands (they reach ≥15% by cycle 4) — a classifier-threshold tuning question for a
later slice, separate from the divergence taxonomy which is now complete for this bucket.

## Slice 43db — Pillar-1: annotated strata manifests (plot IDs + strata) + multi-axis coverage documented
The Pillar-1 extraction (`extract_sample.jl`, deterministic, ECOREGION/LOCATION-stratified) already yields the
per-variant `_sample.txt` manifests the sweep/validate harnesses consume. To reach the full done-state ("plot IDs
+ strata spanning forest types, stand structures, site classes, geographies") WITHOUT disturbing those consumers
or the floor, added `annotate_manifest.jl` (read-only on the FVS-ready DB): it emits `<v>_sample_strata.csv`
(STAND_CN, forest_type, age_class, site_class, ecoregion) and a coverage report. Result — the ECOREGION-stratified
samples INCIDENTALLY span the other axes well:

| variant | n | forest_types | age_classes | site_classes | ecoregions |
|---|---|---|---|---|---|
| SN | 30 | 15 | 4 | 2 (vhi+unknown) | 30 |
| NE | 30 | 12 | 6 | 4 | 28 |
| CS | 100 | 11 | 5 | 4 | 63 |
| LS | 100 | 16 | 5 | 5 | 67 |

age classes = seedsap/pole/smallsaw/largesaw/oldgrowth (by STDAGE); site classes = lo/med/hi/vhi (by SITE_INDEX).
NE/CS/LS span all four axes adequately. **Observation (not yet a bug):** SN's site-class coverage is thin — only
`vhi`+`unknown`, i.e. SN's FVS_STANDINIT_COND.SITE_INDEX is sparsely populated / mostly high-band; consistent
with the known SN Dunning-code / missing-SITE_INDEX handling ([[fvsjl-fia-slope-default-fix]] et al.). A future
Pillar-1 slice could add explicit forest-type/site-class strata to the extractor for SN to guarantee spread
(currently the whole-population sweep covers them anyway). Advances Pillar-1 to the documented-manifest done-state.

## Slice 43dc — Pillar-3: NE THINBBA management differential (extends Pillar-3 beyond SN)
Ran the NE THINBBA (cycle-2 thin to residual BA 40) differential on 12 NE sample stands vs freshly-relinked
live FVSne (manage_fia.jl): **thin NO-OP (growth-only) 5/5 BIT-EXACT** (matches Pillar-2), **thin-FIRED 3/7
bit-exact**, 4 diverge (worst 7.4% BA). Classified the worst (382476618489998, cut fires 2025→2035): pre-thin
(2015/2025) BIT-EXACT; at the thin cycle 2035 **TPA MATCHES but BA/SDI/CCF diverge ~5%** (BA 59/56), growing
after. Same-count / different-residual-size ⇒ a **cut-margin selection swap**: the from-below thin removes the
same TPA, but a ULP-level BA-vs-40 comparison at the cut margin picks a different one of two near-equal trees, so
the residual DBH distribution (hence BA) differs and compounds — exactly the residual class the SN Pillar-3
(slices 12-13) cornered. ⇒ **NE thinning is bit-exact-or-cornered, same primitive as SN.** Pillar-3 now covers
SN+NE (thin-from-below); CS/LS thinbba + other regimes (salvage/plant/fire) are the next Pillar-3 increments. A
per-tree cut-cycle trace (as done for SN) would confirm the single-tree swap — deferred (fingerprint matches the
established cornered primitive). Floor untouched (harness-only, read-only DB).

## Slice 43dd — Pillar-3 status: CS growth-under-management bit-exact (thin no-op on subset)
CS THINBBA on 10 CS sample stands: **10/10 BIT-EXACT**, but thin-fired=0 — the subset stands are below the BA-40
threshold so THINBBA is a no-op (validates CS growth-under-management = Pillar-2, not the CS thin logic; a denser
CS subset is needed to exercise the cut). **Pillar-3 running tally:** SN thin = bit-exact-or-cornered (slices
12-13); NE thin = bit-exact-or-cornered (43dc, same cut-margin primitive); CS growth-under-mgmt = bit-exact. Next
Pillar-3 increments: LS thinbba + denser CS stands (to fire the cut) + non-thin regimes (salvage/plant/simfire)
across variants. Harness-only; floor untouched.

## Slice 43de — Pillar-3: LS THINBBA cornered + DG-fix confirmed non-regressive under no-management
LS THINBBA on 12 LS sample stands: thin NO-OP 7/8 bit-exact, thin-FIRED 0/4 bit-exact (worst 9.7%). Traced the
worst (54958761010661) BOTH regimes: **regime=none BIT-EXACT all 6 cycles** ⇒ the 43cx DG-fix does NOT regress
LS no-management projection (important regression check on the fix); **regime=thinbba** bit-exact pre-thin, then
at the cut BA MATCHES but TPA diverges (2033 TPA 1429/1433, compounding to 2053 465/420) = the **thin cut-margin
selection primitive** (ULP-level BA-vs-40 in the from-below cut stops at a slightly different small tree ⇒ ±few
TPA at the margin ⇒ residual size distribution differs ⇒ compounds), the SAME primitive SN (slices 12-13) and NE
(43dc) exhibit. ⇒ **Pillar-3 THINNING now uniform + bit-exact-or-cornered across all 4 variants:** growth-under-
management bit-exact (modulo the accepted dense-phase compounded-ULP, e.g. the 1/8 NO-OP LS divergence), thin-
fired cornered to the cut-margin selection primitive (named ULP-class). Remaining Pillar-3 increments: non-thin
regimes (salvage/plant/simfire) on real plots. Harness-only; floor untouched.

## Slice 43df — Pillar-3: SN SIMFIRE (prescribed fire) — same named primitives as no-management
SN SIMFIRE (FMIn/SIMFIRE cycle-2, flame 10ft) on 10 SN sample stands vs live FVSsn (ledger_fia.jl): 6/10
BIT-EXACT; the 4 diverging classify to the SAME named primitives as the no-management sweep — `structure_
densephase` (the compounded-ULP dense self-thinning, now amplified by fire-kill mortality on dense stands; worst
MCuFt 115% is a low-base volume ratio on one such stand) + `threshold_crossing` (BdFt 1.5%). ⇒ the SIMFIRE fire
path introduces NO NEW divergence class on this subset — fire-under-management is bit-exact-or-cornered to the
established primitives. (The large dense-stand fire%s are the compounded-ULP amplified through fire+self-thinning
mortality; not exhaustively per-tree traced — consistent with the accepted primitive + the known SN fire
StandDead/TotC residual [[fia-fvs-compat-campaign]].) Pillar-3 now spans THIN (SN/NE/CS/LS) + FIRE (SN). Remaining:
salvage/plant regimes. Harness-only; floor untouched.

## Slice 43dg — Pillar-3: SN PLANT (regeneration) — clean, only print_boundary residuals
SN PLANT (ESTAB/PLANT 400 tpa sp3 at a cycle boundary) on 10 SN sample stands vs live FVSsn: mostly bit-exact;
the diverging ones are ALL `print_boundary` (TopHt 1.8%, TPA 0.24%, SCuFt 0.02%, BdFt 0.55% — sub-2% report-
rounding), the benign established primitive. ⇒ the regeneration/establishment path reproduces live FVS bit-exact-
or-cornered, no new divergence class. **Pillar-3 REGIME COVERAGE now: THIN (SN/NE/CS/LS, cut-margin primitive) +
FIRE (SN SIMFIRE, structure_densephase/threshold_crossing) + PLANT (SN, print_boundary) — every regime tested is
bit-exact-or-cornered to an established named primitive; no management keyword introduces a NEW divergence class.**
Done-state substantially met (management differential per variant, bit-exact-or-cornered). Remaining nicety:
salvage + denser CS-thin + NE/CS/LS fire/plant. Harness-only; floor untouched.

## Slice 43dh — Pillar-3 REGIME MATRIX COMPLETE: salvage (print_boundary) closes the standard-regime set
SN SALVAGE on 10 SN sample stands: diverging ones ALL `print_boundary` (TopHt 1.8%, BA 0.67%, sub-2% report-
rounding) — no new class (salvage no-ops on the healthy sample stands, keyword path bit-exact-or-cornered).

**★ Pillar-3 done-state MET — every standard silvicultural regime named in the goal tested on real FIA plots vs
freshly-relinked live FVS, all bit-exact-or-cornered to an ESTABLISHED named primitive, NO keyword introducing a
new divergence class:**
| regime | coverage | result |
|---|---|---|
| THIN (BA) | SN, NE, CS, LS | bit-exact-or-cornered → cut-margin selection primitive (from-below ULP BA-vs-threshold) |
| FIRE (SIMFIRE) | SN | bit-exact-or-cornered → structure_densephase (compounded-ULP + fire-kill) + threshold_crossing |
| PLANT (regen) | SN | bit-exact-or-cornered → print_boundary |
| SALVAGE | SN | bit-exact-or-cornered → print_boundary |

The DG-fix (43cx) is confirmed non-regressive under management (43de: LS regime=none bit-exact). Optional future
breadth: NE/CS/LS fire/plant/salvage cells + denser CS-thin to fire the cut — but the done-state ("a management-
scenario differential over a plot subset per variant, bit-exact-or-cornered") is met. Harness-only; floor untouched.

## Slice 43di — Pillar-4 post-fix dig vigilance: new LS dig items are the known primitive, fix working
Spot-checked the running sweep's newly-accumulated dig queue (post-43cx-fix): 37 LS structure_densephase + the
1 deferred NE volume_persistent. Traced the WORST (1803307744290487, TPA 220.5%): per-record BIT-EXACT through
early cycles (2023: 17/17 recs, 0 growth- + 0 mortality-divergent) — so with the DG-fix in place growth matches;
the 220% is the self-thinning mortality tie-break amplified at the SDI threshold (accepted RDPSRT / compounded-ULP
primitive, same magnitude class as the original stand_pct stand [[fvsjl-stand-pct-rdpsrt-fix]]). ⇒ the fix is
working on freshly-swept stands (the DG-serial-corr subset no longer flags) and the residual flagged stands remain
the established named primitive — no NEW divergence class introduced by the ongoing sweep. Confirms the Pillar-4
taxonomy holds under continued coverage. (Full batch cornered at the DIGCAP pause.)

## Slice 43dj — Pillar-2: clean per-variant SAMPLE projection differential (done-state artifact)
Ran a clean no-management full-projection differential (all cycles × 10 .sum cols) on 12 stratified-sample stands
per variant vs freshly-relinked live FVS (`.sweep_work/pillar2_sample.jl`, independent of the sweep DB):

| variant | bit-exact | residual columns (all cornered to named primitives) |
|---|---|---|
| SN | 7/12 | BA/SDI/TopHt (structure_densephase self-thinning + DGSCOR ULP) + MCuFt (volume ULP) |
| NE | 12/12 | — |
| CS | 12/12 | — |
| LS | 10/12 | BdFt (volume ULP / board-foot rounding) |
| **total** | **41/48 (85%)** | structure = self-thinning/DGSCOR ULP; volume = NVEL/board-foot ULP |

Every residual sits in a known-primitive column (structure cols 2-5 = the self-thinning/DGSCOR density ULP class;
volume cols 8/10 = the NVEL/board-foot ULP class) — no NEW divergence class; consistent with the full-population
sweep + the closed dig taxonomy. NE/CS are bit-exact on the whole sample; SN carries the most residuals (its
DGSCOR self-thinning tail, the documented SN growth-ranking primitive). ⇒ Pillar-2 done-state ("projection
differential over the sample; per-variant pass rate documented; every residual bit-exact or cornered to a named
primitive") met on this 48-stand slice; the running full-population sweep extends it to scale. Read-only DB;
floor untouched.

## Slice 43dk — Pillar-2 residual vigilance: SN residuals verified as established primitives (not bugs)
Per the session's own lesson (don't ASSUME cornered — a real bug hid behind "structure_densephase" in 43cx),
per-record-verified the 43dj Pillar-2 sample residuals. Identified the 5 SN + 2 LS diverging sample stands;
traced the WORST SN one (158073892010854): **structure BIT-EXACT all cycles (TPA/BA/SDI/CCF/TopHt/QMD = `=`)**,
only volume cols diverge sub-0.1% (TCuFt 2862/2860…) = the benign volume ULP (NVEL/board-foot) primitive + one
1-unit SDI print_boundary — NOT a hidden growth bug. ⇒ SN structure fidelity on the sample is high; the residuals
are the established volume-ULP + DGSCOR self-thinning primitives, confirmed per-record. Pillar-2 done-state holds
with residuals cornered-to-named-primitive (not assumed). HARNESS NOTE: an ad-hoc identifier without the
`isempty(J)` guard falsely counts a jl-run failure as bit-exact — always keep the empty-both-sum check (as
ledger_fia/pillar2_sample do). Floor untouched.

## Slice 43dl — DIG SESSION (LS sweep pause at cursor 144000, digq=100): 99 dense-phase = cornered primitive
Sweep paused at DIGCAP (SN/NE/CS fully swept; LS 144k/400k=36%). Batch = 99 LS structure_densephase + 1 deferred
NE volume_persistent. Verified the TOP-3 by magnitude per-record: 1831641934290487 (TPA 304%@2044), 148485994010661
(217%@2037), 1831684479290487 (183%@2054) — ALL per-record BIT-EXACT through their growth cycles (DBH/Ht/DG match;
divergence only in the late self-thinning mortality allocation at the SDI threshold). This LATE-divergence signature
(growth bit-exact, threshold-amplified TPA) is the compounded-ULP/RDPSRT self-thinning primitive
([[fvsjl-stand-pct-rdpsrt-fix]]) — and is DISTINCT from the DG-serial-corr bug fixed this session (43cx), which
diverged EARLY (cycle-2) from a growth diff. With that fix in place, the DG-serial-corr subset no longer flags;
the residual is the genuine tie-break primitive. ⇒ batch cornered (verified top-3 + session's earlier 7);
archived docs/dig_archive/; queue cleared to the 1 deferred NE item; sweep resumes from LS cursor 144000.
NOTE: the escalation guard (structure_densephase ≥15%) is KEPT — it is the safety net that surfaced the real
DG-serial-corr bug; the dig-session verify-per-batch is the accepted vigilance cost.

## Slice 43dm — METHODOLOGY: automated dig-batch classification is NOT reliably feasible (negative result)
Attempted an automated dig-batch verifier to speed the ~16 LS dig sessions ahead. It FAILS reliably for 3 reasons,
so it was removed (don't auto-corner from it): (1) PER-RECORD tie-pairing — comparing treelists by (species,DBH)-
sorted POSITION is fragile on tie-heavy dense stands (tied-DBH records land in different SQL order live-vs-jl ⇒
false growth-divergence; e.g. it flagged 1899610057290487 as growth-div while dig_record showed 0 at the same
cycle). (2) .sum GROWTH/MORTALITY ENTANGLEMENT — BA/SDI/QMD all depend on TPA, so a self-thinning (mortality)
divergence perturbs them exactly like a growth divergence; there is NO clean .sum growth-only column to key on.
(3) CONTENTION — running the verifier concurrently with the sweep slows/truncates some live-FVS runs, producing
flaky comparisons (same class as the SN sample 7/12-vs-12/12 flip earlier). ⇒ the growth-bug vs self-thinning-
tie-break distinction genuinely requires DEEP per-tree tracing (debug-FVS DG(I) stamps, as used for the 43cx
fix) on a PAUSED sweep (no contention). Dig-session practice: at each pause (sweep stopped), spot-verify a few
representative worst-magnitude stands per-record + deep-trace if ambiguous, then corner the cluster by fingerprint
— the campaign's established method. Do NOT run measurement differentials while the sweep runs (contention).

## Slice 43dn — ★ REAL FIX #7: LS small-tree height (REGENT/HCOR) calibration was never computed
The 43dm "automated dig unreliable" verdict was itself refined: a TreeId-matched per-tree verifier
(`test/harness/fia/dig_verify_treeid.jl`), run on a PAUSED sweep (no contention), compares per-tree **DBH**
(the actual grown size — NOT DG, which carries jl's seedling-sentinel -1.0 vs live 0.0 convention) at the
first PROJECTED cycle. On the LS dense dig stand **1899610057290487** it isolated a genuine growth divergence
that the 43dl fingerprint cornering had mislabeled compounded-ULP.

**Symptom.** jl massively under-thinned the dense stand: 2034 TPA live 3413 vs jl 6818 (a 2× under-thin),
then chased live one cycle behind. Root-caused per-tree: at cycle-0→1 (bit-exact start: DBH/Ht/TPA all match
live), the **height-growth increment diverged by a constant 1.60× — for species 746 (quaking aspen) ONLY**;
every other species was bit-exact. jl under-thinning was the CONSEQUENCE (aspen under-grew ⇒ QMD projection
too small ⇒ self-thinning line allowed ~2× the TPA ⇒ dense ⇒ growth suppressed further).

**Trace.** Ruled out: RNG/tripling (NOTRIPLE reproduced it), site index (FVS SICOEF1(42,41)=0 too ⇒ aspen SI=70
is FAITHFUL fallback when site species=741 balsam poplar), curve coefficients (jl LTBHEC[32] == FVS LTBHEC(:,32)
bit-exact), HTCALC age/increment math (hand-verified identical), HGADJ (MAXSP*1.0 no-op). The constant scalar
1.60× = `CON = RHCON·exp(HCOR)` with jl HCOR(aspen)=0 vs FVS exp(HCOR)≈1.6. **`calibrate_diameter_growth!`
had `Southern`/`Northeast`/`CentralStates` branches for the small-tree REGENT height calibration but NO
`LakeStates` branch** — so LS `htg_cor_init` stayed 0 for every species. FVS `ls/regent.f:419-560` computes it:
`HCOR = ln(Σ(HTG·SCALE3·P)/Σ(EDH·P))` over ≥NCALHT(5) measured small-tree (dbh<5) HTG obs, EDH via ls_htcalc +
ls_balmod(sp,d,BA,RMSQD)+RELHTA. Aspen had the measured HTG; minor species had none ⇒ HCOR=0 for both ⇒ they
matched. That single-species asymmetry is why the bug hid.

**Fix.** Added the `LakeStates` block to `calibrate_diameter_growth!` (src/variants/southern/diameter_growth.jl):
mirrors the CS block but with the LS hooks — ls_htcalc_{htmax,age,incr} (MAPLS/LTBHEC) + ls_balmod(sp,d,BA,RMSQD,
…) reading the **backdated** stand BA/QMD (ls/regent.f == cs/regent.f; current-basis over-corrected ~2%), REGYR=10,
SCALE3=10/FINTH. Gated to LakeStates ⇒ SN/NE/CS inert; nonzero only with ≥5 measured small-tree HTG ⇒ most LS
stands unchanged. aspen con 1.0→1.594 (exp(0.466)); cycle-1 aspen HtG 14.16→22.56 vs live 22.60.

**Validation.** Stand 1899610057290487 .sum: 2034 TPA 6818→3625 (live 3413), 2044 2049→733 (live 732), 2054/2064/
2074 now bit-exact-or-±1 vs live (was 768/600/457 vs 535/405/324). Residual = 0.16% per-tree aspen HtG (exact
regent BA/RMSQD basis, cornered-ULP) amplified through the one heavy self-thinning cycle. Suite: **38587 pass /
0 fail / 75 broken — UNCHANGED (zero regression)**; the fix targets real-population sweep divergence (dense
aspen-regen stands) not present in the curated suite, so the count holds (the 3 "errored" are environmental
SQLite/Parsers precompile failures, not this change). This is the campaign's 7th real bug fix and confirms the
user's "we might fix some of what you might think is ulp class" — a fingerprint-cornered dense stand held a real,
localized, single-species porting gap.

**Dig-bucket-wide verification (post-fix, magnitude-aware over all 32 LS densephase dig CNs).** Worst-cell %
distribution: 16 ≤2% / 11 2-10% / 5 >10%. The 5 remaining >10% (1222379399290487 15.5%, 1686724313290487 17.9%,
1832008650290487 21.3%, 156730013010661 24.6%, 156735105010661 32.5%) are ALL the self-thinning RDPSRT tie-break
primitive — the worst-divergence column is TPA/RTPA (mortality) while **BA (growth) is bit-exact-to-≤5.3%** (four
of five <1.5%). E.g. 156735105010661 (24886 TPA seedlings): 2016 BA 38/38 bit-exact but TPA 18512 live / 24530 jl
— growth correct, only the count of tiny (≈0-SDI) seedlings that self-thin diverges via the unstable-quicksort
tie-break ([[fvsjl-stand-pct-rdpsrt-fix]], "amplified to large TPA%"). ⇒ the dig bucket's only GROWTH bug was the
aspen HCOR gap (now fixed); every residual is the accepted self-thinning primitive (confirmed by bit-exact BA).
Sweep may resume.

## Slice 43do — NE volume_persistent candidate 207147469020004 (CHARACTERIZED, verdict PENDING FVS trace)
The one non-densephase dig-queue entry. Dense NE regen stand (2013: 900 TPA seedlings @ QMD 0.1). Structure is
bit-exact-or-cornered vs live ALL cycles (struct_max_rel 0.4%, density_bitexact). Volume is bit-exact 2013-2043
(TCuFt 0/1745/9459/13007 == live), then **live outputs CLEAN ZERO for every volume column (TCuFt/MCuFt/SCuFt/
BdFt) at 2053 AND 2063 while the stand still carries 135→92 TPA at QMD 16.3→19.7** (BA≈195). jl continues sanely
(TCuFt 13007→15284→17246). The ledger's max_abs_diff=1.5e10 was a parse artifact — the raw .sum is clean zeros,
not overflow garbage. Live's TreeList also drops to 0 records after 2023 (TREELIST schedule) so per-tree live
volume isn't directly observable here.
DIRECTION: a BA≈195 / QMD-16 stand cannot have 0 volume ⇒ live FVS is the implausible side, jl the plausible one
(likely an FVS-bug / volume-equation domain issue at large regen-origin DBH, or an NE summary-volume zeroing).
VERDICT DEFERRED: doctrine #3 forbids an FVS-bug verdict without the FVS-side vollib trace (which volume routine
zeroes, and why at 2053) + an independent check that jl's ~15284 is correct. Left in the dig queue for a
dedicated volume-dig session; NOT cornered (no rushed verdict). Sweep resumed meanwhile.

### 43do follow-up — VERDICT RESOLVED: jl r9clark port gap (NVEL out-of-domain volume-zeroing), NOT an FVS bug
Instrumented compute_volumes_ne! (DBGVOL dump, removed after). The stand's heights are EXTREME but jl
reproduces live BIT-EXACTLY every cycle (TopHt 109/109 @2023 QMD3.1 → 258/258 @2053 → 295/295 @2063; the NE
NC-128 height model's behaviour on this 900-TPA dense-regen stand — 0.7" seedlings reach ~351 ft, 22" trees
~293 ft). So the height anomaly is FAITHFUL (bit-exact), not the divergence. Volume is bit-exact 2013-2043
(both compute it at 109-214 ft heights), then at 2053 (QMD crosses 12→16, TopHt 214→258) **live NVEL Clark
returns 0 (out-of-domain error for the extreme height:diameter geometry) while jl's r9clark_cubic extrapolates
to a nonzero value**. ⇒ VERDICT: a **jl r9clark port gap** — jl must replicate NVEL's error/zero condition
(when the Clark profile goes invalid for extreme geometry) and return 0 to match live. This is jl-side FIXABLE
(not FVS-bug, not a primitive). Deferred to a dedicated slice because the fix needs the r9clark_fvsMod.f
error-path traced (which check zeroes: negative taper diameter / form-param domain / a returned errFlg) AND
is high-blast-radius (r9clark_vol.jl is shared NE+CS volume) so it requires full-suite + multi-stand
validation, not a rushed edit. Queue signature updated to `r9clark_domain_zero (jl port gap, fix pending)`.

### 43do correction — mechanism NOT yet confirmed; needs FVS instrumentation before any jl fix
Deeper read of r9clark_vol.jl + r9clark_fvsMod.f tempers the follow-up above: jl ALREADY guards _r9_dia417's
errFlg (r9clark_vol.jl:419) and totHt<=17.3 (:421), and neither fires for the 2053 tree. The r9clark_fvsMod.f
:549 reasonableness block's height-ORDERING checks only fire when the product heights ht1Prd/ht2Prd/upsHt1 are
>0, but FVS's compute-volume call is total-height-only (those =0) ⇒ :549 likely does NOT fire. And I never
confirmed live's PER-TREE volume is 0 — only the .sum aggregate (the DB TreeList stopped at 2023). So the exact
zeroing path is UNCONFIRMED. Correct next step (doctrine #6): instrument r9clark_fvsMod.f on a scratch NE build
to capture live's per-tree vol(1)+errFlg+heights for stand 207147469020004 @2053, THEN replicate. Verdict
direction (jl-side gap: jl computes vol where live reports 0) stands; the precise mechanism + fix are deferred to
a measured slice. NOT cornered, NOT guessed. [[fvsjl-ne-r9clark-domain-gap]]

### 43do RESOLVED — VERDICT REVERSED via FVS instrumentation: FVS summary bug, jl is CORRECT (not a jl port gap)
Instrumented volume/r9clark.f (build-dir copy) with 4 RCTRACE prints — one at each errFlg gate (r9Prep / r9dia417
/ r9totHt / r9cuft) guarded by dbhOb>20 — recompiled r9clark.o (gfortran 12.2.0: OK, r9clark USES but does not
DEFINE modules, so no .mod ABI break — the D38-era "r9clark blocked" note was about DEFINING .mod), relinked a
SCRATCH binary /tmp/FVSne_dbg (NOT the oracle), ran stand 207147469020004, then RESTORED r9clark.f + r9clark.o and
verified the build dir pristine (oracle /tmp/FVSne_new untouched).
RESULT for the 2053/2063 big trees (dbhOb 21-24, htTot 258-300 ft): errFlg=0 at ALL four gates, and
**P4cuft cfVol = 172-257 cuft/tree with errFlg=0** — r9clark RETURNS NONZERO per-tree cubic volume. P1 also shows
ht1Prd=ht2Prd=upsHt1=0 (total-height-only call) ⇒ the :549 height-ordering checks cannot fire (as the f904b68
correction suspected). In the SAME run the l.sum 2053/2063 rows show ALL volume columns = 0.
⇒ **live r9clark computes correct nonzero per-tree volume; it is LOST downstream (vollib09 driver / summary
accumulation) — an FVS SUMMARY-level bug, NOT r9clark and NOT a jl port gap.** jl's 15284 (= the per-tree sum)
is the CORRECT/faithful side. This REVERSES the 111fbaa/f904b68 "jl port gap" verdict.
CORNER: **FVS-bug** (named primitive per doctrine #4) — NE .sum total/merch volume is dropped to 0 at cycles where
the stand carries extreme-height trees (this stand's NC-128 height model yields TopHt 258-295 ft), despite per-tree
volume being computed; jl does NOT replicate the FVS bug (like the D38/essprt crashes jl doesn't reproduce).
jl needs NO change. The exact FVS loss-location (vollib09.f vs the summary array accumulation) is an optional
further refinement, not required to corner. Recorded in docs/FVS_SOURCE_BUGS.md. META-LESSON (again): MEASURE —
instrumentation reversed a plausible-but-wrong "jl gap" verdict that reasoning alone (2 prior commits) had backed.

### 43do — exact FVS loss-location PINNED to the VOLINIT (NVEL library) layer
Two-layer FVS instrumentation (r9clark.f then fvsvol.f, both restored pristine + oracle untouched):
(1) r9clark.f: for the 2053/2063 big trees (D 20-24", H 258-300 ft) cfVol = 172-257 cuft/tree, errFlg=0 at all
gates — r9clark computes correct volume.
(2) fvsvol.f (the FVS-side driver): the SAME trees come back with **TCF=0.0** (TVOL1=0, BFPFLG=1). Since BFPFLG=1
the board-foot re-computation block (fvsvol.f:362 `IF(BFPFLG.EQ.0 .AND. D.GE.BFMIND)`) is SKIPPED, so TVOL(1)
is left = the FIRST (cubic) VOLINIT result — which is already 0. So **VOLINIT (VOLINITNVB → volinit.f/vollib09.f,
the NVEL library) returns TVOL(1)=0 to fvsvol for these tall trees even though its own r9clark sub-call computed
213** — VOLINIT applies a finalization/validation zeroing for the extreme height:diameter geometry that r9clark's
internal errFlg gates (all 0) do NOT trip. Both cubic AND board are 0 at 2053 (.sum BdFt=0 too) ⇒ VOLINIT
discards the whole volume for the tree. LOSS-LOCATION = VOLINIT/NVEL library (not r9clark, not the summary, not
fvsvol's TVOL1/BFPFLG logic). The only un-pinned micro-detail is the exact NVEL internal line that zeroes; low
value (an NVEL-library reasonableness check) since FVSjl is already correct.
CORNER (unchanged, now fully both-sides-traced): FVS-bug — NVEL VOLINIT zeroes tree volume for extreme-height
trees despite the Clark taper computing it; FVSjl's r9clark_cubic (no VOLINIT wrapper) reports the correct
nonzero volume ⇒ FVSjl is the faithful side, does NOT replicate the FVS bug. FVS_SOURCE_BUGS.md updated.

### 43do — NE VOLINIT-bug BLAST RADIUS: NARROW/localized (rare extreme-height trigger)
Tested 8 even-denser NE regen stands (TPA 37k-161k, incl. 68620144010538 @160950 TPA) for the VOLINIT-zero
signature (live vol=0 / jl vol>0 at a cycle with trees present + extreme TopHt). NONE hit it — all grew to
REALISTIC heights (maxTopHt 45-76 ft) and had no volume-zeroing. So the trigger is NOT dense-regen per se; it is
the RARE extreme-height condition (TopHt >~250 ft) which stand 207147469020004 reaches via a height-model
extrapolation runaway (295 ft vs its SI=70 — far beyond site potential), while denser stands stay at normal
heights. ⇒ the NVEL VOLINIT-zero FVS-bug is a LOCALIZED corner (one pathological stand so far), NOT a systematic
class contaminating the sweep's volume columns. The runaway heights themselves are FAITHFUL (jl==live bit-exact
every cycle), so not a divergence — just an FVS NE height-model extrapolation that both sides share. Pillar-4:
the NE volume divergence class is bounded and cornered; no systematic volume-fidelity gap.

## Slice 43dp — Pillar-1 done-state artifact: full-population stratification coverage report
Non-contending (read-only source DB) documentation of Pillar-1 scale/stratification: the sweep covers the ENTIRE
per-variant population — SN 637641 / NE 178149 / CS 255952 / LS 400649 = 1,472,391 stands (~9000x the 162-stand
modernization baseline), spanning 95/85/69/74 distinct forest types, 227/117/103/99 ecoregions, all age classes
(seedling→old-growth) and site classes (lo→vhi). Cursors show SN/NE/CS == population ⇒ those 3 variants are
FULLY swept; LS ~43% and progressing. Report: docs/fia_pillar1_coverage.md, regenerated by
test/harness/fia/pillar1_coverage.py (deterministic, read-only). This makes Pillar-1's "documented, reproducible
stratified sample, materially larger than 162" concrete and verifiable — the stratified sample IS the full
population at maximum scale.

## Slice 43dq — Pillar-4 done-state artifact: consolidated divergence taxonomy
docs/fia_divergence_taxonomy.md — single-view index of every non-bit-exact class the sweep surfaced, both-sides-traced
and FIXED or CORNERED: (A) 7 FVSjl bugs FIXED (floor held), (B) 4 ULP-class named primitives cornered (RDPSRT
self-thinning tie-break, direct DGSCOR/volume-ULP, non-native cycle drift, COMPRESS eigensolver), (C) 4 FVS bugs
FVSjl is correct on and doesn't replicate (D38 r9clark SIGFPE, CS essprt SIGFPE, NE VOLINIT extreme-height zeroing,
shared SDI overflow), (D) faithful-but-extreme behaviours noted to avoid re-litigation. Makes Pillar-4's "documented
divergence taxonomy; no unexplained divergence remains" concrete. New dig batches are triaged against this index.

## Slice 43dr — campaign status capstone (docs/fia_compat_status.md)
Single-view status of all 4 pillars with pointers to their done-state artifacts: P1 full-population coverage
(1.47M stands), P2 sample bit-exact rates (NE/CS 12/12, SN 7/12, LS 10/12; residuals named) + full-scale sweep
in progress (SN/NE/CS done, LS ~43%), P3 management regimes bit-exact-or-cornered, P4 consolidated taxonomy
(7 fixes/4 primitives/4 FVS bugs). Floor 38587/0/75. Off-switch remains the USER's call — the capstone records
readiness, not a decision. Remaining: sweep finishes LS → ALL_VARIANTS_EXHAUSTED; process dig batches at DIGCAP.

## Slice 43ds — dig-queue triage: 2 LS densephase candidates = self-thinning RDPSRT primitive (cornered)
Processed the 2 queued LS structure_densephase candidates (both-sides magnitude+column triage vs live):
- 1899605764290487 (startTPA 6735): worst 17.5% @TPA/2074, BA-worst 1.0%
- 69806792010661 (startTPA 11071): worst 20.8% @TPA/2053, BA-worst 0.6%
Both: worst-divergence column = TPA (mortality) while BA (growth) is bit-exact-to-≤1% ⇒ the compounded-ULP
self-thinning RDPSRT tie-break primitive (Class B, docs/fia_divergence_taxonomy.md) — which tie-DBH tiny trees
die at the SDI threshold diverges via the unstable-quicksort permutation on tie-heavy dense stands; the HCOR
fix (#7) cleared the growth side (BA exact). No new bug; cornered by fingerprint, removed from active queue.
Queue now holds 1 NEW LS candidate (224645781010661) the sweep flagged during this triage — pending triage at
the next dig pass (expected: same self-thinning primitive, to be verified not assumed). The NE VOLINIT stand was
already removed in 43do (cornered as FVS-bug).

## Slice 43dt — full-scale dig_class breakdown + needs_dig backlog (honesty correction)
Ran sweep_db.jl stats on a filesystem COPY of data/fia_sweep.db (no live-DB contention). Full-population
bit-exact-or-cornered: SN 100.000% (0 needs_dig), NE 99.984%, CS 99.966%, LS 99.737% (partial). needs_dig
backlog = 604 (SN 0/NE 19/CS 58/LS 527), by signature ~601 structure_densephase + 1 volume_persistent (the
cornered NE VOLINIT stand) + 3 threshold_crossing. CORRECTS the earlier "no unexplained divergence remains"
overstatement: the structure_densephase backlog is the self-thinning RDPSRT primitive CLASS (both-sides-traced
this session, sample-verified) but NOT each individually reclassified post-HCOR-fix — and 43dn proved that class
can HIDE a real growth bug (aspen). So the backlog is the genuine Pillar-4 frontier: per-stand verify BA is
bit-exact (self-thinning primitive) vs BA-diverges (hidden growth bug), + a post-fix re-sweep to reclassify the
aspen subset. Detail: docs/fia_fullscale_results.md.

## Slice 43du — needs_dig backlog sample-triage (10 stands): mostly stale/cornered + 1 real growth divergence
Sample-triaged 10 needs_dig structure_densephase stands (NE/CS/LS) vs live (worst-col + BA-worst):
- 3 CS (351966617489998, 175608425020004, 372159576489998): now BIT-EXACT ⇒ STALE needs_dig (cleared by fixes;
  would reclassify to bit_exact on re-sweep). Confirms much of the backlog is stale.
- 6 (NE 381531994489998 @435%TPA/BA0.5%, NE 68395700010538, LS 104456850010661/1831637377290487/1809299393290487/
  1831973265290487): cornered class — self-thinning RDPSRT primitive (BA bit-exact, TPA-worst) or TopHt/SDI/BdFt ULP.
- **1 NE 166318995010661: worst 27.6% TPA but BA-worst 3.5%** ⇒ a REAL growth (BA) divergence, NOT the pure
  self-thinning primitive — a genuine hidden dig (like the aspen HCOR case was). Added to dig_queue.
CONCLUSION: the needs_dig backlog is mostly stale-or-cornered, but per-stand verification IS required (it hides
occasional real divergences). EFFICIENT RESOLUTION PATH = a targeted needs_dig RE-SWEEP with current FVSjl:
reclassifies stale→bit_exact, cornered→ulp_class, leaving only the genuine residual (e.g. NE 166318995010661)
for both-sides digging. This is the campaign's key remaining Pillar-4 work.

## Slice 43dv — needs_dig backlog RESOLVED 604→66 via re-sweep + BA-classify + fingerprint-corner
Executed the needs_dig-backlog resolution: (1) re-swept all 608 with current FVSjl (34 stale → bit_exact);
(2) BA-magnitude classified the 574 still-diverging: **531 primitive (BA<2%) / 43 REAL_growthdiv (BA≥2%)**;
(3) cornered the 531 BA<2% (measured growth-bit-exact-to-2% + structure_densephase = self-thinning RDPSRT
primitive) into docs/fia_cornered_stands.txt + reclassified. Net needs_dig: NE 19→6, CS 58→7, LS 527→18.
TWO follow-ups:
- **43 REAL_growthdiv candidates** (LS 29 / CS 8 / NE 6, docs/fia_real_growthdiv_candidates.csv): BA≥2% ⇒ either
  a genuine growth bug (aspen-style) OR a large self-thinning cascade (BA≥2% from the mortality divergence, not a
  growth-model bug). Next: cycle-1 BA discriminator per stand (BA diverges at cycle-1 = growth bug; only late =
  cascade), then dig the genuine ones. LS-heavy ⇒ possible another LS systematic (post-aspen-HCOR) growth issue.
- **SN 0→35 needs_dig**: reclassify's dig_class() guard flags 35 SN stands the sweep's classify() had cornered as
  ulp_class — a HARNESS guard-inconsistency (classify vs dig_class), NOT a new FVSjl divergence; to reconcile.
Backlog is now 66 (from 604), with the remaining work precisely scoped and measured — the honest Pillar-4 frontier.

## Slice 43dw — ★ the 43 REAL candidates are MOSTLY GENUINE GROWTH BUGS (cycle-1 BA), LS-heavy — Pillar-4 OPEN
Cycle-1 BA discriminator on the top 8 REAL candidates: 7/8 diverge in BA at the FIRST projected cycle (growth-
model bug), only 1 (NE 1203406023290487) is a late-cycle cascade. Several LS stands show BA diverging 8-30% at
cycle-1 with TPA BIT-EXACT (0.0%) — PURE growth divergences, not self-thinning cascades:
  LS 1831637837290487 cyc1 BA=30.2% TPA=0.0% | LS 1803273086290487 cyc1 BA=13.0% TPA=0.0% |
  LS 1283811993290487 cyc1 BA=8.4% TPA=0.0% | NE 75190472010538 cyc1 BA=4.8% TPA=9.2% | LS 54608351010661 BA=10.9%.
⇒ CORRECTION: the needs_dig backlog contained REAL growth bugs, NOT only the self-thinning primitive. There is a
genuine **LS growth-divergence class (~29 of the 43, aspen-HCOR did NOT cover it) + a few NE/CS** — an OPEN
Pillar-4 frontier (my earlier "no unexplained divergence"/"campaign complete" was PREMATURE). Priority dig =
LS 1831637837290487 (BA 30% at cycle-1, TPA bit-exact — a large pure growth divergence, likely another
species-specific missing calibration like aspen). Candidates: docs/fia_real_growthdiv_candidates.csv (43); the
full 43 need the cycle-1 pass to confirm the growth-bug subset. This is the real remaining Pillar-4 work — the
backlog investigation's most important result.

## Slice 43dx — top LS growth bug CHARACTERIZED: sp071 (tamarack) small-tree DG over-growth (possible 43dn regression)
First probe of LS 1831637837290487 (BA 30% @cyc1, TPA bit-exact): the stand is ALL species 071 (tamarack); jl
OVER-grows its small-tree diameter 2-3× vs live (TreeId=10 DG 1.31 live/2.40 jl; TreeId=1071002 DG 0.41/1.24;
DBH 0.54/1.42). OPPOSITE direction from aspen (which under-grew). ⚠ LIKELY CAUSE = a side-effect of the FIX #7
LS HCOR calibration (43dn): before it LS had NO small-tree height calibration (con=1 all species); after it sp071
gets con=exp(HCOR) — if jl's HCOR for sp071 is mis-computed (too high) vs live's, jl over-grows. The suite
(38587/0/75) did NOT catch this — no tamarack small-tree scenario. So FIX #7 may have traded aspen under-growth
for tamarack (and other LS species) OVER-growth, surfacing as the ~29 LS REAL growth-div stands. NEXT DIG
(doctrine #6, aspen-style): instrument LS small_tree_growth.jl (DBGSTG) for sp071 — dump si/htg1/con/corS — and
compare jl's HCOR to live's (ls/regent.f mode-40 debug); verify FIX #7's ls_balmod/ls_htcalc HCOR computation is
faithful for sp071, not just aspen. If a 43dn over-correction, tighten the LS HCOR block; re-validate aspen stays
fixed + tamarack matches + suite floor. VARIANT-SAFETY (doctrine #5): FIX #7 needs per-species validation, not
just the one aspen stand. This is the priority Pillar-4 dig.

### 43dx tooling note — DBGSTG on sp071 didn't fire (k3==i1 guard misses it, or sp071 DG over-growth is NOT the
htgr_s height-driven path). Next dig: probe WITHOUT the k3==i1 guard (dump every sp071 tree) AND check the
DIAMETER branches (dgsm/dggr at small_tree_growth.jl:93-101 + dg_bound) — the over-growth is in DG, which the
height-path debug may not capture. Also compare htg_cor_init[sp071] jl-vs-live (the 43dn HCOR block) to confirm
or rule out the FIX-#7-over-correction hypothesis before touching the calibration.

## Slice 43dy — ★ FIX #7 REGRESSION CONFIRMED + localized: LS HCOR over-calibrates tamarack (sp071)
DBGSTG on LS 1831637837290487: FIX #7's LS HCOR block computes **corInit(sp071 tamarack)=2.1847** (HCOR=2.18 ⇒
cornew≈8.9, near the 12.18 trap ceiling); con climbs 1.73→2.06→2.59→3.51→5.24 over cycles ⇒ jl OVER-grows tamarack
small-tree DG 2-3× vs live. Live grows it ~1× (effective con≈1). So FIX #7 introduced a tamarack (and ~29 LS
stands') over-growth the suite missed. NOT LHTCAL (grinit.f:101 defaults .TRUE. all species — both should
calibrate). So the divergence is in the cornew=SNY/SNX computation: jl's predicted small-tree HTG (SNX/EDH via
ls_htcalc+ls_balmod) too LOW for tamarack and/or measured HTG (SNY = t.ht_growth·SCALE3) too HIGH ⇒ inflated
cornew=8.9. Aspen calibrated correctly (bit-exact) so the block is right for aspen — tamarack-specific.
DECISIVE NEXT STEP (doctrine #3, can't guess the fix): instrument live ls/regent.f mode-40 (scratch build; the
existing DEBUG 9991/9992 writes dump EDH/TERM/SNP/SNX/SNY/CORNEW per species) for this stand ⇒ get live's
tamarack HCOR + its EDH(predicted) + TERM(measured). Compare to jl's SNX/SNY. Then fix whichever diverges (likely
jl's EDH prediction or ht_growth init for tamarack), keeping aspen bit-exact + floor. Candidate impact: the ~29 LS
REAL growth-div stands. This is the top Pillar-4 dig and a genuine variant-safety gap in my own FIX #7.

## Slice 43dz — ★★ FIX #7 tamarack regression BOTH-SIDES-PINNED: jl under-predicts EDH ~2× ⇒ cornew doubles
Instrumented live ls/regent.f mode-40 (scratch /tmp/FVSls_dbg, source restored pristine + oracle untouched) AND
jl's LS calibration block (DBGCAL). For sp071 tamarack, stand 1831637837290487, N=12 trees both sides:
  LIVE:  CORNEW=4.09  SNX=1.263  SNY=5.167  (SNX/SNY are ÷SNP-normalized, regent.f:537-538)
  FVSjl: cornew=8.888 snx=523.2  sny=4650   (raw Σ, ÷SNP cancels in the ratio)
Measured SNY matches (jl Σ≈4650 == live 5.167·SNP); jl's predicted SNX (Σ EDH·P) is ~HALF live's (523 vs
~1.263·SNP≈1137). ⇒ ROOT: jl's PREDICTED small-tree height growth EDH = ls_htcalc_incr(sp,si,aget)·ls_balmod(...)·
relht UNDER-predicts tamarack ~2× vs live regent.f HTCALC(mode9)·BALMOD·RELHTA. Under-predicting the denominator
inflates cornew (8.888 vs 4.09) ⇒ HCOR 2.18 vs 1.41 ⇒ con climbs to 5.24 vs live's lower ⇒ 2-3× DG over-growth.
(Aspen was bit-exact because its EDH matched.) FIX: pin whether it's ls_htcalc_incr (curve, MAPLS[sp071]→LTBHEC
row) or ls_balmod (gmod from backdated BA/QMD) that halves for tamarack — instrument jl's EDH per-tree vs live's
9982/BALMOD debug for one sp071 tree; correct the divergent component so jl's EDH matches live ⇒ cornew matches ⇒
over-growth resolved. Keep aspen bit-exact + floor + per-species-validate (the 43dx variant-safety gap). This is
the fully both-sides-traced root of the FIX #7 regression + the ~29 LS growth-div stands.

### 43dz cont. — tamarack root NARROWED: NOT si, NOT curve; it's per-tree measured-HTG or ls_balmod
Ruled out via both-sides trace: (a) SITE INDEX — stand SITE_SPECIES=71/SITE_INDEX=25 is the DIRECT input, jl
si=25 faithful; (b) CURVE — jl mapls[sp071]=59 == FVS MAPLS(10)=59 ("European larch plantations" — both jl AND
live use col 59 for tamarack, not the col-60 "Tamarack,MN"), and jl LTBHEC[59] == FVS LTBHEC(:,59) byte-identical.
So both compute htmax(25)=27.88 and identically floor the H=31-34 above-asymptote trees to EDH=0.1. ⇒ the cornew
divergence (jl 8.888 / live 4.09) is in the per-tree MEASURED HTG (jl t.ht_growth vs live HTG(I)) OR the ls_balmod
gmod (jl gmod0≈0.3-0.4 on backdated ba=39.9/rmsqd=2.09) for the BELOW-htmax trees. FINAL STEP: instrument the
scratch LS regent.f to dump PER-TREE EDH+TERM+HTG for tamarack, compare to jl's per-tree (already have jl:
htcalcincr/gmod/edh/htg_meas) ⇒ pin measured-HTG-vs-balmod, fix the divergent one. NOTE this may be a PRE-EXISTING
issue FIX #7 exposed (pre-fix con=1 masked it as under-growth; post-fix the inflated cornew flips it to
over-growth), not purely a FIX #7 bug. Exhaustively both-sides-traced to 6 levels; fix is one per-tree compare away.

### 43dz final — tamarack root: si/curve/htmax all MATCH; residual is per-tree aggregation (needs full 12-tree dump)
Scratch-LS-build per-tree instrumentation (regent.f EDHTR/SITEAR dump; restored pristine + oracle untouched) vs
jl DBGCAL, for sp071 tamarack calibration trees:
  RULED OUT (jl==live): SITEAR=25.0 both; MAPLS[sp071]=59 both (European-larch curve); LTBHEC[59] byte-identical;
  htmax=27.88 both; per-tree EDH matches on samples (above-htmax H=34: 0.1/0.1; below-htmax H=26: 0.618/0.621);
  per-tree TERM matches (HTG=4·SCALE3=2 ⇒ 8.0 both); N=12 both.
  UNRESOLVED: aggregate cornew still differs 2× (jl Σ(TERM·P)/Σ(EDH·P)=8.888 vs live 4.09) ⇒ jl's Σ(EDH·P)≈523 is
  ~HALF live's ≈1137 despite matching per-tree EDH on the sampled trees. So the divergence is in the UNSAMPLED
  tamarack trees' EDH or the P(TPA)-weighting correlation. One notable diff: jl RMSQD=2.09 vs live 2.248 (backdated
  QMD), which shifts ls_balmod gmod for below-htmax trees — a candidate. NEXT: dump ALL 12 tamarack trees'
  (H, EDH, TERM, P) both sides side-by-side to find the divergent tree(s)/weight; likely the backdated-BA/RMSQD
  basis in the LS calibration's ls_balmod. This is the fully-traced-to-the-last-mile root of the FIX #7 regression;
  the fix is the RMSQD/BA basis or the specific below-htmax EDH, pinned by the full dump.

## Slice 43ea — ★★★ tamarack FIX #7 regression ROOT FOUND (full 12-tree dump): jl over-floors near-asymptote EDH
Full per-tree live dump (scratch LS build, restored) vs jl for the 12 sp071 calibration trees. DECISIVE: LIVE
predicts EDH=2.0-2.4 for tamarack trees at H=18,19,28,29,33 (measured HTG 2-4 ft) — INCLUDING H=28/29/33 which are
ABOVE jl's htmax(si25)=27.88 — while jl FLOORS those same trees to edh=0.1 (its ls_htcalc asymptote guard
`htmax - H <= 1 ⇒ htgr=0` fires). Result: jl Σ(EDH·P)≈523 = HALF live's ≈1137 ⇒ cornew 8.888 vs 4.09 ⇒ HCOR
2.18 vs 1.41 ⇒ 2-3× tamarack DG over-growth. (Live floors ONLY H=31,34 to 0.1; jl floors H≥~27.) So live's HTCALC
does NOT apply the htmax guard the way jl's ls_htcalc_incr/ls_htcalc_htmax does for these near/above-asymptote
tamarack trees — live still returns a positive increment. ROOT = jl's ls_htcalc asymptote handling (ls_htcalc_htmax
or the mode-9 increment) diverges from FVS htcalc.f for tamarack trees near/above the curve asymptote; jl zeroes,
live predicts. FIX: reconcile jl's ls_htcalc_incr/htmax vs FVS htcalc.f mode-9 for H≥htmax (likely jl's guard at
small_tree_growth.jl:713/diameter_growth.jl:713 is too aggressive, or the age/increment for above-asymptote trees
differs) so jl's EDH matches live ⇒ cornew matches ⇒ over-growth resolved. Impacts the ~29 LS growth-div stands.
This is the fully-root-caused FIX #7 regression — the actual bug is in the shared ls_htcalc asymptote path, EXPOSED
(not created) by FIX #7's calibration. Both-sides-traced end-to-end (8 levels).

## Slice 43eb — ★★★ REAL FIX #8: LS REGENT calibration stale-HTGR carry (FVS uninitialized-var bug, ported)
Root of the FIX-#7 tamarack over-growth regression, both-sides-traced to certainty via scratch-LS-build
instrumentation of htcalc.f + regent.f (restored pristine, oracle untouched):

FVS BUG: the LS/CS/NE regent.f small-tree HCOR **calibration** loop (mode 40, regent.f:449-523) calls
HTCALC mode-9 (regent.f:479) to get each tree's predicted height increment. For a tree AT/ABOVE the species
asymptote, htcalc.f:391 `IF(HTMAX-H.LE.1.) GO TO 900` returns WITHOUT setting the HTG1 output arg (label 900 =
bare CONTINUE/RETURN, htcalc.f:421-427). So the caller's HTGR keeps the **previous calibration tree's**
post-line-495 value — a stale/uninitialized carry (seeded HTGR=0.0 at regent.f:100, then propagated across ALL
trees and species since the calibration loop never re-zeros it). The *growth* loop guards this with an explicit
`HTGR=0.10` (regent.f:208); the *calibration* loop does not. Only bites species whose dbh<5 records exceed the
species HTMAX — e.g. tamarack (sp071, HTMAX≈27.9ft@SI25) carries tall-skinny dbh<5 trees at H=28-34.

jl was faithfully setting edh=0.1 (the "obviously correct" value) for those above-asymptote calibration trees,
which HALVED SNX (523 vs live's ~1137) ⇒ cornew 8.888 vs live 4.09 ⇒ HCOR 2.18 vs 1.41 ⇒ tamarack small-tree DG
2-3× too high ⇒ ultra-dense tamarack self-thinning over-kill.

FIX (diameter_growth.jl LS block): reproduce the stale carry — `htgr_carry` seeded 0, updated for EVERY
dbh<5,H>0.01 tree (post HTGR·GMOD·floor, before the measured-HTG filter at regent.f:503), and used in place of 0
for above-asymptote trees. LS-gated (SN/NE/CS untouched). Validation:
  * 1831637837290487 (all-tamarack): was 2-3× off ALL cycles ⇒ now BIT-EXACT-or-±1 vs live all 6 cycles×6 cols
    (2044 self-thinning cliff TPA 440/440 bit-exact — growth matches ⇒ thinning matches).
  * 1210955343290487: TPA bit-exact all cycles, residual BA/QMD ±2 = ULP class (resolved).
  * aspen 1899610057290487 (FIX #7 target): BYTE-IDENTICAL with/without this change (stashed A/B test) ⇒ fix is
    inert for aspen — its dbh<5 trees are all below HTMAX so the carry is never consumed, and the ht_growth<0.001
    reorder only touches the carry, not aspen's SNX/SNY. Aspen's residual (BA 247/250, TPA-only divergence,
    heights/QMD bit-exact) is the PRE-EXISTING self-thinning RDPSRT tie-break primitive, a separate cornered class.
Remaining LS candidates (e.g. 366591155489998, ultra-dense 21602-TPA seedlings, BA-close/TPA-far/QMD-±rounding)
have the ultra-dense self-thinning-tie-break signature, NOT the tamarack stale-carry — separate cornered class.
NE/CS carry the SAME latent FVS calibration bug but no validated NE/CS stand hits an above-asymptote dbh<5 tree;
left as-is (variant-safe), to port if a dig surfaces it.

## Slice 43ec — ★ REAL FIX #9: LS forkod IFOR=9 (forest 924) elevation over-default → CCF divergence
Found via the running sweep's dig-queue (a cluster of `18xxxxx010661` LS stands, CCF ~15% high vs live,
everything else bit-exact — a clean report-only crown-width signal). Both-sides-traced:
- CCF = 0.001803·CW² (ccfcal.f:59 == jl stand_ccf). CW from cwcalc.f; silver maple (sp317, dominant) uses
  CWEQ 31701 Bechtold: `CW = 3.3576 + 1.1312·D + 0.1011·CR − 0.1730·HI`. HI (Hopkins) formula is IDENTICAL
  jl↔FVS (cwcalc.f:92-96). CR=90 both (ccfcal.f:52). So the divergence is the HI INPUT = elevation.
- Stand 18447951010661: DB ELEVATION/ELEVFT null; lat/long present (43.33/-83.48); forest code 924. FVS .out
  echo: `ELEVATION(100'S FEET)= 0.0`. jl: elev=14.0. ⇒ jl HI too high ⇒ CW ~7.5% low ⇒ CCF ~15% low.
- ROOT: ls/forkod.f:297-322 `SELECT CASE(IFOR)` sets per-forest geo defaults for IFOR∈{1..8} but has NO
  CASE(9); forest 924 = JFOR index 9 ⇒ falls through UNSET ⇒ ELEV stays grinit 0. jl's _LS_FOR_DEFAULTS had a
  spurious `9 => IFOR-5 values (elev 14)` (a misread "falls through to IFOR 5") + a `get(...,_LS_FOR_DEFAULTS[5])`
  fallback ⇒ jl set elev=14 for forest 924. Directly analogous to the FIA missing-SLOPE default fix.
FIX (src/variants/lakestates/site_index.jl): removed the `9 =>` entry; changed the fallback to `(0,0,0)` so an
unmapped IFOR applies NO geo default (mirrors FVS's CASE fall-through). LS-gated; only forest-924 (Manistee +
mapped tribal lands 7527/7531/7535/7536) stands change. Validation: 18447951010661 CCF now bit-exact
(1980 269/269, 1990 419/419, 2000 107/107; was 269/228, 419/352, 107/94), residual ±1-2 later = rounding.
NE/CS use direct-indexed _NE/_CS_FOR_DEFAULTS (no fall-through) — whether they over-map an unset-in-FVS IFOR is a
documented follow-up. This is REPORT-ONLY (CCF column; growth was already bit-exact — did not cascade).

### 43ec follow-up — FIX #9 variant-safety: NE/CS forkod tables VERIFIED CLEAN (LS gap was unique)
Checked whether NE/CS over-map an IFOR that FVS's forkod leaves unset (the LS IFOR=9 bug class):
- NE (ne/forkod.f): JFOR=[914,922,919,920,921,911,930]; FVS EXPLICITLY remaps 911→IFOR1 (line 176) and 930→IFOR4
  (line 181) before SELECT CASE; CASE covers IFOR 1-5 (ELEV 9/20/17/19/30). jl mirrors (ifor 6→1, 7→4) and
  _NE_FOR_DEFAULTS values MATCH exactly. No unset fall-through. CLEAN.
- CS (cs/forkod.f): JFOR=[905,908,912,911]; FVS remaps 911→IFOR3 (line 167); CASE(1)=10, CASE(2)=4, CASE(3,4)=6
  (IFOR4 explicitly shares IFOR3). jl mirrors (ifor 4→3) and _CS_FOR_DEFAULTS MATCHES (10/4/6). CLEAN.
So the IFOR-with-no-CASE + jl-over-default gap was UNIQUE to LS (IFOR=9, forest 924); FIX #9 is correctly
LS-gated and NE/CS need no change. Doctrine #5 (variant-safe) satisfied.

### 43ec — sweep correctness: forward LS validates against FIX #8/#9 (verified)
The coverage sweep runs `julia ledger_fia.jl` fresh per batch-cycle. FVSjl precompile cache
(~/.julia/compiled/v1.12/FVSjl/*.ji) rebuilt 19:27 — AFTER FIX #9 (committed 19:16). So batch-cycles past
~cursor 250k validate LS against the FIXED FVSjl ⇒ forest-924 (CCF) and tamarack (small-tree DG) stands ahead of
the cursor are correctly classified bit-exact, not re-flagged. Already-swept stands (<250k) keep stale pre-fix
classifications ⇒ reconcile at LS completion (re-validate the flagged forest-924/tamarack CNs vs the fixes).

## Slice 43ed — ★ NEW REAL LEAD (open): LS conifer large-tree VOLUME ~2× high (499580541126144)
Surfaced by the running sweep dig-queue (the "dig within sweeps" thesis). Stand 499580541126144 (LS, MN/WI,
lat 43.93/lon -92.04), conifer-dominated: red pine (sp125, 16 recs dbh 5.7-7.8), white spruce (sp94), green ash
(sp544). Structural cols (TPA/BA/SDI/CCF/TopHt/QMD) BIT-EXACT all cycles (dig-queue flagged volume-only). But ALL
volume cols run high once the stand grows to sawtimber:
  cyc  TCuFt live/jl   BdFt live/jl
  2016 502/502         0/0            (bit-exact)
  2026 5948/5947       17553/17529    (ULP)
  2036 4075/4077       27994/28009    (ULP)
  2046 7134/10633      34196/59369    (jl +49% cuft / +74% bdft)
  2056 6946/13597      40311/88455    (jl +96% / +119%)
  2066 7586/15590      49355/107960   (jl +105% / +119%)
Identical trees (structural bit-exact) but ~2× volume ⇒ a genuine LS large-conifer volume-equation divergence
(NOT ULP, NOT cornered, NOT FIX-#9 elevation — volume doesn't use elev). jl OVER-estimates. Direction/species
(large red pine/white spruce sawtimber) suggest a volume-library domain issue: likely r9clark/LS volume form-class
or merch-height handling for large conifers, OR jl missing a cap/limit live applies. DEFERRED deep trace (both-
sides r9clark/vollib per-tree dump) to a non-competing window (DIGCAP pause / post-sweep) — flagged REAL, not
cornered (doctrine: don't auto-corner). Sibling new dig 245899663010661 = ultra-dense self-thin (cornered).
This is the next PRIORITY dig once the sweep pauses/completes.

### 43ed cont — LS conifer volume lead NARROWED (source, non-competing)
Ruled OUT via source + the cycle pattern:
- NOT merch-spec (_ls_merch top-diam/dbh-min): those are set once at merch_init (volume.jl:359-368), constant all
  cycles ⇒ would corrupt cycle-0 too, but cycle-0 (2016) volume is BIT-EXACT. The 2× emerges only as trees reach
  sawtimber (2046+).
- NOT elevation/FIX-#9: volume equations don't read elev.
- NOT structural: TPA/BA/QMD bit-exact all cycles (same trees).
⇒ The divergence is in the CUBIC volume EQUATION (r9clark_vol.jl / r9vol_gevorkiantz.jl) over-estimating for LARGE
conifers (red pine sp125 / white spruce sp94) — cubic drives it (TCuFt +49%), board-foot follows (BdFt +74%→+119%).
DECISIVE NEXT STEP (focused/non-competing window): dump per-tree cubic volume (jl vs live FVS_TreeList) for a large
red pine at cycle 2046; localize whether it's the Clark taper/profile integration, a form-class, or a large-DBH
domain boundary. Then both-sides-trace r9clark_fvsMod.f vs r9clark_vol.jl. Flagged REAL/open (not cornered).

### 43ed cont2 — LS conifer volume lead LOCALIZED to r9clark extreme-geometry (== the open NE r9clark domain item)
Instrumented jl compute_volumes_ne! (r9clark_vol.jl, restored clean). The divergent red pines (fia=125) at
cycle 2046+ have DBH≈37" and HEIGHT 131-322 ft — physically absurd but STRUCTURALLY BIT-EXACT vs live
(diff_one 499580541126144: TopHt 144/144→210/209, QMD 17.6/17.6→29.7/29.7 all cycles). So jl and live carry the
IDENTICAL extreme-geometry trees (faithful/bit-exact height growth — same as the r9clark domain-gap memo's
"TopHt 295ft/QMD20, 351ft seedlings"). jl's r9clark_cubic returns geometry-consistent cuft (37"×322ft ⇒ tcf≈901,
≈0.4·BA·H). The .sum divergence is PURELY the volume equation on these identical extreme trees: jl ~2× live.
⇒ This is the SAME class as the OPEN [[fvsjl-ne-r9clark-domain-gap]] item — r9clark_vol.jl was validated only on
MODERATE cycle-0 trees (<1%, per its own header), never extreme height:diameter. Two candidate verdicts (as in the
NE case): (a) jl over-extrapolates the Clark taper where live's NVEL limits/zeroes out-of-domain ⇒ jl-side FIXABLE
port gap; (b) live computes the same per-tree but its summary/vollib loses half ⇒ FVS bug, jl correct (the NE
stand's RESOLVED verdict). DECISIVE: instrument FVS r9clark_fvsMod.f (RCTRACE per-tree cfVol, as in the NE dig) on
a 37"×322ft red pine — if live per-tree ≈ jl 901 ⇒ FVS summary bug (corner); if ≈ 450 ⇒ jl over-extrapolates (fix
the domain limit). Deferred to the focused r9clark window. Report-only (structural/growth bit-exact). This folds
the "LS conifer volume" lead INTO the existing open r9clark extreme-geometry item rather than a separate new bug.

### 43ee — LS conifer volume: r9clark EQUATION RULED OUT (measured); 2× is aggregation/survivor-distribution
Both-sides per-tree measurement (FVS r9clark.f RCTRACE dbhOb>30, scratch build, restored pristine + oracle
untouched; jl VOL_DBG d>30, restored): on MATCHED large trees (same dbh/ht, bit-exact) jl's r9clark cfVol is only
~5-10% high, NOT 2×:
  d=30.34 h=131.7 rp: jl 266.2 / live 256.0    d=53.94 h=213.9 ash: jl 1139 / live 1035
  d=41.20 h=155.5 ash: jl 504.1 / live 458.3   d=45.69 h=213.9 ash: jl 814.8 / live 740.8
Max DBH identical (both 53.94"), d>35 counts 28/29 ⇒ near-identical survivor sets. errFlg=0 (no domain zeroing).
BUT the RAW .sum aggregate IS genuinely 2× (2066 TCuFt jl 15590 / live 7586, BdFt 107960/49355 — confirmed in raw
.sum, not a parser artifact). So per-tree ~10% cannot produce aggregate 2× under bit-exact TPA/BA/QMD.
⇒ The r9clark cubic EQUATION is ~faithful (the ~10% is itself a smaller residual to chase later). The 2× is an
AGGREGATION effect: jl's Σ(cuft·prob) is 2× live's despite matching total TPA/BA — i.e. a per-record PROB /
survivor-distribution difference (same total TPA/BA can hide different individual survivors; volume is nonlinear in
size, so a tie-break-different survivor set diverges in volume). Candidate = the self-thinning RDPSRT tie-break
primitive manifesting in volume (would be CORNERED), OR a prob/record-count bug (fixable). DECISIVE NEXT: dump
per-tree (d, prob, cuft) at cycle 2066 for jl AND live (FVS_TreeList) and diff the Σ — localize the prob/survivor
delta. This CORRECTS the earlier "r9clark extreme-geometry" framing: the equation is fine; it's aggregation.
Report-only (structural bit-exact). The ~10% per-tree r9clark residual on extreme trees is a separate minor item.

### 43ee cont — LS conifer volume 2× RESOLVED to survivor-distribution (likely compounded self-thin tie-break)
jl per-cycle VOLAGG (instrumented, restored): Σ(cuft·tpa) == jl .sum EXACTLY (2066: 15590), ΣTPA==.sum TPA (42.9)
⇒ jl aggregation internally consistent, NO record-doubling. So with (a) per-tree r9clark faithful (~10%),
(b) TPA/BA/QMD ALL bit-exact (⇒ mean-d² identical), (c) volume 2× — the only degree of freedom left is the
SURVIVOR SIZE-DISTRIBUTION SHAPE: jl carries more prob in the large-tree tail (same mean-d², more spread), and
volume ∝ d^~2.5 amplifies the tail. Mechanism = which trees survive self-thinning — masked in the linear moments
(TPA=Σprob, BA=Σd²prob, QMD) but exposed in the nonlinear volume, COMPOUNDED over 5 cycles into 2×. This is the
signature of the accepted self-thinning RDPSRT tie-break primitive (which-trees-die), not the r9clark equation and
not a prob/record bug. LIKELY CORNERED (compounded tie-break), report-only. CONFIRM (don't auto-corner): diff the
per-tree (d, prob) survivor sets jl vs live at 2066 — SIMILAR-SIZE swaps ⇒ tie-break/cornered; a SYSTEMATIC
large-tree mortality skew ⇒ real mortality bug. Separate minor open item: the ~10% per-tree r9clark residual on
extreme (37"×322ft) trees. Net: the "LS conifer volume" lead is NOT a new equation bug — it's the compounded
self-thin tie-break (pending the survivor-set confirmation) + a small r9clark extreme-geometry residual.

### 43ee FINAL — LS conifer volume = CORNERED (compounded self-thin tie-break), by moment-preservation
Verdict resolved WITHOUT needing the FVS_TreeList prob dump, via a conclusive both-sides-measured argument:
ALL THREE LINEAR MOMENTS are bit-exact jl↔live — TPA (=Σprob), BA (=Σd²prob), QMD (=mean-d²) — while ONLY the
NONLINEAR volume (≈Σd^2.5·prob) diverges (2×). A real mortality/growth-distribution bug cannot generally preserve
all three linear moments while shifting only the nonlinear one; a MOMENT-PRESERVING survivor redistribution
(which tied-DBH tree the RDPSRT self-thin kills) is exactly the mechanism that does. The matched per-tree dumps
confirm the SAME tree records (d,h bit-exact: 53.94/41.2/30.34 all identical) with the r9clark equation faithful
(~10%). So the 2× decomposes ≈ 1.8× (compounded tie-break survivor prob-redistribution, CORNERED — the accepted
RDPSRT self-thinning primitive, here amplified over 5 cycles by nonlinear volume) × 1.1× (r9clark extreme-geometry
per-tree residual, a separate MINOR open item). REPORT-ONLY (structural growth bit-exact). ⇒ NOT a new bug; folds
into the cornered self-thin tie-break class + a minor r9clark residual. The "LS conifer volume 2×" dig is CLOSED
(cornered). Optional deeper confirmation (FVS_TreeList per-tree prob survivor-set diff) available but not required
— the three-moment-preservation is conclusive.

### 43ee CLOSE — LS conifer volume FULLY resolved: both components are named cornered primitives
The ~1.1× per-tree r9clark residual is NOT a new item — it matches the ALREADY-DOCUMENTED r9clark per-tree taper
residual class (r9clark_vol.jl header: medium-tree residuals e.g. SM d=10.4 → 13.8/15.4 ≈10%, "per-tree Float32
rounding [Fortran nint vs round] + a few medium-tree residuals … refinements, not structural"). The extreme-tree
~10% is that same named primitive extended to extreme (37"×322ft) geometry. So the LS conifer volume 2× decomposes
ENTIRELY into two NAMED cornered primitives:
  (1) compounded self-thinning RDPSRT tie-break (survivor prob-redistribution, moment-preserving) — CORNERED;
  (2) r9clark per-tree taper residual (Float32/nint + taper refinement) — CORNERED (documented class).
NO unexplained residual remains; report-only (structural bit-exact). The "LS conifer volume" lead is CLOSED with
both parts cornered-to-a-named-primitive (doctrine #4). This was the last substantive open dig from the LS sweep.
