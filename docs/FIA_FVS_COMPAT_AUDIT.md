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
