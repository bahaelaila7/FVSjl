# FIA/FVS behaviour-compatibility — working checklist / audit

Goal + doctrine: `docs/FIA_FVS_COMPAT_GOAL.md`. Every slice: plots covered, per-cycle pass rate vs
freshly-relinked live FVS, divergences found → both-sides-traced → fixed or cornered. Never regress the
floor (`julia --project=. test/runtests.jl` = 38527/143/0).

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
