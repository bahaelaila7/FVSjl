# Management & disturbance completeness tracker (C1-C5)

**Comprehensive** map of every management / disturbance / keyword-option in C1-C5,
built by classifying the **full 140-keyword master table** (`base/keywds.f`) — not a
hand-picked subset. The earlier version of this file missed ~18 keywords; this one is
exhaustive against keywds.f. Disturbance models that are whole extensions (FFE fire,
insects) are C7; econ is C8 — noted here for completeness but owned there.

> **DO NOT START until [`NATURAL_PROCESS_TODO.md`](NATURAL_PROCESS_TODO.md) is done**
> (it is: dynamics complete, .sum bit-exact 89/90, classification all 90). Same
> discipline: port the real FVS code; tests only catch bugs. C6-coupled items live in
> [`C6_DBS_TODO.md`](C6_DBS_TODO.md).

Legend: ✅ done · 🟡 partial · ⛔ unported · ⚪ N/A in SN · 🧊 C7/C8 extension

## 1. CUTS — thinning / harvest methods (`cuts.f`, ICFLAG)

| keyword | semantics | status |
|---|---|---|
| THINBTA/THINATA/THINBBA/THINABA | from below/above to residual TPA/BA | ✅ |
| THINDBH | proportional DBH-class to residual TPA/BA | ✅ |
| THINAUTO | auto-thin to FULSTK on stocking trigger (recurring) | ✅ (±1-2 TPA) |
| THINPRSC | prescription thin — remove cut-code-marked (KUTKOD≥2) records at cuteff — snt01 stand 3 | ✅ (`_thinprsc!`; stand 3 bit-exact; nps>1 deferred) |
| THINSDI | thin to target SDI (Zeide summation + proportional CUTEFF) | ✅ |
| THINCC | thin to residual crown cover (CCCLS, forest-grown crown width) | ✅ |
| THINHT | thin a height class (label_325 on height) | ✅ |
| THINRDEN | relative-density thin (Curtis RD) | ✅ |
| THINRDSL | relative-density SDI-line (SILVAH RD) thin | ⚪ N/A in SN (RDCLS2 gated VARACD≠NE → no-op) |
| THINMIST | mistletoe (DMR) thin | ⚪ N/A in SN (no dwarf-mistletoe model → IDMR=0 → no-op) |
| THINPT / SETPTHIN | point (plot-specific) thin (per-point + PI/NONSTK) | ✅ |
| THINQFA | Q-factor diameter-dist thin (CUTQFA + 2-record) | ✅ |

## 2. CUTS — modifiers (`cuts.f`, IACTK 201-206 + setup keywords)

| keyword | effect | status |
|---|---|---|
| SPECPREF | per-species cut preference (RDPSRT order) | ✅ |
| SPLEAVE / LEAVESP | leave named species | ⛔ |
| CUTEFF | default cutting efficiency | ⛔ |
| TCONDMLT | thin-condition multiplier | ⛔ |
| YARDLOSS | yarding-loss → scales removed merch/saw/bdft by (1−prlost) **and feeds the FFE down-wood/snag/crown fuel pools** | 🧊 **rolled into C7** (substantive effect is fuel-pool routing; standalone .sum effect nil; its `@test_broken` is the post-thin DGSCOR tail, not yardloss) |
| MINHARV | minimum-harvest threshold (skip cut below it) | ⛔ (was missing) |
| TFIXAREA | treatment fixed-area / pro-rate | ⛔ (was missing) |
| SPGROUP | species groups (vbase/initre.f:4726) referenced by a −N species field | ✅ (kw_spgroup! builds the table; the ISPCC<0 branch is wired into FIXDG/FIXHTG/FIXMORT/HTGSTOP/TOPKILL/MORTMULT/CRNMULT/TREESZCP/SPECPREF via sp_field_matches AND the species-filtering thins THINDBH/SDI/RDEN/CC/PT/QFA via _cut_eligible; bit-exact, test_spgroup.jl) |

## 3. Growth keyword multipliers / overrides (`dgdriv.f`/`htgf.f`/`regent.f`)

| keyword | effect | status |
|---|---|---|
| FIXDG | fix/scale diameter growth | ✅ (one-shot scaler, species×DBH window, scales tripled DG; bit-exact, test_fix_scalers.jl) |
| FIXHTG | fix/scale height growth | ✅ (one-shot scaler, species×DBH window, scales tripled HTG; bit-exact, test_fix_scalers.jl) |
| HTGSTOP / TOPKILL | scale height growth / top-kill (htgstp.f) | ✅ (act 110 HTG×PKIL + act 111 top-kill w/ NORMHT/ITRUNC/Behre/crown; deterministic + stochastic (htgstop_stoch) bit-exact through firing cycle, test_htgstp.jl) |
| BAIMULT | basal-area-increment multiplier (scales DDS) | ✅ (MULTS; bit-exact vs Fortran, test_multipliers.jl) |
| HTGMULT | height-growth multiplier | ✅ (MULTS; bit-exact vs Fortran) |
| CRNMULT | crown-ratio-change multiplier (sn/crown.f:319) | ✅ (active_crn_mult; scales the limited CR change over a DBH window, persistent; bit-exact within ±1 drift, test_crnmult.jl) |
| FIXCW | fix crown width | ⛔ (was missing) |
| REGDMULT / REGHMULT | regen diameter / height growth multiplier | ✅ (MULTS kinds 6/3; regent XRDGRO/XRHGRO; ±1 vs Fortran on regen cycles) |
| NOTRIPLE / NUMTRIP | tripling control (ICL4) | ✅ (NOTRIPLE→icl4=0, NUMTRIP n→icl4=n; bit-exact vs Fortran, test_tripling.jl) |

## 4. Mortality keyword overrides (`morts.f`/`fixmort.f`)

| keyword | effect | status |
|---|---|---|
| FIXMORT | keyword mortality rate override | ✅ normal path (replace/add/max/mult, DBH window, one-shot; bit-exact, test_fixmort.jl) **+ SIZE concentration** PRM(6)=10/20 (KBIG bottom-up/top-down): pools the window's kill into XMORE then re-imposes it whole-record on trees ranked by ∓(DBH+DG/bark) via the faithful `_rdpsrt!`, until XMORE is spent (morts.f:838). fixmort_big.key (pflag=10) bit-exact vs live Fortran on TPA/BA/SDI/QMD across 11 cycles. **KPOINT** multi-plot point concentration (PRM(6)=1/11/21) still deferred (single-plot stands: IPTINV=1, so it's a no-op anyway) |
| MORTMSB / MATUREW | MSB mature-stand break-up mortality (msbmrt.f) | ⚪ EFFECTIVELY INERT for self-thinning/managed stands (verified): fires only when survivors EXCEED the 85% mature self-thinning line, but BAMAX/self-thinning hold the stand AT that line, so TMORE=0 — even a 30-cycle run to QMD 38 doesn't trigger it. Rare-trigger (overmature break-up only); deterministic if ported |
| MORTMULT | mortality-rate multiplier (background only + DBH window, morts.f:518/524) | ✅ (MULTS; DBH window D1≤DBH<D2 via active_mort_mult; bit-exact on bg-mortality cycles, windowed + windowless) |
| TREESZCP | per-species size cap (SIZCAP): DG bound + size-cap mortality + HT cap | ✅ (keyword + morts size-cap floor + htgf HT cap; nomort path bit-exact, see §SIZCAP) |

## 5. Other stand management

| keyword | effect | status |
|---|---|---|
| PRUNE | pruning (option 108, act 249) | 🧊 .sum-INERT in SN (verified: no TPA/BA/TopHt/QMD change) — feeds ecopls.f pruned-log volume (C8 econ), the crown edit doesn't reach SN growth |
| FERTILIZ / FFERT | fertilizer growth response | ⛔ |
| COMPRESS | record compression to a target (comprs.f act=250) | ⛔ |
| ADDFILE | (option 22) redirect input to another unit | ⚪ NOT a tree-add — it just switches IREAD to a Fortran file UNIT (an include-file mechanism); no clean FVSjl mapping (FVSjl uses filenames, not units). 'ADDTREES' is not an SN keyword |
| MANAGED | managed-stand flag (DGF kplant term) | ✅ |
| MGMTID / RESETAGE / SETSITE | mgmt id / reset age / set site mid-run | 🟡 (MGMTID read) |

## 6. Volume / defect keywords (C5 — **.sum-affecting**, keyword-settable)

> NOTE: `DEFECT/BFDEFECT/MCDEFECT` set defect % from the KEY — so G1's defect IS
> reachable from a `.key` (not only DBS). These belong here, not just C6.

| keyword | effect | status |
|---|---|---|
| DEFECT / BFDEFECT / MCDEFECT | per-species/size cull+defect % → reduces cuft/bdft | ⛔ (G1 volume-side; .sum-affecting) |
| BFFDLN / MCFDLN | board/cubic form-class / log-length defaults | ⛔ |
| VOLUME / BFVOLUME | volume-keyword overrides | ⛔ |
| VOLEQNUM / CFVOLEQU / BFVOLEQU | per-species volume-equation selection | 🟡 (R8 Clark default works; explicit override unported) |
| FIAVBC | FIA volume/biomass calc switch | ⛔ |

## 7. Event monitor & activity scheduling (`evmon.f`/`opcycl.f`)

| keyword | effect | status |
|---|---|---|
| IF / THEN / ENDIF | conditional activity scheduling (snt01 stand 2) | ✅ event_monitor.jl (AST evaluator); stand 2 first 2 thins bit-exact; 3rd = class-boundary residual |
| COMPUTE | event-monitor variable assignment | ✅ (kw_compute! parses NAME=expr…END, evaluated each cycle in cuts! before IF conditions read them via compute_vars; ≡ direct ref, bit-exact lead stand, test_compute.jl) |
| TIMEINT | cycle length (period scaling) | 🟡 (kw_timeint! uniform path; DDS·FINT/5 + HTG·FINT/5 + mortality^FINT + year/age step; snt01 bit-exact, TIMEINT-10 TPA ≤8 / volume ≤2% — calibrated-species residual; test_timeint.jl). CYCLEAT + per-cycle GROWTH deferred |
| ESTAB-block (TALLY/PLANT/NATURAL/SPROUT) | establishment scheduling | ✅ PLANT/NATURAL; ⛔ TALLY counts / SPROUT |

## 8. Disturbance models (C7/C8 extensions — owned there)

| keyword | model | status |
|---|---|---|
| FMIN … END | FFE fire (SIMFIRE/SALVAGE/fuels/snags/CWD/carbon) | 🧊 C7 |
| MPB / DFB / DFTM / WSBW / MISTOE / BRUST | mtn pine beetle / DF beetle / DF tussock moth / W spruce budworm / mistletoe / blister rust | 🧊 (no scenario) |
| RDIN / ANIN / RRIN | root-disease model (Western root disease / Annosus) input | 🧊 (no scenario; was missing) |
| PRMFROST / CLIMATE | permafrost / climate-FVS modifiers | 🧊 |
| ECON / CHEAPO | economic analysis | 🧊 C8 (ANNUCST path exists) |

## Validation status — 3-way sweep vs live Fortran (2026-06-22)

The comprehensive 3-way sweep (162 scenarios × with/without management, vs live
Fortran) confirms the **cut logic of every ported method is correct**; the only
management-scenario residuals are *post-cut* tails, not thinning bugs:

| scenario | thin | finding |
|---|---|---|
| `s11`/`s28`/`s29_thinbta` | THINBTA | cut **bit-exact** (536→162 at the thin); post-thin ±4 TPA drift over 7 cycles = the **post-thin DGSCOR/serial-correlation tail** (cut re-ranks the stand → the stochastic DG/mortality responds slightly differently; the increment even flips sign cycle-to-cycle, so it does not propagate). Not a cut bug. |
| `s28_thindbh` | THINDBH | bit-identical (snt01 block 2). |
| `cut_thinprsc` | THINPRSC | Δ2 at the cut because the scenario uses `DESIGN …11.0` = **nps=11 plots** — the deferred multi-plot THINPRSC path (single-plot/snt01 stand 3 is bit-exact). + post-thin tail. |
| `cut_yardloss` | YARDLOSS | removes 0 merch → .sum-neutral; the ±9 TPA is the same post-THINDBH accretion tail, not yardloss. C7-coupled for the fuel pools. |
| `cut_thinsdi` | THINSDI | ✅ **ported** — bit-exact in TPA/BA every cycle (Zeide summation SDI + proportional CUTEFF); ±1 cuft tail only. |
| SPECPREF / IF-THEN | — | cut-preference + event-monitor blocks fire and cut correctly. |

**Conclusion:** ported thinning methods (THINBTA/ATA/BBA/ABA, THINDBH, single-plot
THINPRSC, SPECPREF, IF/THEN event monitor) are cut-exact. Remaining work is the
*unported* methods below + the multi-plot THINPRSC path; the post-thin numeric tail is
the same single-precision/serial-correlation floor seen in the natural-process runs.

## Triage: which ⛔ items are actually APPLIED in SN code (2026-06-23)

Grepped each keyword's effect-variable for READ references in `sn/`+`base/` (beyond
init/keyword-table). This separates real ports from set-but-not-read no-ops:

**Genuinely applied (real ports, each non-trivial):**
- `SIZCAP`/TREESZCP — ✅ DONE. The **TREESZCP** keyword (kw_treeszcp!, keyword_dispatch.jl)
  loads SIZCAP[is,1..4] immediately (no date): field 1 = species (0=all), 2 = cap DBH,
  3 = mortality rate, 4 = IDMFLG flag, 5 = HT cap (field order confirmed empirically vs
  live Fortran). The three effects: (a) DG bound (dg_bound, already present); (b) size-cap
  MORTALITY floor — ported in mortality.jl AFTER _varmrt!, before BAMAX, matching
  sn/morts.f:692: if (D+G)≥SIZCAP[is,1] & IFIX(SIZCAP[is,3])≠1 ⇒ killed=max(killed,
  P·SIZCAP[is,2]·FINT/5)≤P, where **G is OUTSIDE-bark, period-scaled (DG/BARK)·(FINT/5)**
  (the inside-bark diam_growth under-counts which trees reach the cap → too few killed);
  (c) HT cap — htgf.f:286-288 in height_growth!: if HT+HTG>SIZCAP[is,4] ⇒ HTG=max(SIZCAP[is,4]
  −HT, 0.1) (the 0.1 floor: trees already past the cap crawl, never shrink). Validated by
  test/integration/test_treeszcp.jl (3 scenarios vs Fortran, 106 asserts). Residuals: cap
  mid-cycle TPA/BA carry the regen response to cap-driven mortality (QMD bit-exact, endpoint
  matches); htcap TopHt drifts ≤4' as a declining-stand artifact (TPA/BA/QMD bit-exact).
- `FIXMORT` — ✅ normal path DONE (morts.f:1017). apply_fixmort! (keyword_dispatch.jl) overrides
  killed[] AFTER the BA-check (the last word on the kill), one-shot in the date's cycle, over a
  species×DBH window: IP 1 replace (P·rate), 2 add, 3 max, 4 multiply (kill·rate≤P), selected by
  PRM(5) (0/1/2/3), with Fortran's rate clamps. Needed a companion fix to mortality! ordering:
  **TPAMRT (the self-thinning line-reset, morts.f:772) is locked from the BA-check survivors
  BEFORE FIXMORT**, so the forced kill doesn't move next cycle's self-thinning line — without it
  the recovery ran TPA up to ~6% high. Bit-exact every cycle on 3 scenarios (replace, multiply,
  big-tree replace) vs live Fortran (test_fixmort.jl). ✅ **SIZE concentration DONE**
  (PRM(6)=10/20, morts.f:838-935): the window's mortality is pooled into XMORE per the IP, the
  in-window kills are zeroed where the IP replaces, then XMORE is re-imposed WHOLE-RECORD on the
  trees ranked by WORK3=∓(DBH+DG/bark) — KBIG=1 negates so the smallest grown trees go first
  (bottom up), KBIG=2 the largest (top down) — via the faithful `_rdpsrt!` (descending, FVS
  tie-break), killing each record fully (CREDIT+=tpa−killed) until CREDIT reaches XMORE (last
  record partial: killed+=XMORE−CREDIT). fixmort_big.key (pflag=10) is bit-exact vs live Fortran
  on TPA/BA/SDI/TopHt/QMD across all 11 cycles; the only diffs are ±1 in the volume cols (9-12)
  and per-acre growth (24-25) — the documented Float32 DGSCOR/volume noise, not the kill itself.
  DEFERRED: **KPOINT** multi-plot point concentration (PRM(6)=1/11/21, morts.f:937-1015) — a no-op
  on single-plot stands (IPTINV=1) anyway. Species groups (ISPCC<0) ✅ via SPGROUP + sp_field_matches.
- `FIXCW` — cwidth.f (crown-width override). ⚪ **OUTPUT-ONLY for the .sum** (verified): CRWDTH
  is referenced only by the calculator (cwidth.f), record bookkeeping that carries it along
  (comprs/tremov/triple), and OUTPUT consumers (sstage structure-class, svsnad SVS, evldx
  event-monitor var). It never feeds DGF/HTGF/MORTS/DENSE — those use crown RATIO, not width.
  So a FIXCW port changes no .sum growth number; defer until SVS/structure output is in scope.
- `HTGSTP` (HTGSTOP/TOPKILL) — ✅ DONE. htgstp! (keyword_dispatch.jl), called in grow_cycle!
  after TRIPLE/MORTS and before UPDATE (gradd.f:158). act 110 (HTGSTOP) scales HTG by PKIL;
  act 111 (TOPKILL) sets HT=H·(1−PKIL≤0.8), and for tall trees (H≥25, D≥6) whose Behre top
  diameter ≥4 marks a permanent broken top (NORMHT/ITRUNC) and cuts the crown ratio (ICR=−NEW).
  PKIL=BACHLO(AVEPRB,STDPBR), deterministic (=AVEPRB, no RNG) when STDPBR≤0; RANN escape when
  PRB<1; records walked in species-sorted IND1 order for RNG-exactness when stochastic. Needed a
  companion fix to crown_ratio_update! — the **negative-ICR bypass** (sn/crown.f:271): a crown
  already adjusted by topkill/pest models (ICR<0) is restored to +ICR and NOT recomputed that
  cycle; without it the top-killed trees' crown (hence DG/mortality) drifted and TPA ran ~10% high.
  Deterministic scenarios (HTGSTOP 0.5×, TOPKILL 0.5× >30') bit-exact every cycle vs live Fortran
  (test_htgstp.jl) AND a stochastic bare-plant scenario (STDPBR=0.2, htgstop_stoch) bit-exact through the firing cycle — the species-sorted IND1 RNG traversal is confirmed. IMC
  (management code) and ABIRTH (age) are set in Fortran but don't affect the .sum, so skipped.
- FIXDG/FIXHTG — ✅ DONE. grincr.f:451-525: DG/HTG·PRM(2) over a species×DBH window, applied
  in `apply_fix_scalers!` (keyword_dispatch.jl) after all growth / before MORTS. TWO things the
  earlier buggy attempt missed: (1) it is **ONE-SHOT** (OPDONE) — fires only in the cycle whose
  [start, start+period) range holds the keyword date, not every cycle (confirmed empirically:
  0.3× at 1995 drops QMD only in the 1995-cycle, then the gap persists ~constant, not runaway);
  (2) it must scale the **tripled** DG/HTG too — the stash dgU/dgL (htgU/htgL), matching FVS's
  DG(ITFN)/DG(ITFN+1). Reuses the GrowthMultiplier d1/d2 window. Bit-exact every cycle on 3
  scenarios (all/windowed DG, HTG) vs live Fortran (test_fix_scalers.jl). Species groups (ISPCC<0)
  ✅ via SPGROUP + sp_field_matches (test_spgroup.jl).

**Set-but-not-read in SN (0 application refs ⇒ likely NO-OP in SN, or external component):**
- CUTEFF, MINHARV, TCONDMLT — 0 refs.
- ⚠ CORRECTION: CRNMULT and TOPKILL were on this list but are NOT no-ops — both are applied in
  SN (CRNMULT scales the crown-ratio change in sn/crown.f:319; TOPKILL is the htgstp.f act 111
  top-kill). Both are now ✅ ported. The lesson: "0 application refs" from a coarse grep missed
  them because the grep keyed on the keyword name, not the COMMON variable (CRNMLT / IACT 111).
- SPLEAVE/LEAVESP — only `grinit.f:125 LEAVESP(I)=.FALSE.` (init), never checked in the cut logic.
- DEFECT/BFDEFECT/MCDEFECT — CFDEFT/BFDEFT set in sdefet.f/volkey.f, never read in sn/base
  (the per-tree defect reduction is in the NVEL volume LIBRARY, a separate component). Also
  the DEFECT keyword CRASHES this Fortran build on a simple scenario (exit 2). ⇒ verify it is
  even active in SN before porting; FVSjl's R8-Clark volume would need an NVEL-style defect hook.
- ⚠ CAVEAT: "0 refs" used a guessed effect-variable name; some may apply under a different
  COMMON name. Confirm empirically (does the keyword change the .sum?) before declaring no-op.

⇒ The cheap management wins (the 5 MULTS multipliers) are DONE. Every remaining item is a
focused chunk, not a quick port. Several listed-⛔ items are probably SN no-ops.

## Remaining work (current — 2026-06-23, after the growth/mortality surface + SPGROUP)

The entire growth/mortality keyword surface is ported & bit-exact: all CUTS methods (§1),
all MULTS multipliers, FIXDG/FIXHTG, FIXMORT, TREESZCP, HTGSTOP/TOPKILL, CRNMULT, SPGROUP
species groups. The DGSCOR cubic-volume drift is resolved as irreducible Float32 noise
(see DIVERGENCES.md §1). What is left, by kind of work:

**A. Genuinely-applied, .sum-affecting — real ports (ranked most-upstream → downstream):**
- ~~**Cycle calendar** — `TIMEINT`~~ 🟡 RESOLVED & IMPLEMENTED (uniform path). The hidden
  period-scaling mechanism was traced via a Fortran DEBUG dump: the DDS is FINT-independent, but
  dgdriv.f:325/715 applies SCALE=FINT/YR ⇒ DDS·(FINT/5) (YR=5 SN base); a 10-yr cycle ⇒ ×2 DG.
  Ported: diameter_growth! scales DDS by sfint/5; height by FINT/5; mortality^FINT; year/age step.
  Companion fixes: the morts grown-DBH G must NOT re-apply FINT/5 (diam_growth is already cycle-
  scaled); the calibration VMLT uses YR=5 not FINT; the .sum last-row year uses the real period.
  snt01 (period 5) bit-exact; TIMEINT-10 TPA ≤8, volume ≤2% (calibrated-species residual, like the
  DGSCOR tail). REMAINING: full YR-vs-FINT calibration split for bit-exactness; per-cycle GROWTH/
  CYCLEAT boundaries.
- ~~**Tripling control** — `NOTRIPLE` / `NUMTRIP`~~ ✅ DONE — wired through s.control.icl4
  (default 2; NOTRIPLE→0, NUMTRIP n→n). Was a real gap: FVSjl ignored NOTRIPLE (20 cols diverged)
  AND the bare-PLANT scenarios were silently passing for the wrong reason (FVSjl wrongly tripled
  them and the 3× averaging masked the regen DGSCOR tail at ±1). The fix exposed the true ±2 no-trip
  tail; the base+NOTRIPLE stand is bit-exact (test_tripling.jl).
- ~~ADDFILE~~ ⚪ unit-redirect include, not a tree-add (verified); no clean FVSjl mapping.
- ~~**COMPUTE**~~ ✅ DONE — kw_compute! stores NAME=expression defs (parsed with the event-monitor
  expression parser); cuts! evaluates them each cycle before the IF conditions, and _event_var resolves
  them; bit-exact (COMPUTE MYCYC=CYCLE ≡ direct CYCLE; firing thins match Fortran), test_compute.jl.
- **CYCLEAT / TIMEINT** scheduling boundaries (if not covered by the calendar item).
- ~~PRUNE~~ 🧊 .sum-inert in SN (verified) — C8 econ pruned-log volume only.
- **Volume overrides** — `VOLEQNUM`/`CFVOLEQU`/`BFVOLEQU` (🟡), `VOLUME`/`BFVOLUME`,
  `BFFDLN`/`MCFDLN`, `FIAVBC` (C5 volume side).
- **ESTAB TALLY / SPROUT** — tally-count regen + stump sprouting (downstream, C4 regen).
- `COMPRESS` — only fires with the keyword or >~3000 records; comprs.f is a 762-line
  subsystem (incl. EIGEN!) — NOT least-dependent.

**B. Deferred sub-paths inside already-✅ items:**
- FIXMORT point/size concentration reallocation (PRM(6)≥10 — KBIG/KPOINT, morts.f:838-1015).
- THINPRSC multi-plot (nps>1) path.
- ~~Species-group thin-method filter~~ ✅ DONE — `_cut_eligible` (and `_clsstk`/`_sdi_zeide`/
  `_rd_curtis`) take the SPGROUP table; the species-filtering thins (THINDBH/THINSDI/THINRDEN/
  THINCC/THINPT/THINQFA) thread `s.control.sp_groups`, so `THINDBH −1` resolves a group (bit-exact
  vs Fortran, test_spgroup.jl). REMAINING: group-**name** field refs (only −N numeric now).
- IF/THEN snt01 stand-2 3rd-thin class-boundary residual.

**C. Listed-⛔ but likely SN NO-OPS — confirm empirically (does it change the .sum?) first:**
- SPLEAVE/LEAVESP (only grinit init, never read in cut logic), CUTEFF / MINHARV / TCONDMLT
  (0 refs), MORTMSB (QMDMSB=999 inert), FIXCW (verified output-only),
  DEFECT/BFDEFECT/MCDEFECT (external NVEL lib; the keyword crashes this build — verify active).
  ⚠ "0 refs" mis-flagged CRNMULT+TOPKILL once; always confirm with a .sum diff.

**D. Out of scope here (C7/C8):** YARDLOSS (C7 fuel pools), FFE fire, insects, root disease
(RDIN/ANIN/RRIN), PRMFROST/CLIMATE, ECON/CHEAPO.
