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
| SPGROUP (via SPCODES/SPLABEL) | species groups referenced by SPECPREF/LEAVESP | ⛔ |

## 3. Growth keyword multipliers / overrides (`dgdriv.f`/`htgf.f`/`regent.f`)

| keyword | effect | status |
|---|---|---|
| FIXDG | fix/scale diameter growth | ⛔ |
| FIXHTG | fix/scale height growth | ⛔ |
| HTGSTOP / TOPKILL | stop height growth / topkill edits | ⛔ |
| BAIMULT | basal-area-increment multiplier (scales DDS) | ✅ (MULTS; bit-exact vs Fortran, test_multipliers.jl) |
| HTGMULT | height-growth multiplier | ✅ (MULTS; bit-exact vs Fortran) |
| CRNMULT | crown-ratio/width multiplier | ⛔ (was missing) |
| FIXCW | fix crown width | ⛔ (was missing) |
| REGDMULT / REGHMULT | regen diameter / height growth multiplier | ✅ (MULTS kinds 6/3; regent XRDGRO/XRHGRO; ±1 vs Fortran on regen cycles) |
| NOTRIPLE / NUMTRIP | tripling control (LTRIP) | 🟡 (hardcoded limit=2) |

## 4. Mortality keyword overrides (`morts.f`/`fixmort.f`)

| keyword | effect | status |
|---|---|---|
| FIXMORT | keyword mortality rate override | ⛔ |
| MORTMSB | MSB alternate mortality (QMDMSB) | ⛔ (default 999 no-op) |
| MORTMULT | mortality-rate multiplier (background only + DBH window, morts.f:518/524) | ✅ (MULTS; DBH window D1≤DBH<D2 via active_mort_mult; bit-exact on bg-mortality cycles, windowed + windowless) |
| TREESZCP | per-species size cap (SIZCAP): DG bound + size-cap mortality + HT cap | ✅ (keyword + morts size-cap floor + htgf HT cap; nomort path bit-exact, see §SIZCAP) |

## 5. Other stand management

| keyword | effect | status |
|---|---|---|
| PRUNE | pruning (crown/CR edit + pruned-log volume) | ⛔ (was missing) |
| FERTILIZ / FFERT | fertilizer growth response | ⛔ |
| COMPRESS | record compression to a target (comprs.f act=250) | ⛔ |
| ADDFILE / ADDTREES | add tree records mid-run | ⛔ |
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
| COMPUTE | event-monitor variable assignment | ⛔ |
| CYCLEAT / TIMEINT | explicit cycle boundaries / interval | ⛔ |
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
- `FIXMORT` — morts.f subsystem (per-tree rate override + point/size reallocation).
- `FIXCW` — cwidth.f (crown-width override).
- `HTGSTP` (HTGSTOP/TOPKILL) — htgstp.f 200-ln topkill subsystem (reduces HT + volume).
- FIXDG/FIXHTG — grincr.f:481 (DG/HTG·PRM(2) over a DBH window). ⚠ ATTEMPTED a post-growth-
  pass port this session → BUGGY (FIXDG over-scaled QMD 23.8 vs 16.9; FIXHTG corrupted
  TopHt down 58→44 vs ~no-effect). REVERTED. Needs the post-pass + tripling-stash + growth/
  mortality-compensation interaction debugged (isolate on a NOTRIPLE stand first).

**Set-but-not-read in SN (0 application refs ⇒ likely NO-OP in SN, or external component):**
- CRNMULT, TOPKILL, CUTEFF, MINHARV, TCONDMLT — 0 refs.
- SPLEAVE/LEAVESP — only `grinit.f:125 LEAVESP(I)=.FALSE.` (init), never checked in the cut logic.
- DEFECT/BFDEFECT/MCDEFECT — CFDEFT/BFDEFT set in sdefet.f/volkey.f, never read in sn/base
  (the per-tree defect reduction is in the NVEL volume LIBRARY, a separate component). Also
  the DEFECT keyword CRASHES this Fortran build on a simple scenario (exit 2). ⇒ verify it is
  even active in SN before porting; FVSjl's R8-Clark volume would need an NVEL-style defect hook.
- ⚠ CAVEAT: "0 refs" used a guessed effect-variable name; some may apply under a different
  COMMON name. Confirm empirically (does the keyword change the .sum?) before declaring no-op.

⇒ The cheap management wins (the 5 MULTS multipliers) are DONE. Every remaining item is a
focused chunk, not a quick port. Several listed-⛔ items are probably SN no-ops.

## Suggested order (when starting)
1. THINCC/HT/RDEN/RDSL/PT/QFA (remaining label_400 SDI-class thins; THINSDI ✅ done) + LEAVESP/MINHARV
2. Multi-plot (nps>1) THINPRSC path
3. DEFECT/BFDEFECT/MCDEFECT (.sum-affecting volume) + a defect scenario
4. Growth/mort multipliers (FIXDG/HTGMULT/MORTMULT/… — all keyword no-ops at default)
5. PRUNE / FERTILIZ / COMPRESS / ADDFILE
6. FFE fire (C7), insects (no scenario), econ (C8)
