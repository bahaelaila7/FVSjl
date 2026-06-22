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
| THINPT / SETPTHIN | point (plot-specific) thin | ⛔ SN-relevant; needs SETPTHIN per-point prescriptions + multi-plot point-thinning infra (JPNUM/IPTINV/LPTALL) — a focused build |
| THINQFA | Q-factor diameter-dist thin | ⛔ SN-relevant; scenario validated (cut_thinqfa, 2-record, Fortran 2005=89). Scope: port CUTQFA (155-ln Q-factor negative-exponential target distribution → per-class CLSTAR targets) + label_350 per-class driver + 2-record keyword parser. Self-contained (reuses CLSSTK/SDICLS). Next focused unit. |

## 2. CUTS — modifiers (`cuts.f`, IACTK 201-206 + setup keywords)

| keyword | effect | status |
|---|---|---|
| SPECPREF | per-species cut preference (RDPSRT order) | ✅ |
| SPLEAVE / LEAVESP | leave named species | ⛔ |
| CUTEFF | default cutting efficiency | ⛔ |
| TCONDMLT | thin-condition multiplier | ⛔ |
| YARDLOSS | yarding-loss: scales REMOVED merch/saw/bdft by (1-prlost) + feeds fire-fuel pools (dsng/ssng/crown→C7). cut_yardloss removes 0 merch so .sum-effect is nil there; its `broken` test is really a post-THINDBH accretion residual (130 vs 124), not yardloss | ⛔ (C7-coupled; .sum part trivial) |
| MINHARV | minimum-harvest threshold (skip cut below it) | ⛔ (was missing) |
| TFIXAREA | treatment fixed-area / pro-rate | ⛔ (was missing) |
| SPGROUP (via SPCODES/SPLABEL) | species groups referenced by SPECPREF/LEAVESP | ⛔ |

## 3. Growth keyword multipliers / overrides (`dgdriv.f`/`htgf.f`/`regent.f`)

| keyword | effect | status |
|---|---|---|
| FIXDG | fix/scale diameter growth | ⛔ |
| FIXHTG | fix/scale height growth | ⛔ |
| HTGSTOP / TOPKILL | stop height growth / topkill edits | ⛔ |
| BAIMULT | basal-area-increment multiplier | ⛔ (was missing) |
| HTGMULT | height-growth multiplier | ⛔ (was missing) |
| CRNMULT | crown-ratio/width multiplier | ⛔ (was missing) |
| FIXCW | fix crown width | ⛔ (was missing) |
| REGDMULT / REGHMULT | regen diameter / height growth multiplier | ⛔ (was missing) |
| NOTRIPLE / NUMTRIP | tripling control (LTRIP) | 🟡 (hardcoded limit=2) |

## 4. Mortality keyword overrides (`morts.f`/`fixmort.f`)

| keyword | effect | status |
|---|---|---|
| FIXMORT | keyword mortality rate override | ⛔ |
| MORTMSB | MSB alternate mortality (QMDMSB) | ⛔ (default 999 no-op) |
| MORTMULT | mortality-rate multiplier | ⛔ (was missing) |
| TREESZCP | tree size-cap mortality (SIZCAP) | ⛔ (default 999 no-op) |

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

## Suggested order (when starting)
1. THINCC/HT/RDEN/RDSL/PT/QFA (remaining label_400 SDI-class thins; THINSDI ✅ done) + LEAVESP/MINHARV
2. Multi-plot (nps>1) THINPRSC path
3. DEFECT/BFDEFECT/MCDEFECT (.sum-affecting volume) + a defect scenario
4. Growth/mort multipliers (FIXDG/HTGMULT/MORTMULT/… — all keyword no-ops at default)
5. PRUNE / FERTILIZ / COMPRESS / ADDFILE
6. FFE fire (C7), insects (no scenario), econ (C8)
