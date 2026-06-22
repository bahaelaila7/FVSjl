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

Legend: ✅ done · 🟡 partial · ⛔ unported · 🧊 C7/C8 extension

## 1. CUTS — thinning / harvest methods (`cuts.f`, ICFLAG)

| keyword | semantics | status |
|---|---|---|
| THINBTA/THINATA/THINBBA/THINABA | from below/above to residual TPA/BA | ✅ |
| THINDBH | proportional DBH-class to residual TPA/BA | ✅ |
| THINAUTO | auto-thin to FULSTK on stocking trigger | ⛔ |
| THINPRSC | prescription thin — remove cut-code-marked (KUTKOD≥2) records at cuteff — snt01 stand 3 | ✅ (`_thinprsc!`; stand 3 bit-exact; nps>1 deferred) |
| THINSDI | thin to target SDI | ⛔ |
| THINCC | thin to residual **crown competition / cover** | ⛔ (was missing) |
| THINHT | thin a height class | ⛔ |
| THINRDEN | relative-density thin | ⛔ |
| THINRDSL | relative-density **SDI-line / Reineke** thin | ⛔ (was missing) |
| THINMIST | mistletoe (DMR) thin | ⛔ |
| THINPT / SETPTHIN | point (plot-specific) thin | ⛔ |
| THINQFA | Q-factor diameter-dist thin | ⛔ |

## 2. CUTS — modifiers (`cuts.f`, IACTK 201-206 + setup keywords)

| keyword | effect | status |
|---|---|---|
| SPECPREF | per-species cut preference (RDPSRT order) | ✅ |
| SPLEAVE / LEAVESP | leave named species | ⛔ |
| CUTEFF | default cutting efficiency | ⛔ |
| TCONDMLT | thin-condition multiplier | ⛔ |
| YARDLOSS | yarding-loss removed-volume accounting (snt01 stand 4) | ⛔ |
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
| IF / THEN / ENDIF | conditional activity scheduling (snt01 stand 2) | ⛔ |
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

## Suggested order (when starting)
1. THINPRSC + (SALVAGE) + YARDLOSS — snt01 stands 3-4
2. Event monitor (IF/THEN/COMPUTE) — snt01 stand 2
3. DEFECT/BFDEFECT/MCDEFECT (.sum-affecting volume) + a defect scenario
4. Remaining CUTS methods/modifiers (THINSDI/CC/HT/RDEN/RDSL/PT/QFA + LEAVESP/MINHARV/…)
5. Growth/mort multipliers (FIXDG/HTGMULT/MORTMULT/… — all keyword no-ops at default)
6. PRUNE / FERTILIZ / COMPRESS / ADDFILE
7. FFE fire (C7), insects (no scenario), econ (C8)
