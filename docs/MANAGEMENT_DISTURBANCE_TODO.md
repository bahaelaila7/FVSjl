# Management & disturbance completeness tracker

Items in **C3 / C4 / C5** that are management activities, disturbances, or
keyword-gated options (NOT self-contained natural processes). Driven from the FVS
**decision flow** as the completeness oracle.

> **DO NOT START these until [`NATURAL_PROCESS_TODO.md`](NATURAL_PROCESS_TODO.md) is
> complete** (user directive). Tracked here only so nothing is forgotten.
> Same discipline applies: port the real FVS code; tests only catch/fix bugs.

Legend: ⛔ unported · 🟡 partial · ✅ done

## C3 — growth keyword options

| keyword/item | FVS source | effect | status |
|---|---|---|---|
| FIXDG / FIXHTG | `dgdriv.f`/`htgf.f` | per-tree/stand DG·HTG multipliers | ⛔ (no-op at default) |
| HTGSTP (HTGSTOP/TOPKILL) | `htgstp.f` | keyword height-growth stop / topkill edits | ⛔ |
| FFERT | `ffert.f` | fertilizer growth response | ⛔ |
| NOTRIPLE / TRIPLE / ICL4 | `rdmn2.f` | LTRIP tripling control via keyword | 🟡 (hardcoded limit=2) |

## C4 — CUTS methods (9 unported of 14)

| keyword | ICFLAG | semantics | scenario gap | status |
|---|---|---|---|---|
| THINAUTO | 1 | auto-thin to FULSTK on stocking trigger | — | ⛔ |
| THINPRSC | 7 | prescription (per-DBH-class residual table) | snt01 stand 3 | ⛔ |
| xSALVAGE | 9 | salvage dead/damaged | snt01 stand 4 | ⛔ |
| THINSDI | 10 | thin to target SDI | — | ⛔ |
| THINHT | 12 | thin a height class | — | ⛔ |
| THINMIST | 13 | mistletoe (DMR) thin | — | ⛔ |
| THINRDEN | 14 | relative-density thin | — | ⛔ |
| THINPT / SETPTHIN | 15 | point (plot-specific) thin | — | ⛔ |
| THINQFA | 17 | Q-factor diameter-dist thin | — | ⛔ |

(Done: THINBTA/THINATA/THINBBA/THINABA/THINDBH ✅)

## C4 — CUTS modifiers (5 unported of 6)

| keyword | IACTK | effect | status |
|---|---|---|---|
| TCONDMLT | 202 | thin-condition multiplier | ⛔ |
| YARDLOSS | 203 | yarding-loss removed-volume accounting (snt01 stand 4) | ⛔ |
| SPLEAVE / LEAVESP | 206 | leave named species | ⛔ |
| SPGROUP | 125 | define species groups (for SPECPREF/LEAVESP) | ⛔ |
| CUTEFF | 52 | default cutting efficiency | ⛔ |

(Done: SPECPREF ✅)

## C4 — event monitor & scheduling

| item | FVS source | effect | status |
|---|---|---|---|
| EVMON + IF/THEN/ENDIF | `evmon.f`,`algmon.f`,`evtact.f`,`evtstv.f` (~2.4K LOC) | algebraic conditions per cycle gate activities (snt01 stand 2) | ⛔ |
| OPCYCL / OPGET / OPDEL1 | `opcycl.f`/`opnew.f` | full activity scheduler (partial today) | 🟡 |
| COMPRS / COMCUP(keyword) | `comprs.f` (act=250) | COMPRESS-keyword record compression to a target | ⛔ (the zero-PROB-delete half is natural → other tracker) |

## C4 — keyword mortality

| keyword | FVS source | effect | status |
|---|---|---|---|
| MSB alternate mortality | `morts.f` (QMDMSB) | replaces density mortality when QMD>QMDMSB | ⛔ (default 999 ⇒ no-op) |
| SIZE-CAP mortality | `morts.f` (SIZCAP) | kill when d+g ≥ SIZCAP | ⛔ (default 999) |
| FIXMORT | `fixmort.f` | keyword mortality multiplier | ⛔ |

## C4/C7 — disturbance models

| model | FVS source | status |
|---|---|---|
| FFE fire (SIMFIRE/SALVAGE/fuels/snags/CWD/carbon) | `extensions/fire/` (~14K LOC, 69 files) | ⛔ **C7** (snt01 stand 4) |
| insect/disease (MPB/DFB/TM/BWE/MISTOE) | `MPB*/DFB*/TM*/BWE*/MISTOE` | ⛔ (no scenario; 🧊 in plan) |

## C5/C8 — economics & extension output

| item | FVS source | status |
|---|---|---|
| ECON (ECSETP/eccalc beyond ANNUCST) | `econ.f` | 🟡 C8 (ANNUCST-only path exists) |
| extension reports (ESOUT/CVOUT/MPBOUT/…) | various | ⛔/🧊 output |

## Suggested order (when natural is done)
1. THINPRSC + SALVAGE + YARDLOSS (snt01 stands 3-4 cut methods/modifiers)
2. Event monitor (EVMON + IF/THEN) — stand 2; its own chunk
3. Remaining CUTS methods/modifiers (THINSDI/THINHT/… + LEAVESP/SPGROUP/…)
4. Keyword mortality (MSB/SIZECAP/FIXMORT) + growth options (FIXDG/HTGSTP/FFERT)
5. FFE fire (C7), econ (C8), pests
