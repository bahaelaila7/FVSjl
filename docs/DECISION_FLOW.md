# FVS-Southern decision flow (input → output)

> Branch-level companion: [`DECISION_FLOW_DETAILED.md`](DECISION_FLOW_DETAILED.md)
> expands each hot-path routine below into its individual decision branches with
> per-branch port status.
>
> Interactive views (open in a browser; click a node to break it into the functions
> it calls down to atomic math/intrinsics — hover/click shows the **code excerpt**;
> fill = port status, shape = scope ▭ generic / ⬡ Southern-variant / ◆ extension):
> - [`decision_flow.html`](decision_flow.html) — the **FVS** (Fortran oracle) graph,
>   588 routines. The source-of-truth semantics.
> - [`decision_flow_fvsjl.html`](decision_flow_fvsjl.html) — the **FVSjl** graph,
>   126 routines (idiomatic, compressed). Put side by side with the FVS graph to see
>   where the rewrite folded routines together (e.g. `grow_cycle!` = GRINCR+GRADD).
>
> Regenerate either with `tools/gen_callgraph.js <srcRoot> <outHtml> <rootsCSV>`.

A bird's-eye map of how a simulation flows from input to output in the original
Fortran (`FVSsn`), with the corresponding **FVSjl** entry point and **port
status** for each major step. It is intentionally one level deep: it names the
major code snippets (subroutines) and the conditions that gate them, but does not
expand their internals. Use it to (a) onboard, (b) see at a glance what is and is
not ported, and (c) design test scenarios that exercise a specific branch.

Legend: ✅ ported & validated · 🟡 partial · ⛔ not ported (in-scope, planned) ·
🧊 out of scope by the plan (no test scenario exercises it) · 🔁 compressed
(several Fortran routines fold into one FVSjl function).

Oracle source roots: `/workspace/FVSjulia/src` (1:1 transliteration, "Oracle A").
FVSjl roots: `/workspace/FVSjl/src`.

---

## 1. Top level — `fvs.f` main driver

```
INITRE        read keywords (.key) + tree data (.tre); set all options   ✅ initialize()
  ├ keyword dispatch (STDINFO/SITECODE/DESIGN/INVYEAR/NUMCYCLE/thinning…) ✅ engine/keyword_dispatch.jl
  └ HABTYP     decode ecological-unit field → PCOM (table index/alpha)    ✅ variants/southern/habitat.jl
MPBOPS/TMOPS/DFBSCH   insect/disease activity schedules                   🧊 (no pest scenario)
ECSETP        econ setup                                                  🧊 C8 (econ.jl is ANNUCST-only)
OPCYCL/OPLIST schedule activities onto cycles                             🟡 event-monitor wired; full OP* partial
SETUP         per-stand setup (site, RNG seed, density consts)            ✅ setup_growth! / init
NOTRE         build the live tree records (expand, defaults)              ✅ notre!
MPSDLP/DFBINV pest inventory                                              🧊

── PER-CYCLE LOOP (icyc = 1 … NCYC) ───────────────────────────────────
 A. REPORT current stand (before growing it)
    CRATET    dub missing heights, resolve broken-top NORMHT             ✅ dub_missing_heights! (at setup)
    ESFLTR    establishment filter                                        ⛔ C4 regen
    CWIDTH    crown width                                                 ✅ crown width (CSV-driven)
    VOLS      per-tree volume (cuft/merch/sawtimber/bdft)                 ✅ compute_volumes! (R8 Clark etc.)
    PCTILE/DIST/COMP  percentile + distribution + composition stats       🟡 stand_pct!/stats partial
    STATS / DISPLY / SUMOUT  stand tables + the .sum row                  ✅ io/summary.jl (summary_row/write_sum_file)
    DBSREFERENCE / DISPLY    SQLite + .out tables                         🟡 C6 (sndb tables ported separately)
    EXTREE / MISPRT / RDPR / BRPR / SVSTART   tree list, mistletoe,       🧊 / ⛔ (reports & SVS)
              down-wood, snag, SVS reports
    CVGO / CVBROW / CVCNOP    canopy cover                                🧊
 B. GROW the stand one cycle
    TREGRO    = GRINCR (compute increments) + GRADD (apply + update)      ✅ 🔁 grow_cycle!   (see §2)
 C. post-growth reports (EXTREE/DISPLY/RESAGE/…), extension outputs
    (ESOUT/CVOUT/MPBOUT/DFBOUT/TMOUT/BWEOUT/BRROUT), GENPRT               🟡/🧊
```

**FVSjl shape.** The host loop lives in `io/summary.jl:write_sum_file` (and the
test harnesses): per cycle it calls `summary_row` (the REPORT step) then
`grow_cycle!` (the GROW step), exactly mirroring A→B above. Initialization is
`initialize()` + `notre!` + `setup_growth!` + `compute_volumes!`.

---

## 2. The growth cycle — `tregro.f` → `GRINCR` + `GRADD`

This is the heart of C3/C4. In Fortran it is split across two routines; FVSjl
folds the *exercised* parts into the single `grow_cycle!` (🔁).

### 2a. `GRINCR` — compute increments (no records applied yet)

```
RDMN2 / RDTRP        management read; set tripling flag LTRIP            🟡 FVSjl hardcodes TRIPLE_CYCLE_LIMIT=2
FMSDIT               FFE stand SDI for fire                              🧊 C7
SILFTY               silviculture forest type                           🟡 compute_forest_type!
SDICAL → BTSDIX      stand max SDI                                       ✅ stand_sdimax / sdical path
SDICLS → SDIBC       SDI class                                          🟡
SSTAGE               stand structure stage                              ⛔ (stand-class only; not in .sum math)
EVMON(1)             event monitor (pre-thin)                           🟡
save OLD* density (OLDTPA/OLDBA/OLDAVH/ORMSQD/…)                          ✅ implicit in compute_density!
CUTS                 scheduled thinning / harvest                       🟡 THINDBH ported (cuts.jl); other methods pending
CVGO + DENSE(if cut) recompute density post-thin                        ✅ compute_density!
ATSDIX = SDICAL; SDICLS → SDIAC   after-treatment SDI                    🟡
save per-tree vol history (PTOCFV/PDBH/PHT/NCFDEF/NBFDEF)                 ✅ old_cfv/old_tpa snapshots
COMCUP               record compression                                 ⛔ (COMPRS — not the snt blocker)
DFTMGO/MPBGO/DFBGO/BWEGO   insect/disease "go" flags + TMBMAS            🧊
─ growth core ─
DGDRIV               diameter growth (DGF + serial-corr / tripling)     ✅ diameter_growth!
HTGF                 height growth (HTGF/HTCALC/HTDBH)                   ✅ height_growth!
REGENT(false,1)      small-tree (<3") height-driven growth              ✅ small_tree_growth!  (establishment mode ⛔)
FIXDG / FIXHTG       keyword DG/HTG multipliers                         ⛔ (keyword option)
MORTS                periodic mortality                                 ✅ mortality!  (see §3 for branch map)
TRIPLE + REASS       record tripling (cyc1-2), ITRN*=3                  ✅ triple_records!
FFERT                fertilizer response                                ⛔ (keyword option)
```

### 2b. `GRADD` — apply growth & update records

```
MPBCUP/DFBWIN/MISTOE/TMCOUP/BWECUP   insect/disease record updates       🧊
FMMAIN               FFE fire effects + fire-caused mortality            🧊 C7  (this is the s10_fire/fire_* gap)
BRTREG / RDTREG      sprout / planted regeneration                       ⛔ C4 regen
HTGSTP               HTGSTOP / TOPKILL keyword height edits              ⛔ (keyword option; no-op otherwise)
UPDATE               apply DG/HTG → DBH/HT (+ bark, NORMHT for topkill)  ✅ inline in grow_cycle!
RDPSRT / DENSE       sort + recompute density                           ✅ compute_density!
CVGO                 canopy cover                                        🧊
CLAUESTB / ESNUTR    establishment (natural + nutrient hook)            ⛔ C4 regen
DENSE                density (post-establishment)                        ✅
CROWN                crown-ratio update (Weibull, ±change cap)           ✅ crown_ratio_update!
CWIDTH               crown width                                         ✅
PCTILE/DIST/COMP     end-of-cycle percentile/distribution stats         🟡
```

**Order check (FVSjl `grow_cycle!`):** `compute_density!` → `diameter_growth!`
→ `height_growth!` → `small_tree_growth!` → `mortality!` → mort accounting →
`triple_records!` → apply DG/HTG (UPDATE) → `compute_volumes!` (VOLS) → accr
accounting → `crown_ratio_update!` (CROWN) → `cycle++`. This matches the Fortran
order DGDRIV→HTGF→REGENT→MORTS→TRIPLE | UPDATE→DENSE→CROWN→VOLS, with VOLS pulled
inside the cycle (it is the next REPORT step's input).

---

## 3. Worked drill-down — `MORTS` (why module-presence ≠ completeness)

The decision flow above is at *module* granularity: it tells you `MORTS` is
ported. But the bugs found by numeric chasing this project (HABTYP index decode,
BAMAX cap, size-cap) were **missing branches inside a ported module**. A module
map cannot surface those; each ported routine needs its own branch checklist.
`MORTS` (`morts.f`) is the template:

```
RMSQD==0 → reset CEPMRT/SLPMRT                                            ✅
PMSDIL/PMSDIU rescale (>1 ⇒ /100)                                         ✅
compute DQ0/DR0 (Zeide) start QMD, DQ10/DR10 end QMD                      ✅
ICYC>1 & |t-TPAMRT|>1 → reset self-thinning line                          🟡 persisted-once (ok for closed stands)
SDICAL → SDIMAX; const = SDIMAX/0.02483133                               ✅
SDIMAX<5 → background only                                               ✅
t > t85d0  → tn10 = t85d10 (over-dense)                                  ✅
t55d0<t≤t85d0 → solve self-thinning line (iterate treeit)                ✅
t ≤ t55d10 → tn10 = t (no density mortality)                             ✅
per-tree rate: rip = (under limit) ? Hamilton ri : rn ; XMMULT window    ✅
VARMRT  distribute (t-tn10) toward suppressed trees                      ✅
QMD-convergence: recompute d10n from survivors, re-iterate (≤10)         ✅
── after convergence ──
MSB alternate mortality (d10>QMDMSB) → MSBMRT                            ⛔ keyword-only (QMDMSB=999 default)
SIZE-CAP mortality (d+g ≥ SIZCAP)                                        ⛔ keyword-only (SIZCAP=999 default)
BAMAX enforcement: scale kills until BA ≤ BAMAX                          ✅ (was the multi-cycle gap; commit aedecd1)
FIXMORT keyword                                                          ⛔ keyword option
TPAMRT = surviving TPA (for next-cycle line reset)                       🟡
```

The two ⛔ here are keyword-gated (default values make them no-ops), so they are
**correctly absent for the current scenarios** but must be ported before the
HTGSTOP/SIZECAP/MORTMULT/FIXMORT keywords are claimed. The 🟡 line-reset is a
known simplification valid for closed stands (no ingrowth); revisit when
regeneration lands.

---

## 4. What is *not* traced, and why

Grouped by the reason, so "not ported" is a deliberate, visible decision rather
than something discovered by accident mid-debug:

1. **Insect/disease models** (mountain pine beetle, DF beetle, tussock moth,
   western budworm, mistletoe): `MPB*/DFB*/TM*/BWE*/MISTOE`. 🧊 No test scenario
   activates them; not in the C0–C9 plan.
2. **FFE fire** (`FMSDIT`, `FMMAIN`, fire reports): ⛔ **C7**. This is the only
   reason the `s10_fire` / `fire_*` scenarios diverge (no fire-caused mortality).
3. **Regeneration / establishment** (`REGENT` establishment mode, `ESFLTR`,
   `CLAUESTB`, `ESNUTR`, `BRTREG`, `RDTREG`, `ESOUT`): ⛔ **C4** remaining. Not
   exercised by the snt/all/mix scenarios (NOAUTOES).
4. **Economics** (`ECSETP/ECSTATUS`, `eccalc`): 🧊 **C8** (a minimal ANNUCST
   path exists for the econ test only).
5. **Keyword option paths** not yet wired: `FIXDG`, `FIXHTG`, `FFERT`, `HTGSTP`
   (HTGSTOP/TOPKILL), `MORTMULT`, `FIXMORT`, MSB/SIZECAP mortality. ⛔ All are
   no-ops at default settings, so they do not affect current tests.
6. **Canopy cover & misc reports** (`CVGO/CVBROW/CVCNOP/CVOUT`, `EXTREE`,
   `RDPR/BRPR`, `SVS`): 🧊 output-only extensions.

## 5. How to use this for testing

- Each ⛔/🧊 row is a **scenario to add** when that chunk is taken: a key that
  sets the gating keyword (e.g. SIZECAP) turns a no-op branch into a live one and
  gives a Fortran-validatable target.
- Each ✅ row should have at least one scenario in the C10 matrix that makes its
  branch *do something distinguishable*; §3-style branch checklists are the way
  to confirm a ported routine is branch-complete, not just present.
