# BC NEWSPRED (spatial dwarf-mistletoe / NISI) port plan — #196

USER-greenlit 2026-08-13 ("scope then start porting"). Closes the BC YSM multi-cycle
under-mortalization (jl SDI→1586 vs oracle 926) attributed to the missing spatial DM.

## Model overview
`canada/newmist` = the **NISI** (New & Improved Spread & Intensification) *spatial*
dwarf-mistletoe model. 8,761 lines / 50 `.f` files. Distinct from the base `mistoe.f`
(non-spatial, already ported for N-Rockies via `_ie_mis_variant` + CR's own
`dwarf_mistletoe_model.jl`). BC currently has NO DM model wired.

**What it does:** per crown-third of each tree record it maintains 4 life-history
compartments (IMMATURE / LATENT / NON-FLOWERING / FLOWERING). Infections spread
spatially between neighbouring trees (between-tree distance × crown shape × slope ×
spatial autocorrelation), intensify over cycles, and produce a per-tree DMR (dwarf-
mistletoe rating 0-6). DMR drives HEIGHT + DIAMETER growth multipliers (BHTG/DHTG,
mistoe.f) and extra MORTALITY (MISMRT) — the latter is the YSM gap.

## Call graph (entry → payoff)
```
MISIN0 → DMINIT (dminitbc.f)                      # one-time init of DMCOM commons
FVS cycle → MISTOE (mistoe.f)                     # per-cycle DM driver
   ├─ DMTREG (dmtreg.f)                           # NISI spread & intensification core
   │    ├─ DMOPTS (keywords DMCLMP/DMALPH/DMBETA) → DMNAC
   │    ├─ DMMTRX / DMFBRK / DMFSHD / DMADLV      # spread matrix, bark, shade, adult-lives
   │    ├─ DMNB / DMNTRD / DMNDMR / DMSRC         # neighbour build, target reduction, DMR, source
   │    ├─ DMSLST / DMSLOP / DMSAMP / DMOTHR / DMTLST
   │    └─ BNDIST (bndist.f) + DMSHAP + DMCW*     # spatial substrate (distance/shape/crown-width)
   ├─ DMMDMR (dmmdmr.f)  → per-tree DMR
   ├─ DMCYCL (dmcycl.f)  → life-history compartment advance (immature→latent→…→flowering)
   ├─ DMAUTO (dmauto.f)  → SF spatial-autocorrelation reweighting matrix
   ├─ MISINF (infection) + MISMRT (mortality, USEMRT)   # ← PAYOFF: extra mortality
   └─ HTG/DG multipliers (BHTG/DHTG in mistoe.f)        # ← PAYOFF: growth reduction
Reports: DMSUM / DMTLST / DMSLST / misprt.f
```
Keywords (currently unrecognized in jl): MISTOE, NEWSPRED, DMAUTO, DMCLMP, DMALPH, DMBETA.

## Chunk decomposition (dependency-ordered)
- **C0 — Foundation + validation harness.** Port `DMCOM.F77` common → a Julia
  `MistletoeState` (per-crown-third 4-compartment arrays, per-tree DMR, SF matrix,
  RNG seed `DMRNSD`). Build the VALIDATION VEHICLE (see below) so every later chunk
  validates vs a live oracle. NO model logic yet — just state + harness + a loud
  "NEWSPRED active but unported" error so partial states never silently mis-run
  (doctrine #5).
- **C1 — Keywords + DMINIT.** Parse MISTOE/NEWSPRED/DMAUTO/DMCLMP/DMALPH/DMBETA
  (`dmopts.f`, `dmauto.f` keyword halves) + wire `DMINIT` (dminitbc.f) + `DMRNSD`
  (dmrann.f RNG seed). Validate: keyword echo + init state vs FVSbc DEBUG.
- **C2 — Coefficients + BC crown width.** Port `dmblkd.f` (256-line DATA) + the BC
  crown-width routine (dmcw.f dispatch → the BC-applicable dmcw*). Extract to CSV,
  verify line-for-line.
- **C3 — Spatial substrate.** `bndist.f` (between-tree distance, 342), `dmshap.f`
  (crown shape, 378), `dmslop.f` (slope). Validate per-tree distances/shapes vs g16.
- **C4 — NISI spread & intensification core.** `dmtreg.f` (553) + its helpers
  (dmmtrx/dmfbrk/dmfshd/dmadlv/dmnb/dmntrd/dmndmr/dmsrc/dmslst/dmsamp/dmothr/dmtlst).
  Largest chunk; likely split into C4a (matrix+neighbour build) / C4b (spread+DMR).
- **C5 — Life-history + autocorrelation.** `dmcycl.f` (593, compartment advance) +
  `dmauto.f` (SF reweighting). Validate compartment counts per crown-third vs g16.
- **C6 — PAYOFF: stand coupling.** `mistoe.f` driver + `dmmdmr.f` (DMR) + `misinf`
  (infection) + `mismrt` (MISMRT mortality, USEMRT) + BHTG/DHTG growth multipliers.
  Wire into the BC growth+mortality dispatch. THIS closes the YSM gap. Validate the
  YSM .sum TPA/SDI vs oracle.
- **C7 — Reports.** `dmsum.f` / `dmtlst.f` / `dmslst.f` / `misprt.f` DM summary +
  treelist reports (validation-supporting, low-risk).

## Validation strategy — unblocking the crash-blocked A/B
The clean A/B is blocked because the gfortran-16-relinked `FVSbc_clean` SIGSEGVs on
ALL DATABASE reads (isoc23-shim × SQLite C-interop, dbstreesin.f) — NOT a production
bug. Three vehicles, in priority order:
1. **INLINE TREEDATA reproducer (primary).** Extract YSM's trees from the FVS-ready
   DB into an inline `.tre` + a BC keyfile (STDINFO BEC string) with NEWSPRED/MISTOE
   enabled. The INLINE path in FVSbc_clean WORKS (all_BC/all_BC_essf run inline). This
   yields the oracle `.sum` WITH DM mortality → the C6 payoff target. Build in C0.
2. **FVSbc_g16 instrument-replay (per-chunk).** Full gfortran-16 rebuild of
   `bin/FVSbc_buildDir/*.f` (+ newmist) with the isoc23 shim = the durable
   instrumentable oracle (same recipe as FVS{ut,tt,em}_g16). Add DEBUG/WRITE dumps to
   dmtreg/dmcycl/bndist for per-chunk per-tree validation (DMR, distances, compartments).
   Confirm it reads inline stands (avoids the DB path entirely).
3. **DB-capable relink (fallback).** If an inline reproducer can't capture the YSM
   spatial config, relink FVSbc against a non-shim SQLite to unblock the DB path — only
   if 1+2 prove insufficient.

## Integration notes
- Existing jl base DM: `centralrockies/dwarf_mistletoe_model.jl` (CR) +
  `inlandempire/mistoe_coefficients.jl` (`_ie_mis_variant`, N-Rockies). NEWSPRED is a
  SEPARATE spatial model — do NOT overload the base path; add a `bc_newspred!` dispatch
  gated on the NEWSPRED keyword, wired into the BC cycle only (BC-first; other variants
  have their own dmcw* but are out of scope for #196).
- Doctrine: validate bit-exact-or-cornered per chunk vs the inline oracle / g16; MEASURE
  don't infer; never FFI the RNG (DMRNSD/DMRANN is its own stream — expect a realization
  straddle on the spatial draws, same accepted class as ZZRAN).
- Risk: the spatial neighbour model uses per-tree (x,y) positions FVS synthesizes from
  point/plot layout — confirm jl reproduces the same synthetic positions (bndist input)
  before trusting spread counts; this is the likely first "cornered vs real" fork.

## Status
- **C0 validation foundation DONE 2026-08-13.** Key simplification: jl reads the YSM
  metric DB directly (SQLite.jl + the 42eb555 metric-ingest fix) — NO inline reproducer
  needed for the PAYOFF validation. Set up `/workspace/.bcwork/newspred/ysm271.key`
  (YSM029-271 from FVS-BC.YSM-SkyRanch.db, NOTRIPLE, TIMEINT 5/first-4, NUMCYCLE 24).
  jl(DM-free) vs oracle `YSM-SkyRanch.sum.save`:
    - **cyc0 (2018) BIT-EXACT** (jl 2300/10/311 vs oracle 2300/9/311) → stand loads +
      growth correct; the DM effect is the ENTIRE divergence.
    - **C6 target quantified** — NEWSPRED must kill trees AND suppress growth: by 2077
      jl 1713/79/1586 (TPA/BA/SDI) → oracle 1366/42/926 (−20% TPA, −47% BA, −42% SDI).
  Payoff-validation vehicle = jl DB-run vs `.sum.save` (works today). Inline FVSbc
  reproducer for per-chunk internals DEFERRED to C4 (only chunk likely to need it).
- **NEXT: C1** — DMCOM common → `MistletoeState` struct + parse MISTOE/NEWSPRED/DMAUTO/
  DMCLMP/DMALPH/DMBETA/MISTPRT keywords (currently in unrecognized_keywords) + wire
  DMINIT/DMRNSD. Validate: keyword echo + init state.
- Multi-session; each chunk lands + validates before the next. Off-switch untouched (USER's).
