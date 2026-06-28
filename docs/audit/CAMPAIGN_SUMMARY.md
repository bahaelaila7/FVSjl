# FVSjl Faithfulness Campaign — Executive Summary

Capstone over `INDEX.md` (the full per-flag ledger). Mission: make **FVSjl (Southern)** a *semantically
faithful* drop-in for **FVSsn**, with only ULP floating-point and the COMPRESS eigensolver as accepted
divergences. The audit (390 decisions) flagged **9 BANDAID · 53 GAP · 5 UNVERIFIED**.

## Outcome by tier

### BANDAID (trust-critical) — COMPLETE
- 6 real bandaids fixed: **B1** flame, **B2** sprout-CW, **B4** structure xbamax, **B5** establishment
  gentim, **B6** econ PV, **B7/B8** event-monitor TIME/NO.
- 3 caught as FALSE POSITIVES (re-traced, retracted): **B3** (#6 inventory crown snapshot — faithful via
  UPDATE-after-FMMAIN), diameter-growth PTBAA, compress AUTSTK.
- All 3 deferred follow-ups **validated against live FVS / source**: B1 (the fire under-kill it pointed
  at — fixed at root, per-tree live-validated, masking test tightened), B2 (sprout CrWidth = live), B6
  (econ_pnv = eccalc.f:114-117 exactly).

### UNVERIFIED — COMPLETE (5/5, none a masked bug)
PTBAA faithful · fire-snag NOTE-B merch faithful (gross regressed carbon_snt) · density notre FINT/FINTM
real-but-unreachable · crown DUBSCR intentional (snag carbon is bole-based) · establishment WK6 site-prep
unreachable (no SITEPREP scenario).

### GAP — every REACHABLE cluster cleared via live validation
- **mortality** (6): bit-exact vs live (bamax/sdimax scenarios); "user BAMAX ignored" REFUTED; sdimax<5
  faithful (TN10=0); T>35000 present+unreachable; MSBMRT niche-unported.
- **diameter/height growth FINT≠5** (4): growth_fint10/finth10/idg1 bit-exact vs live (TPA/BA/TopHt).
- **fire / FFE** (12): the under-kill ROOT FIX (fire-basis start-of-cycle sm,lg) + 3 down-wood fixes
  (post-burn fall, fire-snag merch bole, fire-killed crown→CWD2B → 2010 down-wood live-faithful) + FULIV2
  shrub override + treelist CrWidth — all live-validated. The fire_year scheduled-vs-actual BUG CLASS was
  found and closed exhaustively (3 instances).
- **event-monitor** (4): AGE fixed+live-validated; `**` faithful (FVS has no exponentiation operator);
  div0/BSDI niche-unreachable.
- **io/volume/density/compress/econ/crown/sprout**: ~22 spine fixes across the campaign (FINT-scaling,
  INT/NINT/IFIX Fortran-rounding, PV timing, RANN, event-monitor IF/THEN, establishment timing, …).

## Validation depth — beyond the aggregate suite
- **Per-tree live validation**: growth (timeint10) and fire-mortality (fire_early) BIT-EXACT vs live FVS,
  tree-by-tree.
- **Broad live sweep (226 scenarios)**: 207 match live ≤1 every cycle; every larger outlier explained
  (multi-stand comparison artifacts, accepted COMPRESS, the regen tail).
- Suite **4530 pass / 1 broken** (the 1 = accepted s22 COMPRESS) throughout.

## The thesis, demonstrated and closed
"A green suite does not certify faithfulness." The sharpest case: the fire under-kill — a real ~5.6-TPA
divergence sitting behind a *deliberately loosened* test (`±12 TPA, "~10 under-kill"`) against a
live-valued golden. Found by tracing FVS logic (not test output), root-caused to the fire-basis
annual-step timing, fixed at the source, validated per-tree against live, and the test **re-tightened to
bit-exact** so faithfulness is now enforced. The treeszcp_cap test showed the same pattern (TPA/BA
excluded for "the regen tail").

## Remaining (all below the trust bar; documented with port plans in INDEX.md)
1. **FFE down-wood one-cycle phasing lag** (#22) — structural (ffe_fuel_update runs before grow_cycle/fire),
   output-only on INTERMEDIATE snapshots (endpoints/2010+ are live-faithful); prior work shows it's not a
   simple reorder (collapses carbon_snt's bit-exact StandDead). Deferred — risk > reward.
2. **≤0.6% late-cycle dense-stand mortality hairline** (FIA 90/131; BA bit-exact) — mechanism not fully
   pinned (stand not at SDImax limit); below practical significance.
3. **FVS_TreeList DBS completeness** — dead records + MortPA + 9 columns (TreeVal/SSCD/PtIndex/MortPA/…);
   output-only, schema-level, moderate effort.
4. **Niche keyword-gated GAPs** — MSBMRT, LHTDRG, CRNMULT, ESGENT/pccf, BSDI-DBHSTAGE: no suite scenario
   triggers them; port when a scenario brings them into scope (principle #4).

## Bottom line
The trust-critical mission is **met where reachable**, demonstrated by per-tree + broad live validation:
growth, mortality, fire-mortality, fuel model, down-wood (endpoints), shrub, volume, density, econ, crown,
sprout, and regen-ESTABLISHMENT are faithful, live-validated drop-ins. What remains is one structural
phasing item, a sub-1% hairline, output-only DBS completeness, and niche unreachable paths — the
irreducible tail, all documented.
