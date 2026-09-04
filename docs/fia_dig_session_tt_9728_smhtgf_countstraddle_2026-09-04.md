# Dig session — TT forward sweep 8025→9728 (population close): 25 escalated structure_densephase (2026-09-04)

Sweep agent, branch `wt/western-sweep`, engine synced to master ab67d239 (src diff vs master EMPTY).
Oracle: `/workspace/.ttwork/FVStt_clean` (ledger default; `FVStt_g16` also relinked).

## Sweep result
TT swept forward **8025 → 9728 = COMPLETE** (cursor at population, one BATCH=2000 cycle, emitted 1703).
**bit_exact 915/1703**; 788 diverging, of which the dig-worthy filter escalated **25** rows (queue never
reached the DIGCAP=100 pause — TT reached population first). Every escalated row is
`signature=structure_densephase`, `worst_col ∈ {CCF, TPA}` (NEVER BA/SDI), `max_rel 15.1–87.5%`,
worst cycle 2033–2073.

Ecoregions of the 25: **M331Dm×4, M331Do×2, M331Du×1, M331Dn×4** (11 in the already-cornered M331D
prefix — surfaced only because the escalation guard never auto-drops a structure col ≥15%), plus
**M331Ja×7 and M332Er×7** (14 in NEW ecoregions not covered by the M331D corner).

## MEASURED VERDICT — CORNERED. Same primitive as the M331D smhtgf-ZRAND corner, extended to
## M331Ja + M332Er. The tt small-tree growth ZRAND draw-order straddle → near-carrying-capacity
## self-thinning count-straddle. NOT a new real bug. No engine change; gate stays 339/11.

### Method
Per-cycle full-.sum A/B (all 10 cols) via `scratchpad/dump_tt.jl` (reuses `ledger_fia.jl`
run_live + parse_sum10 + keytext), regime `none`, on 5 representatives spanning both new ecoregions
and all three stand morphologies (pure-seedling, dense-sapling, mature-pole). Cross-checked initial
DBH distribution from FVS_TREEINIT_COND and the worst_col invariant across all 25.

### Evidence
1. **Divergence enters at the GROWTH cycle, not mortality.** In every trace BA/SDI/CCF diverge one
   or more cycles BEFORE TPA moves; the pure-seedling stands keep TPA bit-identical for the whole run.
   - `750156234290487` (M331Ja, meanDBH 0.64", 95% <3"): **TPA bit-identical every cycle**
     (912→792 both sides); BA/SDI/CCF/TopHt/QMD straddle from 2039 (2069 BA 16/21, CCF 16/30,
     TopHt 14/22) — pure small-tree HEIGHT-growth straddle, the M331D seedling shape exactly.
2. **Density (BA/SDI = the SDImax carrying capacity) is PRESERVED / reconverges; only the
   COUNT×SIZE distribution (TPA/QMD) straddles.** worst_col is CCF or TPA for all 25 — never BA/SDI.
   - `387680576489998` (M332Er, dense sapling TPA 4587, QMD 2.2"): TPA count-straddle
     2216/3467 by 2065, but **BA 250/246, SDI 533/527 reconverge** to the same carrying capacity
     (QMD trades 4.6/3.6). jl thins LESS here.
   - `1856534189290487` (M331Ja, dense TPA 6414, QMD 1.3"): growth straddles first (2033 BA 80/83
     while TPA still bit-identical 6357/6357), TopHt bit-identical, then **BA 241/244, SDI 594/572
     reconverge** at 2073 (TPA 3918/2804, QMD 3.4/4.0). jl thins MORE here.
   - `195350383020004` (M331Dn, MATURE pole meanDBH 8.83", only 15% <3"): sub-material cycle-1 seed
     (2022 BA 240/244, QMD & TopHt bit-identical) compounds through self-thinning to TPA 1050/882 at
     2062, **BA 276/276 converged**, QMD 6.9/7.6. Even in a mature stand the density carrying
     capacity is preserved and the straddle is the count distribution — not a large-tree DG bias.
3. **Direction is stand-specific (jl higher in some, lower in others) ⇒ a BALANCED straddle across
   the cluster**, not a systematic coefficient/model bias. TopHt tracks within ±1–2 ft in the
   density-limited stands.
4. Root mechanism is the one exhaustively TreeId+ZRAND-traced in the M331D dig
   (`fia_dig_session_tt_woodland_selfthin_2026-09-04.md`): `tt/smhtgf.f` (and the small-tree diameter
   path) draws the BACHLO/OLDRN ZRAND deviate per record, FVS draws-once-per-original-then-triples
   while jl triples-then-draws-per-record, desyncing the shared RNG stream (#206 class). The
   compounded small growth seed then trips the near-SDImax self-thinning at a slightly different tree
   COUNT (SN `structure_densephase` global-corner mechanism). Both are named FP/RNG primitives and are
   **ecoregion-independent** — the M331D prefix in the registry was simply narrower than the primitive.

## Corner (durable record on the sweep branch — to be synced to master with the M331D corner)
- Added the 25 CNs to `docs/fia_cornered_stands.txt` (filter drops cornered CNs before the escalation
  guard — the effective mechanism; the cluster prefix alone cannot, because the guard never
  auto-drops a ≥15% structure col).
- Added `M331Ja` and `M332Er` (structure_densephase) rows to `docs/fia_cornered_clusters.tsv`.
- NOTE for the fix-coordinator: `filter_digworthy.jl` reads `/workspace/FVSjl/docs/fia_cornered_stands.txt`
  (main repo). This corner must land in master/main-repo (as the M331D 32 did) to drop these on any re-sweep.

No engine fix. Aligning the small-tree ZRAND draw sequence is a high-regression change to the RNG
path shared by every bit-exact-or-cornered TT stand, with no measured net improvement (balanced
straddle, density preserved) — cornered, not fixed.
